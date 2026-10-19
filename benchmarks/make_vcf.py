"""Write a bgzipped VCF over a region of a reference, the input of the benchmarks of
parse_vcf.

Every REF is taken from the reference, so the file has no REF mismatch. The VCF records
are SNPs, about one every 300 bases, and indels of 1 to 20 bases, one for every nine
SNPs. Half the indels are placed in a repeat run found in the sequence, a run of a unit of
1 to 4 bases repeated at least three times, and a third of those are written at the
right end of the run instead of the left. About 5% of the records have two ALTs, and a
few SNPs inside a deletion carry a `*` ALT. Each sample has a GT, mostly 0/0 and 0/1,
some 1/1 and some ./., a quarter of them phased with a PS, and a GQ. Positions with an N
in the reference get no record.

The same seed gives the same file, byte for byte.

    uv run python benchmarks/make_vcf.py FASTA chr20:1-5000000 --samples 200 \\
        --output benchmarks/data/chr20_5mb_200s.vcf.gz
"""

import argparse
import random
import re
import struct
import zlib
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

from pvc_eval.fasta import IndexedFasta

_MEAN_GAP = 270  # bases between records: SNPs every 300, with one indel per nine SNPs
_INDEL_SHARE = 0.1
_REPEAT_SHARE = 0.5  # of the indels
_RIGHT_ANCHORED_SHARE = 1 / 3  # of the indels in a repeat run
_MULTI_ALT_SHARE = 0.05
_STAR_SHARE = 0.02  # of the deletions of 2 bases or more, which get a SNP inside them
_MAX_INDEL = 20
_REPEAT_SEARCH = 200  # bases after a position in which a repeat run is looked for
_MISSING_SHARE = 0.03  # of the GTs, written ./.
_PHASED_SHARE = 0.25  # of the GTs that are not ./.
_PHASE_BLOCK = 100_000  # bases; a phased GT's PS is the first position of its block
_BASES = "ACGT"
_REPEAT_RUN = re.compile(r"([ACGT]{1,4}?)\1{2,}")
_GQ_VALUES = tuple(range(100))
_GQ_WEIGHTS = tuple(1 + gq * gq for gq in _GQ_VALUES)  # most GQs high, as callers give
_BGZF_BLOCK = 65_280  # the uncompressed bytes of a block, as bgzip writes them
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


@dataclass(frozen=True, slots=True)
class _Event:
    """One VCF record before its sample columns: start is 0-based."""

    start: int
    ref: str
    alts: tuple[str, ...]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("fasta", type=Path, help="the reference, with its .fai beside it")
    parser.add_argument("region", help="contig or contig:start-end, 1-based and inclusive")
    parser.add_argument("--samples", type=int, default=200, help="number of samples")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    contig, start, end = _parse_region(args.region)
    rng = random.Random(args.seed)
    with IndexedFasta(args.fasta) as reference:
        contig_length = reference.get_length(contig)
        end = contig_length if end is None else end
        sequence = reference.get_bases(contig, 0, end)
    samples = [f"sample{index + 1:03d}" for index in range(args.samples)]
    header = _make_header(args.fasta, contig, contig_length, samples)
    events = _make_events(sequence, start, rng)
    lines = (_make_line(contig, event, len(samples), rng) for event in events)
    with args.output.open("wb") as out:
        _write_bgzf(out, _join_lines(header, lines))


def _parse_region(region: str) -> tuple[str, int, int | None]:
    """The contig, the 0-based start and the end of a region, the end None for the whole
    contig."""
    contig, _, stretch = region.partition(":")
    if not stretch:
        return contig, 0, None
    first, last = stretch.replace(",", "").split("-")
    return contig, int(first) - 1, int(last)


def _make_header(fasta: Path, contig: str, contig_length: int, samples: list[str]) -> list[str]:
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    return [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        f"##reference=file://{fasta.resolve()}",
        f"##contig=<ID={contig},length={contig_length}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "\t".join(columns + samples),
    ]


def _make_events(sequence: str, start: int, rng: random.Random) -> Iterator[_Event]:
    """The records over sequence[start:], in order of position, none overlapping another
    except a `*` record inside its deletion."""
    position = start + 1  # an indel needs the base before it
    while True:
        position += 1 + int(rng.expovariate(1 / _MEAN_GAP))
        if position + _REPEAT_SEARCH + 2 * _MAX_INDEL >= len(sequence):
            return
        if rng.random() < _INDEL_SHARE:
            events = _make_indel(sequence, position, rng)
        else:
            events = _make_snp(sequence, position, rng)
        for event in events:
            yield event
            position = max(position, event.start + len(event.ref))


def _make_snp(sequence: str, position: int, rng: random.Random) -> list[_Event]:
    ref = sequence[position]
    if ref not in _BASES:
        return []
    num_alts = 2 if rng.random() < _MULTI_ALT_SHARE else 1
    alts = rng.sample([base for base in _BASES if base != ref], num_alts)
    return [_Event(position, ref, tuple(alts))]


def _make_indel(sequence: str, position: int, rng: random.Random) -> list[_Event]:
    """An indel at position, or in the first repeat run after it, with now and then a
    second ALT or a `*` record inside a deletion."""
    if rng.random() < _REPEAT_SHARE:
        event = _make_repeat_indel(sequence, position, rng)
    else:
        event = _make_plain_indel(sequence, position, rng)
    if event is None or not _is_bases(event.ref + "".join(event.alts)):
        return []
    if rng.random() < _MULTI_ALT_SHARE:
        event = _add_second_alt(event, rng)
    events = [event]
    is_deletion = len(event.ref) >= 3 and len(event.alts) == 1
    if is_deletion and rng.random() < _STAR_SHARE:
        deleted = event.start + 1
        ref = sequence[deleted]
        snp_alt = rng.choice([base for base in _BASES if base != ref])
        events.append(_Event(deleted, ref, (snp_alt, "*")))
    return events


def _make_plain_indel(sequence: str, position: int, rng: random.Random) -> _Event:
    length = _draw_indel_length(rng)
    anchor = sequence[position - 1]
    if rng.random() < 0.5:
        return _Event(position - 1, sequence[position - 1 : position + length], (anchor,))
    inserted = "".join(rng.choices(_BASES, k=length))
    return _Event(position - 1, anchor, (anchor + inserted,))


def _make_repeat_indel(sequence: str, position: int, rng: random.Random) -> _Event | None:
    """A deletion or an insertion of whole units of the first repeat run after position,
    written at the left end of the run, or at its right end."""
    run = _find_repeat_run(sequence, position)
    if run is None:
        return None
    run_start, run_end, unit_length = run
    num_units = (run_end - run_start) // unit_length
    max_units = _MAX_INDEL // unit_length
    is_deletion = rng.random() < 0.5
    num_changed = rng.randint(1, min(max_units, num_units) if is_deletion else max_units)
    changed = num_changed * unit_length
    is_right_anchored = rng.random() < _RIGHT_ANCHORED_SHARE
    if is_deletion:
        anchor = run_end - changed - 1 if is_right_anchored else run_start - 1
        ref = sequence[anchor : anchor + 1 + changed]
        return _Event(anchor, ref, (ref[0],))
    if is_right_anchored:
        anchor, unit = run_end - 1, sequence[run_end - unit_length : run_end]
    else:
        anchor, unit = run_start - 1, sequence[run_start : run_start + unit_length]
    return _Event(anchor, sequence[anchor], (sequence[anchor] + unit * num_changed,))


def _find_repeat_run(sequence: str, position: int) -> tuple[int, int, int] | None:
    """The start, end and unit length of the first repeat run within _REPEAT_SEARCH bases
    after position, 0-based and half-open."""
    match = _REPEAT_RUN.search(sequence, position, position + _REPEAT_SEARCH)
    if match is None:
        return None
    return match.start(), match.end(), len(match.group(1))


def _add_second_alt(event: _Event, rng: random.Random) -> _Event:
    """The event with an insertion after its anchor as a second ALT."""
    inserted = "".join(rng.choices(_BASES, k=_draw_indel_length(rng)))
    return _Event(event.start, event.ref, (*event.alts, event.ref + inserted))


def _draw_indel_length(rng: random.Random) -> int:
    """1 to _MAX_INDEL bases, short ones the most often."""
    return min(_MAX_INDEL, 1 + int(rng.expovariate(1 / 3)))


def _is_bases(text: str) -> bool:
    return all(base in _BASES for base in set(text))


def _make_line(contig: str, event: _Event, num_samples: int, rng: random.Random) -> str:
    qual = f"{rng.uniform(20, 2000):.1f}"
    fields = [contig, str(event.start + 1), ".", event.ref, ",".join(event.alts), qual, "PASS"]
    fields += [".", "GT:GQ:PS"]
    phase_set = str(event.start // _PHASE_BLOCK * _PHASE_BLOCK + 1)
    fields += _make_sample_columns(len(event.alts) + 1, num_samples, phase_set, rng)
    return "\t".join(fields)


def _make_sample_columns(
    num_alleles: int, num_samples: int, phase_set: str, rng: random.Random
) -> list[str]:
    """One GT:GQ:PS per sample, the GTs drawn in Hardy-Weinberg proportions from allele
    frequencies drawn for the record."""
    frequencies = [rng.betavariate(0.4, 2) for _ in range(num_alleles - 1)]
    frequencies.insert(0, max(0.05, 1 - sum(frequencies)))
    pairs = [(a, b) for a in range(num_alleles) for b in range(a, num_alleles)]
    weights = [frequencies[a] * frequencies[b] * (1 if a == b else 2) for a, b in pairs]
    total = sum(weights)
    weights = [weight / total * (1 - _MISSING_SHARE) for weight in weights]
    gts: list[tuple[int, int] | None] = [*pairs, None]
    drawn = rng.choices(gts, [*weights, _MISSING_SHARE], k=num_samples)
    gqs = rng.choices(_GQ_VALUES, _GQ_WEIGHTS, k=num_samples)
    columns = []
    for gt, gq in zip(drawn, gqs, strict=True):
        if gt is None:
            columns.append("./.:.:.")
        elif rng.random() < _PHASED_SHARE:
            first, second = gt if rng.random() < 0.5 else gt[::-1]
            columns.append(f"{first}|{second}:{gq}:{phase_set}")
        else:
            columns.append(f"{gt[0]}/{gt[1]}:{gq}:.")
    return columns


def _join_lines(header: Iterable[str], lines: Iterable[str]) -> Iterator[bytes]:
    for line in header:
        yield (line + "\n").encode("ascii")
    for line in lines:
        yield (line + "\n").encode("ascii")


def _write_bgzf(out: BinaryIO, chunks: Iterable[bytes]) -> None:
    """Write bytes as BGZF, the blocked gzip of bgzip, ending with its empty block."""
    pending = bytearray()
    for chunk in chunks:
        pending += chunk
        while len(pending) >= _BGZF_BLOCK:
            _write_bgzf_block(out, bytes(pending[:_BGZF_BLOCK]))
            del pending[:_BGZF_BLOCK]
    if pending:
        _write_bgzf_block(out, bytes(pending))
    out.write(_BGZF_EOF)


def _write_bgzf_block(out: BinaryIO, data: bytes) -> None:
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    deflated = compressor.compress(data) + compressor.flush()
    # gzip header with the BC extra field, whose value is the block's size less one.
    header = struct.pack(
        "<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, len(deflated) + 25
    )
    out.write(header + deflated + struct.pack("<II", zlib.crc32(data), len(data)))


if __name__ == "__main__":
    main()
