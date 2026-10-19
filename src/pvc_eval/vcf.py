"""The VCF records of one file: stage 2 of doc/spec.md §4, read with the checks of §5.2
and normalised, trimmed then left-aligned. A VCF record is never split or decomposed.

Positions are 0-based and half-open inside pvc_eval. The 1-based POS of a VCF line is
converted here, and nowhere else.
"""

import gzip
import logging
import math
import re
from collections.abc import Generator, Iterator, Sequence
from contextlib import closing
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError

logger = logging.getLogger(__name__)

PLOIDY = 2  # doc/spec.md §2: a global constant

MISSING = sentinel("MISSING")
"""A haplotype the VCF genotype says nothing about, written `.` (doc/spec.md §2)."""

MALFORMED_GT = sentinel("MALFORMED_GT")
"""A GT with a number of alleles other than the ploidy, or an allele index the VCF
record lacks: the sample is invalid at the locus, and the run goes on (§5.4)."""

SPANNING_DELETION = "*"
GZIP_MAGIC = b"\x1f\x8b"
DECIMAL_NUMBER = re.compile(r"[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?")
NUM_FIXED_COLUMNS = 9  # CHROM to FORMAT


@dataclass(frozen=True, slots=True, kw_only=True)
class VcfGenotype:
    """One sample's GT at one VCF record, as written.

    alleles: one entry per haplotype, the index of a VCF allele (0 is REF), or MISSING.
    phased: the GT was written with `|`.
    phase_set: the PS of the sample, or None when the record gives none, which VCF
        reads as one phase set per sample and contig (§5.4).
    quality: the GQ of the sample, or None where absent.
    """

    alleles: tuple[int | MISSING, ...]
    phased: bool
    phase_set: str | None
    quality: int | None


@dataclass(frozen=True, slots=True, kw_only=True)
class VcfRecord:
    """One line of a VCF file.

    start: the first position REF covers, 0-based; `end` is one past the last.
    alleles: REF first, then the ALTs, as written but in upper case; `*` is a spanning
        deletion.
    qual: QUAL, or None for `.`.
    genotypes: one per sample, in the order of the header; a MALFORMED_GT keeps no GQ.
    line_number: the 1-based line of the file, for messages.
    ref_mismatch: REF differs from the reference (§5.2), set by normalising.
    """

    contig: str
    start: int
    alleles: tuple[str, ...]
    qual: float | None
    genotypes: tuple[VcfGenotype | MALFORMED_GT, ...]
    line_number: int
    ref_mismatch: bool = False

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be 0 or more, not {self.start}")
        if not self.alleles or not self.alleles[0]:
            raise ValueError("a VCF record needs a REF of at least one base")

    @property
    def end(self) -> int:
        return self.start + len(self.alleles[0])


def read_samples(path: Path) -> tuple[str, ...]:
    """The sample names of a VCF file, in the order of its columns.

    Raises PvcEvalInputError when the file has no #CHROM line after its ## lines.
    """
    with closing(_numbered_lines(path)) as numbered_lines:
        return _read_header(numbered_lines, path)


def read_vcf_records(path: Path, contigs: Sequence[str]) -> Iterator[VcfRecord]:
    """The VCF records of a file, in the file's order, plain or bgzipped.

    contigs: the contig names of the reference, in the order of its `.fai`.

    Raises PvcEvalInputError, when the stream reaches it, for a VCF record on a contig
    the reference lacks, for records out of the `.fai`'s contig order or out of position
    order, for a file with no #CHROM line after its ## lines, for a line that does not
    hold a VCF record, and for one whose number of sample columns differs from the
    number of samples the #CHROM line names. A blank line is skipped. A record with a symbolic
    allele is skipped, and the number skipped is logged as a warning at the end (§9). A
    QUAL or GQ that is not a number is taken as absent, and those are counted in a
    warning at the end too.
    """
    contig_rank = {contig: rank for rank, contig in enumerate(contigs)}
    # The (contig rank, POS) of the previous record, which this one must not come before,
    # and where that record is, for the message.
    previous_order: tuple[int, int] | None = None
    previous_where = ""
    num_symbolic = 0
    num_not_numbers = 0
    with closing(_numbered_lines(path)) as numbered_lines:
        num_samples = len(_read_header(numbered_lines, path))
        for line_number, line in numbered_lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            contig, pos_1based = _read_position(fields, path, line_number)
            _check_num_samples(fields, num_samples, path, line_number)
            if contig not in contig_rank:
                raise PvcEvalInputError(
                    f"{path}, line {line_number}: contig {contig} is not in the reference."
                    " Compare against the reference the VCF was called on."
                )
            order = (contig_rank[contig], pos_1based)
            where = f"{contig}:{pos_1based} (line {line_number})"
            if previous_order is not None and order < previous_order:
                raise PvcEvalInputError(
                    f"{path} is not sorted: {where} comes after {previous_where}."
                    " Sort it in the contig order of the reference's .fai and by position,"
                    " and run again."
                )
            previous_order, previous_where = order, where
            alts = _read_alts(fields, path, line_number)
            if any(_is_symbolic(alt) for alt in alts):
                num_symbolic += 1
                continue
            record, num_record_not_numbers = _parse_record(
                fields, contig, pos_1based - 1, alts, line_number
            )
            num_not_numbers += num_record_not_numbers
            yield record
    if num_symbolic:
        logger.warning(
            "%s: %d VCF records with symbolic alleles (<DEL>, breakends) were skipped;"
            " pvc_eval does not compare structural variants.",
            path,
            num_symbolic,
        )
    if num_not_numbers:
        logger.warning(
            "%s: %d QUAL or GQ values are not numbers and were taken as absent, as `.` is;"
            " the quality of the rows they would give is `.`.",
            path,
            num_not_numbers,
        )


def parse_genotype(
    gt: str, num_alleles: int, phase_set: str | None, quality: int | None
) -> VcfGenotype | MALFORMED_GT:
    """A GT field as a VcfGenotype, or MALFORMED_GT when it cannot be one (§5.4).

    num_alleles: how many VCF alleles the record has, REF included.
    phase_set, quality: the sample's PS and GQ, kept in the VcfGenotype.
    """
    parts = gt.replace("|", "/").split("/")
    if len(parts) != PLOIDY:
        return MALFORMED_GT
    alleles: list[int | MISSING] = []
    for part in parts:
        if part == ".":
            alleles.append(MISSING)
        elif part.isascii() and part.isdecimal() and int(part) < num_alleles:
            alleles.append(int(part))
        else:
            return MALFORMED_GT
    return VcfGenotype(
        alleles=tuple(alleles), phased="|" in gt, phase_set=phase_set, quality=quality
    )


def _numbered_lines(path: Path) -> Generator[tuple[int, str]]:
    """The lines of a VCF file, plain or bgzipped, each with its 1-based number.

    Raises PvcEvalInputError for a file that is not UTF-8 text, a BCF among them, and
    for a gzip file that is cut short or damaged.
    """
    # Text is decoded and decompressed in blocks, so a fault is found somewhere after the
    # last line read, and only that is known.
    line_number = 0
    try:
        with _open_text(path) as lines:
            for line_number, line in enumerate(lines, start=1):
                yield line_number, line
    except UnicodeDecodeError as err:
        raise PvcEvalInputError(
            f"{path} is not UTF-8 text: a byte at or after line {line_number + 1} cannot be read"
            " as UTF-8. If it is a BCF file, convert it to a VCF with"
            " `bcftools view -Oz -o file.vcf.gz file.bcf`; if it is a VCF in another"
            " encoding, convert it to UTF-8; and run again."
        ) from err
    except (EOFError, gzip.BadGzipFile) as err:
        raise PvcEvalInputError(
            f"{path} is a gzip file cut short or damaged at or after line {line_number + 1}."
            " Compress or copy it again, and run again."
        ) from err


def _open_text(path: Path) -> TextIO:
    with path.open("rb") as raw:
        is_gzipped = raw.read(len(GZIP_MAGIC)) == GZIP_MAGIC
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _read_header(numbered_lines: Iterator[tuple[int, str]], path: Path) -> tuple[str, ...]:
    """The sample names of the #CHROM line, read up to it and no further."""
    for line_number, line in numbered_lines:
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#CHROM"):
            return tuple(line.rstrip("\n").split("\t")[NUM_FIXED_COLUMNS:])
        raise PvcEvalInputError(
            f"{path}, line {line_number}: a VCF file has a #CHROM line after its ## lines,"
            " and this line comes before any. Check that the file is a VCF with its header."
        )
    raise PvcEvalInputError(
        f"{path} has no #CHROM line after its ## lines. Check that the file is a VCF with"
        " its header."
    )


def _check_num_samples(fields: list[str], num_samples: int, path: Path, line_number: int) -> None:
    """Later stages find a sample's genotype by its column, so every line has one per
    sample of the #CHROM line, after FORMAT. A sites-only file has neither."""
    num_found = max(len(fields) - NUM_FIXED_COLUMNS, 0)
    if num_found != num_samples:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has {num_found} sample columns, and the"
            f" #CHROM line names {num_samples} samples. A VCF record has a FORMAT column and"
            " then one column per sample; correct the line or remove it, and run again."
        )


# A line that is not a VCF record stops the run (a provisional choice: doc/spec.md §3 does
# not list it), with the column at fault.
def _read_position(fields: list[str], path: Path, line_number: int) -> tuple[str, int]:
    if len(fields) < NUM_FIXED_COLUMNS - 1:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has {len(fields)} tab-separated columns,"
            " and a VCF record has at least 8, CHROM to INFO. Check that its columns are"
            " separated by tabs, correct the line or remove it, and run again."
        )
    pos = fields[1]
    if not (pos.isascii() and pos.isdecimal() and int(pos) >= 1):
        raise PvcEvalInputError(
            f"{path}, line {line_number}: POS is '{pos}', and it must be a whole number"
            " of 1 or more. Correct the line or remove it, and run again."
        )
    return fields[0], int(pos)


def _read_alts(fields: list[str], path: Path, line_number: int) -> tuple[str, ...]:
    """The ALTs of a line, () for `.`, after checking that REF and every ALT hold
    something."""
    if not fields[3]:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: REF is empty, and it must hold at least one base."
            " Correct the line or remove it, and run again."
        )
    if fields[4] == ".":
        return ()
    alts = tuple(fields[4].split(","))
    if not all(alts):
        raise PvcEvalInputError(
            f"{path}, line {line_number}: ALT '{fields[4]}' has an empty allele; each ALT"
            " must hold at least one base, and a record with no ALT has `.`. Correct the"
            " line or remove it, and run again."
        )
    return alts


def _is_symbolic(alt: str) -> bool:
    """A symbolic ALT (`<DEL>`), a breakend (`A]chr2:20]`) or a single breakend of VCF 4.3
    (`A.`, `.A`); `.` alone is no ALT at all."""
    is_single_breakend = len(alt) > 1 and (alt.startswith(".") or alt.endswith("."))
    return alt.startswith("<") or "[" in alt or "]" in alt or is_single_breakend


def _parse_record(
    fields: list[str], contig: str, start: int, alts: tuple[str, ...], line_number: int
) -> tuple[VcfRecord, int]:
    """The VCF record of a line, and how many of its QUAL and GQ values were taken as
    absent for not being numbers."""
    # VCF bases are case-insensitive, and the reference is compared in upper case.
    alleles = (fields[3].upper(), *(alt.upper() for alt in alts))
    genotypes: list[VcfGenotype | MALFORMED_GT] = []
    qual = _read_qual(fields[5])
    num_not_numbers = _is_not_a_number(fields[5], qual)
    if len(fields) > NUM_FIXED_COLUMNS:
        keys = fields[8].split(":")
        for sample_field in fields[NUM_FIXED_COLUMNS:]:
            # VCF lets trailing sample fields be dropped, so a key may have no value.
            values = dict(zip(keys, sample_field.split(":"), strict=False))
            phase_set = values.get("PS")
            if phase_set == ".":
                phase_set = None
            quality = _read_gq(values.get("GQ"))
            num_not_numbers += _is_not_a_number(values.get("GQ"), quality)
            if "GT" not in values:
                genotypes.append(
                    VcfGenotype(
                        alleles=(MISSING,) * PLOIDY,
                        phased=False,
                        phase_set=phase_set,
                        quality=quality,
                    )
                )
                continue
            genotypes.append(parse_genotype(values["GT"], len(alleles), phase_set, quality))
    record = VcfRecord(
        contig=contig,
        start=start,
        alleles=alleles,
        qual=qual,
        genotypes=tuple(genotypes),
        line_number=line_number,
    )
    return record, num_not_numbers


# QUAL and GQ only annotate a row of the table (§6), where an absent value is written `.`;
# a value that is not a number is taken as absent and counted, rather than stopping the
# run (a provisional choice). Only decimal numbers count, not Python's "nan", "inf" or
# "1_000", and a QUAL must be finite.
def _read_qual(field: str) -> float | None:
    if not DECIMAL_NUMBER.fullmatch(field):
        return None
    qual = float(field)
    return qual if math.isfinite(qual) else None


def _read_gq(field: str | None) -> int | None:
    if field is None or not field.isascii() or not field.isdecimal():
        return None
    return int(field)


def _is_not_a_number(field: str | None, number: float | None) -> bool:
    """A QUAL or GQ was written, and not as `.`, but taken as absent."""
    return number is None and field not in (None, ".")


def normalise(record: VcfRecord, contig_sequence: str) -> VcfRecord:
    """The record trimmed and, when it is an indel, left-aligned.

    contig_sequence: the record's contig in the reference, in upper case, as the
        record's alleles are.

    A record whose REF differs from the reference is returned unchanged and marked
    ref_mismatch, since its alleles cannot be placed on the reference (§5.2).
    """
    if contig_sequence[record.start : record.end] != record.alleles[0]:
        return replace(record, ref_mismatch=True)
    return left_align(trim(record), contig_sequence)


def trim(record: VcfRecord) -> VcfRecord:
    """The record without the bases all its VCF alleles share, at the right end, then
    at the left, keeping at least one base in each (the VCF anchor convention).

    A `*` allele has no bases and takes no part. A record with fewer than two other
    alleles is returned as it is.
    """
    sequences = _sequences(record.alleles)
    if len(sequences) < 2:
        return record
    start, trimmed = _trim(record.start, sequences)
    return replace(record, start=start, alleles=_put_back(record.alleles, trimmed))


def left_align(record: VcfRecord, contig_sequence: str) -> VcfRecord:
    """The record with its indel slid left while the base it uncovers equals the base it
    drops, as `bcftools norm` and vt do, and written with one shared base on its left.

    An ALT is an indel when, without the bases it shares with REF at both ends, one of
    the two is empty, whichever side the shared bases were on. A record moves only when
    every ALT is an indel and every ALT ends at the same place; otherwise it stays where
    it is: a SNP or a complex ALT cannot move, and ALTs that would move apart cannot be
    written as one record (§5.2). An indel at the first base of the contig, with no base
    left of it, stays too.
    """
    sequences = _sequences(record.alleles)
    if len(sequences) < 2:
        return record
    ref = sequences[0]
    alts = sequences[1:]
    indels = [
        indel
        for alt in alts
        if (indel := _leftmost_indel(record.start, ref, alt, contig_sequence)) is not None
    ]
    positions = {indel.position for indel in indels}
    if len(indels) != len(alts) or len(positions) != 1:
        return record
    position = positions.pop()
    if position == 0:
        return record
    # The anchor base is at position - 1, and REF covers the longest deletion.
    end = position + max(indel.num_deleted for indel in indels)
    anchor = contig_sequence[position - 1]
    moved_alts = [
        anchor + indel.inserted + contig_sequence[position + indel.num_deleted : end]
        for indel in indels
    ]
    alleles = [contig_sequence[position - 1 : end], *moved_alts]
    return replace(record, start=position - 1, alleles=_put_back(record.alleles, alleles))


@dataclass(frozen=True, slots=True)
class _Indel:
    """Bases inserted before, or deleted from, a 0-based position of the contig."""

    position: int
    inserted: str
    num_deleted: int


def _leftmost_indel(start: int, ref: str, alt: str, contig_sequence: str) -> _Indel | None:
    """The indel ref -> alt at its leftmost position that keeps a base on its left, or at
    0 when it is at the first base and cannot move; None for a SNP or a complex ALT.

    The shared bases are removed from the left first, so that the indel starts as far
    right as the record allows, and sliding left reaches its leftmost place from there.
    """
    num_prefix = _num_shared(ref, alt)
    ref, alt = ref[num_prefix:], alt[num_prefix:]
    num_suffix = _num_shared(ref[::-1], alt[::-1])
    ref, alt = ref[: len(ref) - num_suffix], alt[: len(alt) - num_suffix]
    if bool(ref) == bool(alt):
        return None
    event = ref or alt  # the inserted or deleted bases
    position = start + num_prefix
    while position > 1 and contig_sequence[position - 1] == event[-1]:
        event = contig_sequence[position - 1] + event[:-1]
        position -= 1
    return _Indel(
        position=position, inserted=event if alt else "", num_deleted=len(event) if ref else 0
    )


def _num_shared(first: str, second: str) -> int:
    """How many bases two sequences share at their start."""
    for index, (base, other_base) in enumerate(zip(first, second, strict=False)):
        if base != other_base:
            return index
    return min(len(first), len(second))


def _trim(start: int, alleles: Sequence[str]) -> tuple[int, list[str]]:
    trimmed = list(alleles)
    while min(map(len, trimmed)) > 1 and len({allele[-1] for allele in trimmed}) == 1:
        trimmed = [allele[:-1] for allele in trimmed]
    while min(map(len, trimmed)) > 1 and len({allele[0] for allele in trimmed}) == 1:
        trimmed = [allele[1:] for allele in trimmed]
        start += 1
    return start, trimmed


def _sequences(alleles: Sequence[str]) -> list[str]:
    """The alleles that are sequences. A `*` has no bases and takes no part in trimming or
    left-aligning; _put_back returns it to its place."""
    return [allele for allele in alleles if allele != SPANNING_DELETION]


def _put_back(alleles: tuple[str, ...], sequences: list[str]) -> tuple[str, ...]:
    """The alleles with the sequences _sequences took out replaced, in order, by these."""
    replacements = iter(sequences)
    return tuple(
        allele if allele == SPANNING_DELETION else next(replacements) for allele in alleles
    )
