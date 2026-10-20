"""The re-spelling program of doc/specs/respell.md: a VCF file written again in other words
that give every sample the same haplotypes, for the check of doc/spec.md §7.4.

What every re-spelling shares is here (§2): the input read as lines and as records, the
records every re-spelling keeps, the touched test, the reach, and the writing of the output,
with its header, its order, its new phase sets and, for a gzipped input, BGZF.

Positions are 0-based and half-open, as everywhere inside pvc_eval, but for POS, which is
1-based where a name says so.

Public: respell_by_decomposing (§3), respell_by_shifting (§4),
respell_with_spanning_deletions (§5), respell_by_merging (§6); make_bgzf, the BGZF writer,
which the tests use to write a bgzipped input.
"""

import bisect
import gzip
import itertools
import random
import struct
import sys
import zlib
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import (
    DEFAULT_PLOIDY,
    DEFAULT_SPAN_LIMIT,
    MISSING,
    Span,
    VcfGenotype,
    VcfRecord,
    parse_vcf,
    trim_alt,
)

_GZIP_MAGIC = b"\x1f\x8b"  # the first two bytes of a gzip file, and so of a bgzipped one
_SPANNING_DELETION = "*"
_PLAIN_BASES = frozenset("ACGT")
# The columns of a VCF line, from 0; the samples' start after FORMAT.
_CHROM_COLUMN = 0
_POS_COLUMN = 1
_REF_COLUMN = 3
_ALT_COLUMN = 4
_QUAL_COLUMN = 5
_FORMAT_COLUMN = 8
_NUM_FIXED_COLUMNS = 9
_GT_LINE_START = "##FORMAT=<ID=GT,"
_GT_LINE = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
_PS_LINE_START = "##FORMAT=<ID=PS,"
_PS_LINE = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
# The largest PS bcftools 1.24 keeps in an Integer; it sets 2,147,483,648 to missing (§2.3).
_MAX_INTEGER_PHASE_SET = 2**31 - 1


@dataclass(frozen=True, slots=True, kw_only=True)
class _InputRecord:
    """A record line of the input.

    line: the line as read, without its end.
    place: the line's place among the input's record lines, from 1.
    fields: the line's tab-separated columns.
    vcf_record: what parse_vcf gives for the line, or None for a record with a symbolic
        allele, which parse_vcf skips.
    """

    line: str
    place: int
    fields: tuple[str, ...]
    vcf_record: VcfRecord | None

    @property
    def contig(self) -> str:
        return self.fields[_CHROM_COLUMN]

    @property
    def pos_1based(self) -> int:
        return int(self.fields[_POS_COLUMN])


@dataclass(frozen=True, slots=True, kw_only=True)
class _RespellingInput:
    """The input of a re-spelling, read as lines and as records (§2.1).

    taken_phase_set_numbers: the whole numbers, sorted, whose text is a PS of the input, of
        those a new phase set can reach (_make_new_phase_set).
    are_phase_sets_integer: every PS of the input is a whole number of at most
        _MAX_INTEGER_PHASE_SET, so that the PS line can declare it Integer.
    """

    header_lines: tuple[str, ...]
    input_records: tuple[_InputRecord, ...]
    is_gzipped: bool
    taken_phase_set_numbers: tuple[int, ...]
    are_phase_sets_integer: bool


@dataclass(frozen=True, slots=True, kw_only=True)
class _OutputRecord:
    """A record line of the output, kept or rewritten, without its end."""

    contig: str
    pos_1based: int
    line: str
    is_rewritten: bool


@dataclass(frozen=True, slots=True)
class _Edit:
    """What one VCF allele does to the reference, trimmed left first
    (doc/specs/reconciliation.md §4.2): deleted, the bases from start that it replaces, none
    for an insertion, which goes in before start; and inserted, the bases it puts there."""

    start: int
    deleted: str
    inserted: str

    @property
    def end(self) -> int:
        return self.start + len(self.deleted)


@dataclass(frozen=True, slots=True)
class _Piece:
    """A record the decomposition writes: one replaced base or one indel, at start, 0-based."""

    start: int
    ref: str
    alt: str


def respell_by_decomposing(
    vcf: Path, reference: IndexedFasta, output: Path, *, ploidy: int = DEFAULT_PLOIDY
) -> None:
    """Respell vcf by decomposing each record that is not kept into records of one replaced
    base or one indel, its pieces, written to output (doc/specs/respell.md §3).

    vcf: a VCF file, plain or gzipped, read with parse_vcf over reference with ploidy.
    output: the file written, bgzipped when vcf is gzipped.

    Raises what parse_vcf raises, before output is written, and ValueError when output is
    vcf itself.
    """
    _check_output(vcf, output)
    respelling_input = _parse_input(vcf, reference, ploidy)
    record_index = _make_record_index(respelling_input.input_records)
    output_records = []
    for input_record in respelling_input.input_records:
        vcf_record = input_record.vcf_record
        if vcf_record is None or _is_kept_whatever_is_written(
            vcf_record, input_record.place, record_index.touched_places
        ):
            output_records.append(_make_kept_record(input_record))
            continue
        pieces = _make_pieces(vcf_record, reference)
        written = [
            (vcf_record.start, vcf_record.alleles),
            *((piece.start, (piece.ref, piece.alt)) for piece in pieces),
        ]
        if not pieces or _is_near_spanning_deletion(
            vcf_record.contig, written, record_index.spanning_deletions, reference
        ):
            output_records.append(_make_kept_record(input_record))
            continue
        new_phase_set = _make_new_phase_set(respelling_input, input_record.place)
        output_records.extend(
            _make_piece_record(input_record, piece, allele_indices, new_phase_set)
            for piece, allele_indices in pieces.items()
        )
    _write_output(output, respelling_input, output_records)


def respell_by_shifting(
    vcf: Path,
    reference: IndexedFasta,
    output: Path,
    *,
    seed: int | None = None,
    ploidy: int = DEFAULT_PLOIDY,
) -> None:
    """Respell vcf by moving each record that is not kept, of one ALT whose edit is an
    insertion or a deletion, further right in its repeat, written to output
    (doc/specs/respell.md §4). A record is moved at most as far as its span, the contig's end
    and the edits of other records that may share a haplotype carrying it allow.

    vcf: a VCF file, plain or gzipped, read with parse_vcf over reference with ploidy.
    output: the file written, bgzipped when vcf is gzipped.
    seed: None to move each record as far as it goes; a whole number to move it by a
        distance drawn evenly from 1 to that, with random.Random(seed), so that the same seed
        writes the same file.

    Raises what parse_vcf raises, before output is written; ValueError when output is vcf
    itself, or seed is neither None nor a whole number.
    """
    draws = _make_draws(seed)
    _check_output(vcf, output)
    respelling_input = _parse_input(vcf, reference, ploidy)
    record_index = _make_record_index(respelling_input.input_records)
    output_records = []
    for input_record in respelling_input.input_records:
        vcf_record = input_record.vcf_record
        if (
            vcf_record is None
            or len(vcf_record.alleles) != 2
            or _is_kept_whatever_is_written(
                vcf_record, input_record.place, record_index.touched_places
            )
        ):
            output_records.append(_make_kept_record(input_record))
            continue
        record_edits = record_index.by_place[input_record.place]
        largest_shift = _count_largest_shift(
            record_edits, record_index.by_contig[vcf_record.contig], reference
        )
        if largest_shift < 1:
            output_records.append(_make_kept_record(input_record))
            continue
        shift = largest_shift if draws is None else draws.randint(1, largest_shift)
        start, alleles = _make_shifted_alleles(record_edits, shift, reference)
        written = [(vcf_record.start, vcf_record.alleles), (start, alleles)]
        if _is_near_spanning_deletion(
            vcf_record.contig, written, record_index.spanning_deletions, reference
        ):
            output_records.append(_make_kept_record(input_record))
            continue
        output_records.append(_make_shifted_record(input_record, start, alleles))
    _write_output(output, respelling_input, output_records)


def respell_with_spanning_deletions(
    vcf: Path, reference: IndexedFasta, output: Path, *, ploidy: int = DEFAULT_PLOIDY
) -> None:
    """Respell vcf by adding a `*` ALT, as GATK writes it, to each record that is not kept
    where a haplotype it gives REF carries, at another record, a deletion of a position of
    its REF, in a phase the file makes known; written to output (doc/specs/respell.md §5).
    The record joins the deletion's phase unit where it is all REF, and an unphased
    heterozygous deletion is then phased in a new phase set of its own.

    vcf: a VCF file, plain or gzipped, read with parse_vcf over reference with ploidy.
    output: the file written, bgzipped when vcf is gzipped.

    Raises what parse_vcf raises, before output is written, and ValueError when output is
    vcf itself.
    """
    _check_output(vcf, output)
    respelling_input = _parse_input(vcf, reference, ploidy)
    record_index = _make_record_index(respelling_input.input_records)
    rewritable_places = {
        place
        for place, record_edits in record_index.by_place.items()
        if not _is_kept_whatever_is_written(
            record_edits.vcf_record, place, record_index.touched_places
        )
        and not _is_near_spanning_deletion(
            record_edits.vcf_record.contig,
            [(record_edits.vcf_record.start, record_edits.vcf_record.alleles)],
            record_index.spanning_deletions,
            reference,
        )
    }
    sample_fields, star_places = _find_spanning_deletion_rewrites(
        record_index, rewritable_places, respelling_input
    )
    output_records = []
    for input_record in respelling_input.input_records:
        place = input_record.place
        if place in sample_fields:
            output_records.append(
                _make_spanned_record(
                    input_record, sample_fields[place], has_star=place in star_places
                )
            )
        else:
            output_records.append(_make_kept_record(input_record))
    _write_output(output, respelling_input, output_records)


def respell_by_merging(
    vcf: Path,
    reference: IndexedFasta,
    output: Path,
    *,
    seed: int | None = None,
    ploidy: int = DEFAULT_PLOIDY,
) -> None:
    """Respell vcf by merging runs of neighbouring records into one complex record each,
    written to output (doc/specs/respell.md §6). A run is two or more consecutive records of
    a contig, each with its REF overlapping the stretch of those before it or starting at
    most 2 positions after it; it is merged when none of its records is kept, has a `*` ALT
    or gives a sample `.`, no other record with an ALT touches its stretch, and each
    sample's heterozygous VCF genotypes in it are in one phase unit.

    vcf: a VCF file, plain or gzipped, read with parse_vcf over reference with ploidy.
    output: the file written, bgzipped when vcf is gzipped.
    seed: None for the longest runs that merge, taken from the left; a whole number to
        bound each run, from its first record, by records that each join with probability
        1/2, drawn with random.Random(seed), so that the same seed writes the same file;
        the longest run within that bound that merges is merged.

    Raises what parse_vcf raises, before output is written; ValueError when output is vcf
    itself, or seed is neither None nor a whole number.
    """
    draws = _make_draws(seed)
    _check_output(vcf, output)
    respelling_input = _parse_input(vcf, reference, ploidy)
    record_index = _make_record_index(respelling_input.input_records)
    input_records = respelling_input.input_records
    mergeable_places = {
        place
        for place, record_edits in record_index.by_place.items()
        if _is_mergeable(record_edits, record_index, reference)
    }
    output_records = []
    first = 0
    while first < len(input_records):
        chained = _find_chained_records(input_records, first, mergeable_places, record_index)
        num_bound = len(chained) if draws is None else _count_drawn_records(len(chained), draws)
        merged = _find_longest_merge(chained[:num_bound], input_records, record_index, reference)
        if merged is None:
            output_records.append(_make_kept_record(input_records[first]))
            first += 1
            continue
        num_records, merged_record = merged
        output_records.append(merged_record)
        first += num_records
    _write_output(output, respelling_input, output_records)


def _make_draws(seed: int | None) -> random.Random | None:
    """The draws of a re-spelling that chooses, from its seed (§7): None for no seed.

    Raises ValueError for a seed that is neither None nor a whole number.
    """
    # bool is an int to Python, and random.Random takes a float or a str as a seed too.
    if seed is not None and (not isinstance(seed, int) or isinstance(seed, bool)):
        raise ValueError(f"the seed must be None or a whole number, not {seed!r}")
    return None if seed is None else random.Random(seed)


# --- Reading (§2.1) ---


def _check_output(vcf: Path, output: Path) -> None:
    """output must not be vcf, under its own name or another: a hard link, or a name that
    differs in case on a file system that ignores it."""
    if output.resolve() == vcf.resolve() or (output.exists() and output.samefile(vcf)):
        raise ValueError(f"the output {output} is the input VCF; give another path to write to")


def _parse_input(vcf: Path, reference: IndexedFasta, ploidy: int) -> _RespellingInput:
    """vcf read with parse_vcf, which raises before anything is written for an input pvc_eval
    cannot read, and as lines, each record line matched to its record by its number."""
    records_by_line = {
        vcf_record.line_number: vcf_record
        for vcf_record in parse_vcf(vcf, reference, span_limit=DEFAULT_SPAN_LIMIT, ploidy=ploidy)
    }
    with vcf.open("rb") as raw:
        is_gzipped = raw.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
    header_lines: list[str] = []
    input_records: list[_InputRecord] = []
    phase_sets: set[str] = set()
    with _open_text(vcf, is_gzipped=is_gzipped) as lines:
        # Numbered as parse_vcf numbers them: universal newlines, a byte-order mark dropped.
        for line_number, line_with_end in enumerate(lines, start=1):
            line = line_with_end.removesuffix("\n")
            if line.startswith("#"):
                header_lines.append(line)
                continue
            if not line.strip():
                continue
            fields = tuple(line.split("\t"))
            phase_sets.update(_find_phase_sets(fields))
            # parse_vcf gives no record for a line with a symbolic allele.
            vcf_record = records_by_line.pop(line_number, None)
            if vcf_record is not None:
                _check_line_matches(fields, line_number, vcf_record)
            input_records.append(
                _InputRecord(
                    line=line, place=len(input_records) + 1, fields=fields, vcf_record=vcf_record
                )
            )
    assert not records_by_line, "every record of parse_vcf is on a line read"
    return _RespellingInput(
        header_lines=tuple(header_lines),
        input_records=tuple(input_records),
        is_gzipped=is_gzipped,
        taken_phase_set_numbers=_find_taken_phase_set_numbers(
            phase_sets, len(input_records) + len(phase_sets)
        ),
        are_phase_sets_integer=all(map(_is_integer_phase_set, phase_sets)),
    )


def _check_line_matches(fields: Sequence[str], line_number: int, vcf_record: VcfRecord) -> None:
    """The record parse_vcf gives for a line is on the line's contig, and its REF, trimmed,
    lies in the line's REF where its start puts it, which is the line's POS once trimmed.

    Raises AssertionError otherwise: the two readers number the lines differently, which is
    a defect of this module, not of the input.
    """
    offset = vcf_record.start - (int(fields[_POS_COLUMN]) - 1)
    trimmed_ref = vcf_record.alleles[0]
    line_ref = fields[_REF_COLUMN].upper()
    if (
        fields[_CHROM_COLUMN] != vcf_record.contig
        or offset < 0
        or line_ref[offset : offset + len(trimmed_ref)] != trimmed_ref
    ):
        raise AssertionError(
            f"line {line_number}, {fields[_CHROM_COLUMN]}:{fields[_POS_COLUMN]}, is matched to"
            f" the record parse_vcf gives at {vcf_record.contig}:{vcf_record.start + 1},"
            " which is not on it"
        )


def _open_text(vcf: Path, *, is_gzipped: bool) -> TextIO:
    if is_gzipped:
        return gzip.open(vcf, "rt", encoding="utf-8-sig")
    return vcf.open(encoding="utf-8-sig")


def _find_phase_sets(fields: Sequence[str]) -> Iterator[str]:
    """The PS values of a record line's samples, but `.`; the last PS key counts, as in
    parse_vcf."""
    if len(fields) <= _NUM_FIXED_COLUMNS:
        return
    key_indices = {key: index for index, key in enumerate(fields[_FORMAT_COLUMN].split(":"))}
    if "PS" not in key_indices:
        return
    ps_index = key_indices["PS"]
    for sample_field in fields[_NUM_FIXED_COLUMNS:]:
        values = sample_field.split(":")
        if ps_index < len(values) and values[ps_index] != ".":
            yield values[ps_index]


def _find_taken_phase_set_numbers(phase_sets: Iterable[str], largest: int) -> tuple[int, ...]:
    """The whole numbers from 1 to largest whose text, as parse_vcf compares it, is one of
    phase_sets, sorted: `7` takes 7, and `007` takes none. Text longer than largest's takes
    none a new phase set can reach, and is not turned into a number."""
    max_digits = len(str(largest))
    return tuple(
        sorted(
            int(phase_set)
            for phase_set in phase_sets
            if phase_set.isascii()
            and phase_set.isdecimal()
            and not phase_set.startswith("0")
            and len(phase_set) <= max_digits
            and int(phase_set) <= largest
        )
    )


def _make_new_phase_set(respelling_input: _RespellingInput, place: int) -> str:
    """The new phase set of the record at place among the input's record lines, from 1: the
    place-th whole number, from 1, that no PS of the input takes (§2.3). It is at most the
    number of record lines plus that of distinct PS values, which _RespellingInput's taken
    numbers reach."""
    taken = respelling_input.taken_phase_set_numbers
    # The place-th free number x is the least with x minus the taken numbers up to x equal
    # to place; counting up from place reaches it.
    number = place
    while (next_number := place + bisect.bisect_right(taken, number)) != number:
        number = next_number
    return str(number)


def _parse_phase_set_number(phase_set: str) -> int | None:
    """A PS as a whole number, or None when it is not one: not decimal digits, or more digits,
    leading zeros aside, than Python turns into a number (§2.3)."""
    if not (phase_set.isascii() and phase_set.isdecimal()):
        return None
    digits = phase_set.lstrip("0") or "0"
    if len(digits) > sys.int_info.default_max_str_digits:
        return None
    return int(digits)


# --- The records every re-spelling keeps (§2.2) ---


def _is_kept_whatever_is_written(
    vcf_record: VcfRecord, place: int, touched_places: set[int]
) -> bool:
    """Whether a record that parse_vcf gives, at place among the record lines, is kept for a
    row of §2.2 that does not depend on what the re-spelling writes: no ALT, a `*` ALT, an
    invalid GT or a touch. The record with a symbolic allele, which parse_vcf does not give,
    and the record near a `*` are left to the caller."""
    alts = vcf_record.alleles[1:]
    if not alts or _SPANNING_DELETION in alts:
        return True
    if any(not isinstance(genotype, VcfGenotype) for genotype in vcf_record.genotypes):
        return True
    return place in touched_places


def _find_edits(vcf_record: VcfRecord) -> tuple[_Edit | None, ...]:
    """The edit of each VCF allele of a record, by its index: None for REF, a `*` and an
    ALT that spells REF."""
    ref = vcf_record.alleles[0]
    edits: list[_Edit | None] = [None]
    for alt in vcf_record.alleles[1:]:
        num_prefix, deleted, inserted = trim_alt(ref, alt)
        is_edit = alt != _SPANNING_DELETION and (deleted or inserted)
        edits.append(_Edit(vcf_record.start + num_prefix, deleted, inserted) if is_edit else None)
    return tuple(edits)


def _is_edit_touching(edit: _Edit, other: _Edit) -> bool:
    """§1: two edits touch when their stretches overlap or lie next to each other; an
    insertion's stretch is the point it goes in at, which touches both positions beside it."""
    return edit.start <= other.end and other.start <= edit.end


@dataclass(frozen=True, slots=True)
class _PhaseUnit:
    """A sample's phase unit at a record (§1): the phase set of a phased GT, None for one
    phased with no PS; or, for an unphased GT, the record's place among the input's record
    lines, a unit of its own."""

    phase_set: str | None
    place: int | None


def _find_phase_unit(genotype: VcfGenotype, place: int) -> _PhaseUnit:
    """The phase unit of a sample's VCF genotype at the record at place. A phased GT with no
    PS is in the one phase set VCF gives those of its sample and contig."""
    if genotype.phased:
        return _PhaseUnit(phase_set=genotype.phase_set, place=None)
    return _PhaseUnit(phase_set=None, place=place)


def _is_haplotype_shareable(
    genotype: VcfGenotype,
    place: int,
    haplotype: int,
    other_genotype: VcfGenotype,
    other_place: int,
    other_haplotype: int,
) -> bool:
    """§1: two VCF alleles of one sample, at the records at place and other_place, may be on
    one haplotype unless both VCF genotypes are in one phase unit, phased in one phase set,
    and the alleles are on different haplotypes."""
    is_one_unit = _find_phase_unit(genotype, place) == _find_phase_unit(other_genotype, other_place)
    return not is_one_unit or haplotype == other_haplotype


@dataclass(frozen=True, slots=True)
class _RecordEdits:
    """A record that parse_vcf gives, with an ALT, its place among the input's record lines,
    from 1, and the edit of each of its VCF alleles, by index (_find_edits)."""

    place: int
    vcf_record: VcfRecord
    edits: tuple[_Edit | None, ...]


@dataclass(frozen=True, slots=True)
class _ContigEdits:
    """The records of a contig with an ALT, with their edits, sorted by REF start; and for
    each, the furthest end of the REFs up to it, so that a search for the REFs that touch a
    stretch stops where no REF before reaches it."""

    record_edits: tuple[_RecordEdits, ...]
    starts: tuple[int, ...]
    furthest_ends: tuple[int, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class _RecordIndex:
    """The records of the input with an ALT, with their edits, found as the rows of §2.2 need
    them: by place; by contig, sorted by REF; those with a `*` ALT, by contig; and the places
    of the touched records."""

    by_place: Mapping[int, _RecordEdits]
    by_contig: Mapping[str, _ContigEdits]
    spanning_deletions: Mapping[str, _ContigEdits]
    touched_places: set[int]


def _make_record_index(input_records: Iterable[_InputRecord]) -> _RecordIndex:
    by_place = {
        input_record.place: _RecordEdits(
            input_record.place, input_record.vcf_record, _find_edits(input_record.vcf_record)
        )
        for input_record in input_records
        if input_record.vcf_record is not None and len(input_record.vcf_record.alleles) > 1
    }
    by_contig = _make_contig_edits(by_place.values())
    return _RecordIndex(
        by_place=by_place,
        by_contig=by_contig,
        spanning_deletions=_make_contig_edits(
            record_edits
            for record_edits in by_place.values()
            if _SPANNING_DELETION in record_edits.vcf_record.alleles[1:]
        ),
        touched_places=_find_touched_places(by_contig),
    )


def _make_contig_edits(record_edits: Iterable[_RecordEdits]) -> dict[str, _ContigEdits]:
    """Records with their edits, by contig, sorted by REF start."""
    by_contig: dict[str, list[_RecordEdits]] = {}
    for one_record_edits in record_edits:
        by_contig.setdefault(one_record_edits.vcf_record.contig, []).append(one_record_edits)
    contig_edits = {}
    for contig, contig_records in by_contig.items():
        contig_records.sort(key=lambda record_edits: record_edits.vcf_record.start)
        furthest_ends = itertools.accumulate(
            (record_edits.vcf_record.end for record_edits in contig_records), max
        )
        contig_edits[contig] = _ContigEdits(
            record_edits=tuple(contig_records),
            starts=tuple(record_edits.vcf_record.start for record_edits in contig_records),
            furthest_ends=tuple(furthest_ends),
        )
    return contig_edits


def _find_ref_touching_records(contig_edits: _ContigEdits, stretch: Span) -> Iterator[_RecordEdits]:
    """The records whose REF overlaps stretch or lies next to it, from the right. Every edit
    that touches stretch is one of theirs, since an edit lies within its REF, ends included."""
    stretch_start, stretch_end = stretch
    index = bisect.bisect_right(contig_edits.starts, stretch_end) - 1
    while index >= 0 and contig_edits.furthest_ends[index] >= stretch_start:
        record_edits = contig_edits.record_edits[index]
        if record_edits.vcf_record.end >= stretch_start:
            yield record_edits
        index -= 1


def _find_touched_places(contig_edits: Mapping[str, _ContigEdits]) -> set[int]:
    """The places of the touched records (§2.2): those at which some sample carries an edit
    that touches an edit it carries at another record of the contig, the two VCF alleles
    being ones that may share a haplotype.

    An edit lies within its record's REF, its ends included, so each record is compared only
    with those whose REF touches its own (§10).
    """
    touched_places: set[int] = set()
    for one_contig_edits in contig_edits.values():
        for record_edits in one_contig_edits.record_edits:
            vcf_record = record_edits.vcf_record
            ref = (vcf_record.start, vcf_record.end)
            if any(
                other.place != record_edits.place and _is_carried_edit_touching(record_edits, other)
                for other in _find_ref_touching_records(one_contig_edits, ref)
            ):
                touched_places.add(record_edits.place)
    return touched_places


def _is_carried_edit_touching(first: _RecordEdits, second: _RecordEdits) -> bool:
    """Whether some sample carries edits at the two records that touch, on haplotypes they
    may share."""
    genotype_pairs = zip(first.vcf_record.genotypes, second.vcf_record.genotypes, strict=True)
    for genotype, other_genotype in genotype_pairs:
        if not (isinstance(genotype, VcfGenotype) and isinstance(other_genotype, VcfGenotype)):
            continue
        for haplotype, edit in _find_carried_edits(genotype, first.edits):
            for other_haplotype, other_edit in _find_carried_edits(other_genotype, second.edits):
                if _is_edit_touching(edit, other_edit) and _is_haplotype_shareable(
                    genotype, first.place, haplotype, other_genotype, second.place, other_haplotype
                ):
                    return True
    return False


def _find_carried_edits(
    genotype: VcfGenotype, edits: Sequence[_Edit | None]
) -> Iterator[tuple[int, _Edit]]:
    """Each haplotype that carries an edit, with it; `.` and REF carry none."""
    for haplotype, allele in enumerate(genotype.alleles):
        if allele is MISSING:
            continue
        edit = edits[allele]
        if edit is not None:
            yield haplotype, edit


def _is_near_spanning_deletion(
    contig: str,
    written: Iterable[tuple[int, Sequence[str]]],
    spanning_deletions: Mapping[str, _ContigEdits],
    reference: IndexedFasta,
) -> bool:
    """§2.2: the reach of a record of contig, or of a record its rewriting would write,
    touches the REF of a record with a `*` ALT.

    written: the input record and those its rewriting would write, each as its start,
        0-based, and its alleles, trimmed as parse_vcf trims them.
    """
    # A contig with no record with a `*` ALT has no REFs to touch.
    contig_spanning_deletions = spanning_deletions.get(contig)
    if contig_spanning_deletions is None:
        return False
    return any(
        next(
            _find_ref_touching_records(
                contig_spanning_deletions, _find_reach(reference, contig, start, alleles)
            ),
            None,
        )
        is not None
        for start, alleles in written
    )


# --- The reach (§1) ---

# The bases a slide reads from the reference at a time.
_SLIDE_BASES = 64


def _find_reach(reference: IndexedFasta, contig: str, start: int, alleles: Sequence[str]) -> Span:
    """The reach of a record (§1): its REF, widened for each ALT whose edit is an insertion
    or a deletion by every position it can slide to (doc/spec.md §5.2), with no limit.

    start: where the REF starts, 0-based. alleles: REF and the ALTs, trimmed as parse_vcf
    trims them.
    """
    ref = alleles[0]
    reach_start, reach_end = start, start + len(ref)
    contig_length = reference.get_length(contig)
    for alt in alleles[1:]:
        if alt == _SPANNING_DELETION:
            continue
        num_prefix, deleted, inserted = trim_alt(ref, alt)
        if bool(deleted) == bool(inserted):
            continue
        event = deleted or inserted
        # The event deletes [position, after), or is inserted before position, where after
        # is position too.
        position = start + num_prefix
        after = position + len(deleted)
        num_left = _count_slide_left(reference, contig, position, event)
        num_right = _count_slide_right(reference, contig, after, event)
        # The anchor of the left-aligned event, if the contig has a base before it; and the
        # last base deleted, or the base after the point of the insertion, slid right.
        reach_start = min(reach_start, max(position - num_left - 1, 0))
        slid_end = after + num_right + (1 if inserted else 0)
        reach_end = max(reach_end, min(slid_end, contig_length))
    return reach_start, reach_end


def _count_slide_left(reference: IndexedFasta, contig: str, position: int, event: str) -> int:
    """How far an event before position slides left: while the base before it equals its
    last base, which then moves to its start."""
    num_left = 0
    while position - num_left > 0:
        stretch_end = position - num_left
        bases = reference.get_bases(contig, max(0, stretch_end - _SLIDE_BASES), stretch_end)
        for base in reversed(bases):
            if base != event[(-1 - num_left) % len(event)]:
                return num_left
            num_left += 1
    return num_left


def _count_slide_right(reference: IndexedFasta, contig: str, after: int, event: str) -> int:
    """How far an event that ends at after slides right: while the base after it equals its
    first base, which then moves to its end."""
    contig_length = reference.get_length(contig)
    num_right = 0
    while after + num_right < contig_length:
        stretch_start = after + num_right
        bases = reference.get_bases(
            contig, stretch_start, min(contig_length, stretch_start + _SLIDE_BASES)
        )
        for base in bases:
            if base != event[num_right % len(event)]:
                return num_right
            num_right += 1
    return num_right


# --- Decomposition (§3) ---


def _make_pieces(vcf_record: VcfRecord, reference: IndexedFasta) -> dict[_Piece, set[int]]:
    """The pieces of a record, in the order of its ALTs, each with the indices of the ALTs
    that make it; a piece two ALTs make is one."""
    allele_indices_of_piece: dict[_Piece, set[int]] = {}
    for allele_index, edit in enumerate(_find_edits(vcf_record)):
        if edit is None:
            continue
        for piece in _make_edit_pieces(vcf_record.contig, edit, reference):
            allele_indices_of_piece.setdefault(piece, set()).add(allele_index)
    return allele_indices_of_piece


def _make_edit_pieces(contig: str, edit: _Edit, reference: IndexedFasta) -> Iterator[_Piece]:
    """The pieces of one edit: its replaced bases and its new ones aligned from the left, a
    SNP where they differ, then the bases left over, inserted after the last aligned
    position or deleted, with the base before as the anchor."""
    deleted, inserted = edit.deleted, edit.inserted
    num_aligned = min(len(deleted), len(inserted))
    for offset in range(num_aligned):
        if deleted[offset] != inserted[offset]:
            yield _Piece(edit.start + offset, deleted[offset], inserted[offset])
    if len(deleted) == len(inserted):
        return
    anchor = edit.start + num_aligned - 1
    left_deleted, left_inserted = deleted[num_aligned:], inserted[num_aligned:]
    if anchor < 0:
        # At the contig's first base, which has none before it, VCF anchors on the base
        # after the indel, at POS 1.
        after = len(left_deleted)
        after_base = reference.get_bases(contig, after, after + 1)
        yield _Piece(0, left_deleted + after_base, left_inserted + after_base)
        return
    anchor_base = reference.get_bases(contig, anchor, anchor + 1)
    yield _Piece(anchor, anchor_base + left_deleted, anchor_base + left_inserted)


def _make_piece_record(
    input_record: _InputRecord, piece: _Piece, allele_indices: set[int], new_phase_set: str
) -> _OutputRecord:
    """A piece written with its VCF genotypes carried over: 1 where the haplotype's ALT makes
    the piece, 0 where it is REF or another ALT, `.` where it is `.`; phased as the record
    was, or in the record's new phase set when it was not, so that the pieces are one phase
    unit as the record was."""
    vcf_record = input_record.vcf_record
    assert vcf_record is not None, "a record with a symbolic allele is kept"
    sample_fields = []
    for genotype in vcf_record.genotypes:
        assert isinstance(genotype, VcfGenotype), "a record with an invalid GT is kept"
        alleles = [
            "." if allele is MISSING else str(int(allele in allele_indices))
            for allele in genotype.alleles
        ]
        if len(alleles) == 1:
            sample_fields.append(_format_sample(alleles, phased=False, phase_set=None))
            continue
        phase_set = genotype.phase_set if genotype.phased else new_phase_set
        sample_fields.append(_format_sample(alleles, phased=True, phase_set=phase_set))
    return _make_rewritten_record(
        input_record,
        piece.start + 1,
        piece.ref,
        [piece.alt],
        input_record.fields[_QUAL_COLUMN],
        sample_fields,
    )


# --- Shift (§4) ---


def _count_largest_shift(
    record_edits: _RecordEdits, contig_edits: _ContigEdits, reference: IndexedFasta
) -> int:
    """How far right the one edit of a record, an insertion or a deletion, can be moved: as
    far as it slides, within the record's span and the contig, and short of every edit of
    another record that a sample carries on a haplotype that may share one carrying the
    indel; 0 when it cannot be moved, or is no indel."""
    vcf_record = record_edits.vcf_record
    edit = record_edits.edits[1]
    if edit is None or bool(edit.deleted) == bool(edit.inserted):
        return 0
    contig_length = reference.get_length(vcf_record.contig)
    num_slide = _count_slide_right(
        reference, vcf_record.contig, edit.end, edit.deleted or edit.inserted
    )
    largest_shift = min(num_slide, _count_span_shift(edit, vcf_record.span, contig_length))
    largest_shift = _count_plain_shift(reference, vcf_record.contig, edit, largest_shift)
    if largest_shift < 1:
        return 0
    # The stretch the indel crosses moved by largest_shift: a to b + k for a deletion of a
    # to b, and the positions its point passes over for an insertion.
    crossed = (edit.start, edit.end + largest_shift)
    for other in _find_ref_touching_records(contig_edits, crossed):
        if other.place == record_edits.place:
            continue
        for other_edit in _find_blocking_edits(record_edits, other):
            # Crossing up to edit.end + k touches an edit that starts at edit.end + k + 1 or
            # before; one that ends before edit.start lies apart whatever k is.
            if other_edit.end >= edit.start:
                largest_shift = min(largest_shift, other_edit.start - edit.end - 1)
    return max(largest_shift, 0)


def _count_plain_shift(
    reference: IndexedFasta, contig: str, edit: _Edit, largest_shift: int
) -> int:
    """How far, up to largest_shift, an indel can move over the reference's plain bases (§4,
    the bases): its deleted bases must be the reference's own and A, C, G or T, and so must
    every position it slides past. parse_vcf accepts a REF A over an M or an N, but
    reconciliation builds a haplotype from the reference's bases."""
    if edit.deleted:
        deleted_bases = reference.get_bases(contig, edit.start, edit.end)
        if deleted_bases != edit.deleted or not _PLAIN_BASES.issuperset(deleted_bases):
            return 0
    slid_past = reference.get_bases(contig, edit.end, edit.end + largest_shift)
    return next(
        (index for index, base in enumerate(slid_past) if base not in _PLAIN_BASES),
        largest_shift,
    )


def _count_span_shift(edit: _Edit, span: Span, contig_length: int) -> int:
    """How far an indel can move with its placement within span (§4): its deleted positions,
    or the base after its point, which the contig's end cuts as it cuts the span. A span the
    span limit cut back to the REF allows none."""
    _, span_end = span
    if edit.deleted:
        return span_end - edit.end
    if span_end == contig_length:
        return contig_length - edit.start
    return span_end - edit.start - 1


def _find_blocking_edits(record_edits: _RecordEdits, other: _RecordEdits) -> Iterator[_Edit]:
    """The edits a sample carries at other on a haplotype that may share one on which it
    carries the indel of record_edits."""
    genotype_pairs = zip(record_edits.vcf_record.genotypes, other.vcf_record.genotypes, strict=True)
    for genotype, other_genotype in genotype_pairs:
        if not (isinstance(genotype, VcfGenotype) and isinstance(other_genotype, VcfGenotype)):
            continue
        carrying = [haplotype for haplotype, _ in _find_carried_edits(genotype, record_edits.edits)]
        for other_haplotype, other_edit in _find_carried_edits(other_genotype, other.edits):
            if any(
                _is_haplotype_shareable(
                    genotype,
                    record_edits.place,
                    haplotype,
                    other_genotype,
                    other.place,
                    other_haplotype,
                )
                for haplotype in carrying
            ):
                yield other_edit


def _make_shifted_alleles(
    record_edits: _RecordEdits, shift: int, reference: IndexedFasta
) -> tuple[int, tuple[str, str]]:
    """The start, 0-based, and the REF and ALT of a record's indel moved right by shift, with
    the base before it as anchor. An inserted sequence turns as it slides: its first base
    moves to its end at each step."""
    vcf_record = record_edits.vcf_record
    edit = record_edits.edits[1]
    assert edit is not None, "only an indel is shifted"
    anchor = edit.start + shift - 1
    if edit.deleted:
        ref = reference.get_bases(vcf_record.contig, anchor, edit.end + shift)
        return anchor, (ref, ref[0])
    anchor_base = reference.get_bases(vcf_record.contig, anchor, anchor + 1)
    turn = shift % len(edit.inserted)
    return anchor, (anchor_base, anchor_base + edit.inserted[turn:] + edit.inserted[:turn])


def _make_shifted_record(
    input_record: _InputRecord, start: int, alleles: tuple[str, str]
) -> _OutputRecord:
    """A shifted record, with its VCF genotypes as read: the shift moves an edit, not a
    haplotype's phase."""
    vcf_record = input_record.vcf_record
    assert vcf_record is not None, "a record with a symbolic allele is kept"
    sample_fields = []
    for genotype in vcf_record.genotypes:
        assert isinstance(genotype, VcfGenotype), "a record with an invalid GT is kept"
        sample_fields.append(
            _format_genotype(genotype.alleles, phased=genotype.phased, phase_set=genotype.phase_set)
        )
    ref, alt = alleles
    return _make_rewritten_record(
        input_record, start + 1, ref, [alt], input_record.fields[_QUAL_COLUMN], sample_fields
    )


# --- Spanning deletions (§5) ---


@dataclass(frozen=True, slots=True)
class _Deletion:
    """A deletion a sample carries over a record's REF: the record that makes it, the
    sample's VCF genotype there and the haplotype that carries it."""

    record_edits: _RecordEdits
    genotype: VcfGenotype
    haplotype: int


@dataclass(frozen=True, slots=True)
class _Star:
    """The haplotypes of a sample's VCF genotype that get a `*`, and the phase its GT is
    written with: phased, and the phase set written; and rephased, the unphased deletion
    phased in its new phase set, phase_set, which the GT joins, or None."""

    haplotypes: frozenset[int]
    phased: bool
    phase_set: str | None
    rephased: _RecordEdits | None = None


def _find_spanning_deletion_rewrites(
    record_index: _RecordIndex, rewritable_places: set[int], respelling_input: _RespellingInput
) -> tuple[dict[int, dict[int, str]], set[int]]:
    """What the `*` re-spelling writes (§5): the records it rewrites, by place, each with the
    sample fields, `GT:PS` by sample index, that are not the record's as read; and the places
    of those that get a `*` ALT. A deletion phased in its new phase set is rewritten with no
    `*`, unless another sample gives it one."""
    sample_fields: dict[int, dict[int, str]] = {}
    star_places: set[int] = set()
    for place in sorted(rewritable_places):
        record_edits = record_index.by_place[place]
        vcf_record = record_edits.vcf_record
        contig_edits = record_index.by_contig[vcf_record.contig]
        star_index = len(vcf_record.alleles)
        deletions_by_sample = _find_deletions_over(record_edits, contig_edits)
        for sample, genotype in enumerate(vcf_record.genotypes):
            assert isinstance(genotype, VcfGenotype), "a record with an invalid GT is kept"
            star = _find_star(
                genotype,
                place,
                deletions_by_sample[sample],
                rewritable_places,
                respelling_input,
            )
            if star is None:
                continue
            if star.rephased is not None:
                rephased_genotype = star.rephased.vcf_record.genotypes[sample]
                assert isinstance(rephased_genotype, VcfGenotype), "it carries the deletion"
                rephased_field = _format_genotype(
                    rephased_genotype.alleles, phased=True, phase_set=star.phase_set
                )
                _add_sample_field(sample_fields, star.rephased.place, sample, rephased_field)
            alleles = [
                star_index if haplotype in star.haplotypes else allele
                for haplotype, allele in enumerate(genotype.alleles)
            ]
            field = _format_genotype(alleles, phased=star.phased, phase_set=star.phase_set)
            _add_sample_field(sample_fields, place, sample, field)
            star_places.add(place)
    return sample_fields, star_places


def _add_sample_field(
    sample_fields: dict[int, dict[int, str]], place: int, sample: int, field: str
) -> None:
    """A sample's field of a rewritten record. A sample whose record gets a `*` gives REF on
    a haplotype, so it is not the unphased deletion carried that a `*` rephases, and one
    record never gets two fields for one sample, though one deletion may be rephased for
    several records under it, in one field."""
    record_fields = sample_fields.setdefault(place, {})
    assert record_fields.get(sample, field) == field, "one field for each sample"
    record_fields[sample] = field


def _find_deletions_over(
    record_edits: _RecordEdits, contig_edits: _ContigEdits
) -> list[list[_Deletion]]:
    """For each sample, by index, the deletions it carries at other records, edits that
    replace positions with no bases, of a position of the REF of record_edits as parse_vcf
    trims it (§5)."""
    vcf_record = record_edits.vcf_record
    deletions_by_sample: list[list[_Deletion]] = [[] for _ in vcf_record.genotypes]
    for other in _find_ref_touching_records(contig_edits, (vcf_record.start, vcf_record.end)):
        if other.place == record_edits.place:
            continue
        for sample, genotype in enumerate(other.vcf_record.genotypes):
            if not isinstance(genotype, VcfGenotype):
                continue
            for haplotype, edit in _find_carried_edits(genotype, other.edits):
                is_deletion = bool(edit.deleted) and not edit.inserted
                if is_deletion and edit.start < vcf_record.end and vcf_record.start < edit.end:
                    deletions_by_sample[sample].append(_Deletion(other, genotype, haplotype))
    return deletions_by_sample


def _find_star(
    genotype: VcfGenotype,
    place: int,
    deletions: Sequence[_Deletion],
    rewritable_places: set[int],
    respelling_input: _RespellingInput,
) -> _Star | None:
    """The haplotypes of a sample's VCF genotype at the record at place that get a `*`
    (§5): those it gives REF and whose deletions over the record are in a phase that is
    known. None when there is none.

    deletions: the deletions the sample carries over the record's REF.
    rewritable_places: the places of the records that are not kept.
    """
    alleles = genotype.alleles
    if not deletions:
        return None
    if len(alleles) == 1:
        # The ploidy is 1, so the deletion is on the record's one haplotype.
        return _Star(frozenset({0}), phased=False, phase_set=None) if alleles[0] == 0 else None
    if all(allele == 0 for allele in alleles):
        return _find_all_ref_star(deletions, genotype, rewritable_places, respelling_input)
    # Both VCF genotypes phased in one phase set, the haplotypes the record gives REF.
    unit = _find_phase_unit(genotype, place)
    haplotypes = frozenset(
        deletion.haplotype
        for deletion in deletions
        if alleles[deletion.haplotype] == 0
        and _find_phase_unit(deletion.genotype, deletion.record_edits.place) == unit
    )
    if not haplotypes:
        return None
    return _Star(haplotypes, phased=True, phase_set=genotype.phase_set)


def _find_all_ref_star(
    deletions: Sequence[_Deletion],
    genotype: VcfGenotype,
    rewritable_places: set[int],
    respelling_input: _RespellingInput,
) -> _Star | None:
    """The `*` of a record that gives REF on every haplotype, which then joins the phase
    unit of the deletions (§5): their phase set when they are phased; its own phase, as
    read, under an unphased homozygous deletion, which every haplotype carries; and the new
    phase set of an unphased heterozygous deletion, which must not be kept. None when the
    deletions come from two or more phase units, which the record cannot join both of."""
    units = {
        _find_phase_unit(deletion.genotype, deletion.record_edits.place) for deletion in deletions
    }
    if len(units) > 1:
        return None
    haplotypes = frozenset(deletion.haplotype for deletion in deletions)
    deletion = deletions[0]
    if deletion.genotype.phased:
        return _Star(haplotypes, phased=True, phase_set=deletion.genotype.phase_set)
    if len(set(deletion.genotype.alleles)) == 1:
        return _Star(haplotypes, phased=genotype.phased, phase_set=genotype.phase_set)
    rephased = deletion.record_edits
    if rephased.place not in rewritable_places:
        return None
    new_phase_set = _make_new_phase_set(respelling_input, rephased.place)
    return _Star(haplotypes, phased=True, phase_set=new_phase_set, rephased=rephased)


def _make_spanned_record(
    input_record: _InputRecord, sample_fields: Mapping[int, str], *, has_star: bool
) -> _OutputRecord:
    """A record of the `*` re-spelling, with the input's POS, REF and ALTs, a `*` after them
    when has_star, and each sample's GT and PS as read but those of sample_fields."""
    vcf_record = input_record.vcf_record
    assert vcf_record is not None, "a record with a symbolic allele is kept"
    fields = []
    for sample, genotype in enumerate(vcf_record.genotypes):
        assert isinstance(genotype, VcfGenotype), "a record with an invalid GT is kept"
        field = sample_fields.get(sample)
        if field is None:
            field = _format_genotype(
                genotype.alleles, phased=genotype.phased, phase_set=genotype.phase_set
            )
        fields.append(field)
    alts = input_record.fields[_ALT_COLUMN].upper().split(",")
    if has_star:
        alts.append(_SPANNING_DELETION)
    return _make_rewritten_record(
        input_record,
        input_record.pos_1based,
        input_record.fields[_REF_COLUMN].upper(),
        alts,
        input_record.fields[_QUAL_COLUMN],
        fields,
    )


# --- Merge (§6) ---

# §6: a record joins a run when its REF starts at most this many positions after the
# stretch of the records before it (M1 has 1 between, M8 has 3).
_MAX_POSITIONS_BETWEEN = 2


def _is_mergeable(
    record_edits: _RecordEdits, record_index: _RecordIndex, reference: IndexedFasta
) -> bool:
    """Whether a record can be in a run (§6): it is not kept for a row of §2.2 of its own,
    has no `*` ALT and gives no sample `.`. The reach of the record a run writes is the
    run's to check."""
    vcf_record = record_edits.vcf_record
    if _is_kept_whatever_is_written(vcf_record, record_edits.place, record_index.touched_places):
        return False
    if any(
        allele is MISSING
        for genotype in vcf_record.genotypes
        if isinstance(genotype, VcfGenotype)
        for allele in genotype.alleles
    ):
        return False
    written = [(vcf_record.start, vcf_record.alleles)]
    return not _is_near_spanning_deletion(
        vcf_record.contig, written, record_index.spanning_deletions, reference
    )


def _find_chained_records(
    input_records: Sequence[_InputRecord],
    first: int,
    mergeable_places: set[int],
    record_index: _RecordIndex,
) -> list[_RecordEdits]:
    """The records a run from input_records[first] can hold: it and those after it,
    consecutive, of its contig and each mergeable, while each one's REF overlaps the stretch
    of those before it or starts at most _MAX_POSITIONS_BETWEEN positions after it. Empty
    when the first is not mergeable."""
    chained: list[_RecordEdits] = []
    stretch_start, stretch_end = 0, 0
    for input_record in itertools.islice(input_records, first, None):
        if input_record.place not in mergeable_places:
            break
        record_edits = record_index.by_place[input_record.place]
        vcf_record = record_edits.vcf_record
        start, end = vcf_record.start, vcf_record.end
        if chained:
            is_overlapping = start < stretch_end and stretch_start < end
            is_after = stretch_end <= start <= stretch_end + _MAX_POSITIONS_BETWEEN
            is_joining = vcf_record.contig == chained[0].vcf_record.contig and (
                is_overlapping or is_after
            )
            if not is_joining:
                break
            stretch_start, stretch_end = min(stretch_start, start), max(stretch_end, end)
        else:
            stretch_start, stretch_end = start, end
        chained.append(record_edits)
    return chained


def _count_drawn_records(num_chained: int, draws: random.Random) -> int:
    """How many of num_chained records, from the first, a seeded run may hold: each after
    the first joins with probability 1/2, and the first that does not ends it."""
    num_records = 1
    while num_records < num_chained and draws.random() < 0.5:
        num_records += 1
    return num_records


def _find_longest_merge(
    chained: Sequence[_RecordEdits],
    input_records: Sequence[_InputRecord],
    record_index: _RecordIndex,
    reference: IndexedFasta,
) -> tuple[int, _OutputRecord] | None:
    """The longest run of the chained records, from the first, that is merged, as how many
    records it holds and the record written for it; None when no run of two or more is."""
    for num_records in range(len(chained), 1, -1):
        merged_record = _make_merged_record(
            chained[:num_records], input_records, record_index, reference
        )
        if merged_record is not None:
            return num_records, merged_record
    return None


def _make_merged_record(
    run: Sequence[_RecordEdits],
    input_records: Sequence[_InputRecord],
    record_index: _RecordIndex,
    reference: IndexedFasta,
) -> _OutputRecord | None:
    """The record a run of mergeable records is written as (§6), or None when it is not
    merged: another record with an ALT touches its stretch; a sample's heterozygous VCF
    genotypes in it are in two phase units; every haplotype gets REF; or the reach of the
    record touches the REF of a record with a `*` ALT (§2.2)."""
    contig = run[0].vcf_record.contig
    stretch = (
        min(record_edits.vcf_record.start for record_edits in run),
        max(record_edits.vcf_record.end for record_edits in run),
    )
    run_places = {record_edits.place for record_edits in run}
    if any(
        other.place not in run_places
        for other in _find_ref_touching_records(record_index.by_contig[contig], stretch)
    ):
        return None
    stretch_start, stretch_end = stretch
    ref = reference.get_bases(contig, stretch_start, stretch_end)
    alt_indices: dict[str, int] = {}
    sample_fields = []
    for sample in range(len(run[0].vcf_record.genotypes)):
        genotypes = [_get_valid_genotype(record_edits, sample) for record_edits in run]
        phase = _find_merged_phase(genotypes, run)
        if phase is None:
            return None
        phased, phase_set = phase
        alleles = [
            0
            if merged_allele == ref
            else alt_indices.setdefault(merged_allele, len(alt_indices) + 1)
            for merged_allele in _make_merged_alleles(genotypes, run, ref, stretch_start)
        ]
        sample_fields.append(_format_genotype(alleles, phased=phased, phase_set=phase_set))
    alts = list(alt_indices)
    if not alts:
        return None
    written = [(stretch_start, (ref, *alts))]
    if _is_near_spanning_deletion(contig, written, record_index.spanning_deletions, reference):
        return None
    first_record = input_records[run[0].place - 1]
    qual = _find_lowest_qual([input_records[record_edits.place - 1] for record_edits in run])
    return _make_rewritten_record(first_record, stretch_start + 1, ref, alts, qual, sample_fields)


def _get_valid_genotype(record_edits: _RecordEdits, sample: int) -> VcfGenotype:
    genotype = record_edits.vcf_record.genotypes[sample]
    assert isinstance(genotype, VcfGenotype), "a record with an invalid GT is kept"
    return genotype


def _get_allele(genotype: VcfGenotype, haplotype: int) -> int:
    allele = genotype.alleles[haplotype]
    assert allele is not MISSING, "a record that gives a sample `.` is not merged"
    return allele


def _find_merged_phase(
    genotypes: Sequence[VcfGenotype], run: Sequence[_RecordEdits]
) -> tuple[bool, str | None] | None:
    """How a sample's GT of a merged run is written (§6): phased, in the phase set of the
    run's heterozygous VCF genotypes, when those are phased, and unphased otherwise; None
    when they are in two or more phase units, which merged would invent a phase. With
    ploidy 1 there is no phase."""
    heterozygous_units = {
        _find_phase_unit(genotype, record_edits.place)
        for genotype, record_edits in zip(genotypes, run, strict=True)
        if len(set(genotype.alleles)) > 1
    }
    if len(heterozygous_units) > 1:
        return None
    if not heterozygous_units:
        return False, None
    (unit,) = heterozygous_units
    if unit.place is not None:
        return False, None
    return True, unit.phase_set


def _make_merged_alleles(
    genotypes: Sequence[VcfGenotype], run: Sequence[_RecordEdits], ref: str, stretch_start: int
) -> Iterator[str]:
    """The allele of each of a sample's haplotypes over a run's stretch, in order: ref, the
    reference over the stretch from stretch_start, with the edits of the VCF alleles the
    haplotype carries at the run's records applied.

    genotypes: the sample's VCF genotypes at the run's records, in the order of run.
    """
    for haplotype in range(len(genotypes[0].alleles)):
        edits = [
            edit
            for genotype, record_edits in zip(genotypes, run, strict=True)
            if (edit := record_edits.edits[_get_allele(genotype, haplotype)]) is not None
        ]
        yield _make_merged_allele(ref, stretch_start, edits)


def _make_merged_allele(ref: str, stretch_start: int, edits: Sequence[_Edit]) -> str:
    """ref, the reference over a stretch from stretch_start, with edits applied in
    position order. The edits a haplotype carries in a run do not touch, or their records
    would be touched (§2.2), so none starts before the one before it ends."""
    pieces = []
    position = stretch_start
    for edit in sorted(edits, key=lambda edit: edit.start):
        assert edit.start >= position, "the edits of one haplotype in a run do not touch"
        pieces += [ref[position - stretch_start : edit.start - stretch_start], edit.inserted]
        position = edit.end
    pieces.append(ref[position - stretch_start :])
    return "".join(pieces)


def _find_lowest_qual(input_records: Sequence[_InputRecord]) -> str:
    """The QUAL of the records with the lowest one, as the input writes it, or `.` when
    none has one (§2.3)."""
    with_qual = [
        input_record
        for input_record in input_records
        if input_record.vcf_record is not None and input_record.vcf_record.qual is not None
    ]
    if not with_qual:
        return "."
    lowest = min(with_qual, key=lambda input_record: input_record.vcf_record.qual)
    return lowest.fields[_QUAL_COLUMN]


# --- Writing (§2.3) ---


def _make_kept_record(input_record: _InputRecord) -> _OutputRecord:
    return _OutputRecord(
        contig=input_record.contig,
        pos_1based=input_record.pos_1based,
        line=input_record.line,
        is_rewritten=False,
    )


def _format_sample(alleles: Sequence[str], *, phased: bool, phase_set: str | None) -> str:
    """A sample's `GT:PS`: the alleles joined by `|` when phased and `/` when not, and the
    phase set of a phased GT, `.` when it has none or the GT is unphased."""
    gt = ("|" if phased else "/").join(alleles)
    return f"{gt}:{phase_set if phased and phase_set is not None else '.'}"


def _format_genotype(
    alleles: Sequence[int | MISSING], *, phased: bool, phase_set: str | None
) -> str:
    """A sample's `GT:PS` from its alleles' indices, `.` for MISSING; a GT of one allele has
    no phase."""
    alleles_written = ["." if allele is MISSING else str(allele) for allele in alleles]
    return _format_sample(alleles_written, phased=phased and len(alleles) > 1, phase_set=phase_set)


def _make_rewritten_record(
    input_record: _InputRecord,
    pos_1based: int,
    ref: str,
    alts: Sequence[str],
    qual: str,
    sample_fields: Sequence[str],
) -> _OutputRecord:
    """A rewritten record: CHROM, POS, REF, ALT, QUAL and each sample's GT and PS, and `.`
    for ID, FILTER and INFO; no FORMAT column when the file has no samples."""
    columns = [input_record.contig, str(pos_1based), ".", ref, ",".join(alts), qual, ".", "."]
    if sample_fields:
        columns += ["GT:PS", *sample_fields]
    return _OutputRecord(
        contig=input_record.contig,
        pos_1based=pos_1based,
        line="\t".join(columns),
        is_rewritten=True,
    )


def _write_output(
    output: Path, respelling_input: _RespellingInput, output_records: Sequence[_OutputRecord]
) -> None:
    """The header and the records, sorted by contig in the input's order and by POS; records
    of one POS keep the order they are given in, which is the input's."""
    contig_ranks = _make_contig_ranks(record.contig for record in respelling_input.input_records)
    sorted_records = sorted(
        output_records, key=lambda record: (contig_ranks[record.contig], record.pos_1based)
    )
    header_lines = _make_header(
        respelling_input, [record for record in output_records if record.is_rewritten]
    )
    lines = [*header_lines, *(record.line for record in sorted_records)]
    data = "".join(f"{line}\n" for line in lines).encode("utf-8")
    output.write_bytes(make_bgzf(data) if respelling_input.is_gzipped else data)


def _make_contig_ranks(contigs: Iterable[str]) -> Mapping[str, int]:
    return {contig: rank for rank, contig in enumerate(dict.fromkeys(contigs))}


def _make_header(
    respelling_input: _RespellingInput, rewritten_records: Sequence[_OutputRecord]
) -> list[str]:
    """The input's header, with the FORMAT lines of GT and PS before #CHROM where the
    rewritten records need them (§2.3).

    GT, when a rewritten record has sample columns, since it writes GT whatever the input's
    FORMAT holds. PS, when every PS of the input and of the rewritten records is a whole
    number of at most _MAX_INTEGER_PHASE_SET: bcftools 1.24 stops on a PS such as `chr1_5`
    declared Integer, warns on one declared String, and sets a larger number to missing; an
    input with such a PS and no PS line already gets its warning.
    """
    header_lines = list(respelling_input.header_lines)
    rewritten_fields = [record.line.split("\t") for record in rewritten_records]
    added_lines = []
    has_samples = any(len(fields) > _NUM_FIXED_COLUMNS for fields in rewritten_fields)
    if has_samples and not any(line.startswith(_GT_LINE_START) for line in header_lines):
        added_lines.append(_GT_LINE)
    written_phase_sets = (
        phase_set for fields in rewritten_fields for phase_set in _find_phase_sets(fields)
    )
    if (
        rewritten_records
        and not any(line.startswith(_PS_LINE_START) for line in header_lines)
        and respelling_input.are_phase_sets_integer
        and all(map(_is_integer_phase_set, written_phase_sets))
    ):
        added_lines.append(_PS_LINE)
    chrom_index = next(
        index for index, line in enumerate(header_lines) if line.startswith("#CHROM")
    )
    header_lines[chrom_index:chrom_index] = added_lines
    return header_lines


def _is_integer_phase_set(phase_set: str) -> bool:
    """Whether a PS is a whole number that bcftools 1.24 keeps in an Integer."""
    number = _parse_phase_set_number(phase_set)
    return number is not None and number <= _MAX_INTEGER_PHASE_SET


# BGZF, the SAM/BAM format specification §4.1: gzip members of at most 64 KiB, each with the
# extra field BC that holds its size, and an empty member at the end. bgzip puts 0xff00
# bytes of text in each, so that a member of text that does not compress still fits.
_BGZF_TEXT_BYTES = 0xFF00
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
# ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN, SI1, SI2, SLEN and BSIZE, the block's size less 1.
_BGZF_HEADER = struct.Struct("<4BI2BH2BHH")
_BGZF_FOOTER = struct.Struct("<II")  # CRC32 and ISIZE, the size of the text


def make_bgzf(data: bytes) -> bytes:
    """data as a BGZF file."""
    blocks = []
    for offset in range(0, len(data), _BGZF_TEXT_BYTES):
        text = data[offset : offset + _BGZF_TEXT_BYTES]
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)  # raw deflate, as gzip holds it
        deflated = compressor.compress(text) + compressor.flush()
        block_size = _BGZF_HEADER.size + len(deflated) + _BGZF_FOOTER.size
        header = _BGZF_HEADER.pack(31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, block_size - 1)
        footer = _BGZF_FOOTER.pack(zlib.crc32(text), len(text))
        blocks.append(header + deflated + footer)
    blocks.append(_BGZF_EOF)
    return b"".join(blocks)
