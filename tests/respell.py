"""The re-spelling program of doc/specs/respell.md: a VCF file written again in other words
that give every sample the same haplotypes, for the check of doc/spec.md §7.4.

What every re-spelling shares is here (§2): the input read as lines and as records, the
records every re-spelling keeps, the touched test, the reach, and the writing of the output,
with its header, its order, its new phase sets and, for a gzipped input, BGZF.

Positions are 0-based and half-open, as everywhere inside pvc_eval, but for POS, which is
1-based where a name says so.

Public: respell_by_decomposing (§3).
"""

import bisect
import gzip
import struct
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
# The columns of a VCF line, from 0; the samples' start after FORMAT.
_CHROM_COLUMN = 0
_POS_COLUMN = 1
_REF_COLUMN = 3
_QUAL_COLUMN = 5
_FORMAT_COLUMN = 8
_NUM_FIXED_COLUMNS = 9
_PS_LINE_START = "##FORMAT=<ID=PS,"
_PS_LINE = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'


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

    max_phase_set: the largest PS of the input that is a whole number, 0 when there is none.
    has_text_phase_set: some PS of the input is not a whole number.
    """

    header_lines: tuple[str, ...]
    input_records: tuple[_InputRecord, ...]
    is_gzipped: bool
    max_phase_set: int
    has_text_phase_set: bool


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
    touched_places = _find_touched_places(respelling_input.input_records)
    spanning_deletion_refs = _find_spanning_deletion_refs(respelling_input.input_records)
    output_records = []
    for input_record in respelling_input.input_records:
        vcf_record = input_record.vcf_record
        if vcf_record is None or _is_kept_whatever_is_written(
            vcf_record, input_record.place, touched_places
        ):
            output_records.append(_make_kept_record(input_record))
            continue
        pieces = _make_pieces(vcf_record, reference)
        written = [
            (vcf_record.start, vcf_record.alleles),
            *((piece.start, (piece.ref, piece.alt)) for piece in pieces),
        ]
        if not pieces or _is_near_spanning_deletion(
            vcf_record.contig, written, spanning_deletion_refs, reference
        ):
            output_records.append(_make_kept_record(input_record))
            continue
        new_phase_set = str(respelling_input.max_phase_set + input_record.place)
        output_records.extend(
            _make_piece_record(input_record, piece, allele_indices, new_phase_set)
            for piece, allele_indices in pieces.items()
        )
    _write_output(output, respelling_input, output_records)


def _is_near_spanning_deletion(
    contig: str,
    written: Iterable[tuple[int, Sequence[str]]],
    spanning_deletion_refs: Mapping[str, _SpanningDeletionRefs],
    reference: IndexedFasta,
) -> bool:
    """§2.2: the reach of a record of contig, or of a record its rewriting would write,
    touches the REF of a record with a `*` ALT.

    written: the input record and those its rewriting would write, each as its start,
        0-based, and its alleles, trimmed as parse_vcf trims them.
    """
    # A contig with no record with a `*` ALT has no REFs to touch.
    refs = spanning_deletion_refs.get(contig)
    if refs is None:
        return False
    return any(
        _is_touching_any_ref(refs, _find_reach(reference, contig, start, alleles))
        for start, alleles in written
    )


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
    whole_phase_sets = [int(phase_set) for phase_set in phase_sets if _is_whole(phase_set)]
    return _RespellingInput(
        header_lines=tuple(header_lines),
        input_records=tuple(input_records),
        is_gzipped=is_gzipped,
        max_phase_set=max(whole_phase_sets, default=0),
        has_text_phase_set=len(whole_phase_sets) < len(phase_sets),
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


def _is_whole(text: str) -> bool:
    return text.isascii() and text.isdecimal()


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


def _is_haplotype_shareable(
    genotype: VcfGenotype, haplotype: int, other_genotype: VcfGenotype, other_haplotype: int
) -> bool:
    """§1: two VCF alleles of one sample may be on one haplotype unless both VCF genotypes
    are phased in one phase set and the alleles are on different haplotypes. A phased GT
    with no PS is in the one phase set VCF gives those of its sample and contig."""
    is_one_phase_set = (
        genotype.phased and other_genotype.phased and genotype.phase_set == other_genotype.phase_set
    )
    return not is_one_phase_set or haplotype == other_haplotype


@dataclass(frozen=True, slots=True)
class _RecordEdits:
    place: int
    vcf_record: VcfRecord
    edits: tuple[_Edit | None, ...]


def _find_touched_places(input_records: Iterable[_InputRecord]) -> set[int]:
    """The places of the touched records (§2.2): those at which some sample carries an edit
    that touches an edit it carries at another record of the contig, the two VCF alleles
    being ones that may share a haplotype.

    An edit lies within its record's REF, its ends included, so the records are swept by REF
    and each is compared only with those before it whose REF reaches its start (§10).
    """
    by_contig: dict[str, list[_RecordEdits]] = {}
    for input_record in input_records:
        vcf_record = input_record.vcf_record
        if vcf_record is None or len(vcf_record.alleles) < 2:
            continue
        record_edits = _RecordEdits(input_record.place, vcf_record, _find_edits(vcf_record))
        by_contig.setdefault(vcf_record.contig, []).append(record_edits)
    touched_places: set[int] = set()
    for contig_records in by_contig.values():
        contig_records.sort(key=lambda record_edits: record_edits.vcf_record.start)
        reaching: list[_RecordEdits] = []
        for record_edits in contig_records:
            start = record_edits.vcf_record.start
            reaching = [earlier for earlier in reaching if earlier.vcf_record.end >= start]
            for earlier in reaching:
                if _is_carried_edit_touching(earlier, record_edits):
                    touched_places.update((earlier.place, record_edits.place))
            reaching.append(record_edits)
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
                    genotype, haplotype, other_genotype, other_haplotype
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


@dataclass(frozen=True, slots=True)
class _SpanningDeletionRefs:
    """The REFs of a contig's records with a `*` ALT, as parse_vcf trims them: their starts,
    sorted, and for each the furthest end of the REFs up to it."""

    starts: tuple[int, ...]
    furthest_ends: tuple[int, ...]


def _make_spanning_deletion_refs(refs: Iterable[Span]) -> _SpanningDeletionRefs:
    sorted_refs = sorted(refs)
    furthest_ends = []
    furthest_end = 0
    for _, end in sorted_refs:
        furthest_end = max(furthest_end, end)
        furthest_ends.append(furthest_end)
    return _SpanningDeletionRefs(tuple(start for start, _ in sorted_refs), tuple(furthest_ends))


def _is_touching_any_ref(refs: _SpanningDeletionRefs, stretch: Span) -> bool:
    """Whether a REF overlaps or lies next to stretch."""
    stretch_start, stretch_end = stretch
    num_starting_by_end = bisect.bisect_right(refs.starts, stretch_end)
    return num_starting_by_end > 0 and refs.furthest_ends[num_starting_by_end - 1] >= stretch_start


def _find_spanning_deletion_refs(
    input_records: Iterable[_InputRecord],
) -> dict[str, _SpanningDeletionRefs]:
    """The REFs of the records with a `*` ALT, by contig; a contig with none is absent."""
    refs: dict[str, list[Span]] = {}
    for input_record in input_records:
        vcf_record = input_record.vcf_record
        if vcf_record is not None and _SPANNING_DELETION in vcf_record.alleles[1:]:
            refs.setdefault(vcf_record.contig, []).append((vcf_record.start, vcf_record.end))
    return {
        contig: _make_spanning_deletion_refs(contig_refs) for contig, contig_refs in refs.items()
    }


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
        respelling_input, has_rewritten=any(record.is_rewritten for record in output_records)
    )
    lines = [*header_lines, *(record.line for record in sorted_records)]
    data = "".join(f"{line}\n" for line in lines).encode("utf-8")
    output.write_bytes(_make_bgzf(data) if respelling_input.is_gzipped else data)


def _make_contig_ranks(contigs: Iterable[str]) -> Mapping[str, int]:
    return {contig: rank for rank, contig in enumerate(dict.fromkeys(contigs))}


def _make_header(respelling_input: _RespellingInput, *, has_rewritten: bool) -> list[str]:
    """The input's header, with a PS line before #CHROM when a rewritten record needs one
    and every PS of the input is a whole number (§2.3): bcftools 1.24 stops on a PS such as
    `chr1_5` declared Integer and warns on one declared String, and an input with such a PS
    and no PS line already gets its warning."""
    header_lines = list(respelling_input.header_lines)
    has_ps_line = any(line.startswith(_PS_LINE_START) for line in header_lines)
    if not has_rewritten or has_ps_line or respelling_input.has_text_phase_set:
        return header_lines
    chrom_index = next(
        index for index, line in enumerate(header_lines) if line.startswith("#CHROM")
    )
    header_lines.insert(chrom_index, _PS_LINE)
    return header_lines


# BGZF, the SAM/BAM format specification §4.1: gzip members of at most 64 KiB, each with the
# extra field BC that holds its size, and an empty member at the end. bgzip puts 0xff00
# bytes of text in each, so that a member of text that does not compress still fits.
_BGZF_TEXT_BYTES = 0xFF00
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
# ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN, SI1, SI2, SLEN and BSIZE, the block's size less 1.
_BGZF_HEADER = struct.Struct("<4BI2BH2BHH")
_BGZF_FOOTER = struct.Struct("<II")  # CRC32 and ISIZE, the size of the text


def _make_bgzf(data: bytes) -> bytes:
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
