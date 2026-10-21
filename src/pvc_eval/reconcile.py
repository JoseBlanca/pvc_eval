"""Reconciliation, stage 4 of doc/spec.md §4, as doc/specs/reconciliation.md §4 gives it:
for each locus, the sequence each haplotype of each paired sample carries over it, in both
files and under every phase assignment the VCF genotypes leave open, named by one shared
list of locus alleles.

Positions are 0-based and half-open, as in Locus and VcfRecord.

Public: reconcile, which reconciles one locus; reconcile_loci, which reconciles a stream of
them; ReconciledLocus, what both give; SampleResult, what one sample gets; TOO_COMPLEX, the
sample result of a sample over the phasing cap; and PHASING_CAP, the default cap.
"""

from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass

from pvc_eval.errors import format_count
from pvc_eval.fasta import FastaWindow, IndexedFasta
from pvc_eval.locus import Locus
from pvc_eval.outcome import InvalidReason, LocusGenotype
from pvc_eval.samples import SamplePair
from pvc_eval.vcf import (
    DEFAULT_PLOIDY,
    MALFORMED_GT,
    MISSING,
    OTHER_PLOIDY_GT,
    Span,
    VcfGenotype,
    VcfRecord,
    trim_alt,
)

PHASING_CAP = 4096  # doc/spec.md §5.4; kept by the owner on 19 October 2026 (Open 2)

TOO_COMPLEX = sentinel("TOO_COMPLEX")
"""The sample result of a sample whose phase assignments would exceed the phasing cap."""

type SampleResult = tuple[LocusGenotype, ...] | InvalidReason | TOO_COMPLEX
"""What reconciliation gives one sample of one file at one locus: its possible genotypes,
one or more, in the order of doc/specs/reconciliation.md §4.3; why it is invalid; or
TOO_COMPLEX."""

_SPANNING_DELETION = "*"
_PLOIDIES = (1, 2)  # doc/specs/reconciliation.md §2: reconciliation handles these alone

# The haplotypes of one phase assignment, each with one VCF allele per record of the
# locus, in the records' order: an index into the record's alleles, or MISSING.
type _Haplotype = tuple[int | MISSING, ...]
# A locus genotype before its locus alleles are given indices: one locus allele per
# haplotype, in the haplotypes' order, or MISSING.
type _SequenceGenotype = tuple[str | MISSING, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class ReconciledLocus:
    """One locus's reconciliation.

    alleles: the locus alleles, the reference over the locus first.
    results_a, results_b: the sample result of each paired sample of file A, and of file
        B, by the sample's index among its file's sample columns.
    """

    locus: Locus
    alleles: tuple[str, ...]
    results_a: Mapping[int, SampleResult]
    results_b: Mapping[int, SampleResult]


@dataclass(frozen=True, slots=True, kw_only=True)
class _Edit:
    """What one VCF allele does to the reference: replace `[start, end)`, 0-based and
    half-open, with bases. An insertion has end equal to start, and goes in before start;
    a deletion has no bases.

    covered: the positions on which the edit covers a `*`, 0-based and half-open (§4.2):
        for a deletion, those its placements delete; for another edit that replaces
        positions, those; None for an insertion, which replaces none.
    """

    start: int
    end: int
    bases: str
    covered: Span | None

    @property
    def is_insertion(self) -> bool:
        return self.start == self.end


# The edit of each VCF allele of a record, by the allele's index: None for REF, for an
# allele that spells REF and for `*`.
type _RecordEdits = tuple[_Edit | None, ...]


def reconcile(
    locus: Locus,
    locus_bases: str,
    sample_pairs: Sequence[SamplePair],
    *,
    ploidy_a: int = DEFAULT_PLOIDY,
    ploidy_b: int = DEFAULT_PLOIDY,
    phasing_cap: int = PHASING_CAP,
) -> ReconciledLocus:
    """The locus alleles of a locus and the sample result of each paired sample of each
    file (doc/specs/reconciliation.md §4). A sample in no pair is not reconciled.

    locus_bases: the reference over the locus, in upper case.
    ploidy_a, ploidy_b: the ploidy each file was read with, 1 or 2.
    phasing_cap: the most phase assignments a sample may have, 1 or more.

    Raises ValueError for a ploidy other than 1 or 2, a phasing cap under 1, bases of
    another length than the locus, and a VCF genotype of a paired sample with another
    number of alleles than its file's ploidy, which means its file was read with another
    ploidy, and for a sample pair whose index is not one of a record's genotypes.
    """
    _check_settings(ploidy_a, ploidy_b, phasing_cap)
    if len(locus_bases) != locus.end - locus.start:
        raise ValueError(
            f"locus_bases has {len(locus_bases)} bases, and the locus, {locus.start}-{locus.end},"
            f" has {locus.end - locus.start}"
        )
    _check_sample_indices(locus, sample_pairs)
    # The locus alleles by their index, the reference first (§4.3); a dict keeps them in
    # order of first appearance, file A's samples first.
    locus_alleles = {locus_bases: 0}
    results_a = _reconcile_file(
        "A",
        locus.records_a,
        sorted({pair.index_a for pair in sample_pairs}),
        ploidy_a,
        locus.start,
        locus_bases,
        phasing_cap,
        locus_alleles,
    )
    results_b = _reconcile_file(
        "B",
        locus.records_b,
        sorted({pair.index_b for pair in sample_pairs}),
        ploidy_b,
        locus.start,
        locus_bases,
        phasing_cap,
        locus_alleles,
    )
    return ReconciledLocus(
        locus=locus, alleles=tuple(locus_alleles), results_a=results_a, results_b=results_b
    )


def reconcile_loci(
    loci: Iterator[Locus],
    reference: IndexedFasta,
    sample_pairs: Sequence[SamplePair],
    *,
    ploidy_a: int = DEFAULT_PLOIDY,
    ploidy_b: int = DEFAULT_PLOIDY,
    phasing_cap: int = PHASING_CAP,
) -> Iterator[ReconciledLocus]:
    """Each locus reconciled, as reconcile gives it, with its bases read from the
    reference; stage 4 of doc/spec.md §4.

    loci: in the order build_loci gives them, contig by contig and by start.

    Raises what reconcile raises; for its settings, when the stream starts.
    """
    _check_settings(ploidy_a, ploidy_b, phasing_cap)
    # A window of the stream's own, as each thread would have.
    reference_window = FastaWindow(reference)
    for locus in loci:
        # The loci of a contig come by start and do not overlap, so no later one needs
        # the bases before this one.
        reference_window.evict_before(locus.start)
        locus_bases = reference_window.get_bases(locus.contig, locus.start, locus.end)
        yield reconcile(
            locus,
            locus_bases,
            sample_pairs,
            ploidy_a=ploidy_a,
            ploidy_b=ploidy_b,
            phasing_cap=phasing_cap,
        )


def _reconcile_file(
    file_name: str,
    records: Sequence[VcfRecord],
    sample_indices: Sequence[int],
    ploidy: int,
    locus_start: int,
    locus_bases: str,
    phasing_cap: int,
    locus_alleles: dict[str, int],
) -> dict[int, SampleResult]:
    """The sample result of each of a file's paired samples, by its index, in the order of
    sample_indices; the locus alleles they carry are added to locus_alleles."""
    _check_genotype_ploidy(records, sample_indices, ploidy, file_name)
    # Made once for the locus, and read by every haplotype of every sample.
    edits = [_make_record_edits(vcf_record, locus_start, locus_bases) for vcf_record in records]
    return {
        sample_index: _make_sample_result(
            locus_start,
            locus_bases,
            records,
            edits,
            sample_index,
            ploidy,
            phasing_cap,
            locus_alleles,
        )
        for sample_index in sample_indices
    }


def _check_settings(ploidy_a: int, ploidy_b: int, phasing_cap: int) -> None:
    for name, ploidy in (("ploidy_a", ploidy_a), ("ploidy_b", ploidy_b)):
        if not _is_whole(ploidy) or ploidy not in _PLOIDIES:
            raise ValueError(f"{name} must be 1 or 2, not {ploidy!r}")
    if not _is_whole(phasing_cap) or phasing_cap < 1:
        raise ValueError(f"phasing_cap must be a whole number of 1 or more, not {phasing_cap!r}")


def _is_whole(number: int) -> bool:
    return isinstance(number, int) and not isinstance(number, bool)


def _check_sample_indices(locus: Locus, sample_pairs: Iterable[SamplePair]) -> None:
    """Raise ValueError when a sample pair names a sample beyond the genotypes of a record
    of its file."""
    for sample_pair in sample_pairs:
        for file_name, records, sample_index in (
            ("A", locus.records_a, sample_pair.index_a),
            ("B", locus.records_b, sample_pair.index_b),
        ):
            for vcf_record in records:
                if sample_index >= len(vcf_record.genotypes):
                    raise ValueError(
                        f"the sample pair {sample_pair.name_a} of file A and"
                        f" {sample_pair.name_b} of file B names sample {sample_index} of file"
                        f" {file_name}, and its VCF record at line {vcf_record.line_number} has"
                        f" genotypes for {format_count(len(vcf_record.genotypes), 'sample')}"
                    )


def _check_genotype_ploidy(
    records: Iterable[VcfRecord], sample_indices: Iterable[int], ploidy: int, file_name: str
) -> None:
    """Raise ValueError when a sample's VCF genotype has another number of alleles than
    the ploidy: the file was read with another ploidy."""
    for vcf_record in records:
        for sample_index in sample_indices:
            vcf_genotype = vcf_record.genotypes[sample_index]
            if isinstance(vcf_genotype, VcfGenotype) and len(vcf_genotype.alleles) != ploidy:
                raise ValueError(
                    f"file {file_name}'s VCF record at line {vcf_record.line_number} has a VCF"
                    f" genotype of {len(vcf_genotype.alleles)} alleles, and file {file_name} is"
                    f" reconciled with ploidy {ploidy}: read it with that ploidy"
                )


def _make_sample_result(
    locus_start: int,
    locus_bases: str,
    records: Sequence[VcfRecord],
    edits: Sequence[_RecordEdits],
    sample_index: int,
    ploidy: int,
    phasing_cap: int,
    locus_alleles: dict[str, int],
) -> SampleResult:
    """A sample's sample result: the first reason of doc/specs/reconciliation.md §4 that
    applies, or its possible genotypes, whose locus alleles are added to locus_alleles."""
    genotypes = [vcf_record.genotypes[sample_index] for vcf_record in records]
    if any(genotype is MALFORMED_GT for genotype in genotypes):
        return InvalidReason.MALFORMED_GT
    vcf_genotypes = [genotype for genotype in genotypes if isinstance(genotype, VcfGenotype)]
    if len(vcf_genotypes) != len(genotypes):
        assert any(genotype is OTHER_PLOIDY_GT for genotype in genotypes)
        return InvalidReason.PLOIDY_MISMATCH
    block_of_record = _find_blocks(vcf_genotypes)
    if _count_assignments(block_of_record) > phasing_cap:
        return TOO_COMPLEX
    locus_genotypes: list[_SequenceGenotype] = []
    for haplotypes in _make_assignments(vcf_genotypes, block_of_record, ploidy):
        locus_genotype = _make_locus_genotype(locus_start, locus_bases, records, edits, haplotypes)
        if locus_genotype is not None:
            locus_genotypes.append(locus_genotype)
    if not locus_genotypes:
        return InvalidReason.CONTRADICTING_RECORDS
    return _make_possible_genotypes(locus_genotypes, locus_alleles)


def _find_blocks(vcf_genotypes: Sequence[VcfGenotype]) -> list[int | None]:
    """The block of each VCF genotype, numbered from 0 in the order of the first record of
    each block, or None for a genotype of a phase unit that is no block (§4.1).

    A phased genotype is in the phase unit of its phase set, None included, which is the
    one phase set of every phased genotype with no PS; any other is a unit of its own. A
    unit is a block when one of its genotypes has two different alleles, `.` counting as
    an allele.
    """
    unit_of_record: list[tuple[str, str | None] | tuple[str, int]] = [
        ("phase set", vcf_genotype.phase_set) if vcf_genotype.phased else ("record", place)
        for place, vcf_genotype in enumerate(vcf_genotypes)
    ]
    heterozygous_units = {
        unit
        for unit, vcf_genotype in zip(unit_of_record, vcf_genotypes, strict=True)
        if len(set(vcf_genotype.alleles)) > 1
    }
    # dict.fromkeys keeps the units in the order of their first record.
    blocks = [unit for unit in dict.fromkeys(unit_of_record) if unit in heterozygous_units]
    block_numbers = {unit: number for number, unit in enumerate(blocks)}
    return [block_numbers.get(unit) for unit in unit_of_record]


def _count_assignments(block_of_record: Sequence[int | None]) -> int:
    """2^max(b - 1, 0) for b blocks (§4.1). With ploidy 1 no genotype has two different
    alleles, so there is no block and one assignment."""
    num_blocks = len({block for block in block_of_record if block is not None})
    return 2 ** max(num_blocks - 1, 0)


def _make_assignments(
    vcf_genotypes: Sequence[VcfGenotype], block_of_record: Sequence[int | None], ploidy: int
) -> Iterator[tuple[_Haplotype, ...]]:
    """The haplotypes of each phase assignment, in the order of §4.1: assignment k swaps
    the two haplotypes of block j, for j of 1 or more, when bit j - 1 of k is set."""
    for assignment in range(_count_assignments(block_of_record)):
        haplotypes = []
        for haplotype in range(ploidy):
            alleles = []
            for vcf_genotype, block in zip(vcf_genotypes, block_of_record, strict=True):
                is_swapped = (
                    block is not None and block >= 1 and bool(assignment >> (block - 1) & 1)
                )
                alleles.append(vcf_genotype.alleles[1 - haplotype if is_swapped else haplotype])
            haplotypes.append(tuple(alleles))
        yield tuple(haplotypes)


def _make_locus_genotype(
    locus_start: int,
    locus_bases: str,
    records: Sequence[VcfRecord],
    edits: Sequence[_RecordEdits],
    haplotypes: Iterable[_Haplotype],
) -> _SequenceGenotype | None:
    """The locus allele of each haplotype of an assignment, or None when one of them is
    impossible."""
    locus_genotype = []
    for haplotype in haplotypes:
        locus_allele = _make_locus_allele(locus_start, locus_bases, records, edits, haplotype)
        if locus_allele is None:
            return None
        locus_genotype.append(locus_allele)
    return tuple(locus_genotype)


def _make_locus_allele(
    locus_start: int,
    locus_bases: str,
    records: Sequence[VcfRecord],
    edits: Sequence[_RecordEdits],
    haplotype: _Haplotype,
) -> str | MISSING | None:
    """A haplotype's locus allele (§4.2): None when it is impossible, MISSING when it is
    missing at a record, and otherwise the reference over the locus with its edits."""
    carried_edits = []
    star_records = []
    for vcf_record, record_edits, allele_index in zip(records, edits, haplotype, strict=True):
        if allele_index is MISSING:
            continue
        if vcf_record.alleles[allele_index] == _SPANNING_DELETION:
            star_records.append(vcf_record)
            continue
        edit = record_edits[allele_index]
        if edit is not None:
            carried_edits.append(edit)
    if _is_colliding(carried_edits):
        return None
    for star_record in star_records:
        if not _is_star_covered(star_record, records, edits, haplotype):
            return None
    if any(allele_index is MISSING for allele_index in haplotype):
        return MISSING
    return _make_edited_bases(locus_start, locus_bases, carried_edits)


def _make_record_edits(vcf_record: VcfRecord, locus_start: int, locus_bases: str) -> _RecordEdits:
    """The edit of each VCF allele of a record, by its index (see _RecordEdits)."""
    return tuple(
        None
        if allele_index == 0 or allele == _SPANNING_DELETION
        else _make_edit(vcf_record, allele_index, locus_start, locus_bases)
        for allele_index, allele in enumerate(vcf_record.alleles)
    )


def _make_edit(
    vcf_record: VcfRecord, allele_index: int, locus_start: int, locus_bases: str
) -> _Edit | None:
    """The edit of a VCF allele that is not `*`, or None when it spells REF: what is left
    of REF and of the allele once the bases they share are removed from the left and then
    from the right (§4.2)."""
    num_prefix, ref, alt = trim_alt(vcf_record.alleles[0], vcf_record.alleles[allele_index])
    if not ref and not alt:
        return None
    start = vcf_record.start + num_prefix
    end = start + len(ref)
    covered: Span | None = (start, end)
    if not ref:
        covered = None
    elif not alt:
        covered = _find_placements(start, end, vcf_record.span, locus_start, locus_bases)
    return _Edit(start=start, end=end, bases=alt, covered=covered)


def _is_colliding(edits: Sequence[_Edit]) -> bool:
    """Whether any two of a haplotype's edits collide."""
    return any(
        _is_collision(edit, other_edit)
        for index, edit in enumerate(edits)
        for other_edit in edits[index + 1 :]
    )


def _is_collision(first: _Edit, second: _Edit) -> bool:
    """Two edits collide when they replace a position in common, when one is an insertion
    and the other replaces the positions on both sides of it, or when both are insertions
    at one point (§4.2)."""
    if first.is_insertion and second.is_insertion:
        return first.start == second.start
    if first.is_insertion:
        return second.start < first.start < second.end
    if second.is_insertion:
        return first.start < second.start < first.end
    return first.start < second.end and second.start < first.end


def _is_star_covered(
    star_record: VcfRecord,
    records: Sequence[VcfRecord],
    edits: Sequence[_RecordEdits],
    haplotype: _Haplotype,
) -> bool:
    """Whether a `*` the haplotype carries is covered (§4.2): a deletion covers it when one
    of its placements deletes a position of the `*` record's REF as the file wrote it, and
    any other edit when, as written, it replaces one. The edits that count are those of the
    alleles the haplotype carries at the other records, and every ALT of a record where it
    is MISSING, since it might carry that one."""
    # §7: a record built without written_ref is measured against its REF.
    star_start, star_end = (
        (star_record.start, star_record.end)
        if star_record.written_ref is None
        else star_record.written_ref
    )
    for vcf_record, record_edits, allele_index in zip(records, edits, haplotype, strict=True):
        if vcf_record is star_record:
            continue
        candidates = record_edits[1:] if allele_index is MISSING else (record_edits[allele_index],)
        for edit in candidates:
            if edit is None or edit.covered is None:
                continue
            start, end = edit.covered
            if start < star_end and star_start < end:
                return True
    return False


def _find_placements(start: int, end: int, span: Span, locus_start: int, locus_bases: str) -> Span:
    """The positions the placements of the deletion of `[start, end)` delete, together,
    0-based and half-open: where it was written and every slide along the reference that deletes only
    positions of its record's span (§4.2). Each placement is a stretch as long as the
    deletion, and the placements are consecutive, so together they are one stretch.

    The span lies within the locus, as a record with an ALT's does (doc/spec.md §4).
    """
    span_start, span_end = span
    deleted_start, deleted_end = start, end
    # A slide left by one base deletes the base before and restores the last one, which
    # gives the same haplotype when the two are equal; likewise to the right.
    while start > span_start and (
        locus_bases[start - 1 - locus_start] == locus_bases[end - 1 - locus_start]
    ):
        start, end = start - 1, end - 1
    first_start = start
    start, end = deleted_start, deleted_end
    while end < span_end and locus_bases[end - locus_start] == locus_bases[start - locus_start]:
        start, end = start + 1, end + 1
    return first_start, end


def _make_edited_bases(locus_start: int, locus_bases: str, edits: Iterable[_Edit]) -> str:
    """The reference over the locus with edits that do not collide, applied in position
    order, an insertion before a replacement that starts at its point."""
    pieces = []
    position = locus_start
    for edit in sorted(edits, key=lambda edit: (edit.start, not edit.is_insertion)):
        pieces.append(locus_bases[position - locus_start : edit.start - locus_start])
        pieces.append(edit.bases)
        position = edit.end
    pieces.append(locus_bases[position - locus_start :])
    return "".join(pieces)


def _make_possible_genotypes(
    locus_genotypes: Iterable[_SequenceGenotype], locus_alleles: dict[str, int]
) -> tuple[LocusGenotype, ...]:
    """A sample's possible genotypes (§4.3): each locus genotype as the sorted indices of
    its locus alleles, MISSING last, without repeats, sorted.

    locus_alleles: the locus alleles met so far, by sequence; a locus allele met for the
        first time is added with the next index, so the order of the genotypes and of the
        haplotypes within each sets the indices.
    """
    possible_genotypes = set()
    for locus_genotype in locus_genotypes:
        indices: list[int | MISSING] = []
        for locus_allele in locus_genotype:
            if locus_allele is MISSING:
                indices.append(MISSING)
            else:
                indices.append(locus_alleles.setdefault(locus_allele, len(locus_alleles)))
        possible_genotypes.add(tuple(sorted(indices, key=_get_sort_key)))
    return tuple(sorted(possible_genotypes, key=lambda genotype: [*map(_get_sort_key, genotype)]))


def _get_sort_key(allele: int | MISSING) -> tuple[bool, int]:
    """The key that sorts a locus allele's index, MISSING after every index."""
    if allele is MISSING:
        return True, 0
    return False, allele
