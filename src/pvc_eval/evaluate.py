"""Evaluation, stage 5 of doc/spec.md §4, as doc/specs/reconciliation.md §5 and §6 give it:
the outcome of each sample pair at a reconciled locus, the pair of possible genotypes that
decided it, and its row of the table.

Public: evaluate_pair, the outcome of two sample results; PairEvaluation, what it gives;
make_outcome_rows, the rows of a reconciled locus.
"""

from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass

from pvc_eval.kind import Kind
from pvc_eval.outcome import InvalidReason, LocusGenotype, Outcome, OutcomeRow
from pvc_eval.reconcile import TOO_COMPLEX, ReconciledLocus, SampleResult
from pvc_eval.samples import SamplePair
from pvc_eval.vcf import MISSING, VcfGenotype, VcfRecord


@dataclass(frozen=True, slots=True)
class PairEvaluation:
    """The outcome of a sample pair at a locus and the pair of possible genotypes it
    reports (doc/specs/reconciliation.md §5); both genotypes are None when the outcome
    reports none: partly outside the regions, too complex and invalid."""

    outcome: Outcome
    genotype_a: LocusGenotype | None
    genotype_b: LocusGenotype | None


# The rows of a locus partly outside the regions: the outcome comes from the locus, not from
# the sample results, and reports no pair (doc/specs/reconciliation.md §5).
_PARTLY_OUTSIDE = PairEvaluation(Outcome.PARTLY_OUTSIDE, None, None)


def evaluate_pair(result_a: SampleResult, result_b: SampleResult) -> PairEvaluation:
    """The outcome of the first rule of doc/specs/reconciliation.md §5 that applies to two
    sample results, file A's and file B's, and the pair it reports."""
    if result_a is TOO_COMPLEX or result_b is TOO_COMPLEX:
        return PairEvaluation(Outcome.TOO_COMPLEX, None, None)
    if isinstance(result_a, InvalidReason) or isinstance(result_b, InvalidReason):
        return PairEvaluation(_get_invalid_outcome(result_a, result_b), None, None)
    first_a, first_b = result_a[0], result_b[0]
    missing_outcome = _get_missing_outcome(result_a, result_b)
    if missing_outcome is not None:
        return PairEvaluation(missing_outcome, first_a, first_b)
    is_reference_a = len(result_a) == 1 and _is_all_ref(first_a)
    is_reference_b = len(result_b) == 1 and _is_all_ref(first_b)
    if is_reference_a and is_reference_b:
        return PairEvaluation(Outcome.REFERENCE_IN_BOTH, first_a, first_b)
    same_genotype = _find_first_pair(result_a, result_b, _is_same_genotype)
    if same_genotype is not None:
        return PairEvaluation(Outcome.SAME_GENOTYPE, *same_genotype)
    same_alleles = _find_first_pair(result_a, result_b, _is_same_alt_set)
    if same_alleles is not None:
        return PairEvaluation(Outcome.SAME_ALLELES, *same_alleles)
    if is_reference_b:
        return PairEvaluation(Outcome.ONLY_IN_A, first_a, first_b)
    if is_reference_a:
        return PairEvaluation(Outcome.ONLY_IN_B, first_a, first_b)
    return PairEvaluation(Outcome.DIFFERENT, first_a, first_b)


def make_outcome_rows(
    reconciled: ReconciledLocus,
    sample_pairs: Sequence[SamplePair],
    *,
    kind: Kind,
    stratum: str | None,
    partly_outside: bool = False,
) -> Iterator[OutcomeRow]:
    """The row of each sample pair at a reconciled locus, in the order of sample_pairs
    (doc/specs/reconciliation.md §6).

    kind: the locus's kind, locus_kind's.
    stratum: the locus's stratum as doc/specs/loci.md gives it, `.` included; None
        without --strata.
    partly_outside: whether --regions is given and the locus has a position outside them;
        every row is then partly outside the regions, whatever evaluate_pair gives, and
        reports no pair (doc/spec.md §5.5).
    """
    locus = reconciled.locus
    records_a = _count_records_with_alt(locus.records_a)
    records_b = _count_records_with_alt(locus.records_b)
    for sample_pair in sample_pairs:
        result_a = reconciled.results_a[sample_pair.index_a]
        result_b = reconciled.results_b[sample_pair.index_b]
        evaluation = _PARTLY_OUTSIDE if partly_outside else evaluate_pair(result_a, result_b)
        quality_a, gq_a = _find_lowest_qual_and_gq(locus.records_a, sample_pair.index_a)
        quality_b, gq_b = _find_lowest_qual_and_gq(locus.records_b, sample_pair.index_b)
        yield OutcomeRow(
            contig=locus.contig,
            start=locus.start,
            end=locus.end,
            kind=kind,
            stratum=stratum,
            str_regions=locus.str_regions,
            sample_a=sample_pair.name_a,
            sample_b=sample_pair.name_b,
            outcome=evaluation.outcome,
            alleles=reconciled.alleles,
            genotype_a=_get_shown_genotype(evaluation.genotype_a, result_a),
            genotype_b=_get_shown_genotype(evaluation.genotype_b, result_b),
            possible_a=_count_possible_genotypes(result_a),
            possible_b=_count_possible_genotypes(result_b),
            quality_a=quality_a,
            quality_b=quality_b,
            gq_a=gq_a,
            gq_b=gq_b,
            records_a=records_a,
            records_b=records_b,
            invalid_reason_a=result_a if isinstance(result_a, InvalidReason) else None,
            invalid_reason_b=result_b if isinstance(result_b, InvalidReason) else None,
        )


def _get_invalid_outcome(result_a: SampleResult, result_b: SampleResult) -> Outcome:
    is_invalid_a = isinstance(result_a, InvalidReason)
    is_invalid_b = isinstance(result_b, InvalidReason)
    if is_invalid_a and is_invalid_b:
        return Outcome.INVALID_IN_BOTH
    return Outcome.INVALID_IN_A if is_invalid_a else Outcome.INVALID_IN_B


def _get_missing_outcome(
    possible_a: Sequence[LocusGenotype], possible_b: Sequence[LocusGenotype]
) -> Outcome | None:
    """The missing outcome of two samples' possible genotypes, or None when neither sample
    is missing: a sample is when every one of its genotypes has a MISSING."""
    is_missing_a = all(_is_partly_missing(genotype) for genotype in possible_a)
    is_missing_b = all(_is_partly_missing(genotype) for genotype in possible_b)
    if is_missing_a and is_missing_b:
        return Outcome.MISSING_IN_BOTH
    if is_missing_a:
        return Outcome.MISSING_IN_A
    if is_missing_b:
        return Outcome.MISSING_IN_B
    return None


def _find_first_pair(
    possible_a: Sequence[LocusGenotype],
    possible_b: Sequence[LocusGenotype],
    is_match: Callable[[LocusGenotype, LocusGenotype], bool],
) -> tuple[LocusGenotype, LocusGenotype] | None:
    """The first pair of possible genotypes that match, A's in order and, for each, B's in
    order (§5), or None when no pair does."""
    for genotype_a in possible_a:
        for genotype_b in possible_b:
            if is_match(genotype_a, genotype_b):
                return genotype_a, genotype_b
    return None


def _is_partly_missing(genotype: LocusGenotype) -> bool:
    return any(allele is MISSING for allele in genotype)


def _is_all_ref(genotype: LocusGenotype) -> bool:
    return all(allele == 0 for allele in genotype)


def _is_same_genotype(genotype_a: LocusGenotype, genotype_b: LocusGenotype) -> bool:
    """Two genotypes are the same when equal, and two genotypes all REF are the same
    whatever their ploidies (decided by the owner on 19 October 2026)."""
    return genotype_a == genotype_b or (_is_all_ref(genotype_a) and _is_all_ref(genotype_b))


def _is_same_alt_set(genotype_a: LocusGenotype, genotype_b: LocusGenotype) -> bool:
    return _get_alt_set(genotype_a) == _get_alt_set(genotype_b)


def _get_alt_set(genotype: LocusGenotype) -> set[int]:
    """The indices of a genotype other than 0; MISSING is no index."""
    return {allele for allele in genotype if allele is not MISSING and allele != 0}


def _get_shown_genotype(
    reported: LocusGenotype | None, sample_result: SampleResult
) -> LocusGenotype | None:
    """The genotype a row shows for a sample: the reported one, or, in a row that reports
    none, the first possible genotype of a sample that has one (§5)."""
    if reported is not None:
        return reported
    if sample_result is TOO_COMPLEX or isinstance(sample_result, InvalidReason):
        return None
    return sample_result[0]


def _count_possible_genotypes(sample_result: SampleResult) -> int | None:
    """How many possible genotypes a sample has: 0 when invalid, None when too complex."""
    if sample_result is TOO_COMPLEX:
        return None
    if isinstance(sample_result, InvalidReason):
        return 0
    return len(sample_result)


def _count_records_with_alt(records: Iterable[VcfRecord]) -> int:
    return sum(1 for vcf_record in records if len(vcf_record.alleles) > 1)


def _find_lowest_qual_and_gq(
    records: Iterable[VcfRecord], sample_index: int
) -> tuple[float | None, int | None]:
    """The lowest QUAL and the lowest GQ among the records at which the sample's VCF
    genotype has an index other than 0, a `*` included (§6), each None when no such record
    has one. A MALFORMED_GT or OTHER_PLOIDY_GT has no alleles and is passed over."""
    quals: list[float] = []
    gqs: list[int] = []
    for vcf_record in records:
        vcf_genotype = vcf_record.genotypes[sample_index]
        if not isinstance(vcf_genotype, VcfGenotype):
            continue
        if not any(allele is not MISSING and allele != 0 for allele in vcf_genotype.alleles):
            continue
        if vcf_record.qual is not None:
            quals.append(vcf_record.qual)
        if vcf_genotype.quality is not None:
            gqs.append(vcf_genotype.quality)
    return min(quals, default=None), min(gqs, default=None)
