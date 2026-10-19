"""Evaluation, doc/specs/reconciliation.md §5 and §6: the cases E1 to E13 at evaluate_pair,
and the cases of the `make_outcome_rows` item of its §9.

Positions are 0-based and half-open, as in Locus and VcfRecord: the VCF records are built
here by hand, with a span that is their REF, which none of them could slide out of.
"""

from collections.abc import Sequence

import pytest

from pvc_eval.evaluate import PairEvaluation, evaluate_pair, make_outcome_rows
from pvc_eval.kind import Kind
from pvc_eval.locus import Locus
from pvc_eval.outcome import InvalidReason, Outcome, OutcomeRow
from pvc_eval.reconcile import TOO_COMPLEX, ReconciledLocus, SampleResult, reconcile
from pvc_eval.samples import SamplePair
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.vcf import MALFORMED_GT, MISSING, OTHER_PLOIDY_GT, VcfGenotype, VcfRecord

CONTIG = "chr1"
INVALID = InvalidReason.MALFORMED_GT


def make_genotype(alleles: tuple[int | MISSING, ...], gq: int | None) -> VcfGenotype:
    """An unphased VCF genotype with its GQ."""
    return VcfGenotype(alleles=alleles, phased=False, phase_set=None, quality=gq)


def make_record(
    start: int,
    alleles: tuple[str, ...],
    qual: float | None,
    genotypes: Sequence[VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT],
) -> VcfRecord:
    """A VCF record at a 0-based start, whose span is its REF."""
    return VcfRecord(
        contig=CONTIG,
        start=start,
        alleles=alleles,
        qual=qual,
        genotypes=tuple(genotypes),
        line_number=start,
        ref_mismatch=False,
        span=(start, start + len(alleles[0])),
    )


EVALUATION_CASES: dict[str, tuple[SampleResult, SampleResult, PairEvaluation]] = {
    "E1": (TOO_COMPLEX, INVALID, PairEvaluation(Outcome.TOO_COMPLEX, None, None)),
    "E2": (INVALID, ((0, 1),), PairEvaluation(Outcome.INVALID_IN_A, None, None)),
    "E3": (
        ((1, MISSING),),
        ((1, MISSING),),
        PairEvaluation(Outcome.MISSING_IN_BOTH, (1, MISSING), (1, MISSING)),
    ),
    "E4": (
        ((MISSING, MISSING),),
        ((1, 2),),
        PairEvaluation(Outcome.MISSING_IN_A, (MISSING, MISSING), (1, 2)),
    ),
    "E5": (
        ((0, 0),),
        ((0, 0),),
        PairEvaluation(Outcome.REFERENCE_IN_BOTH, (0, 0), (0, 0)),
    ),
    "E6": (
        ((0, 1), (1, 2)),
        ((1, 1), (1, 2)),
        PairEvaluation(Outcome.SAME_GENOTYPE, (1, 2), (1, 2)),
    ),
    "E7": (
        ((0, 1), (0, 2)),
        ((1, 1), (2, 2)),
        PairEvaluation(Outcome.SAME_ALLELES, (0, 1), (1, 1)),
    ),
    "E8": (((1,),), ((1, 1),), PairEvaluation(Outcome.SAME_ALLELES, (1,), (1, 1))),
    "E9": (
        ((0, 1), (0, 2)),
        ((0, 0),),
        PairEvaluation(Outcome.ONLY_IN_A, (0, 1), (0, 0)),
    ),
    "E10": (
        ((0, 0),),
        ((0, 0), (1, 2)),
        PairEvaluation(Outcome.SAME_GENOTYPE, (0, 0), (0, 0)),
    ),
    "E11": (
        ((0,),),
        ((0, 0), (1, 2)),
        PairEvaluation(Outcome.SAME_GENOTYPE, (0,), (0, 0)),
    ),
    "E12": (((0, 1),), ((2, 3),), PairEvaluation(Outcome.DIFFERENT, (0, 1), (2, 3))),
    "E13": (((0, 1),), ((1, 2),), PairEvaluation(Outcome.DIFFERENT, (0, 1), (1, 2))),
    "E14": (
        ((0, 1), (0, 2)),
        ((2, 2), (1, 1)),
        PairEvaluation(Outcome.SAME_ALLELES, (0, 1), (1, 1)),
    ),
}


@pytest.mark.parametrize(
    ("result_a", "result_b", "expected"),
    EVALUATION_CASES.values(),
    ids=EVALUATION_CASES.keys(),
)
def test_evaluate_pair_gives_outcome_and_reported_pair_of_spec_case(
    result_a: SampleResult, result_b: SampleResult, expected: PairEvaluation
) -> None:
    assert evaluate_pair(result_a, result_b) == expected


# The cases of §5 that the table of E cases leaves out: the other side of each rule with
# a side, a missing sample against an invalid one, §2's haploid 0 against 0/0, and the
# one genotype that only in A and only in B need.
OTHER_SIDE_CASES: dict[str, tuple[SampleResult, SampleResult, PairEvaluation]] = {
    "too-complex-in-B": (INVALID, TOO_COMPLEX, PairEvaluation(Outcome.TOO_COMPLEX, None, None)),
    "invalid-in-B": (((0, 1),), INVALID, PairEvaluation(Outcome.INVALID_IN_B, None, None)),
    "invalid-in-both": (
        InvalidReason.REF_MISMATCH,
        InvalidReason.CONTRADICTING_RECORDS,
        PairEvaluation(Outcome.INVALID_IN_BOTH, None, None),
    ),
    "invalid-over-missing": (
        ((1, MISSING),),
        INVALID,
        PairEvaluation(Outcome.INVALID_IN_B, None, None),
    ),
    "missing-in-B": (
        ((0, 0),),
        ((0, MISSING), (1, MISSING)),
        PairEvaluation(Outcome.MISSING_IN_B, (0, 0), (0, MISSING)),
    ),
    "only-in-B": (((0, 0),), ((1, 1),), PairEvaluation(Outcome.ONLY_IN_B, (0, 0), (1, 1))),
    "reference-in-both-across-ploidies": (
        ((0,),),
        ((0, 0),),
        PairEvaluation(Outcome.REFERENCE_IN_BOTH, (0,), (0, 0)),
    ),
    # Only in A needs B to hold one genotype, all REF: a first genotype all REF among two
    # is not enough. Likewise for only in B.
    "not-only-in-A-when-B-holds-two": (
        ((1, 1),),
        ((0, 0), (1, 2)),
        PairEvaluation(Outcome.DIFFERENT, (1, 1), (0, 0)),
    ),
    "not-only-in-B-when-A-holds-two": (
        ((0, 0), (1, 2)),
        ((1, 1),),
        PairEvaluation(Outcome.DIFFERENT, (0, 0), (1, 1)),
    ),
}


@pytest.mark.parametrize(
    ("result_a", "result_b", "expected"),
    OTHER_SIDE_CASES.values(),
    ids=OTHER_SIDE_CASES.keys(),
)
def test_evaluate_pair_gives_outcome_of_other_side(
    result_a: SampleResult, result_b: SampleResult, expected: PairEvaluation
) -> None:
    assert evaluate_pair(result_a, result_b) == expected


# The order of §5's "first such pair": A's genotypes in order and, for each, B's in order.
REPORTED_PAIR_CASES: dict[str, tuple[SampleResult, SampleResult, PairEvaluation]] = {
    # Only A's second genotype has B's ALT set.
    "same-alleles-second-of-a": (
        ((0, 1), (0, 2)),
        ((2, 2),),
        PairEvaluation(Outcome.SAME_ALLELES, (0, 2), (2, 2)),
    ),
    # Two pairs match: A's order picks (0, 1), B's order would pick (1, 2).
    "same-genotype-in-order-of-a": (
        ((0, 1), (1, 2)),
        ((1, 2), (0, 1)),
        PairEvaluation(Outcome.SAME_GENOTYPE, (0, 1), (0, 1)),
    ),
    # A rule that names no pair reports the first genotype of each, here where each
    # sample has two.
    "different-first-of-each": (
        ((0, 1), (0, 2)),
        ((1, 2), (2, 3)),
        PairEvaluation(Outcome.DIFFERENT, (0, 1), (1, 2)),
    ),
    "only-in-b-first-of-each": (
        ((0, 0),),
        ((0, 1), (0, 2)),
        PairEvaluation(Outcome.ONLY_IN_B, (0, 0), (0, 1)),
    ),
    "missing-in-a-first-of-each": (
        ((1, MISSING), (2, MISSING)),
        ((1, 2),),
        PairEvaluation(Outcome.MISSING_IN_A, (1, MISSING), (1, 2)),
    ),
}


@pytest.mark.parametrize(
    ("result_a", "result_b", "expected"),
    REPORTED_PAIR_CASES.values(),
    ids=REPORTED_PAIR_CASES.keys(),
)
def test_evaluate_pair_reports_first_pair_in_order_of_a_then_b(
    result_a: SampleResult, result_b: SampleResult, expected: PairEvaluation
) -> None:
    assert evaluate_pair(result_a, result_b) == expected


def test_make_outcome_rows_fills_every_field_for_two_pairs_sharing_a_sample() -> None:
    # Locus 10-16, ATATAT, one STR segment. File A: S0 is 0/1 at the SNP at 10, S1 is
    # 1/1 at the SNP at 13, whose QUAL is absent, and both are 0/0 at a record with no
    # ALT. File B: T0 is 1/1 at the SNP at 10. Both pairs share T0.
    str_segment = GenomeSegment(CONTIG, 10, 16, RegionType.STR, "AT")
    records_a = (
        make_record(10, ("A", "G"), 30.0, [make_genotype((0, 1), 20), make_genotype((0, 0), 50)]),
        make_record(11, ("T",), 99.0, [make_genotype((0, 0), 90), make_genotype((0, 0), 90)]),
        make_record(13, ("T", "C"), None, [make_genotype((0, 0), 40), make_genotype((1, 1), 15)]),
    )
    records_b = (make_record(10, ("A", "G"), 60.0, [make_genotype((1, 1), 99)]),)
    locus = Locus(CONTIG, 10, 16, records_a, records_b, (str_segment,))
    sample_pairs = (
        SamplePair(name_a="S0", index_a=0, name_b="T0", index_b=0),
        SamplePair(name_a="S1", index_a=1, name_b="T0", index_b=0),
    )
    reconciled = reconcile(locus, "ATATAT", sample_pairs)

    rows = list(make_outcome_rows(reconciled, sample_pairs, kind=Kind.STR, stratum="str"))

    alleles = ("ATATAT", "GTATAT", "ATACAT")
    assert rows == [
        OutcomeRow(
            contig=CONTIG,
            start=10,
            end=16,
            kind=Kind.STR,
            stratum="str",
            str_regions=(str_segment,),
            sample_a="S0",
            sample_b="T0",
            outcome=Outcome.SAME_ALLELES,
            alleles=alleles,
            genotype_a=(0, 1),
            genotype_b=(1, 1),
            possible_a=1,
            possible_b=1,
            quality_a=30.0,
            quality_b=60.0,
            gq_a=20,
            gq_b=99,
            records_a=2,
            records_b=1,
            invalid_reason_a=None,
            invalid_reason_b=None,
        ),
        OutcomeRow(
            contig=CONTIG,
            start=10,
            end=16,
            kind=Kind.STR,
            stratum="str",
            str_regions=(str_segment,),
            sample_a="S1",
            sample_b="T0",
            outcome=Outcome.DIFFERENT,
            alleles=alleles,
            genotype_a=(2, 2),
            genotype_b=(1, 1),
            possible_a=1,
            possible_b=1,
            quality_a=None,
            quality_b=60.0,
            gq_a=15,
            gq_b=99,
            records_a=2,
            records_b=1,
            invalid_reason_a=None,
            invalid_reason_b=None,
        ),
    ]


def test_make_outcome_rows_shows_partner_genotype_in_row_invalid_in_a() -> None:
    # Locus 10-12, AC. A's GT is malformed. B is 0/1 at two unphased SNPs, so it has two
    # possible genotypes, AC/GT and AT/GC, and the row shows the first.
    records_a = (make_record(10, ("A", "G"), 30.0, [MALFORMED_GT]),)
    records_b = (
        make_record(10, ("A", "G"), 50.0, [make_genotype((0, 1), 40)]),
        make_record(11, ("C", "T"), 45.0, [make_genotype((0, 1), 41)]),
    )
    locus = Locus(CONTIG, 10, 12, records_a, records_b, ())
    sample_pairs = (SamplePair(name_a="S0", index_a=0, name_b="T0", index_b=0),)
    reconciled = reconcile(locus, "AC", sample_pairs)

    rows = list(make_outcome_rows(reconciled, sample_pairs, kind=Kind.SNP, stratum=None))

    assert rows == [
        OutcomeRow(
            contig=CONTIG,
            start=10,
            end=12,
            kind=Kind.SNP,
            stratum=None,
            str_regions=(),
            sample_a="S0",
            sample_b="T0",
            outcome=Outcome.INVALID_IN_A,
            alleles=("AC", "GT", "AT", "GC"),
            genotype_a=None,
            genotype_b=(0, 1),
            possible_a=0,
            possible_b=2,
            quality_a=None,
            quality_b=45.0,
            gq_a=None,
            gq_b=40,
            records_a=1,
            records_b=2,
            invalid_reason_a=InvalidReason.MALFORMED_GT,
            invalid_reason_b=None,
        )
    ]


def test_make_outcome_rows_gives_reason_and_partner_genotype_in_rows_too_complex() -> None:
    # Locus 10-12, AC, with a phasing cap of 1. A's S0 is 0/1 at two unphased SNPs, two
    # assignments, so too complex. In B, T0's GT is malformed and T1 is 0/1 at the SNP
    # at 10. Both pairs share S0.
    records_a = (
        make_record(10, ("A", "G"), 30.0, [make_genotype((0, 1), 20)]),
        make_record(11, ("C", "T"), 35.0, [make_genotype((0, 1), 25)]),
    )
    records_b = (make_record(10, ("A", "G"), 50.0, [MALFORMED_GT, make_genotype((0, 1), 60)]),)
    locus = Locus(CONTIG, 10, 12, records_a, records_b, ())
    sample_pairs = (
        SamplePair(name_a="S0", index_a=0, name_b="T0", index_b=0),
        SamplePair(name_a="S0", index_a=0, name_b="T1", index_b=1),
    )
    reconciled = reconcile(locus, "AC", sample_pairs, phasing_cap=1)

    rows = list(make_outcome_rows(reconciled, sample_pairs, kind=Kind.SNP, stratum="."))

    common = {
        "contig": CONTIG,
        "start": 10,
        "end": 12,
        "kind": Kind.SNP,
        "stratum": ".",
        "str_regions": (),
        "sample_a": "S0",
        "outcome": Outcome.TOO_COMPLEX,
        "alleles": ("AC", "GC"),
        "genotype_a": None,
        "possible_a": None,
        "quality_a": 30.0,
        "gq_a": 20,
        "records_a": 2,
        "records_b": 1,
        "invalid_reason_a": None,
    }
    assert rows == [
        OutcomeRow(
            **common,
            sample_b="T0",
            genotype_b=None,
            possible_b=0,
            quality_b=None,
            gq_b=None,
            invalid_reason_b=InvalidReason.MALFORMED_GT,
        ),
        OutcomeRow(
            **common,
            sample_b="T1",
            genotype_b=(0, 1),
            possible_b=1,
            quality_b=50.0,
            gq_b=60,
            invalid_reason_b=None,
        ),
    ]


def test_make_outcome_rows_shows_reported_pair_not_first_genotypes() -> None:
    # E6's sample results, given by hand: the reported pair is the second genotype of each.
    locus = Locus(CONTIG, 10, 11, (), (), ())
    reconciled = ReconciledLocus(
        locus=locus,
        alleles=("A", "G", "T"),
        results_a={0: ((0, 1), (1, 2))},
        results_b={0: ((1, 1), (1, 2))},
    )
    sample_pairs = (SamplePair(name_a="S0", index_a=0, name_b="T0", index_b=0),)

    (row,) = make_outcome_rows(reconciled, sample_pairs, kind=Kind.SNP, stratum=None)

    assert (row.outcome, row.genotype_a, row.genotype_b, row.possible_a, row.possible_b) == (
        Outcome.SAME_GENOTYPE,
        (1, 2),
        (1, 2),
        2,
        2,
    )


# The records of file A at locus 10-14, ACGT, each with its QUAL and the one sample's VCF
# genotype, and the quality_a and gq_a of the row.
type _QualityCase = tuple[Sequence[VcfRecord], float | None, int | None]

QUALITY_CASES: dict[str, _QualityCase] = {
    # The `*` at 11 is an ALT of its record; the record at 12, where the sample is 0/0,
    # has a lower QUAL and GQ and is passed over.
    "star-counts-as-alt": (
        (
            make_record(10, ("ACG", "A"), 40.0, [make_genotype((0, 1), 30)]),
            make_record(11, ("C", "T", "*"), 25.0, [make_genotype((0, 2), 12)]),
            make_record(12, ("G", "C"), 5.0, [make_genotype((0, 0), 3)]),
        ),
        25.0,
        12,
    ),
    # QUAL absent at 10 and GQ absent at 13: the lowest of the others.
    "absent-qual-and-gq-passed-over": (
        (
            make_record(10, ("A", "G"), None, [make_genotype((1, 1), 10)]),
            make_record(12, ("G", "C"), 50.0, [make_genotype((0, 1), 35)]),
            make_record(13, ("T", "C"), 45.0, [make_genotype((0, 1), None)]),
        ),
        45.0,
        10,
    ),
    "missing-and-ref-passed-over": (
        (
            make_record(10, ("A", "G"), 10.0, [make_genotype((MISSING, MISSING), 5)]),
            make_record(12, ("G", "C"), 20.0, [make_genotype((0, 0), 6)]),
        ),
        None,
        None,
    ),
    "malformed-and-other-ploidy-passed-over": (
        (
            make_record(10, ("A", "G"), 10.0, [MALFORMED_GT]),
            make_record(12, ("G", "C"), 20.0, [OTHER_PLOIDY_GT]),
        ),
        None,
        None,
    ),
}


@pytest.mark.parametrize(
    ("records_a", "quality", "gq"), QUALITY_CASES.values(), ids=QUALITY_CASES.keys()
)
def test_make_outcome_rows_gives_lowest_quality_and_gq_of_records_with_alt(
    records_a: Sequence[VcfRecord], quality: float | None, gq: int | None
) -> None:
    locus = Locus(CONTIG, 10, 14, tuple(records_a), (), ())
    sample_pairs = (SamplePair(name_a="S0", index_a=0, name_b="T0", index_b=0),)
    reconciled = reconcile(locus, "ACGT", sample_pairs)

    (row,) = make_outcome_rows(reconciled, sample_pairs, kind=Kind.COMPLEX, stratum=None)

    assert (row.quality_a, row.gq_a, row.quality_b, row.gq_b) == (quality, gq, None, None)
