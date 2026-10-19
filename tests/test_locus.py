"""The checks of Locus (doc/specs/reconciliation.md §7, doc/spec.md §4)."""

import pytest

from pvc_eval.locus import Locus
from pvc_eval.vcf import VcfGenotype, VcfRecord


def make_record(
    contig: str, start: int, alleles: tuple[str, ...], span: tuple[int, int]
) -> VcfRecord:
    """A VCF record of one sample, 0/0, 0-based and half-open."""
    genotype = VcfGenotype(alleles=(0, 0), phased=False, phase_set=None, quality=None)
    return VcfRecord(
        contig=contig,
        start=start,
        alleles=alleles,
        qual=None,
        genotypes=(genotype,),
        line_number=3,
        span=span,
    )


def test_locus_holds_a_record_with_no_alt_whose_ref_reaches_beyond_it() -> None:
    record = make_record("chr1", 2, ("ACGTA",), (2, 7))

    locus = Locus("chr1", 4, 6, (record,), (), ())

    assert locus.records_a == (record,)


@pytest.mark.parametrize("end", [5, 4], ids=["end-at-start", "end-before-start"])
def test_locus_raises_for_an_end_not_after_its_start(end: int) -> None:
    with pytest.raises(ValueError, match=f"a locus's end must be after its start: chr1, 5-{end}"):
        Locus("chr1", 5, end, (), (), ())


def test_locus_raises_for_a_record_on_another_contig() -> None:
    record = make_record("chr2", 5, ("A",), (5, 6))

    with pytest.raises(ValueError, match="the locus chr1:4-10 holds a VCF record of chr2, line 3"):
        Locus("chr1", 4, 10, (), (record,), ())


@pytest.mark.parametrize("span", [(3, 6), (5, 11)], ids=["span-starts-before", "span-ends-after"])
def test_locus_raises_for_a_record_with_an_alt_whose_span_is_not_within_it(
    span: tuple[int, int],
) -> None:
    record = make_record("chr1", 5, ("A", "C"), span)

    with pytest.raises(
        ValueError,
        match=f"the locus chr1:4-10 holds the VCF record of line 3, whose span,"
        f" {span[0]}-{span[1]}, is not within it",
    ):
        Locus("chr1", 4, 10, (record,), (), ())
