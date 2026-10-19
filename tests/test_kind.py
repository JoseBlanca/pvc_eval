"""The kind of a locus, locus_kind (doc/specs/loci.md "Kind", doc/spec.md §5.6 and §7.2).

The cases give positions 1-based, as the spec does; the records are built 0-based and
half-open, so `10 A G` starts at 9 and has the span (9, 10).
"""

import pytest

from pvc_eval.kind import Kind, locus_kind
from pvc_eval.locus import Locus
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.vcf import VcfGenotype, VcfRecord

_HETEROZYGOUS = VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=None)


def make_record(start: int, alleles: tuple[str, ...], span: tuple[int, int]) -> VcfRecord:
    """A VCF record on chr1 of one unphased heterozygous sample, 0-based and half-open."""
    return VcfRecord(
        contig="chr1",
        start=start,
        alleles=alleles,
        qual=None,
        genotypes=(_HETEROZYGOUS,),
        line_number=3,
        ref_mismatch=False,
        span=span,
    )


def make_locus(
    start: int,
    end: int,
    records_a: tuple[VcfRecord, ...],
    records_b: tuple[VcfRecord, ...] = (),
    str_regions: tuple[GenomeSegment, ...] = (),
) -> Locus:
    return Locus("chr1", start, end, records_a, records_b, str_regions)


@pytest.mark.parametrize(
    ("locus", "kind"),
    [
        (make_locus(9, 12, (make_record(9, ("CAG", "C", "CAT"), (9, 12)),)), Kind.COMPLEX),
        (make_locus(9, 11, (make_record(9, ("AC", "GT"), (9, 11)),)), Kind.COMPLEX),
        (
            make_locus(
                9, 14, (make_record(9, ("A", "G"), (9, 10)), make_record(13, ("C", "CT"), (13, 14)))
            ),
            Kind.COMPLEX,
        ),
        (make_locus(9, 10, (make_record(9, ("A", "*"), (9, 10)),)), Kind.COMPLEX),
    ],
    ids=["deletion-and-snp-in-one-record", "two-bases-by-two", "records-of-two-kinds", "star-only"],
)
def test_locus_kind_of_the_spec_table(locus: Locus, kind: Kind) -> None:
    assert locus_kind(locus) is kind


def test_alt_equal_to_its_ref_is_left_out() -> None:
    # 10 CA CA,CG: the first ALT makes no edit; the second turns A into G.
    locus = make_locus(9, 11, (make_record(9, ("CA", "CA", "CG"), (9, 11)),))

    assert locus_kind(locus) is Kind.SNP


def test_record_with_no_alt_is_left_out() -> None:
    # 10 ACGT . in file B, whose REF reaches past the SNP at 11 of file A.
    locus = make_locus(
        10, 11, (make_record(10, ("C", "T"), (10, 11)),), (make_record(9, ("ACGT",), (9, 13)),)
    )

    assert locus_kind(locus) is Kind.SNP


def test_star_allele_beside_a_deletion_is_left_out() -> None:
    # 10 AC A deletes 11, and 11 C * in file B says its haplotype carries that deletion.
    locus = make_locus(
        9, 11, (make_record(9, ("AC", "A"), (9, 11)),), (make_record(10, ("C", "*"), (10, 11)),)
    )

    assert locus_kind(locus) is Kind.INDEL


def test_locus_of_snps_of_both_files_is_snp() -> None:
    locus = make_locus(
        9, 16, (make_record(9, ("A", "G"), (9, 10)),), (make_record(14, ("CA", "TA"), (14, 16)),)
    )

    assert locus_kind(locus) is Kind.SNP


def test_locus_of_insertions_and_deletions_is_indel() -> None:
    # 10 A AT,ATT inserts; 14 GCC G in file B deletes two bases.
    locus = make_locus(
        9,
        17,
        (make_record(9, ("A", "AT", "ATT"), (9, 10)),),
        (make_record(13, ("GCC", "G"), (13, 16)),),
    )

    assert locus_kind(locus) is Kind.INDEL


def test_locus_holding_an_str_region_is_str_whatever_its_records() -> None:
    region = GenomeSegment("chr1", 10, 20, RegionType.STR, "AT")
    locus = make_locus(9, 20, (make_record(9, ("A", "G"), (9, 10)),), str_regions=(region,))

    assert locus_kind(locus) is Kind.STR


def test_too_complex_locus_still_has_a_kind() -> None:
    # 14 unphased heterozygous SNPs: 2**13 phasings, over the phasing cap of 4096, so
    # reconciliation gives too complex (doc/spec.md §5.4, §7.2); the kind is still read.
    records = tuple(
        make_record(position, ("A", "C"), (position, position + 1)) for position in range(10, 38, 2)
    )
    locus = make_locus(10, 37, records)

    assert len(records) == 14
    assert locus_kind(locus) is Kind.SNP
