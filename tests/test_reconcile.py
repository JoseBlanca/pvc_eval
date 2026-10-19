"""Reconciliation, doc/specs/reconciliation.md §4: the cases of §4.4 and the other cases of
the `reconcile` item of its §9.

The cases give positions 1-based and inclusive, as the spec does; a Locus is 0-based and
half-open, so a locus 9-12 of the spec is Locus(start=8, end=12).
"""

from collections.abc import Iterator, Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

import pytest

from pvc_eval.fasta import IndexedFasta
from pvc_eval.locus import Locus
from pvc_eval.outcome import InvalidReason
from pvc_eval.reconcile import (
    PHASING_CAP,
    TOO_COMPLEX,
    ReconciledLocus,
    SampleResult,
    reconcile,
    reconcile_loci,
)
from pvc_eval.samples import SamplePair
from pvc_eval.vcf import DEFAULT_SPAN_LIMIT, MISSING, VcfGenotype, VcfRecord, parse_vcf

CONTIG = "chr1"
# The contig of every test: N, which equals no base of a VCF record, so no record slides
# into it, around the bases a case gives.
CONTIG_LENGTH = 80


class Reader:
    """Writes a case's reference and VCF files to a directory and reads them back."""

    def __init__(self, directory: Path, stack: ExitStack) -> None:
        self._directory = directory
        self._stack = stack
        self._num_files = 0

    def make_reference(self, ref_start_1based: int, bases: str) -> IndexedFasta:
        """A FASTA of CONTIG, N but for bases from ref_start_1based on, with its .fai."""
        contig = "N" * (ref_start_1based - 1) + bases
        contig += "N" * (CONTIG_LENGTH - len(contig))
        path = self._make_path("fa")
        path.write_text(f">{CONTIG}\n{contig}\n", encoding="utf-8")
        fai = (
            f"{CONTIG}\t{CONTIG_LENGTH}\t{len(CONTIG) + 2}\t{CONTIG_LENGTH}\t{CONTIG_LENGTH + 1}\n"
        )
        Path(f"{path}.fai").write_text(fai, encoding="utf-8")
        return self._stack.enter_context(IndexedFasta(path))

    def read_records(
        self,
        reference: IndexedFasta,
        lines: Sequence[str],
        *,
        num_samples: int = 1,
        ploidy: int = 2,
        span_limit: int = DEFAULT_SPAN_LIMIT,
    ) -> tuple[VcfRecord, ...]:
        """The records parse_vcf reads from lines `POS REF ALT GT...`, one GT per sample,
        each written `GT` or `GT:PS`."""
        samples = "\t".join(f"S{index}" for index in range(num_samples))
        text = "##fileformat=VCFv4.3\n"
        text += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
        for line in lines:
            pos, ref, alt, *gts = line.split()
            sample_fields = [gt if ":" in gt else f"{gt}:." for gt in gts]
            text += (
                f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:PS\t{'\t'.join(sample_fields)}\n"
            )
        path = self._make_path("vcf")
        path.write_text(text, encoding="utf-8")
        return tuple(parse_vcf(path, reference, span_limit=span_limit, ploidy=ploidy))

    def _make_path(self, suffix: str) -> Path:
        self._num_files += 1
        return self._directory / f"file{self._num_files}.{suffix}"


@pytest.fixture
def reader(tmp_path: Path) -> Iterator[Reader]:
    with ExitStack() as stack:
        yield Reader(tmp_path, stack)


def make_locus(
    locus_1based: tuple[int, int],
    records_a: tuple[VcfRecord, ...],
    records_b: tuple[VcfRecord, ...] = (),
) -> Locus:
    """A Locus over the 1-based, inclusive positions of a case."""
    first, last = locus_1based
    return Locus(CONTIG, first - 1, last, records_a, records_b, ())


def get_locus_bases(reference: IndexedFasta, locus: Locus) -> str:
    return reference.get_bases(locus.contig, locus.start, locus.end)


def pair_samples(num_samples: int) -> tuple[SamplePair, ...]:
    return tuple(
        SamplePair(name_a=f"S{index}", index_a=index, name_b=f"S{index}", index_b=index)
        for index in range(num_samples)
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class Case:
    """A case of §4.4: the reference from a 1-based position, the locus, 1-based and
    inclusive, file A's and file B's VCF lines for one sample, and what reconcile gives."""

    ref_start: int
    bases: str
    locus: tuple[int, int]
    lines_a: Sequence[str]
    lines_b: Sequence[str] = ()
    alleles: tuple[str, ...]
    result_a: SampleResult
    result_b: SampleResult = ((0, 0),)
    span_limit: int = DEFAULT_SPAN_LIMIT


CONTRADICTING = InvalidReason.CONTRADICTING_RECORDS

CASES = {
    # R1, reference 9-13 CAGTC, locus 9-12.
    "R1-deletion-1/1-snp-0/1-is-contradicting": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["9 CAGT C 1/1", "11 G A 0/1"],
        alleles=("CAGT",),
        result_a=CONTRADICTING,
    ),
    "R1-deletion-0/1-snp-0/1": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["9 CAGT C 0/1", "11 G A 0/1"],
        alleles=("CAGT", "CAAT", "C"),
        result_a=((1, 2),),
    ),
    "R1-gatk-star-covered-by-deletion": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["9 CAGT C 0/1", "11 G A,* 1/2"],
        alleles=("CAGT", "CAAT", "C"),
        result_a=((1, 2),),
    ),
    "R1-star-on-anchor-base-not-covered": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["9 CAGT C 0/1", "9 C G,* 1/2"],
        alleles=("CAGT",),
        result_a=CONTRADICTING,
    ),
    "R1-star-with-no-edit-not-covered": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["11 G A,* 1/2"],
        alleles=("CAGT",),
        result_a=CONTRADICTING,
    ),
    "R1-star-covered-by-mnp": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(10, 11),
        lines_a=["10 AG TC 0/1", "11 G A,* 1/2"],
        alleles=("AG", "AA", "TC"),
        result_a=((1, 2),),
    ),
    "R1-star-covered-by-missing-deletion": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["9 CAGT C ./.", "11 G A,* 1/2"],
        alleles=("CAGT",),
        result_a=((MISSING, MISSING),),
    ),
    # Beside R1: an insertion inside the `*` record's REF, 10-11, replaces no position of
    # it, so it does not cover the `*`, carried or where the haplotype is missing.
    "R1-beside-insertion-inside-star-ref-does-not-cover": Case(
        ref_start=10,
        bases="CGT",
        locus=(10, 11),
        lines_a=["10 C CA 1|0", "10 CG C,* 2|1"],
        alleles=("CG",),
        result_a=CONTRADICTING,
    ),
    "R1-beside-missing-insertion-inside-star-ref-does-not-cover": Case(
        ref_start=10,
        bases="CGT",
        locus=(10, 11),
        lines_a=["10 C CA ./.", "10 CG C,* 0/2"],
        alleles=("CG",),
        result_a=CONTRADICTING,
    ),
    # R2, reference 19-24 TAAAAG, locus 19-23.
    # Beside R1: a deletion of 10 that ends before the `*` at 11 does not cover it.
    "R1-beside-deletion-ending-before-star-does-not-cover": Case(
        ref_start=9,
        bases="CAGTC",
        locus=(9, 12),
        lines_a=["9 CA C 0/1", "11 G A,* 1/2"],
        alleles=("CAGT",),
        result_a=CONTRADICTING,
    ),
    "R2-gatk-deletion-and-star": Case(
        ref_start=19,
        bases="TAAAAG",
        locus=(19, 23),
        lines_a=["21 AA A 0/1", "22 A C,* 1/2"],
        alleles=("TAAAA", "TAACA", "TAAA"),
        result_a=((1, 2),),
    ),
    "R2-bcftools-norm-star-covered-by-a-placement": Case(
        ref_start=19,
        bases="TAAAAG",
        locus=(19, 23),
        lines_a=["19 TA T 0/1", "22 A C,* 1/2"],
        alleles=("TAAAA", "TAACA", "TAAA"),
        result_a=((1, 2),),
    ),
    # The deletion's span, 19-23, reaches 3 bases beyond its REF, 19-20, so a limit of 2
    # cuts it to the REF, whose only placement deletes 20.
    "R2-cut-span-star-not-covered": Case(
        ref_start=19,
        bases="TAAAAG",
        locus=(19, 23),
        lines_a=["19 TA T 0/1", "22 A C,* 1/2"],
        alleles=("TAAAA",),
        result_a=CONTRADICTING,
        span_limit=2,
    ),
    # Beside R2: a deletion that slides left onto the `*` at 20, and a limit of 2 that cuts
    # its span, 19-23, to its REF, 22-23, whose placements delete 22 or 23.
    "R2-beside-star-covered-by-a-left-slide": Case(
        ref_start=19,
        bases="TAAAAG",
        locus=(19, 23),
        lines_a=["20 A C,* 1/2", "22 AA A 0/1"],
        alleles=("TAAAA", "TCAAA", "TAAA"),
        result_a=((1, 2),),
    ),
    "R2-beside-cut-span-no-left-slide": Case(
        ref_start=19,
        bases="TAAAAG",
        locus=(19, 23),
        lines_a=["20 A C,* 1/2", "22 AA A 0/1"],
        alleles=("TAAAA",),
        result_a=CONTRADICTING,
        span_limit=2,
    ),
    # R3, reference 21-26 CAGAGT, locus 21-25.
    "R3-as-written": Case(
        ref_start=21,
        bases="CAGAGT",
        locus=(21, 25),
        lines_a=["22 A T 1|0", "23 GAG G 1|0"],
        alleles=("CAGAG", "CTG"),
        result_a=((0, 1),),
    ),
    "R3-bcftools-norm-collides": Case(
        ref_start=21,
        bases="CAGAGT",
        locus=(21, 25),
        lines_a=["21 CAG C 1|0", "22 A T 1|0"],
        alleles=("CAGAG",),
        result_a=CONTRADICTING,
    ),
    # R4, both records 1/1.
    "R4-insertion-between-deleted-bases-collides": Case(
        ref_start=10,
        bases="CAGT",
        locus=(10, 12),
        lines_a=["10 CAG C 1/1", "11 A AT 1/1"],
        alleles=("CAG",),
        result_a=CONTRADICTING,
    ),
    "R4-insertion-at-end-of-deleted-stretch": Case(
        ref_start=10,
        bases="CAGTA",
        locus=(10, 14),
        lines_a=["10 CAG C 1/1", "12 G GT 1/1"],
        alleles=("CAGTA", "CTTA"),
        result_a=((1, 1),),
    ),
    "R4-two-insertions-at-one-point-collide": Case(
        ref_start=10,
        bases="CG",
        locus=(10, 11),
        lines_a=["10 C CA 1/1", "10 C CT 1/1"],
        alleles=("CG",),
        result_a=CONTRADICTING,
    ),
    "R4-replacements-side-by-side": Case(
        ref_start=10,
        bases="CG",
        locus=(10, 11),
        lines_a=["10 C A 1/1", "11 G T 1/1"],
        alleles=("CG", "AT"),
        result_a=((1, 1),),
    ),
    # R5, reference 51-54 TAGC, locus 51-53.
    # Beside R4, the insertion first: the T goes in between 11 and 12, which the MNP
    # replaces, and second: the T goes in at the start of the deleted G.
    "R4-beside-insertion-inside-a-later-mnp-collides": Case(
        ref_start=10,
        bases="CAAAG",
        locus=(10, 14),
        lines_a=["11 A AT 1/1", "11 AAA TTT 1/1"],
        alleles=("CAAAG",),
        result_a=CONTRADICTING,
    ),
    "R4-beside-insertion-at-start-of-an-earlier-deletion": Case(
        ref_start=10,
        bases="CGA",
        locus=(10, 12),
        lines_a=["10 CG C 1/1", "10 C CT 1/1"],
        alleles=("CGA", "CTA"),
        result_a=((1, 1),),
    ),
    "R5-insertion-after-anchor": Case(
        ref_start=51,
        bases="TAGC",
        locus=(51, 53),
        lines_a=["52 A C 1|0", "52 AG A,AAG 2|0"],
        alleles=("TAG", "TCAG"),
        result_a=((0, 1),),
    ),
    # R6, reference 10-12 GCT, locus 10-12.
    "R6-one-insertion-two-spellings": Case(
        ref_start=10,
        bases="GCT",
        locus=(10, 12),
        lines_a=["10 G GCA 1/1"],
        lines_b=["11 C CAC 1/1"],
        alleles=("GCT", "GCACT"),
        result_a=((1, 1),),
        result_b=((1, 1),),
    ),
    # R7, reference 30-32 ACG; file B has 30 A T 0|1 and 32 G C 1|0.
    "R7-unphased": Case(
        ref_start=30,
        bases="ACG",
        locus=(30, 32),
        lines_a=["30 A T 0/1", "32 G C 0/1"],
        lines_b=["30 A T 0|1", "32 G C 1|0"],
        alleles=("ACG", "TCC", "ACC", "TCG"),
        result_a=((0, 1), (2, 3)),
        result_b=((2, 3),),
    ),
    "R7-two-phase-sets": Case(
        ref_start=30,
        bases="ACG",
        locus=(30, 32),
        lines_a=["30 A T 0|1:7", "32 G C 1|0:8"],
        lines_b=["30 A T 0|1", "32 G C 1|0"],
        alleles=("ACG", "ACC", "TCG", "TCC"),
        result_a=((0, 3), (1, 2)),
        result_b=((1, 2),),
    ),
    "R7-one-phase-set": Case(
        ref_start=30,
        bases="ACG",
        locus=(30, 32),
        lines_a=["30 A T 0|1", "32 G C 0|1"],
        lines_b=["30 A T 0|1", "32 G C 1|0"],
        alleles=("ACG", "TCC", "ACC", "TCG"),
        result_a=((0, 1),),
        result_b=((2, 3),),
    ),
    "R7-interleaved-phase-sets": Case(
        ref_start=30,
        bases="ACGTA",
        locus=(30, 34),
        lines_a=["30 A T 0|1:7", "32 G C 0|1:8", "34 A G 1|0:7"],
        alleles=("ACGTA", "ACGTG", "TCCTA", "ACCTG", "TCGTA"),
        result_a=((1, 2), (3, 4)),
    ),
    # Beside R7: `./1` makes a block, and MISSING sorts after every index.
    "R7-beside-dot-one-makes-a-block": Case(
        ref_start=30,
        bases="ACG",
        locus=(30, 32),
        lines_a=["30 A T ./1", "32 G C 0/1"],
        alleles=("ACG", "TCC", "TCG"),
        result_a=((1, MISSING), (2, MISSING)),
    ),
    # Beside R7: the PS of a genotype written with `/` is ignored.
    "R7-beside-ps-of-unphased-genotypes-ignored": Case(
        ref_start=30,
        bases="ACG",
        locus=(30, 32),
        lines_a=["30 A T 0/1:7", "32 G C 0/1:7"],
        alleles=("ACG", "TCC", "ACC", "TCG"),
        result_a=((0, 1), (2, 3)),
    ),
    # Beside R7: three blocks, assignment 1 swapping the second and 2 the third (§4.1).
    "R7-beside-three-blocks-in-order": Case(
        ref_start=30,
        bases="ACGTA",
        locus=(30, 34),
        lines_a=["30 A T 0/1", "32 G C 0/1", "34 A G 0/1"],
        alleles=("ACGTA", "TCCTG", "ACCTA", "TCGTG", "ACGTG", "TCCTA", "ACCTG", "TCGTA"),
        result_a=((0, 1), (2, 3), (4, 5), (6, 7)),
    ),
    # R8, reference 40-44 CAAAG, locus 40-44.
    "R8-edits-that-cancel": Case(
        ref_start=40,
        bases="CAAAG",
        locus=(40, 44),
        lines_a=["40 C CA 0/1", "41 AA A 0/1"],
        alleles=("CAAAG", "CAAG", "CAAAAG"),
        result_a=((0, 0), (1, 2)),
    ),
}


@pytest.mark.parametrize("case", CASES.values(), ids=CASES.keys())
def test_reconcile_gives_the_locus_alleles_and_possible_genotypes_of_each_case(
    reader: Reader, case: Case
) -> None:
    reference = reader.make_reference(case.ref_start, case.bases)
    records_a = reader.read_records(reference, case.lines_a, span_limit=case.span_limit)
    records_b = reader.read_records(reference, case.lines_b)
    locus = make_locus(case.locus, records_a, records_b)

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    assert reconciled == ReconciledLocus(
        locus=locus,
        alleles=case.alleles,
        results_a={0: case.result_a},
        results_b={0: case.result_b},
    )


# R9: 14 SNPs of the 15 positions 30-44, each to another base.
R9_BASES = "ACGTACGTACGTACG"
R9_ALTS = "CATGCATGCATGCAT"


def make_snp_lines(num_snps: int, gt: str) -> list[str]:
    return [f"{30 + index} {R9_BASES[index]} {R9_ALTS[index]} {gt}" for index in range(num_snps)]


def test_reconcile_gives_too_complex_for_fourteen_blocks_R9(reader: Reader) -> None:  # noqa: N802 -- the label of the case of §4.4 in the test's id
    reference = reader.make_reference(30, R9_BASES)
    locus = make_locus((30, 44), reader.read_records(reference, make_snp_lines(14, "0/1")))

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    assert reconciled.results_a == {0: TOO_COMPLEX}
    assert reconciled.alleles == (R9_BASES,)


def test_reconcile_reconciles_thirteen_blocks_R9(reader: Reader) -> None:  # noqa: N802 -- the label of the case of §4.4 in the test's id
    reference = reader.make_reference(30, R9_BASES)
    locus = make_locus((30, 44), reader.read_records(reference, make_snp_lines(13, "0/1")))

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    # 4096 assignments, each a different haplotype and its complement, which make 2^13
    # haplotypes in all.
    possible_genotypes = reconciled.results_a[0]
    assert isinstance(possible_genotypes, tuple)
    assert len(possible_genotypes) == 4096
    assert len(reconciled.alleles) == 8192


def test_reconcile_counts_no_block_for_a_homozygous_genotype(reader: Reader) -> None:
    # 13 blocks and a 1/1: 4096 assignments, within the cap; every haplotype carries the A
    # at 43, so the reference is no haplotype.
    reference = reader.make_reference(30, R9_BASES)
    lines = [*make_snp_lines(13, "0/1"), "43 C A 1/1"]
    locus = make_locus((30, 44), reader.read_records(reference, lines))

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    possible_genotypes = reconciled.results_a[0]
    assert isinstance(possible_genotypes, tuple)
    assert len(possible_genotypes) == 4096
    assert len(reconciled.alleles) == 8193


def test_reconcile_gives_one_assignment_to_one_phase_set_R9(reader: Reader) -> None:  # noqa: N802 -- the label of the case of §4.4 in the test's id
    reference = reader.make_reference(30, R9_BASES)
    locus = make_locus((30, 44), reader.read_records(reference, make_snp_lines(14, "0|1:5")))

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    assert reconciled.alleles == (R9_BASES, "CATGCATGCATGCAG")
    assert reconciled.results_a == {0: ((0, 1),)}


@pytest.mark.parametrize(
    ("phasing_cap", "num_snps", "is_too_complex"),
    [(2, 2, False), (2, 3, True), (1, 1, False), (1, 2, True)],
    ids=["cap-2-two-blocks", "cap-2-three-blocks", "cap-1-one-block", "cap-1-two-blocks"],
)
def test_reconcile_gives_too_complex_over_a_phasing_cap_it_is_given(
    reader: Reader, phasing_cap: int, num_snps: int, is_too_complex: bool
) -> None:
    reference = reader.make_reference(30, R9_BASES)
    locus = make_locus((30, 44), reader.read_records(reference, make_snp_lines(num_snps, "0/1")))

    reconciled = reconcile(
        locus, get_locus_bases(reference, locus), pair_samples(1), phasing_cap=phasing_cap
    )

    assert (reconciled.results_a[0] is TOO_COMPLEX) == is_too_complex


@pytest.mark.parametrize(
    ("lines", "reason"),
    [
        # 44 is G in the reference.
        ([*make_snp_lines(13, "0/1"), "43 C A 3/1", "44 A T 0/0"], InvalidReason.REF_MISMATCH),
        ([*make_snp_lines(14, "0/1"), "44 G T 3/1"], InvalidReason.MALFORMED_GT),
        ([*make_snp_lines(14, "0/1"), "44 G T 1"], InvalidReason.PLOIDY_MISMATCH),
        (["30 A C 3/1", "31 C A 1"], InvalidReason.MALFORMED_GT),
    ],
    ids=[
        "R9-malformed-and-ref-mismatch-is-ref-mismatch",
        "R9-malformed-over-cap-is-malformed",
        "R9-other-ploidy-over-cap-is-ploidy-mismatch",
        "malformed-before-ploidy-mismatch",
    ],
)
def test_reconcile_gives_the_first_reason_that_applies(
    reader: Reader, lines: list[str], reason: InvalidReason
) -> None:
    reference = reader.make_reference(30, R9_BASES)
    locus = make_locus((30, 44), reader.read_records(reference, lines))

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    assert reconciled.results_a == {0: reason}
    assert reconciled.alleles == (R9_BASES,)


def test_reconcile_makes_every_sample_invalid_for_a_no_alt_ref_mismatch_R9(reader: Reader) -> None:  # noqa: N802 -- the label of the case of §4.4 in the test's id
    # 28-31 is NNAC in the reference; the record says GGAC from 28, before the locus.
    reference = reader.make_reference(30, R9_BASES)
    records = reader.read_records(reference, ["28 GGAC . 0/0 0/0", "31 C A 0/1 1/1"], num_samples=2)
    locus = make_locus((31, 31), records)

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(2))

    assert reconciled.results_a == {0: InvalidReason.REF_MISMATCH, 1: InvalidReason.REF_MISMATCH}
    assert reconciled.alleles == ("C",)


@pytest.mark.parametrize(
    ("gt", "ploidy", "result"),
    [("./.", 1, ((MISSING,),)), (".", 2, ((MISSING, MISSING),))],
    ids=["R9-./.-ploidy-1-missing", "R9-.-ploidy-2-missing"],
)
def test_reconcile_gives_missing_for_a_gt_of_only_dots(
    reader: Reader, gt: str, ploidy: int, result: SampleResult
) -> None:
    reference = reader.make_reference(30, "ACG")
    records = reader.read_records(reference, [f"31 C T {gt}"], ploidy=ploidy)
    locus = make_locus((30, 32), records)

    reconciled = reconcile(
        locus, get_locus_bases(reference, locus), pair_samples(1), ploidy_a=ploidy
    )

    assert reconciled.results_a == {0: result}


def test_reconcile_checks_impossible_before_missing(reader: Reader) -> None:
    # The haplotype with the SNP collides with the deletion, and the other is missing at 12:
    # one assignment, impossible.
    reference = reader.make_reference(9, "CAGTC")
    records = reader.read_records(reference, ["9 CAGT C 1/1", "11 G A 0/1", "12 T G ./."])
    locus = make_locus((9, 12), records)

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    assert reconciled.results_a == {0: CONTRADICTING}


def test_reconcile_puts_an_insertion_before_a_replacement_at_its_point(reader: Reader) -> None:
    # doc/spec.md §5.2, the third case: both give CTA.
    reference = reader.make_reference(10, "CG")
    records_a = reader.read_records(reference, ["10 C CT 1/1", "11 G A 1/1"])
    records_b = reader.read_records(reference, ["11 G TA 1/1"])
    locus = make_locus((10, 11), records_a, records_b)

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1))

    assert reconciled.alleles == ("CG", "CTA")
    assert reconciled.results_a == {0: ((1, 1),)}
    assert reconciled.results_b == {0: ((1, 1),)}


def test_reconcile_takes_locus_alleles_from_file_a_then_b_in_header_order(reader: Reader) -> None:
    # Reference 30-32 ACG. A's S0 has ACC, A's S1 TCG, B's S0 TCC on one assignment.
    reference = reader.make_reference(30, "ACG")
    records_a = reader.read_records(reference, ["30 A T 0/0 1/1", "32 G C 0/1 0/0"], num_samples=2)
    records_b = reader.read_records(reference, ["30 A T 0/1 0/0", "32 G C 0/1 0/0"], num_samples=2)
    locus = make_locus((30, 32), records_a, records_b)
    sample_pairs = tuple(reversed(pair_samples(2)))

    reconciled = reconcile(locus, get_locus_bases(reference, locus), sample_pairs)

    assert reconciled.alleles == ("ACG", "ACC", "TCG", "TCC")
    assert reconciled.results_a == {0: ((0, 1),), 1: ((2, 2),)}
    assert reconciled.results_b == {0: ((0, 3), (1, 2)), 1: ((0, 0),)}


def test_reconcile_leaves_out_a_sample_in_no_pair(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")
    records_a = reader.read_records(reference, ["30 A T 0/0 1/1"], num_samples=2)
    records_b = reader.read_records(reference, ["30 A T 0/0 1/1"], num_samples=2)
    locus = make_locus((30, 32), records_a, records_b)
    sample_pairs = (SamplePair(name_a="S0", index_a=0, name_b="S0", index_b=0),)

    reconciled = reconcile(locus, get_locus_bases(reference, locus), sample_pairs)

    assert reconciled.alleles == ("ACG",)
    assert reconciled.results_a == {0: ((0, 0),)}
    assert reconciled.results_b == {0: ((0, 0),)}


def test_reconcile_reads_a_record_with_no_alt_reaching_past_the_locus(reader: Reader) -> None:
    # The record's REF, 30-33, reaches before and after the locus, 31-32.
    reference = reader.make_reference(30, "ACGT")
    records = reader.read_records(reference, ["30 ACGT . 0/0 ./."], num_samples=2)
    locus = make_locus((31, 32), records)

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(2))

    assert reconciled.alleles == ("CG",)
    assert reconciled.results_a == {0: ((0, 0),), 1: ((MISSING, MISSING),)}


def test_reconcile_reads_ploidy_1_with_an_alt_and_a_ref(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")
    records_a = reader.read_records(reference, ["30 A T 1 0"], num_samples=2, ploidy=1)
    locus = make_locus((30, 32), records_a)

    reconciled = reconcile(
        locus, get_locus_bases(reference, locus), pair_samples(2), ploidy_a=1, ploidy_b=1
    )

    assert reconciled.alleles == ("ACG", "TCG")
    assert reconciled.results_a == {0: ((1,),), 1: ((0,),)}
    # File B has no record in the locus: all REF, with its own ploidy.
    assert reconciled.results_b == {0: ((0,),), 1: ((0,),)}


def test_reconcile_gives_one_assignment_at_ploidy_1_with_phase_sets(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")
    records = reader.read_records(reference, ["30 A T 1", "31 C G .|.:7", "32 G C .|.:8"], ploidy=1)
    first_genotype = records[1].genotypes[0]
    assert isinstance(first_genotype, VcfGenotype)
    assert first_genotype.phased
    locus = make_locus((30, 32), records)

    reconciled = reconcile(
        locus, get_locus_bases(reference, locus), pair_samples(1), ploidy_a=1, phasing_cap=1
    )

    assert reconciled.results_a == {0: ((MISSING,),)}


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"ploidy_a": 3}, "ploidy_a must be 1 or 2, not 3"),
        ({"ploidy_b": 0}, "ploidy_b must be 1 or 2, not 0"),
        ({"phasing_cap": 0}, "phasing_cap must be a whole number of 1 or more, not 0"),
        (
            {"ploidy_a": 1},
            "file A's VCF record at line 3 has a VCF genotype of 2 alleles, and file A is"
            " reconciled with ploidy 1",
        ),
    ],
    ids=["ploidy-3", "ploidy-0", "cap-0", "genotype-of-other-ploidy"],
)
def test_reconcile_raises_for_an_argument_out_of_range(
    reader: Reader, arguments: dict[str, int], message: str
) -> None:
    reference = reader.make_reference(30, "ACG")
    locus = make_locus((30, 32), reader.read_records(reference, ["30 A T 0/1"]))

    with pytest.raises(ValueError, match=message):
        reconcile(locus, get_locus_bases(reference, locus), pair_samples(1), **arguments)


def test_reconcile_raises_for_bases_of_another_length_than_the_locus(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")
    locus = make_locus((30, 32), reader.read_records(reference, ["30 A T 0/1"]))

    with pytest.raises(ValueError, match="locus_bases has 2 bases, and the locus, 29-32, has 3"):
        reconcile(locus, "AC", pair_samples(1))


@pytest.mark.parametrize(
    ("phasing_cap", "result"),
    [(1, TOO_COMPLEX), (PHASING_CAP, CONTRADICTING)],
    ids=["over-cap-too-complex", "within-cap-contradicting"],
)
def test_reconcile_gives_too_complex_before_contradicting_records(
    reader: Reader, phasing_cap: int, result: SampleResult
) -> None:
    # Two blocks, and two deletions, of 33-34 and of 34-35, that collide on every haplotype.
    reference = reader.make_reference(30, R9_BASES)
    lines = ["30 A C 0/1", "31 C A 0/1", "32 GTA G 1/1", "33 TAC T 1/1"]
    locus = make_locus((30, 44), reader.read_records(reference, lines))

    reconciled = reconcile(
        locus, get_locus_bases(reference, locus), pair_samples(1), phasing_cap=phasing_cap
    )

    assert reconciled.results_a == {0: result}


def test_reconcile_reads_each_file_with_its_own_ploidy(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")
    records_a = reader.read_records(reference, ["30 A T 0/1"])
    records_b = reader.read_records(reference, ["30 A T 1"], ploidy=1)
    locus = make_locus((30, 32), records_a, records_b)
    expected = ReconciledLocus(
        locus=locus, alleles=("ACG", "TCG"), results_a={0: ((0, 1),)}, results_b={0: ((1,),)}
    )

    reconciled = reconcile(locus, get_locus_bases(reference, locus), pair_samples(1), ploidy_b=1)
    reconciled_loci = list(reconcile_loci(iter([locus]), reference, pair_samples(1), ploidy_b=1))

    assert reconciled == expected
    assert reconciled_loci == [expected]


def test_reconcile_loci_passes_the_phasing_cap_on(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")
    locus = make_locus((30, 32), reader.read_records(reference, ["30 A T 0/1", "32 G C 0/1"]))

    reconciled_loci = list(reconcile_loci(iter([locus]), reference, pair_samples(1), phasing_cap=1))

    assert [reconciled.results_a for reconciled in reconciled_loci] == [{0: TOO_COMPLEX}]


def test_reconcile_loci_raises_for_a_ploidy_when_the_stream_starts(reader: Reader) -> None:
    reference = reader.make_reference(30, "ACG")

    reconciled_loci = reconcile_loci(iter(()), reference, pair_samples(1), ploidy_b=3)

    with pytest.raises(ValueError, match="ploidy_b must be 1 or 2, not 3"):
        next(reconciled_loci)


def test_reconcile_loci_gives_what_reconcile_gives_on_each_locus(reader: Reader) -> None:
    # Two loci of one contig, 30-32 and 40-44 (R8), with the bases read from the FASTA.
    reference = reader.make_reference(30, "ACGNNNNNNNCAAAG")
    records_a = reader.read_records(
        reference, ["30 A T 0/1", "32 G C 0/1", "40 C CA 0/1", "41 AA A 0/1"]
    )
    records_b = reader.read_records(reference, ["30 A T 0|1", "32 G C 1|0"])
    loci = [
        make_locus((30, 32), records_a[:2], records_b),
        make_locus((40, 44), records_a[2:]),
    ]

    reconciled_loci = list(reconcile_loci(iter(loci), reference, pair_samples(1)))

    assert reconciled_loci == [
        ReconciledLocus(
            locus=loci[0],
            alleles=("ACG", "TCC", "ACC", "TCG"),
            results_a={0: ((0, 1), (2, 3))},
            results_b={0: ((2, 3),)},
        ),
        ReconciledLocus(
            locus=loci[1],
            alleles=("CAAAG", "CAAG", "CAAAAG"),
            results_a={0: ((0, 0), (1, 2))},
            results_b={0: ((0, 0),)},
        ),
    ]
    assert reconciled_loci == [
        reconcile(locus, get_locus_bases(reference, locus), pair_samples(1)) for locus in loci
    ]
