"""The loci of build_loci on the spans parse_vcf computes (doc/specs/loci.md, "How it is
verified"): the traded events of doc/spec.md §5.3, the second and third cases of §5.2, a
deletion that only its span joins to the next record, and the STR regions of segment_contig.

The VCF text and the FASTA with its .fai are written in tmp_path, so that a span parse_vcf
computes differently from the one the spec's case assumes shows here. The spec writes
positions 1-based and inclusive; the loci here are 0-based and half-open, so the spec's
"21-25" is (20, 25).
"""

from collections.abc import Iterator
from pathlib import Path

import pytest

from pvc_eval import segments
from pvc_eval.fasta import IndexedFasta
from pvc_eval.loci import build_loci
from pvc_eval.locus import Locus
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.vcf import VcfGenotype, VcfRecord, parse_vcf

HEADER = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
FIRST_RECORD_LINE = 3  # the line of a VCF record that follows the two lines of HEADER
HOMOZYGOUS_ALT = VcfGenotype(alleles=(1, 1), phased=False, phase_set=None, quality=None)


def write_reference(tmp_path: Path, sequence: str) -> Path:
    """Writes a FASTA of one contig, chr1, on one line, as tmp_path/ref.fa, and beside it
    the .fai `samtools faidx` writes for it; returns the FASTA's path."""
    path = tmp_path / "ref.fa"
    path.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
    fai = f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n"
    Path(f"{path}.fai").write_text(fai, encoding="utf-8")
    return path


def write_vcf(tmp_path: Path, name: str, lines: list[str]) -> Path:
    """A VCF file of one sample, HG002, whose records are lines, given as
    `POS REF ALT` on chr1, each with the GT 1/1."""
    body = ""
    for line in lines:
        pos_1based, ref, alt = line.split()
        body += f"chr1\t{pos_1based}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t1/1\n"
    path = tmp_path / name
    path.write_text(HEADER + body, encoding="utf-8")
    return path


def find_loci(
    tmp_path: Path,
    reference: IndexedFasta,
    lines_a: list[str],
    lines_b: list[str],
    *,
    segments_of_chr1: Iterator[GenomeSegment] | None = None,
    chain_distance: int,
) -> list[Locus]:
    """The loci of two VCF files, read by parse_vcf, with no STR region unless given."""
    records_a = parse_vcf(write_vcf(tmp_path, "a.vcf", lines_a), reference)
    records_b = parse_vcf(write_vcf(tmp_path, "b.vcf", lines_b), reference)
    return list(
        build_loci(
            records_a,
            records_b,
            iter(()) if segments_of_chr1 is None else segments_of_chr1,
            ("chr1",),
            chain_distance=chain_distance,
        )
    )


def make_record(
    start: int, alleles: tuple[str, ...], span: tuple[int, int], index_in_file: int
) -> VcfRecord:
    """The record parse_vcf gives for a line of write_vcf, 0-based and half-open;
    index_in_file counts the records of the file from 0."""
    return VcfRecord(
        contig="chr1",
        start=start,
        alleles=alleles,
        qual=None,
        genotypes=(HOMOZYGOUS_ALT,),
        line_number=FIRST_RECORD_LINE + index_in_file,
        span=span,
    )


# The traded events of doc/spec.md §5.3, reference 5-10 AGAACA, with four bases on each side
# that slide nothing. The spans are 5-6 and 9 in file A, and 6 and 8-9 in file B, with one
# position, 7, between them.
TRADED_REFERENCE = "TTTT" + "AGAACA" + "TTTT"
TRADED_LINES_A = ["5 AG A", "9 C A"]
TRADED_LINES_B = ["6 G A", "8 AC A"]
TRADED_A_DELETION = make_record(4, ("AG", "A"), (4, 6), 0)
TRADED_A_SNP = make_record(8, ("C", "A"), (8, 9), 1)
TRADED_B_SNP = make_record(5, ("G", "A"), (5, 6), 0)
TRADED_B_DELETION = make_record(7, ("AC", "A"), (7, 9), 1)


@pytest.mark.parametrize(
    ("chain_distance", "expected"),
    [
        (
            1,
            [
                Locus(
                    "chr1",
                    4,
                    9,
                    (TRADED_A_DELETION, TRADED_A_SNP),
                    (TRADED_B_SNP, TRADED_B_DELETION),
                    (),
                )
            ],
        ),
        (
            0,
            [
                Locus("chr1", 4, 6, (TRADED_A_DELETION,), (TRADED_B_SNP,), ()),
                Locus("chr1", 7, 9, (TRADED_A_SNP,), (TRADED_B_DELETION,), ()),
            ],
        ),
    ],
    ids=["one_locus_at_1", "two_loci_at_0"],
)
def test_traded_events_read_from_vcf_share_a_locus_from_a_chain_distance_of_one(
    tmp_path: Path, chain_distance: int, expected: list[Locus]
) -> None:
    with IndexedFasta(write_reference(tmp_path, TRADED_REFERENCE)) as reference:
        loci = find_loci(
            tmp_path, reference, TRADED_LINES_A, TRADED_LINES_B, chain_distance=chain_distance
        )

    assert loci == expected


def test_deletion_that_can_slide_right_shares_a_locus_with_the_record_it_slides_onto(
    tmp_path: Path,
) -> None:
    # The second case of doc/spec.md §5.2, reference 20-27 GCAGAGTC: A's `21 CAG C` spans
    # 21-25, since the AG at 22-23 repeats at 24-25, and so reaches B's `24 AGT C` at 24-26.
    # Their REFs, 21-23 and 24-26, already touch, so they would share a locus without the
    # span too: the test checks the spans parse_vcf gives, and the next one that only the
    # span joins them.
    fasta = write_reference(tmp_path, "T" * 19 + "GCAGAGTC" + "T" * 4)

    with IndexedFasta(fasta) as reference:
        loci = find_loci(
            tmp_path, reference, ["21 CAG C", "26 T C"], ["24 AGT C"], chain_distance=0
        )

    deletion_a = make_record(20, ("CAG", "C"), (20, 25), 0)
    snp_a = make_record(25, ("T", "C"), (25, 26), 1)
    complex_b = make_record(23, ("AGT", "C"), (23, 26), 0)
    assert loci == [Locus("chr1", 20, 26, (deletion_a, snp_a), (complex_b,), ())]


def test_deletion_joins_a_record_its_ref_leaves_a_gap_to_only_by_its_span(
    tmp_path: Path,
) -> None:
    # doc/specs/loci.md, "How it is verified", 1-based: reference 20-29 GCAGAGAGTC. A's
    # `21 CAG C` deletes the AG at 22-23, which slides right to 26-27, so its span is 21-27
    # and touches B's `28 T C`. By their REFs, 21-23 and 28, a gap of 4 would make two loci
    # at a chain distance of 0.
    fasta = write_reference(tmp_path, "T" * 19 + "GCAGAGAGTC" + "T" * 4)

    with IndexedFasta(fasta) as reference:
        loci = find_loci(tmp_path, reference, ["21 CAG C"], ["28 T C"], chain_distance=0)

    deletion_a = make_record(20, ("CAG", "C"), (20, 27), 0)
    snp_b = make_record(27, ("T", "C"), (27, 28), 0)
    assert loci == [Locus("chr1", 20, 28, (deletion_a,), (snp_b,), ())]


def test_insertion_shares_a_locus_with_a_record_on_the_base_after_its_anchor(
    tmp_path: Path,
) -> None:
    # The third case of doc/spec.md §5.2, reference 10-11 CG: A's `10 C CT` cannot slide and
    # still spans 10-11, the base after its anchor, which B's `11 G TA` at 11 is on. Their
    # REFs, 10 and 11, already touch, so they would share a locus without the span too: the
    # test checks the spans parse_vcf gives.
    fasta = write_reference(tmp_path, "A" * 9 + "CG" + "A" * 4)

    with IndexedFasta(fasta) as reference:
        loci = find_loci(tmp_path, reference, ["10 C CT", "11 G A"], ["11 G TA"], chain_distance=0)

    insertion_a = make_record(9, ("C", "CT"), (9, 11), 0)
    snp_a = make_record(10, ("G", "A"), (10, 11), 1)
    complex_b = make_record(10, ("G", "TA"), (10, 11), 0)
    assert loci == [Locus("chr1", 9, 11, (insertion_a, snp_a), (complex_b,), ())]


# An AT repeat of 10 copies at 26-45, 1-based, with 25 bases on each side that the STR
# finder finds no repeat in. Its segments, found with segment_contig at 6d179b0 on the
# branch str-finder, are, 0-based and half-open, a generic region at (0, 25), the STR at
# (25, 45) and a generic region at (45, 70).
AT_REPEAT_REFERENCE = "GCGTCAGCTGCAGTCCGTGAGCATG" + "AT" * 10 + "CGGTACCTGAGGCTCAGCGTGCAAG"


def test_insertion_before_an_at_repeat_found_by_segment_contig_holds_the_repeat(
    tmp_path: Path,
) -> None:
    # Read with getattr while segment_contig is on the branch str-finder only, so that this
    # module imports without it.
    segment_contig = getattr(segments, "segment_contig", None)
    if segment_contig is None:
        pytest.skip("segment_contig is not in pvc_eval.segments yet")

    # 1-based, `25 G GAT` inserts AT after the G at 25, before the repeat, and slides right
    # through it: its span runs from its anchor to the C after the repeat, 25-46, which is
    # (24, 46) 0-based.
    with IndexedFasta(write_reference(tmp_path, AT_REPEAT_REFERENCE)) as reference:
        loci = find_loci(
            tmp_path,
            reference,
            ["25 G GAT"],
            [],
            segments_of_chr1=segment_contig(reference, "chr1"),
            chain_distance=0,
        )

    insertion_a = make_record(24, ("G", "GAT"), (24, 46), 0)
    at_repeat = GenomeSegment("chr1", 25, 45, RegionType.STR, "AT")
    assert loci == [Locus("chr1", 24, 46, (insertion_a,), (), (at_repeat,))]
