import pytest
from hypothesis import given
from hypothesis import strategies as st

from pvc_eval.vcf import VcfRecord, normalise, trim

HOMOPOLYMER = "GAAAAT"  # G at 0, A at 1 to 4, T at 5
DINUCLEOTIDE = "GCACACAT"  # G at 0, CA repeated from 1 to 6, T at 7


def make_record(start: int, *alleles: str) -> VcfRecord:
    return VcfRecord(
        contig="chr1",
        start=start,
        alleles=alleles,
        qual=None,
        genotypes=(),
        genotype_qualities=(),
        line_number=1,
    )


def placed(record: VcfRecord) -> tuple[int, tuple[str, ...]]:
    return record.start, record.alleles


@pytest.mark.parametrize(
    ("record", "expected"),
    [
        (make_record(0, "ACGT", "AGGT"), (1, ("C", "G"))),
        (make_record(3, "AA", "A"), (3, ("AA", "A"))),
        (make_record(2, "AAT", "AT", "*"), (2, ("AA", "A", "*"))),
        (make_record(2, "A"), (2, ("A",))),
    ],
    ids=["snp_in_padding", "anchor_kept", "spanning_deletion_ignored", "no_alt"],
)
def test_trim_removes_shared_bases_and_keeps_one(
    record: VcfRecord, expected: tuple[int, tuple[str, ...]]
) -> None:
    assert placed(trim(record)) == expected


@pytest.mark.parametrize(
    ("contig", "record", "expected"),
    [
        (HOMOPOLYMER, make_record(3, "AA", "A"), (0, ("GA", "G"))),
        (HOMOPOLYMER, make_record(3, "A", "AA"), (0, ("G", "GA"))),
        (DINUCLEOTIDE, make_record(4, "ACA", "A"), (0, ("GCA", "G"))),
        (HOMOPOLYMER, make_record(2, "AAA", "AA", "A"), (0, ("GAA", "GA", "G"))),
        (HOMOPOLYMER, make_record(3, "AA", "A", "*"), (0, ("GA", "G", "*"))),
    ],
    ids=[
        "deletion_in_homopolymer",
        "insertion_in_homopolymer",
        "deletion_in_dinucleotide_repeat",
        "two_alts_moving_alike",
        "spanning_deletion_kept",
    ],
)
def test_indel_slides_to_the_left_end_of_its_repeat(
    contig: str, record: VcfRecord, expected: tuple[int, tuple[str, ...]]
) -> None:
    assert placed(normalise(record, contig)) == expected


@pytest.mark.parametrize(
    ("contig", "record"),
    [
        (DINUCLEOTIDE, make_record(4, "ACA", "A", "ACG")),
        # An A inserted after position 4 slides 4 bases, a CA inserted there slides 1.
        ("GAAAACACT", make_record(4, "A", "AA", "ACA")),
        (HOMOPOLYMER, make_record(3, "AA", "C")),
        (HOMOPOLYMER, make_record(0, "GA", "G")),
    ],
    ids=["indel_and_snp", "two_indels_moving_apart", "complex", "already_left_aligned"],
)
def test_record_that_cannot_move_stays(contig: str, record: VcfRecord) -> None:
    assert placed(normalise(record, contig)) == placed(trim(record))


def test_ref_mismatch_is_marked_and_not_moved() -> None:
    record = make_record(3, "CA", "C")

    normalised = normalise(record, HOMOPOLYMER)

    assert normalised.ref_mismatch
    assert placed(normalised) == (3, ("CA", "C"))


def haplotype(contig: str, record: VcfRecord, alt_index: int) -> str:
    return contig[: record.start] + record.alleles[alt_index] + contig[record.end :]


@st.composite
def contig_and_indel(draw: st.DrawFn) -> tuple[str, VcfRecord]:
    # Two letters, so that repeats, where indels move, are common.
    contig = draw(st.text(alphabet="AC", min_size=4, max_size=30))
    anchor = draw(st.integers(min_value=0, max_value=len(contig) - 2))
    if draw(st.booleans()):
        length = draw(st.integers(min_value=1, max_value=len(contig) - anchor - 1))
        ref = contig[anchor : anchor + 1 + length]
        alt = contig[anchor]
    else:
        ref = contig[anchor]
        alt = ref + draw(st.text(alphabet="AC", min_size=1, max_size=6))
    return contig, make_record(anchor, ref, alt)


@given(contig_and_indel())
def test_normalising_keeps_the_sequence_the_haplotype_carries(
    case: tuple[str, VcfRecord],
) -> None:
    contig, record = case

    normalised = normalise(record, contig)

    assert haplotype(contig, normalised, 1) == haplotype(contig, record, 1)


@given(contig_and_indel())
def test_normalising_twice_changes_nothing(case: tuple[str, VcfRecord]) -> None:
    contig, record = case

    once = normalise(record, contig)

    assert placed(normalise(once, contig)) == placed(once)
