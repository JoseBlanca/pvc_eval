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
        # The A at 5 is deleted, written with the G on its right as the shared base.
        ("TTGAAAGCT", make_record(5, "AG", "G"), (2, ("GA", "G"))),
    ],
    ids=[
        "deletion_in_homopolymer",
        "insertion_in_homopolymer",
        "deletion_in_dinucleotide_repeat",
        "two_alts_moving_alike",
        "spanning_deletion_kept",
        "deletion_anchored_on_the_right",
    ],
)
def test_indel_slides_to_the_left_end_of_its_repeat(
    contig: str, record: VcfRecord, expected: tuple[int, tuple[str, ...]]
) -> None:
    assert placed(normalise(record, contig)) == expected


@pytest.mark.parametrize(
    ("contig", "record", "expected"),
    [
        (DINUCLEOTIDE, make_record(4, "ACA", "A", "ACG"), (4, ("ACA", "A", "ACG"))),
        # An A inserted after position 4 slides 4 bases, a CA inserted there slides 1.
        ("GAAAACACT", make_record(4, "A", "AA", "ACA"), (4, ("A", "AA", "ACA"))),
        (HOMOPOLYMER, make_record(3, "AA", "C"), (3, ("AA", "C"))),
        (HOMOPOLYMER, make_record(0, "GA", "G"), (0, ("GA", "G"))),
        ("ACGTA", make_record(0, "ACGT", "AGGT"), (1, ("C", "G"))),
        (HOMOPOLYMER, make_record(2, "A"), (2, ("A",))),
        # The A at 0 is deleted, and no base on its left can be the anchor.
        ("ACGT", make_record(0, "AC", "C"), (0, ("AC", "C"))),
    ],
    ids=[
        "indel_and_snp",
        "two_indels_moving_apart",
        "complex",
        "already_left_aligned",
        "snp_in_padding_is_trimmed",
        "no_alt",
        "deletion_of_the_first_base",
    ],
)
def test_record_that_cannot_move_stays(
    contig: str, record: VcfRecord, expected: tuple[int, tuple[str, ...]]
) -> None:
    assert placed(normalise(record, contig)) == expected


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


@st.composite
def contig_and_spelt_indel(draw: st.DrawFn) -> tuple[str, VcfRecord, str]:
    """A contig, one indel spelt with its shared base on the left, on the right, or with
    shared bases on both sides, and the haplotype it gives."""
    contig = draw(st.text(alphabet="AC", min_size=4, max_size=30))
    if draw(st.booleans()):
        # Insert before position, which keeps a base on its left and on its right.
        position = draw(st.integers(min_value=1, max_value=len(contig) - 1))
        num_deleted = 0
        inserted = draw(st.text(alphabet="AC", min_size=1, max_size=6))
    else:
        position = draw(st.integers(min_value=1, max_value=len(contig) - 2))
        num_deleted = draw(st.integers(min_value=1, max_value=len(contig) - 1 - position))
        inserted = ""
    after = position + num_deleted
    spelling = draw(st.sampled_from(["left", "right", "both"]))
    if spelling == "left":
        num_left, num_right = 1, 0
    elif spelling == "right":
        num_left, num_right = 0, 1
    else:
        num_left = draw(st.integers(min_value=1, max_value=position))
        num_right = draw(st.integers(min_value=1, max_value=len(contig) - after))
    start = position - num_left
    ref = contig[start : after + num_right]
    alt = contig[start:position] + inserted + contig[after : after + num_right]
    haplotype = contig[:position] + inserted + contig[after:]
    return contig, make_record(start, ref, alt), haplotype


def leftmost_anchored(contig: str, haplotype: str) -> tuple[int, tuple[str, ...]]:
    """By brute force: the first start whose indel, written with one shared base on the
    left, gives the haplotype."""
    size = len(haplotype) - len(contig)
    for start in range(len(contig)):
        if size < 0:
            ref, alt = contig[start : start + 1 - size], contig[start]
        else:
            ref, alt = contig[start], contig[start] + haplotype[start + 1 : start + 1 + size]
        if contig[:start] + alt + contig[start + len(ref) :] == haplotype:
            return start, (ref, alt)
    raise AssertionError("no start gives the haplotype")


@given(contig_and_spelt_indel())
def test_normalised_indel_is_at_its_leftmost_anchored_place(
    case: tuple[str, VcfRecord, str],
) -> None:
    contig, record, haplotype = case

    assert placed(normalise(record, contig)) == leftmost_anchored(contig, haplotype)


@pytest.mark.parametrize(
    ("start", "alleles", "message"),
    [(-1, ("A", "T"), "start must be 0 or more"), (0, ("", "T"), "REF of at least one base")],
    ids=["negative_start", "empty_ref"],
)
def test_record_that_breaks_its_invariants_is_refused(
    start: int, alleles: tuple[str, ...], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        make_record(start, *alleles)
