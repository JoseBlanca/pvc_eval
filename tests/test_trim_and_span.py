import pytest
from hypothesis import given
from hypothesis import strategies as st

from pvc_eval.vcf import VcfRecord, _find_span, _trim

HOMOPOLYMER = "GAAAAT"  # G at 0, A at 1 to 4, T at 5
DINUCLEOTIDE = "GCACACAT"  # G at 0, CA repeated from 1 to 6, T at 7


@pytest.mark.parametrize(
    ("start", "alleles", "expected"),
    [
        (0, ("ACGT", "AGGT"), (1, ("C", "G"))),
        (3, ("AA", "A"), (3, ("AA", "A"))),
        (2, ("AAT", "AT", "*"), (2, ("AA", "A", "*"))),
        (2, ("A", "*"), (2, ("A", "*"))),
        (2, ("AT", "*"), (2, ("AT", "*"))),
        (2, ("A",), (2, ("A",))),
        (4, ("CAGAG", "CAG"), (4, ("CAG", "C"))),
    ],
    ids=[
        "snp_in_padding",
        "anchor_kept",
        "spanning_deletion_ignored",
        "spanning_deletion_only",
        "spanning_deletion_only_not_trimmed",
        "no_alt",
        "shared_bases_at_both_ends",
    ],
)
def test_trim_removes_shared_bases_and_keeps_one(
    start: int, alleles: tuple[str, ...], expected: tuple[int, tuple[str, ...]]
) -> None:
    assert _trim(start, alleles) == expected


# The three cases of doc/spec.md §5.2, with the 1-based positions of the spec made 0-based
# by padding the contig with one N for each position before the case's first.
CASE_1 = "N" * 20 + "CAGAGT"  # 1-based 21 to 26
CASE_2 = "N" * 19 + "GCAGAGTC"  # 1-based 20 to 27
CASE_3 = "N" * 9 + "CG"  # 1-based 10 to 11


@pytest.mark.parametrize(
    ("contig", "start", "alleles", "expected"),
    [
        # §5.2, the first case: 22 A T and 23 GAG G; the deletion spans 21-25.
        (CASE_1, 21, ("A", "T"), (21, 22)),
        (CASE_1, 22, ("GAG", "G"), (20, 25)),
        # The second case: file A's 21 CAG C spans 21-25, 26 T C 26; file B's 24 AGT C 24-26.
        (CASE_2, 20, ("CAG", "C"), (20, 25)),
        (CASE_2, 25, ("T", "C"), (25, 26)),
        (CASE_2, 23, ("AGT", "C"), (23, 26)),
        # The third case: file A's 10 C CT spans 10-11, 11 G A 11; file B's 11 G TA 11.
        (CASE_3, 9, ("C", "CT"), (9, 11)),
        (CASE_3, 10, ("G", "A"), (10, 11)),
        (CASE_3, 10, ("G", "TA"), (10, 11)),
    ],
    ids=[
        "case_1_snp",
        "case_1_deletion",
        "case_2_deletion",
        "case_2_snp",
        "case_2_complex",
        "case_3_insertion_that_cannot_slide",
        "case_3_snp",
        "case_3_complex",
    ],
)
def test_spans_of_the_three_cases_of_the_spec(
    contig: str, start: int, alleles: tuple[str, ...], expected: tuple[int, int]
) -> None:
    assert _find_span(start, alleles, contig) == expected


@pytest.mark.parametrize(
    ("contig", "start", "alleles", "expected"),
    [
        (HOMOPOLYMER, 3, ("AA", "A"), (0, 5)),
        (HOMOPOLYMER, 1, ("GA", "G"), (0, 5)),
        (HOMOPOLYMER, 3, ("A", "AA"), (0, 6)),
        (DINUCLEOTIDE, 4, ("ACA", "A"), (0, 7)),
        (DINUCLEOTIDE, 2, ("A", "ACA"), (0, 8)),
        # The A at 5 is deleted, written with the G on its right as the shared base; the
        # span holds that G, which is in its REF.
        ("TTGAAAGCT", 5, ("AG", "G"), (2, 7)),
        # An A inserted before the T at 5, written with the T as the shared base.
        (HOMOPOLYMER, 5, ("T", "AT"), (0, 6)),
        # An A inserted after position 4 slides from 0 to 5, a CA inserted there from 3 to 8.
        ("GAAAACACT", 4, ("A", "AA", "ACA"), (0, 9)),
        ("GAAAACACT", 4, ("A", "ACA", "AA"), (0, 9)),
        (DINUCLEOTIDE, 4, ("ACA", "A", "ACG"), (0, 7)),
        (HOMOPOLYMER, 3, ("AA", "A", "*"), (0, 5)),
        (HOMOPOLYMER, 2, ("A", "*"), (2, 3)),
        (HOMOPOLYMER, 2, ("A", "T"), (2, 3)),
        (HOMOPOLYMER, 3, ("AA", "C"), (3, 5)),
        (HOMOPOLYMER, 2, ("A",), (2, 3)),
        # The A at 0 is deleted, and no base on its left can be the anchor.
        ("ACGT", 0, ("AC", "C"), (0, 2)),
        # A T inserted after the last base: there is no base after it to take.
        ("GAT", 2, ("T", "TT"), (1, 3)),
    ],
    ids=[
        "deletion_in_homopolymer",
        "deletion_in_homopolymer_left_aligned",
        "insertion_in_homopolymer",
        "deletion_in_dinucleotide_repeat",
        "insertion_in_dinucleotide_repeat",
        "deletion_anchored_on_the_right",
        "insertion_anchored_on_the_right",
        "union_of_two_insertions",
        "union_of_two_insertions_the_first_reaching_right",
        "union_of_deletion_and_snp",
        "spanning_deletion_beside_a_deletion",
        "spanning_deletion_only",
        "snp",
        "complex",
        "no_alt",
        "deletion_of_the_first_base",
        "insertion_after_the_last_base",
    ],
)
def test_span_covers_every_place_the_record_could_be_written(
    contig: str, start: int, alleles: tuple[str, ...], expected: tuple[int, int]
) -> None:
    assert _find_span(start, alleles, contig) == expected


@pytest.mark.parametrize(
    ("contig", "start", "alleles", "expected"),
    [
        # The deletion of the last A of 10002 slides to the first: its anchor would be
        # 10001 bases left of its REF.
        ("G" + "A" * 10002 + "T", 10001, ("AA", "A"), None),
        # Of 10001 As, 10000 bases: the most a span may reach.
        ("G" + "A" * 10001 + "T", 10000, ("AA", "A"), (0, 10002)),
        # An A inserted after the G slides right past 10001 As, to the T.
        ("G" + "A" * 10001 + "T", 0, ("G", "GA"), None),
        ("G" + "A" * 9999 + "T", 0, ("G", "GA"), (0, 10001)),
        # The A at 10003, deleted after two shared As, slides 10002 bases left to 1, and
        # its anchor, the G at 0, is 10001 bases left of the REF: a slide counted from
        # the REF's start, without the shared bases, would stop one base short of it.
        ("G" + "A" * 10003 + "CT", 10001, ("AAAC", "AAC", "TAAG"), None),
    ],
    ids=[
        "deletion_beyond_10_kb",
        "deletion_at_10_kb",
        "insertion_beyond_10_kb",
        "insertion_at_10_kb",
        "deletion_beyond_10_kb_of_a_long_ref",
    ],
)
def test_span_that_reaches_more_than_10_kb_beyond_its_ref_is_none(
    contig: str, start: int, alleles: tuple[str, ...], expected: tuple[int, int] | None
) -> None:
    assert _find_span(start, alleles, contig) == expected


def make_haplotype(contig: str, start: int, ref: str, alt: str) -> str:
    return contig[:start] + alt + contig[start + len(ref) :]


@st.composite
def make_contig_and_record(draw: st.DrawFn) -> tuple[str, int, tuple[str, ...]]:
    """A contig, and a record on it whose REF agrees with it: up to three ALTs, each of
    them any sequence or `*`."""
    # Two letters, so that repeats, where indels slide, are common.
    contig = draw(st.text(alphabet="AC", min_size=1, max_size=30))
    start = draw(st.integers(min_value=0, max_value=len(contig) - 1))
    end = draw(st.integers(min_value=start + 1, max_value=len(contig)))
    alts = draw(
        st.lists(
            st.text(alphabet="AC", min_size=1, max_size=8) | st.just("*"),
            min_size=0,
            max_size=3,
        )
    )
    return contig, start, (contig[start:end], *alts)


@given(make_contig_and_record())
def test_trimming_keeps_the_sequence_each_alt_gives(
    case: tuple[str, int, tuple[str, ...]],
) -> None:
    contig, start, alleles = case

    trimmed_start, trimmed = _trim(start, alleles)

    for alt, trimmed_alt in zip(alleles[1:], trimmed[1:], strict=True):
        if alt != "*":
            assert make_haplotype(contig, trimmed_start, trimmed[0], trimmed_alt) == make_haplotype(
                contig, start, alleles[0], alt
            )


@given(make_contig_and_record())
def test_span_contains_the_trimmed_ref(case: tuple[str, int, tuple[str, ...]]) -> None:
    contig, start, alleles = case
    trimmed_start, trimmed = _trim(start, alleles)

    span = _find_span(trimmed_start, trimmed, contig)

    assert span is not None
    assert span[0] <= trimmed_start
    assert span[1] >= trimmed_start + len(trimmed[0])


@st.composite
def make_contig_and_indel(draw: st.DrawFn) -> tuple[str, str, int]:
    """A contig, the haplotype one indel gives on it, and the indel's size: the number of
    bases inserted, or minus the number deleted."""
    contig = draw(st.text(alphabet="AC", min_size=2, max_size=30))
    if draw(st.booleans()):
        # Insert before a position, after at least one base.
        position = draw(st.integers(min_value=1, max_value=len(contig)))
        inserted = draw(st.text(alphabet="AC", min_size=1, max_size=6))
        return contig, contig[:position] + inserted + contig[position:], len(inserted)
    position = draw(st.integers(min_value=1, max_value=len(contig) - 1))
    num_deleted = draw(st.integers(min_value=1, max_value=len(contig) - position))
    return contig, contig[:position] + contig[position + num_deleted :], -num_deleted


def find_left_anchored_spellings(
    contig: str, indel_haplotype: str, size: int
) -> list[tuple[int, str, str]]:
    """By brute force: the start, 0-based, REF and ALT of every record that writes the
    indel with one shared base on its left and gives the haplotype."""
    spellings = []
    for start in range(len(contig)):
        if size < 0:
            ref, alt = contig[start : start + 1 - size], contig[start]
        else:
            ref = contig[start]
            alt = ref + indel_haplotype[start + 1 : start + 1 + size]
        if (
            len(ref) == 1 - min(size, 0)
            and make_haplotype(contig, start, ref, alt) == indel_haplotype
        ):
            spellings.append((start, ref, alt))
    return spellings


@given(make_contig_and_indel())
def test_every_spelling_of_one_indel_spans_from_the_first_to_the_last(
    case: tuple[str, str, int],
) -> None:
    contig, indel_haplotype, size = case
    spellings = find_left_anchored_spellings(contig, indel_haplotype, size)
    # The span starts at the anchor of the leftmost spelling. It ends at the last base the
    # rightmost spelling deletes, or at the base after the rightmost insertion's anchor,
    # when the contig has one.
    expected_start = spellings[0][0]
    last_start, last_ref, _ = spellings[-1]
    is_deletion = size < 0
    expected_end = last_start + len(last_ref) if is_deletion else min(last_start + 2, len(contig))

    spans = {_find_span(*_trim(start, (ref, alt)), contig) for start, ref, alt in spellings}

    assert spans == {(expected_start, expected_end)}


@pytest.mark.parametrize(
    ("start", "alleles", "span", "message"),
    [
        (-1, ("A", "T"), (-1, 0), "start must be 0 or more"),
        (0, ("", "T"), (0, 0), "REF of at least one base"),
        (2, ("AA", "A"), (3, 4), "span 3-4 must contain the REF, 2-4"),
        (2, ("AA", "A"), (0, 3), "span 0-3 must contain the REF, 2-4"),
    ],
    ids=["negative_start", "empty_ref", "span_starts_after_ref", "span_ends_before_ref"],
)
def test_record_that_breaks_its_invariants_is_refused(
    start: int, alleles: tuple[str, ...], span: tuple[int, int], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        VcfRecord(
            contig="chr1",
            start=start,
            alleles=alleles,
            qual=None,
            genotypes=(),
            line_number=1,
            ref_mismatch=False,
            span=span,
        )
