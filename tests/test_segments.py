import random
import tracemalloc
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import (
    GenomeSegment,
    RegionType,
    StrFinderParameters,
    _BundleMember,
    _Candidate,
    _drop_near_contig_ends,
    _drop_under_copy_floor,
    _find_candidates,
    _find_cut,
    _find_maximal_scoring_subsequences,
    _find_str_segments,
    _get_copy_floor,
    _get_purity,
    _is_compound,
    _is_in_cluster,
    _make_bundle_member,
    _make_clusters,
    _Rejection,
    _split_bundles,
    _StrSegment,
)


def test_a_segment_holds_its_contig_span_type_and_motif() -> None:
    segment = GenomeSegment("chr1", 58, 80, RegionType.STR, "CA")
    assert (segment.contig, segment.start, segment.end) == ("chr1", 58, 80)
    assert segment.region_type is RegionType.STR
    assert segment.motif == "CA"


def test_region_types_are_the_four_of_the_contract() -> None:
    assert [region_type.value for region_type in RegionType] == [
        "generic",
        "str",
        "str-bundle",
        "satellite",
    ]


@pytest.mark.parametrize(("start", "end"), [(10, 10), (10, 9)], ids=["empty", "reversed"])
def test_a_segment_whose_end_is_not_after_its_start_raises(start: int, end: int) -> None:
    with pytest.raises(ValueError, match="end"):
        GenomeSegment("chr1", start, end, RegionType.GENERIC, None)


# ---- StrFinderParameters --------------------------------------------------------------------

_ZERO_FLOORS = (0, 0, 0, 0, 0, 0)


def test_the_parameters_default_to_the_spec_table() -> None:
    parameters = StrFinderParameters()
    assert (parameters.scan_min_period, parameters.scan_max_period) == (1, 6)
    assert (parameters.min_period, parameters.max_period) == (1, 6)
    assert (parameters.match_reward, parameters.mismatch_penalty) == (2, 7)
    assert parameters.scan_min_copies == 2
    assert parameters.copy_floors == (5, 5, 4, 4, 4, 3)
    assert parameters.min_purity == 0.8
    assert parameters.min_score == 0
    assert parameters.bundle_distance == 15
    assert parameters.contig_end_flank == 15
    assert parameters.max_str_len == 100
    assert parameters.max_bundle_len == 100


# Each field at its limit, as the spec's table gives it.
_AT_THE_LIMIT: list[tuple[str, dict[str, object]]] = [
    ("scan_min_period", {"scan_min_period": 1}),
    ("scan_min_period_up_to_max", {"scan_min_period": 4, "scan_max_period": 4}),
    ("scan_max_period", {"scan_max_period": 6}),
    ("min_period", {"min_period": 1}),
    ("min_period_up_to_max", {"min_period": 4, "max_period": 4}),
    ("max_period", {"max_period": 6}),
    ("match_reward", {"match_reward": 1}),
    ("mismatch_penalty", {"mismatch_penalty": 1}),
    ("scan_min_copies", {"scan_min_copies": 1}),
    ("copy_floors", {"copy_floors": _ZERO_FLOORS}),
    ("min_purity_zero", {"min_purity": 0.0}),
    ("min_purity_one", {"min_purity": 1.0}),
    ("min_score", {"min_score": -5}),
    ("bundle_distance", {"bundle_distance": 1, "max_str_len": 1}),
    ("contig_end_flank", {"contig_end_flank": 0}),
    ("max_str_len", {"bundle_distance": 15, "max_str_len": 15}),
    ("max_bundle_len", {"max_bundle_len": 0}),
]

# Each field one step past its limit, with the field the error must name.
_PAST_THE_LIMIT: list[tuple[str, dict[str, object], str]] = [
    ("scan_min_period", {"scan_min_period": 0}, "scan_min_period"),
    ("scan_min_period_above_max", {"scan_min_period": 5, "scan_max_period": 4}, "scan_min_period"),
    ("scan_max_period", {"scan_max_period": 7}, "scan_max_period"),
    ("min_period", {"min_period": 0}, "min_period"),
    ("min_period_above_max", {"min_period": 5, "max_period": 4}, "min_period"),
    ("max_period", {"max_period": 7}, "max_period"),
    ("match_reward", {"match_reward": 0}, "match_reward"),
    ("mismatch_penalty", {"mismatch_penalty": 0}, "mismatch_penalty"),
    ("scan_min_copies", {"scan_min_copies": 0}, "scan_min_copies"),
    ("copy_floors_negative", {"copy_floors": (5, 5, 4, 4, 4, -1)}, "copy_floors"),
    ("copy_floors_five", {"copy_floors": (5, 5, 4, 4, 4)}, "copy_floors"),
    ("copy_floors_fraction", {"copy_floors": (5, 5, 4.5, 4, 4, 3)}, "copy_floors"),
    ("min_purity_negative", {"min_purity": -0.01}, "min_purity"),
    ("min_purity_above_one", {"min_purity": 1.01}, "min_purity"),
    ("min_purity_nan", {"min_purity": float("nan")}, "min_purity"),
    ("min_purity_huge_int", {"min_purity": 10**400}, "min_purity"),
    ("min_purity_fraction", {"min_purity": Fraction(4, 5)}, "min_purity"),
    ("min_score_fraction", {"min_score": 0.5}, "min_score"),
    ("bundle_distance", {"bundle_distance": 0}, "bundle_distance"),
    ("contig_end_flank", {"contig_end_flank": -1}, "contig_end_flank"),
    ("max_str_len", {"bundle_distance": 15, "max_str_len": 14}, "max_str_len"),
    ("max_bundle_len", {"max_bundle_len": -1}, "max_bundle_len"),
    ("a_bool_is_not_a_whole_number", {"match_reward": True}, "match_reward"),
]


@pytest.mark.parametrize(
    "fields", [fields for _, fields in _AT_THE_LIMIT], ids=[name for name, _ in _AT_THE_LIMIT]
)
def test_a_parameter_at_its_limit_is_accepted(fields: dict[str, object]) -> None:
    parameters = StrFinderParameters(**fields)  # ty: ignore[invalid-argument-type]  # ty cannot match a dict's values to fields
    for field, value in fields.items():
        assert getattr(parameters, field) == value


@pytest.mark.parametrize(
    ("fields", "field"),
    [(fields, field) for _, fields, field in _PAST_THE_LIMIT],
    ids=[name for name, _, _ in _PAST_THE_LIMIT],
)
def test_a_parameter_past_its_limit_raises_naming_the_field(
    fields: dict[str, object], field: str
) -> None:
    with pytest.raises(ValueError, match=rf"\b{field}\b") as raised:
        StrFinderParameters(**fields)  # ty: ignore[invalid-argument-type]  # values of the wrong type, on purpose
    assert repr(fields[field]) in str(raised.value)


def test_min_purity_must_be_a_float_or_an_int() -> None:
    with pytest.raises(
        ValueError, match=r"^min_purity is Fraction\(4, 5\); it must be a float or int from 0 to 1$"
    ):
        StrFinderParameters(min_purity=Fraction(4, 5))  # ty: ignore[invalid-argument-type]  # a Fraction, on purpose


# ---- Step 1, ported from pvc's tandem_repeat.rs at bbcc9f4 ----------------------------------
# pvc's RepeatInterval is 0-based and half-open, as a candidate is, so its numbers carry over.
# pvc's find_tandem_repeats is step 1 with no copy floor; its PeriodRange is the scanned periods.


def _find_like_pvc(bases: str, min_period: int, max_period: int) -> list[_Candidate]:
    parameters = StrFinderParameters(
        scan_min_period=min_period,
        scan_max_period=max_period,
        copy_floors=_ZERO_FLOORS,
    )
    return _find_candidates(bases, parameters)


def _find_rt_subsequences(scores: list[int], block_length: int) -> list[tuple[int, int, int]]:
    """pvc's rt_segments (tandem_repeat.rs 1208), the scores fed in blocks of block_length."""
    blocks = [scores[i : i + block_length] for i in range(0, len(scores), block_length)]
    subsequences = [
        subsequence
        for found, _ in _find_maximal_scoring_subsequences(iter(blocks))
        for subsequence in found
    ]
    return sorted(subsequences)


def _find_brute_subsequences(scores: list[int]) -> list[tuple[int, int, int]]:
    """pvc's brute_segments (tandem_repeat.rs 1220): the maximal scoring subsequences from
    Ruzzo and Tompa's definition, as (start, end, score), `[start, end)` the indices of
    their scores. pvc's are inclusive; the end here is one past pvc's."""
    prefix = [0]
    for score in scores:
        prefix.append(prefix[-1] + score)

    def get_score(i: int, j: int) -> int:
        return prefix[j + 1] - prefix[i]

    optimal = []
    for i in range(len(scores)):
        for j in range(i, len(scores)):
            total = get_score(i, j)
            if total <= 0:
                continue
            if all(
                (i2, j2) == (i, j) or get_score(i2, j2) < total
                for i2 in range(i, j + 1)
                for j2 in range(i2, j + 1)
            ):
                optimal.append((i, j, total))
    maximal = [
        (i, j + 1, total)
        for i, j, total in optimal
        if not any(a <= i and j <= b and (a, b) != (i, j) for a, b, _ in optimal)
    ]
    return sorted(maximal)


# pvc's PeriodRange is both the periods scanned and the periods admitted here.
@pytest.mark.parametrize(
    ("min_period", "max_period"), [(1, 6), (6, 6), (2, 6)], ids=["1-6", "6-6", "2-6"]
)
def test_period_range_accepts_valid_ranges(min_period: int, max_period: int) -> None:
    scanned = StrFinderParameters(scan_min_period=min_period, scan_max_period=max_period)
    assert (scanned.scan_min_period, scanned.scan_max_period) == (min_period, max_period)
    admitted = StrFinderParameters(min_period=min_period, max_period=max_period)
    assert (admitted.min_period, admitted.max_period) == (min_period, max_period)


@pytest.mark.parametrize(("min_period", "max_period"), [(0, 6), (0, 0)], ids=["0-6", "0-0"])
def test_period_range_rejects_zero_min(min_period: int, max_period: int) -> None:
    with pytest.raises(ValueError, match=r"^scan_min_period is 0"):
        StrFinderParameters(scan_min_period=min_period, scan_max_period=max_period)
    with pytest.raises(ValueError, match=r"^min_period is 0"):
        StrFinderParameters(min_period=min_period, max_period=max_period)


def test_period_range_rejects_min_above_max() -> None:
    with pytest.raises(
        ValueError, match=r"^scan_min_period \(3\) must not exceed scan_max_period \(2\)$"
    ):
        StrFinderParameters(scan_min_period=3, scan_max_period=2)
    with pytest.raises(ValueError, match=r"^min_period \(3\) must not exceed max_period \(2\)$"):
        StrFinderParameters(min_period=3, max_period=2)


def test_a_max_period_under_the_default_min_names_both_fields() -> None:
    with pytest.raises(
        ValueError, match=r"^scan_min_period \(1\) must not exceed scan_max_period \(0\)$"
    ):
        StrFinderParameters(scan_max_period=0)


def test_a_limit_set_by_another_field_names_it_and_its_value() -> None:
    with pytest.raises(
        ValueError, match=r"^max_str_len is 14; it must be bundle_distance \(15\) or more$"
    ):
        StrFinderParameters(max_str_len=14)


def test_default_periods_is_one_to_six() -> None:
    parameters = StrFinderParameters()
    assert (parameters.scan_min_period, parameters.scan_max_period) == (1, 6)


def test_scan_params_default_is_two_seven_two() -> None:
    parameters = StrFinderParameters()
    assert parameters.match_reward == 2
    assert parameters.mismatch_penalty == 7
    assert parameters.scan_min_copies == 2


_CRAFTED_SCORES = [
    [],
    [-7, -7, -7],
    [2, 2, 2, 2],
    [2, -7, 2],  # two singletons: the interruption is too costly to bridge
    [2, -3, 2],  # one merged subsequence: the dip is bridged
    [5, -10, 8],  # dip below the left L: two subsequences
    [5, -3, 8],  # dip stays above: one subsequence
    [2, 2, -7, 2, 2],
    [-7, 2, 2, -7, 2, -7, 2, 2, -7],
    [1, 1, 1, -2, 1, 1, 1, -10, 5, 5],
]


@pytest.mark.parametrize("scores", _CRAFTED_SCORES, ids=[str(s) for s in _CRAFTED_SCORES])
def test_ruzzo_tompa_matches_brute_force_on_crafted_cases(scores: list[int]) -> None:
    expected = _find_brute_subsequences(scores)
    assert _find_rt_subsequences(scores, block_length=max(len(scores), 1)) == expected
    assert _find_rt_subsequences(scores, block_length=1) == expected


def test_ruzzo_tompa_matches_brute_force_on_pseudo_random() -> None:
    # pvc's LCG, a u64 that wraps, with the high 32 bits taken as (state >> 33) as u32.
    state = 0x9E3779B97F4A7C15

    def get_next() -> int:
        nonlocal state
        state = (state * 6364136223846793005 + 1442695040888963407) % 2**64
        return (state >> 33) % 2**32

    for _ in range(2000):
        length = get_next() % 160
        scores = [2 if get_next() % 4 == 0 else -7 for _ in range(length)]  # 25% match-like
        expected = _find_brute_subsequences(scores)
        assert _find_rt_subsequences(scores, block_length=max(length, 1)) == expected, scores
        assert _find_rt_subsequences(scores, block_length=7) == expected, scores


def test_finds_a_clean_perfect_tract() -> None:
    seq = "CAGCAGCAGCAGCAGCAGCAGCAG"  # 8 copies, 24 bp
    got = _find_like_pvc(seq, 3, 3)
    assert got == [_Candidate(start=0, end=24, period=3, score=2 * (24 - 3))]


def test_a_substitution_keeps_one_interval() -> None:
    seq = "CAG" * 5 + "CAT" + "CAG" * 5  # the interrupting copy (G to T)
    got = _find_like_pvc(seq, 3, 3)
    period3 = [candidate for candidate in got if candidate.period == 3]
    assert len(period3) == 1, got
    assert (period3[0].start, period3[0].end) == (0, len(seq))


def test_a_single_indel_keeps_one_interval() -> None:
    seq = "CAG" * 6 + "T" + "CAG" * 6  # a 1 bp insertion
    got = _find_like_pvc(seq, 3, 3)
    period3 = [candidate for candidate in got if candidate.period == 3]
    assert len(period3) == 1, got
    assert (period3[0].start, period3[0].end) == (0, len(seq))


def test_an_n_run_emits_nothing() -> None:
    assert _find_like_pvc("NNNNNNNNNNNNNNNNNNNN", 1, 6) == []


def test_a_soft_masked_tract_is_found() -> None:
    seq = "cagcagcagcagcagcagcagcag"  # lowercase (CAG)*8
    got = _find_like_pvc(seq, 3, 3)
    assert len(got) == 1
    assert (got[0].start, got[0].end, got[0].period) == (0, 24, 3)


def test_only_the_primitive_period_is_emitted() -> None:
    seq = "ATATATATATATATATATAT"  # (AT)*10, 20 bp
    got = _find_like_pvc(seq, 2, 6)
    assert {candidate.period for candidate in got} == {2}, got


def test_a_genuine_higher_primitive_period_is_kept() -> None:
    seq = "ACGTACGTACGTACGT"  # (ACGT)*4, primitive period 4
    got = _find_like_pvc(seq, 2, 6)
    assert any(candidate.period == 4 for candidate in got), got
    assert all(candidate.period not in (2, 3) for candidate in got), got


def test_a_homopolymer_is_only_period_one() -> None:
    seq = "AAAAAAAAAAAAAAAAAAAA"  # poly-A, 20 bp
    got = _find_like_pvc(seq, 1, 6)
    assert {candidate.period for candidate in got} == {1}, got
    assert _find_like_pvc(seq, 2, 6) == []


def test_a_period_multiple_re_detection_is_dropped_by_geometry() -> None:
    # A (GAA)*n repeat with stray bases before it: a sloppy period-6 AAGATG that starts
    # earlier is covered by the period-3 one, and goes.
    seq = "ATGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAG"
    got = _find_like_pvc(seq, 2, 6)
    assert any(candidate.period == 3 for candidate in got), got
    assert not any(candidate.period == 6 for candidate in got), got

    # The AAA of every (GAAA) unit covers only a sliver of the period-4 repeat.
    seq = "GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"  # (GAAA)*8
    got = _find_like_pvc(seq, 1, 6)
    assert any(candidate.period == 4 for candidate in got), got


def test_the_min_copies_floor_drops_short_segments() -> None:
    seq = "ACGTACAGTCTGAC"  # one coincidental period-2 match
    got = _find_like_pvc(seq, 2, 2)
    assert all((candidate.end - candidate.start) // 2 >= 2 for candidate in got), got


# ---- Step 1 with copy floors, pvc_eval's own ------------------------------------------------

# Four A's, a C, then 24 A's. Period 6 alone finds [0, 30), motif AAAACA: 4 matches, the
# mismatch of C against A, then 19 matches, 4 * 2 - 7 + 19 * 2 = 39. The period-1 run of
# A's [5, 30) covers 25 of its 30 bases, so with period 1 scanned the period 6 goes.
_PERIOD_SIX_OVER_AN_A_RUN = "AAAACA" + "A" * 24


def test_a_candidate_under_its_copy_floor_still_drops_a_longer_period() -> None:
    alone = _find_like_pvc(_PERIOD_SIX_OVER_AN_A_RUN, 6, 6)
    assert alone == [_Candidate(start=0, end=30, period=6, score=39)]
    floors = (100, 0, 0, 0, 0, 0)
    parameters = StrFinderParameters(copy_floors=floors)
    assert _find_candidates(_PERIOD_SIX_OVER_AN_A_RUN, parameters) == []


# Step 1's drop rule on inputs where a slip in one of its clauses changes the output. The
# expected candidates, (start, end, period, score) with no copy floor, are those of a direct
# transcription of pvc's find_tandem_repeats at bbcc9f4, written by the review of this
# milestone in the session's scratchpad, rev_spec/pvc_ref.py. On 19 October 2026 the command
# was `python3 -c "from pvc_ref import find_tandem_repeats; print(find_tandem_repeats(b'...',
# 1, 6))"` in that directory, with each sequence in place of the dots.
_DROP_RULE_CASES: list[tuple[str, str, list[tuple[int, int, int, int]]]] = [
    (
        "half_covered_is_dropped",  # an overlap of exactly half drops (18, 34) of period 2
        "GGAGGGGAGAAAAAGGAAAGAGAGAGAAAAAAAAGAGAAAGGAAGAGA",
        [
            (0, 2, 1, 2),
            (3, 7, 1, 6),
            (9, 14, 1, 8),
            (14, 16, 1, 2),
            (16, 19, 1, 4),
            (26, 34, 1, 14),
            (37, 40, 1, 4),
            (40, 42, 1, 2),
            (42, 44, 1, 2),
            (6, 10, 2, 4),
            (33, 38, 2, 6),
            (43, 48, 2, 6),
        ],
    ),
    (
        "only_dividing_periods_drop",  # (11, 19) of period 4 stays beside a period 3
        "CCGACCACAAACCACCCACAAGAGGACCGGAAAACCCACGGACGGGCAACCCGCGGGAGGGGCCGCACCGCAACGAAGG",
        [
            (0, 2, 1, 2),
            (4, 6, 1, 2),
            (8, 11, 1, 4),
            (11, 13, 1, 2),
            (14, 17, 1, 4),
            (19, 21, 1, 2),
            (23, 25, 1, 2),
            (26, 28, 1, 2),
            (28, 30, 1, 2),
            (30, 34, 1, 6),
            (34, 37, 1, 4),
            (39, 41, 1, 2),
            (43, 46, 1, 4),
            (47, 49, 1, 2),
            (49, 52, 1, 4),
            (54, 57, 1, 4),
            (58, 62, 1, 6),
            (62, 64, 1, 2),
            (67, 69, 1, 2),
            (71, 73, 1, 2),
            (75, 77, 1, 2),
            (77, 79, 1, 2),
            (5, 9, 2, 4),
            (16, 20, 2, 4),
            (20, 24, 2, 4),
            (51, 55, 2, 4),
            (10, 16, 3, 6),
            (11, 19, 4, 8),
            (37, 45, 4, 8),
            (62, 72, 5, 10),
        ],
    ),
    (
        "a_dropped_candidate_drops_nothing",  # (6, 19) of period 4 stays
        "GAAGGGAAGAAAGAGAGAGGGGGGGAAAGAAAAGAGAGGGAGGAAGGGAAGGGA",
        [
            (1, 3, 1, 2),
            (3, 6, 1, 4),
            (6, 8, 1, 2),
            (9, 12, 1, 4),
            (18, 25, 1, 12),
            (25, 28, 1, 4),
            (29, 33, 1, 6),
            (37, 40, 1, 4),
            (41, 43, 1, 2),
            (43, 45, 1, 2),
            (45, 48, 1, 4),
            (48, 50, 1, 2),
            (50, 53, 1, 4),
            (32, 38, 2, 8),
            (5, 11, 3, 6),
            (38, 44, 3, 6),
            (6, 19, 4, 9),
            (24, 32, 4, 8),
            (35, 43, 4, 8),
            (25, 35, 5, 10),
            (41, 54, 5, 16),
        ],
    ),
    (
        "a_long_shorter_candidate_drops_by_its_end",  # (34, 54) of period 1 drops (38, 54)
        "CAAAACCCACCCAAACAAAACACAACAACACAACAAAAAAAAACAAAAAAAAAACACCACAAACCACAACCCCAAACCC",
        [
            (1, 5, 1, 6),
            (5, 8, 1, 4),
            (9, 12, 1, 4),
            (12, 15, 1, 4),
            (16, 20, 1, 6),
            (23, 25, 1, 2),
            (26, 28, 1, 2),
            (31, 33, 1, 2),
            (34, 54, 1, 20),
            (56, 58, 1, 2),
            (60, 63, 1, 4),
            (63, 65, 1, 2),
            (67, 69, 1, 2),
            (69, 73, 1, 6),
            (73, 76, 1, 4),
            (76, 79, 1, 4),
            (19, 24, 2, 6),
            (27, 32, 2, 6),
            (53, 57, 2, 4),
            (57, 61, 2, 4),
            (64, 68, 2, 4),
            (21, 30, 3, 12),
            (54, 60, 3, 6),
            (4, 13, 4, 10),
            (11, 19, 4, 8),
            (12, 27, 5, 11),
        ],
    ),
    (
        "a_period_waits_for_the_first_untested_start",  # (11, 23) of period 6 goes
        "ACACCCCAAAAACCCCCACCCCCCAAACACAACCCCACACACCACACCACACCCAAAACAC",
        [
            (3, 7, 1, 6),
            (7, 17, 1, 9),
            (18, 24, 1, 10),
            (24, 27, 1, 4),
            (30, 32, 1, 2),
            (32, 36, 1, 6),
            (41, 43, 1, 2),
            (46, 48, 1, 2),
            (51, 54, 1, 4),
            (54, 58, 1, 6),
            (0, 4, 2, 4),
            (26, 31, 2, 6),
            (35, 42, 2, 10),
            (42, 47, 2, 6),
            (47, 52, 2, 6),
            (57, 61, 2, 4),
            (39, 45, 3, 6),
            (44, 50, 3, 6),
            (37, 53, 5, 22),
        ],
    ),
    (
        # After the tail, the periods 5 and 6 over [0, 30) are final some bases before the
        # run of A's that drops them is, and a block of 1 base ends between the two. The
        # tail was found by searching random ones for a case where testing a candidate
        # before its divisors are known lets both through.
        "a_longer_period_waits_for_a_shorter_one_still_open",
        _PERIOD_SIX_OVER_AN_A_RUN + "CAGATTTTCATATTATGCAG",
        [
            (0, 4, 1, 6),
            (5, 30, 1, 48),
            (34, 38, 1, 6),
            (42, 44, 1, 2),
            (39, 43, 2, 4),
            (40, 46, 3, 6),
        ],
    ),
]


@pytest.mark.parametrize(
    ("bases", "expected"),
    [(bases, expected) for _, bases, expected in _DROP_RULE_CASES],
    ids=[name for name, _, _ in _DROP_RULE_CASES],
)
def test_the_drop_rule_gives_pvcs_candidates_in_any_blocks(
    bases: str, expected: list[tuple[int, int, int, int]]
) -> None:
    parameters = StrFinderParameters(copy_floors=_ZERO_FLOORS)
    for block_bases in range(1, 61):
        got = _find_candidates(bases, parameters, block_bases=block_bases)
        found = [(c.start, c.end, c.period, c.score) for c in got]
        assert found == expected, block_bases


def test_a_candidate_under_its_copy_floor_is_not_returned() -> None:
    # Four A's then five, apart: a period-1 floor of 5 keeps only the five.
    bases = "GCTAGAAAACTGTCAGTCCGTAGTAAAAACTGTCAG"
    floors = (5, 0, 0, 0, 0, 0)
    parameters = StrFinderParameters(scan_max_period=1, copy_floors=floors)
    assert _find_candidates(bases, parameters) == [_Candidate(start=24, end=29, period=1, score=8)]


def _make_random_bases(seed: int, length: int, alphabet: str) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice(alphabet) for _ in range(length))


@settings(max_examples=300, deadline=None)
@given(
    bases=st.text(alphabet="ACGTN", max_size=300) | st.text(alphabet="AC", max_size=300),
    floors=st.tuples(*[st.integers(min_value=0, max_value=6)] * 6),
    block_bases=st.integers(min_value=1, max_value=40),
)
def test_the_copy_floors_and_the_blocks_change_nothing_but_what_is_under_the_floor(
    bases: str, floors: tuple[int, int, int, int, int, int], block_bases: int
) -> None:
    every = _find_candidates(bases, StrFinderParameters(copy_floors=_ZERO_FLOORS))
    over_floors = [
        candidate
        for candidate in every
        if (candidate.end - candidate.start) // candidate.period >= floors[candidate.period - 1]
    ]
    parameters = StrFinderParameters(copy_floors=floors)
    assert _find_candidates(bases, parameters, block_bases=block_bases) == over_floors


def test_the_candidates_come_by_period_then_start() -> None:
    bases = _make_random_bases(seed=7, length=3000, alphabet="ACGT")
    got = _find_candidates(bases, StrFinderParameters(copy_floors=_ZERO_FLOORS), block_bases=50)
    keys = [(candidate.period, candidate.start) for candidate in got]
    assert keys == sorted(keys)
    assert {candidate.period for candidate in got} == {1, 2, 3, 4, 5, 6}


@pytest.mark.slow
def test_step_one_holds_far_fewer_than_every_candidate_of_four_million_bases() -> None:
    # The plan's bound. Measured on CPython 3.15t, Apple M5 Pro, 19 October 2026: a peak of
    # 8.0 MB, 4 MB of it the bases as bytes; _find_candidates took 23 to 25 s under
    # tracemalloc, and this test 34 s. The same run in one block of the whole sequence,
    # which holds every candidate and every subsequence of the scores until the end,
    # peaked at 965.6 MB.
    bases = _make_random_bases(seed=19, length=4_000_000, alphabet="ACGT")
    tracemalloc.start()
    try:
        _find_candidates(bases, StrFinderParameters())
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 30_000_000


# ---- Steps 3 to 6, ported from pvc's segment_criteria.rs at bbcc9f4 ------------------------
# pvc's classify is _find_str_segments, steps 4 to 6; its prefilter is _drop_under_copy_floor,
# step 3. pvc's RepeatInterval is 0-based and half-open, as a candidate is, so its numbers
# carry over; pvc's SsrSegment is 1-based and inclusive, so its [s, e] is [s - 1, e) here.

# pvc's params() (segment_criteria.rs 1516), matched_params(0.8, 0, 5): the catalogue's
# periods and copy floors, with a bundle distance of 5 so that the fixtures can be small.
_PARAMS = StrFinderParameters(
    min_period=2,
    max_period=6,
    copy_floors=(10, 5, 4, 3, 3, 3),
    min_purity=0.8,
    min_score=0,
    bundle_distance=5,
)
# pvc's catalog_criteria() (segment_criteria.rs 1525): params() with a bundle distance of 50.
_CATALOG_CRITERIA = replace(_PARAMS, bundle_distance=50)
# What pvc's SsrSegmentCriteria::default() gives at bbcc9f4, which pvc_eval's defaults do not.
_PVC_DEFAULT = StrFinderParameters(
    min_period=1,
    max_period=6,
    copy_floors=(8, 6, 6, 6, 5, 4),
    min_purity=0.8,
    min_score=0,
    bundle_distance=15,
)

# A clean (AT)x8 with five bases of flank each side, the contig most of pvc's tests use.
_CLEAN_AT = "CGCGC" + "ATATATATATATATAT" + "GCGCG"


def _iv(start: int, end: int, period: int, score: int) -> _Candidate:
    """pvc's iv(): a candidate, 0-based and half-open."""
    return _Candidate(start=start, end=end, period=period, score=score)


def _find_whole_contig(
    candidates: list[_Candidate], contig: str, parameters: StrFinderParameters
) -> tuple[_StrSegment, ...]:
    """pvc's classify_whole_contig (segment_criteria.rs 1543): the STR segments alone."""
    return _find_str_segments(candidates, contig, parameters).str_segments


def _find_rejections(
    candidates: list[_Candidate], contig: str, parameters: StrFinderParameters
) -> list[_Rejection]:
    return [reason for _, reason in _find_str_segments(candidates, contig, parameters).rejected]


def test_a_bundle_member_is_trimmed_to_whole_copies_and_carries_its_motif() -> None:
    # pvc reads these bases from contig base 101 and expects the 1-based [103, 126]; from
    # base 0, as the port always reads, that is [2, 26). The trailing CA is cut.
    member = _make_bundle_member(_iv(2, 28, 3, 50), "TTcagcagcagcagcagcagcagcagCATT")
    assert member == _BundleMember(start=2, end=26, motif="CAG")

    # No two clean copies anywhere, so no cut: the detected span stands.
    member = _make_bundle_member(_iv(0, 12, 3, 5), "CAGTTTCAGAAA")
    assert member == _BundleMember(start=0, end=12, motif="CAG")


def test_locus_reports_its_span_period_and_purity() -> None:
    # pvc's [6, 21], 1-based, is 16 bases, 8 copies of AT. Without the contig name, which
    # an STR segment of the port does not hold.
    segment = _StrSegment(start=5, end=21, motif="AT", purity=1.0)
    assert (segment.start, segment.end) == (5, 21)
    assert segment.end - segment.start == 16
    assert len(segment.motif) == 2
    assert segment.purity == 1.0


def test_is_compound_flags_atat_but_not_at_or_atc() -> None:
    assert _is_compound("ATAT")
    assert _is_compound("ATATAT")
    assert not _is_compound("AT")
    assert not _is_compound("ATC")
    assert not _is_compound("ATCG")


def test_is_compound_spares_a_unit_dominated_by_one_base_in_every_rotation() -> None:
    for motif in ["AAAAAC", "AAAAAT", "AAAAAG", "ATTTTT", "CAAAAA", "AAAAC"]:
        for shift in range(len(motif)):
            rotated = motif[shift:] + motif[:shift]
            assert not _is_compound(rotated), rotated
    for motif in ["AA", "AAAA", "ACAC", "AACAAC", "ATATAT", "TATATA"]:
        assert _is_compound(motif), motif


def test_purity_is_one_for_perfect_tiling() -> None:
    assert _get_purity("ATATATAT", "AT") == 1.0
    assert _get_purity("CAGCAGCAG", "CAG") == 1.0


def test_purity_drops_with_one_interruption() -> None:
    assert _get_purity("ATATCTAT", "AT") == 7 / 8


def test_minimal_trim_strips_partial_motifs_at_both_ends() -> None:
    clean = "GATATATATATATC"
    cut = _find_cut(clean, "AT")
    assert cut is not None
    start, end = cut
    assert clean[start:end] == "ATATATATATAT"


def test_minimal_trim_returns_none_without_two_clean_copies() -> None:
    assert _find_cut("ACGT", "AT") is None


def test_split_bundles_sets_aside_whole_clusters_and_keeps_isolated() -> None:
    lone = _iv(5000, 5030, 2, 100)
    candidates = [_iv(100, 130, 2, 100), _iv(150, 180, 3, 100), _iv(200, 230, 2, 100), lone]
    isolated, bundled = _split_bundles(candidates, 50)
    assert isolated == [lone]
    assert [candidate.start for candidate in bundled] == [100, 150, 200]


def test_split_bundles_keeps_all_when_none_close() -> None:
    candidates = [_iv(100, 130, 2, 100), _iv(1000, 1030, 2, 100), _iv(2000, 2030, 2, 100)]
    isolated, bundled = _split_bundles(candidates, 50)
    assert isolated == candidates
    assert bundled == []


def test_classify_emits_a_clean_locus_one_based() -> None:
    # pvc's 1-based [6, 21] is [5, 21) here; without the contig name.
    segments = _find_whole_contig([_iv(5, 21, 2, 100)], _CLEAN_AT, _PARAMS)
    assert segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)


def test_classify_drops_period_over_six_and_compound() -> None:
    contig = "AAAACGATATATATATATATCCCC"
    too_long = _iv(0, 14, 7, 100)
    compound = _iv(6, 20, 4, 100)  # ATAT is (AT)x2
    assert _find_whole_contig([too_long, compound], contig, _PARAMS) == ()


def test_classify_drops_below_copy_number_floor() -> None:
    # (AT)x3 is under the period-2 floor of 5.
    assert _find_whole_contig([_iv(5, 11, 2, 100)], "CCCCCATATATCCCCC", _PARAMS) == ()


def test_classify_keeps_a_locus_whose_flank_clamps_to_one_base() -> None:
    # Step 6's own test of the contig ends: one base before the tract is enough.
    contig = "GATATATATATATATATCGCGCGCGCG"
    segments = _find_whole_contig([_iv(1, 17, 2, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].start == 1  # pvc's 1-based 2


def test_each_rejection_reason_names_the_gate_that_fired() -> None:
    # (AT)x3, under the period-2 floor of 5.
    contig = "CGCGCGCGCGCGCGCGCGCGATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_iv(20, 26, 2, 100)], contig, _PARAMS) == [_Rejection.COPY_FLOOR]

    # ATAT is AT twice.
    contig = "CGCGCGCGCGCGCGCGCGCGATATATATATATATATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_iv(20, 40, 4, 100)], contig, _PARAMS) == [_Rejection.COMPOUND]

    # No ATAT within the reach of the cut: this runs before the copy floor.
    contig = "CGCGCGCGCGCGCGCGCGCGATGGGGGGGGGATATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_iv(20, 39, 2, 100)], contig, _PARAMS) == [_Rejection.NO_CUT]

    # Clean ends and enough copies, but under a purity of 0.8.
    contig = "CGCGCGCGCGCGCGCGCGCGATATATATGGGGGGGGGGATATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_iv(20, 46, 2, 100)], contig, _PARAMS) == [_Rejection.PURITY]

    # The tract starts on the contig's first base.
    contig = "ATATATATATATATATCGCGC"
    assert _find_rejections([_iv(0, 16, 2, 100)], contig, _PARAMS) == [_Rejection.CONTIG_END]


def test_the_walk_reaches_only_one_of_classifications_five_gates() -> None:
    # Step 3 applies the copy floor before step 6 can.
    assert list(_drop_under_copy_floor([_iv(20, 26, 2, 100)], _PARAMS)) == []

    # A compound motif still fires the gate when it is given directly; step 1 never gives one.
    contig = "CGCGCGCGCGCGCGCGCGCGAAAAAAAAAAAAAAAAAAAACGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_iv(20, 40, 2, 100)], contig, _PARAMS) == [_Rejection.COMPOUND]


def test_classify_drops_locus_with_empty_left_flank() -> None:
    contig = "ATATATATATATATATCGCGC"
    assert _find_whole_contig([_iv(0, 16, 2, 100)], contig, _PARAMS) == ()


def test_classify_drops_locus_with_empty_right_flank() -> None:
    contig = "CGCGCATATATATATATATAT"
    assert _find_whole_contig([_iv(5, 21, 2, 100)], contig, _PARAMS) == ()


def test_classify_drops_period_one_homopolymer_and_spares_the_neighbour_ssr() -> None:
    contig = "A" * 20 + "CAG" * 10 + "TTTTT"
    segments = _find_whole_contig([_iv(0, 20, 1, 100), _iv(20, 50, 3, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].motif == "CAG"


def test_prefilter_drops_low_copy_noise() -> None:
    real = _iv(100, 120, 2, 100)
    noise = _iv(500, 504, 2, 10)
    assert list(_drop_under_copy_floor([real, noise], _PARAMS)) == [real]


def test_prefilter_keeps_a_higher_period_tract_that_is_not_a_multiple() -> None:
    two = _iv(100, 120, 2, 100)
    three = _iv(100, 121, 3, 100)
    assert list(_drop_under_copy_floor([two, three], _PARAMS)) == [two, three]


def test_classify_upper_cases_a_soft_masked_tracts_motif() -> None:
    # pvc also compares each result with its former production code (assert_agrees); that
    # code is not ported, so those assertions are dropped.
    contig = "cgcgc" + "atatatatatatatat" + "gcgcg"
    segments = _find_whole_contig([_iv(5, 21, 2, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].motif == "AT"
    assert segments[0].purity == 1.0

    mixed = "CGCGC" + "atatatatatatatat" + "GCGCG"
    assert _find_whole_contig([_iv(5, 21, 2, 100)], mixed, _PARAMS) == segments


def test_classify_reports_the_trimmed_start_not_the_detected_one() -> None:
    # pvc's assert_agrees, against its former production code, is dropped.
    contig = "CGCGC" + "ATC" + "ATATATATATAT" + "GCGCG"
    segments = _find_whole_contig([_iv(5, 20, 2, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].start == 8  # pvc's 1-based 9: the detected 5, plus 3 cut off
    assert segments[0].end == 20
    assert segments[0].end - segments[0].start == 12
    assert segments[0].motif == "AT"


def test_classify_keeps_an_imperfect_tract_above_the_floor_and_drops_one_below() -> None:
    # pvc's assert_agrees_at, against its former production code, is dropped.
    contig = "CGCGC" + "ATATATATATCTATAT" + "GCGCG"  # one C where an A belongs
    detected = _iv(5, 21, 2, 100)
    segments = _find_whole_contig([detected], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].purity == 15 / 16  # pvc's comment: 15/16 = 0.9375

    assert _find_whole_contig([detected], contig, replace(_PARAMS, min_purity=0.99)) == ()


def test_the_copy_floor_counts_the_trimmed_tract_not_the_detected_one() -> None:
    # Two clean CAGGA copies, then five bases that are not one: 3 copies detected, 2 cut.
    contig = "GCTAG" + "CAGGACAGGATTTTT" + "TCGAT"
    detected = _iv(5, 20, 5, 100)
    assert (detected.end - detected.start) // 5 == 3
    assert _find_whole_contig([detected], contig, _PARAMS) == ()

    lowered = replace(_PARAMS, copy_floors=(10, 5, 4, 3, 2, 3))
    kept = _find_whole_contig([detected], contig, lowered)
    assert len(kept) == 1
    assert kept[0].end - kept[0].start == 10


def test_classify_applies_the_score_gate_when_the_floor_actually_bites() -> None:
    # pvc's assert_agrees_at, against its former production code, is dropped.
    parameters = replace(_PARAMS, min_score=50)
    assert len(_find_whole_contig([_iv(5, 21, 2, 100)], _CLEAN_AT, parameters)) == 1
    assert _find_whole_contig([_iv(5, 21, 2, 49)], _CLEAN_AT, parameters) == ()
    assert len(_find_whole_contig([_iv(5, 21, 2, 50)], _CLEAN_AT, parameters)) == 1


def test_classify_at_default_never_rejects_a_score_the_scanner_can_emit() -> None:
    contig = "CGCGC" + "C" * 55 + "ATATATATATATATAT" + "G" * 64

    def count_at_default(score: int) -> int:
        return len(_find_whole_contig([_iv(60, 76, 2, score)], contig, _PVC_DEFAULT))

    assert count_at_default(1) == 1
    assert count_at_default(100) == 1
    assert count_at_default(0) == 1
    assert count_at_default(-1) == 0


def test_classify_rejects_a_zero_flank() -> None:
    # pvc panics in classify; the port refuses the parameters.
    with pytest.raises(ValueError, match=r"\bbundle_distance\b"):
        replace(_PARAMS, bundle_distance=0)


def test_classify_hands_back_the_bundle_members_it_sets_aside() -> None:
    # pvc's assert_agrees_at, against its former production code, is dropped.
    contig = "C" * 100 + "AT" * 10 + "C" * 30 + "AT" * 10 + "C" * 4830 + "AT" * 10 + "C" * 80
    first = _iv(100, 120, 2, 100)
    second = _iv(150, 170, 2, 100)
    lone = _iv(5000, 5020, 2, 100)
    classification = _find_str_segments([first, second, lone], contig, _CATALOG_CRITERIA)
    assert len(classification.str_segments) == 1
    assert classification.str_segments[0].start == 5000  # pvc's 1-based 5001
    assert classification.bundle_members == (first, second)


def test_classify_bundles_only_what_the_flank_test_rejects() -> None:
    classification = _find_str_segments([_iv(5, 21, 2, 100)], _CLEAN_AT, _PARAMS)
    assert len(classification.str_segments) == 1
    assert classification.bundle_members == ()

    classification = _find_str_segments([_iv(5, 11, 2, 100)], "CCCCCATATATCCCCC", _PARAMS)
    assert classification.str_segments == ()
    assert classification.bundle_members == ()


def test_narrowing_the_period_scope_changes_what_is_admitted() -> None:
    contig = "TTTTT" + "CAG" * 10 + "TTTTT"
    cag = _iv(5, 35, 3, 100)
    assert len(_find_whole_contig([cag], contig, _PARAMS)) == 1
    dinucleotides_only = replace(_PARAMS, min_period=2, max_period=2)
    assert _find_whole_contig([cag], contig, dinucleotides_only) == ()

    poly_a = "T" * 5 + "A" * 30 + "TTTTT"
    homopolymer = _iv(5, 35, 1, 100)
    assert _find_whole_contig([homopolymer], poly_a, _PARAMS) == ()
    with_homopolymers = replace(_PARAMS, min_period=1, max_period=6)
    assert len(_find_whole_contig([homopolymer], poly_a, with_homopolymers)) == 1


def test_classify_bundles_at_the_flank_radius() -> None:
    # pvc's two assert_agrees_at, against its former production code, are dropped.
    contig = "C" * 100 + "AT" * 10 + "C" * 30 + "AT" * 10 + "C" * 130
    pair = [_iv(100, 120, 2, 100), _iv(150, 170, 2, 100)]
    assert _find_whole_contig(pair, contig, replace(_PARAMS, bundle_distance=50)) == ()
    assert len(_find_whole_contig(pair, contig, replace(_PARAMS, bundle_distance=5))) == 2


def test_the_minimum_is_looked_up_by_the_tracts_own_period() -> None:
    contig = "C" * 20 + "AT" * 8 + "C" * 164 + "CAG" * 10 + "C" * 170
    dinucleotide = _iv(20, 36, 2, 100)
    trinucleotide = _iv(200, 230, 3, 100)
    picky = replace(_PARAMS, copy_floors=(10, 9, 4, 3, 3, 3))
    segments = _find_whole_contig([dinucleotide, trinucleotide], contig, picky)
    assert len(segments) == 1
    assert segments[0].motif == "CAG"

    mirrored = replace(_PARAMS, copy_floors=(10, 5, 11, 3, 3, 3))
    segments = _find_whole_contig([dinucleotide, trinucleotide], contig, mirrored)
    assert len(segments) == 1
    assert segments[0].motif == "AT"

    assert list(_drop_under_copy_floor([dinucleotide, trinucleotide], picky)) == [trinucleotide]
    assert list(_drop_under_copy_floor([dinucleotide, trinucleotide], mirrored)) == [dinucleotide]


def test_min_copies_uses_the_fallback_only_beyond_the_table() -> None:
    # pvc's table has a value for the periods past 6, which the port does not scan; its
    # assertions on periods 0, 7 and 255 are dropped.
    floors = replace(_PARAMS, copy_floors=(10, 5, 4, 3, 3, 3))
    assert _get_copy_floor(floors, 1) == 10
    assert _get_copy_floor(floors, 6) == 3

    uniform = replace(_PARAMS, copy_floors=(4, 4, 4, 4, 4, 4))
    assert _get_copy_floor(uniform, 1) == 4
    assert _get_copy_floor(uniform, 6) == 4


def test_classify_rejects_a_period_ceiling_no_motif_can_hold() -> None:
    # pvc panics in classify; the port refuses the parameters.
    with pytest.raises(ValueError, match=r"\bmax_period\b"):
        replace(_PARAMS, max_period=7)


def test_classify_rejects_a_nan_purity_floor() -> None:
    # pvc panics in classify; the port refuses the parameters.
    with pytest.raises(ValueError, match=r"\bmin_purity\b"):
        replace(_PARAMS, min_purity=float("nan"))


def test_raising_the_copy_floor_changes_what_is_admitted() -> None:
    tract = _iv(5, 21, 2, 100)
    assert len(_find_whole_contig([tract], _CLEAN_AT, _PARAMS)) == 1
    strict = replace(_PARAMS, copy_floors=(9, 9, 9, 9, 9, 9))
    assert _find_whole_contig([tract], _CLEAN_AT, strict) == ()

    assert list(_drop_under_copy_floor([tract], _PARAMS)) == [tract]
    assert list(_drop_under_copy_floor([tract], strict)) == []


def test_joins_cluster_is_strict_at_the_threshold() -> None:
    assert not _is_in_cluster(130, _iv(180, 210, 2, 100), 50)
    assert _is_in_cluster(130, _iv(179, 209, 2, 100), 50)


def test_a_tract_inside_the_cluster_always_joins_it() -> None:
    assert _is_in_cluster(1000, _iv(500, 520, 2, 2000), 30)
    assert _is_in_cluster(1000, _iv(1000, 1020, 2, 2000), 30)


def test_a_member_spanning_past_its_neighbour_keeps_the_cluster_open() -> None:
    candidates = [
        _iv(0, 43, 2, 100),
        _iv(56, 154, 4, 100),  # the long one
        _iv(62, 105, 2, 100),
        _iv(142, 154, 2, 100),  # 37 bases past the one before, inside the long one
    ]
    isolated, bundled = _split_bundles(candidates, 30)
    assert isolated == []
    assert bundled == candidates
    assert _make_clusters(bundled, 30) == [candidates]


def test_prefilter_drops_an_inverted_interval_instead_of_underflowing() -> None:
    inverted = _iv(120, 100, 2, 100)
    good = _iv(100, 120, 2, 100)
    assert list(_drop_under_copy_floor([inverted, good], _PARAMS)) == [good]
    assert list(_drop_under_copy_floor([_iv(100, 100, 2, 100)], _PARAMS)) == []


_SYNTHETIC_REF = Path(__file__).parent / "data" / "tandem_repeat" / "synthetic_ref.fa"


def test_differential_vs_production_build_loci_on_scanner_output() -> None:
    # pvc compares each contig's STR segments with its former production code
    # (assert_agrees), which is not ported; what is left checks that segments are found.
    contigs = 0
    str_segments = 0
    with IndexedFasta(_SYNTHETIC_REF) as reference:
        for contig in reference.contigs:
            bases = reference.get_bases(contig, 0, reference.get_length(contig))
            # pvc's find_tandem_repeats over periods 1 to 6, then its prefilter.
            scanned = StrFinderParameters(copy_floors=_ZERO_FLOORS)
            candidates = _drop_under_copy_floor(_find_candidates(bases, scanned), _PARAMS)
            str_segments += len(_find_whole_contig(list(candidates), bases, _PARAMS))
            contigs += 1
    assert contigs > 0
    assert str_segments > 0


# ---- Step 2, ported from pvc's repeat_catalog/row.rs at bbcc9f4 -----------------------------


def _make_contig_with_tract(left: int, motif: str, copies: int, right: int) -> str:
    """pvc's contig_with_tract (row.rs 296): G's, the copies of the motif, G's."""
    return "G" * left + motif * copies + "G" * right


@pytest.mark.parametrize(("left", "is_kept"), [(15, True), (14, False)], ids=["15", "14"])
def test_fifteen_bases_of_flank_is_kept_and_fourteen_is_not(left: int, is_kept: bool) -> None:
    bases = _make_contig_with_tract(left, "CAG", 6, 40)
    candidate = _iv(left, left + 18, 3, 42)
    kept = list(_drop_near_contig_ends([candidate], bases, StrFinderParameters()))
    assert kept == ([candidate] if is_kept else [])


@pytest.mark.parametrize(("right", "is_kept"), [(15, True), (14, False)], ids=["15", "14"])
def test_the_right_flank_is_measured_against_the_contigs_own_end(right: int, is_kept: bool) -> None:
    bases = _make_contig_with_tract(40, "CAG", 6, right)
    candidate = _iv(40, 58, 3, 42)
    kept = list(_drop_near_contig_ends([candidate], bases, StrFinderParameters()))
    assert kept == ([candidate] if is_kept else [])


# pvc_eval's own: step 2 measures the cut when there is one, and the detected span when
# there is none. Worked by hand from the spec's step 2 and step 6: ATC(AT)x8 at 13 has the
# motif AT, whose first two clean copies are at offset 3, so its cut starts at 16 while its
# detected span starts at 13; CAGTTTCAGAAA has no two clean copies, so its detected span is
# measured.
_STEP_TWO_CASES = [
    ("cut_clear_of_the_end", "G" * 13 + "ATC" + "AT" * 8 + "G" * 40, _iv(13, 32, 2, 42), True),
    ("no_cut_14_from_the_end", "G" * 14 + "CAGTTTCAGAAA" + "G" * 40, _iv(14, 26, 3, 5), False),
    ("no_cut_15_from_the_end", "G" * 15 + "CAGTTTCAGAAA" + "G" * 40, _iv(15, 27, 3, 5), True),
]


@pytest.mark.parametrize(
    ("bases", "candidate", "is_kept"),
    [case[1:] for case in _STEP_TWO_CASES],
    ids=[case[0] for case in _STEP_TWO_CASES],
)
def test_the_contig_end_flank_is_measured_on_the_cut_or_else_the_detected_span(
    bases: str, candidate: _Candidate, is_kept: bool
) -> None:
    kept = list(_drop_near_contig_ends([candidate], bases, StrFinderParameters()))
    assert kept == ([candidate] if is_kept else [])


# ---- Steps 4 to 6, pvc_eval's own -------------------------------------------------------------
# Boundaries that no ported test reaches. Each expected value is worked by hand from the
# spec's steps, as the comment beside it shows.


def test_a_tract_exactly_two_copies_long_has_no_cut() -> None:
    # Step 6: the search for the start fails before it compares as soon as s + 2m >= L,
    # which at s = 0 already holds for a tract of two copies.
    assert _find_cut("ATAT", "AT") is None


# The cut's start is searched up to offset 3m and its end down to L - 3m: for AT, 6.
_CUT_REACH_CASES = [
    ("start_at_offset_6", "ATCCCC" + "AT" * 6, (6, 18)),
    ("start_past_offset_6", "ATCCCCC" + "AT" * 6, None),
    ("end_at_6_from_the_end", "AT" * 6 + "G" * 6, (0, 12)),
    ("end_past_6_from_the_end", "AT" * 6 + "G" * 7, None),
]


@pytest.mark.parametrize(
    ("tract", "cut"),
    [case[1:] for case in _CUT_REACH_CASES],
    ids=[case[0] for case in _CUT_REACH_CASES],
)
def test_the_cut_reaches_three_motifs_into_the_tract(
    tract: str, cut: tuple[int, int] | None
) -> None:
    assert _find_cut(tract, "AT") == cut


def test_a_cut_at_exactly_the_minimum_purity_is_kept() -> None:
    # The cut ATATGGATAT matches AT at 8 of its 10 bases: 0.8, not under 0.8.
    contig = "CGCGC" + "ATATGGATAT" + "GCGCG"
    segments = _find_whole_contig([_iv(5, 15, 2, 100)], contig, _PARAMS)
    assert segments == (_StrSegment(start=5, end=15, motif="AT", purity=0.8),)


def test_a_candidate_shorter_than_its_period_is_dropped_before_bundling() -> None:
    # Step 4 drops the one-base candidate at 22, so it does not bundle the (AT)x8 at 5.
    short = _iv(22, 23, 2, 1)
    classification = _find_str_segments([_iv(5, 21, 2, 100), short], _CLEAN_AT, _PARAMS)
    assert classification.str_segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)
    assert classification.rejected == ((short, _Rejection.COPY_FLOOR),)


def test_the_candidates_are_bundled_in_order_of_start() -> None:
    # Given by period, as step 1 gives them: the (CAG)x10 at 60 comes after the (AT)x10 at
    # 100 and at 200. Sorted by start, it bundles with the one at 100, 10 bases on; the one
    # at 200 stands alone.
    contig = "C" * 60 + "CAG" * 10 + "C" * 10 + "AT" * 10 + "C" * 80 + "AT" * 10 + "C" * 80
    at_100 = _iv(100, 120, 2, 100)
    at_200 = _iv(200, 220, 2, 100)
    cag_60 = _iv(60, 90, 3, 100)
    classification = _find_str_segments([at_100, at_200, cag_60], contig, _PVC_DEFAULT)
    assert classification.str_segments == (_StrSegment(start=200, end=220, motif="AT", purity=1.0),)
    assert classification.bundle_members == (cag_60, at_100)


# ---- The order of the rules of steps 4 to 6, pvc_eval's own -----------------------------------
# Each case has two rules that apply, so that only the order the spec gives passes. The
# expected values are pvc's, from a transcription of its classify, run as
# uv run python /private/tmp/claude-501/-Users-jose-devel-pvc-eval/2dab6280-84d0-4e89-88dc-7955494ebca9/scratchpad/m2_fix_values.py \
#     /private/tmp/claude-501/-Users-jose-devel-pvc-eval/2dab6280-84d0-4e89-88dc-7955494ebca9/scratchpad/review_m2_spec


def test_a_candidate_under_min_score_bundles_nothing() -> None:
    # The (CAG)x5 at 22 is 1 base from the (AT)x8 at 5, but under min_score 50.
    contig = "CGCGC" + "AT" * 8 + "G" + "CAG" * 5 + "CGCGC"
    candidates = [_iv(5, 21, 2, 100), _iv(22, 37, 3, 10)]
    classification = _find_str_segments(candidates, contig, replace(_PARAMS, min_score=50))
    assert classification.str_segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)
    assert classification.bundle_members == ()
    assert classification.rejected == ()


def test_a_compound_candidate_bundles_nothing() -> None:
    # AAAAAAAA at 22, given as period 2, is compound and is dropped before bundling.
    contig = "CGCGC" + "AT" * 8 + "G" + "A" * 8 + "CGCGC"
    compound = _iv(22, 30, 2, 100)
    classification = _find_str_segments([_iv(5, 21, 2, 100), compound], contig, _PARAMS)
    assert classification.str_segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)
    assert classification.bundle_members == ()
    assert classification.rejected == ((compound, _Rejection.COMPOUND),)


def test_a_candidate_shorter_than_its_period_is_turned_down_before_the_compound_test() -> None:
    # AA given as period 4 is both shorter than its period and, as a motif, compound.
    rejections = _find_rejections([_iv(5, 7, 4, 100)], "CGCGCAAGCGCG", _PARAMS)
    assert rejections == [_Rejection.COPY_FLOOR]


def test_a_candidate_exactly_one_period_long_still_bundles() -> None:
    # CG at 22, one copy of period 2, is not shorter than its period.
    one_copy = _iv(22, 24, 2, 1)
    classification = _find_str_segments([_iv(5, 21, 2, 100), one_copy], _CLEAN_AT, _PARAMS)
    assert classification.str_segments == ()
    assert classification.bundle_members == (_iv(5, 21, 2, 100), one_copy)


def test_step_six_tests_the_copy_floor_then_the_purity_then_the_contig_ends() -> None:
    # The cut ATATGGATAT has 5 copies, under a floor of 6, and a purity of 0.8, under 0.9.
    contig = "CGCGC" + "ATATGGATAT" + "GCGCG"
    parameters = replace(_PARAMS, copy_floors=(10, 6, 4, 3, 3, 3), min_purity=0.9)
    assert _find_rejections([_iv(5, 15, 2, 100)], contig, parameters) == [_Rejection.COPY_FLOOR]

    # The cut ATATGGATATAT has a purity of 10/12, under 0.9, and starts on base 0.
    contig = "ATATGGATATAT" + "GCGCG"
    parameters = replace(_PARAMS, min_purity=0.9)
    assert _find_rejections([_iv(0, 12, 2, 100)], contig, parameters) == [_Rejection.PURITY]
