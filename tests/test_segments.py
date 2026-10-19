import random
import tracemalloc
from fractions import Fraction

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pvc_eval.segments import (
    GenomeSegment,
    RegionType,
    StrFinderParameters,
    _Candidate,
    _find_candidates,
    _find_maximal_scoring_subsequences,
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
    # The plan's bound. Measured with tracemalloc on CPython 3.15t, 19 October 2026: a peak
    # of 8.0 MB, 4 MB of it the bases as bytes, in 14 s. The same run in one block of the
    # whole sequence, which holds every candidate and every subsequence of the scores until
    # the end, peaked at 804.5 MB.
    bases = _make_random_bases(seed=19, length=4_000_000, alphabet="ACGT")
    tracemalloc.start()
    try:
        _find_candidates(bases, StrFinderParameters())
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 30_000_000
