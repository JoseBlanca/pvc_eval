import gzip
import itertools
import random
import sys
import tracemalloc
from collections.abc import Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import (
    _CANONICAL,
    GenomeSegment,
    RegionType,
    StrFinderParameters,
    _BundleMember,
    _Candidate,
    _drop_near_contig_ends,
    _drop_under_copy_floor,
    _find_bundles,
    _find_candidates,
    _find_cut,
    _find_maximal_scoring_subsequences,
    _find_requested_str_regions,
    _find_segments,
    _find_str_segments,
    _get_copy_floor,
    _get_purity,
    _is_compound,
    _is_in_cluster,
    _make_bundle_member,
    _make_clusters,
    _make_grown_satellites,
    _make_repeat_runs,
    _make_runs,
    _make_satellites_absorbing,
    _make_score_blocks,
    _Rejection,
    _Segment,
    _StrSegment,
    segment_contig,
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
        ValueError, match=r"^min_purity is Fraction\(4, 5\); it must be a float or an int$"
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
    """pvc's rt_segments (tandem_repeat.rs 1208), the scores fed in blocks of block_length,
    each block as its runs of scores of one sign, (their sum, how many)."""
    blocks = [
        [
            (sum(run), len(run))
            for run in (list(group) for _, group in itertools.groupby(block, key=lambda s: s > 0))
        ]
        for block in (scores[i : i + block_length] for i in range(0, len(scores), block_length))
    ]
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


# Scored by hand at the defaults, a match 2 and a mismatch -7; score k is the base at k + period.
_SCORE_BLOCK_CASES = [
    # A=A, then N against A, N against N twice, A against N, then A=A.
    ("n_never_matches_n", "AANNNAA", 1, 100, [[(2, 1), (-28, 4), (2, 1)]]),
    # Blocks of 4 bases end after base 3 and base 7: the scores 0 to 2, then 3 to 6.
    ("a_block_end_cuts_a_run_of_matches", "AAAAAAAA", 1, 4, [[(6, 3)], [(8, 4)]]),
    # c=C, a=A, C=c, then T against a.
    ("lower_case_matches_upper_case", "CAcaCT", 2, 100, [[(6, 3), (-7, 1)]]),
]


@pytest.mark.parametrize(
    ("bases", "period", "block_bases", "expected"),
    [case[1:] for case in _SCORE_BLOCK_CASES],
    ids=[case[0] for case in _SCORE_BLOCK_CASES],
)
def test_the_scores_come_as_runs_of_one_sign_cut_at_the_blocks(
    bases: str, period: int, block_bases: int, expected: list[list[tuple[int, int]]]
) -> None:
    canonical = bases.encode("ascii").translate(_CANONICAL)
    got = _make_score_blocks(canonical, period, StrFinderParameters(), block_bases)
    assert list(got) == expected


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
# expected candidates, (start, end, period, score) with no copy floor, are pvc's: its code at
# bbcc9f4 (tandem_repeat.rs `find_tandem_repeats`, `maximal_scoring_subsequences` and
# `is_primitive_motif`), run over periods 1 to 6 through a line-by-line Python transcription
# on 19 October 2026. Milestone 5 checks the whole finder against pvc's own binary.
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


def _make_candidate(start: int, end: int, period: int, score: int) -> _Candidate:
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
    member = _make_bundle_member(_make_candidate(2, 28, 3, 50), "TTcagcagcagcagcagcagcagcagCATT")
    assert member == _BundleMember(start=2, end=26, motif="CAG")

    # No two clean copies anywhere, so no cut: the detected span stands.
    member = _make_bundle_member(_make_candidate(0, 12, 3, 5), "CAGTTTCAGAAA")
    assert member == _BundleMember(start=0, end=12, motif="CAG")


def test_locus_reports_its_span_period_and_purity() -> None:
    # pvc's [6, 21], 1-based, is 16 bases, 8 copies of AT. Without the contig name, which
    # an STR segment of the port does not hold. It checks only that the record gives back
    # the fields it was made with; pvc's test checks its accessors, which the port has not.
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
    lone = _make_candidate(5000, 5030, 2, 100)
    candidates = [
        _make_candidate(100, 130, 2, 100),
        _make_candidate(150, 180, 3, 100),
        _make_candidate(200, 230, 2, 100),
        lone,
    ]
    isolated, bundled = _find_bundles(candidates, 50)
    assert isolated == [lone]
    assert [candidate.start for candidate in bundled] == [100, 150, 200]


def test_split_bundles_keeps_all_when_none_close() -> None:
    candidates = [
        _make_candidate(100, 130, 2, 100),
        _make_candidate(1000, 1030, 2, 100),
        _make_candidate(2000, 2030, 2, 100),
    ]
    isolated, bundled = _find_bundles(candidates, 50)
    assert isolated == candidates
    assert bundled == []


def test_classify_emits_a_clean_locus_one_based() -> None:
    # pvc's 1-based [6, 21] is [5, 21) here; without the contig name.
    segments = _find_whole_contig([_make_candidate(5, 21, 2, 100)], _CLEAN_AT, _PARAMS)
    assert segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)


def test_classify_drops_period_over_six_and_compound() -> None:
    contig = "AAAACGATATATATATATATCCCC"
    too_long = _make_candidate(0, 14, 7, 100)
    compound = _make_candidate(6, 20, 4, 100)  # ATAT is (AT)x2
    assert _find_whole_contig([too_long, compound], contig, _PARAMS) == ()


def test_classify_drops_below_copy_number_floor() -> None:
    # (AT)x3 is under the period-2 floor of 5.
    assert _find_whole_contig([_make_candidate(5, 11, 2, 100)], "CCCCCATATATCCCCC", _PARAMS) == ()


def test_classify_keeps_a_locus_whose_flank_clamps_to_one_base() -> None:
    # Step 6's own test of the contig ends: one base before the repeat is enough.
    contig = "GATATATATATATATATCGCGCGCGCG"
    segments = _find_whole_contig([_make_candidate(1, 17, 2, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].start == 1  # pvc's 1-based 2


def test_each_rejection_reason_names_the_gate_that_fired() -> None:
    # (AT)x3, under the period-2 floor of 5.
    contig = "CGCGCGCGCGCGCGCGCGCGATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_make_candidate(20, 26, 2, 100)], contig, _PARAMS) == [
        _Rejection.COPY_FLOOR
    ]

    # ATAT is AT twice.
    contig = "CGCGCGCGCGCGCGCGCGCGATATATATATATATATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_make_candidate(20, 40, 4, 100)], contig, _PARAMS) == [
        _Rejection.COMPOUND
    ]

    # No ATAT within the reach of the cut: this runs before the copy floor.
    contig = "CGCGCGCGCGCGCGCGCGCGATGGGGGGGGGATATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_make_candidate(20, 39, 2, 100)], contig, _PARAMS) == [
        _Rejection.NO_CUT
    ]

    # Clean ends and enough copies, but under a purity of 0.8.
    contig = "CGCGCGCGCGCGCGCGCGCGATATATATGGGGGGGGGGATATATATCGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_make_candidate(20, 46, 2, 100)], contig, _PARAMS) == [
        _Rejection.PURITY
    ]

    # The repeat starts on the contig's first base.
    contig = "ATATATATATATATATCGCGC"
    assert _find_rejections([_make_candidate(0, 16, 2, 100)], contig, _PARAMS) == [
        _Rejection.CONTIG_END
    ]


def test_the_walk_reaches_only_one_of_classifications_five_gates() -> None:
    # Step 3 applies the copy floor before step 6 can.
    assert list(_drop_under_copy_floor([_make_candidate(20, 26, 2, 100)], _PARAMS)) == []

    # A compound motif still fires the gate when it is given directly; step 1 never gives one.
    contig = "CGCGCGCGCGCGCGCGCGCGAAAAAAAAAAAAAAAAAAAACGCGCGCGCGCGCGCGCGCG"
    assert _find_rejections([_make_candidate(20, 40, 2, 100)], contig, _PARAMS) == [
        _Rejection.COMPOUND
    ]


def test_classify_drops_locus_with_empty_left_flank() -> None:
    contig = "ATATATATATATATATCGCGC"
    assert _find_whole_contig([_make_candidate(0, 16, 2, 100)], contig, _PARAMS) == ()


def test_classify_drops_locus_with_empty_right_flank() -> None:
    contig = "CGCGCATATATATATATATAT"
    assert _find_whole_contig([_make_candidate(5, 21, 2, 100)], contig, _PARAMS) == ()


def test_classify_drops_period_one_homopolymer_and_spares_the_neighbour_ssr() -> None:
    contig = "A" * 20 + "CAG" * 10 + "TTTTT"
    segments = _find_whole_contig(
        [_make_candidate(0, 20, 1, 100), _make_candidate(20, 50, 3, 100)], contig, _PARAMS
    )
    assert len(segments) == 1
    assert segments[0].motif == "CAG"
    assert len(segments[0].motif) == 3  # pvc's period()


def test_prefilter_drops_low_copy_noise() -> None:
    real = _make_candidate(100, 120, 2, 100)
    noise = _make_candidate(500, 504, 2, 10)
    assert list(_drop_under_copy_floor([real, noise], _PARAMS)) == [real]


def test_prefilter_keeps_a_higher_period_tract_that_is_not_a_multiple() -> None:
    two = _make_candidate(100, 120, 2, 100)
    three = _make_candidate(100, 121, 3, 100)
    assert list(_drop_under_copy_floor([two, three], _PARAMS)) == [two, three]


def test_classify_upper_cases_a_soft_masked_tracts_motif() -> None:
    # pvc also compares each result with its former production code (assert_agrees); that
    # code is not ported, so those assertions are dropped.
    contig = "cgcgc" + "atatatatatatatat" + "gcgcg"
    segments = _find_whole_contig([_make_candidate(5, 21, 2, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].motif == "AT"
    assert segments[0].purity == 1.0

    mixed = "CGCGC" + "atatatatatatatat" + "GCGCG"
    assert _find_whole_contig([_make_candidate(5, 21, 2, 100)], mixed, _PARAMS) == segments


def test_classify_reports_the_trimmed_start_not_the_detected_one() -> None:
    # pvc's assert_agrees, against its former production code, is dropped.
    contig = "CGCGC" + "ATC" + "ATATATATATAT" + "GCGCG"
    segments = _find_whole_contig([_make_candidate(5, 20, 2, 100)], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].start == 8  # pvc's 1-based 9: the detected 5, plus 3 cut off
    assert segments[0].end == 20
    assert segments[0].end - segments[0].start == 12
    assert segments[0].motif == "AT"


def test_classify_keeps_an_imperfect_tract_above_the_floor_and_drops_one_below() -> None:
    # pvc's assert_agrees_at, against its former production code, is dropped.
    contig = "CGCGC" + "ATATATATATCTATAT" + "GCGCG"  # one C where an A belongs
    detected = _make_candidate(5, 21, 2, 100)
    segments = _find_whole_contig([detected], contig, _PARAMS)
    assert len(segments) == 1
    assert segments[0].purity == 15 / 16  # pvc's comment: 15/16 = 0.9375

    assert _find_whole_contig([detected], contig, replace(_PARAMS, min_purity=0.99)) == ()


def test_the_copy_floor_counts_the_trimmed_tract_not_the_detected_one() -> None:
    # Two clean CAGGA copies, then five bases that are not one: 3 copies detected, 2 cut.
    contig = "GCTAG" + "CAGGACAGGATTTTT" + "TCGAT"
    detected = _make_candidate(5, 20, 5, 100)
    assert (detected.end - detected.start) // 5 == 3
    assert _find_whole_contig([detected], contig, _PARAMS) == ()

    lowered = replace(_PARAMS, copy_floors=(10, 5, 4, 3, 2, 3))
    kept = _find_whole_contig([detected], contig, lowered)
    assert len(kept) == 1
    assert kept[0].end - kept[0].start == 10


def test_classify_applies_the_score_gate_when_the_floor_actually_bites() -> None:
    # pvc's assert_agrees_at, against its former production code, is dropped.
    parameters = replace(_PARAMS, min_score=50)
    assert len(_find_whole_contig([_make_candidate(5, 21, 2, 100)], _CLEAN_AT, parameters)) == 1
    assert _find_whole_contig([_make_candidate(5, 21, 2, 49)], _CLEAN_AT, parameters) == ()
    assert len(_find_whole_contig([_make_candidate(5, 21, 2, 50)], _CLEAN_AT, parameters)) == 1


def test_classify_at_default_never_rejects_a_score_the_scanner_can_emit() -> None:
    contig = "CGCGC" + "C" * 55 + "ATATATATATATATAT" + "G" * 64

    def count_at_default(score: int) -> int:
        return len(_find_whole_contig([_make_candidate(60, 76, 2, score)], contig, _PVC_DEFAULT))

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
    first = _make_candidate(100, 120, 2, 100)
    second = _make_candidate(150, 170, 2, 100)
    lone = _make_candidate(5000, 5020, 2, 100)
    classification = _find_str_segments([first, second, lone], contig, _CATALOG_CRITERIA)
    assert len(classification.str_segments) == 1
    assert classification.str_segments[0].start == 5000  # pvc's 1-based 5001
    assert classification.bundle_members == (first, second)


def test_classify_bundles_only_what_the_flank_test_rejects() -> None:
    classification = _find_str_segments([_make_candidate(5, 21, 2, 100)], _CLEAN_AT, _PARAMS)
    assert len(classification.str_segments) == 1
    assert classification.bundle_members == ()

    classification = _find_str_segments(
        [_make_candidate(5, 11, 2, 100)], "CCCCCATATATCCCCC", _PARAMS
    )
    assert classification.str_segments == ()
    assert classification.bundle_members == ()


def test_narrowing_the_period_scope_changes_what_is_admitted() -> None:
    contig = "TTTTT" + "CAG" * 10 + "TTTTT"
    cag = _make_candidate(5, 35, 3, 100)
    assert len(_find_whole_contig([cag], contig, _PARAMS)) == 1
    dinucleotides_only = replace(_PARAMS, min_period=2, max_period=2)
    assert _find_whole_contig([cag], contig, dinucleotides_only) == ()

    poly_a = "T" * 5 + "A" * 30 + "TTTTT"
    homopolymer = _make_candidate(5, 35, 1, 100)
    assert _find_whole_contig([homopolymer], poly_a, _PARAMS) == ()
    with_homopolymers = replace(_PARAMS, min_period=1, max_period=6)
    assert len(_find_whole_contig([homopolymer], poly_a, with_homopolymers)) == 1


def test_classify_bundles_at_the_flank_radius() -> None:
    # pvc's two assert_agrees_at, against its former production code, are dropped.
    contig = "C" * 100 + "AT" * 10 + "C" * 30 + "AT" * 10 + "C" * 130
    pair = [_make_candidate(100, 120, 2, 100), _make_candidate(150, 170, 2, 100)]
    assert _find_whole_contig(pair, contig, replace(_PARAMS, bundle_distance=50)) == ()
    assert len(_find_whole_contig(pair, contig, replace(_PARAMS, bundle_distance=5))) == 2


def test_the_minimum_is_looked_up_by_the_tracts_own_period() -> None:
    contig = "C" * 20 + "AT" * 8 + "C" * 164 + "CAG" * 10 + "C" * 170
    dinucleotide = _make_candidate(20, 36, 2, 100)
    trinucleotide = _make_candidate(200, 230, 3, 100)
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
    candidate = _make_candidate(5, 21, 2, 100)
    assert len(_find_whole_contig([candidate], _CLEAN_AT, _PARAMS)) == 1
    strict = replace(_PARAMS, copy_floors=(9, 9, 9, 9, 9, 9))
    assert _find_whole_contig([candidate], _CLEAN_AT, strict) == ()

    assert list(_drop_under_copy_floor([candidate], _PARAMS)) == [candidate]
    assert list(_drop_under_copy_floor([candidate], strict)) == []


def test_joins_cluster_is_strict_at_the_threshold() -> None:
    assert not _is_in_cluster(130, _make_candidate(180, 210, 2, 100), 50)
    assert _is_in_cluster(130, _make_candidate(179, 209, 2, 100), 50)


def test_a_tract_inside_the_cluster_always_joins_it() -> None:
    assert _is_in_cluster(1000, _make_candidate(500, 520, 2, 2000), 30)
    assert _is_in_cluster(1000, _make_candidate(1000, 1020, 2, 2000), 30)


def test_a_member_spanning_past_its_neighbour_keeps_the_cluster_open() -> None:
    candidates = [
        _make_candidate(0, 43, 2, 100),
        _make_candidate(56, 154, 4, 100),  # the long one
        _make_candidate(62, 105, 2, 100),
        _make_candidate(142, 154, 2, 100),  # 37 bases past the one before, inside the long one
    ]
    isolated, bundled = _find_bundles(candidates, 30)
    assert isolated == []
    assert bundled == candidates
    assert _make_clusters(bundled, 30) == [candidates]


def test_prefilter_drops_an_inverted_interval_instead_of_underflowing() -> None:
    inverted = _make_candidate(120, 100, 2, 100)
    good = _make_candidate(100, 120, 2, 100)
    assert list(_drop_under_copy_floor([inverted, good], _PARAMS)) == [good]
    assert list(_drop_under_copy_floor([_make_candidate(100, 100, 2, 100)], _PARAMS)) == []


_SYNTHETIC_REF = Path(__file__).parent / "data" / "tandem_repeat" / "synthetic_ref.fa"


def test_differential_vs_production_build_loci_on_scanner_output() -> None:
    # pvc compares each contig's STR segments with its former production code
    # (assert_agrees), which is not ported. What is left counts the STR segments and the
    # bundle members of each contig, against pvc's: the port's step 1, then pvc's code at
    # bbcc9f4 (segment_criteria.rs `prefilter` and `classify`, with `split_bundles`,
    # `finish_locus`, `minimal_trim`, `recompute_purity` and `is_compound`), run through a
    # line-by-line Python transcription on 19 October 2026. Milestone 5 checks the whole
    # finder against pvc's own binary.
    counts = {}
    with IndexedFasta(_SYNTHETIC_REF) as reference:
        for contig in reference.contigs:
            bases = reference.get_bases(contig, 0, reference.get_length(contig))
            # pvc's find_tandem_repeats over periods 1 to 6, then its prefilter.
            scanned = StrFinderParameters(copy_floors=_ZERO_FLOORS)
            candidates = _drop_under_copy_floor(_find_candidates(bases, scanned), _PARAMS)
            classification = _find_str_segments(candidates, bases, _PARAMS)
            counts[contig] = (
                len(classification.str_segments),
                len(classification.bundle_members),
            )
    assert counts == {"ctg1": (15, 0), "ctg2": (2, 0)}


# ---- Step 2, ported from pvc's repeat_catalog/row.rs at bbcc9f4 -----------------------------


def _make_contig_with_repeat(left: int, motif: str, copies: int, right: int) -> str:
    """pvc's contig_with_tract (row.rs 296): G's, the copies of the motif, G's."""
    return "G" * left + motif * copies + "G" * right


@pytest.mark.parametrize(("left", "is_kept"), [(15, True), (14, False)], ids=["15", "14"])
def test_fifteen_bases_of_flank_is_kept_and_fourteen_is_not(left: int, is_kept: bool) -> None:
    bases = _make_contig_with_repeat(left, "CAG", 6, 40)
    candidate = _make_candidate(left, left + 18, 3, 42)
    kept = list(_drop_near_contig_ends([candidate], bases, StrFinderParameters()))
    assert kept == ([candidate] if is_kept else [])


@pytest.mark.parametrize(("right", "is_kept"), [(15, True), (14, False)], ids=["15", "14"])
def test_the_right_flank_is_measured_against_the_contigs_own_end(right: int, is_kept: bool) -> None:
    bases = _make_contig_with_repeat(40, "CAG", 6, right)
    candidate = _make_candidate(40, 58, 3, 42)
    kept = list(_drop_near_contig_ends([candidate], bases, StrFinderParameters()))
    assert kept == ([candidate] if is_kept else [])


# pvc_eval's own: step 2 measures the cut when there is one, and the detected span when
# there is none. Worked by hand from the spec's step 2 and step 6: ATC(AT)x8 at 13 has the
# motif AT, whose first two clean copies are at offset 3, so its cut starts at 16 while its
# detected span starts at 13; CAGTTTCAGAAA has no two clean copies, so its detected span is
# measured.
_STEP_TWO_CASES = [
    (
        "cut_clear_of_the_end",
        "G" * 13 + "ATC" + "AT" * 8 + "G" * 40,
        _make_candidate(13, 32, 2, 42),
        True,
    ),
    (
        "no_cut_14_from_the_end",
        "G" * 14 + "CAGTTTCAGAAA" + "G" * 40,
        _make_candidate(14, 26, 3, 5),
        False,
    ),
    (
        "no_cut_15_from_the_end",
        "G" * 15 + "CAGTTTCAGAAA" + "G" * 40,
        _make_candidate(15, 27, 3, 5),
        True,
    ),
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
    # which at s = 0 already holds for a span of two copies.
    assert _find_cut("ATAT", "AT") is None


# The cut's start is searched up to offset 3m and its end down to L - 3m: for AT, 6.
_CUT_REACH_CASES = [
    ("start_at_offset_6", "ATCCCC" + "AT" * 6, (6, 18)),
    ("start_past_offset_6", "ATCCCCC" + "AT" * 6, None),
    ("end_at_6_from_the_end", "AT" * 6 + "G" * 6, (0, 12)),
    ("end_past_6_from_the_end", "AT" * 6 + "G" * 7, None),
]


@pytest.mark.parametrize(
    ("span_bases", "cut"),
    [case[1:] for case in _CUT_REACH_CASES],
    ids=[case[0] for case in _CUT_REACH_CASES],
)
def test_the_cut_reaches_three_motifs_into_the_tract(
    span_bases: str, cut: tuple[int, int] | None
) -> None:
    assert _find_cut(span_bases, "AT") == cut


def test_a_cut_at_exactly_the_minimum_purity_is_kept() -> None:
    # The cut ATATGGATAT matches AT at 8 of its 10 bases: 0.8, not under 0.8.
    contig = "CGCGC" + "ATATGGATAT" + "GCGCG"
    segments = _find_whole_contig([_make_candidate(5, 15, 2, 100)], contig, _PARAMS)
    assert segments == (_StrSegment(start=5, end=15, motif="AT", purity=0.8),)


def test_a_candidate_shorter_than_its_period_is_dropped_before_bundling() -> None:
    # Step 4 drops the one-base candidate at 22, so it does not bundle the (AT)x8 at 5.
    short = _make_candidate(22, 23, 2, 1)
    classification = _find_str_segments([_make_candidate(5, 21, 2, 100), short], _CLEAN_AT, _PARAMS)
    assert classification.str_segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)
    assert classification.rejected == ((short, _Rejection.COPY_FLOOR),)


def test_the_candidates_are_bundled_in_order_of_start() -> None:
    # Given by period, as step 1 gives them: the (CAG)x10 at 60 comes after the (AT)x10 at
    # 100 and at 200. Sorted by start, it bundles with the one at 100, 10 bases on; the one
    # at 200 stands alone.
    contig = "C" * 60 + "CAG" * 10 + "C" * 10 + "AT" * 10 + "C" * 80 + "AT" * 10 + "C" * 80
    at_100 = _make_candidate(100, 120, 2, 100)
    at_200 = _make_candidate(200, 220, 2, 100)
    cag_60 = _make_candidate(60, 90, 3, 100)
    classification = _find_str_segments([at_100, at_200, cag_60], contig, _PVC_DEFAULT)
    assert classification.str_segments == (_StrSegment(start=200, end=220, motif="AT", purity=1.0),)
    assert classification.bundle_members == (cag_60, at_100)


# ---- The order of the rules of steps 4 to 6, pvc_eval's own -----------------------------------
# Each case has two rules that apply, so that only the order the spec gives passes. The
# expected values here and in the tests below are pvc's: its code at bbcc9f4
# (segment_criteria.rs `prefilter`, `classify`, `split_bundles`, `joins_cluster`,
# `finish_locus`, `minimal_trim`, `recompute_purity`, `is_compound` and
# `bundle_member_from_bases`, and repeat_catalog/row.rs `row_for_interval`), run through a
# line-by-line Python transcription on 19 October 2026. Milestone 5 checks the whole finder
# against pvc's own binary.


def test_a_candidate_under_min_score_bundles_nothing() -> None:
    # The (CAG)x5 at 22 is 1 base from the (AT)x8 at 5, but under min_score 50.
    contig = "CGCGC" + "AT" * 8 + "G" + "CAG" * 5 + "CGCGC"
    candidates = [_make_candidate(5, 21, 2, 100), _make_candidate(22, 37, 3, 10)]
    classification = _find_str_segments(candidates, contig, replace(_PARAMS, min_score=50))
    assert classification.str_segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)
    assert classification.bundle_members == ()
    assert classification.rejected == ()


def test_a_compound_candidate_bundles_nothing() -> None:
    # AAAAAAAA at 22, given as period 2, is compound and is dropped before bundling.
    contig = "CGCGC" + "AT" * 8 + "G" + "A" * 8 + "CGCGC"
    compound = _make_candidate(22, 30, 2, 100)
    classification = _find_str_segments([_make_candidate(5, 21, 2, 100), compound], contig, _PARAMS)
    assert classification.str_segments == (_StrSegment(start=5, end=21, motif="AT", purity=1.0),)
    assert classification.bundle_members == ()
    assert classification.rejected == ((compound, _Rejection.COMPOUND),)


def test_a_candidate_shorter_than_its_period_is_turned_down_before_the_compound_test() -> None:
    # [5, 7) given as period 4 is shorter than its period, and the four bases from its start,
    # AAAA, are a compound motif. pvc's classify tests the length first.
    rejections = _find_rejections([_make_candidate(5, 7, 4, 100)], "CGCGCAAAAGCG", _PARAMS)
    assert rejections == [_Rejection.COPY_FLOOR]


def test_a_candidate_exactly_one_period_long_still_bundles() -> None:
    # CG at 22, one copy of period 2, is not shorter than its period.
    one_copy = _make_candidate(22, 24, 2, 1)
    classification = _find_str_segments(
        [_make_candidate(5, 21, 2, 100), one_copy], _CLEAN_AT, _PARAMS
    )
    assert classification.str_segments == ()
    assert classification.bundle_members == (_make_candidate(5, 21, 2, 100), one_copy)


def test_step_six_tests_the_copy_floor_then_the_purity_then_the_contig_ends() -> None:
    # The cut ATATGGATAT has 5 copies, under a floor of 6, and a purity of 0.8, under 0.9.
    contig = "CGCGC" + "ATATGGATAT" + "GCGCG"
    parameters = replace(_PARAMS, copy_floors=(10, 6, 4, 3, 3, 3), min_purity=0.9)
    assert _find_rejections([_make_candidate(5, 15, 2, 100)], contig, parameters) == [
        _Rejection.COPY_FLOOR
    ]

    # The cut ATATGGATATAT has a purity of 10/12, under 0.9, and starts on base 0.
    contig = "ATATGGATATAT" + "GCGCG"
    parameters = replace(_PARAMS, min_purity=0.9)
    assert _find_rejections([_make_candidate(0, 12, 2, 100)], contig, parameters) == [
        _Rejection.PURITY
    ]


# ---- The cut's bounds, pvc_eval's own ----------------------------------------------------------
# The expected values are pvc's, from the transcription of its code named above.


def test_the_cut_reaches_no_further_than_half_the_repeat() -> None:
    # 11 bases: the start is searched up to offset 5, so the GCGC at offset 6 is not found.
    assert _find_cut("GCAAAAGCGCG", "GC") is None


def test_the_contig_end_flank_is_measured_on_the_cuts_end() -> None:
    # The detected span ends 12 bases from the contig's end, its cut 15.
    bases = "G" * 40 + "AT" * 8 + "GAT" + "G" * 12
    candidate = _make_candidate(40, 59, 2, 42)
    assert list(_drop_near_contig_ends([candidate], bases, StrFinderParameters())) == [candidate]


def test_a_bundle_member_starts_where_its_cut_starts() -> None:
    member = _make_bundle_member(_make_candidate(0, 11, 2, 5), "ATCATATATAT")
    assert member == _BundleMember(start=3, end=11, motif="AT")


# A contig-end flank of 0 turns step 2 off: a repeat on the first base, or on the last, is
# kept. pvc's values, from the transcription of its row_for_interval with a flank of 0.
_NO_FLANK_CASES = [
    ("on_the_first_base", "AT" * 8 + "G" * 20, _make_candidate(0, 16, 2, 1)),
    ("on_the_last_base", "G" * 20 + "AT" * 8, _make_candidate(20, 36, 2, 1)),
]


@pytest.mark.parametrize(
    ("bases", "candidate"),
    [case[1:] for case in _NO_FLANK_CASES],
    ids=[case[0] for case in _NO_FLANK_CASES],
)
def test_a_contig_end_flank_of_zero_drops_nothing(bases: str, candidate: _Candidate) -> None:
    parameters = StrFinderParameters(contig_end_flank=0)
    assert list(_drop_near_contig_ends([candidate], bases, parameters)) == [candidate]


# pvc's prefilter keeps only a candidate whose end is after its start, and its classify
# passes over one that is not, without a reason. pvc's values, from the transcription of
# its prefilter and classify.
def test_an_empty_candidate_is_dropped_even_at_a_copy_floor_of_zero() -> None:
    no_floors = replace(_PARAMS, copy_floors=_ZERO_FLOORS)
    assert (
        list(
            _drop_under_copy_floor(
                [_make_candidate(5, 5, 2, 1), _make_candidate(9, 7, 2, 1)], no_floors
            )
        )
        == []
    )


def test_an_empty_candidate_is_passed_over_by_step_four_without_a_reason() -> None:
    classification = _find_str_segments([_make_candidate(10, 10, 2, 1)], _CLEAN_AT, _PARAMS)
    assert classification.str_segments == ()
    assert classification.bundle_members == ()
    assert classification.rejected == ()


# ---- Steps 7 and 8, ported from pvc's region_typing/mod.rs at bbcc9f4 ------------------------
#
# pvc's partition tests run its scanning route: a contig-end flank of 0 and the periods
# scanned equal to the periods admitted. pvc's regions are 1-based and inclusive, so its
# [s, e] is [s - 1, e) here; its RepeatInterval is already 0-based and half-open.

# pvc's catalog_config() (mod.rs 1824).
_CATALOG_CONFIG = StrFinderParameters(
    scan_min_period=2,
    scan_max_period=6,
    min_period=2,
    max_period=6,
    match_reward=2,
    mismatch_penalty=7,
    scan_min_copies=2,
    copy_floors=(10, 5, 4, 3, 3, 3),
    min_purity=0.8,
    min_score=0,
    bundle_distance=50,
    contig_end_flank=0,
    max_str_len=1000,
    max_bundle_len=sys.maxsize,  # pvc's u64::MAX: no bundle cap
)
# pvc's TypedRegionConfig::default() at bbcc9f4, which pvc_eval's defaults are not.
_PVC_DEFAULT_CONFIG = StrFinderParameters(
    scan_min_period=1,
    scan_max_period=6,
    min_period=1,
    max_period=6,
    match_reward=2,
    mismatch_penalty=7,
    scan_min_copies=2,
    copy_floors=(8, 6, 6, 6, 5, 4),
    min_purity=0.8,
    min_score=0,
    bundle_distance=15,
    contig_end_flank=0,
    max_str_len=100,
    max_bundle_len=100,
)


def _filler(length: int) -> str:
    """pvc's filler (mod.rs 2203): aperiodic bases, no repeat of any period from 1 to 6."""
    return ("ACGTTGCAAGCTTGCA" * (length // 16 + 1))[:length]


def _put(bases: str, start: int, insert: str) -> str:
    """bases with insert written over it from start, as pvc's copy_from_slice does."""
    return bases[:start] + insert + bases[start + len(insert) :]


def _check_partition(segments: Sequence[_Segment], length: int) -> None:
    """pvc's assert_partitions (mod.rs 1982): the segments run from 0 to the contig's
    length, each starting where the one before ends, with no two of one type in a row."""
    assert segments, "a contig of one base or more has segments"
    expected_start = 0
    previous_type = None
    for segment in segments:
        assert segment.start == expected_start, segments
        assert segment.end > segment.start, segment
        assert segment.region_type is not previous_type, segments
        previous_type = segment.region_type
        expected_start = segment.end
    assert expected_start == length


def _partition(bases: str, parameters: StrFinderParameters) -> list[_Segment]:
    """pvc's partition_resident followed by assert_partitions."""
    segments = list(_find_segments(bases, parameters))
    _check_partition(segments, len(bases))
    return segments


def _get_types(segments: Iterable[_Segment]) -> list[str]:
    """pvc's kinds_of: the types of the segments, in order."""
    return [segment.region_type.value for segment in segments]


def _get_of_type[S: (_Segment, GenomeSegment)](
    segments: Iterable[S], region_type: RegionType
) -> list[S]:
    return [segment for segment in segments if segment.region_type is region_type]


def test_absorption_by_sweep_joins_satellites_that_abut_after_growing_at_radius_zero() -> None:
    # pvc's satellites 22..=86 and 88..=88 and member [86, 88), at a bundle distance of 0,
    # which StrFinderParameters refuses and the absorption step takes as a plain number.
    # The member overlaps only the second satellite; grown, it abuts the first, and the two
    # merge. pvc's half that runs the sweep is dropped: the port has only the fixed point.
    satellites, members, str_segments = _make_grown_satellites(
        [(21, 86), (87, 88)], [_make_candidate(86, 88, 2, 10)], [], 0
    )
    assert satellites == [(21, 88)]
    assert members == []
    assert str_segments == []


def test_a_lone_tract_partitions_as_generic_locus_generic() -> None:
    bases = (
        "ACGTTGCAAGCTCCTAGGATCGATTGCACGGTACCTGAAGCTTGCACTGATCCGTAGGCA"
        + "ATATATATATATATAT"
        + "TGCATTGGACCTAAGCGTTCAGGCTTACGATCCAGGTTACGATCCAAGTGCTTAGCATCG"
    )
    segments = _partition(bases, _CATALOG_CONFIG)
    assert _get_types(segments) == ["generic", "str", "generic"]
    assert segments[1].start == 60
    assert segments[1].motif == "AT"


def test_a_long_array_is_one_satellite_and_swallows_the_locus_inside_it() -> None:
    bases = _filler(60) + "AT" * 1000 + _filler(60)
    segments = _partition(bases, _PVC_DEFAULT_CONFIG)
    satellites = _get_of_type(segments, RegionType.SATELLITE)
    assert len(satellites) == 1
    assert satellites[0].end - satellites[0].start >= 2000
    assert not _get_of_type(segments, RegionType.STR)

    # The control: above the cap, the same repeat is an STR segment.
    segments = _partition(bases, replace(_PVC_DEFAULT_CONFIG, max_str_len=10_000))
    assert len(_get_of_type(segments, RegionType.STR)) == 1
    assert not _get_of_type(segments, RegionType.SATELLITE)


def test_abutting_coverage_runs_merge_into_one() -> None:
    runs = _make_repeat_runs(
        [_make_candidate(0, 10, 2, 1), _make_candidate(10, 20, 3, 1), _make_candidate(29, 40, 2, 1)]
    )
    assert runs == [(0, 20), (29, 40)]
    assert runs[0][1] - runs[0][0] == 20

    # Two 600-base candidates that touch are one run of 1200, over the satellite cap.
    runs = _make_repeat_runs([_make_candidate(0, 600, 2, 1), _make_candidate(600, 1200, 2, 1)])
    assert len(runs) == 1
    assert runs[0][1] - runs[0][0] == 1200


def test_a_run_holds_a_span_inside_another_without_ending_early() -> None:
    assert _make_runs([(0, 50), (10, 20), (20, 30)]) == [(0, 50)]


# A (GA) repeat with an (AAT) repeat after it whose candidates lie inside one another, at
# bundle distance 10 and max_str_len 60. `pvc dev type-regions --min-copies 5,5,4,4,4,3
# --flank-bp 10 --max-str-len 60`, pvc at bbcc9f4, makes it one satellite; a run cut short
# at the inner candidate's end would leave a bundle.
def test_a_candidate_inside_another_does_not_shorten_its_repeat_run(tmp_path: Path) -> None:
    bases = (
        _filler(40) + "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAATAATAATAATAATAAAAATAATAATAAT" + _filler(40)
    )
    parameters = StrFinderParameters(bundle_distance=10, max_str_len=60)
    assert _find_segments_in_fasta("nested", bases, parameters, tmp_path) == [
        (0, 40, "generic", None),
        (40, 103, "satellite", None),
        (103, 142, "generic", None),
    ]


def test_a_repeat_free_contig_is_one_generic_region() -> None:
    bases = "ACGTTGCAAGCTTGCAACGTTGCAAGCTTGCAACGTTGCAAGCTTGCA" * 3
    segments = _partition(bases, _PVC_DEFAULT_CONFIG)
    assert segments == [_Segment(start=0, end=len(bases), region_type=RegionType.GENERIC)]


def test_an_impure_tract_is_generic_when_its_pieces_fall_under_the_copy_floor() -> None:
    bases = _put(_filler(240), 100, "ATATATATGGGGGGGGGGATATATAT")
    segments = _partition(bases, _CATALOG_CONFIG)
    assert not _get_of_type(segments, RegionType.STR)
    assert len(segments) == 1

    # The control: without the interruption the same span is an STR segment.
    bases = _put(_filler(240), 100, "ATATATATATATATATATATATATAT")
    assert _get_types(_partition(bases, _CATALOG_CONFIG)) == ["generic", "str", "generic"]


def test_a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it() -> None:
    # pvc's last assertion, on the counts of rejections, is dropped: the port keeps none.
    bases = _put(_put(_filler(300), 100, "A" * 20), 180, "CAG" * 8)
    segments = _partition(bases, _CATALOG_CONFIG)
    str_segments = _get_of_type(segments, RegionType.STR)
    assert len(str_segments) == 1
    assert str_segments[0].motif == "CAG"
    assert 180 <= str_segments[0].start <= 184
    assert any(
        segment.start <= 109 < segment.end and segment.region_type is RegionType.GENERIC
        for segment in segments
    )


def test_a_homopolymer_of_six_or_more_is_a_period_one_locus_at_default() -> None:
    bases = _put(_filler(300), 100, "A" * 20)
    str_segments = _get_of_type(_partition(bases, _PVC_DEFAULT_CONFIG), RegionType.STR)
    assert len(str_segments) == 1
    assert str_segments[0].motif == "A"
    assert 99 <= str_segments[0].start <= 104


def _make_contig_with_one_guarded_run(run: str, guard: str) -> str:
    """pvc's contig_with_one_guarded_run (mod.rs 2401): the run at 200 in filler, with the
    base on each side set to guard."""
    offset = 200
    bases = _filler(offset * 2 + len(run))
    return bases[: offset - 1] + guard + run + guard + bases[offset + len(run) + 1 :]


def _find_str_segments_at_pvc_default(bases: str) -> list[_Segment]:
    """pvc's loci_at_default (mod.rs 2415)."""
    return _get_of_type(_partition(bases, _PVC_DEFAULT_CONFIG), RegionType.STR)


def test_the_mononucleotide_copy_floor_is_exactly_eight() -> None:
    assert _find_str_segments_at_pvc_default(_make_contig_with_one_guarded_run("A" * 7, "C")) == []
    found = _find_str_segments_at_pvc_default(_make_contig_with_one_guarded_run("A" * 8, "C"))
    assert len(found) == 1
    assert found[0].motif == "A"
    assert found[0].end - found[0].start == 8


def test_the_dinucleotide_copy_floor_is_exactly_six() -> None:
    assert _find_str_segments_at_pvc_default(_make_contig_with_one_guarded_run("AT" * 5, "C")) == []
    found = _find_str_segments_at_pvc_default(_make_contig_with_one_guarded_run("AT" * 6, "C"))
    assert len(found) == 1
    assert found[0].motif == "AT"
    assert found[0].end - found[0].start == 12


def test_a_rejected_repeat_leaves_no_hole() -> None:
    # pvc's two assertions hold of the one satellite the port gives here, as they hold of
    # pvc's own partition: neither checks the type. The candidate of period 1 runs over the
    # whole contig, the ATATAT paid for by the C's and G's, has no cut, and is 126 bases
    # long, over max_str_len.
    bases = "C" * 60 + "ATATAT" + "G" * 60
    segments = _partition(bases, _PVC_DEFAULT_CONFIG)
    assert not _get_of_type(segments, RegionType.STR)
    assert len(segments) == 1


def test_a_tract_at_position_one_is_generic_and_the_partition_still_starts_at_one() -> None:
    segments = _partition("AT" * 8 + "C" * 80, _PVC_DEFAULT_CONFIG)
    assert segments[0].start == 0
    assert not _get_of_type(segments, RegionType.STR)


def test_an_empty_contig_has_no_regions() -> None:
    assert list(_find_segments("", _PVC_DEFAULT_CONFIG)) == []


def _make_contig_with_at_repeats(offsets: Iterable[int], length: int) -> str:
    """pvc's contig_with_tracts_at (mod.rs 2540): an (AT)x10 repeat at each offset, in
    filler."""
    bases = _filler(length)
    for offset in offsets:
        bases = _put(bases, offset, "AT" * 10)
    return bases


def test_two_close_tracts_become_one_bundle_carrying_both() -> None:
    segments = _partition(_make_contig_with_at_repeats([60, 90], 240), _PVC_DEFAULT_CONFIG)
    bundles = _get_of_type(segments, RegionType.STR_BUNDLE)
    assert len(bundles) == 1
    bundle = bundles[0]
    assert len(bundle.members) == 2
    assert bundle.start <= 60
    assert bundle.end >= 110
    assert all(
        bundle.start <= member.start and member.end <= bundle.end for member in bundle.members
    )
    assert all(member.motif in ("AT", "TA") for member in bundle.members)
    assert not _get_of_type(segments, RegionType.STR)


def _partition_two_close_repeats(cap: int) -> list[_Segment]:
    """pvc's two_close_tracts_capped_at (mod.rs 2603), without the counts."""
    bases = _make_contig_with_at_repeats([60, 90], 240)
    return _partition(bases, replace(_PVC_DEFAULT_CONFIG, max_bundle_len=cap))


def test_a_bundle_longer_than_the_cap_is_a_satellite_and_is_counted_as_one() -> None:
    # pvc's four assertions on the counts are dropped: the port keeps none.
    (hull,) = _get_of_type(_partition_two_close_repeats(sys.maxsize), RegionType.STR_BUNDLE)
    hull_length = hull.end - hull.start
    for cap, is_satellite in [
        (hull_length - 1, True),
        (hull_length, False),
        (hull_length + 1, False),
        (0, True),
    ]:
        (at_hull,) = [
            segment
            for segment in _partition_two_close_repeats(cap)
            if (segment.start, segment.end) == (hull.start, hull.end)
        ]
        assert (at_hull.region_type is RegionType.SATELLITE) == is_satellite, cap


def test_a_bundle_at_a_contig_end_is_dropped_unless_capped() -> None:
    def get_types(offsets: list[int], cap: int) -> list[str]:
        bases = _make_contig_with_at_repeats(offsets, 200)
        return _get_types(_partition(bases, replace(_PVC_DEFAULT_CONFIG, max_bundle_len=cap)))

    assert get_types([0, 30], 100) == ["generic"]
    assert get_types([150, 180], 100) == ["generic"]
    assert get_types([1, 31], 100) == ["generic", "str-bundle", "generic"]
    assert get_types([0, 30], 10) == ["satellite", "generic"]

    # One base in from the right end; the C continues neither AT nor TA.
    bases = _make_contig_with_at_repeats([149, 179], 200)[:199] + "C"
    (bundle,) = _get_of_type(_partition(bases, _PVC_DEFAULT_CONFIG), RegionType.STR_BUNDLE)
    assert bundle.end == 199


def test_the_hull_ends_past_a_ragged_last_member() -> None:
    bases = _filler(200) + "TCTC" * 8 + _filler(9) + "CAG" * 8 + "CA" + _filler(200)
    (bundle,) = _get_of_type(_partition(bases, _PVC_DEFAULT_CONFIG), RegionType.STR_BUNDLE)
    # The CAG run is [241, 265) and the CA [265, 267).
    assert max(member.end for member in bundle.members) == 265
    assert bundle.end == 267


def test_three_chained_tracts_become_one_bundle_of_three() -> None:
    segments = _partition(_make_contig_with_at_repeats([60, 110, 160], 300), _CATALOG_CONFIG)
    (bundle,) = _get_of_type(segments, RegionType.STR_BUNDLE)
    assert len(bundle.members) == 3
    assert bundle.start <= 60
    assert bundle.end >= 180


def test_a_tract_clear_of_a_bundle_is_still_admitted() -> None:
    segments = _partition(_make_contig_with_at_repeats([60, 90, 400], 600), _PVC_DEFAULT_CONFIG)
    assert len(_get_of_type(segments, RegionType.STR_BUNDLE)) == 1
    str_segments = _get_of_type(segments, RegionType.STR)
    assert len(str_segments) == 1
    assert str_segments[0].start == 400


def test_a_low_copy_repeat_is_generic_not_a_bundle() -> None:
    segments = _partition(_put(_filler(240), 100, "ATATAT"), _PVC_DEFAULT_CONFIG)
    assert not _get_of_type(segments, RegionType.STR_BUNDLE)
    assert len(segments) == 1


def test_bundle_clusters_regroups_exactly_what_the_split_set_aside() -> None:
    bundled = [
        _make_candidate(100, 130, 2, 100),
        _make_candidate(150, 180, 2, 100),
        _make_candidate(5000, 5030, 2, 100),
        _make_candidate(5040, 5070, 2, 100),
    ]
    clusters = _make_clusters(bundled, 50)
    assert [len(cluster) for cluster in clusters] == [2, 2]
    assert [cluster[0].start for cluster in clusters] == [100, 5000]


_GOLDEN_CATALOG = Path(__file__).parent / "data" / "tandem_repeat" / "golden.ssr_catalog.bed.gz"


def test_the_resident_partition_reproduces_the_golden_catalog() -> None:
    # The settings the golden catalogue was built at: its header gives the purity, the
    # score and the bundle distance; pvc's test pins the rest.
    parameters = replace(
        _CATALOG_CONFIG, max_bundle_len=100, min_purity=0.8, min_score=0, bundle_distance=50
    )
    with gzip.open(_GOLDEN_CATALOG, "rt", encoding="utf-8") as catalog:
        golden = [
            (columns[0], int(columns[1]), int(columns[2]))
            for columns in (line.split("\t") for line in catalog if not line.startswith("#"))
        ]
    assert golden
    found: list[tuple[str, int, int]] = []
    satellites: list[tuple[str, int, int]] = []
    with IndexedFasta(_SYNTHETIC_REF) as reference:
        for contig in reference.contigs:
            segments = list(segment_contig(reference, contig, parameters))
            _check_genome_segments(segments, contig, reference.get_length(contig))
            for segment in segments:
                where = (contig, segment.start, segment.end)
                if segment.region_type is RegionType.STR:
                    found.append(where)
                elif segment.region_type is RegionType.SATELLITE:
                    satellites.append(where)
    assert found

    def overlaps(one: tuple[str, int, int], other: tuple[str, int, int]) -> bool:
        return one[0] == other[0] and one[1] < other[2] and other[1] < one[2]

    missed = [
        golden_str
        for golden_str in golden
        if not any(overlaps(golden_str, other) for other in found)
        and not any(overlaps(golden_str, satellite) for satellite in satellites)
    ]
    assert missed == []


def _make_windowing_fixture() -> str:
    """pvc's windowing_fixture (mod.rs 3071): four (AT)x10 repeats and a 1200-base (AT)
    array at 4000, in 6000 bases of filler."""
    bases = _make_contig_with_at_repeats([990, 1990, 2020, 3500], 6000)
    return _put(bases, 4000, "AT" * 600)


def _make_chained_cluster(span: int, length: int) -> str:
    """pvc's chained_cluster (mod.rs 3089): (AT)x10 repeats every 50 bases from 1000."""
    return _make_contig_with_at_repeats(range(1000, 1000 + span - 20 + 1, 50), length)


def _make_genome_segments(segments: Iterable[_Segment]) -> list[GenomeSegment]:
    """The segments as segment_contig yields them, on a contig named chr1."""
    return [
        GenomeSegment("chr1", segment.start, segment.end, segment.region_type, segment.motif)
        for segment in segments
    ]


def _find_requested(bases: str, requested: list[tuple[int, int]]) -> list[GenomeSegment]:
    """pvc's walk_bed: the STR regions of the whole contig that overlap the requested
    spans, which are pvc's, 1-based and inclusive."""
    spans = [(start - 1, end) for start, end in requested]
    segments = _make_genome_segments(_find_segments(bases, _PVC_DEFAULT_CONFIG))
    return list(_find_requested_str_regions(segments, spans))


def test_a_bed_returns_the_same_findings_the_whole_genome_run_does() -> None:
    # The port segments the whole contig whatever spans are asked for, and a span only
    # chooses which STR regions come back (`str-regions --regions`); generic regions do
    # not come back at all.
    bases = _make_contig_with_at_repeats([], 8000)
    bases = _put(bases, 1000, "AT" * 1500)
    whole = _partition(bases, _CATALOG_CONFIG)
    (truth,) = _get_of_type(whole, RegionType.SATELLITE)
    assert truth.end - truth.start > 3000
    got = list(_find_requested_str_regions(_make_genome_segments(whole), [(1099, 1300)]))
    assert _get_of_type(got, RegionType.SATELLITE) == _make_genome_segments([truth])

    bases = _make_chained_cluster(1500, 6000)
    whole = _partition(bases, _CATALOG_CONFIG)
    (truth,) = _get_of_type(whole, RegionType.STR_BUNDLE)
    assert truth.end - truth.start > 1400
    assert len(truth.members) >= 25
    assert not _get_of_type(whole, RegionType.SATELLITE)
    got = list(_find_requested_str_regions(_make_genome_segments(whole), [(899, 1200)]))
    assert _get_of_type(got, RegionType.STR_BUNDLE) == _make_genome_segments([truth])


def _make_micro_near_satellite(micro_left: bool, gap: int) -> str:
    """pvc's micro_near_satellite (mod.rs 3192): a (CAG)x8 repeat gap bases from a
    1300-base (AT) array, on its left or its right."""
    if micro_left:
        micro_at, array_at = 1000, 1000 + 24 + gap
    else:
        micro_at, array_at = 1000 + 1300 + gap, 1000
    bases = _put(_filler(6000), micro_at, "CAG" * 8)
    return _put(bases, array_at, "AT" * 650)


@pytest.mark.parametrize("micro_left", [True, False], ids=["left", "right"])
def test_a_microsatellite_beside_a_satellite_is_absorbed_into_it(micro_left: bool) -> None:
    near = 15 // 2  # pvc's DEFAULT_BUNDLE_THRESHOLD / 2
    far = 15 * 10
    segments = _partition(_make_micro_near_satellite(micro_left, near), _PVC_DEFAULT_CONFIG)
    assert _get_types(segments) == ["generic", "satellite", "generic"]
    satellite = segments[1]
    assert satellite.end - satellite.start >= 1300 + near + 24

    # The control: far from the array, the same repeat is an STR segment.
    types = _get_types(_partition(_make_micro_near_satellite(micro_left, far), _PVC_DEFAULT_CONFIG))
    assert "str" in types
    assert "satellite" in types


def test_absorption_tracks_the_bundle_radius() -> None:
    def find_first_str_gap(radius: int) -> int | None:
        parameters = replace(_PVC_DEFAULT_CONFIG, bundle_distance=radius)
        for gap in range(1, radius * 3 + 1):
            segments = _find_segments(_make_micro_near_satellite(True, gap), parameters)
            if _get_of_type(segments, RegionType.STR):
                return gap
        return None

    previous = 0
    for radius in [10, 15, 20, 30]:
        gap = find_first_str_gap(radius)
        assert gap is not None, radius
        assert gap > previous, radius
        previous = gap


def test_the_satellite_cap_applies_to_the_cleaned_coverage_not_the_raw() -> None:
    bases = _put(_make_windowing_fixture(), 5200, "CGCG")
    scanned = replace(_PVC_DEFAULT_CONFIG, copy_floors=_ZERO_FLOORS)
    raw = _find_candidates(bases, scanned)
    cleaned = list(_drop_under_copy_floor(raw, _PVC_DEFAULT_CONFIG))

    def find_over_cap(candidates: list[_Candidate]) -> list[tuple[int, int]]:
        return [
            (start, end)
            for start, end in _make_repeat_runs(candidates)
            if end - start > _PVC_DEFAULT_CONFIG.max_str_len
        ]

    from_cleaned = find_over_cap(cleaned)
    assert find_over_cap(raw) != from_cleaned
    satellites = _get_of_type(_partition(bases, _PVC_DEFAULT_CONFIG), RegionType.SATELLITE)
    assert [(satellite.start, satellite.end) for satellite in satellites] == from_cleaned


# The four tests below ask for parts of a contig. pvc's emit_into keeps the findings that
# overlap a requested span whole and cuts generic regions at its edges; the port keeps the
# first rule, for `str-regions --regions`, and has no generic regions to cut. Each test
# names the assertions it drops.


@pytest.mark.parametrize(
    "requested",
    [(900, 1100), (1995, 2100), (4500, 4800), (3000, 3600), (1, 6000)],
    ids=["str", "bundle", "satellite", "interior", "whole"],
)
def test_a_bed_does_not_change_what_things_are(requested: tuple[int, int]) -> None:
    # pvc's assertion on the type at each edge of what comes back holds of the STR regions
    # only, which come back as the whole contig's.
    bases = _make_windowing_fixture()
    whole = _partition(bases, _PVC_DEFAULT_CONFIG)
    got = _find_requested(bases, [requested])
    first, last = requested[0] - 1, requested[1]
    for segment in got:
        for position in (segment.start, segment.end - 1):
            if first <= position < last:
                (truth,) = [other for other in whole if other.start <= position < other.end]
                assert truth.region_type is segment.region_type

    def find_starting_inside(
        segments: Sequence[_Segment] | Sequence[GenomeSegment],
    ) -> list[tuple[int, int]]:
        return [
            (segment.start, segment.end)
            for segment in segments
            if segment.region_type is RegionType.STR and first <= segment.start < last
        ]

    assert find_starting_inside(got) == find_starting_inside(whole)


@pytest.mark.parametrize(
    "requested", [(900, 1100), (1995, 2100), (3000, 3600)], ids=["str", "bundle", "interior"]
)
def test_a_bed_span_is_tiled_by_what_comes_back(requested: tuple[int, int]) -> None:
    # Dropped: the first segment starting at or before the span and the last ending at or
    # after it, and consecutive segments abutting with different types. They hold of the
    # generic regions, which do not come back. What is left: something comes back, in
    # order and without overlaps.
    got = _find_requested(_make_windowing_fixture(), [requested])
    assert got
    assert all(one.end <= other.start for one, other in itertools.pairwise(got))


def test_generic_clips_at_the_edge_and_findings_come_back_whole() -> None:
    # Dropped: the first and last generic regions of (3000, 3600) cut at the span's edges.
    bases = _make_windowing_fixture()
    got = _find_requested(bases, [(4500, 4800)])
    assert len(got) == 1
    satellite = got[0]
    assert satellite.region_type is RegionType.SATELLITE
    assert satellite.start < 4499
    assert satellite.end > 4800
    assert satellite.end - satellite.start > 1000
    whole = _partition(bases, _PVC_DEFAULT_CONFIG)
    assert _make_genome_segments(_get_of_type(whole, RegionType.SATELLITE)) == [satellite]

    (bundle,) = _get_of_type(_find_requested(bases, [(1995, 2100)]), RegionType.STR_BUNDLE)
    assert bundle.start < 1994
    # The members are the whole contig's bundle's: GenomeSegment does not carry them.
    (members,) = [other.members for other in whole if other.start == bundle.start]
    assert all(bundle.start <= member.start and member.end <= bundle.end for member in members)


def test_a_locus_straddling_the_bed_edge_is_emitted_whole() -> None:
    bases = _make_windowing_fixture()
    (str_region,) = _get_of_type(_find_requested(bases, [(900, 1000)]), RegionType.STR)
    assert str_region.end > 1000
    whole = _partition(bases, _PVC_DEFAULT_CONFIG)
    assert [other for other in _make_genome_segments(whole) if other.start == str_region.start] == [
        str_region
    ]


def test_two_spans_sharing_a_scan_span_do_not_leak_the_gap_between_them() -> None:
    # Dropped: the 402 bases the generic regions cover, which do not come back.
    got = _find_requested(_make_windowing_fixture(), [(3000, 3200), (3400, 3600)])
    assert got
    assert all(2999 <= s.start < 3200 or 3399 <= s.start < 3600 for s in got)


# ---- Steps 7 and 8: pvc_eval's own cases -----------------------------------------------------


@pytest.mark.parametrize(
    ("member", "is_absorbed"),
    [
        (_make_candidate(71, 86, 2, 1), True),
        (_make_candidate(70, 85, 2, 1), False),
        (_make_candidate(214, 230, 2, 1), True),
        (_make_candidate(215, 230, 2, 1), False),
    ],
    ids=["left-14", "left-15", "right-14", "right-15"],
)
def test_a_feature_is_absorbed_with_fewer_than_the_bundle_distance_bases_between(
    member: _Candidate, is_absorbed: bool
) -> None:
    # The satellite is [100, 200); the bundle distance 15. doc/specs/str_finder.md step 7.
    satellites, members, _ = _make_grown_satellites([(100, 200)], [member], [], 15)
    if is_absorbed:
        assert (satellites, members) == ([(min(member.start, 100), max(member.end, 200))], [])
    else:
        assert (satellites, members) == ([(100, 200)], [member])


def test_absorption_is_repeated_until_a_pass_absorbs_nothing() -> None:
    # Taken in order, the member at 60 is 20 bases from the satellite at 100 and is not
    # absorbed; the one at 82 is. The grown satellite then reaches the member at 60, which
    # only a second pass absorbs. Likewise the STR segment at 20 is reached only once the
    # member at 60 is absorbed.
    far = _make_candidate(60, 80, 2, 1)
    near = _make_candidate(82, 90, 2, 1)
    str_segment = _StrSegment(start=20, end=50, motif="CA", purity=1.0)
    satellites, members, str_segments = _make_grown_satellites(
        [(100, 200), (400, 500)], [far, near], [str_segment], 15
    )
    assert satellites == [(20, 200), (400, 500)]
    assert members == []
    assert str_segments == []


def test_absorption_merges_satellites_that_a_feature_reaches_on_both_sides() -> None:
    satellites, members, _ = _make_grown_satellites(
        [(0, 120), (140, 300)], [_make_candidate(125, 135, 2, 1)], [], 15
    )
    assert satellites == [(0, 300)]
    assert members == []


_FOUR_SATELLITES = [(0, 50), (100, 200), (300, 400), (500, 600)]
_TWO_CLOSE_SATELLITES = [(0, 50), (100, 200), (230, 330), (500, 600)]


@pytest.mark.parametrize(
    ("satellites", "member", "expected"),
    [
        (_FOUR_SATELLITES, (276, 286), [(0, 50), (100, 200), (276, 400), (500, 600)]),
        (_FOUR_SATELLITES, (275, 285), _FOUR_SATELLITES),
        (_FOUR_SATELLITES, (414, 424), [(0, 50), (100, 200), (300, 424), (500, 600)]),
        (_FOUR_SATELLITES, (415, 425), _FOUR_SATELLITES),
        (_FOUR_SATELLITES[1:], (76, 86), [(76, 200), (300, 400), (500, 600)]),
        (_FOUR_SATELLITES[1:], (75, 85), _FOUR_SATELLITES[1:]),
        (_FOUR_SATELLITES, (614, 620), [(0, 50), (100, 200), (300, 400), (500, 620)]),
        (_FOUR_SATELLITES, (615, 620), _FOUR_SATELLITES),
        (_TWO_CLOSE_SATELLITES, (214, 216), [(0, 50), (100, 330), (500, 600)]),
        (_TWO_CLOSE_SATELLITES, (215, 216), [(0, 50), (100, 200), (215, 330), (500, 600)]),
        (_TWO_CLOSE_SATELLITES, (214, 215), [(0, 50), (100, 215), (230, 330), (500, 600)]),
    ],
    ids=[
        "before-inner-14",
        "before-inner-15",
        "after-inner-14",
        "after-inner-15",
        "before-first-14",
        "before-first-15",
        "after-last-14",
        "after-last-15",
        "between-14-and-14-reaches-both",
        "between-15-and-14-reaches-right",
        "between-14-and-15-reaches-left",
    ],
)
def test_among_several_satellites_a_feature_grows_only_those_at_fewer_than_the_bundle_distance(
    satellites: list[tuple[int, int]], member: tuple[int, int], expected: list[tuple[int, int]]
) -> None:
    # The bundle distance is 15: a gap of 14 bases between the member and a satellite is
    # absorbed, one of 15 is not. doc/specs/str_finder.md step 7, worked by hand.
    candidate = _make_candidate(*member, 2, 1)
    grown, members, _ = _make_grown_satellites(satellites, [candidate], [], 15)
    assert grown == expected
    assert members == ([] if expected != satellites else [candidate])


@st.composite
def _make_satellites_and_span(
    draw: st.DrawFn,
) -> tuple[list[tuple[int, int]], tuple[int, int], int]:
    """Satellites sorted, disjoint and not abutting, a span, and a bundle distance."""
    satellites = []
    position = draw(st.integers(min_value=0, max_value=30))
    for _ in range(draw(st.integers(min_value=0, max_value=8))):
        length = draw(st.integers(min_value=1, max_value=30))
        satellites.append((position, position + length))
        position += length + draw(st.integers(min_value=1, max_value=30))
    start = draw(st.integers(min_value=0, max_value=position + 30))
    span = (start, start + draw(st.integers(min_value=1, max_value=30)))
    return satellites, span, draw(st.integers(min_value=0, max_value=20))


@settings(max_examples=2000, deadline=None)
@given(case=_make_satellites_and_span())
def test_the_satellites_a_span_reaches_are_those_a_scan_of_every_satellite_finds(
    case: tuple[list[tuple[int, int]], tuple[int, int], int],
) -> None:
    # The oracle is step 7's rule written as a scan of every satellite, then the merge of
    # those that overlap or abut, which a bundle distance of 0 can leave.
    satellites, (start, end), distance = case
    reached = [
        (satellite_start, satellite_end)
        for satellite_start, satellite_end in satellites
        if start - satellite_end < distance and satellite_start - end < distance
    ]
    got = _make_satellites_absorbing(satellites, (start, end), distance)
    if not reached:
        assert got is None
        return
    grown = (min(start, reached[0][0]), max(end, reached[-1][1]))
    expected: list[tuple[int, int]] = []
    for satellite_start, satellite_end in sorted(
        [*(satellite for satellite in satellites if satellite not in reached), grown]
    ):
        if expected and satellite_start <= expected[-1][1]:
            expected[-1] = (expected[-1][0], max(expected[-1][1], satellite_end))
        else:
            expected.append((satellite_start, satellite_end))
    assert got == expected


def _write_fasta(directory: Path, contigs: dict[str, str]) -> Path:
    """A FASTA of the contigs, one line each, and its .fai, as samtools faidx writes it,
    except that a contig of no bases keeps its line."""
    fasta = ""
    fai = ""
    for name, bases in contigs.items():
        fasta += f">{name}\n"
        fai += (
            f"{name}\t{len(bases)}\t{len(fasta)}\t{max(len(bases), 1)}\t{max(len(bases), 1) + 1}\n"
        )
        fasta += f"{bases}\n"
    path = directory / "ref.fa"
    path.write_text(fasta, encoding="ascii")
    Path(f"{path}.fai").write_text(fai, encoding="utf-8")
    return path


def _find_segments_in_fasta(
    contig: str, bases: str, parameters: StrFinderParameters, tmp_path: Path
) -> list[tuple[int, int, str, str | None]]:
    """The segments segment_contig gives for a contig of these bases, as (start, end, type,
    motif)."""
    with IndexedFasta(_write_fasta(tmp_path, {contig: bases})) as reference:
        segments = list(segment_contig(reference, contig, parameters))
    _check_genome_segments(segments, contig, len(bases))
    return [
        (segment.start, segment.end, segment.region_type.value, segment.motif)
        for segment in segments
        if segment.region_type is not None
    ]


def _check_genome_segments(segments: Sequence[GenomeSegment], contig: str, length: int) -> None:
    """doc/spec.md §4: one contig's segments, consecutive and non-overlapping, from 0 to
    its length, each with a type."""
    assert all(segment.contig == contig for segment in segments)
    assert all(segment.region_type is not None for segment in segments)
    _check_consecutive([(segment.start, segment.end) for segment in segments], length)


def _check_consecutive(spans: Sequence[tuple[int, int]], length: int) -> None:
    """The spans run from 0 to length, each non-empty and starting where the one before
    ends; none for a length of 0."""
    ends = [0] + [end for _, end in spans]
    assert [start for start, _ in spans] == ends[:-1]
    assert all(end > start for start, end in spans)
    assert ends[-1] == length


# The worked cases of doc/specs/str_finder.md, with the segments `pvc repeat-catalog` and then
# `pvc dev type-regions --min-copies 5,5,4,4,4,3` gave for them, pvc built at bbcc9f4.
_FOUR_A = (
    "GCTAGACATACATACATACAGTCAGCACGACTGTCAGTGATCGCTAGTAGTAGTGATGCACACACACACACACACACACAGTCTGT"
    "AAAACTACGCTACTGCTGTCACTCGACTGCATACATCAGACAGACTCGTATGACAGTCACGCAGA"
)
_WORKED_CASES = [
    (
        "four_a",
        _FOUR_A,
        [(0, 58, "generic", None), (58, 80, "str", "CA"), (80, 151, "generic", None)],
    ),
    (
        "five_a",
        _FOUR_A.replace("AAAA", "AAAAA"),
        [(0, 58, "generic", None), (58, 91, "str-bundle", None), (91, 152, "generic", None)],
    ),
    (
        "near_end",
        "GCTAGACATGAAAAAACTGTCAGTCACACACACACACACACACAGTCTGTAGCTACTGCTGTCACTCGACTGCATACATCAGACAG"
        "ACTCGTATGACAGTCACG",
        [(0, 24, "generic", None), (24, 44, "str", "CA"), (44, 104, "generic", None)],
    ),
    (
        "far_end",
        "GCTAGACATGCATCAGTAAAAAACTGTCAGTCACACACACACACACACACAGTCTGTAGCTACTGCTGTCACTCGACTGCATACATC"
        "AGACAGACTCGTATGACAGTCACG",
        [(0, 17, "generic", None), (17, 51, "str-bundle", None), (51, 111, "generic", None)],
    ),
]


@pytest.mark.parametrize(
    ("contig", "bases", "expected"), _WORKED_CASES, ids=[case[0] for case in _WORKED_CASES]
)
def test_the_worked_cases_give_pvcs_segments(
    contig: str, bases: str, expected: list[tuple[int, int, str, str | None]], tmp_path: Path
) -> None:
    assert _find_segments_in_fasta(contig, bases, StrFinderParameters(), tmp_path) == expected


# 120 A's in pvc's filler. `pvc dev type-regions --min-copies 5,5,4,4,4,3 --min-period 2`,
# pvc at bbcc9f4, makes them a satellite: its catalogue scans periods 1 to 6, and the
# candidate of period 1 makes a repeat run although period 1 is not admitted. pvc's scanning
# route, which scans periods 2 to 6 alone, leaves them generic (doc/specs/str_finder.md,
# step 1).
_A120 = _filler(60) + "A" * 120 + "C" + _filler(59)


def test_a_period_scanned_and_not_admitted_still_makes_a_satellite(tmp_path: Path) -> None:
    parameters = StrFinderParameters(min_period=2)
    assert _find_segments_in_fasta("a120", _A120, parameters, tmp_path) == [
        (0, 60, "generic", None),
        (60, 180, "satellite", None),
        (180, 240, "generic", None),
    ]


def test_a_period_not_scanned_makes_no_repeat_run(tmp_path: Path) -> None:
    parameters = StrFinderParameters(scan_min_period=2, min_period=2)
    assert _find_segments_in_fasta("a120", _A120, parameters, tmp_path) == [
        (0, 240, "generic", None)
    ]


def _make_guarded(run: str, guard: str) -> str:
    """The run with a guard base on each side and 29 bases of pvc's filler beyond it, so
    that 30 bases lie between the run and each contig end."""
    return _filler(29) + guard + run + guard + _filler(29)


# pvc_eval's copy floors of periods 1 and 2 at their boundary, 5 copies each. The segments
# are what `pvc dev type-regions --min-copies 5,5,4,4,4,3` gave, pvc at bbcc9f4.
@pytest.mark.parametrize(
    ("bases", "expected"),
    [
        (_make_guarded("A" * 4, "C"), [(0, 64, "generic", None)]),
        (
            _make_guarded("A" * 5, "C"),
            [(0, 30, "generic", None), (30, 35, "str", "A"), (35, 65, "generic", None)],
        ),
        (_make_guarded("CA" * 4, "G"), [(0, 68, "generic", None)]),
        (
            _make_guarded("CA" * 5, "G"),
            [(0, 30, "generic", None), (30, 40, "str", "CA"), (40, 70, "generic", None)],
        ),
    ],
    ids=["A4", "A5", "CA4", "CA5"],
)
def test_the_default_copy_floors_hold_at_their_boundary(
    bases: str, expected: list[tuple[int, int, str, str | None]], tmp_path: Path
) -> None:
    assert _find_segments_in_fasta("floor", bases, StrFinderParameters(), tmp_path) == expected


# A homopolymer at max_str_len, 100 bases, and one past it, and a run over max_str_len
# whose cut has 10 bases between it and a contig end, so that step 2 drops it. The
# segments are what `pvc dev type-regions --min-copies 5,5,4,4,4,3` gave, pvc at bbcc9f4.
@pytest.mark.parametrize(
    ("bases", "expected"),
    [
        (
            _filler(59) + "C" + "A" * 100 + "C" + _filler(59),
            [(0, 60, "generic", None), (60, 160, "str", "A"), (160, 220, "generic", None)],
        ),
        (
            _filler(59) + "C" + "A" * 101 + "C" + _filler(59),
            [(0, 60, "generic", None), (60, 161, "satellite", None), (161, 221, "generic", None)],
        ),
        (_filler(9) + "C" + "A" * 120 + "C" + _filler(59), [(0, 190, "generic", None)]),
        (_filler(59) + "C" + "A" * 120 + "C" + _filler(9), [(0, 190, "generic", None)]),
    ],
    ids=["100", "101", "near-start", "near-end"],
)
def test_a_repeat_run_longer_than_max_str_len_is_a_satellite_unless_step_2_drops_it(
    bases: str, expected: list[tuple[int, int, str, str | None]], tmp_path: Path
) -> None:
    assert _find_segments_in_fasta("run", bases, StrFinderParameters(), tmp_path) == expected


@pytest.mark.parametrize(
    ("requested", "is_kept"),
    [((0, 10), False), ((0, 11), True), ((19, 30), True), ((20, 30), False)],
    ids=["ends-at-start", "one-base-in", "starts-at-last-base", "starts-at-end"],
)
def test_a_requested_span_keeps_an_str_region_it_overlaps_by_one_base(
    requested: tuple[int, int], is_kept: bool
) -> None:
    str_segment = GenomeSegment("chr1", 10, 20, RegionType.STR, "CA")
    segments = [
        GenomeSegment("chr1", 0, 10, RegionType.GENERIC, None),
        str_segment,
        GenomeSegment("chr1", 20, 30, RegionType.GENERIC, None),
    ]
    got = list(_find_requested_str_regions(segments, [requested]))
    assert got == ([str_segment] if is_kept else [])


def test_requested_spans_in_any_order_keep_the_str_regions_they_overlap() -> None:
    str_segment = GenomeSegment("chr1", 10, 20, RegionType.STR, "CA")
    bundle = GenomeSegment("chr1", 22, 30, RegionType.STR_BUNDLE, None)
    segments = [
        GenomeSegment("chr1", 0, 10, RegionType.GENERIC, None),
        str_segment,
        GenomeSegment("chr1", 20, 22, RegionType.GENERIC, None),
        bundle,
        GenomeSegment("chr1", 30, 40, RegionType.GENERIC, None),
    ]
    requested = [(35, 40), (0, 11), (5, 8), (29, 31)]
    assert list(_find_requested_str_regions(segments, requested)) == [str_segment, bundle]


# The hull ends at the furthest end of any member, not at the end of the member listed
# last: the period-1 G's lie inside the (TGG) repeat and are listed after it. The segments
# are what `pvc dev type-regions --min-copies 5,5,4,4,4,3` gave, pvc at bbcc9f4, with the
# members (59, 101, TGG) and (75, 80, G) its BED lists.
_INNER_LAST_MEMBER = _filler(60) + "GGTGGTGGTGGTGGTGGGGGTGGTAGTGGTGGTGGTGGTGGT" + _filler(60)


def test_a_bundles_hull_ends_at_its_furthest_member_end(tmp_path: Path) -> None:
    expected = [
        (0, 59, "generic", None),
        (59, 102, "str-bundle", None),
        (102, 162, "generic", None),
    ]
    assert _find_segments_in_fasta(
        "inner", _INNER_LAST_MEMBER, StrFinderParameters(), tmp_path
    ) == (expected)
    (bundle,) = _get_of_type(
        _find_segments(_INNER_LAST_MEMBER, StrFinderParameters()), RegionType.STR_BUNDLE
    )
    assert bundle.members == (
        _BundleMember(start=59, end=101, motif="TGG"),
        _BundleMember(start=75, end=80, motif="G"),
    )


# Two bundles of two (AT)x10 repeats, the detected spans of the first ending at 110 and of
# the second starting at 126: 16 bases between them. The segments are what `pvc dev
# type-regions --min-copies 5,5,4,4,4,3 --flank-bp 16` (and `--flank-bp 17`) gave, pvc at
# bbcc9f4.
@pytest.mark.parametrize(
    ("bundle_distance", "expected"),
    [
        (
            16,
            [
                (0, 59, "generic", None),
                (59, 110, "str-bundle", None),
                (110, 126, "generic", None),
                (126, 177, "str-bundle", None),
                (177, 270, "generic", None),
            ],
        ),
        (
            17,
            [(0, 59, "generic", None), (59, 177, "satellite", None), (177, 270, "generic", None)],
        ),
    ],
    ids=["16", "17"],
)
def test_bundles_the_bundle_distance_apart_stay_two(
    bundle_distance: int, expected: list[tuple[int, int, str, str | None]], tmp_path: Path
) -> None:
    bases = _make_contig_with_at_repeats([60, 90, 126, 156], 270)
    parameters = StrFinderParameters(bundle_distance=bundle_distance)
    assert _find_segments_in_fasta("pairs", bases, parameters, tmp_path) == expected


def test_an_empty_contig_has_no_segments(tmp_path: Path) -> None:
    with IndexedFasta(_write_fasta(tmp_path, {"empty": "", "chr1": "ACGT"})) as reference:
        assert list(segment_contig(reference, "empty")) == []


@pytest.mark.parametrize("bases", ["A", "ACG", "AAAAA"], ids=["1", "3", "5"])
def test_a_contig_shorter_than_the_longest_period_is_one_generic_region(
    bases: str, tmp_path: Path
) -> None:
    # Five A's are a candidate of period 1 at both ends of the contig, which step 2
    # drops; at a flank of 0 step 6 turns it down for touching the ends.
    for flank in (15, 0):
        parameters = StrFinderParameters(contig_end_flank=flank)
        assert _find_segments_in_fasta("short", bases, parameters, tmp_path) == [
            (0, len(bases), "generic", None)
        ]


def test_a_contig_the_fai_lacks_raises_value_error(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT"})
    with IndexedFasta(fasta) as reference, pytest.raises(ValueError, match="chr2"):
        segment_contig(reference, "chr2")


def test_a_byte_that_is_not_a_letter_raises_pvc_eval_input_error(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "AC-GT"})
    with (
        IndexedFasta(fasta) as reference,
        pytest.raises(PvcEvalInputError, match="chr1:3 is not a letter"),
    ):
        segment_contig(reference, "chr1")


def test_a_closed_reference_raises_value_error(tmp_path: Path) -> None:
    reference = IndexedFasta(_write_fasta(tmp_path, {"chr1": "ACGT"}))
    reference.close()
    with pytest.raises(ValueError, match="is closed"):
        segment_contig(reference, "chr1")


@pytest.mark.parametrize(
    "parameters",
    [StrFinderParameters(), _PVC_DEFAULT_CONFIG, _CATALOG_CONFIG],
    ids=["defaults", "pvc-default", "catalog"],
)
def test_the_fixtures_segments_run_from_0_to_each_contigs_length(
    parameters: StrFinderParameters,
) -> None:
    with IndexedFasta(_SYNTHETIC_REF) as reference:
        for contig in reference.contigs:
            segments = list(segment_contig(reference, contig, parameters))
            _check_genome_segments(segments, contig, reference.get_length(contig))
            assert {segment.region_type for segment in segments} != {RegionType.GENERIC}


_PLANTED_MOTIFS = ["A", "C", "AT", "CA", "CAG", "GATA", "AACGT", "AACGTT"]


@st.composite
def _make_contigs_with_repeats(draw: st.DrawFn) -> str:
    """Random bases with repeats planted in them, some long, some close to each other,
    some at the contig ends, some with a base changed."""
    pieces = []
    for _ in range(draw(st.integers(0, 8))):
        pieces.append(draw(st.text("ACGTN", min_size=0, max_size=40)))
        motif = draw(st.sampled_from(_PLANTED_MOTIFS))
        run = motif * draw(st.integers(2, 70))
        if draw(st.booleans()):
            at = draw(st.integers(0, len(run) - 1))
            run = run[:at] + draw(st.sampled_from("ACGT")) + run[at + 1 :]
        pieces.append(run)
    pieces.append(draw(st.text("ACGT", min_size=0, max_size=40)))
    return "".join(pieces)


@settings(max_examples=200, deadline=None)
@given(bases=_make_contigs_with_repeats(), flank=st.sampled_from([0, 15]))
def test_the_segments_of_any_contig_run_from_0_to_its_length(bases: str, flank: int) -> None:
    parameters = StrFinderParameters(contig_end_flank=flank, bundle_distance=10, max_str_len=60)
    segments = _find_segments(bases, parameters)
    _check_consecutive([(segment.start, segment.end) for segment in segments], len(bases))


def test_four_threads_give_what_one_thread_gives() -> None:
    with IndexedFasta(_SYNTHETIC_REF) as reference:
        contigs = reference.contigs * 4

        def find_segments(contig: str) -> list[GenomeSegment]:
            return list(segment_contig(reference, contig, _CATALOG_CONFIG))

        one_thread = [find_segments(contig) for contig in contigs]
        assert all(one_thread)
        with ThreadPoolExecutor(max_workers=4) as pool:
            four_threads = list(pool.map(find_segments, contigs))
    assert four_threads == one_thread
