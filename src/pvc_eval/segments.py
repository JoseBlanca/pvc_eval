"""The segmented genome, stage 1 of doc/spec.md §4: what the STR finder yields to the loci.

The STR finder is specified in doc/specs/str_finder.md, whose steps the private functions
here are named after. Positions are 0-based and half-open, as VcfRecord.start and
VcfRecord.span are.

Public: RegionType, the four types of a segment; GenomeSegment, one segment of a contig;
StrFinderParameters, the parameters of the finder; segment_contig, the segmented genome of
one contig; read_str_regions, the STR regions of a --str-regions BED; write_str_regions,
the BED of the `pvc_eval str-regions` command; DEFAULT_NUM_THREADS and MAX_NUM_THREADS,
how many contigs it segments at once when not told and the most it takes.
"""

import bisect
import sys
from collections.abc import Callable, Collection, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, fields
from enum import Enum, StrEnum
from importlib.metadata import version
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from pvc_eval._text_files import check_gzip_whole, open_temporary_output
from pvc_eval.bed import parse_bed
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta


class RegionType(StrEnum):
    """The type of a segment of the segmented genome (doc/spec.md §5.1)."""

    GENERIC = "generic"
    STR = "str"
    STR_BUNDLE = "str-bundle"
    SATELLITE = "satellite"


@dataclass(frozen=True, slots=True)
class GenomeSegment:
    """One segment of a contig, `[start, end)`, 0-based and half-open.

    region_type is None only for a region read from a --str-regions BED without types.
    motif is an STR segment's motif, in upper case, and None for the other three types.

    Raises ValueError when end is not after start.
    """

    contig: str
    start: int
    end: int
    region_type: RegionType | None
    motif: str | None

    def __post_init__(self) -> None:
        if self.end <= self.start:
            raise ValueError(
                f"a segment's end must be after its start: {self.contig}, start {self.start},"
                f" end {self.end}"
            )


_MAX_PERIOD = 6  # the longest motif the finder scans, doc/specs/str_finder.md step 1

# Bases per block of step 1 (_find_candidates): each period takes the scores whose later
# base is in a block at once, as arrays. On 4 Mb of chr1, tracemalloc's peak of
# _find_candidates was 9.2 MiB with blocks of 4,096 bases, 11.8 MiB with 65,536 and 21.0 MiB
# with 262,144, about 49 bytes per base of a block (CPython 3.15t, numpy 2.5.4, Apple M5 Pro,
# 21 October 2026).
_BLOCK_BASES = 1 << 16

# Every byte to its canonical base: A, C, G and T in either case to the upper-case base,
# anything else to 0, which no base equals (pvc tandem_repeat.rs CANONICAL).
_CANONICAL = bytes(
    ord(letter.upper()) if letter in "ACGTacgt" else 0 for letter in map(chr, range(256))
)


@dataclass(frozen=True, slots=True, kw_only=True)
class StrFinderParameters:
    """The parameters of the STR finder, with the defaults of doc/specs/str_finder.md,
    "Parameters".

    copy_floors holds the fewest copies a candidate needs, for periods 1 to 6 in that
    order.

    Raises ValueError naming the field, its value and its limit when a value is outside
    the spec's table.
    """

    scan_min_period: int = 1
    scan_max_period: int = 6
    min_period: int = 1
    max_period: int = 6
    match_reward: int = 2
    mismatch_penalty: int = 7
    scan_min_copies: int = 2
    copy_floors: tuple[int, int, int, int, int, int] = (5, 5, 4, 4, 4, 3)
    min_purity: float = 0.8
    min_score: int = 0
    bundle_distance: int = 15
    contig_end_flank: int = 15
    max_str_len: int = 100
    max_bundle_len: int = 100

    def __post_init__(self) -> None:
        for field in fields(self):
            if field.name not in ("copy_floors", "min_purity"):
                _check_whole_number(field.name, getattr(self, field.name))
        _check_period_range("scan_min_period", "scan_max_period", self)
        _check_period_range("min_period", "max_period", self)
        for name in ("match_reward", "mismatch_penalty", "scan_min_copies", "bundle_distance"):
            _check_at_least(name, getattr(self, name), 1)
        for name in ("contig_end_flank", "max_bundle_len"):
            _check_at_least(name, getattr(self, name), 0)
        _check_at_least("max_str_len", self.max_str_len, self.bundle_distance, "bundle_distance")
        _check_copy_floors(self.copy_floors)
        purity = self.min_purity
        if isinstance(purity, bool) or not isinstance(purity, int | float):
            raise ValueError(f"min_purity is {purity!r}; it must be a float or an int")  # noqa: TRY004 -- the spec gives ValueError for every value outside its table
        # The range alone refuses nan and the infinities, and compares any int exactly.
        if not 0 <= purity <= 1:
            raise ValueError(f"min_purity is {purity!r}; it must be a number from 0 to 1")


def _check_whole_number(name: str, value: object) -> None:
    """Raise ValueError when a parameter is not a whole number. A bool is not one."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} is {value!r}; it must be a whole number")  # noqa: TRY004 -- the spec gives ValueError for every value outside its table


def _check_at_least(name: str, value: int, limit: int, limit_name: str | None = None) -> None:
    """Raise ValueError when a parameter is under its limit. limit_name is the field that
    sets the limit, when one does."""
    if value < limit:
        shown = limit if limit_name is None else f"{limit_name} ({limit})"
        raise ValueError(f"{name} is {value!r}; it must be {shown} or more")


def _check_period_range(min_name: str, max_name: str, parameters: StrFinderParameters) -> None:
    """Raise ValueError unless 1 <= the min period <= the max period <= 6."""
    min_period = getattr(parameters, min_name)
    max_period = getattr(parameters, max_name)
    _check_at_least(min_name, min_period, 1)
    for name, period in ((min_name, min_period), (max_name, max_period)):
        if period > _MAX_PERIOD:
            raise ValueError(f"{name} is {period!r}; it must be {_MAX_PERIOD} or less")
    if min_period > max_period:
        raise ValueError(f"{min_name} ({min_period!r}) must not exceed {max_name} ({max_period!r})")


def _check_copy_floors(copy_floors: object) -> None:
    """Raise ValueError unless copy_floors is a tuple of six whole numbers, 0 or more."""
    if not isinstance(copy_floors, tuple):
        raise ValueError(f"copy_floors is {copy_floors!r}; it must be a tuple")  # noqa: TRY004 -- the spec gives ValueError for every value outside its table
    if len(copy_floors) != _MAX_PERIOD:
        raise ValueError(
            f"copy_floors is {copy_floors!r}; it must be six whole numbers, one for each"
            " period from 1 to 6"
        )
    for floor in copy_floors:
        if isinstance(floor, bool) or not isinstance(floor, int) or floor < 0:
            raise ValueError(
                f"copy_floors is {copy_floors!r}; each must be a whole number, 0 or more"
            )


# A stretch of a contig, `[start, end)`, 0-based and half-open.
type _Span = tuple[int, int]


@dataclass(frozen=True, slots=True, kw_only=True)
class _Candidate:
    """A repeat step 1 finds: the stretch `[start, end)`, 0-based and half-open, its
    period and its Ruzzo-Tompa score. pvc's RepeatInterval."""

    start: int
    end: int
    period: int
    score: int


# Positions of a contig or indices of scores, as an array.
type _Positions = NDArray[np.intp]

# One flag per position, score or candidate.
type _Flags = NDArray[np.bool_]


@dataclass(frozen=True, slots=True, kw_only=True)
class _CandidateArrays:
    """Candidates of one period as arrays, in order of start and so of end: the stretches
    `[starts[i], ends[i])` of the contig, 0-based and half-open, and their scores."""

    starts: _Positions
    ends: _Positions
    scores: NDArray[np.int64]

    def get_where(self, flags: _Flags) -> _CandidateArrays:
        """The candidates whose flag is True."""
        return _CandidateArrays(
            starts=self.starts[flags], ends=self.ends[flags], scores=self.scores[flags]
        )

    def get_slice(self, first: int, stop: int | None = None) -> _CandidateArrays:
        """The candidates `[first, stop)` of the order, 0-based and half-open."""
        return _CandidateArrays(
            starts=self.starts[first:stop],
            ends=self.ends[first:stop],
            scores=self.scores[first:stop],
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class _BlockScores:
    """The scores `[first, first + number of scores)` of one period, 0-based and half-open,
    as arrays. Score k compares the canonical base at k + period with the one at k, and
    boundary k is the place before score k.

    run_starts and run_ends hold the runs of matches, `[run_starts[i], run_ends[i])`, as
    indices from first, in order. totals holds the running total at each boundary from first
    to the end of the scores, totals[k] that at first + k.
    """

    first: int
    run_starts: _Positions
    run_ends: _Positions
    totals: NDArray[np.int64]


@dataclass(frozen=True, slots=True, kw_only=True)
class _OpenSubsequence:
    """A subsequence that later scores can still change, on the stack of Ruzzo and Tompa.

    before and through are the running totals before its first score and after its last;
    `[start, end)` the indices of its scores in the contig, 0-based and half-open. link is
    1 + the index on the stack of the nearest subsequence under it whose total before is
    smaller, or 0 when there is none.
    """

    before: int
    through: int
    start: int
    end: int
    link: int


@dataclass(slots=True, kw_only=True)
class _ScanState:
    """Where the scan of one period stands between two blocks. Not frozen: the call of
    _find_candidates that makes it changes it after each block, and nothing else sees it.

    position is the index of the first score not taken yet, 0-based; total is the running
    total at the boundary before it, and least the least running total at any boundary up
    to there. stack holds the subsequences still open, in Ruzzo and Tompa's order.
    """

    position: int
    total: int
    least: int
    stack: list[_OpenSubsequence]


def _find_candidates(
    bases: str, parameters: StrFinderParameters, *, block_bases: int = _BLOCK_BASES
) -> list[_Candidate]:
    """Step 1 and step 3 of doc/specs/str_finder.md: the candidates of a sequence that
    reach the copy floor of their period, by period and then start, 0-based and half-open.
    With copy floors of 0 it gives what pvc's find_tandem_repeats gives.

    The sequence is read in blocks of block_bases bases, every period in each block before
    the next, and what each period's scan needs from the blocks before is carried in a
    _ScanState: where it stands, its running total and least total, and the subsequences
    still open. A block takes a period's scores whose later base is in it, so that what is
    held for the scores is one block, whatever the scoring. The tests set block_bases to
    make many blocks of a short sequence.

    The candidates made final in a block are tested as main's pure-Python step 1 tests them
    (tests/pure_python_step_one.py): a candidate is tested against the drop rule once every
    candidate of the shorter periods that divide its own and that starts before its end is
    tested; the kept ones are held while a longer period may still need them. What is held
    beyond one block is the subsequences still open and the candidates that wait for them,
    as in main.

    It departs from pvc's tandem_repeat.rs, and from main's pure Python, in these places,
    and gives the same candidates as both:

    - pvc scores the whole contig at once; here, as in main, a block at a time.
    - pvc runs Ruzzo and Tompa one score at a time, and main one run of matches at a time.
      A subsequence never starts or ends inside a run of matches or of mismatches, since no
      score is 0. Here a run that starts on an empty stack, and whose next run also does, is
      a maximal scoring subsequence on its own: the mismatches before that next run bring
      the total to or under its total before, which makes it final as it is, and no later
      run can reach back past an empty stack. Those are found on arrays. Only the other runs
      go through Ruzzo and Tompa's rules 2 to 4 in Python, pvc's loop with its links.
    - pvc's drop rule scans the kept candidates of each shorter period for every candidate,
      and main's scans from a binary search. Here a kept candidate that overlaps at least
      half of a candidate contains its middle, and then covers it by half exactly when it
      is at least half as long; the kept candidates that contain a middle are found by
      binary search, for all the candidates of a block at once.
    - pvc and main test the motifs one by one; here every motif of a block at once.

    <NUMBERS>
    """
    canonical = np.frombuffer(bases.encode("ascii").translate(_CANONICAL), np.uint8)
    periods = [
        period
        for period in range(parameters.scan_min_period, parameters.scan_max_period + 1)
        if period < len(canonical)
    ]
    divisors = {p: [q for q in periods if q < p and p % q == 0] for p in periods}
    multiples = {q: [p for p in periods if p > q and p % q == 0] for q in periods}
    states = {period: _ScanState(position=0, total=0, least=0, stack=[]) for period in periods}
    no_candidates = _make_no_candidates()
    # Per period: the candidates made final and not yet tested against the shorter periods;
    # the kept ones a longer period may still need; the start before which every candidate
    # is tested; and the candidates that reach the copy floor, as objects.
    untested = dict.fromkeys(periods, no_candidates)
    kept = dict.fromkeys(periods, no_candidates)
    tested_to = dict.fromkeys(periods, 0)
    over_floor: dict[int, list[_Candidate]] = {period: [] for period in periods}
    for block_end in range(block_bases, len(canonical) + block_bases, block_bases):
        end = min(block_end, len(canonical))
        # The periods in increasing order, so that a period's divisors are tested first.
        for period in periods:
            state = states[period]
            final = _find_final_candidates(canonical, state, end, period, parameters)
            final = _drop_short_or_repeated(canonical, final, period, parameters)
            queue = _join_candidates(untested[period], final)
            testable_to = min((tested_to[q] for q in divisors[period]), default=len(canonical))
            testable = int(queue.ends.searchsorted(testable_to, side="right"))
            kept_now = queue.get_slice(0, testable)
            untested[period] = queue.get_slice(testable)
            for shorter in divisors[period]:
                kept_now = kept_now.get_where(~_find_covered(kept_now, kept[shorter]))
            if multiples[period]:
                kept[period] = _join_candidates(kept[period], kept_now)
            floored = _drop_under_copy_floor(kept_now, period, parameters)
            over_floor[period].extend(_make_candidate_objects(floored, period))
            # After the last block nothing is open and no candidate starts later; the length
            # is past every end, so the longer periods can test all that is left.
            if untested[period].starts.size:
                tested_to[period] = int(untested[period].starts[0])
            elif end == len(canonical):
                tested_to[period] = len(canonical)
            else:
                tested_to[period] = state.stack[0].start if state.stack else state.position
        for period in periods:
            if multiples[period]:
                still_needed_from = min(tested_to[p] for p in multiples[period])
                unneeded = int(kept[period].ends.searchsorted(still_needed_from, side="right"))
                kept[period] = kept[period].get_slice(unneeded)
    return [candidate for period in periods for candidate in over_floor[period]]


def _make_no_candidates() -> _CandidateArrays:
    """Candidate arrays that hold none."""
    # np.zeros is a function of the module, called once per contig.
    return _CandidateArrays(
        starts=np.zeros(0, np.intp), ends=np.zeros(0, np.intp), scores=np.zeros(0, np.int64)
    )


def _join_candidates(first: _CandidateArrays, second: _CandidateArrays) -> _CandidateArrays:
    """The candidates of first, then those of second, as one set of arrays."""
    # np.concatenate is a function of the module, and the only way to join two arrays.
    return _CandidateArrays(
        starts=np.concatenate((first.starts, second.starts)),
        ends=np.concatenate((first.ends, second.ends)),
        scores=np.concatenate((first.scores, second.scores)),
    )


def _make_candidate_objects(candidates: _CandidateArrays, period: int) -> list[_Candidate]:
    """The candidates of one period, as arrays, made into _Candidate objects, in order."""
    return [
        _Candidate(start=start, end=end, period=period, score=score)
        for start, end, score in zip(
            candidates.starts.tolist(),
            candidates.ends.tolist(),
            candidates.scores.tolist(),
            strict=True,
        )
    ]


def _find_final_candidates(
    canonical: NDArray[np.uint8],
    state: _ScanState,
    end: int,
    period: int,
    parameters: StrFinderParameters,
) -> _CandidateArrays:
    """The maximal scoring subsequences of one period (Ruzzo and Tompa, 1999) that the
    block ending at base end makes final, as candidates in order of start: the scores
    `[k0, k1)` are the bases `[k0, k1 + period)`, 0-based and half-open. It takes the scores
    from state.position to those whose later base is before end, and moves state past them.

    A run starts on an empty stack when its total before is the least so far, the totals of
    the blocks before included. A run that does, and whose next run also does, or which is
    the last run and the scores end on an empty stack, is a subsequence on its own, on
    arrays. The others go through _find_grouped_subsequences with the stack. When the scores
    end on an empty stack, or the contig ends, every subsequence still open is final. A run
    that goes on past end is taken in two pieces; the second, which does not start on an
    empty stack, absorbs the first by rule 4 in the next block.
    """
    scores = _make_block_scores(
        canonical, state.position, max(end - period, state.position), period, state, parameters
    )
    totals = scores.totals
    lowest = np.minimum.accumulate(totals)
    # Whether the stack is empty at each run's start, and at the end of the scores.
    befores = totals[scores.run_starts]
    opens = (befores == lowest[scores.run_starts]) & (befores <= state.least)
    closes = end == len(canonical) or bool(totals[-1] == lowest[-1] and totals[-1] <= state.least)
    on_its_own = opens.copy()  # whether a subsequence starts and ends at the run
    on_its_own[:-1] &= opens[1:]
    if on_its_own.size:
        on_its_own[-1] &= closes
    run_starts = scores.run_starts
    run_ends = scores.run_ends
    grouped = (~on_its_own).nonzero()[0]
    grouped_final = _find_grouped_subsequences(
        state.stack,
        (run_starts[grouped] + scores.first).tolist(),
        (run_ends[grouped] + scores.first).tolist(),
        opens[grouped].tolist(),
        totals[run_starts[grouped]].tolist(),
        totals[run_ends[grouped]].tolist(),
        closes,
    )
    state.position = scores.first + len(totals) - 1
    state.total = int(totals[-1])
    state.least = min(state.least, int(lowest[-1]))
    final = _join_candidates(
        _CandidateArrays(
            starts=run_starts[on_its_own] + scores.first,
            ends=run_ends[on_its_own] + scores.first,
            scores=totals[run_ends[on_its_own]] - totals[run_starts[on_its_own]],
        ),
        grouped_final,
    )
    order = final.starts.argsort(kind="stable")
    return _CandidateArrays(
        starts=final.starts[order], ends=final.ends[order] + period, scores=final.scores[order]
    )


def _make_block_scores(
    canonical: NDArray[np.uint8],
    first: int,
    stop: int,
    period: int,
    state: _ScanState,
    parameters: StrFinderParameters,
) -> _BlockScores:
    """The scores `[first, stop)` of one period, 0-based and half-open: score k is the match
    reward when the canonical bases at k + period and at k are the same A, C, G or T, and
    minus the mismatch penalty otherwise. The running total at first is state.total."""
    later = canonical[first + period : stop + period]
    # A flag per score, with a mismatch before the first and after the last, so that the
    # changes of flag pair up as the starts and ends of the runs of matches. np.zeros is a
    # function of the module, and the only way here to make the array.
    matches = np.zeros(stop - first + 2, np.bool_)
    matches[1:-1] = (later == canonical[first:stop]) & (later != 0)
    changes = (matches[1:] != matches[:-1]).nonzero()[0]
    # The flag before the first score counts as a mismatch, which the penalty added back
    # undoes: totals[0] is state.total.
    steps = matches[:-1] * (parameters.match_reward + parameters.mismatch_penalty)
    totals = (steps - parameters.mismatch_penalty).cumsum() + (
        state.total + parameters.mismatch_penalty
    )
    return _BlockScores(
        first=first, run_starts=changes[0::2], run_ends=changes[1::2], totals=totals
    )


def _find_grouped_subsequences(
    stack: list[_OpenSubsequence],
    starts: Sequence[int],
    ends: Sequence[int],
    opens: Sequence[bool],
    befores: Sequence[int],
    throughs: Sequence[int],
    closes: bool,
) -> _CandidateArrays:
    """The subsequences that runs of matches make final by Ruzzo and Tompa's rules, as
    arrays of their scores' starts and ends, `[start, end)` 0-based and half-open, and
    their scores, in no order. Each run is given by its scores, whether it starts on an
    empty stack, which makes every subsequence on the stack final, and the totals before and
    after it. stack holds the subsequences open before the first run, and is left holding
    those open after the last; when closes is true, those are final too. pvc
    tandem_repeat.rs, maximal_scoring_subsequences, one step per run."""
    final: list[_OpenSubsequence] = []
    for start, end, opens_stack, before, through in zip(
        starts, ends, opens, befores, throughs, strict=True
    ):
        if opens_stack:
            final.extend(stack)
            stack.clear()
        open_one, staying = _make_open_subsequence(stack, before, through, start, end)
        del stack[staying:]
        stack.append(open_one)
    if closes:
        final.extend(stack)
        stack.clear()
    # np.array is a function of the module, and the only way to make an array of a list.
    return _CandidateArrays(
        starts=np.array([one.start for one in final], np.intp),
        ends=np.array([one.end for one in final], np.intp),
        scores=np.array([one.through - one.before for one in final], np.int64),
    )


def _make_open_subsequence(
    stack: Sequence[_OpenSubsequence], before: int, through: int, start: int, end: int
) -> tuple[_OpenSubsequence, int]:
    """The open subsequence that a run of positive scores at the indices `[start, end)`,
    0-based and half-open, makes, and how many of the open subsequences under it stay on
    the stack: the ones above those it absorbs. Rules 2 to 4 of Ruzzo and Tompa, as pvc's
    maximal_scoring_subsequences gives them."""
    staying = len(stack)
    while True:
        reached = _find_smaller_before(stack, before, staying)
        if reached < 0:  # rule 2: no subsequence to reach
            link = 0
            break
        if stack[reached].through >= through:  # rule 3: a separate subsequence
            link = reached + 1
            break
        # Rule 4: it absorbs the subsequence reached and every one above it, and looks again.
        before, start = stack[reached].before, stack[reached].start
        staying = reached
    open_one = _OpenSubsequence(before=before, through=through, start=start, end=end, link=link)
    return open_one, staying


def _find_smaller_before(stack: Sequence[_OpenSubsequence], before: int, top: int) -> int:
    """The index of the highest open subsequence under top whose total before is under
    before, or -1. It hops down the links, past the subsequences whose total before is
    not."""
    link = top
    while link > 0:
        if stack[link - 1].before < before:
            return link - 1
        link = stack[link - 1].link
    return -1


def _drop_short_or_repeated(
    canonical: NDArray[np.uint8],
    candidates: _CandidateArrays,
    period: int,
    parameters: StrFinderParameters,
) -> _CandidateArrays:
    """Drop the candidates of one period with fewer than scan_min_copies copies, and those
    whose motif is not primitive (step 1)."""
    copies = (candidates.ends - candidates.starts) // period
    candidates = candidates.get_where(copies >= parameters.scan_min_copies)
    return candidates.get_where(_find_primitive(canonical, candidates.starts, period))


def _find_primitive(canonical: NDArray[np.uint8], starts: _Positions, period: int) -> _Flags:
    """Whether each motif, the period canonical bases from each start, 0-based, is
    primitive: for no shorter period that divides it does every base equal the one that
    shorter period before it. Two bytes that are not A, C, G or T are both 0, and so equal.
    pvc tandem_repeat.rs, is_primitive_motif."""
    motif = [canonical[starts + offset] for offset in range(period)]
    repeated = np.zeros(len(starts), np.bool_)  # a function of the module, as in _make_block_scores
    for shorter in range(1, period):
        if period % shorter == 0:
            tiles = motif[shorter] == motif[0]
            for offset in range(shorter + 1, period):
                tiles &= motif[offset] == motif[offset - shorter]
            repeated |= tiles
    return ~repeated


def _find_covered(candidates: _CandidateArrays, kept: _CandidateArrays) -> _Flags:
    """Whether a kept candidate overlaps at least half of each candidate, 2 * overlap >=
    length (step 1's drop rule). kept are those of a shorter period that divides the
    candidates' own, in order of start and so of end, as one period's are.

    The overlap of a kept stretch [a, b) and a candidate [s, e) is the least of b - a,
    b - s, e - a and e - s. So it is at least half of e - s exactly when 2a <= s + e <= 2b,
    the kept stretch holds the candidate's middle, and b - a is at least half of e - s.
    The kept stretches that hold a middle are consecutive, found by binary search, and
    there are few: those of one period overlap by fewer bases than the period.
    """
    middles = candidates.starts + candidates.ends  # twice the middle of each candidate
    first = (kept.ends * 2).searchsorted(middles)
    stop = (kept.starts * 2).searchsorted(middles, side="right")
    lengths = candidates.ends - candidates.starts
    kept_lengths = kept.ends - kept.starts
    covered = np.zeros(len(middles), np.bool_)  # a function of the module, as in _make_block_scores
    for offset in range(int((stop - first).max(initial=0))):
        holds = first + offset < stop
        covered[holds] |= 2 * kept_lengths[first[holds] + offset] >= lengths[holds]
    return covered


def _drop_under_copy_floor(
    candidates: _CandidateArrays, period: int, parameters: StrFinderParameters
) -> _CandidateArrays:
    """Step 3 of doc/specs/str_finder.md: drop the candidates of one period whose detected
    span has fewer copies than the copy floor of the period, or whose end is not after
    their start (pvc segment_criteria.rs, prefilter). Positions are 0-based and half-open."""
    copies = (candidates.ends - candidates.starts) // period
    return candidates.get_where(
        (candidates.ends > candidates.starts) & (copies >= _get_copy_floor(parameters, period))
    )


def _get_copy_floor(parameters: StrFinderParameters, period: int) -> int:
    """The fewest copies a candidate of a period needs (pvc MinCopies::for_period)."""
    return parameters.copy_floors[period - 1]


def _get_span_bases(candidate: _Candidate, bases: str) -> str:
    """The bases of a candidate's detected span, in upper case, as pvc's upper gives them."""
    return bases[candidate.start : candidate.end].upper()


def _get_motif(candidate: _Candidate, bases: str) -> str:
    """A candidate's motif: the first period bases of its detected span, in upper case."""
    # Interned, so that the STR segments and bundle members of a motif share one string;
    # there are at most 5,460 motifs of periods 1 to 6, held for the life of the process.
    # Threads share the table safely: whichever interns a motif first, all get equal strings.
    # On chr1 it lowers the peak RSS from 362.0 to 333.5 MB on 64 Mb and from 128.1 to
    # 122.8 MB on 16 Mb, at the same wall time.
    return sys.intern(bases[candidate.start : candidate.start + candidate.period].upper())


def _find_cut_or_span(candidate: _Candidate, bases: str) -> _Span:
    """A candidate's cut, or its detected span when it has no cut, as positions of the
    contig whose bases are bases, 0-based and half-open."""
    cut = _find_cut(_get_span_bases(candidate, bases), _get_motif(candidate, bases))
    if cut is None:
        return (candidate.start, candidate.end)
    return (candidate.start + cut[0], candidate.start + cut[1])


def _drop_near_contig_ends(
    candidates: Iterable[_Candidate], bases: str, parameters: StrFinderParameters
) -> Iterator[_Candidate]:
    """Step 2 of doc/specs/str_finder.md: drop a candidate with fewer than
    contig_end_flank bases between it and either end of the contig (pvc row.rs,
    row_for_interval). bases are the whole contig, and positions 0-based and half-open.

    The distance is measured on the candidate's cut, or on its detected span when it has
    no cut. A contig_end_flank of 0 drops nothing.
    """
    flank = parameters.contig_end_flank
    for candidate in candidates:
        start, end = _find_cut_or_span(candidate, bases)
        if start >= flank and len(bases) - end >= flank:
            yield candidate


class _Rejection(Enum):
    """Why step 4 or step 6 turned a candidate down (pvc RejectionReason)."""

    COPY_FLOOR = "copy floor"  # pvc CopyFloor: under the floor, or shorter than the period
    COMPOUND = "compound"  # pvc Compound
    NO_CUT = "no cut"  # pvc NoCleanTrim
    PURITY = "purity"  # pvc Purity
    CONTIG_END = "contig end"  # pvc FlankClamped: the cut touches an end of the contig


@dataclass(frozen=True, slots=True, kw_only=True)
class _StrSegment:
    """What step 6 makes of a candidate: its cut, `[start, end)`, 0-based and half-open,
    its motif in upper case, and the purity of the cut. pvc's SsrSegment."""

    start: int
    end: int
    motif: str
    purity: float


@dataclass(frozen=True, slots=True, kw_only=True)
class _BundleMember:
    """A bundle member as pvc's partition describes it: its cut, `[start, end)`, 0-based
    and half-open, or its detected span when it has no cut, and its motif in upper case.
    pvc's BundleMember. It holds the cut, where _Classification.bundle_members hold the
    detected spans."""

    start: int
    end: int
    motif: str


@dataclass(frozen=True, slots=True, kw_only=True)
class _Classification:
    """What steps 4 to 6 make of the candidates of a contig, pvc's Classified: the STR
    segments, the bundle members, and each candidate turned down with its reason. The
    bundle members are the candidates themselves, so their detected spans, sorted by
    start and then end. A candidate outside the periods admitted, under min_score, or whose end
    is not after its start is in none of the three, as in pvc."""

    str_segments: tuple[_StrSegment, ...]
    bundle_members: tuple[_Candidate, ...]
    rejected: tuple[tuple[_Candidate, _Rejection], ...]


def _is_compound(motif: str) -> bool:
    """Whether a motif equals itself shifted by a length that divides its own
    (pvc segment_criteria.rs, is_compound)."""
    length = len(motif)
    return any(
        length % shift == 0 and motif[shift:] == motif[: length - shift]
        for shift in range(1, length)
    )


def _find_cut(span_bases: str, motif: str) -> _Span | None:
    """The cut of a detected span, given as its bases, as offsets into them, `[start,
    end)`, or None when it has none (doc/specs/str_finder.md step 6; pvc
    segment_criteria.rs, minimal_trim)."""
    two_copies = motif * 2
    reach = min(3 * len(motif), len(span_bases) // 2)
    start = None
    for offset in range(reach + 1):
        # pvc gives up before it compares, so a span of exactly two copies has no cut.
        if offset + len(two_copies) >= len(span_bases):
            return None
        if span_bases[offset : offset + len(two_copies)] == two_copies:
            start = offset
            break
    if start is None:
        return None
    for offset in range(len(span_bases), len(span_bases) - reach - 1, -1):
        # The two copies found at start end at an offset after start, which this search
        # reaches before any offset under two copies. So neither test below can turn the
        # cut down; both mirror pvc's minimal_trim.
        if offset < len(two_copies):
            return None
        if span_bases[offset - len(two_copies) : offset] == two_copies:
            return (start, offset) if start < offset else None
    return None


def _get_purity(cut_bases: str, motif: str) -> float:
    """The fraction of the bases of a cut that equal the motif repeated from its first
    base (pvc segment_criteria.rs, recompute_purity). The cut is never empty."""
    tiling = (motif * (len(cut_bases) // len(motif) + 1))[: len(cut_bases)]
    matches = sum(base == expected for base, expected in zip(cut_bases, tiling, strict=True))
    return matches / len(cut_bases)


def _is_in_cluster(reach: int, candidate: _Candidate, bundle_distance: int) -> bool:
    """Whether a candidate joins the cluster whose furthest end is reach: fewer than
    bundle_distance bases between them, an overlap counting as none (pvc
    segment_criteria.rs, joins_cluster)."""
    return max(candidate.start - reach, 0) < bundle_distance


def _find_bundles(
    candidates: Sequence[_Candidate], bundle_distance: int
) -> tuple[list[_Candidate], list[_Candidate]]:
    """Step 5 of doc/specs/str_finder.md: the candidates alone in their cluster, and those
    of the clusters of two or more, each in the order given (pvc segment_criteria.rs,
    split_bundles). candidates are sorted by start, then end."""
    isolated: list[_Candidate] = []
    bundled: list[_Candidate] = []
    for cluster in _make_clusters(candidates, bundle_distance):
        if len(cluster) > 1:
            bundled.extend(cluster)
        else:
            isolated.extend(cluster)
    return isolated, bundled


def _make_clusters(
    candidates: Sequence[_Candidate], bundle_distance: int
) -> list[list[_Candidate]]:
    """The clusters of step 5, in order: each candidate joins the cluster before it when
    _is_in_cluster says so. candidates are sorted by start, then end. Over the bundle
    members alone it gives the bundles of step 7 (pvc segment_criteria.rs,
    bundle_clusters)."""
    clusters: list[list[_Candidate]] = []
    reach = 0
    for candidate in candidates:
        if clusters and _is_in_cluster(reach, candidate, bundle_distance):
            clusters[-1].append(candidate)
            reach = max(reach, candidate.end)
        else:
            clusters.append([candidate])
            reach = candidate.end
    return clusters


def _make_bundle_member(candidate: _Candidate, bases: str) -> _BundleMember:
    """A bundle member with its cut, or its detected span when it has no cut, and its
    motif (pvc segment_criteria.rs, bundle_member_from_bases). bases are the whole contig,
    and positions 0-based and half-open."""
    start, end = _find_cut_or_span(candidate, bases)
    return _BundleMember(start=start, end=end, motif=_get_motif(candidate, bases))


def _find_str_segments(
    candidates: Iterable[_Candidate], bases: str, parameters: StrFinderParameters
) -> _Classification:
    """Steps 4 to 6 of doc/specs/str_finder.md over the candidates of a contig (pvc
    segment_criteria.rs, classify, over a whole contig). bases are the whole contig, in
    any case, and positions 0-based and half-open."""
    in_scope: list[_Candidate] = []
    rejected: list[tuple[_Candidate, _Rejection]] = []
    for candidate in candidates:
        if not parameters.min_period <= candidate.period <= parameters.max_period:
            continue
        if candidate.score < parameters.min_score or candidate.end <= candidate.start:
            continue
        if candidate.end - candidate.start < candidate.period:
            rejected.append((candidate, _Rejection.COPY_FLOOR))
        elif _is_compound(_get_motif(candidate, bases)):
            rejected.append((candidate, _Rejection.COMPOUND))
        else:
            in_scope.append(candidate)
    in_scope.sort(key=lambda candidate: (candidate.start, candidate.end))
    isolated, bundle_members = _find_bundles(in_scope, parameters.bundle_distance)
    str_segments: list[_StrSegment] = []
    for candidate in isolated:
        outcome = _make_str_segment(candidate, bases, parameters)
        if isinstance(outcome, _Rejection):
            rejected.append((candidate, outcome))
        else:
            str_segments.append(outcome)
    return _Classification(
        str_segments=tuple(str_segments),
        bundle_members=tuple(bundle_members),
        rejected=tuple(rejected),
    )


def _make_str_segment(
    candidate: _Candidate, bases: str, parameters: StrFinderParameters
) -> _StrSegment | _Rejection:
    """Step 6 of doc/specs/str_finder.md: the STR segment of a candidate alone in its
    cluster, or why there is none (pvc segment_criteria.rs, finish_locus). bases are the
    whole contig, and positions 0-based and half-open."""
    span_bases = _get_span_bases(candidate, bases)
    motif = _get_motif(candidate, bases)
    cut = _find_cut(span_bases, motif)
    if cut is None:
        return _Rejection.NO_CUT
    start, end = candidate.start + cut[0], candidate.start + cut[1]
    if (end - start) // candidate.period < _get_copy_floor(parameters, candidate.period):
        return _Rejection.COPY_FLOOR
    purity = _get_purity(span_bases[cut[0] : cut[1]], motif)
    if purity < parameters.min_purity:
        return _Rejection.PURITY
    if start == 0 or end == len(bases):
        return _Rejection.CONTIG_END
    return _StrSegment(start=start, end=end, motif=motif, purity=purity)


@dataclass(frozen=True, slots=True, kw_only=True)
class _Segment:
    """A segment of a contig as step 8 makes it, `[start, end)`, 0-based and half-open,
    with its type, an STR segment's motif, and an STR bundle's members, which
    GenomeSegment does not carry. pvc's TypedRegion."""

    start: int
    end: int
    region_type: RegionType
    motif: str | None = None
    members: tuple[_BundleMember, ...] = ()


_DEFAULT_PARAMETERS = StrFinderParameters()  # doc/specs/str_finder.md, "Parameters"


def segment_contig(
    reference: IndexedFasta, contig: str, parameters: StrFinderParameters = _DEFAULT_PARAMETERS
) -> Iterator[GenomeSegment]:
    """The segmented genome of one contig: its segments in position order, consecutive and
    non-overlapping from 0 to the contig's length, 0-based and half-open
    (doc/specs/str_finder.md, steps 1 to 8). A contig of no bases has none.

    It reads the whole contig once and resolves all of it when called, before the first
    segment comes out, since step 7 can grow a satellite to either side by any length. It
    keeps no state between calls, so threads can segment contigs of one IndexedFasta at
    once.

    Raises, when called, ValueError for a contig the .fai lacks or a closed reference, and
    PvcEvalInputError for a byte of the FASTA that is not a letter, as
    IndexedFasta.get_bases does.
    """
    bases = reference.get_bases(contig, 0, reference.get_length(contig))
    segments = _find_segments(bases, parameters)
    return (
        GenomeSegment(contig, segment.start, segment.end, segment.region_type, segment.motif)
        for segment in segments
    )


def _find_segments(bases: str, parameters: StrFinderParameters) -> Iterator[_Segment]:
    """Steps 1 to 8 of doc/specs/str_finder.md over the bases of a whole contig: its
    segments, consecutive from 0 to its length, 0-based and half-open (pvc
    region_typing/mod.rs, partition_resident_in, over the whole contig). Steps 1 to 7 run
    before it returns; step 8 makes the segments as they are read."""
    candidates = list(
        _drop_near_contig_ends(_find_candidates(bases, parameters), bases, parameters)
    )
    classification = _find_str_segments(candidates, bases, parameters)
    features = _find_features(candidates, classification, bases, parameters)
    return _make_segments(features, len(bases))


def _find_features(
    candidates: Sequence[_Candidate],
    classification: _Classification,
    bases: str,
    parameters: StrFinderParameters,
) -> list[_Segment]:
    """Step 7 of doc/specs/str_finder.md: the satellites, STR segments and STR bundles of a
    contig, by start (pvc region_typing/mod.rs, resolve_features). candidates are those
    left after step 3, and bases the whole contig; positions are 0-based and half-open."""
    distance = parameters.bundle_distance
    satellites = [
        run for run in _make_repeat_runs(candidates) if run[1] - run[0] > parameters.max_str_len
    ]
    satellites, members, str_segments = _make_grown_satellites(
        satellites, classification.bundle_members, classification.str_segments, distance
    )
    features = [
        _Segment(start=start, end=end, region_type=RegionType.SATELLITE)
        for start, end in satellites
    ]
    features.extend(
        _Segment(
            start=segment.start, end=segment.end, region_type=RegionType.STR, motif=segment.motif
        )
        for segment in str_segments
    )
    for cluster in _make_clusters(members, distance):
        # The hull runs over the detected spans, from the first member's start.
        start, end = cluster[0].start, max(member.end for member in cluster)
        if end - start > parameters.max_bundle_len:
            features.append(_Segment(start=start, end=end, region_type=RegionType.SATELLITE))
        elif start != 0 and end != len(bases):
            members_found = tuple(_make_bundle_member(member, bases) for member in cluster)
            features.append(
                _Segment(
                    start=start, end=end, region_type=RegionType.STR_BUNDLE, members=members_found
                )
            )
    features.sort(key=lambda feature: feature.start)
    return features


def _make_runs(spans: Iterable[_Span]) -> list[_Span]:
    """The spans, those that overlap or abut joined into one, by start; 0-based and
    half-open (pvc region_typing/mod.rs, merge_runs)."""
    runs: list[_Span] = []
    for start, end in sorted(spans):
        if runs and start <= runs[-1][1]:
            runs[-1] = (runs[-1][0], max(runs[-1][1], end))
        else:
            runs.append((start, end))
    return runs


def _make_repeat_runs(candidates: Iterable[_Candidate]) -> list[_Span]:
    """The repeat runs of doc/specs/str_finder.md step 7: the detected spans of the
    candidates, those that overlap or abut joined into one, by start; 0-based and
    half-open (pvc region_typing/mod.rs, coverage_runs)."""
    # After step 3 every candidate ends after its start; the test mirrors pvc's
    # coverage_runs.
    return _make_runs(
        (candidate.start, candidate.end)
        for candidate in candidates
        if candidate.end > candidate.start
    )


def _make_satellites_absorbing(
    satellites: Sequence[_Span], span: _Span, bundle_distance: int
) -> list[_Span] | None:
    """The satellites once a feature's span is absorbed, or None when no satellite reaches
    it: one overlapping it, or with fewer than bundle_distance bases between them. The
    span and the satellites it reaches become one satellite over all of them, which is then
    joined with any it overlaps or abuts (pvc region_typing/mod.rs, absorb_into).
    Positions are 0-based and half-open.

    The satellites are sorted, disjoint and not abutting, as _make_runs gives them, so both
    their starts and their ends ascend. A satellite reaches the span when its end is past
    start - bundle_distance and its start is before end + bundle_distance. The first
    condition holds from some satellite to the last, the second from the first satellite to
    some satellite, so the satellites that meet both are one slice of the list, found by
    two binary searches instead of a scan of every satellite."""
    start, end = span
    first = bisect.bisect_right(
        satellites, start - bundle_distance, key=lambda satellite: satellite[1]
    )
    stop = bisect.bisect_left(satellites, end + bundle_distance, key=lambda satellite: satellite[0])
    if first >= stop:
        return None
    grown = (min(start, satellites[first][0]), max(end, satellites[stop - 1][1]))
    return _make_runs([*satellites[:first], grown, *satellites[stop:]])


def _make_grown_satellites(
    satellites: Sequence[_Span],
    members: Sequence[_Candidate],
    str_segments: Sequence[_StrSegment],
    bundle_distance: int,
) -> tuple[list[_Span], list[_Candidate], list[_StrSegment]]:
    """Absorption, doc/specs/str_finder.md step 7: the satellites grown over every bundle
    member, by its detected span, and every STR segment, by its cut, that one reaches,
    repeated until a pass absorbs nothing; and the members and STR segments left, in their
    order (pvc's plain fixed point, its test helper absorb_by_iteration). Positions are
    0-based and half-open.

    The satellites are taken sorted, disjoint and not abutting, as _make_runs gives them,
    and each absorption gives them back so, since it goes through _make_runs again:
    _make_satellites_absorbing finds the satellites a span reaches by bisection, which
    relies on it."""
    grown = list(satellites)
    members_left = list(members)
    str_segments_left = list(str_segments)
    while True:
        is_absorbed = False
        kept_members: list[_Candidate] = []
        for member in members_left:
            after = _make_satellites_absorbing(grown, (member.start, member.end), bundle_distance)
            if after is None:
                kept_members.append(member)
            else:
                grown, is_absorbed = after, True
        kept_str_segments: list[_StrSegment] = []
        for segment in str_segments_left:
            after = _make_satellites_absorbing(grown, (segment.start, segment.end), bundle_distance)
            if after is None:
                kept_str_segments.append(segment)
            else:
                grown, is_absorbed = after, True
        members_left, str_segments_left = kept_members, kept_str_segments
        if not is_absorbed:
            return grown, members_left, str_segments_left


def _make_segments(features: Iterable[_Segment], length: int) -> Iterator[_Segment]:
    """Step 8 of doc/specs/str_finder.md: the features, sorted by start and not
    overlapping, with one generic region over each stretch between them and before the
    first and after the last, up to the contig's length (pvc region_typing/mod.rs,
    fill_generic_gaps). Positions are 0-based and half-open."""
    covered_to = 0
    for feature in features:
        if feature.start > covered_to:
            yield _Segment(start=covered_to, end=feature.start, region_type=RegionType.GENERIC)
        yield feature
        covered_to = feature.end
    if covered_to < length:
        yield _Segment(start=covered_to, end=length, region_type=RegionType.GENERIC)


def _find_requested_str_regions(
    segments: Iterable[GenomeSegment], requested: Iterable[_Span]
) -> Iterator[GenomeSegment]:
    """The STR regions of one contig's segments, given in position order, that overlap one
    of the contig's requested spans by at least one base, each whole, as `str-regions
    --regions` writes them (doc/specs/str_finder.md, "The str-regions command"; pvc
    region_typing/mod.rs, emit_into). Positions are 0-based and half-open."""
    # The spans joined and sorted, so that one pointer walks them as the segments advance:
    # a span that ends at or before a segment's start overlaps no later segment.
    spans = _make_runs(requested)
    index = 0
    for segment in segments:
        if segment.region_type is RegionType.GENERIC:
            continue
        while index < len(spans) and spans[index][1] <= segment.start:
            index += 1
        if index < len(spans) and spans[index][0] < segment.end:
            yield segment


# The region types a --str-regions BED may give, by the word of its fourth column
# (doc/specs/str_finder.md, read_str_regions). A generic region is not an STR region.
_STR_REGION_TYPES = {
    region_type.value: region_type
    for region_type in RegionType
    if region_type is not RegionType.GENERIC
}
_MOTIF_BASES = frozenset("ACGT")
_BYTE_ORDER_MARK = "\ufeff".encode()
_CORRECT_LINE = " Correct the line or remove it, and run again."


def read_str_regions(path: Path, reference: IndexedFasta) -> Iterator[GenomeSegment]:
    """The STR regions of a BED, one GenomeSegment per line, in the order of the file:
    the BED `str-regions` writes, or a user's (doc/specs/str_finder.md, read_str_regions).
    Positions are 0-based and half-open, as in the BED.

    A line of three columns gives region_type None and motif None. A fourth column is the
    region type, str, str-bundle or satellite, and a fifth the motif, upper-case A, C, G
    and T for str and `.` for the others; later columns are ignored. It reads UTF-8 and
    drops a byte-order mark; lines of only spaces and tabs, lines that start with #, and lines
    whose first field is track or browser are skipped; empty columns at the end of a line
    are no columns.

    Raises PvcEvalInputError, naming the file and the line, when the BED cannot be read, a
    line has fewer than three columns, a start or end is not a whole number, an end is not
    after its start, a contig is not in the reference or an end is past its length, a
    region type is unknown or a motif does not fit it, the BED is not sorted by contig in
    the order of the reference's `.fai` and then by start, or two regions overlap.
    """
    contig_lengths = _get_contig_lengths(reference)
    contig_ranks = {contig: rank for rank, contig in enumerate(reference.contigs)}
    seen_contigs: set[str] = set()
    previous: tuple[int, GenomeSegment] | None = None
    for line_number, columns in _read_bed_lines(path):
        where = f"{path}, line {line_number}"
        contig, start, end = _parse_bed_region(columns, where, contig_lengths)
        region_type, motif = _parse_type_and_motif(columns[3:5], where)
        segment = GenomeSegment(contig, start, end, region_type, motif)
        if previous is not None:
            _check_str_region_order(previous, segment, where, path, seen_contigs, contig_ranks)
        seen_contigs.add(contig)
        previous = (line_number, segment)
        yield segment


def _read_bed_lines(path: Path) -> Iterator[tuple[int, list[str]]]:
    """The columns of each line of a BED that holds a region, with its line number, from 1,
    as _split_bed_line gives them; the lines it gives no columns are skipped.

    Raises PvcEvalInputError when the file is missing, a directory, not readable or not
    UTF-8.
    """
    try:
        bed = path.open("rb")
    except FileNotFoundError:
        raise PvcEvalInputError(f"{path} does not exist. Check its path and run again.") from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{path} is a directory, not a BED file. Check its path and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{path} cannot be read: permission denied. Allow reading it and run again."
        ) from None
    with bed:
        for line_number, line_bytes in enumerate(bed, start=1):
            columns = _split_bed_line(line_bytes, line_number, path)
            if columns:
                yield line_number, columns


def _split_bed_line(line_bytes: bytes, line_number: int, path: Path) -> list[str]:
    """The tab-separated columns of a line of a BED, read as UTF-8, with a byte-order mark
    on the first line dropped and the empty columns at the end of the line taken away; none
    for a line of only spaces and tabs, a line that starts with #, and a line whose first
    field is track or browser (doc/specs/str_finder.md, read_str_regions).

    Raises PvcEvalInputError when the line is not UTF-8.
    """
    if line_number == 1:
        line_bytes = line_bytes.removeprefix(_BYTE_ORDER_MARK)
    try:
        text = line_bytes.decode("utf-8").rstrip("\r\n")
    except UnicodeDecodeError as err:
        raise PvcEvalInputError(
            f"{path} is not UTF-8 text: a byte at or after line {line_number} cannot be read as"
            f" UTF-8. If it is compressed, decompress it with `gunzip -c {path} > file.bed`,"
            " give file.bed instead, and run again."
        ) from err
    if not text.strip(" \t") or text.startswith("#"):
        return []
    if text.split(maxsplit=1)[:1] in (["track"], ["browser"]):
        return []
    columns = text.split("\t")
    while not columns[-1]:
        columns.pop()
    return columns


def _get_contig_lengths(reference: IndexedFasta) -> dict[str, int]:
    """The length of each contig of the reference, by name, in the order of the `.fai`."""
    return {contig: reference.get_length(contig) for contig in reference.contigs}


def _parse_bed_region(
    columns: Sequence[str], where: str, contig_lengths: Mapping[str, int]
) -> tuple[str, int, int]:
    """The contig, start and end of a BED line, 0-based and half-open, checked against
    the reference's contig_lengths. where names the file and the line, for the message.

    Raises PvcEvalInputError for fewer than three columns, a start or end that is not a
    whole number, an end not after its start, a contig not in the reference and an end
    past the contig's length.
    """
    if len(columns) < 3:
        noun = "column" if len(columns) == 1 else "columns"
        raise PvcEvalInputError(
            f"{where} has {len(columns)} tab-separated {noun}, and a BED line needs at least"
            f" three: contig, start and end.{_CORRECT_LINE}"
        )
    contig = columns[0]
    start = _parse_bed_position("start", columns[1], where)
    end = _parse_bed_position("end", columns[2], where)
    if end <= start:
        raise PvcEvalInputError(
            f"{where}: the end, {end}, is not after the start, {start}, and a region holds at"
            f" least one base.{_CORRECT_LINE}"
        )
    length = contig_lengths.get(contig)
    if length is None:
        raise PvcEvalInputError(
            f"{where}: contig {contig} is not in the reference. Use the reference the BED was"
            f" made for, or remove the lines of {contig}, and run again."
        )
    if end > length:
        raise PvcEvalInputError(_make_past_contig_end_message(where, end, contig, length))
    return contig, start, end


def _make_past_contig_end_message(where: str, end: int, contig: str, length: int) -> str:
    """The message for a BED line whose end is past its contig's length."""
    return (
        f"{where}: the end, {end}, is past the end of {contig}, which has {length} bases in the"
        " reference. Use the reference the BED was made for, or correct the line, and run"
        " again."
    )


def _parse_bed_position(name: str, column: str, where: str) -> int:
    """A start or end of a BED line, a whole number of 0 or more."""
    if not (column.isascii() and column.isdecimal()):
        raise PvcEvalInputError(
            f"{where}: the {name} is '{column}', and it must be a whole number of 0 or"
            f" more.{_CORRECT_LINE}"
        )
    return int(column)


def _parse_type_and_motif(
    columns: Sequence[str], where: str
) -> tuple[RegionType | None, str | None]:
    """The region type and the motif of a --str-regions BED line, from its fourth and fifth
    columns, those given of the two: None and None when it has neither, and a motif of None
    for every type but str.

    Raises PvcEvalInputError for an unknown type, a str without a motif of 1 to 6
    upper-case bases, and another type with a fifth column that is not `.`.
    """
    if not columns:
        return None, None
    type_column = columns[0]
    region_type = _STR_REGION_TYPES.get(type_column)
    if region_type is None:
        raise PvcEvalInputError(
            f"{where}: the region type is '{type_column}', and it must be str, str-bundle or"
            " satellite. Make the BED with `pvc_eval str-regions`, or correct the line, and"
            " run again."
        )
    motif = columns[1] if len(columns) > 1 else None
    if region_type is not RegionType.STR:
        if motif is not None and motif != ".":
            raise PvcEvalInputError(
                f"{where}: a region of type {type_column} has no motif, and its fifth column is"
                f" '{motif}'. Put '.' there, and run again."
            )
        return region_type, None
    if motif is None:
        raise PvcEvalInputError(
            f"{where}: an STR segment, of type str, needs its motif in the fifth column: 1 to"
            " 6 bases of A, C, G and T, in upper case. Add it, and run again."
        )
    if not (1 <= len(motif) <= _MAX_PERIOD and set(motif) <= _MOTIF_BASES):
        raise PvcEvalInputError(
            f"{where}: the motif of an STR segment is '{motif}', and it must be 1 to 6 bases"
            " of A, C, G and T, in upper case. Correct it, and run again."
        )
    return region_type, motif


def _check_str_region_order(
    previous: tuple[int, GenomeSegment],
    segment: GenomeSegment,
    where: str,
    path: Path,
    seen_contigs: Collection[str],
    contig_ranks: Mapping[str, int],
) -> None:
    """Raise PvcEvalInputError unless a region of a --str-regions BED comes after the one
    before it: on a later contig of the `.fai` that has had no line yet, or on the same
    contig, starting at or after the end of the one before (doc/spec.md §4)."""
    previous_line, previous_segment = previous
    sort_remedy = (
        " Sort it by contig, in the order of the reference's .fai, and by start, for example"
        f" with `bedtools sort -faidx <the reference's .fai> -i {path}`, and run again."
    )
    if segment.contig != previous_segment.contig:
        if segment.contig in seen_contigs:
            raise PvcEvalInputError(
                f"{where}: the lines of {segment.contig} are not all together; this one comes"
                f" after lines of {previous_segment.contig}.{sort_remedy}"
            )
        if contig_ranks[segment.contig] < contig_ranks[previous_segment.contig]:
            raise PvcEvalInputError(
                f"{where}: contig {segment.contig} comes after {previous_segment.contig}, and the"
                f" reference's .fai puts {segment.contig} before it.{sort_remedy}"
            )
        return
    if segment.start < previous_segment.start:
        raise PvcEvalInputError(
            f"{where}: the start, {segment.start}, is before the start,"
            f" {previous_segment.start}, of line {previous_line} on {segment.contig}."
            f"{sort_remedy}"
        )
    if segment.start < previous_segment.end:
        raise PvcEvalInputError(
            f"{where}: the region {segment.contig} {segment.start} {segment.end} overlaps the"
            f" one of line {previous_line}, {previous_segment.contig}"
            f" {previous_segment.start} {previous_segment.end}, and STR regions must not"
            f" overlap. Merge them with `bedtools merge -i {path}`, and run again. The merged BED"
            " has three columns: its regions still widen the loci, but they have no region"
            " type, which `compare --strata pvc` reads."
        )


def _find_requested_spans(path: Path, reference: IndexedFasta) -> dict[str, list[_Span]]:
    """The regions of a --regions BED, read by pvc_eval.bed.parse_bed, by contig, each
    `(start, end)`, 0-based and half-open, in the order of the file, which need not be
    sorted (doc/specs/str_finder.md, "The str-regions command").

    Raises PvcEvalInputError for what parse_bed refuses, and, naming the file and the line,
    for an end past its contig's length, or for a damaged gzip BED in its place (doc/spec.md
    §5.2).
    """
    contig_lengths = _get_contig_lengths(reference)
    requested: dict[str, list[_Span]] = {}
    for bed_line in parse_bed(path, reference.contigs):
        length = contig_lengths[bed_line.contig]
        if bed_line.end > length:
            where = f"{path}, line {bed_line.line_number}"
            error = PvcEvalInputError(
                _make_past_contig_end_message(where, bed_line.end, bed_line.contig, length)
            )
            check_gzip_whole(error, path, "BED", bed_line.line_number)
            raise error
        requested.setdefault(bed_line.contig, []).append((bed_line.start, bed_line.end))
    return requested


# How many contigs write_str_regions segments at once when not told, and the most it takes:
# the VCF reader's default and limit, DEFAULT_NUM_THREADS and _MAX_NUM_THREADS of
# pvc_eval.vcf (doc/spec.md, "Threads"), as doc/specs/str_finder.md, "The str-regions
# command", gives them.
DEFAULT_NUM_THREADS = 4
MAX_NUM_THREADS = 256


def write_str_regions(
    reference: IndexedFasta,
    output: Path,
    parameters: StrFinderParameters,
    regions: Path | None,
    num_threads: int,
) -> None:
    """Write the BED of `pvc_eval str-regions` (doc/specs/str_finder.md, "The str-regions
    command"): a ## line for pvc_eval's version and one for each parameter, a #chrom line,
    and one row per STR region, contig, start, end, type and motif (`.` but for str),
    0-based and half-open, in the order of the `.fai` and then by start.

    Every contig of the `.fai` is segmented whole, or with a --regions BED only the contigs
    it names, and then only the STR regions that overlap one of its regions by at least
    one base are written, each whole. The BED is written to `<output>.tmp` and renamed to
    output when complete. A run stopped by an error or by Ctrl-C removes the `.tmp` and leaves
    output as it was; one killed by a signal leaves the `.tmp` behind.

    Up to num_threads contigs are segmented at once, each on a thread of its own, or one
    after another in this thread when num_threads is 1; the BED is the same at any number.
    With more than one, an error in a contig, or Ctrl-C, stops the run once the contigs being
    segmented have finished.

    Raises ValueError, before anything is read or written, for a num_threads that is not a
    whole number from 1 to 256. Raises PvcEvalInputError for a --regions BED that
    pvc_eval.bed.parse_bed refuses, an empty one among them, and for one with an end past
    its contig's length; for a byte of the FASTA that is not a letter, the first in the
    order of the `.fai`; and, before any contig is segmented, for an output that is a
    directory or whose `.tmp` cannot be made.
    """
    _check_whole_number("num_threads", num_threads)
    if not 1 <= num_threads <= MAX_NUM_THREADS:
        raise ValueError(
            f"num_threads is {num_threads!r}; it must be a whole number from 1 to {MAX_NUM_THREADS}"
        )
    if regions is None:
        # Without --regions, each whole contig is the one span requested of it.
        requested = {
            contig: [(0, length)] for contig, length in _get_contig_lengths(reference).items()
        }
    else:
        requested = _find_requested_spans(regions, reference)
    temporary = Path(f"{output}.tmp")
    bed = open_temporary_output(output, temporary, "--output", compress=False)
    try:
        with bed:
            bed.writelines(_format_str_regions_header(parameters))
            contigs = [contig for contig in reference.contigs if contig in requested]

            def make_rows(contig: str) -> list[str]:
                return _make_str_region_rows(reference, contig, parameters, requested[contig])

            for rows in _map_in_order(make_rows, contigs, num_threads):
                bed.writelines(rows)
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _make_str_region_rows(
    reference: IndexedFasta,
    contig: str,
    parameters: StrFinderParameters,
    requested: Iterable[_Span],
) -> list[str]:
    """The rows of the `str-regions` BED for one contig, segmented whole: its STR regions
    that overlap one of the requested spans, 0-based and half-open, each row with its line
    end. A thread hands back these rows and not the segments, so that a contig that waits
    for the ones before it to be written holds its rows only."""
    segments = segment_contig(reference, contig, parameters)
    str_regions = _find_requested_str_regions(segments, requested)
    return [_format_str_region(segment) for segment in str_regions]


def _map_in_order[T, R](
    function: Callable[[T], R], items: Iterable[T], num_threads: int
) -> Iterator[R]:
    """function over items, its results in the order of the items: in this thread when
    num_threads is 1, or else on a pool of num_threads threads, with at most two items each
    submitted ahead of the result given (doc/specs/str_finder.md, "Threads"). The same as
    pvc_eval.vcf's, which is private there."""
    if num_threads == 1:
        yield from map(function, items)
        return
    with ThreadPoolExecutor(max_workers=num_threads) as pool:
        yield from pool.map(function, items, buffersize=2 * num_threads)


def _format_str_regions_header(parameters: StrFinderParameters) -> list[str]:
    """The lines before the rows of the `str-regions` BED, each with its line end."""
    lines = [f"##pvc_eval str-regions {version('pvc-eval')}\n"]
    for field in fields(parameters):
        value = getattr(parameters, field.name)
        shown = ",".join(map(str, value)) if isinstance(value, tuple) else str(value)
        lines.append(f"##{field.name}: {shown}\n")
    lines.append("#chrom\tstart\tend\ttype\tmotif\n")
    return lines


def _format_str_region(segment: GenomeSegment) -> str:
    """One row of the `str-regions` BED, with its line end."""
    assert segment.region_type is not None  # segment_contig gives every segment a type
    motif = "." if segment.motif is None else segment.motif
    return (
        f"{segment.contig}\t{segment.start}\t{segment.end}\t{segment.region_type.value}\t{motif}\n"
    )
