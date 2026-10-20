"""The segmented genome, stage 1 of doc/spec.md §4: what the STR finder yields to the loci.

The STR finder is specified in doc/specs/str_finder.md, whose steps the private functions
here are named after. Positions are 0-based and half-open, as VcfRecord.start and
VcfRecord.span are.

Public: RegionType, the four types of a segment; GenomeSegment, one segment of a contig;
StrFinderParameters, the parameters of the finder; segment_contig, the segmented genome of
one contig; read_str_regions, the STR regions of a --str-regions BED; write_str_regions,
the BED of the `pvc_eval str-regions` command.
"""

import bisect
from collections import deque
from collections.abc import Collection, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, fields
from enum import Enum
from importlib.metadata import version
from pathlib import Path
from typing import TextIO

from pvc_eval.bed import parse_bed
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta


class RegionType(Enum):
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

# Bases each period scores before the candidates found so far are tested and the ones
# nothing needs any more are let go (doc/specs/str_finder.md, "What it holds"). On 4 Mb of
# random bases, tracemalloc's peak was 8.0 MB with 4,096, 13.2 MB with 16,384 and 35.1 MB
# with 65,536, in 23 to 26 s each (CPython 3.15t, Apple M5 Pro, 19 October 2026).
_BLOCK_BASES = 1 << 12

# Every byte to its canonical base: A, C, G and T in either case to the upper-case base,
# anything else to 0, which no base equals (pvc tandem_repeat.rs CANONICAL).
_CANONICAL = bytes(
    ord(letter.upper()) if letter in "ACGTacgt" else 0 for letter in map(chr, range(256))
)

# Step 1 compares each canonical base with the one period before it by an XOR, which is 0
# when the two bytes are equal. The earlier base is read through this table, which turns 0,
# a byte that is not A, C, G or T, into 1, and the later one keeps its 0. Read the same way,
# two N would both be 0 and match; read this way they never do, as doc/specs/str_finder.md
# step 1 says, while two equal A, C, G or T still XOR to 0.
_EARLIER_BASES = bytes([1, *range(1, 256)])

# Step 1 reads the XOR of the bases through this table, which turns a mismatch, any byte but
# 0, into 0, and a match into 1: a run of matches ends at the next 0 of what it gives.
_MISMATCHES_AS_ZEROS = bytes([1, *[0] * 255])


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


# A maximal scoring subsequence: the indices of its scores, `[start, end)`, and its score.
type _Subsequence = tuple[int, int, int]

# A run of consecutive scores of one sign: their sum and how many they are.
type _ScoreRun = tuple[int, int]

# A run of positive scores held over an empty stack: the running totals before and after
# it, and the indices of its scores, `[start, end)`.
type _HeldRun = tuple[int, int, int, int]


@dataclass(frozen=True, slots=True, kw_only=True)
class _OpenSubsequence:
    """A subsequence that later scores can still change, on the stack of Ruzzo and Tompa.

    before and through are the running totals before its first score and after its last;
    `[start, end)` the indices of its scores. link is 1 + the index on the stack of the
    nearest subsequence under it whose total before is smaller, or 0 when there is none.
    """

    before: int
    through: int
    start: int
    end: int
    link: int

    def get_subsequence(self) -> _Subsequence:
        """The subsequence as (start, end, score)."""
        return (self.start, self.end, self.through - self.before)


def _find_maximal_scoring_subsequences(
    score_blocks: Iterator[Iterable[_ScoreRun]],
) -> Iterator[tuple[list[_Subsequence], int | None]]:
    """Find the maximal scoring subsequences of a score sequence (Ruzzo and Tompa, 1999),
    given as consecutive blocks of runs of scores of one sign, no score 0. pvc
    tandem_repeat.rs, maximal_scoring_subsequences.

    A subsequence is (start, end, score), `[start, end)` the indices of its scores in the
    whole sequence, 0-based and half-open. After each block it yields the subsequences that
    no later score can change, in order of start, and the index at which the first
    subsequence still open starts, or the next index when none is open: every subsequence
    yielded later starts there or after. After the last block it yields the subsequences
    left and None.

    A maximal scoring subsequence never starts or ends inside a run of scores of one sign,
    since no score is 0: a run of positive scores is taken whole or not at all, and a run
    of negative ones only lowers the total, which is lowest at the run's end. So the loop
    takes one step per run, and gives what one step per score gives.

    pvc holds every open subsequence to the end. Here, when the running total falls to or
    under the total before every open subsequence, those are final and are yielded: a later
    subsequence starts from a total that low, so it can never reach back into them.

    A run of positive scores that comes to an empty stack is held as plain numbers, and
    only while one is held is the stack empty. On the stack it would have no subsequence to
    reach, so it would be open as it is, with link 0, and the least total before would be
    its own. So the next run that brings the total to or under its total before makes it
    final, as the stack would, without an open subsequence being made; any other run, or
    the end of the block, puts it on the stack as it would have been. On 4 Mb of chr1 this
    leaves 112,962 open subsequences made instead of 4.6 million, and str-regions takes 2.31 s
    instead of 4.48 s (doc/reports/perf-str-finder.md, P2).
    """
    stack: list[_OpenSubsequence] = []
    held: _HeldRun | None = None  # a run of positive scores over an empty stack
    lowest_before = 0  # the least total before any open subsequence, when there is one
    total = 0
    index = 0
    for block in score_blocks:
        final: list[_Subsequence] = []
        for score, length in block:
            before = total
            total += score
            index += length
            if held is not None and total <= held[0]:  # a positive run never gets here
                held_before, through, start, end = held
                final.append((start, end, through - held_before))
                held = None
                continue
            if held is not None:
                stack.append(_make_held_open(held))
                lowest_before = held[0]
                held = None
            if score > 0 and not stack:
                held = (before, total, index - length, index)
            elif score > 0:
                open_one, staying = _make_open_subsequence(
                    stack, before, total, index - length, index
                )
                del stack[staying:]
                stack.append(open_one)
                # It takes the total before of a subsequence it absorbs, and absorbs none lower.
                lowest_before = (
                    open_one.before if staying == 0 else min(lowest_before, open_one.before)
                )
            elif stack and total <= lowest_before:
                final.extend(open_one.get_subsequence() for open_one in stack)
                stack.clear()
        if held is not None:  # the next block and the index yielded see it on the stack
            stack.append(_make_held_open(held))
            lowest_before = held[0]
            held = None
        yield final, (stack[0].start if stack else index)
    yield [open_one.get_subsequence() for open_one in stack], None


def _make_held_open(held: _HeldRun) -> _OpenSubsequence:
    """The open subsequence of a held run, as the empty stack it is held over makes it."""
    before, through, start, end = held
    return _OpenSubsequence(before=before, through=through, start=start, end=end, link=0)


def _make_open_subsequence(
    stack: Sequence[_OpenSubsequence], before: int, through: int, start: int, end: int
) -> tuple[_OpenSubsequence, int]:
    """The open subsequence that a run of positive scores at the indices `[start, end)`,
    0-based and half-open, makes, and how many of the open subsequences under it stay on
    the stack: the ones above those it absorbs. Rules 2 to 4 of Ruzzo and Tompa, as pvc's
    maximal_scoring_subsequences gives them."""
    staying = len(stack)
    while True:
        found = _find_smaller_before(stack, before, staying)
        if found < 0:  # rule 2: no subsequence to reach
            link = 0
            break
        if stack[found].through >= through:  # rule 3: a separate subsequence
            link = found + 1
            break
        # Rule 4: it absorbs the subsequence found and every one above it, and looks again.
        before, start = stack[found].before, stack[found].start
        staying = found
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


def _make_score_blocks(
    canonical: bytes, period: int, parameters: StrFinderParameters, block_bases: int
) -> Iterator[list[_ScoreRun]]:
    """The scores of one period as runs of one sign, in blocks that end where the blocks of
    block_bases bases end, so that every period gives the same number of blocks. A run that
    crosses the end of a block is cut there.

    The base at j >= period scores the match reward when it and the base period before it
    are the same A, C, G or T, and minus the mismatch penalty otherwise. Score k is the base
    at k + period. Each block is translated on its own, so what is held is one block.
    """
    reward = parameters.match_reward
    penalty = -parameters.mismatch_penalty
    for block_start in range(0, len(canonical), block_bases):
        first = max(block_start - period, 0)
        last = max(min(block_start + block_bases, len(canonical)) - period, first)
        length = last - first
        earlier = canonical[first:last].translate(_EARLIER_BASES)
        later = canonical[first + period : last + period]
        # A 0 byte where the two bases match, from one XOR of two ints, not a loop per base.
        # With the runs, step 1 takes 9.2 million steps on 4 Mb of chr1 instead of 24 million,
        # and str-regions 4.49 s instead of 5.73 s (doc/reports/perf-str-finder.md, P1).
        compared = (int.from_bytes(earlier) ^ int.from_bytes(later)).to_bytes(length)
        # The runs are found with bytes.find, not with re.finditer: on 3.15t, threads that make
        # match objects at once slow each other, and find makes none. On 4 Mb of chr1, with
        # one segment_contig per thread, eight threads took 5.0 to 5.7 s each with finditer and
        # 2.9 to 3.8 s with find, and one thread 2.3 s and 2.1 s
        # (doc/reports/perf-str-finder.md, P3).
        mismatches = compared.translate(_MISMATCHES_AS_ZEROS)
        runs: list[_ScoreRun] = []
        mismatches_from = 0
        start = compared.find(0)
        while start >= 0:
            end = mismatches.find(0, start)
            if end < 0:
                end = length
            if start > mismatches_from:
                runs.append((penalty * (start - mismatches_from), start - mismatches_from))
            runs.append((reward * (end - start), end - start))
            mismatches_from = end
            start = compared.find(0, end)
        if mismatches_from < length:
            runs.append((penalty * (length - mismatches_from), length - mismatches_from))
        yield runs


def _is_primitive_motif(motif: bytes) -> bool:
    """Whether a motif is not a shorter motif repeated. Two bytes that are not A, C, G or T
    count as equal. pvc tandem_repeat.rs, is_primitive_motif."""
    length = len(motif)
    return not any(
        length % shorter == 0 and motif == motif[:shorter] * (length // shorter)
        for shorter in range(1, length)
    )


def _make_candidates(
    canonical: bytes,
    period: int,
    subsequences: Iterable[_Subsequence],
    parameters: StrFinderParameters,
) -> Iterator[_Candidate]:
    """The candidates of one period's subsequences: the scores `[k0, k1)` are the bases
    `[k0, k1 + period)`, since score k is the base at k + period and compares it with the
    base at k. A candidate is kept when it has scan_min_copies copies and its motif is
    primitive."""
    for start, last_score_end, score in subsequences:
        end = last_score_end + period
        if (end - start) // period < parameters.scan_min_copies:
            continue
        if _is_primitive_motif(canonical[start : start + period]):
            yield _Candidate(start=start, end=end, period=period, score=score)


def _is_covered_by_half(candidate: _Candidate, kept: Sequence[_Candidate]) -> bool:
    """Whether a kept candidate of a shorter period overlaps at least half of a candidate.
    kept is in order of start, and so of end, as one period's candidates are."""
    length = candidate.end - candidate.start
    first = bisect.bisect_right(kept, candidate.start, key=lambda other: other.end)
    for index in range(first, len(kept)):
        other = kept[index]
        if other.start >= candidate.end:
            return False
        overlap = min(other.end, candidate.end) - max(other.start, candidate.start)
        if overlap * 2 >= length:
            return True
    return False


def _drop_covered_candidates(
    queue: deque[_Candidate], testable_to: int, kept_shorter: Sequence[Sequence[_Candidate]]
) -> Iterator[_Candidate]:
    """Take from the queue, in order, the candidates that end at or before testable_to, and
    yield those that no kept candidate of kept_shorter overlaps by at least half (step 1's
    drop rule). kept_shorter holds the kept candidates of each shorter period that divides
    theirs."""
    while queue and queue[0].end <= testable_to:
        candidate = queue.popleft()
        if not any(_is_covered_by_half(candidate, kept) for kept in kept_shorter):
            yield candidate


def _evict_unneeded(kept: list[_Candidate], still_needed_from: int) -> None:
    """Remove from the kept candidates of one period those that end at or before
    still_needed_from, where every candidate of a longer period still to be tested starts.
    kept is in order of end."""
    del kept[: bisect.bisect_right(kept, still_needed_from, key=lambda candidate: candidate.end)]


def _find_candidates(
    bases: str, parameters: StrFinderParameters, *, block_bases: int = _BLOCK_BASES
) -> list[_Candidate]:
    """Step 1 and step 3 of doc/specs/str_finder.md: the candidates of a sequence that
    reach the copy floor of their period, by period and then start, 0-based and half-open.

    With copy floors of 0 it gives what pvc's find_tandem_repeats gives. A candidate is
    dropped when a kept candidate of a shorter period that divides its own overlaps at
    least half of it. A candidate under its floor is still kept for that, and let go once
    every candidate it could drop is known. All periods are scanned side by side, a block
    of block_bases bases at a time, so that what is held does not grow with the sequence;
    the tests set block_bases to make many blocks of a short sequence.
    """
    canonical = bases.encode("ascii").translate(_CANONICAL)
    periods = [
        period
        for period in range(parameters.scan_min_period, parameters.scan_max_period + 1)
        if period < len(canonical)
    ]
    divisors = {p: [q for q in periods if q < p and p % q == 0] for p in periods}
    multiples = {q: [p for p in periods if p > q and p % q == 0] for q in periods}
    # Per period: the candidates found and not yet tested against the shorter periods; the
    # kept ones a longer period may still need; the start before which every candidate is
    # tested; and the ones that reach the copy floor.
    untested: dict[int, deque[_Candidate]] = {period: deque() for period in periods}
    kept: dict[int, list[_Candidate]] = {period: [] for period in periods}
    tested_to = dict.fromkeys(periods, 0)
    found: dict[int, list[_Candidate]] = {period: [] for period in periods}
    streams = [
        _find_maximal_scoring_subsequences(
            _make_score_blocks(canonical, period, parameters, block_bases)
        )
        for period in periods
    ]
    for steps in zip(*streams, strict=True):
        # The periods in increasing order, so that a period's divisors are tested first.
        for period, (subsequences, open_from) in zip(periods, steps, strict=True):
            queue = untested[period]
            queue.extend(_make_candidates(canonical, period, subsequences, parameters))
            # A candidate is tested once every candidate of its divisors that starts before
            # its end is tested.
            testable_to = min((tested_to[q] for q in divisors[period]), default=len(canonical))
            kept_shorter = [kept[q] for q in divisors[period]]
            # The queue's candidates that can be tested now: those of one block, and more
            # only after a block in which a divisor could not be tested.
            kept_now = list(_drop_covered_candidates(queue, testable_to, kept_shorter))
            if multiples[period]:
                kept[period].extend(kept_now)
            found[period].extend(_drop_under_copy_floor(kept_now, parameters))
            # After the last block nothing is open and no candidate starts later; the length
            # is past every end, so the longer periods can test all that is left.
            no_later_start = len(canonical) if open_from is None else open_from
            tested_to[period] = queue[0].start if queue else no_later_start
        for period in periods:
            if multiples[period]:
                _evict_unneeded(kept[period], min(tested_to[p] for p in multiples[period]))
    return [candidate for period in periods for candidate in found[period]]


def _get_copy_floor(parameters: StrFinderParameters, period: int) -> int:
    """The fewest copies a candidate of a period needs (pvc MinCopies::for_period)."""
    return parameters.copy_floors[period - 1]


def _get_span_bases(candidate: _Candidate, bases: str) -> str:
    """The bases of a candidate's detected span, in upper case, as pvc's upper gives them."""
    return bases[candidate.start : candidate.end].upper()


def _get_motif(candidate: _Candidate, bases: str) -> str:
    """A candidate's motif: the first period bases of its detected span, in upper case."""
    return bases[candidate.start : candidate.start + candidate.period].upper()


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


def _drop_under_copy_floor(
    candidates: Iterable[_Candidate], parameters: StrFinderParameters
) -> Iterator[_Candidate]:
    """Step 3 of doc/specs/str_finder.md: drop a candidate whose detected span has fewer
    copies than the copy floor of its period, or whose end is not after its start (pvc
    segment_criteria.rs, prefilter)."""
    for candidate in candidates:
        copies = (candidate.end - candidate.start) // candidate.period
        if candidate.end > candidate.start and copies >= _get_copy_floor(
            parameters, candidate.period
        ):
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


def _find_segments(bases: str, parameters: StrFinderParameters) -> list[_Segment]:
    """Steps 1 to 8 of doc/specs/str_finder.md over the bases of a whole contig: its
    segments, consecutive from 0 to its length, 0-based and half-open (pvc
    region_typing/mod.rs, partition_resident_in, over the whole contig)."""
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


def _make_segments(features: Iterable[_Segment], length: int) -> list[_Segment]:
    """Step 8 of doc/specs/str_finder.md: the features, sorted by start and not
    overlapping, with one generic region over each stretch between them and before the
    first and after the last, up to the contig's length (pvc region_typing/mod.rs,
    fill_generic_gaps). Positions are 0-based and half-open."""
    segments: list[_Segment] = []
    covered_to = 0
    for feature in features:
        if feature.start > covered_to:
            segments.append(
                _Segment(start=covered_to, end=feature.start, region_type=RegionType.GENERIC)
            )
        segments.append(feature)
        covered_to = feature.end
    if covered_to < length:
        segments.append(_Segment(start=covered_to, end=length, region_type=RegionType.GENERIC))
    return segments


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
    for an end past its contig's length.
    """
    contig_lengths = _get_contig_lengths(reference)
    requested: dict[str, list[_Span]] = {}
    for bed_line in parse_bed(path, reference.contigs):
        length = contig_lengths[bed_line.contig]
        if bed_line.end > length:
            where = f"{path}, line {bed_line.line_number}"
            raise PvcEvalInputError(
                _make_past_contig_end_message(where, bed_line.end, bed_line.contig, length)
            )
        requested.setdefault(bed_line.contig, []).append((bed_line.start, bed_line.end))
    return requested


def write_str_regions(
    reference: IndexedFasta,
    output: Path,
    parameters: StrFinderParameters,
    regions: Path | None,
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

    Raises PvcEvalInputError for a --regions BED that pvc_eval.bed.parse_bed refuses, an
    empty one among them, and for one with an end past its contig's length; for a byte of the FASTA that is not a letter; and, before any contig is
    segmented, for an output that is a directory or whose `.tmp` cannot be made.
    """
    if regions is None:
        # Without --regions, each whole contig is the one span requested of it.
        requested = {
            contig: [(0, length)] for contig, length in _get_contig_lengths(reference).items()
        }
    else:
        requested = _find_requested_spans(regions, reference)
    temporary = Path(f"{output}.tmp")
    bed = _open_temporary_output(output, temporary)
    try:
        with bed:
            bed.writelines(_format_str_regions_header(parameters))
            for contig in [contig for contig in reference.contigs if contig in requested]:
                segments = segment_contig(reference, contig, parameters)
                str_regions = _find_requested_str_regions(segments, requested[contig])
                bed.writelines(_format_str_region(segment) for segment in str_regions)
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _open_temporary_output(output: Path, temporary: Path) -> TextIO:
    """Open the file the BED is written to before it is renamed to output.

    Raises PvcEvalInputError, naming --output, when output is a directory or temporary
    cannot be made.
    """
    if output.is_dir():
        raise PvcEvalInputError(
            f"--output {output} is a directory. Give the path of a file to write, and run again."
        )
    cannot = f"--output {output} cannot be written:"
    try:
        return temporary.open("w", encoding="utf-8")
    except FileNotFoundError:
        raise PvcEvalInputError(
            f"{cannot} its directory, {output.parent}, does not exist. Make it, or give another"
            " --output, and run again."
        ) from None
    except NotADirectoryError:
        raise PvcEvalInputError(
            f"{cannot} its path goes through a file as if it were a directory. Correct the"
            " path, and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{cannot} permission denied to write in {output.parent}. Allow writing there, or"
            " give another --output, and run again."
        ) from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{cannot} {temporary}, where it is written first, is a directory. Move it away,"
            " and run again."
        ) from None


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
