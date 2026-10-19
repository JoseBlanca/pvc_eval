"""The segmented genome, stage 1 of doc/spec.md §4: what the STR finder yields to the loci.

The STR finder is specified in doc/specs/str_finder.md, whose steps the private functions
here are named after. Positions are 0-based and half-open, as VcfRecord.start and
VcfRecord.span are.

Public: RegionType, the four types of a segment; GenomeSegment, one segment of a contig;
StrFinderParameters, the parameters of the finder.
"""

import bisect
from collections import deque
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, fields
from enum import Enum


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
# random bases, tracemalloc's peak was 8.0 MB with 4,096, 11.9 MB with 16,384 and 29.4 MB
# with 65,536, at the same time (CPython 3.15t, 19 October 2026).
_BLOCK_BASES = 1 << 12

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
            _check_at_least(name, getattr(self, name), 1, "1")
        for name in ("contig_end_flank", "max_bundle_len"):
            _check_at_least(name, getattr(self, name), 0, "0")
        _check_at_least("max_str_len", self.max_str_len, self.bundle_distance, "bundle_distance")
        _check_copy_floors(self.copy_floors)
        purity = self.min_purity
        is_number = isinstance(purity, int | float) and not isinstance(purity, bool)
        # The range alone refuses nan and the infinities, and compares any int exactly.
        if not (is_number and 0 <= purity <= 1):
            raise ValueError(f"min_purity is {purity!r}; it must be a float or int from 0 to 1")


def _check_whole_number(name: str, value: object) -> None:
    """Raise ValueError when a parameter is not a whole number. A bool is not one."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} is {value!r}; it must be a whole number")  # noqa: TRY004 -- the spec gives ValueError for every value outside its table


def _check_at_least(name: str, value: int, limit: int, limit_name: str) -> None:
    """Raise ValueError when a parameter is under its limit, named as limit_name."""
    if value < limit:
        shown = limit_name if limit_name == str(limit) else f"{limit_name}, {limit}"
        raise ValueError(f"{name} is {value!r}; it must be {shown} or more")


def _check_period_range(min_name: str, max_name: str, parameters: StrFinderParameters) -> None:
    """Raise ValueError unless 1 <= the min period <= the max period <= 6."""
    min_period = getattr(parameters, min_name)
    max_period = getattr(parameters, max_name)
    _check_at_least(min_name, min_period, 1, "1")
    if max_period > _MAX_PERIOD:
        raise ValueError(f"{max_name} is {max_period!r}; it must be {_MAX_PERIOD} or less")
    if min_period > max_period:
        raise ValueError(
            f"{min_name} is {min_period!r}; it must be {max_name}, {max_period}, or less"
        )


def _check_copy_floors(copy_floors: object) -> None:
    """Raise ValueError unless copy_floors is a tuple of six whole numbers, 0 or more."""
    if not isinstance(copy_floors, tuple) or len(copy_floors) != _MAX_PERIOD:
        raise ValueError(
            f"copy_floors is {copy_floors!r}; it must be a tuple of six whole numbers,"
            " one for each period from 1 to 6"
        )
    for floor in copy_floors:
        if isinstance(floor, bool) or not isinstance(floor, int) or floor < 0:
            raise ValueError(
                f"copy_floors is {copy_floors!r}; each must be a whole number, 0 or more"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class _Candidate:
    """A repeat step 1 finds: the stretch `[start, end)`, 0-based and half-open, its
    period and its Ruzzo-Tompa score. pvc's RepeatInterval."""

    start: int
    end: int
    period: int
    score: int


type _Segment = tuple[int, int, int]  # a maximal scoring subsequence: start, end, score


@dataclass(frozen=True, slots=True, kw_only=True)
class _OpenSegment:
    """A segment that later scores can still change, on the stack of Ruzzo and Tompa.

    before and through are the running totals before its first score and after its last;
    start and end the indices of those scores. link is 1 + the index on the stack of the
    nearest segment under it whose total before is smaller, or 0 when there is none.
    """

    before: int
    through: int
    start: int
    end: int
    link: int

    def get_segment(self) -> _Segment:
        """The segment as (start, end, score)."""
        return (self.start, self.end, self.through - self.before)


def _find_maximal_scoring_subsequences(
    score_blocks: Iterator[Iterable[int]],
) -> Iterator[tuple[list[_Segment], int | None]]:
    """Find the maximal scoring subsequences of a score sequence (Ruzzo and Tompa, 1999),
    given as consecutive blocks. pvc tandem_repeat.rs, maximal_scoring_subsequences.

    A segment is (start, end, score), start and end the inclusive indices of its first and
    last score in the whole sequence. After each block it yields the segments that no later
    score can change, in order of start, and the index at which the first segment still
    open starts, or the next index when none is open: every segment yielded later starts
    there or after. After the last block it yields the segments left and None.

    pvc holds every open segment to the end. Here, when the running total falls to or under
    the total before every open segment, those segments are final and are yielded: a later
    segment starts from a total that low, so it can never reach back into them.
    """
    stack: list[_OpenSegment] = []
    lowest_before = 0  # the least total before any open segment, when there is one
    total = 0
    index = 0
    for block in score_blocks:
        final: list[_Segment] = []
        for score in block:
            before = total
            total += score
            index += 1
            if score <= 0:
                if stack and total <= lowest_before:
                    final.extend(open_segment.get_segment() for open_segment in stack)
                    stack.clear()
                continue
            stack.append(_make_open_segment(stack, before, total, index - 1))
            # A merge takes the total before of a segment it cuts, and cuts none lower.
            newest_before = stack[-1].before
            lowest_before = newest_before if len(stack) == 1 else min(lowest_before, newest_before)
        yield final, (stack[0].start if stack else index)
    yield [open_segment.get_segment() for open_segment in stack], None


def _make_open_segment(
    stack: list[_OpenSegment], before: int, through: int, index: int
) -> _OpenSegment:
    """The open segment that a positive score at index makes, after it has absorbed the
    open segments it reaches, which are cut from the stack. Rules 2 to 4 of Ruzzo and Tompa,
    as pvc's maximal_scoring_subsequences gives them."""
    start = index
    while True:
        found = _find_smaller_before(stack, before)
        if found < 0:  # rule 2: no segment to reach
            return _OpenSegment(before=before, through=through, start=start, end=index, link=0)
        if stack[found].through >= through:  # rule 3: a separate segment
            return _OpenSegment(
                before=before, through=through, start=start, end=index, link=found + 1
            )
        # Rule 4: it takes in the segment found and every one above it, and looks again.
        before, start = stack[found].before, stack[found].start
        del stack[found:]


def _find_smaller_before(stack: list[_OpenSegment], before: int) -> int:
    """The index of the highest open segment whose total before is under before, or -1.
    It hops down the links, past the segments whose total before is not."""
    link = len(stack)
    while link > 0:
        if stack[link - 1].before < before:
            return link - 1
        link = stack[link - 1].link
    return -1


def _make_score_blocks(
    canonical: bytes, period: int, parameters: StrFinderParameters, block_bases: int
) -> Iterator[Iterator[int]]:
    """The scores of one period, in blocks that end where the blocks of block_bases bases
    end, so that every period gives the same number of blocks.

    The base at j >= period scores the match reward when it and the base period before it
    are the same A, C, G or T, and minus the mismatch penalty otherwise. Score k is the base
    at k + period.
    """
    reward = parameters.match_reward
    penalty = -parameters.mismatch_penalty
    for block_start in range(0, len(canonical), block_bases):
        first = max(block_start - period, 0)
        last = max(min(block_start + block_bases, len(canonical)) - period, first)
        back = canonical[first:last]
        here = canonical[first + period : last + period]
        yield (
            reward if base and base == earlier else penalty
            for earlier, base in zip(back, here, strict=True)
        )


def _is_primitive_motif(motif: bytes) -> bool:
    """Whether a motif is not a shorter motif repeated. Two bytes that are not A, C, G or T
    count as equal. pvc tandem_repeat.rs, is_primitive_motif."""
    length = len(motif)
    return not any(
        length % shorter == 0 and motif == motif[:shorter] * (length // shorter)
        for shorter in range(1, length)
    )


def _make_candidates(
    canonical: bytes, period: int, segments: Iterable[_Segment], parameters: StrFinderParameters
) -> Iterator[_Candidate]:
    """The candidates of one period's segments: the segment of scores k0 to k1 is the
    stretch `[k0, k1 + period + 1)`, kept when it has scan_min_copies copies and its motif
    is primitive."""
    for first_score, last_score, score in segments:
        start = first_score
        end = last_score + period + 1
        if (end - start) // period < parameters.scan_min_copies:
            continue
        if _is_primitive_motif(canonical[start : start + period]):
            yield _Candidate(start=start, end=end, period=period, score=score)


def _is_covered_by_half(candidate: _Candidate, kept: list[_Candidate]) -> bool:
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
        for period, (segments, open_from) in zip(periods, steps, strict=True):
            queue = untested[period]
            queue.extend(_make_candidates(canonical, period, segments, parameters))
            # A candidate is tested once every candidate of its divisors that starts before
            # its end is tested.
            testable_to = min((tested_to[q] for q in divisors[period]), default=len(canonical))
            while queue and queue[0].end <= testable_to:
                candidate = queue.popleft()
                if any(_is_covered_by_half(candidate, kept[q]) for q in divisors[period]):
                    continue
                if multiples[period]:
                    kept[period].append(candidate)
                copies = (candidate.end - candidate.start) // period
                if copies >= parameters.copy_floors[period - 1]:
                    found[period].append(candidate)
            no_later_start = len(canonical) if open_from is None else open_from
            tested_to[period] = queue[0].start if queue else no_later_start
        # A kept candidate is let go once it ends before every candidate of a multiple of
        # its period that is still to be tested.
        for period in periods:
            if multiples[period]:
                still_needed_from = min(tested_to[p] for p in multiples[period])
                needed = kept[period]
                del needed[: bisect.bisect_right(needed, still_needed_from, key=lambda c: c.end)]
    return [candidate for period in periods for candidate in found[period]]
