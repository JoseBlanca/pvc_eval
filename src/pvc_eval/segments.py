"""The segmented genome, stage 1 of doc/spec.md §4: what the STR finder yields to the loci.

The STR finder is specified in doc/specs/str_finder.md, whose steps the private functions
here are named after. Positions are 0-based and half-open, as VcfRecord.start and
VcfRecord.span are.

Public: RegionType, the four types of a segment; GenomeSegment, one segment of a contig;
StrFinderParameters, the parameters of the finder.
"""

import bisect
from collections import deque
from collections.abc import Iterable, Iterator, Sequence
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
# random bases, tracemalloc's peak was 8.0 MB with 4,096, 13.2 MB with 16,384 and 35.1 MB
# with 65,536, in 23 to 26 s each (CPython 3.15t, Apple M5 Pro, 19 October 2026).
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
            _check_at_least(name, getattr(self, name), 1)
        for name in ("contig_end_flank", "max_bundle_len"):
            _check_at_least(name, getattr(self, name), 0)
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
    if max_period > _MAX_PERIOD:
        raise ValueError(f"{max_name} is {max_period!r}; it must be {_MAX_PERIOD} or less")
    if min_period > max_period:
        raise ValueError(f"{min_name} ({min_period!r}) must not exceed {max_name} ({max_period!r})")


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


# A maximal scoring subsequence: the indices of its scores, `[start, end)`, and its score.
type _Subsequence = tuple[int, int, int]


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
    score_blocks: Iterator[Iterable[int]],
) -> Iterator[tuple[list[_Subsequence], int | None]]:
    """Find the maximal scoring subsequences of a score sequence (Ruzzo and Tompa, 1999),
    given as consecutive blocks. pvc tandem_repeat.rs, maximal_scoring_subsequences.

    A subsequence is (start, end, score), `[start, end)` the indices of its scores in the
    whole sequence, 0-based and half-open. After each block it yields the subsequences that
    no later score can change, in order of start, and the index at which the first
    subsequence still open starts, or the next index when none is open: every subsequence
    yielded later starts there or after. After the last block it yields the subsequences
    left and None.

    pvc holds every open subsequence to the end. Here, when the running total falls to or
    under the total before every open subsequence, those are final and are yielded: a later
    subsequence starts from a total that low, so it can never reach back into them.
    """
    stack: list[_OpenSubsequence] = []
    lowest_before = 0  # the least total before any open subsequence, when there is one
    total = 0
    index = 0
    for block in score_blocks:
        final: list[_Subsequence] = []
        for score in block:
            before = total
            total += score
            index += 1
            if score > 0:
                open_one, staying = _make_open_subsequence(stack, before, total, index - 1)
                del stack[staying:]
                stack.append(open_one)
                # It takes the total before of a subsequence it absorbs, and absorbs none lower.
                lowest_before = (
                    open_one.before if staying == 0 else min(lowest_before, open_one.before)
                )
            elif stack and total <= lowest_before:
                final.extend(open_one.get_subsequence() for open_one in stack)
                stack.clear()
        yield final, (stack[0].start if stack else index)
    yield [open_one.get_subsequence() for open_one in stack], None


def _make_open_subsequence(
    stack: Sequence[_OpenSubsequence], before: int, through: int, index: int
) -> tuple[_OpenSubsequence, int]:
    """The open subsequence that a positive score at index makes, and how many of the open
    subsequences under it stay on the stack: the ones above those it absorbs. Rules 2 to 4
    of Ruzzo and Tompa, as pvc's maximal_scoring_subsequences gives them."""
    start = index
    staying = len(stack)
    while True:
        found = _find_smaller_before(stack, before, staying)
        end = index + 1
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


def _find_cut_or_span(candidate: _Candidate, bases: str) -> tuple[int, int]:
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


def _find_cut(span_bases: str, motif: str) -> tuple[int, int] | None:
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
