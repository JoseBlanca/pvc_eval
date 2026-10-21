"""The pure-Python step 1 of the STR finder, as main had it at 96b3938, kept as the oracle of
the array version in src/pvc_eval/segments.py, as pvc keeps its fixed point for absorption as
the oracle of its sweep.

It is segments.py's step 1 and step 3 of doc/specs/str_finder.md moved here unchanged: the
contig is scored a block of _BLOCK_BASES bases at a time, every period side by side, and
Ruzzo and Tompa's loop takes one step per run of matches or mismatches. Its _find_candidates
gives what the array version's gives, list for list; the tests of segments.py compare the two,
and keep testing the pure functions here (pvc's Ruzzo-Tompa cases, the comparisons, the copy
floor of step 3). Positions are 0-based and half-open.
"""

import bisect
from collections import deque
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

from pvc_eval.segments import _CANONICAL, StrFinderParameters, _Candidate, _get_copy_floor

# Bases each period scores before the candidates found so far are tested and the ones
# nothing needs any more are let go (doc/specs/str_finder.md, "What it holds"). On 4 Mb of
# random bases, tracemalloc's peak was 8.0 MB with 4,096, 13.2 MB with 16,384 and 35.1 MB
# with 65,536, in 23 to 26 s each (CPython 3.15t, Apple M5 Pro, 19 October 2026).
_BLOCK_BASES = 1 << 12

# Step 1 compares each canonical base with the one period before it by an XOR, which is 0
# when the two bytes are equal. The earlier base is read through this table, which turns 0,
# a byte that is not A, C, G or T, into 1, and the later one keeps its 0. Read the same way,
# two N would both be 0 and match; read this way they never do, as doc/specs/str_finder.md
# step 1 says, while two equal A, C, G or T still XOR to 0.
_EARLIER_BASES = bytes([1, *range(1, 256)])

# Step 1 reads the XOR of the bases through this table, which turns a mismatch, any byte but
# 0, into 0, and a match into 1: a run of matches ends at the next 0 of what it gives.
_MISMATCHES_AS_ZEROS = bytes([1, *[0] * 255])


# A maximal scoring subsequence: the indices of its scores, `[start, end)`, and its score.
type _Subsequence = tuple[int, int, int]

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
    comparisons: Iterator[bytes], period: int, parameters: StrFinderParameters
) -> Iterator[tuple[list[_Subsequence], int | None]]:
    """Find the maximal scoring subsequences of one period's scores (Ruzzo and Tompa, 1999),
    given as consecutive blocks of comparisons: score k is the match reward where byte k of
    the blocks together is 0, a match, and minus the mismatch penalty elsewhere. pvc
    tandem_repeat.rs, maximal_scoring_subsequences.

    A subsequence is (start, end, score), `[start, end)` the indices of its scores in the
    whole sequence, 0-based and half-open. After each block it yields the subsequences that
    no later score can change, in order of start, and the index at which the first
    subsequence still open starts, or the next index when none is open: every subsequence
    yielded later starts there or after. After the last block it yields the subsequences
    left and None. A subsequence that cannot be a candidate of step 1 may be left out.

    A maximal scoring subsequence never starts or ends inside a run of matches or of
    mismatches, since no score is 0: a run of matches is taken whole or not at all, and a
    run of mismatches only lowers the total, which is lowest at the run's end. So the loop
    takes one step per run, and gives what one step per score gives.

    pvc holds every open subsequence to the end. Here, when the running total falls to or
    under the total before every open subsequence, those are final and are yielded: a later
    subsequence starts from a total that low, so it can never reach back into them.

    A run of matches that comes to an empty stack is held as plain numbers, and only while
    one is held is the stack empty. On the stack it would have no subsequence to reach, so
    it would be open as it is, with link 0, and the least total before would be its own. So
    the mismatches that bring the total to or under its total before make it final, as the
    stack would, without an open subsequence being made; any other mismatches, or the end
    of the block, put it on the stack as it would have been. On 4 Mb of chr1 this leaves
    112,962 open subsequences made instead of 4.6 million, and str-regions takes 2.31 s
    instead of 4.48 s (doc/reports/perf-str-finder.md, P2).

    On an empty stack, a run of matches shorter than both (scan_min_copies - 1) * period and
    mismatch_penalty // match_reward + 1, with a mismatch after it in its block, changes
    nothing that is a candidate: that one mismatch brings the total to or under the total
    before the run, which then is final on its own, with the period + its length bases,
    under scan_min_copies copies, and the stack is empty again. So on an empty stack the
    loop goes straight to the next run of matches at least that long, or to the run that
    reaches the block's end, which may go on in the next block. The running total is not
    kept over what it skips: on an empty stack, what is yielded later depends only on the
    differences of later totals. On 4 Mb of chr1 the loop takes 1.19 million steps, each a
    run of mismatches and the run of matches after it, instead of 9.2 million runs, and
    str-regions takes 1.23 s instead of 2.06 s (doc/reports/perf-str-finder.md, P8).
    """
    reward = parameters.match_reward
    penalty = parameters.mismatch_penalty
    # The fewest matches of a run that is not skipped on an empty stack, as zero bytes.
    least_run = bytes(min((parameters.scan_min_copies - 1) * period, penalty // reward + 1))
    stack: list[_OpenSubsequence] = []
    held: _HeldRun | None = None  # a run of matches over an empty stack
    lowest_before = 0  # the least total before any open subsequence, when there is one
    total = 0
    block_start = 0  # the index of the block's first score
    for compared in comparisons:
        length = len(compared)
        # Past the block's end the comparisons read as matches, and the mismatches as one
        # mismatch, so that each search below stops at the end at the latest, and a run of
        # matches that reaches the end is never skipped.
        mismatches = compared.translate(_MISMATCHES_AS_ZEROS) + b"\0"
        compared += least_run + b"\0"
        final: list[_Subsequence] = []
        position = 0  # the block's comparisons before it are taken
        while position < length:
            start = compared.find(0, position)
            total -= penalty * (start - position)
            if held is not None and total <= held[0]:
                held_before, through, held_start, held_end = held
                final.append((held_start, held_end, through - held_before))
                held = None
            elif held is not None:
                stack.append(_make_held_open(held))
                lowest_before = held[0]
                held = None
            elif stack and total <= lowest_before:
                final.extend(open_one.get_subsequence() for open_one in stack)
                stack.clear()
            if not stack:  # and nothing is held
                start = compared.find(least_run, start)
            if start == length:
                break
            end = mismatches.find(0, start)
            before = total
            total += reward * (end - start)
            if not stack:
                held = (before, total, block_start + start, block_start + end)
            else:
                open_one, staying = _make_open_subsequence(
                    stack, before, total, block_start + start, block_start + end
                )
                del stack[staying:]
                stack.append(open_one)
                # It takes the total before of a subsequence it absorbs, and absorbs none lower.
                lowest_before = (
                    open_one.before if staying == 0 else min(lowest_before, open_one.before)
                )
            position = end
        if held is not None:  # the next block and the index yielded see it on the stack
            stack.append(_make_held_open(held))
            lowest_before = held[0]
            held = None
        block_start += length
        yield final, (stack[0].start if stack else block_start)
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


def _make_comparisons(canonical: bytes, period: int, block_bases: int) -> Iterator[bytes]:
    """The comparisons of one period, in blocks that end where the blocks of block_bases
    bases end, so that every period gives the same number of blocks. Byte k of the blocks
    together is 0, a match, when the bases at k + period and at k are the same A, C, G or T,
    and not 0 otherwise; score k is the match reward or minus the mismatch penalty. Each
    block is made on its own, so what is held is one block.
    """
    for block_start in range(0, len(canonical), block_bases):
        first = max(block_start - period, 0)
        last = max(min(block_start + block_bases, len(canonical)) - period, first)
        earlier = canonical[first:last].translate(_EARLIER_BASES)
        later = canonical[first + period : last + period]
        # A 0 byte where the two bases match, from one XOR of two ints, not a loop per base.
        # With the runs, step 1 takes 9.2 million steps on 4 Mb of chr1 instead of 24 million,
        # and str-regions 4.49 s instead of 5.73 s (doc/reports/perf-str-finder.md, P1).
        yield (int.from_bytes(earlier) ^ int.from_bytes(later)).to_bytes(last - first)


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
            _make_comparisons(canonical, period, block_bases), period, parameters
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
