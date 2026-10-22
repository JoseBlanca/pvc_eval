"""The selection of --regions (doc/spec.md §3, doc/specs/loci.md "Regions"): a locus is
compared when it shares at least one position with a region, and skipped otherwise; a locus
compared that has a position outside the regions is partly outside them (doc/spec.md §5.5).

Positions are 0-based and half-open, as a BED and a Locus give them.

Public: Regions, the regions of a BED, with overlaps and covers; read_regions, which reads
them; select_loci, which yields the loci they keep.
"""

from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path

from pvc_eval.bed import parse_bed
from pvc_eval.locus import Locus


class Regions:
    """Regions of the genome, by contig, with those that overlap or lie next to each other
    merged. Positions are 0-based and half-open.

    intervals: by contig, the regions as (start, end), in any order.
    """

    __slots__ = ("_ends", "_starts")

    def __init__(self, intervals: Mapping[str, Iterable[tuple[int, int]]]) -> None:
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        for contig, spans in intervals.items():
            merged_spans = _make_merged_spans(spans)
            self._starts[contig] = [start for start, _ in merged_spans]
            self._ends[contig] = [end for _, end in merged_spans]

    def overlaps(self, contig: str, start: int, end: int) -> bool:
        """Whether `[start, end)` of contig, 0-based and half-open, shares at least one
        position with a region."""
        if contig not in self._ends:
            return False
        ends = self._ends[contig]
        first_after_start = bisect_right(ends, start)  # the first region that ends after start
        return first_after_start < len(ends) and self._starts[contig][first_after_start] < end

    def covers(self, contig: str, start: int, end: int) -> bool:
        """Whether every position of `[start, end)` of contig, 0-based and half-open, lies
        in a region. Regions that lie next to each other are one, so a stretch across the
        edge between two of them is covered."""
        if contig not in self._ends:
            return False
        ends = self._ends[contig]
        first_after_start = bisect_right(ends, start)  # the only region that can hold start
        return (
            first_after_start < len(ends)
            and self._starts[contig][first_after_start] <= start
            and end <= ends[first_after_start]
        )


def read_regions(path: Path, contigs: Sequence[str]) -> Regions:
    """The regions of a --regions BED, read whole, since a BED such as GIAB's comes sorted
    by contig name rather than in the `.fai`'s order.

    contigs: the reference's contig names.

    Raises PvcEvalInputError as parse_bed does, a line on a contig the reference lacks
    among them.
    """
    spans_by_contig: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
    for bed_line in parse_bed(path, contigs):
        spans_by_contig[bed_line.contig].append((bed_line.start, bed_line.end))
    return Regions(spans_by_contig)


def select_loci(loci: Iterator[Locus], regions: Regions | None) -> Iterator[Locus]:
    """The loci, in their order, that share at least one position with the regions, and
    every locus when regions is None."""
    if regions is None:
        yield from loci
        return
    for locus in loci:
        if regions.overlaps(locus.contig, locus.start, locus.end):
            yield locus


def _make_merged_spans(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """The spans of one contig sorted, with those that overlap or lie next to each other
    merged into one."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged
