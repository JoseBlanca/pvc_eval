"""The selection of --regions (doc/spec.md §3, doc/specs/loci.md "Regions"): a locus is
compared when it shares at least one position with a region, and skipped otherwise.

Positions are 0-based and half-open, as a BED and a Locus give them.

Public: Regions, the regions of a BED; read_regions, which reads them; select_loci, which
yields the loci they keep.
"""

from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path

from pvc_eval.bed import parse_bed
from pvc_eval.locus import Locus


class Regions:
    """The regions of a --regions BED, by contig, with those that overlap or lie next to
    each other merged. Positions are 0-based and half-open.

    intervals: by contig, the regions as (start, end), sorted, with a gap of at least one
        position between each and the next.
    """

    __slots__ = ("_ends", "_starts")

    def __init__(self, intervals: Mapping[str, Sequence[tuple[int, int]]]) -> None:
        self._starts = {
            contig: [start for start, _ in spans] for contig, spans in intervals.items()
        }
        self._ends = {contig: [end for _, end in spans] for contig, spans in intervals.items()}

    def overlaps(self, contig: str, start: int, end: int) -> bool:
        """Whether `[start, end)` of contig, 0-based and half-open, shares at least one
        position with a region."""
        if contig not in self._ends:
            return False
        ends = self._ends[contig]
        first_after_start = bisect_right(ends, start)  # the first region that ends after start
        return first_after_start < len(ends) and self._starts[contig][first_after_start] < end


def read_regions(path: Path, contigs: Sequence[str]) -> Regions:
    """The regions of a --regions BED, read whole, since a BED such as GIAB's comes sorted
    by contig name rather than in the `.fai`'s order.

    contigs: the reference's contig names, in the `.fai`'s order.

    Raises PvcEvalInputError as parse_bed does, a line on a contig the reference lacks
    among them.
    """
    spans_by_contig: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
    for bed_line in parse_bed(path, contigs):
        spans_by_contig[bed_line.contig].append((bed_line.start, bed_line.end))
    return Regions({contig: _merge(spans) for contig, spans in spans_by_contig.items()})


def select_loci(loci: Iterator[Locus], regions: Regions | None) -> Iterator[Locus]:
    """The loci, in their order, that share at least one position with the regions, and
    every locus when regions is None."""
    if regions is None:
        yield from loci
        return
    for locus in loci:
        if regions.overlaps(locus.contig, locus.start, locus.end):
            yield locus


def _merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The spans of one contig sorted, with those that overlap or lie next to each other
    merged into one."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged
