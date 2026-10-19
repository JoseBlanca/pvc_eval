"""The stratum of a locus (doc/spec.md §3 and §5.7, doc/specs/loci.md "Stratum"): with
--strata pvc, the one region type of the STR regions it holds; with a --strata BED, the
classes of the intervals it shares positions with.

Positions are 0-based and half-open, as a BED and a Locus give them.

Public: pvc_stratum, the stratum of --strata pvc; Strata, the intervals of a --strata BED;
read_strata, which reads them.
"""

from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

from pvc_eval.bed import parse_bed
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.locus import Locus
from pvc_eval.segments import RegionType

# §5.7: the first that applies, from the most complex STR region to the simplest.
_PVC_STRATA = (RegionType.SATELLITE, RegionType.STR_BUNDLE, RegionType.STR)

type _Interval = tuple[int, int, str]  # start, end and class of a line of a --strata BED
type _Piece = tuple[int, int, frozenset[str]]  # a stretch and the classes over all of it


def pvc_stratum(locus: Locus) -> str:
    """The stratum of a locus with --strata pvc (doc/spec.md §5.7): `satellite`,
    `str-bundle` or `str`, the type of its most complex STR region, or `generic` when it
    holds none.

    Raises PvcEvalInputError for an STR region with no region type, read from a
    --str-regions BED without a fourth column. A locus holding no STR region raises
    nothing, so the error comes at the first locus that holds one (decided by the owner
    on 19 October 2026).
    """
    region_types = set()
    for str_region in locus.str_regions:
        if str_region.region_type is None:
            raise PvcEvalInputError(
                f"The STR region {str_region.contig}:{str_region.start + 1}-{str_region.end}"
                " has no region type, since the --str-regions BED has no fourth column there,"
                " and --strata pvc needs the type of every STR region. Give --str-regions a"
                " BED written by `pvc_eval str-regions`, which has the types, or leave"
                " --str-regions out so that the STR finder runs, and run again."
            )
        region_types.add(str_region.region_type)
    for region_type in _PVC_STRATA:
        if region_type in region_types:
            return region_type.value
    return RegionType.GENERIC.value


class Strata:
    """The intervals of a --strata BED, by contig, each with its class. Positions are
    0-based and half-open.

    pieces: by contig, the stretches that one or more intervals cover, as (start, end,
        classes), sorted and not overlapping, each with the classes of the intervals over
        the whole of it.
    """

    __slots__ = ("_classes", "_ends", "_starts")

    def __init__(self, pieces: Mapping[str, Sequence[_Piece]]) -> None:
        self._starts = {
            contig: [piece[0] for piece in by_start] for contig, by_start in pieces.items()
        }
        self._ends = {
            contig: [piece[1] for piece in by_start] for contig, by_start in pieces.items()
        }
        self._classes = {
            contig: [piece[2] for piece in by_start] for contig, by_start in pieces.items()
        }

    def stratum(self, locus: Locus) -> str:
        """The classes of the intervals that share at least one position with the locus,
        each once, in alphabetical order and joined with `,`; `.` when there are none."""
        if locus.contig not in self._ends:
            return "."
        starts = self._starts[locus.contig]
        classes = self._classes[locus.contig]
        locus_classes: set[str] = set()
        piece = bisect_right(self._ends[locus.contig], locus.start)  # the first ending after
        while piece < len(starts) and starts[piece] < locus.end:
            locus_classes |= classes[piece]
            piece += 1
        return ",".join(sorted(locus_classes)) if locus_classes else "."


def read_strata(path: Path, contigs: Sequence[str]) -> Strata:
    """The intervals of a --strata BED, read whole, as read_regions reads a --regions BED;
    its fourth column is the class of each interval.

    contigs: the reference's contig names, in the `.fai`'s order.

    Raises PvcEvalInputError as parse_bed does, a line on a contig the reference lacks
    among them, and for a line with no fourth column.
    """
    intervals_by_contig: defaultdict[str, list[_Interval]] = defaultdict(list)
    for bed_line in parse_bed(path, contigs):
        if bed_line.name is None:
            raise PvcEvalInputError(
                f"{path}, line {bed_line.line_number}: has no fourth column, which names the"
                " class of its interval. Give each line of the --strata BED its class and run"
                " again."
            )
        intervals_by_contig[bed_line.contig].append((bed_line.start, bed_line.end, bed_line.name))
    return Strata(
        {contig: _make_pieces(intervals) for contig, intervals in intervals_by_contig.items()}
    )


def _make_pieces(intervals: Sequence[_Interval]) -> list[_Piece]:
    """The intervals of one contig cut where any of them starts or ends, so that each piece
    has the same classes over the whole of it; stretches with no interval are left out."""
    events = sorted(
        [(start, 1, name) for start, _, name in intervals]
        + [(end, -1, name) for _, end, name in intervals]
    )
    open_classes: Counter[str] = Counter()
    pieces: list[_Piece] = []
    for event_index, (position, change, name) in enumerate(events):
        open_classes[name] += change
        if not open_classes[name]:
            del open_classes[name]
        # The last event at a position: the classes now open hold up to the next one.
        is_last_here = event_index + 1 < len(events) and events[event_index + 1][0] > position
        if is_last_here and open_classes:
            pieces.append((position, events[event_index + 1][0], frozenset(open_classes)))
    return pieces
