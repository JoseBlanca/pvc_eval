"""The stratum of a locus (doc/spec.md §3 and §5.7, doc/specs/loci.md "Stratum"): with
--strata pvc, the one region type of the STR regions it holds; with a --strata BED, the
classes of the intervals it shares positions with.

Positions are 0-based and half-open, as a BED and a Locus give them.

Public: pvc_stratum, the stratum of --strata pvc; Strata, the intervals of a --strata BED;
read_strata, which reads them.
"""

from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

from pvc_eval._text_files import check_gzip_whole
from pvc_eval.bed import BedLine, parse_bed
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.locus import Locus
from pvc_eval.regions import Regions
from pvc_eval.segments import RegionType

# §5.7: the first that applies, from the most complex STR region to the simplest.
_PVC_STRATA = (RegionType.SATELLITE, RegionType.STR_BUNDLE, RegionType.STR)


def pvc_stratum(locus: Locus) -> str:
    """The stratum of a locus with --strata pvc (doc/spec.md §5.7): `satellite`,
    `str-bundle` or `str`, the type of its most complex STR region, or `generic` when it
    holds none.

    Raises PvcEvalInputError for an STR region with no region type, read from a
    --str-regions BED without a fourth column. A locus holding no STR region raises
    nothing, so the error comes at the first locus that holds one (decided by the owner
    on 19 October 2026).
    """
    region_types: set[RegionType] = set()
    for str_region in locus.str_regions:
        if str_region.region_type is None:
            contig, start, end = str_region.contig, str_region.start, str_region.end
            raise PvcEvalInputError(
                f"The STR region {contig}:{start + 1}-{end}, 1-based and inclusive, which the"
                f" --str-regions BED gives as `{contig} {start} {end}`, has no region type,"
                " since that line has no fourth column, and --strata pvc needs the type of every STR region. Give --str-regions a"
                " BED written by `pvc_eval str-regions`, which has the types, or leave"
                " --str-regions out so that the STR finder runs, and run again."
            )
        region_types.add(str_region.region_type)
    for region_type in _PVC_STRATA:
        if region_type in region_types:
            return region_type.value
    return RegionType.GENERIC.value


class Strata:
    """The intervals of a --strata BED, as the regions of each class. Positions are 0-based
    and half-open.

    regions_by_class: the regions each class covers, by class, in any order.
    """

    __slots__ = ("_regions_by_class",)

    def __init__(self, regions_by_class: Mapping[str, Regions]) -> None:
        # Sorted by code point once here, so that stratum gives the classes in that order.
        self._regions_by_class = sorted(regions_by_class.items())

    def stratum(self, locus: Locus) -> str:
        """The classes of the intervals that share at least one position with the locus,
        each once, sorted by code point and joined with `,`; `.` when there are none."""
        locus_classes = [
            name
            for name, regions in self._regions_by_class
            if regions.overlaps(locus.contig, locus.start, locus.end)
        ]
        return ",".join(locus_classes) if locus_classes else "."


def read_strata(path: Path, contigs: Sequence[str]) -> Strata:
    """The intervals of a --strata BED, read whole, as read_regions reads a --regions BED;
    its fourth column is the class of each interval.

    contigs: the reference's contig names.

    Raises PvcEvalInputError as parse_bed does, a line on a contig the reference lacks
    among them, and for a line whose class is missing, empty or holds a `,`; for a gzip
    BED that is damaged, it says so in place of the fault of the line (doc/spec.md §5.2).
    """
    spans_by_class: defaultdict[str, defaultdict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for bed_line in parse_bed(path, contigs):
        try:
            name = _parse_class(bed_line, path)
        except PvcEvalInputError as err:
            check_gzip_whole(err, path, "BED", bed_line.line_number)
            raise
        spans_by_class[name][bed_line.contig].append((bed_line.start, bed_line.end))
    return Strata({name: Regions(spans) for name, spans in spans_by_class.items()})


def _parse_class(bed_line: BedLine, path: Path) -> str:
    """The class of a line of a --strata BED, its fourth column without the spaces at its
    ends."""
    where = f"{path}, line {bed_line.line_number}:"
    name = "" if bed_line.name is None else bed_line.name.strip(" ")
    if not name:
        raise PvcEvalInputError(
            f"{where} has no class: its fourth column, which names the class of its interval,"
            " is missing or empty. Give each line of the --strata BED its class and run again."
        )
    if "," in name:
        raise PvcEvalInputError(
            f"{where} the class {name!r} holds a comma, which the table puts between the"
            " classes of a locus. Replace the commas in the fourth column of the --strata BED,"
            " for instance with `_`, and run again."
        )
    return name
