"""The `pvc_eval` command line (doc/spec.md §3). It has one command so far, `str-regions`
(doc/specs/str_finder.md, "The str-regions command").

A fault of the input stops the run with its message on stderr and exit status 1, and no
traceback; anything else is a defect of pvc_eval and rises with its traceback.

Public: main, the `pvc_eval` entry point.
"""

import argparse
import re
import sys
import textwrap
from collections.abc import Callable, Sequence
from dataclasses import dataclass, fields
from pathlib import Path

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import StrFinderParameters, write_str_regions


@dataclass(frozen=True, slots=True, kw_only=True)
class _FinderOption:
    """An option of `str-regions` that sets a field of StrFinderParameters."""

    option: str
    field: str
    metavar: str
    parse: Callable[[str], object]
    help: str


def _parse_copy_floors(text: str) -> tuple[int, ...]:
    """The copy floors of --min-copies, whole numbers separated by commas."""
    try:
        return tuple(int(floor) for floor in text.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"'{text}' is not whole numbers separated by commas, such as 5,5,4,4,4,3"
        ) from None


# The options of doc/specs/str_finder.md, "The str-regions command", in its order.
_FINDER_OPTIONS = (
    _FinderOption(
        option="--scan-min-period",
        field="scan_min_period",
        metavar="N",
        parse=int,
        help="the shortest period, in bases, the scan for repeats looks for, from 1 to 6",
    ),
    _FinderOption(
        option="--scan-max-period",
        field="scan_max_period",
        metavar="N",
        parse=int,
        help="the longest period, in bases, the scan for repeats looks for, from 1 to 6",
    ),
    _FinderOption(
        option="--min-period",
        field="min_period",
        metavar="N",
        parse=int,
        help="the shortest period, in bases, of a repeat that may be an STR segment or a"
        " member of an STR bundle, from 1 to 6",
    ),
    _FinderOption(
        option="--max-period",
        field="max_period",
        metavar="N",
        parse=int,
        help="the longest period, in bases, of a repeat that may be an STR segment or a"
        " member of an STR bundle, from 1 to 6",
    ),
    _FinderOption(
        option="--min-copies",
        field="copy_floors",
        metavar="N,N,N,N,N,N",
        parse=_parse_copy_floors,
        help="the copy floors: the fewest copies of its motif a repeat needs, for periods 1"
        " to 6 in that order",
    ),
    _FinderOption(
        option="--min-purity",
        field="min_purity",
        metavar="X",
        parse=float,
        help="the smallest fraction, from 0 to 1, of the bases of an STR segment that equal"
        " its motif repeated",
    ),
    _FinderOption(
        option="--min-score",
        field="min_score",
        metavar="N",
        parse=int,
        help="the smallest score the scan may give a repeat",
    ),
    _FinderOption(
        option="--bundle-distance",
        field="bundle_distance",
        metavar="N",
        parse=int,
        help="repeats with fewer than N bases between them make an STR bundle; 1 or more",
    ),
    _FinderOption(
        option="--contig-end-flank",
        field="contig_end_flank",
        metavar="N",
        parse=int,
        help="a repeat with fewer than N bases between it and an end of its contig is"
        " dropped before anything else; 0 keeps every repeat",
    ),
    _FinderOption(
        option="--max-str-len",
        field="max_str_len",
        metavar="N",
        parse=int,
        help="a stretch of repeats longer than N bases is a satellite; at least --bundle-distance",
    ),
    _FinderOption(
        option="--max-str-bundle-len",
        field="max_bundle_len",
        metavar="N",
        parse=int,
        help="an STR bundle longer than N bases is a satellite",
    ),
    _FinderOption(
        option="--scan-match-reward",
        field="match_reward",
        metavar="N",
        parse=int,
        help="what the scan adds to a repeat's score for a base equal to the one a period"
        " before it; 1 or more",
    ),
    _FinderOption(
        option="--scan-mismatch-penalty",
        field="mismatch_penalty",
        metavar="N",
        parse=int,
        help="what the scan takes from a repeat's score for a base that differs from the one"
        " a period before it; 1 or more",
    ),
    _FinderOption(
        option="--scan-min-copies",
        field="scan_min_copies",
        metavar="N",
        parse=int,
        help="the fewest copies of its motif a repeat the scan finds needs to be kept; 1 or more",
    ),
)


# The description of `pvc_eval str-regions --help`, one paragraph per item, in the words of
# doc/specs/str_finder.md.
_STR_REGIONS_ABOUT = (
    "Find the STR regions of a reference, its STR segments, STR bundles and satellites, and"
    " write them as a BED, which `pvc_eval compare --str-regions` reads.",
    "A repeat is a stretch of the reference made of copies of a short sequence, its motif;"
    " the motif's length, 1 to 6 bases, is the repeat's period. An STR segment is a repeat"
    " with at least --bundle-distance bases between it and any other, trimmed to whole"
    " copies of its motif. An STR bundle is a stretch of repeats with fewer than"
    " --bundle-distance bases between each and the next. A satellite is a stretch of"
    " repeats longer than --max-str-len bases, or an STR bundle longer than"
    " --max-str-bundle-len bases.",
    "The scan looks for repeats of periods --scan-min-period to --scan-max-period, and only"
    " those of periods --min-period to --max-period become STR segments or members of STR"
    " bundles; the others still make satellites.",
    "The BED has a ## line for pvc_eval's version and one for each parameter, a #chrom"
    " line, and one row per STR region, in the order of the reference's .fai and then by"
    " start, with five columns: the contig; the start and the end, 0-based and half-open;"
    " the region type, str, str-bundle or satellite; and the motif of an STR segment, or ."
    " for the other two types.",
)


def main(argv: Sequence[str] | None = None) -> None:
    """Run `pvc_eval` with these arguments, or those of the command line. A fault of the
    input prints its message to stderr and exits with status 1."""
    arguments = _make_parser().parse_args(argv)
    try:
        _run_str_regions(arguments)
    except PvcEvalInputError as err:
        print(err, file=sys.stderr)
        sys.exit(1)


def _make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pvc_eval",
        description="Compare two VCF files over a reference genome, haplotype by haplotype.",
    )
    commands = parser.add_subparsers(dest="command", required=True, metavar="command")
    str_regions = commands.add_parser(
        "str-regions",
        help="find the STR regions of a reference and write them as a BED",
        description="\n\n".join(
            textwrap.fill(paragraph, break_on_hyphens=False) for paragraph in _STR_REGIONS_ABOUT
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    str_regions.add_argument(
        "--reference",
        required=True,
        type=Path,
        metavar="FASTA",
        help="the reference, an uncompressed FASTA with its .fai beside it (ref.fa.fai)."
        " Every contig of the .fai is segmented, in its order",
    )
    str_regions.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="BED",
        help="the BED to write. It is written as BED.tmp and renamed when complete",
    )
    str_regions.add_argument(
        "--regions",
        type=Path,
        metavar="BED",
        help="write only the STR regions that overlap these regions by at least one base,"
        " each whole. Only the contigs it names are segmented, each whole, so it changes"
        " which STR regions are written and never where one starts or ends",
    )
    defaults = StrFinderParameters()
    for finder_option in _FINDER_OPTIONS:
        default = getattr(defaults, finder_option.field)
        shown = _format_copy_floors(default) if isinstance(default, tuple) else str(default)
        str_regions.add_argument(
            finder_option.option,
            dest=finder_option.field,
            type=finder_option.parse,
            default=default,
            metavar=finder_option.metavar,
            help=f"{finder_option.help} (default {shown})",
        )
    return parser


def _run_str_regions(arguments: argparse.Namespace) -> None:
    """The `str-regions` command: segment the reference and write its STR regions."""
    parameters = _make_finder_parameters(arguments)
    with IndexedFasta(arguments.reference) as reference:
        write_str_regions(reference, arguments.output, parameters, arguments.regions)


def _make_finder_parameters(arguments: argparse.Namespace) -> StrFinderParameters:
    """The parameters the options give.

    Raises PvcEvalInputError, naming the option, for a value outside its limit.
    """
    values = {field.name: getattr(arguments, field.name) for field in fields(StrFinderParameters)}
    try:
        return StrFinderParameters(**values)
    except ValueError as err:
        copy_floors = arguments.copy_floors
        message = str(err).replace(repr(copy_floors), _format_copy_floors(copy_floors))
        for finder_option in _FINDER_OPTIONS:
            message = re.sub(rf"\b{finder_option.field}\b", finder_option.option, message)
        raise PvcEvalInputError(f"{message}. Correct the option, and run again.") from err


def _format_copy_floors(copy_floors: Sequence[int]) -> str:
    """Copy floors as --min-copies takes them, 5,5,4,4,4,3."""
    return ",".join(map(str, copy_floors))
