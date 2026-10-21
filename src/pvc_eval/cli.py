"""The `pvc_eval` command line (doc/spec.md §3). It has two commands: `compare`
(doc/specs/compare.md, "The `compare` command") and `str-regions` (doc/specs/str_finder.md,
"The str-regions command").

A fault of the input stops the run with its message on stderr and exit status 1, and no
traceback; anything else is a defect of pvc_eval and rises with its traceback. pvc_eval's
log, its warnings and INFO lines, goes to stderr, the message alone on each line.

Public: main, the `pvc_eval` entry point.
"""

import argparse
import logging
import re
import sys
import textwrap
from collections.abc import Callable, Generator, Sequence
from contextlib import closing, contextmanager
from dataclasses import dataclass, fields
from pathlib import Path

from pvc_eval.compare import PVC_STRATA, CompareConfig, compare
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.loci import DEFAULT_CHAIN_DISTANCE
from pvc_eval.reconcile import PHASING_CAP
from pvc_eval.segments import (
    DEFAULT_NUM_THREADS,
    MAX_NUM_THREADS,
    StrFinderParameters,
    write_str_regions,
)
from pvc_eval.table import write_table
from pvc_eval.vcf import DEFAULT_PLOIDY, DEFAULT_SPAN_LIMIT


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


# A paragraph of a --help description that starts with it is a command, printed as it is.
_COMMAND_INDENT = "    "
# The description of `pvc_eval compare --help`, one paragraph per item, in the words of
# doc/spec.md and doc/specs/compare.md.
_COMPARE_ABOUT = (
    "Compare two VCF files, file A and file B, over a reference genome. For each sample"
    " pair, the sequence each haplotype carries over each locus is rebuilt from the VCF"
    " records of each file, so that two spellings of the same genotype agree. A locus is a"
    " stretch of the reference holding VCF records of either file, chained by overlap or by"
    " the chain distance, and widened to every STR region a VCF record touches.",
    "Each sample pair gets one outcome at each locus: same genotype; same alleles, other"
    " genotype (the same non-reference alleles on different numbers of haplotypes); only in"
    " A; only in B; different; reference in both; missing in A, missing in B or missing in"
    " both (a haplotype the genotype says nothing about); invalid in A, invalid in B or"
    " invalid in both (VCF records that claim something impossible for the sample); or too"
    " complex (a sample that needs more phase assignments than --phasing-cap). Each locus"
    " has a kind: SNP, indel, STR (it holds an STR region) or complex.",
    "The STR regions are the tandem repeats of the reference. They come from --str-regions"
    " or, without it, from the STR finder, run with its default parameters on the contigs"
    " that hold a VCF record and, with --regions, only on those --regions names. It takes"
    " from tens of seconds to several minutes for each chromosome of GRCh38. This command"
    " writes them once, for --str-regions on later runs:",
    f"{_COMMAND_INDENT}pvc_eval str-regions --reference ref.fa --output str.bed",
    "The table has a line of column names and one row per sample pair and compared locus,"
    " its fields separated by tabs, positions 1-based and inclusive. A path ending in .gz is"
    " gzipped. Each file is written as PATH.tmp and renamed when complete.",
    "Messages go to stderr: the warnings, a line when the STR finder starts a contig, and a"
    " line when the loci of a contig are done.",
)

# What `pvc_eval` prints when Ctrl-C stops it, and its exit status then
# (doc/specs/compare.md, "The `compare` command").
_STOPPED_MESSAGE = "Stopped; nothing was written."
_STOPPED_STATUS = 130


@dataclass(frozen=True, slots=True, kw_only=True)
class _Limit:
    """The values an option of `compare` takes, checked before anything is read."""

    option: str
    field: str
    is_allowed: Callable[[int | None], bool]
    rule: str


# doc/spec.md §3, in its order.
_COMPARE_LIMITS = (
    _Limit(
        option="--ploidy-a",
        field="ploidy_a",
        is_allowed=lambda value: value in {1, 2},
        rule="a ploidy is 1 or 2",
    ),
    _Limit(
        option="--ploidy-b",
        field="ploidy_b",
        is_allowed=lambda value: value in {1, 2},
        rule="a ploidy is 1 or 2",
    ),
    _Limit(
        option="--phasing-cap",
        field="phasing_cap",
        is_allowed=lambda value: value is not None and value >= 1,
        rule="a phasing cap is 1 or more",
    ),
    _Limit(
        option="--chain-distance",
        field="chain_distance",
        is_allowed=lambda value: value is not None and value >= 0,
        rule="a chain distance is 0 or more",
    ),
    _Limit(
        option="--span-limit",
        field="span_limit",
        is_allowed=lambda value: value is None or value >= 1,
        rule="a span limit is a whole number of 1 or more, or none",
    ),
)

# The word --strata takes instead of a BED, and --span-limit instead of a number.
_PVC_STRATA_WORD = "pvc"
_NO_SPAN_LIMIT_WORD = "none"


def main(argv: Sequence[str] | None = None) -> None:
    """Run `pvc_eval` with these arguments, or those of the command line. A fault of the
    input prints its message to stderr and exits with status 1; Ctrl-C prints that nothing
    was written and exits with status 130, the temporary files removed by their writers."""
    try:
        arguments = _make_parser().parse_args(argv)
        with _log_to_stderr():
            arguments.run(arguments)
    except PvcEvalInputError as err:
        print(err, file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print(_STOPPED_MESSAGE, file=sys.stderr)
        sys.exit(_STOPPED_STATUS)


@contextmanager
def _log_to_stderr() -> Generator[None]:
    """Send pvc_eval's log, INFO and above, to stderr, the message alone on each line, until
    the block ends. It is not passed on to the root logger, so that a handler there, as
    logging.basicConfig() sets up, neither prints a line twice nor gets the INFO lines."""
    logger = logging.getLogger("pvc_eval")
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(message)s"))
    level = logger.level
    propagate = logger.propagate
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    try:
        yield
    finally:
        logger.removeHandler(handler)
        logger.setLevel(level)
        logger.propagate = propagate


def _make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pvc_eval",
        description="Compare two VCF files over a reference genome, haplotype by haplotype.",
    )
    commands = parser.add_subparsers(dest="command", required=True, metavar="command")
    _add_compare_parser(commands)
    str_regions = commands.add_parser(
        "str-regions",
        help="find the STR regions of a reference and write them as a BED",
        description="\n\n".join(
            textwrap.fill(paragraph, break_on_hyphens=False) for paragraph in _STR_REGIONS_ABOUT
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    str_regions.set_defaults(run=_run_str_regions)
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
    str_regions.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_NUM_THREADS,
        metavar="N",
        help="how many contigs are segmented at once, each on a thread of its own, from 1 to"
        f" {MAX_NUM_THREADS}. The peak memory grows with N: a contig being segmented holds 4.6 to"
        " 5.3 MB per Mb of its length, 1.15 GB for GRCh38's chr1. With more than one thread,"
        " Ctrl-C stops the run once the contigs being segmented have finished, which can take"
        f" minutes (default {DEFAULT_NUM_THREADS})",
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


def _add_compare_parser(commands: argparse._SubParsersAction) -> None:
    """The `compare` command and its options, those of doc/spec.md §3 but --evaluator."""
    compare_parser = commands.add_parser(
        "compare",
        help="compare two VCF files over a reference and write the outcome of each locus",
        description="\n\n".join(
            paragraph
            if paragraph.startswith(_COMMAND_INDENT)
            else textwrap.fill(paragraph, break_on_hyphens=False)
            for paragraph in _COMPARE_ABOUT
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    compare_parser.set_defaults(run=_run_compare)
    add = compare_parser.add_argument
    add(
        "--reference",
        required=True,
        type=Path,
        metavar="FASTA",
        help="the reference, an uncompressed FASTA with its .fai beside it (ref.fa.fai). The"
        " contigs' names and order come from the .fai",
    )
    add("--a", required=True, type=Path, metavar="VCF", help="file A, plain or bgzipped")
    add("--b", required=True, type=Path, metavar="VCF", help="file B, plain or bgzipped")
    add(
        "--output",
        required=True,
        type=Path,
        metavar="TSV",
        help="the table to write, one row per sample pair and compared locus",
    )
    add(
        "--regions",
        type=Path,
        metavar="BED",
        help="compare only the loci that share at least one position with these regions. Such"
        " a locus is compared whole, its VCF records outside the regions included; every other"
        " locus is skipped, with no row and no count in the summary (default: every locus)",
    )
    add(
        "--strata",
        metavar="BED|pvc",
        help="a BED whose fourth column names each interval's class; a locus takes the"
        " classes of the intervals it overlaps, joined with ',', or '.' for none. The word"
        " pvc instead gives each locus the type of the STR regions it holds: satellite,"
        " str-bundle, str or generic. A file named pvc is given as ./pvc (default: no"
        " stratum, written '.')",
    )
    add(
        "--str-regions",
        type=Path,
        metavar="BED",
        help="the STR regions to use instead of running the STR finder, as `pvc_eval"
        " str-regions` writes them (default: the STR finder's, with its default parameters)",
    )
    add(
        "--sample-map",
        type=Path,
        metavar="TSV",
        help="two columns: a sample of file A and the sample of file B it is compared with"
        " (default: samples paired by name)",
    )
    for option, file in (("--ploidy-a", "file A"), ("--ploidy-b", "file B")):
        add(
            option,
            type=int,
            default=DEFAULT_PLOIDY,
            metavar="N",
            help=f"the ploidy {file} is read with, the haplotypes of each of its samples: 1 or"
            f" 2 (default {DEFAULT_PLOIDY})",
        )
    add(
        "--phasing-cap",
        type=int,
        default=PHASING_CAP,
        metavar="N",
        help="the most phase assignments pvc_eval tries for one sample at one locus; a"
        f" sample that needs more is too complex. 1 or more (default {PHASING_CAP})",
    )
    add(
        "--span-limit",
        type=_parse_span_limit,
        default=DEFAULT_SPAN_LIMIT,
        metavar="N|none",
        help="how far, in bases, a VCF record's span may reach beyond its REF before it is cut"
        " back to the REF: a whole number of 1 or more, or none for no limit (default"
        f" {DEFAULT_SPAN_LIMIT})",
    )
    add(
        "--chain-distance",
        type=int,
        default=DEFAULT_CHAIN_DISTANCE,
        metavar="N",
        help="VCF records whose spans have at most N reference positions between them are"
        f" compared in one locus: 0 or more (default {DEFAULT_CHAIN_DISTANCE})",
    )
    add(
        "--summary",
        type=Path,
        metavar="TSV",
        help="also write the counts of the rows over all sample pairs, one line per outcome"
        " and stratum and one column per kind (default: no summary)",
    )
    add(
        "--all-rows",
        action="store_true",
        help="also write to the table the rows whose outcome is reference in both; the"
        " summary counts them either way (default: not written)",
    )


def _parse_span_limit(text: str) -> int | None:
    """The value of --span-limit: a whole number, or None for the word none."""
    if text == _NO_SPAN_LIMIT_WORD:
        return None
    try:
        return int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"'{text}' is neither a whole number nor {_NO_SPAN_LIMIT_WORD}"
        ) from None


def _run_compare(arguments: argparse.Namespace) -> None:
    """The `compare` command: compare the two files and write the table and the summary."""
    config = _make_compare_config(arguments)
    # Closing the rows on the way out, whatever happens, stops compare()'s threads and closes
    # their files.
    rows = compare(arguments.a, arguments.b, arguments.reference, config)
    with closing(rows):
        write_table(rows, arguments.output, summary=arguments.summary, all_rows=arguments.all_rows)


def _make_compare_config(arguments: argparse.Namespace) -> CompareConfig:
    """The CompareConfig the options give.

    Raises PvcEvalInputError, naming the option, for a value outside its limit.
    """
    for limit in _COMPARE_LIMITS:
        value = getattr(arguments, limit.field)
        if not limit.is_allowed(value):
            raise PvcEvalInputError(
                f"{limit.option} {value}: {limit.rule}. Correct the option, and run again."
            )
    return CompareConfig(
        regions=arguments.regions,
        strata=_parse_strata(arguments.strata),
        str_regions=arguments.str_regions,
        sample_map=arguments.sample_map,
        ploidy_a=arguments.ploidy_a,
        ploidy_b=arguments.ploidy_b,
        phasing_cap=arguments.phasing_cap,
        span_limit=arguments.span_limit,
        chain_distance=arguments.chain_distance,
    )


def _parse_strata(text: str | None) -> Path | PVC_STRATA | None:
    """The value of CompareConfig.strata for the text of --strata: PVC_STRATA for the word
    pvc, and a path for anything else, ./pvc included. It is told on the text, since
    Path("./pvc") equals Path("pvc")."""
    if text is None:
        return None
    if text == _PVC_STRATA_WORD:
        return PVC_STRATA
    return Path(text)


def _run_str_regions(arguments: argparse.Namespace) -> None:
    """The `str-regions` command: segment the reference and write its STR regions."""
    parameters = _make_finder_parameters(arguments)
    # Checked here and not left to write_str_regions, so that a ValueError raised while the
    # contigs are segmented, a defect, keeps its traceback.
    if not 1 <= arguments.threads <= MAX_NUM_THREADS:
        raise PvcEvalInputError(
            f"--threads is {arguments.threads}; it must be a whole number from 1 to"
            f" {MAX_NUM_THREADS}. Correct the option, and run again."
        )
    with IndexedFasta(arguments.reference) as reference:
        write_str_regions(
            reference, arguments.output, parameters, arguments.regions, arguments.threads
        )


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
