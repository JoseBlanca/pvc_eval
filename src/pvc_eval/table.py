"""The table and the summary of `pvc_eval compare` (doc/specs/compare.md, "The table", "The
summary" and "Writing the files").

Public: TABLE_COLUMNS, the column names of the table; SUMMARY_COLUMNS, those of the summary;
format_table_row, one row's fields as text; write_table, which writes the table and the
summary from the rows.
"""

import gzip
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.kind import Kind
from pvc_eval.outcome import LocusGenotype, Outcome, OutcomeRow
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.vcf import MISSING

# doc/spec.md §6, in its order.
TABLE_COLUMNS: tuple[str, ...] = (
    "chrom",
    "start",
    "end",
    "kind",
    "stratum",
    "str_regions",
    "sample_a",
    "sample_b",
    "outcome",
    "alleles",
    "genotype_a",
    "genotype_b",
    "possible_a",
    "possible_b",
    "quality_a",
    "quality_b",
    "gq_a",
    "gq_b",
    "records_a",
    "records_b",
    "invalid_reason_a",
    "invalid_reason_b",
)
_SUMMARY_KINDS = (Kind.SNP, Kind.INDEL, Kind.STR, Kind.COMPLEX)  # a column each
SUMMARY_COLUMNS: tuple[str, ...] = (
    "outcome",
    "stratum",
    *(kind.value for kind in _SUMMARY_KINDS),
    "total",
)

_NONE = "."  # a field that is None, doc/specs/compare.md "The table"
_GZIP_LEVEL = 6  # doc/specs/compare.md "Writing the files", decided by the owner on 20 October 2026
# The summary's lines come in the order of doc/spec.md §5.5's table, which Outcome keeps.
_OUTCOME_ORDER = {outcome: index for index, outcome in enumerate(Outcome)}

type _SummaryKey = tuple[Outcome, str, Kind]


def format_table_row(row: OutcomeRow) -> tuple[str, ...]:
    """One row's fields as the table writes them, in the order of TABLE_COLUMNS
    (doc/specs/compare.md, "The table").

    start becomes 1-based and end stays as it is, so the locus is 1-based and inclusive, as
    are the STR regions; a field that is None is `.`.
    """
    return (
        row.contig,
        str(row.start + 1),
        str(row.end),
        row.kind.value,
        _NONE if row.stratum is None else row.stratum,
        _format_str_regions(row.str_regions),
        row.sample_a,
        row.sample_b,
        row.outcome.value,
        _NONE if row.alleles is None else ",".join(allele or "-" for allele in row.alleles),
        _format_genotype(row.genotype_a),
        _format_genotype(row.genotype_b),
        _format_optional(row.possible_a),
        _format_optional(row.possible_b),
        _format_quality(row.quality_a),
        _format_quality(row.quality_b),
        _format_optional(row.gq_a),
        _format_optional(row.gq_b),
        str(row.records_a),
        str(row.records_b),
        _NONE if row.invalid_reason_a is None else row.invalid_reason_a.value,
        _NONE if row.invalid_reason_b is None else row.invalid_reason_b.value,
    )


def write_table(
    rows: Iterable[OutcomeRow], output: Path, *, summary: Path | None = None, all_rows: bool = False
) -> None:
    """Write the table to output and, when summary is given, the summary to it
    (doc/specs/compare.md, "The table", "The summary" and "Writing the files").

    The table has every row but, without all_rows, those reference in both; the summary
    counts every row. Each file is written to `<path>.tmp` and renamed to its path when
    every row is written; a path ending in `.gz` is gzipped. When rows raises, or Ctrl-C
    stops the run, the `.tmp` files are removed, the files at the paths are left as they
    were, and the exception goes on.

    Raises PvcEvalInputError, before any row is read, naming --output or --summary, for a
    path that is a directory or whose `.tmp` cannot be made, and for a summary that is the
    same file as output.
    """
    if summary is not None and summary.resolve() == output.resolve():
        raise PvcEvalInputError(
            f"--summary {summary} is the same file as --output {output}. Give the summary a"
            " path of its own, and run again."
        )
    targets = [("--output", output)]
    if summary is not None:
        targets.append(("--summary", summary))
    # Only the temporary files opened here are removed: a `.tmp` that could not be opened,
    # a directory among them, is not pvc_eval's.
    opened: list[tuple[Path, Path, TextIO]] = []
    try:
        for option, path in targets:
            temporary = Path(f"{path}.tmp")
            opened.append((path, temporary, _open_temporary_output(option, path, temporary)))
        files = [file for _, _, file in opened]
        _write_rows(rows, files[0], files[1] if len(files) == 2 else None, all_rows=all_rows)
        for file in files:
            file.close()
        for path, temporary, _ in opened:
            temporary.replace(path)
    except BaseException:
        for _, temporary, file in opened:
            file.close()
            temporary.unlink(missing_ok=True)
        raise


def _write_rows(
    rows: Iterable[OutcomeRow], table: TextIO, summary: TextIO | None, *, all_rows: bool
) -> None:
    """Write the table's lines to table and, when it is given, the summary's to summary."""
    table.write(_format_line(TABLE_COLUMNS))
    counts: Counter[_SummaryKey] = Counter()
    for row in rows:
        counts[row.outcome, _NONE if row.stratum is None else row.stratum, row.kind] += 1
        if all_rows or row.outcome is not Outcome.REFERENCE_IN_BOTH:
            table.write(_format_line(format_table_row(row)))
    if summary is not None:
        summary.writelines(_format_summary(counts))


def _open_temporary_output(option: str, path: Path, temporary: Path) -> TextIO:
    """Open the file that is written before it is renamed to path, gzipped when path ends
    in `.gz`.

    Raises PvcEvalInputError, naming option, when path is a directory or temporary cannot
    be made, with the messages of `str-regions` for --output.
    """
    if path.is_dir():
        raise PvcEvalInputError(
            f"{option} {path} is a directory. Give the path of a file to write, and run again."
        )
    cannot = f"{option} {path} cannot be written:"
    try:
        if path.name.endswith(".gz"):
            return gzip.open(temporary, "wt", compresslevel=_GZIP_LEVEL, encoding="utf-8")
        return temporary.open("w", encoding="utf-8")
    except FileNotFoundError:
        raise PvcEvalInputError(
            f"{cannot} its directory, {path.parent}, does not exist. Make it, or give another"
            f" {option}, and run again."
        ) from None
    except NotADirectoryError:
        raise PvcEvalInputError(
            f"{cannot} its path goes through a file as if it were a directory. Correct the"
            " path, and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{cannot} permission denied to write in {path.parent}. Allow writing there, or"
            f" give another {option}, and run again."
        ) from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{cannot} {temporary}, where it is written first, is a directory. Move it away,"
            " and run again."
        ) from None


def _format_line(fields: Iterable[str]) -> str:
    """Fields separated by tabs, with the line end."""
    return "\t".join(fields) + "\n"


def _format_summary(counts: Counter[_SummaryKey]) -> list[str]:
    """The lines of the summary, its column names first: one line per outcome and stratum
    with a row, by outcome in the order of doc/spec.md §5.5 and then by stratum, by code
    point, a count per kind and their total."""
    lines = [_format_line(SUMMARY_COLUMNS)]
    outcome_strata = sorted(
        {(outcome, stratum) for outcome, stratum, _ in counts},
        key=lambda outcome_stratum: (_OUTCOME_ORDER[outcome_stratum[0]], outcome_stratum[1]),
    )
    for outcome, stratum in outcome_strata:
        kind_counts = [counts[outcome, stratum, kind] for kind in _SUMMARY_KINDS]
        fields = [outcome.value, stratum, *map(str, kind_counts), str(sum(kind_counts))]
        lines.append(_format_line(fields))
    return lines


def _format_str_regions(str_regions: Sequence[GenomeSegment]) -> str:
    """The STR regions, 1-based and inclusive, as `start-end:type`, with `:motif` after an
    STR segment's and `start-end` alone for a region with no type, joined with `,`."""
    if not str_regions:
        return _NONE
    return ",".join(map(_format_str_region, str_regions))


def _format_str_region(segment: GenomeSegment) -> str:
    """One STR region as the table writes it, 1-based and inclusive."""
    text = f"{segment.start + 1}-{segment.end}"
    if segment.region_type is None:
        return text
    text += f":{segment.region_type.value}"
    if segment.region_type is RegionType.STR:
        text += f":{segment.motif}"
    return text


def _format_genotype(genotype: LocusGenotype | None) -> str:
    """The indices joined with `/`, MISSING as `.`; `.` for None."""
    if genotype is None:
        return _NONE
    return "/".join(_NONE if index is MISSING else str(index) for index in genotype)


def _format_quality(quality: float | None) -> str:
    """The shortest decimal that reads back as the same number, which is what repr gives,
    with no `.0` at the end of a whole number; `.` for None."""
    if quality is None:
        return _NONE
    return repr(quality).removesuffix(".0")


def _format_optional(number: int | None) -> str:
    """A whole number, or `.` for None."""
    return _NONE if number is None else str(number)
