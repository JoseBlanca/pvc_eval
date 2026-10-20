"""The table and the summary of `pvc_eval compare` (doc/specs/compare.md, "The table", "The
summary" and "Writing the files").

Public: TABLE_COLUMNS, the column names of the table; SUMMARY_COLUMNS, those of the summary;
format_table_row, one row's fields as text; write_table, which writes the table and the
summary from the rows.
"""

import os
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval._text_files import make_write_error, open_temporary_output
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
_SUMMARY_KINDS = tuple(Kind)  # a column each, in the order of doc/spec.md §5.6
SUMMARY_COLUMNS: tuple[str, ...] = (
    "outcome",
    "stratum",
    *(kind.value for kind in _SUMMARY_KINDS),
    "total",
)

# What the table writes for a field that is None, and the summary for a row with no
# stratum (doc/specs/compare.md, "The table").
_DOT = "."
# The summary's lines come in the order of doc/spec.md §5.5's table, which Outcome keeps.
_OUTCOME_ORDER = {outcome: index for index, outcome in enumerate(Outcome)}

type _SummaryKey = tuple[Outcome, str, Kind]


def format_table_row(row: OutcomeRow) -> tuple[str, ...]:
    """One row's fields as the table writes them, in the order of TABLE_COLUMNS
    (doc/specs/compare.md, "The table").

    The row's positions, those of the locus and of its STR regions, are 0-based and
    half-open, and are written 1-based and inclusive: start becomes start + 1, and end stays
    as it is. A field that is None is `.`.
    """
    return (
        row.contig,
        str(row.start + 1),
        str(row.end),
        row.kind.value,
        _format_stratum(row.stratum),
        _format_str_regions(row.str_regions),
        row.sample_a,
        row.sample_b,
        row.outcome.value,
        _DOT if row.alleles is None else ",".join(allele or "-" for allele in row.alleles),
        _format_genotype(row.genotype_a),
        _format_genotype(row.genotype_b),
        _format_whole_number(row.possible_a),
        _format_whole_number(row.possible_b),
        _format_quality(row.quality_a),
        _format_quality(row.quality_b),
        _format_whole_number(row.gq_a),
        _format_whole_number(row.gq_b),
        str(row.records_a),
        str(row.records_b),
        _DOT if row.invalid_reason_a is None else row.invalid_reason_a.value,
        _DOT if row.invalid_reason_b is None else row.invalid_reason_b.value,
    )


def write_table(
    rows: Iterable[OutcomeRow], output: Path, *, summary: Path | None = None, all_rows: bool = False
) -> None:
    """Write the table to output and, when summary is given, the summary to it
    (doc/specs/compare.md, "The table", "The summary" and "Writing the files").

    The table has every row but, without all_rows, those reference in both; the summary
    counts every row. Each file is written to `<path>.tmp` and renamed to its path when
    every row is written, the summary first; a path ending in `.gz` is gzipped. When rows
    raises, or Ctrl-C stops the run, the `.tmp` files are removed, the files at the paths
    are left as they were, and the exception goes on.

    Raises PvcEvalInputError, before any row is read, naming --output or --summary, for a
    path that is a directory or whose `.tmp` cannot be made, and for a summary that is the
    same file as output, however either is spelled, or one is the other's `.tmp`; and,
    naming the option, for an error of the system while a file is written or renamed, a
    full disk among them.
    """
    if summary is not None:
        _check_distinct_paths(output, summary)
    targets = [("--output", output)]
    if summary is not None:
        targets.append(("--summary", summary))
    # Only the temporary files opened here are removed: a `.tmp` that could not be opened,
    # a directory among them, is not pvc_eval's.
    outputs: list[_Output] = []
    try:
        for option, path in targets:
            temporary = Path(f"{path}.tmp")
            compress = path.name.endswith(".gz")
            file = open_temporary_output(path, temporary, option, compress=compress)
            outputs.append(_Output(option, path, temporary, file))
        if len(outputs) == 2:
            _check_distinct_files(*outputs)
        _write_rows(rows, outputs[0], outputs[1] if len(outputs) == 2 else None, all_rows=all_rows)
        for written in outputs:
            _close_output(written)
        # The summary first and the table last, so a failure between the two renames leaves
        # the earlier table, never a table newer than its summary (doc/specs/compare.md,
        # "Writing the files").
        for written in reversed(outputs):
            _rename_output(written)
    except BaseException:
        _remove_temporary_files(outputs)
        raise


def _check_distinct_paths(output: Path, summary: Path) -> None:
    """Raise PvcEvalInputError when summary and output, or one and the other's `.tmp`, are
    one file by their paths, or are two files that exist and are one.

    Checked before anything is opened, so that a refused run does not empty an earlier
    file.
    """
    pairs = [
        (output, summary),
        (output, Path(f"{summary}.tmp")),
        (Path(f"{output}.tmp"), summary),
    ]
    for first, second in pairs:
        if first.resolve() == second.resolve() or (
            first.exists() and second.exists() and first.samefile(second)
        ):
            raise _make_same_file_error(output, summary)


def _check_distinct_files(table: _Output, summary: _Output) -> None:
    """Raise PvcEvalInputError when the two temporary files, now open, are one file, or one
    is the other's path: names that differ only in letter case on a file system that
    ignores it, which _check_distinct_paths cannot see while neither file exists."""
    table_stat = os.fstat(table.file.fileno())
    summary_stat = os.fstat(summary.file.fileno())
    if (
        os.path.samestat(table_stat, summary_stat)
        or _is_file_at(summary.path, table_stat)
        or _is_file_at(table.path, summary_stat)
    ):
        raise _make_same_file_error(table.path, summary.path)


def _is_file_at(path: Path, file_stat: os.stat_result) -> bool:
    """Whether the file at path, when there is one, is the file of file_stat."""
    try:
        return os.path.samestat(path.stat(), file_stat)
    except FileNotFoundError:
        return False


def _make_same_file_error(output: Path, summary: Path) -> PvcEvalInputError:
    return PvcEvalInputError(
        f"--summary {summary} is the same file as --output {output}. Give the summary a path"
        " of its own, and run again."
    )


@dataclass(frozen=True, slots=True)
class _Output:
    """A file write_table writes: the option that named it, its path, the temporary file it
    is written to first, and that file, open."""

    option: str
    path: Path
    temporary: Path
    file: TextIO


def _write_rows(
    rows: Iterable[OutcomeRow], table: _Output, summary: _Output | None, *, all_rows: bool
) -> None:
    """Write the table's lines to table and, when it is given, the summary's to summary."""
    _write_text(table, _format_line(TABLE_COLUMNS))
    counts: Counter[_SummaryKey] = Counter()
    for row in rows:
        counts[row.outcome, _format_stratum(row.stratum), row.kind] += 1
        if all_rows or row.outcome is not Outcome.REFERENCE_IN_BOTH:
            _write_text(table, _format_line(format_table_row(row)))
    if summary is not None:
        _write_text(summary, "".join(_format_summary(counts)))


def _write_text(output: _Output, text: str) -> None:
    """Write text to output's temporary file.

    Raises PvcEvalInputError, naming output's option, for an error of the system, a full
    disk among them.
    """
    try:
        output.file.write(text)
    except OSError as err:
        raise make_write_error(output.option, output.path, err) from err


def _close_output(output: _Output) -> None:
    """Close output's temporary file, which writes what is left in its buffers.

    Raises PvcEvalInputError, naming output's option, for an error of the system.
    """
    try:
        output.file.close()
    except OSError as err:
        raise make_write_error(output.option, output.path, err) from err


def _rename_output(output: _Output) -> None:
    """Rename output's temporary file to its path.

    Raises PvcEvalInputError, naming output's option, for an error of the system.
    """
    try:
        output.temporary.replace(output.path)
    except OSError as err:
        raise make_write_error(output.option, output.path, err) from err


def _remove_temporary_files(outputs: Iterable[_Output]) -> None:
    """Close and remove the temporary file of each output, after an error or Ctrl-C."""
    for output in outputs:
        # The exception that stopped the run goes on. Closing a file that could not be
        # written raises again, as it writes its buffers; that adds nothing, and must not
        # keep the next temporary file from being removed.
        with suppress(OSError):
            output.file.close()
        output.temporary.unlink(missing_ok=True)


def _format_line(fields: Iterable[str]) -> str:
    """Fields separated by tabs, with the line end."""
    return "\t".join(fields) + "\n"


def _format_summary(counts: Mapping[_SummaryKey, int]) -> list[str]:
    """The lines of the summary, its column names first: one line per outcome and stratum
    with a row, by outcome in the order of doc/spec.md §5.5 and then by stratum, by code
    point, a count per kind and their total."""
    lines = [_format_line(SUMMARY_COLUMNS)]
    outcome_strata = sorted(
        {(outcome, stratum) for outcome, stratum, _ in counts},
        key=lambda outcome_stratum: (_OUTCOME_ORDER[outcome_stratum[0]], outcome_stratum[1]),
    )
    for outcome, stratum in outcome_strata:
        # A kind with no row is written 0 (doc/specs/compare.md, "The summary").
        kind_counts = [counts.get((outcome, stratum, kind), 0) for kind in _SUMMARY_KINDS]
        fields = [outcome.value, stratum, *map(str, kind_counts), str(sum(kind_counts))]
        lines.append(_format_line(fields))
    return lines


def _format_str_regions(str_regions: Sequence[GenomeSegment]) -> str:
    """The STR regions, 1-based and inclusive, as `start-end:type`, with `:motif` after an
    STR segment's and `start-end` alone for a region with no type, joined with `,`."""
    if not str_regions:
        return _DOT
    return ",".join(map(_format_str_region_1based, str_regions))


def _format_str_region_1based(segment: GenomeSegment) -> str:
    """One STR region, 0-based and half-open, as the table writes it, 1-based and inclusive."""
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
        return _DOT
    return "/".join(_DOT if index is MISSING else str(index) for index in genotype)


def _format_quality(quality: float | None) -> str:
    """The shortest decimal that reads back as the same number, which is what repr gives,
    with no `.0` at the end of a whole number; `.` for None."""
    if quality is None:
        return _DOT
    return repr(quality).removesuffix(".0")


def _format_stratum(stratum: str | None) -> str:
    """A stratum as it is, `.` for None."""
    return _DOT if stratum is None else stratum


def _format_whole_number(number: int | None) -> str:
    """A whole number, or `.` for None."""
    return _DOT if number is None else str(number)
