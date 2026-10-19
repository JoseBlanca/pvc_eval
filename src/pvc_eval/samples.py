"""The sample pairs of a comparison (doc/specs/reconciliation.md §3): which sample of file A
is compared with which sample of file B, by name or by a --sample-map file. Only the headers
of the two files are read.

Public: SamplePair and make_sample_pairs.
"""

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from pvc_eval.errors import PvcEvalInputError, format_count
from pvc_eval.vcf import parse_vcf_samples

_logger = logging.getLogger(__name__)

# The column of a VCF line, counted from 1, that holds the first sample: after the nine of
# CHROM to FORMAT (VCF 4.3 §1.5), which vcf.py counts as _NUM_FIXED_COLUMNS.
_FIRST_SAMPLE_COLUMN = 10
# How many samples a message names before it counts the rest (§3).
_NAMES_IN_MESSAGE = 10
_BYTE_ORDER_MARK = "\ufeff"  # U+FEFF, which UTF-8 writes as EF BB BF


@dataclass(frozen=True, slots=True, kw_only=True)
class SamplePair:
    """A sample of file A and the sample of file B it is compared with."""

    name_a: str
    index_a: int  # the sample's place among file A's sample columns, from 0
    name_b: str
    index_b: int  # the sample's place among file B's sample columns, from 0


def make_sample_pairs(
    path_a: Path, path_b: Path, *, sample_map: Path | None = None
) -> tuple[SamplePair, ...]:
    """The sample pairs of files A and B, read from their headers (§3).

    Without sample_map, every name in both headers is a pair, in file A's header order.
    With it, each line of the map, a sample of A and a sample of B separated by a tab, is a
    pair, in the order of the lines; blank lines are skipped.

    The samples of each file that are in no pair are logged in one warning per file.

    Raises PvcEvalInputError for a map that cannot be read or is not UTF-8, a line of the
    map that is not two non-empty fields, a
    sample of the map absent from its file, the same pair on two lines of the map, a name
    in two columns of a header when a pair needs it, and no pair at all.
    """
    samples_a = parse_vcf_samples(path_a)
    samples_b = parse_vcf_samples(path_b)
    columns_a = _get_columns_by_name(samples_a)
    columns_b = _get_columns_by_name(samples_b)
    if sample_map is None:
        sample_pairs = _make_pairs_by_name(path_a, columns_a, path_b, columns_b)
    else:
        sample_pairs = _make_pairs_from_map(sample_map, path_a, columns_a, path_b, columns_b)
    for path, sample_names, paired_indexes in (
        (path_a, samples_a, {pair.index_a for pair in sample_pairs}),
        (path_b, samples_b, {pair.index_b for pair in sample_pairs}),
    ):
        warning = _make_unpaired_warning(path, sample_names, paired_indexes, sample_map)
        if warning is not None:
            _logger.warning("%s", warning)
    return sample_pairs


def _make_pairs_by_name(
    path_a: Path,
    columns_a: Mapping[str, Sequence[int]],
    path_b: Path,
    columns_b: Mapping[str, Sequence[int]],
) -> tuple[SamplePair, ...]:
    sample_pairs = tuple(
        SamplePair(
            name_a=name,
            index_a=_find_column(path_a, columns_a, name),
            name_b=name,
            index_b=_find_column(path_b, columns_b, name),
        )
        for name in columns_a
        if name in columns_b
    )
    if not sample_pairs:
        raise PvcEvalInputError(
            f"{path_a} and {path_b} share no sample name. Pair their samples with"
            " --sample-map and run again."
        )
    return sample_pairs


def _make_pairs_from_map(
    sample_map: Path,
    path_a: Path,
    columns_a: Mapping[str, Sequence[int]],
    path_b: Path,
    columns_b: Mapping[str, Sequence[int]],
) -> tuple[SamplePair, ...]:
    # The line each pair was first given on, to name both lines of a pair given twice.
    line_of_pair: dict[SamplePair, int] = {}
    for line_number, line in enumerate(_read_map_lines(sample_map), start=1):
        if not line.strip():
            continue
        fields = line.rstrip("\r").split("\t")
        if len(fields) != 2 or "" in fields:
            raise _make_fields_error(sample_map, line_number, fields, path_a, path_b)
        name_a, name_b = fields
        if name_a not in columns_a:
            raise _make_absent_error(sample_map, line_number, name_a, path_a, columns_a)
        if name_b not in columns_b:
            raise _make_absent_error(sample_map, line_number, name_b, path_b, columns_b)
        sample_pair = SamplePair(
            name_a=name_a,
            index_a=_find_column(path_a, columns_a, name_a),
            name_b=name_b,
            index_b=_find_column(path_b, columns_b, name_b),
        )
        if sample_pair in line_of_pair:
            raise PvcEvalInputError(
                f"{sample_map}, lines {line_of_pair[sample_pair]} and {line_number}: both"
                f" pair {name_a} of {path_a} with {name_b} of {path_b}, which would give"
                " every row twice. Remove one and run again."
            )
        line_of_pair[sample_pair] = line_number
    if not line_of_pair:
        raise PvcEvalInputError(
            f"{sample_map} holds no pair. Add a line for each pair, a sample of {path_a} and a"
            f" sample of {path_b} separated by a tab, and run again."
        )
    return tuple(line_of_pair)


def _read_map_lines(sample_map: Path) -> list[str]:
    """The lines of a sample map, without their newlines.

    Raises PvcEvalInputError for a map that is missing, a directory, cannot be opened or
    is not UTF-8 text.
    """
    try:
        # Whole, since a map holds a line per pair, and so the bytes that are not UTF-8
        # can be found on their line.
        map_bytes = sample_map.read_bytes()
    except FileNotFoundError:
        raise PvcEvalInputError(
            f"{sample_map} does not exist. Check its path and run again."
        ) from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{sample_map} is a directory, not a sample map. Check its path and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{sample_map} cannot be read: permission denied. Allow reading it and run again."
        ) from None
    try:
        # The byte-order mark that spreadsheet programs write is dropped (§3) here and not
        # by decoding with utf-8-sig, whose error offsets would not count it.
        map_text = map_bytes.decode("utf-8").removeprefix(_BYTE_ORDER_MARK)
    except UnicodeDecodeError as err:
        line_number = map_bytes.count(b"\n", 0, err.start) + 1
        raise PvcEvalInputError(
            f"{sample_map} is not UTF-8 text: line {line_number} cannot be read as UTF-8."
            " Save it as UTF-8 and run again."
        ) from err
    return map_text.split("\n")


def _get_columns_by_name(sample_names: Sequence[str]) -> dict[str, list[int]]:
    """The places of each name among the sample columns, from 0, in the header's order."""
    columns_by_name: dict[str, list[int]] = {}
    for index, name in enumerate(sample_names):
        columns_by_name.setdefault(name, []).append(index)
    return columns_by_name


def _find_column(path: Path, columns_by_name: Mapping[str, Sequence[int]], name: str) -> int:
    """The place of a sample among the sample columns, from 0.

    Raises PvcEvalInputError when the header names it in more than one column, since its
    genotypes could then be either.
    """
    columns = columns_by_name[name]
    if len(columns) == 1:
        return columns[0]
    vcf_columns = [str(index + _FIRST_SAMPLE_COLUMN) for index in columns]
    how_many = "two" if len(columns) == 2 else str(len(columns))
    raise PvcEvalInputError(
        f"{path} names the sample {name} in {how_many} columns,"
        f" {', '.join(vcf_columns[:-1])} and {vcf_columns[-1]}. Rename one of them and"
        " run again."
    )


def _make_absent_error(
    sample_map: Path,
    line_number: int,
    name: str,
    path: Path,
    columns_by_name: Mapping[str, Sequence[int]],
) -> PvcEvalInputError:
    if columns_by_name:
        its_samples = f"whose samples are {_format_names(list(columns_by_name))}"
    else:
        its_samples = "which has no samples"
    return PvcEvalInputError(
        f"{sample_map}, line {line_number}: '{name}' is not a sample of {path}, {its_samples}."
        " Correct the name and run again."
    )


def _make_fields_error(
    sample_map: Path, line_number: int, fields: Sequence[str], path_a: Path, path_b: Path
) -> PvcEvalInputError:
    # Two fields are refused only when one is empty.
    how_many = "an empty field" if len(fields) == 2 else format_count(len(fields), "field")
    return PvcEvalInputError(
        f"{sample_map}, line {line_number}: a line holds a sample of {path_a} and a sample"
        f" of {path_b}, separated by a tab, and this one has {how_many}. Correct the line and"
        " run again."
    )


def _make_unpaired_warning(
    path: Path, sample_names: Sequence[str], paired_indexes: set[int], sample_map: Path | None
) -> str | None:
    """The warning that names the samples of a file in no pair, each name once, or None
    when every sample is in one."""
    # A dict, to name once a name that fills two columns, in the header's order.
    unpaired_names = list(
        dict.fromkeys(
            name for index, name in enumerate(sample_names) if index not in paired_indexes
        )
    )
    if not unpaired_names:
        return None
    if len(unpaired_names) == 1:
        verb, them = "is", "it"
    else:
        verb, them = "are", "them"
    if sample_map is None:
        advice = f"Pair {them} with --sample-map to compare {them}."
    else:
        advice = f"Add {them} to {sample_map} to compare {them}."
    return (
        f"{format_count(len(unpaired_names), 'sample')} of {path} {verb} in no sample pair and {verb} not compared:"
        f" {_format_names(unpaired_names)}. {advice}"
    )


def _format_names(names: Sequence[str]) -> str:
    """The first 10 names, separated by commas, and how many more there are (§3)."""
    names_text = ", ".join(names[:_NAMES_IN_MESSAGE])
    if len(names) > _NAMES_IN_MESSAGE:
        names_text += f" and {len(names) - _NAMES_IN_MESSAGE} more"
    return names_text
