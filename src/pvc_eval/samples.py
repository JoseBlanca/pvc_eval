"""The sample pairs of a comparison (doc/specs/reconciliation.md §3): which sample of file A
is compared with which sample of file B, by name or by a --sample-map file. Only the headers
of the two files are read.

Public: SamplePair and make_sample_pairs.
"""

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.vcf import parse_vcf_samples

_logger = logging.getLogger(__name__)

# The column of a VCF line, counted from 1, that holds the first sample: after the nine of
# CHROM to FORMAT (VCF 4.3 §1.5).
_FIRST_SAMPLE_COLUMN = 10
# How many samples a warning of samples in no pair names before it counts the rest (§3).
_NAMES_IN_WARNING = 10


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

    Raises PvcEvalInputError for a line of the map that is not two non-empty fields, a
    sample of the map absent from its file, the same pair on two lines of the map, a name
    in two columns of a header when a pair needs it, and no pair at all.
    """
    samples_a = parse_vcf_samples(path_a)
    samples_b = parse_vcf_samples(path_b)
    if sample_map is None:
        sample_pairs = _make_pairs_by_name(path_a, samples_a, path_b, samples_b)
    else:
        sample_pairs = _make_pairs_from_map(sample_map, path_a, samples_a, path_b, samples_b)
    for path, sample_names, paired_indexes in (
        (path_a, samples_a, {pair.index_a for pair in sample_pairs}),
        (path_b, samples_b, {pair.index_b for pair in sample_pairs}),
    ):
        warning = _make_unpaired_warning(path, sample_names, paired_indexes, sample_map)
        if warning is not None:
            _logger.warning("%s", warning)
    return sample_pairs


def _make_pairs_by_name(
    path_a: Path, samples_a: Sequence[str], path_b: Path, samples_b: Sequence[str]
) -> tuple[SamplePair, ...]:
    columns_a = _get_columns_by_name(samples_a)
    columns_b = _get_columns_by_name(samples_b)
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
    samples_a: Sequence[str],
    path_b: Path,
    samples_b: Sequence[str],
) -> tuple[SamplePair, ...]:
    columns_a = _get_columns_by_name(samples_a)
    columns_b = _get_columns_by_name(samples_b)
    # The line each pair was first given on, to name both lines of a pair given twice.
    line_of_pair: dict[SamplePair, int] = {}
    with sample_map.open(encoding="utf-8") as map_file:
        for line_number, line in enumerate(map_file, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or "" in fields:
                raise _make_fields_error(sample_map, line_number, fields, path_a, path_b)
            name_a, name_b = fields
            _check_in_file(sample_map, line_number, name_a, path_a, samples_a)
            _check_in_file(sample_map, line_number, name_b, path_b, samples_b)
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
        raise PvcEvalInputError(f"{sample_map} holds no pair.")
    return tuple(line_of_pair)


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


def _check_in_file(
    sample_map: Path, line_number: int, name: str, path: Path, sample_names: Sequence[str]
) -> None:
    if name in sample_names:
        return
    raise PvcEvalInputError(
        f"{sample_map}, line {line_number}: {name} is not a sample of {path}, whose samples"
        f" are {', '.join(sample_names)}. Correct the name and run again."
    )


def _make_fields_error(
    sample_map: Path, line_number: int, fields: Sequence[str], path_a: Path, path_b: Path
) -> PvcEvalInputError:
    if len(fields) == 2:
        how_many = "an empty field"
    else:
        how_many = f"{len(fields)} field" if len(fields) == 1 else f"{len(fields)} fields"
    return PvcEvalInputError(
        f"{sample_map}, line {line_number}: a line holds a sample of {path_a} and a sample"
        f" of {path_b}, separated by a tab, and this one has {how_many}. Correct the line and"
        " run again."
    )


def _make_unpaired_warning(
    path: Path, sample_names: Sequence[str], paired_indexes: set[int], sample_map: Path | None
) -> str | None:
    """The warning that names the samples of a file in no pair, or None when every sample
    is in one. Samples are named by their columns, so a name in two columns that no pair
    needs is named twice."""
    unpaired_names = [
        name for index, name in enumerate(sample_names) if index not in paired_indexes
    ]
    if not unpaired_names:
        return None
    names_text = ", ".join(unpaired_names[:_NAMES_IN_WARNING])
    if len(unpaired_names) > _NAMES_IN_WARNING:
        names_text += f" and {len(unpaired_names) - _NAMES_IN_WARNING} more"
    if len(unpaired_names) == 1:
        samples, verb, them = "1 sample", "is", "it"
    else:
        samples, verb, them = f"{len(unpaired_names)} samples", "are", "them"
    if sample_map is None:
        advice = f"Pair {them} with --sample-map to compare {them}."
    else:
        advice = f"Add {them} to {sample_map} to compare {them}."
    return (
        f"{samples} of {path} {verb} in no sample pair and {verb} not compared:"
        f" {names_text}. {advice}"
    )
