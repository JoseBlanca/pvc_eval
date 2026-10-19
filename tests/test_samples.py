import logging
from collections.abc import Sequence
from pathlib import Path

import pytest

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.samples import SamplePair, make_sample_pairs

_FIXED_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"


def _write_header(path: Path, sample_names: Sequence[str]) -> Path:
    """A VCF file of a header and no records, with these sample columns."""
    line = "\t".join([_FIXED_COLUMNS, *sample_names])
    path.write_text(f"##fileformat=VCFv4.3\n{line}\n", encoding="utf-8")
    return path


def _write_files(
    tmp_path: Path, samples_a: Sequence[str], samples_b: Sequence[str]
) -> tuple[Path, Path]:
    return (
        _write_header(tmp_path / "a.vcf", samples_a),
        _write_header(tmp_path / "b.vcf", samples_b),
    )


def _write_map(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "pairs.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def _get_warnings(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [
        record.getMessage()
        for record in caplog.records
        if record.name == "pvc_eval.samples" and record.levelno == logging.WARNING
    ]


def test_without_map_pairs_shared_names_in_order_of_file_a(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG004", "HG002", "HG003"], ["HG002", "HG004"])

    pairs = make_sample_pairs(path_a, path_b)

    assert pairs == (
        SamplePair(name_a="HG004", index_a=0, name_b="HG004", index_b=1),
        SamplePair(name_a="HG002", index_a=1, name_b="HG002", index_b=0),
    )


def test_map_gives_pairs_in_order_of_its_lines(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003"], ["NA1", "NA2"])
    sample_map = _write_map(tmp_path, "HG003\tNA1\nHG002\tNA2")  # no newline at the end

    pairs = make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert pairs == (
        SamplePair(name_a="HG003", index_a=1, name_b="NA1", index_b=0),
        SamplePair(name_a="HG002", index_a=0, name_b="NA2", index_b=1),
    )


def test_map_pairs_one_sample_of_a_with_two_replicates_of_b(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002_r1", "HG002_r2"])
    sample_map = _write_map(tmp_path, "HG002\tHG002_r1\nHG002\tHG002_r2\n")

    pairs = make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert pairs == (
        SamplePair(name_a="HG002", index_a=0, name_b="HG002_r1", index_b=0),
        SamplePair(name_a="HG002", index_a=0, name_b="HG002_r2", index_b=1),
    )


def test_map_skips_blank_lines(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003"], ["HG002", "HG003"])
    sample_map = _write_map(tmp_path, "\nHG002\tHG003\n  \n\t\n \t \n\nHG003\tHG002\n")

    pairs = make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert pairs == (
        SamplePair(name_a="HG002", index_a=0, name_b="HG003", index_b=1),
        SamplePair(name_a="HG003", index_a=1, name_b="HG002", index_b=0),
    )


@pytest.mark.parametrize(
    ("line", "how_many"),
    [
        ("HG002", "1 field"),
        ("HG002\tHG002\tHG003", "3 fields"),
        ("HG002\t", "an empty field"),
        ("\tHG002", "an empty field"),
    ],
    ids=["one_field", "three_fields", "empty_second_field", "empty_first_field"],
)
def test_map_line_without_two_fields_stops_the_run(
    tmp_path: Path, line: str, how_many: str
) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = _write_map(tmp_path, f"HG002\tHG002\n\n{line}\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 3: a line holds a sample of {path_a} and a sample of {path_b},"
        f" separated by a tab, and this one has {how_many}. Correct the line and run again."
    )


def test_map_sample_absent_from_file_b_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002", "HG003", "HG004"])
    sample_map = _write_map(tmp_path, "HG002\tHG002\nHG002\tHG003\nHG002\tHG005\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 3: 'HG005' is not a sample of {path_b}, whose samples are HG002,"
        " HG003, HG004. Correct the name and run again."
    )


def test_map_sample_absent_from_file_a_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003"], ["HG002"])
    sample_map = _write_map(tmp_path, "HG005\tHG002\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 1: 'HG005' is not a sample of {path_a}, whose samples are HG002,"
        " HG003. Correct the name and run again."
    )


def test_map_line_with_neither_sample_in_its_file_names_the_sample_of_a(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = _write_map(tmp_path, "HG005\tHG006\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 1: 'HG005' is not a sample of {path_a}, whose samples are HG002."
        " Correct the name and run again."
    )


def test_map_name_with_trailing_space_is_quoted_in_the_error(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = _write_map(tmp_path, "HG002\tHG002 \n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 1: 'HG002 ' is not a sample of {path_b}, whose samples are HG002."
        " Correct the name and run again."
    )


def test_map_ignores_a_byte_order_mark(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = tmp_path / "pairs.tsv"
    sample_map.write_bytes(b"\xef\xbb\xbfHG002\tHG002\n")

    pairs = make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert pairs == (SamplePair(name_a="HG002", index_a=0, name_b="HG002", index_b=0),)


def test_map_with_same_pair_on_two_lines_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003"], ["HG002", "HG003"])
    sample_map = _write_map(
        tmp_path, "HG003\tHG003\nHG002\tHG002\nHG002\tHG003\nHG003\tHG002\nHG002\tHG002\n"
    )

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, lines 2 and 5: both pair HG002 of {path_a} with HG002 of {path_b},"
        " which would give every row twice. Remove one and run again."
    )


def test_name_in_two_columns_of_file_a_stops_the_run_when_paired(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003", "HG002"], ["HG002"])

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b)

    assert str(error.value) == (
        f"{path_a} names the sample HG002 in two columns, 10 and 12. Rename one of them and"
        " run again."
    )


def test_name_in_three_columns_of_file_b_names_every_column(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002", "HG003", "HG002", "HG002"])
    sample_map = _write_map(tmp_path, "HG002\tHG002\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{path_b} names the sample HG002 in 3 columns, 10, 12 and 13. Rename one of them and"
        " run again."
    )


def test_name_in_two_columns_that_no_pair_needs_is_warned_of_once(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG003", "HG002", "HG003", "HG004"], ["HG002"])

    pairs = make_sample_pairs(path_a, path_b)

    assert pairs == (SamplePair(name_a="HG002", index_a=1, name_b="HG002", index_b=0),)
    assert _get_warnings(caplog) == [
        f"2 samples of {path_a} are in no sample pair and are not compared: HG003, HG004."
        " Pair them with --sample-map to compare them."
    ]


def test_map_sample_absent_from_file_names_its_first_ten_samples(tmp_path: Path) -> None:
    samples_b = [f"S{number:02}" for number in range(1, 13)]
    path_a, path_b = _write_files(tmp_path, ["HG002"], samples_b)
    sample_map = _write_map(tmp_path, "HG002\tHG005\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 1: 'HG005' is not a sample of {path_b}, whose samples are S01,"
        " S02, S03, S04, S05, S06, S07, S08, S09, S10 and 2 more. Correct the name and run"
        " again."
    )


def test_map_sample_of_file_with_no_samples_says_it_has_none(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], [])
    sample_map = _write_map(tmp_path, "HG002\tHG005\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map}, line 1: 'HG005' is not a sample of {path_b}, which has no samples."
        " Correct the name and run again."
    )


def test_map_that_does_not_exist_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = tmp_path / "pairs.tsv"

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == f"{sample_map} does not exist. Check its path and run again."


def test_map_that_is_a_directory_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = tmp_path / "pairs.tsv"
    sample_map.mkdir()

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map} is a directory, not a sample map. Check its path and run again."
    )


def test_map_that_cannot_be_opened_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = _write_map(tmp_path, "HG002\tHG002\n")
    sample_map.chmod(0o000)

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map} cannot be read: permission denied. Allow reading it and run again."
    )


@pytest.mark.parametrize("start", [b"", b"\xef\xbb\xbf"], ids=["plain", "byte_order_mark"])
def test_map_that_is_not_utf8_names_the_line(tmp_path: Path, start: bytes) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = tmp_path / "pairs.tsv"
    sample_map.write_bytes(start + b"HG002\tHG002\n\nHG\xe9002\tHG002\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == (
        f"{sample_map} is not UTF-8 text: line 3 cannot be read as UTF-8. Save it as UTF-8"
        " and run again."
    )


def test_files_that_share_no_name_stop_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["NA12878"])

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b)

    assert str(error.value) == (
        f"{path_a} and {path_b} share no sample name. Pair their samples with --sample-map"
        " and run again."
    )


def test_map_of_only_blank_lines_stops_the_run(tmp_path: Path) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002"], ["HG002"])
    sample_map = _write_map(tmp_path, "\n\n")

    with pytest.raises(PvcEvalInputError) as error:
        make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert str(error.value) == f"{sample_map} holds no pair."


def test_unpaired_samples_are_warned_of_once_per_file(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003", "HG004"], ["HG002", "NA1"])

    make_sample_pairs(path_a, path_b)

    assert _get_warnings(caplog) == [
        f"2 samples of {path_a} are in no sample pair and are not compared: HG003, HG004."
        " Pair them with --sample-map to compare them.",
        f"1 sample of {path_b} is in no sample pair and is not compared: NA1."
        " Pair it with --sample-map to compare it.",
    ]


def test_unpaired_samples_with_map_are_told_to_add_them_to_it(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003"], ["HG002"])
    sample_map = _write_map(tmp_path, "HG002\tHG002\n")

    make_sample_pairs(path_a, path_b, sample_map=sample_map)

    assert _get_warnings(caplog) == [
        f"1 sample of {path_a} is in no sample pair and is not compared: HG003."
        f" Add it to {sample_map} to compare it."
    ]


def test_warning_names_first_ten_samples_and_counts_the_rest(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    unpaired_names = [f"S{number:02}" for number in range(1, 16)]
    path_a, path_b = _write_files(tmp_path, ["HG002", *unpaired_names], ["HG002"])

    make_sample_pairs(path_a, path_b)

    assert _get_warnings(caplog) == [
        f"15 samples of {path_a} are in no sample pair and are not compared: S01, S02, S03,"
        " S04, S05, S06, S07, S08, S09, S10 and 5 more. Pair them with --sample-map to"
        " compare them."
    ]


@pytest.mark.parametrize(
    ("num_unpaired", "names_text"),
    [
        (10, "S01, S02, S03, S04, S05, S06, S07, S08, S09, S10"),
        (11, "S01, S02, S03, S04, S05, S06, S07, S08, S09, S10 and 1 more"),
    ],
    ids=["ten_named_no_count", "eleven_one_more"],
)
def test_warning_counts_only_samples_past_the_tenth(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, num_unpaired: int, names_text: str
) -> None:
    unpaired_names = [f"S{number:02}" for number in range(1, num_unpaired + 1)]
    path_a, path_b = _write_files(tmp_path, ["HG002", *unpaired_names], ["HG002"])

    make_sample_pairs(path_a, path_b)

    assert _get_warnings(caplog) == [
        f"{num_unpaired} samples of {path_a} are in no sample pair and are not compared:"
        f" {names_text}. Pair them with --sample-map to compare them."
    ]


def test_every_sample_paired_gives_no_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path_a, path_b = _write_files(tmp_path, ["HG002", "HG003"], ["HG003", "HG002"])

    make_sample_pairs(path_a, path_b)

    assert _get_warnings(caplog) == []
