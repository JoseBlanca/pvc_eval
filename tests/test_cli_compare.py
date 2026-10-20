"""The `pvc_eval compare` command (doc/specs/compare.md, "The `compare` command"), on the
four-contig case of tests/compare_case.py: the table and the summary it writes, its
options, the values it refuses, its messages on stderr and its exit statuses."""

import argparse
import inspect
import logging
import re
import subprocess
import sys
import sysconfig
from collections.abc import Generator, Sequence
from pathlib import Path

import pytest
from tests.compare_case import LINES_A, CompareCase, write_bed, write_compare_case

from pvc_eval import cli
from pvc_eval.cli import main
from pvc_eval.compare import PVC_STRATA, CompareConfig, compare
from pvc_eval.outcome import OutcomeRow
from pvc_eval.table import write_table

_TIME = r"(\d+ h )?(\d+ min )?\d+ s since the start\."


def _run(argv: Sequence[str], capsys: pytest.CaptureFixture[str]) -> tuple[int, str, str]:
    """The exit status of `pvc_eval` with these arguments, and what it wrote to stdout and
    to stderr."""
    try:
        main(argv)
    except SystemExit as exit_:
        status = exit_.code
    else:
        status = 0
    assert isinstance(status, int)
    captured = capsys.readouterr()
    return status, captured.out, captured.err


def _make_argv(case: CompareCase, output: Path, *options: str) -> list[str]:
    return [
        "compare",
        "--reference",
        str(case.reference),
        "--a",
        str(case.a),
        "--b",
        str(case.b),
        "--output",
        str(output),
        *options,
    ]


def _write_expected(
    case: CompareCase,
    config: CompareConfig,
    output: Path,
    *,
    summary: Path | None = None,
    all_rows: bool = False,
) -> None:
    """What write_table writes for compare()'s rows."""
    write_table(
        compare(case.a, case.b, case.reference, config),
        output,
        summary=summary,
        all_rows=all_rows,
    )


def _read_column(path: Path, column: str) -> list[str]:
    header, *lines = path.read_text(encoding="utf-8").splitlines()
    index = header.split("\t").index(column)
    return [line.split("\t")[index] for line in lines]


# ---- The table and the summary -----------------------------------------------------------


def test_table_is_what_write_table_writes_for_compares_rows(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    case = write_compare_case(tmp_path)
    output = tmp_path / "loci.tsv"
    status, stdout, _ = _run(_make_argv(case, output), capsys)
    expected = tmp_path / "expected.tsv"
    _write_expected(case, CompareConfig(), expected)
    assert (status, stdout) == (0, "")
    assert output.read_text(encoding="utf-8") == expected.read_text(encoding="utf-8")
    # Seven rows, less chr1:150, reference in both.
    assert _read_column(output, "start") == ["50", "100", "120", "200", "101", "250"]
    assert not Path(f"{output}.tmp").exists()


def test_all_rows_and_summary_are_written(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    case = write_compare_case(tmp_path)
    output = tmp_path / "loci.tsv.gz"
    summary = tmp_path / "summary.tsv"
    status, _, _ = _run(_make_argv(case, output, "--summary", str(summary), "--all-rows"), capsys)
    expected = tmp_path / "expected.tsv.gz"
    expected_summary = tmp_path / "expected_summary.tsv"
    _write_expected(case, CompareConfig(), expected, summary=expected_summary, all_rows=True)
    assert status == 0
    assert output.read_bytes() == expected.read_bytes()
    assert summary.read_text(encoding="utf-8") == expected_summary.read_text(encoding="utf-8")
    # Every row counted, by outcome in the order of doc/spec.md §5.5, then SNP, indel, STR,
    # complex and the total.
    assert summary.read_text(encoding="utf-8").splitlines()[1:] == [
        "reference in both\t.\t1\t0\t0\t0\t1",
        "same genotype\t.\t0\t0\t1\t0\t1",
        "same alleles, other genotype\t.\t1\t0\t0\t0\t1",
        "only in A\t.\t1\t0\t1\t0\t2",
        "only in B\t.\t1\t0\t1\t0\t2",
    ]


@pytest.mark.parametrize(
    ("option", "config", "column", "values"),
    [
        (["--ploidy-a", "1"], CompareConfig(ploidy_a=1), None, None),
        (["--ploidy-b", "1"], CompareConfig(ploidy_b=1), None, None),
        (
            ["--chain-distance", "200", "--phasing-cap", "1"],
            CompareConfig(chain_distance=200, phasing_cap=1),
            None,
            None,
        ),
        (["--span-limit", "1"], CompareConfig(span_limit=1), None, None),
        (["--span-limit", "none"], CompareConfig(span_limit=None), None, None),
        (["--chain-distance", "0"], CompareConfig(chain_distance=0), None, None),
        # HG003 of file A paired with HG004 of file B, and HG002 left out: five rows of
        # seven, those at chr1:100 and chr2:120 being reference in both.
        (
            ["--sample-map", "pairs.tsv"],
            CompareConfig(sample_map=Path("pairs.tsv")),
            "sample_a",
            ["HG003"] * 5,
        ),
        # The BED's one STR region, at chr1 41-60, widens the locus of chr1:50 and replaces
        # the finder's, so the locus of chr3:105 is the SNP alone.
        (
            ["--str-regions", "str.bed"],
            CompareConfig(str_regions=Path("str.bed")),
            "start",
            ["41", "100", "120", "200", "105", "250"],
        ),
    ],
    ids=[
        "ploidy-a",
        "ploidy-b",
        "chain-distance-phasing-cap",
        "span-limit",
        "none",
        "zero",
        "sample-map",
        "str-regions",
    ],
)
def test_options_reach_compare(
    option: list[str],
    config: CompareConfig,
    column: str | None,
    values: list[str] | None,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # chr1:115 AC>A deletes a C of CCC at 116-118, so its span reaches two bases past its
    # REF, which --span-limit 1 cuts.
    lines_a = [*LINES_A[:2], "chr1 115 AC A 20 0/1:30 0/0:30", *LINES_A[2:]]
    case = write_compare_case(tmp_path, lines_a=lines_a)
    (tmp_path / "pairs.tsv").write_text("HG003\tHG004\n", encoding="utf-8")
    write_bed(tmp_path / "str.bed", ["chr1 40 60 str T"])
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "loci.tsv"
    status, _, _ = _run(_make_argv(case, output, *option), capsys)
    expected = tmp_path / "expected.tsv"
    _write_expected(case, config, expected)
    assert status == 0
    assert output.read_text(encoding="utf-8") == expected.read_text(encoding="utf-8")
    if column is not None:
        assert _read_column(output, column) == values


def test_ploidy_a_of_1_makes_file_as_samples_invalid(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    case = write_compare_case(tmp_path)
    output = tmp_path / "loci.tsv"
    status, _, _ = _run(_make_argv(case, output, "--ploidy-a", "1"), capsys)
    assert status == 0
    # File A's genotypes are diploid, so read as haploid each is a ploidy mismatch.
    assert set(_read_column(output, "invalid_reason_a")) == {"ploidy mismatch", "."}
    assert "ploidy mismatch" in _read_column(output, "invalid_reason_a")


def test_strata_pvc_is_the_word_and_dot_slash_pvc_a_bed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    case = write_compare_case(tmp_path)
    write_bed(tmp_path / "pvc", ["chr1 0 300 low_mappability"])
    monkeypatch.chdir(tmp_path)
    word = tmp_path / "word.tsv"
    bed = tmp_path / "bed.tsv"
    assert _run(_make_argv(case, word, "--strata", "pvc"), capsys)[0] == 0
    assert _run(_make_argv(case, bed, "--strata", "./pvc"), capsys)[0] == 0
    expected = tmp_path / "expected.tsv"
    _write_expected(case, CompareConfig(strata=PVC_STRATA), expected)
    assert word.read_text(encoding="utf-8") == expected.read_text(encoding="utf-8")
    assert _read_column(word, "stratum") == ["generic", "str", "str", "generic", "str", "generic"]
    assert _read_column(bed, "stratum") == ["low_mappability", "low_mappability", *["."] * 4]


# ---- Values refused ----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("option", "message"),
    [
        (["--ploidy-a", "3"], "--ploidy-a 3: a ploidy is 1 or 2."),
        (["--ploidy-a", "0"], "--ploidy-a 0: a ploidy is 1 or 2."),
        (["--ploidy-b", "0"], "--ploidy-b 0: a ploidy is 1 or 2."),
        (["--ploidy-b", "3"], "--ploidy-b 3: a ploidy is 1 or 2."),
        (["--phasing-cap", "0"], "--phasing-cap 0: a phasing cap is 1 or more."),
        (["--chain-distance", "-1"], "--chain-distance -1: a chain distance is 0 or more."),
        (
            ["--span-limit", "0"],
            "--span-limit 0: a span limit is a whole number of 1 or more, or none.",
        ),
    ],
    ids=[
        "ploidy-a-3",
        "ploidy-a-0",
        "ploidy-b-0",
        "ploidy-b-3",
        "phasing-cap",
        "chain-distance",
        "span-limit",
    ],
)
def test_value_out_of_range_stops_the_run_before_anything_is_read(
    option: list[str], message: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Not even the reference is read: it does not exist.
    case = CompareCase(reference=tmp_path / "ref.fa", a=tmp_path / "a.vcf", b=tmp_path / "b.vcf")
    output = tmp_path / "loci.tsv"
    status, _, stderr = _run(_make_argv(case, output, *option), capsys)
    assert (status, stderr) == (1, f"{message} Correct the option, and run again.\n")
    assert not output.exists()


@pytest.mark.parametrize(
    "option",
    [["--ploidy-a", "two"], ["--span-limit", "None"], ["--chain-distance", "1.5"]],
    ids=["ploidy-a", "span-limit", "chain-distance"],
)
def test_value_that_is_not_a_number_is_refused_by_argparse(
    option: list[str], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    case = write_compare_case(tmp_path)
    status, _, stderr = _run(_make_argv(case, tmp_path / "loci.tsv", *option), capsys)
    assert status == 2
    assert stderr.startswith("usage: pvc_eval compare")
    assert f"argument {option[0]}: " in stderr


# ---- Errors and Ctrl-C -------------------------------------------------------------------


def test_input_error_from_a_stage_gives_its_message_and_no_table(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # chr1:30 after chr1:50: file A is not sorted, found after the rows of chr1:50 began.
    lines_a = [*LINES_A[:3], "chr1 30 T C 50 0/1:40 1/1:35", *LINES_A[3:]]
    case = write_compare_case(tmp_path, lines_a=lines_a)
    output = tmp_path / "loci.tsv"
    output.write_text("earlier\n", encoding="utf-8")
    status, _, stderr = _run(_make_argv(case, output), capsys)
    assert status == 1
    assert stderr.splitlines()[-1] == (
        f"{case.a} is not sorted: chr1:30 (line 6) comes after chr1:150 (line 5). Sort it in"
        " the contig order of the reference's .fai and by position, and run again."
    )
    assert "Traceback" not in stderr
    assert output.read_text(encoding="utf-8") == "earlier\n"
    assert not Path(f"{output}.tmp").exists()


def _make_stand_in(rows: Sequence[OutcomeRow]) -> Generator[OutcomeRow]:
    """compare()'s rows, then Ctrl-C."""
    yield from rows
    raise KeyboardInterrupt


def test_ctrl_c_stops_with_130_and_writes_nothing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    case = write_compare_case(tmp_path)
    rows = list(compare(case.a, case.b, case.reference))
    monkeypatch.setattr(cli, "compare", lambda *_args, **_kwargs: _make_stand_in(rows))
    output = tmp_path / "loci.tsv"
    summary = tmp_path / "summary.tsv"
    status, _, stderr = _run(_make_argv(case, output, "--summary", str(summary)), capsys)
    assert (status, stderr) == (130, "Stopped; nothing was written.\n")
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "a.vcf",
        "b.vcf",
        "ref.fa",
        "ref.fa.fai",
    ]


def test_ctrl_c_while_the_options_are_parsed_stops_with_130(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    def parse_args(*_args: object, **_kwargs: object) -> argparse.Namespace:
        raise KeyboardInterrupt

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", parse_args)
    case = write_compare_case(tmp_path)
    status, _, stderr = _run(_make_argv(case, tmp_path / "loci.tsv"), capsys)
    assert (status, stderr) == (130, "Stopped; nothing was written.\n")


# A defect of pvc_eval, which the command does not catch: run in a process of its own, so
# that Python's own handler prints the traceback and sets the status.
_DEFECT_SCRIPT = """
import sys
from pvc_eval import cli

def compare(*args, **kwargs):
    raise RuntimeError("a defect of pvc_eval")
    yield

cli.compare = compare
cli.main(sys.argv[1:])
"""


def test_a_defect_rises_with_its_traceback(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path)
    output = tmp_path / "loci.tsv"
    completed = subprocess.run(
        [sys.executable, "-c", _DEFECT_SCRIPT, *_make_argv(case, output)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 1
    assert completed.stderr.startswith("Traceback (most recent call last):")
    assert completed.stderr.endswith("RuntimeError: a defect of pvc_eval\n")
    assert not output.exists()
    assert not Path(f"{output}.tmp").exists()


def test_rows_are_closed_when_the_output_is_refused(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    case = write_compare_case(tmp_path)
    streams: list[Generator[OutcomeRow]] = []

    def make_rows(*_args: object, **_kwargs: object) -> Generator[OutcomeRow]:
        stream = _make_stand_in([])
        streams.append(stream)
        return stream

    monkeypatch.setattr(cli, "compare", make_rows)
    status, _, stderr = _run(_make_argv(case, tmp_path), capsys)
    assert status == 1
    assert stderr.startswith(f"--output {tmp_path} is a directory")
    assert inspect.getgeneratorstate(streams[0]) == inspect.GEN_CLOSED


# ---- Messages ----------------------------------------------------------------------------


def test_warnings_and_info_lines_go_to_stderr(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    case = write_compare_case(tmp_path)
    regions = write_bed(tmp_path / "regions.bed", ["chr1 0 120", "chr3 0 300"])
    status, stdout, stderr = _run(
        _make_argv(case, tmp_path / "loci.tsv", "--regions", str(regions)), capsys
    )
    expected = [
        re.escape(
            f"1 sample of {case.a} is in no sample pair and is not compared: HG003. Pair it"
            " with --sample-map to compare it."
        ),
        re.escape(
            f"1 sample of {case.b} is in no sample pair and is not compared: HG004. Pair it"
            " with --sample-map to compare it."
        ),
        re.escape("chr1: finding STR regions in 300 bases."),
        # chr1's last locus is built only when a record of chr3 comes, after chr3's scan;
        # chr2, which --regions does not name, is not scanned.
        re.escape("chr3: finding STR regions in 300 bases."),
        "chr1: 2 loci compared, 1 skipped by --regions; " + _TIME,
        "chr2: 0 loci compared, 2 skipped by --regions; " + _TIME,
        "chr3: 2 loci compared, 0 skipped by --regions; " + _TIME,
    ]
    assert (status, stdout) == (0, "")
    lines = stderr.splitlines()
    assert len(lines) == len(expected)
    for line, pattern in zip(lines, expected, strict=True):
        assert re.fullmatch(pattern, line), line


class _RecordList(logging.Handler):
    """A handler that keeps the records it is given."""

    def __init__(self) -> None:
        super().__init__()
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def test_a_root_handler_gets_no_line_and_none_is_printed_twice(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # As logging.basicConfig() or pytest's log_cli would set it up.
    case = write_compare_case(tmp_path)
    root_handler = _RecordList()
    logging.getLogger().addHandler(root_handler)
    try:
        status, _, stderr = _run(_make_argv(case, tmp_path / "loci.tsv"), capsys)
    finally:
        logging.getLogger().removeHandler(root_handler)
    assert status == 0
    assert root_handler.records == []
    lines = stderr.splitlines()
    assert len(lines) == 8
    assert len(set(lines)) == len(lines)


def test_the_logger_is_set_back_after_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    case = write_compare_case(tmp_path)
    logger = logging.getLogger("pvc_eval")
    before = (logger.level, logger.propagate, list(logger.handlers))
    assert _run(_make_argv(case, tmp_path / "loci.tsv"), capsys)[0] == 0
    assert (logger.level, logger.propagate, list(logger.handlers)) == before


# ---- --help and the installed command ----------------------------------------------------


def test_help_names_every_option_with_its_default(capsys: pytest.CaptureFixture[str]) -> None:
    status, stdout, _ = _run(["compare", "--help"], capsys)
    assert status == 0
    help_text = " ".join(stdout.split())
    for option in (
        "--reference FASTA",
        "--a VCF",
        "--b VCF",
        "--output TSV",
        "--regions BED",
        "--strata BED|pvc",
        "--str-regions BED",
        "--sample-map TSV",
        "--summary TSV",
        "--all-rows",
    ):
        assert option in help_text
    for option_and_default in (
        "--ploidy-a N",
        "--ploidy-b N",
        "--phasing-cap N",
        "--span-limit N|none",
        "--chain-distance N",
    ):
        assert option_and_default in help_text
    for default in ("(default 2)", "(default 4096)", "(default 10000)", "(default 10)"):
        assert default in help_text
    assert "--evaluator" not in help_text


def test_the_installed_command_runs(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path)
    output = tmp_path / "loci.tsv"
    command = Path(sysconfig.get_path("scripts")) / "pvc_eval"
    completed = subprocess.run(
        [command, *_make_argv(case, output)], check=True, capture_output=True, text=True
    )
    assert completed.stdout == ""
    assert "chr3: finding STR regions in 300 bases.\n" in completed.stderr
    assert output.read_text(encoding="utf-8").splitlines()[-1].startswith("chr3\t250\t250\tSNP")
