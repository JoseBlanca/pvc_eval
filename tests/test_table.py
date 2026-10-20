"""The table and the summary of doc/specs/compare.md, on OutcomeRow values built in the
tests: each field as format_table_row writes it, the worked row of "The table" as its line,
and the files write_table writes, refuses and removes ("Writing the files")."""

import dataclasses
import errno
import gzip
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

import pytest

from pvc_eval._text_files import make_write_error
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.kind import Kind
from pvc_eval.outcome import InvalidReason, Outcome, OutcomeRow
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.table import SUMMARY_COLUMNS, TABLE_COLUMNS, format_table_row, write_table
from pvc_eval.vcf import MISSING

# The worked row of "The table": file A's `100 CAT C 0/1`, QUAL 29.5 and GQ 45, for HG002,
# in the locus chr20:100-110 that the AT repeat at 101-110 widens, with --strata pvc.
_WORKED_ROW = OutcomeRow(
    contig="chr20",
    start=99,
    end=110,
    kind=Kind.STR,
    stratum="str",
    str_regions=(GenomeSegment("chr20", 100, 110, RegionType.STR, "AT"),),
    sample_a="HG002",
    sample_b="HG002",
    outcome=Outcome.ONLY_IN_A,
    alleles=("CATATATATAT", "CATATATAT"),
    genotype_a=(0, 1),
    genotype_b=(0, 0),
    possible_a=1,
    possible_b=1,
    quality_a=29.5,
    quality_b=None,
    gq_a=45,
    gq_b=None,
    records_a=1,
    records_b=0,
    invalid_reason_a=None,
    invalid_reason_b=None,
)
_WORKED_LINE = (
    "chr20\t100\t110\tSTR\tstr\t101-110:str:AT\tHG002\tHG002\tonly in A"
    "\tCATATATATAT,CATATATAT\t0/1\t0/0\t1\t1\t29.5\t.\t45\t.\t1\t0\t.\t.\n"
)
_HEADER = (
    "chrom\tstart\tend\tkind\tstratum\tstr_regions\tsample_a\tsample_b\toutcome\talleles"
    "\tgenotype_a\tgenotype_b\tpossible_a\tpossible_b\tquality_a\tquality_b\tgq_a\tgq_b"
    "\trecords_a\trecords_b\tinvalid_reason_a\tinvalid_reason_b\n"
)
_SUMMARY_HEADER = "outcome\tstratum\tSNP\tindel\tSTR\tcomplex\ttotal\n"


def _make_row(**fields: object) -> OutcomeRow:
    """The worked row with these fields changed."""
    return dataclasses.replace(_WORKED_ROW, **fields)


def _format_fields(row: OutcomeRow) -> dict[str, str]:
    """format_table_row's text, by column name."""
    return dict(zip(TABLE_COLUMNS, format_table_row(row), strict=True))


def _make_segment(
    start: int, end: int, region_type: RegionType | None, motif: str | None = None
) -> GenomeSegment:
    """A segment of chr20, 0-based and half-open."""
    return GenomeSegment("chr20", start, end, region_type, motif)


def test_the_column_names_are_those_of_spec_section_6() -> None:
    assert "\t".join(TABLE_COLUMNS) + "\n" == _HEADER
    assert SUMMARY_COLUMNS == ("outcome", "stratum", "SNP", "indel", "STR", "complex", "total")


def test_the_worked_row_gives_its_line(tmp_path: Path) -> None:
    output = tmp_path / "loci.tsv"
    write_table([_WORKED_ROW], output)
    assert output.read_text(encoding="utf-8") == _HEADER + _WORKED_LINE


@pytest.mark.parametrize(
    ("fields", "expected"),
    [
        ({"start": 0, "end": 1}, {"start": "1", "end": "1"}),
        ({"start": 99, "end": 100}, {"start": "100", "end": "100"}),
        (
            {"stratum": None, "alleles": None, "genotype_a": None, "possible_a": None},
            {"stratum": ".", "alleles": ".", "genotype_a": ".", "possible_a": "."},
        ),
        ({"quality_a": None, "gq_a": None}, {"quality_a": ".", "gq_a": "."}),
        ({"stratum": "."}, {"stratum": "."}),
        ({"genotype_a": (0, MISSING)}, {"genotype_a": "0/."}),
        ({"genotype_a": (MISSING, MISSING)}, {"genotype_a": "./."}),
        ({"genotype_a": (1,), "genotype_b": (0,)}, {"genotype_a": "1", "genotype_b": "0"}),
        ({"genotype_a": (MISSING,)}, {"genotype_a": "."}),
        ({"genotype_b": (1, 2)}, {"genotype_b": "1/2"}),
        ({"alleles": ("AT", "")}, {"alleles": "AT,-"}),
        ({"alleles": ("", "A")}, {"alleles": "-,A"}),
        ({"possible_a": 0, "possible_b": 12}, {"possible_a": "0", "possible_b": "12"}),
        ({"records_a": 3, "records_b": 0}, {"records_a": "3", "records_b": "0"}),
        ({"gq_b": 0}, {"gq_b": "0"}),
    ],
    ids=[
        "first_base_of_contig",
        "one_base_locus",
        "none_as_dot",
        "no_quality_no_gq",
        "stratum_dot_from_bed",
        "missing_haplotype",
        "missing_both_haplotypes",
        "haploid",
        "haploid_missing",
        "two_alt_indices",
        "empty_locus_allele",
        "empty_reference_allele",
        "possible_counts",
        "record_counts",
        "gq_zero",
    ],
)
def test_format_table_row_writes_each_field(
    fields: dict[str, object], expected: dict[str, str]
) -> None:
    formatted = _format_fields(_make_row(**fields))
    assert {column: formatted[column] for column in expected} == expected


@pytest.mark.parametrize(
    ("str_regions", "expected"),
    [
        ((), "."),
        ((_make_segment(100, 110, RegionType.STR, "AT"),), "101-110:str:AT"),
        ((_make_segment(0, 1, RegionType.STR, "A"),), "1-1:str:A"),
        ((_make_segment(100, 140, RegionType.STR_BUNDLE),), "101-140:str-bundle"),
        ((_make_segment(100, 400, RegionType.SATELLITE),), "101-400:satellite"),
        ((_make_segment(100, 110, None),), "101-110"),
        (
            (
                _make_segment(100, 110, RegionType.STR, "AT"),
                _make_segment(110, 130, RegionType.STR, "CAG"),
                _make_segment(140, 150, None),
            ),
            "101-110:str:AT,111-130:str:CAG,141-150",
        ),
    ],
    ids=["none", "str", "str_at_first_base", "bundle", "satellite", "no_type", "three"],
)
def test_format_table_row_writes_the_str_regions_1_based(
    str_regions: tuple[GenomeSegment, ...], expected: str
) -> None:
    assert _format_fields(_make_row(str_regions=str_regions))["str_regions"] == expected


@pytest.mark.parametrize(
    ("qual", "expected"),
    [
        (50.0, "50"),
        (29.53, "29.53"),
        (1e-05, "1e-05"),
        (1e20, "1e+20"),
        (0.0, "0"),
        (3.0e3, "3000"),
        (1234567.0, "1234567"),
        (29.5312345, "29.5312345"),
    ],
    ids=[
        "whole",
        "two_decimals",
        "small",
        "large",
        "zero",
        "whole_thousands",
        "seven_digits",
        "seven_decimals",
    ],
)
def test_format_table_row_writes_the_qual_in_its_shortest_form(qual: float, expected: str) -> None:
    formatted = _format_fields(_make_row(quality_a=qual, quality_b=qual))
    assert (formatted["quality_a"], formatted["quality_b"]) == (expected, expected)


@pytest.mark.parametrize(
    ("kind", "expected"),
    [(Kind.SNP, "SNP"), (Kind.INDEL, "indel"), (Kind.STR, "STR"), (Kind.COMPLEX, "complex")],
    ids=["snp", "indel", "str", "complex"],
)
def test_format_table_row_writes_the_kind_as_its_word(kind: Kind, expected: str) -> None:
    assert _format_fields(_make_row(kind=kind))["kind"] == expected


@pytest.mark.parametrize(
    ("outcome", "expected"),
    [
        (Outcome.TOO_COMPLEX, "too complex"),
        (Outcome.INVALID_IN_A, "invalid in A"),
        (Outcome.INVALID_IN_B, "invalid in B"),
        (Outcome.INVALID_IN_BOTH, "invalid in both"),
        (Outcome.MISSING_IN_A, "missing in A"),
        (Outcome.MISSING_IN_B, "missing in B"),
        (Outcome.MISSING_IN_BOTH, "missing in both"),
        (Outcome.REFERENCE_IN_BOTH, "reference in both"),
        (Outcome.SAME_GENOTYPE, "same genotype"),
        (Outcome.SAME_ALLELES, "same alleles, other genotype"),
        (Outcome.ONLY_IN_A, "only in A"),
        (Outcome.ONLY_IN_B, "only in B"),
        (Outcome.DIFFERENT, "different"),
    ],
    ids=lambda value: str(value).replace(" ", "_").replace(",", ""),
)
def test_format_table_row_writes_the_outcome_as_its_words(outcome: Outcome, expected: str) -> None:
    assert _format_fields(_make_row(outcome=outcome))["outcome"] == expected


def test_format_table_row_takes_each_invalid_reason_from_its_own_file() -> None:
    formatted = _format_fields(
        _make_row(invalid_reason_a=None, invalid_reason_b=InvalidReason.MALFORMED_GT)
    )
    assert (formatted["invalid_reason_a"], formatted["invalid_reason_b"]) == (".", "malformed GT")


@pytest.mark.parametrize(
    ("reason", "expected"),
    [
        (InvalidReason.MALFORMED_GT, "malformed GT"),
        (InvalidReason.PLOIDY_MISMATCH, "ploidy mismatch"),
        (InvalidReason.CONTRADICTING_RECORDS, "contradicting records"),
    ],
    ids=["malformed_gt", "ploidy_mismatch", "contradicting_records"],
)
def test_format_table_row_writes_the_invalid_reason_as_its_words(
    reason: InvalidReason, expected: str
) -> None:
    formatted = _format_fields(_make_row(invalid_reason_a=reason, invalid_reason_b=reason))
    assert (formatted["invalid_reason_a"], formatted["invalid_reason_b"]) == (expected, expected)


# Four rows of one locus, the second reference in both, which the table leaves out
# without all_rows, and the fourth missing in both, which it keeps.
_SHORT_ROWS = [
    _WORKED_ROW,
    _make_row(sample_a="HG003", sample_b="HG003", outcome=Outcome.REFERENCE_IN_BOTH),
    _make_row(sample_a="HG004", sample_b="HG004", outcome=Outcome.DIFFERENT),
    _make_row(sample_a="HG005", sample_b="HG005", outcome=Outcome.MISSING_IN_BOTH),
]
_SHORT_LINES = [
    _WORKED_LINE,
    "chr20\t100\t110\tSTR\tstr\t101-110:str:AT\tHG003\tHG003\treference in both"
    "\tCATATATATAT,CATATATAT\t0/1\t0/0\t1\t1\t29.5\t.\t45\t.\t1\t0\t.\t.\n",
    "chr20\t100\t110\tSTR\tstr\t101-110:str:AT\tHG004\tHG004\tdifferent"
    "\tCATATATATAT,CATATATAT\t0/1\t0/0\t1\t1\t29.5\t.\t45\t.\t1\t0\t.\t.\n",
    "chr20\t100\t110\tSTR\tstr\t101-110:str:AT\tHG005\tHG005\tmissing in both"
    "\tCATATATATAT,CATATATAT\t0/1\t0/0\t1\t1\t29.5\t.\t45\t.\t1\t0\t.\t.\n",
]
# The table without all_rows: every short row but the one reference in both.
_SHORT_TABLE = _HEADER + _SHORT_LINES[0] + _SHORT_LINES[2] + _SHORT_LINES[3]


def test_write_table_leaves_out_reference_in_both_without_all_rows(tmp_path: Path) -> None:
    output = tmp_path / "loci.tsv"
    write_table(iter(_SHORT_ROWS), output)
    assert output.read_text(encoding="utf-8") == _SHORT_TABLE


def test_write_table_writes_every_row_with_all_rows(tmp_path: Path) -> None:
    output = tmp_path / "loci.tsv"
    write_table(iter(_SHORT_ROWS), output, all_rows=True)
    assert output.read_text(encoding="utf-8") == _HEADER + "".join(_SHORT_LINES)


def test_write_table_with_no_rows_writes_the_header_lines(tmp_path: Path) -> None:
    output, summary = tmp_path / "loci.tsv", tmp_path / "summary.tsv"
    write_table(iter([]), output, summary=summary)
    assert output.read_text(encoding="utf-8") == _HEADER
    assert summary.read_text(encoding="utf-8") == _SUMMARY_HEADER


def test_write_table_leaves_no_temporary_file(tmp_path: Path) -> None:
    write_table(iter(_SHORT_ROWS), tmp_path / "loci.tsv", summary=tmp_path / "summary.tsv")
    assert sorted(path.name for path in tmp_path.iterdir()) == ["loci.tsv", "summary.tsv"]


def test_write_table_gzips_a_path_ending_in_gz(tmp_path: Path) -> None:
    output, summary = tmp_path / "loci.tsv.gz", tmp_path / "summary.tsv.gz"
    write_table(iter(_SHORT_ROWS), output, summary=summary)
    assert output.read_bytes()[:2] == b"\x1f\x8b"  # gzip's magic number
    assert gzip.decompress(output.read_bytes()).decode() == _SHORT_TABLE
    assert gzip.decompress(summary.read_bytes()).decode() == (
        _SUMMARY_HEADER
        + "missing in both\tstr\t0\t0\t1\t0\t1\n"
        + "reference in both\tstr\t0\t0\t1\t0\t1\n"
        + "only in A\tstr\t0\t0\t1\t0\t1\n"
        + "different\tstr\t0\t0\t1\t0\t1\n"
    )
    assert sorted(path.name for path in tmp_path.iterdir()) == ["loci.tsv.gz", "summary.tsv.gz"]


def test_write_table_gzips_at_level_6(tmp_path: Path) -> None:
    output = tmp_path / "loci.tsv.gz"
    write_table(iter(_SHORT_ROWS), output)
    # XFL, the ninth byte of the header, is 0 at level 6, 2 at level 9 and 4 at level 1.
    assert output.read_bytes()[8] == 0


def test_write_table_gzips_with_no_file_name_in_the_header(tmp_path: Path) -> None:
    output = tmp_path / "loci.tsv.gz"
    write_table(iter(_SHORT_ROWS), output)
    header = output.read_bytes()
    # FLG, the fourth byte, has FNAME (8) clear, so `gunzip -N` restores no `.tmp` name;
    # MTIME, bytes 4 to 7, is 0.
    assert (header[3] & 8, header[4:8]) == (0, bytes(4))


def test_write_table_writes_utf8_with_no_carriage_return(tmp_path: Path) -> None:
    output, summary = tmp_path / "loci.tsv", tmp_path / "summary.tsv"
    write_table(iter([_make_row(sample_a="Müller")]), output, summary=summary)
    table_bytes = output.read_bytes()
    assert b"\tM\xc3\xbcller\tHG002\t" in table_bytes
    assert b"\r" not in table_bytes + summary.read_bytes()


class _StoppedError(Exception):
    """What the rows raise in the test of a run that stops."""


def _yield_two_rows_then_raise(error: BaseException) -> Iterator[OutcomeRow]:
    yield _SHORT_ROWS[0]
    yield _SHORT_ROWS[2]
    raise error


@pytest.mark.parametrize("earlier", [True, False], ids=["earlier_files", "no_earlier_files"])
@pytest.mark.parametrize(
    "error", [_StoppedError("stopped"), KeyboardInterrupt()], ids=["error", "ctrl_c"]
)
def test_rows_that_raise_leave_the_files_as_they_were(
    tmp_path: Path, earlier: bool, error: BaseException
) -> None:
    output, summary = tmp_path / "loci.tsv", tmp_path / "summary.tsv"
    if earlier:
        output.write_text("earlier table\n", encoding="utf-8")
        summary.write_text("earlier summary\n", encoding="utf-8")
    with pytest.raises(type(error)) as raised:
        write_table(_yield_two_rows_then_raise(error), output, summary=summary)
    assert raised.value is error
    if earlier:
        assert output.read_text(encoding="utf-8") == "earlier table\n"
        assert summary.read_text(encoding="utf-8") == "earlier summary\n"
        assert sorted(path.name for path in tmp_path.iterdir()) == ["loci.tsv", "summary.tsv"]
    else:
        assert list(tmp_path.iterdir()) == []


def test_a_summary_that_cannot_be_renamed_leaves_the_earlier_table(tmp_path: Path) -> None:
    output, summary = tmp_path / "loci.tsv", tmp_path / "summary.tsv"
    output.write_text("earlier table\n", encoding="utf-8")

    def yield_rows_then_take_the_summarys_path() -> Iterator[OutcomeRow]:
        yield _WORKED_ROW
        summary.mkdir()  # after the check before the first row, so the rename fails

    with pytest.raises(PvcEvalInputError) as raised:
        write_table(yield_rows_then_take_the_summarys_path(), output, summary=summary)
    assert str(raised.value) == (
        f"--summary {summary} could not be written: Is a directory. Give another --summary,"
        " and run again."
    )
    # The summary is renamed first, so the table the user reads is the earlier one.
    assert output.read_text(encoding="utf-8") == "earlier table\n"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["loci.tsv", "summary.tsv"]


def _make_path_errors(tmp_path: Path, option: str) -> list[tuple[str, Path, str]]:
    """Each path of option that cannot be written: its id, the path, and the message, as
    `str-regions` gives them for --output. The faults are made by _make_path_fault."""
    missing = tmp_path / "missing" / "out.tsv"
    through_file = tmp_path / "a_file" / "out.tsv"
    read_only = tmp_path / "read_only" / "out.tsv"
    taken = tmp_path / "taken.tsv"
    return [
        (
            "missing_directory",
            missing,
            f"{option} {missing} cannot be written: its directory, {missing.parent}, does not"
            f" exist. Make it, or give another {option}, and run again.",
        ),
        (
            "a_directory",
            tmp_path,
            f"{option} {tmp_path} is a directory. Give the path of a file to write, and run again.",
        ),
        (
            "through_a_file",
            through_file,
            f"{option} {through_file} cannot be written: its path goes through a file as if it"
            " were a directory. Correct the path, and run again.",
        ),
        (
            "permission_denied",
            read_only,
            f"{option} {read_only} cannot be written: permission denied to write in"
            f" {read_only.parent}. Allow writing there, or give another {option}, and run"
            " again.",
        ),
        (
            "temporary_is_a_directory",
            taken,
            f"{option} {taken} cannot be written: {taken}.tmp, where it is written first, is a"
            " directory. Move it away, and run again.",
        ),
    ]


def _make_path_fault(case_id: str, tmp_path: Path) -> None:
    if case_id == "through_a_file":
        (tmp_path / "a_file").write_text("", encoding="utf-8")
    elif case_id == "permission_denied":
        (tmp_path / "read_only").mkdir(mode=0o555)
    elif case_id == "temporary_is_a_directory":
        (tmp_path / "taken.tsv.tmp").mkdir()


_PATH_ERROR_IDS = [case[0] for case in _make_path_errors(Path("out"), "--output")]


def _refuse_to_be_read() -> Iterator[OutcomeRow]:
    """Rows that fail the test when the first is read."""
    pytest.fail("a row was read before the files were opened")
    yield _WORKED_ROW


def _find_path_error(tmp_path: Path, option: str, case_id: str) -> tuple[Path, str]:
    _make_path_fault(case_id, tmp_path)
    return next(
        (path, message)
        for name, path, message in _make_path_errors(tmp_path, option)
        if name == case_id
    )


@pytest.mark.parametrize("case_id", _PATH_ERROR_IDS)
def test_an_output_that_cannot_be_written_is_refused_before_any_row(
    tmp_path: Path, case_id: str
) -> None:
    output, message = _find_path_error(tmp_path, "--output", case_id)
    with pytest.raises(PvcEvalInputError) as raised:
        write_table(_refuse_to_be_read(), output, summary=tmp_path / "summary.tsv")
    assert str(raised.value) == message
    assert not (tmp_path / "summary.tsv.tmp").exists()


@pytest.mark.parametrize("case_id", _PATH_ERROR_IDS)
def test_a_summary_that_cannot_be_written_is_refused_before_any_row(
    tmp_path: Path, case_id: str
) -> None:
    summary, message = _find_path_error(tmp_path, "--summary", case_id)
    output = tmp_path / "loci.tsv"
    with pytest.raises(PvcEvalInputError) as raised:
        write_table(_refuse_to_be_read(), output, summary=summary)
    assert str(raised.value) == message
    assert not output.exists()
    assert not (tmp_path / "loci.tsv.tmp").exists()


@pytest.mark.parametrize("option", ["--output", "--summary"])
def test_a_name_too_long_is_refused_with_the_systems_words(tmp_path: Path, option: str) -> None:
    too_long = tmp_path / ("x" * 300 + ".tsv")
    output, summary = (
        (too_long, tmp_path / "summary.tsv")
        if option == "--output"
        else (tmp_path / "loci.tsv", too_long)
    )
    with pytest.raises(PvcEvalInputError) as raised:
        write_table(_refuse_to_be_read(), output, summary=summary)
    assert str(raised.value) == (
        f"{option} {too_long} could not be written: File name too long. Give another {option},"
        " and run again."
    )
    assert list(tmp_path.iterdir()) == []


# Writes a table larger than the limit on the size of a file the process may write, as a
# full disk would stop it, and prints the message of the error.
_WRITE_PAST_THE_LIMIT = """
import resource, signal, sys
from pathlib import Path
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.kind import Kind
from pvc_eval.outcome import Outcome, OutcomeRow
from pvc_eval.table import write_table

row = OutcomeRow(
    contig="chr20", start=99, end=110, kind=Kind.SNP, stratum=None, str_regions=(),
    sample_a="HG002", sample_b="HG002", outcome=Outcome.SAME_GENOTYPE, alleles=("C", "T"),
    genotype_a=(0, 1), genotype_b=(0, 1), possible_a=1, possible_b=1, quality_a=29.5,
    quality_b=29.5, gq_a=45, gq_b=45, records_a=1, records_b=1, invalid_reason_a=None,
    invalid_reason_b=None,
)
signal.signal(signal.SIGXFSZ, signal.SIG_IGN)  # so a write past the limit raises
resource.setrlimit(resource.RLIMIT_FSIZE, (50_000, 50_000))
directory = Path(sys.argv[1])
try:
    write_table(
        (row for _ in range(10_000)), directory / "loci.tsv", summary=directory / "summary.tsv"
    )
except PvcEvalInputError as err:
    print(err)
"""


def test_a_full_disk_gives_the_specs_message() -> None:
    # The example of doc/specs/compare.md, "Writing the files", which a test cannot make
    # the disk give.
    full_disk = OSError(errno.ENOSPC, "No space left on device")
    assert str(make_write_error("--output", Path("loci.tsv"), full_disk)) == (
        "--output loci.tsv could not be written: No space left on device. Free some space, or"
        " give another --output, and run again."
    )


def test_a_write_that_fails_is_refused_and_leaves_no_temporary_file(tmp_path: Path) -> None:
    completed = subprocess.run(
        [sys.executable, "-c", _WRITE_PAST_THE_LIMIT, str(tmp_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    assert completed.stdout == (
        f"--output {tmp_path / 'loci.tsv'} could not be written: File too large. Give another"
        " --output, and run again.\n"
    )
    assert list(tmp_path.iterdir()) == []


def _is_case_insensitive(directory: Path) -> bool:
    """Whether the file system of directory takes `a` and `A` for one name."""
    probe = directory / "probe"
    probe.write_text("", encoding="utf-8")
    insensitive = (directory / "PROBE").exists()
    probe.unlink()
    return insensitive


# Spellings of --output and --summary that are one file, or one the other's .tmp, and
# whether the test needs a file system that ignores letter case.
_SAME_FILES = [
    ("loci.tsv", "loci.tsv", False),
    ("loci.tsv", "sub/../loci.tsv", False),
    ("loci.tsv", "LOCI.tsv", True),
    ("s.tsv.tmp", "s.tsv", False),
    ("s.tsv", "s.tsv.tmp", False),
    ("S.TSV.TMP", "s.tsv", True),
    ("s.tsv", "S.TSV.TMP", True),
]
_SAME_FILE_IDS = [
    "same_path",
    "through_parent",
    "letter_case",
    "output_is_summarys_tmp",
    "summary_is_outputs_tmp",
    "output_is_summarys_tmp_in_other_case",
    "summary_is_outputs_tmp_in_other_case",
]


@pytest.mark.parametrize("earlier", [False, True], ids=["new_files", "earlier_files"])
@pytest.mark.parametrize(
    ("output_name", "summary_name", "needs_case"), _SAME_FILES, ids=_SAME_FILE_IDS
)
def test_a_summary_that_is_the_output_or_its_tmp_is_refused(
    tmp_path: Path, output_name: str, summary_name: str, needs_case: bool, earlier: bool
) -> None:
    if needs_case and not _is_case_insensitive(tmp_path):
        pytest.skip("the file system tells letter case apart")
    (tmp_path / "sub").mkdir()
    output, summary = tmp_path / output_name, tmp_path / summary_name
    if earlier:
        output.write_text("earlier\n", encoding="utf-8")
    before = sorted(path.name for path in tmp_path.iterdir())
    with pytest.raises(PvcEvalInputError) as raised:
        write_table(_refuse_to_be_read(), output, summary=summary)
    assert str(raised.value) == (
        f"--summary {summary} is the same file as --output {output}. Give the summary a path"
        " of its own, and run again."
    )
    assert sorted(path.name for path in tmp_path.iterdir()) == before
    if earlier:
        assert output.read_text(encoding="utf-8") == "earlier\n"


def _make_summary_row(
    contig: str, kind: Kind, stratum: str | None, sample: str, outcome: Outcome
) -> OutcomeRow:
    return _make_row(
        contig=contig, kind=kind, stratum=stratum, sample_a=sample, sample_b=sample, outcome=outcome
    )


# Five loci, each with the sample pairs HG002 and HG003, over the strata generic and str.
# The first locus's str comes before any generic, and its only in A before any same
# genotype, so that the summary's order is not the order of the rows.
_SUMMARY_ROWS = [
    _make_summary_row("chr1", Kind.STR, "str", "HG002", Outcome.ONLY_IN_A),
    _make_summary_row("chr1", Kind.STR, "str", "HG003", Outcome.SAME_GENOTYPE),
    _make_summary_row("chr1", Kind.SNP, "generic", "HG002", Outcome.SAME_GENOTYPE),
    _make_summary_row("chr1", Kind.SNP, "generic", "HG003", Outcome.REFERENCE_IN_BOTH),
    _make_summary_row("chr1", Kind.INDEL, "generic", "HG002", Outcome.SAME_GENOTYPE),
    _make_summary_row("chr1", Kind.INDEL, "generic", "HG003", Outcome.DIFFERENT),
    _make_summary_row("chr2", Kind.COMPLEX, "str", "HG002", Outcome.SAME_GENOTYPE),
    _make_summary_row("chr2", Kind.COMPLEX, "str", "HG003", Outcome.REFERENCE_IN_BOTH),
    _make_summary_row("chr2", Kind.SNP, "generic", "HG002", Outcome.REFERENCE_IN_BOTH),
    _make_summary_row("chr2", Kind.SNP, "generic", "HG003", Outcome.SAME_GENOTYPE),
]


def test_the_summary_counts_by_outcome_and_stratum_over_every_sample_pair(
    tmp_path: Path,
) -> None:
    output, summary = tmp_path / "loci.tsv", tmp_path / "summary.tsv"
    write_table(iter(_SUMMARY_ROWS), output, summary=summary)
    assert summary.read_text(encoding="utf-8").splitlines() == [
        "outcome\tstratum\tSNP\tindel\tSTR\tcomplex\ttotal",
        "reference in both\tgeneric\t2\t0\t0\t0\t2",
        "reference in both\tstr\t0\t0\t0\t1\t1",
        "same genotype\tgeneric\t2\t1\t0\t0\t3",
        "same genotype\tstr\t0\t0\t1\t1\t2",
        "only in A\tstr\t0\t0\t1\t0\t1",
        "different\tgeneric\t0\t1\t0\t0\t1",
    ]
    # The table, without all_rows, leaves out the three rows reference in both.
    assert len(output.read_text(encoding="utf-8").splitlines()) == 1 + 7


def test_the_summary_sorts_the_strata_by_code_point(tmp_path: Path) -> None:
    summary = tmp_path / "summary.tsv"
    rows = [
        _make_summary_row("chr1", Kind.SNP, "generic", "HG002", Outcome.SAME_GENOTYPE),
        _make_summary_row("chr1", Kind.STR, "Str", "HG002", Outcome.SAME_GENOTYPE),
    ]
    write_table(iter(rows), tmp_path / "loci.tsv", summary=summary)
    assert summary.read_text(encoding="utf-8").splitlines() == [
        "outcome\tstratum\tSNP\tindel\tSTR\tcomplex\ttotal",
        "same genotype\tStr\t0\t0\t1\t0\t1",
        "same genotype\tgeneric\t1\t0\t0\t0\t1",
    ]


def test_the_summary_without_strata_has_a_dot_for_the_stratum(tmp_path: Path) -> None:
    summary = tmp_path / "summary.tsv"
    rows = [
        _make_summary_row("chr1", Kind.SNP, None, "HG002", Outcome.SAME_GENOTYPE),
        _make_summary_row("chr1", Kind.STR, None, "HG002", Outcome.DIFFERENT),
        _make_summary_row("chr2", Kind.SNP, None, "HG002", Outcome.SAME_GENOTYPE),
    ]
    write_table(iter(rows), tmp_path / "loci.tsv", summary=summary, all_rows=True)
    assert summary.read_text(encoding="utf-8").splitlines() == [
        "outcome\tstratum\tSNP\tindel\tSTR\tcomplex\ttotal",
        "same genotype\t.\t2\t0\t0\t0\t2",
        "different\t.\t0\t0\t1\t0\t1",
    ]
