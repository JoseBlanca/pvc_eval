"""The `pvc_eval str-regions` command and read_str_regions (doc/specs/str_finder.md,
"read_str_regions" and "The str-regions command")."""

import os
import re
import subprocess
import sysconfig
from collections.abc import Sequence
from pathlib import Path

import pytest

from pvc_eval.cli import main
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import GenomeSegment, RegionType, read_str_regions

# The worked cases of doc/specs/str_finder.md. Their STR regions are the spec's, which `pvc
# repeat-catalog` and then `pvc dev type-regions --min-copies 5,5,4,4,4,3` gave, pvc built at
# bbcc9f4.
_FOUR_A = (
    "GCTAGACATACATACATACAGTCAGCACGACTGTCAGTGATCGCTAGTAGTAGTGATGCACACACACACACACACACACAGTCTGT"
    "AAAACTACGCTACTGCTGTCACTCGACTGCATACATCAGACAGACTCGTATGACAGTCACGCAGA"
)
_FIVE_A = _FOUR_A.replace("AAAA", "AAAAA")
_NEAR_END = (
    "GCTAGACATGAAAAAACTGTCAGTCACACACACACACACACACAGTCTGTAGCTACTGCTGTCACTCGACTGCATACATCAGACAG"
    "ACTCGTATGACAGTCACG"
)
_FAR_END = (
    "GCTAGACATGCATCAGTAAAAAACTGTCAGTCACACACACACACACACACAGTCTGTAGCTACTGCTGTCACTCGACTGCATACATC"
    "AGACAGACTCGTATGACAGTCACG"
)
# In an order of the .fai that is not the order of the names, so that rows sorted by name
# would differ. `twice` is four_a twice over: its (CA) repeat at 58-80 again at 151 + 58,
# with no repeat across the join that reaches a copy floor (it holds AGAG, two copies of
# period 2, under the floor of 5), so two STR segments, in order of start. `plain` has none.
_WORKED_CONTIGS = {
    "four_a": _FOUR_A,
    "five_a": _FIVE_A,
    "near_end": _NEAR_END,
    "far_end": _FAR_END,
    "twice": _FOUR_A + _FOUR_A,
    "plain": "ACGTTGCA",
}
_WORKED_ROWS = [
    "four_a\t58\t80\tstr\tCA",
    "five_a\t58\t91\tstr-bundle\t.",
    "near_end\t24\t44\tstr\tCA",
    "far_end\t17\t51\tstr-bundle\t.",
    "twice\t58\t80\tstr\tCA",
    "twice\t209\t231\tstr\tCA",
]
_WORKED_SEGMENTS = [
    GenomeSegment("four_a", 58, 80, RegionType.STR, "CA"),
    GenomeSegment("five_a", 58, 91, RegionType.STR_BUNDLE, None),
    GenomeSegment("near_end", 24, 44, RegionType.STR, "CA"),
    GenomeSegment("far_end", 17, 51, RegionType.STR_BUNDLE, None),
    GenomeSegment("twice", 58, 80, RegionType.STR, "CA"),
    GenomeSegment("twice", 209, 231, RegionType.STR, "CA"),
]
_DEFAULT_HEADER = [
    "##pvc_eval str-regions 0.1.0",
    "##scan_min_period: 1",
    "##scan_max_period: 6",
    "##min_period: 1",
    "##max_period: 6",
    "##match_reward: 2",
    "##mismatch_penalty: 7",
    "##scan_min_copies: 2",
    "##copy_floors: 5,5,4,4,4,3",
    "##min_purity: 0.8",
    "##min_score: 0",
    "##bundle_distance: 15",
    "##contig_end_flank: 15",
    "##max_str_len: 100",
    "##max_bundle_len: 100",
    "#chrom\tstart\tend\ttype\tmotif",
]


def _write_fasta(directory: Path, contigs: dict[str, str]) -> Path:
    """A FASTA of the contigs, one line each, and its .fai, as samtools faidx writes it."""
    fasta = ""
    fai = ""
    for name, bases in contigs.items():
        fasta += f">{name}\n"
        fai += f"{name}\t{len(bases)}\t{len(fasta)}\t{len(bases)}\t{len(bases) + 1}\n"
        fasta += f"{bases}\n"
    path = directory / "ref.fa"
    path.write_text(fasta, encoding="ascii")
    Path(f"{path}.fai").write_text(fai, encoding="utf-8")
    return path


def _write_lines(path: Path, lines: Sequence[str]) -> Path:
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return path


def _run(argv: Sequence[str], capsys: pytest.CaptureFixture[str]) -> tuple[int, str]:
    """The exit status of `pvc_eval` with these arguments, and what it wrote to stderr."""
    try:
        main(argv)
    except SystemExit as exit_:
        status = exit_.code
    else:
        status = 0
    assert isinstance(status, int)
    return status, capsys.readouterr().err


def _run_str_regions(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], *options: str
) -> list[str]:
    """The lines of the BED `str-regions` writes for the worked contigs, with these options."""
    fasta = _write_fasta(tmp_path, _WORKED_CONTIGS)
    output = tmp_path / "str.bed"
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(output), *options], capsys
    )
    assert (status, stderr) == (0, "")
    return output.read_text(encoding="utf-8").splitlines()


# ---- The BED str-regions writes ---------------------------------------------------------


def test_the_worked_cases_give_pvcs_str_regions_in_fai_order(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert _run_str_regions(tmp_path, capsys) == _DEFAULT_HEADER + _WORKED_ROWS


def test_the_header_gives_every_parameter_the_options_set(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Every option away from its default, each to a value no other field has, so that an
    # option that sets the wrong field, or none, shows in its ## line.
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT"})
    output = tmp_path / "str.bed"
    options = [
        *("--scan-min-period", "2", "--scan-max-period", "5"),
        *("--min-period", "3", "--max-period", "4"),
        *("--scan-match-reward", "9", "--scan-mismatch-penalty", "8", "--scan-min-copies", "7"),
        *("--min-copies", "8,6,6,6,5,4", "--min-purity", "0.9", "--min-score", "1"),
        *("--bundle-distance", "20", "--contig-end-flank", "10"),
        *("--max-str-len", "150", "--max-str-bundle-len", "120"),
    ]
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(output), *options], capsys
    )
    assert (status, stderr) == (0, "")
    assert output.read_text(encoding="utf-8").splitlines() == [
        "##pvc_eval str-regions 0.1.0",
        "##scan_min_period: 2",
        "##scan_max_period: 5",
        "##min_period: 3",
        "##max_period: 4",
        "##match_reward: 9",
        "##mismatch_penalty: 8",
        "##scan_min_copies: 7",
        "##copy_floors: 8,6,6,6,5,4",
        "##min_purity: 0.9",
        "##min_score: 1",
        "##bundle_distance: 20",
        "##contig_end_flank: 10",
        "##max_str_len: 150",
        "##max_bundle_len: 120",
        "#chrom\tstart\tend\ttype\tmotif",
    ]


def test_the_options_reach_the_finder(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # With a period-1 copy floor of 8, the five A's of five_a are dropped at step 3, so its
    # (CA) repeat stands alone, as in four_a, and is an STR segment.
    lines = _run_str_regions(tmp_path, capsys, "--min-copies", "8,6,6,6,5,4")
    assert "five_a\t58\t80\tstr\tCA" in lines


def test_a_finished_run_leaves_no_temporary_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _run_str_regions(tmp_path, capsys)
    assert sorted(path.name for path in tmp_path.iterdir()) == ["ref.fa", "ref.fa.fai", "str.bed"]


def test_a_run_that_fails_leaves_an_earlier_output_as_it_was(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # chr1 is segmented and its rows written before chr2's byte that is not a letter stops
    # the run; the rows go to the temporary file, which is removed.
    fasta = _write_fasta(tmp_path, {"chr1": _FOUR_A, "chr2": "AC-GT"})
    output = _write_lines(tmp_path / "str.bed", ["earlier"])
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(output)], capsys
    )
    assert status == 1
    assert "chr2:3 is not a letter" in stderr
    assert output.read_text(encoding="utf-8") == "earlier\n"
    assert not Path(f"{output}.tmp").exists()


def test_the_installed_command_runs(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    output = tmp_path / "str.bed"
    command = Path(sysconfig.get_path("scripts")) / "pvc_eval"
    subprocess.run([command, "str-regions", "--reference", fasta, "--output", output], check=True)
    assert output.read_text(encoding="utf-8").splitlines()[-1] == "four_a\t58\t80\tstr\tCA"


# ---- --regions --------------------------------------------------------------------------


def test_a_region_that_cuts_an_str_region_writes_it_whole(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    regions = _write_lines(tmp_path / "regions.bed", ["four_a\t70\t75", "twice\t220\t300"])
    lines = _run_str_regions(tmp_path, capsys, "--regions", str(regions))
    assert lines == [*_DEFAULT_HEADER, "four_a\t58\t80\tstr\tCA", "twice\t209\t231\tstr\tCA"]


@pytest.mark.parametrize(
    ("region", "expected"),
    [
        ("four_a\t79\t80", ["four_a\t58\t80\tstr\tCA"]),
        ("four_a\t80\t90", []),
        ("four_a\t50\t58", []),
        ("four_a\t57\t58\nfour_a\t0\t58", []),
    ],
    ids=["last_base", "just_after", "just_before", "two_before"],
)
def test_a_region_must_overlap_an_str_region_by_a_base(
    region: str, expected: list[str], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    regions = _write_lines(tmp_path / "regions.bed", [region])
    lines = _run_str_regions(tmp_path, capsys, "--regions", str(regions))
    assert lines == _DEFAULT_HEADER + expected


def test_a_contig_the_regions_do_not_name_is_not_segmented(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # chr2 has a byte that is not a letter, which stops the run if chr2 is read.
    fasta = _write_fasta(tmp_path, {"chr2": "AC-GT", "four_a": _FOUR_A})
    regions = _write_lines(
        tmp_path / "regions.bed", ["# a comment", "track name=x", "", "four_a\t0\t151"]
    )
    output = tmp_path / "str.bed"
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(output), "--regions", str(regions)),
        ],
        capsys,
    )
    assert (status, stderr) == (0, "")
    assert output.read_text(encoding="utf-8").splitlines()[-1] == "four_a\t58\t80\tstr\tCA"


def _make_regions_errors(regions: Path) -> list[tuple[str, list[str], str]]:
    """Each refusal of a --regions BED: its id, its lines, and the message."""
    remedy = " Correct the line or remove it, and run again."
    return [
        (
            "contig_not_in_reference",
            ["four_a\t0\t10", "chrX\t0\t10"],
            f"{regions}, line 2: contig chrX is not in the reference. Use the reference the BED"
            " was made for, or remove the lines of chrX, and run again.",
        ),
        (
            "start_not_a_number",
            ["four_a\tten\t20"],
            f"{regions}, line 1: the start is 'ten', and it must be a whole number of 0 or"
            f" more.{remedy}",
        ),
        (
            "start_negative",
            ["four_a\t-1\t20"],
            f"{regions}, line 1: the start is '-1', and it must be a whole number of 0 or"
            f" more.{remedy}",
        ),
        (
            "end_not_a_number",
            ["four_a\t10\t2e1"],
            f"{regions}, line 1: the end is '2e1', and it must be a whole number of 0 or"
            f" more.{remedy}",
        ),
        (
            "end_at_start",
            ["four_a\t10\t10"],
            f"{regions}, line 1: the end, 10, is not after the start, 10, and a region holds at"
            f" least one base.{remedy}",
        ),
        (
            "end_past_contig",
            ["four_a\t10\t152"],
            f"{regions}, line 1: the end, 152, is past the end of four_a, which has 151 bases in"
            " the reference. Use the reference the BED was made for, or correct the line, and"
            " run again.",
        ),
        (
            "two_columns",
            ["four_a\t10"],
            f"{regions}, line 1 has 2 tab-separated columns, and a BED line needs at least"
            f" three: contig, start and end.{remedy}",
        ),
        (
            "one_column",
            ["four_a 10 20"],
            f"{regions}, line 1 has 1 tab-separated column, and a BED line needs at least"
            f" three: contig, start and end.{remedy}",
        ),
    ]


_REGIONS_ERROR_IDS = [case[0] for case in _make_regions_errors(Path("regions.bed"))]


@pytest.mark.parametrize("case_id", _REGIONS_ERROR_IDS)
def test_each_refusal_of_a_regions_bed_names_the_line(
    case_id: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    regions = tmp_path / "regions.bed"
    (lines, message) = next(
        (lines, message)
        for name, lines, message in _make_regions_errors(regions)
        if name == case_id
    )
    _write_lines(regions, lines)
    output = tmp_path / "str.bed"
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(output), "--regions", str(regions)),
        ],
        capsys,
    )
    assert (status, stderr) == (1, f"{message}\n")
    assert not output.exists()
    assert not Path(f"{output}.tmp").exists()


# ---- The errors of the command ----------------------------------------------------------


def test_a_missing_reference_stops_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fasta = tmp_path / "ref.fa"
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(tmp_path / "str.bed")], capsys
    )
    assert (status, stderr) == (1, f"{fasta} does not exist. Check its path and run again.\n")


def test_a_missing_regions_bed_stops_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    regions = tmp_path / "regions.bed"
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(tmp_path / "str.bed"), "--regions", str(regions)),
        ],
        capsys,
    )
    assert (status, stderr) == (1, f"{regions} does not exist. Check its path and run again.\n")


def _make_output_errors(tmp_path: Path) -> list[tuple[str, Path, str]]:
    """Each --output that cannot be written: its id, the path, and the message. The
    paths are made by _make_output_fault."""
    missing = tmp_path / "missing" / "str.bed"
    through_file = tmp_path / "a_file" / "str.bed"
    read_only = tmp_path / "read_only" / "str.bed"
    taken = tmp_path / "taken.bed"
    return [
        (
            "missing_directory",
            missing,
            f"--output {missing} cannot be written: its directory, {missing.parent}, does not"
            " exist. Make it, or give another --output, and run again.",
        ),
        (
            "a_directory",
            tmp_path,
            f"--output {tmp_path} is a directory. Give the path of a file to write, and run again.",
        ),
        (
            "through_a_file",
            through_file,
            f"--output {through_file} cannot be written: its path goes through a file as if it"
            " were a directory. Correct the path, and run again.",
        ),
        (
            "permission_denied",
            read_only,
            f"--output {read_only} cannot be written: permission denied to write in"
            f" {read_only.parent}. Allow writing there, or give another --output, and run"
            " again.",
        ),
        (
            "temporary_is_a_directory",
            taken,
            f"--output {taken} cannot be written: {taken}.tmp, where it is written first, is a"
            " directory. Move it away, and run again.",
        ),
    ]


def _make_output_fault(case_id: str, tmp_path: Path) -> None:
    if case_id == "through_a_file":
        (tmp_path / "a_file").write_text("", encoding="utf-8")
    elif case_id == "permission_denied":
        (tmp_path / "read_only").mkdir(mode=0o555)
    elif case_id == "temporary_is_a_directory":
        (tmp_path / "taken.bed.tmp").mkdir()


_OUTPUT_ERROR_IDS = [case[0] for case in _make_output_errors(Path("out"))]


@pytest.mark.parametrize("case_id", _OUTPUT_ERROR_IDS)
def test_an_output_that_cannot_be_written_stops_the_run_before_segmenting(
    case_id: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # chr1 has a byte that is not a letter, which stops the run if chr1 is segmented first.
    if case_id == "permission_denied" and os.geteuid() == 0:
        pytest.skip("root writes in a directory without write permission")
    fasta = _write_fasta(tmp_path, {"chr1": "AC-GT"})
    (output, message) = next(
        (output, message)
        for name, output, message in _make_output_errors(tmp_path)
        if name == case_id
    )
    _make_output_fault(case_id, tmp_path)
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(output)], capsys
    )
    assert (status, stderr) == (1, f"{message}\n")


def test_an_empty_output_path_stops_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    # An empty --output is the current directory to pathlib.
    monkeypatch.chdir(tmp_path)
    fasta = _write_fasta(tmp_path, {"chr1": "AC-GT"})
    status, stderr = _run(["str-regions", "--reference", str(fasta), "--output", ""], capsys)
    assert (status, stderr) == (
        1,
        "--output . is a directory. Give the path of a file to write, and run again.\n",
    )


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (
            ["--min-purity", "2"],
            "--min-purity is 2.0; it must be a float or int from 0 to 1",
        ),
        (
            ["--min-period", "4", "--max-period", "3"],
            "--min-period (4) must not exceed --max-period (3)",
        ),
        (
            ["--scan-min-period", "0"],
            "--scan-min-period is 0; it must be 1 or more",
        ),
        (
            ["--min-copies", "5,5,4"],
            "--min-copies is 5,5,4; it must be six whole numbers, one for each period from 1 to 6",
        ),
        (
            ["--min-copies", "5,5,4,4,4,-1"],
            "--min-copies is 5,5,4,4,4,-1; each must be a whole number, 0 or more",
        ),
        (
            ["--max-str-len", "10"],
            "--max-str-len is 10; it must be --bundle-distance (15) or more",
        ),
        (
            ["--max-str-bundle-len", "-1"],
            "--max-str-bundle-len is -1; it must be 0 or more",
        ),
        (
            ["--scan-match-reward", "0"],
            "--scan-match-reward is 0; it must be 1 or more",
        ),
    ],
    ids=[
        "min_purity",
        "min_period_above_max",
        "scan_min_period",
        "min_copies_three",
        "min_copies_negative",
        "max_str_len",
        "max_str_bundle_len",
        "scan_match_reward",
    ],
)
def test_an_option_outside_its_limit_stops_the_run_naming_the_option(
    options: list[str], message: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    output = tmp_path / "str.bed"
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(output), *options], capsys
    )
    assert (status, stderr) == (1, f"{message}. Correct the option, and run again.\n")
    assert not output.exists()


def test_help_lists_every_option_with_its_default(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as raised:
        main(["str-regions", "--help"])
    assert raised.value.code == 0
    # argparse starts each option's entry on a line of its own, indented by two spaces.
    entries = {
        entry.split()[0]: " ".join(entry.split())
        for entry in re.split(r"\n  (?=--)", capsys.readouterr().out)[1:]
    }
    assert {
        option: re.findall(r"\(default ([^)]*)\)", entry) for option, entry in entries.items()
    } == {
        "--reference": [],
        "--output": [],
        "--regions": [],
        "--scan-min-period": ["1"],
        "--scan-max-period": ["6"],
        "--min-period": ["1"],
        "--max-period": ["6"],
        "--min-copies": ["5,5,4,4,4,3"],
        "--min-purity": ["0.8"],
        "--min-score": ["0"],
        "--bundle-distance": ["15"],
        "--contig-end-flank": ["15"],
        "--max-str-len": ["100"],
        "--max-str-bundle-len": ["100"],
        "--scan-match-reward": ["2"],
        "--scan-mismatch-penalty": ["7"],
        "--scan-min-copies": ["2"],
    }


# ---- read_str_regions -------------------------------------------------------------------


def test_the_bed_str_regions_writes_reads_back_to_the_same_str_regions(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _run_str_regions(tmp_path, capsys)
    with IndexedFasta(tmp_path / "ref.fa") as reference:
        assert list(read_str_regions(tmp_path / "str.bed", reference)) == _WORKED_SEGMENTS


def test_a_bed_of_three_four_and_more_columns_is_read(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT" * 50, "chr2": "ACGT" * 50})
    bed = _write_lines(
        tmp_path / "str.bed",
        [
            "browser position chr1:1-100",
            "track name=str",
            "# a comment",
            "",
            "chr1\t10\t20",
            "chr1\t20\t30\tstr-bundle",
            "chr1\t30\t40\tsatellite\t.\tmore\tcolumns",
            "chr2\t5\t11\tstr\tACG\t0.9",
            "chr2\t40\t50\tstr\tACGTAC",
        ],
    )
    with IndexedFasta(fasta) as reference:
        assert list(read_str_regions(bed, reference)) == [
            GenomeSegment("chr1", 10, 20, None, None),
            GenomeSegment("chr1", 20, 30, RegionType.STR_BUNDLE, None),
            GenomeSegment("chr1", 30, 40, RegionType.SATELLITE, None),
            GenomeSegment("chr2", 5, 11, RegionType.STR, "ACG"),
            GenomeSegment("chr2", 40, 50, RegionType.STR, "ACGTAC"),
        ]


def _make_str_regions_errors(bed: Path) -> list[tuple[str, list[str], str]]:
    """Each refusal of read_str_regions: its id, the lines of the BED, and the message."""
    remedy = " Correct the line or remove it, and run again."
    sort = (
        " Sort it by contig, in the order of the reference's .fai, and by start, for example"
        f" with `bedtools sort -faidx <the reference's .fai> -i {bed}`, and run again."
    )
    return [
        (
            "two_columns",
            ["chr1\t10\t20", "chr1\t30"],
            f"{bed}, line 2 has 2 tab-separated columns, and a BED line needs at least three:"
            f" contig, start and end.{remedy}",
        ),
        (
            "start_not_a_number",
            ["chr1\t1.5\t20"],
            f"{bed}, line 1: the start is '1.5', and it must be a whole number of 0 or"
            f" more.{remedy}",
        ),
        (
            "end_not_a_number",
            ["chr1\t10\t"],
            f"{bed}, line 1: the end is '', and it must be a whole number of 0 or more.{remedy}",
        ),
        (
            "end_before_start",
            ["chr1\t10\t9"],
            f"{bed}, line 1: the end, 9, is not after the start, 10, and a region holds at"
            f" least one base.{remedy}",
        ),
        (
            "contig_not_in_reference",
            ["chrX\t10\t20"],
            f"{bed}, line 1: contig chrX is not in the reference. Use the reference the BED was"
            " made for, or remove the lines of chrX, and run again.",
        ),
        (
            "end_past_contig",
            ["chr1\t10\t201"],
            f"{bed}, line 1: the end, 201, is past the end of chr1, which has 200 bases in the"
            " reference. Use the reference the BED was made for, or correct the line, and run"
            " again.",
        ),
        (
            "unknown_type",
            ["chr1\t10\t20\tgeneric"],
            f"{bed}, line 1: the region type is 'generic', and it must be str, str-bundle or"
            " satellite. Make the BED with `pvc_eval str-regions`, or correct the line, and run"
            " again.",
        ),
        (
            "pvc_bed",
            ["chr1\t10\t20\tssr_locus\tCA\t2\t5\t0.9\t."],
            f"{bed}, line 1: the region type is 'ssr_locus', and it must be str, str-bundle or"
            " satellite. Make the BED with `pvc_eval str-regions`, or correct the line, and run"
            " again.",
        ),
        (
            "str_without_motif",
            ["chr1\t10\t20\tstr"],
            f"{bed}, line 1: an STR segment, of type str, needs its motif in the fifth column:"
            " 1 to 6 bases of A, C, G and T, in upper case. Add it, and run again.",
        ),
        (
            "str_motif_dot",
            ["chr1\t10\t20\tstr\t."],
            f"{bed}, line 1: the motif of an STR segment is '.', and it must be 1 to 6 bases of"
            " A, C, G and T, in upper case. Correct it, and run again.",
        ),
        (
            "str_motif_lower_case",
            ["chr1\t10\t20\tstr\tca"],
            f"{bed}, line 1: the motif of an STR segment is 'ca', and it must be 1 to 6 bases of"
            " A, C, G and T, in upper case. Correct it, and run again.",
        ),
        (
            "str_motif_seven_bases",
            ["chr1\t10\t20\tstr\tACGTACG"],
            f"{bed}, line 1: the motif of an STR segment is 'ACGTACG', and it must be 1 to 6"
            " bases of A, C, G and T, in upper case. Correct it, and run again.",
        ),
        (
            "str_motif_n",
            ["chr1\t10\t20\tstr\tCN"],
            f"{bed}, line 1: the motif of an STR segment is 'CN', and it must be 1 to 6 bases of"
            " A, C, G and T, in upper case. Correct it, and run again.",
        ),
        (
            "satellite_with_motif",
            ["chr1\t10\t20\tsatellite\tCA"],
            f"{bed}, line 1: a region of type satellite has no motif, and its fifth column is"
            " 'CA'. Put '.' there, and run again.",
        ),
        (
            "bundle_with_motif",
            ["chr1\t10\t20\tstr-bundle\t"],
            f"{bed}, line 1: a region of type str-bundle has no motif, and its fifth column is"
            " ''. Put '.' there, and run again.",
        ),
        (
            "contig_lines_apart",
            ["chr1\t10\t20", "chr2\t10\t20", "chr1\t30\t40"],
            f"{bed}, line 3: the lines of chr1 are not all together; this one comes after lines"
            f" of chr2.{sort}",
        ),
        (
            "contig_before_fai_order",
            ["chr2\t10\t20", "chr1\t10\t20"],
            f"{bed}, line 2: contig chr1 comes after chr2, and the reference's .fai puts chr1"
            f" before it.{sort}",
        ),
        (
            "start_before_previous_start",
            ["chr1\t30\t40", "chr1\t10\t20"],
            f"{bed}, line 2: the start, 10, is before the start, 30, of line 1 on chr1.{sort}",
        ),
        (
            "overlap",
            ["chr1\t10\t20", "chr1\t19\t30"],
            f"{bed}, line 2: the region chr1 19 30 overlaps the one of line 1, chr1 10 20, and"
            " STR regions must not overlap. Merge them with `bedtools merge -i"
            f" {bed}`, which keeps three columns, and run again.",
        ),
        (
            "overlap_after_a_comment",
            ["chr1\t10\t20", "# between", "chr1\t15\t30"],
            f"{bed}, line 3: the region chr1 15 30 overlaps the one of line 1, chr1 10 20, and"
            " STR regions must not overlap. Merge them with `bedtools merge -i"
            f" {bed}`, which keeps three columns, and run again.",
        ),
    ]


_STR_REGIONS_ERROR_IDS = [case[0] for case in _make_str_regions_errors(Path("str.bed"))]


@pytest.mark.parametrize("case_id", _STR_REGIONS_ERROR_IDS)
def test_each_refusal_of_read_str_regions_names_the_line(case_id: str, tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT" * 50, "chr2": "ACGT" * 50})
    bed = tmp_path / "str.bed"
    (lines, message) = next(
        (lines, message)
        for name, lines, message in _make_str_regions_errors(bed)
        if name == case_id
    )
    _write_lines(bed, lines)
    with IndexedFasta(fasta) as reference, pytest.raises(PvcEvalInputError) as raised:
        list(read_str_regions(bed, reference))
    assert str(raised.value) == message


def test_regions_that_abut_are_read(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT" * 50})
    bed = _write_lines(tmp_path / "str.bed", ["chr1\t0\t20", "chr1\t20\t200"])
    with IndexedFasta(fasta) as reference:
        assert list(read_str_regions(bed, reference)) == [
            GenomeSegment("chr1", 0, 20, None, None),
            GenomeSegment("chr1", 20, 200, None, None),
        ]


def test_a_compressed_bed_is_refused(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT" * 50})
    bed = tmp_path / "str.bed.gz"
    bed.write_bytes(b"\x1f\x8b\x08\x00\xff\xfe")
    with IndexedFasta(fasta) as reference, pytest.raises(PvcEvalInputError) as raised:
        list(read_str_regions(bed, reference))
    assert str(raised.value) == (
        f"{bed} is not UTF-8 text: a byte at or after line 1 cannot be read as UTF-8. If it is"
        f" compressed, decompress it with `gunzip -c {bed} > file.bed`, give file.bed instead,"
        " and run again."
    )


def test_a_missing_bed_is_refused(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT" * 50})
    bed = tmp_path / "str.bed"
    with IndexedFasta(fasta) as reference, pytest.raises(PvcEvalInputError) as raised:
        list(read_str_regions(bed, reference))
    assert str(raised.value) == f"{bed} does not exist. Check its path and run again."
