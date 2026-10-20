"""The `pvc_eval str-regions` command and read_str_regions (doc/specs/str_finder.md,
"read_str_regions" and "The str-regions command")."""

import gzip
import os
import re
import subprocess
import sys
import sysconfig
import threading
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import CodeType

import pytest

from pvc_eval.cli import main
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import (
    GenomeSegment,
    RegionType,
    StrFinderParameters,
    read_str_regions,
    segment_contig,
    write_str_regions,
)

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
    """The exit status of `pvc_eval` with these arguments, and what it wrote to stderr. It
    writes nothing to stdout."""
    try:
        main(argv)
    except SystemExit as exit_:
        status = exit_.code
    else:
        status = 0
    assert isinstance(status, int)
    captured = capsys.readouterr()
    assert captured.out == ""
    return status, captured.err


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


# ---- --threads --------------------------------------------------------------------------

# A first contig far longer than the others, 302,000 bases against at most 302, so that with
# several threads the contigs after it are segmented and finished long before it is: rows
# written in the order the contigs finish would put its rows last. It is four_a 2,000 times
# over, so its rows are those of `twice` repeated, one (CA) repeat every 151 bases.
_LONG_COPIES = 2_000
_LONG_CONTIGS = {"long": _FOUR_A * _LONG_COPIES, **_WORKED_CONTIGS}
_LONG_ROWS = [
    f"long\t{58 + 151 * copy}\t{80 + 151 * copy}\tstr\tCA" for copy in range(_LONG_COPIES)
]


def _run_threads(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], contigs: dict[str, str], threads: str
) -> tuple[int, str, Path]:
    """The exit status and stderr of `str-regions --threads threads` over these contigs, and
    the path of the BED it writes."""
    fasta = _write_fasta(tmp_path, contigs)
    output = tmp_path / f"str_{threads}.bed"
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(output), "--threads", threads),
        ],
        capsys,
    )
    return status, stderr, output


def test_threads_write_the_rows_in_fai_order_whichever_contig_finishes_first(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Seven contigs on three threads.
    beds = {}
    for threads in ["1", "3"]:
        status, stderr, output = _run_threads(tmp_path, capsys, _LONG_CONTIGS, threads)
        assert (status, stderr) == (0, "")
        beds[threads] = output.read_bytes()
    assert beds["3"] == beds["1"]
    assert beds["1"].decode().splitlines() == _DEFAULT_HEADER + _LONG_ROWS + _WORKED_ROWS


@pytest.mark.parametrize("threads", ["0", "257"])
def test_threads_outside_1_to_256_stop_the_run_naming_the_option(
    threads: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    status, stderr, output = _run_threads(tmp_path, capsys, {"four_a": _FOUR_A}, threads)
    assert (status, stderr) == (
        1,
        f"--threads is {threads}; it must be a whole number from 1 to 256. Correct the"
        " option, and run again.\n",
    )
    assert not output.exists()
    assert not Path(f"{output}.tmp").exists()


def test_threads_that_are_not_a_number_get_a_usage_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    status, stderr, _ = _run_threads(tmp_path, capsys, {"four_a": _FOUR_A}, "four")
    assert status == 2
    assert stderr.startswith("usage: pvc_eval str-regions")
    assert "argument --threads: invalid int value: 'four'" in stderr


def test_an_error_in_a_later_contig_under_threads_stops_the_run_with_its_message(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Two contigs with a byte that is not a letter, both after the long one: the run stops
    # with the first one's message, as with one thread, whichever finishes first.
    contigs = {"long": _LONG_CONTIGS["long"], "chr2": "AC-GT", "four_a": _FOUR_A, "chr4": "-"}
    output = tmp_path / "str_3.bed"
    output.write_text("earlier\n", encoding="utf-8")
    status, stderr, output = _run_threads(tmp_path, capsys, contigs, "3")
    fasta = tmp_path / "ref.fa"
    assert (status, stderr) == (
        1,
        f"{fasta}: chr2:3 is not a letter where its .fai, line 2, puts a base. Make the .fai"
        f" again with `samtools faidx {fasta}`, check that the FASTA holds only letters"
        " between its header lines, and run again.\n",
    )
    assert output.read_text(encoding="utf-8") == "earlier\n"
    assert not Path(f"{output}.tmp").exists()


@pytest.mark.parametrize("num_threads", [0, 257, True, 2.0], ids=["0", "257", "true", "float"])
def test_write_str_regions_refuses_threads_before_reading_the_regions(
    num_threads: object, tmp_path: Path
) -> None:
    # The --regions BED does not exist, which would stop the run if it were read first.
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    output = tmp_path / "str.bed"
    with (
        IndexedFasta(fasta) as reference,
        pytest.raises(
            ValueError,
            match=re.escape(
                f"num_threads is {num_threads!r}; it must be a whole number from 1 to 256"
            ),
        ),
    ):
        write_str_regions(
            reference,
            output,
            StrFinderParameters(),
            tmp_path / "missing.bed",
            num_threads,  # ty: ignore[invalid-argument-type]  # the wrong types are the test
        )
    assert not Path(f"{output}.tmp").exists()


# ---- Which thread segments which contig -------------------------------------------------

# A tool id of sys.monitoring that nothing else in the tests takes: 0, 1, 2 and 5 are those
# of a debugger, of coverage, of a profiler and of an optimizer.
_MONITORING_TOOL = 4

# The long contig of _LONG_CONTIGS and then twelve short ones, so that with two threads the
# short ones would all be started long before the long one returns, were nothing to hold
# them back.
_LONG_THEN_SHORT = {"long": _LONG_CONTIGS["long"]} | {
    f"short{index}": _FOUR_A for index in range(12)
}

# A call of segment_contig as sys.monitoring sees it: "start" when it starts or "return"
# when it hands back its segments, the contig, and the ident of the thread it runs on.
type _SegmentEvent = tuple[str, str, int]


@pytest.fixture
def segment_contig_events() -> Iterator[list[_SegmentEvent]]:
    """The calls of segment_contig while the test runs, in the order they were recorded.

    sys.monitoring calls back on each start and return of segment_contig's code, and replaces
    nothing, so the code that runs is pvc_eval's own. The tool id and the events are freed
    when the test ends, failed or not.
    """
    # list.append from several threads keeps every item on the free-threaded build, and an
    # event is appended after it happened, so a start recorded before a return happened
    # before it.
    events: list[_SegmentEvent] = []

    def record_start(code: CodeType, offset: int) -> None:  # noqa: ARG001 -- the signature sys.monitoring calls
        # Frame 1 is segment_contig's, which called back.
        events.append(("start", sys._getframe(1).f_locals["contig"], threading.get_ident()))

    def record_return(code: CodeType, offset: int, value: object) -> None:  # noqa: ARG001 -- the signature sys.monitoring calls
        events.append(("return", sys._getframe(1).f_locals["contig"], threading.get_ident()))

    monitoring = sys.monitoring
    monitoring.use_tool_id(_MONITORING_TOOL, "pvc_eval tests")
    try:
        monitoring.register_callback(_MONITORING_TOOL, monitoring.events.PY_START, record_start)
        monitoring.register_callback(_MONITORING_TOOL, monitoring.events.PY_RETURN, record_return)
        monitoring.set_local_events(
            _MONITORING_TOOL,
            segment_contig.__code__,
            monitoring.events.PY_START | monitoring.events.PY_RETURN,
        )
        yield events
    finally:
        monitoring.set_local_events(
            _MONITORING_TOOL, segment_contig.__code__, monitoring.events.NO_EVENTS
        )
        monitoring.free_tool_id(_MONITORING_TOOL)


def _run_long_then_short(tmp_path: Path, capsys: pytest.CaptureFixture[str], threads: str) -> None:
    status, stderr, _ = _run_threads(tmp_path, capsys, _LONG_THEN_SHORT, threads)
    assert (status, stderr) == (0, "")


def test_one_thread_segments_every_contig_in_the_calling_thread(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], segment_contig_events: list[_SegmentEvent]
) -> None:
    _run_long_then_short(tmp_path, capsys, "1")
    # A start and a return for each of the 13 contigs.
    assert len(segment_contig_events) == 26
    assert {thread for _, _, thread in segment_contig_events} == {threading.get_ident()}


def test_two_threads_segment_no_contig_in_the_calling_thread_and_use_at_most_two(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], segment_contig_events: list[_SegmentEvent]
) -> None:
    _run_long_then_short(tmp_path, capsys, "2")
    assert len(segment_contig_events) == 26
    threads = {thread for _, _, thread in segment_contig_events}
    assert threading.get_ident() not in threads
    assert len(threads) <= 2


def test_two_threads_start_at_most_five_contigs_before_the_first_returns(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], segment_contig_events: list[_SegmentEvent]
) -> None:
    # Executor.map with buffersize=2 * 2 submits the long contig and three short ones, and one
    # more when the first result is asked for, before it waits on the long one: 1 + 2N. With
    # no buffersize, the twelve short ones all finish on the second thread while the long one
    # runs.
    _run_long_then_short(tmp_path, capsys, "2")
    long_return = next(
        index
        for index, (kind, contig, _) in enumerate(segment_contig_events)
        if (kind, contig) == ("return", "long")
    )
    started = [contig for kind, contig, _ in segment_contig_events[:long_return] if kind == "start"]
    assert len(started) <= 5, started


# ---- --regions --------------------------------------------------------------------------


def test_a_region_that_cuts_an_str_region_writes_it_whole(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # In an order that is not the .fai's, which the rows keep.
    regions = _write_lines(tmp_path / "regions.bed", ["twice\t220\t300", "four_a\t70\t75"])
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
    """Each refusal of a --regions BED: its id, its lines, and the message. All but the end
    past the contig's length are pvc_eval.bed.parse_bed's."""
    number = (
        " A BED line gives its start and its end as whole numbers of 0 or more; correct the"
        " line and run again."
    )
    return [
        (
            "contig_not_in_reference",
            ["four_a\t0\t10", "chrX\t0\t10"],
            f"{regions}, line 2: the contig 'chrX' is not in the reference. Keep only the"
            " reference's contigs in the BED, for instance with"
            " `awk -F'\\t' 'NR == FNR {keep[$1]; next} $1 in keep'"
            f" ref.fa.fai {regions} > kept.bed`, where ref.fa.fai is the .fai of the"
            " --reference FASTA, and run again with kept.bed.",
        ),
        (
            "start_not_a_number",
            ["four_a\tten\t20"],
            f"{regions}, line 1: the start, 'ten', is not a whole number.{number}",
        ),
        (
            "start_negative",
            ["four_a\t-1\t20"],
            f"{regions}, line 1: the start, '-1', is not a whole number.{number}",
        ),
        (
            "end_not_a_number",
            ["four_a\t10\t2e1"],
            f"{regions}, line 1: the end, '2e1', is not a whole number.{number}",
        ),
        (
            "end_at_start",
            ["four_a\t10\t10"],
            f"{regions}, line 1: the end, 10, is not after the start, 10. A BED line covers the"
            " positions from its start up to its end, 0-based and half-open, so its end is"
            " greater than its start; correct the line and run again.",
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
            f"{regions}, line 1: has 2 columns, and a BED line needs at least three, separated"
            " by tabs: the contig, the start and the end. Correct the line or remove it, and"
            " run again.",
        ),
        (
            "one_column",
            ["four_a"],
            f"{regions}, line 1: has 1 column, and a BED line needs at least three, separated"
            " by tabs: the contig, the start and the end. Correct the line or remove it, and"
            " run again.",
        ),
        (
            "columns_separated_by_spaces",
            ["four_a 10 20"],
            f"{regions}, line 1: has its columns separated by spaces, and a BED separates them"
            " by tabs. Convert it with `awk -v OFS='\\t' '{$1 = $1} 1'"
            f" {regions} > tabs.bed` and run again with tabs.bed.",
        ),
        (
            "no_region",
            ["# only a comment", "track name=x", ""],
            f"{regions} holds no BED line: it is empty, or holds only blank lines, comments, and"
            " track and browser lines. Give a BED with at least one line of contig, start and"
            " end, separated by tabs, and run again.",
        ),
    ]


def test_a_gzipped_regions_bed_is_read(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    regions = tmp_path / "regions.bed.gz"
    regions.write_bytes(gzip.compress(b"twice\t220\t300\n"))
    lines = _run_str_regions(tmp_path, capsys, "--regions", str(regions))
    assert lines == [*_DEFAULT_HEADER, "twice\t209\t231\tstr\tCA"]


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


def test_a_missing_fai_stops_the_run(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    Path(f"{fasta}.fai").unlink()
    status, stderr = _run(
        ["str-regions", "--reference", str(fasta), "--output", str(tmp_path / "str.bed")], capsys
    )
    assert (status, stderr) == (
        1,
        f"{fasta} has no index {fasta}.fai. Make it with `samtools faidx {fasta}` and run again.\n",
    )


def test_a_regions_bed_that_is_a_directory_stops_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    regions = tmp_path / "regions.bed"
    regions.mkdir()
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(tmp_path / "str.bed"), "--regions", str(regions)),
        ],
        capsys,
    )
    assert (status, stderr) == (
        1,
        f"{regions} is a directory, not a BED file. Check its path and run again.\n",
    )


@pytest.mark.skipif(os.geteuid() == 0, reason="root reads a file without read permission")
def test_a_regions_bed_that_cannot_be_read_stops_the_run(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    regions = _write_lines(tmp_path / "regions.bed", ["four_a\t0\t10"])
    regions.chmod(0)
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(tmp_path / "str.bed"), "--regions", str(regions)),
        ],
        capsys,
    )
    assert (status, stderr) == (
        1,
        f"{regions} cannot be read: permission denied. Allow reading it and run again.\n",
    )


def test_copy_floors_that_are_not_numbers_get_a_usage_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # A command line that cannot be parsed keeps argparse's usage message and exit status 2.
    fasta = _write_fasta(tmp_path, {"four_a": _FOUR_A})
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta), "--output", str(tmp_path / "str.bed")),
            *("--min-copies", "5,a,4,4,4,3"),
        ],
        capsys,
    )
    assert status == 2
    assert stderr.startswith("usage: pvc_eval str-regions ")
    assert stderr.splitlines()[-1] == (
        "pvc_eval str-regions: error: argument --min-copies: '5,a,4,4,4,3' is not whole numbers"
        " separated by commas, such as 5,5,4,4,4,3"
    )


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
            "--min-purity is 2.0; it must be a number from 0 to 1",
        ),
        (
            ["--min-period", "4", "--max-period", "3"],
            "--min-period (4) must not exceed --max-period (3)",
        ),
        (
            ["--scan-min-period", "7"],
            "--scan-min-period is 7; it must be 6 or less",
        ),
        (
            ["--min-period", "7"],
            "--min-period is 7; it must be 6 or less",
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
        "scan_min_period_seven",
        "min_period_seven",
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
        "--threads": ["4"],
    }


def test_help_says_what_the_bed_holds_and_what_its_words_mean(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit):
        main(["str-regions", "--help"])
    help_text = " ".join(capsys.readouterr().out.split())
    for sentence in [
        "A repeat is a stretch of the reference made of copies of a short sequence, its motif;"
        " the motif's length, 1 to 6 bases, is the repeat's period.",
        "An STR segment is a repeat with at least --bundle-distance bases between it and any"
        " other, trimmed to whole copies of its motif.",
        "An STR bundle is a stretch of repeats with fewer than --bundle-distance bases between"
        " each and the next.",
        "A satellite is a stretch of repeats longer than --max-str-len bases, or an STR bundle"
        " longer than --max-str-bundle-len bases.",
        "The scan looks for repeats of periods --scan-min-period to --scan-max-period, and only"
        " those of periods --min-period to --max-period become STR segments or members of STR"
        " bundles; the others still make satellites.",
        "one row per STR region, in the order of the reference's .fai and then by start, with"
        " five columns: the contig; the start and the end, 0-based and half-open; the region"
        " type, str, str-bundle or satellite; and the motif of an STR segment, or . for the"
        " other two types.",
        "--max-str-len N a stretch of repeats longer than N bases is a satellite; at least"
        " --bundle-distance (default 100)",
        "--threads N how many contigs are segmented at once, each on a thread of its own, from"
        " 1 to 256. The peak memory grows with N: a contig being segmented holds 4.6 to 5.3 MB"
        " per Mb of its length, 1.15 GB for GRCh38's chr1. With more than one thread, Ctrl-C"
        " stops the run once the contigs being segmented have finished, which can take"
        " minutes (default 4)",
    ]:
        assert sentence in help_text
    assert "by half" not in help_text


# ---- read_str_regions -------------------------------------------------------------------


def test_the_bed_str_regions_writes_reads_back_to_the_same_str_regions(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _run_str_regions(tmp_path, capsys)
    with IndexedFasta(tmp_path / "ref.fa") as reference:
        assert list(read_str_regions(tmp_path / "str.bed", reference)) == _WORKED_SEGMENTS


def test_a_contig_named_like_a_header_word_is_read_and_written(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # A line is a header only when its first field is track or browser, not trackA.
    fasta = _write_fasta(tmp_path, {"trackA": _FOUR_A, "browserB": _FOUR_A})
    regions = _write_lines(tmp_path / "regions.bed", ["trackA\t70\t75", "browserB\t0\t151"])
    output = tmp_path / "str.bed"
    status, stderr = _run(
        [
            *("str-regions", "--reference", str(fasta)),
            *("--output", str(output), "--regions", str(regions)),
        ],
        capsys,
    )
    assert (status, stderr) == (0, "")
    assert output.read_text(encoding="utf-8").splitlines()[-2:] == [
        "trackA\t58\t80\tstr\tCA",
        "browserB\t58\t80\tstr\tCA",
    ]
    with IndexedFasta(fasta) as reference:
        assert list(read_str_regions(output, reference)) == [
            GenomeSegment("trackA", 58, 80, RegionType.STR, "CA"),
            GenomeSegment("browserB", 58, 80, RegionType.STR, "CA"),
        ]


def test_blank_lines_trailing_tabs_crlf_and_a_byte_order_mark_are_read(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path, {"chr1": "ACGT" * 50})
    bed = tmp_path / "str.bed"
    bed.write_bytes(
        "\ufeffchr1\t0\t10\r\n"
        " \n"
        "\t\n"
        " \t \r\n"
        "track\tname=x\n"
        "browser position chr1:1-100\n"
        "chr1\t10\t20\t\n"
        "chr1\t20\t30\tstr-bundle\t\n"
        "chr1\t30\t40\tstr\tCA\t\t\r\n"
        "chr1\t40\t50\tstr\tCA\r\n".encode()
    )
    with IndexedFasta(fasta) as reference:
        assert list(read_str_regions(bed, reference)) == [
            GenomeSegment("chr1", 0, 10, None, None),
            GenomeSegment("chr1", 10, 20, None, None),
            GenomeSegment("chr1", 20, 30, RegionType.STR_BUNDLE, None),
            GenomeSegment("chr1", 30, 40, RegionType.STR, "CA"),
            GenomeSegment("chr1", 40, 50, RegionType.STR, "CA"),
        ]


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
            "start_a_digit_not_ascii",
            ["chr1\t\u0661\t20"],
            f"{bed}, line 1: the start is '\u0661', and it must be a whole number of 0 or"
            f" more.{remedy}",
        ),
        (
            "end_not_a_number",
            ["chr1\t10\tx"],
            f"{bed}, line 1: the end is 'x', and it must be a whole number of 0 or more.{remedy}",
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
            "str_with_an_empty_fifth_column_at_the_end",
            ["chr1\t10\t20\tstr\t"],
            f"{bed}, line 1: an STR segment, of type str, needs its motif in the fifth column:"
            " 1 to 6 bases of A, C, G and T, in upper case. Add it, and run again.",
        ),
        (
            "str_motif_empty",
            ["chr1\t10\t20\tstr\t\t0.9"],
            f"{bed}, line 1: the motif of an STR segment is '', and it must be 1 to 6 bases of"
            " A, C, G and T, in upper case. Correct it, and run again.",
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
            ["chr1\t10\t20\tstr-bundle\tN"],
            f"{bed}, line 1: a region of type str-bundle has no motif, and its fifth column is"
            " 'N'. Put '.' there, and run again.",
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
            f" STR regions must not overlap. Merge them with `bedtools merge -i {bed}`, and run"
            " again. The merged BED has three columns: its regions still widen the loci, but"
            " they have no region type, which `compare --strata pvc` reads.",
        ),
        (
            "overlap_after_a_comment",
            ["chr1\t10\t20", "# between", "chr1\t15\t30"],
            f"{bed}, line 3: the region chr1 15 30 overlaps the one of line 1, chr1 10 20, and"
            f" STR regions must not overlap. Merge them with `bedtools merge -i {bed}`, and run"
            " again. The merged BED has three columns: its regions still widen the loci, but"
            " they have no region type, which `compare --strata pvc` reads.",
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
