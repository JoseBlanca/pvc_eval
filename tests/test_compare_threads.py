"""compare() with its threads (doc/specs/compare.md, "Which thread does what"): on the
four-contig case of tests/compare_case.py, its rows and its first error with 1, 2 and 4 finder
threads are those with none, whatever the indexes beside the VCFs; the contigs it scans; and,
on a file of 50,000 VCF records, the bound of a reader's buffer and the stop.

Every test runs compare() in a thread of its own and waits for it at most _TIMEOUT seconds, so
that a hang fails the test, and fails when an exception reaches threading.excepthook.
"""

import gzip
import logging
import random
import struct
import threading
from collections.abc import Callable, Generator, Iterator, Sequence
from concurrent.futures import Future
from contextlib import closing
from dataclasses import replace
from pathlib import Path
from queue import SimpleQueue

import pytest
from tests.compare_case import CONTIGS, LINES_A, LINES_B, CompareCase, write_compare_case
from tests.test_compare import (
    _CONTIGS_WITH_RECORDS,
    _get_scan_messages,
    _make_config,
    _run_nine_steps,
)

from pvc_eval import compare as compare_module
from pvc_eval.compare import PVC_STRATA, CompareConfig, compare
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.outcome import OutcomeRow
from pvc_eval.segments import GenomeSegment, segment_contig
from pvc_eval.vcf import VcfRecord, parse_vcf

_LOGGER = "pvc_eval.compare"
# Seconds a test waits for compare() or for its threads before it fails as a hang.
_TIMEOUT = 60
_FINDER_THREADS = [1, 2, 4]
# doc/specs/compare.md, "Which thread does what": a reader's buffer holds at most 10,000
# records, handed over in batches of 64, so 156 batches, 9,984 records.
_BATCH_RECORDS = 64
_BUFFER_BATCHES = 10_000 // _BATCH_RECORDS

type _Run = tuple[list[OutcomeRow], Exception | None]


@pytest.fixture(autouse=True)
def _fail_on_thread_exception() -> Iterator[None]:
    """Fails the test when an exception reaches threading.excepthook, after compare()'s threads
    have ended, which they must within _TIMEOUT seconds."""
    hooked: SimpleQueue[threading.ExceptHookArgs] = SimpleQueue()
    threads_before = set(threading.enumerate())
    previous_hook = threading.excepthook
    threading.excepthook = hooked.put
    try:
        yield
        # Only the threads of this test, so that a hang fails it and not the tests after it.
        alive = _join_threads(threads_before, "pvc_eval", _TIMEOUT)
    finally:
        threading.excepthook = previous_hook
    assert not alive, f"threads of compare() still alive after {_TIMEOUT} s: {alive}"
    assert hooked.empty(), f"an exception reached threading.excepthook: {hooked.get()}"


def _get_new_threads(threads_before: set[threading.Thread], prefix: str) -> list[threading.Thread]:
    """The threads alive whose names start with prefix, less threads_before."""
    return [
        thread
        for thread in threading.enumerate()
        if thread.name.startswith(prefix) and thread not in threads_before
    ]


def _join_threads(threads_before: set[threading.Thread], prefix: str, timeout: float) -> list[str]:
    """Waits at most timeout seconds in all for the threads whose names start with prefix,
    less threads_before, and returns the names of those still alive."""
    deadline = threading.Event()
    timer = threading.Timer(timeout, deadline.set)
    timer.start()
    try:
        for thread in _get_new_threads(threads_before, prefix):
            while thread.is_alive() and not deadline.is_set():
                thread.join(0.05)
    finally:
        timer.cancel()
    return [thread.name for thread in _get_new_threads(threads_before, prefix) if thread.is_alive()]


def _wait_until(is_done: Callable[[], bool], what: str) -> None:
    """Returns once is_done() is true, failing the test when it is not within _TIMEOUT
    seconds."""
    deadline = threading.Event()
    timer = threading.Timer(_TIMEOUT, deadline.set)
    timer.start()
    try:
        while not is_done():
            if deadline.wait(0.02):
                pytest.fail(f"{what} did not happen within {_TIMEOUT} s")
    finally:
        timer.cancel()


def _wait_until_steady(get_value: Callable[[], int], minimum: int, what: str) -> int:
    """get_value() once it is minimum or more and has not changed for 0.5 s, failing the
    test when that does not happen within _TIMEOUT seconds."""
    values = [-1]

    def is_steady() -> bool:
        value = get_value()
        if value >= minimum and value == values[-1]:
            return True
        values.append(value)
        threading.Event().wait(0.5)
        return False

    _wait_until(is_steady, what)
    return values[-1]


def _start_call[T](function: Callable[[], T]) -> Future[T]:
    """The future of function's result, or of its exception, called in a thread of the test."""
    future: Future[T] = Future()

    def call() -> None:
        try:
            future.set_result(function())
        except BaseException as error:  # noqa: BLE001 -- raised again by future.result
            future.set_exception(error)

    threading.Thread(target=call, name="test caller", daemon=True).start()
    return future


def _get_result[T](future: Future[T]) -> T:
    """The future's result, failing the test when none comes within _TIMEOUT seconds."""
    try:
        return future.result(timeout=_TIMEOUT)
    except TimeoutError:
        pytest.fail(f"no result within {_TIMEOUT} s: a thread of compare() hangs")


def _call_with_timeout[T](function: Callable[[], T]) -> T:
    """function's result, or its exception, from a thread of the test, failing the test when
    neither comes within _TIMEOUT seconds."""
    return _get_result(_start_call(function))


def _run_compare(case: CompareCase, config: CompareConfig, finder_threads: int) -> _Run:
    """compare()'s rows, and the error that stopped them, from a thread of the test."""
    return _call_with_timeout(lambda: _collect_run(case, config, finder_threads))


def _collect_run(case: CompareCase, config: CompareConfig, finder_threads: int) -> _Run:
    """compare()'s rows, and the error that stopped them."""
    rows: list[OutcomeRow] = []
    try:
        rows.extend(compare(case.a, case.b, case.reference, config, finder_threads=finder_threads))
    except (PvcEvalInputError, ValueError) as error:
        return rows, error
    return rows, None


def _check_same_run(run: _Run, expected: _Run) -> None:
    rows, error = run
    expected_rows, expected_error = expected
    assert rows == expected_rows
    assert type(error) is type(expected_error)
    assert str(error) == str(expected_error)


def _write_index(vcf_path: Path, suffix: str, contigs: Sequence[str]) -> None:
    """Writes beside a VCF a .tbi or a .csi whose header lists contigs, as htslib's tabix and
    CSI formats lay it out (src/pvc_eval/_vcf_index.py), with no bins after it."""
    names = b"".join(contig.encode() + b"\0" for contig in contigs)
    # format 2 (VCF), col_seq 1, col_beg 2, col_end 0, meta "#", skip 0, then l_nm.
    tabix = struct.pack("<7i", 2, 1, 2, 0, ord("#"), 0, len(names)) + names
    if suffix == ".tbi":
        header = b"TBI\x01" + struct.pack("<i", len(contigs)) + tabix
    else:
        # min_shift 14 and depth 5, those of bcftools 1.24's .csi.
        aux = struct.pack("<3i", 14, 5, len(tabix)) + tabix
        header = b"CSI\x01" + aux + struct.pack("<i", len(contigs))
    with gzip.open(vcf_path.with_name(vcf_path.name + suffix), "wb") as index_file:
        index_file.write(header)


def _gzip_vcf(vcf_path: Path) -> Path:
    """The VCF gzipped beside itself, as <name>.gz."""
    gzipped = vcf_path.with_name(vcf_path.name + ".gz")
    gzipped.write_bytes(gzip.compress(vcf_path.read_bytes()))
    return gzipped


def _index_case(
    case: CompareCase,
    *,
    contigs_a: Sequence[str] | None = None,
    contigs_b: Sequence[str] | None = None,
) -> CompareCase:
    """The case with file A gzipped and indexed by a .tbi that lists contigs_a, and file B
    by a .csi that lists contigs_b, each only when its contigs are given."""
    if contigs_a is not None:
        case = replace(case, a=_gzip_vcf(case.a))
        _write_index(case.a, ".tbi", contigs_a)
    if contigs_b is not None:
        case = replace(case, b=_gzip_vcf(case.b))
        _write_index(case.b, ".csi", contigs_b)
    return case


# The contigs with records of the four-contig case: chr1 and chr3 in file A, chr2 and chr3 in B.
_INDEXES = {
    "neither indexed": {},
    "A indexed": {"contigs_a": ["chr1", "chr3"]},
    "B indexed": {"contigs_b": ["chr2", "chr3"]},
    "both indexed": {"contigs_a": ["chr1", "chr3"], "contigs_b": ["chr2", "chr3"]},
}


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
@pytest.mark.parametrize("option", ["no option", "--regions chr3 0 300", "--str-regions"])
@pytest.mark.parametrize("indexes", list(_INDEXES))
def test_rows_with_threads_are_those_of_the_nine_steps(
    tmp_path: Path, indexes: str, option: str, finder_threads: int
) -> None:
    case = _index_case(write_compare_case(tmp_path), **_INDEXES[indexes])
    config = _make_config(option, tmp_path)
    expected = _run_nine_steps(case, config, _CONTIGS_WITH_RECORDS)
    assert len(expected[0]) == {"--regions chr3 0 300": 2}.get(option, 7)
    _check_same_run(_run_compare(case, config, 0), expected)
    _check_same_run(_run_compare(case, config, finder_threads), expected)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
@pytest.mark.parametrize(
    "index_contigs",
    [["chr1", "chr3", "chr4"], ["chr3"], ["chrX", "chr1", "chr3"]],
    ids=["lists chr4, no record on it", "leaves out chr1", "lists chrX, the reference lacks it"],
)
def test_index_out_of_date_gives_the_same_rows(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    index_contigs: list[str],
    finder_threads: int,
) -> None:
    case = write_compare_case(tmp_path)
    expected = _run_compare(case, CompareConfig(), 0)
    indexed_case = _index_case(case, contigs_a=index_contigs)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    _check_same_run(_run_compare(indexed_case, CompareConfig(), finder_threads), expected)
    # A contig with records that the index leaves out is still scanned.
    scanned = {message.split(":")[0] for message in _get_scan_messages(caplog)}
    assert scanned >= set(_CONTIGS_WITH_RECORDS)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
def test_index_that_leaves_out_contig_2_gives_the_same_rows_and_scans_it(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, finder_threads: int
) -> None:
    case = write_compare_case(tmp_path)
    expected = _run_compare(case, CompareConfig(), 0)
    indexed_case = _index_case(case, contigs_a=["chr1", "chr3"], contigs_b=["chr3"])
    caplog.set_level(logging.INFO, logger=_LOGGER)
    _check_same_run(_run_compare(indexed_case, CompareConfig(), finder_threads), expected)
    assert "chr2: finding STR regions in 300 bases." in _get_scan_messages(caplog)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
def test_without_index_contigs_with_records_are_scanned_and_contig_4_is_not(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, finder_threads: int
) -> None:
    case = write_compare_case(tmp_path)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    _run_compare(case, CompareConfig(), finder_threads)
    assert sorted(_get_scan_messages(caplog)) == [
        "chr1: finding STR regions in 300 bases.",
        "chr2: finding STR regions in 300 bases.",
        "chr3: finding STR regions in 300 bases.",
    ]


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
def test_regions_naming_contig_3_scan_contig_3_alone_whatever_the_indexes(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, finder_threads: int
) -> None:
    all_contigs = ["chr1", "chr2", "chr3", "chr4"]
    case = _index_case(write_compare_case(tmp_path), contigs_a=all_contigs, contigs_b=all_contigs)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    config = _make_config("--regions chr3 0 300", tmp_path)
    _run_compare(case, config, finder_threads)
    assert _get_scan_messages(caplog) == ["chr3: finding STR regions in 300 bases."]


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
@pytest.mark.parametrize(
    "last_line_on_chr1",
    ["chr1 250 T C 30 0/1:30 0/0:30", "chr1 300 C T 30 0/1:30 0/0:30"],
    ids=["past the last STR region", "ending on the last base"],
)
@pytest.mark.parametrize(
    "next_line",
    ["chr3 105 C G 40 0/1:33 1/1:30", "chr2 125 A G 40 0/1:33 1/1:30"],
    ids=["then contig 3", "then contig 2"],
)
def test_record_past_last_str_region_then_next_contig_gives_rows_of_the_nine_steps(
    tmp_path: Path, last_line_on_chr1: str, next_line: str, finder_threads: int
) -> None:
    lines_a = ["chr1 50 T C 50 0/1:40 1/1:35", last_line_on_chr1, next_line]
    lines_b = ["chr1 50 T C 50 0/1:40 0/0:35"]
    case = write_compare_case(tmp_path, lines_a=lines_a, lines_b=lines_b)
    expected = _run_nine_steps(case, CompareConfig(), ["chr1", next_line.split()[0]])
    _check_same_run(_run_compare(case, CompareConfig(), finder_threads), expected)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
def test_record_on_last_contig_gives_rows_of_the_nine_steps(
    tmp_path: Path, finder_threads: int
) -> None:
    case = write_compare_case(tmp_path, lines_a=[*LINES_A, "chr4 51 C A 30 0/1:30 0/0:30"])
    expected = _run_nine_steps(case, CompareConfig(), [*_CONTIGS_WITH_RECORDS, "chr4"])
    assert expected[0][-1].contig == "chr4"
    _check_same_run(_run_compare(case, CompareConfig(), finder_threads), expected)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
@pytest.mark.parametrize("indexes", ["neither indexed", "both indexed"])
def test_vcf_out_of_order_on_contig_3_stops_after_the_same_rows(
    tmp_path: Path, indexes: str, finder_threads: int
) -> None:
    lines_a = [*LINES_A[:3], LINES_A[4], LINES_A[3]]
    case = _index_case(write_compare_case(tmp_path, lines_a=lines_a), **_INDEXES[indexes])
    expected = _run_compare(case, CompareConfig(), 0)
    assert expected[0]
    assert isinstance(expected[1], PvcEvalInputError)
    _check_same_run(_run_compare(case, CompareConfig(), finder_threads), expected)


def _write_bad_byte_case(
    directory: Path, bad_contig: str, lines_a: Sequence[str], lines_b: Sequence[str]
) -> CompareCase:
    """The case with bad_contig grown to 140,000 bases, a 0 at 139,001 (1-based), more than
    BLOCK_BASES, 65,536, from its records, so that the finder reads it and a reader does not."""
    bases = CONTIGS[bad_contig] + "N" * 139_700
    bases = bases[:139_000] + "0" + bases[139_001:]
    contigs = {**CONTIGS, bad_contig: bases}
    return write_compare_case(directory, contigs=contigs, lines_a=lines_a, lines_b=lines_b)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
@pytest.mark.parametrize(
    ("bad_contig", "lines_a", "lines_b"),
    [
        ("chr2", LINES_A, LINES_B),
        ("chr3", ["chr1 50 T C 50 0/1:40 1/1:35", "chr1 250 T C 30 0/1:30 0/0:30", LINES_A[3]], []),
        ("chr2", ["chr1 50 T C 50 0/1:40 1/1:35", "chr1 205 G A 30 0/1:30 0/0:30"],
         ["chr2 125 A G 40 0/1:33 1/1:30"]),
    ],
    ids=[
        "chr2, default records",
        "chr3, after a record past chr1's STR regions",
        "chr2 taken early",
    ],
)  # fmt: skip
def test_fasta_byte_not_a_letter_stops_the_finder_after_the_same_rows(
    tmp_path: Path,
    bad_contig: str,
    lines_a: Sequence[str],
    lines_b: Sequence[str],
    finder_threads: int,
) -> None:
    case = _write_bad_byte_case(tmp_path, bad_contig, lines_a, lines_b)
    expected = _run_compare(case, CompareConfig(), 0)
    assert isinstance(expected[1], PvcEvalInputError)
    assert f"{bad_contig}:139001 is not a letter" in str(expected[1])
    _check_same_run(_run_compare(case, CompareConfig(), finder_threads), expected)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
def test_fasta_error_on_contig_with_no_record_listed_by_an_index_is_not_raised(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, finder_threads: int
) -> None:
    case = _write_bad_byte_case(tmp_path, "chr4", LINES_A, LINES_B)
    expected_rows, expected_error = _run_compare(case, CompareConfig(), 0)
    assert expected_error is None
    indexed_case = _index_case(case, contigs_a=["chr1", "chr3", "chr4"])
    scans_ended: list[str] = []
    monkeypatch.setattr(compare_module, "segment_contig", _make_noted_segment_contig(scans_ended))
    rows = compare(
        indexed_case.a, indexed_case.b, indexed_case.reference, finder_threads=finder_threads
    )
    try:
        # Every row but the end of the stream, which stops the threads.
        taken_rows = [_call_with_timeout(lambda: next(rows)) for _ in expected_rows]
        # chr4's scan has raised before the end, so its error is held, and not raised.
        _wait_until(lambda: "chr4" in scans_ended, "the end of chr4's scan")
        assert _call_with_timeout(lambda: list(rows)) == []
    finally:
        _call_with_timeout(rows.close)
    assert taken_rows == expected_rows


def _make_noted_segment_contig(
    scans_ended: list[str],
) -> Callable[[IndexedFasta, str], Iterator[GenomeSegment]]:
    """A stand-in for segment_contig that notes in scans_ended each contig whose scan has
    ended, with its segments or with an exception."""

    def segment_contig_noted(reference: IndexedFasta, contig: str) -> Iterator[GenomeSegment]:
        try:
            return segment_contig(reference, contig)
        finally:
            scans_ended.append(contig)

    return segment_contig_noted


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
def test_strata_pvc_with_untyped_str_region_stops_after_the_same_rows(
    tmp_path: Path, finder_threads: int
) -> None:
    case = write_compare_case(tmp_path)
    str_regions = tmp_path / "str.bed"
    str_regions.write_text(
        "chr1\t100\t116\tstr\tAC\nchr1\t200\t218\tstr\tAGT\nchr2\t120\t140\tstr\tATTC\n"
        "chr3\t100\t120\n",
        encoding="utf-8",
    )
    config = CompareConfig(strata=PVC_STRATA, str_regions=str_regions)
    expected = _run_compare(case, config, 0)
    assert expected[0]
    assert "has no region type" in str(expected[1])
    _check_same_run(_run_compare(case, config, finder_threads), expected)


@pytest.mark.parametrize("finder_threads", _FINDER_THREADS)
@pytest.mark.parametrize(
    "config",
    [
        CompareConfig(ploidy_a=3),
        CompareConfig(ploidy_b=0),
        CompareConfig(phasing_cap=0),
        CompareConfig(span_limit=0),
        CompareConfig(chain_distance=-1),
    ],
    ids=["ploidy_a 3", "ploidy_b 0", "phasing_cap 0", "span_limit 0", "chain_distance -1"],
)
def test_setting_out_of_range_raises_the_same_value_error(
    tmp_path: Path, config: CompareConfig, finder_threads: int
) -> None:
    case = write_compare_case(tmp_path)
    expected = _run_compare(case, config, 0)
    assert isinstance(expected[1], ValueError)
    _check_same_run(_run_compare(case, config, finder_threads), expected)


# The case of the bound and the stop: chr1 of _BIG_CONTIG_BASES random bases, and file A with
# _BIG_NUM_RECORDS SNPs on it, one at every other base, which with a chain distance of 0 do
# not chain, so that the first row comes after a few records.
_BIG_NUM_RECORDS = 50_000
_BIG_CONTIG_BASES = 2 * _BIG_NUM_RECORDS + 10
_BIG_CONFIG = CompareConfig(chain_distance=0)


def _write_big_case(directory: Path) -> CompareCase:
    """The four-contig case with chr1 grown to _BIG_CONTIG_BASES random bases, file A with
    _BIG_NUM_RECORDS SNPs on chr1, and file B with one SNP at its first base."""
    bases = "".join(random.Random(20261020).choices("ACGT", k=_BIG_CONTIG_BASES))
    first_snp = f"chr1 1 {bases[0]} {_make_alt(bases[0])} 50 0/1:40 0/0:40"
    case = write_compare_case(
        directory, contigs={**CONTIGS, "chr1": bases}, lines_a=[], lines_b=[first_snp]
    )
    header = case.a.read_text(encoding="utf-8")
    vcf_lines = [
        f"chr1\t{pos}\t.\t{bases[pos - 1]}\t{_make_alt(bases[pos - 1])}\t50\tPASS\t.\tGT:GQ"
        "\t0/1:40\t0/0:40\n"
        for pos in range(1, _BIG_CONTIG_BASES - 10, 2)
    ]
    assert len(vcf_lines) == _BIG_NUM_RECORDS
    case.a.write_text(header + "".join(vcf_lines), encoding="utf-8")
    return case


def _make_alt(base: str) -> str:
    return "C" if base == "A" else "A"


def _make_counted_parse_vcf(
    counts: dict[Path, list[int]], release: threading.Event | None = None
) -> Callable[..., Iterator[VcfRecord]]:
    """A stand-in for parse_vcf that counts in counts[path][0] the records taken from it, and
    that, given release, waits for it before its first record; each path's counter is changed
    only by the reader of that path."""

    def parse_vcf_counted(path: Path, *args: object, **kwargs: object) -> Iterator[VcfRecord]:
        counter = counts[path]
        if release is not None:
            release.wait(_TIMEOUT)
        with closing(parse_vcf(path, *args, **kwargs)) as vcf_records:  # ty: ignore[invalid-argument-type]  # parse_vcf's own arguments, passed on
            for vcf_record in vcf_records:
                counter[0] += 1
                yield vcf_record

    return parse_vcf_counted


def _make_held_segment_contig(
    held_contigs: frozenset[str] | None, release: threading.Event | None, calls: list[str]
) -> Callable[[IndexedFasta, str], Iterator[GenomeSegment]]:
    """A stand-in for segment_contig that notes each contig in calls and, for the contigs of
    held_contigs, or for every contig when it is None, waits for release before it scans."""

    def segment_contig_held(reference: IndexedFasta, contig: str) -> Iterator[GenomeSegment]:
        calls.append(contig)
        if release is not None and (held_contigs is None or contig in held_contigs):
            release.wait(_TIMEOUT)
        return segment_contig(reference, contig)

    return segment_contig_held


def _start_compare(case: CompareCase, *, finder_threads: int) -> Generator[OutcomeRow]:
    """compare()'s stream over the case with _BIG_CONFIG, not started."""
    return compare(case.a, case.b, case.reference, _BIG_CONFIG, finder_threads=finder_threads)


def test_reader_held_by_the_finder_takes_one_buffer_and_two_batches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = _write_big_case(tmp_path)
    counts = {case.a: [0], case.b: [0]}
    release = threading.Event()
    calls: list[str] = []
    monkeypatch.setattr(compare_module, "parse_vcf", _make_counted_parse_vcf(counts))
    monkeypatch.setattr(
        compare_module, "segment_contig", _make_held_segment_contig(None, release, calls)
    )
    threads_before = set(threading.enumerate())
    rows = _start_compare(case, finder_threads=2)
    try:
        first_row = _start_call(lambda: next(rows))
        try:
            num_taken = _wait_until_steady(
                lambda: counts[case.a][0], 10_000, "the reader's stop on its full buffer"
            )
        finally:
            release.set()
        assert _get_result(first_row).contig == "chr1"
    finally:
        _call_with_timeout(rows.close)
    # The buffer full, one batch the loci side took before it waited on the finder for chr1,
    # and one the reader holds while it waits: 158 batches.
    assert num_taken == (_BUFFER_BATCHES + 2) * _BATCH_RECORDS == 10_112
    # Closed, the stream has waited for its readers.
    assert _get_new_threads(threads_before, "pvc_eval reader") == []


def test_close_after_first_row_ends_readers_and_starts_no_scan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # File A is indexed, so chr2 and chr3, listed by its index and with no record, are
    # planned; with one finder thread, the scan of chr2 is held, or not yet started, at the
    # close, and chr3 waits in line.
    case = _index_case(_write_big_case(tmp_path), contigs_a=["chr1", "chr2", "chr3"])
    release = threading.Event()
    calls: list[str] = []
    held = frozenset({"chr2", "chr3"})
    monkeypatch.setattr(
        compare_module, "segment_contig", _make_held_segment_contig(held, release, calls)
    )
    threads_before = set(threading.enumerate())
    rows = _start_compare(case, finder_threads=1)
    try:
        assert _call_with_timeout(lambda: next(rows)).contig == "chr1"
        _call_with_timeout(rows.close)
        num_calls = len(calls)
        # Closed, the stream has waited for its readers, before it closed the reference.
        assert _get_new_threads(threads_before, "pvc_eval reader") == []
    finally:
        release.set()
    assert _join_threads(threads_before, "pvc_eval", 5) == []
    assert num_calls >= 1
    assert len(calls) == num_calls


def test_planned_contigs_are_scanned_before_a_reader_reaches_them(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = write_compare_case(tmp_path)
    expected = _run_compare(case, CompareConfig(), 0)
    indexed_case = _index_case(case, contigs_a=["chr1", "chr3"])
    counts = {indexed_case.a: [0], indexed_case.b: [0]}
    readers_release = threading.Event()
    calls: list[str] = []
    monkeypatch.setattr(
        compare_module, "parse_vcf", _make_counted_parse_vcf(counts, readers_release)
    )
    monkeypatch.setattr(
        compare_module, "segment_contig", _make_held_segment_contig(frozenset(), None, calls)
    )
    run = _start_call(lambda: _collect_run(indexed_case, CompareConfig(), 2))
    try:
        _wait_until(lambda: len(calls) == 2, "the scans of the planned contigs")
        assert sorted(calls) == ["chr1", "chr3"]
        assert (counts[indexed_case.a][0], counts[indexed_case.b][0]) == (0, 0)
    finally:
        readers_release.set()
    _check_same_run(_get_result(run), expected)


def test_one_finder_thread_scans_at_most_two_contigs_ahead(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = write_compare_case(tmp_path)
    expected_rows, _ = _run_compare(case, CompareConfig(), 0)
    indexed_case = _index_case(case, contigs_a=["chr1", "chr2", "chr3", "chr4"])
    calls: list[str] = []
    monkeypatch.setattr(
        compare_module, "segment_contig", _make_held_segment_contig(frozenset(), None, calls)
    )
    threads_before = set(threading.enumerate())
    rows = compare(indexed_case.a, indexed_case.b, indexed_case.reference, finder_threads=1)
    try:
        # chr1:50's row comes once the loci side has taken chr1's STR regions, and before it
        # needs chr2's: chr2 and chr3 are ahead, and chr4 waits for one of them to be taken.
        first_row = _call_with_timeout(lambda: next(rows))
        _wait_until_steady(lambda: len(calls), 3, "the scans ahead")
        assert calls == ["chr1", "chr2", "chr3"]
        # The finder's thread waits for a contig until the stop, so it is alive here; the
        # readers of this short case may have ended.
        finder_threads = _get_new_threads(threads_before, "pvc_eval finder")
        assert [(thread.name, thread.daemon) for thread in finder_threads] == [
            ("pvc_eval finder 1", True)
        ]
        rest = _call_with_timeout(lambda: list(rows))
    finally:
        _call_with_timeout(rows.close)
    assert [first_row, *rest] == expected_rows


def test_contigs_the_loci_side_passes_leave_the_line_and_stop_counting_as_ahead(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Five contigs, all planned; records on chr1 and chr4 alone. With one finder thread, chr1
    # and chr2 are scanned first; once chr1 is taken, chr3 is started, unless the loci side
    # has passed it by then; chr2 and chr3 have no record and are passed. Once chr4 is taken,
    # nothing is ahead, so chr5 starts.
    contigs = {**CONTIGS, "chr5": CONTIGS["chr4"]}
    base_150 = CONTIGS["chr4"][149]
    lines_a = [
        "chr1 50 T C 50 0/1:40 1/1:35",
        "chr4 51 C A 30 0/1:30 0/0:30",
        f"chr4 150 {base_150} {_make_alt(base_150)} 30 0/1:30 0/0:30",
    ]
    case = write_compare_case(tmp_path, contigs=contigs, lines_a=lines_a, lines_b=[])
    expected_rows, _ = _run_compare(case, CompareConfig(), 0)
    assert [row.contig for row in expected_rows] == ["chr1", "chr4", "chr4"]
    indexed_case = _index_case(case, contigs_a=list(contigs))
    calls: list[str] = []
    monkeypatch.setattr(
        compare_module, "segment_contig", _make_held_segment_contig(frozenset(), None, calls)
    )
    rows = compare(indexed_case.a, indexed_case.b, indexed_case.reference, finder_threads=1)
    try:
        # The row of chr4:51 comes once chr4's STR regions are taken.
        taken_rows = [_call_with_timeout(lambda: next(rows)) for _ in range(2)]
        _wait_until_steady(lambda: len(calls), 4, "the scans ahead")
        assert calls in (
            ["chr1", "chr2", "chr3", "chr4", "chr5"],
            ["chr1", "chr2", "chr4", "chr5"],
        )
        taken_rows += _call_with_timeout(lambda: list(rows))
    finally:
        _call_with_timeout(rows.close)
    assert taken_rows == expected_rows


def test_file_of_50000_records_read_to_its_last_row_gives_the_rows_without_threads(
    tmp_path: Path,
) -> None:
    # The reader's buffer fills and empties some five times, so a reader or a loci side that
    # missed a wake would hang.
    case = _write_big_case(tmp_path)
    expected = _run_compare(case, _BIG_CONFIG, 0)
    assert len(expected[0]) > 10_000
    _check_same_run(_run_compare(case, _BIG_CONFIG, 1), expected)
