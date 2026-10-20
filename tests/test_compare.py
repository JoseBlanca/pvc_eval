"""compare() with no thread, on the four-contig case of tests/compare_case.py: its rows are
those of the nine steps of doc/specs/compare.md, "The stages, joined", run one after
another; it scans the contigs with a record and, with --regions, only those it names; its
errors come after the same rows; and its INFO lines."""

import logging
import re
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path

import pytest
from tests.compare_case import (
    CONTIGS,
    LINES_A,
    STR_REGIONS,
    CompareCase,
    write_bed,
    write_compare_case,
)

from pvc_eval.compare import PVC_STRATA, CompareConfig, _format_duration, compare
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.evaluate import make_outcome_rows
from pvc_eval.fasta import IndexedFasta
from pvc_eval.kind import Kind, locus_kind
from pvc_eval.loci import build_loci
from pvc_eval.locus import Locus
from pvc_eval.outcome import Outcome, OutcomeRow
from pvc_eval.reconcile import reconcile_loci
from pvc_eval.regions import read_regions, select_loci
from pvc_eval.samples import make_sample_pairs
from pvc_eval.segments import GenomeSegment, RegionType, read_str_regions, segment_contig
from pvc_eval.strata import pvc_stratum, read_strata
from pvc_eval.vcf import parse_vcf

_LOGGER = "pvc_eval.compare"
_CONTIGS_WITH_RECORDS = ("chr1", "chr2", "chr3")
_TIME = r"(\d+ h )?(\d+ min )?\d+ s since the start\."

type _Run = tuple[list[OutcomeRow], PvcEvalInputError | None]


def _run_nine_steps(
    case: CompareCase, config: CompareConfig, scanned_contigs: Sequence[str]
) -> _Run:
    """The rows of the nine steps of "The stages, joined", one after another, and the error
    that stopped them. The segments are segment_contig's over scanned_contigs, GENERIC
    segments included, each contig scanned when build_loci first reads its segments."""
    rows: list[OutcomeRow] = []
    try:
        with IndexedFasta(case.reference) as reference:  # step 1
            sample_pairs = make_sample_pairs(case.a, case.b, sample_map=config.sample_map)  # 2
            regions = (  # step 3
                None if config.regions is None else read_regions(config.regions, reference.contigs)
            )
            if config.strata is PVC_STRATA:
                get_stratum = pvc_stratum
            elif config.strata is None:
                get_stratum = _get_no_stratum
            else:
                get_stratum = read_strata(config.strata, reference.contigs).stratum
            span_limit = config.span_limit
            # step 4
            records_a = parse_vcf(case.a, reference, span_limit=span_limit, ploidy=config.ploidy_a)
            records_b = parse_vcf(case.b, reference, span_limit=span_limit, ploidy=config.ploidy_b)
            if config.str_regions is None:  # step 5
                segments = _segment_contigs(reference, scanned_contigs)
            else:
                segments = read_str_regions(config.str_regions, reference)
            loci = build_loci(  # step 6
                records_a,
                records_b,
                segments,
                reference.contigs,
                chain_distance=config.chain_distance,
            )
            reconciled_loci = reconcile_loci(  # steps 7 and 8
                select_loci(loci, regions),
                reference,
                sample_pairs,
                ploidy_a=config.ploidy_a,
                ploidy_b=config.ploidy_b,
                phasing_cap=config.phasing_cap,
            )
            for reconciled in reconciled_loci:  # step 9
                locus = reconciled.locus
                rows.extend(
                    make_outcome_rows(
                        reconciled, sample_pairs, kind=locus_kind(locus), stratum=get_stratum(locus)
                    )
                )
    except PvcEvalInputError as error:
        return rows, error
    return rows, None


def _get_no_stratum(locus: Locus) -> None:  # noqa: ARG001 -- the signature of a stratum
    return None


def _segment_contigs(reference: IndexedFasta, contigs: Iterable[str]) -> Iterator[GenomeSegment]:
    for contig in contigs:
        yield from segment_contig(reference, contig)


def _run_compare(case: CompareCase, config: CompareConfig) -> _Run:
    rows: list[OutcomeRow] = []
    try:
        rows.extend(compare(case.a, case.b, case.reference, config))
    except PvcEvalInputError as error:
        return rows, error
    return rows, None


def _check_same_run(run: _Run, expected: _Run) -> None:
    rows, error = run
    expected_rows, expected_error = expected
    assert rows == expected_rows
    assert type(error) is type(expected_error)
    assert str(error) == str(expected_error)


def _get_messages(caplog: pytest.LogCaptureFixture, text: str) -> list[str]:
    """The messages compare() logged that hold text."""
    return [
        record.getMessage()
        for record in caplog.records
        if record.name == _LOGGER and text in record.getMessage()
    ]


def _get_scan_messages(caplog: pytest.LogCaptureFixture) -> list[str]:
    return _get_messages(caplog, "finding STR regions")


def _make_config(option: str, directory: Path) -> CompareConfig:
    if option == "no option":
        return CompareConfig()
    if option == "--regions chr3":
        return CompareConfig(regions=write_bed(directory / "regions.bed", ["chr3 0 300"]))
    if option == "--strata pvc":
        return CompareConfig(strata=PVC_STRATA)
    if option == "--strata BED":
        strata_lines = ["chr1 0 120 low", "chr3 0 300 all", "chr3 100 130 repeat"]
        return CompareConfig(strata=write_bed(directory / "strata.bed", strata_lines))
    assert option == "--str-regions"
    # Not the finder's: chr2's region left out, and one at chr3 241-260 the finder does
    # not find, which widens the locus of chr3:250.
    str_region_lines = [
        "chr1 100 116 str AC",
        "chr1 200 218 str AGT",
        "chr3 100 120 str CA",
        "chr3 240 260 str-bundle .",
    ]
    return CompareConfig(str_regions=write_bed(directory / "str.bed", str_region_lines))


def test_four_contig_case_str_regions_are_segment_contigs(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path)
    with IndexedFasta(case.reference) as reference:
        str_regions = tuple(
            segment
            for contig in reference.contigs
            for segment in segment_contig(reference, contig)
            if segment.region_type is not RegionType.GENERIC
        )
    assert str_regions == STR_REGIONS


@pytest.mark.parametrize(
    "option", ["no option", "--regions chr3", "--strata pvc", "--strata BED", "--str-regions"]
)
def test_rows_are_those_of_the_nine_steps(tmp_path: Path, option: str) -> None:
    case = write_compare_case(tmp_path)
    config = _make_config(option, tmp_path)
    run = _run_compare(case, config)
    expected = _run_nine_steps(case, config, _CONTIGS_WITH_RECORDS)
    assert expected[0]
    assert expected[1] is None
    _check_same_run(run, expected)


def test_locus_widened_by_str_region_gives_its_row(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path)
    rows = [
        row
        for row in compare(case.a, case.b, case.reference)
        if row.contig == "chr3" and row.start < 200
    ]
    # chr3:105 C>G in both files, 0/1 for HG002 in both, in (CA)10 at 101-120 (1-based).
    assert rows == [
        OutcomeRow(
            contig="chr3",
            start=100,
            end=120,
            kind=Kind.STR,
            stratum=None,
            str_regions=(GenomeSegment("chr3", 100, 120, RegionType.STR, "CA"),),
            sample_a="HG002",
            sample_b="HG002",
            outcome=Outcome.SAME_GENOTYPE,
            alleles=("CACACACACACACACACACA", "CACAGACACACACACACACA"),
            genotype_a=(0, 1),
            genotype_b=(0, 1),
            possible_a=1,
            possible_b=1,
            quality_a=40.0,
            quality_b=41.0,
            gq_a=33,
            gq_b=36,
            records_a=1,
            records_b=1,
            invalid_reason_a=None,
            invalid_reason_b=None,
        )
    ]


def test_contigs_with_records_are_scanned_and_contig_4_is_not(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = write_compare_case(tmp_path)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    list(compare(case.a, case.b, case.reference))
    assert _get_scan_messages(caplog) == [
        "chr1: finding STR regions in 300 bases.",
        "chr2: finding STR regions in 300 bases.",
        "chr3: finding STR regions in 300 bases.",
    ]


def test_regions_naming_contig_3_scan_contig_3_alone(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = write_compare_case(tmp_path)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    list(compare(case.a, case.b, case.reference, _make_config("--regions chr3", tmp_path)))
    assert _get_scan_messages(caplog) == ["chr3: finding STR regions in 300 bases."]


@pytest.mark.parametrize(
    "last_line_on_chr1",
    ["chr1 250 T C 30 0/1:30 0/0:30", "chr1 300 C T 30 0/1:30 0/0:30"],
    ids=["past the last STR region", "ending on the last base"],
)
@pytest.mark.parametrize(
    ("next_line", "widened_locus"),
    [
        ("chr3 105 C G 40 0/1:33 1/1:30", ("chr3", 100, 120)),  # in (CA)10
        ("chr2 125 A G 40 0/1:33 1/1:30", ("chr2", 120, 140)),  # in (ATTC)5
    ],
    ids=["then contig 3", "then contig 2"],
)
def test_record_past_last_str_region_then_next_contig_gives_rows_of_the_nine_steps(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    last_line_on_chr1: str,
    next_line: str,
    widened_locus: tuple[str, int, int],
) -> None:
    # File B has nothing past chr1:50, so when the record at chr1:250 or chr1:300 reads the
    # GENERIC segment at chr2:1, the next record of file A has not been taken. On chr3,
    # chr2 has no record; on chr2, chr2 is decided again at that record.
    lines_a = ["chr1 50 T C 50 0/1:40 1/1:35", last_line_on_chr1, next_line]
    lines_b = ["chr1 50 T C 50 0/1:40 0/0:35"]
    case = write_compare_case(tmp_path, lines_a=lines_a, lines_b=lines_b)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    next_contig = widened_locus[0]
    run = _run_compare(case, CompareConfig())
    expected = _run_nine_steps(case, CompareConfig(), ["chr1", next_contig])
    _check_same_run(run, expected)
    assert [(row.contig, row.start, row.end) for row in run[0]][-1] == widened_locus
    assert _get_scan_messages(caplog) == [
        "chr1: finding STR regions in 300 bases.",
        f"{next_contig}: finding STR regions in 300 bases.",
    ]


def test_vcf_out_of_order_on_contig_3_stops_after_the_same_rows(tmp_path: Path) -> None:
    lines_a = [*LINES_A[:3], LINES_A[4], LINES_A[3]]
    case = write_compare_case(tmp_path, lines_a=lines_a)
    run = _run_compare(case, CompareConfig())
    expected = _run_nine_steps(case, CompareConfig(), _CONTIGS_WITH_RECORDS)
    assert expected[0]
    assert isinstance(expected[1], PvcEvalInputError)
    _check_same_run(run, expected)


def test_fasta_byte_not_a_letter_stops_the_finder_after_the_same_rows(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # chr2 grows to 140,000 bases with a 0 at 139,001 (1-based), more than BLOCK_BASES,
    # 65,536, from its records at 120 and 200, so that parse_vcf's window does not read it.
    chr2 = CONTIGS["chr2"] + "N" * 139_700
    chr2 = chr2[:139_000] + "0" + chr2[139_001:]
    case = write_compare_case(tmp_path, contigs={**CONTIGS, "chr2": chr2})
    caplog.set_level(logging.INFO, logger=_LOGGER)
    run = _run_compare(case, CompareConfig())
    expected = _run_nine_steps(case, CompareConfig(), _CONTIGS_WITH_RECORDS)
    assert expected[0]
    assert "chr2:139001 is not a letter" in str(expected[1])
    _check_same_run(run, expected)
    assert _get_scan_messages(caplog)[-1] == "chr2: finding STR regions in 140,000 bases."


def test_strata_pvc_with_untyped_str_region_stops_after_the_same_rows(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path)
    str_region_lines = [
        "chr1 100 116 str AC",
        "chr1 200 218 str AGT",
        "chr2 120 140 str ATTC",
        "chr3 100 120",
    ]
    str_regions = write_bed(tmp_path / "str.bed", str_region_lines)
    config = CompareConfig(strata=PVC_STRATA, str_regions=str_regions)
    run = _run_compare(case, config)
    expected = _run_nine_steps(case, config, _CONTIGS_WITH_RECORDS)
    assert expected[0]
    assert "has no region type" in str(expected[1])
    _check_same_run(run, expected)


@pytest.mark.parametrize(
    ("config", "finder_threads", "field"),
    [
        (CompareConfig(ploidy_a=3), 2, "ploidy_a"),
        (CompareConfig(ploidy_b=0), 2, "ploidy_b"),
        (CompareConfig(phasing_cap=0), 2, "phasing_cap"),
        (CompareConfig(span_limit=0), 2, "span_limit"),
        (CompareConfig(chain_distance=-1), 2, "chain_distance"),
        (CompareConfig(), -1, "finder_threads"),
        (CompareConfig(), 1.5, "finder_threads"),
    ],
    ids=[
        "ploidy_a 3",
        "ploidy_b 0",
        "phasing_cap 0",
        "span_limit 0",
        "chain_distance -1",
        "finder_threads -1",
        "finder_threads 1.5",
    ],
)
def test_setting_out_of_range_raises_value_error_naming_it(
    tmp_path: Path, config: CompareConfig, finder_threads: int, field: str
) -> None:
    case = write_compare_case(tmp_path)
    rows = compare(case.a, case.b, case.reference, config, finder_threads=finder_threads)
    with pytest.raises(ValueError, match=rf"^{field} "):
        next(rows)


@pytest.mark.parametrize(
    ("option", "expected_lines"),
    [
        (
            "no option",
            [
                "chr1: 3 loci compared; ",
                "chr2: 2 loci compared; ",
                "chr3: 2 loci compared; ",
            ],
        ),
        (
            "--regions chr3",
            [
                "chr1: 0 loci compared, 3 skipped by --regions; ",
                "chr2: 0 loci compared, 2 skipped by --regions; ",
                "chr3: 2 loci compared, 0 skipped by --regions; ",
            ],
        ),
    ],
)
def test_contig_done_line_counts_its_loci(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, option: str, expected_lines: list[str]
) -> None:
    case = write_compare_case(tmp_path)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    list(compare(case.a, case.b, case.reference, _make_config(option, tmp_path)))
    done_lines = _get_messages(caplog, "compared")
    assert len(done_lines) == len(expected_lines)
    for done_line, expected_start in zip(done_lines, expected_lines, strict=True):
        assert re.fullmatch(re.escape(expected_start) + _TIME, done_line), done_line


def test_contig_done_line_for_one_locus_says_locus(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = write_compare_case(tmp_path, lines_a=LINES_A[:1], lines_b=())
    caplog.set_level(logging.INFO, logger=_LOGGER)
    list(compare(case.a, case.b, case.reference))
    done_lines = _get_messages(caplog, "compared")
    assert len(done_lines) == 1
    assert re.fullmatch(r"chr1: 1 locus compared; " + _TIME, done_lines[0])


@pytest.mark.parametrize(
    ("seconds", "expected"),
    [
        (0.4, "0 s"),
        (12.9, "12 s"),
        (60.0, "1 min 0 s"),
        (252.7, "4 min 12 s"),
        (3725.0, "1 h 2 min 5 s"),
    ],
)
def test_duration_is_formatted_in_whole_seconds(seconds: float, expected: str) -> None:
    assert _format_duration(seconds) == expected
