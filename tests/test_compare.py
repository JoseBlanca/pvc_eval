"""compare() with no thread, finder_threads 0, on the four-contig case of
tests/compare_case.py: its rows are those of the nine steps of doc/specs/compare.md, "The
stages, joined", run one after another; it scans the contigs with a record and, with
--regions, only those it names; its errors come after the same rows; and its INFO lines.
tests/test_compare_threads.py runs the same cases with threads."""

import logging
import re
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path

import pytest
from tests.compare_case import (
    CONTIGS,
    LINES_A,
    LINES_B,
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
                partly_outside = regions is not None and not regions.covers(
                    locus.contig, locus.start, locus.end
                )
                rows.extend(
                    make_outcome_rows(
                        reconciled,
                        sample_pairs,
                        kind=locus_kind(locus),
                        stratum=get_stratum(locus),
                        partly_outside=partly_outside,
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
        rows.extend(compare(case.a, case.b, case.reference, config, finder_threads=0))
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
    if option.startswith("--regions "):
        bed_line = option.removeprefix("--regions ")
        return CompareConfig(regions=write_bed(directory / "regions.bed", [bed_line]))
    if option == "--strata pvc":
        return CompareConfig(strata=PVC_STRATA)
    if option == "--strata BED":
        strata_lines = ["chr1 0 120 low", "chr3 0 300 all", "chr3 100 130 repeat"]
        return CompareConfig(strata=write_bed(directory / "strata.bed", strata_lines))
    if option == "--ploidy-a 1":
        return CompareConfig(ploidy_a=1)
    if option == "--ploidy-b 1":
        return CompareConfig(ploidy_b=1)
    if option == "--chain-distance 200 --phasing-cap 1":
        return CompareConfig(chain_distance=200, phasing_cap=1)
    if option == "--span-limit 1":
        return CompareConfig(span_limit=1)
    if option == "--sample-map":
        sample_map = directory / "pairs.tsv"
        sample_map.write_text("HG003\tHG004\n", encoding="utf-8")
        return CompareConfig(sample_map=sample_map)
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
    ("option", "num_rows"),
    [
        ("no option", 7),
        ("--regions chr3 0 300", 2),
        # Within chr3, whose loci are 101-120 and 250 (1-based): none, the first, the second.
        ("--regions chr3 0 50", 0),
        ("--regions chr3 100 101", 1),
        ("--regions chr3 200 300", 1),
        # A region that ends where the locus 101-120 starts shares no position with it.
        ("--regions chr3 90 100", 0),
        ("--strata pvc", 7),
        ("--strata BED", 7),
        ("--str-regions", 7),
        # File A's genotypes are diploid, so each sample of A is invalid; then of B.
        ("--ploidy-a 1", 7),
        ("--ploidy-b 1", 7),
        # chr1's three records chain into one locus, whose two unphased heterozygous
        # records give HG002 of file A two possible genotypes, over the cap; chr2 and chr3
        # chain too.
        ("--chain-distance 200 --phasing-cap 1", 3),
        ("--span-limit 1", 7),
        ("--sample-map", 7),
    ],
)
def test_rows_are_those_of_the_nine_steps(tmp_path: Path, option: str, num_rows: int) -> None:
    # chr1:115 AC>A deletes a C of CCC at 116-118, so its span reaches two bases past its
    # REF, which --span-limit 1 cuts.
    lines_a = [*LINES_A[:2], "chr1 115 AC A 20 0/1:30 0/0:30", *LINES_A[2:]]
    case = write_compare_case(tmp_path, lines_a=lines_a)
    config = _make_config(option, tmp_path)
    run = _run_compare(case, config)
    expected = _run_nine_steps(case, config, _CONTIGS_WITH_RECORDS)
    assert len(expected[0]) == num_rows
    assert expected[1] is None
    _check_same_run(run, expected)


def test_record_on_last_contig_gives_rows_of_the_nine_steps(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # build_loci asks for a segment past chr4:51, the last contig of the .fai.
    lines_a = [*LINES_A, "chr4 51 C A 30 0/1:30 0/0:30"]
    case = write_compare_case(tmp_path, lines_a=lines_a)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    run = _run_compare(case, CompareConfig())
    expected = _run_nine_steps(case, CompareConfig(), [*_CONTIGS_WITH_RECORDS, "chr4"])
    assert [(row.contig, row.start, row.end) for row in expected[0]][-1] == ("chr4", 50, 51)
    _check_same_run(run, expected)
    assert _get_scan_messages(caplog)[-1] == "chr4: finding STR regions in 200 bases."


def test_locus_widened_by_str_region_gives_its_row(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path)
    rows = [
        row
        for row in compare(case.a, case.b, case.reference, finder_threads=0)
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


@pytest.mark.parametrize(
    ("bed_lines", "expected"),
    [
        (["chr3 100 120"], Outcome.SAME_GENOTYPE),
        # Two lines that lie next to each other cover the locus as one region.
        (["chr3 100 110", "chr3 110 120"], Outcome.SAME_GENOTYPE),
        (["chr3 100 101"], Outcome.PARTLY_OUTSIDE),
        (["chr3 119 300"], Outcome.PARTLY_OUTSIDE),
        (["chr3 100 110", "chr3 111 120"], Outcome.PARTLY_OUTSIDE),
        # Only the locus's last base, or only its first, lies outside.
        (["chr3 100 119"], Outcome.PARTLY_OUTSIDE),
        (["chr3 101 120"], Outcome.PARTLY_OUTSIDE),
    ],
    ids=[
        "inside",
        "inside_two_lines",
        "past_the_end",
        "past_the_start",
        "across_a_gap",
        "last_base_outside",
        "first_base_outside",
    ],
)
def test_locus_reaching_past_the_regions_is_partly_outside_them(
    tmp_path: Path, bed_lines: list[str], expected: Outcome
) -> None:
    # The locus at chr3:101-120 (1-based), 0-based (100, 120), same genotype without
    # --regions, as test_locus_widened_by_str_region_gives_its_row shows.
    case = write_compare_case(tmp_path)
    config = CompareConfig(regions=write_bed(tmp_path / "regions.bed", bed_lines))
    rows = [
        (row.contig, row.start, row.end, row.outcome, row.genotype_a, row.genotype_b)
        for row in compare(case.a, case.b, case.reference, config, finder_threads=0)
        if row.contig == "chr3" and row.start < 200
    ]
    assert rows == [("chr3", 100, 120, expected, (0, 1), (0, 1))]


def test_contigs_with_records_are_scanned_and_contig_4_is_not(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = write_compare_case(tmp_path)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    list(compare(case.a, case.b, case.reference, finder_threads=0))
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
    config = _make_config("--regions chr3 0 300", tmp_path)
    list(compare(case.a, case.b, case.reference, config, finder_threads=0))
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


@pytest.mark.parametrize(
    ("bad_contig", "lines_a", "lines_b", "loci_before"),
    [
        # chr2 is taken at the start, as file B's first record; it is scanned at that
        # record, after the loci of chr1:50 and chr1:100.
        ("chr2", LINES_A, LINES_B, [("chr1", 49, 50), ("chr1", 99, 116)]),
        # chr1:250 reads the GENERIC segment at chr2:1, and the locus of chr1:50 is given;
        # chr3 is scanned at chr3:105, before the locus of chr1:250 is given. A chain of
        # segment_contig over chr1 and chr3 scans chr3 at chr1:250, and gives no row.
        (
            "chr3",
            ["chr1 50 T C 50 0/1:40 1/1:35", "chr1 250 T C 30 0/1:30 0/0:30", LINES_A[3]],
            [],
            [("chr1", 49, 50)],
        ),
        # chr1:205, in (AGT)6, reads past chr1's last STR region while chr2:125 of file B
        # is taken, so chr2 is scanned there, before the locus of chr1:50 is given. A chain
        # of segment_contig over chr1 and chr2 scans chr2 at chr2:125, after that locus.
        (
            "chr2",
            ["chr1 50 T C 50 0/1:40 1/1:35", "chr1 205 G A 30 0/1:30 0/0:30"],
            ["chr2 125 A G 40 0/1:33 1/1:30"],
            [],
        ),
    ],
    ids=[
        "chr2, default records",
        "chr3, after a record past chr1's STR regions",
        "chr2 taken early",
    ],
)
def test_fasta_byte_not_a_letter_stops_the_finder_after_its_rows(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    bad_contig: str,
    lines_a: Sequence[str],
    lines_b: Sequence[str],
    loci_before: list[tuple[str, int, int]],
) -> None:
    # The contig grows to 140,000 bases with a 0 at 139,001 (1-based), more than
    # BLOCK_BASES, 65,536, from its records, so that parse_vcf's window does not read it.
    bases = CONTIGS[bad_contig] + "N" * 139_700
    bases = bases[:139_000] + "0" + bases[139_001:]
    contigs = {**CONTIGS, bad_contig: bases}
    case = write_compare_case(tmp_path, contigs=contigs, lines_a=lines_a, lines_b=lines_b)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    rows, error = _run_compare(case, CompareConfig())
    assert [(row.contig, row.start, row.end) for row in rows] == loci_before
    assert str(error) == (
        f"{case.reference}: {bad_contig}:139001 is not a letter where its .fai, line"
        f" {bad_contig[-1]}, puts a base. Make the .fai again with `samtools faidx"
        f" {case.reference}`, check that the FASTA holds only letters between its header"
        " lines, and run again."
    )
    assert _get_scan_messages(caplog)[-1] == f"{bad_contig}: finding STR regions in 140,000 bases."


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
        (CompareConfig(ploidy_a=True), 2, "ploidy_a"),
        (CompareConfig(), True, "finder_threads"),
    ],
    ids=[
        "ploidy_a 3",
        "ploidy_b 0",
        "phasing_cap 0",
        "span_limit 0",
        "chain_distance -1",
        "finder_threads -1",
        "finder_threads 1.5",
        "ploidy_a True",
        "finder_threads True",
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
    ("config", "finder_threads"),
    [
        (CompareConfig(ploidy_a=1, ploidy_b=1), 2),
        (CompareConfig(phasing_cap=1), 2),
        (CompareConfig(chain_distance=0), 2),
        (CompareConfig(span_limit=None), 2),
        (CompareConfig(), 0),
    ],
    ids=["ploidy 1", "phasing_cap 1", "chain_distance 0", "span_limit None", "finder_threads 0"],
)
def test_setting_at_its_edge_gives_rows(
    tmp_path: Path, config: CompareConfig, finder_threads: int
) -> None:
    case = write_compare_case(tmp_path)
    rows = compare(case.a, case.b, case.reference, config, finder_threads=finder_threads)
    assert next(rows).contig == "chr1"


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
        # doc/specs/regions_read.md: a contig --regions does not name is not read, so it
        # has no locus and no line.
        ("--regions chr3 0 300", ["chr3: 2 loci compared, 0 skipped by --regions; "]),
        ("--regions chr3 200 300", ["chr3: 1 locus compared, 1 skipped by --regions; "]),
    ],
)
def test_contig_done_line_counts_its_loci(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, option: str, expected_lines: list[str]
) -> None:
    case = write_compare_case(tmp_path)
    caplog.set_level(logging.INFO, logger=_LOGGER)
    config = _make_config(option, tmp_path)
    list(compare(case.a, case.b, case.reference, config, finder_threads=0))
    done_lines = _get_messages(caplog, "compared")
    assert len(done_lines) == len(expected_lines)
    for done_line, expected_start in zip(done_lines, expected_lines, strict=True):
        assert re.fullmatch(re.escape(expected_start) + _TIME, done_line), done_line


def test_contig_done_line_for_one_locus_says_locus(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    case = write_compare_case(tmp_path, lines_a=LINES_A[:1], lines_b=())
    caplog.set_level(logging.INFO, logger=_LOGGER)
    list(compare(case.a, case.b, case.reference, finder_threads=0))
    done_lines = _get_messages(caplog, "compared")
    assert len(done_lines) == 1
    assert re.fullmatch(r"chr1: 1 locus compared; " + _TIME, done_lines[0])


# doc/specs/regions_read.md: chr1:50 of file A with a REF of G where the reference has T,
# and chr4:51 with G where it has C.
_LINE_REF_DIFFERS_ON_CHR1 = "chr1 50 G C 50 0/1:40 1/1:35"
_LINE_REF_DIFFERS_ON_CHR4 = "chr4 51 G A 30 0/1:30 0/0:30"
# The regions that name contig 3 alone: whole, and cut inside it before, at and after its
# loci, 101-120 and 250 (1-based).
_REGIONS_NAMING_CHR3 = ["chr3 0 300", "chr3 0 50", "chr3 100 101", "chr3 200 300"]


# File A of the four-contig case with chr1:50's REF wrong, and with chr4:51's added.
_LINES_A_REF_DIFFERS_ON_CHR1 = (_LINE_REF_DIFFERS_ON_CHR1, *LINES_A[1:])
_LINES_A_REF_DIFFERS_ON_CHR4 = (*LINES_A, _LINE_REF_DIFFERS_ON_CHR4)


@pytest.mark.parametrize("bed_line", _REGIONS_NAMING_CHR3)
def test_ref_that_differs_on_contig_1_does_not_stop_a_run_with_regions_naming_contig_3(
    tmp_path: Path, bed_line: str
) -> None:
    (tmp_path / "differs").mkdir()
    (tmp_path / "agrees").mkdir()
    case = write_compare_case(tmp_path / "differs", lines_a=_LINES_A_REF_DIFFERS_ON_CHR1)
    agreeing_case = write_compare_case(tmp_path / "agrees")
    config = _make_config(f"--regions {bed_line}", tmp_path)
    expected = _run_nine_steps(agreeing_case, config, _CONTIGS_WITH_RECORDS)
    assert expected[1] is None
    _check_same_run(_run_compare(case, config), expected)


def test_ref_that_differs_on_contig_1_stops_a_run_without_regions(tmp_path: Path) -> None:
    case = write_compare_case(tmp_path, lines_a=_LINES_A_REF_DIFFERS_ON_CHR1)
    rows, error = _run_compare(case, CompareConfig())
    assert rows == []
    assert str(error) == (
        f"{case.a}, line 3: the REF of chr1:50 is 'G', and the reference, read in upper case,"
        " has 'T' there. It is the first record read with no symbolic allele, so the VCF"
        " looks made on another reference: compare against the reference it was called on,"
        " and run again."
    )


def test_bed_line_past_the_end_of_contig_4_reads_no_record_of_it(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # chr4 has 200 bases, so the BED line at 300-400 names no contig, and chr4:51's REF,
    # which differs from the reference, is not read.
    (tmp_path / "differs").mkdir()
    (tmp_path / "agrees").mkdir()
    case = write_compare_case(tmp_path / "differs", lines_a=_LINES_A_REF_DIFFERS_ON_CHR4)
    agreeing_case = write_compare_case(tmp_path / "agrees")
    regions = write_bed(tmp_path / "regions.bed", ["chr3 0 300", "chr4 300 400"])
    config = CompareConfig(regions=regions)
    expected = _run_nine_steps(agreeing_case, config, _CONTIGS_WITH_RECORDS)
    assert len(expected[0]) == 2
    caplog.set_level(logging.INFO, logger=_LOGGER)
    _check_same_run(_run_compare(case, config), expected)
    assert _get_scan_messages(caplog) == ["chr3: finding STR regions in 300 bases."]


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
