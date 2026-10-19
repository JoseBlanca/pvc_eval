import gzip
import itertools
import logging
import random
from collections.abc import Callable, Iterator
from contextlib import ExitStack
from pathlib import Path

import pytest

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import FastaWindow, IndexedFasta
from pvc_eval.vcf import (
    HAPLOID_GT,
    MALFORMED_GT,
    MISSING,
    VcfGenotype,
    VcfRecord,
    _parse_genotype,
    parse_vcf,
    parse_vcf_samples,
)

# The records of these tests agree with this reference: chr1 has A at 0-based 3, 4 and 6,
# and CA at 9 and 10. Its .fai is what `samtools faidx ref.fa` (samtools 1.24) wrote.
REFERENCE_FA = ">chr1\nGGGAAGAGGCAGGGGG\n>chr2\nGGGGGGGG\n"
REFERENCE_FAI = "chr1\t16\t6\t16\t17\nchr2\t8\t29\t8\t9\n"
HEADER = (
    "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG003\n"
)

type MakeWindow = Callable[..., FastaWindow]


@pytest.fixture
def make_window(tmp_path: Path) -> Iterator[MakeWindow]:
    """Makes a FastaWindow over a FASTA and its .fai, written in tmp_path, and closes the
    FASTA at the end of the test."""
    numbers = itertools.count()
    with ExitStack() as stack:

        def make(fasta: str, fai: str, *, block_bases: int = 65_536) -> FastaWindow:
            path = tmp_path / f"ref{next(numbers)}.fa"
            path.write_text(fasta, encoding="utf-8")
            Path(f"{path}.fai").write_text(fai, encoding="utf-8")
            indexed_fasta = stack.enter_context(IndexedFasta(path))
            return FastaWindow(indexed_fasta, block_bases=block_bases)

        yield make


@pytest.fixture
def window(make_window: MakeWindow) -> FastaWindow:
    return make_window(REFERENCE_FA, REFERENCE_FAI)


def make_fasta(contigs: dict[str, str]) -> tuple[str, str]:
    """A FASTA with each contig on one line, and the .fai `samtools faidx` writes for it."""
    fasta = ""
    fai = ""
    for name, sequence in contigs.items():
        fasta += f">{name}\n"
        fai += f"{name}\t{len(sequence)}\t{len(fasta)}\t{len(sequence)}\t{len(sequence) + 1}\n"
        fasta += sequence + "\n"
    return fasta, fai


def get_placement(record: VcfRecord) -> tuple[int, tuple[str, ...], tuple[int, int], bool]:
    return record.start, record.alleles, record.span, record.ref_mismatch


def write_vcf(tmp_path: Path, body: str, *, gzipped: bool = False) -> Path:
    text = HEADER + body
    if gzipped:
        path = tmp_path / "a.vcf.gz"
        with gzip.open(path, "wt", encoding="utf-8") as vcf:
            vcf.write(text)
    else:
        path = tmp_path / "a.vcf"
        path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("gt", "expected"),
    [
        ("0/1", VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=None)),
        ("1|0", VcfGenotype(alleles=(1, 0), phased=True, phase_set=None, quality=None)),
        ("./1", VcfGenotype(alleles=(MISSING, 1), phased=False, phase_set=None, quality=None)),
        (
            "./.",
            VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None),
        ),
        (".", VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None)),
        ("1", HAPLOID_GT),
        ("0", HAPLOID_GT),
        ("3", MALFORMED_GT),
        ("0/1/1", MALFORMED_GT),
        ("./././", MALFORMED_GT),
        ("0/3", MALFORMED_GT),
        ("A/1", MALFORMED_GT),
        ("A", MALFORMED_GT),
    ],
    ids=[
        "unphased",
        "phased",
        "one_missing",
        "both_missing",
        "single_dot_is_both_missing",
        "haploid_alt",
        "haploid_ref",
        "haploid_index_beyond_alts",
        "three_alleles",
        "four_missing",
        "index_beyond_alts",
        "not_a_number",
        "haploid_not_a_number",
    ],
)
def test_gt_gives_its_alleles_or_its_case(
    gt: str, expected: VcfGenotype | MALFORMED_GT | HAPLOID_GT
) -> None:
    assert _parse_genotype(gt, num_alleles=3, phase_set=None, quality=None) == expected


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "bgzipped"])
def test_records_are_read_with_0_based_starts(
    tmp_path: Path, window: FastaWindow, gzipped: bool
) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t10\t.\tCA\tC,CAA\t50.5\tPASS\t.\tGT:GQ:PS\t1|2:30:7\t0/0:.\n"
        "chr2\t3\t.\tG\tT\t.\tPASS\t.\tGT\t./1\t1/1\n",
        gzipped=gzipped,
    )

    first, second = parse_vcf(path, window)

    assert (first.contig, first.start, first.end, first.alleles) == (
        "chr1",
        9,
        11,
        ("CA", "C", "CAA"),
    )
    assert first.span == (9, 12)
    assert first.qual == 50.5
    assert first.genotypes == (
        VcfGenotype(alleles=(1, 2), phased=True, phase_set="7", quality=30),
        VcfGenotype(alleles=(0, 0), phased=False, phase_set=None, quality=None),
    )
    assert (second.contig, second.start, second.qual) == ("chr2", 2, None)
    assert second.genotypes[0] == VcfGenotype(
        alleles=(MISSING, 1), phased=False, phase_set=None, quality=None
    )


def test_samples_are_read_from_the_header(tmp_path: Path) -> None:
    assert parse_vcf_samples(write_vcf(tmp_path, "")) == ("HG002", "HG003")


def test_unsorted_file_stops_the_run_naming_both_lines(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t20\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t15\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(
        PvcEvalInputError, match=r"chr1:15 \(line 4\) comes after chr1:20 \(line 3\)"
    ):
        list(parse_vcf(path, window))


def test_contigs_out_of_the_reference_order_stop_the_run(
    tmp_path: Path, window: FastaWindow
) -> None:
    path = write_vcf(
        tmp_path,
        "chr2\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(PvcEvalInputError, match="not sorted"):
        list(parse_vcf(path, window))


def test_contig_missing_from_the_reference_stops_the_run(
    tmp_path: Path, window: FastaWindow
) -> None:
    path = write_vcf(tmp_path, "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError, match="contig chrUn is not in the reference"):
        list(parse_vcf(path, window))


@pytest.mark.parametrize(
    "alt",
    ["<DEL>", "A]chr2:20]", "[chr2:20[A", "A.", ".A"],
    ids=[
        "symbolic",
        "breakend_right",
        "breakend_left",
        "single_breakend_right",
        "single_breakend_left",
    ],
)
def test_symbolic_alleles_are_skipped_and_counted(
    tmp_path: Path, window: FastaWindow, caplog: pytest.LogCaptureFixture, alt: str
) -> None:
    path = write_vcf(
        tmp_path,
        f"chr1\t5\t.\tA\t{alt}\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, window))

    assert [(record.start, record.alleles) for record in records] == [(6, ("A", "T"))]
    assert "skipped 1 VCF record with symbolic alleles" in caplog.text


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("chr1\t5\t.\tA\tT\t.\t.\n", r"a.vcf, line 3: .*7 tab-separated columns"),
        ("chr1\n", r"a.vcf, line 3: the line has 1 tab-separated column,"),
        ("chr1\tfive\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n", r"a.vcf, line 3: POS is 'five'"),
        ("chr1\t0\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n", r"a.vcf, line 3: POS is '0'"),
        ("chr1\t5\t.\t\tT\t.\t.\t.\tGT\t0/1\t0/1\n", r"a.vcf, line 3: REF is empty"),
        ("chr1\t5\t.\tA\tT,\t.\t.\t.\tGT\t0/1\t0/1\n", r"a.vcf, line 3: ALT 'T,' has an empty"),
        ("chr1\t5\t.\tA\t\t.\t.\t.\tGT\t0/1\t0/1\n", r"a.vcf, line 3: ALT '' has an empty"),
    ],
    ids=[
        "seven_columns",
        "one_column",
        "pos_not_a_number",
        "pos_0",
        "ref_empty",
        "alt_empty",
        "alt_column_empty",
    ],
)
def test_line_that_is_not_a_vcf_record_stops_the_run(
    tmp_path: Path, window: FastaWindow, line: str, message: str
) -> None:
    path = write_vcf(tmp_path, line)

    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, window))


@pytest.mark.parametrize(
    ("line", "found"),
    [
        ("chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\n", "1 sample column,"),
        ("chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\t1/1\n", "3 sample columns,"),
        ("chr1\t5\t.\tA\tT\t.\t.\t.\n", "0 sample columns,"),
    ],
    ids=["one_sample_fewer", "one_sample_more", "no_format_column"],
)
def test_line_with_other_samples_than_the_header_stops_the_run(
    tmp_path: Path, window: FastaWindow, line: str, found: str
) -> None:
    path = write_vcf(tmp_path, line)

    with pytest.raises(
        PvcEvalInputError, match=rf"a.vcf, line 3: the line has {found} .* names 2 samples\."
    ):
        list(parse_vcf(path, window))


def test_sites_only_file_gives_records_without_genotypes(
    tmp_path: Path, window: FastaWindow
) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5\t.\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    (record,) = parse_vcf(path, window)

    assert (record.start, record.alleles, record.genotypes) == (4, ("A", "T"), ())


def test_blank_lines_are_skipped(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
        "\n"
        " \t \n"
        "chr1\t7\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    assert [record.start for record in parse_vcf(path, window)] == [
        4,
        6,
    ]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("##fileformat=VCFv4.3\nchr1\t5\t.\tA\tT\t.\t.\t.\n", r"a.vcf, line 2: .*#CHROM"),
        ("##fileformat=VCFv4.3\n", r"a.vcf has no #CHROM"),
    ],
    ids=["record_before_chrom_line", "only_meta_lines"],
)
def test_file_without_chrom_line_stops_the_run(
    tmp_path: Path, window: FastaWindow, text: str, message: str
) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(PvcEvalInputError, match=message):
        parse_vcf_samples(path)
    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, window))


def test_lower_case_alleles_are_read_in_upper_case_and_normalise_alike(
    tmp_path: Path, window: FastaWindow
) -> None:
    # VCF bases are case-insensitive; the reference is upper case.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\taa\ta\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    lower, upper = parse_vcf(path, window)

    assert get_placement(lower) == get_placement(upper) == (3, ("AA", "A"), (2, 5), False)


def test_file_that_is_not_utf8_text_stops_the_run(tmp_path: Path, window: FastaWindow) -> None:
    path = tmp_path / "a.vcf"
    # A Latin-1 é in a sample name.
    path.write_bytes(HEADER.replace("HG003", "HG003-\xe9").encode("latin-1"))

    with pytest.raises(PvcEvalInputError, match=r"a.vcf .*not UTF-8.*bcftools view -Oz"):
        parse_vcf_samples(path)
    with pytest.raises(PvcEvalInputError, match=r"a.vcf .*not UTF-8.*bcftools view -Oz"):
        list(parse_vcf(path, window))


def test_truncated_gzip_file_stops_the_run(tmp_path: Path, window: FastaWindow) -> None:
    whole = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n" * 50, gzipped=True)
    path = tmp_path / "cut.vcf.gz"
    path.write_bytes(whole.read_bytes()[:-20])

    with pytest.raises(PvcEvalInputError, match=r"cut.vcf.gz .*gzip"):
        list(parse_vcf(path, window))


@pytest.mark.parametrize(
    ("qual", "gq", "expected", "warning"),
    [
        ("abc", "30", (None, 30), "took 1 QUAL or GQ value that"),
        ("nan", "30", (None, 30), "took 1 QUAL or GQ value that"),
        ("inf", "30", (None, 30), "took 1 QUAL or GQ value that"),
        ("1_000", "30", (None, 30), "took 1 QUAL or GQ value that"),
        ("50", "3.5", (50.0, None), "took 1 QUAL or GQ value that"),
        ("nan", "3.5", (None, None), "took 2 QUAL or GQ values that"),
        (".", ".", (None, None), ""),
    ],
    ids=["qual_word", "qual_nan", "qual_inf", "qual_underscore", "gq_fraction", "both", "dots"],
)
def test_qual_or_gq_that_is_not_a_number_is_absent_and_counted(
    tmp_path: Path,
    window: FastaWindow,
    caplog: pytest.LogCaptureFixture,
    qual: str,
    gq: str,
    expected: tuple[float | None, int | None],
    warning: str,
) -> None:
    # The second sample drops its GQ, which VCF allows and which is not counted.
    path = write_vcf(tmp_path, f"chr1\t5\t.\tA\tT\t{qual}\t.\t.\tGT:GQ\t0/1:{gq}\t0/1\n")

    with caplog.at_level(logging.WARNING):
        (record,) = parse_vcf(path, window)

    genotype = record.genotypes[0]
    assert isinstance(genotype, VcfGenotype)
    assert (record.qual, genotype.quality) == expected
    if warning:
        assert warning in caplog.text
    else:
        assert "QUAL or GQ" not in caplog.text


def test_sample_without_gt_has_both_haplotypes_missing(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGQ:PS\t30:7\t.:.\n")

    (record,) = parse_vcf(path, window)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set="7", quality=30),
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None),
    )


def test_phase_set_dot_is_none(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT:PS\t0|1:.\t1|0:12\n")

    (record,) = parse_vcf(path, window)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1), phased=True, phase_set=None, quality=None),
        VcfGenotype(alleles=(1, 0), phased=True, phase_set="12", quality=None),
    )


def test_alt_dot_gives_a_record_with_ref_only(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\n")

    (record,) = parse_vcf(path, window)

    assert (record.start, record.alleles) == (4, ("A",))


def test_records_at_the_same_pos_both_come_out(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\nchr1\t5\t.\tA\tAC\t.\t.\t.\tGT\t0/0\t0/1\n",
    )

    records = list(parse_vcf(path, window))

    assert [get_placement(record) for record in records] == [
        (4, ("A", "T"), (4, 5), False),
        (4, ("A", "AC"), (4, 6), False),
    ]


def test_bgzipped_file_of_several_gzip_members_is_read_whole(
    tmp_path: Path, window: FastaWindow
) -> None:
    # bgzip writes a file as a series of gzip members; gzip reads them as one stream.
    path = tmp_path / "a.vcf.gz"
    path.write_bytes(
        gzip.compress((HEADER + "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n").encode())
        + gzip.compress(b"chr1\t7\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n")
    )

    records = list(parse_vcf(path, window))

    assert [(record.start, record.alleles) for record in records] == [
        (4, ("A", "T")),
        (6, ("A", "G")),
    ]


def test_telomeric_breakend_at_pos_0_is_skipped_and_counted(
    tmp_path: Path, window: FastaWindow, caplog: pytest.LogCaptureFixture
) -> None:
    # VCF 4.3 writes a breakend before the first base of a contig at POS 0.
    path = write_vcf(
        tmp_path,
        "chr1\t0\tbnd1\tN\t.[chr1:500[\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, window))

    assert [(record.start, record.alleles) for record in records] == [(6, ("A", "T"))]
    assert "skipped 1 VCF record with symbolic alleles" in caplog.text


def test_message_names_the_file_by_the_path_given(tmp_path: Path, window: FastaWindow) -> None:
    (tmp_path / "calls").mkdir()
    path = write_vcf(tmp_path / "calls", "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, window))

    assert str(raised.value).startswith(f"{path}, line 3:")


def test_header_of_one_sample_is_named_in_the_singular(tmp_path: Path, window: FastaWindow) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t5\t.\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    with pytest.raises(
        PvcEvalInputError, match=r"0 sample columns, and the #CHROM line names 1 sample\."
    ):
        list(parse_vcf(path, window))


def test_ref_that_differs_from_the_reference_is_marked_and_spans_its_ref(
    tmp_path: Path, window: FastaWindow
) -> None:
    # The reference has AG at 0-based 4 and 5, and the deletion CA -> C cannot be placed.
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tCA\tC\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tAG\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    mismatch, match = parse_vcf(path, window)

    assert get_placement(mismatch) == (4, ("CA", "C"), (4, 6), True)
    assert get_placement(match) == (4, ("AG", "A"), (4, 6), False)


def test_records_are_trimmed_and_never_moved(tmp_path: Path, window: FastaWindow) -> None:
    # chr1 is GGGAAGAGG...: the deletion of the A at 4, written with shared bases on both
    # sides, is trimmed to AA -> A at 3 and not slid to its left-aligned place.
    path = write_vcf(tmp_path, "chr1\t4\t.\tAAGA\tAGA\t.\t.\t.\tGT\t0/1\t0/1\n")

    (record,) = parse_vcf(path, window)

    assert get_placement(record) == (3, ("AA", "A"), (2, 5), False)


def test_records_come_out_in_order_of_span_start(tmp_path: Path, make_window: MakeWindow) -> None:
    # On GAAAAT, the SNP at 0-based 2 spans 2; the deletion of the A at 4, on the line
    # after it, slides to 1 and spans from 0. Two records with one span start keep the
    # order of the file.
    path = write_vcf(
        tmp_path,
        "chr1\t3\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t4\t.\tA\tAA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )
    # The .fai is what `samtools faidx` (samtools 1.24) wrote for this FASTA.
    window = make_window(">chr1\nGAAAAT\n", "chr1\t6\t6\t6\t7\n")

    records = parse_vcf(path, window)

    assert [get_placement(record) for record in records] == [
        (3, ("AA", "A"), (0, 5), False),
        (3, ("A", "AA"), (0, 6), False),
        (2, ("A", "T"), (2, 3), False),
    ]


def make_snp_line(contig: str, pos_1based: int) -> str:
    return f"{contig}\t{pos_1based}\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/1\n"


# The third line stops the run when it is read, so a record that comes out of the reader
# without the error was given before the reader read further than the second line.
NOT_A_RECORD = "chr1\tfive\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
LONG_REFERENCE = make_fasta({"chr1": "C" * 20000, "chr2": "C" * 10})


@pytest.mark.parametrize(
    ("limit", "second_pos_1based"),
    [({}, 10100), ({"span_limit": 5}, 105)],
    ids=["default_10_kb", "limit_5"],
)
def test_record_comes_out_once_no_later_record_can_start_before_it(
    tmp_path: Path, make_window: MakeWindow, limit: dict[str, int], second_pos_1based: int
) -> None:
    # The SNP spans 99. A record at 0-based 10099 or after cannot span from before 99,
    # since a span reaches at most 10 kb left of its REF; with a limit of 5, one at 104.
    path = write_vcf(
        tmp_path,
        make_snp_line("chr1", 100) + make_snp_line("chr1", second_pos_1based) + NOT_A_RECORD,
    )

    records = parse_vcf(path, make_window(*LONG_REFERENCE), **limit)

    assert get_placement(next(records)) == (99, ("C", "T"), (99, 100), False)


@pytest.mark.parametrize(
    ("limit", "second_pos_1based"),
    [({}, 10099), ({"span_limit": 5}, 104), ({"span_limit": None}, 19000)],
    ids=["default_10_kb", "limit_5", "no_limit"],
)
def test_record_is_held_while_a_later_record_could_start_before_it(
    tmp_path: Path,
    make_window: MakeWindow,
    limit: dict[str, int | None],
    second_pos_1based: int,
) -> None:
    # A record at 0-based 10098 could span from 98, before the SNP's 99; with a limit of
    # 5, one at 103; with no limit, any record of the contig.
    path = write_vcf(
        tmp_path,
        make_snp_line("chr1", 100) + make_snp_line("chr1", second_pos_1based) + NOT_A_RECORD,
    )

    records = parse_vcf(path, make_window(*LONG_REFERENCE), **limit)

    with pytest.raises(PvcEvalInputError, match="POS is 'five'"):
        next(records)


@pytest.mark.parametrize("limit", [{}, {"span_limit": None}], ids=["default_10_kb", "no_limit"])
def test_records_held_come_out_at_the_next_contig(
    tmp_path: Path, make_window: MakeWindow, limit: dict[str, None]
) -> None:
    path = write_vcf(tmp_path, make_snp_line("chr1", 100) + make_snp_line("chr2", 1) + NOT_A_RECORD)

    records = parse_vcf(path, make_window(*LONG_REFERENCE), **limit)

    assert get_placement(next(records)) == (99, ("C", "T"), (99, 100), False)


def test_span_beyond_10_kb_is_the_ref_and_counted(
    tmp_path: Path, make_window: MakeWindow, caplog: pytest.LogCaptureFixture
) -> None:
    # The deletion of the last A of 10002 would span from the G before the first, 10001
    # bases left of its REF; the one of 10001 As in chr2 reaches 10000, the most allowed.
    reference = make_fasta({"chr1": "G" + "A" * 10002 + "T", "chr2": "G" + "A" * 10001 + "T"})
    path = write_vcf(
        tmp_path,
        "chr1\t10002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr2\t10001\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        beyond, within = parse_vcf(path, make_window(*reference))

    assert get_placement(beyond) == (10001, ("AA", "A"), (10001, 10003), False)
    assert get_placement(within) == (10000, ("AA", "A"), (0, 10002), False)
    assert (
        "a.vcf: 1 VCF record, at line 3, could slide more than 10,000 bases along the"
        " reference and was given only its REF as its span;"
    ) in caplog.text


def test_spans_beyond_10_kb_are_counted_in_one_warning_from_the_first(
    tmp_path: Path, make_window: MakeWindow, caplog: pytest.LogCaptureFixture
) -> None:
    homopolymer = "G" + "A" * 10002 + "T"
    reference = make_fasta({"chr1": homopolymer, "chr2": homopolymer})
    path = write_vcf(
        tmp_path,
        "chr1\t10002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr2\t10002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, make_window(*reference)))

    assert [get_placement(record)[2] for record in records] == [(10001, 10003), (10001, 10003)]
    assert len(caplog.records) == 1
    assert (
        "a.vcf: 2 VCF records, the first at line 3, could slide more than 10,000 bases along"
        " the reference and each was given only its REF as its span;"
    ) in caplog.text


@pytest.mark.parametrize(
    ("contig", "line", "expected"),
    [
        # An A inserted after the G slides right past 10001 As, to the T.
        ("G" + "A" * 10001 + "T", "chr1\t1\t.\tG\tGA", (0, ("G", "GA"), (0, 1), False)),
        ("G" + "A" * 9999 + "T", "chr1\t1\t.\tG\tGA", (0, ("G", "GA"), (0, 10001), False)),
        # The A at 10003, deleted after two shared As, slides 10002 bases left to 1, and
        # its anchor, the G at 0, is 10001 bases left of the REF: a slide counted from
        # the REF's start, without the shared bases, would stop one base short of it.
        (
            "G" + "A" * 10003 + "CT",
            "chr1\t10002\t.\tAAAC\tAAC,TAAG",
            (10001, ("AAAC", "AAC", "TAAG"), (10001, 10005), False),
        ),
    ],
    ids=["insertion_beyond_10_kb", "insertion_at_10_kb", "deletion_beyond_10_kb_of_a_long_ref"],
)
def test_span_that_reaches_more_than_10_kb_beyond_its_ref_is_the_ref(
    tmp_path: Path,
    make_window: MakeWindow,
    contig: str,
    line: str,
    expected: tuple[int, tuple[str, ...], tuple[int, int], bool],
) -> None:
    path = write_vcf(tmp_path, f"{line}\t.\t.\t.\tGT\t0/1\t0/1\n")

    (record,) = parse_vcf(path, make_window(*make_fasta({"chr1": contig})))

    assert get_placement(record) == expected


# A run of 6 As in chr1 and of 7 in chr2, each after a G and before a T, with 100 Cs on
# either side, so that the bases a span is looked for in start and end inside the contig.
LIMIT_REFERENCE = make_fasta(
    {
        "chr1": "C" * 100 + "G" + "A" * 6 + "T" + "C" * 100,
        "chr2": "C" * 100 + "G" + "A" * 7 + "T" + "C" * 100,
    }
)
# The deletion of the last A of chr1 slides to the G at 100, 5 bases left of its REF; the
# one of chr2 6 bases. An A inserted after the G of chr2 slides right to the T at 108 and
# spans the T, 8 bases right of its REF.
LIMIT_LINES = (
    "chr1\t106\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
    "chr2\t101\t.\tG\tGA\t.\t.\t.\tGT\t0/1\t0/1\n"
    "chr2\t107\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
)
UNCUT_PLACEMENTS = [
    (105, ("AA", "A"), (100, 107), False),
    (100, ("G", "GA"), (100, 109), False),
    (106, ("AA", "A"), (100, 108), False),
]


@pytest.mark.parametrize("limit", [{}, {"span_limit": None}], ids=["default_10_kb", "no_limit"])
def test_spans_within_the_limit_are_not_cut(
    tmp_path: Path, make_window: MakeWindow, limit: dict[str, None]
) -> None:
    path = write_vcf(tmp_path, LIMIT_LINES)

    records = parse_vcf(path, make_window(*LIMIT_REFERENCE), **limit)

    assert [get_placement(record) for record in records] == UNCUT_PLACEMENTS


def test_span_beyond_a_limit_set_by_the_caller_is_the_ref(
    tmp_path: Path, make_window: MakeWindow, caplog: pytest.LogCaptureFixture
) -> None:
    path = write_vcf(tmp_path, LIMIT_LINES)

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, make_window(*LIMIT_REFERENCE), span_limit=5))

    assert [get_placement(record) for record in records] == [
        (105, ("AA", "A"), (100, 107), False),
        (100, ("G", "GA"), (100, 101), False),
        (106, ("AA", "A"), (106, 108), False),
    ]
    assert (
        "a.vcf: 2 VCF records, the first at line 4, could slide more than 5 bases along the"
        " reference and each was given only its REF as its span;"
    ) in caplog.text


def test_with_no_limit_a_span_is_never_cut(
    tmp_path: Path, make_window: MakeWindow, caplog: pytest.LogCaptureFixture
) -> None:
    # The deletion of the last A of 20002 slides 20001 bases left, to the G at 0.
    reference = make_fasta({"chr1": "G" + "A" * 20002 + "T"})
    path = write_vcf(tmp_path, "chr1\t20002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n")

    with caplog.at_level(logging.WARNING):
        (record,) = parse_vcf(path, make_window(*reference), span_limit=None)

    assert get_placement(record) == (20001, ("AA", "A"), (0, 20003), False)
    assert not caplog.text


def test_a_negative_span_limit_is_a_value_error(tmp_path: Path, window: FastaWindow) -> None:
    path = write_vcf(tmp_path, make_snp_line("chr1", 1))

    with pytest.raises(ValueError, match="span_limit must be 0 or more, not -1"):
        list(parse_vcf(path, window, span_limit=-1))


def make_walk(tmp_path: Path) -> tuple[Path, tuple[str, str]]:
    """A contig of 1 Mb, and a VCF of a SNP or a deletion of one base every 997 bases."""
    rng = random.Random(20261019)
    contig = "".join(rng.choices("ACGT", k=1_000_000))
    lines = []
    for pos in range(500, len(contig) - 500, 997):
        if pos % 2:
            lines.append(f"chr1\t{pos + 1}\t.\t{contig[pos : pos + 2]}\t{contig[pos]}")
        else:
            alt = "A" if contig[pos] != "A" else "C"
            lines.append(f"chr1\t{pos + 1}\t.\t{contig[pos]}\t{alt}")
    body = "".join(f"{line}\t.\t.\t.\tGT\t0/1\t0/1\n" for line in lines)
    return write_vcf(tmp_path, body), make_fasta({"chr1": contig})


def test_with_a_limit_the_window_holds_about_the_limit_and_a_block(
    tmp_path: Path, make_window: MakeWindow
) -> None:
    path, reference = make_walk(tmp_path)
    window = make_window(*reference, block_bases=1000)
    most_held = 0
    num_records = 0
    for _ in parse_vcf(path, window, span_limit=2000):
        held_region = window.held_region
        assert held_region is not None
        most_held = max(most_held, held_region[2] - held_region[1])
        num_records += 1

    assert num_records == 1003
    # What the next records can need, 2001 bases behind the last one read, and the block
    # read ahead of it.
    assert most_held <= 2000 + 2 * 1000


def test_with_no_limit_the_reader_never_evicts_and_cuts_no_span(
    tmp_path: Path, make_window: MakeWindow
) -> None:
    path, reference = make_walk(tmp_path)
    window = make_window(*reference)
    # The caller's window holds the first base, and only the caller evicts it.
    window.get_bases("chr1", 0, 1)

    records = list(parse_vcf(path, window, span_limit=None))

    held_region = window.held_region
    assert held_region is not None
    assert held_region[:2] == ("chr1", 0)
    assert records == list(parse_vcf(path, make_window(*reference), span_limit=10**9))


def test_ref_past_the_end_of_its_contig_differs_from_the_reference(
    tmp_path: Path, window: FastaWindow
) -> None:
    # chr2 has 8 bases: the first REF runs one base past its end, the second starts after.
    path = write_vcf(
        tmp_path,
        "chr2\t8\t.\tGG\tG\t.\t.\t.\tGT\t0/1\t0/1\nchr2\t12\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    records = parse_vcf(path, window)

    assert [get_placement(record) for record in records] == [
        (7, ("GG", "G"), (7, 9), True),
        (11, ("G", "T"), (11, 12), True),
    ]


def test_soft_masked_reference_gives_no_mismatch_and_the_same_spans(
    tmp_path: Path, make_window: MakeWindow
) -> None:
    # The window gives soft-masked bases in upper case; §5.2 compares REF without regard
    # to case. The .fai is the one of the reference in upper case.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )
    window = make_window(">chr1\nGGGaagAGGCAGGGGG\n>chr2\nGGGGGGGG\n", REFERENCE_FAI)

    records = parse_vcf(path, window)

    assert [get_placement(record) for record in records] == [
        (3, ("AA", "A"), (2, 5), False),
        (4, ("A", "T"), (4, 5), False),
    ]


def test_gt_with_an_index_the_record_lacks_or_one_allele_is_kept_as_its_case(
    tmp_path: Path,
    window: FastaWindow,
) -> None:
    # The record has one ALT, so 2 is no allele of it.
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/2\t1\n")

    (record,) = parse_vcf(path, window)

    assert record.genotypes == (MALFORMED_GT, HAPLOID_GT)
