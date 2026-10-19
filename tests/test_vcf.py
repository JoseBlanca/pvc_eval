import gzip
import logging
from pathlib import Path

import pytest

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.vcf import (
    HAPLOID_GT,
    MALFORMED_GT,
    MISSING,
    VcfGenotype,
    VcfRecord,
    parse_genotype,
    read_samples,
    read_vcf_records,
)

CONTIGS = ("chr1", "chr2")
# The records of these tests agree with this reference: chr1 has A at 0-based 3, 4 and 6,
# and CA at 9 and 10.
REFERENCE = {"chr1": "GGGAAGAGGCAGGGGG", "chr2": "GGGGGGGG"}
HEADER = (
    "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG003\n"
)


def placed(record: VcfRecord) -> tuple[int, tuple[str, ...], tuple[int, int], bool]:
    return record.start, record.alleles, (record.span_start, record.span_end), record.ref_mismatch


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
def test_parse_genotype(gt: str, expected: VcfGenotype | MALFORMED_GT | HAPLOID_GT) -> None:
    assert parse_genotype(gt, num_alleles=3, phase_set=None, quality=None) == expected


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "bgzipped"])
def test_records_are_read_with_0_based_starts(tmp_path: Path, gzipped: bool) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t10\t.\tCA\tC,CAA\t50.5\tPASS\t.\tGT:GQ:PS\t1|2:30:7\t0/0:.\n"
        "chr2\t3\t.\tG\tT\t.\tPASS\t.\tGT\t./1\t1/1\n",
        gzipped=gzipped,
    )

    first, second = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert (first.contig, first.start, first.end, first.alleles) == (
        "chr1",
        9,
        11,
        ("CA", "C", "CAA"),
    )
    assert (first.span_start, first.span_end) == (9, 12)
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
    assert read_samples(write_vcf(tmp_path, "")) == ("HG002", "HG003")


def test_unsorted_file_stops_the_run_naming_both_lines(tmp_path: Path) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t20\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t15\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(
        PvcEvalInputError, match=r"chr1:15 \(line 4\) comes after chr1:20 \(line 3\)"
    ):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


def test_contigs_out_of_the_reference_order_stop_the_run(tmp_path: Path) -> None:
    path = write_vcf(
        tmp_path,
        "chr2\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(PvcEvalInputError, match="not sorted"):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


def test_contig_missing_from_the_reference_stops_the_run(tmp_path: Path) -> None:
    path = write_vcf(tmp_path, "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError, match="contig chrUn is not in the reference"):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


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
    tmp_path: Path, caplog: pytest.LogCaptureFixture, alt: str
) -> None:
    path = write_vcf(
        tmp_path,
        f"chr1\t5\t.\tA\t{alt}\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))

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
    tmp_path: Path, line: str, message: str
) -> None:
    path = write_vcf(tmp_path, line)

    with pytest.raises(PvcEvalInputError, match=message):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


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
    tmp_path: Path, line: str, found: str
) -> None:
    path = write_vcf(tmp_path, line)

    with pytest.raises(
        PvcEvalInputError, match=rf"a.vcf, line 3: the line has {found} .* names 2 samples\."
    ):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


def test_sites_only_file_gives_records_without_genotypes(tmp_path: Path) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5\t.\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    (record,) = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert (record.start, record.alleles, record.genotypes) == (4, ("A", "T"), ())


def test_blank_lines_are_skipped(tmp_path: Path) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
        "\n"
        " \t \n"
        "chr1\t7\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    assert [record.start for record in read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)] == [
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
def test_file_without_chrom_line_stops_the_run(tmp_path: Path, text: str, message: str) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(PvcEvalInputError, match=message):
        read_samples(path)
    with pytest.raises(PvcEvalInputError, match=message):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


def test_lower_case_alleles_are_read_in_upper_case_and_normalise_alike(tmp_path: Path) -> None:
    # VCF bases are case-insensitive; the reference is upper case.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\taa\ta\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    lower, upper = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert placed(lower) == placed(upper) == (3, ("AA", "A"), (2, 5), False)


def test_file_that_is_not_utf8_text_stops_the_run(tmp_path: Path) -> None:
    path = tmp_path / "a.vcf"
    # A Latin-1 é in a sample name.
    path.write_bytes(HEADER.replace("HG003", "HG003-\xe9").encode("latin-1"))

    with pytest.raises(PvcEvalInputError, match=r"a.vcf .*not UTF-8.*bcftools view -Oz"):
        read_samples(path)
    with pytest.raises(PvcEvalInputError, match=r"a.vcf .*not UTF-8.*bcftools view -Oz"):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


def test_truncated_gzip_file_stops_the_run(tmp_path: Path) -> None:
    whole = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n" * 50, gzipped=True)
    path = tmp_path / "cut.vcf.gz"
    path.write_bytes(whole.read_bytes()[:-20])

    with pytest.raises(PvcEvalInputError, match=r"cut.vcf.gz .*gzip"):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


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
    caplog: pytest.LogCaptureFixture,
    qual: str,
    gq: str,
    expected: tuple[float | None, int | None],
    warning: str,
) -> None:
    # The second sample drops its GQ, which VCF allows and which is not counted.
    path = write_vcf(tmp_path, f"chr1\t5\t.\tA\tT\t{qual}\t.\t.\tGT:GQ\t0/1:{gq}\t0/1\n")

    with caplog.at_level(logging.WARNING):
        (record,) = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    genotype = record.genotypes[0]
    assert isinstance(genotype, VcfGenotype)
    assert (record.qual, genotype.quality) == expected
    if warning:
        assert warning in caplog.text
    else:
        assert "QUAL or GQ" not in caplog.text


def test_sample_without_gt_has_both_haplotypes_missing(tmp_path: Path) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGQ:PS\t30:7\t.:.\n")

    (record,) = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set="7", quality=30),
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None),
    )


def test_phase_set_dot_is_none(tmp_path: Path) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT:PS\t0|1:.\t1|0:12\n")

    (record,) = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1), phased=True, phase_set=None, quality=None),
        VcfGenotype(alleles=(1, 0), phased=True, phase_set="12", quality=None),
    )


def test_alt_dot_gives_a_record_with_ref_only(tmp_path: Path) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\n")

    (record,) = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert (record.start, record.alleles) == (4, ("A",))


def test_records_at_the_same_pos_both_come_out(tmp_path: Path) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\nchr1\t5\t.\tA\tAC\t.\t.\t.\tGT\t0/0\t0/1\n",
    )

    records = list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))

    assert [placed(record) for record in records] == [
        (4, ("A", "T"), (4, 5), False),
        (4, ("A", "AC"), (4, 6), False),
    ]


def test_bgzipped_file_of_several_gzip_members_is_read_whole(tmp_path: Path) -> None:
    # bgzip writes a file as a series of gzip members; gzip reads them as one stream.
    path = tmp_path / "a.vcf.gz"
    path.write_bytes(
        gzip.compress((HEADER + "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n").encode())
        + gzip.compress(b"chr1\t7\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n")
    )

    records = list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))

    assert [(record.start, record.alleles) for record in records] == [
        (4, ("A", "T")),
        (6, ("A", "G")),
    ]


def test_telomeric_breakend_at_pos_0_is_skipped_and_counted(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # VCF 4.3 writes a breakend before the first base of a contig at POS 0.
    path = write_vcf(
        tmp_path,
        "chr1\t0\tbnd1\tN\t.[chr1:500[\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))

    assert [(record.start, record.alleles) for record in records] == [(6, ("A", "T"))]
    assert "skipped 1 VCF record with symbolic alleles" in caplog.text


def test_message_names_the_file_by_the_path_given(tmp_path: Path) -> None:
    (tmp_path / "calls").mkdir()
    path = write_vcf(tmp_path / "calls", "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))

    assert str(raised.value).startswith(f"{path}, line 3:")


def test_header_of_one_sample_is_named_in_the_singular(tmp_path: Path) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t5\t.\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    with pytest.raises(
        PvcEvalInputError, match=r"0 sample columns, and the #CHROM line names 1 sample\."
    ):
        list(read_vcf_records(path, CONTIGS, REFERENCE.__getitem__))


def test_ref_that_differs_from_the_reference_is_marked_and_spans_its_ref(tmp_path: Path) -> None:
    # The reference has AG at 0-based 4 and 5, and the deletion CA -> C cannot be placed.
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tCA\tC\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tAG\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    mismatch, match = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert placed(mismatch) == (4, ("CA", "C"), (4, 6), True)
    assert placed(match) == (4, ("AG", "A"), (4, 6), False)


def test_records_are_trimmed_and_never_moved(tmp_path: Path) -> None:
    # chr1 is GGGAAGAGG...: the deletion of the A at 4, written with shared bases on both
    # sides, is trimmed to AA -> A at 3 and not slid to its left-aligned place.
    path = write_vcf(tmp_path, "chr1\t4\t.\tAAGA\tAGA\t.\t.\t.\tGT\t0/1\t0/1\n")

    (record,) = read_vcf_records(path, CONTIGS, REFERENCE.__getitem__)

    assert placed(record) == (3, ("AA", "A"), (2, 5), False)


def test_records_come_out_in_order_of_span_start(tmp_path: Path) -> None:
    # On GAAAAT, the SNP at 0-based 2 spans 2; the deletion of the A at 4, on the line
    # after it, slides to 1 and spans from 0. Two records with one span start keep the
    # order of the file.
    path = write_vcf(
        tmp_path,
        "chr1\t3\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t4\t.\tA\tAA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    records = read_vcf_records(path, ("chr1",), {"chr1": "GAAAAT"}.__getitem__)

    assert [placed(record) for record in records] == [
        (3, ("AA", "A"), (0, 5), False),
        (3, ("A", "AA"), (0, 6), False),
        (2, ("A", "T"), (2, 3), False),
    ]


def snp_line(contig: str, pos_1based: int) -> str:
    return f"{contig}\t{pos_1based}\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/1\n"


# The third line stops the run when it is read, so a record that comes out of the reader
# without the error was given before the reader read further than the second line.
NOT_A_RECORD = "chr1\tfive\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
LONG_REFERENCE = {"chr1": "C" * 20000, "chr2": "C" * 10}


def test_record_comes_out_once_no_later_record_can_start_before_it(tmp_path: Path) -> None:
    # The SNP spans 99. A record at 0-based 10099 or after cannot span from before 99,
    # since a span reaches at most 10 kb left of its REF.
    path = write_vcf(tmp_path, snp_line("chr1", 100) + snp_line("chr1", 10100) + NOT_A_RECORD)

    records = read_vcf_records(path, CONTIGS, LONG_REFERENCE.__getitem__)

    assert placed(next(records)) == (99, ("C", "T"), (99, 100), False)


def test_record_is_held_while_a_later_record_could_start_before_it(tmp_path: Path) -> None:
    # A record at 0-based 10098 could span from 98, before the SNP's 99.
    path = write_vcf(tmp_path, snp_line("chr1", 100) + snp_line("chr1", 10099) + NOT_A_RECORD)

    records = read_vcf_records(path, CONTIGS, LONG_REFERENCE.__getitem__)

    with pytest.raises(PvcEvalInputError, match="POS is 'five'"):
        next(records)


def test_records_held_come_out_at_the_next_contig(tmp_path: Path) -> None:
    path = write_vcf(tmp_path, snp_line("chr1", 100) + snp_line("chr2", 1) + NOT_A_RECORD)

    records = read_vcf_records(path, CONTIGS, LONG_REFERENCE.__getitem__)

    assert placed(next(records)) == (99, ("C", "T"), (99, 100), False)


def test_span_beyond_10_kb_is_the_ref_and_counted(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # The deletion of the last A of 10002 would span from the G before the first, 10001
    # bases left of its REF; the one of 10001 As in chr2 reaches 10000, the most allowed.
    reference = {"chr1": "G" + "A" * 10002 + "T", "chr2": "G" + "A" * 10001 + "T"}
    path = write_vcf(
        tmp_path,
        "chr1\t10002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr2\t10001\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        beyond, within = read_vcf_records(path, CONTIGS, reference.__getitem__)

    assert placed(beyond) == (10001, ("AA", "A"), (10001, 10003), False)
    assert placed(within) == (10000, ("AA", "A"), (0, 10002), False)
    assert "1 VCF record could slide more than 10 kb" in caplog.text


def test_soft_masked_reference_gives_no_mismatch_and_the_same_spans(tmp_path: Path) -> None:
    # A FASTA reader may give soft-masked bases in lower case; §5.2 compares REF without
    # regard to case.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )
    soft_masked = {"chr1": "GGGaagAGGCAGGGGG", "chr2": "GGGGGGGG"}

    records = read_vcf_records(path, CONTIGS, soft_masked.__getitem__)

    assert [placed(record) for record in records] == [
        (3, ("AA", "A"), (2, 5), False),
        (4, ("A", "T"), (4, 5), False),
    ]
