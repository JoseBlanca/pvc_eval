import gzip
import logging
from pathlib import Path

import pytest

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.vcf import (
    MALFORMED_GT,
    MISSING,
    VcfGenotype,
    parse_genotype,
    read_samples,
    read_vcf_records,
)

CONTIGS = ("chr1", "chr2")
HEADER = (
    "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG003\n"
)


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
        (".", MALFORMED_GT),
        ("0/1/1", MALFORMED_GT),
        ("0/3", MALFORMED_GT),
        ("A/1", MALFORMED_GT),
    ],
    ids=[
        "unphased",
        "phased",
        "one_missing",
        "both_missing",
        "one_allele",
        "three_alleles",
        "index_beyond_alts",
        "not_a_number",
    ],
)
def test_parse_genotype(gt: str, expected: VcfGenotype | MALFORMED_GT) -> None:
    assert parse_genotype(gt, num_alleles=3, phase_set=None, quality=None) == expected


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "bgzipped"])
def test_records_are_read_with_0_based_starts(tmp_path: Path, gzipped: bool) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t10\t.\tCA\tC,CAA\t50.5\tPASS\t.\tGT:GQ:PS\t1|2:30:7\t0/0:.\n"
        "chr2\t3\t.\tG\tT\t.\tPASS\t.\tGT\t./1\t1/1\n",
        gzipped=gzipped,
    )

    first, second = read_vcf_records(path, CONTIGS)

    assert (first.contig, first.start, first.end, first.alleles) == (
        "chr1",
        9,
        11,
        ("CA", "C", "CAA"),
    )
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
        list(read_vcf_records(path, CONTIGS))


def test_contigs_out_of_the_reference_order_stop_the_run(tmp_path: Path) -> None:
    path = write_vcf(
        tmp_path,
        "chr2\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(PvcEvalInputError, match="not sorted"):
        list(read_vcf_records(path, CONTIGS))


def test_contig_missing_from_the_reference_stops_the_run(tmp_path: Path) -> None:
    path = write_vcf(tmp_path, "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError, match="contig chrUn is not in the reference"):
        list(read_vcf_records(path, CONTIGS))


def test_symbolic_alleles_are_skipped_and_counted(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t6\t.\tA\tA]chr2:20]\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(read_vcf_records(path, CONTIGS))

    assert [record.start for record in records] == [6]
    assert "2 VCF records with symbolic alleles" in caplog.text
