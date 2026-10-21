import gzip
import itertools
import logging
import random
import re
import struct
import subprocess
import sys
import threading
import time
import zlib
from collections.abc import Callable, Iterable, Iterator
from contextlib import ExitStack
from pathlib import Path
from typing import TypedDict

import pytest

from pvc_eval import vcf
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import BLOCK_BASES, FastaWindow, IndexedFasta
from pvc_eval.vcf import (
    MALFORMED_GT,
    MISSING,
    OTHER_PLOIDY_GT,
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

type MakeReference = Callable[[str, str], IndexedFasta]


class SpanLimit(TypedDict, total=False):
    """The span_limit a test gives parse_vcf; with none, parse_vcf takes its default."""

    span_limit: int | None


@pytest.fixture
def make_reference(tmp_path: Path) -> Iterator[MakeReference]:
    """Makes an IndexedFasta of a FASTA and its .fai, written in tmp_path, and closes it at
    the end of the test."""
    numbers = itertools.count()
    with ExitStack() as stack:

        def make(fasta: str, fai: str) -> IndexedFasta:
            path = tmp_path / f"ref{next(numbers)}.fa"
            path.write_text(fasta, encoding="utf-8")
            Path(f"{path}.fai").write_text(fai, encoding="utf-8")
            return stack.enter_context(IndexedFasta(path))

        yield make


@pytest.fixture
def reference(make_reference: MakeReference) -> IndexedFasta:
    return make_reference(REFERENCE_FA, REFERENCE_FAI)


# A request to the window that a reader stuck in a loop would reach; the most the tests
# make is about 3,000.
MAX_WINDOW_REQUESTS = 100_000


class RecordingWindow(FastaWindow):
    """A FastaWindow that keeps where it was asked to evict before and the most bases it
    held after a request, and stops a reader that asks for more than MAX_WINDOW_REQUESTS
    regions."""

    def __init__(self, indexed_fasta: IndexedFasta) -> None:
        super().__init__(indexed_fasta)
        self.evicted_before: list[int] = []
        self.most_held = 0
        self.num_requests = 0

    def get_bases(self, contig: str, start: int, end: int) -> str:
        self.num_requests += 1
        if self.num_requests > MAX_WINDOW_REQUESTS:
            raise RuntimeError(f"the reader asked for more than {MAX_WINDOW_REQUESTS} regions")
        bases = super().get_bases(contig, start, end)
        held_region = self.held_region
        if held_region is not None:
            self.most_held = max(self.most_held, held_region[2] - held_region[1])
        return bases

    def evict_before(self, position: int) -> None:
        self.evicted_before.append(position)
        super().evict_before(position)


@pytest.fixture
def windows(monkeypatch: pytest.MonkeyPatch) -> list[RecordingWindow]:
    """The windows parse_vcf makes during the test, each a RecordingWindow."""
    made: list[RecordingWindow] = []

    def make_window(indexed_fasta: IndexedFasta) -> RecordingWindow:
        window = RecordingWindow(indexed_fasta)
        made.append(window)
        return window

    monkeypatch.setattr(vcf, "FastaWindow", make_window)
    return made


def make_fasta(contigs: dict[str, str]) -> tuple[str, str]:
    """A FASTA with each contig on one line, and the .fai `samtools faidx` writes for it."""
    fasta = ""
    fai = ""
    for name, sequence in contigs.items():
        fasta += f">{name}\n"
        fai += f"{name}\t{len(sequence)}\t{len(fasta)}\t{len(sequence)}\t{len(sequence) + 1}\n"
        fasta += sequence + "\n"
    return fasta, fai


def get_placement(record: VcfRecord) -> tuple[int, tuple[str, ...], tuple[int, int]]:
    return record.start, record.alleles, record.span


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
        ("1", OTHER_PLOIDY_GT),
        ("0", OTHER_PLOIDY_GT),
        ("3", MALFORMED_GT),
        ("0/1/1", OTHER_PLOIDY_GT),
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
    gt: str, expected: VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT
) -> None:
    assert (
        _parse_genotype(gt, num_alleles=3, ploidy=2, phase_set=None, quality=None, gts_seen={})
        == expected
    )


def test_same_gt_with_other_number_of_alleles_is_parsed_again(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    # 0/2 is a genotype where the record has two ALTs and malformed where it has one, so a
    # cache that left the number of alleles out of its key would give both the first.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\tA\tG,T\t.\t.\t.\tGT\t0/2\t0/2\nchr1\t6\t.\tG\tA\t.\t.\t.\tGT\t0/2\t0/2\n",
    )

    first, second = parse_vcf(path, reference)

    assert (
        first.genotypes
        == (VcfGenotype(alleles=(0, 2), phased=False, phase_set=None, quality=None),) * 2
    )
    assert second.genotypes == (MALFORMED_GT, MALFORMED_GT)


def test_samples_with_the_same_gt_share_one_alleles_tuple(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    # The same GT and GQ and another PS, so a genotype cache that left PS out of its key
    # would give the second sample the first one's phase set.
    path = write_vcf(tmp_path, "chr1\t4\t.\tA\tG\t.\tPASS\t.\tGT:GQ:PS\t0|1:30:4\t0|1:30:9\n")

    (record,) = parse_vcf(path, reference)

    first, second = record.genotypes
    assert first == VcfGenotype(alleles=(0, 1), phased=True, phase_set="4", quality=30)
    assert second == VcfGenotype(alleles=(0, 1), phased=True, phase_set="9", quality=30)
    assert isinstance(first, VcfGenotype)
    assert isinstance(second, VcfGenotype)
    assert first.alleles is second.alleles


def test_samples_with_the_same_gt_and_ps_keep_their_own_gq(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    # A genotype cache that left GQ out of its key would give the second sample the first
    # one's quality.
    path = write_vcf(tmp_path, "chr1\t4\t.\tA\tG\t.\tPASS\t.\tGT:GQ:PS\t0|1:30:4\t0|1:20:4\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1), phased=True, phase_set="4", quality=30),
        VcfGenotype(alleles=(0, 1), phased=True, phase_set="4", quality=20),
    )


def test_samples_with_the_same_genotype_share_one_genotype(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t4\t.\tA\tG\t.\tPASS\t.\tGT:GQ:PS\t0|1:30:4\t0|1:30:4\n")

    (record,) = parse_vcf(path, reference)

    first, second = record.genotypes
    assert first == VcfGenotype(alleles=(0, 1), phased=True, phase_set="4", quality=30)
    assert first is second


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "bgzipped"])
def test_records_are_read_with_0_based_starts(
    tmp_path: Path, reference: IndexedFasta, gzipped: bool
) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t10\t.\tCA\tC,CAA\t50.5\tPASS\t.\tGT:GQ:PS\t1|2:30:7\t0/0:.\n"
        "chr2\t3\t.\tG\tT\t.\tPASS\t.\tGT\t./1\t1/1\n",
        gzipped=gzipped,
    )

    first, second = parse_vcf(path, reference)

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


def test_unsorted_file_stops_the_run_naming_both_lines(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t12\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t9\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(
        PvcEvalInputError, match=r"chr1:9 \(line 4\) comes after chr1:12 \(line 3\)"
    ):
        list(parse_vcf(path, reference))


def test_contigs_out_of_the_reference_order_stop_the_run(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(
        tmp_path,
        "chr2\t5\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with pytest.raises(PvcEvalInputError, match="not sorted"):
        list(parse_vcf(path, reference))


def test_contig_missing_from_the_reference_stops_the_run(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError, match="contig chrUn is not in the reference"):
        list(parse_vcf(path, reference))


@pytest.mark.parametrize(
    "alt",
    ["<DEL>", "A]chr2:20]", "[chr2:20[A", "A.", ".A", "T,A.", ".A,T"],
    ids=[
        "symbolic",
        "breakend_right",
        "breakend_left",
        "single_breakend_right",
        "single_breakend_left",
        "single_breakend_among_other_alts",
        "single_breakend_before_other_alts",
    ],
)
def test_symbolic_alleles_are_skipped_and_counted(
    tmp_path: Path, reference: IndexedFasta, caplog: pytest.LogCaptureFixture, alt: str
) -> None:
    path = write_vcf(
        tmp_path,
        f"chr1\t5\t.\tA\t{alt}\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, reference))

    assert [(record.start, record.alleles) for record in records] == [(6, ("A", "T"))]
    assert "skipped 1 VCF record with symbolic alleles" in caplog.text


# §5.2: two records with a symbolic allele whose REF the reference does not have, at chr1:5
# (A in the reference) and reaching past the end of chr1's 16 bases.
SYMBOLIC_LINES_THAT_DIFFER = (
    "chr1\t5\t.\tC\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n"
    f"chr1\t6\t.\t{'G' * 20}\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n"
)


def test_symbolic_record_whose_ref_differs_is_skipped_and_counted(
    tmp_path: Path, reference: IndexedFasta, caplog: pytest.LogCaptureFixture
) -> None:
    path = write_vcf(
        tmp_path, SYMBOLIC_LINES_THAT_DIFFER + "chr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, reference))

    assert [(record.start, record.alleles) for record in records] == [(6, ("A", "T"))]
    assert "skipped 2 VCF records with symbolic alleles" in caplog.text


def test_after_symbolic_records_the_first_plain_record_that_differs_stops_the_run(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(
        tmp_path, SYMBOLIC_LINES_THAT_DIFFER + "chr1\t7\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
    )

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, line 5: the REF of chr1:7 is 'C', and the reference, read in upper case, has"
        " 'A' there. It is the first record read with no symbolic allele, so the VCF looks"
        " made on another reference: compare against the reference it was called on, and"
        " run again."
    )


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
    tmp_path: Path, reference: IndexedFasta, line: str, message: str
) -> None:
    path = write_vcf(tmp_path, line)

    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, reference))


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
    tmp_path: Path, reference: IndexedFasta, line: str, found: str
) -> None:
    path = write_vcf(tmp_path, line)

    with pytest.raises(
        PvcEvalInputError, match=rf"a.vcf, line 3: the line has {found} .* names 2 samples\."
    ):
        list(parse_vcf(path, reference))


def test_sites_only_file_gives_records_without_genotypes(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5\t.\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    (record,) = parse_vcf(path, reference)

    assert (record.start, record.alleles, record.genotypes) == (4, ("A", "T"), ())


def test_blank_lines_are_skipped(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
        "\n"
        " \t \n"
        "chr1\t7\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    assert [record.start for record in parse_vcf(path, reference)] == [
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
    tmp_path: Path, reference: IndexedFasta, text: str, message: str
) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(text, encoding="utf-8")

    with pytest.raises(PvcEvalInputError, match=message):
        parse_vcf_samples(path)
    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, reference))


@pytest.mark.parametrize(
    ("kind", "message"),
    [
        ("missing", r"a\.vcf does not exist\. Check its path and run again\.$"),
        ("directory", r"a\.vcf is a directory, not a VCF file\. Check its path and run again\.$"),
        (
            "unreadable",
            r"a\.vcf cannot be read: permission denied\. Allow reading it and run again\.$",
        ),
    ],
    ids=["missing", "directory", "unreadable"],
)
def test_vcf_that_cannot_be_opened_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, kind: str, message: str
) -> None:
    path = tmp_path / "a.vcf"
    if kind == "directory":
        path.mkdir()
    elif kind == "unreadable":
        path.write_text(HEADER, encoding="utf-8")
        path.chmod(0)

    with pytest.raises(PvcEvalInputError, match=message):
        parse_vcf_samples(path)
    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, reference))


@pytest.mark.parametrize(
    ("chrom_line", "column"),
    [
        ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\t\n", 11),
        ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\tHG002\n", 10),
    ],
    ids=["trailing_tab", "two_tabs"],
)
def test_empty_sample_name_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, chrom_line: str, column: int
) -> None:
    path = tmp_path / "a.vcf"
    path.write_text("##fileformat=VCFv4.3\n" + chrom_line, encoding="utf-8")
    message = (
        rf"a\.vcf, line 2: the #CHROM line has an empty sample name in column {column}; each"
        r" sample needs a name, and the columns are separated by one tab\."
    )

    with pytest.raises(PvcEvalInputError, match=message):
        parse_vcf_samples(path)
    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, reference))


def test_lower_case_alleles_are_read_in_upper_case_and_normalise_alike(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    # VCF bases are case-insensitive; the reference is upper case.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\taa\ta\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    lower, upper = parse_vcf(path, reference)

    assert get_placement(lower) == get_placement(upper) == (3, ("AA", "A"), (2, 5))


def test_file_that_is_not_utf8_text_stops_the_run(tmp_path: Path, reference: IndexedFasta) -> None:
    path = tmp_path / "a.vcf"
    # A Latin-1 é in a sample name.
    path.write_bytes(HEADER.replace("HG003", "HG003-\xe9").encode("latin-1"))

    with pytest.raises(PvcEvalInputError, match=r"a.vcf .*not UTF-8.*bcftools view -Oz"):
        parse_vcf_samples(path)
    with pytest.raises(PvcEvalInputError, match=r"a.vcf .*not UTF-8.*bcftools view -Oz"):
        list(parse_vcf(path, reference))


def test_truncated_gzip_file_stops_the_run(tmp_path: Path, reference: IndexedFasta) -> None:
    whole = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n" * 50, gzipped=True)
    path = tmp_path / "cut.vcf.gz"
    path.write_bytes(whole.read_bytes()[:-20])

    with pytest.raises(PvcEvalInputError, match=r"cut.vcf.gz .*gzip"):
        list(parse_vcf(path, reference))


def test_gzip_file_with_damaged_compressed_data_stops_the_run(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    whole = gzip.compress((HEADER + "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n" * 50).encode())
    path = tmp_path / "damaged.vcf.gz"
    # The first byte of the deflate stream, after the gzip header of 10, given the block type
    # 3, which deflate does not have (RFC 1951, 3.2.3), so that zlib cannot decompress it.
    path.write_bytes(whole[:10] + bytes([whole[10] | 0b110]) + whole[11:])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path} is a gzip file cut short or damaged at or after line 1. Compress or copy it"
        " again, and run again."
    )


# Two records of the reference above, and a change to each that makes one of them a fault of
# the file's content: the second record out of order, or the first with its two sample
# columns made one.
TWO_RECORDS = "chr1\t9\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t12\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
# Each change with the fault it makes, and the line of that fault.
CONTENT_FAULTS = [
    (b"chr1\t12\t", b"chr1\t02\t", r"is not sorted: chr1:2 \(line 4\) comes after chr1:9", 4),
    (
        b"0/1\t0/1\nchr1\t12",
        b"0/1 0/1\nchr1\t12",
        r"line 3: the line has 1 sample column, and the #CHROM line names 2 samples",
        3,
    ),
]
CONTENT_FAULT_IDS = ["not_sorted", "sample_columns"]


def write_damaged_gzip(path: Path, text: str, old: bytes, new: bytes) -> Path:
    """Write text gzipped in stored deflate blocks, which hold the text as it is, with old,
    found once, replaced by new of the same length after gzip has taken its CRC-32: gzip
    gives the changed text, and finds the damage only at the end of the file."""
    compressed = gzip.compress(text.encode(), compresslevel=0, mtime=0)
    assert compressed.count(old) == 1
    assert len(old) == len(new)
    path.write_bytes(compressed.replace(old, new))
    return path


def make_damaged_message(path: Path, at_line: str) -> str:
    """The message of a damaged gzip file; at_line is " at line N", or "" with no line."""
    return (
        f"{path} is damaged: gzip's own check of the whole file fails (`gzip -t {path}` shows"
        f" it), so the fault pvc_eval found{at_line} may come from the damaged bytes rather than"
        " from the file's content. Copy or compress the file again from its source, and run"
        " again."
    )


@pytest.mark.parametrize("num_threads", [1, 4])
@pytest.mark.parametrize(
    ("old", "new", "message", "line_number"), CONTENT_FAULTS, ids=CONTENT_FAULT_IDS
)
def test_damaged_gzip_read_as_a_fault_of_a_line_is_reported_damaged(
    tmp_path: Path,
    reference: IndexedFasta,
    old: bytes,
    new: bytes,
    message: str,
    line_number: int,
    num_threads: int,
) -> None:
    path = write_damaged_gzip(tmp_path / "a.vcf.gz", HEADER + TWO_RECORDS, old, new)

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference, num_threads=num_threads))

    assert str(raised.value) == make_damaged_message(path, f" at line {line_number}")
    # The fault of the line the damage made is kept as the cause.
    assert re.search(message, str(raised.value.__cause__))


@pytest.mark.parametrize("gzipped", [False, True], ids=["plain", "gzipped_whole"])
@pytest.mark.parametrize(
    ("old", "new", "message"),
    [(old, new, message) for old, new, message, _ in CONTENT_FAULTS],
    ids=CONTENT_FAULT_IDS,
)
def test_undamaged_file_with_a_fault_of_a_line_keeps_its_message(
    tmp_path: Path, reference: IndexedFasta, old: bytes, new: bytes, message: str, gzipped: bool
) -> None:
    path = write_vcf(tmp_path, TWO_RECORDS.replace(old.decode(), new.decode()), gzipped=gzipped)

    with pytest.raises(PvcEvalInputError, match=message):
        list(parse_vcf(path, reference))


@pytest.mark.parametrize("num_threads", [1, 4])
def test_damaged_gzip_found_while_threads_parse_samples_is_reported_damaged(
    tmp_path: Path, make_reference: MakeReference, num_threads: int
) -> None:
    fasta, lines = make_threaded_input(random.Random(1))
    # The POS of line 2,503 moved back by 1,000 bases, behind the line before it, once
    # several batches have gone to the threads.
    path = write_damaged_gzip(
        tmp_path / "a.vcf.gz", HEADER + lines, b"\nchr1\t25011\t", b"\nchr1\t24011\t"
    )

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, make_reference(*fasta), num_threads=num_threads))

    assert str(raised.value) == make_damaged_message(path, " at line 2503")


def test_damaged_gzip_read_as_a_fault_of_the_header_is_reported_damaged(
    tmp_path: Path,
) -> None:
    path = write_damaged_gzip(
        tmp_path / "a.vcf.gz", HEADER + TWO_RECORDS, b"##fileformat", b"#Xfileformat"
    )

    with pytest.raises(PvcEvalInputError) as raised:
        parse_vcf_samples(path)

    assert str(raised.value) == make_damaged_message(path, " at line 1")
    assert "line 1: a VCF file has a #CHROM line" in str(raised.value.__cause__)


def test_damaged_gzip_read_as_text_that_is_not_utf8_is_reported_damaged(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_damaged_gzip(tmp_path / "a.vcf.gz", HEADER + TWO_RECORDS, b"HG003", b"HG\xff03")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    # A byte that is not UTF-8 is found somewhere after the last line read, on no line known.
    assert str(raised.value) == make_damaged_message(path, "")
    assert isinstance(raised.value.__cause__, UnicodeDecodeError)


# A VCF of 60 SNPs on a reference of 64 As, gzipped as gzip writes it.
SMALL_VCF = HEADER + "".join(
    f"chr1\t{pos_1based}\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1|0\n" for pos_1based in range(1, 61)
)


def test_every_damaged_byte_of_a_gzip_vcf_that_stops_the_read_is_reported_damaged(
    tmp_path: Path, make_reference: MakeReference
) -> None:
    reference = make_reference(*make_fasta({"chr1": "A" * 64}))
    compressed = gzip.compress(SMALL_VCF.encode(), mtime=0)
    path = tmp_path / "a.vcf.gz"
    messages: set[str] = set()
    num_read_whole = 0
    # Every byte after the gzip header of 10, which gzip does not check whole: a change to
    # its first two makes the file no gzip file, and one to its time or its last two bytes
    # changes nothing that is read.
    for offset in range(10, len(compressed)):
        damaged = bytearray(compressed)
        damaged[offset] ^= 0x40
        path.write_bytes(damaged)
        try:
            records = list(parse_vcf(path, reference))
        except PvcEvalInputError as err:
            # The line numbers are dropped, so that each message counts once.
            messages.add(re.sub(r"line \d+", "line N", str(err)))
            continue
        # A read with no error has had gzip check every byte, so it gives the records.
        assert len(records) == 60
        num_read_whole += 1

    assert messages <= {
        make_damaged_message(path, " at line N"),
        make_damaged_message(path, ""),
        f"{path} is a gzip file cut short or damaged at or after line N. Compress or copy it"
        " again, and run again.",
    }
    # The padding bits after the end of the last deflate block, which zlib does not read.
    assert num_read_whole <= 1


def make_bgzf_block(text: bytes) -> bytes:
    """A BGZF block of text, as bgzip writes one: a gzip member whose header holds the BC
    subfield, which gives the size of the block less 1 (SAM/BAM format specification,
    4.1)."""
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    deflated = compressor.compress(text) + compressor.flush()
    block_size = 18 + len(deflated) + 8
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
    trailer = struct.pack("<II", zlib.crc32(text), len(text))
    return header + struct.pack("<H", block_size - 1) + deflated + trailer


# The empty block that ends a BGZF file, the 28 bytes the specification gives.
BGZF_END = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
BGZF_BLOCKS = [
    make_bgzf_block((HEADER + "chr1\t9\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n").encode()),
    make_bgzf_block(b"chr1\t12\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n"),
]


def test_empty_bgzf_block_is_the_end_block_of_the_specification() -> None:
    assert make_bgzf_block(b"") == BGZF_END


def test_bgzf_file_with_its_end_block_is_read(tmp_path: Path, reference: IndexedFasta) -> None:
    path = tmp_path / "a.vcf.gz"
    path.write_bytes(b"".join(BGZF_BLOCKS) + BGZF_END)

    assert [record.start for record in parse_vcf(path, reference)] == [8, 11]


@pytest.mark.parametrize("num_blocks", [1, 2], ids=["cut_after_block_1", "cut_after_block_2"])
def test_bgzf_file_cut_between_blocks_is_cut_short(
    tmp_path: Path, reference: IndexedFasta, num_blocks: int
) -> None:
    path = tmp_path / "a.vcf.gz"
    path.write_bytes(b"".join(BGZF_BLOCKS[:num_blocks]))
    message = (
        f"{path} is a gzip file cut short or damaged: it is in BGZF, as bgzip and bcftools"
        " write, and lacks the empty block that ends such a file. Copy or compress the file"
        " again from its source, and run again."
    )

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))
    assert str(raised.value) == message
    with pytest.raises(PvcEvalInputError) as raised:
        parse_vcf_samples(path)
    assert str(raised.value) == message


_BYTE_ORDER_MARK_VCF = ("﻿" + HEADER + "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n").encode()


@pytest.mark.parametrize(
    "contents",
    [_BYTE_ORDER_MARK_VCF, gzip.compress(_BYTE_ORDER_MARK_VCF)],
    ids=["plain", "gzipped"],
)
def test_byte_order_mark_at_the_start_is_dropped(
    tmp_path: Path, reference: IndexedFasta, contents: bytes
) -> None:
    path = tmp_path / "a.vcf"
    path.write_bytes(contents)

    assert parse_vcf_samples(path) == ("HG002", "HG003")
    assert [(record.start, record.alleles) for record in parse_vcf(path, reference)] == [
        (4, ("A", "T"))
    ]


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
    reference: IndexedFasta,
    caplog: pytest.LogCaptureFixture,
    qual: str,
    gq: str,
    expected: tuple[float | None, int | None],
    warning: str,
) -> None:
    # The second sample drops its GQ, which VCF allows and which is not counted.
    path = write_vcf(tmp_path, f"chr1\t5\t.\tA\tT\t{qual}\t.\t.\tGT:GQ\t0/1:{gq}\t0/1\n")

    with caplog.at_level(logging.WARNING):
        (record,) = parse_vcf(path, reference)

    genotype = record.genotypes[0]
    assert isinstance(genotype, VcfGenotype)
    assert (record.qual, genotype.quality) == expected
    if warning:
        assert warning in caplog.text
    else:
        assert "QUAL or GQ" not in caplog.text


def test_sample_without_gt_has_both_haplotypes_missing(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGQ:PS\t30:7\t.:.\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set="7", quality=30),
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None),
    )


def test_phase_set_dot_is_none(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT:PS\t0|1:.\t1|0:12\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1), phased=True, phase_set=None, quality=None),
        VcfGenotype(alleles=(1, 0), phased=True, phase_set="12", quality=None),
    )


def test_dropped_trailing_sample_fields_are_absent(tmp_path: Path, reference: IndexedFasta) -> None:
    # VCF lets a sample drop trailing fields: HG002 drops PS, HG003 drops GQ and PS too.
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT:GQ:PS\t1|0:20\t0|1\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(1, 0), phased=True, phase_set=None, quality=20),
        VcfGenotype(alleles=(0, 1), phased=True, phase_set=None, quality=None),
    )


def test_dropped_gt_after_another_key_is_missing(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGQ:GT\t30\t30:0/1\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=30),
        VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=30),
    )


def test_sample_values_beyond_the_format_keys_are_ignored(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0|1:30:4\t0|1\n")

    (record,) = parse_vcf(path, reference)

    assert (
        record.genotypes
        == (VcfGenotype(alleles=(0, 1), phased=True, phase_set=None, quality=None),) * 2
    )


def test_key_named_twice_in_format_takes_its_last_copy(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    # HG003 drops the last copy of GQ, so it has no GQ, although it has a value at the first.
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGQ:GT:GQ\t10:0/1:20\t10:1/1\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=20),
        VcfGenotype(alleles=(1, 1), phased=False, phase_set=None, quality=None),
    )


def test_alt_dot_gives_a_record_with_ref_only(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\n")

    (record,) = parse_vcf(path, reference)

    assert (record.start, record.alleles) == (4, ("A",))


@pytest.mark.parametrize(
    ("alts", "alts_without_dots"),
    [
        ("T,.", "T"),
        (".,T", "T"),
        ("<DEL>,.", "<DEL>"),
        ("T,..", "T"),
        ("T,. ", "T"),
        ("T,.,G", "T,G"),
    ],
    ids=["dot_last", "dot_first", "symbolic_and_dot", "two_dots", "dot_and_space", "dot_between"],
)
def test_alt_dot_among_other_alts_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, alts: str, alts_without_dots: str
) -> None:
    path = write_vcf(tmp_path, f"chr1\t5\t.\tA\t{alts}\t.\t.\t.\tGT\t0/1\t0/0\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, line 3: ALT '{alts}' holds `.` beside other alleles, which VCF does not"
        " allow: `.` stands alone, for a record with no ALT. Write the ALT without the `.`,"
        f" here '{alts_without_dots}', and run again."
    )


@pytest.mark.parametrize("alt", ["..", ". "], ids=["two_dots", "dot_and_space"])
def test_alt_dot_written_with_more_than_one_character_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, alt: str
) -> None:
    path = write_vcf(tmp_path, f"chr1\t5\t.\tA\t{alt}\t.\t.\t.\tGT\t0/1\t0/0\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, line 3: ALT '{alt}' is neither an allele nor `.`, which VCF writes alone,"
        " as one character, for a record with no ALT. Write the ALT as `.`, or as the"
        " alleles it stands for, and run again."
    )


def test_records_at_the_same_pos_both_come_out(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\nchr1\t5\t.\tA\tAC\t.\t.\t.\tGT\t0/0\t0/1\n",
    )

    records = list(parse_vcf(path, reference))

    assert [get_placement(record) for record in records] == [
        (4, ("A", "T"), (4, 5)),
        (4, ("A", "AC"), (4, 6)),
    ]


def test_bgzipped_file_of_several_gzip_members_is_read_whole(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    # bgzip writes a file as a series of gzip members; gzip reads them as one stream.
    path = tmp_path / "a.vcf.gz"
    path.write_bytes(
        gzip.compress((HEADER + "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n").encode())
        + gzip.compress(b"chr1\t7\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n")
    )

    records = list(parse_vcf(path, reference))

    assert [(record.start, record.alleles) for record in records] == [
        (4, ("A", "T")),
        (6, ("A", "G")),
    ]


def test_telomeric_breakend_at_pos_0_is_skipped_and_counted(
    tmp_path: Path, reference: IndexedFasta, caplog: pytest.LogCaptureFixture
) -> None:
    # VCF 4.3 writes a breakend before the first base of a contig at POS 0.
    path = write_vcf(
        tmp_path,
        "chr1\t0\tbnd1\tN\t.[chr1:500[\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr1\t7\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, reference))

    assert [(record.start, record.alleles) for record in records] == [(6, ("A", "T"))]
    assert "skipped 1 VCF record with symbolic alleles" in caplog.text


def test_message_names_the_file_by_the_path_given(tmp_path: Path, reference: IndexedFasta) -> None:
    (tmp_path / "calls").mkdir()
    path = write_vcf(tmp_path / "calls", "chrUn\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value).startswith(f"{path}, line 3:")


def test_header_of_one_sample_is_named_in_the_singular(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = tmp_path / "a.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t5\t.\tA\tT\t.\t.\t.\n",
        encoding="utf-8",
    )

    with pytest.raises(
        PvcEvalInputError, match=r"0 sample columns, and the #CHROM line names 1 sample\."
    ):
        list(parse_vcf(path, reference))


@pytest.mark.parametrize("num_threads", [1, 4])
def test_ref_that_differs_from_the_reference_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, num_threads: int
) -> None:
    # The reference has AG at 0-based 4 and 5; the second line agrees with it.
    path = write_vcf(
        tmp_path,
        "chr1\t5\t.\tCa\tC\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tAG\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )
    records = parse_vcf(path, reference, num_threads=num_threads)
    given: list[VcfRecord] = []

    with pytest.raises(PvcEvalInputError) as raised:
        given.extend(records)

    assert str(raised.value) == (
        f"{path}, line 3: the REF of chr1:5 is 'Ca', and the reference, read in upper case,"
        " has 'AG' there. It is"
        " the first record read with no symbolic allele, so the VCF looks made on another"
        " reference: compare against the reference it was called on, and run again."
    )
    assert given == []


@pytest.mark.parametrize(
    ("reference_base", "ref"),
    [
        ("R", "A"),
        ("R", "g"),
        ("M", "C"),
        ("n", "T"),
        ("A", "N"),
        ("G", "r"),
        ("T", "U"),
        ("K", "B"),
    ],
    ids=[
        "reference_r_against_a",
        "reference_r_against_lower_case_g",
        "reference_m_against_c",
        "reference_lower_case_n_against_t",
        "ref_n_against_a",
        "ref_lower_case_r_against_g",
        "ref_u_against_t",
        "two_codes_that_share_a_base",
    ],
)
def test_iupac_codes_that_share_a_base_agree(
    tmp_path: Path, make_reference: MakeReference, reference_base: str, ref: str
) -> None:
    path = write_vcf(tmp_path, f"chr1\t3\t.\t{ref}\t.\t.\t.\t.\tGT\t0/0\t0/0\n")
    reference = make_reference(*make_fasta({"chr1": f"GG{reference_base}GG"}))

    (record,) = parse_vcf(path, reference)

    assert record.alleles == (ref.upper(),)


@pytest.mark.parametrize(
    ("reference_base", "ref", "found"),
    [
        ("A", "c", "is 'c', and the reference, read in upper case, has 'A' there."),
        ("R", "C", "is 'C', and the reference, read in upper case, has 'R' there."),
        (
            "N",
            "x",
            "is 'x', and the reference, read in upper case, has 'N' there ('X' is no IUPAC code).",
        ),
        (
            "x",
            "X",
            "is 'X', and the reference, read in upper case, has 'X' there ('X' is no IUPAC code).",
        ),
    ],
    ids=["a_against_c", "reference_r_against_c", "no_code_against_n", "no_code_against_itself"],
)
def test_iupac_codes_that_share_no_base_differ(
    tmp_path: Path, make_reference: MakeReference, reference_base: str, ref: str, found: str
) -> None:
    # A letter that is no IUPAC code stands for no base, so it agrees with nothing.
    path = write_vcf(tmp_path, f"chr1\t3\t.\t{ref}\t.\t.\t.\t.\tGT\t0/0\t0/0\n")
    reference = make_reference(*make_fasta({"chr1": f"GG{reference_base}GG"}))

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, line 3: the REF of chr1:3 {found} It is the first record read with no"
        " symbolic allele, so the VCF looks made on another reference: compare against the"
        " reference it was called on, and run again."
    )


@pytest.mark.parametrize(
    ("line", "found"),
    [
        ("chr1\t5\t.\tAC\tA", "is 'AC', and the reference, read in upper case, has 'AG' there."),
        ("chr1\t5\t.\tAGC\tA", "is 'AGC', and the reference, read in upper case, has 'AGA' there."),
    ],
    ids=["second_of_two_bases", "third_of_three_bases"],
)
def test_ref_that_agrees_at_its_first_base_and_differs_later_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, line: str, found: str
) -> None:
    # The reference has AGA at 1-based 5 to 7.
    path = write_vcf(tmp_path, f"{line}\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, line 3: the REF of chr1:5 {found} It is the first record read with no"
        " symbolic allele, so the VCF looks made on another reference: compare against the"
        " reference it was called on, and run again."
    )


# 20 bases of REF are written whole in the message, and one more only where they differ.
LONG_REF_SEQUENCE = "G" + "ACGT" * 6 + "G"


@pytest.mark.parametrize(
    ("ref", "found"),
    [
        (
            "ACGTACGTACGTACGTACGA",
            "is 'ACGTACGTACGTACGTACGA', and the reference, read in upper case, has"
            " 'ACGTACGTACGTACGTACGT' there.",
        ),
        (
            "ACGTACGTACGTACGTACGAA",
            "(21 bases) has 'A' at chr1:21, where the reference, read in upper case, has 'T'.",
        ),
        (
            "ACGTACGTAXGTACGTACGTA",
            "(21 bases) has 'X' at chr1:11, where the reference, read in upper case, has 'C'"
            " ('X' is no IUPAC code).",
        ),
    ],
    ids=["20_bases_written_whole", "21_bases_where_they_differ", "21_bases_no_code"],
)
def test_ref_of_more_than_20_bases_is_named_where_it_first_differs(
    tmp_path: Path, make_reference: MakeReference, ref: str, found: str
) -> None:
    path = write_vcf(tmp_path, f"chr1\t2\t.\t{ref}\t.\t.\t.\t.\tGT\t0/0\t0/0\n")
    reference = make_reference(*make_fasta({"chr1": LONG_REF_SEQUENCE}))

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, line 3: the REF of chr1:2 {found} It is the first record read with no"
        " symbolic allele, so the VCF looks made on another reference: compare against the"
        " reference it was called on, and run again."
    )


def test_records_are_trimmed_and_never_moved(tmp_path: Path, reference: IndexedFasta) -> None:
    # chr1 is GGGAAGAGG...: the deletion of the A at 4, written with shared bases on both
    # sides, is trimmed to AA -> A at 3 and not slid to its left-aligned place.
    path = write_vcf(tmp_path, "chr1\t4\t.\tAAGA\tAGA\t.\t.\t.\tGT\t0/1\t0/1\n")

    (record,) = parse_vcf(path, reference)

    assert get_placement(record) == (3, ("AA", "A"), (2, 5))


def test_records_come_out_in_order_of_span_start(
    tmp_path: Path, make_reference: MakeReference
) -> None:
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
    reference = make_reference(">chr1\nGAAAAT\n", "chr1\t6\t6\t6\t7\n")

    records = parse_vcf(path, reference)

    assert [get_placement(record) for record in records] == [
        (3, ("AA", "A"), (0, 5)),
        (3, ("A", "AA"), (0, 6)),
        (2, ("A", "T"), (2, 3)),
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
    tmp_path: Path, make_reference: MakeReference, limit: SpanLimit, second_pos_1based: int
) -> None:
    # The SNP spans 99. A record at 0-based 10099 or after cannot span from before 99,
    # since a span reaches at most 10 kb left of its REF; with a limit of 5, one at 104.
    path = write_vcf(
        tmp_path,
        make_snp_line("chr1", 100) + make_snp_line("chr1", second_pos_1based) + NOT_A_RECORD,
    )

    records = parse_vcf(path, make_reference(*LONG_REFERENCE), **limit)

    assert get_placement(next(records)) == (99, ("C", "T"), (99, 100))


@pytest.mark.parametrize(
    ("limit", "second_pos_1based"),
    [({}, 10099), ({"span_limit": 5}, 104), ({"span_limit": None}, 19000)],
    ids=["default_10_kb", "limit_5", "no_limit"],
)
def test_record_is_held_while_a_later_record_could_start_before_it(
    tmp_path: Path,
    make_reference: MakeReference,
    limit: SpanLimit,
    second_pos_1based: int,
) -> None:
    # A record at 0-based 10098 could span from 98, before the SNP's 99; with a limit of
    # 5, one at 103; with no limit, any record of the contig.
    path = write_vcf(
        tmp_path,
        make_snp_line("chr1", 100) + make_snp_line("chr1", second_pos_1based) + NOT_A_RECORD,
    )

    records = parse_vcf(path, make_reference(*LONG_REFERENCE), **limit)

    with pytest.raises(PvcEvalInputError, match="POS is 'five'"):
        next(records)


@pytest.mark.parametrize("limit", [{}, {"span_limit": None}], ids=["default_10_kb", "no_limit"])
def test_records_held_come_out_at_the_next_contig(
    tmp_path: Path, make_reference: MakeReference, limit: SpanLimit
) -> None:
    path = write_vcf(tmp_path, make_snp_line("chr1", 100) + make_snp_line("chr2", 1) + NOT_A_RECORD)

    records = parse_vcf(path, make_reference(*LONG_REFERENCE), **limit)

    assert get_placement(next(records)) == (99, ("C", "T"), (99, 100))


def test_span_beyond_10_kb_is_the_ref_and_counted(
    tmp_path: Path, make_reference: MakeReference, caplog: pytest.LogCaptureFixture
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
        beyond, within = parse_vcf(path, make_reference(*reference))

    assert get_placement(beyond) == (10001, ("AA", "A"), (10001, 10003))
    assert get_placement(within) == (10000, ("AA", "A"), (0, 10002))
    assert (
        "a.vcf: 1 VCF record, at line 3, could slide more than 10,000 bases along the"
        " reference and was given only its REF as its span;"
    ) in caplog.text


def test_spans_beyond_10_kb_are_counted_in_one_warning_from_the_first(
    tmp_path: Path, make_reference: MakeReference, caplog: pytest.LogCaptureFixture
) -> None:
    homopolymer = "G" + "A" * 10002 + "T"
    reference = make_fasta({"chr1": homopolymer, "chr2": homopolymer})
    path = write_vcf(
        tmp_path,
        "chr1\t10002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr2\t10002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, make_reference(*reference)))

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
        ("G" + "A" * 10001 + "T", "chr1\t1\t.\tG\tGA", (0, ("G", "GA"), (0, 1))),
        ("G" + "A" * 9999 + "T", "chr1\t1\t.\tG\tGA", (0, ("G", "GA"), (0, 10001))),
        # The A at 10003, deleted after two shared As, slides 10002 bases left to 1, and
        # its anchor, the G at 0, is 10001 bases left of the REF: a slide counted from
        # the REF's start, without the shared bases, would stop one base short of it.
        (
            "G" + "A" * 10003 + "CT",
            "chr1\t10002\t.\tAAAC\tAAC,TAAG",
            (10001, ("AAAC", "AAC", "TAAG"), (10001, 10005)),
        ),
    ],
    ids=["insertion_beyond_10_kb", "insertion_at_10_kb", "deletion_beyond_10_kb_of_a_long_ref"],
)
def test_span_that_reaches_more_than_10_kb_beyond_its_ref_is_the_ref(
    tmp_path: Path,
    make_reference: MakeReference,
    contig: str,
    line: str,
    expected: tuple[int, tuple[str, ...], tuple[int, int], bool],
) -> None:
    path = write_vcf(tmp_path, f"{line}\t.\t.\t.\tGT\t0/1\t0/1\n")

    (record,) = parse_vcf(path, make_reference(*make_fasta({"chr1": contig})))

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
    (105, ("AA", "A"), (100, 107)),
    (100, ("G", "GA"), (100, 109)),
    (106, ("AA", "A"), (100, 108)),
]


@pytest.mark.parametrize("limit", [{}, {"span_limit": None}], ids=["default_10_kb", "no_limit"])
def test_spans_within_the_limit_are_not_cut(
    tmp_path: Path, make_reference: MakeReference, limit: SpanLimit
) -> None:
    path = write_vcf(tmp_path, LIMIT_LINES)

    records = parse_vcf(path, make_reference(*LIMIT_REFERENCE), **limit)

    assert [get_placement(record) for record in records] == UNCUT_PLACEMENTS


def test_span_beyond_a_limit_set_by_the_caller_is_the_ref(
    tmp_path: Path, make_reference: MakeReference, caplog: pytest.LogCaptureFixture
) -> None:
    path = write_vcf(tmp_path, LIMIT_LINES)

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, make_reference(*LIMIT_REFERENCE), span_limit=5))

    assert [get_placement(record) for record in records] == [
        (105, ("AA", "A"), (100, 107)),
        (100, ("G", "GA"), (100, 101)),
        (106, ("AA", "A"), (106, 108)),
    ]
    assert (
        "a.vcf: 2 VCF records, the first at line 4, could slide more than 5 bases along the"
        " reference and each was given only its REF as its span;"
    ) in caplog.text


def test_with_no_limit_a_span_is_never_cut(
    tmp_path: Path, make_reference: MakeReference, caplog: pytest.LogCaptureFixture
) -> None:
    # The deletion of the last A of 20002 slides 20001 bases left, to the G at 0.
    reference = make_fasta({"chr1": "G" + "A" * 20002 + "T"})
    path = write_vcf(tmp_path, "chr1\t20002\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n")

    with caplog.at_level(logging.WARNING):
        (record,) = parse_vcf(path, make_reference(*reference), span_limit=None)

    assert get_placement(record) == (20001, ("AA", "A"), (0, 20003))
    assert not caplog.text


@pytest.mark.parametrize(
    "span_limit", [-1, 0, 1.5, "10", True], ids=["negative", "zero", "fraction", "text", "bool"]
)
def test_a_span_limit_that_is_not_a_whole_number_of_1_or_more_is_a_value_error(
    tmp_path: Path, reference: IndexedFasta, span_limit: object
) -> None:
    path = write_vcf(tmp_path, make_snp_line("chr1", 1))

    with pytest.raises(
        ValueError,
        match=rf"^span_limit must be a whole number of 1 or more, or None, not {span_limit!r}$",
    ):
        list(parse_vcf(path, reference, span_limit=span_limit))  # ty: ignore[invalid-argument-type]  # a wrong type on purpose


def test_span_limit_of_1_is_named_in_the_singular(
    tmp_path: Path, make_reference: MakeReference, caplog: pytest.LogCaptureFixture
) -> None:
    path = write_vcf(tmp_path, LIMIT_LINES)

    with caplog.at_level(logging.WARNING):
        list(parse_vcf(path, make_reference(*LIMIT_REFERENCE), span_limit=1))

    assert (
        "a.vcf: 3 VCF records, the first at line 3, could slide more than 1 base along the"
        " reference and each was given only its REF as its span;"
    ) in caplog.text


# A run of 200 As after a G: the deletion of the A at 2 slides left to the G and right to
# the last A, past the first bases its span is looked for in.
LONG_SLIDE_REFERENCE = make_fasta({"chr1": "G" + "A" * 200 + "T" + "C" * 300})
LONG_SLIDE_LINE = "chr1\t2\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\n"


@pytest.mark.parametrize("limit", [{}, {"span_limit": None}], ids=["default_10_kb", "no_limit"])
def test_deletion_that_slides_far_right_spans_its_whole_slide(
    tmp_path: Path, make_reference: MakeReference, limit: SpanLimit
) -> None:
    path = write_vcf(tmp_path, LONG_SLIDE_LINE)

    (record,) = parse_vcf(path, make_reference(*LONG_SLIDE_REFERENCE), **limit)

    assert get_placement(record) == (1, ("AA", "A"), (0, 201))


def test_span_reaching_one_base_beyond_a_limit_of_64_is_the_ref(
    tmp_path: Path, make_reference: MakeReference
) -> None:
    # The deletion's span ends 198 bases after its REF, and is first looked for 64 bases
    # on each side of it: the most a limit of 64 allows.
    path = write_vcf(tmp_path, LONG_SLIDE_LINE)

    (record,) = parse_vcf(path, make_reference(*LONG_SLIDE_REFERENCE), span_limit=64)

    assert get_placement(record) == (1, ("AA", "A"), (1, 3))


def test_with_no_limit_an_insertion_slides_to_the_end_of_its_contig(
    tmp_path: Path, make_reference: MakeReference, windows: list[RecordingWindow]
) -> None:
    # The window stops a reader that keeps asking for bases past the contig's end.
    path = write_vcf(tmp_path, "chr1\t1\t.\tG\tGA\t.\t.\t.\tGT\t0/1\t0/1\n")
    reference = make_reference(*make_fasta({"chr1": "G" + "A" * 100}))

    (record,) = parse_vcf(path, reference, span_limit=None)

    assert get_placement(record) == (0, ("G", "GA"), (0, 101))
    assert len(windows) == 1


def test_a_byte_that_is_not_a_base_where_a_record_reads_stops_the_run(
    tmp_path: Path, make_reference: MakeReference
) -> None:
    path = write_vcf(tmp_path, "chr1\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n")
    reference = make_reference(">chr1\nGGGA1GAGG\n", "chr1\t9\t6\t9\t10\n")

    with pytest.raises(PvcEvalInputError, match=r"ref\d+\.fa: chr1:5 is not a letter"):
        list(parse_vcf(path, reference))


def make_walk(tmp_path: Path) -> tuple[Path, tuple[str, str], list[int]]:
    """A contig of 1 Mb, a VCF of a SNP or a deletion of one base every 997 bases, and the
    0-based POS of each record."""
    rng = random.Random(20261019)
    contig = "".join(rng.choices("ACGT", k=1_000_000))
    positions = list(range(500, len(contig) - 500, 997))
    lines = []
    for pos in positions:
        if pos % 2:
            lines.append(f"chr1\t{pos + 1}\t.\t{contig[pos : pos + 2]}\t{contig[pos]}")
        else:
            alt = "A" if contig[pos] != "A" else "C"
            lines.append(f"chr1\t{pos + 1}\t.\t{contig[pos]}\t{alt}")
    body = "".join(f"{line}\t.\t.\t.\tGT\t0/1\t0/1\n" for line in lines)
    return write_vcf(tmp_path, body), make_fasta({"chr1": contig}), positions


def test_with_a_limit_the_window_holds_the_limit_behind_and_a_block_ahead(
    tmp_path: Path, make_reference: MakeReference, windows: list[RecordingWindow]
) -> None:
    path, fasta, positions = make_walk(tmp_path)

    records = list(parse_vcf(path, make_reference(*fasta), span_limit=2000))

    assert len(records) == len(positions) == 1003
    (window,) = windows
    # Before each record, the bases no record from its POS on can need: its span is
    # looked for at most span_limit + 1 bases before its POS.
    assert window.evicted_before == [pos - 2001 for pos in positions]
    # The 2001 bases behind the last record read, the record, and a block read ahead.
    assert window.most_held <= 2001 + 100 + BLOCK_BASES


def test_with_no_limit_the_window_holds_a_block_behind_and_a_block_ahead(
    tmp_path: Path, make_reference: MakeReference, windows: list[RecordingWindow]
) -> None:
    path, fasta, positions = make_walk(tmp_path)
    reference = make_reference(*fasta)

    records = list(parse_vcf(path, reference, span_limit=None))

    window = windows[0]
    # Before each record, a block behind its POS; a slide that reaches further is read
    # again.
    assert window.evicted_before == [pos - BLOCK_BASES for pos in positions]
    assert window.most_held <= 2 * BLOCK_BASES + 100
    assert records == list(parse_vcf(path, reference, span_limit=10**9))


def test_each_stream_reads_through_a_window_of_its_own(
    tmp_path: Path, make_reference: MakeReference, windows: list[RecordingWindow]
) -> None:
    path, fasta, _ = make_walk(tmp_path)
    reference = make_reference(*fasta)

    pairs = list(zip(parse_vcf(path, reference), parse_vcf(path, reference), strict=True))

    assert len(windows) == 2
    assert all(first == second for first, second in pairs)


@pytest.mark.parametrize(
    ("line", "where"),
    [
        ("chr2\t8\t.\tGG\tG", "line 3: the REF of chr2:8 reaches past"),
        ("chr2\t12\t.\tG\tT", "line 3: the REF of chr2:12 reaches past"),
    ],
    ids=["runs_over_the_end", "starts_past_the_end"],
)
def test_ref_past_the_end_of_its_contig_stops_the_run(
    tmp_path: Path, reference: IndexedFasta, line: str, where: str
) -> None:
    # chr2 has 8 bases.
    path = write_vcf(tmp_path, f"{line}\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, reference))

    assert str(raised.value) == (
        f"{path}, {where} the end of chr2, which has 8 bases in the reference, so the VCF"
        " looks made on another reference. Compare against the reference the VCF was"
        " called on."
    )


def test_ref_on_the_last_base_of_its_contig_is_read(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr2\t7\t.\tGG\tG\t.\t.\t.\tGT\t0/1\t0/1\n")

    (record,) = parse_vcf(path, reference)

    assert get_placement(record) == (6, ("GG", "G"), (0, 8))


def test_soft_masked_reference_differs_from_no_ref_and_gives_the_same_spans(
    tmp_path: Path, make_reference: MakeReference
) -> None:
    # The window gives soft-masked bases in upper case; §5.2 compares REF without regard
    # to case. The .fai is the one of the reference in upper case.
    path = write_vcf(
        tmp_path,
        "chr1\t4\t.\tAA\tA\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
    )
    reference = make_reference(">chr1\nGGGaagAGGCAGGGGG\n>chr2\nGGGGGGGG\n", REFERENCE_FAI)

    records = parse_vcf(path, reference)

    assert [get_placement(record) for record in records] == [
        (3, ("AA", "A"), (2, 5)),
        (4, ("A", "T"), (4, 5)),
    ]


def test_gt_with_an_index_the_record_lacks_or_one_allele_is_kept_as_its_case(
    tmp_path: Path,
    reference: IndexedFasta,
) -> None:
    # The record has one ALT, so 2 is no allele of it.
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/2\t1\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (MALFORMED_GT, OTHER_PLOIDY_GT)


def test_one_allele_read_with_ploidy_2_is_other_ploidy(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t1\t0/1\n")

    (record,) = parse_vcf(path, reference, ploidy=2)

    assert record.genotypes == (
        OTHER_PLOIDY_GT,
        VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=None),
    )


def test_one_allele_read_with_ploidy_1_is_a_genotype(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t1\t0\n")

    (record,) = parse_vcf(path, reference, ploidy=1)

    assert record.genotypes == (
        VcfGenotype(alleles=(1,), phased=False, phase_set=None, quality=None),
        VcfGenotype(alleles=(0,), phased=False, phase_set=None, quality=None),
    )


def test_two_alleles_read_with_ploidy_1_are_other_ploidy(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1\n")

    (record,) = parse_vcf(path, reference, ploidy=1)

    assert record.genotypes == (
        OTHER_PLOIDY_GT,
        VcfGenotype(alleles=(1,), phased=False, phase_set=None, quality=None),
    )


def test_two_dots_read_with_ploidy_1_are_missing(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t./.\t.\n")

    (record,) = parse_vcf(path, reference, ploidy=1)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING,), phased=False, phase_set=None, quality=None),
        VcfGenotype(alleles=(MISSING,), phased=False, phase_set=None, quality=None),
    )


def test_sample_without_gt_read_with_ploidy_1_has_its_haplotype_missing(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGQ:PS\t30:7\t.:.\n")

    (record,) = parse_vcf(path, reference, ploidy=1)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING,), phased=False, phase_set="7", quality=30),
        VcfGenotype(alleles=(MISSING,), phased=False, phase_set=None, quality=None),
    )


def test_one_dot_read_with_ploidy_2_is_missing_on_both_haplotypes(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t.\t./.\n")

    (record,) = parse_vcf(path, reference, ploidy=2)

    assert record.genotypes == (
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None),
        VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None),
    )


def test_three_alleles_read_with_ploidy_3_are_a_genotype(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1/1\t./1/0\n")

    (record,) = parse_vcf(path, reference, ploidy=3)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1, 1), phased=False, phase_set=None, quality=None),
        VcfGenotype(alleles=(MISSING, 1, 0), phased=False, phase_set=None, quality=None),
    )


def test_index_beyond_two_alts_is_malformed(tmp_path: Path, reference: IndexedFasta) -> None:
    path = write_vcf(tmp_path, "chr1\t5\t.\tA\tT,C\t.\t.\t.\tGT\t3/1\t2/1\n")

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        MALFORMED_GT,
        VcfGenotype(alleles=(2, 1), phased=False, phase_set=None, quality=None),
    )


@pytest.mark.parametrize("ploidy", [0, -1, 1.5, True], ids=["zero", "negative", "fraction", "bool"])
def test_a_ploidy_that_is_not_a_whole_number_of_1_or_more_is_a_value_error(
    tmp_path: Path, reference: IndexedFasta, ploidy: object
) -> None:
    path = write_vcf(tmp_path, make_snp_line("chr1", 1))
    # The call gives the stream without checking; the check is made when it starts.
    records = parse_vcf(path, reference, ploidy=ploidy)  # ty: ignore[invalid-argument-type]  # a wrong type on purpose

    with pytest.raises(
        ValueError, match=rf"^ploidy must be a whole number of 1 or more, not {ploidy!r}$"
    ):
        next(records)


# More lines than 4 threads take in their first batches (2 each of vcf._BATCH_SIZE), so
# that the pool has several batches for every thread.
NUM_THREADED_LINES = 3000


def make_threaded_input(rng: random.Random) -> tuple[tuple[str, str], str]:
    """A reference of random bases, and VCF lines on it, one every 10 bases: SNPs,
    insertions, deletions and records of two ALTs, a symbolic record every 97 lines, and
    GTs of each kind of §5.4 with a GQ, sometimes not a number, and a PS."""
    sequence = "".join(rng.choice("ACGT") for _ in range(NUM_THREADED_LINES * 10 + 20))
    gts = ["0/0", "0/1", "1|0", "1/1", "./.", ".", "0/2", "1", "0/1/1"]
    lines = ""
    for index in range(NUM_THREADED_LINES):
        pos_1based = 11 + index * 10
        ref = sequence[pos_1based - 1]
        other = "T" if ref != "T" else "G"
        alt = [other, ref + "A", f"{other},{ref}C"][index % 3]
        if index % 7 == 0:
            ref, alt = sequence[pos_1based - 1 : pos_1based + 1], ref
        if index % 97 == 0:
            alt = "<DEL>"
        samples = "\t".join(
            f"{rng.choice(gts)}:{rng.choice([str(rng.randrange(100)), 'x', '.'])}:"
            f"{rng.randrange(5)}"
            for _ in range(2)
        )
        lines += f"chr1\t{pos_1based}\t.\t{ref}\t{alt}\t30\t.\t.\tGT:GQ:PS\t{samples}\n"
    return make_fasta({"chr1": sequence}), lines


def test_records_and_warnings_are_the_same_with_one_thread_or_several(
    tmp_path: Path,
    make_reference: MakeReference,
    caplog: pytest.LogCaptureFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta, lines = make_threaded_input(random.Random(20261019))
    path = write_vcf(tmp_path, lines)
    reference = make_reference(*fasta)
    # The thread that parsed each batch, by the number of threads.
    batch_threads: dict[int, list[threading.Thread]] = {1: [], 4: []}
    parse_batch = vcf._parse_batch

    def record_thread(
        placed_lines: list[vcf._PlacedLine], ploidy: int, seen_by_thread: threading.local
    ) -> tuple[list[tuple[int, VcfRecord]], int]:
        batch_threads[num_threads].append(threading.current_thread())
        return parse_batch(placed_lines, ploidy, seen_by_thread)

    monkeypatch.setattr(vcf, "_parse_batch", record_thread)
    records: dict[int, list[VcfRecord]] = {}
    warnings: dict[int, list[str]] = {}
    for num_threads in (1, 4):
        caplog.clear()
        with caplog.at_level(logging.WARNING):
            records[num_threads] = list(parse_vcf(path, reference, num_threads=num_threads))
        warnings[num_threads] = caplog.messages

    assert records[4] == records[1]
    assert warnings[4] == warnings[1]
    assert len(records[1]) == NUM_THREADED_LINES - len(range(0, NUM_THREADED_LINES, 97))
    assert len(warnings[1]) == 2
    # One thread parses every batch itself; four leave every batch to the pool.
    assert set(batch_threads[1]) == {threading.current_thread()}
    assert threading.current_thread() not in batch_threads[4]
    assert len(batch_threads[4]) > 4 * 4


# More lines than 4 threads take in their first batches, as NUM_THREADED_LINES is.
NUM_PLOIDY_LINES = 1000


@pytest.mark.parametrize("num_threads", [1, 4])
def test_two_streams_of_other_ploidies_read_together_keep_their_own_genotypes(
    tmp_path: Path, make_reference: MakeReference, num_threads: int
) -> None:
    # The same GTs read with ploidy 1 and 2 by two streams open at once in one process,
    # taken in turns on this thread, each with its own pool, so a GT or genotype kept for
    # one ploidy in anything the streams share and given to the other would show here.
    lines = "".join(
        f"chr1\t{pos_1based}\t.\tC\tT\t.\t.\t.\tGT\t1\t0/1\n"
        if pos_1based % 2
        else f"chr1\t{pos_1based}\t.\tC\tT\t.\t.\t.\tGT\t.\t./.\n"
        for pos_1based in range(1, NUM_PLOIDY_LINES + 1)
    )
    path = write_vcf(tmp_path, lines)
    reference = make_reference(*make_fasta({"chr1": "C" * NUM_PLOIDY_LINES}))

    pairs = list(
        zip(
            parse_vcf(path, reference, num_threads=num_threads, ploidy=1),
            parse_vcf(path, reference, num_threads=num_threads, ploidy=2),
            strict=True,
        )
    )

    haploid = VcfGenotype(alleles=(1,), phased=False, phase_set=None, quality=None)
    diploid = VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=None)
    missing_1 = VcfGenotype(alleles=(MISSING,), phased=False, phase_set=None, quality=None)
    missing_2 = VcfGenotype(alleles=(MISSING, MISSING), phased=False, phase_set=None, quality=None)
    assert len(pairs) == NUM_PLOIDY_LINES
    for record_1, record_2 in pairs:
        if record_1.start % 2 == 0:
            assert record_1.genotypes == (haploid, OTHER_PLOIDY_GT)
            assert record_2.genotypes == (OTHER_PLOIDY_GT, diploid)
        else:
            assert record_1.genotypes == (missing_1, missing_1)
            assert record_2.genotypes == (missing_2, missing_2)


@pytest.mark.parametrize("num_threads", [1, 4])
def test_with_threads_an_error_comes_after_the_records_before_it(
    tmp_path: Path, make_reference: MakeReference, num_threads: int
) -> None:
    fasta, lines = make_threaded_input(random.Random(1))
    path = write_vcf(tmp_path, lines + NOT_A_RECORD)
    records = parse_vcf(path, make_reference(*fasta), span_limit=5, num_threads=num_threads)
    given: list[VcfRecord] = []

    with pytest.raises(PvcEvalInputError, match="POS is 'five'"):
        given.extend(records)

    # Every record but the last comes out before the error: the next record, 10 bases on,
    # shows that no later record can start before it with a limit of 5. 31 lines are
    # symbolic.
    assert len(given) == NUM_THREADED_LINES - 31 - 1
    assert given[-1].line_number == NUM_THREADED_LINES + 1


def make_snps_with_one_ref_wrong(index_wrong: int) -> tuple[tuple[str, str], str]:
    """A reference of random bases, and a SNP every 10 bases on it, NUM_THREADED_LINES of
    them, from 1-based 11; the REF of the SNP of 0-based index index_wrong is T where the
    reference has A, and every other REF agrees with the reference."""
    rng = random.Random(20261019)
    sequence = [rng.choice("ACGT") for _ in range(NUM_THREADED_LINES * 10 + 20)]
    lines = ""
    for index in range(NUM_THREADED_LINES):
        pos_1based = 11 + index * 10
        if index == index_wrong:
            sequence[pos_1based - 1] = "A"
        ref = "T" if index == index_wrong else sequence[pos_1based - 1]
        alt = "G" if ref != "G" else "C"
        lines += f"chr1\t{pos_1based}\t.\t{ref}\t{alt}\t30\t.\t.\tGT\t0/1\t1/1\n"
    return make_fasta({"chr1": "".join(sequence)}), lines


@pytest.mark.parametrize("num_threads", [1, 4])
@pytest.mark.parametrize(
    "index_wrong",
    [0, 1, 100, vcf._BATCH_SIZE - 1, vcf._BATCH_SIZE, 5 * vcf._BATCH_SIZE],
    ids=[
        "first_record",
        "second_record",
        "middle_of_a_batch",
        "last_of_a_batch",
        "first_of_a_batch",
        "first_of_a_batch_past_the_first_of_each_thread",
    ],
)
def test_ref_that_differs_comes_after_the_records_before_it(
    tmp_path: Path, make_reference: MakeReference, index_wrong: int, num_threads: int
) -> None:
    fasta, lines = make_snps_with_one_ref_wrong(index_wrong)
    path = write_vcf(tmp_path, lines)
    records = parse_vcf(path, make_reference(*fasta), span_limit=5, num_threads=num_threads)
    given: list[VcfRecord] = []

    with pytest.raises(PvcEvalInputError) as raised:
        given.extend(records)

    # The header takes two lines. Every record before the wrong one but the last comes out
    # before the error: the record after it, 10 bases on, would show that no later record
    # can start before it with a limit of 5, as for every other input error.
    line_number = index_wrong + 3
    pos_1based = 11 + index_wrong * 10
    if index_wrong == 0:
        remedy = (
            "It is the first record read with no symbolic allele, so the VCF looks made on"
            " another reference: compare against the reference it was called on, and run"
            " again."
        )
    elif index_wrong == 1:
        remedy = (
            "The 1 record before it with no symbolic allele agrees with the reference, so"
            " this line looks damaged: correct it or remove it, and run again."
        )
    else:
        remedy = (
            f"The {index_wrong} records before it with no symbolic allele agree with the"
            " reference, so this line looks damaged: correct it or remove it, and run again."
        )
    assert str(raised.value) == (
        f"{path}, line {line_number}: the REF of chr1:{pos_1based} is 'T', and the"
        f" reference, read in upper case, has 'A' there. {remedy}"
    )
    assert [record.line_number for record in given] == list(range(3, line_number - 1))
    assert [record.start for record in given] == [10 + 10 * index for index in range(len(given))]


@pytest.mark.parametrize(
    "num_threads",
    [0, -1, 257, 2**62, 1.5, "4", True],
    ids=["zero", "negative", "above_256", "huge", "fraction", "text", "bool"],
)
def test_a_num_threads_that_is_not_a_whole_number_from_1_to_256_is_a_value_error(
    tmp_path: Path, reference: IndexedFasta, num_threads: object
) -> None:
    path = write_vcf(tmp_path, make_snp_line("chr1", 1))

    with pytest.raises(
        ValueError,
        match=rf"^num_threads must be a whole number from 1 to 256, not {num_threads!r}$",
    ):
        list(parse_vcf(path, reference, num_threads=num_threads))  # ty: ignore[invalid-argument-type]  # a wrong type on purpose


def test_a_pool_of_4_threads_runs_4_items_at_once_and_gives_them_in_order() -> None:
    # Each item waits until 4 are running, so a pool of fewer threads breaks the barrier.
    barrier = threading.Barrier(4, timeout=5)

    def wait_for_4(item: int) -> int:
        barrier.wait()
        return item

    assert list(vcf._map_in_order(wait_for_4, range(8), 4)) == list(range(8))


def make_placed_line(gt: str, gq: int) -> vcf._PlacedLine:
    """A SNP at the first base of chr1 with one sample, of this GT and GQ."""
    fields = ["chr1", "1", ".", "A", "C", "30", ".", ".", "GT:GQ", f"{gt}:{gq}"]
    return vcf._PlacedLine(
        pos=0,
        contig="chr1",
        start=0,
        alleles=("A", "C"),
        qual=30.0,
        line_number=gq,
        span=(0, 1),
        fields=fields,
    )


def test_genotypes_and_gts_a_thread_keeps_are_bounded() -> None:
    # Every genotype has a GQ of its own, so none is seen twice, and the thread passes its
    # bound of 4,096 twice in 10,000 lines. The bound is checked when a batch starts, so a
    # thread holds at most one batch more. The first batch's GT, 1/1, is not seen again.
    seen_by_thread = threading.local()
    batch_size = 64
    vcf._parse_batch([make_placed_line("1/1", gq) for gq in range(batch_size)], 2, seen_by_thread)
    most_held = 0
    for start in range(batch_size, 10_000, batch_size):
        lines = [make_placed_line("0/1", gq) for gq in range(start, start + batch_size)]
        vcf._parse_batch(lines, 2, seen_by_thread)
        most_held = max(most_held, len(seen_by_thread.genotypes))

    assert most_held == 4096 + batch_size
    assert len(seen_by_thread.genotypes) < 4096
    assert set(seen_by_thread.gts) == {("0/1", 2, 2)}


def test_genotypes_a_thread_keeps_for_one_ploidy_are_not_given_for_another() -> None:
    # One thread's store, given batches of one GT with ploidy 1 and then 2: a key without
    # the ploidy would give the second batch the first one's genotype.
    seen_by_thread = threading.local()

    (_, haploid_record), *_ = vcf._parse_batch([make_placed_line("1", 30)], 1, seen_by_thread)[0]
    (_, diploid_record), *_ = vcf._parse_batch([make_placed_line("1", 30)], 2, seen_by_thread)[0]

    assert haploid_record.genotypes == (
        VcfGenotype(alleles=(1,), phased=False, phase_set=None, quality=30),
    )
    assert diploid_record.genotypes == (OTHER_PLOIDY_GT,)


# More digits than Python's int() reads from text by default (4,300).
TOO_MANY_DIGITS = "1" * 5000


def test_pos_of_more_digits_than_python_reads_stops_the_run(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(tmp_path, f"chr1\t{TOO_MANY_DIGITS}\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n")

    with pytest.raises(
        PvcEvalInputError, match=r"a\.vcf, line 3: POS has 5,000 digits, which is past the end"
    ):
        list(parse_vcf(path, reference))


def test_allele_index_of_more_digits_than_python_reads_is_malformed(
    tmp_path: Path, reference: IndexedFasta
) -> None:
    path = write_vcf(
        tmp_path, f"chr1\t4\t.\tA\tG\t.\t.\t.\tGT\t0/{TOO_MANY_DIGITS}\t{TOO_MANY_DIGITS}\n"
    )

    (record,) = parse_vcf(path, reference)

    assert record.genotypes == (MALFORMED_GT, MALFORMED_GT)


def test_gq_of_more_digits_than_python_reads_is_absent_and_counted(
    tmp_path: Path, reference: IndexedFasta, caplog: pytest.LogCaptureFixture
) -> None:
    path = write_vcf(tmp_path, f"chr1\t4\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:{TOO_MANY_DIGITS}\t0/1:20\n")

    with caplog.at_level(logging.WARNING):
        (record,) = parse_vcf(path, reference)

    assert record.genotypes == (
        VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=None),
        VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=20),
    )
    assert "took 1 QUAL or GQ value that could not be read as a number as absent" in caplog.text


# doc/specs/regions_read.md: parse_vcf reads only the records of the contigs it is given.
# A reference of four contigs, in this .fai order; chr10's name begins with chr1's. Every
# line of CONTIG_LINES agrees with it.
FOUR_CONTIGS = {"chr1": "ACGTACGTAC", "chr2": "GGGGGGGG", "chr10": "TTTTTTTT", "chr3": "CCCCCCCC"}
CONTIG_LINES = (
    "chr1\t2\t.\tC\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
    "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n"
    "chr2\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n"
    "chr10\t4\t.\tT\tC\t.\t.\t.\tGT\t0/1\t0/1\n"
    "chr3\t6\t.\tC\tG\t.\t.\t.\tGT\t0/1\t0/1\n"
)


@pytest.fixture
def four_contigs(make_reference: MakeReference) -> IndexedFasta:
    return make_reference(*make_fasta(FOUR_CONTIGS))


def get_places(records: Iterable[VcfRecord]) -> list[tuple[str, int, int]]:
    """The contig, 0-based start and line number of each record."""
    return [(record.contig, record.start, record.line_number) for record in records]


@pytest.mark.parametrize(
    ("contigs", "places"),
    [
        ({"chr1"}, [("chr1", 1, 3), ("chr1", 4, 4)]),
        ({"chr10"}, [("chr10", 3, 6)]),
        (["chr2", "chr3"], [("chr2", 2, 5), ("chr3", 5, 7)]),
        (
            None,
            [("chr1", 1, 3), ("chr1", 4, 4), ("chr2", 2, 5), ("chr10", 3, 6), ("chr3", 5, 7)],
        ),
        (
            {"chr1", "chr2", "chr10", "chr3"},
            [("chr1", 1, 3), ("chr1", 4, 4), ("chr2", 2, 5), ("chr10", 3, 6), ("chr3", 5, 7)],
        ),
        ({"chrX"}, []),
        ({"chrX", "chr3"}, [("chr3", 5, 7)]),
        (set(), []),
    ],
    ids=[
        "chr1_not_chr10",
        "chr10_not_chr1",
        "chr2_and_chr3_as_a_list",
        "none_reads_every_record",
        "every_contig",
        "a_name_the_reference_lacks",
        "a_name_the_reference_lacks_and_chr3",
        "no_contig",
    ],
)
def test_only_records_of_contigs_given_are_read(
    tmp_path: Path,
    four_contigs: IndexedFasta,
    contigs: set[str] | list[str] | None,
    places: list[tuple[str, int, int]],
) -> None:
    path = write_vcf(tmp_path, CONTIG_LINES)

    records = parse_vcf(path, four_contigs, contigs=contigs)

    assert get_places(records) == places


@pytest.mark.parametrize("num_threads", [1, 4])
@pytest.mark.parametrize(
    "chr2_line",
    [
        "chr2\t3\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
        "chr2\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1\n",
        "chr2\t3\t.\tG\tT\t.\t.\t.\n",
        "chr2\t0\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
        "chr2\t3\t.\tG\tT\t.\t.\t.\tGT\t0/x\t0/1\n",
        "chr2\t3\t.\tG\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n",
        "chr2\t3\t.\tG\tT\thigh\t.\t.\tGT:GQ\t0/1:high\t0/1:20\n",
        "chr2\t3\t.\tGGGGGGGGGG\tG\t.\t.\t.\tGT\t0/1\t0/1\n",
        "chr2\t3\t.\t\t\t.\t.\t.\tGT\t0/1\t0/1\n",
    ],
    ids=[
        "ref_differs",
        "too_few_sample_columns",
        "sites_only",
        "pos_0",
        "malformed_gt",
        "symbolic_alt",
        "qual_and_gq_not_numbers",
        "ref_past_the_contig_end",
        "ref_and_alt_empty",
    ],
)
def test_line_on_a_contig_not_read_gives_no_error_and_no_warning(
    tmp_path: Path,
    four_contigs: IndexedFasta,
    caplog: pytest.LogCaptureFixture,
    chr2_line: str,
    num_threads: int,
) -> None:
    lines = CONTIG_LINES.split("\n")
    path = write_vcf(tmp_path, "\n".join([*lines[:2], chr2_line.rstrip("\n"), *lines[3:]]))

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, four_contigs, contigs={"chr1"}, num_threads=num_threads))

    assert get_places(records) == [("chr1", 1, 3), ("chr1", 4, 4)]
    assert caplog.text == ""


def test_symbolic_record_on_a_contig_read_is_still_counted(
    tmp_path: Path, four_contigs: IndexedFasta, caplog: pytest.LogCaptureFixture
) -> None:
    # One symbolic record on chr1, read, and one on chr2, not read: the warning counts one.
    path = write_vcf(
        tmp_path,
        "chr1\t2\t.\tC\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n"
        "chr2\t3\t.\tG\t<DEL>\t.\t.\t.\tGT\t0/1\t0/1\n",
    )

    with caplog.at_level(logging.WARNING):
        records = list(parse_vcf(path, four_contigs, contigs={"chr1"}))

    assert records == []
    assert "skipped 1 VCF record with symbolic alleles" in caplog.text


# The message parse_vcf gives for a line out of order.
SORT_ADVICE = " Sort it in the contig order of the reference's .fai and by position, and run again."


@pytest.mark.parametrize(
    ("body", "contigs", "message"),
    [
        (
            "chr2\t3\t.\tG\tT\t.\t.\n",
            {"chr1"},
            "line 3: the line has 7 tab-separated columns, and a VCF record has at least 8,"
            " CHROM to INFO. Check that its columns are separated by tabs, correct the line or"
            " remove it, and run again.",
        ),
        (
            "chr2\n",
            {"chr1"},
            "line 3: the line has 1 tab-separated column, and a VCF record has at least 8,"
            " CHROM to INFO. Check that its columns are separated by tabs, correct the line or"
            " remove it, and run again.",
        ),
        (
            # Two faults: too few columns comes first.
            "chrUn\tfive\t.\tG\tT\t.\t.\n",
            {"chr1"},
            "line 3: the line has 7 tab-separated columns, and a VCF record has at least 8,"
            " CHROM to INFO. Check that its columns are separated by tabs, correct the line or"
            " remove it, and run again.",
        ),
        (
            "chr2\tfive\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
            {"chr1"},
            "line 3: POS is 'five', and it must be a whole number of 1 or more. Correct the"
            " line or remove it, and run again.",
        ),
        (
            # Two faults: a POS that is not a number comes before the contig.
            "chrUn\t-3\t.\tG\tT\t.\t.\t.\n",
            {"chr1"},
            "line 3: POS is '-3', and it must be a whole number of 1 or more. Correct the"
            " line or remove it, and run again.",
        ),
        (
            f"chr2\t{TOO_MANY_DIGITS}\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
            {"chr1"},
            "line 3: POS has 5,000 digits, which is past the end of any contig. Correct the"
            " line or remove it, and run again.",
        ),
        (
            "chrUn\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1\n",
            {"chr1"},
            "line 3: contig chrUn is not in the reference. Compare against the reference the"
            " VCF was called on.",
        ),
        (
            # Two faults: a contig the reference lacks comes before the order.
            "chr2\t5\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\nchrUn\t3\t.\tG\tT\t.\t.\t.\n",
            {"chr1"},
            "line 4: contig chrUn is not in the reference. Compare against the reference the"
            " VCF was called on.",
        ),
    ],
    ids=[
        "seven_columns",
        "one_column",
        "seven_columns_and_contig_lacking",
        "pos_not_a_number",
        "pos_not_a_number_and_contig_lacking",
        "pos_of_too_many_digits",
        "contig_lacking",
        "contig_lacking_and_out_of_order",
    ],
)
def test_line_on_a_contig_not_read_gives_the_error_of_its_first_fault(
    tmp_path: Path, four_contigs: IndexedFasta, body: str, contigs: set[str], message: str
) -> None:
    path = write_vcf(tmp_path, body)

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, four_contigs, contigs=contigs))

    assert str(raised.value) == f"{path}, {message}"


@pytest.mark.parametrize(
    ("body", "contigs", "out_of_order"),
    [
        (
            # Both lines on chr2, not read; POS 0 is taken for its order.
            "chr2\t5\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr2\t0\t.\tG\tT\t.\t.\t.\n",
            {"chr1"},
            "chr2:0 (line 4) comes after chr2:5 (line 3).",
        ),
        (
            # A line of chr2, not read, before chr1 by the .fai's order.
            "chr2\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t2\t.\tC\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
            {"chr1"},
            "chr1:2 (line 4) comes after chr2:3 (line 3).",
        ),
        (
            # A line of chr2, not read, after one of chr10, not read, which the .fai puts
            # after chr2.
            "chr1\t2\t.\tC\tA\t.\t.\t.\tGT\t0/1\t0/1\n"
            "chr10\t4\t.\tT\tC\t.\t.\t.\tGT\t0/1\t0/1\n"
            "chr2\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
            {"chr1"},
            "chr2:3 (line 5) comes after chr10:4 (line 4).",
        ),
        (
            # A line of chr2, read, after one of chr10, not read.
            "chr10\t4\t.\tT\tC\t.\t.\t.\tGT\t0/1\t0/1\nchr2\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0/1\n",
            {"chr2"},
            "chr2:3 (line 4) comes after chr10:4 (line 3).",
        ),
        (
            # A line of chr1, not read, after one of chr1, read, further on.
            "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\nchr1\t2\t.\tC\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
            {"chr2"},
            "chr1:2 (line 4) comes after chr1:5 (line 3).",
        ),
        (
            # The line before is not the line of the record before: a blank line and a
            # comment line come between.
            "chr2\t5\t.\tG\tT\t.\t.\t.\n\n# a comment\nchr2\t4\t.\tG\tT\t.\t.\t.\n",
            {"chr1"},
            "chr2:4 (line 6) comes after chr2:5 (line 3).",
        ),
        (
            "chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n\n# a comment\n"
            "chr1\t2\t.\tC\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
            None,
            "chr1:2 (line 6) comes after chr1:5 (line 3).",
        ),
    ],
    ids=[
        "not_read_after_not_read",
        "read_after_not_read_of_a_later_contig",
        "not_read_after_not_read_of_a_later_contig",
        "read_after_not_read",
        "same_contig_neither_read",
        "blank_and_comment_lines_between_neither_read",
        "blank_and_comment_lines_between_both_read",
    ],
)
def test_line_out_of_order_with_a_contig_not_read_stops_the_run(
    tmp_path: Path,
    four_contigs: IndexedFasta,
    body: str,
    contigs: set[str] | None,
    out_of_order: str,
) -> None:
    path = write_vcf(tmp_path, body)

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, four_contigs, contigs=contigs))

    assert str(raised.value) == f"{path} is not sorted: {out_of_order}{SORT_ADVICE}"


def test_ref_that_differs_on_the_first_record_read_says_so(
    tmp_path: Path, four_contigs: IndexedFasta
) -> None:
    # chr1's two lines agree with the reference and are not read; chr2's REF differs.
    lines = CONTIG_LINES.split("\n")
    bad_chr2_line = "chr2\t3\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1"
    path = write_vcf(tmp_path, "\n".join([*lines[:2], bad_chr2_line, *lines[3:]]))

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, four_contigs, contigs={"chr2", "chr3"}))

    assert str(raised.value) == (
        f"{path}, line 5: the REF of chr2:3 is 'A', and the reference, read in upper case, has"
        " 'G' there. It is the first record read with no symbolic allele, so the VCF looks"
        " made on another reference: compare against the reference it was called on, and"
        " run again."
    )


def test_ref_that_differs_counts_only_the_records_read_before_it(
    tmp_path: Path, four_contigs: IndexedFasta
) -> None:
    # chr1's two lines are read and agree, chr2's and chr10's agree and are not read, and
    # chr3's REF differs: 2 records, not 4.
    lines = CONTIG_LINES.split("\n")
    bad_chr3_line = "chr3\t6\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\n"
    path = write_vcf(tmp_path, "\n".join([*lines[:4], bad_chr3_line]))

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_vcf(path, four_contigs, contigs={"chr1", "chr3"}))

    assert str(raised.value) == (
        f"{path}, line 7: the REF of chr3:6 is 'A', and the reference, read in upper case, has"
        " 'C' there. The 2 records before it with no symbolic allele agree with the"
        " reference, so this line looks damaged: correct it or remove it, and run again."
    )


def test_contigs_given_as_one_string_raises_type_error(
    tmp_path: Path, four_contigs: IndexedFasta
) -> None:
    # A str is a collection of its letters, which would read nothing in silence.
    path = write_vcf(tmp_path, CONTIG_LINES)

    with pytest.raises(
        TypeError,
        match=r"^contigs must be a collection of contig names or None, not the string 'chr1'$",
    ):
        list(parse_vcf(path, four_contigs, contigs="chr1"))


# A parse_vcf stream with a pool, dropped unclosed in a reference cycle and finalized by the
# garbage collector in a thread that holds threading's _active_limbo_lock, as
# threading.enumerate() and Thread.start() hold it. threading has no public way to hold that
# lock, and holding it is what makes a finalizer that joins the pool's workers wait forever:
# a worker that ends takes the lock to leave threading's list of threads. It runs in a
# process of its own, so that a hang fails the test at its timeout and leaves the lock of
# this process free.
_FINALIZED_STREAM_SCRIPT = """
import gc, sys, threading
from pathlib import Path
from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import parse_vcf

directory = Path(sys.argv[1])
gc.disable()
with IndexedFasta(directory / "ref.fa") as reference:
    records = parse_vcf(directory / "a.vcf", reference, num_threads=2, span_limit=1)
    next(records)
    cycle = [records]
    cycle.append(cycle)
    del records, cycle
    with threading._active_limbo_lock:
        gc.collect()
print("collected")
"""
# Seconds the process is given; it takes well under one when nothing waits.
_FINALIZED_STREAM_TIMEOUT = 30


def test_stream_finalized_while_threads_are_listed_does_not_hang(tmp_path: Path) -> None:
    # 5,000 records one base apart, so that with a span limit of 1 the first record comes
    # out while the pool still has batches to parse.
    fasta, fai = make_fasta({"chr1": "A" * 5000})
    (tmp_path / "ref.fa").write_text(fasta, encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text(fai, encoding="utf-8")
    lines = (f"chr1\t{pos}\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/1\n" for pos in range(1, 5000))
    (tmp_path / "a.vcf").write_text(HEADER + "".join(lines), encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-c", _FINALIZED_STREAM_SCRIPT, str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
        timeout=_FINALIZED_STREAM_TIMEOUT,
    )

    assert (completed.returncode, completed.stdout, completed.stderr) == (0, "collected\n", "")


# Seconds a worker of the test below parses its item for; close() must not wait for it.
_SLOW_ITEM_SECONDS = 5


def test_stream_closed_early_does_not_wait_for_a_worker_still_parsing() -> None:
    release = threading.Event()

    def parse_item(item: int) -> int:
        # Item 1 is still being parsed when the stream is closed, until released or for
        # _SLOW_ITEM_SECONDS at most, so that a close() that waits for it ends.
        if item == 1:
            release.wait(_SLOW_ITEM_SECONDS)
        return item

    results = vcf._map_in_order(parse_item, range(10), 2)
    assert next(results) == 0
    started = time.monotonic()
    results.close()
    close_seconds = time.monotonic() - started
    release.set()

    assert close_seconds < 1
