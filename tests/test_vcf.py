import gzip
import itertools
import logging
import random
import threading
from collections.abc import Callable, Iterator
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
        f"{path}, line 3: the REF of chr1:5 is Ca, and the reference has AG there, so the"
        " VCF looks made on another reference, or damaged. Compare against the reference"
        " the VCF was called on."
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
    ("reference_base", "ref"),
    [("A", "C"), ("R", "C"), ("N", "X"), ("X", "X")],
    ids=["a_against_c", "reference_r_against_c", "no_code_against_n", "no_code_against_itself"],
)
def test_iupac_codes_that_share_no_base_differ(
    tmp_path: Path, make_reference: MakeReference, reference_base: str, ref: str
) -> None:
    # A letter that is no IUPAC code stands for no base, so it agrees with nothing.
    path = write_vcf(tmp_path, f"chr1\t3\t.\t{ref}\t.\t.\t.\t.\tGT\t0/0\t0/0\n")
    reference = make_reference(*make_fasta({"chr1": f"GG{reference_base}GG"}))

    with pytest.raises(PvcEvalInputError, match=r"line 3: the REF of chr1:3 "):
        list(parse_vcf(path, reference))


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
    [0, 100, vcf._BATCH_SIZE - 1, vcf._BATCH_SIZE, 5 * vcf._BATCH_SIZE],
    ids=[
        "first_record",
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
    assert str(raised.value) == (
        f"{path}, line {line_number}: the REF of chr1:{pos_1based} is T, and the reference"
        " has A there, so the VCF looks made on another reference, or damaged. Compare"
        " against the reference the VCF was called on."
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
