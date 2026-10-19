import gzip
import os
import random
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta

# Lines of five bases, a last line shorter than the others, lower-case bases, a contig
# whose last line has no line end. Its .fai is what `samtools faidx one.fa` (samtools
# 1.24) wrote for it.
ONE_FA = b">chr1 first contig\nACGTa\ncgtNN\nAC\n>chr2\nGGGG\nTT\n>chr3\nacgtACG"
ONE_FAI = "chr1\t12\t19\t5\t6\nchr2\t6\t40\t4\t5\nchr3\t7\t54\t7\t8\n"

# Windows line ends, so bytes per line is bases per line + 2, and lines of three bases.
# Its .fai is what `samtools faidx crlf.fa` (samtools 1.24) wrote for it.
CRLF_FA = b">c1\r\nACG\r\nTTA\r\nG\r\n>c2\r\nAAAA\r\n"
CRLF_FAI = "c1\t7\t5\t3\t5\nc2\t4\t23\t4\t6\n"


def write_fasta(directory: Path, name: str, fasta: bytes, fai: str | None) -> Path:
    path = directory / name
    path.write_bytes(fasta)
    if fai is not None:
        Path(f"{path}.fai").write_text(fai, encoding="utf-8")
    return path


def fasta_and_fai(contigs: dict[str, str], line_width: int, line_end: str) -> tuple[bytes, str]:
    """A FASTA of the contigs, wrapped at line_width, and the .fai samtools would write."""
    fasta = ""
    fai_lines = []
    for name, sequence in contigs.items():
        fasta += f">{name}\n"
        lines = [sequence[i : i + line_width] for i in range(0, len(sequence), line_width)]
        fai_lines.append(
            f"{name}\t{len(sequence)}\t{len(fasta)}\t{line_width}\t{line_width + len(line_end)}\n"
        )
        fasta += "".join(line + line_end for line in lines)
    return fasta.encode("ascii"), "".join(fai_lines)


@pytest.fixture
def one_fa(tmp_path: Path) -> Path:
    return write_fasta(tmp_path, "one.fa", ONE_FA, ONE_FAI)


@pytest.fixture
def crlf_fa(tmp_path: Path) -> Path:
    return write_fasta(tmp_path, "crlf.fa", CRLF_FA, CRLF_FAI)


def test_contigs_and_lengths_come_from_the_fai(one_fa: Path) -> None:
    with IndexedFasta(one_fa) as indexed_fasta:
        assert indexed_fasta.contigs == ("chr1", "chr2", "chr3")
        assert [indexed_fasta.length(contig) for contig in indexed_fasta.contigs] == [12, 6, 7]


def test_contig_order_is_the_fai_order_not_the_fasta_order(tmp_path: Path) -> None:
    fai = "chr3\t7\t54\t7\t8\nchr1\t12\t19\t5\t6\nchr2\t6\t40\t4\t5\n"
    path = write_fasta(tmp_path, "one.fa", ONE_FA, fai)
    with IndexedFasta(path) as indexed_fasta:
        assert indexed_fasta.contigs == ("chr3", "chr1", "chr2")
        assert indexed_fasta.bases("chr3", 0, 7) == "ACGTACG"


@pytest.mark.parametrize(
    ("contig", "start", "end", "expected"),
    [
        ("chr1", 0, 3, "ACG"),
        ("chr1", 0, 1, "A"),
        ("chr1", 9, 12, "NAC"),
        ("chr1", 11, 12, "C"),
        ("chr1", 4, 5, "A"),
        ("chr1", 4, 6, "AC"),
        ("chr1", 3, 11, "TACGTNNA"),
        ("chr1", 5, 5, ""),
        ("chr1", 0, 0, ""),
        ("chr1", 12, 12, ""),
        ("chr1", 0, 12, "ACGTACGTNNAC"),
        ("chr2", 0, 6, "GGGGTT"),
        ("chr2", 3, 5, "GT"),
        ("chr3", 0, 7, "ACGTACG"),
        ("chr3", 6, 7, "G"),
    ],
    ids=[
        "first_bases",
        "first_base",
        "last_bases",
        "last_base",
        "last_base_of_a_line",
        "across_one_line_end",
        "across_two_line_ends",
        "empty_inside",
        "empty_at_start",
        "empty_at_end",
        "whole_contig",
        "whole_second_contig",
        "second_contig_across_line_end",
        "whole_contig_without_final_line_end",
        "last_base_without_final_line_end",
    ],
)
def test_bases_of_a_region_in_upper_case(
    one_fa: Path, contig: str, start: int, end: int, expected: str
) -> None:
    with IndexedFasta(one_fa) as indexed_fasta:
        assert indexed_fasta.bases(contig, start, end) == expected


@pytest.mark.parametrize(
    ("contig", "start", "end", "expected"),
    [
        ("c1", 0, 3, "ACG"),
        ("c1", 2, 5, "GTT"),
        ("c1", 2, 7, "GTTAG"),
        ("c1", 6, 7, "G"),
        ("c1", 0, 7, "ACGTTAG"),
        ("c2", 0, 4, "AAAA"),
        ("c2", 4, 4, ""),
    ],
    ids=[
        "first_line",
        "across_one_line_end",
        "across_two_line_ends_to_the_end",
        "last_base",
        "whole_contig",
        "whole_second_contig",
        "empty_at_end",
    ],
)
def test_bases_with_windows_line_ends(
    crlf_fa: Path, contig: str, start: int, end: int, expected: str
) -> None:
    with IndexedFasta(crlf_fa) as indexed_fasta:
        assert indexed_fasta.bases(contig, start, end) == expected


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools is not installed")
def test_bases_agree_with_a_fai_made_by_samtools(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "one.fa", ONE_FA, None)
    subprocess.run(["samtools", "faidx", str(path)], check=True)
    with IndexedFasta(path) as indexed_fasta:
        assert indexed_fasta.contigs == ("chr1", "chr2", "chr3")
        assert indexed_fasta.bases("chr1", 3, 11) == "TACGTNNA"


@given(
    sequences=st.lists(
        st.text(alphabet="ACGTacgtN", min_size=1, max_size=30), min_size=1, max_size=4
    ),
    line_width=st.integers(min_value=1, max_value=12),
    line_end=st.sampled_from(["\n", "\r\n"]),
)
def test_every_region_is_the_upper_case_slice_of_its_contig(
    tmp_path_factory: pytest.TempPathFactory, sequences: list[str], line_width: int, line_end: str
) -> None:
    # tmp_path_factory and not tmp_path, which hypothesis refuses: tmp_path would be one
    # directory for all the examples.
    contigs = {f"contig{i}": sequence for i, sequence in enumerate(sequences)}
    fasta, fai = fasta_and_fai(contigs, line_width, line_end)
    path = write_fasta(tmp_path_factory.mktemp("example"), "ref.fa", fasta, fai)
    with IndexedFasta(path) as indexed_fasta:
        assert indexed_fasta.contigs == tuple(contigs)
        for name, sequence in contigs.items():
            assert indexed_fasta.length(name) == len(sequence)
            for start in range(len(sequence) + 1):
                for end in range(start, len(sequence) + 1):
                    assert indexed_fasta.bases(name, start, end) == sequence[start:end].upper()


def test_threads_reading_one_indexed_fasta_get_the_bases_one_thread_gets(tmp_path: Path) -> None:
    rng = random.Random(20261019)
    contigs = {
        f"chr{i}": "".join(rng.choices("ACGTacgtN", k=rng.randint(2000, 5000))) for i in range(4)
    }
    fasta, fai = fasta_and_fai(contigs, 60, "\n")
    path = write_fasta(tmp_path, "ref.fa", fasta, fai)
    num_threads, regions_per_thread = 8, 200
    regions = []
    for _ in range(num_threads * regions_per_thread):
        name = rng.choice(list(contigs))
        start = rng.randint(0, len(contigs[name]))
        end = rng.randint(start, min(start + 300, len(contigs[name])))
        regions.append((name, start, end))
    with IndexedFasta(path) as indexed_fasta:
        one_thread = [indexed_fasta.bases(*region) for region in regions]
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            several_threads = list(
                executor.map(lambda region: indexed_fasta.bases(*region), regions)
            )
    assert several_threads == one_thread
    assert one_thread == [contigs[name][start:end].upper() for name, start, end in regions]


def test_a_missing_fai_says_how_to_make_one(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "one.fa", ONE_FA, None)
    expected = (
        f"{path} has no index {path}.fai. Make it with `samtools faidx {path}` and run again."
    )
    with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
        IndexedFasta(path)


def test_a_missing_fasta_says_to_check_the_path(tmp_path: Path) -> None:
    path = tmp_path / "ref.fa"
    expected = f"{path} does not exist. Check the path given to --reference, and run again."
    with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
        IndexedFasta(path)


def test_an_empty_fai_says_to_make_it_again(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "one.fa", ONE_FA, "")
    expected = (
        f"{path}.fai names no contig."
        f" Make the .fai again with `samtools faidx {path}` and run again."
    )
    with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
        IndexedFasta(path)


def test_a_fai_that_is_not_utf8_is_named_with_its_line(tmp_path: Path) -> None:
    # The .fai is what `samtools faidx` (samtools 1.24) wrote for this FASTA, whose second
    # contig name is in Latin-1: it copies the name as it is.
    path = tmp_path / "lat.fa"
    path.write_bytes(b">chr1\nACGT\n>chr\xe92\nACGT\n")
    Path(f"{path}.fai").write_bytes(b"chr1\t4\t6\t4\t5\nchr\xe92\t4\t18\t4\t5\n")
    expected = (
        f"{path}.fai, line 2, is not UTF-8 text."
        f" Make the .fai again with `samtools faidx {path}` and run again."
    )
    with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
        IndexedFasta(path)


def test_a_compressed_fasta_says_to_decompress_it(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "one.fa.gz", gzip.compress(ONE_FA), ONE_FAI)
    expected = (
        f"{path} is compressed with gzip or bgzip, and pvc_eval reads only an uncompressed"
        f" FASTA. Decompress it with `gunzip -c {path} > file.fa`, index it with"
        " `samtools faidx file.fa`, and give file.fa instead."
    )
    with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
        IndexedFasta(path)


@pytest.mark.parametrize(
    ("fai", "line_number", "fault"),
    [
        (
            "chr1\t12\t19\t5\t6\nchr2\t6\t40\t4\n",
            2,
            "has 4 columns where a .fai has 5: name, length, offset, bases per line,"
            " bytes per line",
        ),
        (
            "chr1\t12\t19\t5\t6\t7\n",
            1,
            "has 6 columns where a .fai has 5: name, length, offset, bases per line,"
            " bytes per line",
        ),
        (
            "chr1\ttwelve\t19\t5\t6\n",
            1,
            "has twelve as its length, which is not a whole number of 0 or more",
        ),
        (
            "chr1\t12\t-19\t5\t6\n",
            1,
            "has -19 as its offset, which is not a whole number of 0 or more",
        ),
        ("\t12\t19\t5\t6\n", 1, "has no contig name"),
        ("chr1\t12\t19\t0\t1\n", 1, "has 0 bases per line, where a line holds at least one base"),
        (
            "chr1\t12\t19\t5\t5\n",
            1,
            "has 5 bytes per line for 5 bases per line, which leaves no byte for the line end",
        ),
        (
            "chr1\t12\t19\t5\t6\nchr1\t6\t40\t4\t5\n",
            2,
            "names contig chr1 again, which line 1 named first",
        ),
        (
            "chr1\t12\t19\t5\t6\nchr2\t6\t40\t4\t5\nchr3\t8\t54\t7\t8\n",
            3,
            "puts the bases of chr3 past the end of the FASTA, which has 61 bytes",
        ),
        (
            "chr1\t12\t19\t5\t6\nchr2\t6\t70\t4\t5\n",
            2,
            "puts the bases of chr2 past the end of the FASTA, which has 61 bytes",
        ),
        (
            "chr1\t12\t19\t5\t6\nchr2\t6\t40\t4\t5\nchr3\t7\t55\t7\t8\n",
            3,
            "puts the bases of chr3 past the end of the FASTA, which has 61 bytes",
        ),
    ],
    ids=[
        "too_few_columns",
        "too_many_columns",
        "length_not_a_number",
        "negative_offset",
        "no_name",
        "no_bases_per_line",
        "no_byte_for_line_end",
        "contig_twice",
        "length_past_end_of_file",
        "offset_past_end_of_file",
        "last_base_one_byte_past_end_of_file",
    ],
)
def test_a_wrong_fai_line_is_named(tmp_path: Path, fai: str, line_number: int, fault: str) -> None:
    path = write_fasta(tmp_path, "one.fa", ONE_FA, fai)
    expected = (
        f"{path}.fai, line {line_number}: {fault}."
        f" Make the .fai again with `samtools faidx {path}` and run again."
    )
    with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
        IndexedFasta(path)


def not_bases_message(path: Path, region: str, fai_line_number: int) -> str:
    return (
        f"{path}: the bytes of {region} that its .fai, line {fai_line_number}, gives are not"
        f" all letters. Make the .fai again with `samtools faidx {path}`, check that the"
        " FASTA holds only letters between its header lines, and run again."
    )


# Each FASTA differs from the one its .fai was made for.
@pytest.mark.parametrize(
    ("fasta", "fai", "contig", "start", "end", "region", "fai_line_number"),
    [
        # The lines have four bases, not the five the .fai says, so the region read
        # crosses a line end the .fai does not know of.
        (b">chr1\nACGT\nACGT\nACGT\n", "chr1\t12\t6\t5\t6\n", "chr1", 0, 12, "chr1:1-12", 1),
        # The .fai is what `samtools faidx` (samtools 1.24) wrote before the header of chr1
        # was lengthened by four bytes, ">chr1 first contig" to ">chr1 the first contig".
        (
            b">chr1 the first contig\nACGTACGTAC\nGG\n>chr2\nTTTT\n",
            "chr1\t12\t19\t10\t11\nchr2\t4\t39\t4\t5\n",
            "chr2",
            0,
            3,
            "chr2:1-3",
            2,
        ),
        (
            b">chr1 the first contig\nACGTACGTAC\nGG\n>chr2\nTTTT\n",
            "chr1\t12\t19\t10\t11\nchr2\t4\t39\t4\t5\n",
            "chr1",
            0,
            12,
            "chr1:1-12",
            1,
        ),
        # Windows line ends, with a .fai that says one byte of line end.
        (b">c1\r\nACG\r\nTTA\r\nG\r\n", "c1\t4\t5\t4\t5\n", "c1", 0, 4, "c1:1-4", 1),
        (b">chr1\nAC T\n", "chr1\t4\t6\t4\t5\n", "chr1", 0, 4, "chr1:1-4", 1),
        (">chr1\nACÑT\n".encode(), "chr1\t5\t6\t5\t6\n", "chr1", 0, 5, "chr1:1-5", 1),
    ],
    ids=[
        "shorter_lines",
        "longer_header_digit_in_region",
        "longer_header_line_end_in_region",
        "windows_line_end_in_region",
        "space_in_sequence_line",
        "byte_that_is_not_ascii",
    ],
)
def test_bytes_that_are_not_bases_where_the_fai_says_are_named(
    tmp_path: Path,
    fasta: bytes,
    fai: str,
    contig: str,
    start: int,
    end: int,
    region: str,
    fai_line_number: int,
) -> None:
    path = write_fasta(tmp_path, "ref.fa", fasta, fai)
    expected = not_bases_message(path, region, fai_line_number)
    with (
        IndexedFasta(path) as indexed_fasta,
        pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"),
    ):
        indexed_fasta.bases(contig, start, end)


def test_a_fasta_cut_short_after_opening_is_named(tmp_path: Path) -> None:
    path = write_fasta(tmp_path, "one.fa", ONE_FA, ONE_FAI)
    expected = not_bases_message(path, "chr3:1-7", 3)
    with IndexedFasta(path) as indexed_fasta:
        os.truncate(path, len(ONE_FA) - 2)
        with pytest.raises(PvcEvalInputError, match=f"^{re.escape(expected)}$"):
            indexed_fasta.bases("chr3", 0, 7)


@pytest.mark.parametrize(
    ("contig", "start", "end", "message"),
    [
        ("chr9", 0, 1, "contig chr9 is not in {path}"),
        ("chr1", -1, 3, "start must be 0 or more, not -1"),
        ("chr1", 5, 4, "end 4 is before start 5"),
        ("chr1", 0, 13, "end 13 is past the length of chr1, 12"),
        ("chr1", 13, 13, "end 13 is past the length of chr1, 12"),
    ],
    ids=[
        "unknown_contig",
        "negative_start",
        "end_before_start",
        "end_past_length",
        "start_past_length",
    ],
)
def test_a_region_outside_its_contig_is_a_value_error(
    one_fa: Path, contig: str, start: int, end: int, message: str
) -> None:
    with (
        IndexedFasta(one_fa) as indexed_fasta,
        pytest.raises(ValueError, match=f"^{re.escape(message.format(path=one_fa))}$"),
    ):
        indexed_fasta.bases(contig, start, end)


def test_length_of_an_unknown_contig_is_a_value_error(one_fa: Path) -> None:
    with (
        IndexedFasta(one_fa) as indexed_fasta,
        pytest.raises(ValueError, match=f"^{re.escape(f'contig chr9 is not in {one_fa}')}$"),
    ):
        indexed_fasta.length("chr9")


def test_a_closed_indexed_fasta_reads_nothing(one_fa: Path) -> None:
    with IndexedFasta(one_fa) as indexed_fasta:
        pass
    with pytest.raises(OSError, match="Bad file descriptor"):
        indexed_fasta.bases("chr1", 0, 3)


def test_closing_twice_leaves_a_file_opened_between_alone(one_fa: Path, tmp_path: Path) -> None:
    indexed_fasta = IndexedFasta(one_fa)
    indexed_fasta.close()
    other_path = tmp_path / "other.txt"
    other_path.write_bytes(b"other file")
    # Likely to get the descriptor number the IndexedFasta had.
    other_fd = os.open(other_path, os.O_RDONLY)
    try:
        indexed_fasta.close()
        assert os.pread(other_fd, 10, 0) == b"other file"
    finally:
        os.close(other_fd)
