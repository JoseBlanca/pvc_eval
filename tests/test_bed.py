"""The lines of a --regions or --strata BED, parse_bed, and the refusals that read_regions
and read_strata share (doc/specs/loci.md "Regions" and "Stratum").

A BED line is 0-based and half-open, as the intervals parse_bed yields are.
"""

import gzip
from pathlib import Path

import pytest

from pvc_eval.bed import BedLine, parse_bed
from pvc_eval.errors import PvcEvalInputError

_CONTIGS = ("chr1", "chr2")


def write_bed(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "regions.bed"
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return path


def test_parse_bed_skips_blank_and_header_lines_and_reads_the_fourth_column(
    tmp_path: Path,
) -> None:
    lines = [
        "# a comment",
        "track name=regions",
        "browser position chr1:1-1000",
        "",
        "chr2\t0\t5",
        "chr1\t300\t400\tsegdup\textra",
        "chr1\t10\t20\t",
    ]

    assert list(parse_bed(write_bed(tmp_path, lines), _CONTIGS)) == [
        BedLine(contig="chr2", start=0, end=5, name=None, line_number=5),
        BedLine(contig="chr1", start=300, end=400, name="segdup", line_number=6),
        BedLine(contig="chr1", start=10, end=20, name=None, line_number=7),
    ]


def test_parse_bed_reads_lines_that_end_in_carriage_return(tmp_path: Path) -> None:
    path = tmp_path / "regions.bed"
    path.write_bytes(b"chr1\t0\t10\r\n")

    assert list(parse_bed(path, _CONTIGS)) == [
        BedLine(contig="chr1", start=0, end=10, name=None, line_number=1)
    ]


@pytest.mark.parametrize(
    ("line", "message"),
    [
        (
            "chr1\t100",
            "has 2 columns, and a BED line needs at least three, separated by tabs: the"
            " contig, the start and the end. Correct the line or remove it, and run again.",
        ),
        (
            "chr1\t1e3\t2000",
            "the start, 1e3, is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t\u00b2\t2000",
            "the start, \u00b2, is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t100\t-20",
            "the end, -20, is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t100\t100",
            "the end, 100, is not after the start, 100. A BED line covers the positions from"
            " its start up to its end, 0-based and half-open, so its end is greater than its"
            " start; correct the line and run again.",
        ),
        (
            "chrUn\t100\t200",
            "the contig chrUn is not in the reference. Keep only the reference's contigs in the"
            " BED, for instance with `awk -F'\\t' 'NR == FNR {{keep[$1]; next}} $1 in keep'"
            " ref.fa.fai {path} > kept.bed`, where ref.fa.fai is the .fai of the --reference"
            " FASTA, and run again with kept.bed.",
        ),
    ],
    ids=[
        "two-columns",
        "start-not-whole",
        "start-superscript",
        "end-not-whole",
        "end-not-after-start",
        "contig",
    ],
)
def test_parse_bed_refuses_a_line(tmp_path: Path, line: str, message: str) -> None:
    path = write_bed(tmp_path, ["chr1\t0\t10", line])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == f"{path}, line 2: {message.format(path=path)}"


def test_parse_bed_refuses_a_missing_file(tmp_path: Path) -> None:
    path = tmp_path / "absent.bed"

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == f"{path} does not exist. Check its path and run again."


def test_parse_bed_refuses_a_gzipped_bed(tmp_path: Path) -> None:
    path = tmp_path / "regions.bed.gz"
    path.write_bytes(gzip.compress(b"chr1\t0\t10\n"))

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == (
        f"{path} is not UTF-8 text: a byte at or after line 1 cannot be read as UTF-8. If it"
        f" is compressed, decompress it with `gunzip -c {path} > file.bed` and give file.bed"
        " instead; if it is text in another encoding, convert it to UTF-8; and run again."
    )


def test_parse_bed_refuses_a_line_separated_by_spaces(tmp_path: Path) -> None:
    path = write_bed(tmp_path, ["chr1\t0\t10", "chr1 100  200"])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == (
        f"{path}, line 2: has its columns separated by spaces, and a BED separates them by"
        f" tabs. Convert it with `awk -v OFS='\\t' '{{$1 = $1}} 1' {path} > tabs.bed` and run"
        " again with tabs.bed."
    )


def test_parse_bed_refuses_a_directory(tmp_path: Path) -> None:
    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(tmp_path, _CONTIGS))

    assert str(raised.value) == (
        f"{tmp_path} is a directory, not a BED file. Check its path and run again."
    )


def test_parse_bed_refuses_a_file_it_cannot_read(tmp_path: Path) -> None:
    path = write_bed(tmp_path, ["chr1\t0\t10"])
    path.chmod(0)
    try:
        with pytest.raises(PvcEvalInputError) as raised:
            list(parse_bed(path, _CONTIGS))
    finally:
        path.chmod(0o644)

    assert str(raised.value) == (
        f"{path} cannot be read: permission denied. Allow reading it and run again."
    )
