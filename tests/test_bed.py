"""The lines of a --regions or --strata BED, parse_bed, and the refusals that read_regions
and read_strata share (doc/specs/loci.md "Regions" and "Stratum").

A BED line is 0-based and half-open, as the intervals parse_bed yields are.
"""

import gzip
from pathlib import Path

import pytest

from pvc_eval.bed import BedLine, parse_bed
from pvc_eval.errors import PvcEvalInputError
from pvc_eval.locus import Locus
from pvc_eval.regions import read_regions
from pvc_eval.strata import read_strata

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
            "the start, '1e3', is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t\u00b2\t2000",
            "the start, '\u00b2', is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t100\t-20",
            "the end, '-20', is not a whole number. A BED line gives its start and its end as"
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
            "the contig 'chrUn' is not in the reference. Keep only the reference's contigs in the"
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


# The BED of _COMPRESSED_LINES, as bgzip writes it: three gzip members, the last one empty.
# Made with `printf '# GIAB\nchr1\t100\t200\tsegdup\nchr2\t0\t50\thomopolymer\n' | bgzip -c`,
# bgzip of htslib 1.24.
_BGZIPPED_BED = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00%\x00\x01\x07\x00\xf8\xff"
    b"# GIAB\n\xc9I\x98;\x07\x00\x00\x00\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC"
    b'\x02\x00C\x00K\xce(2\xe4440\xe04\x02\xe2\xe2\xd4\xf4\x94\xd2\x02\xae\xe4\x8c"#N\x03NS'
    b'\x03\xce\x8c\xfc\xdc\xfc\x82\xfc\x9c\xca\xdc\xd4".\x00\x0c\xd6]\xb2*\x00\x00\x00\x1f\x8b'
    b"\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00"
    b"\x00\x00\x00"
)
_COMPRESSED_LINES = b"# GIAB\nchr1\t100\t200\tsegdup\nchr2\t0\t50\thomopolymer\n"
_COMPRESSED_BEDS = [gzip.compress(_COMPRESSED_LINES), _BGZIPPED_BED]
_COMPRESSED_IDS = ["gzipped", "bgzipped"]


def write_compressed_bed(tmp_path: Path, compressed: bytes) -> Path:
    path = tmp_path / "regions.bed.gz"
    path.write_bytes(compressed)
    return path


@pytest.mark.parametrize("compressed", _COMPRESSED_BEDS, ids=_COMPRESSED_IDS)
def test_parse_bed_reads_a_gzipped_bed(tmp_path: Path, compressed: bytes) -> None:
    path = write_compressed_bed(tmp_path, compressed)

    assert list(parse_bed(path, _CONTIGS)) == [
        BedLine(contig="chr1", start=100, end=200, name="segdup", line_number=2),
        BedLine(contig="chr2", start=0, end=50, name="homopolymer", line_number=3),
    ]


@pytest.mark.parametrize("compressed", _COMPRESSED_BEDS, ids=_COMPRESSED_IDS)
def test_read_regions_reads_a_gzipped_bed(tmp_path: Path, compressed: bytes) -> None:
    regions = read_regions(write_compressed_bed(tmp_path, compressed), _CONTIGS)

    assert [
        regions.overlaps("chr1", 99, 100),
        regions.overlaps("chr1", 199, 200),
        regions.overlaps("chr1", 200, 201),
        regions.overlaps("chr2", 49, 50),
    ] == [False, True, False, True]


@pytest.mark.parametrize("compressed", _COMPRESSED_BEDS, ids=_COMPRESSED_IDS)
def test_read_strata_reads_a_gzipped_bed(tmp_path: Path, compressed: bytes) -> None:
    strata = read_strata(write_compressed_bed(tmp_path, compressed), _CONTIGS)

    assert [
        strata.stratum(Locus("chr1", 150, 160, (), (), ())),
        strata.stratum(Locus("chr2", 10, 20, (), (), ())),
        strata.stratum(Locus("chr2", 50, 60, (), (), ())),
    ] == ["segdup", "homopolymer", "."]


@pytest.mark.parametrize(
    "compressed",
    [b"", gzip.compress(b""), gzip.compress(b"#chrom\tstart\tend\n\n")],
    ids=["plain", "gzipped", "gzipped-header-only"],
)
def test_parse_bed_refuses_an_empty_bed(tmp_path: Path, compressed: bytes) -> None:
    path = write_compressed_bed(tmp_path, compressed)

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == _make_empty_message(path)


def test_parse_bed_refuses_a_bed_of_only_blank_comment_and_header_lines(
    tmp_path: Path,
) -> None:
    path = write_bed(tmp_path, ["# a comment", "track name=regions", "browser hide all", ""])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == _make_empty_message(path)


def test_read_regions_and_read_strata_refuse_an_empty_bed(tmp_path: Path) -> None:
    path = write_bed(tmp_path, [])

    with pytest.raises(PvcEvalInputError) as raised_by_regions:
        read_regions(path, _CONTIGS)
    with pytest.raises(PvcEvalInputError) as raised_by_strata:
        read_strata(path, _CONTIGS)

    assert [str(raised_by_regions.value), str(raised_by_strata.value)] == [
        _make_empty_message(path),
        _make_empty_message(path),
    ]


def _make_empty_message(path: Path) -> str:
    return (
        f"{path} holds no BED line: it is empty, or holds only blank lines, comments, and"
        " track and browser lines. Give a BED with at least one line of contig, start and"
        " end, separated by tabs, and run again."
    )


def test_parse_bed_refuses_a_gzipped_bed_that_is_not_utf8(tmp_path: Path) -> None:
    # The gunzip hint the message gave when a gzipped BED was refused went with the
    # refusal: a BED that is still not UTF-8 is compressed by another program, or in
    # another encoding.
    path = write_compressed_bed(tmp_path, gzip.compress(b"chr1\t0\t10\tsegdup\xff\n"))

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == (
        f"{path} is not UTF-8 text: a byte at or after line 1 cannot be read as UTF-8. If it"
        " is compressed by a program other than gzip or bgzip, decompress it; if it is text"
        " in another encoding, convert it to UTF-8; and run again."
    )


# 200 lines, so that a fault found at the end of the member is after line 200.
_DAMAGED_GZIP = gzip.compress(b"chr1\t0\t10\n" * 200)


@pytest.mark.parametrize(
    ("compressed", "line_number"),
    [
        # Its last 8 bytes, the CRC-32 and the length that end a gzip member, cut.
        (_DAMAGED_GZIP[:-8], 201),
        # The first bytes of the deflate stream, after the header of 10, flipped, so that
        # zlib cannot decompress it.
        (
            _DAMAGED_GZIP[:12]
            + bytes(byte ^ 0xFF for byte in _DAMAGED_GZIP[12:20])
            + _DAMAGED_GZIP[20:],
            1,
        ),
        # Its CRC-32 zeroed.
        (_DAMAGED_GZIP[:-8] + bytes(4) + _DAMAGED_GZIP[-4:], 201),
    ],
    ids=["cut-short", "deflate-damaged", "crc-wrong"],
)
def test_parse_bed_refuses_a_damaged_gzip(
    tmp_path: Path, compressed: bytes, line_number: int
) -> None:
    path = write_compressed_bed(tmp_path, compressed)

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value) == (
        f"{path} is a gzip file cut short or damaged at or after line {line_number}. Compress"
        " or copy it again, and run again."
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


@pytest.mark.parametrize(
    ("contig", "contigs", "named", "command"),
    [
        ("1", ("chr1", "chr2"), "chr1", '$1 !~ /^(#|track|browser|chr)/ {{$1 = "chr" $1}} 1'),
        ("chr1", ("1", "2"), "1", '$1 !~ /^(#|track|browser)/ {{sub(/^chr/, "", $1)}} 1'),
    ],
    ids=["chr-missing", "chr-extra"],
)
def test_parse_bed_says_to_rename_a_contig_named_in_another_style(
    tmp_path: Path, contig: str, contigs: tuple[str, ...], named: str, command: str
) -> None:
    path = write_bed(tmp_path, [f"{contig}\t100\t200"])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, contigs))

    assert str(raised.value) == (
        f"{path}, line 1: the contig '{contig}' is not in the reference, which names it"
        f" '{named}'. Name the BED's contigs as the reference does, for instance with"
        f" `awk -F'\\t' -v OFS='\\t' '{command.format()}' {path} > renamed.bed`, and run again"
        " with renamed.bed."
    )


def test_parse_bed_quotes_a_contig_so_that_a_byte_order_mark_shows(tmp_path: Path) -> None:
    # A byte-order mark at the start of a line other than the first is not dropped.
    path = write_bed(tmp_path, ["chr1\t0\t10", "\ufeffchr2\t0\t10"])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value).startswith(
        f"{path}, line 2: the contig '\\ufeffchr2' is not in the reference."
    )


_BYTE_ORDER_MARK_BED = "\ufeffchr1\t190\t300\n".encode()


@pytest.mark.parametrize(
    "contents",
    [_BYTE_ORDER_MARK_BED, gzip.compress(_BYTE_ORDER_MARK_BED)],
    ids=["plain", "gzipped"],
)
def test_parse_bed_drops_a_byte_order_mark_at_the_start(tmp_path: Path, contents: bytes) -> None:
    path = tmp_path / "regions.bed"
    path.write_bytes(contents)

    assert list(parse_bed(path, _CONTIGS)) == [
        BedLine(contig="chr1", start=190, end=300, name=None, line_number=1)
    ]


def test_parse_bed_drops_only_the_first_of_two_byte_order_marks(tmp_path: Path) -> None:
    path = write_bed(tmp_path, ["\ufeff\ufeffchr1\t0\t10"])

    with pytest.raises(PvcEvalInputError) as raised:
        list(parse_bed(path, _CONTIGS))

    assert str(raised.value).startswith(
        f"{path}, line 1: the contig '\\ufeffchr1' is not in the reference."
    )


def test_parse_bed_reads_a_contig_whose_name_starts_with_track_or_browser(
    tmp_path: Path,
) -> None:
    lines = ["track name=x", "browser hide all", "tracks_1\t0\t10", "browser2\t5\t15"]

    assert list(parse_bed(write_bed(tmp_path, lines), ("tracks_1", "browser2"))) == [
        BedLine(contig="tracks_1", start=0, end=10, name=None, line_number=3),
        BedLine(contig="browser2", start=5, end=15, name=None, line_number=4),
    ]
