"""The selection of --regions: read_regions, Regions and select_loci (doc/specs/loci.md
"Regions", doc/spec.md §3 and §7.2).

The spec's cases give positions 1-based and inclusive; loci are built 0-based and
half-open, and a BED line is 0-based and half-open, so the locus at 100 to 120 is (99, 120)
and the region at 110 to 200 is the BED line `chr1 109 200`.
"""

import gzip
from pathlib import Path

import pytest

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.locus import Locus
from pvc_eval.regions import read_regions, select_loci

_CONTIGS = ("chr1", "chr2")


def write_bed(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "regions.bed"
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return path


def make_locus(contig: str, start: int, end: int) -> Locus:
    return Locus(contig, start, end, (), (), ())


def test_region_sharing_a_position_keeps_the_locus(tmp_path: Path) -> None:
    regions = read_regions(write_bed(tmp_path, ["chr1\t109\t200"]), _CONTIGS)
    locus = make_locus("chr1", 99, 120)

    assert list(select_loci(iter([locus]), regions)) == [locus]


def test_region_next_to_the_locus_skips_it(tmp_path: Path) -> None:
    regions = read_regions(write_bed(tmp_path, ["chr1\t120\t200"]), _CONTIGS)

    assert list(select_loci(iter([make_locus("chr1", 99, 120)]), regions)) == []


def test_no_regions_keeps_every_locus() -> None:
    loci = [make_locus("chr1", 99, 120), make_locus("chr2", 5, 6)]

    assert list(select_loci(iter(loci), None)) == loci


def test_select_loci_keeps_the_order_of_the_loci(tmp_path: Path) -> None:
    regions = read_regions(write_bed(tmp_path, ["chr2\t0\t10", "chr1\t0\t10"]), _CONTIGS)
    loci = [
        make_locus("chr1", 2, 3),
        make_locus("chr1", 50, 60),
        make_locus("chr2", 0, 1),
        make_locus("chr2", 9, 12),
    ]

    assert list(select_loci(iter(loci), regions)) == [loci[0], loci[2], loci[3]]


def test_unsorted_bed_with_overlapping_and_adjacent_lines(tmp_path: Path) -> None:
    # chr1: 300-400 and 100-200, 150-250 overlapping, 250-260 next to it: 100-260 and 300-400.
    lines = [
        "# a comment",
        "track name=regions",
        "browser position chr1:1-1000",
        "",
        "chr2\t0\t5",
        "chr1\t300\t400\tname",
        "chr1\t150\t250",
        "chr1\t100\t200",
        "chr1\t250\t260",
    ]
    regions = read_regions(write_bed(tmp_path, lines), _CONTIGS)

    overlaps = {
        (contig, start, end): regions.overlaps(contig, start, end)
        for contig, start, end in [
            ("chr1", 99, 100),
            ("chr1", 99, 101),
            ("chr1", 259, 260),
            ("chr1", 260, 300),
            ("chr1", 399, 500),
            ("chr1", 400, 500),
            ("chr2", 4, 5),
            ("chr2", 5, 6),
        ]
    }
    assert overlaps == {
        ("chr1", 99, 100): False,
        ("chr1", 99, 101): True,
        ("chr1", 259, 260): True,
        ("chr1", 260, 300): False,
        ("chr1", 399, 500): True,
        ("chr1", 400, 500): False,
        ("chr2", 4, 5): True,
        ("chr2", 5, 6): False,
    }


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
    ids=["two-columns", "start-not-whole", "end-not-whole", "end-not-after-start", "contig"],
)
def test_read_regions_refuses_a_line(tmp_path: Path, line: str, message: str) -> None:
    path = write_bed(tmp_path, ["chr1\t0\t10", line])

    with pytest.raises(PvcEvalInputError) as raised:
        read_regions(path, _CONTIGS)

    assert str(raised.value) == f"{path}, line 2: {message.format(path=path)}"


def test_read_regions_refuses_a_missing_file(tmp_path: Path) -> None:
    path = tmp_path / "absent.bed"

    with pytest.raises(PvcEvalInputError) as raised:
        read_regions(path, _CONTIGS)

    assert str(raised.value) == f"{path} does not exist. Check its path and run again."


def test_read_regions_refuses_a_gzipped_bed(tmp_path: Path) -> None:
    path = tmp_path / "regions.bed.gz"
    path.write_bytes(gzip.compress(b"chr1\t0\t10\n"))

    with pytest.raises(PvcEvalInputError) as raised:
        read_regions(path, _CONTIGS)

    assert str(raised.value) == (
        f"{path} is not UTF-8 text: a byte at or after line 1 cannot be read as UTF-8. If it"
        f" is compressed, decompress it with `gunzip -c {path} > file.bed` and give file.bed"
        " instead; if it is text in another encoding, convert it to UTF-8; and run again."
    )
