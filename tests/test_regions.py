"""The selection of --regions: read_regions, Regions and select_loci (doc/specs/loci.md
"Regions", doc/spec.md §3 and §7.2).

The spec's cases give positions 1-based and inclusive; loci are built 0-based and
half-open, and a BED line is 0-based and half-open, so the locus at 100 to 120 is (99, 120)
and the region at 110 to 200 is the BED line `chr1 109 200`.
"""

from pathlib import Path

from pvc_eval.locus import Locus
from pvc_eval.regions import Regions, read_regions, select_loci
from pvc_eval.vcf import VcfGenotype, VcfRecord

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


def test_record_with_no_alt_reaching_past_the_locus_does_not_count(tmp_path: Path) -> None:
    # The locus at 100 to 120 holds a record with no ALT whose REF covers 100 to 130; the
    # region at 125 to 200 shares positions with that REF only.
    regions = read_regions(write_bed(tmp_path, ["chr1\t124\t200"]), _CONTIGS)
    genotype = VcfGenotype(alleles=(0, 0), phased=False, phase_set=None, quality=None)
    no_alt = VcfRecord(
        contig="chr1",
        start=99,
        alleles=("A" * 31,),
        qual=None,
        genotypes=(genotype,),
        line_number=3,
        ref_mismatch=False,
        span=(99, 130),
    )
    locus = Locus("chr1", 99, 120, (no_alt,), (), ())

    assert list(select_loci(iter([locus]), regions)) == []


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


def test_region_nested_inside_another(tmp_path: Path) -> None:
    regions = read_regions(write_bed(tmp_path, ["chr1\t100\t300", "chr1\t150\t160"]), _CONTIGS)

    assert regions.overlaps("chr1", 200, 210) is True


def test_contig_with_no_bed_line_has_no_region(tmp_path: Path) -> None:
    regions = read_regions(write_bed(tmp_path, ["chr1\t0\t10"]), _CONTIGS)

    assert regions.overlaps("chr2", 0, 50) is False


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


def test_regions_sorts_and_merges_what_it_is_given() -> None:
    regions = Regions({"chr1": [(300, 400), (100, 300), (150, 160)]})

    assert [regions.overlaps("chr1", 150, 160), regions.overlaps("chr1", 200, 210)] == [True, True]
