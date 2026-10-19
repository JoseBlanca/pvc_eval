"""The stratum of a locus: pvc_stratum for --strata pvc, and read_strata and Strata for a
--strata BED (doc/specs/loci.md "Stratum", doc/spec.md §5.7 and §7.2).

The spec's cases give positions 1-based and inclusive; loci and STR regions are built
0-based and half-open, and a BED line is 0-based and half-open, so the locus at 100 to 120
is (99, 120) and an interval at 90 to 101 is the BED line `chr1 89 101`.
"""

from pathlib import Path

import pytest

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.locus import Locus
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.strata import pvc_stratum, read_strata
from pvc_eval.vcf import VcfGenotype, VcfRecord

_CONTIGS = ("chr1", "chr2")

# An AT repeat at 11 to 20, the STR segment of §5.7, and the insertion of AT after 10.
_STR = GenomeSegment("chr1", 10, 20, RegionType.STR, "AT")
_BUNDLE = GenomeSegment("chr1", 40, 80, RegionType.STR_BUNDLE, None)
_SATELLITE = GenomeSegment("chr1", 100, 400, RegionType.SATELLITE, None)


def make_record(start: int, alleles: tuple[str, ...], span: tuple[int, int]) -> VcfRecord:
    """A VCF record on chr1 of one heterozygous sample, 0-based and half-open."""
    genotype = VcfGenotype(alleles=(0, 1), phased=False, phase_set=None, quality=None)
    return VcfRecord(
        contig="chr1",
        start=start,
        alleles=alleles,
        qual=None,
        genotypes=(genotype,),
        line_number=3,
        ref_mismatch=False,
        span=span,
    )


_INSERTION = make_record(9, ("A", "AAT"), (9, 20))


def make_locus(start: int, end: int, str_regions: tuple[GenomeSegment, ...]) -> Locus:
    return Locus("chr1", start, end, (_INSERTION,), (), str_regions)


def write_bed(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "strata.bed"
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("locus", "stratum"),
    [
        (make_locus(9, 20, (_STR,)), "str"),
        (make_locus(9, 20, ()), "generic"),
        (make_locus(9, 400, (_STR, _SATELLITE)), "satellite"),
        (make_locus(9, 80, (_STR, _BUNDLE)), "str-bundle"),
        (Locus("chr1", 40, 400, (), (), (_BUNDLE, _SATELLITE)), "satellite"),
    ],
    ids=[
        "insertion-before-str",
        "no-str-region",
        "str-chained-to-satellite",
        "str-chained-to-bundle",
        "bundle-chained-to-satellite",
    ],
)
def test_pvc_stratum(locus: Locus, stratum: str) -> None:
    assert pvc_stratum(locus) == stratum


@pytest.mark.parametrize(
    "str_regions",
    [
        (GenomeSegment("chr1", 10, 20, None, None),),
        (GenomeSegment("chr1", 10, 20, None, None), _SATELLITE),
        (_SATELLITE, GenomeSegment("chr1", 400, 410, None, None)),
    ],
    ids=["alone", "before-a-satellite", "after-a-satellite"],
)
def test_pvc_stratum_raises_for_an_str_region_with_no_type(
    str_regions: tuple[GenomeSegment, ...],
) -> None:
    locus = make_locus(9, str_regions[-1].end, str_regions)
    untyped = next(region for region in str_regions if region.region_type is None)

    with pytest.raises(PvcEvalInputError) as raised:
        pvc_stratum(locus)

    assert str(raised.value) == (
        f"The STR region chr1:{untyped.start + 1}-{untyped.end} has no region type, since the --str-regions BED has no"
        " fourth column there, and --strata pvc needs the type of every STR region. Give"
        " --str-regions a BED written by `pvc_eval str-regions`, which has the types, or"
        " leave --str-regions out so that the STR finder runs, and run again."
    )


def test_strata_joins_the_classes_of_the_intervals_the_locus_shares_positions_with(
    tmp_path: Path,
) -> None:
    lines = ["chr1\t89\t101\thomopolymer", "chr1\t114\t500\tsegdup"]
    strata = read_strata(write_bed(tmp_path, lines), _CONTIGS)

    assert strata.stratum(Locus("chr1", 99, 120, (), (), ())) == "homopolymer,segdup"


def test_strata_gives_a_dot_for_a_locus_under_no_interval(tmp_path: Path) -> None:
    lines = ["chr1\t89\t99\thomopolymer", "chr1\t120\t500\tsegdup", "chr2\t99\t120\tsegdup"]
    strata = read_strata(write_bed(tmp_path, lines), _CONTIGS)

    assert strata.stratum(Locus("chr1", 99, 120, (), (), ())) == "."


def test_strata_sorts_the_classes_and_gives_each_once(tmp_path: Path) -> None:
    lines = ["chr1\t0\t10\tzeta", "chr1\t20\t30\talpha", "chr1\t40\t50\tzeta"]
    strata = read_strata(write_bed(tmp_path, lines), _CONTIGS)

    assert strata.stratum(Locus("chr1", 0, 50, (), (), ())) == "alpha,zeta"


def test_strata_gives_a_dot_on_a_contig_with_no_bed_line(tmp_path: Path) -> None:
    strata = read_strata(write_bed(tmp_path, ["chr1\t0\t10\tzeta"]), _CONTIGS)

    assert strata.stratum(Locus("chr2", 0, 50, (), (), ())) == "."


def test_strata_from_an_unsorted_bed_with_overlapping_lines(tmp_path: Path) -> None:
    lines = [
        "# a comment",
        "track name=strata",
        "chr2\t0\t50\tsegdup",
        "chr1\t300\t400\tsegdup",
        "chr1\t0\t200\tlowmap",
        "chr1\t100\t250\tsegdup",
        "chr1\t150\t160\tlowmap",
        "chr1\t190\t260\tgc_high\textra",
    ]
    strata = read_strata(write_bed(tmp_path, lines), _CONTIGS)

    stratum_by_span = {
        span: strata.stratum(Locus(contig, *span, (), (), ()))
        for contig, span in [
            ("chr1", (0, 100)),
            ("chr1", (150, 155)),
            ("chr1", (199, 200)),
            ("chr1", (200, 201)),
            ("chr1", (255, 300)),
            ("chr1", (260, 300)),
            ("chr1", (250, 350)),
            ("chr2", (49, 60)),
        ]
    }
    assert stratum_by_span == {
        (0, 100): "lowmap",
        (150, 155): "lowmap,segdup",
        (199, 200): "gc_high,lowmap,segdup",
        (200, 201): "gc_high,segdup",
        (255, 300): "gc_high",
        (260, 300): ".",
        (250, 350): "gc_high,segdup",
        (49, 60): "segdup",
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
            "chr1\tx\t2000\tsegdup",
            "the start, x, is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t100\t2000.5\tsegdup",
            "the end, 2000.5, is not a whole number. A BED line gives its start and its end as"
            " whole numbers of 0 or more; correct the line and run again.",
        ),
        (
            "chr1\t100\t50\tsegdup",
            "the end, 50, is not after the start, 100. A BED line covers the positions from"
            " its start up to its end, 0-based and half-open, so its end is greater than its"
            " start; correct the line and run again.",
        ),
        (
            "chrUn\t100\t200\tsegdup",
            "the contig chrUn is not in the reference. Keep only the reference's contigs in the"
            " BED, for instance with `awk -F'\\t' 'NR == FNR {{keep[$1]; next}} $1 in keep'"
            " ref.fa.fai {path} > kept.bed`, where ref.fa.fai is the .fai of the --reference"
            " FASTA, and run again with kept.bed.",
        ),
        (
            "chr1\t100\t200",
            "has no fourth column, which names the class of its interval. Give each line of"
            " the --strata BED its class and run again.",
        ),
        (
            "chr1\t100\t200\t",
            "has no fourth column, which names the class of its interval. Give each line of"
            " the --strata BED its class and run again.",
        ),
    ],
    ids=[
        "two-columns",
        "start-not-whole",
        "end-not-whole",
        "end-not-after-start",
        "contig",
        "no-fourth-column",
        "empty-fourth-column",
    ],
)
def test_read_strata_refuses_a_line(tmp_path: Path, line: str, message: str) -> None:
    path = write_bed(tmp_path, ["chr1\t0\t10\tsegdup", line])

    with pytest.raises(PvcEvalInputError) as raised:
        read_strata(path, _CONTIGS)

    assert str(raised.value) == f"{path}, line 2: {message.format(path=path)}"
