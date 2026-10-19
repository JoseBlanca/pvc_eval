"""The loci of build_loci (doc/specs/loci.md, "Building the loci"; doc/spec.md §5.3).

The cases are those of doc/specs/loci.md, "How it is verified". The spec writes positions
1-based and inclusive; here they are 0-based and half-open, so the spec's "SNP at 100" is a
span (99, 100) and its "11-20" is (10, 20).
"""

import tracemalloc
from collections.abc import Iterable, Iterator, Sequence
from itertools import pairwise

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pvc_eval.loci import DEFAULT_CHAIN_DISTANCE, build_loci
from pvc_eval.locus import Locus
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.vcf import VcfRecord

CONTIGS = ("chr1", "chr2", "chr3")


def make_record(
    start: int,
    ref: str,
    alts: tuple[str, ...],
    span: tuple[int, int],
    *,
    contig: str = "chr1",
    line_number: int = 1,
) -> VcfRecord:
    """A VCF record with no sample, 0-based and half-open."""
    return VcfRecord(
        contig=contig,
        start=start,
        alleles=(ref, *alts),
        qual=None,
        genotypes=(),
        line_number=line_number,
        ref_mismatch=False,
        span=span,
    )


def make_snp(start: int, *, contig: str = "chr1", line_number: int = 1) -> VcfRecord:
    """A SNP at the 0-based position start, whose span is its REF."""
    return make_record(
        start, "A", ("G",), (start, start + 1), contig=contig, line_number=line_number
    )


def make_no_alt(
    start: int, length: int, *, contig: str = "chr1", line_number: int = 1
) -> VcfRecord:
    """A record with no ALT whose REF, and so its span, is [start, start + length)."""
    return make_record(
        start, "A" * length, (), (start, start + length), contig=contig, line_number=line_number
    )


def make_region(
    start: int, end: int, *, contig: str = "chr1", region_type: RegionType = RegionType.STR
) -> GenomeSegment:
    """An STR region, or a segment of another type, [start, end)."""
    motif = "AT" if region_type is RegionType.STR else None
    return GenomeSegment(contig, start, end, region_type, motif)


def make_locus(
    start: int,
    end: int,
    records_a: tuple[VcfRecord, ...] = (),
    records_b: tuple[VcfRecord, ...] = (),
    str_regions: tuple[GenomeSegment, ...] = (),
    *,
    contig: str = "chr1",
) -> Locus:
    return Locus(contig, start, end, records_a, records_b, str_regions)


def run_build_loci(
    records_a: Iterable[VcfRecord],
    records_b: Iterable[VcfRecord] = (),
    segments: Iterable[GenomeSegment] = (),
    *,
    chain_distance: int = DEFAULT_CHAIN_DISTANCE,
    contigs: Sequence[str] = CONTIGS,
) -> list[Locus]:
    return list(
        build_loci(
            iter(records_a),
            iter(records_b),
            iter(segments),
            contigs,
            chain_distance=chain_distance,
        )
    )


def test_default_chain_distance_is_ten() -> None:
    assert DEFAULT_CHAIN_DISTANCE == 10


# The chain distance at its boundaries. The spec's locus ending at 100 takes a SNP at 111
# (gap 10) and leaves one at 112 (gap 11); with 0, SNPs at 100 and 101 chain and at 100 and
# 102 do not.
SNP_100 = make_snp(99, line_number=1)
SNP_101 = make_snp(100, line_number=2)
SNP_102 = make_snp(101, line_number=2)
SNP_111 = make_snp(110, line_number=2)
SNP_112 = make_snp(111, line_number=2)


@pytest.mark.parametrize(
    ("records_a", "chain_distance", "expected"),
    [
        ((SNP_100, SNP_111), 10, [make_locus(99, 111, (SNP_100, SNP_111))]),
        (
            (SNP_100, SNP_112),
            10,
            [make_locus(99, 100, (SNP_100,)), make_locus(111, 112, (SNP_112,))],
        ),
        ((SNP_100, SNP_101), 0, [make_locus(99, 101, (SNP_100, SNP_101))]),
        (
            (SNP_100, SNP_102),
            0,
            [make_locus(99, 100, (SNP_100,)), make_locus(101, 102, (SNP_102,))],
        ),
    ],
    ids=["gap_10_chains", "gap_11_does_not", "adjacent_chain_at_0", "gap_1_does_not_at_0"],
)
def test_records_chain_up_to_the_chain_distance(
    records_a: tuple[VcfRecord, ...], chain_distance: int, expected: list[Locus]
) -> None:
    assert run_build_loci(records_a, chain_distance=chain_distance) == expected


# The STR table: an AT repeat at 11-20.
AT_REPEAT = make_region(10, 20)
SECOND_REPEAT = make_region(20, 30)
SNP_9 = make_snp(8)
SNP_10 = make_snp(9)
SNP_21 = make_snp(20)


@pytest.mark.parametrize(
    ("records_a", "regions", "expected"),
    [
        ((SNP_10,), (AT_REPEAT,), [make_locus(9, 20, (SNP_10,), (), (AT_REPEAT,))]),
        ((SNP_21,), (AT_REPEAT,), [make_locus(10, 21, (SNP_21,), (), (AT_REPEAT,))]),
        ((SNP_9,), (AT_REPEAT,), [make_locus(8, 9, (SNP_9,))]),
        (
            (SNP_10,),
            (AT_REPEAT, SECOND_REPEAT),
            [make_locus(9, 20, (SNP_10,), (), (AT_REPEAT,))],
        ),
    ],
    ids=[
        "touched_at_its_start",
        "touched_at_its_end",
        "one_position_between",
        "region_next_to_a_touched_one_stays_out",
    ],
)
def test_a_record_widens_its_locus_to_the_str_regions_it_touches(
    records_a: tuple[VcfRecord, ...], regions: tuple[GenomeSegment, ...], expected: list[Locus]
) -> None:
    assert run_build_loci(records_a, (), regions) == expected


def test_a_locus_holds_an_str_region_in_a_gap_between_its_records() -> None:
    first = make_snp(9, line_number=1)
    second = make_snp(17, line_number=2)
    homopolymer = make_region(11, 16)

    loci = run_build_loci((first, second), (), (homopolymer,))

    assert loci == [make_locus(9, 18, (first, second), (), (homopolymer,))]


# The table of a locus that reaches back: a SNP at 100 and a satellite at 105-3000 that it
# does not touch, with a chain distance of 10.
def make_satellite(start: int, end: int) -> GenomeSegment:
    return make_region(start, end, region_type=RegionType.SATELLITE)


SATELLITE = make_satellite(104, 3000)
SNP_2990 = make_snp(2989, line_number=2)
SNP_3001 = make_snp(3000, line_number=2)
SNP_3005 = make_snp(3004, line_number=2)
SATELLITE_AT_112 = make_satellite(111, 3000)
SATELLITE_AT_205 = make_satellite(204, 3000)
SATELLITE_AT_101 = make_satellite(100, 200)
HOMOPOLYMER_AT_102 = make_region(101, 103, region_type=RegionType.STR)


@pytest.mark.parametrize(
    ("records_a", "regions", "expected"),
    [
        (
            (SNP_100, SNP_2990),
            (SATELLITE,),
            [make_locus(99, 3000, (SNP_100, SNP_2990), (), (SATELLITE,))],
        ),
        (
            (SNP_100, SNP_3001),
            (SATELLITE,),
            [make_locus(99, 3001, (SNP_100, SNP_3001), (), (SATELLITE,))],
        ),
        (
            (SNP_100, SNP_3005),
            (SATELLITE,),
            [make_locus(99, 100, (SNP_100,)), make_locus(3004, 3005, (SNP_3005,))],
        ),
        (
            (SNP_100, SNP_2990),
            (SATELLITE_AT_112,),
            [
                make_locus(99, 100, (SNP_100,)),
                make_locus(111, 3000, (SNP_2990,), (), (SATELLITE_AT_112,)),
            ],
        ),
        (
            (SNP_100, SNP_2990),
            (SATELLITE_AT_101, SATELLITE_AT_205),
            [make_locus(99, 3000, (SNP_100, SNP_2990), (), (SATELLITE_AT_101, SATELLITE_AT_205))],
        ),
        (
            (SNP_100, SNP_2990),
            (HOMOPOLYMER_AT_102, SATELLITE),
            [make_locus(99, 3000, (SNP_100, SNP_2990), (), (HOMOPOLYMER_AT_102, SATELLITE))],
        ),
    ],
    ids=[
        "record_inside_the_satellite",
        "record_touching_its_end",
        "no_record_touches_it",
        "satellite_beyond_the_chain_distance",
        "distance_measured_from_the_locus_not_the_snp",
        "region_read_only_with_the_next_record",
    ],
)
def test_a_locus_that_reaches_back_joins_the_one_before_it(
    records_a: tuple[VcfRecord, ...], regions: tuple[GenomeSegment, ...], expected: list[Locus]
) -> None:
    assert run_build_loci(records_a, (), regions) == expected


# The table of records with no ALT: `100 CAGTCAGTCAGTC .`, 100-112, in file B.
NO_ALT_100 = make_no_alt(99, 13, line_number=5)
SNP_99 = make_snp(98, line_number=1)
SNP_113 = make_snp(112, line_number=2)
SNP_110 = make_snp(109, line_number=2)


@pytest.mark.parametrize(
    ("records_a", "expected"),
    [
        ((SNP_99, SNP_113), [make_locus(98, 99, (SNP_99,)), make_locus(112, 113, (SNP_113,))]),
        (
            (SNP_100, SNP_112),
            [
                make_locus(99, 100, (SNP_100,), (NO_ALT_100,)),
                make_locus(111, 112, (SNP_112,), (NO_ALT_100,)),
            ],
        ),
        ((SNP_100, SNP_110), [make_locus(99, 110, (SNP_100, SNP_110), (NO_ALT_100,))]),
        ((), []),
    ],
    ids=["overlaps_neither", "in_both", "reaches_past_the_locus", "no_record_with_an_alt"],
)
def test_a_record_with_no_alt_is_in_every_locus_its_ref_overlaps(
    records_a: tuple[VcfRecord, ...], expected: list[Locus]
) -> None:
    assert run_build_loci(records_a, (NO_ALT_100,)) == expected


def test_a_record_with_no_alt_in_the_gap_between_chained_snps_is_in_their_locus() -> None:
    first = make_snp(9, line_number=1)
    no_alt = make_no_alt(11, 1, line_number=2)
    second = make_snp(14, line_number=3)

    loci = run_build_loci((first, no_alt, second))

    assert loci == [make_locus(9, 15, (first, no_alt, second))]


def test_a_record_whose_only_alt_is_star_starts_a_locus() -> None:
    star = make_record(9, "A", ("*",), (9, 10))

    assert run_build_loci((star,)) == [make_locus(9, 10, (star,))]


def test_an_str_region_at_the_start_of_the_next_contig_is_not_dropped() -> None:
    snp_chr1 = make_snp(29, contig="chr1", line_number=1)
    region_chr2 = make_region(0, 20, contig="chr2")
    snp_chr2 = make_snp(4, contig="chr2", line_number=2)

    loci = run_build_loci((snp_chr1, snp_chr2), (), (region_chr2,))

    assert loci == [
        make_locus(29, 30, (snp_chr1,)),
        make_locus(0, 20, (snp_chr2,), (), (region_chr2,), contig="chr2"),
    ]


def test_an_str_region_on_the_next_contig_is_read_past_the_ones_of_the_first() -> None:
    snp_chr1 = make_snp(9, contig="chr1", line_number=1)
    far_region_chr1 = make_region(4999, 5010, contig="chr1")
    region_chr2 = make_region(89, 99, contig="chr2")
    snp_chr2 = make_snp(99, contig="chr2", line_number=2)

    loci = run_build_loci((snp_chr1, snp_chr2), (), (far_region_chr1, region_chr2))

    assert loci == [
        make_locus(9, 10, (snp_chr1,)),
        make_locus(89, 100, (snp_chr2,), (), (region_chr2,), contig="chr2"),
    ]


def test_a_locus_never_spans_two_contigs() -> None:
    last_of_chr1 = make_snp(999, contig="chr1", line_number=1)
    first_of_chr2 = make_snp(0, contig="chr2", line_number=2)

    loci = run_build_loci((last_of_chr1, first_of_chr2), chain_distance=1000)

    assert loci == [
        make_locus(999, 1000, (last_of_chr1,)),
        make_locus(0, 1, (first_of_chr2,), contig="chr2"),
    ]


def test_a_contig_with_records_of_one_file_gives_loci_with_none_of_the_other() -> None:
    snp_a_chr1 = make_snp(9, contig="chr1", line_number=1)
    snp_a_chr2 = make_snp(9, contig="chr2", line_number=2)
    snp_b_chr1 = make_snp(14, contig="chr1", line_number=1)

    loci = run_build_loci((snp_a_chr1, snp_a_chr2), (snp_b_chr1,))

    assert loci == [
        make_locus(9, 15, (snp_a_chr1,), (snp_b_chr1,)),
        make_locus(9, 10, (snp_a_chr2,), (), contig="chr2"),
    ]


class CountedIterator[T]:
    """An iterator that counts what is read from it, and keeps the last item read."""

    def __init__(self, items: Iterable[T]) -> None:
        self._items = iter(items)
        self.count = 0
        self.last: T | None = None

    def __iter__(self) -> Iterator[T]:
        return self

    def __next__(self) -> T:
        item = next(self._items)
        self.count += 1
        self.last = item
        return item


@pytest.mark.parametrize(
    ("records_a", "records_b", "expected", "segments_read"),
    [
        ((), (), [], 0),
        (
            (SNP_100, SNP_112),
            (),
            [make_locus(99, 100, (SNP_100,)), make_locus(111, 112, (SNP_112,))],
            1,
        ),
        ((), (SNP_100,), [make_locus(99, 100, (), (SNP_100,))], 1),
        ((NO_ALT_100, make_no_alt(200, 3, line_number=6)), (), [], 1),
    ],
    ids=["both_empty", "file_b_empty", "file_a_empty", "only_records_with_no_alt"],
)
def test_empty_inputs_read_no_segment_they_do_not_need(
    records_a: tuple[VcfRecord, ...],
    records_b: tuple[VcfRecord, ...],
    expected: list[Locus],
    segments_read: int,
) -> None:
    # One STR region on chr3, after every record: it is read only to know that it is past.
    segments = CountedIterator([make_region(0, 10, contig="chr3")])

    loci = list(build_loci(iter(records_a), iter(records_b), segments, CONTIGS))

    assert loci == expected
    assert segments.count == segments_read


def test_generic_segments_between_the_str_regions_change_no_locus() -> None:
    records_a = (SNP_100, SNP_2990)
    str_regions = (HOMOPOLYMER_AT_102, SATELLITE)
    segments = (
        make_region(0, 101, region_type=RegionType.GENERIC),
        HOMOPOLYMER_AT_102,
        make_region(103, 104, region_type=RegionType.GENERIC),
        SATELLITE,
        make_region(3000, 5000, region_type=RegionType.GENERIC),
    )
    expected = [make_locus(99, 3000, records_a, (), str_regions)]

    assert run_build_loci(records_a, (), str_regions) == expected
    assert run_build_loci(records_a, (), segments) == expected


# What build_loci checks: each fault raises ValueError when the stream reaches it.
@pytest.mark.parametrize(
    ("records_a", "records_b", "segments", "message"),
    [
        ((make_snp(9, contig="chrX"),), (), (), "file A.*chrX"),
        ((make_snp(9),), (), (make_region(0, 5, contig="chrX"),), "segments.*chrX"),
        (
            (make_snp(9, contig="chr2"), make_snp(19, contig="chr1")),
            (),
            (),
            "file A.*chr1",
        ),
        ((make_snp(9),), (make_snp(29), make_snp(19)), (), "file B.*19"),
        (
            (make_snp(9, contig="chr2"),),
            (),
            (make_region(0, 5, contig="chr2"), make_region(30, 35, contig="chr1")),
            "segments.*chr1 at start 30 comes after",
        ),
        (
            (make_snp(99),),
            (),
            (make_region(30, 35), make_region(20, 25)),
            "segments.*start 20 comes after",
        ),
        (
            (make_snp(99),),
            (),
            (make_region(20, 35), make_region(30, 40)),
            "segments.*30-40 overlaps",
        ),
    ],
    ids=[
        "record_on_an_unknown_contig",
        "segment_on_an_unknown_contig",
        "record_on_an_earlier_contig",
        "record_with_an_earlier_span_start",
        "str_region_on_an_earlier_contig",
        "str_region_with_an_earlier_start",
        "str_region_overlapping_the_one_before",
    ],
)
def test_inputs_out_of_order_or_on_an_unknown_contig_raise(
    records_a: tuple[VcfRecord, ...],
    records_b: tuple[VcfRecord, ...],
    segments: tuple[GenomeSegment, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        run_build_loci(records_a, records_b, segments)


@pytest.mark.parametrize("chain_distance", [-1, 1.5, True], ids=["negative", "float", "bool"])
def test_a_chain_distance_that_is_not_a_whole_number_raises_at_the_first_next(
    chain_distance: int,
) -> None:
    loci = build_loci(iter(()), iter(()), iter(()), CONTIGS, chain_distance=chain_distance)

    with pytest.raises(ValueError, match="chain_distance"):
        next(loci)


# The cases doc/spec.md §7.2 lists for loci.
def test_chaining_runs_across_both_files() -> None:
    first_a = make_snp(9, line_number=1)
    middle_b = make_snp(15, line_number=1)
    last_a = make_snp(22, line_number=2)

    loci = run_build_loci((first_a, last_a), (middle_b,))

    # The two SNPs of file A leave a gap of 12; each leaves less than 10 to file B's.
    assert loci == [make_locus(9, 23, (first_a, last_a), (middle_b,))]


# The traded events of doc/spec.md §5.3, reference 5-10 AGAACA: spans 5-6 and 9 in file A,
# 6 and 8-9 in file B.
TRADED_A_DELETION = make_record(4, "AG", ("A",), (4, 6), line_number=1)
TRADED_A_SNP = make_snp(8, line_number=2)
TRADED_B_SNP = make_snp(5, line_number=1)
TRADED_B_DELETION = make_record(7, "AC", ("A",), (7, 9), line_number=2)


@pytest.mark.parametrize(
    ("chain_distance", "expected"),
    [
        (
            1,
            [
                make_locus(
                    4,
                    9,
                    (TRADED_A_DELETION, TRADED_A_SNP),
                    (TRADED_B_SNP, TRADED_B_DELETION),
                )
            ],
        ),
        (
            0,
            [
                make_locus(4, 6, (TRADED_A_DELETION,), (TRADED_B_SNP,)),
                make_locus(7, 9, (TRADED_A_SNP,), (TRADED_B_DELETION,)),
            ],
        ),
    ],
    ids=["one_locus_at_1", "two_loci_at_0"],
)
def test_traded_events_share_a_locus_from_a_chain_distance_of_one(
    chain_distance: int, expected: list[Locus]
) -> None:
    loci = run_build_loci(
        (TRADED_A_DELETION, TRADED_A_SNP),
        (TRADED_B_SNP, TRADED_B_DELETION),
        chain_distance=chain_distance,
    )

    assert loci == expected


def test_adjacent_snps_of_the_two_files_chain_at_0() -> None:
    snp_a = make_snp(9)
    snp_b = make_snp(10)

    loci = run_build_loci((snp_a,), (snp_b,), chain_distance=0)

    assert loci == [make_locus(9, 11, (snp_a,), (snp_b,))]


def test_a_record_with_no_alt_starts_no_locus() -> None:
    no_alt = make_no_alt(9, 1, line_number=1)
    snp = make_snp(30, line_number=2)

    loci = run_build_loci((no_alt, snp))

    assert loci == [make_locus(30, 31, (snp,))]


def test_an_insertion_next_to_a_repeat_pulls_the_repeat_into_its_locus() -> None:
    # `9 G GC` spans 9-10, the base after its anchor included; the repeat is at 11-20.
    insertion = make_record(8, "G", ("GC",), (8, 10))

    loci = run_build_loci((insertion,), (), (AT_REPEAT,))

    assert loci == [make_locus(8, 20, (insertion,), (), (AT_REPEAT,))]


def test_complex_records_at_the_two_ends_of_a_repeat_share_its_locus() -> None:
    # `9 GA TT` (9-10) touches the repeat at 11-20 at its start, `21 TC GG` (21-22) at its
    # end. At a chain distance of 0 only the repeat joins them.
    complex_a = make_record(8, "GA", ("TT",), (8, 10))
    complex_b = make_record(20, "TC", ("GG",), (20, 22))

    loci = run_build_loci((complex_a,), (complex_b,), (AT_REPEAT,), chain_distance=0)

    assert loci == [make_locus(8, 22, (complex_a,), (complex_b,), (AT_REPEAT,))]


# What build_loci reads.
def test_build_loci_reads_one_record_past_the_locus_it_yields() -> None:
    snps = CountedIterator(make_snp(100 * index) for index in range(100_000))
    loci_count = 0

    for index, _locus in enumerate(build_loci(snps, iter(()), iter(()), CONTIGS)):
        loci_count = index + 1
        assert snps.count <= loci_count + 1

    assert loci_count == 100_000


def make_guarded_segments(
    segments: Iterable[GenomeSegment], records: CountedIterator[VcfRecord]
) -> Iterator[GenomeSegment]:
    """The segments, raising when one is asked for after the one read before it already
    came after the end of the span of the last record read."""
    previous: GenomeSegment | None = None
    for segment in segments:
        last_record = records.last
        if last_record is None:
            raise AssertionError("a segment was read before any record")
        if previous is not None and (CONTIGS.index(previous.contig), previous.start) > (
            CONTIGS.index(last_record.contig),
            last_record.span[1],
        ):
            raise AssertionError(f"read past {previous} for {last_record.span}")
        yield segment
        previous = segment


def test_build_loci_reads_the_segments_only_as_far_as_the_records_need() -> None:
    # Each contig holds a GENERIC segment and then an STR region every 20 positions.
    segments = [
        make_region(start, start + length, contig=contig, region_type=region_type)
        for contig in CONTIGS
        for block in range(0, 1000, 20)
        for start, length, region_type in (
            (block, 15, RegionType.GENERIC),
            (block + 15, 5, RegionType.STR),
        )
    ]
    records = CountedIterator(
        [
            make_snp(10, contig="chr1", line_number=1),
            make_snp(300, contig="chr1", line_number=2),
            make_snp(40, contig="chr2", line_number=3),
        ]
    )
    read_segments = CountedIterator(make_guarded_segments(segments, records))

    loci = list(build_loci(records, iter(()), read_segments, CONTIGS))

    assert len(loci) == 3
    # All 100 of chr1, and those of chr2 up to the first that starts after 41: at 0, 15,
    # 20, 35, 40 and 55.
    assert read_segments.count == 106


def make_adjacent_records_with_no_alt(count: int) -> Iterator[VcfRecord]:
    for position in range(count):
        yield make_no_alt(position, 1, line_number=position + 1)


def find_peak_memory(count: int) -> tuple[int, list[Locus]]:
    """The peak tracemalloc traces while build_loci reads count records with no ALT, one
    position apart, and then a SNP of file B over the last of them."""
    snp = make_snp(count - 1)
    tracemalloc.start()
    try:
        loci = list(
            build_loci(make_adjacent_records_with_no_alt(count), iter((snp,)), iter(()), CONTIGS)
        )
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak, loci


def test_records_with_no_alt_are_not_kept_once_passed() -> None:
    peak_10_000, _ = find_peak_memory(10_000)
    peak_100_000, loci = find_peak_memory(100_000)

    assert abs(peak_100_000 - peak_10_000) < 0.1 * peak_10_000, (peak_10_000, peak_100_000)
    assert loci == [
        make_locus(
            99_999, 100_000, (make_no_alt(99_999, 1, line_number=100_000),), (make_snp(99_999),)
        )
    ]


# A property: build_loci against a slow version of doc/spec.md §5.3 written from its rules.
CONTIG_LENGTH = 60


def find_gap(first: tuple[int, int], second: tuple[int, int]) -> int:
    """The positions between two stretches, 0-based and half-open: 0 when they are next to
    each other, and below 0 when they overlap."""
    return max(first[0], second[0]) - min(first[1], second[1])


def is_touching(record: VcfRecord, region: GenomeSegment) -> bool:
    span_start, span_end = record.span
    return span_start <= region.end and region.start <= span_end


def is_overlapping(first: tuple[int, int], second: tuple[int, int]) -> bool:
    return first[0] < second[1] and second[0] < first[1]


def has_alt(record: VcfRecord) -> bool:
    return len(record.alleles) > 1


def build_loci_slowly(
    records_a: Sequence[VcfRecord],
    records_b: Sequence[VcfRecord],
    segments: Sequence[GenomeSegment],
    contigs: Sequence[str],
    chain_distance: int,
) -> list[Locus]:
    """doc/spec.md §5.3, steps 1 to 4, applied until nothing changes. Every record has a
    line number of its own."""
    loci = []
    for contig in contigs:
        regions = [
            segment
            for segment in segments
            if segment.contig == contig and segment.region_type is not RegionType.GENERIC
        ]
        left = [
            record
            for record in (*records_a, *records_b)
            if record.contig == contig and has_alt(record)
        ]
        built: list[tuple[int, int, list[VcfRecord]]] = []
        while left:
            # Step 1: the leftmost record not yet in a locus.
            first = min(left, key=lambda record: record.span[0])
            left.remove(first)
            (start, end), members = first.span, [first]
            changed = True
            while changed:
                changed = False
                # Step 2: every record that overlaps the locus or chains to it.
                for record in list(left):
                    if find_gap(record.span, (start, end)) <= chain_distance:
                        left.remove(record)
                        members.append(record)
                        start, end = min(start, record.span[0]), max(end, record.span[1])
                        changed = True
                # Step 3: every STR region a record of the locus touches.
                for region in regions:
                    touched = any(is_touching(record, region) for record in members)
                    if touched and (region.start < start or region.end > end):
                        start, end = min(start, region.start), max(end, region.end)
                        changed = True
                # Step 4: a locus within the chain distance of one built before it.
                for before in list(built):
                    if find_gap((before[0], before[1]), (start, end)) <= chain_distance:
                        built.remove(before)
                        members.extend(before[2])
                        start, end = min(start, before[0]), max(end, before[1])
                        changed = True
            built.append((start, end, members))
        for start, end, members in sorted(built, key=lambda locus: locus[0]):
            in_locus = [
                record
                for record in (*records_a, *records_b)
                if record in members
                or (
                    record.contig == contig
                    and not has_alt(record)
                    and is_overlapping(record.span, (start, end))
                )
            ]
            loci.append(
                Locus(
                    contig,
                    start,
                    end,
                    tuple(record for record in records_a if record in in_locus),
                    tuple(record for record in records_b if record in in_locus),
                    tuple(
                        region
                        for region in regions
                        if is_overlapping((region.start, region.end), (start, end))
                    ),
                )
            )
    return loci


@st.composite
def draw_segments(draw: st.DrawFn, contig: str) -> list[GenomeSegment]:
    """STR regions of one contig, in order and apart, with GENERIC segments between them
    or not."""
    with_generic = draw(st.booleans())
    segments = []
    position = 0
    for gap, length, region_type in draw(
        st.lists(
            st.tuples(
                st.integers(0, 15),
                st.integers(1, 15),
                st.sampled_from([RegionType.STR, RegionType.STR_BUNDLE, RegionType.SATELLITE]),
            ),
            max_size=8,
        )
    ):
        start = position + gap
        if start + length > CONTIG_LENGTH:
            break
        if with_generic and gap > 0:
            segments.append(
                make_region(position, start, contig=contig, region_type=RegionType.GENERIC)
            )
        segments.append(make_region(start, start + length, contig=contig, region_type=region_type))
        position = start + length
    return segments


@st.composite
def draw_records(draw: st.DrawFn, first_line_number: int) -> list[VcfRecord]:
    """Records of one file on chr1 and chr2, with an ALT and without, in the order of
    parse_vcf: by contig and then by span start."""
    records = []
    drawn = draw(
        st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(0, CONTIG_LENGTH - 1),
                st.sampled_from([1, 1, 1, 2, 3, 5, 12, 25]),
                st.booleans(),
            ),
            max_size=10,
        )
    )
    for offset, (contig, span_start, length, with_alt) in enumerate(drawn):
        span = (span_start, min(CONTIG_LENGTH, span_start + length))
        line_number = first_line_number + offset
        if with_alt:
            records.append(
                make_record(span_start, "A", ("G",), span, contig=contig, line_number=line_number)
            )
        else:
            records.append(
                make_no_alt(
                    span_start, span[1] - span_start, contig=contig, line_number=line_number
                )
            )
    return sorted(records, key=lambda record: (CONTIGS.index(record.contig), record.span[0]))


def check_invariants(
    loci: Sequence[Locus],
    records: Sequence[VcfRecord],
    chain_distance: int,
) -> None:
    """The invariants of doc/specs/loci.md, "A property"."""
    for record in records:
        holding = [locus for locus in loci if record in (*locus.records_a, *locus.records_b)]
        if has_alt(record):
            assert len(holding) == 1, record
            assert holding[0].start <= record.span[0]
            assert record.span[1] <= holding[0].end
        else:
            overlapped = [
                locus
                for locus in loci
                if locus.contig == record.contig
                and is_overlapping(record.span, (locus.start, locus.end))
            ]
            assert holding == overlapped, record
    for before, after in pairwise(loci):
        assert CONTIGS.index(before.contig) <= CONTIGS.index(after.contig)
        if before.contig == after.contig:
            assert after.start - before.end > chain_distance


@settings(max_examples=500, deadline=None)
@given(
    records_a=draw_records(1),
    records_b=draw_records(101),
    segments_chr1=draw_segments("chr1"),
    segments_chr2=draw_segments("chr2"),
    chain_distance=st.integers(0, 20),
)
def test_build_loci_gives_the_loci_of_the_rules(
    records_a: list[VcfRecord],
    records_b: list[VcfRecord],
    segments_chr1: list[GenomeSegment],
    segments_chr2: list[GenomeSegment],
    chain_distance: int,
) -> None:
    segments = segments_chr1 + segments_chr2

    loci = run_build_loci(records_a, records_b, segments, chain_distance=chain_distance)

    assert loci == build_loci_slowly(records_a, records_b, segments, CONTIGS, chain_distance)
    check_invariants(loci, records_a + records_b, chain_distance)
