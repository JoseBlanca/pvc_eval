"""The loci, stage 3 of doc/spec.md §4: the VCF records of file A and file B and the STR
regions of the segmented genome, merged into the stretches of reference that are compared
(doc/spec.md §5.3, doc/specs/loci.md "Building the loci").

Positions are 0-based and half-open, as VcfRecord.span and GenomeSegment are. Across
contigs a position is compared as the pair of its contig's rank in `contigs` and the
position, so that every position of one contig comes before those of the next.

Public: build_loci, which yields the loci; DEFAULT_CHAIN_DISTANCE.
"""

import heapq
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from enum import IntEnum

from pvc_eval.locus import Locus
from pvc_eval.segments import GenomeSegment, RegionType
from pvc_eval.vcf import VcfRecord

DEFAULT_CHAIN_DISTANCE = 10  # doc/spec.md §5.3, kept by the owner on 19 October 2026

type _Position = tuple[int, int]  # the rank of a contig in `contigs`, and a position on it
type _RankedSegment = tuple[int, GenomeSegment]  # a segment and the rank of its contig


class _SourceFile(IntEnum):
    """The file a record comes from; file A is read before file B where both are equal."""

    A = 0
    B = 1


@dataclass(frozen=True, slots=True)
class _RankedRecord:
    """A VCF record as build_loci reads it: with the rank of its contig, its file and its
    place in that file's stream."""

    rank: int
    source_file: _SourceFile
    index_in_file: int
    vcf_record: VcfRecord

    @property
    def start_position(self) -> _Position:
        """The rank and the 0-based start of the span."""
        return (self.rank, self.vcf_record.span[0])

    @property
    def end_position(self) -> _Position:
        """The rank and the end of the span, one past its last position."""
        return (self.rank, self.vcf_record.span[1])


@dataclass(slots=True)
class _OpenLocus:
    """The locus being built, `[start, end)`, with its records with an ALT. It is changed
    as it widens, and turned into a Locus once, when it is yielded."""

    rank: int
    contig: str
    start: int
    end: int
    records: list[_RankedRecord]


def build_loci(
    records_a: Iterator[VcfRecord],
    records_b: Iterator[VcfRecord],
    segments: Iterator[GenomeSegment],
    contigs: Sequence[str],
    *,
    chain_distance: int = DEFAULT_CHAIN_DISTANCE,
) -> Iterator[Locus]:
    """The loci of the VCF records of both files (doc/spec.md §5.3), contig by contig in
    the order of contigs and by start, yielded as the records are read.

    records_a, records_b: each file's records in the order parse_vcf yields them, by
        contig in the order of contigs and then by span start, records with no ALT among
        them.
    segments: the segmented genome, by contig in the order of contigs and then by start,
        with or without GENERIC segments, which are read and dropped. It is read only as
        far as the records need: up to the first segment after the span of the last record.
    contigs: the reference's contig names, in the order of its `.fai`.
    chain_distance: the most positions that may lie between a record and a locus for the
        record to chain to it, a whole number of 0 or more.

    Positions are 0-based and half-open. A locus holds its records with an ALT, the records
    with no ALT whose REF overlaps it, and the STR regions that lie in it.

    Raises ValueError, when the stream starts, for a chain_distance that is not a whole
    number of 0 or more and for a contig named twice in contigs; and when the stream
    reaches it, for a record or a segment on a contig that contigs lacks, a record that
    comes before the one read before it in its file, a segment that comes before the one
    read before it, and an STR region that overlaps the STR region read before it.
    """
    _check_chain_distance(chain_distance)
    ranks = {contig: rank for rank, contig in enumerate(contigs)}
    _check_contigs(contigs, ranks)
    ranked_records = heapq.merge(
        _make_ranked_records(records_a, _SourceFile.A, ranks),
        _make_ranked_records(records_b, _SourceFile.B, ranks),
        key=lambda ranked_record: ranked_record.start_position,
    )
    ranked_segments = _make_ranked_segments(segments, ranks)
    last_segment_position: _Position | None = None
    held_regions: list[_RankedSegment] = []
    # The records with no ALT that the reading has not passed, by end, and those it has
    # passed that the open locus or an STR region still holds.
    no_alt_ahead: list[tuple[_Position, _SourceFile, int, _RankedRecord]] = []
    no_alt_passed: list[_RankedRecord] = []
    locus: _OpenLocus | None = None

    ranked_record = next(ranked_records, None)
    while ranked_record is not None:
        last_segment_position = _read_segments_up_to(
            ranked_record.end_position, ranked_segments, last_segment_position, held_regions
        )
        # The reach is found anew for each record, after its segments are read: an STR
        # region that starts within the chain distance of the locus may be read only now.
        is_yielded = False
        if locus is not None and not _is_taken(ranked_record, locus, held_regions, chain_distance):
            yield _make_locus(locus, held_regions, _get_no_alt(no_alt_ahead, no_alt_passed))
            locus = None
            is_yielded = True
        if _is_with_alt(ranked_record.vcf_record):
            locus = _make_widened_locus(locus, ranked_record, held_regions)
        else:
            entry = (
                ranked_record.end_position,
                ranked_record.source_file,
                ranked_record.index_in_file,
                ranked_record,
            )
            heapq.heappush(no_alt_ahead, entry)

        ranked_record = next(ranked_records, None)
        next_start = None if ranked_record is None else ranked_record.start_position
        region_count = len(held_regions)
        held_regions = [
            ranked_region
            for ranked_region in held_regions
            if _is_region_held(ranked_region, locus, next_start)
        ]
        # A record with no ALT the reading has passed is held only by the open locus or by
        # an STR region, so it is checked again only when one of them goes: checking every
        # one at each record took time in the square of the records over a satellite.
        if is_yielded or len(held_regions) < region_count:
            no_alt_passed = [
                no_alt
                for no_alt in no_alt_passed
                if _is_no_alt_held(no_alt, locus, next_start, held_regions)
            ]
        while no_alt_ahead and (next_start is None or no_alt_ahead[0][0] <= next_start):
            no_alt = heapq.heappop(no_alt_ahead)[-1]
            if _is_no_alt_held(no_alt, locus, next_start, held_regions):
                no_alt_passed.append(no_alt)
    if locus is not None:
        yield _make_locus(locus, held_regions, no_alt_passed)


def _get_no_alt(
    no_alt_ahead: Iterable[tuple[_Position, _SourceFile, int, _RankedRecord]],
    no_alt_passed: Iterable[_RankedRecord],
) -> list[_RankedRecord]:
    """The records with no ALT held, ahead of the reading and passed by it."""
    return [*no_alt_passed, *(entry[-1] for entry in no_alt_ahead)]


def _check_chain_distance(chain_distance: int) -> None:
    is_whole = isinstance(chain_distance, int) and not isinstance(chain_distance, bool)
    if not (is_whole and chain_distance >= 0):
        raise ValueError(
            f"chain_distance must be a whole number of 0 or more, not {chain_distance!r}"
        )


def _check_contigs(contigs: Sequence[str], ranks: Mapping[str, int]) -> None:
    if len(ranks) < len(contigs):
        repeated = next(contig for rank, contig in enumerate(contigs) if ranks[contig] != rank)
        raise ValueError(f"contigs names {repeated} twice; each contig must be named once")


def _is_with_alt(vcf_record: VcfRecord) -> bool:
    return len(vcf_record.alleles) > 1


def _make_ranked_records(
    vcf_records: Iterable[VcfRecord], source_file: _SourceFile, ranks: Mapping[str, int]
) -> Iterator[_RankedRecord]:
    """One file's records with the rank of their contigs, checked to be on a contig of
    ranks and in order of span start."""
    file_name = source_file.name
    previous: _RankedRecord | None = None
    for index_in_file, vcf_record in enumerate(vcf_records):
        span_start = vcf_record.span[0]
        if vcf_record.contig not in ranks:
            raise ValueError(
                f"file {file_name}: a VCF record on {vcf_record.contig}, which contigs lacks,"
                f" at span start {span_start} (0-based)"
            )
        ranked_record = _RankedRecord(
            ranks[vcf_record.contig], source_file, index_in_file, vcf_record
        )
        if previous is not None and ranked_record.start_position < previous.start_position:
            raise ValueError(
                f"file {file_name}: a VCF record on {vcf_record.contig} at span start"
                f" {span_start} comes after one on {previous.vcf_record.contig} at span start"
                f" {previous.vcf_record.span[0]} (0-based); build_loci needs each file's"
                " records by contig, in the order of contigs, and then by span start"
            )
        yield ranked_record
        previous = ranked_record


def _make_ranked_segments(
    segments: Iterable[GenomeSegment], ranks: Mapping[str, int]
) -> Iterator[_RankedSegment]:
    """The segments with the rank of their contigs. Each is checked to be on a contig of
    ranks and to start no earlier than the one before it, and each STR region not to
    overlap the STR region before it."""
    previous_segment: _RankedSegment | None = None
    previous_region: _RankedSegment | None = None
    for segment in segments:
        if segment.contig not in ranks:
            raise ValueError(
                f"segments: a segment on {segment.contig}, which contigs lacks, at start"
                f" {segment.start} (0-based)"
            )
        ranked_segment = (ranks[segment.contig], segment)
        # A GENERIC segment out of order would stop the reading of segments before an STR
        # region a record touches, and the locus would come out too narrow.
        _check_segment_order(ranked_segment, previous_segment)
        previous_segment = ranked_segment
        if segment.region_type is not RegionType.GENERIC:
            _check_region_overlap(ranked_segment, previous_region)
            previous_region = ranked_segment
        yield ranked_segment


def _check_segment_order(ranked_segment: _RankedSegment, previous: _RankedSegment | None) -> None:
    if previous is None:
        return
    rank, segment = ranked_segment
    previous_rank, previous_segment = previous
    if (rank, segment.start) < (previous_rank, previous_segment.start):
        raise ValueError(
            f"segments: a segment on {segment.contig} at start {segment.start} comes after"
            f" one on {previous_segment.contig} at start {previous_segment.start} (0-based);"
            " build_loci needs the segments by contig, in the order of contigs, and then by"
            " start"
        )


def _check_region_overlap(ranked_region: _RankedSegment, previous: _RankedSegment | None) -> None:
    if previous is None:
        return
    rank, region = ranked_region
    previous_rank, previous_region = previous
    if rank == previous_rank and region.start < previous_region.end:
        raise ValueError(
            f"segments: the STR region {region.contig} {region.start}-{region.end} overlaps"
            f" the one before it, {previous_region.start}-{previous_region.end} (0-based,"
            " half-open); STR regions do not overlap"
        )


def _read_segments_up_to(
    position: _Position,
    ranked_segments: Iterator[_RankedSegment],
    last_segment_position: _Position | None,
    held_regions: list[_RankedSegment],
) -> _Position | None:
    """Reads the segments up to the first whose 0-based start comes after position, adds
    the STR regions among them to held_regions, and returns the rank and the start of the
    last segment read. position is a rank and the end of a span, one past its last
    position."""
    while last_segment_position is None or last_segment_position <= position:
        ranked_segment = next(ranked_segments, None)
        if ranked_segment is None:
            return last_segment_position
        rank, segment = ranked_segment
        last_segment_position = (rank, segment.start)
        if segment.region_type is not RegionType.GENERIC:
            held_regions.append(ranked_segment)
    return last_segment_position


def _find_reach(
    locus: _OpenLocus, held_regions: Iterable[_RankedSegment], chain_distance: int
) -> int:
    """The reach of the open locus: the last 0-based span start at which a later record is
    still taken by it, an inclusive bound. It is the locus's end, one past its last
    position, plus the chain distance, or the end of a held STR region on its contig that
    starts no further, whichever is further.

    A record that starts past end + chain_distance and no further than such a region's end
    touches the region, which joins the record's locus to this one (§5.3 step 4).
    """
    chained_end = locus.end + chain_distance
    reach = chained_end
    for rank, region in held_regions:
        if rank == locus.rank and region.start <= chained_end:
            reach = max(reach, region.end)
    return reach


def _is_taken(
    ranked_record: _RankedRecord,
    locus: _OpenLocus,
    held_regions: Iterable[_RankedSegment],
    chain_distance: int,
) -> bool:
    """Whether a record, of either kind, is taken by the open locus: on its contig, and
    starting no further than its reach."""
    if ranked_record.rank != locus.rank:
        return False
    span_start = ranked_record.vcf_record.span[0]
    return span_start <= _find_reach(locus, held_regions, chain_distance)


def _is_touching(vcf_record: VcfRecord, region: GenomeSegment) -> bool:
    """Whether a record's span overlaps an STR region or lies next to it (§5.3), both on
    one contig, 0-based and half-open: the span ends on the region's first position or
    starts on the position after its last."""
    span_start, span_end = vcf_record.span
    return span_start <= region.end and region.start <= span_end


def _is_overlapping(start: int, end: int, other_start: int, other_end: int) -> bool:
    """Whether [start, end) and [other_start, other_end), 0-based and half-open on one
    contig, share a position."""
    return start < other_end and other_start < end


def _make_widened_locus(
    locus: _OpenLocus | None, ranked_record: _RankedRecord, held_regions: Iterable[_RankedSegment]
) -> _OpenLocus:
    """The open locus with a record with an ALT added, widened over the record's span and
    over the STR regions it touches; the locus given is changed and returned, and a new
    one is made when none is open."""
    vcf_record = ranked_record.vcf_record
    span_start, span_end = vcf_record.span
    if locus is None:
        locus = _OpenLocus(ranked_record.rank, vcf_record.contig, span_start, span_end, [])
    locus.records.append(ranked_record)
    locus.start = min(locus.start, span_start)
    locus.end = max(locus.end, span_end)
    for rank, region in held_regions:
        if rank == locus.rank and _is_touching(vcf_record, region):
            locus.start = min(locus.start, region.start)
            locus.end = max(locus.end, region.end)
    return locus


def _is_region_held(
    ranked_region: _RankedSegment, locus: _OpenLocus | None, next_start: _Position | None
) -> bool:
    """Whether an STR region may still be in a locus: in a gap of the open locus, or
    touched by a later record."""
    rank, region = ranked_region
    if locus is not None and rank == locus.rank and region.end > locus.start:
        return True
    return next_start is not None and (rank, region.end) >= next_start


def _is_no_alt_held(
    no_alt: _RankedRecord,
    locus: _OpenLocus | None,
    next_start: _Position | None,
    held_regions: Iterable[_RankedSegment],
) -> bool:
    """Whether a record with no ALT may still be in a locus: one that overlaps the open
    locus, a later one, or a held STR region that may still widen a later locus over it."""
    span_start, span_end = no_alt.vcf_record.span
    if locus is not None and no_alt.rank == locus.rank and span_end > locus.start:
        return True
    if next_start is not None and no_alt.end_position > next_start:
        return True
    return any(
        rank == no_alt.rank and _is_overlapping(span_start, span_end, region.start, region.end)
        for rank, region in held_regions
    )


def _make_locus(
    locus: _OpenLocus,
    held_regions: Iterable[_RankedSegment],
    held_no_alt: Iterable[_RankedRecord],
) -> Locus:
    """The Locus of the open locus, with the held STR regions and records with no ALT that
    overlap it; each file's records in the order of its stream."""
    overlapping_no_alt = [
        no_alt
        for no_alt in held_no_alt
        if no_alt.rank == locus.rank
        and _is_overlapping(*no_alt.vcf_record.span, locus.start, locus.end)
    ]
    ranked_records = sorted(
        [*locus.records, *overlapping_no_alt],
        key=lambda ranked_record: (ranked_record.source_file, ranked_record.index_in_file),
    )
    return Locus(
        locus.contig,
        locus.start,
        locus.end,
        tuple(
            ranked.vcf_record for ranked in ranked_records if ranked.source_file is _SourceFile.A
        ),
        tuple(
            ranked.vcf_record for ranked in ranked_records if ranked.source_file is _SourceFile.B
        ),
        tuple(
            region
            for rank, region in held_regions
            if rank == locus.rank
            and _is_overlapping(region.start, region.end, locus.start, locus.end)
        ),
    )
