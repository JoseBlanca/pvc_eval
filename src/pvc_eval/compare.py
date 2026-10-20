"""compare(), the stages of doc/spec.md §4 joined, from two VCF files to the rows of the
table (doc/specs/compare.md, "What compare() does").

Positions are 0-based and half-open, as everywhere inside pvc_eval.

This version runs every stage in the caller's thread, whatever finder_threads is: the
readers are read as build_loci asks for their records, and the STR finder scans a contig
when build_loci first needs its segments.

Public: compare, which yields the rows; CompareConfig, its options that change the rows;
PVC_STRATA, the value of CompareConfig.strata for --strata pvc; DEFAULT_FINDER_THREADS.
"""

import logging
import time
from collections.abc import Callable, Container, Generator, Iterator, Sequence
from contextlib import ExitStack, closing
from dataclasses import dataclass
from itertools import groupby
from operator import attrgetter
from pathlib import Path
from typing import cast

from pvc_eval.evaluate import make_outcome_rows
from pvc_eval.fasta import IndexedFasta
from pvc_eval.kind import locus_kind
from pvc_eval.loci import DEFAULT_CHAIN_DISTANCE, build_loci
from pvc_eval.locus import Locus
from pvc_eval.outcome import OutcomeRow
from pvc_eval.reconcile import PHASING_CAP, reconcile_loci
from pvc_eval.regions import Regions, read_regions
from pvc_eval.samples import make_sample_pairs
from pvc_eval.segments import GenomeSegment, RegionType, read_str_regions, segment_contig
from pvc_eval.strata import pvc_stratum, read_strata
from pvc_eval.vcf import DEFAULT_PLOIDY, DEFAULT_SPAN_LIMIT, VcfRecord, parse_vcf

_logger = logging.getLogger(__name__)

PVC_STRATA = sentinel("PVC_STRATA")
"""The value of CompareConfig.strata for --strata pvc, apart from a BED whose path is pvc."""

DEFAULT_FINDER_THREADS = 2  # doc/specs/compare.md, "Which thread does what"


@dataclass(frozen=True, slots=True, kw_only=True)
class CompareConfig:
    """The options of doc/spec.md §3 that change the rows, with the defaults of the stages.

    regions: the --regions BED; a locus that shares no position with it is skipped.
    strata: the --strata BED, PVC_STRATA for --strata pvc, or None for no stratum.
    str_regions: the --str-regions BED, read instead of running the STR finder.
    sample_map: the --sample-map file of sample pairs; without it, samples pair by name.
    ploidy_a, ploidy_b: the haplotypes of each sample of file A and of file B, 1 or 2.
    phasing_cap: the most possible genotypes a sample may have, 1 or more.
    span_limit: the most bases a span may reach beyond its REF, 1 or more, or None for no
        limit.
    chain_distance: the most positions between a record and a locus for it to chain to
        the locus, 0 or more.

    The stages check the numbers before compare() gives its first row.
    """

    regions: Path | None = None
    strata: Path | PVC_STRATA | None = None
    str_regions: Path | None = None
    sample_map: Path | None = None
    ploidy_a: int = DEFAULT_PLOIDY
    ploidy_b: int = DEFAULT_PLOIDY
    phasing_cap: int = PHASING_CAP
    span_limit: int | None = DEFAULT_SPAN_LIMIT
    chain_distance: int = DEFAULT_CHAIN_DISTANCE


def compare(
    a: Path,
    b: Path,
    reference: Path,
    config: CompareConfig = CompareConfig(),  # noqa: B008 -- frozen, so one instance is safe
    *,
    finder_threads: int = DEFAULT_FINDER_THREADS,
) -> Iterator[OutcomeRow]:
    """The rows of file A against file B over the reference, every row included: locus by
    locus, contig by contig in the order of the reference's .fai, and in each locus in the
    order of the sample pairs. They are the rows of the stages of doc/spec.md §4 run one
    after another (doc/specs/compare.md, "The stages, joined").

    a, b: file A and file B, plain or bgzipped VCF.
    reference: the FASTA, with its .fai beside it.
    finder_threads: how many threads the STR finder may use, 0 for none. This version
        uses none whatever the number.

    Without config.str_regions, the STR finder scans a contig only once a record on it
    has reached the loci and, with config.regions, the regions name it. compare() logs at
    INFO "chr1: finding STR regions in 248,956,422 bases." when a scan starts, and
    "chr20: 31,204 loci compared, 1,880 skipped by --regions; 4 min 12 s since the
    start." when the loci of a contig are done, the clause about --regions only with
    config.regions.

    Raises, when the stream starts, ValueError for a finder_threads that is not a whole
    number of 0 or more, and, before the first row, the ValueError a stage raises for a
    setting of config out of range, naming its field; and PvcEvalInputError for the faults
    of the input the stages raise it for, when the stream reaches them.
    """
    _check_finder_threads(finder_threads)
    start_time = time.monotonic()
    with IndexedFasta(reference) as indexed_fasta, ExitStack() as stack:
        contigs = indexed_fasta.contigs
        sample_pairs = make_sample_pairs(a, b, sample_map=config.sample_map)
        regions = None if config.regions is None else read_regions(config.regions, contigs)
        get_stratum = _make_get_stratum(config.strata, contigs)
        records_a = _parse_vcf_until_exit(stack, a, indexed_fasta, config, config.ploidy_a)
        records_b = _parse_vcf_until_exit(stack, b, indexed_fasta, config, config.ploidy_b)
        # The contigs of the records build_loci has taken, from which the segments are decided.
        taken_contigs: set[str] = set()
        if config.str_regions is None:
            segments = _make_segments(indexed_fasta, taken_contigs, regions)
        else:
            segments = read_str_regions(config.str_regions, indexed_fasta)
        loci = build_loci(
            _note_taken_contigs(records_a, taken_contigs),
            _note_taken_contigs(records_b, taken_contigs),
            segments,
            contigs,
            chain_distance=config.chain_distance,
        )
        reconciled_loci = reconcile_loci(
            _drop_skipped_loci(loci, regions, start_time),
            indexed_fasta,
            sample_pairs,
            ploidy_a=config.ploidy_a,
            ploidy_b=config.ploidy_b,
            phasing_cap=config.phasing_cap,
        )
        for reconciled in reconciled_loci:
            locus = reconciled.locus
            yield from make_outcome_rows(
                reconciled, sample_pairs, kind=locus_kind(locus), stratum=get_stratum(locus)
            )


def _check_finder_threads(finder_threads: int) -> None:
    if not (_is_whole(finder_threads) and finder_threads >= 0):
        raise ValueError(
            f"finder_threads must be a whole number of 0 or more, not {finder_threads!r}"
        )


def _is_whole(number: object) -> bool:
    return isinstance(number, int) and not isinstance(number, bool)


def _parse_vcf_until_exit(
    stack: ExitStack, path: Path, reference: IndexedFasta, config: CompareConfig, ploidy: int
) -> Iterator[VcfRecord]:
    """parse_vcf's records of a file, its stream closed when stack exits, so that a run
    stopped early closes the file and parse_vcf's pool at once."""
    # parse_vcf is a generator function, so its stream has close().
    vcf_records = cast(
        Generator[VcfRecord],
        parse_vcf(path, reference, span_limit=config.span_limit, ploidy=ploidy),
    )
    return stack.enter_context(closing(vcf_records))


def _make_get_stratum(
    strata: Path | PVC_STRATA | None, contigs: Sequence[str]
) -> Callable[[Locus], str | None]:
    """The function that gives a locus its stratum: pvc_stratum for --strata pvc, the
    classes of a --strata BED, read whole now, and None without --strata."""
    if strata is PVC_STRATA:
        return pvc_stratum
    if strata is None:
        return _get_no_stratum
    return read_strata(strata, contigs).stratum


def _get_no_stratum(locus: Locus) -> None:  # noqa: ARG001 -- the signature of a stratum
    return None


def _note_taken_contigs(
    vcf_records: Iterator[VcfRecord], taken_contigs: set[str]
) -> Iterator[VcfRecord]:
    """The records, each one's contig added to taken_contigs as it is handed over, so that
    a record counts as taken when build_loci takes it (doc/specs/compare.md, "How the
    segments reach build_loci")."""
    for vcf_record in vcf_records:
        taken_contigs.add(vcf_record.contig)
        yield vcf_record


def _make_segments(
    reference: IndexedFasta, taken_contigs: Container[str], regions: Regions | None
) -> Iterator[GenomeSegment]:
    """The segments build_loci reads, contig by contig in the order of the .fai: the STR
    regions of each contig that needs them, and for each other contig a GENERIC segment at
    its first base, 0 to 1, which is not part of the segmented genome
    (doc/specs/compare.md, "How the segments reach build_loci").

    A contig that does not need its STR regions when it is reached is decided again when
    build_loci asks for the segment after its GENERIC one. build_loci asks that only for a
    record on that contig or a later one, which has been taken by then; so if the contig
    still does not need them, it has no record, or the regions do not name it.
    """
    for contig in reference.contigs:
        if not _is_needed(contig, reference, taken_contigs, regions):
            # It tells build_loci that no segment before position 0 of contig is to come.
            yield GenomeSegment(contig, 0, 1, RegionType.GENERIC, None)
            if not _is_needed(contig, reference, taken_contigs, regions):
                continue
        yield from _find_str_regions(reference, contig)


def _is_needed(
    contig: str, reference: IndexedFasta, taken_contigs: Container[str], regions: Regions | None
) -> bool:
    """Whether a contig needs its STR regions: a record on it has been taken, and regions,
    when given, name it."""
    if contig not in taken_contigs:
        return False
    return regions is None or regions.overlaps(contig, 0, reference.get_length(contig))


def _find_str_regions(reference: IndexedFasta, contig: str) -> Iterator[GenomeSegment]:
    """The STR regions of a contig, from the STR finder run on it when the first is asked
    for; build_loci drops the GENERIC segments, so they are not handed to it."""
    length = reference.get_length(contig)
    _logger.info("%s: finding STR regions in %s bases.", contig, f"{length:,}")
    for segment in segment_contig(reference, contig):
        if segment.region_type is not RegionType.GENERIC:
            yield segment


def _drop_skipped_loci(
    loci: Iterator[Locus], regions: Regions | None, start_time: float
) -> Iterator[Locus]:
    """The loci select_loci keeps, with the counts of each contig logged when its loci are
    done: when the first locus of a later contig is built, or when the loci end.

    start_time: the time.monotonic() of the start of the run.
    """
    for contig, contig_loci in groupby(loci, key=attrgetter("contig")):
        num_built = 0
        num_compared = 0
        for locus in contig_loci:
            num_built += 1
            # select_loci's rule (doc/specs/loci.md, "Regions"), here so that it is counted.
            if regions is None or regions.overlaps(locus.contig, locus.start, locus.end):
                num_compared += 1
                yield locus
        _log_contig_done(contig, num_compared, num_built - num_compared, regions, start_time)


def _log_contig_done(
    contig: str, num_compared: int, num_skipped: int, regions: Regions | None, start_time: float
) -> None:
    noun = "locus" if num_compared == 1 else "loci"
    skipped = "" if regions is None else f", {num_skipped:,} skipped by --regions"
    duration = _format_duration(time.monotonic() - start_time)
    _logger.info(
        "%s: %s %s compared%s; %s since the start.",
        contig,
        f"{num_compared:,}",
        noun,
        skipped,
        duration,
    )


def _format_duration(seconds: float) -> str:
    """A duration in whole seconds, rounded down: "12 s", "4 min 12 s", "1 h 2 min 5 s"."""
    minutes, whole_seconds = divmod(int(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours} h {minutes} min {whole_seconds} s"
    if minutes:
        return f"{minutes} min {whole_seconds} s"
    return f"{whole_seconds} s"
