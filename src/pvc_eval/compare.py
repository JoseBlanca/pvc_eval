"""compare(), the stages of doc/spec.md §4 joined, from two VCF files to the rows of the
table (doc/specs/compare.md, "What compare() does").

Positions are 0-based and half-open, as everywhere inside pvc_eval.

With finder_threads 0 every stage runs in the caller's thread: the readers are read as
build_loci asks for their records, and the STR finder scans a contig when build_loci first
needs its segments. With 1 or more, each VCF is read in a thread of its own into a bounded
buffer, and the STR finder runs on a pool of finder_threads threads, a few contigs ahead of
the loci (doc/specs/compare.md, "Which thread does what"). The rows are the same.

Public: compare, which yields the rows; CompareConfig, its options that change the rows;
PVC_STRATA, the value of CompareConfig.strata for --strata pvc; DEFAULT_FINDER_THREADS.
"""

import functools
import logging
import sys
import threading
import time
from collections.abc import Callable, Container, Generator, Iterable, Iterator, Sequence
from contextlib import ExitStack, closing, suppress
from dataclasses import dataclass
from itertools import groupby
from operator import attrgetter
from pathlib import Path
from queue import Queue, ShutDown
from typing import cast

from pvc_eval._vcf_index import parse_index_contigs
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

# doc/specs/compare.md, "Which thread does what": a reader hands its records over in batches
# of _BATCH_RECORDS, so that its buffer is locked once per batch and not per record, and its
# buffer holds at most 10,000 records, so _BUFFER_BATCHES batches, 9,984 records.
_BATCH_RECORDS = 64
_BUFFER_BATCHES = 10_000 // _BATCH_RECORDS

_END = sentinel("_END")
"""What a reader puts in its buffer after its last record, when its reading ends with no
error."""

type _Buffer = Queue[list[VcfRecord] | BaseException | _END]
"""A reader's buffer: its batches of records, and then how its reading ended, _END or the
exception that ended it."""


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
) -> Generator[OutcomeRow]:
    """The rows of file A against file B over the reference, every row included: locus by
    locus, contig by contig in the order of the reference's .fai, and in each locus in the
    order of the sample pairs. They are the rows of the stages of doc/spec.md §4 run one
    after another (doc/specs/compare.md, "The stages, joined").

    a, b: file A and file B, plain or bgzipped VCF.
    reference: the FASTA, with its .fai beside it.
    finder_threads: how many threads the STR finder uses, a whole number of 0 or more. With
        0, compare() makes no thread: every stage runs in the caller's thread. With 1 or
        more, each VCF is also read in a thread of its own, and the indexes beside the VCFs
        plan the contigs to scan. Closing the stream stops the threads; it waits for the
        readers, and not for a scan in progress.

    With config.regions, the VCF records of a contig the regions do not name are not read,
    and their lines are checked only for their place (doc/specs/regions_read.md).

    Without config.str_regions, the STR regions of a contig reach the loci only once a
    record on it has reached them and, with config.regions, the regions name it. A thread
    of the finder may scan a contig earlier, or one that turns out not to be needed, and
    drop its STR regions. compare() logs at
    INFO "chr1: finding STR regions in 248,956,422 bases." when a scan starts, and
    "chr20: 31,204 loci compared, 1,880 skipped by --regions; 4 min 12 s since the
    start." when the loci of a contig are done, the clause about --regions only with
    config.regions.

    Raises, when the stream starts, ValueError for a finder_threads that is not a whole
    number of 0 or more; it checks no field of config. The stages check those: before the
    first row, the stage that takes a field out of range raises ValueError naming it. Raises
    PvcEvalInputError for the faults of the input the stages raise it for, when the stream
    reaches them, with or without threads after the same rows.
    """
    _check_finder_threads(finder_threads)
    start_time = time.monotonic()
    with IndexedFasta(reference) as indexed_fasta, ExitStack() as stack:
        contigs = indexed_fasta.contigs
        sample_pairs = make_sample_pairs(a, b, sample_map=config.sample_map)
        regions = None if config.regions is None else read_regions(config.regions, contigs)
        get_stratum = _make_get_stratum(config.strata, contigs)
        contigs_read = _find_contigs_read(indexed_fasta, regions)
        if finder_threads == 0:
            records_a = _parse_vcf_until_exit(
                stack, a, indexed_fasta, config, config.ploidy_a, contigs_read
            )
            records_b = _parse_vcf_until_exit(
                stack, b, indexed_fasta, config, config.ploidy_b, contigs_read
            )
            finder = None
        else:
            records_a, records_b, finder = _start_threads(
                stack,
                (a, b),
                reference,
                indexed_fasta,
                config,
                regions,
                contigs_read,
                finder_threads,
            )
        # The contigs of the records build_loci has taken, from which the segments are decided.
        taken_contigs: set[str] = set()
        if config.str_regions is None:
            segments = _make_segments(indexed_fasta, taken_contigs, regions, finder)
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
            # It also logs the counts of each contig when its loci are done.
            _drop_skipped_loci(loci, regions, start_time),
            indexed_fasta,
            sample_pairs,
            ploidy_a=config.ploidy_a,
            ploidy_b=config.ploidy_b,
            phasing_cap=config.phasing_cap,
        )
        for reconciled in reconciled_loci:
            locus = reconciled.locus
            partly_outside = regions is not None and not regions.covers(
                locus.contig, locus.start, locus.end
            )
            yield from make_outcome_rows(
                reconciled,
                sample_pairs,
                kind=locus_kind(locus),
                stratum=get_stratum(locus),
                partly_outside=partly_outside,
            )


def _check_finder_threads(finder_threads: int) -> None:
    if not (_is_whole(finder_threads) and finder_threads >= 0):
        raise ValueError(
            f"finder_threads must be a whole number of 0 or more, not {finder_threads!r}"
        )


def _is_whole(number: object) -> bool:
    return isinstance(number, int) and not isinstance(number, bool)


def _find_contigs_read(reference: IndexedFasta, regions: Regions | None) -> frozenset[str] | None:
    """The contigs whose VCF records are read: those regions names when given, or None for
    every contig (doc/specs/regions_read.md)."""
    if regions is None:
        return None
    return frozenset(
        contig for contig in reference.contigs if _is_named(contig, reference, regions)
    )


def _parse_vcf_until_exit(
    stack: ExitStack,
    path: Path,
    reference: IndexedFasta,
    config: CompareConfig,
    ploidy: int,
    contigs_read: frozenset[str] | None,
) -> Iterator[VcfRecord]:
    """parse_vcf's records of a file on contigs_read, or on every contig when it is None,
    its stream closed when stack exits, so that a run stopped early closes the file and
    parse_vcf's pool at once."""
    # parse_vcf is a generator function, so its stream has close().
    vcf_records = cast(
        Generator[VcfRecord],
        parse_vcf(
            path, reference, span_limit=config.span_limit, ploidy=ploidy, contigs=contigs_read
        ),
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
    reference: IndexedFasta,
    taken_contigs: Container[str],
    regions: Regions | None,
    finder: _Finder | None,
) -> Iterator[GenomeSegment]:
    """The segments build_loci reads, contig by contig in the order of the .fai: the STR
    regions of each contig that needs them, and for each other contig a GENERIC segment at
    its first base, 0 to 1, which is not part of the segmented genome
    (doc/specs/compare.md, "How the segments reach build_loci").

    A contig that does not need its STR regions when it is reached is decided again when
    build_loci asks for the segment after its GENERIC one. build_loci asks that only for a
    record on that contig or a later one, which has been taken by then; so if the contig
    still does not need them, it has no record, or the regions do not name it.

    finder: the pool that scans the contigs, or None to scan each in this thread.
    """
    for rank, contig in enumerate(reference.contigs):
        if finder is not None:
            finder.drop_passed(rank)
        if not _is_needed(contig, reference, taken_contigs, regions):
            # It tells build_loci that no segment before position 0 of contig is to come.
            yield GenomeSegment(contig, 0, 1, RegionType.GENERIC, None)
            if not _is_needed(contig, reference, taken_contigs, regions):
                continue
        if finder is None:
            yield from _find_str_regions(reference, contig)
        else:
            yield from finder.take_str_regions(contig)


def _is_needed(
    contig: str, reference: IndexedFasta, taken_contigs: Container[str], regions: Regions | None
) -> bool:
    """Whether a contig needs its STR regions: a record on it has been taken, and regions,
    when given, name it."""
    return contig in taken_contigs and _is_named(contig, reference, regions)


def _is_named(contig: str, reference: IndexedFasta, regions: Regions | None) -> bool:
    """Whether regions, when given, name a contig."""
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


def _start_threads(
    stack: ExitStack,
    paths: tuple[Path, Path],
    reference_path: Path,
    reference: IndexedFasta,
    config: CompareConfig,
    regions: Regions | None,
    contigs_read: frozenset[str] | None,
    finder_threads: int,
) -> tuple[Iterator[VcfRecord], Iterator[VcfRecord], _Finder | None]:
    """Starts a reader for each of file A and file B and, without config.str_regions, the
    finder's pool with the contigs the indexes plan in line; returns the records of each
    reader, as the loci side takes them, and the finder (doc/specs/compare.md, "Which thread
    does what").

    reference_path: the FASTA that reference reads, which each reader and each thread of
        the finder opens for itself.
    contigs_read: the contigs whose records the readers read, or None for every contig.

    When stack exits, the stop is set; no thread is waited for (doc/specs/compare.md, "The
    stop").
    """
    stop = threading.Event()
    finder = None
    if config.str_regions is None:
        finder = _Finder(reference, regions, finder_threads, stop)
    buffers: tuple[_Buffer, _Buffer] = (
        Queue(maxsize=_BUFFER_BATCHES),
        Queue(maxsize=_BUFFER_BATCHES),
    )
    stack.callback(_stop_threads, stop, finder, buffers)
    if finder is not None:
        # The reference lacks a contig an index lists: a record on it stops the run when a
        # reader reaches it, so its scan is not planned.
        index_contigs = parse_index_contigs(paths[0]) | parse_index_contigs(paths[1])
        finder.put(contig for contig in reference.contigs if contig in index_contigs)
        finder.start(reference_path)
    for path, ploidy, buffer, file_name in zip(
        paths, (config.ploidy_a, config.ploidy_b), buffers, "AB", strict=True
    ):
        open_vcf_records = functools.partial(
            _open_vcf_records,
            path,
            span_limit=config.span_limit,
            ploidy=ploidy,
            contigs=contigs_read,
        )
        threading.Thread(
            target=_parse_vcf_into_buffer,
            args=(open_vcf_records, reference_path, buffer, finder),
            name=f"pvc_eval reader {file_name}",
            daemon=True,
        ).start()
    return _take_records(buffers[0]), _take_records(buffers[1]), finder


def _stop_threads(
    stop: threading.Event, finder: _Finder | None, buffers: Iterable[_Buffer]
) -> None:
    """Sets the stop, wakes the finder's threads that wait for a contig, and shuts the
    buffers down, which wakes a reader waiting for room; each reader then ends within one
    batch of parse_vcf's work. No thread is waited for: a stop run by the garbage collector,
    in a thread that may hold threading's lock on its list of threads, would wait forever
    for a thread that ends, which takes that lock (doc/specs/compare.md, "The stop").

    While Python is exiting it does nothing: shutting a buffer down takes a lock, which
    Python then refuses, and every thread of compare() is a daemon thread, which ends with
    the process."""
    if sys.is_finalizing():
        return
    stop.set()
    if finder is not None:
        finder.stop()
    for buffer in buffers:
        buffer.shutdown(immediate=True)


def _open_vcf_records(
    path: Path,
    reference: IndexedFasta,
    *,
    span_limit: int | None,
    ploidy: int,
    contigs: frozenset[str] | None,
) -> Generator[VcfRecord]:
    """parse_vcf's stream over reference, which a reader thread has opened for itself."""
    # parse_vcf is a generator function, so its stream has close().
    return cast(
        Generator[VcfRecord],
        parse_vcf(path, reference, span_limit=span_limit, ploidy=ploidy, contigs=contigs),
    )


def _parse_vcf_into_buffer(
    open_vcf_records: Callable[[IndexedFasta], Generator[VcfRecord]],
    reference_path: Path,
    buffer: _Buffer,
    finder: _Finder | None,
) -> None:
    """A reader thread: parse_vcf's records put in the buffer in batches of _BATCH_RECORDS,
    the last perhaps shorter, followed by how the reading ended, _END or the exception that
    ended it, which the loci side raises where it comes. Each contig the reader reaches is
    put in the finder's line. At the stop the buffer is shut down, and the reading ends with
    its stream closed.

    open_vcf_records: parse_vcf with its arguments but the reference, called in this thread
        with the reader's own handle on the reference.
    reference_path: the FASTA, which the reader opens for itself and closes when its reading
        ends, so that compare() closes its own without waiting for it, and nothing is opened
        for a reader that does not start.
    """
    # After the stop, nothing more is taken from the buffer.
    with suppress(ShutDown):
        _put_vcf_records(open_vcf_records, reference_path, buffer, finder)


def _put_vcf_records(
    open_vcf_records: Callable[[IndexedFasta], Generator[VcfRecord]],
    reference_path: Path,
    buffer: _Buffer,
    finder: _Finder | None,
) -> None:
    """What _parse_vcf_into_buffer does, raising ShutDown at the stop."""
    batch: list[VcfRecord] = []
    ending: BaseException | _END = _END
    try:
        with (
            IndexedFasta(reference_path) as reference,
            closing(open_vcf_records(reference)) as vcf_records,
        ):
            contig = None
            for vcf_record in vcf_records:
                if finder is not None and vcf_record.contig != contig:
                    finder.put([vcf_record.contig])
                contig = vcf_record.contig
                batch.append(vcf_record)
                if len(batch) == _BATCH_RECORDS:
                    buffer.put(batch)
                    batch = []
    except ShutDown:
        raise
    except BaseException as error:  # noqa: BLE001 -- raised by the loci side, or dropped at the stop
        ending = error
    if batch:
        buffer.put(batch)
    buffer.put(ending)


def _take_records(buffer: _Buffer) -> Iterator[VcfRecord]:
    """The records of a reader's buffer, one at a time, and then the exception that ended
    its reading, raised after them."""
    while True:
        item = buffer.get()
        if item is _END:
            return
        if isinstance(item, BaseException):
            raise item
        yield from item


class _Finder:
    """The STR finder's pool: threads that scan, each one contig at a time, the contigs in
    line, a few contigs ahead of the loci side (doc/specs/compare.md, "Which contigs are
    scanned", "How far ahead").

    A contig is put in line when planned from an index or when a reader reaches it, once,
    and only when the regions, when given, name it. A thread takes the contig of lowest rank
    in line while fewer than num_threads + 1 contigs are ahead, from the start of their scan
    until the loci side takes their STR regions; the contig the loci side waits for, it
    takes whatever the number ahead. A contig the loci side passes leaves the line, and the
    STR regions of its scan are dropped.
    """

    def __init__(
        self,
        reference: IndexedFasta,
        regions: Regions | None,
        num_threads: int,
        stop: threading.Event,
    ) -> None:
        self._reference = reference
        self._regions = regions
        self._num_threads = num_threads
        self._stop = stop
        self._ranks = {contig: rank for rank, contig in enumerate(reference.contigs)}
        # What follows is changed only with self._condition held, by the readers, the
        # finder's threads and the loci side: the line is shared by design.
        self._condition = threading.Condition()
        self._line: set[str] = set()
        self._ever_in_line: set[str] = set()
        self._ahead: set[str] = set()
        self._results: dict[str, tuple[GenomeSegment, ...] | BaseException] = {}
        # The rank of the contig the loci side is on; it has passed those of lower rank.
        self._loci_rank = 0
        self._waited_contig: str | None = None

    def start(self, reference_path: Path) -> None:
        """Starts the finder's threads. Each opens the FASTA at reference_path for itself and
        closes it when it ends, so that a scan still running after compare() has closed its
        reference reads a descriptor of its own, and nothing is opened for a thread that does
        not start. They are daemon threads, so that a scan in progress does not keep Python
        from exiting.
        """
        for number in range(1, self._num_threads + 1):
            thread = threading.Thread(
                target=self._find_in_line,
                args=(reference_path,),
                name=f"pvc_eval finder {number}",
                daemon=True,
            )
            thread.start()

    def put(self, contigs: Iterable[str]) -> None:
        """Puts contigs in line: each that the regions, when given, name, has not been in
        line before, and the loci side has not passed."""
        with self._condition:
            for contig in contigs:
                if contig in self._ever_in_line or self._ranks[contig] < self._loci_rank:
                    continue
                if _is_named(contig, self._reference, self._regions):
                    self._line.add(contig)
                    self._ever_in_line.add(contig)
            self._condition.notify_all()

    def drop_passed(self, rank: int) -> None:
        """Notes that the loci side is on the contig of this rank, and drops from the line,
        from the contigs ahead and from the results those of lower rank."""
        with self._condition:
            self._loci_rank = rank
            for contig in [*self._line, *self._ahead]:
                if self._ranks[contig] < rank:
                    self._line.discard(contig)
                    self._ahead.discard(contig)
                    self._results.pop(contig, None)
            self._condition.notify_all()

    def take_str_regions(self, contig: str) -> tuple[GenomeSegment, ...]:
        """The STR regions of a contig, for the loci side, waiting for its scan, which starts
        as soon as a thread is free; raises the exception its scan raised.

        The contig has been put in line: a record on it has been taken, so its reader has
        reached it, and the regions, when given, name it.
        """
        with self._condition:
            self._waited_contig = contig
            self._condition.notify_all()
            while contig not in self._results:
                self._condition.wait()
            result = self._results.pop(contig)
            self._ahead.discard(contig)
            self._waited_contig = None
            self._condition.notify_all()
        if isinstance(result, BaseException):
            raise result
        return result

    def stop(self) -> None:
        """Wakes the threads waiting for a contig, once the stop is set, so that they end."""
        with self._condition:
            self._condition.notify_all()

    def _find_in_line(self, reference_path: Path) -> None:
        """A finder thread: the contigs taken from the line, scanned one after another in the
        FASTA at reference_path, which it opens for itself, until the stop. A scan's
        exception, or the one that opening the FASTA raised, is passed to the loci side for
        each contig it takes, or dropped after the stop."""
        try:
            reference = IndexedFasta(reference_path)
        except BaseException as error:  # noqa: BLE001 -- raised by the loci side, or dropped at the stop
            while (contig := self._take_contig()) is not None:
                self._put_result(contig, error)
            return
        with reference:
            while (contig := self._take_contig()) is not None:
                try:
                    result: tuple[GenomeSegment, ...] | BaseException = tuple(
                        _find_str_regions(reference, contig)
                    )
                except BaseException as error:  # noqa: BLE001 -- raised by the loci side, or dropped at the stop
                    result = error
                self._put_result(contig, result)

    def _take_contig(self) -> str | None:
        """The next contig to scan, taken from the line, waiting while there is none to
        start; None at the stop."""
        with self._condition:
            while not self._stop.is_set():
                contig = self._find_next_contig()
                if contig is not None:
                    self._line.remove(contig)
                    self._ahead.add(contig)
                    return contig
                self._condition.wait()
            return None

    def _find_next_contig(self) -> str | None:
        """The contig in line to start now, with self._condition held: the one the loci side
        waits for, or the lowest ranked while fewer than num_threads + 1 are ahead."""
        if self._waited_contig is not None and self._waited_contig in self._line:
            return self._waited_contig
        if self._line and len(self._ahead) < self._num_threads + 1:
            return min(self._line, key=self._ranks.__getitem__)
        return None

    def _put_result(self, contig: str, result: tuple[GenomeSegment, ...] | BaseException) -> None:
        """Keeps a scan's STR regions, or its exception, for the loci side; drops them when
        the loci side has passed the contig or the stop is set."""
        with self._condition:
            if self._stop.is_set() or contig not in self._ahead:
                return
            self._results[contig] = result
            self._condition.notify_all()
