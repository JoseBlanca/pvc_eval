"""The VCF records of one file: stage 2 of doc/spec.md §4, read with the checks of §5.2,
normalised, given their spans and put in order of span start. Normalising trims the bases
a record's alleles share and marks a REF that differs from the reference; a VCF record is
never moved, split or decomposed.

Positions are 0-based and half-open inside pvc_eval. The 1-based POS of a VCF line is
converted here, and nowhere else.

Public: parse_vcf, which gives the VCF records; parse_vcf_samples, which gives the sample
names; VcfRecord and VcfGenotype; Span, the type of VcfRecord.span; the sentinels a
genotype is compared with, MISSING, MALFORMED_GT and HAPLOID_GT; and the constants PLOIDY,
DEFAULT_SPAN_LIMIT and DEFAULT_NUM_THREADS.
"""

import functools
import gzip
import heapq
import logging
import math
import re
import sys
import threading
from collections.abc import Callable, Generator, Iterable, Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import BLOCK_BASES, FastaWindow, IndexedFasta

_logger = logging.getLogger(__name__)

PLOIDY = 2  # doc/spec.md §2: a global constant

MISSING = sentinel("MISSING")
"""A haplotype the VCF genotype says nothing about, written `.` (doc/spec.md §2)."""

MALFORMED_GT = sentinel("MALFORMED_GT")
"""A GT with a number of alleles other than the ploidy, or an allele index the VCF
record lacks: the sample is invalid at the locus, and the run goes on (§5.4)."""

HAPLOID_GT = sentinel("HAPLOID_GT")
"""A GT of one allele index, such as `1` on chrX of a male sample: the sample is invalid
at the locus with its own reason, `haploid GT` (§5.4, §6, Open 1)."""

_SPANNING_DELETION = "*"
_GZIP_MAGIC = b"\x1f\x8b"
_DECIMAL_NUMBER = re.compile(r"[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?")
_NUM_FIXED_COLUMNS = 9  # CHROM to FORMAT
# doc/spec.md §5.2: the span limit when the caller sets none. A span that would reach
# further than the limit beyond its REF, on either side, is cut to the REF, which bounds
# how long the reader holds a record and how many bases it holds of the reference.
DEFAULT_SPAN_LIMIT = 10_000
# The threads that parse the sample columns when the caller sets no number: the owner's
# choice of 19 October 2026, where the speed stopped growing. On the 64 Mb, 200-sample
# benchmark input at the default span limit, parse_vcf took 13.1 s with 1 thread, 3.9 s
# with 4 and 4.4 s with 8, the rest being the calling thread's own work; at 8 threads the
# peak memory was 64 MB against 37 MB, the same as on a 5 Mb input, so it does not grow
# with the file. The 4 was measured on one stream. Two streams read together from one
# thread, as the loci stage will merge files A and B, read their lines one after the other:
# on the 5 Mb GATK-like input against itself, one stream took 0.35 s, two zipped 0.75 s
# with 4 threads each and 0.76 s with 2 each, and each stream in a thread of its own
# 0.57 s. How two streams share the cores is for the loci or the command-line milestone.
DEFAULT_NUM_THREADS = 4
# The most threads parse_vcf takes. Past 4 the measured time stopped falling, and each
# thread adds its own batches and genotypes to the memory. 256 is far above the cores of
# the machines pvc_eval is meant for; a larger number is a mistake, which 2**62 showed by
# reaching the pool as itertools.islice's own ValueError.
_MAX_NUM_THREADS = 256
# How many bases on each side of a REF its span is first looked for in. A slide that runs
# into the end of them looks again in twice as many, up to one more than the span limit, so
# the bases a record asks the window for are a few times its slide, whatever the limit,
# and a record with no limit does not ask for its whole contig.
_FIRST_REACH = 64

type Span = tuple[int, int]
"""A stretch of a contig, 0-based and half-open."""


@dataclass(frozen=True, slots=True, kw_only=True)
class VcfGenotype:
    """One sample's GT at one VCF record, as written.

    alleles: one entry per haplotype, the index of a VCF allele (0 is REF), or MISSING.
    phased: the GT was written with `|`.
    phase_set: the PS of the sample, or None when the record gives none, which VCF
        reads as one phase set per sample and contig (§5.4).
    quality: the GQ of the sample, or None where absent.
    """

    alleles: tuple[int | MISSING, ...]
    phased: bool
    phase_set: str | None
    quality: int | None


@dataclass(frozen=True, slots=True, kw_only=True)
class VcfRecord:
    """One line of a VCF file.

    start: the first position REF covers, 0-based; `end` is one past the last.
    alleles: REF first, then the ALTs, in upper case and trimmed (§5.2); `*` is a
        spanning deletion.
    qual: QUAL, or None for `.`.
    genotypes: one per sample, in the order of the header; a MALFORMED_GT or HAPLOID_GT
        keeps no GQ.
    line_number: the 1-based line of the file, for messages.
    ref_mismatch: REF, as written, differs from the reference (§5.2).
    span: the record's span (§5.2), 0-based and half-open; it always contains the REF.
    """

    contig: str
    start: int
    alleles: tuple[str, ...]
    qual: float | None
    genotypes: tuple[VcfGenotype | MALFORMED_GT | HAPLOID_GT, ...]
    line_number: int
    ref_mismatch: bool
    span: Span

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be 0 or more, not {self.start}")
        if not self.alleles or not self.alleles[0]:
            raise ValueError("a VCF record needs a REF of at least one base")
        span_start, span_end = self.span
        if span_start > self.start or span_end < self.end:
            raise ValueError(
                f"span {span_start}-{span_end} must contain the REF, {self.start}-{self.end}"
            )

    @property
    def end(self) -> int:
        return self.start + len(self.alleles[0])


def parse_vcf_samples(path: Path) -> tuple[str, ...]:
    """The sample names of a VCF file, in the order of its columns.

    Raises PvcEvalInputError for a path that does not exist, is a directory or cannot be
    read, when the file has no #CHROM line after its ## lines, and for an empty sample name.
    """
    with closing(_read_numbered_lines(path)) as numbered_lines:
        return _parse_samples(numbered_lines, path)


def parse_vcf(
    path: Path,
    reference: IndexedFasta,
    *,
    span_limit: int | None = DEFAULT_SPAN_LIMIT,
    num_threads: int = DEFAULT_NUM_THREADS,
) -> Iterator[VcfRecord]:
    """The VCF records of a file, plain or bgzipped, normalised and given their spans
    (§5.2), in order of span start, and in the file's order where two spans start at one
    position.

    reference: the reference. The contigs, their order and their lengths are those of its
        `.fai`. The reader reads it through a FastaWindow of its own, which gives the bases
        in upper case, so a soft-masked reference gives no REF mismatch (§5.2).
    span_limit: the most bases a span may reach beyond its REF on either side, a whole
        number of 1 or more, or None for no limit (§5.2). A record whose span would reach further is given its REF as
        its span, and the number of those is logged as a warning at the end.
    num_threads: how many threads parse the sample columns, a whole number from 1 to 256.
        The calling thread reads the lines, checks them and gives them their spans. With 1
        it parses the sample columns too, and no pool is made: the run is serial, as a
        test or a profile may want. The records are the same
        whatever the number.

    With a limit, a record is held only until no later record of its contig can have a
    span that starts before its own: a span starts at most span_limit bases before its
    REF, and the file is sorted by POS. With no limit, the records of a contig are held
    until the next contig begins.

    The window holds the bases from span_limit + 1 before the last record read, or with
    no limit from BLOCK_BASES before it, to a block after the furthest base read, and
    once a slide reaches further back, up to a block more behind.

    A record whose REF differs from the reference is marked ref_mismatch and given its
    trimmed REF as its span: such a record makes its locus invalid for every sample of
    its file (§5.4), so its slide is not looked for.

    Raises PvcEvalInputError, when the stream reaches it, for a VCF record on a contig
    the reference lacks or whose REF reaches past the end of its contig, for records out
    of the `.fai`'s contig order or out of position order, for a path that does not exist,
    is a directory or cannot be read, for a file with no #CHROM line after its ## lines, for
    an empty sample name in it, for a line that does not hold a VCF record, and for one whose
    number of sample columns differs from the number of samples the #CHROM line names. A
    blank line is skipped. A record with a symbolic allele is skipped, and the number
    skipped is logged as a warning at the end (§9). A QUAL or GQ that is not a number is
    taken as absent, and those are counted in a warning at the end too. Raises
    ValueError, when the stream starts, for a span_limit
    that is neither None nor a whole number of 1 or more, and for a num_threads that is not
    a whole number from 1 to 256.
    """
    # §5.2: a limit of 0 would cut every insertion, whose span holds the base after it.
    if span_limit is not None and not _is_whole_number_from_1(span_limit):
        raise ValueError(
            f"span_limit must be a whole number of 1 or more, or None, not {span_limit!r}"
        )
    if not (_is_whole_number_from_1(num_threads) and num_threads <= _MAX_NUM_THREADS):
        raise ValueError(
            f"num_threads must be a whole number from 1 to {_MAX_NUM_THREADS}, not {num_threads!r}"
        )
    with closing(_parse_records_in_file_order(path, reference, span_limit, num_threads)) as records:
        yield from _sort_by_span_start(records, span_limit)


def _is_whole_number_from_1(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 1


def _sort_by_span_start(
    records: Iterable[tuple[int, VcfRecord]], span_limit: int | None
) -> Iterator[VcfRecord]:
    """The records in order of span start, and of the file where two spans start at one
    position.

    records: each record with the 0-based position of its POS, in the file's order.
    """
    # A heap of (span start, place in the file, record); the place keeps ties in the
    # file's order and is never equal, so records are never compared.
    held: list[tuple[int, int, VcfRecord]] = []
    contig = ""
    for place_in_file, (pos, vcf_record) in enumerate(records):
        if vcf_record.contig != contig:
            while held:
                yield heapq.heappop(held)[2]
            contig = vcf_record.contig
        heapq.heappush(held, (vcf_record.span[0], place_in_file, vcf_record))
        # With no limit, a later record of the contig can span from its first base.
        if span_limit is None:
            continue
        # Every later record of the contig has its POS at pos or after, and its REF,
        # trimmed, too, so its span starts at pos - span_limit or after.
        while held and held[0][0] <= pos - span_limit:
            yield heapq.heappop(held)[2]
    while held:
        yield heapq.heappop(held)[2]


# The sample columns of this many lines are parsed together, by one thread. With 8 threads
# on the 5 Mb, 200-sample GATK-like benchmark input, batches of 64 lines took 0.45 s and
# 64 MB, and of 256 lines 0.46 s and 108 MB.
_BATCH_SIZE = 64


def _parse_records_in_file_order(
    path: Path, reference: IndexedFasta, span_limit: int | None, num_threads: int
) -> Generator[tuple[int, VcfRecord]]:
    """The records of parse_vcf in the file's order, each with the 0-based
    position of its POS.

    This thread reads the lines, checks them and gives them their spans; the sample columns
    of each batch of lines are parsed by one of num_threads threads, or by this one when
    num_threads is 1, and the batches are taken back in the file's order.
    """
    status = _ReadStatus()
    with (
        closing(_parse_lines_without_samples(path, reference, span_limit, status)) as placed_lines,
        closing(
            # Each thread keeps what it has parsed in seen_by_thread, for this file alone.
            _map_in_order(
                functools.partial(_parse_batch, seen_by_thread=threading.local()),
                _batch_lines(placed_lines, status),
                num_threads,
            )
        ) as parsed_batches,
    ):
        for records, num_gq_not_numbers in parsed_batches:
            # The batches are taken back in this thread, so only this thread changes status.
            status.num_not_numbers += num_gq_not_numbers
            yield from records
    if status.error is not None:
        raise status.error
    _log_warnings(path, span_limit, status)


def _batch_lines(
    placed_lines: Iterator[_PlacedLine], status: _ReadStatus
) -> Generator[list[_PlacedLine]]:
    """The placed lines in batches of _BATCH_SIZE, the last one perhaps shorter.

    An input error ends the batches with the lines before it, and is kept in status, so
    that it is raised after their records, as it would be with no batches.
    """
    batch: list[_PlacedLine] = []
    try:
        for placed_line in placed_lines:
            batch.append(placed_line)
            if len(batch) == _BATCH_SIZE:
                yield batch
                batch = []
    except PvcEvalInputError as err:
        status.error = err
    if batch:
        yield batch


def _map_in_order[T, R](
    function: Callable[[T], R], items: Iterable[T], num_threads: int
) -> Generator[R]:
    """function over items, its results in the order of the items: in this thread when
    num_threads is 1, or else on a pool of num_threads threads, which take at most two items
    each ahead of the result given. items is read in this thread."""
    if num_threads == 1:
        yield from map(function, items)
        return
    with ThreadPoolExecutor(max_workers=num_threads) as pool:
        yield from pool.map(function, items, buffersize=2 * num_threads)


@dataclass(slots=True)
class _ReadStatus:
    """What the reader of a file warns about at its end, and the input error that stopped
    it; changed only by the thread that reads the lines."""

    num_symbolic: int = 0
    num_not_numbers: int = 0
    num_beyond_reach: int = 0
    first_beyond_reach_line: int = 0
    error: PvcEvalInputError | None = None


@dataclass(frozen=True, slots=True, kw_only=True)
class _PlacedLine:
    """A VCF line read, checked, normalised and given its span, whose sample columns are
    left to parse.

    pos: the 0-based position of POS. fields: the columns of the line. The others are the
    VcfRecord's.
    """

    pos: int
    contig: str
    start: int
    alleles: tuple[str, ...]
    qual: float | None
    line_number: int
    ref_mismatch: bool
    span: Span
    fields: list[str]


def _parse_lines_without_samples(
    path: Path, reference: IndexedFasta, span_limit: int | None, status: _ReadStatus
) -> Generator[_PlacedLine]:
    """The lines of a VCF file that hold a record, in the file's order, checked, normalised
    and given their spans, their sample columns left to parse; what is skipped and what is
    taken as absent is added to status."""
    contig_rank = {contig: rank for rank, contig in enumerate(reference.contigs)}
    # A window of the stream's own, so that two streams over one reference do not drop
    # each other's buffer.
    reference_window = FastaWindow(reference)
    # The (contig rank, POS) of the previous record, which this one must not come before,
    # and where that record is, for the message.
    previous_order: tuple[int, int] | None = None
    previous_where = ""
    with closing(_read_numbered_lines(path)) as numbered_lines:
        num_samples = len(_parse_samples(numbered_lines, path))
        for line_number, line in numbered_lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            contig, pos_1based = _parse_position(fields, path, line_number)
            _check_num_samples(fields, num_samples, path, line_number)
            if contig not in contig_rank:
                raise PvcEvalInputError(
                    f"{path}, line {line_number}: contig {contig} is not in the reference."
                    " Compare against the reference the VCF was called on."
                )
            order = (contig_rank[contig], pos_1based)
            where = f"{contig}:{pos_1based} (line {line_number})"
            if previous_order is not None and order < previous_order:
                raise PvcEvalInputError(
                    f"{path} is not sorted: {where} comes after {previous_where}."
                    " Sort it in the contig order of the reference's .fai and by position,"
                    " and run again."
                )
            previous_order, previous_where = order, where
            alts = _parse_alts(fields, path, line_number)
            if any(_is_symbolic(alt) for alt in alts):
                status.num_symbolic += 1
                continue
            if pos_1based == 0:
                raise _make_pos_error(fields[1], path, line_number)
            pos = pos_1based - 1
            contig_length = reference.get_length(contig)
            # §5.2: a REF past the end of its contig means the VCF was made on another
            # reference, as a contig the reference lacks does.
            if pos + len(fields[3]) > contig_length:
                raise PvcEvalInputError(
                    f"{path}, line {line_number}: the REF of {contig}:{pos_1based} reaches past"
                    f" the end of {contig}, which has {contig_length} bases in the reference,"
                    " so the VCF looks made on another reference. Compare against the"
                    " reference the VCF was called on."
                )
            # §5.2: with a limit, no base before pos - span_limit - 1 is needed again. This
            # record and every later one of the contig have their POS at pos or after (the
            # file is sorted), their trimmed REF too, and their span is looked for at most
            # span_limit + 1 bases before it; the records held for their order need no more
            # bases. With no limit no distance bounds a slide, so the window keeps a block
            # behind, and reads again what a longer slide asks for. A record on the next
            # contig makes the window drop the buffer.
            evict_distance = BLOCK_BASES if span_limit is None else span_limit + 1
            reference_window.evict_before(pos - evict_distance)
            placed_line, is_qual_absent, is_beyond_reach = _parse_line_without_samples(
                fields, contig, pos, alts, line_number, reference_window, contig_length, span_limit
            )
            status.num_not_numbers += is_qual_absent
            if is_beyond_reach:
                status.num_beyond_reach += 1
                status.first_beyond_reach_line = status.first_beyond_reach_line or line_number
            yield placed_line


def _log_warnings(path: Path, span_limit: int | None, status: _ReadStatus) -> None:
    """The warnings of the end of a file, for what its reader skipped, took as absent or
    cut."""
    if status.num_symbolic:
        _logger.warning(
            "%s: skipped %s with symbolic alleles (<DEL>, breakends); pvc_eval does not"
            " compare structural variants.",
            path,
            _format_count(status.num_symbolic, "VCF record"),
        )
    if status.num_not_numbers:
        _logger.warning(
            "%s: took %s that could not be read as a number as absent, as `.` is, in the"
            " quality of the table.",
            path,
            _format_count(status.num_not_numbers, "QUAL or GQ value"),
        )
    if status.num_beyond_reach:
        assert span_limit is not None, "only a limit cuts a span"
        if status.num_beyond_reach == 1:
            which, given = f"at line {status.first_beyond_reach_line}", "was given"
        else:
            which, given = f"the first at line {status.first_beyond_reach_line}", "each was given"
        _logger.warning(
            "%s: %s, %s, could slide more than %s along the reference and %s only"
            " its REF as its span; another spelling of the same event further along may"
            " fall into another locus and be compared apart.",
            path,
            _format_count(status.num_beyond_reach, "VCF record"),
            which,
            f"{span_limit:,} {'base' if span_limit == 1 else 'bases'}",
            given,
        )


type _ParsedGt = tuple[tuple[int | MISSING, ...], bool]
"""A GT's alleles, and whether it was written with `|`."""

_PARSED_MISSING_GT: _ParsedGt = ((MISSING,) * PLOIDY, False)


type _Genotype = VcfGenotype | MALFORMED_GT | HAPLOID_GT
"""What a sample's GT parses to."""

type _GenotypesSeen = dict[tuple[str, int, str | None, int | None], _Genotype]
"""The genotypes parsed, by GT, number of alleles, PS and GQ."""

type _GtsSeen = dict[tuple[str, int], _ParsedGt | MALFORMED_GT | HAPLOID_GT]
"""What _parse_gt gave, by GT and number of alleles."""


def _parse_genotype(
    gt: str, num_alleles: int, phase_set: str | None, quality: int | None, gts_seen: _GtsSeen
) -> _Genotype:
    """A GT field as a VcfGenotype; HAPLOID_GT for one allele index, such as `1`; or
    MALFORMED_GT when it cannot be either (§5.4).

    num_alleles: how many VCF alleles the record has, REF included.
    phase_set, quality: the sample's PS and GQ, kept in the VcfGenotype.
    gts_seen: the GTs parsed so far, which this adds to.
    """
    gt_key = (gt, num_alleles)
    parsed_gt = gts_seen.get(gt_key)
    if parsed_gt is None:
        parsed_gt = gts_seen[gt_key] = _parse_gt(gt, num_alleles)
    if parsed_gt is MALFORMED_GT or parsed_gt is HAPLOID_GT:
        return parsed_gt
    alleles, phased = parsed_gt
    return VcfGenotype(alleles=alleles, phased=phased, phase_set=phase_set, quality=quality)


def _parse_gt(gt: str, num_alleles: int) -> _ParsedGt | MALFORMED_GT | HAPLOID_GT:
    """What a GT field says whatever the sample's PS and GQ: its alleles and whether it is
    phased, HAPLOID_GT or MALFORMED_GT (§5.4).

    A single `.` says nothing about either haplotype, and is both MISSING.
    """
    if gt == ".":
        return _PARSED_MISSING_GT
    parts = gt.replace("|", "/").split("/")
    if len(parts) == 1 and _is_allele_index(gt, num_alleles):
        return HAPLOID_GT
    if len(parts) != PLOIDY:
        return MALFORMED_GT
    alleles: list[int | MISSING] = []
    for part in parts:
        if part == ".":
            alleles.append(MISSING)
        elif _is_allele_index(part, num_alleles):
            alleles.append(int(part))
        else:
            return MALFORMED_GT
    return tuple(alleles), "|" in gt


def _is_allele_index(part: str, num_alleles: int) -> bool:
    index = _parse_decimal(part)
    return index is not None and index < num_alleles


# int() reads from text every number of fewer digits than this, whatever
# sys.set_int_max_str_digits sets, so only longer text is checked against its limit.
_DIGITS_ALWAYS_READ = sys.int_info.str_digits_check_threshold


def _parse_decimal(text: str) -> int | None:
    """The whole number text writes in ASCII decimal digits, or None when it is not one or
    has more digits than int() reads from text (4,300 by default)."""
    if not (text.isascii() and text.isdecimal()):
        return None
    if len(text) >= _DIGITS_ALWAYS_READ and _has_too_many_digits(text):
        return None
    return int(text)


def _has_too_many_digits(text: str) -> bool:
    """text has more digits than int() reads (sys.get_int_max_str_digits(), 0 for none)."""
    max_digits = sys.get_int_max_str_digits()
    return max_digits != 0 and len(text) > max_digits


def _read_numbered_lines(path: Path) -> Generator[tuple[int, str]]:
    """The lines of a VCF file, plain or bgzipped, each with its 1-based number.

    Raises PvcEvalInputError for a file that is not UTF-8 text, a BCF among them, and
    for a gzip file that is cut short or damaged.
    """
    # Text is decoded and decompressed in blocks, so a fault is found somewhere after the
    # last line read, and only that is known.
    line_number = 0
    try:
        with _open_text(path) as lines:
            for line_number, line in enumerate(lines, start=1):
                yield line_number, line
    except UnicodeDecodeError as err:
        raise PvcEvalInputError(
            f"{path} is not UTF-8 text: a byte at or after line {line_number + 1} cannot be read"
            " as UTF-8. If it is a BCF file, convert it to a VCF with"
            " `bcftools view -Oz -o file.vcf.gz file.bcf`; if it is a VCF in another"
            " encoding, convert it to UTF-8; and run again."
        ) from err
    except (EOFError, gzip.BadGzipFile) as err:
        raise PvcEvalInputError(
            f"{path} is a gzip file cut short or damaged at or after line {line_number + 1}."
            " Compress or copy it again, and run again."
        ) from err


def _open_text(path: Path) -> TextIO:
    """Raises PvcEvalInputError for a path that does not exist, is a directory or cannot
    be read."""
    try:
        with path.open("rb") as raw:
            is_gzipped = raw.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
    except FileNotFoundError:
        raise PvcEvalInputError(f"{path} does not exist. Check its path and run again.") from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{path} is a directory, not a VCF file. Check its path and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{path} cannot be read: permission denied. Allow reading it and run again."
        ) from None
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _parse_samples(numbered_lines: Iterator[tuple[int, str]], path: Path) -> tuple[str, ...]:
    """The sample names of the #CHROM line, read up to it and no further."""
    for line_number, line in numbered_lines:
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#CHROM"):
            return _parse_sample_names(line, path, line_number)
        raise PvcEvalInputError(
            f"{path}, line {line_number}: a VCF file has a #CHROM line after its ## lines,"
            " and this line comes before any. Check that the file is a VCF with its header."
        )
    raise PvcEvalInputError(
        f"{path} has no #CHROM line after its ## lines. Check that the file is a VCF with"
        " its header."
    )


def _parse_sample_names(chrom_line: str, path: Path, line_number: int) -> tuple[str, ...]:
    """The sample names of a #CHROM line; raises PvcEvalInputError for an empty one, which
    a trailing tab or two tabs in a row make."""
    names = tuple(chrom_line.rstrip("\n").split("\t")[_NUM_FIXED_COLUMNS:])
    if "" in names:
        column = _NUM_FIXED_COLUMNS + names.index("") + 1
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the #CHROM line has an empty sample name in column"
            f" {column}; each sample needs a name, and the columns are separated by one tab."
            " Correct the line and run again."
        )
    return names


def _check_num_samples(fields: list[str], num_samples: int, path: Path, line_number: int) -> None:
    """Later stages find a sample's genotype by its column, so every line has one per
    sample of the #CHROM line, after FORMAT. A sites-only file has neither."""
    num_found = max(len(fields) - _NUM_FIXED_COLUMNS, 0)
    if num_found != num_samples:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has {_format_count(num_found, 'sample column')},"
            f" and the #CHROM line names {_format_count(num_samples, 'sample')}. A VCF record has a FORMAT column and"
            " then one column per sample; correct the line or remove it, and run again."
        )


# A line that is not a VCF record stops the run (a provisional choice: doc/spec.md §3 does
# not list it), with the column at fault.
def _parse_position(fields: list[str], path: Path, line_number: int) -> tuple[str, int]:
    if len(fields) < _NUM_FIXED_COLUMNS - 1:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has"
            f" {_format_count(len(fields), 'tab-separated column')},"
            " and a VCF record has at least 8, CHROM to INFO. Check that its columns are"
            " separated by tabs, correct the line or remove it, and run again."
        )
    pos = fields[1]
    # VCF 4.3 writes a telomeric breakend at POS 0, so 0 is refused only once the record
    # is known not to be symbolic.
    if not (pos.isascii() and pos.isdecimal()):
        raise _make_pos_error(pos, path, line_number)
    pos_1based = _parse_decimal(pos)
    if pos_1based is None:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: POS has {len(pos):,} digits, which is past the end of"
            " any contig. Correct the line or remove it, and run again."
        )
    return fields[0], pos_1based


def _format_count(number: int, noun: str) -> str:
    """A number and its noun, in the singular for 1: "1 sample", "3 samples"."""
    return f"{number} {noun}" if number == 1 else f"{number} {noun}s"


def _make_pos_error(pos: str, path: Path, line_number: int) -> PvcEvalInputError:
    return PvcEvalInputError(
        f"{path}, line {line_number}: POS is '{pos}', and it must be a whole number"
        " of 1 or more. Correct the line or remove it, and run again."
    )


def _parse_alts(fields: list[str], path: Path, line_number: int) -> tuple[str, ...]:
    """The ALTs of a line, () for `.`, after checking that REF and every ALT hold
    something."""
    if not fields[3]:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: REF is empty, and it must hold at least one base."
            " Correct the line or remove it, and run again."
        )
    if fields[4] == ".":
        return ()
    alts = tuple(fields[4].split(","))
    if not all(alts):
        raise PvcEvalInputError(
            f"{path}, line {line_number}: ALT '{fields[4]}' has an empty allele; each ALT"
            " must hold at least one base, and a record with no ALT has `.`. Correct the"
            " line or remove it, and run again."
        )
    return alts


def _is_symbolic(alt: str) -> bool:
    """A symbolic ALT (`<DEL>`), a breakend (`A]chr2:20]`) or a single breakend of VCF 4.3
    (`A.`, `.A`); `.` alone is no ALT at all."""
    is_single_breakend = len(alt) > 1 and (alt.startswith(".") or alt.endswith("."))
    return alt.startswith("<") or "[" in alt or "]" in alt or is_single_breakend


def _parse_line_without_samples(
    fields: list[str],
    contig: str,
    pos: int,
    alts: tuple[str, ...],
    line_number: int,
    reference_window: FastaWindow,
    contig_length: int,
    span_limit: int | None,
) -> tuple[_PlacedLine, bool, bool]:
    """A line normalised and with its span, its sample columns left to parse; whether its
    QUAL was taken as absent for not being a number; and whether its span was cut to its
    REF for reaching more than span_limit beyond it.

    pos: the 0-based position of POS; the REF lies within the contig.
    """
    # VCF bases are case-insensitive, and the window gives the reference in upper case.
    alleles = (fields[3].upper(), *(alt.upper() for alt in alts))
    ref_bases = reference_window.get_bases(contig, pos, pos + len(alleles[0]))
    ref_mismatch = ref_bases != alleles[0]
    start, trimmed = _trim(pos, alleles)
    ref_span = (start, start + len(trimmed[0]))
    span = (
        ref_span
        if ref_mismatch
        else _find_span_in_window(
            start, trimmed, reference_window, contig, contig_length, span_limit
        )
    )
    is_beyond_reach = span is None
    if span is None:
        span = ref_span
    qual = _parse_qual(fields[5])
    placed_line = _PlacedLine(
        pos=pos,
        contig=contig,
        start=start,
        alleles=trimmed,
        qual=qual,
        line_number=line_number,
        ref_mismatch=ref_mismatch,
        span=span,
        fields=fields,
    )
    return placed_line, _is_taken_as_absent(fields[5], qual), is_beyond_reach


# GTs and whole genotypes repeat, even with a phase set per sample, so each thread keeps
# those it has parsed: equal genotypes then are one frozen object, and genotypes with one GT
# share one alleles tuple. On the 5 Mb, 200-sample GATK-like benchmark input, of 3,625,400
# genotypes one thread builds 515,559 and eight threads 604,124; a cache shared by the
# threads, functools.lru_cache, would take a lock on every sample.
# A thread starts them again, genotypes and GTs, when a batch starts and it holds more
# than this many genotypes, so it holds at most this many and one batch's, 64 lines times
# the samples of a line; a GT is kept only beside a genotype, so the GTs are bounded too.
_MAX_GENOTYPES_SEEN = 4096


def _parse_batch(
    placed_lines: Iterable[_PlacedLine], seen_by_thread: threading.local
) -> tuple[list[tuple[int, VcfRecord]], int]:
    """The VCF records of a batch of placed lines, each with the 0-based position of its
    POS, and how many of their GQ values were taken as absent for not being numbers.

    seen_by_thread: where the thread that runs this keeps the genotypes and the GTs it has
        parsed, as `genotypes` and `gts`, so that no two threads share them.
    """
    if (
        not hasattr(seen_by_thread, "genotypes")
        or len(seen_by_thread.genotypes) > _MAX_GENOTYPES_SEEN
    ):
        seen_by_thread.genotypes = {}
        seen_by_thread.gts = {}
    genotypes_seen: _GenotypesSeen = seen_by_thread.genotypes
    gts_seen: _GtsSeen = seen_by_thread.gts
    records: list[tuple[int, VcfRecord]] = []
    num_not_numbers = 0
    for line in placed_lines:
        genotypes, num_line_not_numbers = _parse_genotypes(
            line.fields, len(line.alleles), genotypes_seen, gts_seen
        )
        num_not_numbers += num_line_not_numbers
        record = VcfRecord(
            contig=line.contig,
            start=line.start,
            alleles=line.alleles,
            qual=line.qual,
            genotypes=genotypes,
            line_number=line.line_number,
            ref_mismatch=line.ref_mismatch,
            span=line.span,
        )
        records.append((line.pos, record))
    return records, num_not_numbers


# An index past the end of every sample field, for a key that FORMAT lacks.
_ABSENT_KEY_INDEX = sys.maxsize


def _parse_genotypes(
    fields: list[str], num_alleles: int, genotypes_seen: _GenotypesSeen, gts_seen: _GtsSeen
) -> tuple[tuple[_Genotype, ...], int]:
    """The genotypes of a line's sample columns, one per sample, and how many of their GQ
    values were taken as absent for not being numbers.

    num_alleles: how many VCF alleles the record has, REF included.
    genotypes_seen, gts_seen: the genotypes and GTs parsed so far, which this adds to.
    """
    if len(fields) <= _NUM_FIXED_COLUMNS:
        return (), 0
    genotypes: list[_Genotype] = []
    num_not_numbers = 0
    keys = fields[8].split(":")
    # The keys are found once per record and each sample read by position, rather than
    # into a dict of its own: parse_vcf is 39% faster on the 5 Mb, 200-sample GATK-like
    # benchmark input, and 42% on the 64 Mb one. VCF does not allow a key twice; when it
    # is, the last copy counts, as it did when each sample was read into a dict. Unlike that
    # dict, a sample that drops the last copy has no value for the key, even when it has one
    # at an earlier copy; so a GQ that is not a number at that earlier copy is not counted
    # as taken as absent.
    key_indices = {key: index for index, key in enumerate(keys)}
    gt_index = key_indices.get("GT", _ABSENT_KEY_INDEX)
    gq_index = key_indices.get("GQ", _ABSENT_KEY_INDEX)
    ps_index = key_indices.get("PS", _ABSENT_KEY_INDEX)
    for sample_field in fields[_NUM_FIXED_COLUMNS:]:
        # VCF lets trailing sample fields be dropped, so an index past the last value, or
        # the index of a key FORMAT lacks, has no value; a GT with none says what `.` says.
        values = sample_field.split(":")
        num_values = len(values)
        gt = values[gt_index] if gt_index < num_values else "."
        gq = values[gq_index] if gq_index < num_values else None
        phase_set = values[ps_index] if ps_index < num_values else None
        if phase_set == ".":
            phase_set = None
        quality = _parse_gq(gq)
        num_not_numbers += _is_taken_as_absent(gq, quality)
        # A genotype seen before is looked up here, and not in _parse_genotype, to save a
        # call per sample: most samples are seen before.
        genotype_key = (gt, num_alleles, phase_set, quality)
        genotype = genotypes_seen.get(genotype_key)
        if genotype is None:
            genotype = genotypes_seen[genotype_key] = _parse_genotype(
                gt, num_alleles, phase_set, quality, gts_seen
            )
        genotypes.append(genotype)
    return tuple(genotypes), num_not_numbers


# QUAL and GQ only annotate a row of the table (§6), where an absent value is written `.`;
# a value that is not a number is taken as absent and counted, rather than stopping the
# run (a provisional choice). Only decimal numbers count, not Python's "nan", "inf" or
# "1_000", and a QUAL must be finite.
def _parse_qual(field: str) -> float | None:
    if not _DECIMAL_NUMBER.fullmatch(field):
        return None
    qual = float(field)
    return qual if math.isfinite(qual) else None


def _parse_gq(field: str | None) -> int | None:
    # _parse_decimal's body, written out since this runs once per sample.
    if field is None or not field.isascii() or not field.isdecimal():
        return None
    if len(field) >= _DIGITS_ALWAYS_READ and _has_too_many_digits(field):
        return None
    return int(field)


def _is_taken_as_absent(field: str | None, number: float | None) -> bool:
    """A QUAL or GQ was written, and not as `.`, but taken as absent."""
    return number is None and field not in (None, ".")


def _trim(start: int, alleles: Sequence[str]) -> tuple[int, tuple[str, ...]]:
    """A record's start and alleles without the bases all its VCF alleles share, at the
    right end, then at the left, keeping at least one base in each (the VCF anchor
    convention); start is 0-based.

    A `*` allele has no bases and takes no part. A record with fewer than two other
    alleles is returned as it is.
    """
    sequences = _drop_spanning_deletions(alleles)
    if len(sequences) < 2:
        return start, tuple(alleles)
    start, trimmed = _trim_sequences(start, sequences)
    return start, _make_trimmed_alleles(alleles, trimmed)


def _find_span_in_window(
    start: int,
    alleles: Sequence[str],
    reference_window: FastaWindow,
    contig: str,
    contig_length: int,
    span_limit: int | None,
) -> Span | None:
    """The span of a trimmed record whose REF agrees with the reference (§5.2), or None
    when it would reach more than span_limit bases beyond the REF on either side; with no
    limit, never None.

    start: where REF starts, 0-based.

    The span is looked for in the bases within a reach of the REF, and again in twice as
    many while a slide runs into their end short of the contig's. With a limit the reach
    stops at span_limit + 1: a slide that runs into the end of those bases reaches beyond
    the limit, and how far it goes is not needed.
    """
    end = start + len(alleles[0])
    reach = _FIRST_REACH
    while True:
        if span_limit is not None:
            reach = min(reach, span_limit + 1)
        stretch_start = max(0, start - reach)
        stretch_end = min(contig_length, end + reach)
        stretch = reference_window.get_bases(contig, stretch_start, stretch_end)
        span_start, span_end = _find_span(start, alleles, stretch, stretch_start, contig_length)
        # A slide that ran into an end of the stretch short of the contig's gives a span
        # that starts before the stretch or ends at its end or after (see _find_span). An
        # insertion that stops one base short of the end gives that too, and is only
        # looked for again.
        is_cut_by_stretch = span_start < stretch_start or (
            span_end >= stretch_end and stretch_end < contig_length
        )
        if not is_cut_by_stretch or (span_limit is not None and reach > span_limit):
            break
        reach *= 2
    if span_limit is not None and (start - span_start > span_limit or span_end - end > span_limit):
        return None
    return span_start, span_end


def _find_span(
    start: int, alleles: Sequence[str], stretch: str, stretch_start: int, contig_length: int
) -> Span:
    """The span of a trimmed record whose REF agrees with the reference (§5.2), as far as a
    stretch of its contig shows it.

    start: where REF starts, 0-based. stretch: the bases of the contig from stretch_start,
    0-based, in upper case, the REF among them. contig_length: the number of bases of the
    contig, which tells the end of the stretch from the end of the contig.

    The span is the union of the REF and of the span of each ALT. An ALT that is an
    insertion or a deletion spans from the anchor of its left-aligned form to the last
    base it deletes, or the base after the point it is inserted at, when slid as far
    right as it goes. A SNP, a complex ALT and a `*` add nothing to the REF. A slide stops
    at the ends of the stretch; one that stops there, short of the ends of the contig,
    gives a span that starts one base before the stretch or ends at its end or after.
    """
    end = start + len(alleles[0])
    span_start, span_end = start, end
    for alt in alleles[1:]:
        if alt == _SPANNING_DELETION:
            continue
        indel_span = _find_indel_span(start, alleles[0], alt, stretch, stretch_start, contig_length)
        if indel_span is not None:
            span_start = min(span_start, indel_span[0])
            span_end = max(span_end, indel_span[1])
    return span_start, span_end


def _find_indel_span(
    start: int, ref: str, alt: str, stretch: str, stretch_start: int, contig_length: int
) -> Span | None:
    """The span of one ALT when it is an insertion or a deletion, 0-based and half-open,
    as far as the stretch shows it (see _find_span); None for a SNP or a complex ALT.

    The ALT is an indel when, without the bases it shares with REF at both ends, one of
    the two is empty, whichever side the shared bases were on.
    """
    num_prefix = _count_shared(ref, alt)
    deleted, inserted = ref[num_prefix:], alt[num_prefix:]
    num_suffix = _count_shared(deleted[::-1], inserted[::-1])
    deleted = deleted[: len(deleted) - num_suffix]
    inserted = inserted[: len(inserted) - num_suffix]
    if bool(deleted) == bool(inserted):
        return None
    event = deleted or inserted
    # The event is the deletion of [position, after), or an insertion before position,
    # where after is position too.
    # Positions in the stretch, not in the contig.
    position = start + num_prefix - stretch_start
    after = position + len(deleted)
    # The event slides left while the base before it equals its last base, and right
    # while the base after it equals its first; each slide turns the event round by one
    # base, so its bases are read round from the end it slid from.
    num_left = 0
    while (
        position - num_left > 0
        and stretch[position - num_left - 1] == event[(-1 - num_left) % len(event)]
    ):
        num_left += 1
    num_right = 0
    while (
        after + num_right < len(stretch)
        and stretch[after + num_right] == event[num_right % len(event)]
    ):
        num_right += 1
    # The anchor is the base before the left-aligned event; an event slid to the first
    # base of the contig has none, and its span starts there. An insertion's span holds
    # the base after its point (§5.2: so that C -> CT at 10 overlaps a record at 11),
    # unless it is inserted after the last base of the contig.
    span_start = max(stretch_start + position - num_left - 1, 0)
    span_end = stretch_start + after + num_right + (1 if inserted else 0)
    return span_start, min(span_end, contig_length)


def _count_shared(first: str, second: str) -> int:
    """How many bases two sequences share at their start."""
    for index, (base, other_base) in enumerate(zip(first, second, strict=False)):
        if base != other_base:
            return index
    return min(len(first), len(second))


def _trim_sequences(start: int, alleles: Sequence[str]) -> tuple[int, list[str]]:
    trimmed = list(alleles)
    while min(map(len, trimmed)) > 1 and len({allele[-1] for allele in trimmed}) == 1:
        trimmed = [allele[:-1] for allele in trimmed]
    while min(map(len, trimmed)) > 1 and len({allele[0] for allele in trimmed}) == 1:
        trimmed = [allele[1:] for allele in trimmed]
        start += 1
    return start, trimmed


def _drop_spanning_deletions(alleles: Sequence[str]) -> list[str]:
    """The alleles that are sequences. A `*` has no bases and takes no part in trimming;
    _make_trimmed_alleles returns it to its place."""
    return [allele for allele in alleles if allele != _SPANNING_DELETION]


def _make_trimmed_alleles(alleles: Sequence[str], sequences: list[str]) -> tuple[str, ...]:
    """The alleles with the sequences _drop_spanning_deletions took out replaced, in order, by these."""
    replacements = iter(sequences)
    return tuple(
        allele if allele == _SPANNING_DELETION else next(replacements) for allele in alleles
    )
