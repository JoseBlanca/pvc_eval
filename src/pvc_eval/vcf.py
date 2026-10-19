"""The VCF records of one file: stage 2 of doc/spec.md §4, read with the checks of §5.2,
normalised, given their spans and put in order of span start. Normalising trims the bases
a record's alleles share and marks a REF that differs from the reference; a VCF record is
never moved, split or decomposed.

Positions are 0-based and half-open inside pvc_eval. The 1-based POS of a VCF line is
converted here, and nowhere else.
"""

import gzip
import heapq
import logging
import math
import re
from collections.abc import Callable, Generator, Iterable, Iterator, Sequence
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError

logger = logging.getLogger(__name__)

PLOIDY = 2  # doc/spec.md §2: a global constant

MISSING = sentinel("MISSING")
"""A haplotype the VCF genotype says nothing about, written `.` (doc/spec.md §2)."""

MALFORMED_GT = sentinel("MALFORMED_GT")
"""A GT with a number of alleles other than the ploidy, or an allele index the VCF
record lacks: the sample is invalid at the locus, and the run goes on (§5.4)."""

HAPLOID_GT = sentinel("HAPLOID_GT")
"""A GT of one allele index, such as `1` on chrX of a male sample: the sample is invalid
at the locus with its own reason, `haploid GT` (§5.4, §6, Open 1)."""

SPANNING_DELETION = "*"
GZIP_MAGIC = b"\x1f\x8b"
DECIMAL_NUMBER = re.compile(r"[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?")
NUM_FIXED_COLUMNS = 9  # CHROM to FORMAT
# doc/spec.md §5.2: a span that would reach further than this beyond its REF, on either
# side, is cut to the REF, which bounds how long the reader holds a record.
MAX_SPAN_REACH = 10_000

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
    span_start, span_end: the record's span (§5.2), 0-based and half-open; it always
        contains the REF.
    """

    contig: str
    start: int
    alleles: tuple[str, ...]
    qual: float | None
    genotypes: tuple[VcfGenotype | MALFORMED_GT | HAPLOID_GT, ...]
    line_number: int
    ref_mismatch: bool
    span_start: int
    span_end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be 0 or more, not {self.start}")
        if not self.alleles or not self.alleles[0]:
            raise ValueError("a VCF record needs a REF of at least one base")
        if self.span_start > self.start or self.span_end < self.end:
            raise ValueError(
                f"span {self.span_start}-{self.span_end} must contain the REF,"
                f" {self.start}-{self.end}"
            )

    @property
    def end(self) -> int:
        return self.start + len(self.alleles[0])


def read_samples(path: Path) -> tuple[str, ...]:
    """The sample names of a VCF file, in the order of its columns.

    Raises PvcEvalInputError when the file has no #CHROM line after its ## lines.
    """
    with closing(_numbered_lines(path)) as numbered_lines:
        return _read_header(numbered_lines, path)


def read_vcf_records(
    path: Path, contigs: Sequence[str], contig_sequence: Callable[[str], str]
) -> Iterator[VcfRecord]:
    """The VCF records of a file, plain or bgzipped, normalised and given their spans
    (§5.2), in order of span start, and in the file's order where two spans start at one
    position.

    contigs: the contig names of the reference, in the order of its `.fai`.
    contig_sequence: the sequence of a contig of the reference, by name, in upper or
        lower case; it is compared in upper case, so a soft-masked reference gives no
        REF mismatch (§5.2). It is called once for each contig that has a VCF record.

    A record is held only until no later record of its contig can have a span that
    starts before its own: a span starts at most MAX_SPAN_REACH bases before its REF,
    and the file is sorted by POS. A record whose span would reach further than that
    beyond its REF is given its REF as its span, and the number of those is logged as a
    warning at the end. A record whose REF differs from the reference is marked
    ref_mismatch and given its REF as its span, since it cannot slide along a reference
    it disagrees with.

    Raises PvcEvalInputError, when the stream reaches it, for a VCF record on a contig
    the reference lacks, for records out of the `.fai`'s contig order or out of position
    order, for a file with no #CHROM line after its ## lines, for a line that does not
    hold a VCF record, and for one whose number of sample columns differs from the
    number of samples the #CHROM line names. A blank line is skipped. A record with a symbolic
    allele is skipped, and the number skipped is logged as a warning at the end (§9). A
    QUAL or GQ that is not a number is taken as absent, and those are counted in a
    warning at the end too.
    """
    with closing(_records_in_file_order(path, contigs, contig_sequence)) as records:
        yield from _in_span_order(records)


def _in_span_order(records: Iterable[tuple[int, VcfRecord]]) -> Iterator[VcfRecord]:
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
        heapq.heappush(held, (vcf_record.span_start, place_in_file, vcf_record))
        # Every later record of the contig has its POS at pos or after, and its REF,
        # trimmed, too, so its span starts at pos - MAX_SPAN_REACH or after.
        while held and held[0][0] <= pos - MAX_SPAN_REACH:
            yield heapq.heappop(held)[2]
    while held:
        yield heapq.heappop(held)[2]


def _records_in_file_order(
    path: Path, contigs: Sequence[str], contig_sequence: Callable[[str], str]
) -> Generator[tuple[int, VcfRecord]]:
    """The records of read_vcf_records in the file's order, each with the 0-based
    position of its POS."""
    contig_rank = {contig: rank for rank, contig in enumerate(contigs)}
    # The (contig rank, POS) of the previous record, which this one must not come before,
    # and where that record is, for the message.
    previous_order: tuple[int, int] | None = None
    previous_where = ""
    num_symbolic = 0
    num_not_numbers = 0
    num_beyond_reach = 0
    sequence_contig = ""
    sequence = ""
    with closing(_numbered_lines(path)) as numbered_lines:
        num_samples = len(_read_header(numbered_lines, path))
        for line_number, line in numbered_lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            contig, pos_1based = _read_position(fields, path, line_number)
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
            alts = _read_alts(fields, path, line_number)
            if any(_is_symbolic(alt) for alt in alts):
                num_symbolic += 1
                continue
            if pos_1based == 0:
                raise _pos_error(fields[1], path, line_number)
            if contig != sequence_contig:
                sequence_contig, sequence = contig, contig_sequence(contig).upper()
            record, num_record_not_numbers, is_beyond_reach = _parse_record(
                fields, contig, pos_1based - 1, alts, line_number, sequence
            )
            num_not_numbers += num_record_not_numbers
            num_beyond_reach += is_beyond_reach
            yield pos_1based - 1, record
    if num_symbolic:
        logger.warning(
            "%s: skipped %s with symbolic alleles (<DEL>, breakends); pvc_eval does not"
            " compare structural variants.",
            path,
            _count(num_symbolic, "VCF record"),
        )
    if num_not_numbers:
        logger.warning(
            "%s: took %s that could not be read as a number as absent, as `.` is, in the"
            " quality of the table.",
            path,
            _count(num_not_numbers, "QUAL or GQ value"),
        )
    if num_beyond_reach:
        logger.warning(
            "%s: %s could slide more than 10 kb along the reference, and %s given only"
            " its REF as its span; another spelling of the same event further along may"
            " fall into another locus and be compared apart.",
            path,
            _count(num_beyond_reach, "VCF record"),
            "was" if num_beyond_reach == 1 else "were each",
        )


def parse_genotype(
    gt: str, num_alleles: int, phase_set: str | None, quality: int | None
) -> VcfGenotype | MALFORMED_GT | HAPLOID_GT:
    """A GT field as a VcfGenotype; HAPLOID_GT for one allele index, such as `1`; or
    MALFORMED_GT when it cannot be either (§5.4).

    A single `.` says nothing about either haplotype, and is both MISSING.

    num_alleles: how many VCF alleles the record has, REF included.
    phase_set, quality: the sample's PS and GQ, kept in the VcfGenotype.
    """
    if gt == ".":
        return VcfGenotype(
            alleles=(MISSING,) * PLOIDY, phased=False, phase_set=phase_set, quality=quality
        )
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
    return VcfGenotype(
        alleles=tuple(alleles), phased="|" in gt, phase_set=phase_set, quality=quality
    )


def _is_allele_index(part: str, num_alleles: int) -> bool:
    return part.isascii() and part.isdecimal() and int(part) < num_alleles


def _numbered_lines(path: Path) -> Generator[tuple[int, str]]:
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
    with path.open("rb") as raw:
        is_gzipped = raw.read(len(GZIP_MAGIC)) == GZIP_MAGIC
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _read_header(numbered_lines: Iterator[tuple[int, str]], path: Path) -> tuple[str, ...]:
    """The sample names of the #CHROM line, read up to it and no further."""
    for line_number, line in numbered_lines:
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#CHROM"):
            return tuple(line.rstrip("\n").split("\t")[NUM_FIXED_COLUMNS:])
        raise PvcEvalInputError(
            f"{path}, line {line_number}: a VCF file has a #CHROM line after its ## lines,"
            " and this line comes before any. Check that the file is a VCF with its header."
        )
    raise PvcEvalInputError(
        f"{path} has no #CHROM line after its ## lines. Check that the file is a VCF with"
        " its header."
    )


def _check_num_samples(fields: list[str], num_samples: int, path: Path, line_number: int) -> None:
    """Later stages find a sample's genotype by its column, so every line has one per
    sample of the #CHROM line, after FORMAT. A sites-only file has neither."""
    num_found = max(len(fields) - NUM_FIXED_COLUMNS, 0)
    if num_found != num_samples:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has {_count(num_found, 'sample column')},"
            f" and the #CHROM line names {_count(num_samples, 'sample')}. A VCF record has a FORMAT column and"
            " then one column per sample; correct the line or remove it, and run again."
        )


# A line that is not a VCF record stops the run (a provisional choice: doc/spec.md §3 does
# not list it), with the column at fault.
def _read_position(fields: list[str], path: Path, line_number: int) -> tuple[str, int]:
    if len(fields) < NUM_FIXED_COLUMNS - 1:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has"
            f" {_count(len(fields), 'tab-separated column')},"
            " and a VCF record has at least 8, CHROM to INFO. Check that its columns are"
            " separated by tabs, correct the line or remove it, and run again."
        )
    pos = fields[1]
    # VCF 4.3 writes a telomeric breakend at POS 0, so 0 is refused only once the record
    # is known not to be symbolic.
    if not (pos.isascii() and pos.isdecimal()):
        raise _pos_error(pos, path, line_number)
    return fields[0], int(pos)


def _count(number: int, noun: str) -> str:
    """A number and its noun, in the singular for 1: "1 sample", "3 samples"."""
    return f"{number} {noun}" if number == 1 else f"{number} {noun}s"


def _pos_error(pos: str, path: Path, line_number: int) -> PvcEvalInputError:
    return PvcEvalInputError(
        f"{path}, line {line_number}: POS is '{pos}', and it must be a whole number"
        " of 1 or more. Correct the line or remove it, and run again."
    )


def _read_alts(fields: list[str], path: Path, line_number: int) -> tuple[str, ...]:
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


def _parse_record(
    fields: list[str],
    contig: str,
    pos: int,
    alts: tuple[str, ...],
    line_number: int,
    contig_sequence: str,
) -> tuple[VcfRecord, int, bool]:
    """The VCF record of a line, normalised and with its span; how many of its QUAL and
    GQ values were taken as absent for not being numbers; and whether its span was cut
    to its REF for reaching more than MAX_SPAN_REACH beyond it.

    pos: the 0-based position of POS.
    """
    # VCF bases are case-insensitive, and the reference is compared in upper case.
    alleles = (fields[3].upper(), *(alt.upper() for alt in alts))
    ref_mismatch = contig_sequence[pos : pos + len(alleles[0])] != alleles[0]
    start, trimmed = trim(pos, alleles)
    ref_span = (start, start + len(trimmed[0]))
    span = ref_span if ref_mismatch else vcf_record_span(start, trimmed, contig_sequence)
    is_beyond_reach = span is None
    if span is None:
        span = ref_span
    genotypes: list[VcfGenotype | MALFORMED_GT | HAPLOID_GT] = []
    qual = _read_qual(fields[5])
    num_not_numbers = _is_not_a_number(fields[5], qual)
    if len(fields) > NUM_FIXED_COLUMNS:
        keys = fields[8].split(":")
        for sample_field in fields[NUM_FIXED_COLUMNS:]:
            # VCF lets trailing sample fields be dropped, so a key may have no value.
            values = dict(zip(keys, sample_field.split(":"), strict=False))
            phase_set = values.get("PS")
            if phase_set == ".":
                phase_set = None
            quality = _read_gq(values.get("GQ"))
            num_not_numbers += _is_not_a_number(values.get("GQ"), quality)
            if "GT" not in values:
                genotypes.append(
                    VcfGenotype(
                        alleles=(MISSING,) * PLOIDY,
                        phased=False,
                        phase_set=phase_set,
                        quality=quality,
                    )
                )
                continue
            genotypes.append(parse_genotype(values["GT"], len(alleles), phase_set, quality))
    record = VcfRecord(
        contig=contig,
        start=start,
        alleles=trimmed,
        qual=qual,
        genotypes=tuple(genotypes),
        line_number=line_number,
        ref_mismatch=ref_mismatch,
        span_start=span[0],
        span_end=span[1],
    )
    return record, num_not_numbers, is_beyond_reach


# QUAL and GQ only annotate a row of the table (§6), where an absent value is written `.`;
# a value that is not a number is taken as absent and counted, rather than stopping the
# run (a provisional choice). Only decimal numbers count, not Python's "nan", "inf" or
# "1_000", and a QUAL must be finite.
def _read_qual(field: str) -> float | None:
    if not DECIMAL_NUMBER.fullmatch(field):
        return None
    qual = float(field)
    return qual if math.isfinite(qual) else None


def _read_gq(field: str | None) -> int | None:
    if field is None or not field.isascii() or not field.isdecimal():
        return None
    return int(field)


def _is_not_a_number(field: str | None, number: float | None) -> bool:
    """A QUAL or GQ was written, and not as `.`, but taken as absent."""
    return number is None and field not in (None, ".")


def trim(start: int, alleles: Sequence[str]) -> tuple[int, tuple[str, ...]]:
    """A record's start and alleles without the bases all its VCF alleles share, at the
    right end, then at the left, keeping at least one base in each (the VCF anchor
    convention); start is 0-based.

    A `*` allele has no bases and takes no part. A record with fewer than two other
    alleles is returned as it is.
    """
    sequences = _sequences(alleles)
    if len(sequences) < 2:
        return start, tuple(alleles)
    start, trimmed = _trim(start, sequences)
    return start, _put_back(alleles, trimmed)


def vcf_record_span(start: int, alleles: Sequence[str], contig_sequence: str) -> Span | None:
    """The span of a trimmed record whose REF agrees with the reference (§5.2), or None
    when it would reach more than MAX_SPAN_REACH bases beyond the REF.

    start: where REF starts, 0-based. contig_sequence: the record's contig, in upper case.

    The span is the union of the REF and of the span of each ALT. An ALT that is an
    insertion or a deletion spans from the anchor of its left-aligned form to the last
    base it deletes, or the base after the point it is inserted at, when slid as far
    right as it goes. A SNP, a complex ALT and a `*` add nothing to the REF.
    """
    end = start + len(alleles[0])
    span_start, span_end = start, end
    for alt in alleles[1:]:
        if alt == SPANNING_DELETION:
            continue
        indel_span = _indel_span(start, alleles[0], alt, contig_sequence)
        if indel_span is not None:
            span_start = min(span_start, indel_span[0])
            span_end = max(span_end, indel_span[1])
    if start - span_start > MAX_SPAN_REACH or span_end - end > MAX_SPAN_REACH:
        return None
    return span_start, span_end


def _indel_span(start: int, ref: str, alt: str, contig_sequence: str) -> Span | None:
    """The span of one ALT when it is an insertion or a deletion, 0-based and half-open;
    None for a SNP or a complex ALT.

    The ALT is an indel when, without the bases it shares with REF at both ends, one of
    the two is empty, whichever side the shared bases were on.
    """
    num_prefix = _num_shared(ref, alt)
    deleted, inserted = ref[num_prefix:], alt[num_prefix:]
    num_suffix = _num_shared(deleted[::-1], inserted[::-1])
    deleted = deleted[: len(deleted) - num_suffix]
    inserted = inserted[: len(inserted) - num_suffix]
    if bool(deleted) == bool(inserted):
        return None
    event = deleted or inserted
    # The event is the deletion of [position, after), or an insertion before position,
    # where after is position too.
    position = start + num_prefix
    after = position + len(deleted)
    # A slide this long reaches more than MAX_SPAN_REACH beyond the REF, whatever the REF,
    # so sliding stops there and the caller cuts the span; a satellite of megabases is
    # not walked for each record in it.
    max_slide = MAX_SPAN_REACH + len(ref) + 1
    # The event slides left while the base before it equals its last base, and right
    # while the base after it equals its first; each slide turns the event round by one
    # base, so its bases are read round from the end it slid from.
    num_left = 0
    while (
        num_left < max_slide
        and position - num_left > 0
        and contig_sequence[position - num_left - 1] == event[(-1 - num_left) % len(event)]
    ):
        num_left += 1
    num_right = 0
    while (
        num_right < max_slide
        and after + num_right < len(contig_sequence)
        and contig_sequence[after + num_right] == event[num_right % len(event)]
    ):
        num_right += 1
    # The anchor is the base before the left-aligned event; an event slid to the first
    # base of the contig has none, and its span starts there. An insertion's span holds
    # the base after its point (§5.2: so that C -> CT at 10 overlaps a record at 11),
    # unless it is inserted after the last base of the contig.
    span_start = max(position - num_left - 1, 0)
    span_end = after + num_right + (1 if inserted else 0)
    return span_start, min(span_end, len(contig_sequence))


def _num_shared(first: str, second: str) -> int:
    """How many bases two sequences share at their start."""
    for index, (base, other_base) in enumerate(zip(first, second, strict=False)):
        if base != other_base:
            return index
    return min(len(first), len(second))


def _trim(start: int, alleles: Sequence[str]) -> tuple[int, list[str]]:
    trimmed = list(alleles)
    while min(map(len, trimmed)) > 1 and len({allele[-1] for allele in trimmed}) == 1:
        trimmed = [allele[:-1] for allele in trimmed]
    while min(map(len, trimmed)) > 1 and len({allele[0] for allele in trimmed}) == 1:
        trimmed = [allele[1:] for allele in trimmed]
        start += 1
    return start, trimmed


def _sequences(alleles: Sequence[str]) -> list[str]:
    """The alleles that are sequences. A `*` has no bases and takes no part in trimming;
    _put_back returns it to its place."""
    return [allele for allele in alleles if allele != SPANNING_DELETION]


def _put_back(alleles: Sequence[str], sequences: list[str]) -> tuple[str, ...]:
    """The alleles with the sequences _sequences took out replaced, in order, by these."""
    replacements = iter(sequences)
    return tuple(
        allele if allele == SPANNING_DELETION else next(replacements) for allele in alleles
    )
