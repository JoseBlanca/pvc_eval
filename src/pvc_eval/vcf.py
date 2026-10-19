"""The VCF records of one file: stage 2 of doc/spec.md §4, read with the checks of §5.2,
normalised, given their spans and put in order of span start. Normalising trims the bases
a record's alleles share and marks a REF that differs from the reference; a VCF record is
never moved, split or decomposed.

Positions are 0-based and half-open inside pvc_eval. The 1-based POS of a VCF line is
converted here, and nowhere else.

Public: parse_vcf, which gives the VCF records; parse_vcf_samples, which gives the sample
names; VcfRecord and VcfGenotype; Span, the type of VcfRecord.span; the sentinels a
genotype is compared with, MISSING, MALFORMED_GT and OTHER_PLOIDY_GT; and the constants
DEFAULT_PLOIDY and DEFAULT_SPAN_LIMIT.
"""

import gzip
import heapq
import logging
import math
import re
from collections.abc import Generator, Iterable, Iterator, Sequence
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError, format_count
from pvc_eval.fasta import BLOCK_BASES, FastaWindow, IndexedFasta

_logger = logging.getLogger(__name__)

# doc/specs/reconciliation.md §2: the ploidy a file is read with when the user gives none,
# decided by the owner on 19 October 2026.
DEFAULT_PLOIDY = 2

MISSING = sentinel("MISSING")
"""A haplotype the VCF genotype says nothing about, written `.` (doc/spec.md §2)."""

MALFORMED_GT = sentinel("MALFORMED_GT")
"""A GT with an allele that is neither `.` nor an index of one of the VCF record's
alleles: the sample is invalid at the locus, and the run goes on (§5.4)."""

OTHER_PLOIDY_GT = sentinel("OTHER_PLOIDY_GT")
"""A GT whose alleles are each `.` or an index of one of the VCF record's alleles, not
all `.`, and whose number of alleles is not the ploidy the file is read with, such as `1`
on chrX of a male sample in a file read as diploid: the sample is invalid at the locus
with its own reason, `ploidy mismatch` (doc/specs/reconciliation.md §2). A GT of only
`.` is MISSING instead, whatever its number of alleles."""

_SPANNING_DELETION = "*"
_GZIP_MAGIC = b"\x1f\x8b"
_DECIMAL_NUMBER = re.compile(r"[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?")
_NUM_FIXED_COLUMNS = 9  # CHROM to FORMAT
# doc/spec.md §5.2: the span limit when the caller sets none. A span that would reach
# further than the limit beyond its REF, on either side, is cut to the REF, which bounds
# how long the reader holds a record and how many bases it holds of the reference.
DEFAULT_SPAN_LIMIT = 10_000
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
    genotypes: one per sample, in the order of the header; a MALFORMED_GT or OTHER_PLOIDY_GT
        keeps no GQ.
    line_number: the 1-based line of the file, for messages.
    ref_mismatch: REF, as written, differs from the reference (§5.2).
    span: the record's span (§5.2), 0-based and half-open; it always contains the REF.
    """

    contig: str
    start: int
    alleles: tuple[str, ...]
    qual: float | None
    genotypes: tuple[VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT, ...]
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

    Raises PvcEvalInputError when the file has no #CHROM line after its ## lines.
    """
    with closing(_read_numbered_lines(path)) as numbered_lines:
        return _parse_samples(numbered_lines, path)


def parse_vcf(
    path: Path,
    reference: IndexedFasta,
    *,
    span_limit: int | None = DEFAULT_SPAN_LIMIT,
    ploidy: int = DEFAULT_PLOIDY,
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
    ploidy: how many haplotypes each sample has, a whole number of 1 or more
        (doc/specs/reconciliation.md §2). A GT with another number of alleles is
        OTHER_PLOIDY_GT, and one of only `.` is MISSING on every haplotype.

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
    of the `.fai`'s contig order or out of position order, for a file with no #CHROM line
    after its ## lines, for a line that does not hold a VCF record, and for one whose
    number of sample columns differs from the number of samples the #CHROM line names. A
    blank line is skipped. A record with a symbolic allele is skipped, and the number
    skipped is logged as a warning at the end (§9). A QUAL or GQ that is not a number is
    taken as absent, and those are counted in a warning at the end too. Raises
    ValueError, when the stream starts, for a span_limit
    that is neither None nor a whole number of 1 or more, and for a ploidy that is not a
    whole number of 1 or more.
    """
    # §5.2: a limit of 0 would cut every insertion, whose span holds the base after it.
    is_whole = isinstance(span_limit, int) and not isinstance(span_limit, bool)
    if span_limit is not None and not (is_whole and span_limit >= 1):
        raise ValueError(
            f"span_limit must be a whole number of 1 or more, or None, not {span_limit!r}"
        )
    if not (isinstance(ploidy, int) and not isinstance(ploidy, bool) and ploidy >= 1):
        raise ValueError(f"ploidy must be a whole number of 1 or more, not {ploidy!r}")
    records_in_file_order = _parse_records_in_file_order(path, reference, span_limit, ploidy)
    with closing(records_in_file_order) as records:
        yield from _sort_by_span_start(records, span_limit)


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


def _parse_records_in_file_order(
    path: Path, reference: IndexedFasta, span_limit: int | None, ploidy: int
) -> Generator[tuple[int, VcfRecord]]:
    """The records of parse_vcf in the file's order, each with the 0-based
    position of its POS."""
    contig_rank = {contig: rank for rank, contig in enumerate(reference.contigs)}
    # A window of the stream's own, so that two streams over one reference do not drop
    # each other's buffer.
    reference_window = FastaWindow(reference)
    # The (contig rank, POS) of the previous record, which this one must not come before,
    # and where that record is, for the message.
    previous_order: tuple[int, int] | None = None
    previous_where = ""
    num_symbolic = 0
    num_not_numbers = 0
    num_beyond_reach = 0
    first_beyond_reach_line = 0
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
                num_symbolic += 1
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
            record, num_record_not_numbers, is_beyond_reach = _parse_record(
                fields,
                contig,
                pos,
                alts,
                line_number,
                reference_window,
                contig_length,
                span_limit,
                ploidy,
            )
            num_not_numbers += num_record_not_numbers
            if is_beyond_reach:
                num_beyond_reach += 1
                first_beyond_reach_line = first_beyond_reach_line or line_number
            yield pos, record
    if num_symbolic:
        _logger.warning(
            "%s: skipped %s with symbolic alleles (<DEL>, breakends); pvc_eval does not"
            " compare structural variants.",
            path,
            format_count(num_symbolic, "VCF record"),
        )
    if num_not_numbers:
        _logger.warning(
            "%s: took %s that could not be read as a number as absent, as `.` is, in the"
            " quality of the table.",
            path,
            format_count(num_not_numbers, "QUAL or GQ value"),
        )
    if num_beyond_reach:
        assert span_limit is not None, "only a limit cuts a span"
        if num_beyond_reach == 1:
            which, given = f"at line {first_beyond_reach_line}", "was given"
        else:
            which, given = f"the first at line {first_beyond_reach_line}", "each was given"
        _logger.warning(
            "%s: %s, %s, could slide more than %s along the reference and %s only"
            " its REF as its span; another spelling of the same event further along may"
            " fall into another locus and be compared apart.",
            path,
            format_count(num_beyond_reach, "VCF record"),
            which,
            f"{span_limit:,} {'base' if span_limit == 1 else 'bases'}",
            given,
        )


def _parse_genotype(
    gt: str, num_alleles: int, ploidy: int, phase_set: str | None, quality: int | None
) -> VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT:
    """A GT field as a VcfGenotype; MALFORMED_GT when an allele is neither `.` nor an
    index of one of the record's VCF alleles; or OTHER_PLOIDY_GT when its number of
    alleles is not the ploidy (doc/specs/reconciliation.md §2).

    A GT of only `.`, such as `.` or `./.`, says nothing about any haplotype, and is
    MISSING on every one of them, whatever its number of alleles.

    num_alleles: how many VCF alleles the record has, REF included.
    ploidy: how many haplotypes the sample has.
    phase_set, quality: the sample's PS and GQ, kept in the VcfGenotype.
    """
    parts = gt.replace("|", "/").split("/")
    alleles: list[int | MISSING] = []
    for part in parts:
        if part == ".":
            alleles.append(MISSING)
        elif _is_allele_index(part, num_alleles):
            alleles.append(int(part))
        else:
            return MALFORMED_GT
    if all(allele is MISSING for allele in alleles):
        alleles = [MISSING] * ploidy
    elif len(alleles) != ploidy:
        return OTHER_PLOIDY_GT
    return VcfGenotype(
        alleles=tuple(alleles), phased="|" in gt, phase_set=phase_set, quality=quality
    )


def _is_allele_index(part: str, num_alleles: int) -> bool:
    return part.isascii() and part.isdecimal() and int(part) < num_alleles


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
    with path.open("rb") as raw:
        is_gzipped = raw.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _parse_samples(numbered_lines: Iterator[tuple[int, str]], path: Path) -> tuple[str, ...]:
    """The sample names of the #CHROM line, read up to it and no further."""
    for line_number, line in numbered_lines:
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#CHROM"):
            return tuple(line.rstrip("\n").split("\t")[_NUM_FIXED_COLUMNS:])
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
    num_found = max(len(fields) - _NUM_FIXED_COLUMNS, 0)
    if num_found != num_samples:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has {format_count(num_found, 'sample column')},"
            f" and the #CHROM line names {format_count(num_samples, 'sample')}. A VCF record has a FORMAT column and"
            " then one column per sample; correct the line or remove it, and run again."
        )


# A line that is not a VCF record stops the run (a provisional choice: doc/spec.md §3 does
# not list it), with the column at fault.
def _parse_position(fields: list[str], path: Path, line_number: int) -> tuple[str, int]:
    if len(fields) < _NUM_FIXED_COLUMNS - 1:
        raise PvcEvalInputError(
            f"{path}, line {line_number}: the line has"
            f" {format_count(len(fields), 'tab-separated column')},"
            " and a VCF record has at least 8, CHROM to INFO. Check that its columns are"
            " separated by tabs, correct the line or remove it, and run again."
        )
    pos = fields[1]
    # VCF 4.3 writes a telomeric breakend at POS 0, so 0 is refused only once the record
    # is known not to be symbolic.
    if not (pos.isascii() and pos.isdecimal()):
        raise _make_pos_error(pos, path, line_number)
    return fields[0], int(pos)


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


def _parse_record(
    fields: list[str],
    contig: str,
    pos: int,
    alts: tuple[str, ...],
    line_number: int,
    reference_window: FastaWindow,
    contig_length: int,
    span_limit: int | None,
    ploidy: int,
) -> tuple[VcfRecord, int, bool]:
    """The VCF record of a line, normalised and with its span; how many of its QUAL and
    GQ values were taken as absent for not being numbers; and whether its span was cut
    to its REF for reaching more than span_limit beyond it.

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
    genotypes, num_gq_not_numbers = _parse_genotypes(fields, len(alleles), ploidy)
    record = VcfRecord(
        contig=contig,
        start=start,
        alleles=trimmed,
        qual=qual,
        genotypes=genotypes,
        line_number=line_number,
        ref_mismatch=ref_mismatch,
        span=span,
    )
    num_not_numbers = _is_taken_as_absent(fields[5], qual) + num_gq_not_numbers
    return record, num_not_numbers, is_beyond_reach


def _parse_genotypes(
    fields: list[str], num_alleles: int, ploidy: int
) -> tuple[tuple[VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT, ...], int]:
    """The genotypes of a line's sample columns, one per sample, and how many of their GQ
    values were taken as absent for not being numbers.

    num_alleles: how many VCF alleles the record has, REF included.
    ploidy: how many haplotypes each sample has.
    """
    if len(fields) <= _NUM_FIXED_COLUMNS:
        return (), 0
    genotypes: list[VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT] = []
    num_not_numbers = 0
    keys = fields[8].split(":")
    for sample_field in fields[_NUM_FIXED_COLUMNS:]:
        # VCF lets trailing sample fields be dropped, so a key may have no value.
        values = dict(zip(keys, sample_field.split(":"), strict=False))
        phase_set = values.get("PS")
        if phase_set == ".":
            phase_set = None
        quality = _parse_gq(values.get("GQ"))
        num_not_numbers += _is_taken_as_absent(values.get("GQ"), quality)
        # A sample with no GT is MISSING on every haplotype, as a GT of `.` is (§2).
        gt = values.get("GT", ".")
        genotypes.append(_parse_genotype(gt, num_alleles, ploidy, phase_set, quality))
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
    if field is None or not field.isascii() or not field.isdecimal():
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
