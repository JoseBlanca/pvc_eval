"""The lines of a --regions or --strata BED (doc/specs/loci.md "Regions" and "Stratum").

Positions are 0-based and half-open, as a BED gives them.

Public: BedLine, one interval of a BED; parse_bed, which yields them.
"""

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from pvc_eval._text_files import open_numbered_lines
from pvc_eval.errors import PvcEvalInputError

# The first word of the header lines of a BED that is not a comment; the STR finder's spec
# has read_str_regions skip the same lines.
_HEADER_WORDS = ("track", "browser")

# How to put right a BED that is not UTF-8 text, the end of the message that says so. A
# BED is read as parse_vcf reads a VCF, gzipped or not, decided by the owner on
# 19 October 2026, since GIAB's BEDs come gzipped.
_NOT_UTF8_ADVICE = (
    "If it is compressed by a program other than gzip or bgzip, decompress it; if it is text"
    " in another encoding, convert it to UTF-8; and run again."
)


@dataclass(frozen=True, slots=True, kw_only=True)
class BedLine:
    """One interval of a BED, `[start, end)`, 0-based and half-open.

    name: the fourth column, or None when the line has none or it is empty.
    line_number: the 1-based line of the file, for messages.
    """

    contig: str
    start: int
    end: int
    name: str | None
    line_number: int


def parse_bed(path: Path, contigs: Sequence[str]) -> Iterator[BedLine]:
    """The intervals of a BED, plain, gzipped or bgzipped, in the order of its lines, from
    its first three columns and its fourth. A byte-order mark at the start of the text is
    dropped. Blank lines, lines that start with `#`, and
    lines whose first word is `track` or `browser` are skipped.

    contigs: the reference's contig names.

    Raises PvcEvalInputError, naming the file and the line and how to put it right, for a
    file that does not exist, is a directory, cannot be read, is not UTF-8 text once
    decompressed, or is a gzip file cut short or damaged; a BED with no line once the
    skipped ones are left out; a line with fewer than three columns separated by tabs; a
    start or end that is not a whole number, an end not after its start; and a contig the
    reference lacks.
    """
    reference_contigs = frozenset(contigs)
    has_bed_line = False
    with open_numbered_lines(path, "BED", _NOT_UTF8_ADVICE) as numbered_lines:
        for line_number, line in numbered_lines:
            if _is_skipped(line):
                continue
            has_bed_line = True
            yield _parse_bed_line(line, path, line_number, reference_contigs)
    # Decided by the owner on 19 October 2026: an empty --regions BED would skip every
    # locus without a word.
    if not has_bed_line:
        raise PvcEvalInputError(
            f"{path} holds no BED line: it is empty, or holds only blank lines, comments, and"
            " track and browser lines. Give a BED with at least one line of contig, start and"
            " end, separated by tabs, and run again."
        )


def _is_skipped(line: str) -> bool:
    """Whether a line of a BED is blank, a comment or a header line."""
    words = line.split(maxsplit=1)
    return not words or line.startswith("#") or words[0] in _HEADER_WORDS


def _parse_bed_line(
    line: str, path: Path, line_number: int, reference_contigs: frozenset[str]
) -> BedLine:
    where = f"{path}, line {line_number}:"
    columns = line.rstrip("\r\n").split("\t")
    if len(columns) < 3 and len(line.split()) >= 3:
        raise PvcEvalInputError(
            f"{where} has its columns separated by spaces, and a BED separates them by tabs."
            f" Convert it with `awk -v OFS='\\t' '{{$1 = $1}} 1' {path} > tabs.bed` and run"
            " again with tabs.bed."
        )
    if len(columns) < 3:
        raise PvcEvalInputError(
            f"{where} has {len(columns)} column{'s' if len(columns) > 1 else ''}, and a BED"
            " line needs at least three, separated by tabs: the contig, the start and the end."
            " Correct the line or remove it, and run again."
        )
    contig, start_text, end_text = columns[:3]
    start = _parse_position(start_text, "start", where)
    end = _parse_position(end_text, "end", where)
    if end <= start:
        raise PvcEvalInputError(
            f"{where} the end, {end}, is not after the start, {start}. A BED line covers the"
            " positions from its start up to its end, 0-based and half-open, so its end is"
            " greater than its start; correct the line and run again."
        )
    if contig not in reference_contigs:
        raise PvcEvalInputError(_make_contig_message(contig, where, path, reference_contigs))
    name = columns[3] if len(columns) > 3 and columns[3] else None
    return BedLine(contig=contig, start=start, end=end, name=name, line_number=line_number)


def _parse_position(text: str, column: str, where: str) -> int:
    if not (text.isascii() and text.isdigit()):
        raise PvcEvalInputError(
            f"{where} the {column}, {text!r}, is not a whole number. A BED line gives its start"
            " and its end as whole numbers of 0 or more; correct the line and run again."
        )
    return int(text)


def _make_contig_message(
    contig: str, where: str, path: Path, reference_contigs: frozenset[str]
) -> str:
    """The message for a line on a contig the reference lacks."""
    # A BED that names chr1 as 1, or 1 as chr1: keeping only the reference's contigs
    # would keep no line, so the message says to rename them.
    if f"chr{contig}" in reference_contigs:
        named = f"chr{contig}"
        command = '$1 !~ /^(#|track|browser|chr)/ {$1 = "chr" $1} 1'
    elif contig.startswith("chr") and contig.removeprefix("chr") in reference_contigs:
        named = contig.removeprefix("chr")
        command = '$1 !~ /^(#|track|browser)/ {sub(/^chr/, "", $1)} 1'
    else:
        # Decided by the owner on 19 October 2026: the run stops, over skipping such lines
        # with a warning, also for GIAB's BED against a reference cut to one chromosome.
        return (
            f"{where} the contig {contig!r} is not in the reference. Keep only the reference's"
            " contigs in the BED, for instance with"
            " `awk -F'\\t' 'NR == FNR {keep[$1]; next} $1 in keep'"
            f" ref.fa.fai {path} > kept.bed`, where ref.fa.fai is the .fai"
            " of the --reference FASTA, and run again with kept.bed."
        )
    return (
        f"{where} the contig {contig!r} is not in the reference, which names it {named!r}."
        " Name the BED's contigs as the reference does, for instance with"
        f" `awk -F'\\t' -v OFS='\\t' '{command}' {path} > renamed.bed`, and run again with"
        " renamed.bed."
    )
