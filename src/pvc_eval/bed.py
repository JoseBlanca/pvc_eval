"""The lines of a --regions or --strata BED (doc/specs/loci.md "Regions" and "Stratum").

Positions are 0-based and half-open, as a BED gives them.

Public: BedLine, one interval of a BED; parse_bed, which yields them.
"""

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from pvc_eval.errors import PvcEvalInputError

_SKIPPED_PREFIXES = ("#", "track", "browser")  # as read_str_regions skips them


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
    """The intervals of a BED, in the order of its lines, from its first three columns and
    its fourth. Blank lines and lines that start with `#`, `track` or `browser` are skipped.

    contigs: the reference's contig names.

    Raises PvcEvalInputError, naming the file and the line and how to put it right, for a
    file that does not exist, is a directory, cannot be read or is not UTF-8 text; a line
    with fewer than three columns separated by tabs; a start or end that is not a whole
    number, an end not after its start; and a contig the reference lacks.
    """
    reference_contigs = frozenset(contigs)
    line_number = 0
    try:
        with path.open(encoding="utf-8") as bed:
            for line_number, line in enumerate(bed, start=1):
                if not line.strip() or line.startswith(_SKIPPED_PREFIXES):
                    continue
                yield _parse_bed_line(line, path, line_number, reference_contigs)
    except FileNotFoundError:
        raise PvcEvalInputError(f"{path} does not exist. Check its path and run again.") from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{path} is a directory, not a BED file. Check its path and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{path} cannot be read: permission denied. Allow reading it and run again."
        ) from None
    except UnicodeDecodeError as err:
        raise PvcEvalInputError(
            f"{path} is not UTF-8 text: a byte at or after line {line_number + 1} cannot be read"
            f" as UTF-8. If it is compressed, decompress it with `gunzip -c {path} > file.bed`"
            " and give file.bed instead; if it is text in another encoding, convert it to"
            " UTF-8; and run again."
        ) from err


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
        # Decided by the owner on 19 October 2026: the run stops, over skipping such lines
        # with a warning, also for GIAB's BED against a reference cut to one chromosome.
        raise PvcEvalInputError(
            f"{where} the contig {contig} is not in the reference. Keep only the reference's"
            " contigs in the BED, for instance with"
            " `awk -F'\\t' 'NR == FNR {keep[$1]; next} $1 in keep'"
            f" ref.fa.fai {path} > kept.bed`, where ref.fa.fai is the .fai"
            " of the --reference FASTA, and run again with kept.bed."
        )
    name = columns[3] if len(columns) > 3 and columns[3] else None
    return BedLine(contig=contig, start=start, end=end, name=name, line_number=line_number)


def _parse_position(text: str, column: str, where: str) -> int:
    if not (text.isascii() and text.isdigit()):
        raise PvcEvalInputError(
            f"{where} the {column}, {text}, is not a whole number. A BED line gives its start"
            " and its end as whole numbers of 0 or more; correct the line and run again."
        )
    return int(text)
