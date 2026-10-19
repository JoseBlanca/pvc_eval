"""A FASTA read through its `.fai`, one region at a time. The reference (doc/spec.md §3) is
read this way.

Only the bytes of the region asked for are read from disk, so the FASTA is never held in
memory. Positions are 0-based and half-open, as in the `.fai`.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import Self

from pvc_eval.errors import PvcEvalInputError

GZIP_MAGIC = b"\x1f\x8b"  # the first two bytes of a gzip or bgzip file
NUM_FAI_COLUMNS = 5  # name, length, offset, bases per line, bytes per line
LINE_END_BYTES = b"\r\n"


@dataclass(frozen=True, slots=True, kw_only=True)
class _FaiLine:
    """Where the bases of one contig are in the FASTA, as one line of its `.fai` says.

    offset: the byte of the FASTA that holds the first base.
    line_number: the 1-based line of the `.fai`, for messages.
    """

    length: int
    offset: int
    line_bases: int
    line_bytes: int
    line_number: int

    def byte_of(self, position: int) -> int:
        """The byte of the FASTA that holds the base at a 0-based position."""
        line, column = divmod(position, self.line_bases)
        return self.offset + line * self.line_bytes + column


class IndexedFasta:
    """A FASTA with its `.fai` beside it (`ref.fa.fai`).

    The contig names and their order are those of the `.fai`. An IndexedFasta is a context
    manager that holds one open file descriptor and closes it on exit. Its reads do not
    depend on a file position, so several threads can read one IndexedFasta at once.

    Raises PvcEvalInputError when the FASTA is compressed, when the `.fai` is missing, and
    for a line of the `.fai` that is malformed or whose bases fall past the end of the
    FASTA; FileNotFoundError when the FASTA itself is missing.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        with path.open("rb") as fasta:
            magic = fasta.read(len(GZIP_MAGIC))
            fasta_size = os.fstat(fasta.fileno()).st_size
        if magic == GZIP_MAGIC:
            raise PvcEvalInputError(
                f"{path} is compressed with gzip or bgzip, and pvc_eval reads only an"
                f" uncompressed FASTA. Decompress it with `gunzip -c {path} > file.fa`,"
                " index it with `samtools faidx file.fa`, and give file.fa instead."
            )
        self._fai_lines = _read_fai(path, fasta_size)
        # One descriptor read with os.pread, which takes the byte to start at and so
        # leaves no file position for threads to share.
        self._fd = os.open(path, os.O_RDONLY)
        self.contigs = tuple(self._fai_lines)
        """The contig names, in the order of the `.fai`."""

    def length(self, contig: str) -> int:
        """The number of bases of a contig.

        Raises ValueError for a contig the `.fai` lacks.
        """
        return self._fai_line(contig).length

    def bases(self, contig: str, start: int, end: int) -> str:
        """The bases of a region of a contig, 0-based and half-open, in upper case.

        Raises ValueError for a contig the `.fai` lacks, a negative start, an end before
        the start or past the end of the contig. Raises PvcEvalInputError when the bytes
        of the region are not bases where the `.fai` says they are.
        """
        fai_line = self._fai_line(contig)
        if start < 0:
            raise ValueError(f"start must be 0 or more, not {start}")
        if end < start:
            raise ValueError(f"end {end} is before start {start}")
        if end > fai_line.length:
            raise ValueError(f"end {end} is past the length of {contig}, {fai_line.length}")
        if start == end:
            return ""
        first_byte = fai_line.byte_of(start)
        raw = os.pread(self._fd, fai_line.byte_of(end - 1) + 1 - first_byte, first_byte)
        # Pad the first line to a whole one, so that every line of the region starts at a
        # multiple of the bytes per line, and its line end follows its bases.
        padding = start % fai_line.line_bases
        padded = b"N" * padding + raw
        lines = (
            padded[line_start : line_start + fai_line.line_bases]
            for line_start in range(0, len(padded), fai_line.line_bytes)
        )
        region = b"".join(lines)[padding:]
        where = f"{self.path}: the bases of {contig}:{start + 1}-{end}"
        if len(region) != end - start or any(byte in region for byte in LINE_END_BYTES):
            raise PvcEvalInputError(
                f"{where} are not where its .fai, line {fai_line.line_number}, says."
                f" Make the .fai again with `samtools faidx {self.path}` and run again."
            )
        if not region.isascii():
            raise PvcEvalInputError(
                f"{where} hold a byte that is not ASCII, where a FASTA holds only letters."
                " Check that it is a FASTA file, and run again."
            )
        return region.decode("ascii").upper()

    def close(self) -> None:
        """Close the FASTA. A read after this raises OSError."""
        if self._fd >= 0:
            os.close(self._fd)
            self._fd = -1

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def _fai_line(self, contig: str) -> _FaiLine:
        try:
            return self._fai_lines[contig]
        except KeyError:
            raise ValueError(f"contig {contig} is not in {self.path}") from None


def _read_fai(fasta: Path, fasta_size: int) -> dict[str, _FaiLine]:
    """The lines of the `.fai` of a FASTA of fasta_size bytes, by contig, in its order."""
    fai = Path(f"{fasta}.fai")
    remedy = f" Make the .fai again with `samtools faidx {fasta}` and run again."
    try:
        text = fai.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise PvcEvalInputError(
            f"{fasta} has no index {fai}. Make it with `samtools faidx {fasta}` and run again."
        ) from None
    fai_lines: dict[str, _FaiLine] = {}
    for line_number, line in enumerate(text.splitlines(), start=1):
        where = f"{fai}, line {line_number}:"
        name, fai_line = _parse_fai_line(line, line_number, where, remedy)
        if name in fai_lines:
            raise PvcEvalInputError(
                f"{where} names contig {name} again, which line"
                f" {fai_lines[name].line_number} named first.{remedy}"
            )
        if fai_line.length and fai_line.byte_of(fai_line.length - 1) >= fasta_size:
            raise PvcEvalInputError(
                f"{where} puts the bases of {name} past the end of the FASTA, which has"
                f" {fasta_size} bytes.{remedy}"
            )
        fai_lines[name] = fai_line
    return fai_lines


def _parse_fai_line(line: str, line_number: int, where: str, remedy: str) -> tuple[str, _FaiLine]:
    """The contig name of a `.fai` line and where its bases are."""
    columns = line.split("\t")
    if len(columns) != NUM_FAI_COLUMNS:
        raise PvcEvalInputError(
            f"{where} has {len(columns)} columns where a .fai has {NUM_FAI_COLUMNS}:"
            f" name, length, offset, bases per line, bytes per line.{remedy}"
        )
    name = columns[0]
    if not name:
        raise PvcEvalInputError(f"{where} has no contig name.{remedy}")
    numbers = []
    for column, value in zip(
        ("length", "offset", "bases per line", "bytes per line"), columns[1:], strict=True
    ):
        if not value.isdecimal():
            raise PvcEvalInputError(
                f"{where} has {value} as its {column}, which is not a whole number of 0 or"
                f" more.{remedy}"
            )
        numbers.append(int(value))
    length, offset, line_bases, line_bytes = numbers
    if line_bases == 0:
        raise PvcEvalInputError(
            f"{where} has 0 bases per line, where a line holds at least one base.{remedy}"
        )
    if line_bytes <= line_bases:
        raise PvcEvalInputError(
            f"{where} has {line_bytes} bytes per line for {line_bases} bases per line, which"
            f" leaves no byte for the line end.{remedy}"
        )
    return name, _FaiLine(
        length=length,
        offset=offset,
        line_bases=line_bases,
        line_bytes=line_bytes,
        line_number=line_number,
    )
