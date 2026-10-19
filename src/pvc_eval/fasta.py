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
FAI_COLUMNS = ("name", "length", "offset", "bases per line", "bytes per line")


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

    Raises PvcEvalInputError when the FASTA is missing or compressed, when the `.fai` is
    missing or names no contig, and for a line of the `.fai` that is not UTF-8, is
    malformed or puts bases past the end of the FASTA.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        try:
            with path.open("rb") as fasta:
                magic = fasta.read(len(GZIP_MAGIC))
                fasta_size = os.fstat(fasta.fileno()).st_size
        except FileNotFoundError:
            raise PvcEvalInputError(
                f"{path} does not exist. Check the path given to --reference, and run again."
            ) from None
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
        num_bytes = fai_line.byte_of(end - 1) + 1 - first_byte
        fasta_bytes = os.pread(self._fd, num_bytes, first_byte)
        # fasta_bytes is shorter only when the FASTA was cut short after it was opened,
        # and then the region stays empty, which isalpha takes as not letters.
        region = bytearray()
        if len(fasta_bytes) == num_bytes:
            region = _drop_line_ends(
                fasta_bytes, start % fai_line.line_bases, end - start, fai_line
            )
        del fasta_bytes  # so that at most two copies of the region are held at once
        # A byte that is not a letter, or a line end where the .fai puts a base: the .fai
        # was made for another FASTA. isalpha is true only for ASCII letters.
        if not region.isalpha():
            raise PvcEvalInputError(
                f"{self.path}: the bytes of {contig}:{start + 1}-{end} that its .fai, line"
                f" {fai_line.line_number}, gives are not all letters. Make the .fai again with"
                f" `samtools faidx {self.path}`, check that the FASTA holds only letters"
                " between its header lines, and run again."
            )
        upper_region = region.upper()
        del region
        return upper_region.decode("ascii")

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
        fai_bytes = fai.read_bytes()
    except FileNotFoundError:
        raise PvcEvalInputError(
            f"{fasta} has no index {fai}. Make it with `samtools faidx {fasta}` and run again."
        ) from None
    fai_lines: dict[str, _FaiLine] = {}
    for line_number, line_bytes in enumerate(fai_bytes.splitlines(), start=1):
        try:
            line = line_bytes.decode("utf-8")
        except UnicodeDecodeError as err:
            raise PvcEvalInputError(
                f"{fai}, line {line_number}, is not UTF-8 text.{remedy}"
            ) from err
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
    if not fai_lines:
        raise PvcEvalInputError(f"{fai} names no contig.{remedy}")
    return fai_lines


def _parse_fai_line(line: str, line_number: int, where: str, remedy: str) -> tuple[str, _FaiLine]:
    """The contig name of a `.fai` line and where its bases are."""
    columns = line.split("\t")
    if len(columns) != len(FAI_COLUMNS):
        raise PvcEvalInputError(
            f"{where} has {len(columns)} columns where a .fai has {len(FAI_COLUMNS)}:"
            f" {', '.join(FAI_COLUMNS)}.{remedy}"
        )
    name = columns[0]
    if not name:
        raise PvcEvalInputError(f"{where} has no contig name.{remedy}")
    numbers = []
    for column, value in zip(FAI_COLUMNS[1:], columns[1:], strict=True):
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


def _drop_line_ends(
    fasta_bytes: bytes, first_column: int, num_bases: int, fai_line: _FaiLine
) -> bytearray:
    """The bases of fasta_bytes, which start at first_column of a line, without line ends.

    The bytes the .fai says hold line ends are dropped without being looked at; a line end
    elsewhere is left in the bases, for the caller to find.
    """
    region = bytearray(num_bases)
    line_end_bytes = fai_line.line_bytes - fai_line.line_bases
    with memoryview(fasta_bytes) as source, memoryview(region) as target:
        source_at = 0
        target_at = 0
        line_bases = fai_line.line_bases - first_column
        while target_at < num_bases:
            num_copied = min(line_bases, num_bases - target_at)
            target[target_at : target_at + num_copied] = source[source_at : source_at + num_copied]
            target_at += num_copied
            source_at += num_copied + line_end_bytes
            line_bases = fai_line.line_bases
    return region
