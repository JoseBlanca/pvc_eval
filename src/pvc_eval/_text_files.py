"""The lines of a text file, plain, gzipped or bgzipped, for the readers of VCFs and BEDs,
and the temporary file an output is written to before it is renamed.

Public, to pvc_eval's own modules: open_numbered_lines, which gives a file's lines with
their numbers; open_temporary_output, which opens the temporary file of an output.
"""

import errno
import gzip
import inspect
import io
import zlib
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError

# The first two bytes of a gzip file, and so of a bgzipped one.
_GZIP_MAGIC = b"\x1f\x8b"
# The level an output is gzipped at (doc/specs/compare.md, "Writing the files", decided by
# the owner on 20 October 2026).
_GZIP_LEVEL = 6
# How many bytes a gzip file checked whole is decompressed by at a time.
_CHECK_CHUNK_BYTES = 1 << 20


@contextmanager
def open_numbered_lines(
    path: Path, file_kind: str, not_utf8_advice: str
) -> Generator[Generator[tuple[int, str]]]:
    """The lines of a file, plain, gzipped or bgzipped, each with its 1-based number, read
    as UTF-8 text with a byte-order mark at the start of the text dropped; the file is
    closed when the with block ends.

    file_kind: what the file is, "VCF" or "BED", for the message of a path that is a
        directory.
    not_utf8_advice: the sentences that end the message of a file that is not UTF-8 text,
        saying how to put it right.

    Raises PvcEvalInputError, naming the file, for a path that does not exist, is a
    directory or cannot be read, for a file that is not UTF-8 text once decompressed, and
    for a gzip file cut short or damaged. When reading a gzip file stops at a fault of its
    content, text that is not UTF-8 or a PvcEvalInputError raised in the with block before
    the last line, the file is checked whole, and if gzip finds it cut short or damaged the
    error says so in place of the fault, which is kept as its cause (doc/spec.md §5.2).
    """
    place = _ReadPlace()
    numbered_lines = _read_numbered_lines(path, file_kind, not_utf8_advice, place)
    try:
        yield numbered_lines
    except PvcEvalInputError as err:
        # gzip checks its data only at the end of each member, so a line read before the
        # end may be made of damaged bytes. A read that reached the end, and an error of
        # the reader's own, which closes it, have had every byte checked.
        is_suspended = inspect.getgeneratorstate(numbered_lines) == inspect.GEN_SUSPENDED
        if is_suspended and _is_damaged_gzip(path):
            # The fault is of the line last read.
            raise _make_damaged_error(path, place.line_number) from err
        raise
    finally:
        numbered_lines.close()


@dataclass(slots=True)
class _ReadPlace:
    """The number of the last line a reader gave, 0 before the first."""

    line_number: int = 0


def _read_numbered_lines(
    path: Path, file_kind: str, not_utf8_advice: str, place: _ReadPlace
) -> Generator[tuple[int, str]]:
    """The lines of open_numbered_lines, which raises what it says for the file itself,
    and for text that is not UTF-8 checks a gzip file whole."""
    # Text is decoded and decompressed in blocks, so a fault is found somewhere after the
    # last line read, and only that is known.
    try:
        with _open_text(path) as lines:
            for line_number, line in enumerate(lines, start=1):
                place.line_number = line_number
                yield line_number, line
    except FileNotFoundError:
        raise PvcEvalInputError(f"{path} does not exist. Check its path and run again.") from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{path} is a directory, not a {file_kind} file. Check its path and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{path} cannot be read: permission denied. Allow reading it and run again."
        ) from None
    except UnicodeDecodeError as err:
        # The byte is somewhere after the last line read, on no line known.
        if _is_damaged_gzip(path):
            raise _make_damaged_error(path, None) from err
        raise PvcEvalInputError(
            f"{path} is not UTF-8 text: a byte at or after line {place.line_number + 1} cannot be read"
            f" as UTF-8. {not_utf8_advice}"
        ) from err
    except (EOFError, gzip.BadGzipFile, zlib.error) as err:
        raise PvcEvalInputError(
            f"{path} is a gzip file cut short or damaged at or after line {place.line_number + 1}."
            " Compress or copy it again, and run again."
        ) from err


def _is_damaged_gzip(path: Path) -> bool:
    """Whether a file is a gzip file that gzip, decompressing it whole and checking the
    CRC-32 and the length of each member, finds cut short or damaged."""
    if not _is_gzipped(path):
        return False
    try:
        with gzip.open(path, "rb") as decompressed:
            while decompressed.read(_CHECK_CHUNK_BYTES):
                pass
    except EOFError, gzip.BadGzipFile, zlib.error:
        return True
    return False


def _make_damaged_error(path: Path, line_number: int | None) -> PvcEvalInputError:
    """The error of a gzip file whose damage was read as text before gzip found it, at the
    line of the fault the text made, or None when that is not known."""
    at_line = "" if line_number is None else f" at line {line_number}"
    return PvcEvalInputError(
        f"{path} is damaged: gzip's own check of the whole file fails (`gzip -t {path}` shows"
        f" it), so the fault pvc_eval found{at_line} may come from the damaged bytes rather than"
        " from the file's content. Copy or compress the file again from its source, and run"
        " again."
    )


def _is_gzipped(path: Path) -> bool:
    """Whether a file is gzipped, known by its first two bytes."""
    with path.open("rb") as raw:
        return raw.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC


def _open_text(path: Path) -> TextIO:
    """Open a file, plain or gzipped, known by its first two bytes, as UTF-8 text, with a
    byte-order mark at the start of the text dropped and one anywhere else kept.
    """
    # utf-8-sig drops the byte-order mark only at the start of the text; kept, it would
    # become part of the first line, a BED's first contig or a VCF's ##fileformat
    # (doc/specs/loci.md, Regions).
    if _is_gzipped(path):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open(encoding="utf-8-sig")


def open_temporary_output(path: Path, temporary: Path, option: str, *, compress: bool) -> TextIO:
    """Open temporary, the file an output is written to before it is renamed to path, as
    UTF-8 text, gzipped when compress is true.

    option: the command-line option that gave path, "--output" or "--summary", for the
    messages.

    Lines end in "\n" on every system. A gzipped file's header holds no file name and a
    time of 0, so `gunzip -N` does not restore the `.tmp` name (doc/specs/compare.md,
    "Writing the files").

    Raises PvcEvalInputError, naming option, when path is a directory or temporary cannot
    be made, with the system's words for a fault other than a missing directory, a path
    through a file, a permission denied, or a temporary that is a directory.
    """
    if path.is_dir():
        raise PvcEvalInputError(
            f"{option} {path} is a directory. Give the path of a file to write, and run again."
        )
    cannot = f"{option} {path} cannot be written:"
    try:
        if compress:
            return _open_gzipped_output(temporary)
        return temporary.open("w", encoding="utf-8", newline="\n")
    except FileNotFoundError as err:
        raise PvcEvalInputError(
            f"{cannot} its directory, {path.parent}, does not exist. Make it, or give another"
            f" {option}, and run again."
        ) from err
    except NotADirectoryError as err:
        raise PvcEvalInputError(
            f"{cannot} its path goes through a file as if it were a directory. Correct the"
            " path, and run again."
        ) from err
    except PermissionError as err:
        raise PvcEvalInputError(
            f"{cannot} permission denied to write in {path.parent}. Allow writing there, or"
            f" give another {option}, and run again."
        ) from err
    except IsADirectoryError as err:
        raise PvcEvalInputError(
            f"{cannot} {temporary}, where it is written first, is a directory. Move it away,"
            " and run again."
        ) from err
    except OSError as err:
        raise make_write_error(option, path, err) from err


def make_write_error(option: str, path: Path, err: OSError) -> PvcEvalInputError:
    """The error for a fault of the system while an output is opened, written or renamed,
    a full disk, a read-only file system or a name too long, naming option and giving the
    system's words (doc/specs/compare.md, "Writing the files")."""
    advice = "Free some space, or give another" if err.errno == errno.ENOSPC else "Give another"
    return PvcEvalInputError(
        f"{option} {path} could not be written: {err.strerror}. {advice} {option}, and run again."
    )


def _open_gzipped_output(temporary: Path) -> TextIO:
    """Open temporary to be written as gzipped UTF-8 text, with no name and no time in the
    gzip header."""
    raw = temporary.open("wb", buffering=0)
    try:
        compressed = gzip.GzipFile(
            filename="", mode="wb", compresslevel=_GZIP_LEVEL, fileobj=raw, mtime=0
        )
    except BaseException:
        raw.close()
        raise
    # GzipFile closes, when it is closed, the file it holds as myfileobj, which it sets only
    # for a file it opened itself; given this one, it would leave it open.
    compressed.myfileobj = raw
    return io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
