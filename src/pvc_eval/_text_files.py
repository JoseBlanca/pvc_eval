"""The lines of a text file, plain, gzipped or bgzipped, for the readers of VCFs and BEDs,
and the temporary file an output is written to before it is renamed.

Public, to pvc_eval's own modules: read_numbered_lines, which yields a file's lines with
their numbers; open_temporary_output, which opens the temporary file of an output.
"""

import gzip
import zlib
from collections.abc import Generator
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError

# The first two bytes of a gzip file, and so of a bgzipped one.
_GZIP_MAGIC = b"\x1f\x8b"
# The level an output is gzipped at (doc/specs/compare.md, "Writing the files", decided by
# the owner on 20 October 2026).
_GZIP_LEVEL = 6


def read_numbered_lines(
    path: Path, file_kind: str, not_utf8_advice: str
) -> Generator[tuple[int, str]]:
    """The lines of a file, plain, gzipped or bgzipped, each with its 1-based number, read
    as UTF-8 text with a byte-order mark at the start of the text dropped.

    file_kind: what the file is, "VCF" or "BED", for the message of a path that is a
        directory.
    not_utf8_advice: the sentences that end the message of a file that is not UTF-8 text,
        saying how to put it right.

    Raises PvcEvalInputError, naming the file, for a path that does not exist, is a
    directory or cannot be read, for a file that is not UTF-8 text once decompressed, and
    for a gzip file cut short or damaged.
    """
    # Text is decoded and decompressed in blocks, so a fault is found somewhere after the
    # last line read, and only that is known.
    line_number = 0
    try:
        with _open_text(path) as lines:
            for line_number, line in enumerate(lines, start=1):
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
        raise PvcEvalInputError(
            f"{path} is not UTF-8 text: a byte at or after line {line_number + 1} cannot be read"
            f" as UTF-8. {not_utf8_advice}"
        ) from err
    except (EOFError, gzip.BadGzipFile, zlib.error) as err:
        raise PvcEvalInputError(
            f"{path} is a gzip file cut short or damaged at or after line {line_number + 1}."
            " Compress or copy it again, and run again."
        ) from err


def _open_text(path: Path) -> TextIO:
    """Open a file, plain or gzipped, known by its first two bytes, as UTF-8 text, with a
    byte-order mark at the start of the text dropped and one anywhere else kept.
    """
    with path.open("rb") as raw:
        is_gzipped = raw.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
    # utf-8-sig drops the byte-order mark only at the start of the text; kept, it would
    # become part of the first line, a BED's first contig or a VCF's ##fileformat
    # (doc/specs/loci.md, Regions).
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open(encoding="utf-8-sig")


def open_temporary_output(path: Path, temporary: Path, option: str, *, compress: bool) -> TextIO:
    """Open temporary, the file an output is written to before it is renamed to path, as
    UTF-8 text, gzipped when compress is true.

    option: the command-line option that gave path, "--output" or "--summary", for the
    messages.

    Raises PvcEvalInputError, naming option, when path is a directory or temporary cannot
    be made.
    """
    if path.is_dir():
        raise PvcEvalInputError(
            f"{option} {path} is a directory. Give the path of a file to write, and run again."
        )
    cannot = f"{option} {path} cannot be written:"
    try:
        if compress:
            return gzip.open(temporary, "wt", compresslevel=_GZIP_LEVEL, encoding="utf-8")
        return temporary.open("w", encoding="utf-8")
    except FileNotFoundError:
        raise PvcEvalInputError(
            f"{cannot} its directory, {path.parent}, does not exist. Make it, or give another"
            f" {option}, and run again."
        ) from None
    except NotADirectoryError:
        raise PvcEvalInputError(
            f"{cannot} its path goes through a file as if it were a directory. Correct the"
            " path, and run again."
        ) from None
    except PermissionError:
        raise PvcEvalInputError(
            f"{cannot} permission denied to write in {path.parent}. Allow writing there, or"
            f" give another {option}, and run again."
        ) from None
    except IsADirectoryError:
        raise PvcEvalInputError(
            f"{cannot} {temporary}, where it is written first, is a directory. Move it away,"
            " and run again."
        ) from None
