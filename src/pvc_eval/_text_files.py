"""The lines of a text file, plain, gzipped or bgzipped, for the readers of VCFs and BEDs.

Public, to pvc_eval's own modules: read_numbered_lines, which yields a file's lines with
their numbers.
"""

import gzip
import zlib
from collections.abc import Generator
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError

# The first two bytes of a gzip file, and so of a bgzipped one.
_GZIP_MAGIC = b"\x1f\x8b"


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
