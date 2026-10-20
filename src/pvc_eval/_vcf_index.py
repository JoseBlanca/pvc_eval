"""The contigs that the indexes of a bgzipped VCF list, so that their scans for STR regions can be
planned before any record is read (doc/specs/compare.md, "Which contigs are scanned").

An index, a .tbi or a .csi, is BGZF, which gzip reads as one stream, and lists the contigs with
a record in its VCF in a header at its start; only that header is read. The layouts, all
integers int32 little-endian, as htslib's tabix and CSI formats give them and as bcftools 1.24
writes them for a VCF:

- .tbi: b"TBI\\x01", n_ref, six integers (format, col_seq, col_beg, col_end, meta, skip), l_nm,
  then l_nm bytes of NUL-terminated contig names.
- .csi: b"CSI\\x01", min_shift, depth, l_aux, then l_aux bytes of aux, then n_ref. The aux of a
  VCF's .csi holds the .tbi's header from its six integers on: l_nm at byte 24 of aux and the
  names from byte 28.

Public, to pvc_eval's own modules: parse_index_contigs, which gives the contigs of the .tbi and
the .csi beside a VCF.
"""

import gzip
import logging
import struct
import zlib
from io import BufferedIOBase
from pathlib import Path

_logger = logging.getLogger(__name__)

_TBI_MAGIC = b"TBI\x01"
_CSI_MAGIC = b"CSI\x01"
_INDEX_SUFFIXES = (".tbi", ".csi")
# The format field of the tabix header: 2 is VCF, in htslib's tbx.h (TBX_VCF).
_FORMAT_VCF = 2
# Bytes of the tabix header in a .csi's aux before the names: the six integers and l_nm.
_AUX_BEFORE_NAMES = 28


class _UnreadableIndexError(Exception):
    """An index whose header is cut short or is not that of a VCF's .tbi or .csi."""


def parse_index_contigs(vcf_path: Path) -> list[str]:
    """The contigs listed by the .tbi and the .csi beside a VCF, <vcf_path>.tbi and
    <vcf_path>.csi, each read when it exists, their contigs joined: those of the .tbi first,
    in its order, then those only the .csi lists. A VCF with neither, plain text ones among
    them, gives no contigs.

    An index that exists and cannot be read, its header cut short, or not an index of a VCF,
    gives no contigs, and a warning is logged that names it. Contigs that the reference lacks
    are not dropped here.
    """
    contigs: dict[str, None] = {}
    for suffix in _INDEX_SUFFIXES:
        index_path = vcf_path.with_name(vcf_path.name + suffix)
        try:
            index_contigs = _parse_index(index_path)
        except FileNotFoundError:
            continue
        except OSError, EOFError, zlib.error, _UnreadableIndexError, UnicodeDecodeError:
            # OSError includes gzip.BadGzipFile, a file that is not gzip; EOFError and
            # zlib.error, a gzip stream cut short or damaged.
            _logger.warning(
                "%s could not be read as the index of a VCF; the contigs of %s are scanned for"
                " STR regions as the reading reaches them, which can take longer.",
                index_path,
                vcf_path,
            )
            continue
        contigs.update(dict.fromkeys(index_contigs))
    return list(contigs)


def _parse_index(index_path: Path) -> list[str]:
    """The contigs that the header of a .tbi or a .csi lists, in its order.

    Raises _UnreadableIndexError for a header cut short or not that of a VCF's index, and the
    errors of gzip for a file that is not gzip or is cut short inside the header.
    """
    with gzip.open(index_path, "rb") as index_file:
        magic = _read_exactly(index_file, len(_TBI_MAGIC))
        if magic == _TBI_MAGIC:
            num_contigs, names = _read_tbi_header(index_file)
        elif magic == _CSI_MAGIC:
            num_contigs, names = _read_csi_header(index_file)
        else:
            raise _UnreadableIndexError
    return _parse_names(names, num_contigs)


def _read_tbi_header(index_file: BufferedIOBase) -> tuple[int, bytes]:
    """The number of contigs and the bytes of their names, from the .tbi header after the magic."""
    num_contigs, file_format, *_, names_length = struct.unpack("<8i", _read_exactly(index_file, 32))
    _check_format(file_format)
    return num_contigs, _read_exactly(index_file, names_length)


def _read_csi_header(index_file: BufferedIOBase) -> tuple[int, bytes]:
    """The number of contigs and the bytes of their names, from the .csi header after the magic,
    whose aux holds a tabix header."""
    *_, aux_length = struct.unpack("<3i", _read_exactly(index_file, 12))
    if aux_length < _AUX_BEFORE_NAMES:
        raise _UnreadableIndexError
    aux = _read_exactly(index_file, aux_length)
    file_format, *_, names_length = struct.unpack_from("<7i", aux)
    _check_format(file_format)
    names = aux[_AUX_BEFORE_NAMES : _AUX_BEFORE_NAMES + names_length]
    if len(names) != names_length:
        raise _UnreadableIndexError
    (num_contigs,) = struct.unpack("<i", _read_exactly(index_file, 4))
    return num_contigs, names


def _check_format(file_format: int) -> None:
    """Raise _UnreadableIndexError unless the tabix header's format is VCF."""
    if file_format != _FORMAT_VCF:
        raise _UnreadableIndexError


def _read_exactly(index_file: BufferedIOBase, size: int) -> bytes:
    """The next size bytes of the file; raises _UnreadableIndexError when it ends first, or
    when size, a length from the header, is negative."""
    # Checked before the read: gzip reads the rest of the file for -1 and raises ValueError
    # for less.
    if size < 0:
        raise _UnreadableIndexError
    data = index_file.read(size)
    if len(data) != size:
        raise _UnreadableIndexError
    return data


def _parse_names(names: bytes, num_contigs: int) -> list[str]:
    """The contig names of NUL-terminated UTF-8 bytes, which must be num_contigs of them.

    Raises _UnreadableIndexError for names of another number, which a last name that is not
    NUL-terminated gives too, and UnicodeDecodeError for names that are not UTF-8.
    """
    contigs = [name.decode("utf-8") for name in names.split(b"\0")[:-1]]
    if len(contigs) != num_contigs:
        raise _UnreadableIndexError
    return contigs
