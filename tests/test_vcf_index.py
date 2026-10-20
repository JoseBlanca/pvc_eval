"""The contigs that the .tbi and .csi of a VCF list (doc/specs/compare.md, "Which contigs are
scanned" and "The index reader" of "How it is verified").

The index files are bcftools 1.24's, made as tests/data/README.md says: those of
records_on_b_d.vcf.gz list ctgB and ctgD, those of records_on_a_b.vcf.gz ctgA and ctgB, and the
header of both VCFs names ctgA to ctgD.
"""

import gzip
import logging
import shutil
import struct
from pathlib import Path

import pytest

from pvc_eval._vcf_index import parse_index_contigs

_DATA = Path(__file__).parent / "data" / "vcf_index"
_VCF_B_D = "records_on_b_d.vcf.gz"
_VCF_A_B = "records_on_a_b.vcf.gz"


def _copy_vcf(tmp_path: Path, indexes: list[tuple[str, str]]) -> Path:
    """Copy records_on_b_d.vcf.gz to tmp_path, with each (source, suffix) of indexes copied
    beside it as records_on_b_d.vcf.gz<suffix>."""
    vcf_path = tmp_path / _VCF_B_D
    shutil.copyfile(_DATA / _VCF_B_D, vcf_path)
    for source, suffix in indexes:
        shutil.copyfile(_DATA / source, tmp_path / f"{_VCF_B_D}{suffix}")
    return vcf_path


def _warnings(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [
        record.getMessage()
        for record in caplog.records
        if record.levelno == logging.WARNING and record.name == "pvc_eval._vcf_index"
    ]


def _unreadable_warning(index_path: Path, vcf_path: Path, option: str) -> str:
    """The warning of doc/specs/compare.md, "Which contigs are scanned", with the option of
    `bcftools index` that makes the index: -t for a .tbi, -c for a .csi."""
    return (
        f"{index_path} could not be read as the index of a VCF and is ignored; a contig of"
        f" {vcf_path} that no other index lists is scanned for STR regions only once the reading"
        f" reaches it, which can take longer. To use the index, make it again with"
        f" `bcftools index {option} {vcf_path}`."
    )


_INDEX_OPTIONS = {".tbi": "-t", ".csi": "-c"}


@pytest.mark.parametrize(
    "index_name",
    [f"{_VCF_B_D}.tbi", f"{_VCF_B_D}.csi"],
    ids=["tbi", "csi"],
)
def test_one_index_gives_the_contigs_with_records_and_not_the_header_ones(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, index_name: str
) -> None:
    suffix = index_name.removeprefix(_VCF_B_D)
    vcf_path = _copy_vcf(tmp_path, [(index_name, suffix)])
    with caplog.at_level(logging.WARNING):
        contigs = parse_index_contigs(vcf_path)
    assert contigs == ["ctgB", "ctgD"]
    assert _warnings(caplog) == []


def test_both_indexes_give_their_contigs_joined(tmp_path: Path) -> None:
    # The .csi is that of another VCF, as an index out of date would be.
    vcf_path = _copy_vcf(tmp_path, [(f"{_VCF_B_D}.tbi", ".tbi"), (f"{_VCF_A_B}.csi", ".csi")])
    assert parse_index_contigs(vcf_path) == ["ctgB", "ctgD", "ctgA"]


def test_vcf_with_no_index_gives_no_contigs(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    vcf_path = _copy_vcf(tmp_path, [])
    with caplog.at_level(logging.WARNING):
        contigs = parse_index_contigs(vcf_path)
    assert contigs == []
    assert _warnings(caplog) == []


def test_plain_text_vcf_gives_no_contigs(caplog: pytest.LogCaptureFixture) -> None:
    # records_on_b_d.vcf has no index beside it; its bgzipped copy has two.
    with caplog.at_level(logging.WARNING):
        contigs = parse_index_contigs(_DATA / "records_on_b_d.vcf")
    assert contigs == []
    assert _warnings(caplog) == []


def _decompressed(name: str) -> bytes:
    return gzip.decompress((_DATA / name).read_bytes())


def _with_int32(name: str, offset: int, value: int) -> bytes:
    """The index of tests/data/vcf_index named, with the int32 at offset of its decompressed
    bytes set to value, gzipped again."""
    index_bytes = bytearray(_decompressed(name))
    struct.pack_into("<i", index_bytes, offset, value)
    return gzip.compress(bytes(index_bytes))


# Each case: the suffix of the index, and its bytes. The decompressed .tbi of records_on_b_d is
# 214 bytes; its header is the magic, n_ref and six int32 to byte 32, l_nm, and the 10 bytes of
# names "ctgB\0ctgD\0" from byte 36. The .csi's is the magic, min_shift, depth and l_aux to byte
# 16, the 38 bytes of aux, which end with the same names, and n_ref at byte 54.
_UNREADABLE_INDEXES = {
    # The compressed .tbi is 124 bytes and the .csi 134; cut at 30, their headers are cut, and
    # their first 46 and 52 bytes hold the whole of it.
    "tbi_compressed_cut": (".tbi", (_DATA / f"{_VCF_B_D}.tbi").read_bytes()[:30]),
    "csi_compressed_cut": (".csi", (_DATA / f"{_VCF_B_D}.csi").read_bytes()[:30]),
    "tbi_cut_in_fixed_fields": (".tbi", gzip.compress(_decompressed(f"{_VCF_B_D}.tbi")[:30])),
    "tbi_cut_in_names": (".tbi", gzip.compress(_decompressed(f"{_VCF_B_D}.tbi")[:41])),
    "csi_cut_in_aux": (".csi", gzip.compress(_decompressed(f"{_VCF_B_D}.csi")[:45])),
    "csi_cut_before_n_ref": (".csi", gzip.compress(_decompressed(f"{_VCF_B_D}.csi")[:54])),
    # The format field at byte 8, from 2, VCF, to 0, generic (htslib's tbx.h).
    "tbi_not_of_a_vcf": (".tbi", _with_int32(f"{_VCF_B_D}.tbi", 8, 0)),
    # n_ref at byte 4, from 2 to 3, more than the names.
    "tbi_n_ref_not_its_names": (".tbi", _with_int32(f"{_VCF_B_D}.tbi", 4, 3)),
    # l_nm at byte 32, from 10 to -1, -2 and the lowest int32: a gzip read of -1 bytes reads the
    # rest of the file, and one of less raises ValueError.
    "tbi_names_length_minus_1": (".tbi", _with_int32(f"{_VCF_B_D}.tbi", 32, -1)),
    "tbi_names_length_minus_2": (".tbi", _with_int32(f"{_VCF_B_D}.tbi", 32, -2)),
    "tbi_names_length_lowest": (".tbi", _with_int32(f"{_VCF_B_D}.tbi", 32, -(2**31))),
    # The first byte of the names, at byte 36, from "c" to 0xff, which is not UTF-8.
    "tbi_names_not_utf8": (
        ".tbi",
        gzip.compress(
            _decompressed(f"{_VCF_B_D}.tbi")[:36] + b"\xff" + _decompressed(f"{_VCF_B_D}.tbi")[37:]
        ),
    ),
    # l_aux at byte 12, from 38 to 0, as a .csi with no tabix header has.
    "csi_without_tabix_header": (".csi", _with_int32(f"{_VCF_B_D}.csi", 12, 0)),
    # l_nm at byte 24 of aux, from 10 to 11, one past the end of aux.
    "csi_names_past_aux": (".csi", _with_int32(f"{_VCF_B_D}.csi", 16 + 24, 11)),
    # A valid .tbi with its magic alone changed, to a BAM index's.
    "tbi_magic_of_bai": (
        ".tbi",
        gzip.compress(b"BAI\x01" + _decompressed(f"{_VCF_B_D}.tbi")[4:]),
    ),
    # The format field of the .csi, at byte 0 of aux and 16 of the file, from 2, VCF, to 0.
    "csi_not_of_a_vcf": (".csi", _with_int32(f"{_VCF_B_D}.csi", 16, 0)),
    # l_aux at byte 12, from 38 to 27, one short of the tabix header before the names.
    "csi_aux_27_bytes": (".csi", _with_int32(f"{_VCF_B_D}.csi", 12, 27)),
    # l_nm from 10 to 9, and the NUL that ends the names removed: b"ctgB\0ctgD".
    "tbi_last_name_not_nul_terminated": (
        ".tbi",
        gzip.compress(
            _decompressed(f"{_VCF_B_D}.tbi")[:32]
            + struct.pack("<i", 9)
            + _decompressed(f"{_VCF_B_D}.tbi")[36:45]
            + _decompressed(f"{_VCF_B_D}.tbi")[46:]
        ),
    ),
    # The compressed bytes 18 to 25, the start of the deflate data after the BGZF block's
    # 18-byte gzip header, set to 0xff.
    "tbi_deflate_damaged": (
        ".tbi",
        (_DATA / f"{_VCF_B_D}.tbi").read_bytes()[:18]
        + b"\xff" * 8
        + (_DATA / f"{_VCF_B_D}.tbi").read_bytes()[26:],
    ),
    "wrong_magic": (".tbi", (_DATA / _VCF_B_D).read_bytes()),
    "not_gzip": (".tbi", b"TBI\x01 is not compressed\n"),
    "empty": (".csi", b""),
}


@pytest.mark.parametrize(
    ("suffix", "index_bytes"),
    list(_UNREADABLE_INDEXES.values()),
    ids=list(_UNREADABLE_INDEXES),
)
def test_index_that_cannot_be_read_gives_the_warning_and_no_contigs(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, suffix: str, index_bytes: bytes
) -> None:
    vcf_path = _copy_vcf(tmp_path, [])
    index_path = tmp_path / f"{_VCF_B_D}{suffix}"
    index_path.write_bytes(index_bytes)
    with caplog.at_level(logging.WARNING):
        contigs = parse_index_contigs(vcf_path)
    assert contigs == []
    assert _warnings(caplog) == [_unreadable_warning(index_path, vcf_path, _INDEX_OPTIONS[suffix])]


def test_index_that_cannot_be_read_leaves_the_other_ones_contigs(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    vcf_path = _copy_vcf(tmp_path, [(f"{_VCF_A_B}.csi", ".csi")])
    tbi_path = tmp_path / f"{_VCF_B_D}.tbi"
    tbi_path.write_bytes(b"not an index")
    with caplog.at_level(logging.WARNING):
        contigs = parse_index_contigs(vcf_path)
    assert contigs == ["ctgA", "ctgB"]
    assert _warnings(caplog) == [_unreadable_warning(tbi_path, vcf_path, "-t")]


# Only the header is read, so an index cut after it is read whole: the bins that follow are cut
# away. The .tbi's header is 36 bytes and the 10 of names, the .csi's 16 bytes, the 38 of aux and
# the 4 of n_ref.
_INDEXES_CUT_AFTER_HEADER = {
    "tbi_decompressed": (".tbi", gzip.compress(_decompressed(f"{_VCF_B_D}.tbi")[:46])),
    "csi_decompressed": (".csi", gzip.compress(_decompressed(f"{_VCF_B_D}.csi")[:58])),
    "tbi_compressed": (".tbi", (_DATA / f"{_VCF_B_D}.tbi").read_bytes()[:60]),
    "csi_compressed": (".csi", (_DATA / f"{_VCF_B_D}.csi").read_bytes()[:60]),
}


@pytest.mark.parametrize(
    ("suffix", "index_bytes"),
    list(_INDEXES_CUT_AFTER_HEADER.values()),
    ids=list(_INDEXES_CUT_AFTER_HEADER),
)
def test_index_cut_after_its_header_still_gives_its_contigs(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, suffix: str, index_bytes: bytes
) -> None:
    vcf_path = _copy_vcf(tmp_path, [])
    (tmp_path / f"{_VCF_B_D}{suffix}").write_bytes(index_bytes)
    with caplog.at_level(logging.WARNING):
        contigs = parse_index_contigs(vcf_path)
    assert contigs == ["ctgB", "ctgD"]
    assert _warnings(caplog) == []
