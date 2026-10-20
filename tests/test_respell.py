"""The re-spelling program of doc/specs/respell.md: what every re-spelling does (§2) and the
decomposition (§3), at respell_by_decomposing, with the expected records as literals in the
notation of the cases, which tests.respell_cases describes; each case's output against its
input through the oracle of §8 and through bcftools; and the properties of §9.2 over the
files of tests.respell_generator.
"""

import gzip
import re
import shutil
import struct
import subprocess
import time
from collections import Counter
from collections.abc import Sequence
from contextlib import ExitStack
from pathlib import Path

import pytest
from hypothesis import given, settings
from tests.conftest import CONTIG, FILEFORMAT_LINE, CaseFiles, make_chrom_line, write_case_vcf
from tests.haplotype_oracle import check_same_haplotypes
from tests.respell import _find_reach, _make_bgzf, respell_by_decomposing
from tests.respell_cases import (
    DECOMPOSITION_CASES,
    RespellingCase,
    read_case_lines,
    read_header,
    read_vcf_lines,
)
from tests.respell_generator import (
    RESPELLING_CASE_PARAMS,
    RESPELLING_PARAMS,
    GeneratedVcf,
    Respelling,
    check_generated_respellings,
    check_respelled_haplotypes,
    make_generated_vcf,
    write_and_respell_case,
    write_generated_vcf,
)

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import VcfGenotype, parse_vcf

_PS_INTEGER = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
# The 28 bytes that end a BGZF file, an empty block (the SAM/BAM format specification, §4.1.2).
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _write_vcf(
    path: Path,
    case_lines: Sequence[str],
    *,
    header_lines: Sequence[str] = (FILEFORMAT_LINE,),
    has_chrom: bool = False,
) -> Path:
    return write_case_vcf(
        path, case_lines, header_lines=header_lines, has_chrom=has_chrom, fills_ps=False
    )


def _decompose(
    case_files: CaseFiles,
    directory: Path,
    bases: str,
    case_lines: Sequence[str],
    *,
    ploidy: int = 2,
) -> list[str]:
    """The decomposition's records, of a case on a reference that is N but for bases from
    position 1."""
    reference = case_files.write_reference(1, bases)
    vcf = _write_vcf(directory / "a.vcf", case_lines)
    output = directory / "b.vcf"
    respell_by_decomposing(vcf, reference, output, ploidy=ploidy)
    return read_case_lines(output)


@pytest.mark.parametrize("case", [pytest.param(case, id=case.name) for case in DECOMPOSITION_CASES])
def test_decomposition_gives_the_case_of_the_spec(
    case_files: CaseFiles, tmp_path: Path, case: RespellingCase
) -> None:
    output_lines = _decompose(
        case_files, tmp_path, case.bases, case.input_lines, ploidy=case.ploidy
    )
    assert output_lines == list(case.output_lines)


# The cases of the tests of reconciliation that the D cases do not repeat.
@pytest.mark.parametrize(
    ("case_lines", "expected"),
    [
        pytest.param(
            ["2 A AGG 0/1", "4 GTA G 1/0"],
            ["2 A AGG 0|1:1", "4 GTA G 1|0:2"],
            id="a-phase-set-for-each-record",
        ),
        pytest.param(
            ["2 AC GT 0/1", "5 T G 0/1"],
            ["2 A G 0|1:1", "3 C T 0|1:1", "5 T G 0|1:2"],
            id="apart-decomposed",
        ),
        pytest.param(
            ["2 AC GT 0/1", "4 G A 0/0"],
            ["2 A G 0|1:1", "3 C T 0|1:1", "4 G A 0|0:2"],
            id="next-to-but-not-carried-decomposed",
        ),
        pytest.param(
            ["2 AC GT 0|1"], ["2 A G 0|1:.", "3 C T 0|1:."], id="phased-with-no-ps-keeps-none"
        ),
    ],
)
def test_decomposition_carries_the_phase_of_each_record(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str], expected: Sequence[str]
) -> None:
    assert _decompose(case_files, tmp_path, "GACGTA", case_lines) == expected


# §1, "may share a haplotype" and "touch", on reference GACGTA.
@pytest.mark.parametrize(
    ("case_lines", "ploidy", "expected"),
    [
        pytest.param(
            ["2 AC GT 1|0", "4 G A 0/1"],
            2,
            ["2 AC GT 1|0", "4 G A 0/1"],
            id="phased-with-no-ps-and-unphased-kept",
        ),
        pytest.param(
            ["2 AC GT 0/0 0/1", "4 G A 0/0 0/1"],
            2,
            ["2 AC GT 0/0 0/1", "4 G A 0/0 0/1"],
            id="touching-on-the-second-sample-kept",
        ),
        # A malformed GT carries no edit, so it touches nothing; its record is kept for itself.
        pytest.param(
            ["2 AC GT 0/5", "4 G A 0/1"],
            2,
            ["2 AC GT 0/5", "4 G A 0|1:2"],
            id="an-invalid-gt-touches-nothing",
        ),
        pytest.param(
            ["2 AC GT 1|0", "4 G A 0|1"],
            2,
            ["2 A G 1|0:.", "3 C T 1|0:.", "4 G A 0|1:."],
            id="phased-with-no-ps-is-one-phase-set",
        ),
        pytest.param(
            ["2 AC GT 1|0", "4 G A 0|1:5"],
            2,
            ["2 AC GT 1|0", "4 G A 0|1:5"],
            id="no-ps-and-a-ps-are-two-phase-sets",
        ),
        pytest.param(
            ["2 AC GT 1|0:5", "4 G A 1|0:5"],
            2,
            ["2 AC GT 1|0:5", "4 G A 1|0:5"],
            id="phased-on-one-haplotype-kept",
        ),
        pytest.param(
            ["2 AC GT 1|0:5", "4 G A 0/1"],
            2,
            ["2 AC GT 1|0:5", "4 G A 0/1"],
            id="phased-and-unphased-kept",
        ),
        pytest.param(["2 AC GT 1", "4 G A 1"], 1, ["2 AC GT 1", "4 G A 1"], id="ploidy-1-kept"),
        pytest.param(
            ["2 A AT 0/1", "3 C G 0/1"],
            2,
            ["2 A AT 0/1", "3 C G 0/1"],
            id="insertion-touches-the-position-after",
        ),
        pytest.param(
            ["2 A G 0/1", "2 A AT 0/1"],
            2,
            ["2 A G 0/1", "2 A AT 0/1"],
            id="insertion-touches-the-position-before",
        ),
        pytest.param(
            ["2 A AT 0/1", "2 A AG 0/1"],
            2,
            ["2 A AT 0/1", "2 A AG 0/1"],
            id="two-insertions-at-one-point-kept",
        ),
        pytest.param(
            ["2 A AT 0/1", "4 G C 0/1"],
            2,
            ["2 A AT 0|1:1", "4 G C 0|1:2"],
            id="insertion-one-position-apart-decomposed",
        ),
    ],
)
def test_records_whose_edits_touch_on_a_haplotype_they_may_share_are_kept(
    case_files: CaseFiles,
    tmp_path: Path,
    case_lines: Sequence[str],
    ploidy: int,
    expected: Sequence[str],
) -> None:
    assert _decompose(case_files, tmp_path, "GACGTA", case_lines, ploidy=ploidy) == expected


# §2.2: each record every re-spelling keeps, written as read beside a decomposed neighbour.
@pytest.mark.parametrize(
    ("bases", "case_lines", "expected"),
    [
        pytest.param(
            "GACGTA",
            ["2 A . ./.", "5 T G 0/1"],
            ["2 A . ./.", "5 T G 0|1:2"],
            id="no-alt",
        ),
        pytest.param(
            "GACGTACGTA",
            ["2 A T,* 1/2", "8 G C 0/1"],
            ["2 A T,* 1/2", "8 G C 0|1:2"],
            id="a-star-alt",
        ),
        # The insertion of an A after 2 slides to the point after 6, so its reach holds 7.
        pytest.param(
            "GCAAAATCGTA",
            ["2 C CA 0/1", "7 T G,* 0/1", "10 T A 0/1"],
            ["2 C CA 0/1", "7 T G,* 0/1", "10 T A 0|1:3"],
            id="reach-touches-a-star-ref",
        ),
        # AC -> TCT is complex and reaches its REF alone, but its piece `3 C CT` slides over
        # the Ts to the point after 6, so the piece's reach holds 7.
        pytest.param(
            "GACTTTGCGTA",
            ["2 AC TCT 0/1", "7 G A,* 0/1", "10 T A 0/1"],
            ["2 AC TCT 0/1", "7 G A,* 0/1", "10 T A 0|1:3"],
            id="a-piece-reach-touches-a-star-ref",
        ),
        # ACGTA -> TCGTC reaches 2-6, which holds the `*` at 4; its pieces, SNPs at 2 and 6,
        # do not reach it.
        pytest.param(
            "GACGTA",
            ["2 ACGTA TCGTC 0/1", "4 G T,* 0/2"],
            ["2 ACGTA TCGTC 0/1", "4 G T,* 0/2"],
            id="the-record-reach-holds-a-star-ref",
        ),
        pytest.param(
            "GACGTA",
            ["2 A G,* 0/2", "3 C T 0/1"],
            ["2 A G,* 0/2", "3 C T 0/1"],
            id="next-to-a-star-ref-on-the-left",
        ),
        pytest.param(
            "GACGTA",
            ["2 A T 0/1", "3 C G,* 0/2"],
            ["2 A T 0/1", "3 C G,* 0/2"],
            id="next-to-a-star-ref-on-the-right",
        ),
        pytest.param(
            "GACGTA",
            ["2 A <DEL> 0/1", "5 T G 0/1"],
            ["2 A <DEL> 0/1", "5 T G 0|1:2"],
            id="symbolic",
        ),
        pytest.param(
            "GACGTA",
            ["2 A T 0/5", "5 T G 0/1"],
            ["2 A T 0/5", "5 T G 0|1:2"],
            id="malformed-gt",
        ),
        pytest.param(
            "GACGTA",
            ["2 A T 1", "5 T G 0/1"],
            ["2 A T 1", "5 T G 0|1:2"],
            id="other-ploidy",
        ),
        pytest.param(
            "GACGTA",
            ["2 AC GT 0/1", "4 G A 0/1", "6 A C 0/1"],
            ["2 AC GT 0/1", "4 G A 0/1", "6 A C 0|1:3"],
            id="touched",
        ),
        # §3: an ALT that spells REF makes no piece.
        pytest.param(
            "GACGTA",
            ["2 A A 0/1", "5 T G 0/1"],
            ["2 A A 0/1", "5 T G 0|1:2"],
            id="no-alt-makes-a-piece",
        ),
    ],
)
def test_a_record_every_respelling_keeps_is_written_as_read(
    case_files: CaseFiles,
    tmp_path: Path,
    bases: str,
    case_lines: Sequence[str],
    expected: Sequence[str],
) -> None:
    assert _decompose(case_files, tmp_path, bases, case_lines) == expected


# §2.2, last paragraph: a record at the contig's last base is decomposed as any other.
def test_a_record_at_the_last_base_of_its_contig_is_decomposed(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_contigs({CONTIG: "ACGT"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["3 GT G 0/1"])
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf") == ["3 GT G 0|1:1"]


def test_a_ref_that_differs_from_the_reference_stops_before_any_output(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1", "4 T C 0/1"])
    output = tmp_path / "b.vcf"
    with pytest.raises(PvcEvalInputError, match="line 4: the REF of chr1:4 is 'T'"):
        respell_by_decomposing(vcf, reference, output)
    assert not output.exists()


def test_an_output_that_is_the_input_raises(case_files: CaseFiles, tmp_path: Path) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    text = vcf.read_text()
    (tmp_path / "sub").mkdir()
    with pytest.raises(ValueError, match="is the input"):
        respell_by_decomposing(vcf, reference, tmp_path / "sub" / ".." / "a.vcf")
    assert vcf.read_text() == text


def test_an_output_that_is_a_hard_link_to_the_input_raises(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    text = vcf.read_text()
    output = tmp_path / "b.vcf"
    output.hardlink_to(vcf)
    with pytest.raises(ValueError, match="is the input"):
        respell_by_decomposing(vcf, reference, output)
    assert vcf.read_text() == text


def test_an_output_named_as_the_input_in_other_case_raises_where_case_is_ignored(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # macOS's file system, APFS as it comes, ignores case, so A.VCF there is a.vcf; on a file
    # system that does not, A.VCF is another file and is written.
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    text = vcf.read_text()
    output = tmp_path / "A.VCF"
    if not output.exists():
        respell_by_decomposing(vcf, reference, output)
        assert read_case_lines(output) == ["2 A T 0|1:1"]
        return
    with pytest.raises(ValueError, match="is the input"):
        respell_by_decomposing(vcf, reference, output)
    assert vcf.read_text() == text


def test_a_star_on_another_contig_keeps_nothing(case_files: CaseFiles, tmp_path: Path) -> None:
    reference = case_files.write_contigs({"chr1": "GACGTA", "chr2": "GACGTA"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["chr1 2 A T 0/1", "chr2 2 A G,* 0/2"], has_chrom=True)
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf", has_chrom=True) == [
        "chr1 2 A T 0|1:1",
        "chr2 2 A G,* 0/2",
    ]


# §2.3, the order: records of one POS keep the order of the input records they come from.
@pytest.mark.parametrize(
    ("case_lines", "expected"),
    [
        pytest.param(
            ["2 A . ./.", "2 A T 0/1"], ["2 A . ./.", "2 A T 0|1:2"], id="no-alt-then-rewritten"
        ),
        pytest.param(
            ["2 A T 0/1", "2 A <DEL> 0/1"],
            ["2 A T 0|1:1", "2 A <DEL> 0/1"],
            id="rewritten-then-symbolic",
        ),
        # The deletion of the A at 2 is written anchored on the G at 1, at the POS of the SNP.
        pytest.param(
            ["1 G T 0/1 0/0", "2 AC C 0/0 0/1"],
            ["1 G T 0|1:1 0|0:1", "1 GA G 0|0:2 0|1:2"],
            id="a-piece-anchored-before-its-pos",
        ),
    ],
)
def test_records_of_one_pos_keep_the_order_of_the_input(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str], expected: Sequence[str]
) -> None:
    assert _decompose(case_files, tmp_path, "GACGTA", case_lines) == expected


def test_contigs_keep_the_order_of_the_input_and_not_of_their_names(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_contigs({"chr2": "AAAAAAAA", "chr1": "CCGTCCCC"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["chr2 3 A G 0/1", "chr1 2 CGT C 0/1"], has_chrom=True)
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf", has_chrom=True) == [
        "chr2 3 A G 0|1:1",
        "chr1 2 CGT C 0|1:2",
    ]


def test_records_are_sorted_by_the_contigs_of_the_input(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # By POS alone chr2's record would come first; the new phase sets count the records of
    # both contigs.
    reference = case_files.write_contigs({"chr1": "AAAAAAAA", "chr2": "CCGTCCCC"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["chr1 3 A G 0/1", "chr2 2 CGT C 0/1"], has_chrom=True)
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf", has_chrom=True) == [
        "chr1 3 A G 0|1:1",
        "chr2 2 CGT C 0|1:2",
    ]


# §2.3, a new phase set: M + n, M the largest whole-number PS of the input, kept and skipped
# records included, and n the record's place among the record lines.
@pytest.mark.parametrize(
    ("case_lines", "expected"),
    [
        pytest.param(
            ["2 A T 0|1:7", "5 T G 0/1"], ["2 A T 0|1:7", "5 T G 0|1:9"], id="after-a-ps-of-7"
        ),
        pytest.param(
            ["2 A T 0|1:7", "4 G <DEL> 0|1:30", "6 A C 0/1"],
            ["2 A T 0|1:7", "4 G <DEL> 0|1:30", "6 A C 0|1:33"],
            id="after-the-ps-of-a-skipped-record",
        ),
        pytest.param(
            ["2 A T 0|1:chr1_5", "5 T G 0/1"],
            ["2 A T 0|1:chr1_5", "5 T G 0|1:2"],
            id="a-ps-that-is-no-number-counts-as-none",
        ),
    ],
)
def test_a_new_phase_set_is_one_no_record_of_the_input_uses(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str], expected: Sequence[str]
) -> None:
    assert _decompose(case_files, tmp_path, "GACGTA", case_lines) == expected


@pytest.mark.parametrize(
    ("case_lines", "header_lines", "expected_header"),
    [
        pytest.param(
            ["2 A T 0/1", "5 T G 0/1"],
            [FILEFORMAT_LINE, "##contig=<ID=chr1>"],
            [FILEFORMAT_LINE, "##contig=<ID=chr1>", _PS_INTEGER],
            id="added-once",
        ),
        pytest.param(
            ["2 A T 0/1", "5 T G 0/1"],
            [FILEFORMAT_LINE, _PS_INTEGER, "##contig=<ID=chr1>"],
            [FILEFORMAT_LINE, _PS_INTEGER, "##contig=<ID=chr1>"],
            id="not-added-when-present",
        ),
        pytest.param(
            ["2 AC GT 0/1", "4 G A 0/1"],
            [FILEFORMAT_LINE],
            [FILEFORMAT_LINE],
            id="not-added-when-nothing-is-rewritten",
        ),
        pytest.param(
            ["2 A T 0|1:chr1_5", "5 T G 0/1"],
            [FILEFORMAT_LINE],
            [FILEFORMAT_LINE],
            id="not-added-when-a-ps-is-no-number",
        ),
    ],
)
def test_the_header_gains_a_ps_line_where_a_rewritten_record_needs_one(
    case_files: CaseFiles,
    tmp_path: Path,
    case_lines: Sequence[str],
    header_lines: Sequence[str],
    expected_header: Sequence[str],
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", case_lines, header_lines=header_lines)
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_header(tmp_path / "b.vcf") == [*expected_header, make_chrom_line(1)]


def test_kept_records_are_written_byte_for_byte_and_rewritten_ones_with_their_columns(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # The records at 2 and 4 touch and are kept, with their CRLF line ends made LF; the one
    # at 7 keeps its QUAL, and its ALT is written in upper case.
    reference = case_files.write_reference(1, "GACGTAG")
    vcf = tmp_path / "a.vcf"
    vcf.write_bytes(
        b"##fileformat=VCFv4.3\r\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\r\n"
        b"chr1\t2\trs1\tAC\tGT\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\r\n"
        b"chr1\t4\trs2\tG\tA\t40\tPASS\tDP=12\tGT:GQ\t0/1:20\r\n"
        b"chr1\t7\trs3\tG\tc\t35.5\tq10\tDP=9\tGT:GQ\t1/1:40\r\n"
    )
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert (tmp_path / "b.vcf").read_bytes() == (
        b"##fileformat=VCFv4.3\n"
        + _PS_INTEGER.encode()
        + b"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n"
        b"chr1\t2\trs1\tAC\tGT\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\n"
        b"chr1\t4\trs2\tG\tA\t40\tPASS\tDP=12\tGT:GQ\t0/1:20\n"
        b"chr1\t7\t.\tG\tC\t35.5\t.\t.\tGT:PS\t1|1:3\n"
    )


def test_a_file_with_no_samples_is_written_with_no_format_column(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 AC GT"])
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_vcf_lines(tmp_path / "b.vcf")[-2:] == [
        "chr1\t2\t.\tA\tG\t.\t.\t.",
        "chr1\t3\t.\tC\tT\t.\t.\t.",
    ]


def test_a_gzipped_input_gives_a_bgzipped_output_that_parse_vcf_reads(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    plain = _write_vcf(tmp_path / "a.vcf", ["2 AC GT 0/1"])
    vcf = tmp_path / "a.vcf.gz"
    vcf.write_bytes(gzip.compress(plain.read_bytes()))
    output = tmp_path / "b.vcf.gz"
    respell_by_decomposing(vcf, reference, output)
    compressed = output.read_bytes()
    # A BGZF block is a gzip member with the extra field BC (SAM/BAM specification §4.1).
    assert compressed[:4] == b"\x1f\x8b\x08\x04"
    assert compressed[12:16] == b"BC\x02\x00"
    assert compressed.endswith(_BGZF_EOF)
    piece_genotype = VcfGenotype(alleles=(0, 1), phased=True, phase_set="1", quality=None)
    assert [
        (vcf_record.start, vcf_record.alleles, vcf_record.genotypes)
        for vcf_record in parse_vcf(output, reference)
    ] == [(1, ("A", "G"), (piece_genotype,)), (2, ("C", "T"), (piece_genotype,))]


def test_bgzf_blocks_hold_at_most_65280_bytes_of_text_each() -> None:
    text = b"".join(b"chr1\t%d\t.\tA\tT\t.\t.\t.\n" % position for position in range(10_000))
    compressed = _make_bgzf(text)
    block_text_sizes = []
    offset = 0
    while offset < len(compressed):
        assert compressed[offset + 12 : offset + 16] == b"BC\x02\x00"
        (block_size_minus_1,) = struct.unpack_from("<H", compressed, offset + 16)
        (text_size,) = struct.unpack_from("<I", compressed, offset + block_size_minus_1 - 3)
        block_text_sizes.append(text_size)
        offset += block_size_minus_1 + 1
    assert offset == len(compressed)
    assert block_text_sizes == [65280, 65280, 65280, len(text) - 3 * 65280, 0]
    assert gzip.decompress(compressed) == text


# §1, the reach: the span of parse_vcf, which these spans are not cut, for the records of
# doc/spec.md §5.2 and of doc/specs/respell.md §3 and §4, each on a reference that is N but
# for its bases from ref_start_1based.
@pytest.mark.parametrize(
    ("ref_start_1based", "bases", "line"),
    [
        pytest.param(6, "CTAAAAGT", "9 A G 1|0", id="spec-5.2-moved-snp"),
        pytest.param(6, "CTAAAAGT", "10 AA A 1|0", id="spec-5.2-moved-deletion"),
        pytest.param(21, "CAGAGT", "22 A T 1|0", id="spec-5.2-case-1-snp"),
        pytest.param(21, "CAGAGT", "23 GAG G 1|0", id="spec-5.2-case-1-deletion"),
        pytest.param(20, "GCAGAGTC", "21 CAG C 1/1", id="spec-5.2-case-2-deletion"),
        pytest.param(20, "GCAGAGTC", "26 T C 1/1", id="spec-5.2-case-2-snp"),
        pytest.param(20, "GCAGAGTC", "24 AGT C 1/1", id="spec-5.2-case-2-complex"),
        pytest.param(10, "CG", "10 C CT 1/1", id="spec-5.2-case-3-insertion"),
        pytest.param(10, "CG", "11 G TA 1/1", id="spec-5.2-case-3-complex"),
        pytest.param(1, "GCAAAAT", "2 CA C 0/1", id="S1-deletion-in-a-homopolymer"),
        pytest.param(1, "TCACACAG", "2 C CAC 1/1", id="S2-insertion-rotated"),
        pytest.param(1, "GACGTA", "5 T TA 0/1", id="D1-piece"),
        pytest.param(1, "GCAAAAT", "2 CA C,CAA 1/2", id="S9-two-alts"),
        pytest.param(1, "ACGT", "1 AC C 0/1", id="D10-at-the-first-base"),
    ],
)
def test_the_reach_is_the_span_where_the_span_is_not_cut(
    case_files: CaseFiles, ref_start_1based: int, bases: str, line: str
) -> None:
    reference = case_files.write_reference(ref_start_1based, bases)
    (vcf_record,) = case_files.write_and_parse_records(reference, [line])
    assert _find_reach(reference, CONTIG, vcf_record.start, vcf_record.alleles) == vcf_record.span


# S7 and S8: the slides stop at the contig's end.
@pytest.mark.parametrize(
    "line", [pytest.param("2 CA C 0/1", id="S7"), pytest.param("2 C CA 0/1", id="S8")]
)
def test_the_reach_stops_at_the_end_of_the_contig(case_files: CaseFiles, line: str) -> None:
    reference = case_files.write_contigs({CONTIG: "GCAAAA"})
    (vcf_record,) = case_files.write_and_parse_records(reference, [line])
    assert vcf_record.span == (1, 6)
    assert _find_reach(reference, CONTIG, vcf_record.start, vcf_record.alleles) == (1, 6)


# A slide of 99 bases, longer than the 64 the reach reads at a time, to the left and to the
# right, in a homopolymer of 100 As at 3-102.
@pytest.mark.parametrize(
    "line",
    [pytest.param("101 AA A 0/1", id="left"), pytest.param("2 CA C 0/1", id="right")],
)
def test_the_reach_slides_past_one_read_of_the_reference(case_files: CaseFiles, line: str) -> None:
    reference = case_files.write_contigs({CONTIG: "GC" + "A" * 100 + "TG"})
    (vcf_record,) = case_files.write_and_parse_records(reference, [line])
    assert vcf_record.span == (1, 102)
    assert _find_reach(reference, CONTIG, vcf_record.start, vcf_record.alleles) == (1, 102)


def test_the_reach_has_no_limit_where_the_span_is_cut(case_files: CaseFiles) -> None:
    # The deletion of an A at 3 slides to 9: its span would reach 6 bases past its REF, 2-3.
    reference = case_files.write_reference(1, "GCAAAAAAAT")
    (vcf_record,) = case_files.write_and_parse_records(reference, ["2 CA C 0/1"], span_limit=2)
    assert vcf_record.span == (1, 3)
    assert _find_reach(reference, CONTIG, vcf_record.start, vcf_record.alleles) == (1, 9)


def _write_timing_vcf(directory: Path, reference: IndexedFasta, num_snps: int) -> Path:
    """A VCF of one sample, all 0/1: num_snps SNPs ten bases apart, at 5, 15, 25...; beside
    every tenth, at 6, 16..., a second SNP it touches; and before every hundredth from the
    fiftieth, at 501, 1501..., a deletion of the 40 bases after it, which touches the five
    SNPs of 505 to 535."""
    bases = reference.get_bases(CONTIG, 0, 10 * num_snps)
    lines = [FILEFORMAT_LINE, make_chrom_line(1)]

    def add_record(pos: int, ref: str, alt: str) -> None:
        lines.append(f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t0/1")

    def add_snp(pos: int) -> None:
        ref = bases[pos - 1]
        add_record(pos, ref, "A" if ref != "A" else "C")

    for index in range(num_snps):
        if index % 100 == 50:
            pos = 10 * index + 1
            add_record(pos, bases[pos - 1 : pos + 40], bases[pos - 1])
        add_snp(10 * index + 5)
        if index % 10 == 0:
            add_snp(10 * index + 6)
    vcf = directory / "timing.vcf"
    vcf.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return vcf


# §10: the touched test sweeps the records by REF. Comparing every pair of 20,000 records
# would take some 96 s (doc/specs/respell.md §10).
@pytest.mark.slow
def test_decomposing_20000_snps_with_touching_records_takes_under_10_seconds(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    num_snps = 20_000
    reference = case_files.write_contigs({CONTIG: "GATTACA" * (10 * num_snps // 7 + 1)})
    vcf = _write_timing_vcf(tmp_path, reference, num_snps)
    started = time.perf_counter()
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    elapsed = time.perf_counter() - started
    output_lines = read_case_lines(tmp_path / "b.vcf")
    # 22,200 records: of the 20,000 SNPs, 2,000 touch a second SNP and 600 more a deletion,
    # and are kept with the 2,000 second SNPs and the 200 deletions.
    assert len(output_lines) == 22_200
    assert sum(line.endswith("0/1") for line in output_lines) == 4_800
    assert elapsed < 10, f"{elapsed:.1f} s"


# --- The oracle and bcftools (§2.3, §8, §9) ---

_BCFTOOLS = shutil.which("bcftools")
_needs_bcftools = pytest.mark.skipif(
    _BCFTOOLS is None, reason="bcftools is not on the path, so its reading is not checked"
)
# A position in a message of bcftools, `chr1:5`, which differs between a file and its
# re-spelling since a re-spelling moves records (§2.3).
_MESSAGE_POSITION = re.compile(r"[^\s:]+:\d+")
_ORACLE_EXAMPLES = 300  # §9.2
_BCFTOOLS_EXAMPLES = 50  # §9.2


@pytest.mark.parametrize(("respelling", "case"), RESPELLING_CASE_PARAMS)
def test_the_respelling_of_a_case_passes_the_oracle(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling, case: RespellingCase
) -> None:
    reference, vcf, output = write_and_respell_case(case_files, tmp_path, respelling, case)
    comparison = check_respelled_haplotypes(vcf, output, reference, ploidy=case.ploidy)
    assert comparison.num_windows >= 1


def test_an_oracle_check_of_a_respelling_fails_on_a_window_over_the_cap(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # A re-spelling of a 14-base MNP as 14 unphased SNPs loses its phase: the oracle cannot
    # compare the window, with 2^13 phase assignments in B, and the check must not pass it.
    bases = "GACGTACGTACGTACGT"
    reference = case_files.write_reference(1, bases)
    ref = bases[1:15]
    alt = "".join("C" if base == "A" else "A" for base in ref)
    vcf = _write_vcf(tmp_path / "a.vcf", [f"2 {ref} {alt} 0/1"])
    output = _write_vcf(
        tmp_path / "b.vcf",
        [
            f"{pos} {ref_base} {alt_base} 0/1"
            for pos, ref_base, alt_base in zip(range(2, 16), ref, alt, strict=True)
        ],
    )
    with pytest.raises(AssertionError, match="1 window or pair over the oracle's cap"):
        check_respelled_haplotypes(vcf, output, reference, ploidy=2)


def _read_bcftools_messages(vcf: Path) -> set[str]:
    """What `bcftools view` writes to standard error for vcf, each line with its positions
    removed; it must exit with 0."""
    assert _BCFTOOLS is not None, "the test is skipped without bcftools"
    result = subprocess.run(
        [_BCFTOOLS, "view", str(vcf)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (vcf, result.stderr)
    return {_MESSAGE_POSITION.sub("<position>", line) for line in result.stderr.splitlines()}


def _check_bcftools_reads_as_input(vcf: Path, output: Path) -> None:
    """§2.3, valid: bcftools 1.24 reads output with exit status 0 and no message it does not
    also give for vcf."""
    new_messages = _read_bcftools_messages(output) - _read_bcftools_messages(vcf)
    assert not new_messages, (read_vcf_lines(vcf), read_vcf_lines(output), new_messages)


@_needs_bcftools
@pytest.mark.parametrize(("respelling", "case"), RESPELLING_CASE_PARAMS)
def test_bcftools_reads_the_respelling_of_a_case_as_its_input(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling, case: RespellingCase
) -> None:
    _, vcf, output = write_and_respell_case(case_files, tmp_path, respelling, case)
    _check_bcftools_reads_as_input(vcf, output)


@settings(max_examples=_ORACLE_EXAMPLES, deadline=None)
@given(generated=make_generated_vcf())
def test_the_gts_of_a_generated_record_with_no_alt_are_of_0_and_dot(
    generated: GeneratedVcf,
) -> None:
    # §9.2: an allele 1 at a record with no ALT would make the GT malformed.
    for line in generated.lines:
        _, _, _, alt, *sample_fields = line.split()
        if alt == ".":
            gts = [sample_field.split(":")[0] for sample_field in sample_fields]
            alleles = {allele for gt in gts for allele in re.split(r"[/|]", gt)}
            assert alleles <= {"0", "."}, line


# §9.2: at least one generated file in five has a pair of windows a block joins, so that a
# re-spelling that loses or invents a phase set between windows is likely to be caught. The
# count is of the generator alone, a file against itself, on fixed files: over random runs
# of 300 it ranged from 57 to 102 files on 20 October 2026, around the 60 asserted, as
# hypothesis varies its draws between runs, and the fixed run gives 73; the oracle property
# stays random.
def test_a_fifth_of_the_generated_files_have_a_pair_of_windows_a_block_joins(
    tmp_path_factory: pytest.TempPathFactory,
) -> None:
    tally: Counter[str] = Counter()

    @settings(max_examples=_ORACLE_EXAMPLES, deadline=None, derandomize=True)
    @given(generated=make_generated_vcf())
    def count_linked_pairs(generated: GeneratedVcf) -> None:
        with ExitStack() as stack:
            reference, vcf = write_generated_vcf(
                generated, tmp_path_factory.mktemp("generated"), stack
            )
            comparison = check_same_haplotypes(vcf, vcf, reference, ploidy=generated.ploidy)
        tally["files"] += 1
        tally["with a linked pair"] += comparison.num_linked_pairs > 0

    count_linked_pairs()
    assert tally["with a linked pair"] * 5 >= tally["files"], tally


# §9.2: each property runs over generated files, and asserts that at least one file in ten
# had a record rewritten.
@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_the_respelling_of_a_generated_file_passes_the_oracle(
    tmp_path_factory: pytest.TempPathFactory, respelling: Respelling
) -> None:
    def check(
        generated: GeneratedVcf, reference: IndexedFasta, vcf: Path, output: Path
    ) -> dict[str, int]:
        check_respelled_haplotypes(vcf, output, reference, ploidy=generated.ploidy)
        return {}

    check_generated_respellings(respelling, check, tmp_path_factory, max_examples=_ORACLE_EXAMPLES)


@_needs_bcftools
@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_bcftools_reads_the_respelling_of_a_generated_file_as_its_input(
    tmp_path_factory: pytest.TempPathFactory, respelling: Respelling
) -> None:
    def check(
        generated: GeneratedVcf, reference: IndexedFasta, vcf: Path, output: Path
    ) -> dict[str, int]:
        del generated, reference  # bcftools reads the files alone
        _check_bcftools_reads_as_input(vcf, output)
        return {}

    check_generated_respellings(
        respelling, check, tmp_path_factory, max_examples=_BCFTOOLS_EXAMPLES
    )
