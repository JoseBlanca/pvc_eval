"""The re-spelling program of doc/specs/respell.md: what every re-spelling does (§2), the
decomposition (§3), at respell_by_decomposing, the shift (§4), at respell_by_shifting, and
the spanning deletions (§5), at respell_with_spanning_deletions, with the expected records as literals in the notation of the cases, which
tests.respell_cases describes; each case's output against its input through the oracle of §8
and through bcftools; and the properties of §9.2 over the files of tests.respell_generator.
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
from typing import cast

import pytest
from hypothesis import given, settings
from tests.conftest import CONTIG, FILEFORMAT_LINE, CaseFiles, make_chrom_line, write_case_vcf
from tests.haplotype_oracle import check_same_haplotypes
from tests.respell import (
    _find_reach,
    make_bgzf,
    respell_by_decomposing,
    respell_by_shifting,
    respell_with_spanning_deletions,
)
from tests.respell_cases import (
    RespellingCase,
    read_case_lines,
    read_header,
    read_vcf_lines,
)
from tests.respell_generator import (
    RESPELLING_CASE_PARAMS,
    RESPELLING_PARAMS,
    RESPELLINGS,
    GeneratedVcf,
    Respelling,
    check_generated_respellings,
    check_respelled_haplotypes,
    is_rewritten,
    make_generated_vcf,
    write_and_respell,
    write_and_respell_case,
    write_generated_vcf,
)

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import VcfGenotype, parse_vcf

_PS_INTEGER = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
_GT_STRING = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
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


_DECOMPOSITION, _SHIFT, _SPANNING_DELETIONS = (
    next(respelling for respelling in RESPELLINGS if respelling.name == name)
    for name in ("decomposition", "shift", "spanning deletions")
)


def _respell(
    case_files: CaseFiles,
    directory: Path,
    respelling: Respelling,
    bases: str,
    case_lines: Sequence[str],
    *,
    ploidy: int = 2,
) -> list[str]:
    """The records a re-spelling writes, with no seed, of a case on a reference that is N but
    for bases from position 1."""
    case = RespellingCase(
        name="a case of this module",
        bases=bases,
        input_lines=tuple(case_lines),
        output_lines=(),
        ploidy=ploidy,
    )
    _, _, output = write_and_respell_case(case_files, directory, respelling, case)
    return read_case_lines(output)


@pytest.mark.parametrize(("respelling", "case"), RESPELLING_CASE_PARAMS)
def test_the_respelling_gives_the_case_of_the_spec(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling, case: RespellingCase
) -> None:
    _, _, output = write_and_respell_case(case_files, tmp_path, respelling, case)
    assert read_case_lines(output) == list(case.output_lines)


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
    assert _respell(case_files, tmp_path, _DECOMPOSITION, "GACGTA", case_lines) == expected


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
    assert (
        _respell(case_files, tmp_path, _DECOMPOSITION, "GACGTA", case_lines, ploidy=ploidy)
        == expected
    )


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
    assert _respell(case_files, tmp_path, _DECOMPOSITION, bases, case_lines) == expected


# §2.2, last paragraph: a record at the contig's last base is decomposed as any other.
def test_a_record_at_the_last_base_of_its_contig_is_decomposed(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_contigs({CONTIG: "ACGT"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["3 GT G 0/1"])
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf") == ["3 GT G 0|1:1"]


@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_a_ref_that_differs_from_the_reference_stops_before_any_output(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1", "4 T C 0/1"])
    output = tmp_path / "b.vcf"
    with pytest.raises(PvcEvalInputError, match="line 4: the REF of chr1:4 is 'T'"):
        respelling.respell(vcf, reference, output, 2, None)
    assert not output.exists()


@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_an_output_that_is_the_input_raises(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    text = vcf.read_text()
    (tmp_path / "sub").mkdir()
    with pytest.raises(ValueError, match="is the input"):
        respelling.respell(vcf, reference, tmp_path / "sub" / ".." / "a.vcf", 2, None)
    assert vcf.read_text() == text


@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_an_output_that_is_a_hard_link_to_the_input_raises(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling
) -> None:
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    text = vcf.read_text()
    output = tmp_path / "b.vcf"
    output.hardlink_to(vcf)
    with pytest.raises(ValueError, match="is the input"):
        respelling.respell(vcf, reference, output, 2, None)
    assert vcf.read_text() == text


@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_an_output_named_as_the_input_in_other_case_raises_where_case_is_ignored(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling
) -> None:
    # macOS's file system, APFS as it comes, ignores case, so A.VCF there is a.vcf; on a file
    # system that does not, A.VCF is another file and is written.
    reference = case_files.write_reference(1, "GACGTA")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    text = vcf.read_text()
    output = tmp_path / "A.VCF"
    if not output.exists():
        respelling.respell(vcf, reference, output, 2, None)
        respelling.respell(vcf, reference, tmp_path / "b.vcf", 2, None)
        assert output.read_bytes() == (tmp_path / "b.vcf").read_bytes()
        return
    with pytest.raises(ValueError, match="is the input"):
        respelling.respell(vcf, reference, output, 2, None)
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
    assert _respell(case_files, tmp_path, _DECOMPOSITION, "GACGTA", case_lines) == expected


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
        # Python turns at most 4,300 digits into a number.
        pytest.param(
            ["2 A T 0|1:" + "1" * 4301, "5 T G 0/1"],
            ["2 A T 0|1:" + "1" * 4301, "5 T G 0|1:2"],
            id="a-ps-of-4301-digits-counts-as-none",
        ),
        pytest.param(
            ["2 A T 0|1:" + "0" * 10 + "1" * 4300, "5 T G 0/1"],
            ["2 A T 0|1:" + "0" * 10 + "1" * 4300, "5 T G 0|1:" + "1" * 4299 + "3"],
            id="a-ps-of-4300-digits-after-leading-zeros-counts",
        ),
    ],
)
def test_a_new_phase_set_is_one_no_record_of_the_input_uses(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str], expected: Sequence[str]
) -> None:
    assert _respell(case_files, tmp_path, _DECOMPOSITION, "GACGTA", case_lines) == expected


@pytest.mark.parametrize(
    ("case_lines", "header_lines", "expected_header"),
    [
        pytest.param(
            ["2 A T 0/1", "5 T G 0/1"],
            [FILEFORMAT_LINE, "##contig=<ID=chr1>"],
            [FILEFORMAT_LINE, "##contig=<ID=chr1>", _GT_STRING, _PS_INTEGER],
            id="added-once",
        ),
        pytest.param(
            ["2 A T 0/1", "5 T G 0/1"],
            [FILEFORMAT_LINE, _GT_STRING, _PS_INTEGER, "##contig=<ID=chr1>"],
            [FILEFORMAT_LINE, _GT_STRING, _PS_INTEGER, "##contig=<ID=chr1>"],
            id="not-added-when-present",
        ),
        # §2.3: a rewritten record writes GT, so a header with no GT line gains one.
        pytest.param(
            ["2 A T 0/1", "5 T G 0/1"],
            [FILEFORMAT_LINE, _PS_INTEGER],
            [FILEFORMAT_LINE, _PS_INTEGER, _GT_STRING],
            id="gt-added-when-ps-is-present",
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
            [FILEFORMAT_LINE, _GT_STRING],
            id="not-added-when-a-ps-is-no-number",
        ),
        # bcftools 1.24 keeps a PS of 2,147,483,647 in an Integer and sets one more to missing.
        pytest.param(
            ["2 A T 0|1:2147483645", "5 T G 0/1"],
            [FILEFORMAT_LINE],
            [FILEFORMAT_LINE, _GT_STRING, _PS_INTEGER],
            id="added-when-the-new-phase-set-is-2147483647",
        ),
        pytest.param(
            ["2 A T 0|1:2147483646", "5 T G 0/1"],
            [FILEFORMAT_LINE],
            [FILEFORMAT_LINE, _GT_STRING],
            id="not-added-when-the-new-phase-set-is-2147483648",
        ),
        pytest.param(
            ["2 A T 0|1:2147483648", "4 G <DEL> 0/1"],
            [FILEFORMAT_LINE],
            [FILEFORMAT_LINE, _GT_STRING],
            id="not-added-when-a-ps-of-the-input-is-2147483648",
        ),
    ],
)
def test_the_header_gains_the_gt_and_ps_lines_a_rewritten_record_needs(
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
        + _GT_STRING.encode()
        + b"\n"
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
    # With no sample column, no GT is written and the header gains no GT line.
    assert read_header(tmp_path / "b.vcf") == [FILEFORMAT_LINE, _PS_INTEGER, make_chrom_line(0)]
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
    compressed = make_bgzf(text)
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


# --- The shift (§4) ---


# §2.2 for the shift, on reference GCAAAATCGTTTA: each record every re-spelling keeps, and
# the records the shift keeps for itself, written as read beside the deletion of the T at 10,
# which the shift moves to 12.
_SHIFTED_NEIGHBOUR = "9 GT G 0/1"
_NEIGHBOUR_SHIFTED = "11 TT T 0/1:."


@pytest.mark.parametrize(
    "case_lines",
    [
        pytest.param(["2 C . ./."], id="no-alt"),
        pytest.param(["2 C G,* 0/2"], id="a-star-alt"),
        # The deletion of the A at 3 slides to 6, so its reach, 2-6, lies next to the `*`
        # record's REF at 7.
        pytest.param(["2 CA C 0/1", "7 T G,* 0/2"], id="reach-touches-a-star-ref"),
        pytest.param(["2 C <DEL> 0/1"], id="symbolic"),
        pytest.param(["2 CA C 0/5"], id="malformed-gt"),
        pytest.param(["2 CA C 1"], id="other-ploidy"),
        pytest.param(["2 CA C 0/1", "3 A G 0/1"], id="touched"),
        pytest.param(["2 C G 0/1"], id="a-snp"),
        pytest.param(["2 CA GT 0/1"], id="a-complex-record"),
        pytest.param(["2 CA C,CAA 1/2"], id="two-alts"),
        pytest.param(["5 AA A 0/1"], id="already-at-its-last-placement"),
    ],
)
def test_a_record_the_shift_keeps_is_written_as_read(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str]
) -> None:
    lines = [*case_lines, _SHIFTED_NEIGHBOUR]
    assert _respell(case_files, tmp_path, _SHIFT, "GCAAAATCGTTTA", lines) == [
        *case_lines,
        _NEIGHBOUR_SHIFTED,
    ]


# §4, the phase: what stops a shift and what does not, on reference GCAAAAT.
@pytest.mark.parametrize(
    ("case_lines", "ploidy", "expected"),
    [
        # The point of the insertion passes over 3 and 4, and a third position, 5, would lie
        # next to the SNP at 6.
        pytest.param(
            ["2 C CA 0/1", "6 A G 0/1"],
            2,
            ["4 A AA 0/1:.", "6 A G 0/1"],
            id="insertion-moved-short-of-a-snp",
        ),
        pytest.param(
            ["2 CA C 1|0:9", "5 A G 0|1:8"],
            2,
            ["2 CA C 1|0:9", "5 A G 0|1:8"],
            id="two-phase-sets-kept",
        ),
        pytest.param(
            ["2 CA C 1|0", "5 A G 0|1"],
            2,
            ["5 AA A 1|0:.", "5 A G 0|1"],
            id="phased-with-no-ps-on-the-other-haplotype",
        ),
        pytest.param(
            ["2 CA C 1|0", "5 A G 0|1:9"],
            2,
            ["2 CA C 1|0", "5 A G 0|1:9"],
            id="no-ps-and-a-ps-are-two-phase-sets",
        ),
        pytest.param(
            ["2 CA C 1|0:9", "5 A G 1|0:9"],
            2,
            ["2 CA C 1|0:9", "5 A G 1|0:9"],
            id="phased-on-one-haplotype-kept",
        ),
        pytest.param(["2 CA C 1", "5 A G 1"], 1, ["2 CA C 1", "5 A G 1"], id="ploidy-1-kept"),
        pytest.param(["2 CA C 1"], 1, ["5 AA A 1:."], id="ploidy-1-moved"),
        pytest.param(["2 CA C ./1"], 2, ["5 AA A ./1:."], id="missing-carried"),
        # A `.` and REF carry no edit, and a record with no ALT makes none.
        pytest.param(
            ["2 CA C 0/1", "5 A G ./0"],
            2,
            ["5 AA A 0/1:.", "5 A G ./0"],
            id="a-missing-allele-does-not-stop-it",
        ),
        pytest.param(
            ["2 CA C 0/1", "5 A . ./."],
            2,
            ["5 AA A 0/1:.", "5 A . ./."],
            id="a-record-with-no-alt-does-not-stop-it",
        ),
        # The record at 1 has a REF that starts before the deletion and two ALTs, the first
        # of which, carried, is a SNP at 7: the deletion moves two positions and not three.
        pytest.param(
            ["1 GCAAAAT GCAAAAC,G 0/1", "2 CA C 0/1"],
            2,
            ["1 GCAAAAT GCAAAAC,G 0/1", "4 AA A 0/1:."],
            id="an-edit-of-a-record-whose-ref-starts-before",
        ),
        # The same REF, whose carried ALT is now a SNP at 1, apart from the deletion at 3:
        # the REF reaches the stretch the deletion crosses, but the edit does not.
        pytest.param(
            ["1 GCAAAAT TCAAAAT,G 0/1", "2 CA C 0/1"],
            2,
            ["1 GCAAAAT TCAAAAT,G 0/1", "5 AA A 0/1:."],
            id="an-edit-apart-on-the-left-in-a-ref-that-reaches-it",
        ),
        # A SNP just past the last placement, 6: the deletion moves short of it, to 5, and
        # the insertion's point to the one between 5 and 6.
        pytest.param(
            ["2 CA C 0/1", "7 T G 0/1"],
            2,
            ["4 AA A 0/1:.", "7 T G 0/1"],
            id="deletion-short-of-a-snp-past-its-last-placement",
        ),
        pytest.param(
            ["2 C CA 0/1", "7 T G 0/1"],
            2,
            ["5 A AA 0/1:.", "7 T G 0/1"],
            id="insertion-short-of-a-snp-past-its-last-placement",
        ),
        # The SNP at 1 has a REF that ends before the deletion, sorted after the long REF of
        # the record at 1, whose carried SNP at 7 stops the deletion at 5: the search by REF
        # must go past the short REF to the long one.
        pytest.param(
            ["1 GCAAAAT GCAAAAC,G 0/1", "1 G T 0/1", "2 CA C 0/1"],
            2,
            ["1 GCAAAAT GCAAAAC,G 0/1", "1 G T 0/1", "4 AA A 0/1:."],
            id="a-long-ref-behind-a-short-one-that-does-not-reach",
        ),
        # Deletions of 3 and of 5 on the two haplotypes of one phase set: neither stops the
        # other, and both move to 6.
        pytest.param(
            ["2 CA C 1|0:9", "4 AA A 0|1:9"],
            2,
            ["5 AA A 1|0:9", "5 AA A 0|1:9"],
            id="two-deletions-on-two-haplotypes-both-moved",
        ),
    ],
)
def test_the_shift_stops_short_of_an_edit_that_may_share_the_haplotype(
    case_files: CaseFiles,
    tmp_path: Path,
    case_lines: Sequence[str],
    ploidy: int,
    expected: Sequence[str],
) -> None:
    assert _respell(case_files, tmp_path, _SHIFT, "GCAAAAT", case_lines, ploidy=ploidy) == expected


# §4, the bases: parse_vcf accepts a REF A where the reference has M or N, and the shift
# slides only over the reference's own A, C, G and T, which reconciliation builds from.
@pytest.mark.parametrize(
    ("bases", "case_lines", "expected"),
    [
        pytest.param("GCMAAAT", ["2 CA C 0/1"], ["2 CA C 0/1"], id="deleting-an-m-kept"),
        pytest.param("GCNAAAT", ["2 CA C 0/1"], ["2 CA C 0/1"], id="deleting-an-n-kept"),
        # The As slid past are the reference's own up to the N at 6: the deletion moves to 5.
        pytest.param(
            "GCAAANT",
            ["2 CA C 0/1"],
            ["4 AA A 0/1:."],
            id="slides-up-to-an-n",
        ),
        # The inserted N would slide past the reference's N at 3.
        pytest.param("GCNAT", ["2 C CN 0/1"], ["2 C CN 0/1"], id="inserting-an-n-kept"),
    ],
)
def test_the_shift_slides_only_over_plain_bases_of_the_reference(
    case_files: CaseFiles,
    tmp_path: Path,
    bases: str,
    case_lines: Sequence[str],
    expected: Sequence[str],
) -> None:
    assert _respell(case_files, tmp_path, _SHIFT, bases, case_lines) == expected


def test_the_shift_stops_only_for_edits_of_its_own_contig(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_contigs({"chr1": "GCAAAAT", "chr2": "GCAAAAT"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["chr1 2 CA C 0/1", "chr2 5 A G 0/1"], has_chrom=True)
    respell_by_shifting(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf", has_chrom=True) == [
        "chr1 5 AA A 0/1:.",
        "chr2 5 A G 0/1",
    ]


# §4, the span: with the span limit of 10,000, the deletion of an A at 3 in a run of As
# moves when its span reaches 10,000 positions past its REF and stays when the limit cuts it.
@pytest.mark.parametrize(
    ("num_as", "expected"),
    [
        pytest.param(10_001, ["10002 AA A 0/1:."], id="span-at-the-limit-moved"),
        pytest.param(10_002, ["2 CA C 0/1"], id="span-cut-kept"),
    ],
)
def test_a_span_cut_by_the_span_limit_allows_no_shift(
    case_files: CaseFiles, tmp_path: Path, num_as: int, expected: Sequence[str]
) -> None:
    reference = case_files.write_contigs({CONTIG: "GC" + "A" * num_as + "T"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 CA C 0/1"])
    respell_by_shifting(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf") == expected


def test_a_shifted_record_is_written_with_its_qual_and_its_gt_and_ps(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, "GCAAAAT")
    vcf = tmp_path / "a.vcf"
    vcf.write_text(
        f"{FILEFORMAT_LINE}\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n"
        "chr1\t2\trs1\tCa\tc\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\n",
        encoding="utf-8",
    )
    respell_by_shifting(vcf, reference, tmp_path / "b.vcf")
    assert read_vcf_lines(tmp_path / "b.vcf") == [
        FILEFORMAT_LINE,
        _GT_STRING,
        _PS_INTEGER,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0",
        "chr1\t5\t.\tAA\tA\t50\t.\t.\tGT:PS\t0/1:.",
    ]


# §4: with a seed, each shift is drawn from 1 to the largest, here 3 for the deletion and 5
# for the insertion of S2.
_SEEDED_LINES = ["2 CA C 0/1", "11 C CAC 1/1"]
_SEEDED_BASES = "GCAAAATNNTCACACAG"


def test_the_same_seed_writes_the_same_bytes(case_files: CaseFiles, tmp_path: Path) -> None:
    reference = case_files.write_reference(1, _SEEDED_BASES)
    vcf = _write_vcf(tmp_path / "a.vcf", _SEEDED_LINES)
    outputs = [tmp_path / "b1.vcf", tmp_path / "b2.vcf"]
    for output in outputs:
        respell_by_shifting(vcf, reference, output, seed=20261020)
    assert outputs[0].read_bytes() == outputs[1].read_bytes()


@pytest.mark.parametrize(
    "seed",
    [pytest.param(True, id="bool"), pytest.param(1.5, id="float"), pytest.param("7", id="str")],
)
def test_a_seed_that_is_not_a_whole_number_raises_before_any_output(
    case_files: CaseFiles, tmp_path: Path, seed: object
) -> None:
    reference = case_files.write_reference(1, _SEEDED_BASES)
    vcf = _write_vcf(tmp_path / "a.vcf", _SEEDED_LINES)
    output = tmp_path / "b.vcf"
    # A caller without a type checker can pass any of these; cast lets the test do so.
    wrong_seed = cast(int, seed)
    message = f"^the seed must be None or a whole number, not {seed!r}$"
    with pytest.raises(ValueError, match=message):
        respell_by_shifting(vcf, reference, output, seed=wrong_seed)
    assert not output.exists()


def test_seeds_write_every_placement_from_one_to_the_largest_shift(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    placements: set[str] = set()
    for seed in range(40):
        directory = tmp_path / str(seed)
        directory.mkdir()
        reference = case_files.write_reference(1, _SEEDED_BASES)
        vcf = _write_vcf(directory / "a.vcf", _SEEDED_LINES)
        respell_by_shifting(vcf, reference, directory / "b.vcf", seed=seed)
        placements.update(read_case_lines(directory / "b.vcf"))
    assert placements == {
        "3 AA A 0/1:.",
        "4 AA A 0/1:.",
        "5 AA A 0/1:.",
        "12 A ACA 1/1:.",
        "13 C CAC 1/1:.",
        "14 A ACA 1/1:.",
        "15 C CAC 1/1:.",
        "16 A ACA 1/1:.",
    }


# --- The spanning deletions (§5) ---


# §2.2 for the `*` re-spelling, on reference GCAAAATNCAGTC: each record every re-spelling
# keeps, written as read beside a record at 11 that the deletion at 9 gives a `*`.
_STAR_NEIGHBOURS = ["9 CAGT C 0|1:5", "11 G A 0/0"]
_STAR_NEIGHBOURS_WRITTEN = ["9 CAGT C 0|1:5", "11 G A,* 0|2:5"]


@pytest.mark.parametrize(
    "case_lines",
    [
        pytest.param(["2 C . ./."], id="no-alt"),
        pytest.param(["2 C G,* 0/2"], id="a-star-alt"),
        # The deletion of the A at 3 slides to 6, so its reach, 2-6, lies next to the `*`
        # record's REF at 7.
        pytest.param(["2 CA C 0/1", "7 T G,* 0/2"], id="reach-touches-a-star-ref"),
        pytest.param(["2 C <DEL> 0/1"], id="symbolic"),
        pytest.param(["2 CA C 0/5"], id="malformed-gt"),
        pytest.param(["2 CA C 1"], id="other-ploidy"),
        pytest.param(["2 CA C 0/1", "3 A G 0/1"], id="touched"),
    ],
)
def test_a_record_the_star_respelling_keeps_is_written_as_read(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str]
) -> None:
    lines = [*case_lines, *_STAR_NEIGHBOURS]
    assert _respell(case_files, tmp_path, _SPANNING_DELETIONS, "GCAAAATNCAGTC", lines) == [
        *case_lines,
        *_STAR_NEIGHBOURS_WRITTEN,
    ]


# §2.2 for the records under a deletion, each of which would get a `*` but for its row, on
# the reference of §5: CAGTC at 9-13, with G at 8 and N before it.
@pytest.mark.parametrize(
    "case_lines",
    [
        # The REF of the record at 10 lies next to that of the `*` record at 9.
        pytest.param(["8 GCAGT G 0|1:5", "9 C A,* 0/2", "10 A G 0|0:5"], id="next-to-a-star-ref"),
        pytest.param(["8 GCAGT G 0|1:5", "11 G A 0/5"], id="malformed-gt"),
        pytest.param(["8 GCAGT G 0|1:5", "11 G A 0"], id="other-ploidy"),
        # The SNP at 11, on the haplotype without the deletion, lies next to the SNP at 12,
        # which may share it.
        pytest.param(["8 GCAGT G 0|1:5", "11 G A 1|0:5", "12 T C 0/1"], id="touched"),
    ],
)
def test_a_record_under_a_deletion_that_is_kept_gets_no_star(
    case_files: CaseFiles, tmp_path: Path, case_lines: Sequence[str]
) -> None:
    assert (
        _respell(case_files, tmp_path, _SPANNING_DELETIONS, "NNNNNNNGCAGTC", case_lines)
        == case_lines
    )


# §5, which records get a `*`, on the reference of §5: CAGTC at 9-13, with N before it.
@pytest.mark.parametrize(
    ("case_lines", "ploidy", "expected"),
    [
        # The first condition: an edit that replaces positions with no bases.
        pytest.param(
            ["9 CAGT CC 0|1:5", "11 G A 0|0:5"],
            2,
            ["9 CAGT CC 0|1:5", "11 G A 0|0:5"],
            id="a-replacement-over-the-record-gives-no-star",
        ),
        # CAGT -> CT is the deletion of AG at 10-11, once T is trimmed at the right.
        pytest.param(
            ["9 CAGT CT 0|1:5", "11 G A 0|0:5"],
            2,
            ["9 CAGT CT 0|1:5", "11 G A,* 0|2:5"],
            id="a-deletion-trimmed-at-the-right",
        ),
        # The deletion of 10 lies next to the REF at 11, and deletes none of it.
        pytest.param(
            ["9 CA C 0|1:5", "11 G A 0|0:5"],
            2,
            ["9 CA C 0|1:5", "11 G A 0|0:5"],
            id="a-deletion-next-to-the-ref-gives-no-star",
        ),
        # The second condition: a `.` is left as it is, and REF on one haplotype is not REF
        # on every one.
        pytest.param(
            ["9 CAGT C 1/1", "11 G A 0/."],
            2,
            ["9 CAGT C 1/1", "11 G A 0/."],
            id="unphased-with-a-missing-allele-gives-no-star",
        ),
        # The third condition, the phase.
        pytest.param(
            ["9 CAGT C 0|1:5", "11 G A 0|0:7"],
            2,
            ["9 CAGT C 0|1:5", "11 G A,* 0|2:5"],
            id="all-ref-leaves-its-phase-set-for-the-deletions",
        ),
        pytest.param(
            ["9 CAGT C 0|1", "11 G A 0/0"],
            2,
            ["9 CAGT C 0|1", "11 G A,* 0|2:."],
            id="all-ref-joins-a-deletion-phased-with-no-ps",
        ),
        pytest.param(
            ["9 CAGT C 1/1", "11 G A 0|0:7"],
            2,
            ["9 CAGT C 1/1", "11 G A,* 2|2:7"],
            id="under-a-homozygous-deletion-the-phase-is-as-read",
        ),
        pytest.param(
            ["9 CAGT C 1/.", "11 G A 0/0"],
            2,
            ["9 CAGT C 1|.:1", "11 G A,* 2|0:1"],
            id="a-deletion-with-a-missing-allele-is-heterozygous",
        ),
        pytest.param(
            ["9 CAGT C 0/1", "10 A G 0/0", "11 G A 0/0"],
            2,
            ["9 CAGT C 0|1:1", "10 A G,* 0|2:1", "11 G A,* 0|2:1"],
            id="two-records-under-one-deletion-share-its-new-phase-set",
        ),
        pytest.param(
            ["9 CAGT C 0|1:5", "11 G A 1|0:6"],
            2,
            ["9 CAGT C 0|1:5", "11 G A 1|0:6"],
            id="two-phase-sets-kept",
        ),
        pytest.param(
            ["9 CAGT C 0|1:5", "11 G A 1/0"],
            2,
            ["9 CAGT C 0|1:5", "11 G A 1/0"],
            id="phased-and-unphased-kept",
        ),
        pytest.param(
            ["9 CAGT C 1|0:5", "11 G A 0|0:5"],
            2,
            ["9 CAGT C 1|0:5", "11 G A,* 2|0:5"],
            id="the-star-on-the-haplotype-of-the-deletion",
        ),
        # The deletions over 11 come from phase sets 5 and 6, two phase units: the record at
        # 11 cannot join both.
        pytest.param(
            ["9 CAGT C 0|1:5", "10 AG A 1|0:6", "11 G A 0/0"],
            2,
            ["9 CAGT C 0|1:5", "10 AG A 1|0:6", "11 G A 0/0"],
            id="deletions-of-two-phase-units-give-no-star",
        ),
        # An unphased heterozygous deletion that is kept, here for the SNP next to it on a
        # haplotype they may share, cannot be given a new phase set.
        pytest.param(
            ["9 CAGT C 0/1", "11 G A 0/0", "13 C T 0/1"],
            2,
            ["9 CAGT C 0/1", "11 G A 0/0", "13 C T 0/1"],
            id="a-kept-unphased-deletion-gives-no-star",
        ),
        # A phased deletion that is kept still gives one: it is not rewritten.
        pytest.param(
            ["9 CAGT C 0|1:5", "11 G A 0/0", "13 C T 0|1:5"],
            2,
            ["9 CAGT C 0|1:5", "11 G A,* 0|2:5", "13 C T 0|1:5"],
            id="a-kept-phased-deletion-gives-a-star",
        ),
        pytest.param(
            ["9 CAGT C 1 0", "11 G A 0 1"],
            1,
            ["9 CAGT C 1 0", "11 G A,* 2:. 1:."],
            id="ploidy-1-of-two-samples",
        ),
    ],
)
def test_a_record_gets_a_star_where_its_haplotype_carries_a_deletion_in_known_phase(
    case_files: CaseFiles,
    tmp_path: Path,
    case_lines: Sequence[str],
    ploidy: int,
    expected: Sequence[str],
) -> None:
    assert (
        _respell(
            case_files, tmp_path, _SPANNING_DELETIONS, "NNNNNNNNCAGTC", case_lines, ploidy=ploidy
        )
        == expected
    )


def test_a_record_given_a_star_is_written_with_its_qual_and_its_gt_and_ps(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, "NNNNNNNNCAGTC")
    vcf = tmp_path / "a.vcf"
    vcf.write_text(
        f"{FILEFORMAT_LINE}\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n"
        "chr1\t9\trs1\tCAGT\tC\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\n"
        "chr1\t11\trs2\tg\ta\t40\tPASS\tDP=12\tGT:GQ\t0/0:20\n",
        encoding="utf-8",
    )
    respell_with_spanning_deletions(vcf, reference, tmp_path / "b.vcf")
    assert read_vcf_lines(tmp_path / "b.vcf") == [
        FILEFORMAT_LINE,
        _GT_STRING,
        _PS_INTEGER,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0",
        "chr1\t9\t.\tCAGT\tC\t50\t.\t.\tGT:PS\t0|1:1",
        "chr1\t11\t.\tG\tA,*\t40\t.\t.\tGT:PS\t0|2:1",
    ]


def test_a_star_is_given_only_by_a_deletion_of_its_own_contig(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    bases = "NNNNNNNNCAGTC"
    reference = case_files.write_contigs({"chr1": bases, "chr2": bases})
    vcf = _write_vcf(
        tmp_path / "a.vcf", ["chr1 9 CAGT C 0|1:5", "chr2 11 G A 0|0:5"], has_chrom=True
    )
    respell_with_spanning_deletions(vcf, reference, tmp_path / "b.vcf")
    assert read_case_lines(tmp_path / "b.vcf", has_chrom=True) == [
        "chr1 9 CAGT C 0|1:5",
        "chr2 11 G A 0|0:5",
    ]


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


# §9.2: at least one generated file in ten has a record each re-spelling rewrites, so that a
# re-spelling that kept everything fails. As for the pairs below, the count is of fixed
# files: the files with a shifted record ranged from 22 to 77 over 30 random runs of 300 on
# 20 October 2026, around the 30 asserted, since hypothesis draws the files of one run
# alike, and the fixed run gives 40 for the shift, 238 for the decomposition and 47 for the
# spanning deletions; the properties stay random. Over 2,000 random files the spanning
# deletions wrote a `*` in 180, 9%.
@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_a_tenth_of_the_generated_files_have_a_record_each_respelling_rewrites(
    tmp_path_factory: pytest.TempPathFactory, respelling: Respelling
) -> None:
    tally: Counter[str] = Counter()

    @settings(max_examples=_ORACLE_EXAMPLES, deadline=None, derandomize=True)
    @given(generated=make_generated_vcf())
    def count_rewritten(generated: GeneratedVcf) -> None:
        with ExitStack() as stack:
            directory = tmp_path_factory.mktemp("generated")
            _, vcf, output = write_and_respell(generated, respelling, directory, stack)
        tally["files"] += 1
        tally["rewritten"] += is_rewritten(vcf, output)

    count_rewritten()
    assert tally["rewritten"] * 10 >= tally["files"], tally


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


# §9.2: each property runs over generated files.
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
