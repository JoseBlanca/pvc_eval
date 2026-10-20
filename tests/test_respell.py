"""The re-spelling program of doc/specs/respell.md: what every re-spelling does (§2) and the
decomposition (§3), at respell_by_decomposing, with the expected records as literals.

A case is written as the spec writes it, one record a line, `POS REF ALT` and one sample
field per sample, on chr1, and the output's records are read back the same way: a record
written as read keeps its input's sample fields, `0/1`, and a rewritten one has `GT:PS`,
`0|1:1`.
"""

import gzip
import struct
import time
from collections.abc import Sequence
from pathlib import Path

import pytest
from tests.conftest import CONTIG, CaseFiles
from tests.respell import _compress_bgzf, _find_reach, respell_by_decomposing

from pvc_eval.errors import PvcEvalInputError
from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import VcfGenotype, parse_vcf

_FILEFORMAT = "##fileformat=VCFv4.3"
_PS_INTEGER = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">'
# The 28 bytes that end a BGZF file, an empty block (the SAM/BAM format specification, §4.1.2).
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _make_chrom_line(num_samples: int) -> str:
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if num_samples:
        columns += ["FORMAT", *(f"S{index}" for index in range(num_samples))]
    return "\t".join(columns)


def _format_case_line(case_line: str, *, has_chrom: bool) -> str:
    """A VCF line from a case's `[CHROM] POS REF ALT sample...`; FORMAT is GT:PS when a
    sample field has a PS, and GT otherwise."""
    fields = case_line.split()
    chrom = fields.pop(0) if has_chrom else CONTIG
    pos, ref, alt, *sample_fields = fields
    columns = [chrom, pos, ".", ref, alt, ".", ".", "."]
    if sample_fields:
        has_ps = any(":" in sample_field for sample_field in sample_fields)
        columns += ["GT:PS" if has_ps else "GT", *sample_fields]
    return "\t".join(columns)


def _write_vcf(
    path: Path,
    case_lines: Sequence[str],
    *,
    header_lines: Sequence[str] = (_FILEFORMAT,),
    has_chrom: bool = False,
) -> Path:
    num_samples = len(case_lines[0].split()) - (4 if has_chrom else 3)
    lines = [
        *header_lines,
        _make_chrom_line(num_samples),
        *(_format_case_line(line, has_chrom=has_chrom) for line in case_lines),
    ]
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return path


def _read_lines(path: Path) -> list[str]:
    with path.open("rb") as raw:
        is_gzipped = raw.read(2) == b"\x1f\x8b"
    text = gzip.decompress(path.read_bytes()).decode() if is_gzipped else path.read_text()
    return text.splitlines()


def _read_case_lines(path: Path, *, has_chrom: bool = False) -> list[str]:
    """The output's records in the notation of the cases."""
    case_lines = []
    for line in _read_lines(path):
        if line.startswith("#"):
            continue
        chrom, pos, _, ref, alt, *columns = line.split("\t")
        sample_fields = columns[4:]
        case_line = " ".join([pos, ref, alt, *sample_fields])
        case_lines.append(f"{chrom} {case_line}" if has_chrom else case_line)
    return case_lines


def _read_header(path: Path) -> list[str]:
    return [line for line in _read_lines(path) if line.startswith("#")]


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
    return _read_case_lines(output)


# doc/specs/respell.md §3: D1 to D12, each with its reference from position 1, its input,
# its output and its ploidy.
@pytest.mark.parametrize(
    ("bases", "case_lines", "expected", "ploidy"),
    [
        pytest.param(
            "GACGTA",
            ["2 ACGT TCGAA 0/1"],
            ["2 A T 0|1:1", "5 T A 0|1:1", "5 T TA 0|1:1"],
            2,
            id="D1-aligned-from-the-left",
        ),
        pytest.param(
            "GACGT",
            ["2 ACG TCA,TCG 1/2"],
            ["2 A T 1|1:1", "4 G A 1|0:1"],
            2,
            id="D2-a-piece-two-alts-share-is-one-record",
        ),
        pytest.param(
            "GACGTA",
            ["2 AC GT 1|0:7", "5 T G 0|1:7"],
            ["2 A G 1|0:7", "3 C T 1|0:7", "5 T G 0|1:7"],
            2,
            id="D3-pieces-stay-in-their-phase-set",
        ),
        pytest.param(
            "GACGTA",
            ["2 AC GT 0/1", "4 G A 0/1"],
            ["2 AC GT 0/1", "4 G A 0/1"],
            2,
            id="D4-next-to-and-unphased-kept",
        ),
        pytest.param(
            "GACGTA",
            ["2 AC GT 1|0:5", "4 G A 0|1:5"],
            ["2 A G 1|0:5", "3 C T 1|0:5", "4 G A 0|1:5"],
            2,
            id="D5-phased-on-different-haplotypes-decomposed",
        ),
        pytest.param(
            "GACGTA",
            ["2 AC GT 1|0:5", "4 G A 0|1:6"],
            ["2 AC GT 1|0:5", "4 G A 0|1:6"],
            2,
            id="D6-two-phase-sets-kept",
        ),
        pytest.param(
            "GACGTA",
            ["2 AC GT 0/1 0/0", "4 G A 0/0 1/1"],
            ["2 A G 0|1:1 0|0:1", "3 C T 0|1:1 0|0:1", "4 G A 0|0:2 1|1:2"],
            2,
            id="D7-touching-only-across-samples",
        ),
        pytest.param(
            "GACT", ["2 AC GT ./1"], ["2 A G .|1:1", "3 C T .|1:1"], 2, id="D8-missing-carried"
        ),
        pytest.param("GACT", ["2 AC GT 1"], ["2 A G 1:.", "3 C T 1:."], 1, id="D9-ploidy-1"),
        pytest.param(
            "ACGT", ["1 AC C 0/1"], ["1 AC C 0|1:1"], 2, id="D10-deletion-of-the-first-base"
        ),
        pytest.param(
            "ACGT", ["1 A CA 0/1"], ["1 A CA 0|1:1"], 2, id="D11-insertion-before-the-first-base"
        ),
        pytest.param(
            "TAAC",
            ["3 AC G,AA 0/1", "3 A G,AA 1/2"],
            ["3 AC G,AA 0/1", "3 A G,AA 1/2"],
            2,
            id="D12-insertion-inside-a-replaced-stretch-kept",
        ),
    ],
)
def test_decomposition_gives_the_case_of_the_spec(
    case_files: CaseFiles,
    tmp_path: Path,
    bases: str,
    case_lines: Sequence[str],
    expected: Sequence[str],
    ploidy: int,
) -> None:
    assert _decompose(case_files, tmp_path, bases, case_lines, ploidy=ploidy) == expected


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
    assert _read_case_lines(tmp_path / "b.vcf") == ["3 GT G 0|1:1"]


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
        assert _read_case_lines(output) == ["2 A T 0|1:1"]
        return
    with pytest.raises(ValueError, match="is the input"):
        respell_by_decomposing(vcf, reference, output)
    assert vcf.read_text() == text


def test_records_are_sorted_by_the_contigs_of_the_input(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # By POS alone chr2's record would come first; the new phase sets count the records of
    # both contigs.
    reference = case_files.write_contigs({"chr1": "AAAAAAAA", "chr2": "CCGTCCCC"})
    vcf = _write_vcf(tmp_path / "a.vcf", ["chr1 3 A G 0/1", "chr2 2 CGT C 0/1"], has_chrom=True)
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    assert _read_case_lines(tmp_path / "b.vcf", has_chrom=True) == [
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
            [_FILEFORMAT, "##contig=<ID=chr1>"],
            [_FILEFORMAT, "##contig=<ID=chr1>", _PS_INTEGER],
            id="added-once",
        ),
        pytest.param(
            ["2 A T 0/1", "5 T G 0/1"],
            [_FILEFORMAT, _PS_INTEGER, "##contig=<ID=chr1>"],
            [_FILEFORMAT, _PS_INTEGER, "##contig=<ID=chr1>"],
            id="not-added-when-present",
        ),
        pytest.param(
            ["2 AC GT 0/1", "4 G A 0/1"],
            [_FILEFORMAT],
            [_FILEFORMAT],
            id="not-added-when-nothing-is-rewritten",
        ),
        pytest.param(
            ["2 A T 0|1:chr1_5", "5 T G 0/1"],
            [_FILEFORMAT],
            [_FILEFORMAT],
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
    assert _read_header(tmp_path / "b.vcf") == [*expected_header, _make_chrom_line(1)]


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
    assert _read_lines(tmp_path / "b.vcf")[-2:] == [
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
    compressed = _compress_bgzf(text)
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


def test_the_reach_has_no_limit_where_the_span_is_cut(case_files: CaseFiles) -> None:
    # The deletion of an A at 3 slides to 9: its span would reach 6 bases past its REF, 2-3.
    reference = case_files.write_reference(1, "GCAAAAAAAT")
    (vcf_record,) = case_files.write_and_parse_records(reference, ["2 CA C 0/1"], span_limit=2)
    assert vcf_record.span == (1, 3)
    assert _find_reach(reference, CONTIG, vcf_record.start, vcf_record.alleles) == (1, 9)


def _write_snps(directory: Path, reference: IndexedFasta, num_snps: int) -> Path:
    bases = reference.get_bases(CONTIG, 0, 10 * num_snps)
    lines = [_FILEFORMAT, _make_chrom_line(1)]
    for index in range(num_snps):
        pos = 10 * index + 5
        ref = bases[pos - 1]
        alt = "A" if ref != "A" else "C"
        lines.append(f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t0/1")
    vcf = directory / "snps.vcf"
    vcf.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return vcf


# §10: the touched test sweeps the records by REF. Comparing every pair of 20,000 records
# would take some 96 s (doc/specs/respell.md §10).
@pytest.mark.slow
def test_decomposing_20000_snps_takes_under_10_seconds(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    num_snps = 20_000
    reference = case_files.write_contigs({CONTIG: "GATTACA" * (10 * num_snps // 7 + 1)})
    vcf = _write_snps(tmp_path, reference, num_snps)
    started = time.perf_counter()
    respell_by_decomposing(vcf, reference, tmp_path / "b.vcf")
    elapsed = time.perf_counter() - started
    assert len(_read_case_lines(tmp_path / "b.vcf")) == num_snps
    assert elapsed < 10, f"{elapsed:.1f} s"
