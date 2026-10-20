"""The cases of doc/specs/respell.md in its notation, for the tests of the re-spelling
program and of what checks it: RespellingCase, the D cases of §3 as DECOMPOSITION_CASES,
and the readers of a VCF file in the notation, read_case_lines, read_header and
read_vcf_lines. A case's input is written with write_case_vcf of tests.conftest, with
fills_ps=False.

A case is written as the spec writes it, one record a line, `POS REF ALT` and one sample
field per sample, on chr1, and an output's records are read back the same way: a record
written as read keeps its input's sample fields, `0/1`, and a rewritten one has `GT:PS`,
`0|1:1`.
"""

import gzip
from dataclasses import dataclass
from pathlib import Path

_GZIP_MAGIC = b"\x1f\x8b"


@dataclass(frozen=True, slots=True, kw_only=True)
class RespellingCase:
    """A case of the spec: its name, the reference's bases from position 1 (N after them),
    the input's lines, the output's lines and the ploidy."""

    name: str
    bases: str
    input_lines: tuple[str, ...]
    output_lines: tuple[str, ...]
    ploidy: int = 2


DECOMPOSITION_CASES = (
    RespellingCase(
        name="D1-aligned-from-the-left",
        bases="GACGTA",
        input_lines=("2 ACGT TCGAA 0/1",),
        output_lines=("2 A T 0|1:1", "5 T A 0|1:1", "5 T TA 0|1:1"),
    ),
    RespellingCase(
        name="D2-a-piece-two-alts-share-is-one-record",
        bases="GACGT",
        input_lines=("2 ACG TCA,TCG 1/2",),
        output_lines=("2 A T 1|1:1", "4 G A 1|0:1"),
    ),
    RespellingCase(
        name="D3-pieces-stay-in-their-phase-set",
        bases="GACGTA",
        input_lines=("2 AC GT 1|0:7", "5 T G 0|1:7"),
        output_lines=("2 A G 1|0:7", "3 C T 1|0:7", "5 T G 0|1:7"),
    ),
    RespellingCase(
        name="D4-next-to-and-unphased-kept",
        bases="GACGTA",
        input_lines=("2 AC GT 0/1", "4 G A 0/1"),
        output_lines=("2 AC GT 0/1", "4 G A 0/1"),
    ),
    RespellingCase(
        name="D5-phased-on-different-haplotypes-decomposed",
        bases="GACGTA",
        input_lines=("2 AC GT 1|0:5", "4 G A 0|1:5"),
        output_lines=("2 A G 1|0:5", "3 C T 1|0:5", "4 G A 0|1:5"),
    ),
    RespellingCase(
        name="D6-two-phase-sets-kept",
        bases="GACGTA",
        input_lines=("2 AC GT 1|0:5", "4 G A 0|1:6"),
        output_lines=("2 AC GT 1|0:5", "4 G A 0|1:6"),
    ),
    RespellingCase(
        name="D7-touching-only-across-samples",
        bases="GACGTA",
        input_lines=("2 AC GT 0/1 0/0", "4 G A 0/0 1/1"),
        output_lines=("2 A G 0|1:1 0|0:1", "3 C T 0|1:1 0|0:1", "4 G A 0|0:2 1|1:2"),
    ),
    RespellingCase(
        name="D8-missing-carried",
        bases="GACT",
        input_lines=("2 AC GT ./1",),
        output_lines=("2 A G .|1:1", "3 C T .|1:1"),
    ),
    RespellingCase(
        name="D9-ploidy-1",
        bases="GACT",
        input_lines=("2 AC GT 1",),
        output_lines=("2 A G 1:.", "3 C T 1:."),
        ploidy=1,
    ),
    RespellingCase(
        name="D10-deletion-of-the-first-base",
        bases="ACGT",
        input_lines=("1 AC C 0/1",),
        output_lines=("1 AC C 0|1:1",),
    ),
    RespellingCase(
        name="D11-insertion-before-the-first-base",
        bases="ACGT",
        input_lines=("1 A CA 0/1",),
        output_lines=("1 A CA 0|1:1",),
    ),
    RespellingCase(
        name="D12-insertion-inside-a-replaced-stretch-kept",
        bases="TAAC",
        input_lines=("3 AC G,AA 0/1", "3 A G,AA 1/2"),
        output_lines=("3 AC G,AA 0/1", "3 A G,AA 1/2"),
    ),
)


def read_vcf_lines(path: Path) -> list[str]:
    """The lines of a VCF file, plain or gzipped, without their ends."""
    with path.open("rb") as raw:
        is_gzipped = raw.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
    data = path.read_bytes()
    text = (gzip.decompress(data) if is_gzipped else data).decode("utf-8")
    return text.splitlines()


def read_case_lines(path: Path, *, has_chrom: bool = False) -> list[str]:
    """The records of a VCF file in the notation of the cases, with CHROM first when
    has_chrom."""
    case_lines = []
    for line in read_vcf_lines(path):
        if line.startswith("#"):
            continue
        chrom, pos, _, ref, alt, *columns = line.split("\t")
        sample_fields = columns[4:]
        case_line = " ".join([pos, ref, alt, *sample_fields])
        case_lines.append(f"{chrom} {case_line}" if has_chrom else case_line)
    return case_lines


def read_header(path: Path) -> list[str]:
    """The header lines of a VCF file, #CHROM included."""
    return [line for line in read_vcf_lines(path) if line.startswith("#")]
