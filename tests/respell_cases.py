"""The cases of doc/specs/respell.md in its notation, for the tests of the re-spelling
program and of what checks it: RespellingCase, the D cases of §3 as DECOMPOSITION_CASES,
the S cases of §4 as SHIFT_CASES, the T cases of §5 as SPANNING_DELETION_CASES, the M cases
of §6 as MERGE_CASES, and the readers of a VCF file in the notation, read_case_lines,
read_header and read_vcf_lines. A case's input is written with
write_case_vcf of tests.conftest, with fills_ps=False.

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
    """A case of the spec: its name, the reference's bases from position 1, the input's
    lines, the output's lines and the ploidy. The bases are followed by N, but when
    is_whole_contig, where they are the whole contig, as in S7 and S8."""

    name: str
    bases: str
    input_lines: tuple[str, ...]
    output_lines: tuple[str, ...]
    ploidy: int = 2
    is_whole_contig: bool = False


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

SHIFT_CASES = (
    RespellingCase(
        name="S1-deletion-moved-to-the-last-a",
        bases="GCAAAAT",
        input_lines=("2 CA C 0/1",),
        output_lines=("5 AA A 0/1:.",),
    ),
    RespellingCase(
        name="S2-insertion-rotated-as-it-slides",
        bases="TCACACAG",
        input_lines=("2 C CAC 1/1",),
        output_lines=("7 A ACA 1/1:.",),
    ),
    RespellingCase(
        name="S3-next-to-a-snp-that-may-share-its-haplotype-kept",
        bases="GCAAAAT",
        input_lines=("2 CA C 0/1", "5 A G 0/1"),
        output_lines=("2 CA C 0/1", "5 A G 0/1"),
    ),
    RespellingCase(
        name="S4-moved-short-of-a-snp",
        bases="GCAAAAAT",
        input_lines=("2 CA C 0/1", "7 A G 0/1"),
        output_lines=("4 AA A 0/1:.", "7 A G 0/1"),
    ),
    RespellingCase(
        name="S5-snp-on-the-other-haplotype",
        bases="GCAAAAT",
        input_lines=("2 CA C 1|0:9", "5 A G 0|1:9"),
        output_lines=("5 AA A 1|0:9", "5 A G 0|1:9"),
    ),
    RespellingCase(
        name="S6-snp-of-another-sample",
        bases="GCAAAAT",
        input_lines=("2 CA C 0/1 0/0", "5 A G 0/0 0/1"),
        output_lines=("5 AA A 0/1:. 0/0:.", "5 A G 0/0 0/1"),
    ),
    RespellingCase(
        name="S7-deletion-of-the-last-base-of-the-contig",
        bases="GCAAAA",
        input_lines=("2 CA C 0/1",),
        output_lines=("5 AA A 0/1:.",),
        is_whole_contig=True,
    ),
    RespellingCase(
        name="S8-insertion-after-the-last-base-of-the-contig",
        bases="GCAAAA",
        input_lines=("2 C CA 0/1",),
        output_lines=("6 A AA 0/1:.",),
        is_whole_contig=True,
    ),
    RespellingCase(
        name="S9-two-alts-kept",
        bases="GCAAAAT",
        input_lines=("2 CA C,CAA 1/2",),
        output_lines=("2 CA C,CAA 1/2",),
    ),
)

# §5: the reference of doc/spec.md §5.4, CAGTC at 9-13 with N before it, so that `9 CAGT C`
# deletes 10-12 and slides nowhere.
_SPANNING_DELETION_BASES = "NNNNNNNNCAGTC"

SPANNING_DELETION_CASES = (
    RespellingCase(
        name="T1-phased-in-one-phase-set",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 0|1:5", "11 G A 1|0:5"),
        output_lines=("9 CAGT C 0|1:5", "11 G A,* 1|2:5"),
    ),
    RespellingCase(
        name="T2-homozygous-deletion-over-a-record-of-another-sample",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 1/1 0/0", "11 G A 0/0 0/1"),
        output_lines=("9 CAGT C 1/1 0/0", "11 G A,* 2/2:. 0/1:."),
    ),
    RespellingCase(
        name="T3-all-ref-shares-the-new-phase-set-of-the-deletion",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 0/1 0/0", "11 G A 0/0 0/1"),
        output_lines=("9 CAGT C 0|1:1 0/0:.", "11 G A,* 0|2:1 0/1:."),
    ),
    RespellingCase(
        name="T4-all-ref-joins-the-phase-set-of-the-deletion",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 0|1:5", "11 G A 0/0"),
        output_lines=("9 CAGT C 0|1:5", "11 G A,* 0|2:5"),
    ),
    RespellingCase(
        name="T5-ploidy-1",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 1", "11 G A 0"),
        output_lines=("9 CAGT C 1", "11 G A,* 2:."),
        ploidy=1,
    ),
    RespellingCase(
        name="T6-missing-under-the-deletion-kept",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 1|0:5", "11 G A .|1:5"),
        output_lines=("9 CAGT C 1|0:5", "11 G A .|1:5"),
    ),
    RespellingCase(
        name="T7-unphased-snp-under-the-deletion-kept",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 0/1", "11 G A 0/1"),
        output_lines=("9 CAGT C 0/1", "11 G A 0/1"),
    ),
    RespellingCase(
        name="T8-a-star-already-there-kept",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("9 CAGT C 0/1", "11 G A,* 1/2"),
        output_lines=("9 CAGT C 0/1", "11 G A,* 1/2"),
    ),
    # §5, the second way: the deletions over 11-13, of 11 and of 13, come from phase sets 5
    # and 6, and the record joins neither, though the first way alone would give a `*` on
    # the first haplotype.
    RespellingCase(
        name="all-ref-under-two-phase-units-kept",
        bases="NNNNNNNNCAGTCAGTCA",
        input_lines=("10 AG A 1|0:5", "11 GTC ACA 0|0:5", "12 TC T 0|1:6"),
        output_lines=("10 AG A 1|0:5", "11 GTC ACA 0|0:5", "12 TC T 0|1:6"),
    ),
    # §5, the first condition: R's REF as parse_vcf trims it. `12 TC TA` is `13 C A`, which a
    # deletion of 12 alone does not delete, so its haplotype gets no `*`.
    RespellingCase(
        name="trimmed-ref-past-the-deletion-kept",
        bases=_SPANNING_DELETION_BASES,
        input_lines=("11 GT G 0|1:5", "12 TC TA 1|0:5"),
        output_lines=("11 GT G 0|1:5", "12 TC TA 1|0:5"),
    ),
)


MERGE_CASES = (
    RespellingCase(
        name="M1-one-position-between",
        bases="GACGTA",
        input_lines=("2 A T 1/1", "4 G C 0/1"),
        output_lines=("2 ACG TCG,TCC 1/2:.",),
    ),
    RespellingCase(
        name="M2-phased-in-one-phase-set",
        bases="GACGTA",
        input_lines=("2 A T 1|0:5", "4 G C 0|1:5"),
        output_lines=("2 ACG TCG,ACC 1|2:5",),
    ),
    RespellingCase(
        name="M3-two-unphased-heterozygous-records-kept",
        bases="GACGTA",
        input_lines=("2 A T 0/1", "4 G C 0/1"),
        output_lines=("2 A T 0/1", "4 G C 0/1"),
    ),
    RespellingCase(
        name="M4-two-phase-sets-kept",
        bases="GACGTA",
        input_lines=("2 A T 1|0:5", "4 G C 0|1:6"),
        output_lines=("2 A T 1|0:5", "4 G C 0|1:6"),
    ),
    RespellingCase(
        name="M5-an-alt-for-each-distinct-haplotype",
        bases="GACGTA",
        input_lines=("2 A T 1|0:5 1/1", "4 G C 0|1:5 1/1"),
        output_lines=("2 ACG TCG,ACC,TCC 1|2:5 3/3:.",),
    ),
    RespellingCase(
        name="M6-a-haplotype-of-ref",
        bases="GACGTA",
        input_lines=("2 A T 0|1:5", "4 G C 0|1:5"),
        output_lines=("2 ACG TCC 0|1:5",),
    ),
    RespellingCase(
        name="M7-a-missing-allele-kept",
        bases="GACGTA",
        input_lines=("2 A T ./1", "4 G C 1/1"),
        output_lines=("2 A T ./1", "4 G C 1/1"),
    ),
    RespellingCase(
        name="M8-three-positions-between-kept",
        bases="GACGTA",
        input_lines=("2 A T 1/1", "6 A C 1/1"),
        output_lines=("2 A T 1/1", "6 A C 1/1"),
    ),
    RespellingCase(
        name="M9-a-run-of-three",
        bases="GACGTA",
        input_lines=("2 A T 1/1", "4 G C 1/1", "6 A G 0/1"),
        output_lines=("2 ACGTA TCCTA,TCCTG 1/2:.",),
    ),
    RespellingCase(
        name="M10-ploidy-1",
        bases="GACGTA",
        input_lines=("2 A T 1", "4 G C 0"),
        output_lines=("2 ACG TCG 1:.",),
        ploidy=1,
    ),
    RespellingCase(
        name="M11-an-insertion-in-the-run",
        bases="GACGTA",
        input_lines=("2 A AT 1/1", "4 G C 0/1"),
        output_lines=("2 ACG ATCG,ATCC 1/2:.",),
    ),
    RespellingCase(
        name="M12-a-kept-record-touching-both-stretches",
        bases="GACGTA",
        input_lines=("2 A T 1/1", "3 C A 2/0", "4 G C 0/1"),
        output_lines=("2 A T 1/1", "3 C A 2/0", "4 G C 0/1"),
    ),
    # §6: a run whose every haplotype gets REF has no ALT to write.
    RespellingCase(
        name="a-run-all-ref-kept",
        bases="GACGTA",
        input_lines=("2 A T 0/0", "4 G C 0/0"),
        output_lines=("2 A T 0/0", "4 G C 0/0"),
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
