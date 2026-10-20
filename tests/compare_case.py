"""The four-contig case of doc/specs/compare.md, "How it is verified", which the tests of
compare() share: a reference of four contigs, chr1 to chr3 of 300 bases each with STR
regions the STR finder admits and chr4 of 200 bases with none, and the VCF text of file A, with records on chr1 and chr3, and of
file B, with records on chr2 and chr3, two samples in each and one sample pair by name,
HG002.

Positions in the VCF lines are 1-based, as in a VCF; the STR regions are 0-based and
half-open, as GenomeSegment gives them.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from pvc_eval.segments import GenomeSegment, RegionType

# Filler with no run of one base longer than two and no pair of bases twice in a row, with
# perfect repeats planted in it: chr1, (AC)8 at 101-116 and (AGT)6 at 201-218; chr2,
# (ATTC)5 at 121-140; chr3, (CA)10 at 101-120; chr4, none (1-based, inclusive). Where a
# repeat meets the filler, a C makes a run of three: CCC at 116-118 of chr1 and at 99-101
# of chr3.
CONTIGS = {
    "chr1": (
        "AAGCGGCACTTGTAAGTGCCACGCCGCTTGGTCTTCTGTTGTTCGCTGGTGCTGAGCAAGCACGCCATAAGGCCAAGGCCATAC"
        "CAAGATAGTAGTCTAGACACACACACACACACCCATCACCTAGACGGTGACATTCAACAACCACTTGTCCTTAATCATGAAGG"
        "ATAAGCATACAAGAACTCAGTTCGTAGAAGTCAAGTAGTAGTAGTAGTAGTAAGCCTAACGTCGTCGACTAGCGTCTGCTTAT"
        "CTATGTTGGACCTTAGTTCAATCTATCGCTCATTGCTCAGATAGTGAAGC"
    ),
    "chr2": (
        "TGCACTTGCAGTAGATTCGTCTGAGGTACTCAGACTCGAATGCGGAGTGCTTGTCTGGCACTCGCCGTTGGTGAGGTTCGGTT"
        "ACGTCAAGCGATAGCTGTCGGCTACCGGCTGGAGCCAATTCATTCATTCATTCATTCTTCTTAATCACTGTAGACCACTAGTA"
        "TCATCACCAGCCGTACAATCACTGTCACCTCGGTCTTGGAATGGTGCTCAACCTACAGTACCGACAATGCCGGATTATGAGCT"
        "GGTCTTGTTGCTTCTGGACGTCCGCAACGAGGTATTAGCCTATGATTCCGC"
    ),
    "chr3": (
        "CGTTCCAGCCTTATTGCCAATTCGAGGTATCGAATACCGCACGAACTCAGGTAGGAGGTGCAAGTAGAATTCCAAGCGAACCT"
        "AGAACCAATAGCATTCCCACACACACACACACACACATTGCGATAGATGGCTTGTCCTGGTACTATTATTGGCCTTCTGGTGG"
        "ATACTAAGGTCGATTCTAAGATCAAGTTATCCGCTTGACGCCTCGCCATTGCCTACCAATCCGTAAGAAGTTAATCCTAGCTA"
        "GGACATCCGTCAGTACCGGACCAGAAGACGCTCGAAGCAACTTGTACAACG"
    ),
    "chr4": (
        "CGCACCGACTCAGTTGCAACTCCGAACCAGCCTTAGCAGATAAGGCGTCACCTCAGTTAATAACTACTGCCGGCGGTTGTCTG"
        "TTGAAGTTATGCCGACCTCCTCAGTCAGCCATAGCCTCCGGCATAATCGGATGCTACGGTGGAGTCCTTCTGACATACAAGCT"
        "TGAACAACAGGAAGGATCTACCTAGACCACCACC"
    ),
}

# The STR regions segment_contig finds in CONTIGS with its default parameters, checked by
# test_compare.py.
STR_REGIONS = (
    GenomeSegment("chr1", 100, 116, RegionType.STR, "AC"),
    GenomeSegment("chr1", 200, 218, RegionType.STR, "AGT"),
    GenomeSegment("chr2", 120, 140, RegionType.STR, "ATTC"),
    GenomeSegment("chr3", 100, 120, RegionType.STR, "CA"),
)

SAMPLES_A = ("HG002", "HG003")
SAMPLES_B = ("HG002", "HG004")

# Each line is CHROM POS REF ALT QUAL and a GT:GQ per sample.
LINES_A = (
    "chr1 50 T C 50 0/1:40 1/1:35",
    "chr1 100 GAC G 29.5 0/1:45 0/0:50",  # a deletion of AC in (AC)8
    "chr1 150 T A 12 0/0:20 0/1:22",
    "chr3 105 C G 40 0/1:33 1/1:30",  # a SNP in (CA)10, which widens its locus
    "chr3 250 G T 60 1/1:50 0/1:41",
)
LINES_B = (
    "chr2 120 A AATTC 35 0/1:38 0/0:40",  # an insertion of ATTC before (ATTC)5
    "chr2 200 T G 45 1/1:42 0/1:39",
    "chr3 105 C G 41 0/1:36 0/0:37",
    "chr3 250 G T 58 0/1:47 1/1:44",
)

_FASTA_LINE_BASES = 60


@dataclass(frozen=True, slots=True, kw_only=True)
class CompareCase:
    """The paths of a case: the reference, with its .fai beside it, and files A and B."""

    reference: Path
    a: Path
    b: Path


def write_compare_case(
    directory: Path,
    *,
    contigs: Mapping[str, str] = CONTIGS,
    lines_a: Sequence[str] = LINES_A,
    lines_b: Sequence[str] = LINES_B,
) -> CompareCase:
    """Writes the reference, its .fai and the two VCF files of a case to directory."""
    reference = directory / "ref.fa"
    _write_fasta(reference, contigs)
    path_a = directory / "a.vcf"
    path_b = directory / "b.vcf"
    _write_vcf(path_a, SAMPLES_A, lines_a)
    _write_vcf(path_b, SAMPLES_B, lines_b)
    return CompareCase(reference=reference, a=path_a, b=path_b)


def write_bed(path: Path, lines: Sequence[str]) -> Path:
    """Writes lines, their columns separated by spaces, as a BED separated by tabs."""
    path.write_text("".join("\t".join(line.split()) + "\n" for line in lines), encoding="utf-8")
    return path


def _write_fasta(path: Path, contigs: Mapping[str, str]) -> None:
    """A FASTA of lines of _FASTA_LINE_BASES bases, and its .fai."""
    fasta_text = ""
    fai_lines = []
    for name, bases in contigs.items():
        fasta_text += f">{name}\n"
        fai_lines.append(
            f"{name}\t{len(bases)}\t{len(fasta_text.encode())}\t{_FASTA_LINE_BASES}"
            f"\t{_FASTA_LINE_BASES + 1}\n"
        )
        for start in range(0, len(bases), _FASTA_LINE_BASES):
            fasta_text += bases[start : start + _FASTA_LINE_BASES] + "\n"
    path.write_text(fasta_text, encoding="utf-8")
    Path(f"{path}.fai").write_text("".join(fai_lines), encoding="utf-8")


def _write_vcf(path: Path, samples: Sequence[str], lines: Sequence[str]) -> None:
    text = "##fileformat=VCFv4.3\n"
    text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    for line in lines:
        chrom, pos, ref, alt, qual, *sample_fields = line.split()
        fixed = (chrom, pos, ".", ref, alt, qual, "PASS", ".", "GT:GQ")
        text += "\t".join((*fixed, *sample_fields)) + "\n"
    path.write_text(text, encoding="utf-8")
