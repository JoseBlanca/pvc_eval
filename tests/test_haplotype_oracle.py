"""The oracle of doc/specs/respell.md §8, at check_same_haplotypes: the cases O1 to O6, a
window over ORACLE_CAP, the rules of one haplotype over one window, and the modules it
imports. Files are written in the notation of the cases, which tests.respell_cases
describes.
"""

import re
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

import pytest
from tests.conftest import CaseFiles, make_chrom_line, write_case_vcf
from tests.haplotype_oracle import ORACLE_CAP, HaplotypeComparison, check_same_haplotypes

from pvc_eval.fasta import IndexedFasta

_REPOSITORY = Path(__file__).parents[1]
_O_BASES = "GACGTACGTAC"  # the reference of the O cases, from position 1


def _check(
    directory: Path,
    reference: IndexedFasta,
    lines_a: Sequence[str],
    lines_b: Sequence[str],
    *,
    ploidy: int = 2,
) -> HaplotypeComparison:
    vcf_a = write_case_vcf(directory / "a.vcf", lines_a, fills_ps=False)
    vcf_b = write_case_vcf(directory / "b.vcf", lines_b, fills_ps=False)
    return check_same_haplotypes(vcf_a, vcf_b, reference, ploidy=ploidy)


# O1, O2, O3 and O6, and wrong re-spellings the oracle must reject besides them.
@pytest.mark.parametrize(
    ("ref_start_1based", "bases", "lines_a", "lines_b", "message"),
    [
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 0/1", "4 G C 0/1"],
            ["2 ACG TCG,ACC 1|2:1"],
            "chr1:2-4, sample S0: file A gives 2 genotypes and file B 1 genotype;"
            " only in A: {ACG, TCC}",
            id="O1-a-merge-that-invents-a-phase",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 0|1:5", "10 A G 0|1:5"],
            ["2 A T 0|1:5", "10 A G 0|1:6"],
            "chr1:2-2 and 10-10, sample S0: file A gives 1 genotype and file B 2 genotypes;"
            " only in B: {(A, G), (T, A)}",
            id="O2-a-phase-set-lost-between-two-windows",
        ),
        pytest.param(
            9,
            "CAGTC",
            ["9 CAGT C 0|1:5", "11 G A 1|0:5"],
            ["9 CAGT C 0|1:5", "11 G A,* 2|1:5"],
            "chr1:9-12, sample S0: file A gives 1 genotype and file B no genotype;"
            " only in A: {C, CAAT}",
            id="O3-a-star-on-the-haplotype-without-the-deletion",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 0|1:5", "6 A C 1|1:5", "10 A G 0|1:5"],
            ["2 A T 0|1:5", "6 A C 1|1:5", "10 A G 0|1:6"],
            "chr1:2-2 and 10-10, sample S0: file A gives 1 genotype and file B 2 genotypes;"
            " only in B: {(A, G), (T, A)}",
            id="O6-a-homozygous-window-joins-nothing",
        ),
        # The pieces of a phased record put in a new phase set lose its phase with the
        # record at 6.
        pytest.param(
            1,
            _O_BASES,
            ["2 AC GT 0|1:5", "6 A C 0|1:5"],
            ["2 A G 0|1:9", "3 C T 0|1:9", "6 A C 0|1:5"],
            "chr1:2-3 and 6-6, sample S0: file A gives 1 genotype and file B 2 genotypes;"
            " only in B: {(AC, C), (GT, A)}",
            id="pieces-of-a-phased-record-in-a-new-phase-set",
        ),
        # Two unphased records given one phase set: B gives one genotype where A gives two.
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 0/1", "10 A G 0/1"],
            ["2 A T 0|1:3", "10 A G 0|1:3"],
            "chr1:2-2 and 10-10, sample S0: file A gives 2 genotypes and file B 1 genotype;"
            " only in A: {(A, G), (T, A)}",
            id="two-unphased-records-joined",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 0/1", "2 AC . ./."],
            ["2 A T 0/1"],
            "chr1:2-2, sample S0: file A gives 1 genotype and file B 1 genotype;"
            " only in A: {., .}; only in B: {A, T}",
            id="a-record-with-no-alt-makes-its-window-missing",
        ),
        # A `.` might be the deletion that covers the `*`, so A is missing, not impossible.
        pytest.param(
            9,
            "CAGTC",
            ["9 CAGT C ./.", "11 G A,* 1/2"],
            ["9 CAGT C 0/0", "11 G A,* 1/2"],
            "chr1:9-12, sample S0: file A gives 1 genotype and file B no genotype;"
            " only in A: {., .}",
            id="a-star-beside-a-missing-haplotype-is-missing",
        ),
        # §8: a `.` covers a `*` only through an edit of an ALT of its own record, and the
        # T at 2 replaces no base of the `*` record's REF, 3.
        pytest.param(
            1,
            _O_BASES,
            ["2 A T .|0:1", "3 C G,* 2|0:1"],
            ["2 A T .|0:1", "3 C G,* 1|0:1"],
            "chr1:2-3, sample S0: file A gives no genotype and file B 1 genotype;"
            " only in B: {AC, .}",
            id="a-missing-haplotype-whose-alts-do-not-cover-the-star",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 0/5"],
            ["2 A T 0/1"],
            "chr1:2-2, sample S0: file A gives malformed GT and file B 1 genotype;"
            " only in B: {A, T}",
            id="a-malformed-gt-in-one-file",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 A T 1"],
            ["2 A T 0/5"],
            "chr1:2-2, sample S0: file A gives ploidy mismatch and file B malformed GT",
            id="ploidy-mismatch-against-malformed-gt",
        ),
        # §5.4: an insertion between two positions one edit replaces collides with it.
        pytest.param(
            1,
            _O_BASES,
            ["2 ACG A 1/1", "3 C CT 1/1"],
            ["2 ACG A 1/1", "3 C CT 0/0"],
            "chr1:2-4, sample S0: file A gives no genotype and file B 1 genotype;"
            " only in B: {A, A}",
            id="an-insertion-inside-a-deletion-collides",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 A AT 1/1", "2 A AG 1/1"],
            ["2 A AT 1/1", "2 A AG 0/0"],
            "chr1:2-3, sample S0: file A gives no genotype and file B 1 genotype;"
            " only in B: {ATC, ATC}",
            id="two-insertions-at-one-point-collide",
        ),
        pytest.param(
            1,
            _O_BASES,
            ["2 AC GT 1/1", "3 C A 1/1"],
            ["2 AC GT 1/1", "3 C A 0/0"],
            "chr1:2-3, sample S0: file A gives no genotype and file B 1 genotype;"
            " only in B: {GT, GT}",
            id="two-replacements-of-one-position-collide",
        ),
    ],
)
def test_a_file_whose_haplotypes_differ_raises_with_the_window_and_both_sets(
    case_files: CaseFiles,
    tmp_path: Path,
    ref_start_1based: int,
    bases: str,
    lines_a: Sequence[str],
    lines_b: Sequence[str],
    message: str,
) -> None:
    reference = case_files.write_reference(ref_start_1based, bases)
    with pytest.raises(AssertionError, match=f"^{re.escape(message)}$"):
        _check(tmp_path, reference, lines_a, lines_b)


# O4, O5 and spellings the oracle must take as the same.
@pytest.mark.parametrize(
    ("lines_a", "lines_b", "comparison"),
    [
        pytest.param(
            ["2 AC GT ./1"],
            ["2 A G .|1:1", "3 C T .|1:1"],
            HaplotypeComparison(num_windows=1, num_linked_pairs=0, num_over_cap=0),
            id="O4-a-missing-haplotype-decomposed",
        ),
        pytest.param(
            ["2 A T 0/1"],
            ["2 A T 0/1"],
            HaplotypeComparison(num_windows=1, num_linked_pairs=0, num_over_cap=0),
            id="O5-one-window-and-no-pair",
        ),
        pytest.param(
            ["2 A T 0|1:5 0/1", "10 A G 0|1:5 1/1"],
            ["2 A T 0|1:5 0/1", "10 A G 0|1:5 1/1"],
            HaplotypeComparison(num_windows=4, num_linked_pairs=1, num_over_cap=0),
            id="windows-and-pairs-counted-over-every-sample",
        ),
        # Replacements next to each other, and an insertion at the end of a deleted
        # stretch, do not collide.
        pytest.param(
            ["2 A T 1/1", "3 C G 1/1"],
            ["2 AC TG 1/1"],
            HaplotypeComparison(num_windows=1, num_linked_pairs=0, num_over_cap=0),
            id="replacements-next-to-each-other",
        ),
        pytest.param(
            ["2 ACG A 1/1", "4 G GT 1/1"],
            ["2 ACG AT 1/1"],
            HaplotypeComparison(num_windows=1, num_linked_pairs=0, num_over_cap=0),
            id="an-insertion-at-the-end-of-a-deletion",
        ),
        pytest.param(
            ["2 A T 1", "4 G C 1"],
            ["2 A T 1", "4 G C 1"],
            HaplotypeComparison(num_windows=2, num_linked_pairs=0, num_over_cap=0),
            id="ploidy-mismatch-in-both",
        ),
    ],
)
def test_two_spellings_of_the_same_haplotypes_pass(
    case_files: CaseFiles,
    tmp_path: Path,
    lines_a: Sequence[str],
    lines_b: Sequence[str],
    comparison: HaplotypeComparison,
) -> None:
    reference = case_files.write_reference(1, _O_BASES)
    assert _check(tmp_path, reference, lines_a, lines_b) == comparison


def test_ploidy_1_gives_one_genotype_a_window(case_files: CaseFiles, tmp_path: Path) -> None:
    reference = case_files.write_reference(1, _O_BASES)
    # The records at 2-3 and 4 lie next to each other: one window, 2-4.
    lines_a = ["2 AC GT 1", "4 G C 0"]
    with pytest.raises(AssertionError, match=re.escape("only in A: {GTG}; only in B: {ACG}")):
        _check(tmp_path, reference, lines_a, ["2 AC GT 0", "4 G C 0"], ploidy=1)
    assert _check(
        tmp_path, reference, lines_a, ["2 A G 1", "3 C T 1"], ploidy=1
    ) == HaplotypeComparison(num_windows=1, num_linked_pairs=0, num_over_cap=0)


def test_a_missing_haplotype_whose_alts_do_not_cover_a_star_is_impossible_with_ploidy_1(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    reference = case_files.write_reference(1, _O_BASES)
    message = "chr1:2-3, sample S0: file A gives no genotype and file B 1 genotype; only in B: {.}"
    with pytest.raises(AssertionError, match=f"^{re.escape(message)}$"):
        _check(
            tmp_path,
            reference,
            ["2 A T .", "3 C G,* 2"],
            ["2 A T .", "3 C G,* 1"],
            ploidy=1,
        )


def test_a_window_over_the_cap_is_counted_and_not_compared(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # Fourteen unphased heterozygous SNPs that touch, one window: 2^13 = 8192 phase
    # assignments in A, over the cap of 4096. B, all in one phase set, would differ.
    assert ORACLE_CAP == 4096
    bases = "GACGTACGTACGTACGT"
    reference = case_files.write_reference(1, bases)
    positions = range(2, 16)
    snps = [f"{pos} {bases[pos - 1]} {'C' if bases[pos - 1] == 'A' else 'A'}" for pos in positions]
    lines_a = [f"{snp} 0/1" for snp in snps]
    lines_b = [f"{snp} 0|1:1" for snp in snps]
    assert _check(tmp_path, reference, lines_a, lines_b) == HaplotypeComparison(
        num_windows=0, num_linked_pairs=0, num_over_cap=1
    )
    # Thirteen give 4096 assignments, within the cap, and are compared.
    with pytest.raises(AssertionError, match=r"^chr1:2-14, sample S0: file A gives "):
        _check(tmp_path, reference, lines_a[:13], lines_b[:13])


def test_the_message_of_many_genotypes_gives_ten_and_how_many_more(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # Twelve heterozygous SNPs that touch: unphased, 2^11 = 2048 possible genotypes in A; in
    # one phase set, 1 in B, which A also has.
    bases = "GACGTACGTACGTACGT"
    reference = case_files.write_reference(1, bases)
    snps = [
        f"{pos} {bases[pos - 1]} {'C' if bases[pos - 1] == 'A' else 'A'}" for pos in range(2, 14)
    ]
    with pytest.raises(AssertionError) as raised:
        _check(
            tmp_path, reference, [f"{snp} 0/1" for snp in snps], [f"{snp} 0|1:1" for snp in snps]
        )
    message = str(raised.value)
    assert message.startswith(
        "chr1:2-13, sample S0: file A gives 2048 genotypes and file B 1 genotype; only in A: {"
    )
    assert message.endswith("} and 2037 more")
    assert message.count("{") == 10
    assert "only in B" not in message


@pytest.mark.parametrize(
    ("line_b", "chrom_line_b", "samples_b"),
    [
        pytest.param("2 A T 0/1", make_chrom_line(1), "S0", id="another-number"),
        pytest.param(
            "2 A T 0/1 0/1",
            make_chrom_line(2).replace("S0\tS1", "S1\tS0"),
            "S1, S0",
            id="the-same-names-in-another-order",
        ),
    ],
)
def test_files_of_different_samples_raise_value_error_naming_both(
    case_files: CaseFiles, tmp_path: Path, line_b: str, chrom_line_b: str, samples_b: str
) -> None:
    reference = case_files.write_reference(1, _O_BASES)
    vcf_a = write_case_vcf(tmp_path / "a.vcf", ["2 A T 0/1 0/1"], fills_ps=False)
    vcf_b = write_case_vcf(tmp_path / "b.vcf", [line_b], fills_ps=False)
    num_samples_b = len(line_b.split()) - 3
    vcf_b.write_text(vcf_b.read_text().replace(make_chrom_line(num_samples_b), chrom_line_b))
    message = f"{vcf_a} has the samples S0, S1 and {vcf_b} the samples {samples_b}"
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        check_same_haplotypes(vcf_a, vcf_b, reference)


def test_the_oracle_imports_nothing_of_what_the_check_of_7_4_tests() -> None:
    # In a process of its own, so that what the other tests imported is not counted.
    code = (
        "import sys, tests.haplotype_oracle;"
        " print(sorted(name for name in sys.modules if name.startswith('pvc_eval')))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=_REPOSITORY, capture_output=True, text=True, check=True
    )
    assert result.stdout.split("\n")[0] == str(
        ["pvc_eval", "pvc_eval._text_files", "pvc_eval.errors", "pvc_eval.fasta", "pvc_eval.vcf"]
    )
