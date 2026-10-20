"""The check of doc/spec.md §7.4 at check_respelling_outcomes (doc/specs/respell.md §9.3):
every case of each re-spelling, input against output; the AT repeat through
segment_contig; hand-made files B that are not re-spellings, each of which it must reject;
the property over generated files, with their planted repeats as a --str-regions BED; and
each re-spelling of chromosome 20 of GIAB's HG002 benchmark (§9.4), which needs two paths
from the environment and is skipped, naming them, without them:

- PVC_EVAL_HG002_VCF, GIAB's HG002 v4.2.1 benchmark VCF for GRCh38, bgzipped;
- PVC_EVAL_GRCH38, GRCh38, the analysis set, with its .fai beside it.
"""

import gzip
import os
import re
from pathlib import Path

import pytest
from tests.conftest import CaseFiles, write_case_vcf
from tests.respell import _make_bgzf, respell_by_decomposing
from tests.respell_cases import RespellingCase
from tests.respell_check import _make_checked_rows, check_respelling_outcomes
from tests.respell_generator import (
    RESPELLING_CASE_PARAMS,
    RESPELLING_PARAMS,
    GeneratedVcf,
    Respelling,
    check_generated_respellings,
    check_respelled_haplotypes,
    is_rewritten,
    write_and_respell_case,
    write_planted_repeats,
)
from tests.test_loci_from_vcf import AT_REPEAT_REFERENCE

from pvc_eval.fasta import IndexedFasta
from pvc_eval.kind import Kind
from pvc_eval.outcome import Outcome
from pvc_eval.segments import GenomeSegment, RegionType, StrFinderParameters, write_str_regions

_PROPERTY_EXAMPLES = 300  # doc/specs/respell.md §9.2
_GIAB_CONTIG = "chr20"
_GIAB_SEED = 20261020  # §9.4: the seed of the re-spellings that draw, fixed


def _write_vcf(path: Path, case_lines: list[str]) -> Path:
    return write_case_vcf(path, case_lines, fills_ps=False)


@pytest.mark.parametrize(("respelling", "case"), RESPELLING_CASE_PARAMS)
def test_the_respelling_of_a_case_gives_the_outcomes_of_its_input_against_itself(
    case_files: CaseFiles, tmp_path: Path, respelling: Respelling, case: RespellingCase
) -> None:
    reference, vcf, output = write_and_respell_case(case_files, tmp_path, respelling, case)

    check_respelling_outcomes(vcf, output, reference, ploidy=case.ploidy)


def test_an_insertion_shifted_through_an_at_repeat_is_one_str_locus_of_same_genotype(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # §9.1: the AT repeat at 26-45, 1-based, found by segment_contig. `25 G GAT` inserts AT
    # before the repeat, and `45 T TAT` after it, the same haplotype: both span 25-46, so
    # they share one locus, which holds the repeat, (25, 45) 0-based and half-open.
    reference = case_files.write_contigs({"chr1": AT_REPEAT_REFERENCE})
    vcf = _write_vcf(tmp_path / "a.vcf", ["25 G GAT 0/1"])
    output = _write_vcf(tmp_path / "b.vcf", ["45 T TAT 0/1"])

    rows = list(_make_checked_rows(vcf, output, reference, str_regions=None, ploidy=2))

    at_repeat = GenomeSegment("chr1", 25, 45, RegionType.STR, "AT")
    assert [(row.start, row.end, row.kind, row.str_regions, row.outcome) for row in rows] == [
        (24, 46, Kind.STR, (at_repeat,), Outcome.SAME_GENOTYPE)
    ]


def test_a_file_whose_own_result_is_too_complex_passes_against_its_decomposition(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # 14 unphased heterozygous SNPs two bases apart, one locus at the chain distance of 10:
    # 2^13 phase assignments in A, over the cap of 4096, and as many in B, where each SNP is
    # rewritten in a phase set of its own. Too complex comes from A's own result.
    bases = "GACGTACGTACGTACGTACGTACGTACGTAC"
    reference = case_files.write_reference(1, bases)
    lines = [
        f"{pos} {bases[pos - 1]} {'T' if bases[pos - 1] != 'T' else 'A'} 0/1"
        for pos in range(2, 30, 2)
    ]
    vcf = _write_vcf(tmp_path / "a.vcf", lines)
    output = tmp_path / "b.vcf"
    respell_by_decomposing(vcf, reference, output)

    rows = list(_make_checked_rows(vcf, output, reference, str_regions=None, ploidy=2))

    assert is_rewritten(vcf, output)
    assert [(row.start, row.end, row.outcome) for row in rows] == [(1, 28, Outcome.TOO_COMPLEX)]


# Files B that are not re-spellings of A, on reference GACGTACGTACGTACGT from position 1:
# the case, A's lines, B's lines and the message, whose locus is 1-based and inclusive.
_MNP_REF = "ACGTACGTACGTAC"  # 2-15
_MNP_ALT = "CAAACAAACAAACA"
_NOT_RESPELLINGS = [
    pytest.param(
        ["2 A T 0/1"],
        ["2 A C 0/1"],
        "chr1:2-2, sample S0: different, an outcome a file does not give against itself;"
        " file A gives 1 genotype and file B 1 genotype; only in A: A/T; only in B: A/C",
        id="wrong-sequence",
    ),
    pytest.param(
        # B's GT names an allele past its ALTs, so B is invalid where A is not.
        ["2 A T 0/1"],
        ["2 A T 0/3"],
        "chr1:2-2, sample S0: invalid in B, an outcome a file does not give against itself;"
        " file A gives 1 genotype and file B invalid, malformed GT; only in A: A/T",
        id="invalid-in-b-alone",
    ),
    pytest.param(
        # Only the second sample differs, so a check of the first sample alone passes it.
        ["2 A T 0/1 0/1"],
        ["2 A T 0/1 1/1"],
        "chr1:2-2, sample S1: same alleles, other genotype, an outcome a file does not give"
        " against itself; file A gives 1 genotype and file B 1 genotype; only in A: A/T;"
        " only in B: T/T",
        id="wrong-genotype-in-the-second-sample",
    ),
    pytest.param(
        # The phase between the two SNPs is lost: B has A's genotype and one more, so the
        # outcome is still same genotype, and only the results tell them apart.
        ["2 A T 0|1:5", "4 G C 0|1:5"],
        ["2 A T 0/1", "4 G C 0/1"],
        "chr1:2-4, sample S0: same genotype, and file B's result is not file A's;"
        " file A gives 1 genotype and file B 2 genotypes; only in B: ACC/TCG",
        id="lost-phase-assignment",
    ),
    pytest.param(
        # A 14-base MNP written as 14 unphased SNPs: 2^13 phase assignments, over the cap.
        [f"2 {_MNP_REF} {_MNP_ALT} 0/1"],
        [
            f"{pos} {ref_base} {alt_base} 0/1"
            for pos, ref_base, alt_base in zip(range(2, 16), _MNP_REF, _MNP_ALT, strict=True)
        ],
        "chr1:2-15, sample S0: too complex, where file A's own result is not too complex;"
        f" file A gives 1 genotype and file B too complex; only in A: {_MNP_REF}/{_MNP_ALT}",
        id="too-complex-in-b-alone",
    ),
]


@pytest.mark.parametrize(("lines_a", "lines_b", "message"), _NOT_RESPELLINGS)
def test_a_file_b_that_is_not_a_respelling_raises(
    case_files: CaseFiles, tmp_path: Path, lines_a: list[str], lines_b: list[str], message: str
) -> None:
    reference = case_files.write_reference(1, "GACGTACGTACGTACGT")
    vcf = _write_vcf(tmp_path / "a.vcf", lines_a)
    output = _write_vcf(tmp_path / "b.vcf", lines_b)

    with pytest.raises(AssertionError, match=f"^{re.escape(message)}$"):
        check_respelling_outcomes(vcf, output, reference)


def test_a_message_lists_ten_genotypes_of_a_file_at_most(
    case_files: CaseFiles, tmp_path: Path
) -> None:
    # 12 SNPs of one phase set in A, unphased in B: B has 2^11 possible genotypes, of which
    # the message lists those A lacks, ten of them, and how many more there are.
    bases = "GACGTACGTACGTACGTACGTACGT"
    reference = case_files.write_reference(1, bases)
    alts = {pos: "T" if bases[pos - 1] != "T" else "A" for pos in range(2, 25, 2)}
    vcf = _write_vcf(
        tmp_path / "a.vcf", [f"{pos} {bases[pos - 1]} {alt} 0|1:5" for pos, alt in alts.items()]
    )
    output = _write_vcf(
        tmp_path / "b.vcf", [f"{pos} {bases[pos - 1]} {alt} 0/1" for pos, alt in alts.items()]
    )

    with pytest.raises(AssertionError) as raised:
        check_respelling_outcomes(vcf, output, reference)

    head, listed = str(raised.value).split("; only in B: ")
    assert head == (
        "chr1:2-24, sample S0: same genotype, and file B's result is not file A's;"
        " file A gives 1 genotype and file B 2048 genotypes"
    )
    assert listed.endswith(" and 2037 more")
    assert listed.count(" or ") == 9


def test_files_of_different_samples_raise(case_files: CaseFiles, tmp_path: Path) -> None:
    reference = case_files.write_reference(1, "GACGT")
    vcf = _write_vcf(tmp_path / "a.vcf", ["2 A T 0/1"])
    output = _write_vcf(tmp_path / "b.vcf", ["2 A T 0/1 0/1"])

    message = f"{vcf} has the samples S0 and {output} the samples S0, S1"
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        check_respelling_outcomes(vcf, output, reference)


def test_files_with_no_samples_raise(case_files: CaseFiles, tmp_path: Path) -> None:
    # With no sample there would be no row, and any file B would pass.
    reference = case_files.write_reference(1, "GACGT")
    vcf = write_case_vcf(tmp_path / "a.vcf", ["2 A T"], num_samples=0)
    output = write_case_vcf(tmp_path / "b.vcf", ["2 A C"], num_samples=0)

    message = f"{vcf} names no sample, so the check would compare nothing"
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        check_respelling_outcomes(vcf, output, reference)


# §9.2: the property asserts that at least one file in ten had a record rewritten, and that
# some row is in a locus with an STR region, so that a check that did not read the BED of
# the planted repeats fails: 9,415 of 12,152 rows were over 5,000 files on 20 October 2026.
@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_a_generated_file_against_its_respelling_gives_the_outcomes_of_itself(
    tmp_path_factory: pytest.TempPathFactory, respelling: Respelling
) -> None:
    def check(
        generated: GeneratedVcf, reference: IndexedFasta, vcf: Path, output: Path
    ) -> dict[str, int]:
        rows = _make_checked_rows(
            vcf,
            output,
            reference,
            str_regions=write_planted_repeats(generated, vcf.parent),
            ploidy=generated.ploidy,
        )
        return {"rows in an STR locus": sum(bool(row.str_regions) for row in rows)}

    tally = check_generated_respellings(
        respelling, check, tmp_path_factory, max_examples=_PROPERTY_EXAMPLES
    )
    assert tally["rows in an STR locus"] > 0, tally


# --- GIAB's HG002, §9.4 ---


_HG002_WHAT = "GIAB's HG002 v4.2.1 benchmark VCF for GRCh38"
_GRCH38_WHAT = "GRCh38, the analysis set, with its .fai beside it"


def _get_giab_paths() -> tuple[Path, Path]:
    """The HG002 VCF and GRCh38 the environment names. The test is skipped, naming each
    variable unset, when one is (§9.4); it fails, saying what is wrong, when a variable names
    a path that is not there, is a directory, or has no .fai beside it where one is needed."""
    needed = (
        ("PVC_EVAL_HG002_VCF", _HG002_WHAT, None),
        ("PVC_EVAL_GRCH38", _GRCH38_WHAT, ".fai"),
    )
    paths = []
    unset = []
    wrong = []
    for variable, what, beside in needed:
        value = os.environ.get(variable)
        if not value:
            unset.append(f"{variable} is not set; it names {what}")
            continue
        path = Path(value)
        problem = _find_path_problem(path, beside)
        if problem is None:
            paths.append(path)
        else:
            wrong.append(f"{variable} names {path}, which {problem}")
    if wrong:
        pytest.fail(". ".join(wrong))
    if unset:
        pytest.skip(". ".join(unset))
    hg002_vcf, grch38 = paths
    return hg002_vcf, grch38


def _find_path_problem(path: Path, beside: str | None) -> str | None:
    """What is wrong with a path the environment names, or None."""
    if path.is_dir():
        return "is a directory"
    if not path.is_file():
        return "is not there"
    if beside is not None and not Path(f"{path}{beside}").is_file():
        return f"has no {beside} beside it"
    return None


def test_the_giab_test_is_skipped_naming_each_variable_unset(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("PVC_EVAL_HG002_VCF", raising=False)
    monkeypatch.delenv("PVC_EVAL_GRCH38", raising=False)

    message = (
        f"PVC_EVAL_HG002_VCF is not set; it names {_HG002_WHAT}."
        f" PVC_EVAL_GRCH38 is not set; it names {_GRCH38_WHAT}"
    )
    with pytest.raises(pytest.skip.Exception, match=f"^{re.escape(message)}$"):
        _get_giab_paths()


@pytest.mark.parametrize(
    ("vcf_name", "fasta_name", "problems"),
    [
        pytest.param(
            "absent.vcf.gz",
            "ref.fa",
            ["PVC_EVAL_HG002_VCF names {vcf}, which is not there"],
            id="absent",
        ),
        pytest.param(
            "directory",
            "ref.fa",
            ["PVC_EVAL_HG002_VCF names {vcf}, which is a directory"],
            id="directory",
        ),
        pytest.param(
            "a.vcf.gz",
            "no_fai.fa",
            ["PVC_EVAL_GRCH38 names {fasta}, which has no .fai beside it"],
            id="no-fai",
        ),
        pytest.param(
            "directory",
            "absent.fa",
            [
                "PVC_EVAL_HG002_VCF names {vcf}, which is a directory",
                "PVC_EVAL_GRCH38 names {fasta}, which is not there",
            ],
            id="both",
        ),
    ],
)
def test_the_giab_test_fails_naming_what_is_wrong_with_a_path_set(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    vcf_name: str,
    fasta_name: str,
    problems: list[str],
) -> None:
    (tmp_path / "directory").mkdir()
    (tmp_path / "a.vcf.gz").write_bytes(b"")
    (tmp_path / "ref.fa").write_text(">chr20\nA\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("chr20\t1\t7\t1\t2\n", encoding="utf-8")
    (tmp_path / "no_fai.fa").write_text(">chr20\nA\n", encoding="utf-8")
    vcf, fasta = tmp_path / vcf_name, tmp_path / fasta_name
    monkeypatch.setenv("PVC_EVAL_HG002_VCF", str(vcf))
    monkeypatch.setenv("PVC_EVAL_GRCH38", str(fasta))

    message = ". ".join(problems).format(vcf=vcf, fasta=fasta)
    with pytest.raises(pytest.fail.Exception, match=f"^{re.escape(message)}$"):
        _get_giab_paths()


def _write_contig_vcf(vcf: Path, contig: str, output: Path) -> int:
    """Write the header of a bgzipped VCF and the records of one contig, plain, to output,
    and return how many records it wrote. The VCF is sorted, so its reading stops where the
    contig's records end."""
    prefix = f"{contig}\t".encode()
    num_records = 0
    with gzip.open(vcf, "rb") as lines, output.open("wb") as written:
        for line in lines:
            if line.startswith(prefix):
                num_records += 1
            elif num_records:
                break
            elif not line.startswith(b"#"):
                continue
            written.write(line)
    return num_records


def test_the_records_of_one_contig_are_written_with_the_header(tmp_path: Path) -> None:
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr2>",
        "##contig=<ID=chr20>",
        "##contig=<ID=chr21>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002",
    ]
    records = {
        contig: [f"{contig}\t{pos}\t.\tA\tT\t.\t.\t.\tGT\t0/1" for pos in (5, 9)]
        for contig in ("chr2", "chr20", "chr21")
    }
    lines = [*header, *records["chr2"], *records["chr20"], *records["chr21"]]
    vcf = tmp_path / "all.vcf.gz"
    vcf.write_bytes(_make_bgzf("".join(f"{line}\n" for line in lines).encode()))
    output = tmp_path / "chr20.vcf"

    num_records = _write_contig_vcf(vcf, "chr20", output)

    assert num_records == 2
    assert output.read_text(encoding="utf-8").splitlines() == [*header, *records["chr20"]]


@pytest.mark.slow
@pytest.mark.parametrize("respelling", RESPELLING_PARAMS)
def test_each_respelling_of_hg002_chr20_passes_the_oracle_and_gives_the_outcomes_of_itself(
    tmp_path: Path, respelling: Respelling
) -> None:
    hg002_vcf, grch38 = _get_giab_paths()
    vcf = tmp_path / "a.vcf"
    num_records = _write_contig_vcf(hg002_vcf, _GIAB_CONTIG, vcf)
    assert num_records > 0, f"{hg002_vcf} has no record on {_GIAB_CONTIG}"
    contig_bed = tmp_path / "contig.bed"
    with IndexedFasta(grch38) as reference:
        length = reference.get_length(_GIAB_CONTIG)
        contig_bed.write_text(f"{_GIAB_CONTIG}\t0\t{length}\n", encoding="utf-8")
        str_regions = tmp_path / "str_regions.bed"
        write_str_regions(reference, str_regions, StrFinderParameters(), contig_bed)
        output = tmp_path / "b.vcf"
        respelling.respell(vcf, reference, output, 2, _GIAB_SEED)

        comparison = check_respelled_haplotypes(vcf, output, reference, ploidy=2)
        check_respelling_outcomes(vcf, output, reference, str_regions=str_regions)

    assert comparison.num_windows > 0
