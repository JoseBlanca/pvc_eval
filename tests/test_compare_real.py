"""The check of doc/spec.md §7.3 over real data, as doc/specs/compare.md, "A file against
itself over real data", gives it: GIAB's HG002 v4.2.1 benchmark VCF for GRCh38 compared
with itself by compare(), over chr20 whole.

It needs two paths from the environment, and is skipped, with a message that names each
variable that is unset or whose file is missing (doc/specs/compare.md, "Decided by the
owner", point 6):

- PVC_EVAL_GRCH38, the GRCh38 analysis set without alts, with its .fai beside it;
- PVC_EVAL_HG002_VCF, the benchmark VCF, bgzipped, with its .tbi beside it.

File A and file B are one file, so every sample result of B is the one of A, and a row
can only be same genotype, reference in both, missing in both, invalid in both or too
complex. The last three are the ones where the file itself has, for the sample at the
locus, a missing haplotype, an invalid sample or a sample over the phasing cap; the test
checks that from what the row holds.
"""

import os
from pathlib import Path

import pytest

from pvc_eval.compare import CompareConfig, compare
from pvc_eval.fasta import IndexedFasta
from pvc_eval.outcome import Outcome, OutcomeRow
from pvc_eval.vcf import MISSING

pytestmark = pytest.mark.slow

_CONTIG = "chr20"
# How many rows the test shows when it fails, so that a broken stage does not print
# tens of thousands of them.
_ROWS_SHOWN = 20
# The fields of OutcomeRow that come in pairs, one for each file, without their _a and _b.
_SAMPLE_FIELDS = ("genotype", "possible", "quality", "gq", "records", "invalid_reason")


def _get_path_from_environment(variable: str, what: str, beside: str) -> tuple[Path | None, str]:
    """The path a variable names, and "" when it and the file beside it exist; otherwise
    None and why it cannot be used."""
    value = os.environ.get(variable)
    if not value:
        return None, f"{variable} is not set; it names {what}"
    path = Path(value)
    if not path.is_file():
        return None, f"{variable} names {path}, which does not exist; it names {what}"
    if not Path(f"{path}{beside}").is_file():
        return None, f"{variable} names {path}, which has no {path.name}{beside} beside it"
    return path, ""


def _find_why_broken(row: OutcomeRow) -> str:
    """Why a row breaks the property of a file against itself, or "" when it does not.

    Every row: its two samples have the same value in each field of _SAMPLE_FIELDS, since
    they are one sample result of one file's records. Then, by outcome:
    - same genotype and reference in both: nothing more;
    - missing in both: the genotype shown has a MISSING haplotype;
    - invalid in both: the invalid reason is not None;
    - too complex: the count of possible genotypes is None, a sample over the phasing cap.
    Any other outcome, those "in A" or "in B" alone, only in A or B, same alleles, other
    genotype, and different, cannot come from one file.
    """
    for field in _SAMPLE_FIELDS:
        if getattr(row, f"{field}_a") != getattr(row, f"{field}_b"):
            return f"the two samples have different {field}_a and {field}_b"
    match row.outcome:
        case Outcome.SAME_GENOTYPE | Outcome.REFERENCE_IN_BOTH:
            return ""
        case Outcome.MISSING_IN_BOTH:
            if row.genotype_a is None or MISSING not in row.genotype_a:
                return "missing in both, but the genotype has no missing haplotype"
            return ""
        case Outcome.INVALID_IN_BOTH:
            if row.invalid_reason_a is None:
                return "invalid in both, but the sample has no invalid reason"
            return ""
        case Outcome.TOO_COMPLEX:
            if row.possible_a is not None:
                return "too complex, but the sample is not over the phasing cap"
            return ""
        case _:
            return f"the outcome {row.outcome} cannot come from a file against itself"


def test_hg002_against_itself_over_chr20_gives_only_outcomes_of_one_file(
    tmp_path: Path,
) -> None:
    reference, reference_problem = _get_path_from_environment(
        "PVC_EVAL_GRCH38", "the GRCh38 FASTA, with its .fai beside it", ".fai"
    )
    vcf, vcf_problem = _get_path_from_environment(
        "PVC_EVAL_HG002_VCF",
        "GIAB's HG002 v4.2.1 benchmark VCF for GRCh38, bgzipped, with its .tbi beside it",
        ".tbi",
    )
    if reference is None or vcf is None:
        pytest.skip(". ".join(problem for problem in (reference_problem, vcf_problem) if problem))
    with IndexedFasta(reference) as indexed_fasta:
        contig_length = indexed_fasta.get_length(_CONTIG)
    regions = tmp_path / "chr20.bed"
    regions.write_text(f"{_CONTIG}\t0\t{contig_length}\n", encoding="utf-8")

    rows = compare(vcf, vcf, reference, CompareConfig(regions=regions))
    row_count = 0
    contigs = set()
    broken = []
    for row in rows:
        row_count += 1
        contigs.add(row.contig)
        reason = _find_why_broken(row)
        if reason:
            broken.append(f"{reason}: {row}")

    # The rows of HG002 v4.2.1 against itself over chr20, counted by the command
    # `pvc_eval compare` at 0c43267 and by this test at 5bafe18, on 20 October 2026: a
    # stage that drops or adds loci changes the count.
    assert row_count == 81_257
    assert contigs == {_CONTIG}
    assert not broken, f"{len(broken)} rows break the property:\n" + "\n".join(broken[:_ROWS_SHOWN])
