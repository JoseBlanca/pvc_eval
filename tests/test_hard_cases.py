"""The 50 cases of doc/specs/hard_cases.md through the oracle and through compare(), as its
§5 gives them, and the excerpts of tests/hard_cases.py checked against the references.

test_hard_case writes a case's excerpt and its two files, checks that the oracle of
tests/haplotype_oracle.py gives each file the case's possible genotypes over the locus, and
that compare() gives the case's row and no other.

The two excerpt tests need a reference each from the environment, and are skipped, with a
message that names the variable, when it is unset or its file or the .fai beside it is
missing: PVC_EVAL_GRCH38 for the 45 cases on GRCh38, PVC_EVAL_TOMATO for the five on SL4.0,
E12 to E16.
"""

import os
from collections.abc import Sequence
from pathlib import Path

import pytest
from tests.compare_case import write_fasta
from tests.haplotype_oracle import find_possible_haplotypes
from tests.hard_cases import HARD_CASES, HardCase, Reference

from pvc_eval.compare import compare
from pvc_eval.fasta import IndexedFasta
from pvc_eval.outcome import OutcomeRow
from pvc_eval.segments import GenomeSegment, RegionType, segment_contig

# doc/specs/hard_cases.md §1: the bases the STR regions of an excerpt are compared over,
# beyond each end of it.
_WIDENING = 5_000
# The fields of OutcomeRow §5 checks: QUAL and GQ are not, since the cases set them.
_CHECKED_FIELDS = (
    "contig",
    "start",
    "end",
    "kind",
    "str_regions",
    "outcome",
    "alleles",
    "genotype_a",
    "genotype_b",
    "possible_a",
    "possible_b",
    "invalid_reason_a",
    "invalid_reason_b",
    "records_a",
    "records_b",
)


def _write_vcf(path: Path, contig: str, sample: str, lines: Sequence[str]) -> Path:
    """Writes a case's lines, `POS REF ALT GT`, as a VCF of one sample with QUAL 50, FILTER
    PASS and INFO `.`; FORMAT is GT:PS for a GT with its PS, GT otherwise."""
    text_lines = [
        "##fileformat=VCFv4.3",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]
    for line in lines:
        pos, ref, alt, sample_field = line.split()
        vcf_format = "GT:PS" if ":" in sample_field else "GT"
        text_lines.append(
            "\t".join((contig, pos, ".", ref, alt, "50", "PASS", ".", vcf_format, sample_field))
        )
    path.write_text("".join(f"{text_line}\n" for text_line in text_lines), encoding="utf-8")
    return path


def _make_expected_row(case: HardCase) -> dict[str, object]:
    """The fields of _CHECKED_FIELDS the case expects, with its excerpt positions made
    0-based and half-open, as OutcomeRow holds them."""
    start_1based, end = case.locus_1based
    str_regions = tuple(
        GenomeSegment(case.contig, region_start_1based - 1, region_end, region_type, motif)
        for region_start_1based, region_end, region_type, motif in case.str_regions_1based
    )
    return {
        "contig": case.contig,
        "start": start_1based - 1,
        "end": end,
        "kind": case.kind,
        "str_regions": str_regions,
        "outcome": case.outcome,
        "alleles": case.alleles,
        "genotype_a": case.genotype_a,
        "genotype_b": case.genotype_b,
        "possible_a": case.possible_a,
        "possible_b": case.possible_b,
        "invalid_reason_a": case.invalid_reason_a,
        "invalid_reason_b": case.invalid_reason_b,
        "records_a": case.records_a,
        "records_b": case.records_b,
    }


def _get_checked_fields(row: OutcomeRow) -> dict[str, object]:
    return {field: getattr(row, field) for field in _CHECKED_FIELDS}


@pytest.mark.parametrize("case", HARD_CASES, ids=[case.name for case in HARD_CASES])
def test_hard_case(case: HardCase, tmp_path: Path) -> None:
    reference = tmp_path / "excerpt.fa"
    write_fasta(reference, {case.contig: case.bases})
    vcf_a = _write_vcf(tmp_path / "a.vcf", case.contig, case.sample, case.lines_a)
    vcf_b = _write_vcf(tmp_path / "b.vcf", case.contig, case.sample, case.lines_b)
    start_1based, end = case.locus_1based

    with IndexedFasta(reference) as indexed_fasta:
        oracle_genotypes = [
            find_possible_haplotypes(vcf, indexed_fasta, case.contig, start_1based - 1, end, 0)
            for vcf in (vcf_a, vcf_b)
        ]
    rows = list(compare(vcf_a, vcf_b, reference, finder_threads=0))

    assert oracle_genotypes == [case.oracle_genotypes_a, case.oracle_genotypes_b]
    # One row for each of these cases, compared field by field, so that a failure names the
    # fields that differ.
    assert len(rows) == 1, rows
    assert _get_checked_fields(rows[0]) == _make_expected_row(case)


def _get_reference(variable: str, what: str) -> Path:
    """The FASTA a variable names; skips the test when it is unset, or when the file or
    the .fai beside it is missing."""
    value = os.environ.get(variable)
    if not value:
        pytest.skip(f"{variable} is not set; it names {what}")
    path = Path(value)
    if not path.is_file():
        pytest.skip(f"{variable} names {path}, which does not exist; it names {what}")
    if not Path(f"{path}.fai").is_file():
        pytest.skip(f"{variable} names {path}, which has no {path.name}.fai beside it")
    return path


def _find_excerpt_problems(
    reference_path: Path, reference: Reference, directory: Path
) -> list[str]:
    """What differs between the cases of a reference and the reference itself: an
    excerpt that is not the reference's bases at its stretch, or STR regions of a locus
    other than those segment_contig finds over the excerpt widened by _WIDENING bases on
    each side."""
    problems = []
    with IndexedFasta(reference_path) as indexed_fasta:
        for case in HARD_CASES:
            if case.reference is not reference:
                continue
            excerpt_end = case.offset + len(case.bases)
            bases = indexed_fasta.get_bases(case.contig, case.offset, excerpt_end).upper()
            if bases != case.bases:
                problems.append(f"{case.name}: the excerpt is not the reference's bases")
            widened_start = case.offset - _WIDENING
            widened = indexed_fasta.get_bases(
                case.contig, widened_start, excerpt_end + _WIDENING
            ).upper()
            str_regions = _find_locus_str_regions(case, widened, widened_start, directory)
            if str_regions != case.str_regions_1based:
                problems.append(
                    f"{case.name}: the widened excerpt gives the STR regions {str_regions}"
                )
    return problems


def _find_locus_str_regions(
    case: HardCase, widened: str, widened_start: int, directory: Path
) -> tuple[tuple[int, int, RegionType | None, str | None], ...]:
    """The STR regions segment_contig finds inside a case's locus over the widened bases,
    which start at widened_start of the contig, 0-based; in excerpt positions, 1-based and
    inclusive, as the case gives them."""
    fasta = directory / f"{case.name}.fa"
    write_fasta(fasta, {case.contig: widened})
    # An excerpt position p is the widened sequence's 0-based position p - 1 + shift.
    shift = case.offset - widened_start
    start_1based, end = case.locus_1based
    with IndexedFasta(fasta) as indexed_fasta:
        return tuple(
            (segment.start - shift + 1, segment.end - shift, segment.region_type, segment.motif)
            for segment in segment_contig(indexed_fasta, case.contig)
            if segment.region_type is not RegionType.GENERIC
            and start_1based - 1 + shift <= segment.start
            and segment.end <= end + shift
        )


def test_grch38_excerpts_are_the_references(tmp_path: Path) -> None:
    path = _get_reference("PVC_EVAL_GRCH38", "the GRCh38 FASTA, with its .fai beside it")
    assert _find_excerpt_problems(path, Reference.GRCH38, tmp_path) == []


def test_tomato_excerpts_are_the_references(tmp_path: Path) -> None:
    path = _get_reference("PVC_EVAL_TOMATO", "the SL4.0 tomato FASTA, with its .fai beside it")
    assert _find_excerpt_problems(path, Reference.SL4_0, tmp_path) == []
