"""The check of doc/spec.md §7.4, as doc/specs/respell.md §9.3 gives it: file A against a
re-spelling of itself, file B, through the stages of pvc_eval, reading, loci,
reconciliation and evaluation, joined here as compare() will join them. Every row must be
one that a file gives against itself (§7.3), and every sample's result in B must be its
result in A.

Positions are 0-based and half-open, but in the messages, which give a locus 1-based and
inclusive as a VCF does.

Public: check_respelling_outcomes.
"""

import itertools
from collections.abc import Iterator, Sequence
from pathlib import Path

from pvc_eval.evaluate import make_outcome_rows
from pvc_eval.fasta import IndexedFasta
from pvc_eval.kind import locus_kind
from pvc_eval.loci import build_loci
from pvc_eval.outcome import InvalidReason, Outcome, OutcomeRow
from pvc_eval.reconcile import TOO_COMPLEX, ReconciledLocus, SampleResult, reconcile_loci
from pvc_eval.samples import SamplePair
from pvc_eval.segments import GenomeSegment, read_str_regions, segment_contig
from pvc_eval.vcf import DEFAULT_PLOIDY, MISSING, parse_vcf, parse_vcf_samples

# The outcomes a file gives against itself (doc/spec.md §7.3). Of the three that come from a
# cause in the sample, missing in both and invalid in both need A's own result checked no
# further: evaluate_pair gives them only when both results are missing, or both invalid. Too
# complex comes when either result is, so for it A's own result is checked.
_ALLOWED_OUTCOMES = frozenset(
    {
        Outcome.SAME_GENOTYPE,
        Outcome.REFERENCE_IN_BOTH,
        Outcome.MISSING_IN_BOTH,
        Outcome.INVALID_IN_BOTH,
        Outcome.TOO_COMPLEX,
    }
)


def check_respelling_outcomes(
    vcf_a: Path,
    vcf_b: Path,
    reference: IndexedFasta,
    *,
    str_regions: Path | None = None,
    ploidy: int = DEFAULT_PLOIDY,
) -> None:
    """Check file A against file B, a re-spelling of it, through parse_vcf, build_loci,
    reconcile_loci and make_outcome_rows, each sample paired with itself (doc/specs/respell.md
    §9.3): at every locus, every row's outcome is one a file gives against itself; missing in
    both, invalid in both and too complex come only where A's own result is so; and B's
    result is A's, the same possible genotypes with the same indices, or the same reason.

    str_regions: a --str-regions BED; without one, the STR regions are found by
        segment_contig over every contig of the reference, which on a genome takes long.
    ploidy: what both files are read and reconciled with, 1 or 2.

    Raises AssertionError at the first row that fails, with the locus, the sample, the
    outcome and both results; ValueError when the headers name different samples; and what
    the stages raise.
    """
    for _ in _make_checked_rows(vcf_a, vcf_b, reference, str_regions=str_regions, ploidy=ploidy):
        pass


def _make_checked_rows(
    vcf_a: Path,
    vcf_b: Path,
    reference: IndexedFasta,
    *,
    str_regions: Path | None,
    ploidy: int,
) -> Iterator[OutcomeRow]:
    """The rows of check_respelling_outcomes, each yielded once checked, so that a test can
    see them."""
    sample_pairs = _make_sample_pairs(vcf_a, vcf_b)
    segments: Iterator[GenomeSegment]
    if str_regions is None:
        segments = itertools.chain.from_iterable(
            segment_contig(reference, contig) for contig in reference.contigs
        )
    else:
        segments = read_str_regions(str_regions, reference)
    loci = build_loci(
        parse_vcf(vcf_a, reference, ploidy=ploidy),
        parse_vcf(vcf_b, reference, ploidy=ploidy),
        segments,
        reference.contigs,
    )
    reconciled_loci = reconcile_loci(
        loci, reference, sample_pairs, ploidy_a=ploidy, ploidy_b=ploidy
    )
    for reconciled in reconciled_loci:
        rows = make_outcome_rows(
            reconciled, sample_pairs, kind=locus_kind(reconciled.locus), stratum=None
        )
        for sample_pair, row in zip(sample_pairs, rows, strict=True):
            _check_row(reconciled, sample_pair, row)
            yield row


def _make_sample_pairs(vcf_a: Path, vcf_b: Path) -> tuple[SamplePair, ...]:
    """Each sample paired with itself; ValueError when the headers name different samples."""
    samples = parse_vcf_samples(vcf_a)
    samples_b = parse_vcf_samples(vcf_b)
    if samples_b != samples:
        raise ValueError(
            f"{vcf_a} has the samples {', '.join(samples)}"
            f" and {vcf_b} the samples {', '.join(samples_b)}"
        )
    return tuple(
        SamplePair(name_a=name, index_a=index, name_b=name, index_b=index)
        for index, name in enumerate(samples)
    )


def _check_row(reconciled: ReconciledLocus, sample_pair: SamplePair, row: OutcomeRow) -> None:
    """Raise AssertionError when a row is not what file A gives against itself, or when the
    sample's result in B is not its result in A."""
    result_a = reconciled.results_a[sample_pair.index_a]
    result_b = reconciled.results_b[sample_pair.index_b]
    outcome = row.outcome
    problem = None
    if outcome not in _ALLOWED_OUTCOMES:
        problem = f"{outcome}, an outcome a file does not give against itself"
    elif outcome is Outcome.TOO_COMPLEX and result_a is not TOO_COMPLEX:
        problem = f"{outcome}, where file A's own result is not too complex"
    elif result_b != result_a:
        problem = f"{outcome}, and file B's result is not file A's"
    if problem is None:
        return
    locus = reconciled.locus
    alleles = reconciled.alleles
    raise AssertionError(
        f"{locus.contig}:{locus.start + 1}-{locus.end}, sample {sample_pair.name_a}: {problem};"
        f" file A gives {_format_result(result_a, alleles)}"
        f" and file B {_format_result(result_b, alleles)}"
    )


def _format_result(result: SampleResult, alleles: Sequence[str]) -> str:
    """A sample result for a message: its possible genotypes, each as its locus alleles'
    sequences, `-` for an empty one and `.` for a missing one, joined by `/`; or why it has
    none."""
    if result is TOO_COMPLEX:
        return "too complex"
    if isinstance(result, InvalidReason):
        return f"invalid, {result}"
    return " or ".join(
        "/".join("." if allele is MISSING else alleles[allele] or "-" for allele in genotype)
        for genotype in result
    )
