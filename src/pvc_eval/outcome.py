"""The types every evaluator shares (doc/specs/reconciliation.md §7), in a module that holds
only them, as pvc_eval.locus holds Locus.

Public: InvalidReason, why a sample is invalid at a locus; LocusGenotype, a locus genotype
as indices into a locus's locus alleles.
"""

from enum import StrEnum

from pvc_eval.vcf import MISSING


class InvalidReason(StrEnum):
    """Why a sample is invalid at a locus, in the words the table writes (doc/spec.md §6)."""

    REF_MISMATCH = "REF mismatch"
    MALFORMED_GT = "malformed GT"
    PLOIDY_MISMATCH = "ploidy mismatch"
    CONTRADICTING_RECORDS = "contradicting records"


type LocusGenotype = tuple[int | MISSING, ...]
"""One locus allele per haplotype, as its index into the locus alleles, or MISSING; sorted,
with MISSING after every index (doc/specs/reconciliation.md §4.3)."""
