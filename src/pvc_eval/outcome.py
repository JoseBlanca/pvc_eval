"""The types every evaluator shares (doc/specs/reconciliation.md §7), in a module that holds
only them, as pvc_eval.locus holds Locus.

Public: Outcome, the outcome of a sample pair at a locus; InvalidReason, why a sample is
invalid at a locus; LocusGenotype, a locus genotype as indices into a locus's locus alleles;
OutcomeRow, one row of the table.
"""

from dataclasses import dataclass
from enum import StrEnum

from pvc_eval.kind import Kind
from pvc_eval.segments import GenomeSegment
from pvc_eval.vcf import MISSING


class Outcome(StrEnum):
    """The outcome of a sample pair at a locus, in the words the table writes
    (doc/spec.md §5.5)."""

    TOO_COMPLEX = "too complex"
    INVALID_IN_A = "invalid in A"
    INVALID_IN_B = "invalid in B"
    INVALID_IN_BOTH = "invalid in both"
    MISSING_IN_A = "missing in A"
    MISSING_IN_B = "missing in B"
    MISSING_IN_BOTH = "missing in both"
    REFERENCE_IN_BOTH = "reference in both"
    SAME_GENOTYPE = "same genotype"
    SAME_ALLELES = "same alleles, other genotype"
    ONLY_IN_A = "only in A"
    ONLY_IN_B = "only in B"
    DIFFERENT = "different"


class InvalidReason(StrEnum):
    """Why a sample is invalid at a locus, in the words the table writes (doc/spec.md §6)."""

    MALFORMED_GT = "malformed GT"
    PLOIDY_MISMATCH = "ploidy mismatch"
    CONTRADICTING_RECORDS = "contradicting records"


type LocusGenotype = tuple[int | MISSING, ...]
"""One locus allele per haplotype, as its index into the locus alleles, or MISSING; sorted,
with MISSING after every index (doc/specs/reconciliation.md §4.3)."""


@dataclass(frozen=True, slots=True, kw_only=True)
class OutcomeRow:
    """One row of the table (doc/spec.md §6): one sample pair at one locus, with the fields
    of doc/specs/reconciliation.md §6, which every evaluator fills.

    Positions are 0-based and half-open, as in Locus; the table writer makes them 1-based.
    A value the table writes as `.` is None, except stratum.

    contig, start, end: the locus.
    stratum: the locus's stratum as doc/specs/loci.md gives it, `.` for a locus no interval
        of the BED overlaps; None without --strata.
    str_regions: the STR regions the locus contains.
    alleles: the locus alleles, REF first; None from an evaluator that builds none.
    genotype_a, genotype_b: the reported genotype; in a row too complex or invalid, the
        first possible genotype of a sample that has one; None for a sample too complex or
        invalid, or from another evaluator.
    possible_a, possible_b: how many possible genotypes the sample has; 0 for an invalid
        sample; None for a sample too complex, or from another evaluator.
    quality_a, quality_b: the lowest QUAL among the file's records in the locus at which the
        sample's VCF genotype has an index other than 0, `*` included; None when no such
        record has a QUAL.
    gq_a, gq_b: the lowest GQ of the sample among the same records; None when none has one.
    records_a, records_b: how many of the file's records in the locus have an ALT.
    invalid_reason_a, invalid_reason_b: the sample's reason when it is invalid, whatever
        the outcome; None otherwise.
    """

    contig: str
    start: int
    end: int
    kind: Kind
    stratum: str | None
    str_regions: tuple[GenomeSegment, ...]
    sample_a: str
    sample_b: str
    outcome: Outcome
    alleles: tuple[str, ...] | None
    genotype_a: LocusGenotype | None
    genotype_b: LocusGenotype | None
    possible_a: int | None
    possible_b: int | None
    quality_a: float | None
    quality_b: float | None
    gq_a: int | None
    gq_b: int | None
    records_a: int
    records_b: int
    invalid_reason_a: InvalidReason | None
    invalid_reason_b: InvalidReason | None
