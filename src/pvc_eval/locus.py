"""A locus, what stage 3 of doc/spec.md §4 yields to the reconciliation: the type that
passes between the loci and the evaluators, fixed by doc/spec.md §4 and changed only there.

Positions are 0-based and half-open, as VcfRecord.start and VcfRecord.span are.

Public: Locus.
"""

from dataclasses import dataclass

from pvc_eval.segments import GenomeSegment
from pvc_eval.vcf import VcfRecord


@dataclass(frozen=True, slots=True)
class Locus:
    """A stretch of one contig, `[start, end)`, 0-based and half-open, with the VCF records
    of both files that fall in it (doc/spec.md §5.3).

    It covers the spans of all its records with an ALT and every STR region it holds, and
    so their REFs. Each file's records are in the order parse_vcf yields them, by span
    start; a record with no ALT is in every locus its REF overlaps, and its REF may reach
    beyond the locus. Each record's genotypes are in the order of its file's header.

    records_a, records_b: file A's records in the locus, and file B's.
    str_regions: the STR regions the locus contains, in position order, with no GENERIC
        segment.
    """

    contig: str
    start: int
    end: int
    records_a: tuple[VcfRecord, ...]
    records_b: tuple[VcfRecord, ...]
    str_regions: tuple[GenomeSegment, ...]
