"""The kind of a locus (doc/spec.md §5.6), as doc/specs/loci.md gives it.

Public: Kind, the four kinds, in the words the table writes; locus_kind, the kind of a
locus.
"""

from collections.abc import Iterator
from enum import Enum

from pvc_eval.locus import Locus
from pvc_eval.vcf import trim_alt


class Kind(Enum):
    """The kind of a locus, read from its VCF records and its STR regions (doc/spec.md §5.6)."""

    SNP = "SNP"
    INDEL = "indel"
    STR = "STR"
    COMPLEX = "complex"


def locus_kind(locus: Locus) -> Kind:
    """The kind of a locus (doc/spec.md §5.6, doc/specs/loci.md "Kind").

    A locus holding an STR region is STR. Otherwise the kind is read from every ALT of
    every record with an ALT of both files, each trimmed against its REF, first at the left
    and then at the right, as an edit is (§5.4); `*` alleles and an ALT equal to its REF
    are left out. Every ALT a substitution of one base gives SNP, every ALT an insertion or
    a deletion gives indel, and anything else complex, also a locus with no ALT left.
    """
    if locus.str_regions:
        return Kind.STR
    edit_kinds = set(_get_edit_kinds(locus))
    if edit_kinds == {Kind.SNP}:
        return Kind.SNP
    if edit_kinds == {Kind.INDEL}:
        return Kind.INDEL
    # §5.6's "every ALT" would make a locus with no ALT left both SNP and indel; the loci
    # spec decides it is complex.
    return Kind.COMPLEX


def _get_edit_kinds(locus: Locus) -> Iterator[Kind]:
    """The kind of each edit of the locus: SNP, INDEL or COMPLEX."""
    for vcf_record in (*locus.records_a, *locus.records_b):
        ref = vcf_record.alleles[0]
        for alt in vcf_record.alleles[1:]:
            if alt in ("*", ref):
                continue
            yield _get_edit_kind(ref, alt)


def _get_edit_kind(ref: str, alt: str) -> Kind:
    """The kind of the edit ALT makes on REF, once the bases they share are trimmed."""
    _, ref, alt = trim_alt(ref, alt)
    if len(ref) == 1 and len(alt) == 1:
        return Kind.SNP
    if not ref or not alt:
        return Kind.INDEL
    return Kind.COMPLEX
