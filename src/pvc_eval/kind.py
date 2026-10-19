"""The kind of a locus (doc/spec.md §5.6), as doc/specs/loci.md gives it.

Public: Kind, the four kinds, in the words the table writes; locus_kind, the kind of a
locus.
"""

from collections.abc import Iterator
from enum import Enum

from pvc_eval.locus import Locus


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
    edit_kinds = set(_find_edit_kinds(locus))
    if edit_kinds == {Kind.SNP}:
        return Kind.SNP
    if edit_kinds == {Kind.INDEL}:
        return Kind.INDEL
    # §5.6's "every ALT" would make a locus with no ALT left both SNP and indel; the loci
    # spec decides it is complex.
    return Kind.COMPLEX


def _find_edit_kinds(locus: Locus) -> Iterator[Kind]:
    """The kind of each edit of the locus: SNP, INDEL or COMPLEX."""
    for vcf_record in (*locus.records_a, *locus.records_b):
        ref = vcf_record.alleles[0]
        for alt in vcf_record.alleles[1:]:
            if alt in ("*", ref):
                continue
            yield _find_edit_kind(ref, alt)


def _find_edit_kind(ref: str, alt: str) -> Kind:
    """The kind of the edit ALT makes on REF, once the bases they share are trimmed."""
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref, alt = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(ref), len(alt)) and ref[-1 - suffix] == alt[-1 - suffix]:
        suffix += 1
    ref, alt = ref[: len(ref) - suffix], alt[: len(alt) - suffix]
    if len(ref) == 1 and len(alt) == 1:
        return Kind.SNP
    if not ref or not alt:
        return Kind.INDEL
    return Kind.COMPLEX
