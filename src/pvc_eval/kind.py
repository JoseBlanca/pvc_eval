"""The kind of a locus (doc/spec.md §5.6), as doc/specs/loci.md gives it.

Public: Kind, the four kinds, in the words the table writes.
"""

from enum import Enum


class Kind(Enum):
    """The kind of a locus, read from its VCF records and its STR regions (doc/spec.md §5.6)."""

    SNP = "SNP"
    INDEL = "indel"
    STR = "STR"
    COMPLEX = "complex"
