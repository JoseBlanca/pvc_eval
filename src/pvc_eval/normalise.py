"""Normalising a VCF record: trim, then left-align (doc/spec.md §5.2).

A VCF record is never split or decomposed. Positions are 0-based and half-open.
"""

from collections.abc import Sequence
from dataclasses import replace

from pvc_eval.vcf import SPANNING_DELETION, VcfRecord


def normalise(record: VcfRecord, contig_sequence: str) -> VcfRecord:
    """The record trimmed and, when it is an indel, left-aligned.

    contig_sequence: the record's contig in the reference, in upper case.

    A record whose REF differs from the reference is returned unchanged and marked
    ref_mismatch, since its alleles cannot be placed on the reference (§5.2).
    """
    if contig_sequence[record.start : record.end] != record.alleles[0].upper():
        return replace(record, ref_mismatch=True)
    return left_align(trim(record), contig_sequence)


def trim(record: VcfRecord) -> VcfRecord:
    """The record without the bases all its VCF alleles share, at the right end, then
    at the left, keeping at least one base in each (the VCF anchor convention).

    A `*` allele has no bases and takes no part. A record with fewer than two other
    alleles is returned as it is.
    """
    sequences = [allele for allele in record.alleles if allele != SPANNING_DELETION]
    if len(sequences) < 2:
        return record
    start, trimmed = _trim(record.start, sequences)
    return replace(record, start=start, alleles=_put_back(record.alleles, trimmed))


def left_align(record: VcfRecord, contig_sequence: str) -> VcfRecord:
    """The trimmed record slid left while the base it uncovers equals the base it drops,
    as `bcftools norm` and vt do.

    A record moves only when every ALT is an indel and every ALT moves by the same
    number of bases; otherwise it stays where it is: a SNP or a complex ALT cannot move,
    and ALTs that would move apart cannot be written as one record (§5.2).
    """
    ref = record.alleles[0]
    alts = [allele for allele in record.alleles[1:] if allele != SPANNING_DELETION]
    if not alts:
        return record
    shifts = {_indel_shift(record.start, ref, alt, contig_sequence) for alt in alts}
    if len(shifts) != 1:
        return record
    shift = shifts.pop()
    if not shift:
        return record
    new_start = record.start - shift
    flank = contig_sequence[new_start : record.start]
    moved = [
        allele if allele == SPANNING_DELETION else (flank + allele)[: len(allele)]
        for allele in record.alleles
    ]
    return trim(replace(record, start=new_start, alleles=tuple(moved)))


def _indel_shift(start: int, ref: str, alt: str, contig_sequence: str) -> int | None:
    """How many bases left the indel ref -> alt can slide; None for a SNP or a complex
    ALT, which cannot move."""
    start, (ref, alt) = _trim(start, [ref, alt])
    is_indel = len(ref) != len(alt) and min(len(ref), len(alt)) == 1 and ref[0] == alt[0]
    if not is_indel:
        return None
    event = max(ref, alt, key=len)[1:]  # the inserted or deleted bases, after the anchor
    event_start = start + 1
    # Slide while an anchor base stays to the left (event_start - 2 >= 0).
    while event_start > 1 and contig_sequence[event_start - 1] == event[-1]:
        event = contig_sequence[event_start - 1] + event[:-1]
        event_start -= 1
    return start - (event_start - 1)


def _trim(start: int, alleles: Sequence[str]) -> tuple[int, list[str]]:
    trimmed = list(alleles)
    while min(map(len, trimmed)) > 1 and len({allele[-1] for allele in trimmed}) == 1:
        trimmed = [allele[:-1] for allele in trimmed]
    while min(map(len, trimmed)) > 1 and len({allele[0] for allele in trimmed}) == 1:
        trimmed = [allele[1:] for allele in trimmed]
        start += 1
    return start, trimmed


def _put_back(alleles: tuple[str, ...], trimmed: list[str]) -> tuple[str, ...]:
    """The alleles with each sequence replaced by its trimmed form, `*` left in place."""
    replacements = iter(trimmed)
    return tuple(
        allele if allele == SPANNING_DELETION else next(replacements) for allele in alleles
    )
