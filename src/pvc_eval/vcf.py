"""Reading the VCF records of one file: stage 2 of doc/spec.md §4, with the checks of §5.2.

Positions are 0-based and half-open inside pvc_eval. The 1-based POS of a VCF line is
converted here, and nowhere else.
"""

import gzip
import logging
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from pvc_eval.errors import PvcEvalInputError

logger = logging.getLogger(__name__)

PLOIDY = 2  # doc/spec.md §2: a global constant

MISSING = sentinel("MISSING")
"""A haplotype the VCF genotype says nothing about, written `.` (doc/spec.md §2)."""

MALFORMED_GT = sentinel("MALFORMED_GT")
"""A GT with a number of alleles other than the ploidy, or an allele index the VCF
record lacks: the sample is invalid at the locus, and the run goes on (§5.4)."""

SPANNING_DELETION = "*"
GZIP_MAGIC = b"\x1f\x8b"
NUM_FIXED_COLUMNS = 9  # CHROM to FORMAT


@dataclass(frozen=True, slots=True)
class VcfGenotype:
    """One sample's GT at one VCF record, as written.

    alleles: one entry per haplotype, the index of a VCF allele (0 is REF), or MISSING.
    phased: the GT was written with `|`.
    phase_set: the PS of the sample, or None when the record gives none, which VCF
        reads as one phase set per sample and contig (§5.4).
    """

    alleles: tuple[int | MISSING, ...]
    phased: bool
    phase_set: str | None


@dataclass(frozen=True, slots=True, kw_only=True)
class VcfRecord:
    """One line of a VCF file.

    start: the first position REF covers, 0-based; `end` is one past the last.
    alleles: REF first, then the ALTs, as written; `*` is a spanning deletion.
    qual: QUAL, or None for `.`.
    genotypes: one per sample, in the order of the header.
    genotype_qualities: the GQ of each sample, None where absent.
    line_number: the 1-based line of the file, for messages.
    ref_mismatch: REF differs from the reference (§5.2), set by normalising.
    """

    contig: str
    start: int
    alleles: tuple[str, ...]
    qual: float | None
    genotypes: tuple[VcfGenotype | MALFORMED_GT, ...]
    genotype_qualities: tuple[int | None, ...]
    line_number: int
    ref_mismatch: bool = False

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"start must be 0 or more, not {self.start}")
        if not self.alleles or not self.alleles[0]:
            raise ValueError("a VCF record needs a REF of at least one base")

    @property
    def end(self) -> int:
        return self.start + len(self.alleles[0])


def read_samples(path: Path) -> tuple[str, ...]:
    """The sample names of a VCF file, in the order of its columns."""
    with _open_text(path) as lines:
        for line in lines:
            if line.startswith("#CHROM"):
                return tuple(line.rstrip("\n").split("\t")[NUM_FIXED_COLUMNS:])
            if not line.startswith("##"):
                break
    raise PvcEvalInputError(f"{path.name} has no #CHROM header line; is it a VCF file?")


def read_vcf_records(path: Path, contigs: Sequence[str]) -> Iterator[VcfRecord]:
    """The VCF records of a file, in the file's order, plain or bgzipped.

    contigs: the contig names of the reference, in the order of its `.fai`.

    Raises PvcEvalInputError, when the stream reaches it, for a VCF record on a contig
    the reference lacks, for records out of the `.fai`'s contig order or out of position
    order, and for a line that does not hold a VCF record. A record with a symbolic
    allele is skipped, and the number skipped is logged as a warning at the end (§9).
    """
    contig_rank = {contig: rank for rank, contig in enumerate(contigs)}
    previous: tuple[int, int, str, int] | None = None  # rank, pos, contig, line
    num_symbolic = 0
    with _open_text(path) as lines:
        for line_number, line in enumerate(lines, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            contig, pos_1based = _read_position(fields, path, line_number)
            if contig not in contig_rank:
                raise PvcEvalInputError(
                    f"{path.name}, line {line_number}: contig {contig} is not in the reference."
                    " Compare against the reference the VCF was called on."
                )
            current = (contig_rank[contig], pos_1based, contig, line_number)
            if previous is not None and current[:2] < previous[:2]:
                raise PvcEvalInputError(
                    f"{path.name} is not sorted: {contig}:{pos_1based} (line {line_number})"
                    f" comes after {previous[2]}:{previous[1]} (line {previous[3]})."
                    " Sort it in the contig order of the reference's .fai and by position,"
                    " and run again."
                )
            previous = current
            alts = _read_alts(fields[4])
            if any(_is_symbolic(alt) for alt in alts):
                num_symbolic += 1
                continue
            yield _parse_record(fields, contig, pos_1based - 1, alts, line_number)
    if num_symbolic:
        logger.warning(
            "%s: %d VCF records with symbolic alleles (<DEL>, breakends) were skipped;"
            " pvc_eval does not compare structural variants.",
            path.name,
            num_symbolic,
        )


def parse_genotype(gt: str, num_alleles: int, phase_set: str | None) -> VcfGenotype | MALFORMED_GT:
    """A GT field as a VcfGenotype, or MALFORMED_GT when it cannot be one (§5.4).

    num_alleles: how many VCF alleles the record has, REF included.
    """
    parts = gt.replace("|", "/").split("/")
    if len(parts) != PLOIDY:
        return MALFORMED_GT
    alleles: list[int | MISSING] = []
    for part in parts:
        if part == ".":
            alleles.append(MISSING)
        elif part.isascii() and part.isdecimal() and int(part) < num_alleles:
            alleles.append(int(part))
        else:
            return MALFORMED_GT
    return VcfGenotype(alleles=tuple(alleles), phased="|" in gt, phase_set=phase_set)


def _open_text(path: Path) -> TextIO:
    with path.open("rb") as raw:
        is_gzipped = raw.read(len(GZIP_MAGIC)) == GZIP_MAGIC
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _read_position(fields: list[str], path: Path, line_number: int) -> tuple[str, int]:
    if len(fields) < NUM_FIXED_COLUMNS - 1 or not fields[1].isdecimal():
        raise PvcEvalInputError(
            f"{path.name}, line {line_number}: not a VCF record; it needs at least eight"
            " tab-separated columns and a whole number in POS."
        )
    return fields[0], int(fields[1])


def _read_alts(alt_field: str) -> tuple[str, ...]:
    if alt_field == ".":
        return ()
    return tuple(alt_field.split(","))


def _is_symbolic(alt: str) -> bool:
    return alt.startswith("<") or "[" in alt or "]" in alt


def _parse_record(
    fields: list[str], contig: str, start: int, alts: tuple[str, ...], line_number: int
) -> VcfRecord:
    alleles = (fields[3], *alts)
    genotypes: list[VcfGenotype | MALFORMED_GT] = []
    genotype_qualities: list[int | None] = []
    if len(fields) > NUM_FIXED_COLUMNS:
        keys = fields[8].split(":")
        for sample_field in fields[NUM_FIXED_COLUMNS:]:
            # VCF lets trailing sample fields be dropped, so a key may have no value.
            values = dict(zip(keys, sample_field.split(":"), strict=False))
            phase_set = values.get("PS")
            genotypes.append(
                parse_genotype(
                    values.get("GT", "./."),
                    len(alleles),
                    None if phase_set in (None, ".") else phase_set,
                )
            )
            genotype_qualities.append(_int_or_none(values.get("GQ")))
    return VcfRecord(
        contig=contig,
        start=start,
        alleles=alleles,
        qual=_float_or_none(fields[5]),
        genotypes=tuple(genotypes),
        genotype_qualities=tuple(genotype_qualities),
        line_number=line_number,
    )


# QUAL and GQ only annotate a row of the table (§6), where an absent value is written `.`;
# a value that is not a number is taken as absent rather than stopping the run.
def _float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _int_or_none(value: str | None) -> int | None:
    if value is None or not value.isascii() or not value.isdecimal():
        return None
    return int(value)
