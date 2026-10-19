"""The decomposition of doc/specs/reconciliation.md §9's property and of doc/spec.md §7.4:
file B made from file A by decomposing each VCF record into records of one replaced base or
one indel, phased in a phase set of their own for each record of A. It is a function of the
tests, make_decomposed_lines, checked here on cases worked by hand.
"""

from collections.abc import Iterator, Sequence
from contextlib import ExitStack
from pathlib import Path

import pytest
from tests.test_reconcile import CONTIG, CaseFiles

from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import MALFORMED_GT, MISSING, OTHER_PLOIDY_GT, VcfGenotype, VcfRecord


def make_decomposed_lines(vcf_records: Sequence[VcfRecord], reference: IndexedFasta) -> list[str]:
    """The VCF lines of file B, `POS REF ALT GT:PS...` as CaseFiles reads them, made from
    file A's records by decomposing each into records of one replaced base or one indel.

    Each ALT's edit is found as doc/specs/reconciliation.md §4.2 finds it, the bases it
    shares with REF removed from the left and then from the right, and the replaced
    stretch and its new bases are aligned from the left: a SNP at each position where they
    differ, then the bases left over as an insertion after the stretch's last aligned
    position or a deletion of the positions left over, each with the base before it as
    its anchor. A piece two ALTs share is one record. Every haplotype carries the pieces
    of its own ALT, REF (0) at the others' and `.` where A has `.`, so a record of A that
    is `1/2` gives `1|0` at the pieces of ALT 1 alone. All the pieces of one record of A
    are phased in a phase set of their own, numbered from 1, so no block is added
    (doc/spec.md §5.4). The lines are sorted by POS, as parse_vcf asks.
    """
    numbered_lines: list[tuple[int, str]] = []
    for phase_set, vcf_record in enumerate(vcf_records, start=1):
        alts_of_piece: dict[tuple[int, str, str], set[int]] = {}
        for allele_index in range(1, len(vcf_record.alleles)):
            for piece in _make_pieces(vcf_record, allele_index, reference):
                alts_of_piece.setdefault(piece, set()).add(allele_index)
        for (pos_1based, ref, alt), allele_indices in alts_of_piece.items():
            gts = [
                _format_piece_genotype(vcf_genotype, allele_indices, phase_set)
                for vcf_genotype in vcf_record.genotypes
            ]
            numbered_lines.append((pos_1based, f"{pos_1based} {ref} {alt} {' '.join(gts)}"))
    return [line for _, line in sorted(numbered_lines, key=lambda numbered: numbered[0])]


def _make_pieces(
    vcf_record: VcfRecord, allele_index: int, reference: IndexedFasta
) -> Iterator[tuple[int, str, str]]:
    """The pieces of one ALT's edit, each as (POS, REF, ALT), POS 1-based."""
    ref, alt = vcf_record.alleles[0], vcf_record.alleles[allele_index]
    assert alt != "*", "the generated records have no `*`"
    num_prefix = _count_shared(ref, alt)
    ref, alt = ref[num_prefix:], alt[num_prefix:]
    num_suffix = _count_shared(ref[::-1], alt[::-1])
    ref, alt = ref[: len(ref) - num_suffix], alt[: len(alt) - num_suffix]
    start = vcf_record.start + num_prefix  # 0-based, the first replaced position
    num_aligned = min(len(ref), len(alt))
    for offset in range(num_aligned):
        if ref[offset] != alt[offset]:
            yield start + offset + 1, ref[offset], alt[offset]
    # The anchor is the last aligned position, or the base before the edit when none is.
    anchor = start + num_aligned - 1  # 0-based
    anchor_base = reference.get_bases(CONTIG, anchor, anchor + 1)
    if len(alt) > num_aligned:
        yield anchor + 1, anchor_base, anchor_base + alt[num_aligned:]
    elif len(ref) > num_aligned:
        yield anchor + 1, anchor_base + ref[num_aligned:], anchor_base


def _count_shared(first: str, second: str) -> int:
    """How many bases two sequences share at their start."""
    for index, (base, other_base) in enumerate(zip(first, second, strict=False)):
        if base != other_base:
            return index
    return min(len(first), len(second))


def _format_piece_genotype(
    vcf_genotype: VcfGenotype | MALFORMED_GT | OTHER_PLOIDY_GT,
    allele_indices: set[int],
    phase_set: int,
) -> str:
    """A sample's GT:PS at a piece carried by the ALTs allele_indices."""
    assert isinstance(vcf_genotype, VcfGenotype), "the generated GTs are all well formed"
    alleles = [
        "." if allele is MISSING else str(int(allele in allele_indices))
        for allele in vcf_genotype.alleles
    ]
    return f"{'|'.join(alleles)}:{phase_set}"


@pytest.fixture
def case_files(tmp_path: Path) -> Iterator[CaseFiles]:
    with ExitStack() as stack:
        yield CaseFiles(tmp_path, stack)


@pytest.mark.parametrize(
    ("bases", "lines_a", "expected_lines_b"),
    [
        # ACGT -> TCGAA replaces 2-5: SNPs at 2 and 5, and the A left over goes in after 5.
        (
            "GACGTA",
            ["2 ACGT TCGAA 0/1"],
            ["2 A T 0|1:1", "5 T A 0|1:1", "5 T TA 0|1:1"],
        ),
        # TCA replaces 2-4, TCG only 2 (its G is REF's): the SNP at 2 is one record, carried
        # by both haplotypes; the SNP at 4 by the first alone.
        ("GACGT", ["2 ACG TCA,TCG 1/2"], ["2 A T 1|1:1", "4 G A 1|0:1"]),
        # AC -> GT, with the first haplotype missing.
        ("GACT", ["2 AC GT ./1"], ["2 A G .|1:1", "3 C T .|1:1"]),
        # Each record of A has a phase set of its own, and REF where a haplotype carries it.
        (
            "GACGTA",
            ["2 A AGG 0/1", "4 GTA G 1/0"],
            ["2 A AGG 0|1:1", "4 GTA G 1|0:2"],
        ),
    ],
    ids=["complex", "multi-alt-1/2", "missing-./1", "phase-set-per-record"],
)
def test_make_decomposed_lines_decomposes_each_record(
    case_files: CaseFiles,
    bases: str,
    lines_a: Sequence[str],
    expected_lines_b: Sequence[str],
) -> None:
    reference = case_files.write_reference(1, bases)
    records_a = case_files.write_and_parse_records(reference, lines_a)
    assert make_decomposed_lines(records_a, reference) == expected_lines_b
