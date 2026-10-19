"""A file against its decomposition, the property of doc/specs/reconciliation.md §9 and of
doc/spec.md §7.3 and §7.4: file B is file A with each VCF record that no edit of another
record of the same haplotype touches decomposed into records of one replaced base or one
indel, phased in a phase set of its own for each record of A, and the records that touch
kept whole. Every row is then same genotype, reference in both, missing in both, invalid in
both or too complex, each of the last three only where file A has its cause.

The decomposition is a function of these tests, make_decomposed_lines, checked first on
cases worked by hand.
"""

from collections.abc import Iterator, Sequence
from contextlib import ExitStack
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from tests.test_reconcile import CONTIG, CaseFiles, get_locus_bases, pair_samples

from pvc_eval.evaluate import make_outcome_rows
from pvc_eval.fasta import IndexedFasta
from pvc_eval.kind import Kind
from pvc_eval.locus import Locus
from pvc_eval.outcome import InvalidReason, Outcome
from pvc_eval.reconcile import PHASING_CAP, TOO_COMPLEX, reconcile
from pvc_eval.vcf import MALFORMED_GT, MISSING, OTHER_PLOIDY_GT, VcfGenotype, VcfRecord

_BASES = "ACGT"

# An edit, as doc/specs/reconciliation.md §4.2 finds it: the first position it replaces,
# 0-based, the bases of REF it replaces, none for an insertion, and the bases it puts there.
type _Edit = tuple[int, str, str]


def make_decomposed_lines(vcf_records: Sequence[VcfRecord], reference: IndexedFasta) -> list[str]:
    """The VCF lines of file B, `POS REF ALT GT...` as CaseFiles reads them, made from file
    A's records by decomposing each that nothing touches into records of one replaced base
    or one indel (doc/specs/reconciliation.md §9).

    A record is kept whole, as A wrote it, when an edit of an ALT a sample carries at it
    touches, overlaps or lies next to (doc/spec.md §2), an edit of an ALT the same sample
    carries at another record. The records of A are unphased, so any two such ALTs may be
    on one haplotype.

    Each ALT's edit is found as §4.2 finds it, the bases it shares with REF removed from
    the left and then from the right, and the replaced stretch and its new bases are
    aligned from the left: a SNP at each position where they differ, then the bases left
    over as an insertion after the stretch's last aligned position or a deletion of the
    positions left over, each with the base before it as its anchor. A piece two ALTs
    share is one record. Every haplotype carries the pieces of its own ALT, REF (0) at the
    others' and `.` where A has `.`, so a record of A that is `1/2` gives `1|0` at the
    pieces of ALT 1 alone. All the pieces of one record of A are phased, `GT:PS`, in a
    phase set of their own, numbered by the record's place from 1, so no block is added
    (doc/spec.md §5.4). The lines are sorted by POS, as parse_vcf asks.
    """
    edits = [
        [_find_edit(vcf_record, allele_index) for allele_index in range(1, len(vcf_record.alleles))]
        for vcf_record in vcf_records
    ]
    numbered_lines: list[tuple[int, str]] = []
    for place, vcf_record in enumerate(vcf_records):
        if _is_touched(vcf_records, edits, place):
            numbered_lines.append((vcf_record.start + 1, _format_whole_line(vcf_record)))
            continue
        alts_of_piece: dict[tuple[int, str, str], set[int]] = {}
        for allele_index, edit in enumerate(edits[place], start=1):
            for piece in _make_pieces(edit, reference):
                alts_of_piece.setdefault(piece, set()).add(allele_index)
        for (pos_1based, ref, alt), allele_indices in alts_of_piece.items():
            gts = [
                _format_piece_genotype(vcf_genotype, allele_indices, place + 1)
                for vcf_genotype in vcf_record.genotypes
            ]
            numbered_lines.append((pos_1based, f"{pos_1based} {ref} {alt} {' '.join(gts)}"))
    return [line for _, line in sorted(numbered_lines, key=lambda numbered: numbered[0])]


def _find_edit(vcf_record: VcfRecord, allele_index: int) -> _Edit:
    """The edit of one ALT, which the generated records never make `*` or REF."""
    ref, alt = vcf_record.alleles[0], vcf_record.alleles[allele_index]
    assert alt != "*", "the generated records have no `*`"
    num_prefix = _count_shared(ref, alt)
    ref, alt = ref[num_prefix:], alt[num_prefix:]
    num_suffix = _count_shared(ref[::-1], alt[::-1])
    return vcf_record.start + num_prefix, ref[: len(ref) - num_suffix], alt[: len(alt) - num_suffix]


def _is_touched(
    vcf_records: Sequence[VcfRecord], edits: Sequence[Sequence[_Edit]], place: int
) -> bool:
    """Whether an edit a sample carries at the record at place touches one it carries at
    another record."""
    for sample_index in range(len(vcf_records[place].genotypes)):
        carried = _get_carried_edits(vcf_records[place], edits[place], sample_index)
        for other_place, other_record in enumerate(vcf_records):
            if other_place == place:
                continue
            other_carried = _get_carried_edits(other_record, edits[other_place], sample_index)
            if any(_is_touching(edit, other) for edit in carried for other in other_carried):
                return True
    return False


def _get_carried_edits(
    vcf_record: VcfRecord, record_edits: Sequence[_Edit], sample_index: int
) -> list[_Edit]:
    """The edits of the ALTs a sample carries at a record, on either haplotype."""
    vcf_genotype = vcf_record.genotypes[sample_index]
    assert isinstance(vcf_genotype, VcfGenotype), "the generated GTs are all well formed"
    assert not vcf_genotype.phased, "the generated GTs are all unphased"
    return [
        record_edits[allele - 1]
        for allele in vcf_genotype.alleles
        if allele is not MISSING and allele != 0
    ]


def _is_touching(edit: _Edit, other: _Edit) -> bool:
    """Whether two edits overlap or lie next to each other: their replaced stretches,
    an insertion's being the point it goes in at, share or meet at a position."""
    start, ref, _ = edit
    other_start, other_ref, _ = other
    return start <= other_start + len(other_ref) and other_start <= start + len(ref)


def _format_whole_line(vcf_record: VcfRecord) -> str:
    """A record of A as it was written, with its GTs unphased."""
    gts = []
    for vcf_genotype in vcf_record.genotypes:
        assert isinstance(vcf_genotype, VcfGenotype), "the generated GTs are all well formed"
        gts.append(
            "/".join("." if allele is MISSING else str(allele) for allele in vcf_genotype.alleles)
        )
    ref, *alts = vcf_record.alleles
    return f"{vcf_record.start + 1} {ref} {','.join(alts)} {' '.join(gts)}"


def _make_pieces(edit: _Edit, reference: IndexedFasta) -> Iterator[tuple[int, str, str]]:
    """The pieces of one edit, each as (POS, REF, ALT), POS 1-based."""
    start, ref, alt = edit
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
        # A replacement of 2-3 and a SNP at 5 lie apart: both are decomposed.
        ("GACGTA", ["2 AC GT 0/1", "5 T G 0/1"], ["2 A G 0|1:1", "3 C T 0|1:1", "5 T G 0|1:2"]),
        # A replacement of 2-3 and a SNP at 4 lie next to each other: both are kept whole.
        ("GACGTA", ["2 AC GT 0/1", "4 G A 0/1"], ["2 AC GT 0/1", "4 G A 0/1"]),
        # The SNP at 4 is carried by no haplotype, so it touches nothing.
        ("GACGTA", ["2 AC GT 0/1", "4 G A 0/0"], ["2 A G 0|1:1", "3 C T 0|1:1", "4 G A 0|0:2"]),
        # Two records whose edits touch, carried by different samples, are decomposed.
        (
            "GACGTA",
            ["2 AC GT 0/1 0/0", "4 G A 0/0 1/1"],
            ["2 A G 0|1:1 0|0:1", "3 C T 0|1:1 0|0:1", "4 G A 0|0:2 1|1:2"],
        ),
        # An insertion inside the stretch a deletion-and-SNP replaces, in a homopolymer:
        # doc/specs/reconciliation.md §9's first example, kept whole.
        ("TAAC", ["3 AC G,AA 0/1", "3 A G,AA 1/2"], ["3 A G,AA 1/2", "3 AC G,AA 0/1"]),
    ],
    ids=[
        "complex",
        "multi-alt-1/2",
        "missing-./1",
        "phase-set-per-record",
        "apart-decomposed",
        "next-to-kept-whole",
        "touching-but-not-carried",
        "touching-on-other-samples",
        "touching-insertion-kept-whole",
    ],
)
def test_make_decomposed_lines_decomposes_each_record_nothing_touches(
    case_files: CaseFiles,
    bases: str,
    lines_a: Sequence[str],
    expected_lines_b: Sequence[str],
) -> None:
    reference = case_files.write_reference(1, bases)
    num_samples = len(lines_a[0].split()) - 3
    records_a = case_files.write_and_parse_records(reference, lines_a, num_samples=num_samples)
    assert make_decomposed_lines(records_a, reference) == expected_lines_b


@st.composite
def _make_case(draw: st.DrawFn) -> tuple[str, list[str], int, int]:
    """A reference of 4 to 12 bases, from position 1; 1 to 4 VCF lines of file A, sorted by
    POS, each with 1 or 2 ALTs and an unphased GT per sample; the number of samples, 1 or 2;
    and a phasing cap. Every POS is 2 or later, so the base before every edit, the anchor
    of its pieces, is one of the reference's."""
    bases = draw(st.text(_BASES, min_size=4, max_size=12))
    num_samples = draw(st.integers(1, 2))
    lines = []
    for _ in range(draw(st.integers(1, 4))):
        pos = draw(st.integers(2, len(bases)))
        ref_length = draw(st.integers(1, min(4, len(bases) - pos + 1)))
        ref = bases[pos - 1 : pos - 1 + ref_length]
        alt_strategy = st.one_of(
            st.text(_BASES, min_size=1, max_size=4),
            st.text(_BASES, max_size=3).map(lambda tail, anchor=ref[0]: anchor + tail),
        ).filter(lambda alt, ref=ref: alt != ref)
        alts = draw(st.lists(alt_strategy, min_size=1, max_size=2, unique=True))
        allele = st.sampled_from([".", *map(str, range(len(alts) + 1))])
        gts = [f"{draw(allele)}/{draw(allele)}" for _ in range(num_samples)]
        lines.append((pos, f"{pos} {ref} {','.join(alts)} {' '.join(gts)}"))
    lines.sort(key=lambda numbered: numbered[0])
    phasing_cap = draw(st.sampled_from([1, 2, PHASING_CAP]))
    return bases, [line for _, line in lines], num_samples, phasing_cap


_ALLOWED = {
    Outcome.SAME_GENOTYPE,
    Outcome.REFERENCE_IN_BOTH,
    Outcome.MISSING_IN_BOTH,
    Outcome.INVALID_IN_BOTH,
    Outcome.TOO_COMPLEX,
}


@settings(max_examples=300, deadline=None)
@given(case=_make_case())
def test_a_file_against_its_decomposition_gives_only_the_outcomes_of_itself(
    tmp_path_factory: pytest.TempPathFactory, case: tuple[str, list[str], int, int]
) -> None:
    bases, lines_a, num_samples, phasing_cap = case
    with ExitStack() as stack:
        case_files = CaseFiles(tmp_path_factory.mktemp("case"), stack)
        reference = case_files.write_reference(1, bases)
        records_a = case_files.write_and_parse_records(reference, lines_a, num_samples=num_samples)
        lines_b = make_decomposed_lines(records_a, reference)
        records_b = case_files.write_and_parse_records(reference, lines_b, num_samples=num_samples)
        spans = [vcf_record.span for vcf_record in (*records_a, *records_b)]
        locus = Locus(
            CONTIG,
            min(start for start, _ in spans),
            max(end for _, end in spans),
            records_a,
            records_b,
            (),
        )
        sample_pairs = pair_samples(num_samples)
        reconciled = reconcile(
            locus, get_locus_bases(reference, locus), sample_pairs, phasing_cap=phasing_cap
        )
    rows = make_outcome_rows(reconciled, sample_pairs, kind=Kind.SNP, stratum=None)
    for sample_pair, row in zip(sample_pairs, rows, strict=True):
        assert row.outcome in _ALLOWED, (lines_a, lines_b, row)
        # A's own sample result, which does not depend on file B: what A gives against itself.
        result_a = reconciled.results_a[sample_pair.index_a]
        if row.outcome is Outcome.TOO_COMPLEX:
            assert result_a is TOO_COMPLEX, (lines_a, lines_b, row)
        elif row.outcome is Outcome.INVALID_IN_BOTH:
            assert isinstance(result_a, InvalidReason), (lines_a, lines_b, row)
        elif row.outcome is Outcome.MISSING_IN_BOTH:
            assert isinstance(result_a, tuple), (lines_a, lines_b, row)
            assert all(any(allele is MISSING for allele in genotype) for genotype in result_a), (
                lines_a,
                lines_b,
                row,
            )
