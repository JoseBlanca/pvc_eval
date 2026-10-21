"""A file against its decomposition, the property of doc/specs/reconciliation.md §9 and of
doc/spec.md §7.3 and §7.4: file B is file A respelled by respell_by_decomposing
(doc/specs/respell.md §3), each VCF record that no edit of another record of the same
haplotype touches decomposed into records of one replaced base or one indel in the record's
phase unit, and the others written as read. Every row is then same genotype, reference in
both, missing in both, invalid in both or too complex, each of the last three only where
file A has its cause.
"""

from contextlib import ExitStack

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from tests.conftest import CONTIG, CaseFiles, get_locus_bases, pair_samples
from tests.respell import respell_by_decomposing
from tests.respell_check import check_outcome_row

from pvc_eval.evaluate import make_outcome_rows
from pvc_eval.kind import Kind
from pvc_eval.locus import Locus
from pvc_eval.reconcile import PHASING_CAP, reconcile
from pvc_eval.vcf import parse_vcf

_BASES = "ACGT"


def _make_allele_strategy(num_alts: int) -> st.SearchStrategy[str]:
    """An allele of a GT: `.` one draw in ten, an index of the record's alleles otherwise."""
    return st.integers(0, 9).flatmap(
        lambda draw_number: st.just(".") if draw_number == 0 else st.integers(0, num_alts).map(str)
    )


@st.composite
def _make_case(draw: st.DrawFn) -> tuple[str, list[str], int, int]:
    """A reference of 4 to 12 bases, from position 1; 1 to 4 VCF lines of file A, sorted by
    POS, each with 1 or 2 ALTs and an unphased GT per sample, `.` in about one allele of
    ten; the number of samples, 1 or 2; and a phasing cap, 4096 in three cases of five, 1
    or 2 in the others, so that most rows compare sequences. POS may be 1, so an
    indel may start at the contig's first base."""
    bases = draw(st.text(_BASES, min_size=4, max_size=12))
    num_samples = draw(st.integers(1, 2))
    lines = []
    for _ in range(draw(st.integers(1, 4))):
        pos = draw(st.integers(1, len(bases)))
        ref_length = draw(st.integers(1, min(4, len(bases) - pos + 1)))
        ref = bases[pos - 1 : pos - 1 + ref_length]
        alt_strategy = st.one_of(
            st.text(_BASES, min_size=1, max_size=4),
            st.text(_BASES, max_size=3).map(lambda tail, anchor=ref[0]: anchor + tail),
        ).filter(lambda alt, ref=ref: alt != ref)
        alts = draw(st.lists(alt_strategy, min_size=1, max_size=2, unique=True))
        allele = _make_allele_strategy(len(alts))
        gts = [f"{draw(allele)}/{draw(allele)}" for _ in range(num_samples)]
        lines.append((pos, f"{pos} {ref} {','.join(alts)} {' '.join(gts)}"))
    lines.sort(key=lambda numbered: numbered[0])
    phasing_cap = draw(st.sampled_from([PHASING_CAP, PHASING_CAP, PHASING_CAP, 1, 2]))
    return bases, [line for _, line in lines], num_samples, phasing_cap


@settings(max_examples=300, deadline=None)
@given(case=_make_case())
def test_a_file_against_its_decomposition_gives_only_the_outcomes_of_itself(
    tmp_path_factory: pytest.TempPathFactory, case: tuple[str, list[str], int, int]
) -> None:
    bases, lines_a, num_samples, phasing_cap = case
    with ExitStack() as stack:
        case_files = CaseFiles(tmp_path_factory.mktemp("case"), stack)
        reference = case_files.write_reference(1, bases)
        vcf_a = case_files.write_vcf(lines_a, num_samples=num_samples)
        vcf_b = vcf_a.with_name("b.vcf")
        respell_by_decomposing(vcf_a, reference, vcf_b)
        records_a = tuple(parse_vcf(vcf_a, reference))
        records_b = tuple(parse_vcf(vcf_b, reference))
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
    # B has A's blocks and its haplotypes A's sequences, so each row is one a file gives
    # against itself and each sample's result is A's (doc/specs/respell.md §9.3).
    rows = make_outcome_rows(reconciled, sample_pairs, kind=Kind.SNP, stratum=None)
    for sample_pair, row in zip(sample_pairs, rows, strict=True):
        check_outcome_row(reconciled, sample_pair, row)
