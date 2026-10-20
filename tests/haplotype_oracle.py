"""The oracle of doc/specs/respell.md §8: whether two VCF files of the same samples give every
sample the same haplotypes over one reference, worked out again without pvc_eval's
reconciliation, so that a re-spelling is shown right without the code that the check of
doc/spec.md §7.4 tests. It reads the files with parse_vcf and imports nothing else of
pvc_eval but the reference reader.

The oracle compares, for each sample, its possible genotypes over each window, a stretch made
by joining the spans of the records with an ALT of both files while they touch; and over
each pair of windows that a block joins, the possible genotypes over both together.

Positions are 0-based and half-open, but in the messages, which give windows 1-based and
inclusive as a VCF does.

Public: check_same_haplotypes, HaplotypeComparison, ORACLE_CAP.
"""

import bisect
import itertools
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import (
    DEFAULT_PLOIDY,
    MALFORMED_GT,
    MISSING,
    OTHER_PLOIDY_GT,
    VcfGenotype,
    VcfRecord,
    parse_vcf,
    parse_vcf_samples,
)

# The most phase assignments a window, or a pair, is compared with in either file:
# doc/specs/reconciliation.md §4.1's default phasing cap.
ORACLE_CAP = 4096

_SPANNING_DELETION = "*"
_MALFORMED_GT_REASON = "malformed GT"
_PLOIDY_MISMATCH_REASON = "ploidy mismatch"
_MAX_LISTED = 10  # the genotypes a message lists of each file, so that it stays readable

_IMPOSSIBLE = sentinel("_IMPOSSIBLE")
"""A haplotype whose edits collide, or that carries a `*` no edit covers (§8)."""
_OVER_CAP = sentinel("_OVER_CAP")
"""A window or a pair with more than ORACLE_CAP phase assignments."""

# What one haplotype gives over one window: a sequence, or MISSING.
type _HaplotypeValue = str | MISSING
# A possible genotype over one window or a pair of them: for each haplotype its values, one
# per window, the haplotypes sorted so that the tuple stands for their multiset.
type _Genotype = tuple[tuple[_HaplotypeValue, ...], ...]
# What a sample gives over a window or a pair: its possible genotypes, the reason it is
# invalid there, or _OVER_CAP.
type _SampleResult = frozenset[_Genotype] | str | _OVER_CAP
# A phase unit of a sample (doc/specs/reconciliation.md §4.1): the phase set of a phased GT,
# None for a phased GT with no PS, or the line of an unphased one, a unit of its own.
type _UnitKey = tuple[str, str | int | None]


@dataclass(frozen=True, slots=True, kw_only=True)
class HaplotypeComparison:
    """What check_same_haplotypes compared.

    num_windows: windows compared, over every sample.
    num_linked_pairs: pairs of windows compared together, over every sample.
    num_over_cap: windows and pairs not compared, over ORACLE_CAP in either file.
    """

    num_windows: int
    num_linked_pairs: int
    num_over_cap: int


@dataclass(frozen=True, slots=True)
class _Edit:
    """What one VCF allele does to the reference, trimmed left first: it replaces the bases
    of deleted from start, and puts inserted there; an insertion goes in before start."""

    start: int
    deleted: str
    inserted: str

    @property
    def end(self) -> int:
        return self.start + len(self.deleted)


@dataclass(frozen=True, slots=True)
class _Window:
    """A window (§1), 0-based and half-open, with the reference's bases over it."""

    start: int
    end: int
    bases: str


def check_same_haplotypes(
    vcf_a: Path, vcf_b: Path, reference: IndexedFasta, *, ploidy: int = DEFAULT_PLOIDY
) -> HaplotypeComparison:
    """Check that two VCF files give every sample the same haplotypes (doc/specs/respell.md
    §8): over each window the same possible genotypes, and over each pair of windows a block
    of either file joins, the same possible genotypes over both.

    Both files are read with parse_vcf over reference with ploidy.

    Raises AssertionError at the first difference, with the contig, the window or the pair,
    1-based and inclusive, the sample and what each file gives there; and ValueError when the
    two headers name different samples.
    """
    samples = parse_vcf_samples(vcf_a)
    samples_b = parse_vcf_samples(vcf_b)
    if samples_b != samples:
        raise ValueError(
            f"{vcf_a} has the samples {', '.join(samples)}"
            f" and {vcf_b} the samples {', '.join(samples_b)}"
        )
    records_a = _read_records_by_contig(vcf_a, reference, ploidy)
    records_b = _read_records_by_contig(vcf_b, reference, ploidy)
    num_windows = num_linked_pairs = num_over_cap = 0
    for contig in dict.fromkeys([*records_a, *records_b]):
        contig_records = (records_a.get(contig, ()), records_b.get(contig, ()))
        windows = _make_windows(reference, contig, contig_records)
        windows_records = [_find_window_records(windows, records) for records in contig_records]
        for sample_index, sample in enumerate(samples):
            pairs = sorted(
                {
                    pair
                    for window_records in windows_records
                    for pair in _find_linked_pairs(window_records, sample_index)
                }
            )
            comparisons: list[tuple[tuple[int, ...], bool]] = [
                ((index,), False) for index in range(len(windows))
            ]
            comparisons += [(pair, True) for pair in pairs]
            for window_indices, is_pair in comparisons:
                results = [
                    _find_sample_result(
                        [windows[index] for index in window_indices],
                        [window_records[index] for index in window_indices],
                        sample_index,
                        ploidy,
                    )
                    for window_records in windows_records
                ]
                # A reason needs no phase assignment, and is compared whatever the cap.
                has_reason = any(isinstance(result, str) for result in results)
                if any(result is _OVER_CAP for result in results) and not has_reason:
                    num_over_cap += 1
                    continue
                result_a, result_b = results
                if result_a != result_b:
                    where = _format_windows(contig, [windows[index] for index in window_indices])
                    raise AssertionError(
                        f"{where}, sample {sample}: {_format_difference(result_a, result_b)}"
                    )
                if is_pair:
                    num_linked_pairs += 1
                else:
                    num_windows += 1
    return HaplotypeComparison(
        num_windows=num_windows, num_linked_pairs=num_linked_pairs, num_over_cap=num_over_cap
    )


def _read_records_by_contig(
    vcf: Path, reference: IndexedFasta, ploidy: int
) -> dict[str, list[VcfRecord]]:
    records: dict[str, list[VcfRecord]] = {}
    for vcf_record in parse_vcf(vcf, reference, ploidy=ploidy):
        records.setdefault(vcf_record.contig, []).append(vcf_record)
    return records


# --- Windows and pairs ---


def _has_alt(vcf_record: VcfRecord) -> bool:
    return len(vcf_record.alleles) > 1


def _find_window_span(vcf_record: VcfRecord) -> tuple[int, int]:
    """The stretch a record with an ALT adds to the windows: its span, or its REF when its
    only ALT is `*`."""
    if all(alt == _SPANNING_DELETION for alt in vcf_record.alleles[1:]):
        return vcf_record.start, vcf_record.end
    return vcf_record.span


def _make_windows(
    reference: IndexedFasta, contig: str, contig_records: Sequence[Sequence[VcfRecord]]
) -> list[_Window]:
    """The windows of a contig: the spans of the records with an ALT of both files, sorted,
    joined while they overlap or lie next to each other."""
    spans = sorted(
        _find_window_span(vcf_record)
        for records in contig_records
        for vcf_record in records
        if _has_alt(vcf_record)
    )
    joined: list[list[int]] = []
    for start, end in spans:
        if joined and start <= joined[-1][1]:
            joined[-1][1] = max(joined[-1][1], end)
        else:
            joined.append([start, end])
    return [_Window(start, end, reference.get_bases(contig, start, end)) for start, end in joined]


def _find_window_records(
    windows: Sequence[_Window], records: Sequence[VcfRecord]
) -> list[list[VcfRecord]]:
    """The records of one file in each window: those with an ALT whose span it holds, and
    those with no ALT whose REF overlaps it, where only their `.` matters."""
    starts = [window.start for window in windows]
    window_records: list[list[VcfRecord]] = [[] for _ in windows]
    for vcf_record in records:
        if _has_alt(vcf_record):
            start, _ = _find_window_span(vcf_record)
            window_records[bisect.bisect_right(starts, start) - 1].append(vcf_record)
            continue
        for index, window in enumerate(windows):
            if window.start < vcf_record.end and vcf_record.start < window.end:
                window_records[index].append(vcf_record)
    return window_records


def _find_unit_key(vcf_record: VcfRecord, genotype: VcfGenotype) -> _UnitKey:
    if genotype.phased:
        return "phase set", genotype.phase_set
    return "line", vcf_record.line_number


def _is_heterozygous(genotype: VcfGenotype) -> bool:
    """Two different alleles, `.` counting as one (doc/specs/reconciliation.md §4.1)."""
    return len(set(genotype.alleles)) > 1


def _find_blocks(windows_records: Sequence[Sequence[VcfRecord]], sample: int) -> list[_UnitKey]:
    """The blocks of a sample in some windows: its phase units with a heterozygous VCF
    genotype there, in the order first met."""
    blocks: dict[_UnitKey, None] = {}
    for records in windows_records:
        for vcf_record in records:
            genotype = vcf_record.genotypes[sample]
            if isinstance(genotype, VcfGenotype) and _is_heterozygous(genotype):
                blocks[_find_unit_key(vcf_record, genotype)] = None
    return list(blocks)


def _find_linked_pairs(
    window_records: Sequence[Sequence[VcfRecord]], sample: int
) -> Iterator[tuple[int, int]]:
    """The pairs of windows a block of one file joins for a sample: of the windows where it
    has a heterozygous VCF genotype, each with the next. A window where it is homozygous
    joins nothing (§8)."""
    windows_of_block: dict[_UnitKey, list[int]] = {}
    for index, records in enumerate(window_records):
        for block in _find_blocks([records], sample):
            windows_of_block.setdefault(block, []).append(index)
    for indices in windows_of_block.values():
        yield from itertools.pairwise(indices)


# --- One sample over one window or a pair ---


def _find_sample_result(
    windows: Sequence[_Window],
    windows_records: Sequence[Sequence[VcfRecord]],
    sample: int,
    ploidy: int,
) -> _SampleResult:
    """A sample's possible genotypes over some windows together: for every phase assignment
    of its blocks there, the multiset of its haplotypes' values over the windows; or the
    reason it is invalid there; or _OVER_CAP."""
    genotypes = [
        vcf_record.genotypes[sample] for records in windows_records for vcf_record in records
    ]
    if any(genotype is MALFORMED_GT for genotype in genotypes):
        return _MALFORMED_GT_REASON
    if any(genotype is OTHER_PLOIDY_GT for genotype in genotypes):
        return _PLOIDY_MISMATCH_REASON
    blocks = _find_blocks(windows_records, sample)
    # The first block keeps the order written, since a genotype is a multiset.
    orders = list(itertools.permutations(range(ploidy)))
    if len(orders) ** max(len(blocks) - 1, 0) > ORACLE_CAP:
        return _OVER_CAP
    possible_genotypes: set[_Genotype] = set()
    for block_orders in itertools.product(orders, repeat=max(len(blocks) - 1, 0)):
        order_of_unit = dict(zip(blocks[1:], block_orders, strict=True))
        haplotypes = []
        for haplotype in range(ploidy):
            values = _find_haplotype_values(
                windows, windows_records, sample, haplotype, order_of_unit
            )
            if values is _IMPOSSIBLE:
                break
            haplotypes.append(values)
        else:
            possible_genotypes.add(tuple(sorted(haplotypes, key=_make_sort_key)))
    return frozenset(possible_genotypes)


def _find_haplotype_values(
    windows: Sequence[_Window],
    windows_records: Sequence[Sequence[VcfRecord]],
    sample: int,
    haplotype: int,
    order_of_unit: Mapping[_UnitKey, tuple[int, ...]],
) -> tuple[_HaplotypeValue, ...] | _IMPOSSIBLE:
    """What one haplotype gives over each of the windows, or _IMPOSSIBLE when it is
    impossible over any."""
    values = []
    for window, records in zip(windows, windows_records, strict=True):
        value = _find_haplotype_value(window, records, sample, haplotype, order_of_unit)
        if value is _IMPOSSIBLE:
            return _IMPOSSIBLE
        values.append(value)
    return tuple(values)


def _find_haplotype_value(
    window: _Window,
    records: Sequence[VcfRecord],
    sample: int,
    haplotype: int,
    order_of_unit: Mapping[_UnitKey, tuple[int, ...]],
) -> _HaplotypeValue | _IMPOSSIBLE:
    """What one haplotype gives over one window (§8): _IMPOSSIBLE when two of its edits
    collide, or when it carries a `*` that neither an edit it carries nor an edit of an ALT
    of a record where it is `.` covers; MISSING when it is `.` at a record; otherwise the
    reference over the window with its edits applied."""
    edits: list[_Edit] = []
    # The edits of every ALT of the records where the haplotype is `.`, which it might carry.
    missing_edits: list[_Edit] = []
    star_refs: list[tuple[int, int]] = []
    is_missing = False
    for vcf_record in records:
        genotype = vcf_record.genotypes[sample]
        assert isinstance(genotype, VcfGenotype), "an invalid GT is the sample's reason"
        order = order_of_unit.get(_find_unit_key(vcf_record, genotype))
        allele = genotype.alleles[haplotype if order is None else order[haplotype]]
        if allele is MISSING:
            is_missing = True
            missing_edits.extend(_make_alt_edits(vcf_record))
        elif vcf_record.alleles[allele] == _SPANNING_DELETION:
            star_refs.append((vcf_record.start, vcf_record.end))
        elif allele != 0:
            edit = _make_edit(vcf_record, vcf_record.alleles[allele])
            if edit is not None:
                edits.append(edit)
    if any(_is_colliding(edit, other) for edit, other in itertools.combinations(edits, 2)):
        return _IMPOSSIBLE
    covering_edits = [*edits, *missing_edits]
    if any(not _is_covered(star_ref, covering_edits) for star_ref in star_refs):
        return _IMPOSSIBLE
    if is_missing:
        return MISSING
    return _apply_edits(window, edits)


def _make_alt_edits(vcf_record: VcfRecord) -> Iterator[_Edit]:
    for alt in vcf_record.alleles[1:]:
        if alt != _SPANNING_DELETION:
            edit = _make_edit(vcf_record, alt)
            if edit is not None:
                yield edit


def _make_edit(vcf_record: VcfRecord, alt: str) -> _Edit | None:
    """The edit of an ALT: the bases it shares with REF removed, first at the left and then
    at the right; None when it spells REF."""
    ref = vcf_record.alleles[0]
    num_left = 0
    while num_left < min(len(ref), len(alt)) and ref[num_left] == alt[num_left]:
        num_left += 1
    deleted, inserted = ref[num_left:], alt[num_left:]
    while deleted and inserted and deleted[-1] == inserted[-1]:
        deleted, inserted = deleted[:-1], inserted[:-1]
    if not deleted and not inserted:
        return None
    return _Edit(vcf_record.start + num_left, deleted, inserted)


def _is_colliding(edit: _Edit, other: _Edit) -> bool:
    """doc/spec.md §5.4: two edits that replace a common position, an insertion between two
    positions the other replaces, or two insertions at one point."""
    if edit.deleted and other.deleted:
        return edit.start < other.end and other.start < edit.end
    if edit.deleted:
        return edit.start < other.start < edit.end
    if other.deleted:
        return other.start < edit.start < other.end
    return edit.start == other.start


def _is_covered(star_ref: tuple[int, int], edits: Sequence[_Edit]) -> bool:
    """Whether an edit, as written, replaces a position of a `*` record's REF."""
    star_start, star_end = star_ref
    return any(edit.deleted and edit.start < star_end and star_start < edit.end for edit in edits)


def _apply_edits(window: _Window, edits: Sequence[_Edit]) -> str:
    """The reference over the window with the edits applied in position order, an insertion
    before a replacement that starts at its point."""
    pieces = []
    position = window.start
    for edit in sorted(edits, key=lambda edit: (edit.start, bool(edit.deleted))):
        assert window.start <= edit.start, "an edit is in its window"
        assert edit.end <= window.end, "an edit is in its window"
        pieces.append(window.bases[position - window.start : edit.start - window.start])
        pieces.append(edit.inserted)
        position = edit.end
    pieces.append(window.bases[position - window.start :])
    return "".join(pieces)


# --- Messages ---


def _make_sort_key(values: Sequence[_HaplotypeValue]) -> tuple[tuple[bool, str], ...]:
    return tuple((value is MISSING, "" if value is MISSING else value) for value in values)


def _format_windows(contig: str, windows: Sequence[_Window]) -> str:
    return f"{contig}:" + " and ".join(f"{window.start + 1}-{window.end}" for window in windows)


def _format_difference(result_a: _SampleResult, result_b: _SampleResult) -> str:
    """What each file gives, its reason or how many possible genotypes, and then the
    genotypes each file gives that the other does not, ten at most for each (§8)."""
    texts = [f"file A gives {_format_count(result_a)} and file B {_format_count(result_b)}"]
    for name, result, other in (("A", result_a, result_b), ("B", result_b, result_a)):
        if isinstance(result, frozenset):
            other_genotypes = other if isinstance(other, frozenset) else frozenset()
            only = sorted(_format_genotype(genotype) for genotype in result - other_genotypes)
            if only:
                texts.append(f"only in {name}: {_format_list(only)}")
    return "; ".join(texts)


def _format_count(result: _SampleResult) -> str:
    if result is _OVER_CAP:
        return f"more than {ORACLE_CAP} phase assignments"
    if isinstance(result, str):
        return result
    if not result:
        return "no genotype"
    return "1 genotype" if len(result) == 1 else f"{len(result)} genotypes"


def _format_list(genotypes: Sequence[str]) -> str:
    """The first ten genotypes, and how many more there are."""
    text = " or ".join(genotypes[:_MAX_LISTED])
    num_more = len(genotypes) - _MAX_LISTED
    return f"{text} and {num_more} more" if num_more > 0 else text


def _format_genotype(genotype: _Genotype) -> str:
    return "{" + ", ".join(_format_haplotype(values) for values in genotype) + "}"


def _format_value(value: _HaplotypeValue) -> str:
    if value is MISSING:
        return "."
    return value or "-"


def _format_haplotype(values: Sequence[_HaplotypeValue]) -> str:
    texts = [_format_value(value) for value in values]
    return texts[0] if len(texts) == 1 else f"({', '.join(texts)})"
