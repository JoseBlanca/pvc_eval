"""The VCF files the properties of doc/specs/respell.md §9.2 generate with hypothesis: small
contigs of random bases and planted repeats, and a few records on them dense enough to
interact. make_generated_vcf is the strategy, GeneratedVcf what it gives, PlantedRepeat
one repeat of it, and write_generated_vcf writes one to a directory; and
check_respelled_haplotypes, the oracle's check of a re-spelling as its tests make it.

hypothesis draws the ends of a range of whole numbers more often than its middle, so a rare
kind is drawn at the top of its range, where it is also not what shrinking leads to. Over
2,000 files on 20 October 2026 the shares were: two contigs in 18% of files (§9.2 asks one
in five), ploidy 1 in 21% (one in five), a contig with no N in 27% (one in four); no ALT in
4.1% of records (one in twenty), a `*` ALT in 4.8% (one in twenty); `.` in 8.2% of alleles
(one in ten); a malformed GT in 4.4% of GTs and one of the other ploidy in 1.6% (one in
fifty each).

Positions are 0-based and half-open, but for POS in the lines, which is 1-based.
"""

from collections.abc import Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from hypothesis import strategies as st
from tests.conftest import CaseFiles, write_case_vcf
from tests.haplotype_oracle import ORACLE_CAP, HaplotypeComparison, check_same_haplotypes

from pvc_eval.fasta import IndexedFasta

_BASES = "ACGT"
_CONTIG_NAMES = ("chr1", "chr2")
# §9.2: a contig's bases, before any N, and the pieces they are made of.
_MIN_CONTIG_BASES = 12
_MAX_CONTIG_BASES = 60
_PADDING_BASES = 4  # the most N on each side of a contig that is padded


class _RecordKind(Enum):
    """What a record with one ALT is drawn as."""

    SNP = "SNP"
    REPLACEMENT = "replacement"
    REPEAT_INDEL = "indel of a planted repeat"
    RANDOM_INDEL = "random indel"
    DELETION_OVER = "deletion over another record"


@dataclass(frozen=True, slots=True)
class PlantedRepeat:
    """A repeat the generator planted: copies of unit from start to end, 0-based, half-open."""

    contig: str
    start: int
    end: int
    unit: str


@dataclass(frozen=True, slots=True, kw_only=True)
class GeneratedVcf:
    """A generated file.

    contigs: each contig's name and bases, in the order of the file.
    lines: the records in the notation of the cases with CHROM first, `chr1 POS REF ALT`
        and one sample field per sample, `GT` or `GT:PS`, sorted by contig and POS.
    num_samples: the samples, S0, S1...
    ploidy: what the file is read with; a GT of the other ploidy is drawn too.
    repeats: the repeats planted, for a --str-regions BED.
    seed: a seed for the re-spellings that draw.
    """

    contigs: tuple[tuple[str, str], ...]
    lines: tuple[str, ...]
    num_samples: int
    ploidy: int
    repeats: tuple[PlantedRepeat, ...]
    seed: int


@dataclass(frozen=True, slots=True)
class _Contig:
    """A contig's bases, the stretch of them that are not N, and its planted repeats."""

    name: str
    bases: str
    core_start: int
    core_end: int
    repeats: tuple[PlantedRepeat, ...]


@dataclass(frozen=True, slots=True)
class _Record:
    """A record before its genotypes: start, 0-based, REF and the ALTs."""

    start: int
    ref: str
    alts: tuple[str, ...]


def check_respelled_haplotypes(
    vcf: Path, output: Path, reference: IndexedFasta, *, ploidy: int
) -> HaplotypeComparison:
    """Check with the oracle that output, a re-spelling of vcf, gives every sample the same
    haplotypes, and that no window or pair was over the oracle's cap (doc/specs/respell.md
    §9.2), since one left uncompared could hide a lost phase set."""
    comparison = check_same_haplotypes(vcf, output, reference, ploidy=ploidy)
    num_over_cap = comparison.num_over_cap
    if num_over_cap:
        raise AssertionError(
            f"{num_over_cap} window or pair over the oracle's cap of {ORACLE_CAP} phase"
            f" assignments, not compared, between {vcf} and {output}"
        )
    return comparison


def write_generated_vcf(
    generated: GeneratedVcf, directory: Path, stack: ExitStack
) -> tuple[IndexedFasta, Path]:
    """The reference, opened in stack, and the VCF file of a generated file, in directory."""
    reference = CaseFiles(directory, stack).write_contigs(dict(generated.contigs))
    vcf = write_case_vcf(
        directory / "a.vcf",
        generated.lines,
        num_samples=generated.num_samples,
        has_chrom=True,
        fills_ps=False,
    )
    return reference, vcf


@st.composite
def make_generated_vcf(draw: st.DrawFn) -> GeneratedVcf:
    """A file of §9.2: one contig, or two in one file of five; 1 to 3 samples, of ploidy 2
    in four files of five and 1 otherwise; and on each contig 1 to 6 records."""
    num_contigs = 2 if draw(st.integers(0, 4)) == 4 else 1
    num_samples = draw(st.integers(1, 3))
    ploidy = 2 if draw(st.integers(0, 4)) < 4 else 1
    contigs = [draw(_make_contig(name)) for name in _CONTIG_NAMES[:num_contigs]]
    lines = []
    for contig in contigs:
        for record in draw(_make_records(contig)):
            num_alts = 0 if record.alts == (".",) else len(record.alts)
            sample_fields = [draw(_make_sample_field(num_alts, ploidy)) for _ in range(num_samples)]
            fields = [contig.name, str(record.start + 1), record.ref, ",".join(record.alts)]
            lines.append(" ".join([*fields, *sample_fields]))
    return GeneratedVcf(
        contigs=tuple((contig.name, contig.bases) for contig in contigs),
        lines=tuple(lines),
        num_samples=num_samples,
        ploidy=ploidy,
        repeats=tuple(repeat for contig in contigs for repeat in contig.repeats),
        seed=draw(st.integers(0, 2**32 - 1)),
    )


# --- The contigs ---


@st.composite
def _make_contig(draw: st.DrawFn, name: str) -> _Contig:
    """12 to 60 bases made of 2 to 5 pieces, each 3 to 8 random bases or a repeat of a unit
    of 1 to 3 bases in 3 to 6 copies; padded with N on both sides but in one contig of four,
    which is those bases alone, so that records reach its ends. Pieces past 60 bases are
    left out, and random bases are added to reach 12."""
    padding = 0 if draw(st.integers(0, 3)) == 3 else draw(st.integers(1, _PADDING_BASES))
    core = ""
    repeats = []
    for _ in range(draw(st.integers(2, 5))):
        if draw(st.booleans()):
            unit = draw(st.text(_BASES, min_size=1, max_size=3))
            piece = unit * draw(st.integers(3, 6))
        else:
            unit = ""
            piece = draw(st.text(_BASES, min_size=3, max_size=8))
        if len(core) + len(piece) > _MAX_CONTIG_BASES:
            break
        if unit:
            start = padding + len(core)
            repeats.append(PlantedRepeat(name, start, start + len(piece), unit))
        core += piece
    if len(core) < _MIN_CONTIG_BASES:
        core += draw(st.text(_BASES, min_size=_MIN_CONTIG_BASES - len(core), max_size=8))
    return _Contig(
        name=name,
        bases="N" * padding + core + "N" * padding,
        core_start=padding,
        core_end=padding + len(core),
        repeats=tuple(repeats),
    )


# --- The records ---


@st.composite
def _make_records(draw: st.DrawFn, contig: _Contig) -> list[_Record]:
    """1 to 6 records, sorted by POS: one in twenty with no ALT, one in twenty with a `*` ALT,
    and the others a SNP, a replacement, the insertion or the deletion of whole units of a
    planted repeat at its edge, a random indel, or a deletion over a record drawn before it,
    each with a second ALT one time in two."""
    records: list[_Record] = []
    for _ in range(draw(st.integers(1, 6))):
        kind = draw(st.integers(0, 19))
        if kind == 19:
            records.append(draw(_make_no_alt_record(contig)))
            continue
        record = draw(_make_one_alt_record(contig, records))
        if kind == 18:
            records.append(_Record(record.start, record.ref, (*record.alts, "*")))
            continue
        if draw(st.booleans()):
            second_alt = draw(
                _make_alt(record.ref).filter(lambda alt, record=record: alt not in record.alts)
            )
            record = _Record(record.start, record.ref, (*record.alts, second_alt))
        records.append(record)
    return sorted(records, key=lambda record: record.start)


@st.composite
def _make_no_alt_record(draw: st.DrawFn, contig: _Contig) -> _Record:
    """A record with no ALT, `.`, over 1 to 3 bases."""
    start = draw(st.integers(contig.core_start, contig.core_end - 1))
    length = draw(st.integers(1, min(3, contig.core_end - start)))
    return _Record(start, contig.bases[start : start + length], (".",))


@st.composite
def _make_one_alt_record(draw: st.DrawFn, contig: _Contig, records: Sequence[_Record]) -> _Record:
    """A record of one ALT, of a kind drawn evenly. The indel of a planted repeat on a contig
    with none, and the deletion over another record where no record is drawn before it or
    the one drawn is at the first base of the contig's bases, are a SNP instead."""
    kind = draw(st.sampled_from(_RecordKind))
    if kind is _RecordKind.REPEAT_INDEL and contig.repeats:
        return draw(_make_repeat_indel(contig))
    if kind is _RecordKind.DELETION_OVER and records:
        record = draw(_make_deletion_over(contig, draw(st.sampled_from(records))))
        if record is not None:
            return record
    if kind is _RecordKind.REPLACEMENT:
        return draw(_make_replacement(contig))
    if kind is _RecordKind.RANDOM_INDEL:
        return draw(_make_random_indel(contig))
    return draw(_make_snp(contig))


@st.composite
def _make_snp(draw: st.DrawFn, contig: _Contig) -> _Record:
    start = draw(st.integers(contig.core_start, contig.core_end - 1))
    ref = contig.bases[start]
    return _Record(start, ref, (draw(st.sampled_from([base for base in _BASES if base != ref])),))


@st.composite
def _make_replacement(draw: st.DrawFn, contig: _Contig) -> _Record:
    """REF of 2 to 4 bases, and an ALT of 1 to 4 random bases, so complex as often as not."""
    start = draw(st.integers(contig.core_start, contig.core_end - 2))
    length = draw(st.integers(2, min(4, contig.core_end - start)))
    ref = contig.bases[start : start + length]
    alt = draw(st.text(_BASES, min_size=1, max_size=4).filter(lambda alt: alt != ref))
    return _Record(start, ref, (alt,))


@st.composite
def _make_random_indel(draw: st.DrawFn, contig: _Contig) -> _Record:
    """The insertion of 1 to 3 random bases after an anchor, or the deletion of 1 to 3 bases
    after it."""
    start = draw(st.integers(contig.core_start, contig.core_end - 1))
    anchor = contig.bases[start]
    num_after = contig.core_end - start - 1
    if num_after == 0 or draw(st.booleans()):
        return _Record(start, anchor, (anchor + draw(st.text(_BASES, min_size=1, max_size=3)),))
    length = draw(st.integers(1, min(3, num_after)))
    return _Record(start, contig.bases[start : start + 1 + length], (anchor,))


@st.composite
def _make_repeat_indel(draw: st.DrawFn, contig: _Contig) -> _Record:
    """The insertion or the deletion of whole units of a planted repeat, at its left edge,
    anchored on the base before it, or at its right edge, anchored on its last base or on the
    base before the units deleted."""
    repeat = draw(st.sampled_from(contig.repeats))
    unit = repeat.unit
    num_copies = (repeat.end - repeat.start) // len(unit)
    num_units = draw(st.integers(1, num_copies))
    has_left_anchor = repeat.start > contig.core_start
    if has_left_anchor and draw(st.booleans()):
        anchor = repeat.start - 1
        anchor_base = contig.bases[anchor]
        if draw(st.booleans()):
            return _Record(anchor, anchor_base, (anchor_base + unit * num_units,))
        deleted = unit * num_units
        return _Record(anchor, anchor_base + deleted, (anchor_base,))
    last = repeat.end - 1
    if draw(st.booleans()):
        last_base = contig.bases[last]
        return _Record(last, last_base, (last_base + unit * num_units,))
    if num_units == num_copies and not has_left_anchor:
        # Deleting every copy needs the base before the repeat as the anchor.
        num_units -= 1
    anchor = repeat.end - num_units * len(unit) - 1
    return _Record(anchor, contig.bases[anchor : repeat.end], (contig.bases[anchor],))


@st.composite
def _make_deletion_over(draw: st.DrawFn, contig: _Contig, record: _Record) -> _Record | None:
    """A deletion that deletes the first base of record's REF, anchored 1 to 3 bases before
    it, so that the `*` re-spelling has records under deletions; None when the record is at
    the first base of the contig's bases."""
    if record.start == contig.core_start:
        return None
    anchor = draw(st.integers(max(contig.core_start, record.start - 3), record.start - 1))
    end = draw(st.integers(record.start + 1, min(contig.core_end, record.start + 3)))
    return _Record(anchor, contig.bases[anchor:end], (contig.bases[anchor],))


def _make_alt(ref: str) -> st.SearchStrategy[str]:
    """A second ALT for a REF: a SNP of its first base, its first base with 0 to 3 random
    bases after it, or 1 to 4 random bases; never REF itself."""
    snp = st.sampled_from([base for base in _BASES if base != ref[0]]).map(
        lambda first: first + ref[1:]
    )
    anchored = st.text(_BASES, max_size=3).map(lambda tail: ref[0] + tail)
    random_bases = st.text(_BASES, min_size=1, max_size=4)
    return st.one_of(snp, anchored, random_bases).filter(lambda alt: alt != ref)


# --- The genotypes ---


@st.composite
def _make_sample_field(draw: st.DrawFn, num_alts: int, ploidy: int) -> str:
    """A sample's field: a GT phased in one draw of two, a phased one with a PS of `1`, `2`
    or none; `.` in one allele of ten; malformed, with an allele past the ALTs, in one GT of
    fifty, but at a record with no ALT, whose GTs are of `0` and `.` alone (§9.2); and of the
    other ploidy in one GT of fifty."""
    kind = draw(st.integers(0, 49))
    num_alleles = (3 - ploidy) if kind == 48 else ploidy
    alleles = [draw(_make_allele(num_alts)) for _ in range(num_alleles)]
    if kind == 49 and num_alts > 0:
        alleles[0] = str(num_alts + 1)
    if num_alleles == 1:
        return alleles[0]
    if not draw(st.booleans()):
        return "/".join(alleles)
    phase_set = draw(st.sampled_from(["1", "2", None]))
    gt = "|".join(alleles)
    return gt if phase_set is None else f"{gt}:{phase_set}"


def _make_allele(num_alts: int) -> st.SearchStrategy[str]:
    """`.` one draw in ten, an index of the record's alleles otherwise."""
    return st.integers(0, 9).flatmap(
        lambda draw_number: st.just(".") if draw_number == 9 else st.integers(0, num_alts).map(str)
    )
