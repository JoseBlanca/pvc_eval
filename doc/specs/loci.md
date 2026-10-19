# The loci

*Status: approved by the owner on 19 October 2026, with kind and stratum in this spec rather than
in reconciliation's. No code yet. It
depends on `doc/spec.md` §2 (vocabulary), §4 (the two contracts, "What stage 1 yields" and "What
stage 3 yields"), §5.2 (spans), §5.3 (loci), §5.6 (kind) and §5.7 (strata), and on
`doc/specs/str_finder.md` for the segments.*

The loci are stage 3 of pvc_eval's pipeline. `build_loci` takes the VCF records of file A and of
file B and the STR regions of the reference, and yields the loci: the stretches of reference that
are compared, each with the records of both files and the STR regions that fall in it. A locus that
is too narrow splits one event between two rows, and each of them comes out *different*. A locus
that is too wide joins unrelated events in one row, which is *same genotype* only when all of them
agree. This spec gives `build_loci` in full: how it merges its inputs,
what it holds while it works, what it does at each boundary, and what it raises. After it come the
three things read from a locus once it is built: whether `--regions` selects it, its kind and its
stratum. No point is open.

`doc/spec.md` §3, §4 and §5.3 already give the rules of the loci, the `Locus` type and the
signature of `build_loci`. This spec decides, beyond them:

- the order in which the two files and the segments are read, and which records with no ALT and
  which STR regions are held, and for how long;
- that a locus holds every STR region that lies inside it, also one that no record touches;
- the `ValueError`s of `build_loci`, for inputs out of order or on an unknown contig;
- how `--regions` and `--strata` BEDs are read, and the text of a stratum from a BED;
- that a locus whose only ALTs are `*` is of kind complex.

Positions in the cases below are 1-based and inclusive, as in a VCF and in `doc/spec.md`: "span
100–103" is the four positions 100 to 103. In the code and in the tests they are 0-based and
half-open, as `VcfRecord.span` is, so span 100–103 is `(99, 103)`.

## Vocabulary

The terms of `doc/spec.md` §2 keep their meaning. These are added:

| term | meaning | not |
|---|---|---|
| **record with an ALT** | a VCF record with at least one ALT, `*` included. Only these start, chain and widen loci (§5.3); a record whose ALT is `.` is a **record with no ALT** | |
| **gap** | the number of reference positions that lie between two stretches of one contig: 0 when they are next to each other, and none when they overlap | "distance" |
| **chains** | a record chains to a locus when its span overlaps the locus or leaves a gap of at most the chain distance | "joins", "links" |
| **held** | kept by `build_loci` after it was read from an input and before it is yielded in a locus or dropped | "buffered", "cached" |

## Building the loci

### The rules

`doc/spec.md` §5.3 gives the rules, and they are not repeated here beyond this summary. A locus
starts from the leftmost record with an ALT that is not yet in a locus. It takes in every record
with an ALT that chains to it, and widens to every STR region that a record with an ALT in it
touches, until neither adds anything. A locus that then lies within the chain distance of the
locus before it on its contig is joined to it (§5.3 step 4). These sections say how each rule
applies at its boundaries.

**Chaining.** A record whose span starts *g* positions after the end of the locus leaves a gap of
*g* − 1. It chains when that gap is at most the chain distance. With the default of 10, a locus
ending at 100 takes a SNP at 111 (gap 10, positions 101 to 110) and leaves a SNP at 112 (gap 11).
With a chain distance of 0, SNPs at 100 and 101 are one locus (gap 0) and SNPs at 100 and 102 are
two. A record's span can start before the locus ends; it then overlaps the locus and chains whatever
the chain distance. The chain distance is measured from the locus, which by then may hold STR
regions, and not from the nearest record in it. In the code, with the locus from `start` to `end`
and a span from `span_start` to `span_end`, all 0-based and half-open, the record chains when
`span_start <= end + chain_distance`; records come by span start, so their `span_start` is never
before the locus's first record.

**Touching an STR region.** A record touches an STR region when its span overlaps the region or
ends on the position just before the region's first, or starts on the position just after its
last: in 0-based half-open terms, `span_start <= region.end and region.start <= span_end`. Only a
record with an ALT touches. The locus itself does not: an STR region next to the
locus but touched by no record of it stays out. With an AT repeat as an STR region at 11–20:

| records with an ALT | locus | STR regions |
|---|---|---|
| SNP at 10 | 10–20 | 11–20, touched at its start |
| SNP at 21 | 11–21 | 11–20, touched at its end |
| SNP at 9 | 9 | none: position 10 lies between |
| SNP at 10, and a second STR region at 21–30 | 10–20 | 11–20 only: no record touches 21–30, which lies next to the first region |

**The STR regions a locus holds** are those that overlap it, and each of them lies wholly inside it.
Most were touched by a record. The others lie in a gap between two records that chained without
touching them. SNPs at 10 and 18 chain (gap 7). A homopolymer at 12–16 touches neither, as one
position lies between it and each SNP. The locus is 10–18, and its `str_regions` holds 12–16, so
its kind is STR (§5.6). The locus allele covers the repeat whichever record touched it, so
`str_regions` describes what the row compares. An STR region that overlaps a locus always lies
wholly inside it: one that reached past the locus's edge would overlap the record or the STR region
at that edge, and STR regions do not overlap.

**A locus that reaches back to the one before it.** A locus can widen to an STR region that
starts before its first record. When that brings it to within the chain distance of the locus
before it, the two are one locus (§5.3 step 4, decided by the owner on 19 October 2026, over
keeping them apart). The distance is measured between the two loci, whose edges can be STR regions,
and not to their records. With a SNP at 100 and a satellite at 105–3000 that it does not touch:

| other records with an ALT | loci | why |
|---|---|---|
| SNP at 2990 | 100–3000 | the SNP at 2990 pulls in the satellite, which leaves a gap of 4 after the locus at 100 |
| SNP at 3001 | 100–3001 | the SNP touches the satellite at its end, and pulls it in the same way |
| SNP at 3005 | 100 and 3005 | no record touches the satellite, which lies inside neither locus and is in neither |
| SNP at 2990, satellite at 112–3000 instead | 100 and 112–3000 | the satellite leaves a gap of 11, more than the chain distance |
| SNP at 2990, satellite at 205–3000 instead, and a satellite at 101–200 | 100–3000 | the SNP at 100 touches 101–200, so its locus is 100–200, and 205–3000 leaves a gap of 4 after it, although the SNP is 104 positions away |
| SNP at 2990, and a homopolymer at 102–103 besides the satellite | 100–3000 | as in the first row; the homopolymer touches no record and lies inside the joined locus. The reach (below) must be computed after the satellite is read |

Two loci therefore never lie within the chain distance of each other. They never overlap either:
an STR region that overlapped a locus would lie inside it, as above.

**Records with no ALT.** A record with no ALT is in every locus its REF overlaps, and in no other
(`doc/spec.md` §5.2). It does not chain, touch or widen, and its REF can reach past the locus.

| records | loci | the record with no ALT, `100 CAGTCAGTCAGTC .` (100–112) |
|---|---|---|
| SNPs at 99 and 113 | 99 and 113 | in neither: it overlaps neither locus |
| SNPs at 100 and 112 | 100 and 112 (gap 11) | in both |
| SNPs at 100 and 110 | 100–110 (gap 9) | in it, reaching 2 positions past its end |
| none with an ALT | none | in no locus, and nothing is yielded |

`bcftools view -a` writes such records: keeping a subset of samples and trimming the ALTs none of
them carries turned `10 CAGT C` into `10 CAGT .` (bcftools 1.24, run on 19 October 2026). A record
whose only ALT is `*` is a record with an ALT, and starts a locus like any other.

**Contigs.** A locus never spans two contigs, whatever the chain distance. SNPs on the last
position of chr1 and on the first of chr2 are two loci. A locus can start on a contig's first
position or end on its last; spans already stop at the contig's ends (§5.2), so nothing here needs
the contig's length. A contig where only one file has records gives loci whose other tuple of
records is empty. A contig where neither file has a record with an ALT gives no locus, and its
segments are not needed.

**Empty inputs.** Two empty streams of records give no loci. One empty stream gives the other file's
loci, each with an empty tuple for the empty file. In both cases `build_loci` reads no segment it
does not need (see "What `build_loci` holds").

### The order of the loci

The two streams of records come each in the order `parse_vcf` yields: by contig in the `.fai`'s
order, then by span start. `build_loci` reads them as one stream, ordered by the rank of the contig
in `contigs`, then by span start, and file A before file B where both are equal. Which of two equal
records is read first changes no locus, since the rules make a locus a set. Inside a locus each
file's records stay in their stream's order, records with no ALT among them (§4).

The loci come contig by contig in the order of `contigs`, and within a contig by start, as §4
requires. On a contig they do not overlap. The segments come by contig in the `.fai`'s order and
then by start. `build_loci` needs `contigs` to compare positions across contigs. When file A's next
record is on chr2 and file B's on chr3, it decides which comes first. When the next STR region is
on chr2 and the records are on chr3, it decides whether that region is behind the records and
dropped, or ahead and kept.

### What `build_loci` holds

`doc/spec.md` §4 says nothing holds more than the current locus and a short look-ahead. In what
follows, positions are 0-based and half-open, and every position is compared with its contig: as
the pair of the contig's rank in `contigs` and the position, so that any position of chr2 comes
after every position of chr1. Without the contig, an STR region read ahead from the next contig
would be dropped or never read. **The next start** is the contig and `span_start` of the next
record not yet handled, of either file and of either kind, whatever its contig. When no record is
left, there is no next start, and every condition that uses it is false. Every record handled
later starts there or after, since both files come in that order.

`build_loci` reads the records one at a time, in the order of the section above. Before it handles
a record, of either kind, it reads the segments up to the first one that comes after the end of
its span: until `(segment rank, segment.start) > (record rank, span_end)`, the rank of each being
that of its contig. Every STR region that the record's span could touch has then been read. After
the last record it reads no more segments. GENERIC segments, the generic regions of `segment_contig`, whose
`region_type` is `RegionType.GENERIC`, are read and dropped.

A locus is open from the record with an ALT that starts it until it is yielded. Its **reach** is
the furthest position a later record can start at and still join it: `end + chain_distance`, or
the `end` of a held STR region on its contig with `region.start <= end + chain_distance`, whichever
is further. While the locus is open, each record read on its contig with `span_start <= reach` is
taken: a record with an ALT joins the locus and widens it, with the STR regions it touches, and a
record with no ALT is held. A record with an ALT that starts past `end + chain_distance` and no
further than the end of such a region always touches the region, since the region starts before
the record and ends at or after its start; by §5.3 that region joins the record's locus, which is
then within the chain distance of this one, so step 4 makes the two one locus. For each next
record, `build_loci` first reads the segments up to it, as above, then computes the reach from the
locus and every STR region held at that moment, and only then decides whether the record is taken
or the locus is yielded. The reach must not be kept from the last time the locus widened: an STR
region that starts within the chain distance is often read only when the next record's segments
are, as in the last row of the table above. The locus is yielded when the next record is on a
later contig, starts past the reach, or does not exist. No later record can then join it. It is yielded with its records
with an ALT, and with the held STR regions and records with no ALT that overlap it.

What is held besides the open locus, and when each is dropped:

- **An STR region** is held while a locus is open on its contig and `region.end > start` of the
  locus, since it may lie in a gap of it, or while its end is at or after the next start, since a
  later record could touch it. A region on a later contig than the open locus is held by the second
  condition. Otherwise it is dropped. A region can be held far ahead of the open locus, as the
  satellite of the case above is held while the locus at 100 is built.
- **A record with no ALT** is held while a locus is open on its contig and `record.end > start` of
  the locus, while its end is after the next start, or while it overlaps a held STR region, which could still
  widen a later locus over it. Otherwise it is dropped.

So the records with no ALT of an all-sites VCF are dropped as soon as the reading passes them,
except over an STR region that may still join a locus. Over a satellite of 5 kb with no record with
an ALT in it, that holds about 5,000 records until the reading passes the satellite's end. When
such a satellite starts within the chain distance after the open locus, the locus is held too,
until the reading passes the satellite's end: the third row of the table of a locus that reaches
back holds the locus at 100 until the SNP at 3005 is read. The open locus itself grows with the
data, as §5.3's simulation shows, and with the satellites it holds; nothing here bounds it, which is what the phasing cap (§5.4) is for.

Because the segments are read only as far as the records need, a segment iterator that runs the
STR finder contig by contig, lazily, is asked for no contig past the last record. `compare` builds
such an iterator over the contigs that hold records (see "Not in this spec").

### What `build_loci` checks

The inputs come from `parse_vcf` and from `segment_contig` or `read_str_regions`, which check the
files and give the user the errors of `doc/spec.md` §3. What `build_loci` finds wrong is a fault of
the program that called it, such as a filter between the stages that reorders records. So it raises
`ValueError`, naming the file (A or B) or the segments and the position, when the stream reaches it:

- a record or a segment on a contig that `contigs` lacks;
- a record of one file that comes before the one read before it, by contig rank or by span start;
- a segment, GENERIC or not, that comes before the one read before it, by contig rank or by start.
  GENERIC segments are checked too: the reading of segments stops at the first one past a record,
  of any type, so a GENERIC segment out of order would stop it before an STR region the record
  touches, and the locus would come out too narrow with no error;
- an STR region that overlaps the STR region read before it. GENERIC segments are not checked
  for overlap, since none of them becomes part of a locus;
- a contig named twice in `contigs`, when the stream starts.

A chain distance that is not a whole number of 0 or more raises `ValueError` when the stream starts,
that is, at the first `next()`, as `parse_vcf` does for its span limit. `True` and `False` are not
whole numbers here.

## The Python interface

`pvc_eval.locus` holds `Locus`, and `pvc_eval.segments` holds `GenomeSegment` and `RegionType`,
both exactly as `doc/spec.md` §4 gives them, fields included.

`build_loci`, in `pvc_eval.loci`, yields the loci of the records of both files. Its signature is the
one of `doc/spec.md` §4. `records_a` and `records_b` are in the order of `parse_vcf`. `segments` is
in the order of `read_str_regions`, with or without GENERIC segments. `contigs` holds the
reference's contig names in the `.fai`'s order. Positions are 0-based and half-open. It raises
`ValueError` as the section above says.

```python
def build_loci(records_a: Iterator[VcfRecord], records_b: Iterator[VcfRecord],
               segments: Iterator[GenomeSegment], contigs: Sequence[str], *,
               chain_distance: int = 10) -> Iterator[Locus]: ...

DEFAULT_CHAIN_DISTANCE = 10
```

## What is read from a locus

These do not change which loci there are. `compare` applies them to each locus that `build_loci`
yields.

### Regions

`--regions` keeps a locus when it shares at least one position with a region, and skips it
otherwise (`doc/spec.md` §3). A locus is the stretch from `start` to `end`; a record with no ALT
whose REF reaches past it does not count. A locus at 100–120 is kept by a region at 110–200, and
compared whole, its records at 100–109 included. It is skipped by a region at 121–200, which only
lies next to it.

`read_regions` reads the BED whole, since BEDs such as GIAB's come sorted by contig name rather
than in the `.fai`'s order. GIAB's benchmark regions for HG002 on GRCh38 hold about 300,000 lines.
It skips blank lines and lines that start with `#`, `track` or `browser`, as `read_str_regions`
does, and uses the first three columns of the others. Regions that overlap or lie next to each other
are merged. It raises `PvcEvalInputError`, naming the file and line and how to put it right, for a
line with fewer than three columns, a start or end that is not a whole number, or an end not after
its start, and for a line on a contig the reference lacks. The last stops the run at the first such
line even when the reference was cut to one chromosome and the BED is GIAB's, which covers chr1–22,
X and Y: the message names the contig and says to keep only the reference's contigs in the BED, and
how (decided by the owner on 19 October 2026, over skipping such lines with a warning).

`select_loci` yields, in their order, the loci that share a position with the regions, and every
locus when `regions` is None.

```python
class Regions:
    def overlaps(self, contig: str, start: int, end: int) -> bool: ...

def read_regions(path: Path, contigs: Sequence[str]) -> Regions: ...
def select_loci(loci: Iterator[Locus], regions: Regions | None) -> Iterator[Locus]: ...
```

### Kind

`locus_kind`, in `pvc_eval.kind`, gives the kind of `doc/spec.md` §5.6. A locus holding an STR
region is STR. Otherwise the kind is read from every ALT of every record with an ALT in the locus.
Records with no ALT, `*` alleles and an ALT equal to its REF, which makes no edit, are left out. Each ALT is trimmed against its record's REF,
first at the left and then at the right, as an edit is (§5.4). It is then a substitution of one
base, an insertion, a deletion, or anything else. When every ALT is a substitution of one base the
locus is SNP, when every ALT is an insertion or a deletion it is indel, and otherwise it is complex.
These cases decide it:

| records in the locus | kind |
|---|---|
| `10 CAG C,CAT` | complex: the first ALT deletes AG and the second turns G into T |
| `10 AC GT` | complex: two bases replaced by two |
| `10 A G` and `14 C CT` | complex: records of two kinds |
| `10 A *` alone, with no deletion in the locus | complex: no ALT is left to be a SNP or an indel |

The last case is decided here: §5.6's "every ALT" would make a locus with no ALT left a SNP and an
indel at once. Such a locus is invalid for every sample that carries the `*` (§5.4) whichever kind
it gets.

```python
class Kind(Enum):
    SNP = "SNP"
    INDEL = "indel"
    STR = "STR"
    COMPLEX = "complex"

def locus_kind(locus: Locus) -> Kind: ...
```

### Stratum

With `--strata pvc`, `pvc_stratum`, in `pvc_eval.strata`, gives the one class of `doc/spec.md` §5.7
from the region types of `locus.str_regions`: `satellite`, `str-bundle`, `str` or `generic`. An STR
region with no region type, from a `--str-regions` BED without a fourth column, raises
`PvcEvalInputError`. Its message names the region and says that `--strata pvc` needs the types that
`pvc_eval str-regions` writes. A locus holding no STR region is `generic` and raises nothing, so
the error comes at the first locus that holds an untyped region, after the rows of earlier loci are
written (decided by the owner on 19 October 2026, over reading the BED once before the comparison
to find an untyped line before any row is written).

With a BED, `read_strata` reads it whole, as `read_regions` does, with the same refusals, plus one
for a line with no fourth column. `Strata.stratum` gives the classes of the intervals that share at
least one position with the locus, each once, in alphabetical order and joined with `,`, so the same
classes always give the same text. It gives `.` when the locus shares no position with any of them.
A locus at 100–120 under `homopolymer` at 90–101 and `segdup` at 115–500 is `homopolymer,segdup`.

```python
def pvc_stratum(locus: Locus) -> str: ...

class Strata:
    def stratum(self, locus: Locus) -> str: ...

def read_strata(path: Path, contigs: Sequence[str]) -> Strata: ...
```

## How it is verified

Each check is made at the function named, with records built by hand as `VcfRecord` values, their
spans given as literals, and segments as lists of `GenomeSegment`. So the tests need no FASTA and do
not wait for the STR finder. Each case's loci are compared with a literal list of `Locus` values.

**At `build_loci`, the cases of this spec:** the chain distance at a gap of 10 and of 11, and at 0
with adjacent SNPs and with a gap of 1; the four rows of the STR table; the STR region in a gap; the
six rows of the table of a locus that reaches back; the four rows of the table of records with no
ALT; a record with no ALT in the gap between two chained SNPs (SNP at 10, `12 A .`, SNP at 15: one
locus 10–15 holding it); a `*`-only record starting a locus; two cases across contigs, 0-based with
a chain distance of 10: a satellite at (110, 3000) with SNPs at (99, 100) and (2989, 2990), one
locus (99, 3000), the satellite starting exactly at the chain distance; a SNP at (9, 10) and one at
(399, 400) on chr1 with a satellite at (0, 500) on chr2 and no segment on chr1, two loci on chr1
and the satellite in neither; a record with no ALT over an STR region that no record touches yet,
kept until a later record pulls the region in; a SNP at (29, 30) on chr1, then on chr2 an STR region at (0, 20) and a SNP at (4, 5), which give a chr2 locus (0, 20) holding the
region; and a SNP at (9, 10) with a region at (4999, 5010) on chr1, then on chr2 a region at
(89, 99) and a SNP at (99, 100), which give a chr2 locus (89, 100) holding it; SNPs at the last
position of chr1 and the first of chr2 with a chain distance of
1,000; a contig with records in one file only; both streams empty, one empty, and one holding only
records with no ALT; segments with and without GENERIC segments between the STR regions, giving the
same loci; each `ValueError` of "What `build_loci` checks".

**At `build_loci`, the cases `doc/spec.md` §7.2 lists for loci:** chaining across both files; the
traded events of §5.3, one locus at chain distance 1 and two at 0; adjacent SNPs, chained at 0; a
record with no ALT, which starts no locus; an insertion next to a repeat pulling it in; two complex
records at the two ends of a repeat. The traded events and the second and third cases of §5.2 are
also run from VCF text and a FASTA in `tmp_path` through `parse_vcf`, so that the spans `parse_vcf`
computes are the ones the loci rely on. In those two cases of §5.2 the records' REFs already touch,
so they share a locus even without the span. One more case is joined by the span alone: reference
20–29 `GCAGAGAGTC`, file A `21 CAG C`, file B `28 T C`, chain distance 0. A's deletion of 22–23
slides right to 26–27, so its span is 21–27 and the two records are one locus, 21–28; by their REFs
alone, 21–23 and 28, they would leave a gap of 4 and make two.

**What it reads, at `build_loci`.** The record streams and the segments are iterators that count
what is read from them. When the *k*-th locus is yielded, counting from 1, from 100,000 SNPs
of file A 100 positions apart and an empty file B, at most *k* + 1 records have been read from file
A. The segment iterator raises if it is read past the first segment that comes after
the end of the span of the last record handled, by contig rank and then position. Records with no ALT are not kept once passed: with 10,000 and with
100,000 of them one position apart, made by a generator as each is read, followed by one SNP, the peak that tracemalloc traces while
`build_loci` runs differs by less than 10%, and the SNP's locus holds only the one that overlaps
it. STR regions are not kept once passed either: with 10,000 and with 100,000 SNPs 100 positions
apart, each followed by an STR region it does not touch, the peak differs by less than 10%. Each
memory test runs `build_loci` once on a small input before it measures, so that what the first
call allocates is not counted. With no records the segment iterator is never read.

**A property, at `build_loci`, with hypothesis.** A slow version of `doc/spec.md` §5.3, written in
the test from the rules alone, holds every record and STR region in lists and applies the rules
until nothing changes. On generated records with an ALT and without, of both files, on two contigs,
with STR regions, chain distances from 0 to 20, records and STR regions near the end of the first
contig and the start of the second, and STR regions within the chain distance of each other, its
loci equal `build_loci`'s. The same cases
also check that every record with an ALT is in exactly one locus, inside it; that no two loci of a
contig lie within the chain distance of each other; that each record with
no ALT is in exactly the loci its REF overlaps; and that the loci of a contig come by start and do
not overlap.

**With the STR finder.** One test feeds `build_loci` from `segment_contig` over a small FASTA in
`tmp_path` that holds an AT repeat and an insertion of AT before it. It is skipped, with a message,
until `segment_contig` exists on the branch.

**At `select_loci`, `locus_kind`, `pvc_stratum` and `Strata.stratum`:** the cases of their
sections; a region list read from an unsorted BED with overlapping lines; each refusal of
`read_regions` and `read_strata`; the strata cases of `doc/spec.md` §7.2 at `pvc_stratum`; and
one locus of each kind, including a *too complex* one, which is §7.2's.

## Not in this spec

- How `compare` finds the contigs that hold records and chains `segment_contig` over them: the
  spec of `compare` (`doc/spec.md` §3).
- Reading `--str-regions` and the STR finder: `doc/specs/str_finder.md`.
- What a locus gives the comparison, its locus alleles and possible genotypes: reconciliation,
  `doc/spec.md` §5.4.
- How `str_regions`, the kind and the stratum are written in the table: `doc/spec.md` §6.
