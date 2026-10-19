# pvc_eval — specification

*Status: draft, 2026-10-19, for review before any code. pvc_eval is its own project: it imports
nothing from pvc. The one thing it takes from pvc is a port of its STR finder, with that finder's
tests (§5.1).*

## 1. What pvc_eval does

pvc_eval compares two VCF files over the same reference genome and says, for every sample and every
place where either file has a variant, whether the two files give that sample the same genotype.

Comparing VCF files is hard because one genotype can be written in many ways. A deletion inside a
repeat can sit at any of several positions. A deletion and a nearby SNP can be one VCF record or two.
One caller writes a spanning deletion with a `*` allele, another does not. A comparison that matches
VCF records by position and allele calls these different when they are the same; inside tandem
repeats that is most of the disagreement it reports.

pvc_eval avoids this by never comparing VCF records. It groups the VCF records of both files into
**loci** — stretches of reference that contain every way the files could have written the same
events — rebuilds, for each sample, the sequence each haplotype carries over the locus, and compares
those sequences. Two files agree when the sequences agree, however each spelled them.

It treats the two files alike. Neither is the truth: if the user wants one to be, words like
*missed*, *false call*, *recall* and *precision* are the user's reading of pvc_eval's table, not
something pvc_eval computes. It does not filter: a user who wants only PASS records, or a quality
threshold, filters the VCFs first. Its Python API exposes each stage (§4) for anything the command
line does not do.

## 2. Vocabulary

Each term has one meaning in this document and in the code; the right-hand column lists words not to
use for it.

**Inputs**

| term | meaning | not |
|---|---|---|
| **reference** | the genome FASTA, and only that | "reference" for either VCF |
| **file A, file B** | the two VCF files compared | "truth", "query", "baseline", "calls" |
| **VCF record** | one line of a VCF file | "site", "variant", "locus" |
| **VCF allele** | REF or one ALT of a VCF record, as written | |
| **VCF genotype** | a VCF record's GT for one sample, as written | |
| **normalise** | trim the bases a VCF record's alleles share at either end (§5.2). A VCF record is never moved, split or decomposed. | "decompose", "atomise" |
| **VCF record span** | every reference position at which the VCF record could have been written: its REF after trimming, widened for an indel to the whole stretch it can slide along (§5.2) | |
| **sample pair** | a sample of file A matched to a sample of file B: same name unless the user maps them | |
| **ploidy** | the number of haplotypes a sample has at every locus; a global constant, 2 | |

**From the reference**

| term | meaning | not |
|---|---|---|
| **segmented genome** | each contig of the reference cut by the STR finder (§5.1) into consecutive, non-overlapping regions, each a generic region or an STR region | |
| **generic region** | a stretch of the segmented genome that holds no tandem repeat the STR finder admits; pvc's `generic` | "ordinary sequence" |
| **STR region** | a stretch the STR finder (§5.1) reports as a tandem repeat: an **STR segment** (one repeat), an **STR bundle** (two or more repeats within 15 bases of each other) or a **satellite** (a repeat run longer than 100 bases). pvc's terms. | "tract", "STR locus", "repeat locus" |
| **region type** | which of the four a region of the segmented genome is: `generic`, `str` (an STR segment), `str-bundle` or `satellite` | "kind", which is §5.6's |

**Building what is compared**

| term | meaning | not |
|---|---|---|
| **locus** | the unit compared: a stretch of reference holding VCF records of either file, chained by overlap and widened to every STR region a VCF record touches (§5.3) | any other use of "locus" |
| **touches** | overlaps, or lies immediately next to. The span of an insertion of a repeat unit can end on the base just before the repeat; without "next to" it would not pull the repeat in. | |
| **phasing cap** | the most phase assignments (§5.4) pvc_eval tries for one sample at one locus | "record cap" |
| **haplotype** | one of a sample's chromosomes; a sample has *ploidy* haplotypes | "copy", "strand", "phase" |
| **locus allele** | the sequence one haplotype carries over the locus, built from the reference and VCF alleles. The locus's REF is the reference over the locus. | "haplotype" for the sequence |
| **locus genotype** | which locus allele each haplotype carries, unordered | |
| **possible genotypes** | the locus genotypes a sample's VCF genotypes allow: one when they are phased or unambiguous, several when unphased VCF records leave the phase open, none when they contradict each other | "candidate genotypes" |
| **missing** | a VCF genotype says nothing about a haplotype (`.`); the locus allele built from it is missing | |
| **invalid** | a file's VCF records claim something impossible for a sample at a locus (§5.4) | "missing" |

**Annotations on each row of the output**

| term | meaning |
|---|---|
| **outcome** | the result of comparing one sample pair at one locus (§5.5) |
| **kind** | SNP, indel, STR or complex (§5.6) |
| **stratum** | the class `--strata` gives the locus: from a stratification BED (homopolymer, segmental duplication, …), or, with `--strata pvc`, the region type of the STR regions the locus holds (§5.7) |
| **quality** | for one file and sample, the lowest QUAL among the VCF records giving the sample a non-reference VCF allele at the locus |

## 3. Interface

```
pvc_eval compare --reference ref.fa --a a.vcf.gz --b b.vcf.gz --output loci.tsv
                 [--regions regions.bed] [--strata strata.bed|pvc] [--str-regions str.bed]
                 [--sample-map pairs.tsv] [--phasing-cap N] [--span-limit N|none]
                 [--summary summary.tsv] [--all-rows]
                 [--evaluator pvc|vcfeval|happy|naive]
```

| option | what it does |
|---|---|
| `--reference` | FASTA with its `.fai`. Contig names and order come from the `.fai`. |
| `--a`, `--b` | the two VCFs, plain or bgzipped. |
| `--regions` | compare only the loci that overlap these regions, by at least one position; every other locus is skipped and gets no row, nor a count in the summary. A locus that overlaps is compared whole, its VCF records outside the regions included, since loci are built from all VCF records (§5.3). Without it, every locus is compared. |
| `--strata` | BED with a fourth column naming each interval's class; a locus takes the class of the intervals it overlaps (several are joined with `,`; none gives `.`). The word `pvc` instead of a path gives each locus one region type from the STR regions it holds — `satellite`, `str-bundle`, `str` or `generic` — with no BED needed (§5.7); a file named `pvc` is given as `./pvc`. |
| `--str-regions` | BED of STR regions to use instead of running the STR finder. |
| `--sample-map` | two columns: a sample of A, a sample of B. Without it, samples are paired by name. |
| `--phasing-cap` | default 4096 (§5.4, **Open 2**). |
| `--span-limit` | how far, in bases, a VCF record's span may reach beyond its REF before it is cut back to the REF (§5.2): a whole number of 1 or more, default 10000, or `none` for no limit. |
| `--output` | the table of §6, one row per sample pair and locus. |
| `--summary` | counts of outcomes by kind and stratum. |
| `--all-rows` | also write the rows whose outcome is *reference in both* (§5.5). |
| `--evaluator` | what decides each locus's outcome: `pvc`, pvc_eval's own rules (§5.4, §5.5), the default; or `vcfeval`, `happy` or `naive`, another evaluator wrapped behind the same interface (§8). |

**Errors, which stop the run:** a VCF not sorted in the `.fai`'s contig order and by position; a
contig in a VCF but not in the reference; a VCF record whose REF reaches past the end of its
contig; no sample pairs; a sample named in `--sample-map` absent
from its file. Samples present in only one file and not mapped are listed as a warning and not
compared.

**Not errors:** everything about the content of a VCF record — a REF that disagrees with the
reference, a malformed or haploid GT, contradicting VCF records. Those make the affected locus
*invalid* for the affected sample (§5.4) and the run goes on.

The Python API offers each stage of §4 as a function over iterators, so a user can, for example,
filter VCF records between reading and grouping, or read loci and genotypes without evaluating them.
Above the stages, one function does what `compare` does, whichever the evaluator: it takes the
two VCFs, the reference and the options of the table above, and yields the rows of §6 in position
order, one per sample pair and compared locus, *reference in both* included; `--all-rows` only
decides which of them the command line writes.

```python
def compare(a: Path, b: Path, reference: Path, config: CompareConfig,
            evaluator: Evaluator = PVC) -> Iterator[OutcomeRow]: ...
```

`CompareConfig` holds `--regions`, `--strata`, `--str-regions`, `--sample-map` and
`--phasing-cap`.

## 4. The pipeline

Five stages, each a function from one stream to the next. Nothing holds more than the current locus
and a short look-ahead, so memory does not grow with the genome.

```
reference ──► segmented genome ────────────┐
file A ──► VCF records, normalised ──┐      │
                                     ├──► loci ──┬──► reconciled loci ──► outcomes
file B ──► VCF records, normalised ──┘           └──► another evaluator's decisions ──► outcomes (§8)
```

1. **Segmented genome** (§5.1): the STR finder, run per contig over the reference, yields the
   contig's generic regions and STR regions in position order. Only contigs that appear in a VCF
   are scanned.
2. **VCF records** (§5.2): each file is read as a stream of VCF records, checked against the
   reference, normalised, given its span, and kept in order of span start.
3. **Loci** (§5.3): the two streams of VCF records and the STR regions of the segmented genome are
   merged into a stream of loci, each holding the VCF records of both files that fall in it.
4. **Reconciliation** (§5.4): each locus gets its list of locus alleles, and each sample in each file
   its possible genotypes as indices into that list.
5. **Evaluation** (§5.5): each sample pair at each locus gets an outcome.

Stages 1 to 3 are the same for every evaluator, so every evaluator reports on the same loci and
their tables line up row by row. With another evaluator, stages 4 and 5 are replaced by its wrapper
(§8).

## 5. The rules

### 5.1 STR regions

pvc_eval ports the STR finder of pvc (repository `pop_var_caller`, at the commit recorded in
`doc/port_provenance.md`). pvc's own rules are in its `doc/devel/ng/spec/typed_regions.md` (§2.1–§2.5,
§4, §5), `repeat_catalog.md` (§3–§5), `ssr_repeat_scanner.md` (§3) and `str_bundle_loci.md` (§5);
where those documents and pvc's code disagree, the code is followed (pvc's `typed_regions_cli.md`
header example, for one, is out of date).

**What the finder does, in pvc's order.** For one contig's uppercased sequence:

1. **Find candidate repeats** (pvc `tandem_repeat::find_tandem_repeats`). For every period from 1 to
   6, each base scores +2 if it equals the base one period earlier and both are A, C, G or T, and −7
   otherwise (N never matches). The stretches of maximal total score (Ruzzo–Tompa's algorithm) are
   the candidates. A candidate is kept if it holds at least 2 copies of its period and its first
   period's bases are not themselves a repeat of a shorter unit; a candidate is dropped when a kept
   candidate of a shorter period that divides its period covers at least half of it.
2. **Bundles, then STR segments** (pvc `region_typing/mod.rs` `BundleMember`, its
   `typed_regions.md` §2.4; `repeat_catalog/row.rs`, `segment_criteria.rs`). A candidate with
   another within 15 bases on either side is an STR bundle member; the bundle is the hull of its
   members' detected spans, and its members go through none of the tests below. Any other candidate
   is admitted as an **STR segment**, trimmed to whole motifs, unless it has fewer copies than the
   floor for its period, its motif is compound, it has no clean whole-motif cut at both ends, its
   purity after the cut is below the minimum, or it lies within 15 bases of a contig end; a rejected
   candidate's bases stay in a generic region.
3. **Satellites** (pvc `region_typing::resolve_features`). A run of repeat longer than 100 bases is a
   satellite, and swallows every STR segment or bundle member within 15 bases of it, growing to cover
   them. A bundle whose span exceeds 100 bases becomes a satellite. A bundle touching a contig end is
   dropped.

The result is the segmented genome: pvc's partition of the contig into STR segments, STR bundles,
satellites and generic regions. Loci are built from its STR regions (§5.3), and `--strata pvc`
reads its region types (§5.7).

**Satellites are STR regions in pvc_eval**, although pvc does not call them STRs. An indel in a
satellite can be written in as many places as in a short repeat, so a VCF record touching a
satellite must pull it in like any other repeat. A satellite can be kilobases long and hold many
VCF records; what bounds the work is the phasing cap (§5.4), not the length.

**Parameters and defaults**, as in pvc unless said otherwise:

| parameter | default | pvc's source |
|---|---|---|
| periods | 1 to 6 | `segment_criteria.rs` `DEFAULT_MIN_PERIOD`, `DEFAULT_MAX_PERIOD` |
| match reward, mismatch penalty | 2, 7 | `tandem_repeat.rs` |
| copies a candidate needs | 2 | `cli/repeat_catalog.rs` `SCANNER_MIN_COPIES` |
| copy floor, periods 1 to 6 | **5, 5, 4, 4, 4, 3** (see below) | catalog floors, `repeat_catalog/criteria.rs`; pvc's calling floors are 8, 6, 6, 6, 5, 4 |
| minimum purity | 0.8 | `segment_criteria.rs` |
| bundle and flank distance | 15 bases | `DEFAULT_BUNDLE_THRESHOLD`, `CATALOG_MIN_FLANK_BP` |
| longest STR segment or bundle | 100 bases | `region_typing/mod.rs` `DEFAULT_MAX_STR_LEN`, `DEFAULT_MAX_STR_BUNDLE_LEN` |

**The copy floor is pvc's catalog floor, not its calling floor.** pvc's calling floors (8 copies of a
homopolymer, 6 of a dinucleotide, …) answer "does this repeat stutter enough to need the STR model?".
pvc_eval's question is "could an event here be written in more than one place?", which a repeat of
5 copies already raises. The lower floor admits more repeats, so more loci are widened to an STR
region and more have kind STR (§5.6); it changes no outcome of a locus that was already whole. Every
parameter can be set on the command line; the check against pvc (below) is run
with pvc's calling floors, to compare like with like.

**What is ported, and from where.** pvc runs the finder in two ways that give the same partition:
from a stored catalogue (a Parquet file built once per reference), and by scanning the sequence
directly (`region_typing::partition_resident_in`). pvc_eval ports the second, which needs no catalogue
file, plus the catalogue route's one extra rule: nothing within 15 bases of a contig end. About 900
lines of pvc's Rust do this: the candidate finder, the admission rules with their helpers
(`minimal_trim`, `recompute_purity`, `is_compound`, `split_bundles`), and the satellite and bundle
resolution (ported from its plain fixed-point form, `absorb_into`, which pvc's tests show equal to the
faster sweep it uses).

**Where the finder runs.** Scanning a whole human genome in Python would take far longer than a
benchmark comparison. So:

- `pvc_eval str-regions --reference ref.fa [--regions r.bed] --output str.bed` runs the finder once
  and writes STR regions as a BED, for `compare --str-regions` to reuse. Its fourth column is
  each region's type, `str`, `str-bundle` or `satellite`, which `--strata pvc` and the
  `str_regions` column (§6) read; its fifth is the motif of an STR segment, `.` for the other two.
  A `--str-regions` BED without those columns can be used for loci, and with `--strata pvc` it is
  an error. pvc gives the same answer whether or not a region list is given (its `typed_regions.md`
  §2.5), because it always scans whole contigs; the port does the same, and a region list only
  limits what is written.
- Without `--str-regions`, `compare` runs the finder itself over the contigs its VCFs touch.

How long the port takes on GRCh38 is measured first, before deciding whether the finder needs
vectorising (numpy) to be usable without the cached BED.

**Tests that come with the port.** pvc has about 140 tests of these functions (35 for the candidate
finder, 58 for the admission rules, 49 and two property tests for the partition). Each is rewritten
in Python with its inputs and expected outputs unchanged; their inline test sequences are copied as
they are. pvc's data-driven fixtures (`tests/data/tandem_repeat/synthetic_ref.fa` and the 16-locus
`golden.ssr_catalog.bed.gz`) are copied into pvc_eval's tests. Then the whole-port check of §7:
pvc_eval's STR regions and `pvc dev type-regions`'s, over the GIAB benchmark regions of GRCh38, must
be identical. `type-regions` does not scan: it reads the repeat catalogue pvc builds beside the
reference (pvc `cli/typed_regions.rs`, `RepeatCatalog::open_beside_reference`). So the check needs
a pvc binary and that catalogue; it is skipped, with a message naming what is missing, where either
is absent, and pvc_eval never needs pvc to run.

### 5.2 Reading and normalising VCF records

**Checks.** Each VCF record's contig must be in the reference, its REF must lie within the contig's
length, and records must come in the `.fai`'s contig order and by position; otherwise the run stops
with an error naming the file and line. A position past a contig's end means the VCF was made on
another reference, as a contig the reference lacks does (decided by the owner on 19 October 2026). A
VCF record whose REF differs from the reference at its position, compared without regard to case
(references are often soft-masked in lower case), is kept and marked *REF mismatch*; every locus
containing it is invalid for every sample of that file (§5.4).

**Normalising.** Remove the bases all of a VCF record's alleles but `*` share at the right end,
then at the left end, keeping at least one base in every allele (the VCF anchor convention). A
record with no ALT, or whose only ALTs are `*`, is not trimmed. The record is never moved: its locus alleles are
built from its VCF alleles where it was written (§5.4). Moving it, as `bcftools norm` does, would
break valid files: each record is left-aligned on its own, so a deletion can slide onto a SNP the
same haplotype carries (the example below).

**The span.** A VCF record's span is every position at which it could have been written, so that
two spellings of one event always overlap (§5.3):

- A record that is not an indel spans its REF after trimming.
- A deletion spans from its anchor base when slid as far left as it goes (its left-aligned form,
  as `bcftools norm` writes it) to the last base it deletes when slid as far right as it goes. It
  slides left while the base before it equals its last deleted base, and right while the base
  after it equals its first.
- An insertion spans likewise from the anchor of its left-aligned form to the base after the point
  it is inserted at when slid as far right as it goes. That last base is in the span even when the
  insertion cannot slide, so that `C → CT` at 10 overlaps a record at 11 (the third case below).
- A record with several ALTs spans the union of its ALTs' spans.
- A span always holds the record's trimmed REF, since its locus alleles are built where it was
  written.
- At the ends of a contig a span stops at the contig: an insertion after its last base has no base
  after it in its span, and a deletion slid to its first base has no anchor before it.

Every spelling of an indel with its shared base on the left has the same span. Any other spelling,
such as a deletion written with its shared base on the right, spans that span and its own trimmed
REF, which can reach one base further. So all spellings overlap and share a locus, and none is
moved (decided by the owner on 19 October 2026).

**Order.** A span can start before the record's POS, by as much as the deletion or insertion can
slide left, so the stream holds VCF records until the reader is further right than the longest
such distance and emits them in order of span start. A record whose span would reach more than
the span limit beyond its REF is given its trimmed REF as its span and reported in a warning. The
span limit is 10 kb unless the user sets another whole number of 1 or more, or none, with
`--span-limit` or in the Python API (decided by the owner on 19 October 2026; a limit of 0 would
cut every insertion, whose span always holds the base after it). With no limit no span is cut, and since no distance then bounds how far
left a later record can start, the stream holds the records of a contig until the next contig
begins.

**Threads.** The reader parses the sample columns of batches of VCF records on a pool of threads,
4 by default, and 1 parses them in the calling thread with no pool; whatever the number, it gives
the same records, in the same order, with the same errors and warnings (decided by the owner on
19 October 2026, for a measured gain of about three times on 200 samples). The number is set in
the Python API; whether the command line sets it too, and how two files read together share the
cores, is not decided yet.

**Three cases.** The first is valid, and would be made invalid by moving its deletion; the other
two are one haplotype written twice, which only the span puts in one locus.

```
reference 21–26  CAGAGT                    the AG at 22–23 repeats at 24–25
22  A    T    1|0                          one haplotype: SNP at 22 …
23  GAG  G    1|0                          … and loss of 24–25: CTGT
```

`bcftools norm` writes the deletion as `21 CAG C`, deleting 22–23, which collides with the SNP at
22. Kept where it was written, the haplotype is `CTGT` and the sample is valid. The deletion's span
is 21–25, so both records are in one locus.

```
reference 20–27  GCAGAGTC
file A   21  CAG  C   1/1                  span 21–25: the deletion of 22–23 can slide to 24–25
file A   26  T    C   1/1
file B   24  AGT  C   1/1                  span 24–26
```

Both files give `GCAGCC`. Without the span, A's deletion (21–23) and B's record (24–26) would fall
into different loci, giving *only in A* at one and *different* at the other.

```
reference 10–11  CG
file A   10  C    CT  1/1                  span 10–11
file A   11  G    A   1/1
file B   11  G    TA  1/1                  span 11
```

Both files give `CTA`, in one locus because A's insertion spans the base after its anchor.

### 5.3 Loci

**Building a locus.** Starting from the leftmost VCF record not yet in a locus, of either file:

1. The locus is that record's span (§5.2).
2. Add every VCF record, of either file, whose span overlaps the locus, and widen the locus to cover
   its span. Repeat until no more records overlap. Overlap chains: if A overlaps B and B overlaps C,
   all three are in the locus, whichever files they come from.
3. Widen the locus to cover every STR region a VCF record in it **touches** — overlaps or lies
   immediately next to. Go back to step 2, since the wider locus may now overlap more records. Stop
   when neither step adds anything.

A locus is a run of reference positions, written 1-based and inclusive, as in VCF.

**Why STR regions widen loci too.** The span (§5.2) makes the spellings of one indel overlap. A
complex record — one that changes more than one base in an allele — has no span beyond its REF,
and two callers can write the same events in a repeat as complex records at different places in
it. Taking the whole repeat whenever a VCF record touches it puts every spelling into one locus.

### 5.4 Reconciliation

For each locus, each file and each sample, pvc_eval builds the sample's possible genotypes, then
expresses the possible genotypes of both files over one shared list of locus alleles.

**Building a locus allele.** A locus allele is built for one haplotype from the reference over the
locus and the VCF alleles that haplotype carries, each at the position it was written at (§5.2):

- Each VCF allele other than REF is an **edit**: after trimming the bases it shares with its record's
  REF, first at the left end and then at the right, it replaces some reference positions (none, for
  an insertion) with some bases (none, for a deletion). Trimming first matters: a VCF record widened
  over flank bases for one allele's sake spells the reference there in its other alleles, and that
  is no edit. The order matters too: `52 AG → A,AAG` trimmed from the left inserts its A after
  position 52, and from the right before it, which differ when the haplotype changes position 52.
  Left first puts the insertion after the anchor, as VCF does.
- A `*` allele makes no edit: it says the haplotype's bases here are removed by another VCF record,
  which carries that edit. The haplotype must carry a deletion, from another VCF record of the
  locus, covering a position of the `*` record's REF; if it does not, the haplotype is
  **impossible**.
- The edits of one haplotype are applied to the reference over the locus in position order. Two
  edits that change the same reference position, or two insertions at the same point, cannot both
  apply: the haplotype is **impossible**.
- If any VCF genotype leaves the haplotype missing (`.`), the locus allele is missing, however many
  other VCF records give it a definite allele.

**Possible genotypes.** A sample's VCF genotypes across the locus's VCF records say which VCF allele
each haplotype carries, up to phase.

- VCF genotypes written with `|` and the same `PS` (or no `PS`, which VCF reads as one phase set per
  sample and contig) fix which VCF alleles go together on one haplotype.
- Between unphased VCF records, or different phase sets, every assignment of VCF alleles to
  haplotypes is tried. A heterozygous unphased VCF record, or a phase set holding a heterozygous
  VCF genotype, is one **block**; with ploidy 2 and *b* blocks there are 2^(b−1) assignments.
  Homozygous VCF genotypes, and records of one phase set, add none.

**The phasing cap.** A sample whose assignments at a locus would exceed the phasing cap (default
4096, that is 13 blocks; **Open 2**) is not reconciled there, and its sample pairs get the outcome
*too complex* (§5.5). The count is made before any assignment is built. This bounds the work,
which doubles with every block, and only the work: a satellite holding hundreds of phased VCF
records is reconciled, and a locus of a few bases with 14 unphased heterozygous records is not.

Each assignment gives one locus genotype, unless one of its haplotypes is impossible. The sample's
possible genotypes are the distinct locus genotypes that result. A sample whose VCF genotypes allow
none is **invalid** at the locus.

**A sample is also invalid** at a locus if any of its VCF genotypes there is malformed — a number of
alleles other than the ploidy, or an allele index the VCF record lacks — or if any VCF record of its
file in the locus has a REF mismatch (§5.2). A haploid VCF genotype, such as `1` on chrX of a male
sample, is not malformed but makes the sample invalid too, with its own reason (**Open 1**).

**The example of an invalid sample.** Reference positions 9 to 13 read `CAGTC`. One file has

```
POS REF  ALT  GT
9   CAGT C    1/1     both haplotypes lose AGT (positions 10–12)
11  G    A    0/1     one haplotype has A at position 11
```

Position 11 no longer exists on either haplotype, yet one is said to carry an A there: no assignment
works, so the sample is invalid at this locus. The locus is 9 to 12, the deletion's span. With the
deletion `0/1` the same records are fine — one haplotype `C`, the other `CAAT` — and GATK writes
the second record as `G → A,*` with `1/2`, which says exactly that.

**Shared locus alleles.** The locus's alleles are REF followed by every distinct non-missing locus
allele appearing in any possible genotype of any sample of either file, in order of first appearance.
Each possible genotype becomes a sorted tuple of indices into that list, with missing haplotypes
written as missing. This is what lets §5.5 compare genotypes by number: two files that spelled one
sequence differently now name it by the same index.

### 5.5 Evaluation

Each sample pair at each compared locus gets exactly one outcome. The first rule that applies wins:

| outcome | when |
|---|---|
| **too complex** | either file's sample exceeds the phasing cap at the locus and was not reconciled (§5.4) |
| **invalid in A**, **invalid in B**, **invalid in both** | that file's sample is invalid at the locus |
| **missing in A**, **missing in B**, **missing in both** | every possible genotype of that file's sample has a missing haplotype |
| **reference in both** | in both files, the only possible genotype is all REF |
| **same genotype** | some possible genotype of A equals some possible genotype of B |
| **same alleles, other genotype** | some possible genotype of A and some of B carry the same set of non-reference locus alleles, on different numbers of haplotypes (one in A, both in B, say) |
| **only in A**, **only in B** | one file gives the sample a non-reference locus allele, and the other's only possible genotype is all REF |
| **different** | both give the sample non-reference locus alleles, and no pair of possible genotypes matches as above |

When several possible genotypes exist, the pair reported in the table is the one that decided the
outcome; for *different*, the first in order.

*Reference in both* is not one of the outcomes agreed beforehand: it is needed because a locus exists
wherever any sample of either file has a variant, so in a cohort most samples at most loci are
reference in both files. It is written to the table only with `--all-rows`; the summary always
counts it.

### 5.6 Kind

The kind of a locus is read from the VCF records of both files in it, after normalising (§5.2),
and not from its locus alleles, so that it is defined at every locus: one that is *too complex*,
one compared by another evaluator (§8), and one where no sample carries an ALT. `*` alleles are
left out.

| kind | when |
|---|---|
| **STR** | the locus contains an STR region |
| **SNP** | every ALT of every record is one base replacing a one-base REF |
| **indel** | every ALT of every record is one insertion or one deletion after the anchor base |
| **complex** | anything else: a complex record, or records of different kinds |

### 5.7 Strata from pvc's region types

`--strata pvc` sorts loci by the part of pvc's partition (§5.1) they fall in, which is the
stratification most comparisons of pvc want. The four classes are the region types (§2), pvc's
region kinds with pvc_eval's names: `generic` (a generic region; pvc `RegionKind::Generic`), `str`
(an STR segment; `ssr_locus` in pvc's typed-regions BED), `str-bundle` (`ssr_bundle`) and
`satellite`.

A locus has exactly one of the four classes, the first that applies, from the most complex STR
region to the simplest:

| stratum | when the locus contains |
|---|---|
| `satellite` | a satellite |
| `str-bundle` | an STR bundle, and no satellite |
| `str` | an STR segment, and no satellite or STR bundle |
| `generic` | no STR region |

The generic regions a locus overlaps beside its STR regions do not count. The rule for a
stratification BED, every class the locus overlaps, would make nearly every STR locus
`generic,str`: a locus widened to an STR region (§5.3) usually keeps the base before it, where an
insertion's span starts.

For example, an AT repeat is an STR segment at positions 11 to 20, and a VCF record inserts `AT`
after position 10. The locus is 10 to 20; its stratum is `str`, not `generic,str`. A locus whose
records chain that STR segment to a satellite 20 bases away is `satellite` (within 15 bases the
satellite would have swallowed the segment, §5.1), and one that chains it to an STR bundle is
`str-bundle`.

The STR regions are the ones the loci were built from, the finder's or `--str-regions`'s, so the
stratum always agrees with the `str_regions` column (§6). It refines the kind of §5.6: a locus of
kind STR is never `generic`, and every other kind is.

## 6. The output table

Tab-separated, one row per sample pair and locus (without `--all-rows`, rows whose outcome is
*reference in both* are left out). Columns:

| column | content |
|---|---|
| `chrom`, `start`, `end` | the locus, 1-based and inclusive |
| `kind`, `stratum` | §5.6; `.` when `--strata` is not given |
| `str_regions` | the STR regions in the locus, as `start-end:type`, with `:motif` after an STR segment's, joined with `,`; `.` for none |
| `sample_a`, `sample_b` | the sample pair |
| `outcome` | §5.5 |
| `alleles` | the locus alleles, REF first, joined with `,` (a deletion of the whole locus is written `-`) |
| `genotype_a`, `genotype_b` | the reported possible genotypes, as indices into `alleles` joined with `/`, `.` for a missing haplotype |
| `possible_a`, `possible_b` | how many possible genotypes the sample had in each file |
| `quality_a`, `quality_b` | §2, *quality*; `.` when the file gives the sample no non-reference VCF allele |
| `gq_a`, `gq_b` | the lowest GQ among the sample's VCF records in the locus that give it a non-reference VCF allele, `.` when absent |
| `records_a`, `records_b` | how many VCF records of each file the locus holds |
| `invalid_reason_a`, `invalid_reason_b` | for an invalid outcome: `contradicting records`, `REF mismatch`, `malformed GT` or `haploid GT`; `.` otherwise |

The summary (`--summary`) counts rows by outcome, kind and stratum.

## 7. Tests

Nothing about pvc_eval is trusted until these pass.

1. **The STR finder against pvc:** the port's own copy of pvc's tests (§5.1), and identical STR
   regions, types included, from the port and from `pvc dev type-regions` over the GIAB benchmark
   regions of GRCh38 (with pvc's repeat catalogue, §5.1).
2. **Each stage on hand-written cases**, including:
   - normalising and spans: trimming at both ends; a `*`-only record, not trimmed; deletion and
     insertion spans in homopolymers and longer repeats; a multi-ALT record; a complex record; the
     three cases of §5.2; a span cut at 10 kb;
   - loci: chaining across both files; an insertion next to a repeat pulling it in; two
     complex records at the two ends of a repeat; with `--regions`, a locus
     sticking out of a region compared whole, and one sharing no position with any region skipped;
   - reconciliation: the invalid example of §5.4 and its valid variant; GATK's `*` form, and a `*`
     on a haplotype with no deletion (impossible); the insertion of `52 AG → A,AAG` after the
     anchor; phase sets; the phasing cap at 13 blocks (reconciled) and 14 (*too complex*), and 14
     heterozygous records in one phase set (one block, reconciled); `./1`; a REF mismatch, and a
     lower-case reference that is no mismatch; malformed GT; haploid GT;
   - kind: one locus of each kind, and a *too complex* locus, which still has one;
   - evaluation: one case per outcome;
   - strata with `--strata pvc`: the insertion of §5.7 (`str`), a locus with no STR region
     (`generic`), the STR segment chained to a satellite (`satellite`) and to a bundle
     (`str-bundle`), a bundle chained to a satellite (`satellite`), and a `--str-regions` BED
     without a type column (an error).
3. **A file against itself:** every row *same genotype*, *reference in both*, *missing in both*,
   *invalid in both* or *too complex*; the last three only where the file itself has a missing
   haplotype, an invalid sample or a sample over the phasing cap.
4. **A file against re-spellings of itself** — every VCF record decomposed into single-base and
   single-indel records with genotypes carried over; neighbouring records merged into complex ones;
   indels shifted right within their repeats; GATK-style `*` records added — every row as for the
   file against itself, with the same outcome at the same locus. The re-spelling program is part of
   the tests.
5. **The hard cases from pvc's STR work:** the VCF records behind the cases recorded in pvc's
   reports (GIAB's split indel-plus-SNP truth records against pvc's single repeat record; GATK's
   overlapping records at STR edge deletions; HipSTR's whole-repeat records). They are copied into
   this spec, each with the outcome it must give, before it is approved (**Open 3**).

## 8. Comparing evaluators

Other evaluators run behind `compare`'s interface (§3): the same two VCFs, reference and options in,
the table of §6 out. There are three: RTG vcfeval, hap.py, and naive matching (`bcftools norm`, then
matching VCF records by position and alleles).

**A wrapper yields outcomes, not reconciled loci.** vcfeval and hap.py do not build locus alleles
or possible genotypes; they give a decision for each VCF record (matched, or unmatched in one
file), with some grades of partial match. A wrapper that filled reconciled loci from those would be
inventing what the tool never computed. So a wrapper takes the loci of §4 stage 3, runs its tool
once per sample pair (vcfeval and hap.py compare one sample at a time), and turns the decisions of
the VCF records in each locus into one outcome of §5.5. The columns the tool has nothing for —
`alleles`, `genotype_a`, `genotype_b`, `possible_a`, `possible_b` — are `.`.

**What the wrappers are for.** Once pvc_eval passes §7, the same comparisons are run with every
evaluator, and every locus where they disagree is classified by hand. That classification is how
pvc_eval's rules are checked against the established tools'.

hap.py ships only for x86-64; on an arm64 machine it runs from its Docker image under emulation,
which is slow but enough for benchmark regions. How each wrapper turns a locus's record decisions
into an outcome — a locus whose records got different decisions, the outcomes a tool cannot tell
apart, the A/B asymmetry of tools that call one file the truth — is its own short spec under
`doc/specs/`, written after this one is settled.

## 9. Not covered

- **Speed at cohort scale.** pvc_eval is Python and meant for benchmarks of megabases and up to a
  few hundred samples. Whole genomes of thousands of samples would need another implementation,
  written against the same tests.
- **Structural variants and symbolic alleles** (`<DEL>`, breakends): a VCF record with one is
  skipped and counted in a warning.
- **Ploidy other than 2**: functions are written for any ploidy where that costs nothing; where one
  is limited to 2 (the phase enumeration of §5.4 is the likely one), it raises an error for any
  other ploidy.

## 10. Open points

**Open 1. Haploid genotypes.** A haploid VCF genotype (`1`, `0`) is what callers write on chrX and
chrY of male samples, and on chrM. The options:

- *Invalid*, as now: the sample is invalid at the locus, with the reason `haploid GT`. Nothing is
  compared on those contigs, and the table says why.
- *Read as homozygous*: `1` is compared as `1/1`. A file that writes male chrX haploid then agrees
  with one that writes it diploid and homozygous, which is what both mean; a heterozygous call on
  male chrX shows as *same alleles, other genotype*.
- *Support ploidy 1*: a locus genotype of one haplotype, compared only with another of one. Most
  faithful, and the most work; §9 leaves ploidy other than 2 out.

I recommend reading it as homozygous: it compares what the files mean, at the cost of one rule.
Meanwhile the implementer makes it invalid, which the other two can replace without touching
anything else.

**Open 2. The phasing cap's default.** 4096 assignments is 13 blocks; how long they take per sample
and locus is not measured yet. A higher cap compares more loci in unphased cohort data and
costs time exponentially; a lower one gives more *too complex* rows. I recommend 4096 until §8's
comparisons show how many rows it costs. Meanwhile the implementer uses 4096.

**Open 3. The hard cases of §7.5.** The VCF records of the cases in pvc's reports, and the outcome
each must give, are to be copied into §7 from pvc's reports. Until then §7.5 checks nothing; the
implementer leaves it out of the first plan.
