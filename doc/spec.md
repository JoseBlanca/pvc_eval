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
| **ploidy** | the number of haplotypes each sample of a file has at every locus: 1 or 2, given per file, 2 by default | |

**From the reference**

| term | meaning | not |
|---|---|---|
| **segmented genome** | each contig of the reference cut by the STR finder (§5.1) into consecutive, non-overlapping regions, each a generic region or an STR region | |
| **generic region** | a stretch of the segmented genome that holds no tandem repeat the STR finder admits; pvc's `generic` | "ordinary sequence" |
| **STR region** | a stretch the STR finder (§5.1) reports as a tandem repeat: an **STR segment** (one repeat), an **STR bundle** (two or more repeats with fewer than 15 bases between each and the next) or a **satellite** (a repeat run longer than 100 bases). pvc's terms. | "tract", "STR locus", "repeat locus" |
| **region type** | which of the four a region of the segmented genome is: `generic`, `str` (an STR segment), `str-bundle` or `satellite` | "kind", which is §5.6's |

**Building what is compared**

| term | meaning | not |
|---|---|---|
| **locus** | the unit compared: a stretch of reference holding VCF records of either file, chained by overlap or by the chain distance and widened to every STR region a VCF record touches (§5.3) | any other use of "locus" |
| **chain distance** | the most reference positions that may lie between a locus and the span of a VCF record it still takes in (§5.3) | |
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
                 [--sample-map pairs.tsv] [--ploidy-a N] [--ploidy-b N] [--phasing-cap N]
                 [--span-limit N|none] [--chain-distance N]
                 [--summary summary.tsv] [--all-rows]
```

`compare()`, the command, the table writer and the summary are given in full in
`doc/specs/compare.md`, approved by the owner on 20 October 2026, with how its threads scan only
the contigs whose loci need STR regions.

| option | what it does |
|---|---|
| `--reference` | FASTA with its `.fai`. Contig names and order come from the `.fai`. |
| `--a`, `--b` | the two VCFs, plain or bgzipped. |
| `--regions` | compare only the loci that overlap these regions, by at least one position; every other locus is skipped and gets no row, nor a count in the summary. A locus that overlaps is compared whole, its VCF records outside the regions included, since loci are built from all VCF records (§5.3). The VCF records of a contig it does not name are not read beyond their contig and position (`doc/specs/regions_read.md`). Without it, every locus is compared. |
| `--strata` | BED with a fourth column naming each interval's class; a locus takes the class of the intervals it overlaps (several are joined with `,`; none gives `.`). The word `pvc` instead of a path gives each locus one region type from the STR regions it holds — `satellite`, `str-bundle`, `str` or `generic` — with no BED needed (§5.7); a file named `pvc` is given as `./pvc`. |
| `--str-regions` | BED of STR regions to use instead of running the STR finder. |
| `--sample-map` | two columns: a sample of A, a sample of B. Without it, samples are paired by name. |
| `--ploidy-a`, `--ploidy-b` | the ploidy each file is read with, 1 or 2, default 2 (§5.4). |
| `--phasing-cap` | default 4096 (§5.4). |
| `--chain-distance` | default 10 (§5.3): a whole number of 0 or more. |
| `--span-limit` | how far, in bases, a VCF record's span may reach beyond its REF before it is cut back to the REF (§5.2): a whole number of 1 or more, default 10000, or `none` for no limit. |
| `--output` | the table of §6, one row per sample pair and locus. |
| `--summary` | counts of the rows over all sample pairs, one line per outcome and stratum, one column per kind (§6). |
| `--all-rows` | also write the rows whose outcome is *reference in both* (§5.5). |

**Errors, which stop the run:** a VCF not sorted in the `.fai`'s contig order and by position; a
contig in a VCF but not in the reference; a VCF record whose REF reaches past the end of its
contig, or whose REF differs from the reference (§5.2); an ALT `.` among other ALTs (§5.2); a
gzipped VCF or BED that is damaged (§5.2); no sample pairs; a sample named in `--sample-map` absent
from its file; a line of `--sample-map` that is not two fields; a pair given twice; a sample name
given twice in a header, when a pair needs it. Samples present in only one file and not mapped are listed as a warning and not
compared.

**Not errors:** everything else about the content of a VCF record — a malformed GT or one of
another ploidy, contradicting VCF records. Those make the affected locus
*invalid* for the affected sample (§5.4) and the run goes on.

The Python API offers each stage of §4 as a function over iterators, so a user can, for example,
filter VCF records between reading and grouping, or read loci and genotypes without evaluating them.
Above the stages, one function does what `compare` does: it takes the
two VCFs, the reference and the options of the table above, and yields the rows of §6 in position
order, one per sample pair and compared locus, *reference in both* included; `--all-rows` only
decides which of them the command line writes.

```python
def compare(a: Path, b: Path, reference: Path, config: CompareConfig = CompareConfig(), *,
            finder_threads: int = DEFAULT_FINDER_THREADS) -> Iterator[OutcomeRow]: ...
```

`--evaluator`, which would choose another evaluator of §8, and an `evaluator` parameter of
`compare()` are left out until a wrapper of §8 exists (decided by the owner on 20 October 2026):
adding them then breaks nobody, since their default is pvc_eval's own rules.

`CompareConfig` holds `--regions`, `--strata`, `--str-regions`, `--sample-map`, `--ploidy-a`,
`--ploidy-b`, `--phasing-cap`,
`--span-limit` and `--chain-distance`.

## 4. The pipeline

Five stages, each a function from one stream to the next. Nothing holds more than the current locus
and a short look-ahead, so memory does not grow with the genome; with `--span-limit none` the
look-ahead is a whole contig (§5.2).

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

**What stage 1 yields.** The STR finder and the loci are built apart, so the type that passes
between them is fixed here, in `pvc_eval.segments`, and changed only in this spec:

```python
class RegionType(StrEnum):
    GENERIC = "generic"
    STR = "str"
    STR_BUNDLE = "str-bundle"
    SATELLITE = "satellite"

@dataclass(frozen=True, slots=True)
class GenomeSegment:
    contig: str
    start: int                # 0-based
    end: int                  # one past the last position; end > start
    region_type: RegionType | None   # None only from a --str-regions BED without types
    motif: str | None         # an STR segment's motif, upper case; None for the other three

def segment_contig(reference: IndexedFasta, contig: str,
                   parameters: StrFinderParameters = ...) -> Iterator[GenomeSegment]: ...
def read_str_regions(path: Path, reference: IndexedFasta) -> Iterator[GenomeSegment]: ...
```

Coordinates are 0-based and half-open, as `VcfRecord.start` and `VcfRecord.span` are.
`segment_contig` yields one contig's segments in position order, consecutive and non-overlapping,
from 0 to the contig's length. `read_str_regions` reads a `--str-regions` BED (§5.1) and yields its
STR regions only, sorted by contig in the reference's `.fai` order and by start. A BED whose
regions overlap, or name a contig the reference lacks, stops the run with an error naming the
line. The loci (§5.3) take an iterator of either: they ignore `GENERIC` segments and never rely
on the segments covering a contig. `StrFinderParameters` holds the parameters of §5.1, and the finder's own spec gives it.

**What stage 3 yields.** The loci and the reconciliation are built apart too, so the type that
passes between them is fixed here, in `pvc_eval.locus`, a module that holds only this type, and
changed only in this spec:

```python
@dataclass(frozen=True, slots=True)
class Locus:
    contig: str
    start: int                              # 0-based
    end: int                                # one past the last position
    records_a: tuple[VcfRecord, ...]        # file A's records in the locus
    records_b: tuple[VcfRecord, ...]        # file B's
    str_regions: tuple[GenomeSegment, ...]  # the STR regions the locus contains

def build_loci(records_a: Iterator[VcfRecord], records_b: Iterator[VcfRecord],
               segments: Iterator[GenomeSegment], contigs: Sequence[str], *,
               chain_distance: int = 10) -> Iterator[Locus]: ...
```

`build_loci` is in `pvc_eval.loci`, and the loci's own spec, `doc/specs/loci.md`, gives it in
full. `contigs` is the reference's contig names in the `.fai`'s order (`IndexedFasta.contigs`),
which orders the two files' contigs and the segments' between them. A `Locus` covers the spans of
all its records with an ALT and every STR region it holds (§5.3), so `start` and `end` hold their
REFs. Each file's records are in the order `parse_vcf` yields them, by span start. Records with no
ALT are among them: one is in every locus its REF overlaps, and its REF may reach beyond the locus
(§5.2). `str_regions` holds no `GENERIC` segment, in position order. Loci come
contig by contig in the `.fai`'s order and by start, and do not overlap. The locus knows nothing
of samples: each record's genotypes are in the order of its file's header (`parse_vcf_samples`),
and the reference bases over the locus are read from the reference.

Stages 1 to 3 are the same for every evaluator, so every evaluator reports on the same loci and
their tables line up row by row. With another evaluator, stages 4 and 5 are replaced by its wrapper
(§8).

## 5. The rules

### 5.1 STR regions

pvc_eval ports the STR finder of pvc (repository `pop_var_caller`, at the commit recorded in
`doc/port_provenance.md`). pvc's own rules are in its `doc/devel/ng/spec/typed_regions.md` (§2.1–§2.5,
§4, §5), `repeat_catalog.md` (§3–§5), `ssr_repeat_scanner.md` (§3) and `str_bundle_loci.md` (§5);
where those documents and pvc's code disagree, the code is followed (pvc's `typed_regions_cli.md`
header example, for one, is out of date). The finder's own spec, `doc/specs/str_finder.md`, gives
each step exactly, with pvc's function for each, its parameters, the `str-regions` command and how
the port is checked.

**What the finder does, in pvc's order.** For one contig's uppercased sequence:

1. **Find candidate repeats** (pvc `tandem_repeat::find_tandem_repeats`). For every period from 1 to
   6, each base scores +2 if it equals the base one period earlier and both are A, C, G or T, and −7
   otherwise (N never matches). The stretches of maximal total score (Ruzzo–Tompa's algorithm) are
   the candidates. A candidate is kept if it holds at least 2 copies of its period and its first
   period's bases are not themselves a repeat of a shorter unit; a candidate is dropped when a kept
   candidate of a shorter period that divides its period covers at least half of it.
2. **Filter, bundle, then admit STR segments** (pvc `region_typing/segment_criteria.rs`:
   `prefilter`, the classification before `split_bundles`, then `finish_locus`). In this order:
   - A candidate with fewer than 15 bases between it and either end of the contig, measured on its
     whole-motif cut or, when it has none, on its detected span, is dropped (pvc's catalogue route,
     `repeat_catalog/row.rs` `row_for_interval`). It bundles nothing and counts towards no
     satellite.
   - A candidate is dropped if its period is outside the range, its score is under the floor, it
     has fewer copies than the floor for its period as detected, or its motif is compound.
   - Of the candidates left, one with fewer than 15 bases between it and another, on either side,
     is an STR bundle member. The bundle is the hull of its members' detected spans. Members go through none of
     the tests below.
   - Any other candidate is admitted as an **STR segment**, trimmed to whole motifs, unless it has
     no clean whole-motif cut at both ends, fewer copies than the floor after the cut, purity after
     the cut below the minimum, or a cut that starts on the contig's first base or ends on its
     last.

   A dropped or rejected candidate is no STR segment and no bundle member, but it still counts
   towards satellites (step 3): pvc builds satellites from the runs of every candidate that passed
   the copy floor as detected (`region_typing/mod.rs`, `coverage_runs`), whatever later tests
   rejected it. Its bases are in a generic region only if no satellite covers them.
3. **Satellites** (pvc `region_typing::resolve_features`). A run of repeat longer than 100 bases is a
   satellite, and swallows every STR segment or bundle member that overlaps it or has fewer than 15
   bases between them, growing to cover them. A bundle whose span exceeds 100 bases becomes a satellite. A bundle touching a contig end is
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

**What is ported, and from where.** pvc runs the finder in two ways: from a stored catalogue (a
Parquet file built once per reference), which `pvc dev type-regions` reads, and by scanning the
sequence directly (`region_typing::partition_resident_in`). They differ in two rules: the catalogue
route drops candidates near a contig end (step 2), and it always scans periods 1 to 6, whatever
periods it then admits. pvc_eval ports the scanning route, which needs no catalogue file, with both
of the catalogue route's rules, so that it gives what `dev type-regions` gives. About 900
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
  an error; its regions are written in the `str_regions` column (§6) as `start-end` alone. The port always segments whole contigs,
  so a region list only limits what is written and never changes a region; pvc's
  `dev type-regions --regions` can return a shorter satellite than a whole-genome run
  (`doc/specs/str_finder.md`, "The same check over GRCh38").
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

**Reading.** A VCF, plain or gzipped (known by its first two bytes), is read as UTF-8 text; a
byte-order mark at its start is dropped, as for a BED, and a gzip file cut short or damaged stops
the run with a message naming the line it was reached at (decided by the owner on 19 October
2026). gzip checks its data only at the end of each compressed block, so damaged bytes reach the
reader as text first, and a line made of them would be reported as a fault of the file's content,
a column missing or records out of order: on 21 October 2026, of 5,087 copies of one gzipped VCF
each with one byte changed, 4,202 were reported so. So when reading a gzipped VCF or BED stops at
a fault of its content, the file is first checked whole, and if gzip finds it damaged, the run
stops with a message saying that the file is damaged and must be compressed or copied again, in
place of the fault of the line (decided by the owner on 21 October 2026). This holds also where
the fault is found after the line is read, by sample pairing, `--strata` or `--regions`. A good
file is read once; the whole check costs only a read that already failed. A BGZF file, as `bgzip`
and `bcftools` write, ends with an empty block that marks its end; one cut between two blocks
reads as a valid shorter file, so a BGZF file without that block is cut short, as htslib also
says (found by the review of 21 October 2026).

**The ALT `.`.** `.` alone as the ALT is a record with no ALT (below). `.` as one of several ALTs
(`G A,.`) means nothing in VCF, and stops the run with an error naming the file and the line
(decided by the owner on 21 October 2026); it was read as an allele. The rule holds also for a
record with a symbolic ALT, `<DEL>,.`, which is otherwise skipped.

**Checks.** Each VCF record's contig must be in the reference, its REF must lie within the contig's
length, and records must come in the `.fai`'s contig order and by position; otherwise the run stops
with an error naming the file and line. A position past a contig's end means the VCF was made on
another reference, as a contig the reference lacks does, and so do contigs out of the `.fai`'s
order: the VCFs, the FASTA and its `.fai` must agree (both decided by the owner on 19 October
2026). A
VCF record whose REF differs from the reference at its position, compared without regard to case
(references are often soft-masked in lower case), stops the run too, with an error naming the file,
the line, the position, the REF and the reference's bases there: a VCF that disagrees with its
reference was made on another one, or is damaged, and no comparison over it can be trusted (decided
by the owner on 19 October 2026, in place of marking the record and making its loci invalid).
Two bases agree when the bases their IUPAC codes stand for overlap: `N` agrees with any base, `R`
with `A` or `G`, `A` only with codes that include `A`; so a reference that writes `M` or `N` at a
site where the VCF has a plain base, as some references and callers do, is no difference, while a
VCF made on another build still is (decided by the owner on 19 October 2026). Only the letters of
ASCII are upper-cased for the comparison: a letter outside ASCII is no IUPAC code and agrees with
no base, so a REF holding one differs from the reference. Upper-cased by Python's rules, `ß` became
`SS`, one base longer, and stopped the run with a traceback, and `ſ` became `S` and agreed with `G`
(found by the review of milestone 1 of `doc/plans/hard-cases.md`, 21 October 2026). The error says how
many records of the file agreed with the reference before the failing one, which tells a wrong
reference (none) from a damaged line (many); with `--regions`, only the records of the contigs it
names are compared with the reference and counted (`doc/specs/regions_read.md`). A record with a symbolic allele (§9) is checked only
for its contig and its order, not for its REF or its length, since it is skipped and nothing is
compared from it (decided by the owner on 19 October 2026); a file made on another reference still
stops at its first record with plain alleles.

**Normalising.** Remove the bases all of a VCF record's alleles but `*` share at the right end,
then at the left end, keeping at least one base in every allele (the VCF anchor convention). A
record with no ALT, or whose only ALTs are `*`, is not trimmed. The record is never moved: its
locus alleles are built from its trimmed VCF alleles where it stands (§5.4). Moving it, as `bcftools norm` does, would
break valid files: each record is left-aligned on its own, so a deletion can slide onto a SNP the
same haplotype carries (the first case below), or past one, which changes the haplotype without
any sign of it:

```
reference 6–13  CTAAAAGT
9   A    G    1|0                          the haplotype is CTAGAGT
10  AA   A    1|0
```

`bcftools norm` (1.24) writes the deletion as `7 TA T`, and the haplotype becomes `CTGAAGT`.
So pvc_eval is best given the VCFs as their callers wrote them. A file already run through
`bcftools norm` is compared as it now reads: where the move made records collide, the sample is
invalid there (§5.4); where it changed a haplotype silently, pvc_eval cannot know (decided by the
owner on 19 October 2026, over trying every placement of colliding indels, which guessed among
haplotypes that differ).

**The span.** A VCF record's span is every position at which it could have been written, so that
two spellings of one event always overlap (§5.3):

- A record that is not an indel spans its REF after trimming.
- A deletion spans from its anchor base when slid as far left as it goes (its left-aligned form,
  as `bcftools norm` writes it) to the last base it deletes when slid as far right as it goes. It
  slides left while the base before it equals its last deleted base, and right while the base
  after it equals its first.
- An insertion spans likewise from the anchor of its left-aligned form to the base after the point
  it is inserted at when slid as far right as it goes. An insertion slides right while the base
  after it equals its first inserted base, which then moves to its end (`CA` inserted before `C`
  becomes `AC` inserted after it), and left while the base before it equals its last inserted
  base, which moves to its start. The last base of the span is in it even when the insertion
  cannot slide, so that `C → CT` at 10 overlaps a record at 11 (the third case below).
- A record with several ALTs spans the union of its ALTs' spans.
- A span always holds the record's trimmed REF, since its locus alleles are built where it was
  written. The span of a record with a `*` allele also holds its REF as the file wrote it, before
  trimming, also when the span limit cuts the span back: a `*` names those bases (§5.4), and the
  deletion that covers it must fall in the same locus (decided by the owner on 21 October 2026,
  over leaving such a deletion to another locus, where the `*` was judged uncovered).
- At the ends of a contig a span stops at the contig: an insertion after its last base has no base
  after it in its span, and a deletion slid to its first base has no anchor before it.

Every spelling of an indel with its shared base on the left has the same span. Any other spelling,
such as a deletion written with its shared base on the right, spans that span and its own trimmed
REF, which can reach one base further. So all spellings overlap and share a locus, and none is
moved (decided by the owner on 19 October 2026).

**Records with no ALT.** A VCF record whose ALT is `.`, as all-sites VCFs write, starts no locus
and widens none: it takes no part in §5.3. It is in every locus its REF overlaps, and there it
gives its samples' VCF genotypes like any record: `0/0` is REF, `./.` is missing. It is left out of
the kind (§5.6) and of the `records_a`, `records_b` counts (§6).

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
two are one haplotype written twice, whose spellings overlap only through the span.

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

Both files give `GCAGCC`. By overlap alone, without the span, A's deletion (21–23) and B's record
(24–26) would fall into different loci, giving *only in A* at one and *different* at the other.
They touch, so the chain distance (§5.3) joins them too; `doc/specs/loci.md` has a case that only
the span joins.

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
2. Add every VCF record not yet in a locus, of either file, whose span overlaps the locus or has
   at most the chain distance of reference positions between it and the locus (default 10, so
   spans that touch always chain), and widen the locus to cover its span. Repeat until no record
   is added. Chaining is transitive: if A chains with B and B with C, all three are in the locus,
   whichever files they come from.
3. Widen the locus to cover every STR region a VCF record in it **touches** — overlaps or lies
   immediately next to. Go back to step 2, since the wider locus may now overlap more records. Stop
   when neither step adds anything.
4. If the locus has at most the chain distance of reference positions between it and the locus
   before it on the contig, the two are one locus, covering both. This can happen only when an STR
   region has widened the locus to the left.

A locus is a run of reference positions, written 1-based and inclusive, as in VCF.

With a SNP at 100, a satellite at 105–3000 that it does not touch, and a record at 2990 inside the
satellite, the record makes a locus 105–3000, and step 4 joins it to the locus at 100: one locus
100–3000. The distance is measured between the loci, not to their records: if the SNP at 100
touches a satellite at 101–200, its locus is 100–200, which the locus 205–3000 joins although the
SNP itself is 104 positions away (decided by the owner on 19 October 2026, over keeping such loci
apart, which would lose the chain distance next to STR regions; the look-ahead then holds a locus
until the reading has passed every STR region that starts within the chain distance after it). The
loci's own spec, `doc/specs/loci.md`, gives the rules in full, with their cases.

**Why a chain distance** (decided by the owner on 19 October 2026). Two events can be traded
between spellings: one file deletes a base and changes the next one, the other changes the first
and deletes the second. Each event's spellings overlap, but the two files' records need not:

```
reference 5–10  AGAACA
file A   5  AG  A   1/1                    deletes the G at 6
file A   9  C   A   1/1
file B   6  G   A   1/1
file B   8  AC  A   1/1                    deletes the C at 9
```

Both files give `AAAAA`. The spans are 5–6 and 9 in A, and 6 and 8–9 in B. By overlap alone they
make two loci, 5–6 and 8–9, each *different*. With one position between them, a chain distance of
1 or more makes one locus, *same genotype*. No finite distance catches every such trade, since
inside a repeat the events can lie any distance apart; 10 catches those within a few bases, at the
cost of wider loci and, in dense unphased data, more blocks per locus (§5.4). The comparison with
other evaluators (§8) shows what is left.

How far 10 chains was simulated, with variant sites of both files together placed at random over
2 Mb and 30% of them unphased heterozygous in each sample:

| sites | longest locus | sample-loci over the phasing cap |
|---|---|---|
| 1 per 30 bases | 56 bases | 0 |
| 1 per 10 bases | 145 bases | 0 |
| 1 per 5 bases | 402 bases | 0.9% |
| 1 per 3 bases | 1,235 bases | 32% |

A long locus is one row, and all its events must agree for *same genotype*. The owner kept 10 on
19 October 2026 with this table in front of them, as real data seldom reaches the dense end.

**Why STR regions widen loci too.** The span (§5.2) makes the spellings of one indel overlap. A
complex record — one that changes more than one base in an allele — has no span beyond its REF,
and two callers can write the same events in a repeat as complex records at different places in
it. Taking the whole repeat whenever a VCF record touches it puts every spelling into one locus.

### 5.4 Reconciliation

The rules below are given in full, with their cases and the Python interface, in
`doc/specs/reconciliation.md`.

For each locus, each file and each sample, pvc_eval builds the sample's possible genotypes, then
expresses the possible genotypes of both files over one shared list of locus alleles.

**Building a locus allele.** A locus allele is built for one haplotype from the reference over the
locus and the VCF alleles that haplotype carries, each where its record stands after normalising
(§5.2): trimmed, never moved, which is how `VcfRecord` holds it. "As written" below means that.
So `p CAA → CA` is normalised to `p CA → C` and deletes position p+1, whichever A its caller
meant; the two are the same haplotype unless another edit lies at p+1 or p+2.

- Each VCF allele of the normalised record other than REF is an **edit**: after trimming the bases
  it still shares with REF, first at the left end and then at the right, it replaces some reference positions (none, for
  an insertion) with some bases (none, for a deletion). Trimming first matters: a VCF record widened
  over flank bases for one allele's sake spells the reference there in its other alleles, and that
  is no edit. The order matters too: `52 AG → A,AAG` trimmed from the left inserts its A after
  position 52, and from the right before it, which differ when the haplotype changes position 52.
  Left first puts the insertion after the anchor, as VCF does.
- A `*` allele makes no edit: it says the haplotype's bases here are removed by another VCF record,
  which carries that edit. The haplotype must carry a deletion, from another VCF record of the
  locus, that deletes a position of the `*` record's REF as the file wrote it, before normalising,
  where the deletion was written or at any other placement its span allows (§5.2), or another edit
  that replaces such a position as written, as GATK writes `*` under an MNP; a haplotype missing
  at a VCF record might carry any of its ALTs. If it carries none, the haplotype is **impossible**.
  The placements count because `bcftools norm` moves the deletion and leaves the `*` record where
  it was. A deletion's anchor base is in its span but deleted by no placement, so a `*` on it
  alone is not covered. The REF counts as written because a `*` names the bases of the REF its
  caller wrote: pvc pads a record with a `*` over bases its alleles share, and normalising would
  trim away the very bases the deletion reaches (decided by the owner on 21 October 2026, over
  the normalised REF; `doc/specs/hard_cases.md` §6, cases S11 and E5).
- The edits of one haplotype are applied to the reference over the locus in position order. Two
  edits **collide** if they change the same reference position, if one is an insertion between
  two positions the other replaces, or if they are two insertions at the same point, as written. Colliding edits cannot both apply: the haplotype is **impossible**.
- If any VCF genotype leaves the haplotype missing (`.`), the locus allele is missing, however many
  other VCF records give it a definite allele. Impossibility is checked first: a haplotype whose
  definite edits collide, or whose `*` is not covered, is impossible even when it is also
  missing.

**Possible genotypes.** A sample's VCF genotypes across the locus's VCF records say which VCF allele
each haplotype carries, up to phase.

- VCF genotypes written with `|` and the same `PS` (or no `PS`, which VCF reads as one phase set per
  sample and contig) fix which VCF alleles go together on one haplotype.
- Between unphased VCF records, or different phase sets, every assignment of VCF alleles to
  haplotypes is tried. A heterozygous unphased VCF record, or a phase set holding a heterozygous
  VCF genotype, is one **block**; with ploidy 2 and *b* blocks there are 2^max(b−1, 0)
  assignments. A VCF genotype is heterozygous here when its two alleles differ, `.` counting as an
  allele, so `./1` is a block. Homozygous VCF genotypes, and records of one phase set, add none.

**The phasing cap.** A sample whose assignments at a locus would exceed the phasing cap (default
4096, that is 13 blocks; kept by the owner on 19 October 2026, to be lowered if it costs too much
time) is not reconciled there, and its sample pairs get the outcome
*too complex* (§5.5). The count is made before any assignment is built, and a sample over the cap
adds nothing to the shared locus alleles. This bounds the work,
which doubles with every block, and only the work: a satellite holding hundreds of phased VCF
records is reconciled, and a locus of a few bases with 14 unphased heterozygous records is not.

Each assignment gives one locus genotype, unless one of its haplotypes is impossible. The sample's
possible genotypes are the distinct locus genotypes that result. A sample whose VCF genotypes allow
none is **invalid** at the locus.

**A sample is also invalid** at a locus if any of its VCF genotypes there is malformed — an allele
that is neither `.` nor an index of the VCF record's alleles. A VCF genotype with a number of alleles other than its file's
ploidy, such as `1` on chrX of a male sample in a file read with ploidy 2, is not malformed but
makes the sample invalid too, with the reason `ploidy mismatch` (decided by the owner on 19
October 2026). A VCF genotype whose alleles are all `.` is missing, whatever their number.

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
allele appearing in any possible genotype of any paired sample of either file, in order of first
appearance, samples in header order, file A before file B. A locus allele whose sequence equals REF, because its
edits cancel or an ALT spells the reference, is REF, index 0. Each possible genotype becomes a
sorted tuple of indices into that list, with missing haplotypes written as missing and sorted after
every index, and a sample's possible genotypes are in the order of those tuples. This is what lets §5.5 compare genotypes by number: two files that spelled one
sequence differently now name it by the same index.

### 5.5 Evaluation

The rules below are given in full, with their cases, in `doc/specs/reconciliation.md`.

Each sample pair at each compared locus gets exactly one outcome. The first rule that applies wins:

| outcome | when |
|---|---|
| **too complex** | either file's sample exceeds the phasing cap at the locus and was not reconciled (§5.4) |
| **invalid in A**, **invalid in B**, **invalid in both** | that file's sample is invalid at the locus |
| **missing in A**, **missing in B**, **missing in both** | every possible genotype of that file's sample has a missing haplotype |
| **reference in both** | in both files, the only possible genotype is all REF |
| **same genotype** | some possible genotype of A equals some possible genotype of B; two genotypes all REF are equal whatever their ploidies |
| **same alleles, other genotype** | some possible genotype of A and some of B carry the same set of non-reference locus alleles, on different numbers of haplotypes (one in A, both in B, say) |
| **only in A**, **only in B** | one file gives the sample a non-reference locus allele, and the other's only possible genotype is all REF |
| **different** | both give the sample non-reference locus alleles, and no pair of possible genotypes matches as above |

When several possible genotypes exist, the pair reported in the table is the first, taking A's
possible genotypes in order and, for each, B's in order, that meets the rule of the outcome. Where
the rule names no pair (*missing*, *reference in both*, *only in*, *different*), it is the first
possible genotype of each file. *Too complex* and *invalid* report no pair; in their rows a sample
that was reconciled shows its first possible genotype (§6). Rows come locus by
locus, and within a locus in the order of the sample pairs: file A's header order, or the order of
`--sample-map`.

`./1` in both files is *missing in both*, not *same genotype*: half of each genotype is unknown,
so nothing says the files agree there (decided by the owner on 19 October 2026).

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
| **SNP** | every ALT of every record, trimmed against REF as an edit is (§5.4), replaces one base with one base |
| **indel** | every ALT of every record, trimmed so, is one insertion or one deletion |
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
| `str_regions` | the STR regions in the locus, 1-based and inclusive, as `start-end:type`, with `:motif` after an STR segment's and `start-end` alone for a region with no type, joined with `,`; `.` for none |
| `sample_a`, `sample_b` | the sample pair |
| `outcome` | §5.5 |
| `alleles` | the locus alleles, REF first, joined with `,` (a deletion of the whole locus is written `-`) |
| `genotype_a`, `genotype_b` | the reported possible genotypes, as indices into `alleles` joined with `/`, `.` for a missing haplotype; `.` for a *too complex* or invalid sample |
| `possible_a`, `possible_b` | how many possible genotypes the sample had in each file; `.` for a *too complex* sample, 0 for an invalid one |
| `quality_a`, `quality_b` | §2, *quality*; `.` when the file gives the sample no non-reference VCF allele |
| `gq_a`, `gq_b` | the lowest GQ among the sample's VCF records in the locus that give it a non-reference VCF allele, `.` when absent |
| `records_a`, `records_b` | how many VCF records with an ALT of each file the locus holds |
| `invalid_reason_a`, `invalid_reason_b` | for an invalid sample, whatever the outcome, the first that applies: `malformed GT`, `ploidy mismatch`, `contradicting records`; `.` otherwise |

The summary (`--summary`) counts the rows, *reference in both* included, summed over all sample
pairs: one line per outcome and stratum, one column per kind and a total (decided by the owner on
20 October 2026, to be tried first). `doc/specs/compare.md` gives how each field of a row and of
the summary is written.

## 7. Tests

Nothing about pvc_eval is trusted until these pass.

1. **The STR finder against pvc:** the port's own copy of pvc's tests (§5.1), and identical STR
   regions, types included, from the port and from `pvc dev type-regions` over the GIAB benchmark
   regions of GRCh38 (with pvc's repeat catalogue, §5.1).
2. **Each stage on hand-written cases**, including:
   - normalising and spans: trimming at both ends; a `*`-only record, not trimmed; deletion and
     insertion spans in homopolymers and longer repeats; a multi-ALT record; a complex record; the
     three cases of §5.2; a span cut at 10 kb; a REF that differs from the reference, which stops
     the run, and a lower-case reference, which is no difference;
   - loci: chaining across both files; the traded events of §5.3, one locus at chain distance 1
     and two at 0; adjacent SNPs, chained at 0; a record with no ALT, which starts no locus; an insertion next to a repeat pulling it in; two
     complex records at the two ends of a repeat; with `--regions`, a locus
     sticking out of a region compared whole, and one sharing no position with any region skipped;
   - reconciliation: the invalid example of §5.4 and its valid variant; GATK's `*` form, and a `*`
     on a haplotype with no deletion, and one on a deletion's anchor base only (impossible); the
     `bcftools norm` output of GATK's `*` form (valid), and of the first case of §5.2 (invalid); an insertion that slides with its bases rotated; the insertion of `52 AG → A,AAG` after the
     anchor; phase sets; the phasing cap at 13 blocks (reconciled) and 14 (*too complex*), and 14
     heterozygous records in one phase set (one block, reconciled); `./1`; malformed GT; ploidy
     mismatch, and a file of ploidy 1;
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
   indels shifted right within their repeats, where no other edit of the same haplotype lies in
   the stretch they cross; GATK-style `*` records added — every row in the set allowed for
   the file against itself. Records decomposed from one record are written phased, in one phase
   set, so the re-spelling adds no blocks (§5.4). Each re-spelling rewrites only records none of
   whose edits on a haplotype touches an edit of another record on the same haplotype, edits
   being judged as §5.4 makes them: collisions are judged on
   edits as written (§5.4), so where the edits of one haplotype touch, the outcome can follow the
   spelling (decided by the owner on 19 October 2026, over collision rules that do not depend on
   it). The re-spelling program is part of the tests; its spec is `doc/specs/respell.md`.
   `bcftools norm` is not one of the re-spellings, since it can change a haplotype (§5.2).
5. **The hard cases from pvc's STR work:** 50 cases, each the VCF records of two real files at one
   locus from pvc's reports, with a few hundred bases of the reference and the row `compare()` must
   give: GIAB's split records against pvc's one repeat record, deletions at a repeat's edge written
   with overlapping records, HipSTR's whole-repeat records, and other spellings. They are given,
   with their outcomes, in `doc/specs/hard_cases.md`, approved by the owner on 21 October 2026.

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
- **Ploidy other than 1 or 2**: functions are written for any ploidy where that costs nothing; where
  one is limited to 1 and 2 (the phase enumeration of §5.4), it raises an error for any other
  ploidy.

## 10. Open points

Open 1, haploid genotypes, and Open 2, the phasing cap's default, were answered by the owner on
19 October 2026: each file is read with one ploidy, 1 or 2, and a GT of another ploidy makes the
sample invalid (§5.4); the cap stays at 4096, to be lowered if it costs too much time.

Open 3, the hard cases of §7.5, was answered by the owner on 21 October 2026: the 50 cases of
`doc/specs/hard_cases.md`, each with the outcome it must give.
