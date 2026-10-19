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
| **normalise** | trim the bases a VCF record's alleles share at either end, then left-align an indel (§5.2). Never split or decompose a VCF record. | "decompose", "atomise" |
| **VCF record span** | the reference positions a VCF record's REF covers, after normalising | |
| **sample pair** | a sample of file A matched to a sample of file B: same name unless the user maps them | |
| **ploidy** | the number of haplotypes a sample has at every locus; a global constant, 2 | |

**From the reference**

| term | meaning | not |
|---|---|---|
| **segmented genome** | each contig of the reference cut by the STR finder (§5.1) into consecutive, non-overlapping regions, each a generic region or an STR region | |
| **generic region** | a stretch of the segmented genome that holds no tandem repeat the STR finder admits; pvc's `generic` | "ordinary sequence" |
| **STR region** | a stretch the STR finder (§5.1) reports as a tandem repeat: an **STR segment** (one repeat), an **STR bundle** (two or more repeats within 15 bases of each other) or a **satellite** (a repeat run longer than 100 bases). pvc's terms. | "tract", "STR locus", "repeat locus" |

**Building what is compared**

| term | meaning | not |
|---|---|---|
| **locus** | the unit compared: a stretch of reference holding VCF records of either file, chained by overlap and widened to every STR region a VCF record touches (§5.3) | any other use of "locus" |
| **touches** | overlaps, or lies immediately next to. A left-aligned insertion of a repeat unit sits on the base just before the repeat; without "next to" it would never pull the repeat in. | |
| **record cap** | the most VCF records, of both files together, a locus may hold | |
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
| **stratum** | the class `--strata` gives the locus: from a stratification BED (homopolymer, segmental duplication, …), or, with `--strata pvc`, the kind of STR region the locus holds (§5.7) |
| **quality** | for one file and sample, the lowest QUAL among the VCF records giving the sample a non-reference VCF allele at the locus |

## 3. Interface

```
pvc_eval compare --reference ref.fa --a a.vcf.gz --b b.vcf.gz --output loci.tsv
                 [--regions regions.bed] [--strata strata.bed|pvc] [--str-regions str.bed]
                 [--sample-map pairs.tsv] [--record-cap N] [--summary summary.tsv] [--all-rows]
                 [--evaluator pvc|vcfeval|happy|naive]
```

| option | what it does |
|---|---|
| `--reference` | FASTA with its `.fai`. Contig names and order come from the `.fai`. |
| `--a`, `--b` | the two VCFs, plain or bgzipped. |
| `--regions` | compare only the loci that overlap these regions, by at least one position; every other locus is skipped and gets no row, nor a count in the summary. A locus that overlaps is compared whole, its VCF records outside the regions included, since loci are built from all VCF records (§5.3). Without it, every locus is compared. |
| `--strata` | BED with a fourth column naming each interval's class; a locus takes the class of the intervals it overlaps (several are joined with `,`; none gives `.`). The word `pvc` instead of a path gives each locus one class from the STR regions it holds — `satellite`, `str-bundle`, `str` or `generic` — with no BED needed (§5.7); a file named `pvc` is given as `./pvc`. |
| `--str-regions` | BED of STR regions to use instead of running the STR finder. |
| `--sample-map` | two columns: a sample of A, a sample of B. Without it, samples are paired by name. |
| `--record-cap` | default 30 (§5.3). |
| `--output` | the table of §6, one row per sample pair and locus. |
| `--summary` | counts of outcomes by kind and stratum. |
| `--all-rows` | also write the rows whose outcome is *reference in both* (§5.5). |
| `--evaluator` | what decides each locus's outcome: `pvc`, pvc_eval's own rules (§5.4, §5.5), the default; or `vcfeval`, `happy` or `naive`, another evaluator wrapped behind the same interface (§8). |

**Errors, which stop the run:** a VCF not sorted in the `.fai`'s contig order and by position; a
contig in a VCF but not in the reference; no sample pairs; a sample named in `--sample-map` absent
from its file. Samples present in only one file and not mapped are listed as a warning and not
compared.

**Not errors:** everything about the content of a VCF record — a REF that disagrees with the
reference, a malformed GT, contradicting VCF records. Those make the affected locus *invalid* for the
affected sample (§5.4) and the run goes on.

The Python API offers each stage of §4 as a function over iterators, so a user can, for example,
filter VCF records between reading and grouping, or read loci and genotypes without evaluating them.
Above the stages, one function does what `compare` does, whichever the evaluator: it takes the
two VCFs, the reference and the options of the table above, and yields the rows of §6 in position
order, one per sample pair and locus.

```python
def compare(a: Path, b: Path, reference: Path, config: CompareConfig,
            evaluator: Evaluator = PVC) -> Iterator[OutcomeRow]: ...
```

`CompareConfig` holds the options of the table above but `--output`, `--summary` and
`--evaluator`.

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
   reference, normalised, and kept in position order.
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
2. **Admit or reject each candidate as one STR** (pvc `repeat_catalog/row.rs`, `segment_criteria.rs`).
   A candidate is rejected — and its bases stay in a generic region — if it has fewer copies than the
   floor for its period, its motif is compound, it has no clean whole-motif cut at both ends, its
   purity after the cut is below the minimum, or it lies within 15 bases of a contig end. Candidates
   chained within 15 bases of each other become an **STR bundle**; a lone admitted candidate is an
   **STR segment**, trimmed to whole motifs.
3. **Satellites** (pvc `region_typing::resolve_features`). A run of repeat longer than 100 bases is a
   satellite, and swallows every STR segment or bundle member within 15 bases of it, growing to cover
   them. A bundle whose span exceeds 100 bases becomes a satellite. A bundle touching a contig end is
   dropped.

The result is the segmented genome: pvc's partition of the contig into STR segments, STR bundles,
satellites and generic regions. Loci are built from its STR regions (§5.3), and `--strata pvc`
reads its kinds (§5.7).

**Satellites are STR regions in pvc_eval.** The vocabulary of §2 names two kinds of STR region; this
section adds a third. An indel in a satellite can be written in as many places as in a short repeat,
so a VCF record touching a satellite must pull it in like any other repeat. A satellite can be
kilobases long; the record cap (§5.3) keeps such loci bounded.

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
5 copies already raises. The lower floor only makes some loci wider, which costs nothing in
correctness. Every parameter can be set on the command line; the check against pvc (below) is run
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
  and writes STR regions as a BED, for `compare --str-regions` to reuse; its fourth column is each
  region's kind, `str`, `str-bundle` or `satellite`, which `--strata pvc` reads (§5.7). A
  `--str-regions` BED without that column can be used for loci, and with `--strata pvc` it is an
  error. pvc gives the same answer whether or not a region list is given (its `typed_regions.md`
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
be identical. That check needs a pvc binary; it is skipped, with a message, where none is found,
and pvc_eval never needs pvc to run.

### 5.2 Reading and normalising VCF records

**Checks.** Each VCF record's contig must be in the reference, and records must come in the `.fai`'s
contig order and by position; otherwise the run stops with an error naming the file and line. A
VCF record whose REF differs from the reference at its position is kept and marked *REF mismatch*;
every locus containing it is invalid for every sample of that file (§5.4).

**Normalising.** For each VCF record:

1. **Trim.** Remove the bases all its VCF alleles share at the right end, then at the left end,
   keeping at least one base in every allele (the VCF anchor convention). The VCF record span shrinks
   to what remains.
2. **Left-align.** If the record is an indel — one VCF allele empty after removing the anchor base —
   slide it left while the base it would uncover equals the base it would drop, as `bcftools norm`
   and vt do. A VCF record with several ALTs is left-aligned only if every ALT moves by the same
   amount; otherwise it stays where it is. A VCF record that changes more than one base in an allele
   (a complex record) has no left-aligned form and stays where its caller put it; §5.3 is what makes
   that harmless inside repeats.

**Order after normalising.** Left-alignment can move a VCF record before ones already read. The
stream holds VCF records until the reader is further right than the longest distance any record has
moved, and emits them in order. A record that would move more than 10 kb is not moved and is
reported in a warning.

### 5.3 Loci

**Building a locus.** Starting from the leftmost VCF record not yet in a locus, of either file:

1. The locus is that record's span.
2. Add every VCF record, of either file, whose span overlaps the locus, and widen the locus to cover
   it. Repeat until no more records overlap. Overlap chains: if A overlaps B and B overlaps C, all
   three are in the locus, whichever files they come from.
3. Widen the locus to cover every STR region a VCF record in it **touches** — overlaps or lies
   immediately next to. Go back to step 2, since the wider locus may now overlap more records. Stop
   when neither step adds anything.

A locus is a run of reference positions, written 1-based and inclusive, as in VCF.

**Why STR regions widen loci and nothing else does.** After normalising, two spellings of one indel
start at the same position, so they overlap. Two spellings of several events do not always: a caller
that writes "delete one AT, then change the G after the repeat" as one complex VCF record leaves it
at the end of the repeat, while another caller's separate deletion left-aligns to the start of the
repeat. Taking the whole repeat whenever a VCF record touches it puts both into one locus.

**The record cap.** A locus holding more VCF records than the record cap (default 30, both files
together) is not reconciled; every sample pair gets the outcome *too many records*. This bounds the
work in §5.4, whose number of possible genotypes can double with every unphased heterozygous VCF
record.

### 5.4 Reconciliation

For each locus, each file and each sample, pvc_eval builds the sample's possible genotypes, then
expresses the possible genotypes of both files over one shared list of locus alleles.

**Building a locus allele.** A locus allele is built for one haplotype from the reference over the
locus and the VCF alleles that haplotype carries:

- Each VCF allele other than REF is an **edit**: after trimming the bases it shares with its record's
  REF at either end, it replaces some reference positions (none, for an insertion) with some bases
  (none, for a deletion). Trimming first matters: a VCF record widened over flank bases for one
  allele's sake spells the reference there in its other alleles, and that is no edit.
- A `*` allele makes no edit: it says the haplotype's bases here are removed by another VCF record,
  which carries that edit.
- The edits of one haplotype are applied to the reference over the locus in position order. Two
  edits that change the same reference position, or two insertions at the same point, cannot both
  apply: the haplotype is **impossible**.
- If any VCF genotype leaves the haplotype missing (`.`), the locus allele is missing, however many
  other VCF records give it a definite allele.

**Possible genotypes.** A sample's VCF genotypes across the locus's VCF records say which VCF allele
each haplotype carries, up to phase:

- VCF genotypes written with `|` and the same `PS` (or no `PS`, which VCF reads as one phase set per
  sample and contig) fix which VCF alleles go together on one haplotype.
- Between unphased VCF records, or different phase sets, every assignment of VCF alleles to
  haplotypes is tried. With ploidy 2 and *h* heterozygous unphased VCF records, that is 2^(h−1)
  assignments.

Each assignment gives one locus genotype, unless one of its haplotypes is impossible. The sample's
possible genotypes are the distinct locus genotypes that result. A sample whose VCF genotypes allow
none is **invalid** at the locus.

**A sample is also invalid** at a locus if any of its VCF genotypes there is malformed — a number of
alleles other than the ploidy, or an allele index the VCF record lacks — or if any VCF record of its
file in the locus has a REF mismatch (§5.2).

**The example of an invalid sample.** Reference positions 9 to 13 read `CAGTC`. One file has

```
POS REF  ALT  GT
9   CAGT C    1/1     both haplotypes lose AGT (positions 10–12)
11  G    A    0/1     one haplotype has A at position 11
```

Position 11 no longer exists on either haplotype, yet one is said to carry an A there: no assignment
works, so the sample is invalid at this locus. With the deletion `0/1` the same records are fine —
one haplotype `CC`, the other `CAATC` — and GATK writes the second record as `G → A,*` with `1/2`,
which says exactly that.

**Shared locus alleles.** The locus's alleles are REF followed by every distinct non-missing locus
allele appearing in any possible genotype of any sample of either file, in order of first appearance.
Each possible genotype becomes a sorted tuple of indices into that list, with missing haplotypes
written as missing. This is what lets §5.5 compare genotypes by number: two files that spelled one
sequence differently now name it by the same index.

### 5.5 Evaluation

Each sample pair at each compared locus gets exactly one outcome. The first rule that applies wins:

| outcome | when |
|---|---|
| **too many records** | the locus holds more VCF records than the record cap and was not reconciled (§5.3) |
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

The kind of a locus is read from its locus alleles:

| kind | when |
|---|---|
| **STR** | the locus contains an STR region |
| **SNP** | every non-reference locus allele has the length of REF and differs from it at one position |
| **indel** | every non-reference locus allele differs from REF by one insertion or one deletion, and nothing else |
| **complex** | anything else: several changed positions, or a change of length together with a substitution |

### 5.7 Strata from pvc's region kinds

`--strata pvc` sorts loci by the part of pvc's partition (§5.1) they fall in, which is the
stratification most comparisons of pvc want. The four classes are pvc's region kinds, with
pvc_eval's names: `generic` (a generic region; pvc `RegionKind::Generic`), `str` (an STR segment;
`ssr_locus` in pvc's typed-regions BED), `str-bundle` (`ssr_bundle`) and `satellite`.

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
`generic,str`: a locus widened to an STR region (§5.3) usually keeps the base before it, where a
left-aligned insertion sits.

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
| `str_regions` | the STR regions in the locus, as `start-end:motif` joined with `,`; `.` for none |
| `sample_a`, `sample_b` | the sample pair |
| `outcome` | §5.5 |
| `alleles` | the locus alleles, REF first, joined with `,` (a deletion of the whole locus is written `-`) |
| `genotype_a`, `genotype_b` | the reported possible genotypes, as indices into `alleles` joined with `/`, `.` for a missing haplotype |
| `possible_a`, `possible_b` | how many possible genotypes the sample had in each file |
| `quality_a`, `quality_b` | §2, *quality*; `.` when the file gives the sample no non-reference VCF allele |
| `gq_a`, `gq_b` | the lowest GQ among the sample's VCF records in the locus that give it a non-reference VCF allele, `.` when absent |
| `records_a`, `records_b` | how many VCF records of each file the locus holds |
| `invalid_reason_a`, `invalid_reason_b` | for an invalid outcome: `contradicting records`, `REF mismatch` or `malformed GT` |

The summary (`--summary`) counts rows by outcome, kind and stratum.

## 7. Tests

Nothing about pvc_eval is trusted until these pass.

1. **The STR finder against pvc:** the port's own copy of pvc's tests (§5.1), and identical STR
   regions, kinds included, from the port and from `pvc dev type-regions` over the GIAB benchmark
   regions of GRCh38.
2. **Each stage on hand-written cases**, including:
   - normalising: trimming at both ends; left-alignment in homopolymers and longer repeats; a
     multi-ALT record whose ALTs would move by different amounts; a complex record, not moved;
   - loci: chaining across both files; a left-aligned insertion next to a repeat pulling it in; the
     complex-record-plus-separate-records example of §5.3; the record cap; with `--regions`, a locus
     sticking out of a region compared whole, and one sharing no position with any region skipped;
   - reconciliation: the invalid example of §5.4 and its valid variant; GATK's `*` form; phase
     sets; `./1`; a REF mismatch; malformed GT;
   - evaluation: one case per outcome;
   - strata with `--strata pvc`: the insertion of §5.7 (`str`), a locus with no STR region
     (`generic`), the STR segment chained to a satellite (`satellite`) and to a bundle
     (`str-bundle`), a bundle chained to a satellite (`satellite`), and a `--str-regions` BED
     without a kind column (an error).
3. **A file against itself:** every row *same genotype* or *reference in both*.
4. **A file against re-spellings of itself** — every VCF record decomposed into single-base and
   single-indel records with genotypes carried over; neighbouring records merged into complex ones;
   indels shifted right within their repeats; GATK-style `*` records added — every row *same
   genotype* or *reference in both*. The re-spelling program is part of the tests.
5. **The hard cases from pvc's STR work:** the VCF records behind the cases recorded in pvc's
   reports (GIAB's split indel-plus-SNP truth records against pvc's single repeat record; GATK's
   overlapping records at STR edge deletions; HipSTR's whole-repeat records), each with the outcome
   it must give.

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
