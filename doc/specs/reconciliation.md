# Reconciliation and evaluation — specification

*Status: approved by the owner on 19 October 2026. No code of these stages exists yet. It
builds stages 4 and 5 of `doc/spec.md` §4 and stands on `doc/spec.md` §2, §4, §5.2, §5.4, §5.5 and
§6, on the VCF reader (`src/pvc_eval/vcf.py`) and on the loci's spec, `doc/specs/loci.md`, which
gives the kind and the stratum of a locus.*

Reconciliation takes one locus at a time. It rebuilds the sequence each haplotype of each paired
sample carries over the locus, in both files, under every phase the VCF genotypes leave open. It
then names those sequences by one shared list, so that two files that spelled a haplotype
differently give it the same index. Evaluation compares the two samples of each sample pair at the
locus and gives them one outcome of `doc/spec.md` §5.5, with the pair of possible genotypes that
decided it. Together they fill the columns of `doc/spec.md` §6 that depend on genotypes. Every rule below
has a case with its result, and each case is a test.

This spec also gives what goes in and what comes out: the sample pairs, the ploidy each file is
read with, and `OutcomeRow`, one row of the table, which every evaluator yields (`doc/spec.md` §8).

The owner decided, on 19 October 2026:

- `OutcomeRow` follows the Python API, 0-based and with `None` and enums, over mirroring the TSV;
- an insertion collides with an edit that replaces the bases on both sides of it, over letting
  only shared positions collide;
- a `*` is covered by a deletion its haplotype might carry where its genotype is `.`, over making
  such a haplotype impossible;
- each file is read with one ploidy, given by the user, 2 by default, which settles `doc/spec.md`
  Open 1; a GT of another ploidy makes the sample invalid;
- the phasing cap stays at 4096 by default, to be lowered if it costs too much time (`doc/spec.md`
  Open 2);
- a `*` is covered by any edit that replaces a position of its record's REF, as GATK writes `*`
  under an MNP, over deletions only (Open 4 of the draft);
- two genotypes all REF count as equal whatever their ploidies (Open 5 of the draft).

The decisions change `doc/spec.md` beyond its §5.4 and §5.5, as §10 lists. Four smaller points
are decided here, each where it stands: a GT of only `.` is missing whatever its
  number of alleles (§2); a sample in no pair is not reconciled (§4); a reason is written whenever a sample is invalid, also in a row *too complex*
  because of its partner (§6); and `quality` and `gq` count a `*` as an ALT (§6).

Positions in the cases are 1-based and inclusive, as in VCF and in `doc/spec.md`; in the Python
interface they are 0-based and half-open, as everywhere inside pvc_eval.

## 1. Vocabulary

The terms of `doc/spec.md` §2 keep their meaning: locus, locus allele, locus genotype, possible
genotypes, missing, invalid, haplotype, phasing cap. §5.4 of it defines an edit, a block and an
impossible haplotype in bold, and this spec makes each of them exact. Names from the code are those of `doc/spec.md` §4 and of the VCF reader. A `Locus` holds a
locus's records of each file; a `VcfRecord` is one VCF record, normalised, with its span
(`doc/spec.md` §5.2: every position the record could have been written at, cut back to its REF
when it would reach further than the span limit, 10 kb by default); a `VcfGenotype` is one
sample's GT at it. Where a GT cannot be a genotype the reader puts a marker in its place,
`MALFORMED_GT` or, after §2, `OTHER_PLOIDY_GT`, and `MISSING` stands for a `.` haplotype. These markers are
sentinels, objects made by Python 3.15's `sentinel()` that are compared with `is` and can stand in
a type hint (`.claude/skills/coding/python.md`).

It adds:

| term | meaning | not |
|---|---|---|
| **edit** | what one VCF allele of a normalised record does to the reference: replace a stretch of positions (none, for an insertion) with some bases (none, for a deletion) | "variant", "event" |
| **phase unit** | the VCF genotypes of one sample whose order on the haplotypes is fixed together: one unphased VCF genotype, or every phased VCF genotype of one phase set | "phase block" |
| **block** | a phase unit with a heterozygous VCF genotype in it (`doc/spec.md` §5.4) | |
| **phase assignment** | one choice of order for every block of a sample at a locus; it gives each haplotype one VCF allele of each VCF record | "phasing" |
| **placement** | a stretch a deletion could delete instead of the one it was written over and give the same haplotype: one of its slides of `doc/spec.md` §5.2 | "shift" |
| **sample result** | what reconciliation gives one sample of one file at one locus: its possible genotypes, *too complex*, or *invalid* with its reason | |

## 2. Ploidy

Each file is read with one ploidy, 1 or 2, which the user gives for each file and which is 2 when
they do not (decided by the owner on 19 October 2026). A GT whose number of alleles is not that
ploidy, such as `1` on chrX of a male sample in a file read as diploid, is not malformed. It makes
the sample invalid at the locus with its own reason, `ploidy mismatch`, so the table says why
nothing was compared there. A file whose GTs all have one allele, a set of haploid calls, is read
with ploidy 1, and its samples are compared like any other.

So the VCF reader takes the ploidy:

- `parse_vcf(path, reference, *, span_limit=..., ploidy=DEFAULT_PLOIDY)`. `ploidy` is a whole
  number of 1 or more, and anything else raises `ValueError` when the stream starts, as
  `span_limit` does.
- A GT is `MALFORMED_GT` when one of its alleles is neither `.` nor an index of one of the
  record's VCF alleles.
- A GT whose alleles are all `.` or such indices, and whose number of alleles differs from the
  ploidy, is `OTHER_PLOIDY_GT`, which replaces `HAPLOID_GT`, unless every allele is `.`.
- A GT whose alleles are all `.`, a single `.` or `./.` in a file of either ploidy, is `MISSING`
  on every haplotype, ploidy times (decided here: it says nothing about any haplotype, and
  callers write `./.` for a no-call in haploid files too).
- A sample with no GT in a record's FORMAT is `MISSING` on every haplotype, ploidy times, as now.
- `PLOIDY` becomes `DEFAULT_PLOIDY = 2`.

Reconciliation handles ploidy 1 and 2. A ploidy of 1 has no phase to choose: each sample has one
phase assignment. A higher ploidy raises `ValueError` in `reconcile`, as `doc/spec.md` §9 foresees
for the phase enumeration, and the command line refuses `--ploidy-a` or `--ploidy-b` other than 1
or 2 before it reads anything (the spec of `compare`).

Two files read with different ploidies are compared by the rules of §5, where a genotype of one
haplotype never equals one of two, with one rule of their own (decided by the owner on 19 October 2026): two genotypes
all REF count as equal whatever their ploidies. A haploid `1` against a diploid `1/1` is *same
alleles, other genotype* (E8 of §5); a haploid `0` against `0/0` is *reference in both*.

## 3. Sample pairs

A sample pair names a sample of file A and a sample of file B by their names and by their places
among the sample columns of their files' headers, which is where each VCF record holds their
genotypes.

- **Without `--sample-map`**, every name in both headers is a pair, in file A's header order.
- **With `--sample-map`**, each line of the map is a pair, in the order of its lines. A line has
  two fields separated by a tab, a sample of A and a sample of B, and no header. A line of nothing
  but spaces and tabs is blank and skipped. The map is read as UTF-8, and a byte-order mark at its
  start, which spreadsheet programs write, is ignored. Names are compared exactly, spaces
  included. A sample may be in several pairs, as one sample of A compared with two replicates in B.

These stop the run with a `PvcEvalInputError`, whose message names the file and the line:

- a map that cannot be read: "pairs.tsv does not exist. Check its path and run again.", and
  likewise for a directory or a file that cannot be opened; a map that is not UTF-8: "pairs.tsv
  is not UTF-8 text: line 3 cannot be read as UTF-8. Save it as UTF-8 and run again.";

- a line of the map that is not two non-empty fields: "pairs.tsv, line 3: a line holds a sample
  of a.vcf.gz and a sample of b.vcf.gz, separated by a tab, and this one has 1 field. Correct the
  line and run again.", or, for two fields with one empty, "and this one has an empty field";
- a sample of the map absent from its file, the sample of A first when neither is in its file,
  with the name in quotes so that a space in it shows: "pairs.tsv, line 3: 'HG005' is not a sample of
  b.vcf.gz, whose samples are HG002, HG003, HG004. Correct the name and run again." The list
  names the first 10 samples and then how many more there are; a file with no sample columns
  gives "'HG005' is not a sample of b.vcf.gz, which has no samples.";
- the same pair on two lines: "pairs.tsv, lines 2 and 5: both pair HG002 of a.vcf.gz with HG002
  of b.vcf.gz, which would give every row twice. Remove one and run again.";
- a name that appears twice in the header of a file, when a pair needs it, since its genotypes
  could be either column (VCF 4.3 §1.5 does not allow duplicate sample names): "a.vcf.gz names the
  sample HG002 in two columns, 10 and 12. Rename one of them and run again." The columns are
  those of the VCF line, the first sample in column 10, and a name in more than two columns has
  them all: "in 3 columns, 10, 12 and 13";
- no pair at all: "a.vcf.gz and b.vcf.gz share no sample name. Pair their samples with
  --sample-map and run again." With a map, "pairs.tsv holds no pair. Add a line for each
  pair, a sample of a.vcf.gz and a sample of b.vcf.gz separated by a tab, and run again."

The samples of each file that are in no pair are given in one warning per file and are not
compared: "2 samples of a.vcf.gz are in no sample pair and are not compared: HG003, HG004. Pair
them with --sample-map to compare them." With a map, the advice is to add them to it. The warning
names each sample once, a name that fills two columns too, the first 10 of them, and then says
how many more there are.

## 4. Reconciliation

For each locus, reconciliation works on each sample of each file that is in a sample pair, file
A's in its header order and then file B's, each sample once however many pairs it is in. A sample
in no pair is not reconciled and adds nothing to the locus alleles (decided here, since its rows
are never written). The sample is given one sample result by the first of these that applies:

1. **Invalid, `malformed GT`**, when any of its VCF genotypes in the locus is `MALFORMED_GT`.
2. **Invalid, `ploidy mismatch`**, when any is `OTHER_PLOIDY_GT`.
3. **Too complex**, when its phase assignments would exceed the phasing cap (§4.1).
4. **Invalid, `contradicting records`**, when no phase assignment gives a locus genotype (§4.2).
5. Otherwise, **its possible genotypes**, one or more.

A REF that differs from the reference stops the run when the VCF is read (`doc/spec.md` §5.2,
decided by the owner on 19 October 2026, over making its loci invalid), so no such record
reaches reconciliation.

The order is that of the invalid reasons in `doc/spec.md` §6, with the cap placed after the reasons that need no
phase assignment. A malformed GT has no alleles, so its blocks could not be counted. And a sample
whose GTs are wrong is invalid however many blocks it has.

### 4.1 Phase assignments and the cap

A sample's VCF genotypes at the locus's records of its file fall into phase units. `Locus` holds
each file's records by span start, and in the file's order where two spans start together
(`doc/spec.md` §4); that order numbers the blocks below.

- A VCF genotype written with `|` belongs to the phase unit of its phase set: its `PS`, or, with
  no `PS`, the one phase set VCF gives every phased genotype of the sample that has none (VCF 4.3
  §1.6.2, `PS`).
- A VCF genotype written with `/`, or a single `.`, is a phase unit of its own; its `PS`, if any,
  is ignored, as VCF says.

A phase unit is a block when one of its VCF genotypes has two different alleles, `.` counting as
an allele, so `./1` makes a block and `./.` does not. With ploidy 2 and *b* blocks there are
2^max(*b*−1, 0) phase assignments: the first block keeps the order written, since a locus genotype
is unordered, and each other block keeps it or swaps its two haplotypes. With ploidy 1 there is
one.

The count is made from the blocks, before any locus allele is built. A sample whose count exceeds
the phasing cap is *too complex*: 13 blocks give 4096 assignments, within the default cap, and 14
give 8192, over it. Fourteen heterozygous VCF genotypes phased in one phase set are one block and
one assignment. The cap is a whole number of 1 or more, `PHASING_CAP = 4096` by default.

The assignments are taken in a fixed order, because the order of the locus alleles (§4.3) follows
it. The blocks are numbered from 0 in the order of their first VCF record in the locus, and
assignment *k*, for *k* from 0 to 2^(*b*−1)−1, swaps block *j* (for *j* ≥ 1) when bit *j*−1 of *k*
is set. So assignment 0 is every VCF genotype as written, and assignment 1 swaps the second block
alone.

### 4.2 Building a locus allele

An assignment gives each haplotype one VCF allele of each record of the file in the locus: an
index into the record's alleles, or `MISSING`. The haplotype's locus allele is built from those,
on the reference over the locus, in three steps.

**The edits.** Each VCF allele the haplotype carries that is neither REF (index 0), `MISSING`
nor `*` gives one edit. It is found from the record as `VcfRecord` holds it, normalised: the bases
the allele shares with REF are removed first from the left and then from the right, and what is
left of REF is the stretch of positions replaced, and what is left of the allele the bases put in
their place. Left first puts an insertion after the anchor base, as VCF writes it. An allele that
leaves nothing of either, one that spells REF, makes no edit.

```
51–54 TAGC   52  AG  A,AAG    A:   deletes G at 53
                              AAG: inserts A after 52 (from the right first, it would go before 52)
```

**Impossible.** The haplotype is impossible when either of these holds:

- Two of its edits **collide** as written (`doc/spec.md` §5.2 says why they are never moved).
  Two edits that replace positions collide when they share one. An insertion collides with an edit
  that replaces both the position before the insertion and the one after it. Two insertions collide
  when they go in at the same point. An insertion at either end of a replaced stretch, or a
  replacement next to another, does not collide.
- It carries a `*` that no deletion **covers**. A deletion covers a `*` when one of its placements
  deletes a position of the `*` record's REF; any other edit covers it when, as written, it
  replaces such a position, as GATK writes `*` under an MNP (decided by the owner on 19 October
  2026). An edit of a record whose allele the haplotype carries counts. So does an edit among the
  alleles of a record where the haplotype is `MISSING`, since it might carry that one. The placements of a deletion are where it was
  written and every slide of `doc/spec.md` §5.2 that deletes only positions of its record's span:
  a record whose span the span limit cut back to its REF allows only the slides within its REF. A deletion is an
  edit that replaces one or more positions with no bases. Where the span is not cut, the base
  before the deletion's left-aligned form, its anchor in `bcftools norm`'s spelling, is deleted
  by none of its placements; the base before it as written may be, when it slides left.

**Missing.** A haplotype that is not impossible and carries `MISSING` at any record has a missing
locus allele, however many other records give it a definite allele. Impossibility comes first: a
haplotype missing at one record and with colliding edits at others is impossible.

**The sequence.** Otherwise the edits are applied to the reference over the locus in position
order, an insertion before a replacement that starts at its point. Every edit lies within the
locus: the locus covers the span, and so the REF, of every record with an ALT (`doc/spec.md` §4), and a record with
no ALT makes no edit, however far its REF reaches.

An assignment in which any haplotype is impossible gives no locus genotype. The others each give
one, the multiset of their haplotypes' locus alleles.

### 4.3 The shared locus alleles and the possible genotypes

The locus alleles are the reference over the locus, at index 0, followed by every distinct locus
allele, not missing, of every locus genotype of §4.2, in order of first appearance. They are taken
sample by sample in the order of §4, and within a sample assignment by assignment in the order of
§4.1, haplotype by haplotype. A locus allele whose sequence is the reference's is index 0. Samples
that are invalid or too complex add none. A locus allele can be empty, a deletion of the whole
locus.

Each locus genotype becomes the tuple of its locus alleles' indices, sorted, with `MISSING` after
every index. A sample's possible genotypes are the distinct tuples, sorted in that order. The
indices let §5 compare two files by number, since one sequence has one index whichever file
spelled it.

When any VCF genotype of a sample in the locus has a `.`, every possible genotype of the sample
has a missing haplotype, since each assignment puts the `.` on one haplotype or the other. So a
sample's possible genotypes are either all missing somewhere or none is, and §5 never meets a
mix.

### 4.4 Cases

Each case gives the reference, file A's records for one sample, and the result, which the tests
assert as literals. All of them were worked with the rules above in a throwaway script, and the
two `bcftools norm` outputs were run with bcftools 1.24 (`bcftools norm -f ref.fa in.vcf`).

**R1. The invalid example of `doc/spec.md` §5.4.** Reference 9–13 `CAGTC`, locus 9–12.

| records | result |
|---|---|
| `9 CAGT C 1/1`, `11 G A 0/1` | invalid, `contradicting records`: the SNP collides with the deletion on the haplotype that carries it |
| `9 CAGT C 0/1`, `11 G A 0/1` | alleles `CAGT,CAAT,C`; possible genotypes `(1, 2)`. Assignment 0 puts both on one haplotype, which is impossible |
| `9 CAGT C 0/1`, `11 G A,* 1/2` (GATK) | alleles `CAGT,CAAT,C`; possible genotypes `(1, 2)`: against the second row, *same genotype* |
| `9 CAGT C 0/1`, `9 C G,* 1/2` | invalid, `contradicting records`: a `*` on the deletion's anchor base alone is not covered |
| `11 G A,* 1/2` | invalid, `contradicting records`: no edit covers the `*` |
| `10 AG TC 0/1`, `11 G A,* 1/2` | locus 10–11, alleles `AG,AA,TC`, `(1, 2)`: the MNP replaces 11 and covers the `*` |
| `9 CAGT C ./.`, `11 G A,* 1/2` | possible genotypes `(MISSING, MISSING)`: the `*` might be covered |

**R2. A `*` after `bcftools norm`.** Reference 19–24 `TAAAAG`, locus 19–23. GATK's form,
`21 AA A 0/1` with `22 A C,* 1/2`, gives alleles `TAAAA,TAACA,TAAA` and `(1, 2)`. `bcftools norm`
writes the deletion as `19 TA T` and leaves the `*` at 22. That deletion deletes 20, but its
placements reach 23, so it covers the `*`, and the result is the same: *same genotype*. If the
span limit had cut the deletion's span to its REF, the only placement would delete 20 and the
sample would be invalid.

**R3. The first case of `doc/spec.md` §5.2.** Reference 21–26 `CAGAGT`, locus 21–25. As written,
`22 A T 1|0` with `23 GAG G 1|0`, the alleles are `CAGAG,CTG` and the possible genotypes `(0, 1)`.
`bcftools norm` writes the deletion as `21 CAG C 1|0`, which collides with the SNP at 22 on the
phased haplotype: invalid, `contradicting records`. Against each other, *invalid in B*.

**R4. Insertions and deleted bases.** Both records are `1/1`.

| reference, locus | records | result |
|---|---|---|
| 10–12 `CAG`, followed by `T` | `10 CAG C`, `11 A AT` | invalid: the T goes in between 11 and 12, both deleted |
| 10–14 `CAGTA` | `10 CAG C`, `12 G GT` | alleles `CAGTA,CTTA`, `(1, 1)`: the T goes in after 12, at the end of the deleted stretch |
| 10–11 `CG` | `10 C CA`, `10 C CT` | invalid: two insertions after 10 |
| 10–11 `CG` | `10 C A`, `11 G T` | alleles `CG,AT`, `(1, 1)`: two replacements side by side |

**R5. The insertion after the anchor.** Reference 51–54 `TAGC`, locus 51–53, `52 A C 1|0` and
`52 AG A,AAG 2|0`: the first haplotype has C at 52 and an A inserted after it, `TCAG`, and the
possible genotypes are `(0, 1)` with alleles `TAG,TCAG`. Trimmed from the right first, the
insertion would go before 52 and give `TACG`.

**R6. One insertion, two spellings.** Reference 10–12 `GCT`. `10 G GCA 1/1` in A and
`11 C CAC 1/1` in B both span 10–12, since the insertion slides from after 10 to after 11 and an
insertion's span holds the base after it; the locus is 10–12. Both give `GCACT`, alleles `GCT,GCACT`, `(1, 1)` in both: *same genotype*.

**R7. Phase.** Reference 30–32 `ACG`, SNPs `30 A T` and `32 G C`.

File B has the two SNPs phased on different haplotypes, `30 A T 0|1` and `32 G C 1|0`, which give
it one possible genotype, `ACC` and `TCG`.

| file A | possible genotypes of A, alleles | against B |
|---|---|---|
| `0/1`, `0/1` | `(0, 1)`, `(2, 3)`; `ACG,TCC,ACC,TCG` | *same genotype*, reported `(2, 3)` and `(2, 3)` |
| `0|1` PS 7, `1|0` PS 8 | `(0, 3)`, `(1, 2)`; `ACG,ACC,TCG,TCC` | *same genotype*, reported `(1, 2)` and `(1, 2)` |
| `0|1`, `0|1` | `(0, 1)`; `ACG,TCC,ACC,TCG` | *different*, reported `(0, 1)` and `(2, 3)` |

The alleles are listed in the order A's assignments meet them, so they differ between the first
two rows, and B's genotype has other indices in each. In the third row A has both SNPs on one
haplotype and B on two, so no genotype matches, and their ALT sets differ too.

Phase sets can interleave. Reference 30–34 `ACGTA`, `30 A T 0|1` PS 7, `32 G C 0|1` PS 8,
`34 A G 1|0` PS 7: PS 7 is block 0, by its record at 30, and PS 8 block 1. Assignment 0 gives
`ACGTG` and `TCCTA`, assignment 1 `ACCTG` and `TCGTA`: alleles `ACGTA,ACGTG,TCCTA,ACCTG,TCGTA`,
possible genotypes `(1, 2)`, `(3, 4)`.

**R8. Edits that cancel.** Reference 40–44 `CAAAG`, `40 C CA 0/1` and `41 AA A 0/1`, unphased.
The insertion of an A slides along the homopolymer and spans 40–44; the locus is 40–44. One
assignment puts the insertion and the deletion of an A on one haplotype, which is then the
reference, and gives `(0, 0)`; the other gives `(1, 2)`, alleles `CAAAG,CAAG,CAAAAG`. Against a
file with no records there, *same genotype*, since some possible genotype of A matches. The two
records cannot tell which, and `doc/spec.md` §5.5 asks only whether the files can agree.

**R9. The cap and the order of reasons.** Fourteen unphased `0/1` SNPs in one locus: *too
complex*; thirteen: reconciled. Fourteen `0|1` SNPs with one `PS`: reconciled, one assignment. A
sample with a malformed GT at one record and a GT of another ploidy at another: invalid,
`malformed GT`. A malformed GT, or a GT of another ploidy, in a sample with fourteen blocks: invalid, `malformed
GT` or `ploidy mismatch`. `./.` read with ploidy 1, and
`.` with ploidy 2: missing.

## 5. Evaluation

Evaluation takes the sample results of the two samples of a pair and gives the outcome of
`doc/spec.md` §5.5, the first rule that applies, and the pair of possible genotypes the row
reports. With *P*<sub>A</sub> and *P*<sub>B</sub> the two samples' possible genotypes, in their
order, where a genotype is **all REF** when every index in it is 0 and its **ALT set** is the set
of its indices other than 0:

| outcome | when | reported |
|---|---|---|
| partly outside the regions | `make_outcome_rows` is told the locus is (`doc/spec.md` §5.5) | none |
| too complex | either sample is too complex | none |
| invalid in A, in B, in both | that sample is invalid | none |
| missing in A, in B, in both | every genotype of that sample's *P* has a `MISSING` | the first of each |
| reference in both | *P*<sub>A</sub> and *P*<sub>B</sub> each hold one genotype, all REF | the first of each |
| same genotype | some genotype of *P*<sub>A</sub> equals one of *P*<sub>B</sub>, two genotypes all REF counting as equal | the first such pair |
| same alleles, other genotype | some genotype of *P*<sub>A</sub> and one of *P*<sub>B</sub> have the same ALT set | the first such pair |
| only in A | *P*<sub>B</sub> holds one genotype, all REF | the first of each |
| only in B | *P*<sub>A</sub> holds one genotype, all REF | the first of each |
| different | anything else | the first of each |

"The first such pair" takes the genotypes of *P*<sub>A</sub> in order and, for each, those of
*P*<sub>B</sub> in order, as `doc/spec.md` §5.5 says. The rules leave no pair without an outcome.
*Different* takes what is left, and in what is left both samples have a genotype with an ALT.
*Only in A* needs no condition on A: had A a genotype all REF, *same genotype* would have
applied. Within one ploidy two genotypes all REF are equal anyway; the words "counting as equal"
matter only across ploidies (decided by the owner on 19 October 2026).

*Partly outside the regions* depends on the locus and `--regions`, not on the sample results, so
`evaluate_pair` never gives it: `make_outcome_rows` gives it to every pair of a locus it is told
lies partly outside the regions, whatever `evaluate_pair` would give (decided by the owner on 22
October 2026). *Partly outside the regions*, *too complex* and *invalid* report no pair. In their
rows a sample that is too complex or invalid has no genotype, and a sample that was reconciled
shows its first possible genotype, as `doc/spec.md` §6 says column by column. So a row *invalid in A* still shows what B called. `PairEvaluation` holds only the reported pair;
`make_outcome_rows` takes the partner's genotype from its own sample result.

**Cases**, as `evaluate_pair` sees them, with `.` for `MISSING`:

| E | A | B | outcome | reported |
|---|---|---|---|---|
| 1 | too complex | invalid | too complex | none |
| 2 | invalid | `(0, 1)` | invalid in A | none |
| 3 | `(1, .)` | `(1, .)` | missing in both | `(1, .)`, `(1, .)` |
| 4 | `(., .)` | `(1, 2)` | missing in A | `(., .)`, `(1, 2)` |
| 5 | `(0, 0)` | `(0, 0)` | reference in both | `(0, 0)`, `(0, 0)` |
| 6 | `(0, 1)`, `(1, 2)` | `(1, 1)`, `(1, 2)` | same genotype | `(1, 2)`, `(1, 2)` |
| 7 | `(0, 1)`, `(0, 2)` | `(1, 1)`, `(2, 2)` | same alleles, other genotype | `(0, 1)`, `(1, 1)` |
| 8 | `(1,)` | `(1, 1)` | same alleles, other genotype | `(1,)`, `(1, 1)` |
| 9 | `(0, 1)`, `(0, 2)` | `(0, 0)` | only in A | `(0, 1)`, `(0, 0)` |
| 10 | `(0, 0)` | `(0, 0)`, `(1, 2)` | same genotype | `(0, 0)`, `(0, 0)` |
| 11 | `(0,)` | `(0, 0)`, `(1, 2)` | same genotype | `(0,)`, `(0, 0)` |
| 12 | `(0, 1)` | `(2, 3)` | different | `(0, 1)`, `(2, 3)` |
| 13 | `(0, 1)` | `(1, 2)` | different | `(0, 1)`, `(1, 2)` |
| 14 | `(0, 1)`, `(0, 2)` | `(2, 2)`, `(1, 1)` | same alleles, other genotype | `(0, 1)`, `(1, 1)` |

E3 is the owner's `./1` against `./1` (`doc/spec.md` §5.5). E10 is R8's sample, in file B, against a file A
with no records there. E11 is the same across two ploidies. In E13 both carry
locus allele 1, on one haplotype each, but B carries a second ALT, so the ALT sets differ. In E14 both `(0, 1)` with `(1, 1)` and `(0, 2)` with
`(2, 2)` share their ALT set; A's order decides, so the pair is `(0, 1)` and `(1, 1)`, where
taking B's genotypes first would give `(0, 2)` and `(2, 2)`.

## 6. The row of the table

`OutcomeRow` is one row of `doc/spec.md` §6, for one sample pair at one locus. Every evaluator
yields it (`doc/spec.md` §8), so its fields are fixed here, for all of them. It follows the Python
API (decided by the owner on 19 October 2026). Positions are 0-based and half-open, as in
`Locus`, and the table writer converts them to the 1-based ones of `doc/spec.md` §6. A value the table writes as `.`
is `None`, except a stratum, which is the loci spec's text, and the closed sets are enums whose values are the words the table writes.

| field | type | value |
|---|---|---|
| `contig`, `start`, `end` | `str`, `int`, `int` | the locus |
| `kind` | `Kind` | `locus_kind(locus)` (`doc/specs/loci.md`) |
| `stratum` | `str \| None` | the stratum as the loci spec gives it, `.` included for a locus no interval of the BED overlaps; `None` without `--strata` |
| `str_regions` | `tuple[GenomeSegment, ...]` | `locus.str_regions` |
| `sample_a`, `sample_b` | `str` | the names of the pair |
| `outcome` | `Outcome` | §5 |
| `alleles` | `tuple[str, ...] \| None` | the locus alleles, REF first; `None` from an evaluator that builds none |
| `genotype_a`, `genotype_b` | `LocusGenotype \| None` | the reported genotype; in a row *partly outside the regions*, *too complex* or *invalid*, the first possible genotype of a sample that has one (§5); `None` for a sample too complex or invalid, or from another evaluator |
| `possible_a`, `possible_b` | `int \| None` | how many possible genotypes; 0 for an invalid sample; `None` for a sample too complex, or from another evaluator |
| `quality_a`, `quality_b` | `float \| None` | the lowest QUAL among the records of the file in the locus at which the sample's VCF genotype has an index other than 0; records with no QUAL are passed over; `None` when no record is left |
| `gq_a`, `gq_b` | `int \| None` | the lowest GQ among the same records, likewise |
| `records_a`, `records_b` | `int` | how many records of the file in the locus have an ALT, whatever the sample carries at them |
| `invalid_reason_a`, `invalid_reason_b` | `InvalidReason \| None` | the sample's reason when it is invalid, whatever the outcome; `None` otherwise |

`quality` and `gq` count a `*` as an index other than 0, since it is an ALT of the record (§2 of
`doc/spec.md` says "non-reference VCF allele"). They are filled for a sample too complex or
invalid as for any other, from the records; a malformed or other-ploidy GT has no alleles and is
passed over. The reason is filled whenever the sample is invalid, also in a row *too complex*
because of its partner, so the table says what was wrong with each file. `doc/spec.md` §6 writes
it only for an invalid outcome, and is changed to match (§10).

## 7. The Python interface

The sample pairs come from the two files' headers and the map. `make_sample_pairs` raises the
errors of §3 and logs its warnings; it reads only the headers.

```python
# pvc_eval.samples
@dataclass(frozen=True, slots=True, kw_only=True)
class SamplePair:
    name_a: str
    index_a: int    # the sample's place among file A's sample columns, from 0
    name_b: str
    index_b: int

def make_sample_pairs(path_a: Path, path_b: Path, *,
                      sample_map: Path | None = None) -> tuple[SamplePair, ...]: ...
```

The types every evaluator shares are in a module that holds only them, as `pvc_eval.locus` holds
`Locus`:

```python
# pvc_eval.outcome
class Outcome(StrEnum):
    PARTLY_OUTSIDE = "partly outside the regions"
    TOO_COMPLEX = "too complex"
    INVALID_IN_A = "invalid in A"
    INVALID_IN_B = "invalid in B"
    INVALID_IN_BOTH = "invalid in both"
    MISSING_IN_A = "missing in A"
    MISSING_IN_B = "missing in B"
    MISSING_IN_BOTH = "missing in both"
    REFERENCE_IN_BOTH = "reference in both"
    SAME_GENOTYPE = "same genotype"
    SAME_ALLELES = "same alleles, other genotype"
    ONLY_IN_A = "only in A"
    ONLY_IN_B = "only in B"
    DIFFERENT = "different"

class InvalidReason(StrEnum):
    MALFORMED_GT = "malformed GT"
    PLOIDY_MISMATCH = "ploidy mismatch"
    CONTRADICTING_RECORDS = "contradicting records"

type LocusGenotype = tuple[int | MISSING, ...]   # sorted, MISSING last (§4.3)

@dataclass(frozen=True, slots=True, kw_only=True)
class OutcomeRow:
    ...   # the fields of §6, in its order
```

`reconcile` gives one locus's reconciliation. `locus_bases` is the reference over the locus, in
upper case, so its length is `locus.end - locus.start`; a test passes it as a string, and
`reconcile_loci`, stage 4 of `doc/spec.md` §4, reads it through a `FastaWindow` of its own. Each
sample result is keyed by the sample's index in its file's header. Both raise `ValueError` for a
ploidy other than 1 or 2, a phasing cap under 1, bases of another length than the locus, or a VCF
genotype with another number of alleles than its file's ploidy, which means the records were
read with another ploidy, and for a sample pair whose index is not one of a record's genotypes;
each is a defect of the caller, and the caller reconciles each file with the ploidy it read it
with. `Locus` raises `ValueError` when it is made with an end not after its start, a record on
another contig, or a record with an ALT whose span is not within the locus (`doc/spec.md` §4);
the loci build no such `Locus`, and a `Locus` built by hand is checked where it is made.

```python
# pvc_eval.reconcile
PHASING_CAP = 4096  # doc/spec.md §5.4; kept by the owner on 19 October 2026 (Open 2)

TOO_COMPLEX = sentinel("TOO_COMPLEX")
type SampleResult = tuple[LocusGenotype, ...] | InvalidReason | TOO_COMPLEX
# the possible genotypes, one or more, in the order of §4.3

@dataclass(frozen=True, slots=True, kw_only=True)
class ReconciledLocus:
    locus: Locus
    alleles: tuple[str, ...]                  # the locus alleles, REF first
    results_a: Mapping[int, SampleResult]     # the paired samples of file A, by header index
    results_b: Mapping[int, SampleResult]

def reconcile(locus: Locus, locus_bases: str, sample_pairs: Sequence[SamplePair], *,
              ploidy_a: int = DEFAULT_PLOIDY, ploidy_b: int = DEFAULT_PLOIDY,
              phasing_cap: int = PHASING_CAP) -> ReconciledLocus: ...

def reconcile_loci(loci: Iterator[Locus], reference: IndexedFasta,
                   sample_pairs: Sequence[SamplePair], *,
                   ploidy_a: int = DEFAULT_PLOIDY, ploidy_b: int = DEFAULT_PLOIDY,
                   phasing_cap: int = PHASING_CAP) -> Iterator[ReconciledLocus]: ...
```

`evaluate_pair` gives the outcome of §5 from two sample results. `make_outcome_rows`, stage 5,
gives a reconciled locus's rows, one per sample pair in the order of `sample_pairs`. It takes the
kind and the stratum as values, since they come from the locus alone (`doc/specs/loci.md`), and
`partly_outside`, true when `--regions` is given and the locus has a position outside them
(`Regions.covers`, `doc/specs/loci.md`), which makes every row *partly outside the regions*.

```python
# pvc_eval.evaluate
@dataclass(frozen=True, slots=True)
class PairEvaluation:
    outcome: Outcome
    genotype_a: LocusGenotype | None   # the reported pair (§5); None when none is reported
    genotype_b: LocusGenotype | None

def evaluate_pair(result_a: SampleResult, result_b: SampleResult) -> PairEvaluation: ...

def make_outcome_rows(reconciled: ReconciledLocus, sample_pairs: Sequence[SamplePair], *,
                      kind: Kind, stratum: str | None,
                      partly_outside: bool = False) -> Iterator[OutcomeRow]: ...
```

`reconcile` and `evaluate` are verbs the verb table of `.claude/skills/coding/python.md` lacks.
They are added to it in the commit of the code that first uses them, meaning "build a locus's
locus alleles and its samples' sample results (§4)" and "give the outcome of a sample pair (§5)".

## 8. How it runs

Reconciliation holds one locus, its bases and its sample results; nothing grows with the genome.
The work for one sample is its assignments, at most the phasing cap, times its haplotypes, times
the edits of a haplotype. A locus is reconciled apart from every other, so loci can later be
spread over threads without changing a result. How long 4096 assignments take is not measured
yet; the performance work measures it on a locus of 13 blocks, and the owner lowers the default
cap if it costs too much (`doc/spec.md`, Open 2).

## 9. How it is verified

Each check is made at the highest function at which it can be seen, with its expected result as
a literal.

- **At `make_sample_pairs`**, over headers written to `tmp_path`: pairing by name in A's order;
  a map's order, a sample in two pairs, a blank line; every error and the warning of §3, with
  their text; the first 10 names and the count after them.
- **At `reconcile`**, on `Locus` values built in the test, with their records read by `parse_vcf`
  from a few lines written to `tmp_path`, so that each record has the span the reader gives it:
  every case of §4.4, R1 to R9, with its locus alleles and its possible genotypes. Then the
  order of the locus alleles across two samples and two files; a sample in no pair, which adds
  nothing; a file with no record in the locus, all REF; a record with no ALT whose REF reaches
  past the locus, `0/0` (all REF) and `./.` (missing); and ploidy 1, with a `1` and a `0`.
- **At `evaluate_pair`**: every case E1 to E14 of §5.
- **At `make_outcome_rows`**: every field of §6 on one locus with two pairs, one of them sharing a
  sample; the partner's genotype in a row *invalid in A*; the reason in a row *too complex*;
  with `partly_outside`, a row *partly outside the regions* for a pair that would be *same
  genotype*, one that would be *reference in both*, one with an invalid sample and one with a
  sample too complex, each showing each sample's first possible genotype, or none for the
  invalid or too complex one, and the reason;
  `quality` and `gq` with a QUAL absent, and with a `*`.
- **At `parse_vcf`**: `1` read with ploidy 2 is `OTHER_PLOIDY_GT`, with ploidy 1 a genotype;
  `1/1` read with ploidy 1 is `OTHER_PLOIDY_GT`; `./.` read with ploidy 1 and `.` with ploidy 2
  are `MISSING` on every haplotype; `3/1` with two ALTs is `MALFORMED_GT`; a ploidy of 0 raises.
- **A property**, with hypothesis, at `reconcile` and `make_outcome_rows`: file B is file A
  re-spelled, each record decomposed into records of one replaced base or one indel, phased,
  with a phase set of its own for each record of A so that no block is added (`doc/spec.md`
  §7.4). A record is decomposed only where none of the edits a haplotype carries at it touches
  an edit the same haplotype carries at another record, as `doc/spec.md` §7.4 says, edits as
  §4.2 makes them, trimmed and where written: two records of one haplotype that touch
  can give a file and its decomposition different outcomes under the collision rules of §4.2,
  which judge edits as written. Reference `TAAC` with `3 AC G,AA 0/1` and `3 A G,AA 1/2` is
  invalid, and decomposed is not; reference `GTTT` with `2 TT GCTA 1/1` and `3 T G,TC 2/2` is
  valid, and decomposed is invalid (decided by the owner on 19 October 2026). Records that touch
  are kept whole in B. The records of both are read by `parse_vcf` over a generated reference, and the `Locus`
  is built in the test over the spans of both, with no `--regions`. Every row is then *same genotype*, *reference in
  both*, *missing in both*, *invalid in both* or *too complex*, the last three only where A
  itself has a missing haplotype, an invalid sample or a sample over the cap (`doc/spec.md`
  §7.3). Each sample's result in B is also the same as in A, the same possible genotypes with
  the same indices, or the same reason: B has A's blocks, one phase set for each record of A,
  and its haplotypes have A's sequences. The outcome alone would let a B that lost phase
  assignments pass, since *same genotype* needs only one pair to match. B being A itself would check little, since both then get the same results by
  construction. The other re-spellings of §7.4 move records, and are tested over the whole
  pipeline once `build_loci` exists.

## 10. Changes to `doc/spec.md`

The owner's decisions change `doc/spec.md` in places outside its §5.4 and §5.5; the sections
below are its own. They go in a commit of
their own, after this spec is approved:

- §2, *ploidy*: the number of haplotypes each sample of a file has at every locus, 1 or 2, given
  per file, 2 by default.
- §3: the options `--ploidy-a N` and `--ploidy-b N`, 1 or 2, default 2, and `CompareConfig`
  holds them; and three more errors that stop the run, from §3 of this spec: a line of
  `--sample-map` that is not two fields, a pair given twice, and a sample name given twice in a
  header that a pair needs.
- §5.4: an insertion collides with an edit that replaces both bases beside it; a `*` is covered
  by a deletion the haplotype might carry where it is missing, and by any other edit that
  replaces a position of its record's REF; a GT of another ploidy, in place
  of a haploid GT, with the reason `ploidy mismatch`, and a GT of only `.` missing whatever its
  number of alleles; the shared locus alleles taken from the paired samples, not from "any sample
  of either file"; a pointer to this spec.
- §5.5: a row *too complex* or *invalid* shows the reconciled partner's first genotype; two
  genotypes all REF of different ploidies count as equal; a pointer to this
  spec.
- §6: `ploidy mismatch` in place of `haploid GT`, and the reason written whenever the sample is
  invalid.
- §7.2: "ploidy mismatch, and a file of ploidy 1" in place of "haploid GT".
- §9: ploidy other than 1 or 2.
- §10: Open 1 closed with the owner's answer.

The VCF reader changes as §2 says, in a milestone of its own before reconciliation.

## 11. Open points

No point is open. The two the review raised were answered by the owner on 19 October 2026:
**Open 4**, a `*` under an edit that is not a deletion, which GATK writes with MNPs turned on
(`AssemblyBasedCallerUtils.java` at GATK commit e8c49f6), is covered by that edit (§4.2); and
**Open 5**, two genotypes all REF of different ploidies count as equal (§5).

## 12. Not in this spec

- `compare()`, its evaluator interface, the command line and its options: the spec of `compare`.
- Writing `OutcomeRow` as a TSV line, and the summary: the table's own milestone, from
  `doc/spec.md` §6.
- The kind, the stratum and `--regions`: `doc/specs/loci.md`.
- The wrappers of other evaluators, which yield `OutcomeRow` with `alleles`, `genotype_*` and
  `possible_*` set to `None`: their specs under `doc/specs/` (`doc/spec.md` §8).
- Running loci on threads: the performance work, once a baseline is measured.
