# The re-spelling program — specification

*Status: approved by the owner on 20 October 2026, who answered Open 1 the same day: keep the
rule. No code of it exists yet. The tests
of reconciliation already decompose records, in `_make_decomposed_lines` of
`tests/test_reconcile_property.py`, and the decomposition here replaces it. It builds the program
`doc/spec.md` §7.4 asks for, and stands on `doc/spec.md` §2, §5.2, §5.4 and §7, on
`doc/specs/reconciliation.md` §4 for the edits, phase units and blocks, and on pvc_eval's VCF
reader, `parse_vcf`, which gives each record normalised, with its span and its VCF genotypes.*

`doc/spec.md` §7.4 checks pvc_eval's central claim: that two spellings of the same genotypes
agree. It compares a file with re-spellings of itself, and every row must be in the set allowed
for a file against itself (§7.3). The re-spelling program makes those re-spellings. It writes a
VCF file again in other words, in four ways: each record decomposed into one-base and one-indel
records, indels shifted right within their repeats, `*` alleles added under deletions as GATK
writes them, and neighbouring records merged into complex ones. A re-spelling that changed a
haplotype would make §7.4 fail where pvc_eval is right, so the program is checked by an oracle of
its own, written for it, that does not use the code §7.4 tests. This spec gives the four re-spellings, which records each one rewrites and which it keeps
as they are, the files they write, the oracle, and how all of it is checked. It all lives under
`tests/`, as §7.4 says, and adds nothing to `src/`.

The owner decided on 19 October 2026 that collisions are judged on edits as written, so where
two edits of one haplotype touch, the outcome can follow the spelling (`doc/spec.md` §7.4). Every
re-spelling therefore rewrites only records none of whose edits touches an edit of another record
on the same haplotype. Beyond that rule, this spec decides:

- that the pieces of a phased VCF genotype stay in its phase set, and only an unphased one gets
  a phase set of its own, so that no re-spelling adds or removes a block;
- when two records "are on the same haplotype" although the phase is not known: unless both are
  phased in one phase set on different haplotypes, they may be, and that counts;
- that a record near a `*` is kept as it is by every re-spelling, since its spellings can cover
  the `*` differently;
- that a rewritten record keeps only CHROM, POS, REF, ALT, QUAL and the GT and PS of each sample;
- what the oracle compares where the phase is not known (§8).

The owner decided on 20 October 2026 (**Open 1**, at the end) that the `*` re-spelling does not
write GATK's commonest case, a deletion and a SNP under it, both heterozygous and unphased.

Positions in the cases are 1-based and inclusive, as in a VCF. In the Python interface they are
0-based and half-open, as everywhere inside pvc_eval.

## 1. Vocabulary

The terms of `doc/spec.md` §2 keep their meaning, and so do those of
`doc/specs/reconciliation.md` §1: an **edit** is what one VCF allele does to the reference once
the bases it shares with REF are trimmed, left first; a **phase unit** is the VCF genotypes of one
sample whose order on the haplotypes is fixed together, one unphased VCF genotype or every phased
VCF genotype of one phase set; a **block** is a phase unit with a heterozygous VCF genotype in
it; a **placement** is one of the stretches a deletion can slide to (`doc/spec.md` §5.2). This
spec adds:

| term | meaning | not |
|---|---|---|
| **re-spelling** | a VCF file written from another so that every sample's haplotypes have the same sequences and the same phase; also each of the four ways of making one | "normalisation", which is §5.2's trimming |
| **rewritten record** | a record of the input that a re-spelling replaces by records it writes itself | |
| **kept record** | a record of the input that a re-spelling writes as the input has it, byte for byte | "skipped" |
| **carries** | a haplotype carries the edit of the VCF allele its VCF genotype gives it at a record; `.` and REF carry none | |
| **may share a haplotype** | two VCF alleles of one sample, at two records, may be on one haplotype unless both VCF genotypes are phased in one phase set and the alleles are on different haplotypes. Two records in different phase sets, or unphased, may. A sample of ploidy 1 has one haplotype, which all its alleles share. | |
| **touch** | as `doc/spec.md` §2: two stretches touch when they overlap or lie next to each other. The stretch of an insertion is the point it goes in at, between two positions, so it touches both. | |
| **reach** | the stretch a record's spellings can take: its REF, widened for each ALT whose edit is an insertion or a deletion by every position it can slide to (`doc/spec.md` §5.2), with no limit. The re-spellings work it out with a slide of their own, for records of the input and for the records they would write | "span", which `parse_vcf` gives and which the span limit, 10 kb by default, cuts back to the REF |
| **window** | the oracle's unit: a stretch of a contig made by joining the spans of the records with an ALT of both files while they overlap or lie next to each other (§8) | "locus", which is §5.3's |

## 2. What every re-spelling does

### 2.1 Reading

A re-spelling reads its input twice: as text, line by line, to keep the header and the kept
records as they are; and with `parse_vcf`, with the reference and the ploidy it is given, for the
records normalised, their spans and their VCF genotypes. Each record of `parse_vcf` is matched to
its line by the line number it carries. Both read with `DEFAULT_SPAN_LIMIT`, the default of
`doc/spec.md` §5.2, as the check of §7.4 does. `parse_vcf` yields records in order of span start, so the
re-spelling puts them back in the order of the file first.

Whatever stops `parse_vcf` stops the re-spelling, with the same `PvcEvalInputError`, pvc_eval's
error for a fault of the input, before it
writes anything: a REF that differs from the reference, a contig the reference lacks, a REF past
its contig's end, records out of order. A file that pvc_eval cannot read has no re-spelling.

### 2.2 The records every re-spelling keeps

Some records are kept by all four re-spellings, because rewriting them would change what the
file says or would test nothing:

| record | why it is kept |
|---|---|
| a record with no ALT (`.`) | it makes no edit, and its `./.` makes a haplotype missing wherever its REF reaches (`doc/spec.md` §5.2) |
| a record with a `*` ALT | a `*` makes no edit and says another record's edit covers it |
| a record whose reach, or the reach of any record its rewriting would write, overlaps or touches the REF as written of a record with a `*` ALT, the stretch its file's REF covers before trimming (`doc/specs/reconciliation.md` §4.2, decided by the owner on 21 October 2026) | reconciliation covers a `*` by every placement of a deletion and by the edits of a missing haplotype (`doc/specs/reconciliation.md` §4.2), so a rewritten neighbour could cover a `*` that was not covered, or the reverse |
| a record with a symbolic allele | `parse_vcf` skips it, and so does pvc_eval (`doc/spec.md` §9); it is written where it stood |
| a record at which some sample's GT is malformed, or of another ploidy than the one the file is read with | the sample is invalid wherever the record is (`doc/specs/reconciliation.md` §4), and its GT has no alleles to carry into pieces |
| a **touched** record: one at which some sample carries an edit that touches an edit the same sample carries at another record of the contig, the two VCF alleles being ones that may share a haplotype | where edits of one haplotype touch, the outcome can follow the spelling (owner, 19 October 2026) |

The edits compared are those of `doc/specs/reconciliation.md` §4.2, as written: the edit of each
ALT trimmed left first, where the record stands. A `.` carries no edit, as in reconciliation's
collisions, whichever haplotype it is on. The `*` rows above take care of the one place where reconciliation reads the edits
of a missing haplotype.

A missing haplotype, `.`, is otherwise no reason to keep a record: the re-spellings carry it into
what they write, as each section says. So is a record at the first or the last base of its
contig: a VCF allele there is written as VCF writes it (§3), and an indel slides up to the
contig's end and no further (§4).

Each re-spelling has its own conditions too, in its section. A record that does not meet them is
kept.

### 2.3 Writing

**Kept records** are written as the input has them, byte for byte, except the line end, which is
`\n`.

**Rewritten records** are written with these columns: CHROM; POS; ID `.`; REF and ALT, upper
case; QUAL as the input record has it, or for a merged record the lowest QUAL of the records
merged, `.` when none has one; FILTER `.`; INFO `.`; FORMAT `GT:PS`; and for each sample its GT
and its PS. A file with no sample columns has no FORMAT column either. The GT is written with `|` when the VCF genotype is phased and with `/` when it is
not, unless the re-spelling's section changes its phase; a GT of ploidy 1 has no separator. The PS is the phase set of a phased GT, `.` when it has
none or the GT is unphased. Nothing else of the input record is kept: pvc_eval reads no ID,
FILTER or INFO, and of the FORMAT fields only GT, PS and GQ, the genotype quality, which only
fills the `gq` columns of `doc/spec.md` §6, which §7.4 does not check.

**A new phase set** is a PS that no record of the input uses: the n-th whole number, counting
from 1, that is not the PS of any record of the input, where n is the place of the record it is
made for among the input's record lines, from 1. Each record gets its own, and a sample's new PS
is the record's. A PS of the input is compared as `parse_vcf` compares it, as text, so `007`
does not take 7. A new phase set is at most the number of record lines plus the number of
distinct PS values of the input, whatever those values are; numbering from the input's largest
PS, as the spec first said, went past bcftools' 32-bit limit and Python's 4,300 digits on inputs
with a PS near either (the review of milestone 5). Where the input has no PS, the n-th number is n.

**The header** is the input's, line by line, with up to two lines added before `#CHROM` when the
output writes a rewritten record. The first, when a rewritten record has sample columns and the
header has no `##FORMAT=<ID=GT,` line, since a rewritten record writes GT whatever the input's
FORMAT holds:

```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
```

The second, when the header has no `##FORMAT=<ID=PS,` line, and every PS of the input and every
new phase set the output writes is a whole number of at most 2,147,483,647, the largest that
bcftools 1.24 keeps in an `Integer`; it sets 2,147,483,648 to missing, with a warning (seen on 20
October 2026):

```
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
```

bcftools 1.24 warns about a FORMAT field the header lacks, even with only `.` as its values. It
stops with exit status 8 on a PS such as `chr1_5` declared `Integer`, and warns that PS "should be
declared as Type=Integer" when it is declared `String`. An input with such a PS, or one above
2,147,483,647, and no PS line already gets the warning that PS is not in its header, so the
output adds no line and gets the same warnings (seen by the reviews on 20 October 2026, which found that the first version of this
rule, a `Type=String` line, broke the rule of validity below).

**The order.** Records are sorted by contig, in the input's order of contigs, and by POS, which is
what `parse_vcf` asks (`doc/spec.md` §5.2). Records of one POS keep the order of the input
records they come from, and the records written for one input record keep the order of their
section.

**Compression.** When the input starts with the two bytes of gzip, the output is bgzipped, written
by a BGZF writer of the tests, so that `parse_vcf` and `bcftools` read it as they read the input;
otherwise it is plain text. The name of the output is the caller's.

**Valid.** `parse_vcf` reads the output with the same reference and ploidy without an error, and
`bcftools view` (1.24) reads it with exit status 0 and no message on standard error that it does
not also give for the input, the messages compared with their positions removed, since bcftools
writes one in each (`FORMAT 'GT' at chr1:5 is not defined in the header`) and a re-spelling moves
records. bcftools 1.24 warns, and goes on, when the header lacks a FORMAT line or a contig line
(seen on 20 October 2026 with a header of two lines); those warnings come from the input's
header, and the output adds no other.

## 3. Decomposition

**What it does.** Each record that is not kept is replaced by records of one replaced base or one
indel, which are its **pieces**, with the VCF genotypes carried over.

**The pieces.** Each ALT's edit is found as reconciliation finds it, and its replaced stretch and
its new bases are aligned from the left: a SNP at each position where they differ, then the
bases left over, as an insertion after the stretch's last aligned position or as a deletion of
the positions left over. An indel is written with the base before it as its anchor; at the
contig's first base, which has none before it, with the base after it, at POS 1, as VCF writes it.
A piece two ALTs share is one record. An ALT that spells REF makes no piece, and a record none of
whose ALTs makes a piece is kept.

**The genotypes.** At each piece, each haplotype gets `1` when the ALT its VCF genotype gives it
makes that piece, `0` when it is REF or another ALT, and `.` when it is `.`. So a record of A that
is `1/2` gives `1|0` at the pieces of ALT 1 alone.

**The phase.** The pieces of one record are one phase unit, as the record was, so the
decomposition adds no block. For each sample:

- a VCF genotype phased in a phase set keeps it, `|` and the same PS, so its pieces stay in the
  phase unit of its other records;
- a phased one with no PS stays phased with none;
- an unphased one is written phased, `|`, in the record's new phase set (§2.3), which no other
  record uses: one unit, as before;
- with ploidy 1 there is no phase, and the GT is written as one allele.

**The cases**, D1 to D12, each a test with its output as a literal. Each gives the reference from
position 1, the input, and the output, on `chr1`, as `POS REF ALT` and one sample field per
sample, `GT:PS` for a rewritten record. A record written *as read* is kept. The cases of the
other sections, S, T, M and O, are written the same way.

| D | reference | input | output |
|---|---|---|---|
| 1 | `GACGTA` | `2 ACGT TCGAA 0/1` | `2 A T 0\|1:1`, `5 T A 0\|1:1`, `5 T TA 0\|1:1`: aligned from the left, `ACGT` and `TCGA` differ at 2 and 5, and the A left over is inserted after 5 |
| 2 | `GACGT` | `2 ACG TCA,TCG 1/2` | `2 A T 1\|1:1`, `4 G A 1\|0:1`: the SNP at 2 is one record, which both ALTs make |
| 3 | `GACGTA` | `2 AC GT 1\|0:7`, `5 T G 0\|1:7` | `2 A G 1\|0:7`, `3 C T 1\|0:7`, `5 T G 0\|1:7`: the pieces stay in phase set 7 |
| 4 | `GACGTA` | `2 AC GT 0/1`, `4 G A 0/1` | both as read: the edits of 2–3 and 4 lie next to each other, and unphased they may share a haplotype |
| 5 | `GACGTA` | `2 AC GT 1\|0:5`, `4 G A 0\|1:5` | `2 A G 1\|0:5`, `3 C T 1\|0:5`, `4 G A 0\|1:5`: the same edits, phased on different haplotypes, touch on none |
| 6 | `GACGTA` | `2 AC GT 1\|0:5`, `4 G A 0\|1:6` | both as read: in two phase sets the phase between them is not known |
| 7 | `GACGTA` | `2 AC GT 0/1 0/0`, `4 G A 0/0 1/1` | `2 A G 0\|1:1 0\|0:1`, `3 C T 0\|1:1 0\|0:1`, `4 G A 0\|0:2 1\|1:2`: the edits touch only across samples |
| 8 | `GACT` | `2 AC GT ./1` | `2 A G .\|1:1`, `3 C T .\|1:1` |
| 9 | `GACT`, ploidy 1 | `2 AC GT 1` | `2 A G 1:.`, `3 C T 1:.` |
| 10 | `ACGT` | `1 AC C 0/1` | `1 AC C 0\|1:1`: the deletion of the first base is anchored on the base after it |
| 11 | `ACGT` | `1 A CA 0/1` | `1 A CA 0\|1:1`: so is an insertion before it |
| 12 | `TAAC` | `3 AC G,AA 0/1`, `3 A G,AA 1/2` | both as read: the insertion of an A after 3 falls between the two positions `AC → G` replaces, which collides, while decomposed it would not (`doc/specs/reconciliation.md` §9) |

D1, D2, D7, D8, D10, D11 and D12 are the cases of `_make_decomposed_lines`, whose output differs in
two ways: a kept record is now written as read, not with its GT rewritten, and the lines of one
POS come in the order of the input, not in order of span start.

## 4. Shift

**What it does.** An indel inside a repeat can be written at any of its placements, and callers
choose differently: `bcftools norm` and most callers write the leftmost, some the rightmost. The
shift writes indels further right than the input has them.

**Which records.** A record that is not kept, with one ALT, whose edit is an insertion or a
deletion that can slide right. A SNP, a complex record and a record of several ALTs are kept:
the ALTs of one record would slide by different distances, and one record has one POS.

**How far.** A deletion of positions *a* to *b* slides right one position while the base after
*b* equals the base at *a*; an insertion slides right while the base after its point equals its
first inserted base, which then moves to its end (`doc/spec.md` §5.2). The record is moved by a
distance *k* of 1 or more, at most the largest one meeting all of these:

- **the span**: the placement lies within the record's span as `parse_vcf` gives it. A span the
  span limit cut back to the REF allows no shift, since the shifted record's span would be its
  new REF and no longer meet the original;
- **the contig**: no slide goes past the contig's last base;
- **the bases**: the reference over the deleted bases and over every position the indel slides
  past is A, C, G or T, and the deleted bases are the reference's own. `parse_vcf` accepts a REF
  `A` where the reference has `M` or `N` (`doc/spec.md` §5.2), while reconciliation builds each
  haplotype from the reference's bases; sliding on the VCF's letters would delete an `A` in place
  of the `M` and change the haplotype (`GCMAAAT` with `2 CA C`, found by the review of
  milestone 4);
- **the phase**: for every sample, and every haplotype carrying the indel, no edit of another
  record that may share that haplotype touches the stretch the indel crosses: positions *a* to
  *b* + *k* for a deletion, and the positions the point passes over for an insertion. An edit
  carried only by another sample, or phased on the other haplotype, does not stop it; one whose
  phase is not known does.

Without a seed the record is moved by the largest *k*; with a seed, by a *k* drawn evenly from 1
to the largest, so that the placements in between are written too.

**What is written.** The indel at its new placement, with the base before it as anchor, and the
record's VCF genotypes as read: the shift moves an edit, not a haplotype's phase.

| S | reference | input | output |
|---|---|---|---|
| 1 | `GCAAAAT` | `2 CA C 0/1` | `5 AA A 0/1:.`: the deletion of the A at 3 moves to the A at 6 |
| 2 | `TCACACAG` | `2 C CAC 1/1` | `7 A ACA 1/1:.`: `AC` inserted after 2 is `CA` inserted after 7, rotated as it slides |
| 3 | `GCAAAAT` | `2 CA C 0/1`, `5 A G 0/1` | both as read: deleting 4 would lie next to the SNP at 5, which may share its haplotype |
| 4 | `GCAAAAAT` | `2 CA C 0/1`, `7 A G 0/1` | `4 AA A 0/1:.`, `7 A G 0/1` as read: the deletion moves two positions, and a third would lie next to the SNP |
| 5 | `GCAAAAT` | `2 CA C 1\|0:9`, `5 A G 0\|1:9` | `5 AA A 1\|0:9`, then `5 A G 0\|1:9` as read: the SNP is on the other haplotype |
| 6 | `GCAAAAT` | `2 CA C 0/1 0/0`, `5 A G 0/0 0/1` | `5 AA A 0/1:. 0/0:.`, `5 A G 0/0 0/1` as read: the SNP is another sample's |
| 7 | `GCAAAA`, the whole contig | `2 CA C 0/1` | `5 AA A 0/1:.`: the deletion of the contig's last base |
| 8 | `GCAAAA`, the whole contig | `2 C CA 0/1` | `6 A AA 0/1:.`: the insertion after the contig's last base |
| 9 | `GCAAAAT` | `2 CA C,CAA 1/2` | as read: two ALTs |

Shifting changes no span, since every spelling of an indel with its anchor on the left has the
same span (`doc/spec.md` §5.2), so the loci of a file and its shift are the same. What the shift
tests is the reconciliation of edits that stand at different placements in the two files.

## 5. Spanning deletions

**What it does.** Where a haplotype carries a deletion, a record whose REF it deletes cannot give
that haplotype a base there. A file that writes the record with REF on that haplotype, as most
callers do, and one that writes `*`, as GATK does, say the same thing (`doc/spec.md` §5.4). This
re-spelling writes the `*`.

**Which records.** A record R that is not kept, with an ALT, gets a `*` on haplotype *h* of a
sample when all of these hold:

- *h* carries, at another record D, a deletion: an edit that replaces one or more positions with
  no bases, which deletes, as written, a position of R's REF as `parse_vcf` trims it. A `*` is
  measured against its record's REF as written (`doc/specs/reconciliation.md` §4.2), which holds
  the trimmed REF, so such a deletion covers it however the record is written. The trimmed REF is
  asked for so that the re-spelling adds a `*` only where it would have been covered before the
  rule of 21 October 2026 too: `12 TC TA` is `13 C A` trimmed, and a `*` there for a deletion of
  12 alone is now covered, but is not added;
- R's VCF genotype gives *h* REF, `0`. A `.` is left as it is. An ALT on *h* would be an edit
  touching the deletion, and R would be kept;
- the phase between R and D on *h* is known, in one of three ways:
  1. both VCF genotypes are phased in one phase set;
  2. R gives REF on every haplotype of the sample. R is then written in D's phase unit, which
     ties it to D without adding a block, since an all-REF VCF genotype is no block: in D's
     phase set when D is phased; in a new phase set made for D, given to both, when D is
     unphased and heterozygous; and with its phase as read when D is unphased and homozygous,
     since every haplotype then gets a `*`. When the deletions over R of that sample come from
     records of two or more phase units, R would have to join both, and that sample gets no `*`
     at R, even on a haplotype where the first way alone would give one: `10 AG A 1|0:5`,
     `11 GTC ACA 0|0:5`, `12 TC T 0|1:6` is written as read (the review of milestone 5 found
     that the other reading gives no wrong outcome either, and this one is the simpler);
  3. the ploidy is 1.

  When D is unphased and heterozygous, D is rewritten as well, with its VCF genotype phased in
  its new phase set, and so D must not be a kept record.

**What is written.** R with `*` added as its last ALT, and the `*`'s index on every haplotype
that gets one; D, when rewritten, with only its phase changed. Neither adds a block: R joins D's
phase unit only where R is all REF, which is no block, and D keeps its own unit under a new name.

The cases use the reference of `doc/spec.md` §5.4, `CAGTC` at 9–13, with `N` before it, so
`9 CAGT C` deletes 10–12 and slides nowhere.

| T | input | output |
|---|---|---|
| 1 | `9 CAGT C 0\|1:5`, `11 G A 1\|0:5` | `9 CAGT C 0\|1:5` as read, `11 G A,* 1\|2:5` |
| 2 | `9 CAGT C 1/1 0/0`, `11 G A 0/0 0/1` | `9 CAGT C 1/1 0/0` as read, `11 G A,* 2/2:. 0/1:.`: the record under the deletion is there for the second sample |
| 3 | `9 CAGT C 0/1 0/0`, `11 G A 0/0 0/1` | `9 CAGT C 0\|1:1 0/0:.`, `11 G A,* 0\|2:1 0/1:.`: the deletion and the record under it share the deletion's new phase set |
| 4 | `9 CAGT C 0\|1:5`, `11 G A 0/0` | `9 CAGT C 0\|1:5` as read, `11 G A,* 0\|2:5` |
| 5 | `9 CAGT C 1`, `11 G A 0`, ploidy 1 | `9 CAGT C 1` as read, `11 G A,* 2:.` |
| 6 | `9 CAGT C 1\|0:5`, `11 G A .\|1:5` | both as read: the haplotype with the deletion is `.` at 11 |
| 7 | `9 CAGT C 0/1`, `11 G A 0/1` | both as read: unphased, the SNP may share the deletion's haplotype, where it collides with it (Open 1, answered) |
| 8 | `9 CAGT C 0/1`, `11 G A,* 1/2` | both as read: a `*` is there already |

Only records already in the file get a `*`. GATK writes a `*` only at a record some sample needs
for its own allele, so a multi-sample file is where most of them appear, as in T2.

## 6. Merge

**What it does.** Two callers can write the same events as one complex record or as several
neighbouring ones. The merge joins runs of neighbouring records into one complex record.

**A run** is two or more records, consecutive in the order of §2.3, of one contig, each one's REF
overlapping the stretch of the records before it, or starting after it with at most 2 positions
between them (M1 has 1, M8 has 3). The run's
stretch runs from the first position of its REFs to the last; the REFs are the normalised ones of
`parse_vcf`. A run is merged only when all of these hold:

- no record of it is kept, has a `*` ALT, or gives any sample a `.`. A `.` would make the merged
  allele of its haplotype missing and drop what the other records say;
- for every sample with a heterozygous VCF genotype in the run, the merged alleles of its
  haplotypes differ. Phased records whose edits spell the same sequence on both haplotypes, such
  as one deletion of a repeat written at two placements, `2 ACA A 1|0:5` and `4 ACA A 0|1:5` on
  `GACACAT`, would merge into a homozygous record and remove a block, which turns a sample over
  the phasing cap in file A into one under it in file B (the review of milestone 6);
- the reference over the run's stretch is A, C, G or T. The merged REF holds the bases between
  the records' REFs, which no input REF covered and `parse_vcf` never checked; a letter that is no
  IUPAC code there, such as `X`, makes `parse_vcf` reject the output (the review of milestone 6);
- no other record with an ALT, symbolic ones aside, since `parse_vcf` does not read them, has a
  REF overlapping or touching the run's stretch, since the
  merged record replaces the whole stretch, gaps included, and its edit would touch it;
- for every sample, the heterozygous VCF genotypes of the run are all in one phase unit: one
  phase set, or a single unphased record. Two unphased heterozygous records give two possible
  phases, and merged they would give one: the merge would invent a phase. Homozygous VCF genotypes
  go with any, and with ploidy 1 there is no phase.

**What is written.** One record, at the first position of the stretch, whose REF is the reference
over the stretch. Each haplotype's allele is that REF with the edits of the VCF alleles it
carries at the run's records applied, which the second condition makes possible, since no two
of them touch. The ALTs are the distinct alleles that differ from REF, in order of first
appearance, samples in the order of the header and haplotypes in order; a haplotype whose allele
is REF gets `0`. The GT is phased, in the phase set of the run's heterozygous VCF genotypes, when
those are phased, and unphased otherwise. A merge across samples can need as many ALTs as there
are distinct haplotypes, as M5 shows.

A run whose every haplotype gets REF, all `0/0`, has no ALT to write and is not merged.

Without a seed, the runs are as long as the conditions allow, taken from the left. With a seed,
each record that could join the run before it does so with probability 1/2, so that runs of
every length are written.

| M | reference | input | output |
|---|---|---|---|
| 1 | `GACGTA` | `2 A T 1/1`, `4 G C 0/1` | `2 ACG TCG,TCC 1/2:.` |
| 2 | `GACGTA` | `2 A T 1\|0:5`, `4 G C 0\|1:5` | `2 ACG TCG,ACC 1\|2:5` |
| 3 | `GACGTA` | `2 A T 0/1`, `4 G C 0/1` | both as read: two unphased heterozygous records |
| 4 | `GACGTA` | `2 A T 1\|0:5`, `4 G C 0\|1:6` | both as read: two phase sets |
| 5 | `GACGTA` | `2 A T 1\|0:5 1/1`, `4 G C 0\|1:5 1/1` | `2 ACG TCG,ACC,TCC 1\|2:5 3/3:.` |
| 6 | `GACGTA` | `2 A T 0\|1:5`, `4 G C 0\|1:5` | `2 ACG TCC 0\|1:5`: the first haplotype is REF |
| 7 | `GACGTA` | `2 A T ./1`, `4 G C 1/1` | both as read: a `.` |
| 8 | `GACGTA` | `2 A T 1/1`, `6 A C 1/1` | both as read: 3 positions lie between them |
| 9 | `GACGTA` | `2 A T 1/1`, `4 G C 1/1`, `6 A G 0/1` | `2 ACGTA TCCTA,TCCTG 1/2:.`: a run of three |
| 10 | `GACGTA`, ploidy 1 | `2 A T 1`, `4 G C 0` | `2 ACG TCG 1:.` |
| 11 | `GACGTA` | `2 A AT 1/1`, `4 G C 0/1` | `2 ACG ATCG,ATCC 1/2:.`: an insertion in the run |
| 12 | `GACGTA` | `2 A T 1/1`, `3 C A 2/0`, `4 G C 0/1` | all as read: the malformed GT at 3 keeps its record, which touches both stretches |

The merged record's span is its REF, or wider where an ALT is a pure indel (`doc/spec.md` §5.2),
so a merge can change which records share a locus. The loci are built from both files together
and hold both spellings, which is what §7.4 tests.

## 7. The Python interface

The program is in `tests/respell.py`, the oracle in `tests/haplotype_oracle.py`, and the check of
§7.4 in `tests/respell_check.py`; the tests import them as `tests.conftest` is imported now.
`respell` joins the table of `.claude/skills/coding/python.md` that gives each verb of a function
name one meaning, in the commit of the code that
first uses it, meaning "write a VCF file again in other words that give every sample the same
haplotypes (`doc/specs/respell.md`)".

Each re-spelling reads `vcf`, a VCF file, plain or gzipped, with `parse_vcf` over `reference` and
`ploidy`, and writes the re-spelling to `output`, bgzipped when `vcf` is gzipped (§2.3). It
raises what `parse_vcf` raises (§2.1), and `ValueError` for an `output` that is `vcf` itself. The
two that choose take a `seed`; with `None` they choose the largest shift and the longest runs,
which is what the cases of §4 and §6 give, and with a whole number they draw, with
`random.Random(seed)`, so that the same seed writes the same file.

```python
# tests.respell
def respell_by_decomposing(vcf: Path, reference: IndexedFasta, output: Path, *,
                           ploidy: int = DEFAULT_PLOIDY) -> None: ...
def respell_by_shifting(vcf: Path, reference: IndexedFasta, output: Path, *,
                        seed: int | None = None, ploidy: int = DEFAULT_PLOIDY) -> None: ...
def respell_with_spanning_deletions(vcf: Path, reference: IndexedFasta, output: Path, *,
                                    ploidy: int = DEFAULT_PLOIDY) -> None: ...
def respell_by_merging(vcf: Path, reference: IndexedFasta, output: Path, *,
                       seed: int | None = None, ploidy: int = DEFAULT_PLOIDY) -> None: ...
```

The oracle compares two files of the same samples, in the same order, over one reference (§8).
It returns what it compared, and raises `AssertionError` at the first difference, with the
contig, the window, 1-based, the sample, and what each file gives there; and `ValueError` when the
two headers name different samples.

```python
# tests.haplotype_oracle
@dataclass(frozen=True, slots=True, kw_only=True)
class HaplotypeComparison:
    num_windows: int            # windows compared, over every sample
    num_linked_pairs: int       # pairs of windows compared together (§8)
    num_over_cap: int           # windows and pairs not compared, over ORACLE_CAP in either file

ORACLE_CAP = 4096   # phase assignments; doc/specs/reconciliation.md §4.1's default cap

def check_same_haplotypes(vcf_a: Path, vcf_b: Path, reference: IndexedFasta, *,
                          ploidy: int = DEFAULT_PLOIDY) -> HaplotypeComparison: ...
```

The check of §7.4 runs file A against file B through the stages of pvc_eval, reading,
loci, reconciliation and evaluation (§9.3), and raises
`AssertionError` at the first row out of the allowed set, or the first sample whose result in B
differs from A's, with the locus, the sample, the outcome and both results. `str_regions` is a
`--str-regions` BED; without one, the STR regions come from `segment_contig` over the reference.

```python
# tests.respell_check
def check_respelling_outcomes(vcf_a: Path, vcf_b: Path, reference: IndexedFasta, *,
                              str_regions: Path | None = None,
                              ploidy: int = DEFAULT_PLOIDY) -> None: ...
```

Two helpers are public so that the test modules share them: `check_same_samples`, in
`tests.haplotype_oracle`, which raises `ValueError` when two files' headers do not name the same
samples in the same order, or name none; and `check_outcome_row`, in `tests.respell_check`, the
rule of §9.3 for one row, which the property of `tests/test_reconcile_property.py` also uses.

`tests/test_reconcile_property.py` keeps its property, with file B made by
`respell_by_decomposing` from file A written to its directory, and its hand cases of
`_make_decomposed_lines` become cases of `respell_by_decomposing` (§3), so the tests have one
decomposition, not two.

## 8. The oracle

The oracle rebuilds, for every sample, the haplotype sequences each file gives, and checks they
are equal. It is a short function of the tests, written for this purpose. It reads both files with
`parse_vcf`, which gives the records, their trimmed alleles, their spans and their VCF genotypes,
and imports nothing else of pvc_eval but the reference reader: not `pvc_eval.reconcile`,
`pvc_eval.loci` or `pvc_eval.evaluate`, which are what §7.4 tests. A test checks that its
module imports none of them.

**What is determined.** Where every VCF genotype is phased in one phase set, a sample's two
haplotypes are fixed, and the oracle could compare them base by base over a whole contig. Where
the phase is not known, the file gives a set of possibilities instead, and the oracle compares
the sets. Two kinds of facts are determined:

- **within a window**, the sample's possible genotypes over it: for every phase assignment of the
  blocks with a heterozygous VCF genotype in the window, the multiset of the haplotypes' sequences over the window;
- **between two windows a block joins**, the possible genotypes over both together: for every
  assignment of the blocks with a heterozygous VCF genotype in either, the multiset of the
  haplotypes' pairs of sequences, one over each window. With the windows of a contig numbered
  from the left, a phase set with heterozygous VCF genotypes in windows 1, 4 and 9 joins 1 with 4
  and 4 with 9. A window where it has only homozygous ones joins nothing: phase does not pass
  through it, and joining through it would let a B that moved window 9's record to another phase
  set pass (found by the review).

A **window**, of §1, is built from both files, so that every spelling of an event in either file
falls in one: the spans of all records with an ALT, of both files, sorted, joined while they
overlap or lie next to each other. A record whose only ALT is `*` has its REF as its span;
records with a symbolic allele are not read. Windows, and the phase sets that join them, are
per contig. The span of a record holds every edit it makes and every
placement of its indels, so a haplotype's edits in a window change nothing outside it, and
the sequence of a window is the reference over it with the edits applied. A record with no ALT
belongs to every window its REF overlaps, where only its `.` matters.

Phase between windows is carried only by blocks, and two windows where one block is
heterozygous are linked through every such window of it in between, so consecutive pairs are the
ones compared. A file B that
loses a phase set, or that gives two records of A a phase set they did not share, changes the
possible genotypes over some such pair, of A's blocks or of B's, and the oracle sees it.

**One haplotype over one window**, for one phase assignment, is one of three things:

- **impossible**, when two edits it carries collide as `doc/spec.md` §5.4 says (they replace a
  common position, an insertion lies between two positions one edit replaces, or two insertions
  go in at one point), or when it carries a `*` and no edit it carries replaces, as written, a
  position of that record's REF as written, the stretch its file's REF covers before trimming; the assignment then gives no genotype;
- **missing**, when it is not impossible and is `.` at some record of the window. A `*` is also
  covered, as in reconciliation, by an edit of any ALT of a record at which the haplotype is
  `.`, since the haplotype might carry it; a `.` at a record none of whose ALTs covers the `*`
  leaves the haplotype impossible (found by the review of milestone 2);
- otherwise **a sequence**: the reference over the window with its edits applied in position
  order, an insertion before a replacement that starts at its point.

These are reconciliation's rules, written again more simply: the `*` is covered only by an edit
as written, not by a deletion's other placements. The oracle compares two files with the same
rule, so where a file already in the input relies on a placement, both files are impossible
there alike; every `*` a re-spelling writes stands under its deletion as written.

A sample with a malformed GT, or a GT of another ploidy, at a record of the window is, there,
`malformed GT` or `ploidy mismatch`, which must be the same in both files; a pair with such a
window is compared the same way, by that reason. With ploidy 1 there
are no assignments and a window has one possible genotype.

A window, or a pair, whose blocks give more than `ORACLE_CAP` assignments in either file is not
compared and is counted in `num_over_cap`, unless a sample is `malformed GT` or `ploidy mismatch`
there, which needs no assignment and is compared whatever the cap. The tests of the re-spellings
assert that the count is 0 (§9.2), since a re-spelling that lost a phase set over 14 pieces was
otherwise passed (the review of milestone 2); the check over a real VCF reports it. A message
gives, for each file, how many genotypes it has and those the other lacks, ten at most.

**The cases**, at `check_same_haplotypes`, on reference `GACGTACGTAC` from position 1, or on
`CAGTC` at 9–13 for O3:

| O | file A | file B | result |
|---|---|---|---|
| 1 | `2 A T 0/1`, `4 G C 0/1` | `2 ACG TCG,ACC 1\|2:1` | raises: over the window 2–4, A gives `{ACG, TCC}` and `{ACC, TCG}`, B only `{ACC, TCG}`. A merge of two unphased heterozygous records invents a phase |
| 2 | `2 A T 0\|1:5`, `10 A G 0\|1:5` | `2 A T 0\|1:5`, `10 A G 0\|1:6` | raises: each window is the same, but over the pair 2 and 10, A gives `{(A, A), (T, G)}` alone and B also `{(A, G), (T, A)}` |
| 3 | `9 CAGT C 0\|1:5`, `11 G A 1\|0:5` | `9 CAGT C 0\|1:5`, `11 G A,* 2\|1:5` | raises: in B the `*` is on the haplotype without the deletion, which is impossible, so B gives no genotype over 9–12 and A gives `{CAAT, C}` |
| 4 | `2 AC GT ./1` | `2 A G .\|1:1`, `3 C T .\|1:1` | passes: both give `{., GT}` over 2–3, one window |
| 5 | `2 A T 0/1` | `2 A T 0/1` | passes, with 1 window and no pair |
| 6 | `2 A T 0\|1:5`, `6 A C 1\|1:5`, `10 A G 0\|1:5` | the same, `10 A G 0\|1:6` | raises: phase set 5 is homozygous at 6, so it joins 2 with 10, where A gives `{(A, A), (T, G)}` alone and B also `{(A, G), (T, A)}` |

## 9. How it is verified

Each check is made at the function that shows it, with its expected result as a literal.

### 9.1 The cases

- **At each re-spelling**, over files written to `tmp_path`: every case of §3, §4, §5 and §6,
  with the output's record lines as literals. For each row of §2.2 and each re-spelling, one input
  with that record beside a record the re-spelling rewrites: the row's record is written as read
  and its neighbour is rewritten; for a REF mismatch, `PvcEvalInputError` and
  no output file. The header: the PS line added once, and not added when the input has one. A
  bgzipped input gives an output that starts with the gzip bytes and that `parse_vcf` reads. The
  same seed gives the same bytes, for the shift and the merge.
- **At `check_same_haplotypes`**: O1 to O6; and every case of §3 to §6, input against output,
  passes.
- **At `check_respelling_outcomes`**: every case of §3 to §6, input against output, passes; and
  one case through `segment_contig`, the AT repeat of `tests/test_loci_from_vcf.py` with
  `25 G GAT 0/1` shifted through the repeat, whose single locus holds the STR region and gives
  *same genotype*.
- **`bcftools view`** reads the output of every case, as §2.3 says, when bcftools is on the
  path; the test is skipped, with a message, when it is not.

### 9.2 The properties

The files are generated by hypothesis, small enough to shrink well and dense enough that the
records interact:

- a reference of 12 to 60 bases, made of 2 to 5 pieces, each either 3 to 8 random bases or a
  repeat of a unit of 1 to 3 bases in 3 to 6 copies; in one file of four the contig is exactly
  those bases, so that records reach its ends, and otherwise it is padded with `N`;
- 1 to 6 records on each contig, sorted by POS, each with 1 or 2 ALTs: a SNP, a replacement of 2 to 4 bases, the
  insertion or the deletion of whole units of a repeat at its edge, or a random indel; one record
  in twenty with no ALT, whose GTs are drawn from `0` and `.` alone, and one in twenty with a
  `*` ALT;
- 1 to 3 samples, ploidy 2 in four files of five and 1 otherwise; three samples in four with
  every GT phased in phase set 1, so that phase sets join windows, and the others with each GT
  phased in one draw of two, with a PS of `1`, `2` or none, `.` in one allele of ten, malformed in one GT of fifty, and
  of the other ploidy in one GT of fifty;
- in one file of five, a second contig, so that the contigs' order and new phase sets across
  contigs are exercised.

Every property runs a fixed 300 examples, with `deadline=None`, and the seed of the shift and the
merge is drawn too:

- **the oracle**, for each re-spelling: `check_same_haplotypes` of the file against its
  re-spelling passes;
- **§7.4**, for each re-spelling: `check_respelling_outcomes` passes;
- **bcftools**, for each re-spelling, 50 examples: `bcftools view` reads the output as §2.3 says;
  skipped without bcftools;
- **the re-spellings do something**: of the 300 files of each property, at least one in ten has
  a record its re-spelling rewrites, so that a re-spelling that kept everything would fail. The
  count is made, for each re-spelling, on a fixed set of 1,000 generated files, in a test of its
  own, and the properties stay random: in random runs of 300 the shifted files ranged from 22 to
  77, around the 30 asserted (the review of milestone 4), and a fixed set of 300 left the shift and
  the `*` 1 and 2 files above it after the generator changed in milestone 6, so the set is 1,000. The generator is tuned until each
  re-spelling reaches it, the `*` one by drawing deletions over other records;
- **the oracle sees phase between windows**: of the 300 files, at least one in five has a pair
  of windows a block joins (§8), so that a re-spelling that loses or invents a phase set between
  windows is likely to be caught in one run. With phase drawn as above, one file in fourteen had
  one, and a defect seen only over a pair was missed by a run of 300 about one time in nine (the
  review of milestone 2). Drawing PS 1 for most phased GTs reached one file in seven; the
  samples phased whole reach one in four (25% over 2,000 files). The share varies by a third
  between runs of 300, so it is asserted by a test of the generator alone, which draws the same
  files every run, and the properties stay random.
  Every test that checks a re-spelling with the oracle asserts that no window and no pair was
  over `ORACLE_CAP`.

A failure the properties find, once shrunk, is written into §3 to §6 as a case, and its test is
kept with its literals.

### 9.3 The check of §7.4

`check_respelling_outcomes` takes file A and its re-spelling B through the stages of
`doc/spec.md` §4: `parse_vcf` for each; `build_loci` over both, with segments from
`read_str_regions` over the BED, or from `segment_contig`; `reconcile_loci`; and
`make_outcome_rows`, with `locus_kind` for the kind and no `--strata` class. The sample pairs pair each
sample with itself. It asserts, at every locus and for every pair:

- the outcome is *same genotype*, *reference in both*, *missing in both*, *invalid in both* or
  *too complex* (`doc/spec.md` §7.3);
- *missing in both*, *invalid in both* and *too complex* only where A's own sample result at the
  locus is missing, invalid or too complex, which is what A gives against itself;
- the sample's result in B equals A's: the same possible genotypes with the same indices, or the
  same reason. The re-spellings keep A's blocks and A's sequences, so this holds, and it is
  stronger than the outcome, which would let a B that lost a phase assignment pass, since *same
  genotype* needs only one pair to match (`doc/specs/reconciliation.md` §9).

The property gives the BED of the repeats it planted, which `read_str_regions` reads without a
type column.

`compare()` joins these stages; another session is building it, and until it is on `main` the
check joins them itself, as `tests/test_reconcile_property.py` does. When it arrives, the rows
are taken from `compare()`, and the results, which `compare()` does not yield, still from
`reconcile_loci`.

### 9.4 A real VCF

Each re-spelling of one chromosome of GIAB's HG002 benchmark is checked with the oracle and with
`check_respelling_outcomes`, its seed fixed, marked `slow`. Where the VCF and its reference are
is the owner's to say; the session building `compare()` asks the same, and these tests read them
from where that answer puts them. Until then they are skipped, with a message that names what is
missing.

## 10. How it runs

A re-spelling holds the input's records and lines in memory, which for a GIAB chromosome is some
10^5 records. Whether a record is touched is found by a sweep over the records in order of REF,
comparing each only with those whose REF touches it, since an edit lies within its REF; the
shift's crossed stretch and the records near a `*` are found the same way over the reach. The
square of the records, which `_make_decomposed_lines` took, would take some 40 minutes on a
chromosome: 0.96 s for 2,000 records in the reconciliation work (`doc/reports/reconciliation.md`,
19 October 2026, on the owner's Mac), times 2,500. The oracle holds the records of both files, as
the re-spellings hold a file, and works on one window, or one pair, at a time,
with at most `ORACLE_CAP` assignments for each sample; on 20,000 one-sample records it took
34 MB and 3.3 s (the review of milestone 2). Neither is measured yet;
the real VCF's run is the first measurement.

## 11. Open points

No point is open. **Open 1. GATK's commonest `*`**, answered by the owner on 20 October 2026: keep
the rule. The question was this. GATK writes a deletion `9 CAGT C 0/1` and a SNP under it
`11 G A 0/1` as `11 G A,* 1/2`, both unphased. Unphased, the SNP may share the deletion's
haplotype, where the two collide, so by the rule of 19 October 2026 the `*` re-spelling keeps
both records (T7) and never writes this case. Reconciliation already gives the two spellings the
same result: its case R1 (`doc/specs/reconciliation.md` §4.4, third row) is this one, tested at
`reconcile`, the function that reconciles one locus, since the assignment that puts
them on one haplotype is impossible in both.

- **Keep the rule (recommended).** The `*` re-spelling writes a `*` only where the phase is known
  (§5). The unphased GATK case stays tested at `reconcile` by R1, not through the whole pipeline.
  Nothing more to build.
- **Allow it.** R gets a `*` on the haplotype that gives it REF when D is heterozygous and every
  ALT R gives the sample collides with D's deletion: the only assignments left put the deletion
  on R's REF haplotype. It tests the commonest real case through the pipeline, at the cost of a
  rule that departs from "edits that touch are kept" and that the oracle has to follow too: one
  more condition in §5, its cases, and a property that draws the case.

The owner chose to keep the rule, on 20 October 2026.

## 12. Not in this spec

- `bcftools norm`, which can change a haplotype (`doc/spec.md` §5.2), is not a re-spelling.
- Chaining re-spellings, such as a shift of a decomposition: each reads what another writes, so a
  later plan can add it as a property without changing the program.
- Re-spelling file B of a comparison, or offering re-spellings on the command line: nothing in
  `src/` or `pvc_eval`'s commands, which the owner would decide.
- The comparison of evaluators of `doc/spec.md` §8, which may run on these re-spellings later.
- `doc/spec.md` §7.4 gains a pointer to this spec in a commit of its own once the owner approves
  it.
