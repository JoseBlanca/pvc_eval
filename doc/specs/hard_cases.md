# pvc's hard cases — specification

*Status: approved by the owner on 21 October 2026, with each expected outcome, and with its Open 1
answered the same day: a `*` is measured against its record's REF as written (§6). No test of it
exists yet. It answers `doc/spec.md` Open 3 and replaces the text of its §7.5. It stands on
`doc/spec.md` §5.2 to §5.5, on `doc/specs/loci.md`, `doc/specs/reconciliation.md` and
`doc/specs/compare.md`, and on the oracle of `doc/specs/respell.md` §8, `tests/haplotype_oracle.py`.*

pvc's work on tandem repeats met, in its reports, the disagreements between VCF files that pvc_eval
was made to judge: GIAB writing one haplotype as an indel and a SNP where pvc writes one record over
the repeat, deletions across the edge of a repeat written as overlapping records with a `*` allele,
and HipSTR writing the whole repeat in one record. This spec takes 50 of them, each with the VCF
records of two real files at one locus. They are numbered by kind: S1 to S14 for GIAB's split
records against pvc's repeat record, E1 to E16 for deletions at a repeat's edge, H1 to H15 for
HipSTR's whole-repeat records, and O1 to O5 for other spellings. Each case gives the row of the
table that `compare()` must give for it: the outcome, the locus, its locus alleles and the genotype
of each file. Each becomes a test that runs `compare()` on a few hundred bases of the real
reference.

The owner approved each expected outcome on 21 October 2026, since the check stands on them. 34
of the 50 are *same genotype*: the files spell one haplotype differently and pvc_eval must see
through it. 11 are *different* and 4 *same alleles, other genotype*, where the files truly
disagree. One, E11, is *invalid in B*, file B being pvc's: pvc wrote what it meant to be TR's
haplotypes, and the collision rule the owner chose on 19 October 2026 makes them impossible.

pvc_eval at commit 96b3938 gives 48 of the 50 rows below. It gives *invalid in B* for S11 and E5,
where pvc pads a record with a `*` allele, because it measures the `*` against the record's REF
trimmed. The owner decided on 21 October 2026 that it is measured against the REF as written (§6),
which makes both *same genotype*; `doc/spec.md` §5.4 and the code change first (§7).

Positions are 1-based and inclusive, as in a VCF, unless they are called excerpt positions.

## 1. What a case is

A case holds:

- **The excerpt**: 301 bases of the reference, from 150 bases before a position of the locus to 150
  after, upper case, written into the test as a literal. The test writes it as a FASTA of one
  contig, named as the reference's contig, so that it needs no download. An excerpt position is the
  position minus the offset each case gives. The excerpts are GRCh38's, from the analysis set that
  `PVC_EVAL_GRCH38` names, and for five cases SL4.0, the tomato reference of pvc's tomato runs.
- **The records of file A and of file B**, each as its file wrote it, with its file and line: POS,
  REF, ALT and the sample's GT, with its PS where the file gives one. The test writes them shifted
  to excerpt positions, QUAL 50, FILTER PASS, INFO `.`, for one sample named as in the file. Every
  other field is left out, since none changes an outcome.
- **The haplotypes each file says**, as the locus alleles of the row and the genotype each file is
  given, and, for a file with more than one possible genotype, how many.
- **The expected row**: the locus, its kind and STR regions, the outcome, the locus alleles and the
  reported genotypes.

File A is GIAB's file wherever a case has one, so that *invalid in B* or *only in A* always names
the caller. A case with no GIAB file takes the files in the order the case names them.

### How the expected rows were made

Each expected row was worked out twice, under the rules of `doc/spec.md` at commit 96b3938, and
the two agree on all 50:

1. **From the haplotypes.** For each file, the oracle of `doc/specs/respell.md` §8, which rebuilds
   haplotypes from the records without pvc_eval's reconciliation, gave the sample's possible
   genotypes over the stretches its records span. Joined with the reference between them, these
   are the possible genotypes over the locus. The outcome came from them by the rules of
   `doc/spec.md` §5.5, in a script written for this and apart from `pvc_eval.evaluate`.
2. **From `compare()`**, on the excerpt and the two VCF files, with the default `CompareConfig` and
   `finder_threads=0`.

Then the owner changed the rules for a `*` in a padded record (§6), which change the outcome of S11
and E5 and the locus of E3, E15 and E16, and no other row. The expected rows are those the oracle
and `compare()` both give under those rules, at commit 5590a35 of this branch, where the code
follows them; S11's and E5's were first worked with a copy of the oracle changed so, by this
spec's writer and again by its reviewer.

The oracle and pvc_eval read records by the rules of §5.4, with one difference: the oracle counts a
`*` as covered only by an edit as written, where pvc_eval also counts every placement of a deletion.
No case here turns on it: with coverage as written only, all 50 rows are the same (worked by the
spec's reviewer). So the oracle checks the reconciliation's code, not its rules. Where a rule gives
a result that the file's caller did not mean, both give it, and the case says so (E11).

Every record was read back from its file at the line named and compared with the case: POS, REF, ALT
and the sample's GT, all 50 cases, no difference. An excerpt is a contig of its own, and the STR
finder drops a repeat within 15 bases of a contig's end (`doc/spec.md` §5.1), so an excerpt could
give other STR regions than the whole reference. The STR regions of each excerpt were compared with
those the finder gives on the excerpt widened by 5,000 bases on each side. Over every locus they are
the same, in all 50 cases; they differ only in repeats at an excerpt's edges, more than 100 bases
from any locus.

## 2. The files

| short name | file | what it is |
|---|---|---|
| TR | `pop_var_caller/benchmarks/ssr_hg002/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz` | GIAB's tandem-repeat benchmark v1.0.1 for HG002, the 50,000-record subset pvc's tract reports scored against. Every GT phased, no PS |
| GIAB-HG002, -HG003, -HG004 | `pop_var_caller/benchmarks/giab/all_bench_regions/vcfs/HG00x_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` | GIAB's small-variant benchmark v4.2.1. GTs mostly unphased |
| pvc-ng30 | `pop_var_caller/benchmarks/ssr_hg002/results/ng/HG002_30x.raw.vcf` | pvc on HG002 at 30x, 2 September 2026, the run `tract_genotype_accuracy_2026-09-03.md` scored |
| pvc-main30 | `pop_var_caller-str-group-b/tmp/tr/main/30x/HG002.vcf` | a later pvc run, October 2026 |
| pvc-touch30, pvc-base30 | `pop_var_caller-str-group-b/tmp/tr/{d_touch,d_base}/30x/HG002.sorted.vcf.gz` | pvc with and without the change of 18 October 2026 that widens a repeat record over a deletion touching the repeat (`doc/devel/reports/ng_str_edge_touching_2026-10-18.md`) |
| pvc-m07-30 | `pop_var_caller-str-group-b/tmp/tr/m07/30x/HG002.sorted.vcf.gz` | pvc 0.1.0 with `--defaults` on HG002 at 30x |
| pvc-300x-HG003, -HG004 | `pop_var_caller-coverage-depth-scale/tmp/depth_scale/after2/300x/HG00x/calls.vcf` | pvc on HG003 and HG004 at 300x, 17 October 2026 |
| pvc-explain30-HG004 | `pop_var_caller-str-group-b/tmp/explain/30x.HG004.vcf` | pvc on HG004 at 30x |
| pvc-ng30-HG004 | `pop_var_caller/benchmarks/giab/results/per_sample/30x/ng/HG004.raw.vcf` | an older pvc run on HG004 at 30x |
| HipSTR30 | `pop_var_caller/benchmarks/ssr_hg002/results/hipstr/HG002_30x.str.vcf.gz` | HipSTR 0.6.2 on HG002 at 30x |
| freebayes30 | `pop_var_caller/benchmarks/ssr_hg002/results/freebayes/HG002_30x.fb.vcf.gz` | freebayes on HG002 at 30x |
| GATK-tomato, pvc-tomato | `pop_var_caller-ref-reads-in-repeats/tmp/tomato/{gatk_N40,main40}.vcf` | GATK HaplotypeCaller and pvc, jointly on 40 tomato samples, SL4.0 |

The paths are under `/Users/jose/devel/`. Line numbers count the lines of the file decompressed, its
header included. These files stay in pvc's directories and the tests do not read them.

No human GATK call set on this machine has a `*` allele or two overlapping records: GATK's are
slices of 29 to 78 records, and the runs behind the edge-deletion report were deleted. The
overlapping records of humans here are pvc's own. Since 10 October 2026, when a deletion crosses the
edge of a repeat, pvc widens its record of the repeat over the deleted flank bases, the **widened
repeat record**, and writes a `*` for that haplotype in its record of those flank bases, the **flank
record**, "how GATK writes it" (`doc/devel/ng/spec/str_edge_deletion_alleles.md` §6). GATK's own are
the tomato cases, E12 to E16, which have no GIAB file.

## 3. What the cases show

**Spellings that agree.** pvc_eval gives *same genotype* where the two files write one haplotype in
these ways, each against the other:

- an indel and a SNP as two records against one record over the repeat whose ALT carries both (S7 to
  S10, S13, H3, H5);
- two indels at one position, phased on different haplotypes, against one record of two ALTs (S1,
  H10);
- an indel anchored before the repeat against a record that starts on its first base (S6), or an
  indel at another place in the repeat (E1, H1, H2);
- a record whose REF reaches into the flank to carry a SNP against the SNP written apart (S12, H7,
  H9, H10);
- a deletion across a repeat's edge, written in a widened repeat record beside a flank record with a
  `*` allele, against the deletion written plainly (E2, E3, E6);
- four heterozygous phased records that give one sequence on both haplotypes, against one record
  `1/1` (O1);
- one insertion placed at three points of a compound repeat (O3 to O5).

Where a file leaves the phase open, one of its possible genotypes is enough (S7, S8, S10, E4, E9,
H6, H8, H11).

**The locus takes the whole stretch.** In S3 the two files give the same repeat, A×22, and pvc adds
a SNP on the base after it. The SNP lies next to the STR region, which it therefore touches, and in
the span of file A's insertion, every position the insertion could be written at (`doc/spec.md`
§5.2), which ends on the base after the repeat. So the locus holds it, and the row is *different*:
pvc_eval compares what each file says over the whole locus, not over the repeat alone (`doc/spec.md`
§5.3). pvc's tract scorer, which compares the repeat alone, counts it right.

**A `*` in a padded record.** In S11 and E5 pvc writes a flank record with a `*` allele and with
REF and ALT padded over bases they share. Normalising trims those bases (§5.2): in S11 the REF
shrinks from 24760292–24760298 to 24760294–24760296, in E5 to the anchor base 104588961, and the
deletion pvc meant to cover the `*` reaches neither. Measured against the trimmed REF, as
`doc/spec.md` §5.4 says until now, the `*` is not covered and the sample is *invalid in B*.
Measured against the REF as written, as the owner decided on 21 October 2026, the deletion covers
it and both are *same genotype*. The cases expect the second.

**A collision as written.** In E11 pvc's two records, trimmed from the left first, give edits that
collide as written, and the sample is *invalid in B*. pvc meant TR's haplotypes. The collision rule
is the owner's of 19 October 2026 (`doc/spec.md` §7.4), and this case stands on it.

Of pvc-touch30's 53 records with a `*` allele, compared with TR over the same excerpts, E5 and E11
are the only ones that make the sample invalid under the rules until now, and E11 the only one
under the rule decided.

**Where the files truly disagree.** pvc's old repeat record keeps the wrong base (S4) or puts an
insertion on one haplotype (S14); the touching-deletion change put a `*` on the wrong haplotype (E8)
or a deletion twice (E10); HipSTR carries a flank base no other file has (H12) or calls a
heterozygote homozygous (H13). Where the two of GIAB's files disagree with each other, the cases
show both against the same caller: E6 and E7, H14 and H15. A haplotype on one copy in one file and
on both in the other is *same alleles, other genotype* (E13, E14, H13, H15).

**HipSTR's phase.** HipSTR writes every GT with `|` and no PS, yet it phases no record with another:
its FORMAT field PQ, the probability it gives the phase, is 0.50 on every heterozygote, and DSNP,
the reads it phased by SNPs, is 0 everywhere. pvc_eval reads them as one phase set per sample and
contig (§5.4), so a locus holding two HipSTR records would be compared under a phase HipSTR never
claimed. None of the 15 HipSTR cases holds two HipSTR records; the files hold one record per repeat.

## 4. The cases

Each table gives a case's locus, its two files, the outcome and the genotype each file is given, as
indices into the locus alleles, which §4.5 lists with each case's records and source. A case's locus
is its only row: no case here gives more than one locus. The kind of every locus is STR but in S13,
S14, H5, H6 and H9, where the STR finder admits no repeat that the records touch, and the kind is
complex.

### 4.1 GIAB's split records against pvc's repeat record

The kind `doc/spec.md` §7.5 names first. Sources:
`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §3 "The genotype-coding trap", §3.1,
§3.4b and §3.4c (S1 to S6); `giab_unexplained_fn_2026-10-03.md` §4, "The 'same variant written
another way' cases" (S7 to S11); `doc/devel/ng/spec/str_bundle_loci.md` §10 step 4 (S12);
`doc/devel/reports/reviews/ng_indel_genotypes_vs_giab_2026-09-11.md`, "Two things it cost" (S13,
S14). All in `/Users/jose/devel/pop_var_caller/`.

| case | locus | file A | file B | outcome | A, B |
|---|---|---|---|---|---|
| S1 | chr1:4416244–4416275 | TR | pvc-ng30 | same genotype | 1/2, 1/2 |
| S2 | chr1:150329032–150329048 | TR | pvc-ng30 | same genotype | 1/1, 1/1 |
| S3 | chr1:9955403–9955423 | TR | pvc-ng30 | different | 1/1, 1/2 |
| S4 | chr1:14722151–14722179 | TR | pvc-ng30 | different | 0/1, 0/2 |
| S5 | chr1:14722151–14722179 | TR | pvc-main30 | same genotype | 0/1, 0/1 |
| S6 | chr1:69233430–69233465 | TR | pvc-ng30 | same genotype | 1/1, 1/1 |
| S7 | chr11:123865140–123865192 | GIAB-HG003 | pvc-300x-HG003 | same genotype | 1/2, 1/2 |
| S8 | chr14:90749397–90749422 | GIAB-HG004 | pvc-300x-HG004 | same genotype | 3/4, 3/4 |
| S9 | chr10:4853442–4853458 | GIAB-HG004 | pvc-300x-HG004 | same genotype | 1/2, 1/2 |
| S10 | chr20:48207407–48207436 | GIAB-HG004 | pvc-300x-HG004 | same genotype | 3/4, 3/4 |
| S11 | chr8:24760287–24760317 | GIAB-HG004 | pvc-300x-HG004 | same genotype | 2/3, 2/3 |
| S12 | chr4:28575007–28575040 | TR | pvc-main30 | same genotype | 1/2, 1/2 |
| S13 | chr9:11024723–11024725 | GIAB-HG004 | pvc-explain30-HG004 | same genotype | 1/1, 1/1 |
| S14 | chr9:11024723–11024725 | GIAB-HG004 | pvc-ng30-HG004 | different | 1/1, 1/2 |

### 4.2 Deletions at a repeat's edge, written with overlapping records

`doc/spec.md` §7.5's "GATK's overlapping records at STR edge deletions". E1 is the locus
`doc/devel/reports/ng_str_edge_touching_2026-10-18.md` names in "What it does". E2 to E11 were found
among pvc-touch30's 53 records with a `*` allele, which `str_edge_deletion_alleles.md` §6 makes; no
report names them. E12 to E16 were found among GATK's 642 `*` records in GATK-tomato, which has no
GIAB file: file A is GATK, file B pvc.

| case | locus | file A | file B | outcome | A, B |
|---|---|---|---|---|---|
| E1 | chr5:7613158–7613171 | TR | pvc-touch30 | same genotype | 0/1, 0/1 |
| E2 | chr13:31685219–31685246 | TR | pvc-touch30 | same genotype | 1/1, 1/1 |
| E3 | chr20:3029115–3029159 | TR | pvc-touch30 | same genotype | 1/1, 1/1 |
| E4 | chr20:3029115–3029159 | TR | pvc-base30 | same genotype | 1/1, 1/1 |
| E5 | chr8:104588961–104588997 | TR | pvc-touch30 | same genotype | 1/2, 1/2 |
| E6 | chr12:13880994–13881016 | TR | pvc-touch30 | same genotype | 1/2, 1/2 |
| E7 | chr12:13880994–13881016 | GIAB-HG002 | pvc-touch30 | different | 0/1, 2/4 |
| E8 | chr4:98053533–98053546 | TR | pvc-touch30 | different | 0/1, 1/2 |
| E9 | chr4:98053533–98053546 | TR | pvc-base30 | same genotype | 0/1, 0/1 |
| E10 | chr16:5228832–5228881 | TR | pvc-touch30 | different | 1/1, 1/2 |
| E11 | chr1:64356008–64356056 | TR | pvc-touch30 | invalid in B | 1/2, . |
| E12 | SL4.0ch05:576140–576154 | GATK-tomato | pvc-tomato | same genotype | 1/1, 1/1 |
| E13 | SL4.0ch05:576140–576154 | GATK-tomato | pvc-tomato | same alleles, other genotype | 1/1, 0/1 |
| E14 | SL4.0ch05:576140–576154 | GATK-tomato | pvc-tomato | same alleles, other genotype | 1/1, 0/1 |
| E15 | SL4.0ch05:58429419–58429436 | GATK-tomato | pvc-tomato | different | 1/1, 1/2 |
| E16 | SL4.0ch05:58429419–58429436 | GATK-tomato | pvc-tomato | different | 1/1, 2/2 |

### 4.3 HipSTR's whole-repeat records

No report of pvc names a locus where HipSTR and another file spell one haplotype differently:
`doc/devel/reports/ssr_hipstr_concordance_2026-06-25.md` is on tomato, with no GIAB file, and
`tract_genotype_investigation/hipstr.md` gives HipSTR's scores and its traits (records reaching into
the flank, `|` without phasing). These cases were found by comparing HipSTR30's 3,163 records with
an ALT against TR, GIAB-HG002 and pvc-m07-30, by the haplotypes each gives over the repeat; 2,447 of
the 3,163 give the same as all three files. File A is GIAB's where the case has one, so HipSTR is
file B but in H2 and H11.

| case | locus | file A | file B | outcome | A, B |
|---|---|---|---|---|---|
| H1 | chr1:113339221–113339260 | TR | HipSTR30 | same genotype | 0/1, 0/1 |
| H2 | chr1:113339221–113339260 | HipSTR30 | pvc-m07-30 | same genotype | 0/1, 0/1 |
| H3 | chr1:72941245–72941319 | TR | HipSTR30 | same genotype | 1/1, 1/1 |
| H4 | chr1:72941245–72941319 | TR | GIAB-HG002 | same genotype | 1/1, 1/1 |
| H5 | chr1:80541526–80541557 | TR | HipSTR30 | same genotype | 0/1, 0/1 |
| H6 | chr1:80541526–80541557 | GIAB-HG002 | HipSTR30 | same genotype | 0/1, 0/1 |
| H7 | chr1:14451931–14451957 | TR | HipSTR30 | same genotype | 0/1, 0/1 |
| H8 | chr1:14451930–14451957 | GIAB-HG002 | HipSTR30 | same genotype | 0/1, 0/1 |
| H9 | chr1:61522963–61522979 | TR | HipSTR30 | same genotype | 0/1, 0/1 |
| H10 | chr1:102475773–102475817 | TR | HipSTR30 | same genotype | 1/2, 1/2 |
| H11 | chr1:102475773–102475817 | HipSTR30 | pvc-m07-30 | same genotype | 1/2, 1/2 |
| H12 | chr1:172778363–172778407 | TR | HipSTR30 | different | 1/2, 2/3 |
| H13 | chr10:12391976–12392017 | TR | HipSTR30 | same alleles, other genotype | 0/1, 1/1 |
| H14 | chr1:93319248–93319304 | TR | HipSTR30 | different | 0/1, 2/2 |
| H15 | chr1:93319249–93319304 | GIAB-HG002 | HipSTR30 | same alleles, other genotype | 0/1, 1/1 |

### 4.4 Other spellings

O1 and O2 are in `doc/devel/reports/tract_genotype_investigation/unseen.md`, among "spelling and
rotation disagreements at the same length"; O3 to O5 in its §(b).

| case | locus | file A | file B | outcome | A, B |
|---|---|---|---|---|---|
| O1 | chr6:41399584–41399642 | TR | pvc-touch30 | same genotype | 1/1, 1/1 |
| O2 | chr6:41399584–41399642 | TR | freebayes30 | different | 1/1, 1/6 |
| O3 | chr3:37126860–37126872 | TR | pvc-touch30 | same genotype | 1/1, 1/1 |
| O4 | chr3:37126860–37126872 | TR | HipSTR30 | same genotype | 1/1, 1/1 |
| O5 | chr3:37126860–37126872 | TR | freebayes30 | same genotype | 1/1, 1/1 |

### 4.5 The cases in full

Each case gives its excerpt and offset; its records, each with the short name of its file and its
line (§2) and as the file wrote it; the locus with its kind and STR regions; the locus alleles,
numbered, which are the haplotypes the two files say; and the genotype of each file, the number of
possible genotypes of each, and the outcome. A sample with more than one possible genotype has its
others among the alleles listed.

**S1.** TR writes a heterozygote of two insertions as two phased records at one position; pvc one
record of two ALTs over the repeat. Two insertions at one point on different haplotypes do not
collide.

```
excerpt  chr1:4416100-4416400, 301 bases; excerpt position = position - 4416099; sample HG002
A  TR:180  4416244 C CGTGT 0|1
A  TR:181  4416244 C CGTGTGT 1|0
B  pvc-ng30:2677  4416245 GTGTGTGTGTGTGTGTGT GTGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGT 1/2
locus    chr1:4416244-4416275 (excerpt 145-176), kind STR; STR regions: 4416245-4416275 str-bundle
  0  CGTGTGTGTGTGTGTGTGTGAAAGAGAGAGAG
  1  CGTGTGTGTGTGTGTGTGTGTGTGTGAAAGAGAGAGAG
  2  CGTGTGTGTGTGTGTGTGTGTGTGAAAGAGAGAGAG
possible genotypes  A: 1/2   B: 1/2
row  A 1/2  B 1/2  possible 1, 1  records 2, 1  -> same genotype
```

**S2.** TR writes the stretch as four records, an insertion followed by SNPs on the next two
positions; pvc as a deletion and a repeat record. The insertion goes in between 150329035 and
150329036 and the SNP replaces 150329036, which is no collision.

```
excerpt  chr1:150328890-150329190, 301 bases; excerpt position = position - 150328889; sample HG002
A  TR:1583  150329033 T A 1|1
A  TR:1584  150329035 A ACTTTTT 1|1
A  TR:1585  150329036 A T 1|1
A  TR:1586  150329037 C T 1|1
B  pvc-ng30:3583  150329032 GT G 1/1
B  pvc-ng30:3584  150329038 TTTTTTTTTT TTTTTTTTTTTTTTTTT 1/1
locus    chr1:150329032-150329048 (excerpt 143-159), kind STR; STR regions: 150329038-150329047 str T
  0  GTAAACTTTTTTTTTTG
  1  GAAACTTTTTTTTTTTTTTTTTG
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 4, 2  -> same genotype
```

**S3.** Over the repeat both files give A×22. pvc's SNP on 9955423, the base after the repeat,
touches the STR region and lies in the insertion's span, so the locus takes it in, and B carries an
allele A does not.

```
excerpt  chr1:9955260-9955560, 301 bases; excerpt position = position - 9955259; sample HG002
A  TR:283  9955403 C CAAA 1|1
B  pvc-ng30:2747  9955404 AAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAA 1/1
B  pvc-ng30:2748  9955423 G A 0/1
locus    chr1:9955403-9955423 (excerpt 144-164), kind STR; STR regions: 9955404-9955422 str A
  0  CAAAAAAAAAAAAAAAAAAAG
  1  CAAAAAAAAAAAAAAAAAAAAAAG
  2  CAAAAAAAAAAAAAAAAAAAAAAA
possible genotypes  A: 1/1   B: 1/2
row  A 1/1  B 1/2  possible 1, 1  records 1, 2  -> different
```

**S4.** Same length, other sequence: TR's haplotype deletes the repeat's last A and the G after it,
with a SNP inside; pvc's old record deletes two A's inside and keeps the G.

```
excerpt  chr1:14722006-14722306, 301 bases; excerpt position = position - 14722005; sample HG002
A  TR:317  14722158 A T 0|1
A  TR:318  14722161 AAG A 0|1
B  pvc-ng30:2776  14722151 AAAAAAAAAAAA AAAAAAATAA 0/1
locus    chr1:14722151-14722179 (excerpt 146-174), kind STR; STR regions: 14722151-14722179 str-bundle
  0  AAAAAAAAAAAAGAATAAAAAATAAAAAA
  1  AAAAAAATAAAAATAAAAAATAAAAAA
  2  AAAAAAATAAGAATAAAAAATAAAAAA
possible genotypes  A: 0/1   B: 0/2
row  A 0/1  B 0/2  possible 1, 1  records 2, 1  -> different
```

**S5.** The same locus with pvc's later record, widened over the G, which spells TR's sequence.

```
excerpt  chr1:14722006-14722306, 301 bases; excerpt position = position - 14722005; sample HG002
A  TR:317  14722158 A T 0|1
A  TR:318  14722161 AAG A 0|1
B  pvc-main30:2770  14722151 AAAAAAAAAAAAG AAAAAAATAAA 0/1
locus    chr1:14722151-14722179 (excerpt 146-174), kind STR; STR regions: 14722151-14722179 str-bundle
  0  AAAAAAAAAAAAGAATAAAAAATAAAAAA
  1  AAAAAAATAAAAATAAAAAATAAAAAA
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 2, 1  -> same genotype
```

**S6.** TR anchors the deletion on the base before the repeat, pvc writes the whole repeat from its
first base.

```
excerpt  chr1:69233295-69233595, 301 bases; excerpt position = position - 69233294; sample HG002
A  TR:966  69233430 TATAATAATA T 1|1
B  pvc-ng30:3258  69233431 ATAATAATAATAATAATAATAATAATAATAATA ATAATAATAATAATAATAATAATA 1/1
locus    chr1:69233430-69233465 (excerpt 136-171), kind STR; STR regions: 69233431-69233463 str ATA
  0  TATAATAATAATAATAATAATAATAATAATAATAAT
  1  TATAATAATAATAATAATAATAATAAT
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
```

**S7.** HG003. GIAB writes a record of two ALTs and an unphased SNP; pvc one record whose ALT
carries the SNP. GIAB's two phase assignments give two possible genotypes, and one of them is pvc's.

```
excerpt  chr11:123865020-123865320, 301 bases; excerpt position = position - 123865019; sample HG003
A  GIAB-HG003:2756918  123865154 CGTGTGTGT CGTGTGT,C 2/1
A  GIAB-HG003:2756919  123865171 G A 0/1
B  pvc-300x-HG003:3112  123865155 GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTGTGTGTATGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 1/2
locus    chr11:123865140-123865192 (excerpt 121-173), kind STR; STR regions: 123865140-123865192 str-bundle
  0  AAAAATTCTGTAACCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  1  AAAAATTCTGTAACCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  2  AAAAATTCTGTAACCGTGTGTGTGTGTGTATGTGTGTGTGTGTGTGTGTGT
  3  AAAAATTCTGTAACCGTGTGTGTATGTGTGTGTGTGTGTGTGTGT
  4  AAAAATTCTGTAACCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
possible genotypes  A: 1/2, 3/4   B: 1/2
row  A 1/2  B 1/2  possible 2, 1  records 2, 1  -> same genotype
```

**S8.** HG004. The pattern of S7 on a T repeat, the SNP on the repeat's last base.

```
excerpt  chr14:90749262-90749562, 301 bases; excerpt position = position - 90749261; sample HG004
A  GIAB-HG004:3265346  90749402 CTT CT,C 2/1
A  GIAB-HG004:3265347  90749422 T G 0/1
B  pvc-300x-HG004:3106  90749403 TTTTTTTTTTTTTTTTTTTT TTTTTTTTTTTTTTTTTG,TTTTTTTTTTTTTTTTTTT 1/2
locus    chr14:90749397-90749422 (excerpt 136-161), kind STR; STR regions: 90749397-90749422 str-bundle
  0  TTTTTCTTTTTTTTTTTTTTTTTTTT
  1  TTTTTCTTTTTTTTTTTTTTTTTT
  2  TTTTTCTTTTTTTTTTTTTTTTTTG
  3  TTTTTCTTTTTTTTTTTTTTTTTG
  4  TTTTTCTTTTTTTTTTTTTTTTTTT
possible genotypes  A: 1/2, 3/4   B: 3/4
row  A 3/4  B 3/4  possible 2, 1  records 2, 1  -> same genotype
```

**S9.** HG004. GIAB's `0/1` deletion and `1/1` SNP are pvc's one `1/2` record.

```
excerpt  chr10:4853300-4853600, 301 bases; excerpt position = position - 4853299; sample HG004
A  GIAB-HG004:2407486  4853442 CA C 0/1
A  GIAB-HG004:2407487  4853458 A C 1/1
B  pvc-300x-HG004:3028  4853443 AAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAC,AAAAAAAAAAAAAAC 1/2
locus    chr10:4853442-4853458 (excerpt 143-159), kind STR; STR regions: 4853443-4853458 str A
  0  CAAAAAAAAAAAAAAAA
  1  CAAAAAAAAAAAAAAAC
  2  CAAAAAAAAAAAAAAC
possible genotypes  A: 1/2   B: 1/2
row  A 1/2  B 1/2  possible 1, 1  records 2, 1  -> same genotype
```

**S10.** HG004. The pattern of S7 on a CA repeat.

```
excerpt  chr20:48207270-48207570, 301 bases; excerpt position = position - 48207269; sample HG004
A  GIAB-HG004:3904085  48207407 GCACA GCACACA,G 2/1
A  GIAB-HG004:3904086  48207435 A G 0/1
B  pvc-300x-HG004:3235  48207408 CACACACACACACACACACACACACACA CACACACACACACACACACACACACACACA,CACACACACACACACACACACACG 1/2
locus    chr20:48207407-48207436 (excerpt 138-167), kind STR; STR regions: 48207408-48207435 str CA
  0  GCACACACACACACACACACACACACACAT
  1  GCACACACACACACACACACACACAT
  2  GCACACACACACACACACACACACACACACGT
  3  GCACACACACACACACACACACACGT
  4  GCACACACACACACACACACACACACACACAT
possible genotypes  A: 1/2, 3/4   B: 3/4
row  A 3/4  B 3/4  possible 2, 1  records 2, 1  -> same genotype
```

**S11.** HG004. pvc writes the 2-base deletion as `TATACAC → TATAC,*` and the 4-base deletion in the
repeat as a record of its own. The `*` record's REF as written is 24760292–24760298, and the 4-base
deletion, as written at 24760298–24760301 or slid to 24760297–24760300, deletes a position of it,
so it covers the `*`: pvc's one possible genotype is GIAB's second. Trimmed, the REF would be
24760294–24760296, which no placement deletes, and the sample would be invalid (§6).

```
excerpt  chr8:24760150-24760450, 301 bases; excerpt position = position - 24760149; sample HG004
A  GIAB-HG004:2066390  24760294 TAC T 0/1
A  GIAB-HG004:2066391  24760296 CACAT C 0/1
B  pvc-300x-HG004:2886  24760292 TATACAC TATAC,* 1/2
B  pvc-300x-HG004:2887  24760297 ACATATATATATATATATAT ATATATATATATATAT 0/1
locus    chr8:24760287-24760317 (excerpt 138-168), kind STR; STR regions: 24760287-24760317 str-bundle
  0  AAAAATATACACATATATATATATATATATA
  1  AAAAATATATATATATATATATATA
  2  AAAAATATACATATATATATATATATA
  3  AAAAATATACATATATATATATATATATA
possible genotypes  A: 0/1, 2/3   B: 2/3
row  A 2/3  B 2/3  possible 2, 1  records 2, 2  -> same genotype
```

**S12.** TR writes three phased records across an STR bundle, pvc one record over the bundle,
widened over the flank SNP at 28575040.

```
excerpt  chr4:28574875-28575175, 301 bases; excerpt position = position - 28574874; sample HG002
A  TR:7780  28575011 ACC A 1|0
A  TR:7781  28575017 CCCCCCA C 0|1
A  TR:7782  28575040 G C 1|0
B  pvc-main30:7107  28575007 ACACACCCCCCCCCCCACACACACACACACACAG ACACACCCCCCCACACACACACACACAG,ACACACCCCCCCCCACACACACACACACACAC 1/2
locus    chr4:28575007-28575040 (excerpt 133-166), kind STR; STR regions: 28575007-28575039 str-bundle
  0  ACACACCCCCCCCCCCACACACACACACACACAG
  1  ACACACCCCCCCCCACACACACACACACACAC
  2  ACACACCCCCCCACACACACACACACAG
possible genotypes  A: 1/2   B: 1/2
row  A 1/2  B 1/2  possible 1, 1  records 3, 1  -> same genotype
```

**S13.** HG004. GIAB writes an insertion and a SNP two positions apart, pvc one complex record.

```
excerpt  chr9:11024574-11024874, 301 bases; excerpt position = position - 11024573; sample HG004
A  GIAB-HG004:2245326  11024723 C CAG 1/1
A  GIAB-HG004:2245327  11024725 C G 1/1
B  pvc-explain30-HG004:2917  11024723 CAC CAGAG 1/1
locus    chr9:11024723-11024725 (excerpt 150-152), kind complex; STR regions: none
  0  CAC
  1  CAGAG
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 2, 1  -> same genotype
```

**S14.** HG004. pvc's older run puts the insertion on one haplotype only.

```
excerpt  chr9:11024574-11024874, 301 bases; excerpt position = position - 11024573; sample HG004
A  GIAB-HG004:2245326  11024723 C CAG 1/1
A  GIAB-HG004:2245327  11024725 C G 1/1
B  pvc-ng30-HG004:2916  11024723 C CAG,G 0/1
B  pvc-ng30-HG004:2917  11024725 C G 1/1
locus    chr9:11024723-11024725 (excerpt 150-152), kind complex; STR regions: none
  0  CAC
  1  CAGAG
  2  CAG
possible genotypes  A: 1/1   B: 1/2
row  A 1/1  B 1/2  possible 1, 1  records 2, 2  -> different
```

**E1.** The locus of pvc's report on deletions that touch an STR segment. TR anchors the deletion
before the repeat; pvc's widened repeat record starts one base later.

```
excerpt  chr5:7613015-7613315, 301 bases; excerpt position = position - 7613014; sample HG002
A  TR:10175  7613158 AAT A 1|0
B  pvc-touch30:8329  7613159 ATAAAAAAAAAAA AAAAAAAAAAA 0/1
locus    chr5:7613158-7613171 (excerpt 144-157), kind STR; STR regions: 7613161-7613171 str A
  0  AATAAAAAAAAAAA
  1  AAAAAAAAAAAA
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 1, 1  -> same genotype
```

**E2.** pvc writes the deletion across the repeat's left edge in a widened repeat record and turns
the flank record into `AC → AG,*`. Trimmed, the `*` record's REF is 31685224, which STR allele 1
replaces, so it covers the `*`.

```
excerpt  chr13:31685085-31685385, 301 bases; excerpt position = position - 31685084; sample HG002
A  TR:26438  31685223 ACTAT A 1|1
A  TR:26439  31685228 A G 1|1
B  pvc-touch30:17007  31685223 AC AG,* 1/2
B  pvc-touch30:17008  31685224 CTATATATATATATATATATATA GTATATATATATATATATA,CTATATATATATATATATA 1/2
locus    chr13:31685219-31685246 (excerpt 135-162), kind STR; STR regions: 31685219-31685246 str-bundle
  0  AAAAACTATATATATATATATATATATA
  1  AAAAAGTATATATATATATATATA
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 2, 2  -> same genotype
```

**E3.** pvc's run with the touching-deletion change: `GA → GG,*` beside a repeat record over the
bundle.

```
excerpt  chr20:3028985-3029285, 301 bases; excerpt position = position - 3028984; sample HG002
A  TR:34237  3029116 A G 1|1
A  TR:34238  3029132 T G 1|1
A  TR:34239  3029134 G T 1|1
A  TR:34240  3029135 T TTTG 1|1
B  pvc-touch30:21736  3029115 GA GG,* 1/2
B  pvc-touch30:21737  3029116 ATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG,ATTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG 1/2
locus    chr20:3029115-3029159 (excerpt 131-175), kind STR; STR regions: 3029117-3029159 str-bundle
  0  GATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG
  1  GGTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 4, 2  -> same genotype
```

**E4.** The same locus in pvc's run before that change, `GA → G,GG` unphased beside the repeat
record: two possible genotypes, one of them TR's.

```
excerpt  chr20:3028985-3029285, 301 bases; excerpt position = position - 3028984; sample HG002
A  TR:34237  3029116 A G 1|1
A  TR:34238  3029132 T G 1|1
A  TR:34239  3029134 G T 1|1
A  TR:34240  3029135 T TTTG 1|1
B  pvc-base30:21766  3029115 GA G,GG 1/2
B  pvc-base30:21767  3029117 TTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG,TTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG 1/2
locus    chr20:3029115-3029159 (excerpt 131-175), kind STR; STR regions: 3029117-3029159 str-bundle
  0  GATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG
  1  GGTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG
  2  GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG
  3  GGGTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG
possible genotypes  A: 1/1   B: 1/1, 2/3
row  A 1/1  B 1/1  possible 1, 2  records 4, 2  -> same genotype
```

**E5.** pvc writes the 18-base deletion in the repeat record and turns the insertion record into
`ACGACGCCGA → ACGACGCCGACGCCGA,*`. The `*` record's REF as written is 104588961–104588970, which
the 18-base deletion, at 104588966–104588983 as written, deletes in part, so it covers the `*` and
pvc gives TR's haplotypes. Trimmed, the record would be `104588961 A → ACGACGC,*`, whose REF is
its anchor base alone, which no placement deletes, and the sample would be invalid (§6).

```
excerpt  chr8:104588825-104589125, 301 bases; excerpt position = position - 104588824; sample HG002
A  TR:17773  104588961 A ACGACGC 0|1
A  TR:17774  104588964 ACGCCGACGCCGCCGCCGC A 1|0
B  pvc-touch30:12607  104588961 ACGACGCCGA ACGACGCCGACGCCGA,* 1/2
B  pvc-touch30:12608  104588965 CGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC CGCCGCCGCCGCCGC 0/1
locus    chr8:104588961-104588997 (excerpt 137-173), kind STR; STR regions: 104588971-104588997 str CGC
  0  ACGACGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC
  1  ACGACGCCGCCGCCGCCGC
  2  ACGACGCCGACGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC
possible genotypes  A: 1/2   B: 1/2
row  A 1/2  B 1/2  possible 1, 1  records 2, 2  -> same genotype
```

**E6.** pvc writes the deletion across the repeat's right edge in the repeat record and the flank
SNP as `C → T,*`; TR writes three phased records.

```
excerpt  chr12:13880855-13881155, 301 bases; excerpt position = position - 13880854; sample HG002
A  TR:23634  13880994 T TTGTGTG 0|1
A  TR:23635  13881009 TGTGTGC T 1|0
A  TR:23636  13881015 C T 0|1
B  pvc-touch30:15740  13880995 TGTGTGTGTGTGTGTGTGTGC TGTGTGTGTGTGTGT,TGTGTGTGTGTGTGTGTGTGTGTGTGC 1/2
B  pvc-touch30:15741  13881015 C T,* 1/2
locus    chr12:13880994-13881016 (excerpt 140-162), kind STR; STR regions: 13880995-13881014 str TG
  0  TTGTGTGTGTGTGTGTGTGTGCG
  1  TTGTGTGTGTGTGTGTG
  2  TTGTGTGTGTGTGTGTGTGTGTGTGTGTG
possible genotypes  A: 1/2   B: 1/2
row  A 1/2  B 1/2  possible 1, 1  records 3, 2  -> same genotype
```

**E7.** The same pvc records against GIAB v4.2.1, which lacks the SNP at 13881015. TR and GIAB
v4.2.1 disagree here, and the locus is outside GIAB v4.2.1's benchmark regions.

```
excerpt  chr12:13880855-13881155, 301 bases; excerpt position = position - 13880854; sample HG002
A  GIAB-HG002:2834424  13880994 T TTGTGTG 0/1
A  GIAB-HG002:2834425  13881009 TGTGTGC T 0/1
B  pvc-touch30:15740  13880995 TGTGTGTGTGTGTGTGTGTGC TGTGTGTGTGTGTGT,TGTGTGTGTGTGTGTGTGTGTGTGTGC 1/2
B  pvc-touch30:15741  13881015 C T,* 1/2
locus    chr12:13880994-13881016 (excerpt 140-162), kind STR; STR regions: 13880995-13881014 str TG
  0  TTGTGTGTGTGTGTGTGTGTGCG
  1  TTGTGTGTGTGTGTGTGTGTGTG
  2  TTGTGTGTGTGTGTGTG
  3  TTGTGTGTGTGTGTGTGTGTGTGTGTGCG
  4  TTGTGTGTGTGTGTGTGTGTGTGTGTGTG
possible genotypes  A: 0/1, 2/3   B: 2/4
row  A 0/1  B 2/4  possible 2, 1  records 2, 2  -> different
```

**E8.** pvc's widened repeat allele already spells the deletion and the SNP, and the `*` went to the
flank record's other haplotype, so pvc deletes the three bases on both.

```
excerpt  chr4:98053390-98053690, 301 bases; excerpt position = position - 98053389; sample HG002
A  TR:8408  98053533 TAAC T 1|0
A  TR:8409  98053539 A T 1|0
B  pvc-touch30:7554  98053533 TAAC T,* 1/2
B  pvc-touch30:7555  98053536 CAAAAAAAAAA TAAAAAAA 0/1
locus    chr4:98053533-98053546 (excerpt 144-157), kind STR; STR regions: 98053537-98053546 str A
  0  TAACAAAAAAAAAA
  1  TAATAAAAAAA
  2  TAAAAAAAAAA
possible genotypes  A: 0/1   B: 1/2
row  A 0/1  B 1/2  possible 1, 1  records 2, 2  -> different
```

**E9.** The same locus in pvc's run before the change, unphased: one of its two possible genotypes
is TR's.

```
excerpt  chr4:98053390-98053690, 301 bases; excerpt position = position - 98053389; sample HG002
A  TR:8408  98053533 TAAC T 1|0
A  TR:8409  98053539 A T 1|0
B  pvc-base30:7565  98053533 TAAC T 0/1
B  pvc-base30:7566  98053537 AAAAAAAAAA AATAAAAAAA 0/1
locus    chr4:98053533-98053546 (excerpt 144-157), kind STR; STR regions: 98053537-98053546 str A
  0  TAACAAAAAAAAAA
  1  TAATAAAAAAA
  2  TAACAATAAAAAAA
  3  TAAAAAAAAAA
possible genotypes  A: 0/1   B: 0/1, 2/3
row  A 0/1  B 0/1  possible 1, 2  records 2, 2  -> same genotype
```

**E10.** pvc carries a 7-base and a 20-base deletion on one haplotype, where TR has the 20-base
deletion on both.

```
excerpt  chr16:5228700-5229000, 301 bases; excerpt position = position - 5228699; sample HG002
A  TR:29998  5228832 AAGAGAGCGAGAGAGAGAGAG A 1|1
A  TR:29999  5228859 A T 1|1
B  pvc-touch30:18978  5228832 AAGAGAGC A,* 1/2
B  pvc-touch30:18979  5228839 CGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA CGAGAGAGAGAGAGAGAGAGAGA,TGAGAGAGAGAGAGAGAGAGAGA 1/2
locus    chr16:5228832-5228881 (excerpt 133-182), kind STR; STR regions: 5228840-5228881 str GA
  0  AAGAGAGCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA
  1  AAGAGAGTGAGAGAGAGAGAGAGAGAGAGA
  2  AGAGAGAGAGAGAGAGAGAGAGA
possible genotypes  A: 1/1   B: 1/2
row  A 1/1  B 1/2  possible 1, 1  records 2, 2  -> different
```

**E11.** pvc means TR's haplotypes: STR allele 1 (−22) with the flank deletion `TAA → T`, and STR
allele 2 (−14) with the `*`. Trimmed from the left first, allele 1 deletes 64356032–64356053 and
allele 2 64356041–64356054, and `TAA → T` deletes 64356053–64356054: it collides with both as
written, and the sample is invalid. The collision rule judges edits as written (decided by the owner
on 19 October 2026).

```
excerpt  chr1:64355880-64356180, 301 bases; excerpt position = position - 64355879; sample HG002
A  TR:932  64356030 TATATATATATATATATATATATAA T 0|1
A  TR:933  64356040 TATATATATATATAA T 1|0
B  pvc-touch30:3169  64356009 ATATATATATATATATATATATATATATATATATATATATATATAA ATATATATATATATATATATATAA,ATATATATATATATATATATATATATATATAT 1/2
B  pvc-touch30:3170  64356052 TAA T,* 1/2
locus    chr1:64356008-64356056 (excerpt 129-177), kind STR; STR regions: 64356009-64356052 str AT
  0  CATATATATATATATATATATATATATATATATATATATATATATAAAA
  1  CATATATATATATATATATATATATATATATATAA
  2  CATATATATATATATATATATATAA
possible genotypes  A: 1/2   B: none
row  A 1/2  B .  possible 1, 0  records 2, 2  -> invalid in B, contradicting records
```

**E12.** Tomato, sample SRS3394689. GATK writes the deletion of the flank A and a record `A → T,*`
over it; pvc a repeat record and a record of three ALTs and `*`. Both delete the A on both
haplotypes.

```
excerpt  SL4.0ch05:576000-576300, 301 bases; excerpt position = position - 575999; sample SRS3394689
A  GATK-tomato:65098  576151 TA T 1/1
A  GATK-tomato:65099  576152 A T,* 2/2
B  pvc-tomato:71841  576140 TTTTTTTTTTTTA TTTTTTTTTTT 0/0
B  pvc-tomato:71842  576152 ATA TA,TTA,* 1/1
locus    SL4.0ch05:576140-576154 (excerpt 141-155), kind STR; STR regions: 576140-576151 str T
  0  TTTTTTTTTTTTATA
  1  TTTTTTTTTTTTTA
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 2, 2  -> same genotype
```

**E13.** Tomato, sample SRS3394663: GATK deletes the A on both haplotypes, pvc on one.

```
excerpt  SL4.0ch05:576000-576300, 301 bases; excerpt position = position - 575999; sample SRS3394663
A  GATK-tomato:65098  576151 TA T 1/1
A  GATK-tomato:65099  576152 A T,* 2/2
B  pvc-tomato:71841  576140 TTTTTTTTTTTTA TTTTTTTTTTT 0/0
B  pvc-tomato:71842  576152 ATA TA,TTA,* 0/1
locus    SL4.0ch05:576140-576154 (excerpt 141-155), kind STR; STR regions: 576140-576151 str T
  0  TTTTTTTTTTTTATA
  1  TTTTTTTTTTTTTA
possible genotypes  A: 1/1   B: 0/1
row  A 1/1  B 0/1  possible 1, 1  records 2, 2  -> same alleles, other genotype
```

**E14.** Tomato, sample SRS3394701: GATK changes the A to T on both haplotypes, pvc on one.

```
excerpt  SL4.0ch05:576000-576300, 301 bases; excerpt position = position - 575999; sample SRS3394701
A  GATK-tomato:65098  576151 TA T 0/0
A  GATK-tomato:65099  576152 A T,* 1/1
B  pvc-tomato:71841  576140 TTTTTTTTTTTTA TTTTTTTTTTT 0/0
B  pvc-tomato:71842  576152 ATA TA,TTA,* 0/2
locus    SL4.0ch05:576140-576154 (excerpt 141-155), kind STR; STR regions: 576140-576151 str T
  0  TTTTTTTTTTTTATA
  1  TTTTTTTTTTTTTTA
possible genotypes  A: 1/1   B: 0/1
row  A 1/1  B 0/1  possible 1, 1  records 2, 2  -> same alleles, other genotype
```

**E15.** Tomato, sample SRS3394663. GATK writes a 5-base deletion and two records whose `*` it
covers; pvc the same deletion on both haplotypes and a contraction of the repeat on one.

```
excerpt  SL4.0ch05:58429278-58429578, 301 bases; excerpt position = position - 58429277; sample SRS3394663
A  GATK-tomato:74913  58429428 ATAAAT A 1/1
A  GATK-tomato:74914  58429429 T A,* 2/2
A  GATK-tomato:74915  58429433 T A,* 2/2
B  pvc-tomato:83029  58429420 AAAAAAAAAT AAAAAAAA,AAAAAT,AAAAAAA 0/2
B  pvc-tomato:83030  58429429 TAAATAAA AAAAAAA,AAA,* 2/2
locus    SL4.0ch05:58429419-58429436 (excerpt 142-159), kind STR; STR regions: 58429420-58429428 str A
  0  CAAAAAAAAATAAATAAA
  1  CAAAAAAAAAAAA
  2  CAAAAAAAA
possible genotypes  A: 1/1   B: 1/2
row  A 1/1  B 1/2  possible 1, 1  records 3, 2  -> different
```

**E16.** Tomato, sample SRS3394698: GATK gives two phased SNPs on both haplotypes, pvc the 5-base
deletion on both.

```
excerpt  SL4.0ch05:58429278-58429578, 301 bases; excerpt position = position - 58429277; sample SRS3394698
A  GATK-tomato:74913  58429428 ATAAAT A 0/0
A  GATK-tomato:74914  58429429 T A,* 1|1:58429429
A  GATK-tomato:74915  58429433 T A,* 1|1:58429429
B  pvc-tomato:83029  58429420 AAAAAAAAAT AAAAAAAA,AAAAAT,AAAAAAA 0/0
B  pvc-tomato:83030  58429429 TAAATAAA AAAAAAA,AAA,* 2/2
locus    SL4.0ch05:58429419-58429436 (excerpt 142-159), kind STR; STR regions: 58429420-58429428 str A
  0  CAAAAAAAAATAAATAAA
  1  CAAAAAAAAAAAAAAAAA
  2  CAAAAAAAAAAAA
possible genotypes  A: 1/1   B: 2/2
row  A 1/1  B 2/2  possible 1, 1  records 3, 2  -> different
```

**H1.** TR places the 8-base deletion after the repeat's two interruptions; HipSTR writes the whole
repeat, from its first base.

```
excerpt  chr1:113339083-113339383, 301 bases; excerpt position = position - 113339082; sample HG002
A  TR:1458  113339221 ATGTGTGTG A 1|0
B  HipSTR30:3213  113339207 GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 0|1
locus    chr1:113339221-113339260 (excerpt 139-178), kind STR; STR regions: 113339222-113339259 str TG
  0  ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  1  ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 1, 1  -> same genotype
```

**H2.** HipSTR against pvc, whose record covers the uninterrupted part of the repeat only.

```
excerpt  chr1:113339083-113339383, 301 bases; excerpt position = position - 113339082; sample HG002
A  HipSTR30:3213  113339207 GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 0|1
B  pvc-m07-30:3467  113339222 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 0/1
locus    chr1:113339221-113339260 (excerpt 139-178), kind STR; STR regions: 113339222-113339259 str TG
  0  ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  1  ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 1, 1  -> same genotype
```

**H3.** HipSTR writes the whole repeat with the inserted copy; TR writes an insertion and a SNP.

```
excerpt  chr1:72941135-72941435, 301 bases; excerpt position = position - 72941134; sample HG002
A  TR:1010  72941294 T TATCTATCC 1|1
A  TR:1011  72941298 C T 1|1
B  HipSTR30:3016  72941260 TCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTA TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTA 1|1
locus    chr1:72941245-72941319 (excerpt 111-185), kind STR; STR regions: 72941245-72941319 str-bundle
  0  AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCATTTTT
  1  AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTATCATTTTT
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 2, 1  -> same genotype
```

**H4.** TR against GIAB v4.2.1, which writes the same haplotype as two insertions, after 72941297
and after 72941298.

```
excerpt  chr1:72941135-72941435, 301 bases; excerpt position = position - 72941134; sample HG002
A  TR:1010  72941294 T TATCTATCC 1|1
A  TR:1011  72941298 C T 1|1
B  GIAB-HG002:98215  72941297 C CTATCCAT 1/1
B  GIAB-HG002:98216  72941298 C CT 1/1
locus    chr1:72941245-72941319 (excerpt 111-185), kind STR; STR regions: 72941245-72941319 str-bundle
  0  AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCATTTTT
  1  AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTATCATTTTT
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 2, 2  -> same genotype
```

**H5.** HipSTR writes the whole repeat with the insertion and the SNP on one haplotype; TR writes
them as two records phased on one haplotype.

```
excerpt  chr1:80541386-80541686, 301 bases; excerpt position = position - 80541385; sample HG002
A  TR:1070  80541526 C CAGATAGATAGATAGATAGATAGATAGAT 1|0
A  TR:1071  80541557 G A 1|0
B  HipSTR30:3058  80541499 AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGGTAGATAGATAGATAGAT AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT 0|1
locus    chr1:80541526-80541557 (excerpt 141-172), kind complex; STR regions: none
  0  CAGATAGATAGATAGATAGATAGATAGATAGG
  1  CAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 2, 1  -> same genotype
```

**H6.** The same HipSTR record against GIAB v4.2.1's two unphased records: two possible genotypes,
one of them HipSTR's.

```
excerpt  chr1:80541386-80541686, 301 bases; excerpt position = position - 80541385; sample HG002
A  GIAB-HG002:109203  80541526 C CAGATAGATAGATAGATAGATAGATAGAT 0/1
A  GIAB-HG002:109204  80541557 G A 0/1
B  HipSTR30:3058  80541499 AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGGTAGATAGATAGATAGAT AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT 0|1
locus    chr1:80541526-80541557 (excerpt 141-172), kind complex; STR regions: none
  0  CAGATAGATAGATAGATAGATAGATAGATAGG
  1  CAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA
  2  CAGATAGATAGATAGATAGATAGATAGATAGA
  3  CAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGG
possible genotypes  A: 0/1, 2/3   B: 0/1
row  A 0/1  B 0/1  possible 2, 1  records 2, 1  -> same genotype
```

**H7.** HipSTR's REF starts 4 bases before the repeat, and its ALT has the same length as REF; TR
writes four phased SNPs, one of them in that flank.

```
excerpt  chr1:14451792-14452092, 301 bases; excerpt position = position - 14451791; sample HG002
A  TR:311  14451931 C T 1|0
A  TR:312  14451937 T G 1|0
A  TR:313  14451938 G T 1|0
A  TR:314  14451940 T G 1|0
B  HipSTR30:2683  14451930 TCTTTTTTGTTTTGTTTTGTTTTG TTTTTTTGTTGTTGTTTTGTTTTG 0|1
locus    chr1:14451931-14451957 (excerpt 140-166), kind STR; STR regions: 14451932-14451957 str-bundle
  0  CTTTTTTGTTTTGTTTTGTTTTGTTTT
  1  TTTTTTGTTGTTGTTTTGTTTTGTTTT
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 4, 1  -> same genotype
```

**H8.** The same HipSTR record against GIAB v4.2.1, which writes the same haplotype as a deletion
and an insertion, unphased.

```
excerpt  chr1:14451792-14452092, 301 bases; excerpt position = position - 14451791; sample HG002
A  GIAB-HG002:19368  14451930 TC T 0/1
A  GIAB-HG002:19369  14451940 T TG 0/1
B  HipSTR30:2683  14451930 TCTTTTTTGTTTTGTTTTGTTTTG TTTTTTTGTTGTTGTTTTGTTTTG 0|1
locus    chr1:14451930-14451957 (excerpt 139-166), kind STR; STR regions: 14451932-14451957 str-bundle
  0  TCTTTTTTGTTTTGTTTTGTTTTGTTTT
  1  TTTTTTTGTTGTTGTTTTGTTTTGTTTT
  2  TCTTTTTTGTTGTTGTTTTGTTTTGTTTT
  3  TTTTTTTGTTTTGTTTTGTTTTGTTTT
possible genotypes  A: 0/1, 2/3   B: 0/1
row  A 0/1  B 0/1  possible 2, 1  records 2, 1  -> same genotype
```

**H9.** HipSTR's REF reaches 2 bases into the flank and its ALT carries TR's flank SNP at 61522963.

```
excerpt  chr1:61522820-61523120, 301 bases; excerpt position = position - 61522819; sample HG002
A  TR:910  61522963 G A 1|0
A  TR:911  61522975 A ATATATATATATATATATATATATATATATG 1|0
B  HipSTR30:2961  61522962 TGTATATATATATATA TATATATATATATATATATATATATATATATATATATATATATGTA 0|1
locus    chr1:61522963-61522979 (excerpt 144-160), kind complex; STR regions: none
  0  GTATATATATATATATG
  1  ATATATATATATATATATATATATATATATATATATATATATGTATG
possible genotypes  A: 0/1   B: 0/1
row  A 0/1  B 0/1  possible 1, 1  records 2, 1  -> same genotype
```

**H10.** HipSTR writes one phased record of two ALTs whose REF reaches past the repeat to carry the
flank SNP at 102475817; TR writes an insertion and a deletion at one position, on different
haplotypes, and the SNP.

```
excerpt  chr1:102475645-102475945, 301 bases; excerpt position = position - 102475644; sample HG002
A  TR:1332  102475773 T TAC 0|1
A  TR:1333  102475773 TACAC T 1|0
A  TR:1334  102475817 T C 0|1
B  HipSTR30:3163  102475774 ACACACACACACACACACACACACACACACACACACACACACAT ACACACACACACACACACACACACACACACACACACACAT,ACACACACACACACACACACACACACACACACACACACACACACAC 1|2
locus    chr1:102475773-102475817 (excerpt 129-173), kind STR; STR regions: 102475774-102475815 str AC
  0  TACACACACACACACACACACACACACACACACACACACACACAT
  1  TACACACACACACACACACACACACACACACACACACACAT
  2  TACACACACACACACACACACACACACACACACACACACACACACAC
possible genotypes  A: 1/2   B: 1/2
row  A 1/2  B 1/2  possible 1, 1  records 3, 1  -> same genotype
```

**H11.** The same HipSTR record against pvc's record of two ALTs and its `AT → AC` beside it,
unphased.

```
excerpt  chr1:102475645-102475945, 301 bases; excerpt position = position - 102475644; sample HG002
A  HipSTR30:3163  102475774 ACACACACACACACACACACACACACACACACACACACACACAT ACACACACACACACACACACACACACACACACACACACAT,ACACACACACACACACACACACACACACACACACACACACACACAC 1|2
B  pvc-m07-30:3413  102475774 ACACACACACACACACACACACACACACACACACACACACAC ACACACACACACACACACACACACACACACACACACAC,ACACACACACACACACACACACACACACACACACACACACACAC 1/2
B  pvc-m07-30:3414  102475816 AT AC 0/1
locus    chr1:102475773-102475817 (excerpt 129-173), kind STR; STR regions: 102475774-102475815 str AC
  0  TACACACACACACACACACACACACACACACACACACACACACAT
  1  TACACACACACACACACACACACACACACACACACACACAT
  2  TACACACACACACACACACACACACACACACACACACACACACACAC
  3  TACACACACACACACACACACACACACACACACACACACAC
  4  TACACACACACACACACACACACACACACACACACACACACACACAT
possible genotypes  A: 1/2   B: 1/2, 3/4
row  A 1/2  B 1/2  possible 1, 2  records 1, 2  -> same genotype
```

**H12.** HipSTR's second ALT carries G → C at 172778363, a flank base TR and pvc leave as it is; the
lengths agree.

```
excerpt  chr1:172778235-172778535, 301 bases; excerpt position = position - 172778234; sample HG002
A  TR:1784  172778363 G GCT 0|1
A  TR:1785  172778363 G GCTCTCTCTCTCTCTCTCT 1|0
B  HipSTR30:3357  172778362 AGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT AGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT,ACCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT 1|2
locus    chr1:172778363-172778407 (excerpt 129-173), kind STR; STR regions: 172778364-172778405 str CT
  0  GCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA
  1  GCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA
  2  GCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA
  3  CCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA
possible genotypes  A: 1/2   B: 2/3
row  A 1/2  B 2/3  possible 1, 1  records 2, 1  -> different
```

**H13.** The same 2-base deletion, on both haplotypes in HipSTR and on one in TR. HipSTR gives its
genotype a posterior probability, its FORMAT field Q, of 0.34.

```
excerpt  chr10:12391846-12392146, 301 bases; excerpt position = position - 12391845; sample HG002
A  TR:20498  12391976 AAC A 1|0
B  HipSTR30:3816  12391977 ACACACACACACACACACACACACACACACACACACACAC ACACACACACACACACACACACACACACACACACACAC 1|1
locus    chr10:12391976-12392017 (excerpt 131-172), kind STR; STR regions: 12391977-12392016 str AC
  0  AACACACACACACACACACACACACACACACACACACACACA
  1  AACACACACACACACACACACACACACACACACACACACA
possible genotypes  A: 0/1   B: 1/1
row  A 0/1  B 1/1  possible 1, 1  records 1, 1  -> same alleles, other genotype
```

**H14.** TR deletes the flank T at 93319249 on one haplotype; HipSTR deletes 4 bases of the repeat
on both.

```
excerpt  chr1:93319126-93319426, 301 bases; excerpt position = position - 93319125; sample HG002
A  TR:1161  93319248 AT A 1|0
B  HipSTR30:3117  93319250 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 1|1
locus    chr1:93319248-93319304 (excerpt 123-179), kind STR; STR regions: 93319250-93319303 str TG
  0  ATTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  1  ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  2  ATTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
possible genotypes  A: 0/1   B: 2/2
row  A 0/1  B 2/2  possible 1, 1  records 1, 1  -> different
```

**H15.** GIAB v4.2.1 against the same HipSTR record: the 4-base deletion on one haplotype against
both. TR and GIAB v4.2.1 disagree here.

```
excerpt  chr1:93319126-93319426, 301 bases; excerpt position = position - 93319125; sample HG002
A  GIAB-HG002:128803  93319249 TTGTG T 0/1
B  HipSTR30:3117  93319250 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 1|1
locus    chr1:93319249-93319304 (excerpt 124-179), kind STR; STR regions: 93319250-93319303 str TG
  0  TTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
  1  TTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
possible genotypes  A: 0/1   B: 1/1
row  A 0/1  B 1/1  possible 1, 1  records 1, 1  -> same alleles, other genotype
```

**O1.** TR writes four records, each heterozygous and phased, whose two haplotypes have the same
sequence; pvc writes one record, `1/1`.

```
excerpt  chr6:41399460-41399760, 301 bases; excerpt position = position - 41399459; sample HG002
A  TR:12602  41399584 GCCT G 0|1
A  TR:12603  41399613 CTCCTCT C 1|0
A  TR:12604  41399625 TTCC T 0|1
A  TR:12605  41399628 C T 1|0
B  pvc-touch30:9746  41399585 CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCT CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCT 1/1
locus    chr6:41399584-41399642 (excerpt 125-183), kind STR; STR regions: 41399585-41399641 str CCT
  0  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCTC
  1  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTC
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 4, 1  -> same genotype
```

**O2.** The same TR records against freebayes, which writes the haplotypes otherwise.

```
excerpt  chr6:41399460-41399760, 301 bases; excerpt position = position - 41399459; sample HG002
A  TR:12602  41399584 GCCT G 0|1
A  TR:12603  41399613 CTCCTCT C 1|0
A  TR:12604  41399625 TTCC T 0|1
A  TR:12605  41399628 C T 1|0
B  freebayes30:31443  41399613 CTCCTCTT CT,CTCCTCCT 1/2
B  freebayes30:31444  41399625 T C 0/1
B  freebayes30:31445  41399628 C T 0/1
locus    chr6:41399584-41399642 (excerpt 125-183), kind STR; STR regions: 41399585-41399641 str CCT
  0  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCTC
  1  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTC
  2  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCCTCCTCCTCCTCCTC
  3  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCTTCCTCCTCCTCCTC
  4  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCCTCCTCCTCCTCCTC
  5  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTC
  6  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCCTCCTCCTCCTCCTC
  7  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCTTCCTCCTCCTCCTC
  8  GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCCTCCTCCTCCTCCTC
possible genotypes  A: 1/1   B: 1/6, 2/3, 4/5, 7/8
row  A 1/1  B 1/6  possible 1, 4  records 4, 3  -> different
```

**O3.** One insertion of AAAAG written at three places in a run of A's followed by an AAAG repeat.
TR inserts it inside the A run; pvc writes the A run whole. pvc_eval's STR finder makes the A run
alone an STR region.

```
excerpt  chr3:37126720-37127020, 301 bases; excerpt position = position - 37126719; sample HG002
A  TR:5893  37126868 A AAAAAG 1|1
B  pvc-touch30:5857  37126860 AAAAAAAAAAAA AAAAAAAAAAAAAGAAA 1/1
locus    chr3:37126860-37126872 (excerpt 141-153), kind STR; STR regions: 37126860-37126871 str A
  0  AAAAAAAAAAAAG
  1  AAAAAAAAAAAAAGAAAG
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
```

**O4.** TR against HipSTR, which writes the AAAG repeat whole. Its record trims to an insertion,
whose span lies in the locus of O3.

```
excerpt  chr3:37126720-37127020, 301 bases; excerpt position = position - 37126719; sample HG002
A  TR:5893  37126868 A AAAAAG 1|1
B  HipSTR30:10361  37126869 AAAGAAAGAAAGAAAGAAAGAAAG AAAAGAAAGAAAGAAAGAAAGAAAGAAAG 1|1
locus    chr3:37126860-37126872 (excerpt 141-153), kind STR; STR regions: 37126860-37126871 str A
  0  AAAAAAAAAAAAG
  1  AAAAAAAAAAAAAGAAAG
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
```

**O5.** TR against freebayes, which writes both repeats and their flank bases.

```
excerpt  chr3:37126720-37127020, 301 bases; excerpt position = position - 37126719; sample HG002
A  TR:5893  37126868 A AAAAAG 1|1
B  freebayes30:16445  37126859 CAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAAGT CAAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAAGAAAGT 1/1
locus    chr3:37126860-37126872 (excerpt 141-153), kind STR; STR regions: 37126860-37126871 str A
  0  AAAAAAAAAAAAG
  1  AAAAAAAAAAAAAGAAAG
possible genotypes  A: 1/1   B: 1/1
row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
```

## 5. How it is verified

`tests/test_hard_cases.py` holds one test, `test_hard_case`, parametrized by the 50 cases, each
named by its case, S1 to O5, in pytest's report. The cases themselves are literals in
`tests/hard_cases.py`, one value per case: the excerpt's contig, offset and bases; the records of
each file as VCF lines in excerpt positions, each with its original line beside it as a comment; the
sample's name; and the expected row. Each test:

1. writes the excerpt as a FASTA with its `.fai`, and the two VCF files, to `tmp_path`;
2. checks the haplotypes with the oracle: the possible genotypes `tests/haplotype_oracle.py` gives
   each file over the locus are the case's, as sequences, or none for a file whose sample is
   invalid;
3. runs `compare()` with the default `CompareConfig` and `finder_threads=0`, and checks that it
   yields the case's rows, one for these cases, field by field: the locus (0-based and half-open in
   `OutcomeRow`), the kind, the STR regions with their types and motifs, the outcome, the locus
   alleles, both reported genotypes, both counts of possible genotypes, both invalid reasons and
   both counts of records. QUAL and GQ are not checked, since the cases set them.

Step 2 needs the oracle to give its possible genotypes, which it now uses only inside
`check_same_haplotypes`. It gains one public function, in the tests and not in `pvc_eval`:

```python
def find_possible_haplotypes(vcf: Path, reference: IndexedFasta, contig: str, start: int,
                             end: int, sample: int, *, ploidy: int = DEFAULT_PLOIDY
                             ) -> frozenset[tuple[str, ...]] | None: ...
```

It reads the records of one file with `parse_vcf`, takes those of its windows
(`doc/specs/respell.md` §8) inside the stretch from `start` to `end`, 0-based and half-open, and
gives the sample's possible genotypes over that stretch, with the phase assignments of the blocks of
all those windows taken together, as the oracle's check of a pair of windows takes them. Taken
window by window and joined, they would give genotypes no assignment gives: in S4, a second one for
TR, whose two windows hold one phase set. Each genotype is a sorted tuple of sequences, one per
haplotype, with the reference's bases between the windows. It gives `None` when the sample is
invalid, and raises `ValueError` when a window crosses the stretch's ends or the sample is over
`ORACLE_CAP`. No case here has a missing haplotype, so a missing one raises `ValueError` too.

**Each test can fail.** The plan shows the test failing before it is trusted, at each of its two
checks: the oracle's, with an ALT of file B changed by one base in a case that is *same genotype*;
and `compare()`'s, with the oracle's literal left as it is and one field of the expected row
changed, the outcome of one case and one locus allele of another.

**The excerpts against the reference.** Two more tests read each excerpt's stretch from a reference
and check that it is the case's bases, and that the STR regions of each locus are those the STR
finder gives over the excerpt widened by 5,000 bases on each side.
`test_grch38_excerpts_are_the_references` does it for the 45 GRCh38 cases, from the FASTA that
`PVC_EVAL_GRCH38` names, and `test_tomato_excerpts_are_the_references` for the five tomato cases,
E12 to E16, from the FASTA that `PVC_EVAL_TOMATO` names, the SL4.0 chromosomes at
`/Users/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa` on this machine. Each
is skipped on its own, with a message naming its variable, when the variable is not set or names no
file, as `tests/test_compare_real.py` is. The 50 tests of the cases need neither.

## 6. A `*` is measured against its record's REF as written

The owner decided on 21 October 2026, answering Open 1 of this spec's draft, that a deletion covers
a `*` when one of its placements deletes a position of the `*` record's REF as the file wrote it,
before trimming (§5.2), and no longer of the trimmed REF. A `*` names the bases of the REF its
caller wrote, whichever side the caller padded, and pvc pads such records on purpose
(`doc/devel/ng/spec/str_edge_deletion_alleles.md` §6). The rest of the rule stands: the placements
of a deletion, the other edits that cover a `*` as written, and a haplotype `.` at a record, which
might carry any of its ALTs.

The options not taken were to keep measuring against the trimmed REF, which leaves S11 and E5
*invalid in B*, and to trim a `*` record from the left first, as an edit is trimmed, which gives
the same rows here but works only for a record padded on the left.

What the decision changes in the cases: S11 and E5 become *same genotype*, with the same locus
alleles as before, S11 reporting `2/3` and `2/3` and E5 `1/2` and `1/2`. No other case changes,
and neither do `doc/specs/reconciliation.md` R1 and R2, whose `*` records are not padded: worked
with a copy of the oracle changed so, by this spec's writer on S11, E5, E2, E3, E6 and S1, and by
its reviewer on every case. A file can now cover a `*` with a deletion of a flank base it padded
the REF over; no case here does.

The review of the code that follows the rule found that a deletion in the padding could still fall
in another locus than the `*` record, which is chained by its span, built from its trimmed REF;
the `*` was then judged uncovered there. So the owner decided, the same day, that the span of a
record with a `*` holds its REF as written (`doc/spec.md` §5.2). That moves the locus of three
cases, whose padded `*` record reaches past it, and changes no outcome or genotype: E3 starts at
3029115 in place of 3029116, and E15 and E16 end at 58429436 in place of 58429434; their locus
alleles hold the bases added.

## 7. Changes to other specs and to the code

In this order, each in a commit of its own, the specs before the code:

1. `doc/spec.md` §5.4: the `*` rule above, in place of "a position of the `*` record's normalised
   REF". `doc/specs/reconciliation.md` §4.2 likewise, with a case of a padded `*` record made from
   S11's records over their own few bases.
2. `doc/specs/respell.md`: where it judges which records a re-spelling keeps whole by "the REF of a
   record with a `*` ALT", that REF is the REF as written, so that no re-spelling changes whether a
   padded `*` is covered. The `*` re-spelling keeps adding a `*` only where a deletion reaches the
   trimmed REF, which the REF as written holds.
3. The code: `VcfRecord` keeps the stretch its REF covered as written, which it does not now, and
   `pvc_eval.reconcile`, `tests/haplotype_oracle.py` and the re-spelling program measure a `*`
   against it. The tests of R1 and R2 pass unchanged, and the new case of `reconciliation.md` is
   added.
4. Once the owner approves the cases, `doc/spec.md` §7.5 becomes a pointer to this spec, and §10
   marks Open 3 answered, with the date.

The tests of S11 and E5 fail until step 3 is done; the plan puts it before the tests of the cases.

## 8. Not in this spec

- **Cases reconstructed from a report's words.** Every case has its records on disk. One report
  describes a record that is no longer there: `giab_unexplained_fn_2026-10-03.md` §4 gives S11's pvc
  record as `TATACAC → TATAC` `1/1` from a run whose worktree is gone. With that record the row
  would be *different*, since the 2-base deletion would then be on both haplotypes. The later run of
  pvc-300x-HG004 is used instead.
- **Two HipSTR records in one locus** (§3): a test of it would be made up, and belongs with the
  owner's decision on how to read `|` without PS, if they take one.
- **vcfeval's verdicts** on the same records, which pvc's reports quote: the comparison of
  evaluators is `doc/spec.md` §8.
- **pvc's own errors** that the cases show, S4, S14, E8 and E10: they are pvc's to fix, and pvc_eval
  reports them as *different*, as it should.

