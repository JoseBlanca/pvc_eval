# Work report: the re-spelling program

*Done, 21 October 2026. The plan is `doc/plans/respell.md`, built from `doc/specs/respell.md`,
on the branch `respell`; the last commit of code is 630d173. On the owner's order of 21 October
2026, `main` was merged into the branch and the GIAB test given the `num_threads` that `main`'s
`write_str_regions` asks for (496ebe3): the five checks gave "1886 passed, 9 skipped" and the four
GIAB tests passed in 132 s; then the branch was merged into `main`. Nothing is pushed.*

## For the owner

The re-spelling program of `doc/spec.md` §7.4 is built, and pvc_eval passes the check it exists
for. A VCF file compared with each of its four re-spellings gives, at every locus and for every
sample, the same result as the file against itself:

- on GIAB's HG002 v4.2.1 benchmark VCF, chr20, against each re-spelling: all 81,257 rows *same
  genotype*, for all four (8 min 30 s together on the owner's Mac, 21 October 2026);
- on 5,000 random small files with repeats for each re-spelling, run outside the suite at 630d173:
  no failure of the oracle described below, and none of the §7.4 check.

What exists now, all under `tests/`, nothing in `src/`:

- `tests/respell.py`: four functions, each writing a VCF again in other words. Each record is
  decomposed into one-base and one-indel records; indels are shifted right within their repeats;
  `*` is added under deletions where the phase is known; and neighbouring records are merged
  into complex ones. Each leaves unchanged a record whose edits overlap, or lie next to, the edits
  of another record that may be on the same haplotype, as you decided on 19 October 2026.
- `tests/haplotype_oracle.py`: an oracle, a check written apart from pvc_eval's reconciliation,
  loci and evaluation, that rebuilds each sample's haplotypes from two files and compares them.
  Where the phase is not known it compares the possible genotypes over each stretch where the two
  files' records chain together, and over each two such stretches one phase set joins.
- `tests/respell_check.py`: the check of §7.4 through `parse_vcf`, `build_loci`,
  `reconcile_loci` and `make_outcome_rows`.
- A generator of small files with repeats, and the properties over it: the oracle, bcftools and
  §7.4 for each re-spelling. The suite went from 914 to 1,385 passing tests and takes 96 s.

**What is asked of you: the merge into `main`.** `main` has moved since this branch began. It now
holds `compare()` and the STR finder's speed-up. A trial merge has no textual conflict and passes
pytest (1,777 passed, 9 skipped), but `ty check` finds one call that no longer fits:
`write_str_regions` on `main` takes a required `num_threads`, and the GIAB test of this branch
calls it without one. The test is skipped unless the GIAB variables are set, so pytest does not
see it.

- **I merge `main` into `respell` first, add the argument, run the five checks (`uv sync
  --locked`, `ruff format --check`, `ruff check`, `ty check`, `pytest`) and the GIAB tests, and
  hand the branch back to you for its merge into `main` (recommended).** A merge commit and one
  small commit on this branch, about 15 minutes with the GIAB tests; `main` then receives a
  branch that passes as it stands.
- **You merge `respell` into `main` as it is**, and the call is fixed on `main` afterwards. `main`
  then has a failing `ty check` until the fix.

With `compare()`, the function that runs a whole comparison, now on `main`, the plan's one
deferred change can be made: the check of §7.4 taking its rows from `compare()`. It changes no
result. I can write it as a plan of one milestone when you want it; nothing waits on it.

Five things you may want to know, each in the milestone that found it:

- The `*` re-spelling writes no `*` on HG002 chr20, which has no phased genotype and no `0/0`
  record; it is checked on generated files only. A phased or a multi-sample VCF would check it
  on real data (milestone 5).
- In pvc_eval itself, outside this plan, `parse_vcf` computes a record's span by sliding on the
  VCF's own letters, while reconciliation builds each haplotype from the reference's bases.
  Where a REF matches the reference only by an IUPAC code, `A` over `M`, the two can disagree
  (milestone 4).
- Also outside this plan: `parse_vcf` accepts `.` as one ALT among others (`G A,.`) and reads it
  as an allele (milestone 5); and a gzip file with one byte changed is reported as having no
  `#CHROM` line, not as damaged (milestone 1).
- The spec changed during the work in a dozen places. Each change is a commit of its own, made
  before the code, and none changes what pvc_eval gives. The largest: a new phase set is the n-th
  whole number the input does not use; the shift and the merge work only where the reference is
  A, C, G or T; and a merge that would make a heterozygous sample homozygous is not made. Such a
  merge would remove a block, a set of genotypes whose phase pvc_eval has to try both ways, so a
  sample over the phasing cap (4,096 phase assignments) in file A would fall under it in file B
  and the two would get different outcomes (milestones 4 to 6).
- How sure this is: every case of the spec is a test with its output as a literal; in every
  milestone, reviewers broke the code on purpose, and each break they made now fails a test;
  four defects that changed a haplotype or a block were found by review and fixed, each with the
  case that showed it; over the six milestones the reviewers ran more than 100,000 generated
  files through pvc_eval's reconciliation besides the suite's.

The checks on the commit the work starts from, 9dbec00: `uv sync --locked` clean, `ruff format
--check`, `ruff check` and `ty check` pass, `uv run pytest -q` gives "914 passed, 4 skipped"; the
skipped tests need pvc and GRCh38.

## Milestone 1: decomposition, and what the four re-spellings share

Done in b957dea and 52afc56, and fixed in a3b208f to 26f8b51. `respell_by_decomposing` writes a
VCF with each record nothing touches decomposed into one-base and one-indel records, and
`tests/respell.py` holds what the other three re-spellings will use: the reading, the records
kept, the touched test with the phase, the reach, the writing, new phase sets, the header, the
order and a BGZF writer. The reconciliation property makes its file B with it, and
`_make_decomposed_lines` is gone.

- `uv run pytest tests/test_respell.py` gives "82 passed": D1 to D12, every row of the spec's
  §2.2, the header, bgzip, the output that is the input, the order, and the reach against
  `parse_vcf`'s spans. The plan asked for at least 25.
- `uv run pytest tests/test_reconcile_property.py` passes, with file B made by
  `respell_by_decomposing`; `grep -rn _make_decomposed_lines tests` finds nothing.
- The 20,000-SNP timing test, with 2,600 touching records, takes 0.23 s against a limit of 10 s;
  a sweep that keeps every earlier record takes 133 s and fails it.
- The five checks pass on 26f8b51; the whole suite gives "984 passed, 4 skipped".

Changed in the spec: §2.3 asked for a `Type=String` PS line when an input PS is not a whole
number, and also that the output give bcftools no message the input does not. bcftools 1.24
warns about a PS declared `String`, so the two could not both hold; two reviewers found it. The
spec now adds no PS line then, which leaves the output with the input's own warnings (6294b8e).
It changes only the header of the test program's output, nothing pvc_eval gives.

Changed in the plan: the case notation and the D cases are in `tests/respell_cases.py`, so that
milestones 2 and 3 can run them; `tests/conftest.py` has one VCF writer for the cases.

The review found one blocker: an output named as the input in another case, on macOS, or a hard
link to it, overwrote the input instead of raising. Fixed, with two tests. The correctness and
spec reviewers each compared input and output through `reconcile` over generated files, phased
and not, ploidy 1 and 2, records at both ends of a contig: 16,000 files in all, with no
difference. Fixed besides: seven plausible wrong rules no test could fail on (the touch on a
second sample, a phased GT with no PS against an unphased one, the records next to a `*` on
either side, a slide past one read of the reference, the order of a kept and a rewritten record
at one POS, contigs not in name order), each now a case; the test of the records near a `*` now
takes what any re-spelling writes; each line read is checked against the record matched to it;
names from the verb table. Not taken: a PS above 2^31, which bcftools sets to missing, and one
of more than 4,300 digits, which Python will not convert; no real file has either.

For the owner: a gzip file with one byte changed is reported by `parse_vcf` as having no
`#CHROM` line after its `##` lines, not as damaged. It is outside this plan.

## Milestone 2: the oracle and the generated files

Done in 5691144 and d954d50, fixed in a first round, 3d62c69 to 705187d, and a second, 6bbb7bb to
bc18ac6. `check_same_haplotypes` rebuilds each sample's haplotypes from two files and compares
them window by window and over pairs of windows a phase set joins, without pvc_eval's
reconciliation; `tests/respell_generator.py` makes the small files with repeats that the
properties of every re-spelling will use.

- `uv run pytest tests/test_haplotype_oracle.py` passes: O1 to O6 with their messages as
  literals, a case for each rule of the spec's §8, the cap, the samples, and that the oracle
  loads none of pvc_eval's reconciliation, loci or evaluation.
- In `tests/test_respell.py`, every D case passes the oracle and bcftools, and the oracle
  property over 300 random decomposed files passes, with at least one in ten rewritten and none
  over the cap; the bcftools property over 50 files passes.
- Over one run of 300 files, 246 had a rewritten record. Over 5,000 random files outside the
  suite, the decomposition differed from its input in none: 3,633 rewritten, 1,375 with a pair
  of windows a phase set joins, none over the cap.
- The five checks pass on bc18ac6; the whole suite gives "1055 passed, 4 skipped".

Changed in the spec (3d62c69, 3ec59d3, bfbe79d): a `.` covers a `*` only through an ALT of its
own record, as in reconciliation, where the spec said any `.` might; a sample's malformed GT or
ploidy mismatch is compared whatever the cap, and the tests assert no window over the cap; the
oracle holds both files whole, which on 20,000 records took 34 MB and 3.3 s; and the generator
gives at least one file in five a pair of windows a phase set joins, by phasing every GT of three
samples in four in one phase set. With phase drawn as the spec first said, one file in fourteen
had a pair, and a defect seen only between windows went unseen in about one run of 300 in nine.

The review found three blockers. Records with no ALT were given GTs whose `1` is malformed, so
half of them made their windows compare nothing. Four rules of the oracle, among them that a
collision comes before a missing haplotype and that `.` counts as an allele of a heterozygous GT,
could be broken with no test failing. And the oracle held both files whole, against the spec's
§10, which now says what the code does. Fixed besides: a `*` beside an unrelated `.` taken as
missing, which would have let a wrong `*` re-spelling pass; a window over the cap passed
silently, as a decomposition that lost the phase of 14 pieces did; failure messages of up to 145
KB, now the difference and its counts; records with no ALT placed by scanning every window, 24 s
on 8,000 records, now 0.6 s on 20,000; and one table of re-spellings that the case tests and the
properties run over, so the later milestones add a row. The second looks confirmed each blocker's
fix and found one more case, a record with no ALT just after a window, now tested. Not taken: a
record with no ALT and `0/.` across two windows links them, which changes nothing, since its `.`
haplotype is missing in both.

## Milestone 3: the check of §7.4, for decomposition

Done in d4cd386 and fixed in 9ab3d38 to e9b6e84. `check_respelling_outcomes` takes a file and its
re-spelling through `parse_vcf`, `build_loci`, `reconcile_loci` and `make_outcome_rows`, and
asserts at every locus and for every sample what `doc/spec.md` §7.4 asks: an allowed outcome,
*too complex* only where file A's own result is, and file B's result equal to A's.

- `uv run pytest tests/test_respell_check.py` passes: every D case; the AT repeat through the STR
  finder, one STR locus of *same genotype*; six files B that are not re-spellings, each raising
  with its message (a wrong sequence, a lost phase, a wrong second sample, a 14-base MNP *too
  complex* in B alone, and two of the messages); a file A that is itself *too complex*, which
  passes; the property over 300 decomposed files with their repeats as a `--str-regions` BED,
  with at least one row in an STR region; and the GIAB test with its path checks.
- **GIAB HG002 v4.2.1, chr20, against its decomposition: it passes.** 85,951 records decomposed
  into 86,896; the oracle compared 84,717 windows and no pair, since the VCF is unphased; all
  81,257 rows are *same genotype*. 125 s and 455 MB, of which 116 s are the STR finder, which on
  this branch is the version before its speed-up reached `main`. The test reads the VCF and
  GRCh38 from `PVC_EVAL_HG002_VCF` and `PVC_EVAL_GRCH38`, as the compare session's spec settles.
- Over 5,000 generated files outside the suite, no failure; of 12,152 rows, *same genotype*
  4,171, *invalid in both* 3,351, *missing in both* 2,580, *reference in both* 2,050, *too
  complex* 0. *Invalid in both* is 28% of the rows because a contig of 12 to 60 bases is often one
  locus, and one malformed GT makes a sample invalid over all of it.
- The five checks pass on e9b6e84; the whole suite gives "1084 passed, 5 skipped".

Changed in the plan (5d44fa5): milestone 4 runs the AT repeat case with its file B made by the
shift, which this milestone wrote by hand; and when `compare()` reaches `main`, the check takes
the reference as a path, as `compare()` does.

The review found no row that the check passes while pvc_eval treats B differently from A. It
found one blocker, the outcome missing from one message, which `doc/specs/respell.md` §7 asks
for; fixed, confirmed by the reviewer's second look. Fixed besides: a check of the first sample
alone, or one that rejected every *too complex* row, passed every test, and each now has a case;
messages of up to 104,571 characters, now the counts and the difference; a GIAB variable set to a
wrong path skipped the test, and now fails it saying what is wrong; a file with no samples passed
against any other; the cause in A is checked for *too complex* only, since `evaluate_pair` gives
*missing in both* and *invalid in both* only when A has the cause; and the harness of the case
tests and the properties is in one place, so milestones 4 to 6 add a row of the table of
re-spellings. Not taken now: the reference as a path, deferred to `compare()` as above.

## Milestone 4: shift

Done in 5e28f5d and fixed in f3afcea to 92df6bb. `respell_by_shifting` moves each indel that
nothing touches right within its repeat, as far as the spec's §4 allows, or a distance drawn
from a seed, with the VCF genotypes as read.

- `uv run pytest tests/test_respell.py tests/test_respell_check.py` passes: S1 to S9 with their
  outputs as literals; the records the shift keeps; the phase cases; the blocking edit just past
  the last placement; a long REF behind a short one; the span limit at 10,000 and 10,001 bases;
  IUPAC bases in the reference; the seed; the AT repeat with its file B now made by the shift;
  the oracle, bcftools and §7.4 properties.
- **GIAB HG002 v4.2.1, chr20, against its shift: it passes.** 10,458 of 85,951 records shifted,
  all 81,257 rows *same genotype*, 126 s, 557 MB, of which the shift itself takes 1.1 s.
- Over 5,000 generated files outside the suite, 1,174 had a shifted record, with no failure of
  the oracle or of §7.4. The reviewers ran 16,000 more files and 2,700 of a generator of their
  own through `reconcile` and the §7.4 check, seeded and not: no difference.
- The five checks pass on 92df6bb; the whole suite gives "1175 passed, 6 skipped".

Changed in the spec: the shift slides only over the reference's own A, C, G or T (abaf79d); a PS
line is declared `Integer` only when every PS fits in 32 bits, and a GT line is added when the
header lacks one (33a1c41); the rewritten share is counted on a fixed set of 300 files, since in
random runs the shifted files ranged from 22 to 77 around the 30 asserted (2617cca).

The review found one wrong shift: where the VCF's REF matches the reference only by an IUPAC code,
`A` over `M`, the shift slid on the VCF's letters and changed the haplotype (reference `GCMAAAT`,
`2 CA C`); fixed, confirmed by the reviewer's second look, which also notes that one of its four
cases, the slide up to an `N`, would pass without the fix. Fixed besides: two outputs bcftools
warned about, a PS of 2^31 or more and an input with no GT in its FORMAT; the tests of a bad input
for every re-spelling, not only the decomposition; four plausible wrong rules no test could fail
on, at the end of the slide and in the search by REF; one search by REF where there were three;
the seed checked as a whole number; and the generator draws a repeat's indel at either edge again.
Not taken: splitting `tests/respell.py`, which reads by its sections at 900 lines.

For the owner, outside this plan: `parse_vcf` gives a record its span by sliding on the VCF's own
letters, while reconciliation builds each haplotype from the reference's bases. Where a REF
matches the reference only by an IUPAC code, as above, the two can disagree. It is in pvc_eval's
code, and only a reference with IUPAC codes under a deletion reaches it.

## Milestone 5: spanning deletions

Done in 249c68d and fixed in c41246d to 483cd2f. `respell_with_spanning_deletions` writes a `*`
on each haplotype that carries a deletion over a record whose REF it deletes, where the spec's
§5 says the phase between the two is known; GATK's unphased case is kept as read, as the owner
decided for Open 1.

- `uv run pytest tests/test_respell.py tests/test_respell_check.py` passes: T1 to T8 and the
  trimmed-REF case `12 TC TA`, with their outputs as literals, through the oracle, bcftools and
  §7.4; the cases of each of §5's three ways the phase is known and of what is kept; the oracle,
  bcftools and §7.4 properties.
- **GIAB HG002 v4.2.1, chr20: it passes, and writes no `*`.** chr20 has no phased genotype and
  no `0/0` record, so none of §5's three ways holds; all 81,257 rows are *same genotype*. A
  check of this re-spelling on real data needs a phased or a multi-sample VCF.
- Over 5,000 random generated files, 782 got a `*` (15.6%), and every one passed the oracle and
  §7.4. Two reviewers ran 28,000 and 7,500 more files, one with a generator of its own that
  reached every branch of the rule, through `reconcile`, the oracle and a reading of §5 written
  apart from the code: no difference.
- The five checks pass on 483cd2f; the whole suite gives "1268 passed, 7 skipped".

Changed in the spec: a new phase set is the n-th whole number, from 1, that no PS of the input
uses (d9c746f), where it was numbered from the input's largest PS, which went past bcftools'
32-bit limit and Python's 4,300 digits on inputs with a PS near either; and an all-REF record
under deletions of two phase units gets no `*`, also where one of them shares its own phase set
(bf23b11).

The review found one blocker, found by two reviewers from two sides: condition 1 of §5, both
genotypes phased in one phase set, had no test that could fail, and each of three broken
versions wrote a `*` that changed reconciliation's result while the whole suite passed. Fixed
with cases a broken rule fails, confirmed by the second look. Fixed besides: two unphased
deletions counted as one phase unit; the share of files with a `*` met its bound only on the
fixed files, and the generator now draws more deletions over other records; one function gives
the phase unit of a genotype, which the merge will use; and the search for deletions runs once
per record, not once per sample.

For the owner, outside this plan: `parse_vcf` accepts `.` as one ALT among others (`G A,.`) and
reads it as an allele.

## Milestone 6: merge

Done in f4450ac and fixed in 02c1c51 to 630d173. `respell_by_merging` joins runs of neighbouring
records into one complex record, where the spec's §6 allows it, with as many ALTs as the samples
need.

- `uv run pytest tests/test_respell.py tests/test_respell_check.py` passes: M1 to M12 and the run
  all REF with their outputs as literals; the cases of each condition of §6; the merged phase,
  QUAL and order; the join probability, each of four outputs within 400 to 600 of 2,000 seeds;
  a sample at the phasing cap through the §7.4 check; the oracle, bcftools and §7.4 properties.
- **GIAB HG002 v4.2.1, chr20, against its merge: it passes.** Without a seed 594 records are
  merged into 290, with seed 20261020 293 into 145; all 81,257 rows *same genotype*; 130 s, of
  which the merge takes 6.5 s.
- On 1,000 fixed generated files each re-spelling rewrites a record in at least 110 (decomposition
  736, shift 142, `*` 110, merge 156); the generated files now leave the seed out one time in
  three, so the properties check the longest runs and the largest shifts too.
- The five checks pass on 630d173; the whole suite gives "1385 passed, 8 skipped".

Changed in the spec (e18d686, 623f6c3, 51effc9): a run is merged only over plain bases; symbolic
records do not stop it; a run that would give a sample with a heterozygous genotype equal merged
alleles is not merged; and the rewritten share is counted on 1,000 fixed files, where 300 left the
shift and the `*` 1 and 2 files above their bound.

The review found one blocker: phased records whose edits spell one sequence on both haplotypes,
as one deletion of a repeat written at two placements, merged into a homozygous record and
removed a block. At the phasing cap the §7.4 check then failed where pvc_eval is right: file A
*too complex*, file B 4,096 genotypes. The oracle could not see it, since the haplotypes are the
same. Fixed, confirmed by the reviewer's second look over 12,000 more re-spellings. Fixed
besides: a merged REF over a reference letter that is no IUPAC code, which `parse_vcf` rejected;
four plausible wrong rules no test could fail on; a merge whose time grew with the cube of a
chain of records it could not join, 43 s for 800 records, now 1.2 s; a duplicate header line
when the input puts `ID=` after another key; and four repetitions in the code. Not taken: runs
formed leftward too, which only records written untrimmed show.

## The check of the whole plan

On 630d173, 21 October 2026: the four GIAB tests, with `PVC_EVAL_HG002_VCF` and
`PVC_EVAL_GRCH38` set, "4 passed" in 510 s. Outside the suite, 5,000 random generated files for
each re-spelling, each with the generator's seed:

| re-spelling | rewritten | failures | same genotype | invalid in both | missing in both | reference in both |
|---|---|---|---|---|---|---|
| decomposition | 3,575 | 0 | 4,159 | 3,261 | 2,637 | 2,509 |
| shift | 601 | 0 | 4,267 | 3,119 | 2,763 | 2,322 |
| `*` | 593 | 0 | 4,129 | 3,121 | 2,548 | 2,245 |
| merge | 702 | 0 | 3,830 | 3,166 | 2,472 | 2,536 |

No row was *too complex*: the generated files are too small to reach the cap, which only the hand
cases reach.

## For whoever next revises the skills

The owner can stop here.

- **A reviewer's worktree is needed until its second look.** In milestone 2 the session removed
  the spec and tests reviewers' worktrees as soon as they reported, and neither could be resumed
  for its second look; fresh reviewers had to be briefed instead. The coding skill could say to
  keep each reviewer's worktree until the milestone is ticked.
- **A deliberate break can hang.** In milestone 4 one of the tests reviewer's breaks made a slide
  loop forever, and its run held the review for 47 minutes until the session stopped it.
  `review.md` could ask every break to run with a time limit, as the session asked from
  milestone 5 on.
- **Generated shares around a bound are a lottery.** Asserting "one file in ten" on random runs
  of 300 failed about one run in ten; on 300 fixed files the margins fell to 1 or 2 when the
  generator changed. Counting on 1,000 fixed files held. A rule for the planning skill: a share
  asserted over generated files is counted on fixed files, with a margin measured over many runs.
- **The oracle needed its own review as much as the code.** It accepted a wrong `*` beside a `.`
  and passed silently over windows over its cap, and its first generator gave pairs of windows to
  one file in fourteen. An oracle is code under test.
- **Merging `main` back** was not in the plan; a trial merge, made in a worktree of its own at
  the end, found an argument added on `main` that the branch's slow test lacked. A plan that runs
  beside other branches could end with that step.
- **Tokens per milestone**, from the `Agent` tool's results. The implementer's count includes its
  fixes, and each reviewer's its second look.

| milestone | implementation and fixes | review (5 categories) |
|---|---|---|
| 1, decomposition | 250,000 | 489,000 |
| 2, oracle and generator | 349,000 | 582,000, with two fresh second looks |
| 3, check of §7.4 | 213,000 | 312,000 |
| 4, shift | 321,000 | 543,000 |
| 5, spanning deletions | 268,000 | 486,000 |
| 6, merge | 373,000, with the final check | 523,000 |
