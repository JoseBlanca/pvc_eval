# Work report: the re-spelling program

*Under way, 20 October 2026. The plan is `doc/plans/respell.md`, built from
`doc/specs/respell.md`, on the branch `respell`. Nothing is merged into `main` or pushed.*

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
