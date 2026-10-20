# Plan: reading only the contigs `--regions` names

*Status: under way, 20 October 2026. The owner asked for the spec, `doc/specs/regions_read.md`,
to be reviewed and built without waiting for their approval. On the branch `regions-read`, in
`.claude/worktrees/regions-read`; the report is `doc/reports/regions-read.md`.*

One milestone, since the change is one keyword argument of `parse_vcf` and one call in
`compare()`, with their tests.

## Before the milestone

On c90733e: the five checks of `python.md` pass, with "1306 passed, 5 skipped".

## Milestones

- [ ] 1. `parse_vcf` reads only the contigs it is given, and `compare()` gives it those
  `--regions` names

**What it gives.** A run with `--regions` reads the records of the contigs the regions name, and
skips the lines of the others after the four checks of the spec.

**The spec it builds.** `doc/specs/regions_read.md`, whole.

**Deliverables.**

- `contigs` in `parse_vcf` (`src/pvc_eval/vcf.py`), and its use in `src/pvc_eval/compare.py`;
  the REF error's wording for the first record read.
- Every test of the spec's "How it is verified" but the run on real data, in
  `tests/test_vcf.py`, `tests/test_compare.py` and `tests/test_compare_threads.py`: at least 14
  new tests, all passing; the five checks pass.
- The run of HG002 against itself over chr20, with the command and the variables of
  `doc/specs/compare.md`, "Decided by the owner" point 6: the same 81,257 rows and summary as
  before, its time and peak in the report; and `tests/test_compare_real.py` passing.

**What could go wrong.** The order of the checks on a skipped line, which decides which message a
line with two faults gives, and a name that begins another (chr1 and chr10). The spec's tests name
both.
