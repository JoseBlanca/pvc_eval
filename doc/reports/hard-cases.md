# Report: pvc's hard cases

*Done, 21 October 2026. The work of `doc/plans/hard-cases.md` on the branch `hard-cases`, which
asks the owner for its merge into `main`.*

**What exists now.** The 50 hard cases of `doc/specs/hard_cases.md`, which answer `doc/spec.md`
Open 3, are tests: `uv run pytest tests/test_hard_cases.py` runs each case through the haplotype
oracle and through `compare()` on 301 bases of GRCh38 or SL4.0, and all 50 give the row the owner
approved. Two more tests check the excerpts and their STR regions against the references when
`PVC_EVAL_GRCH38` and `PVC_EVAL_TOMATO` name them, and a third checks that they would catch an STR
region across a locus's end. pvc_eval measures a `*` against its record's REF as written, and the
span of a record with a `*` holds that REF, as the owner decided on 21 October 2026; `doc/spec.md`
§5.2, §5.4, §7.5 and §10, `doc/specs/reconciliation.md` and `doc/specs/respell.md` say so.

**Checks at the end**, at 2e7e3f3: `ruff format --check`, `ruff check` and `ty check` pass;
`uv run pytest -q` gives "1968 passed, 11 skipped" (1886 passed, 9 skipped at a050b57); the hard
cases give "51 passed, 2 skipped" without the two variables and "53 passed" with them.

**What is asked of the owner.** The merge of `hard-cases` into `main`, and nothing is pushed. Two
things found on the way are not fixed here: a REF of non-ASCII letters stops the run with a
traceback (milestone 1, below), and the span-limit warning's text, which no longer fits a record
with a `*`.

## Milestone 1: a `*` measured against its record's REF as written

Done at 5590a35, from 35d6eef. `uv run pytest -q` gives "1905 passed, 9 skipped": 19 new tests,
each shown failing on the code before its change. `ruff format --check`, `ruff check` and `ty
check` pass.

- `VcfRecord.written_ref`, set by `parse_vcf` where the REF as written differs from the trimmed
  REF: `tests/test_vcf.py`, padded and unpadded records, `None`, and the `ValueError`.
- R10 at `reconcile`, `(1, 2)`, which 35d6eef gave as invalid; a deletion in the left padding,
  at `reconcile` and at the oracle; the re-spelling program keeping a record next to a padded
  `*` record's REF as written, on the right and on the left.
- The span of a record with a `*` holds its REF as written: `parse_vcf` with and without the span
  limit, `compare()` on the reviewer's case, which gave *invalid in A* and *same genotype* at two
  loci where it now gives one locus, *same genotype*, and the oracle on a deletion next to the
  padding.
- The 50 hard cases, run through `compare()` and the oracle at 5590a35 by a script of the
  session's: both agree on all 50. S11 and E5 give *same genotype*, as the spec says.

**What the review found, and the owner decided.** Five reviewers, one per category.

- Two reviewers, correctness and spec, found apart that a deletion in the padding of a `*` record
  could fall in another locus, since spans were built from the trimmed REF; the owner decided on
  21 October 2026 that the span of a record with a `*` holds its REF as written (`doc/spec.md`
  §5.2). It moved R10's locus to 10–19, a case this branch added, and the loci of E3, E15 and E16
  by one or two bases, with no outcome changed; the specs were corrected first.
- The design reviewer measured the field on every record at 721 bytes per record against 601
  before, on 100,000 one-sample SNPs; the owner decided it is set only where it differs, 610.
- Three reviewers found no test of the padding on the left of the trimmed REF; the tests above
  were added.
- Not taken: the errors reviewer's question whether a `written_ref` with a negative start should
  be refused; the spec asks only that it hold the trimmed REF, as it asks of a span.
- The correctness reviewer's second look at the span fix found it settled, with its tests failing
  without it, and parse_vcf's order kept on 1,500 random records. Its nit, left for the owner: the
  span-limit warning still says a cut record "was given only its REF as its span", where a record
  with a `*` now keeps its REF as written too.
- The tests reviewer noted that no property test compares `reconcile` and the oracle on records
  with a `*`; the two read one field, and the R10 tests guard it. Left as it is.

**For the owner.** Outside this plan, the errors reviewer found that a REF of non-ASCII letters,
`ß` against a reference `G`, stops the run with a Python traceback from `_check_ref` in
`pvc_eval.vcf` instead of an error naming the line, and that `ſ` is read as the IUPAC code `S`.
It is on `main` too, and is not fixed here.

## Milestone 2: the 50 hard cases through `compare()`

Done at 2e7e3f3, from 4d8d44d.

- `tests/hard_cases.py`: one `HardCase` per case, generated from §4.5 of the spec by a script kept
  outside the repository, with each record's and each row's line of §4.5 beside it as a comment and
  the possible genotypes as §4.5's indices. Two reviewers parsed §4.5 with scripts of their own and
  found every field of all 50 equal to it.
- `tests/test_hard_cases.py`: `test_hard_case`, parametrized by the 50 cases, checks the oracle's
  possible genotypes of each file and the one row of `compare()`, field by field. The two excerpt
  tests, and a test that they report an STR region across a locus's end.
- `find_possible_genotypes` in `tests/haplotype_oracle.py`, with tests on R10, on S4's file A,
  whose two windows must be taken together, on windows outside the stretch, on invalid samples and
  on each `ValueError`.
- Each check of the test shown failing on a working copy: S1 with an ALT base of file B changed
  fails the oracle's check; S3's outcome and one locus allele of H1, changed, fail `compare()`'s.
  The tests reviewer also reverted the `*` rule (S11 and E5 fail), the span rule (E3, E15 and E16
  fail), and parts of the evaluation (S8, S10, S11, E13, E14, H13 and H15 fail).

**What the review found.** No wrong outcome. Taken: the excerpt tests kept only the STR regions
wholly inside a locus, and now keep those that overlap or touch it; the oracle function's windows
at the stretch's ends and its invalid samples were untested; the literals gave the possible
genotypes as sequences, which a reader could not check against §4.5; the function was named
`find_possible_haplotypes` for what are possible genotypes, so the spec's §5 and the code changed;
the environment skip was a fourth copy, now one helper in `tests/conftest.py`. Not taken:
`tests/test_respell_check.py` keeps its own helper, which fails rather than skips on a wrong path,
as one of its tests asserts.

**What changed in the plan.** The hard cases give 51 tests without the variables, not 50: the
review added the test of an STR region across a locus's end. `tests/compare_case.py`'s FASTA writer
became public, `write_fasta`, for the cases to use.

## For whoever next revises the skills

The owner can stop here.

- The spec of the cases was written from research by four subagents over pvc's reports, and each
  case was then run by a scratch harness, through the oracle and through `compare()`, before the
  spec was written. Running the cases first is what found S11 and E5, where the spec's rule and the
  caller's meaning parted; a spec that only copied the reports' cases with their reported outcomes
  would have missed both. The spec skill says to work cases by hand; for cases from real files, a
  harness that runs both the oracle and the code is cheaper and found more.
- One answer of the owner, the `*` against the REF as written, led the review of milestone 1 to two
  more decisions of the owner, the span and the field's memory. A rule change made before the plan,
  in the session's chat, was specified more thinly than one made in a spec reviewed by the
  `spec-reviewer`; the next such change could go through that reviewer before the code.
- Tokens, from the results of the `Agent` tool: milestone 1, implementation and two rounds of fixes
  192,802; review 43,000 to 60,000 for each of five categories, and 76,810 for the correctness
  reviewer with its second look. Milestone 2, implementation and one round of fixes 210,362;
  review 64,523 to 95,483 for each of four categories.
