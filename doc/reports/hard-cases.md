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

**Checks at the end**, at 6ae1f1a: `ruff format --check`, `ruff check` and `ty check` pass;
`uv run pytest -q` gives "1987 passed, 11 skipped" (1886 passed, 9 skipped at a050b57); the hard
cases give "51 passed, 2 skipped" without the two variables and "53 passed" with them: the 50
cases, the test that the excerpt tests catch an STR region across a locus's end, and the two
excerpt tests, which need the variables. `PVC_EVAL_GRCH38` and `PVC_EVAL_TOMATO` name the FASTA of
GRCh38 and of SL4.0, the tomato reference of the five tomato cases.

**What is asked of the owner.** The merge of `hard-cases` into `main`; nothing is pushed.

- A REF of letters outside ASCII, `ß` against a reference `G`, stopped the run with a Python
  traceback, and `ſ` was read as the IUPAC code `S`, on `main` too. At the owner's word of 21
  October 2026 it is fixed on this branch by `doc/plans/non-ascii-ref.md` (below). An ALT holding a character outside ASCII, which was read as `SS` for `ß`, now stops the run too,
  as the owner decided the same day.
- The warning for a record whose span the span limit cut says it "was given only its REF as its
  span". A reviewer thought it no longer fits a record with a `*`, which keeps its REF as written.
  It does fit: a user reads "its REF" as the REF the file wrote, which is now exactly the span of
  such a record. It is left as it is.

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
  §5.2). It moved the locus of R10, the padded `*` case this branch added to `reconciliation.md`, to 10–19,
  with the same genotype, and the loci of E3, E15 and E16
  by one or two bases, with no outcome changed; the specs were corrected first.
- The design reviewer measured the field on every record at 721 bytes per record against 601
  before, on 100,000 one-sample SNPs; the owner decided it is set only where it differs, 610 bytes per record, as tracemalloc measured
  the records `parse_vcf` held.
- Three reviewers found no test of the padding on the left of the trimmed REF; the tests above
  were added.
- Not taken: the errors reviewer's question whether a `written_ref` with a negative start should
  be refused; the spec asks only that it hold the trimmed REF, as it asks of a span.
- The correctness reviewer's second look at the span fix found it settled, with its tests failing
  without it, and parse_vcf's order kept on 1,500 random records. Its nit, that the span-limit warning's "only its REF as its span" no longer fits a record with a
  `*`, does not hold: the REF the file wrote is that record's span.
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

## The plan of a REF letter outside ASCII

`doc/plans/non-ascii-ref.md`, one milestone, done at a0b9dbc from 178789c. `doc/spec.md` §5.2 says
a letter outside ASCII is no base; `_check_ref` upper-cases only ASCII letters, so such a REF
differs from the reference and stops the run with the error of §5.2, which now writes each letter
outside ASCII with its code point, `'K' (U+212A)`, so that a look-alike of a base cannot read as one.
`tests/test_vcf.py` has six cases: `ß` and `ſ` at a `G`, `ß` inside a longer REF and ending a REF
of 21 bases, the Kelvin sign, and a combining mark after a `G`; each failed at 178789c. Two
reviewers, errors and tests; the errors reviewer tried 15 characters outside ASCII, every one
stopped with the error naming the file, line and position, and its finding, the look-alikes, is
the code points above; the tests reviewer's missing case, the REF of 21 bases, was added.

Its milestone 2, an ALT holding a character outside ASCII, done at 6ae1f1a from f84572f, at the
owner's decision of 21 October 2026. `doc/spec.md` §5.2 and §3 say such an ALT stops the run, and
that a record with a symbolic allele is skipped whatever its other ALTs hold. The error names the
file, the line, the position and each character with its code point; an ALT column over 20
characters is not quoted, and the message says which ALT holds the character and where.
`test_alt_holding_a_character_outside_ascii_stops_the_run` has 12 cases, each accepted before.
Two reviewers, errors and design. Taken: a `str.isascii()` test first, since the check runs on
every record (19 ms on the 98,149 records of chr20 in Python's loop, about 3 ms so); the message
capped as the REF's is, after a 100,000-base insertion gave one of 100 kB; "character" for "letter".
Whole `parse_vcf` on chr20 against GRCh38, best of 5, two runs each: 0.945 s and 0.964 s at
e839edc, 0.961 s and 0.976 s at 6ae1f1a, within the spread of one run, 0.1 to 0.3 s. Not taken:
the design reviewer's wish to report a REF that differs before a damaged ALT on the same line; the
errors reviewer found the order of a line's faults sensible, and such a line is damaged either way.
