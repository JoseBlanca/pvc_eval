# Report: pvc's hard cases

*Under way, 21 October 2026. The work of `doc/plans/hard-cases.md` on the branch `hard-cases`.*

The plan starts on commit a050b57, where `ruff format --check`, `ruff check` and `ty check` pass
and `uv run pytest -q` gives "1886 passed, 9 skipped".

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
- The tests reviewer noted that no property test compares `reconcile` and the oracle on records
  with a `*`; the two read one field, and the R10 tests guard it. Left as it is.

**For the owner.** Outside this plan, the errors reviewer found that a REF of non-ASCII letters,
`ß` against a reference `G`, stops the run with a Python traceback from `_check_ref` in
`pvc_eval.vcf` instead of an error naming the line, and that `ſ` is read as the IUPAC code `S`.
It is on `main` too, and is not fixed here.
