# Report: compare

*Under way, 20 October 2026. The work of `doc/plans/compare.md`, which builds
`doc/specs/compare.md`, on the branch `compare`; nothing is merged into `main` or pushed.*

## Before the first milestone

On 1d4a1f8, with `Kind` and `RegionType` already `StrEnum` (ca0463f):
`uv run python -c "import sys; print(sys._is_gil_enabled())"` printed `False`; `uv sync --locked`
was clean; `ruff format --check` gave "40 files already formatted"; `ruff check` and `ty check`
gave "All checks passed!"; `pytest -q` gave "914 passed, 4 skipped".

## Milestone 1: `compare()` with no thread

Built in 37f033b and fixed in 02cb2b9 to da3f025. `compare()`, `CompareConfig` and `PVC_STRATA`
are in `src/pvc_eval/compare.py`; the four-contig case is in `tests/compare_case.py`, its STR
regions checked against `segment_contig` (chr1 100–116 AC and 200–218 AGT, chr2 120–140 ATTC, chr3
100–120 CA, 0-based; none on chr4). On da3f025 `uv run pytest tests/test_compare.py -q` gives "50
passed", the plan asked for at least 15; the five checks pass, with "964 passed, 4 skipped".

- The rows equal those of the nine steps, written out in the test with `segment_contig` over every
  contig with records, with no option and with each option and setting, `--regions` cut inside a
  contig among them. Breaking the segment iterator seven ways fails tests each time (the
  implementer's five and two of the reviewers').
- Only the contigs with records, and with `--regions` only those it names, are scanned, by the
  INFO lines.
- The errors: a VCF out of order and an untyped `--str-regions` region after the nine steps' rows;
  a FASTA byte that is not a letter, raised by the finder, after literal rows.

**The spec changed twice, by the session, after the review.** Both are in commits of their own,
0ad2d19 and eef6533, and neither changes a row of the table.

- Before an error of the STR finder, `compare()` can give one locus more or fewer than a chain of
  `segment_contig` over a list of contigs known before the reading. That chain reads a contig's
  segments at another record than the spec's iterator, which decides a contig only from the
  records `build_loci` has taken. Three reviewers found it, with cases. The spec now names
  `compare()` with `finder_threads=0` as the stages run one after another, which the threads of
  milestone 5 must match. The command line writes no table on an error, so its user sees the same
  message either way.
- `compare()` no longer repeats the stages' checks of the ploidies, the phasing cap, the span limit
  and the chain distance, which were copies with the same messages. The stages raise the same
  `ValueError` before any row, but after an error of a header or a BED.

Also added to the spec: a contig whose records all have no ALT gets no progress line, and the
loci counted as skipped on a contig `--regions` does not name are built without STR regions, so
there can be more of them (be57c54).

**What the review found that mattered, all fixed:** `--regions` cut inside a contig had no test,
so a copy of `select_loci`'s rule could drift (blocker; a second look confirmed the fix); six
settings that no test set; no test that a value at the edge of its range is accepted; no record
on the last contig. Not taken: none.

The implementer added the verb "note", for a pass-through that records something about the items
it hands on, to the table of `python.md`, as that file asks.
