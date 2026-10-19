# Report: the STR finder

*19 October 2026. Under way. The work of `doc/plans/str-finder.md`, which builds
`doc/specs/str_finder.md`, on the branch `str-finder`.*

The plan started from commit bec2800, on which `uv sync --locked`, `ruff format --check`,
`ruff check` and `ty check` pass and `pytest -q` gives "213 passed in 1.38s".

## Milestone 1: the segment contract and the candidates

Commits 90a3c39 to 932d773. `src/pvc_eval/segments.py` and `tests/test_segments.py`.

- **The contract**, `RegionType` and `GenomeSegment` as `doc/spec.md` §4 gives them, in commit
  90a3c39 alone, sent to the loci session on 19 October 2026.
- **`StrFinderParameters`**: each field tested at its limit and one past it; each raises
  `ValueError` naming the field, its value and the limit.
- **Step 1**, with the 17 ported tests of pvc's `tandem_repeat.rs`, all passing.
- **The memory bound.** On 4 million random bases, at the default copy floors, step 1 peaks at
  8.0 MB under `tracemalloc`; holding every candidate, in one block, it peaks at 965.6 MB (Apple
  M5 Pro, CPython 3.15t). The test takes about 34 s and is marked `slow`.
- Checks on 932d773: ruff and ty clean, `pytest -q` "301 passed in 37.56s" (213 before).

**Changed in the plan.** Step 2, the contig-end flank, moved to milestone 2: it needs the cut of
step 6, and it only drops candidates, so it can run after step 1.

**What the review found and is fixed.** Five reviewers, one per category.
- *Tests, a blocker.* Four mutants of step 1's drop rule, and one that tests a candidate before its
  divisors are final, changed the output and passed every test, because two tests took their
  expected values from step 1 itself. Six sequences now carry pvc's candidates as literals, from a
  reviewer's transcription of pvc's code, over block sizes 1 to 60, and each mutant fails them.
- *Design.* The Ruzzo–Tompa "segments" were renamed subsequences, so that the word keeps the
  spec's meaning; two functions that nested four levels were split; internal ends made half-open.
- *Errors.* Three messages of the parameter checks corrected.
- The spec and correctness reviewers found nothing. Each compared step 1 with its own
  transcription of pvc's `find_tandem_repeats`: 7,800 random and adversarial sequences, every block
  size from 1 to 4,096, no difference.

**Not taken.** More checks on `GenomeSegment`: §4 states only that the end is after the start, and
the contract is shared. The duplicate default tests: they are pvc's, ported by name.

**What the owner should know.**
- pvc's test `a_period_multiple_re_detection_is_dropped_by_geometry` never reaches the rule it is
  named for, at bbcc9f4: its period-6 stretches are one copy long or have a motif that is not
  primitive. A test that does reach it was added.
- On input made of long runs of one base, the queue of candidates waiting for their test grows by
  about one per 500 bases, because period 1 keeps one subsequence open. The output stays right, and
  pvc has the same bound. Whether real sequence does this is for the performance work.
