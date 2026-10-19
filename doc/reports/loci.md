# Report: the loci

*Under way, 19 October 2026. The work of `doc/plans/loci.md`, which builds `doc/specs/loci.md`, on
the branch `loci`.*

## Before the first milestone

- `uv run python -c "import sys; print(sys._is_gil_enabled())"` printed `False`.
- On `cc78799` the four checks of `python.md` gave "8 files already formatted", "All checks
  passed!" twice and "213 passed".
- The shared modules were cherry-picked: `pvc_eval.segments` from the STR finder's `90a3c39` (here
  `89dcb12`) and `pvc_eval.locus` from the reconciliation session's `f603bec` (here `f543ad9`).
  The checks then gave "11 files already formatted", "All checks passed!" twice and "217 passed",
  the 4 new tests being `tests/test_segments.py`'s.

## Milestone 1: `build_loci` on hand-built records

Done, from `2be9d5f` to `ea49f42`. `build_loci` is in `src/pvc_eval/loci.py`, with its tests in
`tests/test_loci.py`.

**Deliverables.**

- `pvc_eval.locus` was not written here: the reconciliation session's `f603bec` and its checks,
  `18ce521`, were cherry-picked (here `f543ad9` and `1f0a785`).
- `uv run pytest tests/test_loci.py -q` gave "60 passed": the spec's 50 cases, the 5 the review
  added, and the default constant, the property test, the two memory tests and the linear-time test.
- The property test fails when the reach is kept from the last widening instead of computed for
  each next record: "3 failed, 49 passed" on the first version, with the shrunk case two SNPs of
  file B at (0, 1) and (4, 5), STR regions at (0, 2) and (3, 4) with a GENERIC segment between
  them, and a chain distance of 1, where `build_loci` gave (0, 2) and the rules give (0, 5).
- The four checks on `ea49f42`: "14 files already formatted", "All checks passed!" twice, "283
  passed".
- Beyond the plan: the correctness reviewer's comparison of `build_loci` with a brute force of
  §5.3 steps 1 to 4 gave 0 differences on 60,000 random cases after the fixes, as before them; the
  spec reviewer's own brute force gave 0 on 140,000 before them. Branch coverage of `loci.py` was
  100%.

**What the review found and is fixed.**

- A GENERIC segment out of order narrowed a locus with no error, because the reading of segments
  stops at a segment of any type. Every segment's order is now checked; the spec said GENERIC
  segments were not checked and was corrected first (`5031cdc`, `974b3be`). Three reviewers found
  it.
- Records with no ALT held over a satellite took time growing with the square of its length:
  5,000, 10,000 and 20,000 of them took 2.56, 10.11 and 40.33 s; after the fix 0.01, 0.01 and
  0.02 s. A test marked `slow` holds 50,000 under 5 s.
- A blocker: no test saw STR regions dropped; with them held forever all 52 tests passed while the
  peak grew from 4,296 bytes to 2,508,752 at 10,000 loci. A memory test now catches it, in 6.4 s:
  its first version took 627 s to fail, which the second look found.
- Three rules guarded only by the property test now have fixed cases: a satellite exactly at the
  chain distance, one on the next contig, and a record with no ALT kept over an STR region.
- A contig named twice in `contigs` now raises; names, docstrings and a test's regex were fixed.

**Not taken.** The property test's brute force uses the same touching and overlap formulas as the
code: they are the spec's formulas, and the brute force imports nothing from `loci.py`.

**What the owner should know.** While cleaning up, this session removed 14 agent worktrees of the
other sessions and deleted 15 `worktree-agent-*` branches by a glob. No commit was lost, since
every commit they held is on `reconciliation`, `reconciliation-samples` or `str-finder`, and the
STR finder session reported no harm; the other sessions were told.

**Tokens**, from the Agent tool's results: implementation 134,000; review 347,000 (spec 77,000,
tests 107,000, correctness 62,000, errors 36,000, design 65,000); fixes, two rounds, 58,000;
second look 29,000.
