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
