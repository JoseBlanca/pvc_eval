# Report: compare

*Under way, 20 October 2026. The work of `doc/plans/compare.md`, which builds
`doc/specs/compare.md`, on the branch `compare`; nothing is merged into `main` or pushed.*

## Before the first milestone

On 1d4a1f8, with `Kind` and `RegionType` already `StrEnum` (ca0463f):
`uv run python -c "import sys; print(sys._is_gil_enabled())"` printed `False`; `uv sync --locked`
was clean; `ruff format --check` gave "40 files already formatted"; `ruff check` and `ty check`
gave "All checks passed!"; `pytest -q` gave "914 passed, 4 skipped".
