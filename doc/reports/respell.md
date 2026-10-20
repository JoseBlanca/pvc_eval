# Work report: the re-spelling program

*Under way, 20 October 2026. The plan is `doc/plans/respell.md`, built from
`doc/specs/respell.md`, on the branch `respell`. Nothing is merged into `main` or pushed.*

The checks on the commit the work starts from, 9dbec00: `uv sync --locked` clean, `ruff format
--check`, `ruff check` and `ty check` pass, `uv run pytest -q` gives "914 passed, 4 skipped"; the
skipped tests need pvc and GRCh38.
