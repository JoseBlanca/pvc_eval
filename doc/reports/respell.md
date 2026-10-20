# Work report: the re-spelling program

*Under way, 20 October 2026. The plan is `doc/plans/respell.md`, built from
`doc/specs/respell.md`, on the branch `respell`. Nothing is merged into `main` or pushed.*

The checks on the commit the work starts from, 9dbec00: `uv sync --locked` clean, `ruff format
--check`, `ruff check` and `ty check` pass, `uv run pytest -q` gives "914 passed, 4 skipped"; the
skipped tests need pvc and GRCh38.

## Milestone 1: decomposition, and what the four re-spellings share

Done in b957dea and 52afc56, and fixed in a3b208f to 26f8b51. `respell_by_decomposing` writes a
VCF with each record nothing touches decomposed into one-base and one-indel records, and
`tests/respell.py` holds what the other three re-spellings will use: the reading, the records
kept, the touched test with the phase, the reach, the writing, new phase sets, the header, the
order and a BGZF writer. The reconciliation property makes its file B with it, and
`_make_decomposed_lines` is gone.

- `uv run pytest tests/test_respell.py` gives "82 passed": D1 to D12, every row of the spec's
  §2.2, the header, bgzip, the output that is the input, the order, and the reach against
  `parse_vcf`'s spans. The plan asked for at least 25.
- `uv run pytest tests/test_reconcile_property.py` passes, with file B made by
  `respell_by_decomposing`; `grep -rn _make_decomposed_lines tests` finds nothing.
- The 20,000-SNP timing test, with 2,600 touching records, takes 0.23 s against a limit of 10 s;
  a sweep that keeps every earlier record takes 133 s and fails it.
- The five checks pass on 26f8b51; the whole suite gives "984 passed, 4 skipped".

Changed in the spec: §2.3 asked for a `Type=String` PS line when an input PS is not a whole
number, and also that the output give bcftools no message the input does not. bcftools 1.24
warns about a PS declared `String`, so the two could not both hold; two reviewers found it. The
spec now adds no PS line then, which leaves the output with the input's own warnings (6294b8e).
It changes only the header of the test program's output, nothing pvc_eval gives.

Changed in the plan: the case notation and the D cases are in `tests/respell_cases.py`, so that
milestones 2 and 3 can run them; `tests/conftest.py` has one VCF writer for the cases.

The review found one blocker: an output named as the input in another case, on macOS, or a hard
link to it, overwrote the input instead of raising. Fixed, with two tests. The correctness and
spec reviewers each compared input and output through `reconcile` over generated files, phased
and not, ploidy 1 and 2, records at both ends of a contig: 16,000 files in all, with no
difference. Fixed besides: seven plausible wrong rules no test could fail on (the touch on a
second sample, a phased GT with no PS against an unphased one, the records next to a `*` on
either side, a slide past one read of the reference, the order of a kept and a rewritten record
at one POS, contigs not in name order), each now a case; the test of the records near a `*` now
takes what any re-spelling writes; each line read is checked against the record matched to it;
names from the verb table. Not taken: a PS above 2^31, which bcftools sets to missing, and one
of more than 4,300 digits, which Python will not convert; no real file has either.

For the owner: a gzip file with one byte changed is reported by `parse_vcf` as having no
`#CHROM` line after its `##` lines, not as damaged. It is outside this plan.
