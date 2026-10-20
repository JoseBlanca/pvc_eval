# Plan: compare

*Status: approved by the owner on 20 October 2026, who asked for it to be written and run; done on 20 October 2026. It
builds `doc/specs/compare.md`, approved by the owner on 20 October 2026, on the branch `compare`
in `.claude/worktrees/compare`; the owner's handoff named the branch `compare`, so it is not
`plan/compare`. The work report is `doc/reports/compare.md`.*

The plan joins the stages that exist on `main` into `compare()`, writes the table and the
summary, adds the `pvc_eval compare` command, and then runs the reading of the VCFs and the STR
finder in threads. Six milestones, in the order the owner agreed on 20 October 2026: `compare()`
with no thread first, so that the threads of milestone 5 have rows to be compared with; the
table, the command and the index reader next; the threads; and the check over real data last.

## In and out

In: everything `doc/specs/compare.md` gives, in a new `pvc_eval.compare` and `pvc_eval.table`,
the `compare` command in `pvc_eval.cli`, and their tests. The spec's seven open points were
answered by the owner on 20 October 2026, and its section "Decided by the owner" lists them. The
changes to `doc/spec.md`, `doc/specs/loci.md` and `Kind`, `RegionType` as `StrEnum` are done,
in c81328d, 12df04f and ca0463f.

Out, and where it goes:

- the re-spellings of `doc/spec.md` §7.4, and pvc's hard cases of §7.5: the next pieces;
- the other evaluators of §8, and `--evaluator`;
- the STR finder's speed: a performance review on the branch `perf/str-finder`, not merged.
  This branch has `main`'s finder, so milestone 6's times are for that one, five to nine times
  slower than the reviewed one (the spec's "How it runs").

## Before the first milestone

- The branch `compare` holds the approved spec, at 1d4a1f8 or after.
- `uv run python -c "import sys; print(sys._is_gil_enabled())"` prints `False`.
- The five checks of `python.md` pass. On ca0463f: `uv sync --locked` clean, "40 files already
  formatted", "All checks passed!" from ruff and from ty, and "914 passed, 4 skipped".
- bcftools 1.24 is at `/opt/homebrew/bin/bcftools`, for the index files of milestone 4.

## The test case shared by the milestones

Milestone 1 writes, in `tests/conftest.py` or a helper module of the tests, a function that
writes to `tmp_path` the four-contig case of the spec's "How it is verified": a FASTA with its
`.fai`, four contigs of a few hundred bases, with STR regions that `segment_contig` admits on
contigs 1 to 3 and none on contig 4, and VCF text for file A (records on contigs 1 and 3) and
file B (2 and 3), with two samples in each, one sample pair by name. Milestones 3, 4 and 5 use it.
Its STR regions are checked once, in a test, against `segment_contig`'s output written as
literals.

## Milestones

Milestones 1, 2 and 4 need nothing of each other. They are run one after another all the same,
since their subagents would commit to the one branch of this worktree at once.

- [x] 1. `compare()` with no thread

**What it gives.** `compare()` and `CompareConfig`, which a Python user can call: the rows of two
VCFs, from the stages run one after another in the caller's thread.

**The spec it builds.** `doc/specs/compare.md`: "The stages, joined"; "Which contigs are scanned",
the rule of which contigs need scanning, not the planned contigs; "How the segments reach
`build_loci`", whole; "Errors, and stopping early", the first paragraph; "The Python interface"
for `PVC_STRATA`, `CompareConfig`, `DEFAULT_FINDER_THREADS` and `compare`; and from "The `compare`
command", the two INFO lines that `compare()` logs. In this milestone every value of
`finder_threads` runs with no thread, as 0 does, so that the default, 2, works before the threads
exist; milestone 5 gives the values from 1 up their threads, and its tests compare them with 0.

**Deliverables.**

- `src/pvc_eval/compare.py` and `tests/test_compare.py`.
- At `compare()`, on the four-contig case: its rows equal those of the nine steps of "The stages,
  joined" written out in the test, with no option, with `--regions` naming contig 3 alone, with
  `--strata pvc`, with a `--strata` BED, and with `--str-regions`; one locus's rows as literal
  `OutcomeRow` values.
- At `compare()`, which contigs are scanned, from the INFO lines "finding STR regions": the
  contigs with records and not contig 4; with `--regions` naming contig 3, contig 3 alone.
- The one-base GENERIC segment: a record past the last STR region of contig 1, and one whose
  span ends on contig 1's last base, each followed by a record on contig 3, give the nine steps'
  rows.
- The errors: a VCF out of order on contig 3; a FASTA byte that is not a letter on a contig of
  more than 65,536 bases, more than 65,536 bases from any record, raised by the finder; a
  `--strata pvc` with an untyped `--str-regions` region; each with the same type and message as
  the nine steps. The rows before the first and the third are the nine steps'; before the finder's
  error, they are literal rows, since a list of contigs known before the reading can read a
  contig's segments one locus earlier or later (the spec's "The stages, joined", changed after the
  review of this milestone). A `ValueError` for each setting out of range and for a
  `finder_threads` under 0.
- The INFO line when a contig's loci are done, with its counts, matched by pattern.
- `uv run pytest tests/test_compare.py -q` runs at least 15 tests and all pass; the five checks
  pass. On 1d4a1f8 the file does not exist and pytest exits with 4.

**What it stands on.** The stages on `main`; nothing of this plan.

**What could go wrong.** The segment iterator's rule is the part a wrong reading breaks silently:
a contig that should get its STR regions gets none, and its loci come out narrower with plausible
rows. The rows-equal-to-the-nine-steps tests guard it, so in them the nine steps feed `build_loci`
the segments of every contig with records, by `segment_contig` over each, not by the iterator.

- [x] 2. The table and the summary

**What it gives.** `pvc_eval.table`: `TABLE_COLUMNS`, `SUMMARY_COLUMNS`, `format_table_row` and
`write_table`, which writes the table and the wide summary from any iterator of `OutcomeRow`.

**The spec it builds.** `doc/specs/compare.md`: "The table", "The summary", "Writing the files",
and "The Python interface" for these four names.

**Deliverables.**

- `src/pvc_eval/table.py` and `tests/test_table.py`.
- At `format_table_row`, on `OutcomeRow` values built in the test, each field as a literal: the
  cases listed in the spec's "The table, at `format_table_row`", and the worked row of "The table"
  as its exact line.
- At `write_table`: the header and the rows of a short list, with and without `all_rows`; a `.gz`
  path read back with `gzip`; rows whose iterator raises after two rows, leaving no output, no
  `.tmp`, and an older file at the path unchanged; each refusal of "Writing the files" with its
  message; the summary over two sample pairs and two strata, lines, order and counts as literals,
  a 0 among them, *reference in both* counted without `all_rows`; a summary without `--strata`.
- `uv run pytest tests/test_table.py -q` runs at least 12 tests and all pass; the five checks pass.

**What it stands on.** `OutcomeRow`, on `main`. Not on milestone 1.

**What could go wrong.** The 1-based coordinates of `start` and of the STR regions, and the
QUAL's shortest form: a test of each at a boundary, a QUAL of 50.0, 29.53, 1e-05 and 1e+20.

- [x] 3. The `pvc_eval compare` command

**What it gives.** pvc_eval runs end to end from the command line: two VCFs in, a table and a
summary out, with no thread yet.

**The spec it builds.** `doc/specs/compare.md`, "The `compare` command", whole.

**Deliverables.**

- `compare` in `src/pvc_eval/cli.py`, its tests in `tests/test_cli_compare.py`.
- At `main()`, on the four-contig case: the table equals what `write_table` writes for
  `compare()`'s rows; `--summary` and `--all-rows`; `--strata pvc` and `--strata ./pvc`, a BED
  named `pvc`, told apart; each option's value out of range, with its message and status 1; a
  value that is not a number, status 2; a `PvcEvalInputError` from a stage, its message on stderr
  with no traceback, status 1, and no table left; the INFO lines on stderr, by pattern.
- Ctrl-C: `KeyboardInterrupt` raised from inside the rows, through a stand-in for `compare` in
  `pvc_eval.cli`, gives status 130, "Stopped; nothing was written." and no `.tmp`.
- `uv run pvc_eval compare --help` names every option with its default.
- `uv run pytest tests/test_cli_compare.py -q` runs at least 10 tests and all pass; the five checks
  pass.

**What it stands on.** Milestones 1 and 2.

- [x] 4. The index reader

**What it gives.** The contigs a `.tbi` or `.csi` lists, read with `gzip` and `struct`, for
milestone 5 to plan its scans.

**The spec it builds.** `doc/specs/compare.md`, "Which contigs are scanned": the planned contigs,
an index that cannot be read and its warning, both a `.tbi` and a `.csi`, contigs the reference
lacks; and "The index reader" of "How it is verified".

**Deliverables.**

- The reader in a module of its own, `src/pvc_eval/_vcf_index.py`, public only to pvc_eval's
  modules as `src/pvc_eval/_text_files.py` is. The spec says `pvc_eval.compare`; a file of its own
  keeps it testable apart and out of the file milestones 1 and 5 change.
- `tests/data/`: a small bgzipped VCF with its `.tbi` and `.csi`, made with `bcftools view -Oz`,
  `bcftools index -t` and `bcftools index -c` (1.24), the commands and what the files hold in
  `tests/data/README.md`, its header naming contigs it has no record on.
- Tests: the contigs of the `.tbi`, of the `.csi`, and of both together as literal lists; a VCF
  with no index and a plain-text VCF give no contigs; an index cut short and one that is not an
  index give the warning with its text, and no contigs.
- `uv run pytest tests/test_vcf_index.py -q` runs at least 6 tests and all pass; the five checks
  pass.

**What it stands on.** Nothing of this plan.

**What could go wrong.** The `.csi` layout: its contig names sit inside its auxiliary block, after
a 28-byte header, which the spec's reviewer and the spec's writer checked only on bcftools 1.24's
files.

- [x] 5. The threads

**What it gives.** The default run: each VCF read in a thread of its own, the STR finder on a pool
of `finder_threads` threads, planned from the indexes, a few contigs ahead of the loci.

**The spec it builds.** `doc/specs/compare.md`: "Which thread does what"; "Which contigs are
scanned", the planned contigs and "How far ahead"; "Errors, and stopping early", whole; and "How
it runs", the bounds of "What is held".

**Deliverables.**

- The threads in `src/pvc_eval/compare.py`; their tests in `tests/test_compare_threads.py`.
- At `compare()`, on the four-contig case: the rows with `finder_threads` 0, 1, 2 and 4 equal each
  other and the nine steps', with both files indexed, one, and neither, with `--regions`, and with
  `--str-regions`; an index that lists contig 4 with no record on it, and one that leaves out
  contig 2, give the same rows.
- The contigs an index lists that the reference lacks are dropped from the planned contigs here,
  since the reader of milestone 4 does not take the reference: an index that lists a contig the
  FASTA lacks gives the same rows as no index, and no traceback.
- The errors of milestone 1 again with 1, 2 and 4 threads: the same type and message after the same
  rows. A FASTA error on a contig no record is on, listed by an index, is not raised.
- The bound and the stop, as the spec's "The bound and the stop" gives them: between 10,000 and
  10,000 plus two batches taken from a held reader; after a close, no reader thread alive within
  5 seconds and no scan started after it. No exception reaches `threading.excepthook` in any of
  these tests: a test sets a hook that fails it.
- From the review of milestone 3: `compare()` is annotated as returning a `Generator`, which the
  command closes, and the `cast` in `pvc_eval.cli` goes; a test at `main()` shows the command runs
  `compare()` with two finder threads.
- Each test that waits on a thread has a timeout, so that a hang fails it.
- `uv run pytest tests/test_compare_threads.py -q` runs at least 12 tests and all pass, in each of
  ten runs in a row, made by a shell loop; the five checks pass.

**What it stands on.** Milestones 1 and 4.

**What could go wrong.** A hang or a race. The spec's reviewer found no cycle among the readers'
buffers, the pool's limit and the loci side; the tests of this milestone are where one would
show, which is why each has a timeout and the file is run ten times. The other risk is a row that
differs only with threads: the rows-equal tests guard it.

- [x] 6. A file against itself over real data

**What it gives.** The check of `doc/spec.md` §7.3 over chr20 of GIAB's HG002, and the time and
peak memory of a real run, for the report.

**The spec it builds.** `doc/specs/compare.md`, "A file against itself over real data", and the
owner's answer to point 6.

**Deliverables.**

- `tests/test_compare_real.py`, marked slow, reading `PVC_EVAL_GRCH38` and `PVC_EVAL_HG002_VCF` and
  skipped with a message naming each that is missing. It runs `compare()` on the VCF against
  itself with `--regions` naming chr20 and checks that every row's outcome is *same genotype*,
  *reference in both*, *missing in both*, *invalid in both* or *too complex*, the last three only
  where the file has a missing haplotype, an invalid sample or a sample over the phasing cap.
- Run with the two variables set to the paths of the spec's point 6: it passes, and the report
  gives its wall time, its peak memory (`/usr/bin/time -l`) and the count of rows by outcome.
- Without the variables, `uv run pytest -q` skips it with its message.

**What it stands on.** Milestones 3 and 5.

**What could go wrong.** A row that breaks the property is a finding about the stages, not about
`compare`. It changes an outcome, so it is a stop: the session asks the owner, with the locus,
its records and the row, after finishing what does not depend on the answer.

## The check of the whole plan

Milestone 6 is that check, with the command `pvc_eval compare` run once by the session on HG002
against itself over chr20, `--summary` included, and its summary in the report.
