# Report: reading only the contigs `--regions` names

*Done, 20 October 2026. The work of `doc/plans/regions-read.md`, which builds
`doc/specs/regions_read.md`, on the branch `regions-read`, last commit of code 79e968b. Nothing is
merged into `main` or pushed. The owner asked for the spec to be reviewed and built without their
approval; nothing here waited on them.*

## For the owner

**With `--regions`, `pvc_eval compare` now reads only the records of the chromosomes the regions
name.** GIAB's HG002 v4.2.1 against itself over chr20 of GRCh38, `--regions` naming chr20, takes
25.7 s, against 39.1 s before, and gives the same 81,257 rows and the same summary, all *same
genotype* (Apple M5 Pro, CPython 3.15t, 20 October 2026, at 79e968b). Its peak was 503 MB
resident, against 400 to 464 MB measured before the change on different runs; nobody has looked
into why it is not lower. The 25.7 s are now the STR finder's scan of chr20 and the loci of chr20;
the other 21 chromosomes cost the 4 s of skipping their lines.

The checks on 79e968b: `uv sync --locked` clean, ruff and ty pass, `pytest -q` "1372 passed, 5
skipped" three runs in a row (1306 before), and the HG002 test "19 passed".

**Decided in your name, against your answer.** You chose to read whole chromosomes by seeking
through the VCF's index, and to read a VCF with no index whole, with a message. Measured first,
skipping the lines of other chromosomes after looking at their first columns costs 4 s on HG002
and needs no index: it reads whole chromosomes, as you chose, keeps exact line numbers in error
messages, which seeking would lose, and speeds up a VCF without an index too. So there is no index
reader of compressed blocks and no message. Each skipped line is still checked for its columns, its
POS, its contig and its order; its REF, sample columns, ALT, QUAL and GT are not. So a record on
chr1 whose REF shows a VCF made on another build no longer stops a chr20 run; the first record of
chr20 still would.

**Asked of you:** whether to keep this design, and the merge of `regions-read` into `main`.
`main` has not moved since the branch was made from it (7207b8e).

**Two defects of the `compare` work found and fixed here**, since this branch changes the same
files:

- **A hang.** A `parse_vcf` stream closed by the garbage collector waited for its thread pool from
  inside the collector, and when the collector ran in a thread that held Python's list of threads,
  as it does while starting or listing threads, both waited forever. It hung one test run for 10
  minutes, and is the likely cause of the failure a reviewer saw once in 15 runs during the
  `compare` work. `compare()` had the same pattern in its stop. Now a stream closed early does not
  wait for its pool, and `compare()`'s stop waits for no thread: each reader and each scanning
  thread opens its own handle on the reference (spec of `compare` changed in 4cdb41d).
- **An error printed at exit**, 2 runs in 120, when a `compare()` stream was left open as Python
  exited. The stop now does nothing while Python exits.

After these, a stress test of 8,884 threaded runs with random early closes and streams dropped
unclosed gave no hang, no wrong row, nothing printed and no handle left open; the suite passed 20
runs in a row.

**Left open:** with `--regions`, the last records of the last chromosome read reach the comparison
about 3.6 s late, while the rest of the file is skipped; the run's total time and the table do not
change. And the time to the first row is not measured.

## The milestone

Built in b917c7f, fixed in f710ce0 to 79e968b. Reviewed in spec, tests, and correctness with
errors and design. The spec reviewer found the four checks of a skipped line, their order, their
messages and their line numbers equal to today's on 23 cases. The correctness reviewer found the
records read equal to a whole read filtered to the same contigs on 240 random cases. Fixed:
`contigs` given as one string read nothing silently (now a `TypeError` when the stream starts); no
test checked the earlier line number of the out-of-order message; then the two defects above, and
handles on the reference left open when a thread failed to start. The spec's case of a malformed
GT on a skipped line could not fail, since a malformed GT is no error on a contig read either; the
spec says so now (f619ad8). Not taken: none.

## For whoever next revises the skills

- **Measure the cheap design before writing the agreed one.** The owner's answers assumed a reader
  that seeks through the index. A two-minute measurement of skipping lines by their first columns
  showed it simpler and as fast in practice, and changed the spec before any code.
- **A defect a reviewer saw once and could not name was real.** The one failure in 15 runs during
  the `compare` work, left as unexplained, was a deadlock between the garbage collector and
  threads. A failure seen once in code with threads is worth a loop with `-rf` and
  `faulthandler` before it is written off; `python.md` could say that a finalizer, or anything the
  garbage collector may run, never joins a thread.
