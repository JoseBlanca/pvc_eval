# Report: str-regions on threads

*20 October 2026. Done. The work of `doc/plans/str-regions-threads.md`, on the branch
`plan/str-regions-threads`, from main at 38642f1; nothing is merged into main or pushed.*

**`pvc_eval str-regions --threads N` segments up to N contigs at once, and its BED is the same at
any N.** On GRCh38's chromosomes 13 to 22 and Y (855 Mb, as one FASTA), on an Apple M5 Pro with
CPython 3.15t, one run each, every BED equal to the 1-thread one:

| threads | wall time | peak RSS |
|---|---|---|
| 1 | 250.9 s | 919 MB |
| 2 | 142.1 s | 1,404 MB |
| 4 (the default) | 84.1 s | 2,351 MB |
| 8 | 60.5 s | 3,354 MB |

Other sessions' test runs used up to two cores during the 1-thread run, and the load average rose
from 1.3 to 6.7 over the four runs, so these locate the gain rather than measure it to the second.

What is asked of the owner:
- **The merge** of `plan/str-regions-threads` into main.
- **The default.** It is 4, as the spec says, the same as the VCF reader's. By the plan's rule, 8
  is proposed: on the input above it is 28% faster than 4, at a peak of 3.35 GB, under the plan's
  4 GB. But the input holds the smaller chromosomes. Over all of GRCh38, eight threads would hold
  chr1 to chr8 together, each 1.0 to 1.2 GB at its peak, so the peak would be near 9 GB, against
  near 5 GB at 4; neither is measured. The session recommends keeping 4, which a user raises with
  `--threads` when the memory is there.

What exists now that did not:
- `write_str_regions(reference, output, parameters, regions, num_threads)`: the contigs, in the
  order of the `.fai`, go to a `ThreadPoolExecutor` through `Executor.map(..., buffersize=2 * N)`,
  each thread returning the formatted rows of its contig; with 1 there is no pool. `num_threads`
  outside 1 to 256 raises `ValueError` before anything is read or written.
- `--threads N` in `pvc_eval str-regions`, with its `--help` line (the range, the default, the
  memory per Mb, the wait on Ctrl-C). A value outside 1 to 256 stops with exit status 1 and a
  message naming `--threads`; a value that is not a number gets argparse's message and status 2.
- `DEFAULT_NUM_THREADS = 4` and `MAX_NUM_THREADS = 256`, public in `pvc_eval.segments`.
- 17 tests, 934 to 951; ruff and ty clean.

## Milestone 1

Commits ea0a1b2 to 7d01ea0, with the spec and plan in a5dec6c, c86b80e and 30d2ac9.

**Changed in the spec before the code.** The spec reviewer and the first reader found the first
draft left `buffersize`, the exit status of `--threads 0`, the test of the order and the timing of
an error under threads to be guessed; the spec now says each. After the code review, the spec
says that Ctrl-C can take minutes on human chromosomes, not a minute, and names a test for each
rule of the pool.

**What the review found and is fixed.** Five reviewers, one per category. The correctness and
errors reviewers found no defect: the BED was the same at 1, 2, 3, 8 and 256 threads with and
without `--regions`, every `--threads` value gave the spec's message and status, and an error, a
failed write and Ctrl-C each left no `.tmp` and no threads behind.
- *Tests, two blockers from two reviewers.* Nothing could see N: five mutations (no pool ever, a
  pool at 1, a pool of 4 whatever N, the command line ignoring `--threads`, no `buffersize`) left
  every test passing. Three tests now watch `segment_contig` with `sys.monitoring`, which records
  its calls and replaces nothing: at 1 thread every contig is segmented in the calling thread; at
  2, none is, on at most 2 threads; at 2, at most 1 + 2N = 5 contigs have started when a long
  first contig returns. Each mutation fails one of them.
- *Tests.* The test of an error in a later contig could not tell the `.fai`'s order from the order
  contigs fail; its first faulty contig is now 30,000 bases with the bad byte at its end, so that it
  fails last (10,000 already caught the mutation in 230 runs of 230). No test accepted 1 or 256;
  both are now accepted, at the command line and in the function.
- *Design.* The limit of 256 was written in `cli.py` and `segments.py`; it is now one public name,
  imported by the command line, which still checks `--threads` itself so that a `ValueError` from
  the function stays a defect. The new whole-number check was replaced by the module's own.

The 17 new tests passed 10 runs out of 10 after the fixes.

**What the owner should know.**
- Ctrl-C prints Python's traceback, with or without threads, as before this change. A one-line
  message would be a change of the spec.
- After a second Ctrl-C the command leaves `write_str_regions` while worker threads still run,
  and the reference is closed under them. `IndexedFasta._read_upper_bytes` reads `self._fd` twice,
  so a worker could read a file descriptor that has been freed; the correctness reviewer could not
  make a wrong read happen, since the process is shutting down by then. `fasta.py` is outside this
  plan.
- `vcf.py` keeps its own copy of the default, the limit and the ordered map. The design reviewer
  proposes one private module for them when a third user comes, `compare`'s threads.
