# Plan: str-regions on threads

*20 October 2026. Approved by the owner's order of that day ("put the threads in"). One milestone.*

**What and why.** `pvc_eval str-regions` segments the contigs of a reference one after another;
over GRCh38 that takes 17.1 minutes on one core (`doc/reports/perf-str-finder.md`). The contigs
share nothing, and `segment_contig` is safe to call from several threads, so this change segments
up to `--threads N` contigs at once and writes their rows in the order of the `.fai`, as
`doc/specs/str_finder.md`, "The `str-regions` command", **Threads**, now says. That section and the
spec's "How it is verified" are what this milestone builds.

## Milestone 1: `--threads` in `str-regions`

Builds: `write_str_regions(reference, output, parameters, regions, num_threads)` in
`src/pvc_eval/segments.py`, and the option `--threads N` in `src/pvc_eval/cli.py`, with its
`--help` line (what it does, default 4, and that memory grows with N), and their tests in
`tests/test_str_regions.py`.

How: a worker function that segments one contig with `segment_contig` and returns the formatted
rows of its STR regions, as `write_str_regions` makes them today (`_find_requested_str_regions`
over the contig's requested spans, then `_format_str_region`). The contigs are mapped over it in
the order of the `.fai` (only those `--regions` names): in this thread when `num_threads` is 1, and
otherwise with `ThreadPoolExecutor(max_workers=num_threads)` and
`pool.map(worker, contigs, buffersize=2 * num_threads)`. `buffersize` is an argument of
`concurrent.futures.Executor.map` since Python 3.14: map then submits at most that many tasks ahead
of the result being read, instead of all of them. The rows are written as each contig's list comes
back. `src/pvc_eval/vcf.py` has the same pattern in `_map_in_order`; it is private there and
`vcf.py` is not changed, so `segments.py` has its own few lines. `num_threads` is a required
parameter of `write_str_regions`, after `regions`; the existing callers in the tests pass it.
`segments.py` gets its own `DEFAULT_NUM_THREADS = 4` (public, for the command line's default) and a
private limit of 256, each with a comment that they are the VCF reader's, `doc/spec.md` and
`vcf.py`'s `DEFAULT_NUM_THREADS`. `write_str_regions` raises `ValueError` for a `num_threads` that is
not a whole number from 1 to 256, before it reads the `--regions` BED or opens the `.tmp`; the
command line turns `--threads` outside 1 to 256 into a `PvcEvalInputError` with exit status 1, as
it does for the other parameters, naming `--threads` (read how `cli.py` does that for them).

Deliverables, each checked:
- the tests of the spec's "How it is verified" for `--threads`, which fail before the change, and the
  whole suite passing, with the checks of `.claude/skills/coding/python.md` (with `PVC_EVAL_PVC` and
  `PVC_EVAL_GRCH38` set, so that `tests/test_whole_port.py` runs);
- the BEDs of `/Users/jose/devel/pvc_eval/.claude/worktrees/perf-str-finder/benchmarks/data/chr1_10mb_4mb.fa`
  and `chr1_10mb_16mb.fa` beside it equal, under `cmp`, to `base/chr1_10mb_4mb.bed` and
  `base/chr1_10mb_16mb.bed` there, at 1 and at 4 threads (one contig each: this checks the code path,
  not the pool);
- a measurement, reported, not committed: wall time and peak RSS of `str-regions` over GRCh38's
  chromosomes 13 to 22 and Y (about 1 Gb, written as a FASTA with `benchmarks/make_str_slices.py`)
  at 1, 2, 4 and 8 threads, the BED compared with the 1-thread one each time. This also shows that
  N > 1 runs contigs at once. The default of 4 stays unless another N is at least 25% faster than 4
  at a peak under 4 GB, in which case it is proposed to the owner with both numbers.

What could go wrong: the output depending on which thread finished first (the test with contigs of
different lengths guards it); a finished contig holding its segments while it waits (the worker
returns rows); Ctrl-C waiting for the running contigs (said in the spec and the help).
