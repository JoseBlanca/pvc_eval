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

How: a worker that segments one contig and returns the formatted rows of its STR regions, those the
requested spans overlap; the contigs mapped over it in the order of the `.fai`, in this thread when
`num_threads` is 1, and otherwise on a `ThreadPoolExecutor` whose `map` has `buffersize` so that at
most about N results wait (as `vcf.py`'s `_map_in_order` does; `vcf.py` is not changed). The rows are
written as they come. The checks of `num_threads` follow `parse_vcf`'s: a whole number from 1 to 256,
`ValueError` from the function, and an argparse error naming `--threads` from the command line.

Deliverables, each checked:
- the tests of the spec's "How it is verified" for `--threads`, which fail before the change;
- the BEDs of `benchmarks/data/chr1_10mb_4mb.fa` and `chr1_10mb_16mb.fa` (in the
  `perf-str-finder` worktree's `benchmarks/data/`) unchanged at 1 and 4 threads;
- a measurement, reported, not committed: wall time and peak RSS of `str-regions` over GRCh38's
  chromosomes 13 to 22 and Y (about 1 Gb) at 1, 2, 4 and 8 threads, the BED compared with the
  1-thread one each time, from which the default of 4 is confirmed or a different one proposed to the
  owner.

What could go wrong: the output depending on which thread finished first (the test with contigs of
different lengths guards it); a finished contig holding its segments while it waits (the worker
returns rows); Ctrl-C waiting for the running contigs (said in the spec and the help).
