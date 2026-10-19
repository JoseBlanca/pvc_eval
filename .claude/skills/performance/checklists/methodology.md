# Methodology

**Purpose.** Before any finding about code is acted on, the workload is
one a user runs, the measurements can be repeated, and the profile shows
the run that matters. A finding here outranks the others, because every
other finding is judged with these numbers.

**Sent.** Always.

## Rules

- **The workload is the user's.** A command as `doc/spec.md` §3 gives it,
  on an input of the size §9 names: megabases, up to a few hundred samples.
  A micro-benchmark of one function stands in for it only after the
  profile of the workload has named the function. A workload with no input
  of the right size is a finding, and making the input is its fix.
- **Two sizes.** One input to profile in seconds and one at the size of the
  spec, or ten times the first. A cost that grows faster than the input
  shows only on the second: a quadratic step looks free on a hundred
  records.
- **The output is part of the measurement.** The baseline's output is
  kept, and every experiment's output is compared with it byte for byte.
  A workload whose output is not kept cannot tell a faster run from a
  wrong one.
- **The noise is measured.** The baseline is run twice, apart; the
  difference is the smallest change that can be claimed. hyperfine's
  spread inside one call is not that number: a busy machine moves the
  mean between calls by more than the spread within one.
- **Same conditions.** The same machine, the same Python (`uv run python
  -VV`), the same input, the files warm in the page cache on both sides,
  nothing else running. A baseline and a candidate are timed in one
  hyperfine call and again with the order reversed, since hyperfine runs
  one command's runs and then the other's (`profiling_environment.md`).
- **Sampling before tracing.** Tachyon for where the time goes; cProfile
  for how many calls, never for how long a small function takes, since its
  overhead falls on every call. A finding whose time comes only from
  cProfile is at most Likely.
- **Enough samples.** A Tachyon profile with fewer than about 10000
  samples gives shares that move by points from run to run
  (`tachyon.md`). A profile under that is a finding of this category.
- **The profiler measures the program that runs.** A tool that turns the
  GIL back on, pyinstrument, scalene, line_profiler, profiles a program
  whose threads take turns, which pvc_eval on 3.15t is not. Its output is
  not evidence, and a finding that rests on it is filed again on Tachyon's.
- **Wall and CPU, apart.** A profile in `wall` mode and one in `cpu` mode
  over all threads. Where they differ, the time is spent waiting, and the
  question is what is waited on, not how fast the code is.
- **Peak, not total.** The memory target is the peak, the maximum resident
  set size of `/usr/bin/time -l`. A total of allocations says how much was
  allocated over the run, which matters to time and not to the peak.
- **Start-up is its own question.** For a run under a second, the import
  of pvc_eval can be most of it: `uv run python -X importtime -c "import
  pvc_eval" 2> imports.txt`. For the workloads of the spec, it is not.
- **One change per measurement.** A finding that bundles two changes is
  split.

## Output

Prepend `[workload]`, `[measure]` or `[profile]` to the title of a
finding of this category, so that the session puts it first in the report.
