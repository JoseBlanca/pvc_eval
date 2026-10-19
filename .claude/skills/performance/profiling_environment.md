# The profiling environment: what runs on this machine

This is inventory, not method: which tools were run on this machine with
the project's Python, which work, the commands that worked, and what does
not. Every entry was checked on 19 October 2026 on a scratch package made
with `uv init --package`; entries drift, so a tool far from that date is
checked again before a plan names it.

## The machine and the Python

An Apple M5 Pro, 18 cores, 64 GB, macOS 27.0.1. uv 0.13.1. The project's
Python is CPython 3.15.0, free-threaded, pinned as `3.15+freethreaded`
(`uv python pin 3.15t` writes it so). The owner chose it on 19 October
2026 to run threads in parallel. A timing is valid for the Python it was
taken on, and the report says which.

## The first trap: tools that turn the GIL back on

On 3.15t a C extension that has not declared that it runs without the GIL
turns the GIL back on for the whole process when it is imported, and
Python prints a `RuntimeWarning` that says so. A profiler that does this
measures another program, one whose threads take turns. Every tool below
was checked with `sys._is_gil_enabled()` printed from inside the run:

| tool | GIL during the run | use on 3.15t |
|---|---|---|
| Tachyon, `profiling.sampling` | not loaded into the process | yes, with the permission `tachyon.md` describes |
| cProfile, `profiling.tracing` | off | yes |
| tracemalloc | off | yes |
| pyperf | off | yes |
| pytest, hypothesis | off | yes |
| pyinstrument 5.1.3 | **turned on** | no |
| scalene 2.3.0 | **turned on** | no |
| line_profiler | **turned on** | no |
| memray 1.20.0 | off, but measures nothing | no: see below |

A tool not in the table is checked the same way before it is used.

## How a tool is run: `uv run --with`, not `uvx`

`uvx <tool>` runs the tool in an environment of its own, which does not
have pvc_eval and does not read the project's pin. A tool is run inside the
project's environment, with `uv run --with <tool>`, which adds it for one
command and leaves `pyproject.toml` and `uv.lock` as they were, or is added
with `uv add --dev` when every review uses it.

## Present and checked

- **hyperfine 2.0.0** (`/opt/homebrew/bin/hyperfine`, checked 19 October
  2026; it was 1.20.0 when this file was first written): the wall time of a
  command, with its spread.

      hyperfine -N --warmup 2 --runs 10 --export-json base.json \
        'uv run --no-sync pvc_eval compare ...'

  `--no-sync` keeps uv from syncing in the middle of a benchmark. uv's own
  start costs about 5 ms with it and 10 ms without, measured on 3.14 on 19
  October 2026: nothing against a run of seconds, and the same on both
  sides. Two commands in one call are compared, but not interleaved:
  hyperfine runs all the runs of the first, then all of the second, and has
  no option to alternate them (checked on 2.0.0, 19 October 2026). So a
  baseline and a candidate in two worktrees are timed in one call and then
  in a second call with the order reversed; a drift of the machine during
  a call then favours each side once, and a gain is believed only when both
  calls show it.

- **`/usr/bin/time -l`** (macOS): the peak memory. The line to read is
  "maximum resident set size", in bytes. Through `uv run` it is the
  Python process's, and the line "peak memory footprint" is uv's own: a
  script that held 300 MB gave 317308928 and 8847888 on those two lines
  through `uv run --no-sync`, and 317308928 and 307970768 run directly.
  Take the median of three runs; differences of a few MB are the
  allocator's noise.

- **Tachyon**, the sampling profiler: where the time goes, thread by
  thread. Its own file, `tachyon.md`, says how to use it. On this Mac it
  needs root, which a session does not have: the owner takes the captures
  from a block the session writes, and the session reads them with
  `replay`, which needs no permission.

- **cProfile**, `profiling.tracing` under its 3.15 name (the same class):
  exact call counts, which say when a function is called a million times.
  Its overhead falls on every call, so the times it shows for small
  functions are inflated; the time is Tachyon's.

      uv run python -m cProfile -o prof.out -m pvc_eval compare ...
      uv run python -c "import pstats; pstats.Stats('prof.out').sort_stats('ncalls').print_stats(25)"

  A count is the same from run to run once the modules are compiled: five
  runs of one script gave 698 calls each, and the run just after the
  `__pycache__` was removed gave 763, the compiling counted with it. So run
  once before counting, and gate on the line of the function the finding
  names, not on the total. One run that day, the first on that script, gave
  600697 in total for the same 698; why was not found, and is one more
  reason to read the function's own line.

- **tracemalloc** (standard library): where the peak memory comes from, on
  3.15t, now that memray cannot. On a script that held a 300 MB bytearray
  and 200000 short strings it gave a peak of 314.3 MB and named both lines,
  286 MiB and 13.5 MiB. It slows the run, so its time is not looked at.

      uv run python -X tracemalloc=10 scripts/peak.py

  where the script calls the stage under study and prints
  `tracemalloc.get_traced_memory()` and the top lines of
  `tracemalloc.take_snapshot().statistics("lineno")`. Such a script is kept
  in `benchmarks/` with the workload.

- **pyperf**: a micro-benchmark of one function, with its own warm-up:
  `uv run --with pyperf python -m pyperf timeit -s 'setup' 'call'`. For the
  rare finding about one function; the workload is what decides.

- **ruff 0.17.1 and ty 0.0.89** run on 3.15t, and both accept the
  `sentinel` builtin in type hints; ty narrows `x is MISSING` as it does
  `x is None`.

## Not usable

- **memray 1.20.0**: on 3.15.0 and on 3.15t it reported a peak of 0 B for
  the 300 MB script, with and without `--native`, while on 3.14 it reported
  300.164 MB. Check a newer release against that script before using it.
- **pyinstrument, scalene, line_profiler**: they turn the GIL on, as above.
- **py-spy**: needs root on macOS, and Tachyon does its work.
- **`uvx`** for anything that imports pvc_eval, as above.

## Datasets

None yet. The first review makes them, as step 1 of `SKILL.md` says, and
lists them here with where they came from and their sizes.
