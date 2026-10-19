# Tachyon, the sampling profiler of Python 3.15

Tachyon is `profiling.sampling` of the standard library, new in Python
3.15. It reads the stacks of a running Python process from outside, a
thousand times a second by default, so the program runs at its own speed
and every thread, and the state of each, can be seen. It is the profiler
of choice for pvc_eval, decided by the owner on 19 October 2026, for two
reasons: it does not load anything into the process, so it profiles the
free-threaded 3.15t that pvc_eval runs on as it is (pyinstrument, scalene
and line_profiler load C modules that turn the GIL back on, which profiles
another program); and it sees threads, which is what 3.15t was chosen for.

This file is read by the session and by every `perf-reviewer` before a
profile is taken or a measurement plan names Tachyon. It is written from
the documentation of Python 3.15
(https://docs.python.org/3/library/profiling.sampling.html), read on 19
October 2026, and from what was run that day, which the last section
lists. What was not run is marked so.

## Permission: the owner takes the captures

Tachyon reads another process's memory, which the system guards. On
macOS that takes root, or a profiler signed with the
`com.apple.security.cs.debugger` entitlement; without either it stops with
"Permission Error: Unable to access process memory on macOS", as it did
here on 19 October 2026. On Linux it takes root, `CAP_SYS_PTRACE`, or
`kernel.yama.ptrace_scope` set to 0.

Neither the session nor its subagents have root. So the owner decided on 19
October 2026 to take the captures themselves, from blocks of commands the
session writes for them to copy and paste. Only the capture needs the
permission. Everything after it, the tables, the flame graphs, the
heatmaps and the diff against a baseline, is made by the session with
`replay` from the capture file, which needs no permission. So the owner
pastes one block per capture and nothing else.

### The block the session writes for the owner

The owner will paste it as it is, so the session fills in every part, with
absolute paths, and checks it before handing it over:

    cd /Users/jose/devel/pvc_eval/.claude/worktrees/perf-loci
    uv sync --locked
    mkdir -p .profiles
    sudo -v
    .venv/bin/python -m pvc_eval compare --reference /abs/ref.fa --a /abs/a.vcf.gz --b /abs/b.vcf.gz --output .profiles/base-wall.tsv &
    PID=$!
    sleep 1
    sudo .venv/bin/python -m profiling.sampling attach -a -r 10khz --binary -o .profiles/base-wall.bin $PID
    wait $PID
    sudo chown "$(id -un)" .profiles/base-wall.bin

Why each part is so, since each was a trap:

- **`cd` to the worktree** of the review, absolute, as the first line: the
  owner's shell can be anywhere.
- **`uv sync --locked`**, so the environment is the commit's.
- **`sudo -v` before the program starts**, so that the password is asked
  while nothing runs, and the attach then starts without waiting for it.
- **The program is started by `.venv/bin/python`, not `uv run`.** `uv run`
  stays as the parent of Python, so `$!` would be uv's process and Tachyon
  would attach to uv; checked on 19 October 2026 (uv 0.13.1: `$!` was
  76577, Python's own id 76578). This is the one place where `python.md`'s
  rule of `uv run` for everything gives way.
- **Only the profiler runs as root.** The program runs as the owner, so
  everything it writes, and the `.venv`, stay theirs. `sudo uv run` is
  never used: it would leave files owned by root in `.venv` and uv's cache.
  The profiler is the same `.venv/bin/python`, a link to the project's
  3.15t, which is what the rule of matching versions below asks.
- **`sleep 1`** gives the interpreter time to start before the attach; it
  loses the first second of the run, which on a workload of tens of seconds
  does not matter.
- **`chown`**, because the capture was written by root.
- **No comments in the block**: zsh, the owner's shell, does not take `#`
  comments typed at its prompt unless `interactivecomments` is set. What
  the block does is said in the message, before it.
- **One capture per block**, its name saying what it is: `base-wall`,
  `base-cpu`, `p3-wall` for candidate P3. The output of the run goes beside
  it, for the comparison with the baseline's.

`.profiles/` is in `.gitignore`; captures are not committed. The message
around the block says what it measures, about how long it runs, and what
the owner sends back: only that it ended, or the last lines it printed
when it did not.

The block was run on 19 October 2026 in zsh on a scratch package, with
`sudo` replaced by a function that runs its command without root: the
attach found the program and stopped only at the permission error, which
`sudo` lifts. A program that ended in under a second was gone before the
attach ("Process with PID ... does not exist"), so the workload of a block
runs for at least a few seconds. The block with the real `sudo` has not
been run; the session of the first review that has a capture taken with it
corrects this section with what happened.

## Versions must match exactly

The profiler and the program must be the same Python: the same minor
version, the same release when it is a pre-release, and the same build,
since a free-threaded profiler cannot read a standard program or the
reverse. So the profiler is the project's own interpreter,
`.venv/bin/python` in the owner's block and `uv run python` for a
`replay`. The `python3.15` of another install does not do.

## The four commands

| command | what it does | for |
|---|---|---|
| `run [-m] target args` | starts the program and samples it to the end | a workload from its start |
| `attach PID` | samples a running process for `-d` seconds, or to its end | the large workload, started apart; the only form when only the profiler has root |
| `dump PID` | prints one snapshot of every stack and exits; each thread marked on CPU, idle, waiting, has exception | a run that hangs, or threads that wait on each other |
| `replay file.bin` | turns a binary profile into any other format, without running again | one capture, several views |

## Modes: what a sample counts

`--mode` decides which samples are kept.

- **`wall`**, the default: every sample, waiting included. The wall time of
  a run is what the owner's target is, so the first profile is in this
  mode.
- **`cpu`**: only samples where the thread is on a core. A function high in
  `wall` and low in `cpu` is waiting, on a file, a lock or another thread,
  and making its code faster gains nothing. The pair of profiles is the
  first thing to compare.
- **`gil`**: samples where the thread holds the GIL. 3.15t has no GIL, so
  this mode says nothing for pvc_eval; it is for a comparison against the
  standard build, if the owner asks for one.
- **`exception`**: samples while an exception is in flight or an `except`
  block runs. pvc_eval raises only on a fault of the input, so time here on
  a clean workload means exceptions are used as control flow somewhere hot.

## The options that matter here

- **`-a`, `--all-threads`.** Only the main thread is sampled without it.
  Every profile of pvc_eval once it uses threads takes `-a`; a profile
  without it misses the workers and shows the main thread waiting.
- **`-r`, the rate**, 1 kHz by default; `10khz` or `20khz` for a run of a
  few seconds. The documentation gives the error of a share as about ±0.5
  points for a function at 5% with 100000 samples, and anywhere from 3% to
  7% with 1000. So a profile has at least 10000 samples before a share is
  believed: ten seconds at 1 kHz, one second at 10 kHz; under that, raise
  the rate or profile a larger input. `--realtime-stats` shows the rate
  actually reached, which can fall below the one asked for.
- **`-d`, the duration**, for `attach`, or to stop a long run: 10 to 30
  seconds of a steady workload is usually enough.
- **`--native`** adds `<native>` frames for time in C, the standard
  library's gzip or `str` methods among it, which is otherwise given to the
  Python function that called it.
- **`--no-gc`** hides `<GC>` frames. They are left on: much time in `<GC>`
  means many objects made and dropped, a finding for
  `checklists/memory.md`.
- **`--opcodes`** shows the bytecode being run, and whether the
  interpreter specialized it (`LOAD_ATTR_INSTANCE_VALUE` against
  `LOAD_ATTR`). For `checklists/interpreter.md`, on one hot function, never
  as a first profile; it costs memory.
- **`--blocking`** stops the program at each sample for a consistent stack.
  Only when the stacks of a normal profile do not make sense, which the
  documentation says happens most with generators, as pvc_eval's stages
  are; then at 1 ms or slower, since it slows the program.
- **`--subprocesses`** follows child processes; pvc_eval has none unless
  the owner allows processes.

## Outputs

| option | what it gives | when |
|---|---|---|
| `--pstats`, the default | a table to stdout: samples, share, self and total time per function, and a summary of hot spots; `--sort`, `-l` (15 rows by default), `--no-summary` | the text a reviewer reads and quotes |
| `--binary -o prof.bin` | a compact capture, for `replay` and for `--diff-flamegraph` | every profile that is kept: one capture, any view later |
| `--flamegraph -o x.html` | an interactive HTML flame graph, with a filter by thread | the owner's view of where the time goes |
| `--diff-flamegraph base.bin -o diff.html` | the flame graph of this run coloured against a baseline capture: red worse, blue better, purple new | an experiment of the `performance` skill: what the change moved |
| `--heatmap -o dir` | one HTML page per source file, each line coloured by its time | the lines inside a hot function, in place of line_profiler, which turns the GIL on |
| `--collapsed` | one line per stack with its count | a diff by script, or another tool |
| `--gecko -o x.json` | a file for the Firefox Profiler, with markers for CPU and GC over time | the timeline of the threads; `--mode` cannot be used with it |
| `--live` | a view like `top` in the terminal | a person watching; not for a subagent |

A text output is what goes into a report or a reviewer's prompt; an HTML
file is given to the owner by its path.

## The recipes

What the owner captures, each with the block above, changing only the
attach line and the names:

| capture | the attach line takes |
|---|---|
| the baseline, where the wall time goes | `-a -r 10khz --binary -o .profiles/base-wall.bin` |
| the baseline, waiting told from working | `-a --mode cpu -r 10khz --binary -o .profiles/base-cpu.bin` |
| a candidate of an experiment | the same as the baseline it is compared with, `-o .profiles/p3-wall.bin`, in the candidate's worktree |
| bytecode of a hot function | `-a -r 10khz --opcodes --binary -o .profiles/base-opcodes.bin` |

What the session then makes from them, with no permission, in the
worktree:

    uv run python -m profiling.sampling replay --sort cumul-pct -l 30 .profiles/base-wall.bin
    uv run python -m profiling.sampling replay --sort tottime -l 30 .profiles/base-cpu.bin
    uv run python -m profiling.sampling replay --flamegraph -o .profiles/base-wall.html .profiles/base-wall.bin
    uv run python -m profiling.sampling replay --heatmap -o .profiles/base-heat .profiles/base-wall.bin
    uv run python -m profiling.sampling replay --diff-flamegraph /abs/base-wall.bin -o .profiles/p3-diff.html .profiles/p3-wall.bin

`cumul-pct` sorts by the share of samples in a function and everything it
calls, so that a stage of the pipeline and its helpers read as one;
`tottime` by the time in the function itself, which is where a change of
code acts. `replay` takes `--diff-flamegraph` with the baseline's capture,
as its `--help` showed on 19 October 2026, so an experiment needs one
capture of the candidate and no second run. The tables go into the report
and to the reviewers as printed; the HTML files go to the owner by path.

A run that hangs, or threads that wait on each other, is a different
block: the owner runs `sudo .venv/bin/python -m profiling.sampling dump -a
<PID>` from the worktree, with the process id the session gives, and sends
back what it printed.


## What Tachyon does not answer

- **How many times a function was called.** Sampling counts time, not
  calls. The count is cProfile's, `profiling.tracing` in 3.15, which is the
  same class under its new name.
- **A difference of 1% or 2%.** Under the noise of sampling. hyperfine on
  the workload answers it, or nothing does.
- **Memory.** Tachyon samples stacks, not allocations. The peak is
  `/usr/bin/time -l`, and where it comes from is tracemalloc's, as
  `profiling_environment.md` says.
- **A run under a second**, which gives too few samples: profile a larger
  input.
- **A function shorter than the interval**, which can fall between samples
  every time.

## What was run on 19 October 2026

- `python -m profiling.sampling run --help` on 3.15.0rc2, and `replay
  --help` on 3.15.0: the options above are there, `--diff-flamegraph`,
  `--sort` and `-l` among those of `replay`.
- `run` on a script, in the default and in `cpu` mode, on this Mac without
  root: the permission error of the first section.
- `run` inside a Linux container of Apple's `container` (image
  `python:3.15-rc-slim`, 3.15.0rc3, as root, with `--cap-add
  CAP_SYS_PTRACE`): it printed nothing and had not ended after ten
  minutes, on a script that runs in 0.2 s, and was stopped. A container is
  not a way around the permission until someone finds why.
- Signing a copy of the interpreter with the debugger entitlement was not
  tried: the session was refused it as a weakening of a protection, and
  the owner chose to take the captures instead.

So every recipe above comes from the documentation and is to be checked
on its first use, and this section brought up to date.
