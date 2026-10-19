---
name: performance
description: How the wall time and the peak memory of pvc_eval are reviewed and improved. Use it when the owner asks for a performance review, a profile, where the time or the memory goes, or to make something faster or smaller; and when a spec or a plan has a target of time or memory. The session measures a baseline and profiles it, sends reviewer subagents in parallel, one per category, to find candidates, then tries each candidate as an experiment, one at a time, and keeps only what the measurement and the code justify.
---

# Performance

The targets are two: the wall time of a run and its peak memory. A
change is judged on both, and on what it does to the code.

This skill comes from the performance review of popnei, which came from
pop_var_caller's, ported on 19 October 2026 to pure Python, with the rules
the owner gave that day for accepting a change. pvc_eval had no code then,
so the categories and the thresholds are to be revised after the first
review that runs on it.

## The principles

- **Measure, do not infer.** Where the time and the memory go is found by
  a profile, and whether a change helps is found by running it. A reviewer
  who reasons that a line is slow is right less often than a profile is:
  in pop_var_caller a first round applied seven findings matched by
  pattern and the benchmark could not see a change; the second round, with
  a profile, found gains in the same code.
- **Right before fast.** No change alters an outcome. Every experiment runs
  the tests first, and its output on the benchmark input is compared with
  the baseline's, byte for byte, before any timing is looked at. A faster
  wrong answer is not a result.
- **Simple code is worth something.** The owner's rule, given on 19
  October 2026: a change that makes the code clearer, simpler or better
  structured is accepted even when its gain is small, or when it costs a
  little time or memory. A change that makes the code more complicated is
  rejected when its gain is small or nothing. The table under "Deciding"
  puts numbers on that.
- **Hot paths only.** A change to code that runs once, at start up or in
  the handling of an error, is noise. What is hot is what the profile says.
- **One change per measurement.** Two changes measured together give a
  result nobody can read.
- **Timings are taken one at a time.** Two runs at once on this machine
  measure each other. Reviewers read, count and profile on small inputs in
  parallel; wall times on the benchmark input are taken by one agent at a
  time.
- **The hard trade-offs are the owner's.** When a change buys time with
  memory, or memory with time, and the table does not settle it, the
  session asks the owner, with both numbers, as the `writing` skill asks
  for a decision.

The format of a finding and its severities are in `finding_format.md`, the
tools that work on this machine in `profiling_environment.md`, Tachyon,
the profiler, in `tachyon.md`, and the categories in `checklists/`. The
session and every reviewer read the first three before they start.

## 0. What is already known

A report of an earlier review, `doc/reports/perf-*.md`, says which
experiments were run, what they gave and what was left open. The review
starts where it ended. The spec may hold a target: `doc/spec.md` §9 says
pvc_eval is meant for benchmarks of megabases and up to a few hundred
samples.

## 1. The workload and the baseline

A review is of a workload, not of a file: a command as a user runs it, on
an input of a size a user has. Choose one to three, write each down with
its command and its input, and keep them in `benchmarks/` so that the next
review runs the same. One input small enough to profile in seconds, and
one at the size the spec names, because a cost that grows faster than the
input does not show on a small one. When there is no input of the right
size, making one is the first deliverable, and the review says so.

Then the baseline, on the commit under review, with nothing else running:

- the wall time, with `hyperfine`, as `profiling_environment.md` gives it;
- the peak memory, the maximum resident set size of `/usr/bin/time -l`;
- the output of the run, kept, to compare every experiment against.

And the noise: the baseline measured twice, apart. The difference between
the two is the smallest change any experiment can claim.

## 2. Where the time and the memory go

Before any reviewer is sent: a profile of the workload with Tachyon, in
`wall` and in `cpu` mode, over all threads, kept as a capture, as
`tachyon.md` says; a profile of the peak memory with tracemalloc; and
cProfile when the question is how many times a function is called. Their
output, as printed, is given to every reviewer, so that each starts from
the same evidence and none runs them again. A review without a profile
gives a long list of candidates matched by pattern, most of which show
nothing when tried.

Tachyon needs root on this Mac, which a session does not have, so the
owner takes each capture from a block of commands the session writes for
them to paste, and the session makes every table and graph from the
capture file itself (`tachyon.md`, first section). The session asks for
all the captures it needs at once, the baseline in `wall` and in `cpu`
mode, and waits for them before sending the reviewers. When the owner
cannot take them, the review goes on with cProfile's counts and
tracemalloc alone, says so at the top of the report, and files no finding
above Likely. It does not turn to pyinstrument or scalene, which turn the
GIL back on and profile another program (`profiling_environment.md`).

When the time grows faster than the input, a second profile on an input
ten times larger shows where.

## 3. The reviewers

| category | sent when |
|---|---|
| `methodology` | always: that the workload, the measurements and the profile can be trusted before any finding about code is acted on |
| `algorithms` | always: complexity, data structures, work done twice |
| `memory` | always: what the run holds, what grows with the input, copies |
| `interpreter` | the profile shows time in a few Python functions called millions of times |
| `io` | the profile shows time in reading or writing, in gzip, in parsing lines |
| `parallelism` | always once pvc_eval uses threads; before, when the profile shows independent work that takes most of the time |

One `perf-reviewer` per category, all in one message. Each gets: its
category; the commit; the workload and its command; the baseline numbers
and the noise; the output of the profiles, as printed; and what the code
does not show. A reviewer that runs code gets `isolation: "worktree"` and
is told to check out the commit first, since such a worktree starts on
`main`, and then `uv sync`.

## 4. Collecting, and the report

Read every report. Promote the findings on which several categories agree
that the same place is hot; demote those the profile shows to be cold. A
reviewer that reported what it could not show gets its category again,
with the gap named.

`doc/reports/perf-<scope>-<date>.md`, written as the `writing` skill says,
for the owner. It is kept, because the next review starts from its
numbers. In this order:

1. The workload, the commit, the machine, and the baseline with its noise.
2. Where the time and the memory go, from the profiles, in a few lines.
3. The candidates in the order they will be tried, each in the form of
   `finding_format.md`, numbered P1, P2, so that they can be referred to.
4. What the code already does well, when there is something specific.

The experiments then run, and each candidate is closed in the report as
they do.

## 5. The experiments

A candidate is a hypothesis, and a measurement decides it. The session
takes the candidates in the order of the report and gives each to a
subagent, `general-purpose`, in the worktree of the review, on a branch
`perf/<scope>`, one at a time. The subagent:

1. makes the one change;
2. runs the checks of `.claude/skills/coding/python.md`, and compares the
   output of the workload with the baseline's; when they differ, it stops
   and reports;
3. measures the wall time and the peak memory as the baseline was
   measured, the count first when the mechanism can be counted, the calls
   of a function in cProfile, the bytes in tracemalloc, since a count is the
   same on every run; and, when Tachyon can be run, a `--diff-flamegraph`
   against the baseline's capture, which shows what the change moved;
4. when the difference in wall time is under twice the noise, or under
   about 5%, measures both sides again, the baseline included, before
   trusting it;
5. reports the numbers as printed, the lines it added and removed, and
   what the change does to how the code reads;
6. commits the change on the branch when the session accepts it, and
   reverts it otherwise.

After several changes are kept, the workload is measured once more against
the first baseline, since gains do not always add up.

## Deciding

What the session does with a measured candidate. Time and memory are
compared with the baseline; "same" is within the noise.

| what the change does to the code | time | memory | verdict |
|---|---|---|---|
| clearer, simpler or better structured | any gain, same, or a small loss | same or better, or a small loss | accept; a loss above about 5% of either is said to the owner |
| about the same | better by more than the noise | same or better | accept |
| about the same | better | worse | accept when the time gained is much larger than the memory lost, in relative terms, and the memory does not grow with the input: half the time for 1% more memory is a clear win. Otherwise the owner's |
| about the same | worse | better | accept when memory is what limits the users and the time lost is small. Otherwise the owner's |
| more complicated | better by 20% or more, or a peak that stops growing with the input | not worse | accept, with a comment at the change that says what it buys, with its numbers |
| more complicated | better by less than 20% | any | reject |
| any | same | same | reject, unless it is simpler |

The 20% and the 5% are where to start, not limits; the report says when
the session decided against them, and why. A change of complexity, a
cost that grew as the square of the input and now grows with it, or a
peak that grew with the genome and now does not, can be accepted on a
smaller gain at today's size, when the larger input shows the gap
widening.

These are never accepted by the session, only proposed to the owner with
their numbers: a dependency, compiled code of any kind, numpy among them,
a change of Python version or build, processes, a cache with no bound, and
a change to any output. Threads are not on this list: the owner chose the
free-threaded 3.15t on 19 October 2026 so that pvc_eval can use them, and
a change that adds them is judged by the table like any other, with the
cost of its code counted in full (`checklists/parallelism.md`).

A candidate is closed in the report as one of: kept, with its commit and
its numbers; no gain, closed, which is an expected outcome and what the
experiment is for; rejected for its complexity, with its numbers; not run,
with the reason; or the owner's.

## What the owner gets

A reply in chat as CLAUDE.md describes them: the workload and the
baseline; where the time and the memory go; what was tried and what each
experiment gave, in seconds and megabytes on the stated input and machine;
what was kept; what is theirs to decide; and the path of the report.
Nothing is merged into `main` without the owner's order.
