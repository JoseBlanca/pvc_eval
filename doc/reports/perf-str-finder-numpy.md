# Performance review of the STR finder with numpy, 21 October 2026

This report is about step 1 of the STR finder, the search for candidate repeats, rewritten with
numpy. The STR finder is `segment_contig` in `src/pvc_eval/segments.py`, and the command
`pvc_eval str-regions` runs it over every contig of a reference (`doc/specs/str_finder.md`). The
owner decided on 21 October 2026 that numpy may be a runtime dependency, and this review measures
what it gains. It is the second review of the finder. The first, `doc/reports/perf-str-finder.md`
of 20 October 2026, recommended against numpy; this review starts where it ended and tests its
reading. The work is on the branch `perf/str-finder-numpy`, from main at 96b3938, on an Apple M5
Pro (18 cores, 64 GB) with CPython 3.15.0 free-threaded (3.15t, GIL off) and numpy 2.5.4.

**The result.** With step 1 on numpy arrays, `str-regions` takes about 0.3 of main's time, and
gives the same BED and the same segments everywhere it was run:

| input | main | numpy | peak RSS, main → numpy |
|---|---|---|---|
| chr1, 1 thread | 79.2 s | 24.9 s | 1,157.6 → 1,176.8 MB |
| GRCh38, 1 thread | 1,169.6 s | 376.1 s | 1,855.9 → 1,977.3 MB |
| GRCh38, 4 threads | 389.4 s | 172.6 s | 5,213.4 → 5,093.4 MB |
| GRCh38, 8 threads | 391.5 s | 98.4 s | 8,326.0 → 8,213.6 MB |

These were measured between 15:26 and 16:36 while other sessions' jobs held the load average
between 3 and 14, so the times are slower than on a quiet machine on both sides (main took
1,020.8 s over GRCh38 at 12:44); chr1 is the median of 2 pairs run alternately, GRCh38 one run
per side, back to back. The first review's reading was half right. The loop at the heart of step
1, Ruzzo and Tompa's algorithm, still runs in Python, but array code now settles 97.6% of the runs
of matches without it and does all the work on the candidates. What numpy costs: a dependency, about 10 MB and 24 ms at start, about 20
MB more at the peak of one contig, and 175 more lines in `segments.py`, with main's pure-Python
step 1 kept in the tests as the oracle.

**What is the owner's**, with the numbers at the end of this report:

- the merge of `perf/str-finder-numpy` into main, which this report recommends;
- the peak memory at 8 threads, which was 11% to 14% higher than main's on chromosomes 13 to 22
  and Y in two runs, and 1% lower over GRCh38 in one; it does not need to hold up the merge;
- pvc's single sweep for absorption (step 7, which grows the satellites over the STR segments and
  bundle members near them), now 47% of a run over GRCh38 on one thread;
- a compiled extension for step 1, which is not recommended now.

The limit on the match reward and the mismatch penalty, which the owner approved during the
review, is in the spec and the code (below).

## How step 1 works now

Step 1 scores, for each period p from 1 to 6, every base against the base p before it: the match
reward (2) when both are the same A, C, G or T, minus the mismatch penalty (7) otherwise. Its
candidates are the maximal scoring subsequences of those scores, found by Ruzzo and Tompa's
algorithm, which keeps a stack of the subsequences still open. Each must hold at least two copies
of its motif (the copies test), and its motif must not be a shorter motif repeated (the motif
test). A candidate is then dropped when a kept
candidate of a shorter period that divides its own covers at least half of it (the drop rule),
and step 3 drops those under the copy floor of their period. Main does this in pure Python, a
block of 4,096 bases at a time, one step per run of matches.

The new code reads the contig in blocks of 65,536 bases and carries each period's state from one
block to the next: how far it has scored, its running total, the least running total so far,
and the stack. Per block and period:

- The runs of matches come from one comparison of two arrays, and the running total before each
  run from the runs' starts and lengths. The stack is empty at a run's start exactly when its
  total there is the least so far.
- A run that starts on an empty stack, and whose next run also does, is a maximal scoring
  subsequence on its own: 97.6% of the runs on 4 Mb of chr1, all found on arrays. The other
  runs, 110,889 on 4 Mb where main's loop takes 1,186,811 steps, go through Ruzzo and Tompa's
  rules in Python, as main and pvc do.
- The copies test, the motif test, the drop rule and the copy floor run on arrays, with main's
  bookkeeping of which candidates can be tested yet and which kept ones a longer period may still
  need. `_Candidate` objects are made only for the candidates found, about 30,000 on 4 Mb
  instead of a million.

The code departs from pvc's `tandem_repeat.rs` in how the work is done, not in what it finds, and
the docstring of `_find_candidates` lists each departure with the reason the candidates are the
same. Main's pure-Python step 1 is kept unchanged in `tests/pure_python_step_one.py` as the
oracle, as pvc keeps its plain fixed point of absorption as a test helper. A hypothesis test
compares the two over generated sequences, parameters and block sizes, and pvc's crafted cases of
Ruzzo and Tompa run against both. Outside the suite, 84,000 more random cases on the final code
gave no difference from the oracle.

## The workload and the baseline

The workload is the first review's (`benchmarks/README.md`):

    uv run --no-sync pvc_eval str-regions --reference <stretch>.fa --output out.bed

on stretches of GRCh38 chr1 from position 10,000,000, of 4, 16 and 64 Mb, at the default
parameters; then all of chr1, chromosomes 13 to 22 and Y as one FASTA (855 Mb, the input of
`doc/reports/str-regions-threads.md`), and all of GRCh38. Every experiment had to give main's BED,
byte for byte, and main's sha256 of every segment, generic regions included
(`benchmarks/digest_str_finder.py`), on the four stretches; what was kept was also compared on
chr1, on the 855 Mb input at 1, 4 and 8 threads, and on GRCh38 at 1, 4 and 8 threads, and passed
the whole-port test (`tests/test_whole_port.py`), which compares pvc_eval's STR regions with those
of a pvc built at bbcc9f4 on 600 generated contigs and three slices of GRCh38.

The baseline is main at 96b3938, measured on 21 October 2026 between 12:09 and 12:14, when the
load average was 1.2 to 1.9:

| input | wall time | peak RSS |
|---|---|---|
| 4 Mb | 1.212 s | 59.0 MB |
| 16 Mb | 4.941 s | 113.1 MB |
| 64 Mb | 20.246 s | 339.7 MB |
| chr1 | 76.8 s | 1,157.6 MB |
| GRCh38 | 1,020.8 s | 1,855.2 MB |

The times are means of alternating pairs (below); chr1 is 2 runs and GRCh38 1. Python with
`pvc_eval.segments` imported holds 28.5 MB; with numpy, 38.9 MB.

**The noise, and how each candidate was judged.** Other sessions ran pvc's benchmarks on this
machine through the day, and the load average went from 1.2 to 34. Main against itself gave a
median ratio of 1.11 over 2 pairs on 64 Mb at load 8 to 34, and 1.000 over 10 pairs on 4 Mb at
load 1.4 to 1.7. So every candidate was timed against its base in alternating pairs from two
worktrees (`benchmarks/paired.sh`), and judged by the median of the ratios of its pairs, with
main against itself in the same window. `paired.sh` was changed for this review, on the
methodology reviewer's finding that it always ran A first: it now runs A first in odd pairs and
B first in even ones, records the CPU time and the load before each run, and compares the outputs
of every pair. A gain is claimed beyond twice the same-code deviation of its window, and at least
3%.

Threads are measured with a new script, `benchmarks/threads_str_finder.py`: N threads each run
`segment_contig` on the 4 Mb stretch at once, and it prints the time of each call and the CPU
time of all of them against N times one thread's. CPU time that grows with N is contention
between the threads, and the load of other jobs hardly moves it; the wall time of an 855 Mb run
with `--threads` measures, under load, mostly which jobs of others overlap it.

**No Tachyon capture was taken.** Tachyon needs root, which the owner gives by pasting a block of
commands, and the session ran without the owner at the terminal. The evidence is cProfile's
counts, the step timer (`benchmarks/stages_str_finder.py`), timers around single functions, and
tracemalloc. The first review's Tachyon profiles of nearly the same code found nothing waiting.

## Where the time and the memory went at the start

On 16 Mb, step 1 (with step 3, which it holds) took 5.40 s of 6.25 s and absorption (step 7)
0.49 s, at load 8 to 10. Over GRCh38 the first review gave step 1 about 80% of the time and
absorption about 220 s. cProfile on 4 Mb split step 1 into the Ruzzo–Tompa loop, 1.04 s of 3.54
s under the profiler, and the work on the candidates, about 2.3 s: 999,207 calls building
candidates, 1,092,462 motif tests and 998,562 passes of the drop rule. cProfile inflates the
second part, whose functions are called a million times, so these are counts more than shares.

The peak memory of one contig came at step 7, at about 4.6 to 5.3 MB per Mb of contig, as the
first review found; step 1 held about 0.7 to 1.0 MB per Mb less than the peak.

## The candidates

Five reviewers, one for each of the methodology, algorithms, memory, interpreter and parallelism
categories, read the code, the profiles and the first review, with numpy in scope. They agreed on
one place, step 1, and three of them prototyped it. Their findings that matter here:

- **The algorithms reviewer**: array code can decide most of step 1 exactly. On 4 Mb, of 4.6
  million runs of matches over the six periods, 4.49 million start on an empty stack and are
  followed by a run that also does, so each is a subsequence on its own; 110,889 runs need the
  stack. Its prototype gave `_find_candidates`'s candidates on 1 Mb, 4 Mb and 23,000 random
  sequences, and took step 1 on 4 Mb from 2.34 s to 0.53 s.
- **The interpreter reviewer**: numpy pays 0.16 to 0.66 µs per call, so it wins on operations over
  a block of thousands of bases and loses on any call per run or per candidate. Today's per-base
  work in bytes (about 42 ms per 4 Mb) was never the cost; the steps in Python were.
- **The memory reviewer**: an array version holds about 60 bytes per base of the stretch it works
  on, for six periods, so it must work in blocks: 3.7 MiB at 65,536 bases, 14 GB if it held all of
  chr1.
- **The parallelism reviewer**: numpy's functions called through the module contend across
  threads on 3.15t, and its array methods do not. 100,000 calls of `np.cumsum` on 64 bytes took
  0.18 to 0.20 s on one thread and 0.78 to 1.81 s on each of eight threads; the same through
  `np.minimum.accumulate`, 0.04 s and 0.17 to 0.23 s; `a.nonzero()` and `a.searchsorted()` did
  not contend. So array code for threads makes few numpy calls per block, through methods.
- **The methodology reviewer**: the fixes to `paired.sh` and the thread script above.

The candidates, in the order they were tried:

- **P1. Step 1 on numpy arrays.** *Kept, in its second form (commits 8831a58, e6aa98a).*

  The first form (8831a58) cut each contig where every period's stack was empty at once, so that
  each piece could be worked whole, and doubled the piece while no such place was found. Its
  outputs equalled main's everywhere, and against main in alternating pairs at load 3 to 7 it gave
  median ratios of 0.365 on 4 Mb (20 pairs), 0.338 on 16 Mb (10) and 0.336 on 64 Mb (5), and chr1
  in 29.0 s against 84.9 s. Then the correctness reviewer of the code review that followed (three
  reviewers, of correctness, tests and design) found that a match reward high against
  the penalty, which the spec accepts, leaves some stack open everywhere, so the piece grows to
  the whole contig:

  | `--scan-match-reward` / `--scan-mismatch-penalty` | main, 4 Mb | first form, 4 Mb |
  |---|---|---|
  | 2 / 7, the defaults | 59 MB | 85 MB |
  | 2 / 3 | not measured | 126 MB |
  | 1 / 1 | 68.2 MB | 178.8 MB |
  | 7 / 1 | 50.7 MB | 717.7 MB |

  At 7 / 1 that is about 170 bytes per base, tens of GB on chr1, where main holds 51 MB. Memory
  that grows with the contig is a defect, so the second form (e6aa98a) reads fixed blocks and
  carries each period's state across them, as described above. It gave 67.9 MB at 7 / 1 and
  89.2 MB at 1 / 1 on 4 Mb, and against main, at load 1.2 to 1.9: median ratios 0.367 on 4 Mb (20
  pairs), 0.345 on 16 Mb (10) and 0.340 on 64 Mb (5); chr1 0.373 (2 pairs, 76.8 → 28.6 s); peak RSS
  20 to 21 MB above main's on each stretch, not growing with it.

  It is more complicated than main's code, and its gain is far over the 20% that the performance skill asks of a change that makes
  the code more complicated: kept, with
  its numbers in the docstring of `_find_candidates`. The two code reviews each found real
  defects, now fixed: step 3 had moved into the oracle, so pvc's ported tests of it checked test
  code (now `_drop_under_copy_floor` in `segments.py` again); six likely mistakes in the array code, three found by each
  review, passed every fixed test and were caught only by the property test, and only under some
  of its random seeds (each now has a fixed case that fails on it); and the
  running totals overflowed int64 at a mismatch penalty of 2^40 on 16 Mb (they now restart at 0
  wherever a block ends on an empty stack, and the parameters have a limit, below).

- **P2. numpy methods in place of numpy's module functions.** *Folded into P1.* After the
  parallelism reviewer's measurement, P1 was written with methods and a few module calls per
  block, each marked in the code, and the rule is stated with its numbers beside `_BLOCK_BASES`.
  It was not measured on its own.

- **P3. The block size.** *No gain, closed.* Against 65,536 bases, on 16 Mb in 8 pairs: blocks of
  16,384 took 8% longer (median ratio 1.082) for 4.5 MB less peak RSS, and blocks of 262,144 2%
  less (0.980) for 11.6 MB more per thread. The comment on `_BLOCK_BASES` gives these numbers.

- **P4. The running totals only at the runs' boundaries.** *Kept, commit 9cc05d4.* The total
  before a run is worked from the runs' starts and lengths, and the least total so far is the
  least over the runs' starts, since mismatches only lower it; three passes over every score go.
  Against e6aa98a at load 1.7 to 2.4: median ratios 0.844 on 4 Mb (20 pairs), 0.837 on 16 Mb (10)
  and 0.839 on 64 Mb (3), and 3 to 4 MB less peak. The code reads about as plainly, so a gain over
  the noise keeps it.

- **M1. Let go of the candidates before absorption.** *Not run.* The memory reviewer found that
  step 7 keeps the list of candidates only to make the repeat runs, and that freeing it lowers
  tracemalloc's peak on 4 Mb from 19.1 to 17.8 MiB. It is outside step 1 and small beside the
  noise of RSS (about 4 MB on 4 Mb); it can be tried with the sweep for absorption.

- **M2. Candidates as arrays through steps 4 to 7.** *Not run.* It could save up to 4 MiB of 17
  on 4 Mb, but only by rewriting steps 4 to 7 on indices, away from pvc's structure, and the peak
  is not what limits a run.

**The limit on the reward and the penalty.** The spec's table allowed any match reward or mismatch
penalty of 1 or more, and said pvc leaves them unbounded; pvc holds them as i32. The owner
approved pvc's limit, 2,147,483,647, on 21 October 2026. It is in `doc/specs/str_finder.md`
(419c85b), in `StrFinderParameters`, which refuses a larger value naming the field, and in the
`--help` of `--scan-match-reward` and `--scan-mismatch-penalty` (83c3dff). At that limit the
64-bit totals cannot overflow on a contig shorter than 2^31 bases, and a slow test runs 1 Mb at
the largest values against the oracle.

## The result on chr1 and GRCh38

Measured on 21 October 2026 between 15:26 and 16:36, at the final commit, 419c85b, against main;
other sessions held the load average between 3 and 14.5. Main against itself in this window: a
median ratio of 1.010 over 10 pairs on 4 Mb.

| input | main | numpy | median ratio | peak RSS, main → numpy |
|---|---|---|---|---|
| 4 Mb, 20 pairs | 2.395 s | 0.744 s | 0.304 | 59.1 → 76.0 MB |
| 16 Mb, 10 pairs | 9.380 s | 2.636 s | 0.278 | 113.2 → 130.1 MB |
| 64 Mb, 5 pairs | 23.352 s | 7.408 s | 0.291 | 339.7 → 357.4 MB |
| chr1, 2 pairs | 79.200 s | 24.945 s | 0.315 | 1,157.6 → 1,176.8 MB |

The means are slower than on a quiet machine on both sides (main took 1.212 s on 4 Mb at 12:09);
the ratios are what was judged. Every pair gave equal outputs.

**Threads.** One `segment_contig` of 4 Mb per thread, three rounds per side, load 2.8 to 4.3. The
CPU time of all the calls against N times one thread's, and the slowest call:

| threads | main | numpy |
|---|---|---|
| 1 | 1.21 to 1.22 s | 0.31 to 0.33 s |
| 4 | 1.17 to 1.22 ×; slowest call 1.50 to 1.55 s | 1.26 to 1.36 ×; slowest call 0.45 to 0.47 s |
| 8 | 1.63 to 1.71 ×; slowest call 2.23 to 2.30 s | 1.68 to 1.82 ×; slowest call 0.64 to 0.70 s |

So the numpy version loses a little more to running beside other threads, 0.1 or so on that
ratio, and still runs each call 3.3 to 3.9 times faster at 4 and 8 threads. On the second form of
P1, before P4, it had lost less than main at 8 threads (1.45 to 1.52 against 1.55 to 1.66); the
difference between the two is within what the load moves.

**`str-regions --threads` over chromosomes 13 to 22 and Y and over GRCh38**, one run per side,
main then numpy, outputs equal at every count:

| input | threads | main | numpy | peak RSS, main → numpy | load before |
|---|---|---|---|---|---|
| chr13–22, Y | 1 | 268.5 s | 94.3 s | 917.7 → 929.0 MB | 4.3, 5.1 |
| chr13–22, Y | 4 | 103.6 s | 34.8 s | 2,352.7 → 2,445.4 MB | 6.4, 7.1 |
| chr13–22, Y | 8 | 80.0 s | 25.2 s | 3,356.3 → 3,834.5 MB | 7.3, 8.3 |
| GRCh38 | 1 | 1,169.6 s | 376.1 s | 1,855.9 → 1,977.3 MB | 8.0, 9.4 |
| GRCh38 | 4 | 389.4 s | 172.6 s | 5,213.4 → 5,093.4 MB | 3.5, 6.9 |
| GRCh38 | 8 | 391.5 s | 98.4 s | 8,326.0 → 8,213.6 MB | 9.0, 8.6 |

Main's GRCh38 run at 8 threads took as long as its run at 4, at a load average of 9, and was not
run again, so the factor of 4 at 8 threads over GRCh38 is likely too high; on the 855 Mb input,
where main's runs scaled as expected, the numpy version is 3.0 times faster at 4 threads and 3.2
at 8. Over GRCh38 it is 3.1 times faster at 1 thread and 2.3 at 4.

**Memory.** On one contig the peak is about 20 MB above main's, at every length from 4 Mb to
chr1: 10.4 MB is numpy's import, and the rest one block's arrays and the stack. Over a whole
reference at 1 thread it was 121 MB above (1,977 against 1,856 MB): the peak of a run over several
contigs depends on how the allocator reuses the memory of the contigs before (the first review),
and numpy's arrays come from the system's allocator, not Python's; this was not looked into
further. At 8 threads the peak over chromosomes 13 to 22 and Y was 478 MB above main's here and
373 MB above in a run of the second form of P1 at 12:27 (3,888.7 against 3,515.7 MB), 11% to 14%;
over GRCh38 at 8 threads it was 112 MB below. Which contigs are segmented at once depends on
timing, and the faster finder changes it; one run per side cannot say more.

## What numpy cannot take

Over GRCh38 on one thread, after the change (`stages_str_finder.py`, 377.5 s in all):

| step | seconds | share |
|---|---|---|
| step 1 and step 3 | 142.3 | 38% |
| step 7, absorption | 178.2 | 47% |
| steps 4 to 6 | 32.6 | 9% |
| reading, step 2, step 8 | 24.4 | 6% |

In main it was about 880 s, in the first review's run of the step timer over GRCh38, which
took 1,104 s in all. Within step 1, on 16 Mb (timers around single functions,
two runs, 0.76 s): the Ruzzo–Tompa loop in Python takes 0.23 s, 30%, for the runs that need the
stack; the array work per base about 0.33 s; the candidate work about 0.16 s; the joining and
sorting of each block's subsequences the remaining 0.04 s. The loop runs once
per run that needs the stack, in order, with its stack, so numpy cannot take it; nor does the
per-base work get cheaper with numpy than one pass of each operation per block.

**A compiled extension** (Rust through PyO3, or C) would take all of step 1, the loop included, as
pvc's own `tandem_repeat.rs` runs it; pvc_eval would then build and ship compiled code for each
platform, which the owner has kept out so far. By this review's numbers it would take at most
step 1's 142 s of a GRCh38 run; absorption's 178 s would remain. Not recommended before the sweep.

**Absorption** is now the largest step. Its fixed point repeats passes until none absorbs anything,
up to about 30 passes on some chromosomes: chr4 spends 27.3 s in it, of 46.5 million calls, against
9.2 s in step 1. pvc's single sweep (`absorb_by_sweep`), which the first review measured at chr4
23.8 → 0.91 s, and chr18 8.03 → 0.52 s (15 to 26 times less), would take it to about a twentieth, so a GRCh38
run on one thread from about 376 s to about 210 s. It was not built here: the brief of this
review kept it for a change of its own.

## The owner's decisions

**The merge of `perf/str-finder-numpy` into main: recommended.** It makes `str-regions` about 3
times faster on one thread and 2.3 to 3.2 times faster on 4 and 8 threads, gives the same outputs on
GRCh38, chr1, the 855 Mb input and every random case tried, and costs:

- numpy 2.5.4 as a runtime dependency, below 3 (c8ad848), with the GIL off, which
  `tests/test_free_threading.py` checks on every run of the tests; 10.4 MB and 24 ms at start;
- about 20 MB more peak per contig;
- 464 lines added and 289 removed in `segments.py` (1,370 → 1,545 lines), a step 1 that no longer
  reads line by line as pvc's, and main's step 1 kept as a 367-line oracle in the tests; 375 lines
  of tests added;
- a rule to keep: numpy's module functions contend across threads on 3.15t, so per-block code
  uses array methods (the comment beside `_BLOCK_BASES`).

The branch also holds the spec's limit on the reward and the penalty, the two benchmark scripts and
this report. Its commits change only the finder, its tests, `cli.py`'s two `--help` lines,
`doc/specs/str_finder.md`'s table, `benchmarks/` and this report.

**The peak at 8 threads.** On chromosomes 13 to 22 and Y it was 11% to 14% higher than main's in
two runs, and over GRCh38 1% lower in one. The cause was not found, and one or two runs per side cannot
separate it from which contigs happened to be segmented at once. It does not need to hold up the
merge: 8 threads is not the default, and over the whole genome the peak was lower. If users run
with 8 threads where memory is tight, several runs per side would settle it. The default of 4 threads showed no such rise (2,445 against 2,353 MB on the 855 Mb input,
5,093 against 5,213 MB over GRCh38).

**pvc's single sweep for absorption: recommended next.** About 165 s of the 376 s left over GRCh38
on one thread; about 60 lines and the fixed point kept as the oracle, as the first review
describes.

**A compiled step 1: not recommended now.** At most 142 s of 376 s, for compiled code to build and
ship.

**Also for the owner to know.**

- `doc/specs/str_finder.md`, "What it holds", still says that a measurement decides whether the
  finder needs numpy. A sentence that says it uses numpy, in blocks of 65,536 bases, holding about
  50 bytes per base of a block per thread, would bring it up to date; it was not changed here.
- `tests/test_cli_compare.py::test_warnings_and_info_lines_go_to_stderr` fails about one run in
  three on main too (2 of 6 runs at 96b3938): log lines from threads arrive in another order. It is
  outside the finder.

## For whoever next revises the skills

The owner can stop here.

- Without the owner at the terminal there were no Tachyon captures. Timers around single functions
  (`time.perf_counter` in a wrapper set from a script) gave the shares within step 1 without
  cProfile's inflation of functions called a million times; they are cheap when the wrapped
  functions are called per block.
- The first array design passed every test, every benchmark and a code review on default
  parameters, and needed 170 bytes per base under a scoring the spec allows. The memory reviewer's
  bound was right for the defaults only. A check of memory at the edges of the parameter table is
  worth a line in the memory checklist.
- `paired.sh` now alternates the order and records CPU time and load; under load the CPU-time ratio
  matched the wall-time ratio within 2% in every call, which is a check that a pair was not
  distorted by descheduling.
- numpy's module functions contend across threads on 3.15t and its methods do not: worth a line in
  `python.md` or the parallelism checklist once numpy is merged.
