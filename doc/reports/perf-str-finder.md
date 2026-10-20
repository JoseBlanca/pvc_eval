# Performance review of the STR finder, 20 October 2026

This report is about the wall time and the peak memory of the STR finder: `segment_contig`
(`src/pvc_eval/segments.py`), which cuts one contig of a reference into generic regions, STR
segments, STR bundles and satellites, and the command `pvc_eval str-regions`, which runs it over
every contig and writes the STR regions as a BED (`doc/specs/str_finder.md`). It is the first
review of the finder. Its measurement decides whether the finder needs numpy, which is a runtime
dependency and the owner's decision. It was made on the branch `perf/str-finder`, from main at
b71c83c, on an Apple M5 Pro (18 cores, 64 GB) with CPython 3.15.0 free-threaded (3.15t, GIL off).
The experiments are closed below, one by one, as they are run.

## The workload and the baseline

The inputs are stretches of GRCh38 chr1 (the analysis set without alts), each written as a FASTA
of one contig by `benchmarks/make_str_slices.py`, all starting at chr1:10,000,000 so that each
longer one holds the shorter ones: 1, 4, 16 and 64 Mb. The finder's cost per contig is not linear
in the contig's length (step 7 below), so the stretches are contiguous rather than several short
slices. `benchmarks/README.md` gives the commands. The workload is

    uv run --no-sync pvc_eval str-regions --reference <stretch>.fa --output out.bed

at the default parameters. Two outputs are kept from the baseline and every experiment must give
them unchanged: the BED, and a sha256 of every segment `segment_contig` yields, generic regions
included (`benchmarks/digest_str_finder.py`), since the BED leaves the generic regions out.

The baseline, at commit fcd71eb (the finder's code of b71c83c), measured on 20 October 2026
between 14:37 and 14:46:

| input | wall time | peak RSS | STR regions (str / str-bundle / satellite) |
|---|---|---|---|
| 4 Mb | 5.956 ± 0.174 s | 57.1 MB | 17,213 / 5,488 / 34 |
| 16 Mb | 24.274 ± 0.213 s | 128.7 MB | 67,662 / 20,284 / 131 |
| 64 Mb | 105.75 s | 393.9 MB | 374,113 rows |

The wall times are the mean ± sd of hyperfine (10 runs on 4 Mb, 5 on 16 Mb, 1 run of
`/usr/bin/time` on 64 Mb); the peak RSS is the median of 3 maximum resident set sizes (1 on 64
Mb). Python with pvc_eval imported takes 27 MB.

**The noise.** On 4 Mb the hyperfine runs fell into two levels 5% apart within one call, in CPU
time as well as wall time, so the core's speed changed under them. So each experiment is timed in
pairs, baseline and candidate run alternately from two worktrees, and judged by the median of
the ratios of its pairs. The same code against itself, at 15:15 with a load average of 2.4 to
3.6 from other sessions, gave a median ratio of 1.0017 over 10 pairs on 4 Mb (single pairs from
0.942 to 1.012) and 0.9954 over 5 pairs on 16 Mb. A gain is claimed when the median is more than
3% from 1 on 4 Mb, or 2% on 16 Mb.

## Where the time and the memory go

**Time.** Two steps of the finder's eight take nearly all of it:

| input | step 1, the candidates | step 7, absorption | all other steps |
|---|---|---|---|
| 4 Mb | 5.56 s | 0.10 s | 0.06 s |
| 64 Mb | 90.34 s | 14.62 s | 1.12 s |

(`benchmarks/stages_str_finder.py`, one run each.) Step 1 grows with the length of the contig.
It scores every base once per period, 6 times, and passes each of the 24 million scores of 4 Mb
through Python, building a frozen dataclass for each of the 6.76 million positive ones. Step 7
grows with the number of passes times the features times the satellites of the contig, which is
close to the square of its length: 2.4 million comparisons of a feature with a satellite on 4 Mb,
56.5 million on 16 Mb and 690 million on 64 Mb, in 3, 5 and 5 passes.

The Tachyon profiles (16 Mb in wall and cpu mode, 64 Mb in wall mode) name the same places, and
the two modes agree, so nothing waits. Read by the time between samples rather than by the count
of samples, they give step 1 75% to 85% and step 7 6% (16 Mb) to 14% (64 Mb). The count of samples
alone overstates step 7 1.6 to 1.7 times, because the sampler fails to read about half of its
attempts during step 1 and almost none during step 7; this was found by the methodology reviewer
from the timestamps of the samples, and matches the step timer.

**Forecast for chr1** (249 Mb, 18.5 Mb of it N), from the 64 Mb stretch, which has the same
features and satellites per Mb of sequence as chr1 (counted on pvc's typing of chr1 at
pvc_eval's floors): step 1 about 330 s, step 7 about 180 to 225 s, so 9 to 10 minutes. GRCh38 one
contig after another: about 100 minutes.

**Memory.** The peak comes at the end of a contig, when it holds at once the bases (1 byte per
base), the candidates, the STR segments and their motifs, the features and the segments: about
5.1 MB per Mb traced by tracemalloc, against 5.5 MB per Mb of RSS. Forecast: about 1.4 GB on chr1,
the largest contig of GRCh38. The garbage collector takes under 1% of the time.

**Threads.** `str-regions` runs one contig after another. `segment_contig` is safe to call from
several threads, but on today's code eight threads each segmenting 1 Mb take 2.5 to 3.3 s each
against 1.47 s for one thread alone, while eight processes take 1.8 s each. The building of
frozen dataclasses does not scale across threads on 3.15t: in a micro-benchmark, a million builds
take 0.30 s on one thread and 1.10 s on each of eight, against 0.15 s and 0.19 s without `frozen`.

## The candidates, in the order they are tried

Every candidate gave the baseline's outputs when its reviewer tried it. P1 to P3 work on step 1,
P4 on step 7, P5 and P6 on the peak memory.

- **P1. Step 1: one step per run of equal-sign scores.** A maximal scoring subsequence never
  starts or ends inside a run of scores of one sign, since no score is 0. So the runs of matches
  of a period are found per block with C-level operations on bytes (an XOR of two translated
  slices and `re.finditer`), and the Ruzzo–Tompa loop takes one step per run: 9.2 million runs
  instead of 24 million scores on 4 Mb. Tried by the algorithms reviewer: step 1 from 5.6 s to 4.3
  s on 4 Mb. Code: about the same size.

  *Kept, commit 7c3294a.* The code adds 65 lines to segments.py and removes 23: the score
  function finds the runs (about 15 lines, with the translation table and a comment on why the
  earlier base is read through it), and the loop takes a run's length. Three tests of the runs
  were added; a deliberate break, N matching N, failed two of them. The BEDs and the digests of
  every segment equal the baseline's on 4, 16 and 64 Mb, and 912 tests pass. On 4 Mb, cProfile
  counts 4,598,463 calls of `_make_open_subsequence` instead of 6,762,676, and the loop takes
  9,199,504 steps instead of 23,999,979. In alternating pairs against the baseline:

  | input | baseline | P1 | median ratio | peak RSS |
  |---|---|---|---|---|
  | 4 Mb, 10 pairs | 5.732 ± 0.017 s | 4.487 ± 0.013 s | 0.783 | 57.1 → 57.0 MB |
  | 16 Mb, 5 pairs | 24.042 ± 0.120 s | 18.942 ± 0.045 s | 0.787 | 128.7 → 128.6 MB |

  The code is somewhat more complicated and the gain is over 20%, so it is kept, with a comment
  at the change that gives these numbers.
- **P2. Step 1: a run alone on an empty stack becomes final without an open subsequence.** It is
  held as plain numbers until the gap after it decides it. The open subsequences built on 4 Mb
  fall from 6.76 million to 113,000. Tried with P1: step 1 2.4 s on 4 Mb. Code: about 15 lines and
  a second path in the loop.

  *Kept, commit 587a4df.* The code adds 35 lines and removes 1: the held run, its type, a helper
  that puts it on the stack, and a docstring paragraph that states the invariant (a run is held
  only while the stack is empty) and why it is made final as the stack would make it. No test
  was added: making a held run final one unit too early failed 6 tests, and not putting it on the
  stack at a block's end failed 21. Outputs equal on 4, 16 and 64 Mb; 912 tests pass. On 4 Mb,
  112,962 open subsequences are built instead of 4,598,463, and step 1 takes 2.13 s instead of
  4.25 s. In alternating pairs against P1:

  | input | P1 | P2 | median ratio | peak RSS |
  |---|---|---|---|---|
  | 4 Mb, 10 pairs | 4.479 ± 0.010 s | 2.312 ± 0.006 s | 0.516 | 57.1 → 57.1 MB |
  | 16 Mb, 5 pairs | 19.156 ± 0.165 s | 10.298 ± 0.097 s | 0.539 | 128.6 → 128.5 MB |

  More complicated, and a gain far over 20%: kept, with the numbers in the docstring.
- **P3. Step 1: `_OpenSubsequence` without `frozen`.** It is a private working object that never
  leaves the function. The frozen form costs 283 ns per build against 135 ns, and is what stops
  threads from scaling. Measured after P1 and P2, which leave few builds, and at one and eight
  threads.

  *Changed into the search for what makes threads contend, and kept as commit 0f6c2d7.* After P2 only 112,962 open subsequences are built on 4 Mb, so unfreezing them could not
  matter, but threads still lost: on 4 Mb, one `segment_contig` per thread, eight threads took
  4.9 to 5.7 s each and eight processes 3.1 s. Trying one change at a time, with outputs equal,
  the cause was the match objects of `re.finditer` in P1's score function: in a micro-benchmark
  of step 1's pieces, `finditer` took 0.36 s on one thread and 4.35 to 4.45 s on each of eight,
  while `translate`, the XOR and `split` scaled. The shared parameters, the tables, generator
  expressions and lambdas were not causes. The runs are now found with two `bytes.find` calls
  over the XOR and a second translation of it, which make no objects; 15 lines added, 5
  removed, about as readable, with a comment that gives these numbers.

  | 4 Mb, seconds per call | 1 thread | 4 threads | 8 threads | 8 processes |
  |---|---|---|---|---|
  | before (P2) | 2.28 | 2.90 to 3.04 | 4.98 to 5.71 | 3.08 to 3.26 |
  | after | 2.06 | 2.44 to 2.64 | 2.85 to 3.81 | 2.83 to 2.98 |

  One thread, in alternating pairs against P2: 2.114 ± 0.007 s instead of 2.317 ± 0.009 s on
  4 Mb (median ratio 0.912), 9.446 ± 0.127 s instead of 10.294 ± 0.223 s on 16 Mb (0.921); peak
  unchanged. Outputs equal on 4, 16 and 64 Mb; 912 tests pass.

  A smaller cause remains. A frozen dataclass's generated `__init__` reads `object` from a
  closure cell that every thread shares, and on 3.15t a shared cell contends: 4 million calls
  through it took 0.18 s on one thread and 3.5 to 3.8 s on each of eight, against 0.21 to 0.27 s
  through a module global. `_Candidate`, built about a million times per 4 Mb, stays frozen:
  unfrozen it gave about 5% on one thread and an uneven gain at eight threads (median of the
  slowest thread 3.89 s against 4.76 s over nine runs), it departs from python.md's rule for
  records, and it would need `unsafe_hash=True`, since a test puts candidates in a set. The
  shared cell is a behaviour of CPython 3.15t that could be reported upstream.
- **P4. Step 7: find the satellites a feature reaches by bisection.** The satellites stay sorted
  and disjoint, so those within the bundle distance of a span are one slice of the list, found
  with two binary searches instead of a scan of all of them. The fixed point and its passes stay
  as they are. On 4 Mb, 0.018 s against 0.059 s; forecast 0.5 s against 14.6 s on 64 Mb. Code:
  about the same, with the invariant of sorted, disjoint satellites stated.

  *Kept, commit d9601a8.* Two binary searches, on the satellites' ends
  and on their starts, replace the scan; the rest of absorption is unchanged, and on 64 Mb it
  still makes 2,377,895 calls in 5 passes. The function is shorter (34 lines changed in
  segments.py), and its docstring states the invariant and why the satellites reached are one
  slice. Twelve tests were added: eleven hand-worked cases at gaps of exactly 14 and 15 bases
  around several satellites, and a hypothesis comparison with a scan written in the test from
  the spec's rule. Each of four deliberate breaks (either search one base off, the satellite
  before the slice not looked at, only one satellite of the slice grown) failed 3 to 11 tests.
  Outputs equal on 4, 16 and 64 Mb; 924 tests pass. Step 7 on 64 Mb takes 1.84 s instead of
  14.62 s. In alternating pairs against P3:

  | input | P3 | P4 | median ratio | peak RSS |
  |---|---|---|---|---|
  | 64 Mb, 3 pairs | 46.677 ± 0.060 s | 34.280 ± 0.050 s | 0.735 | 393.6 → 393.6 MB |
  | 16 Mb, 5 pairs | 9.360 ± 0.024 s | 8.398 ± 0.023 s | 0.898 | 128.5 → 128.5 MB |
- **P5. Step 8 as a generator.** `_make_segments` yields instead of building the list of all
  segments, and `_find_segments` returns before the first segment, which frees the candidates
  while the segments are written. tracemalloc on 4 Mb: the peak from 20.9 to 19.6 MB. Code:
  simpler, 4 lines fewer.

  *Kept, commit 1bc0f7b.* `_make_segments` yields, and `_find_segments` returns its iterator
  once steps 1 to 7 are done, so `segment_contig` still resolves the whole contig when called.
  10 lines added, 13 removed. tracemalloc's peak on 4 Mb: 21.8 MiB instead of 23.1. In
  alternating pairs against P4, the peak RSS on 64 Mb is 362.0 MB instead of 393.6 MB, and on
  16 Mb 128.1 instead of 128.5; the wall time is the same (median ratios 0.998 on 64 Mb, 3
  pairs, and 0.996 on 16 Mb, 5 pairs). Outputs equal; 924 tests pass.
- **P6. Intern the motifs.** Each STR segment and bundle member holds its own copy of its motif,
  of at most 5,460 distinct ones. With P5: the peak 17.9 MB on 4 Mb. Code: one line.

  *Kept, commit ea99a42.* `_get_motif` returns `sys.intern` of the motif, with a comment that
  says the table is bounded and safe for threads. In alternating pairs against P5, the peak RSS
  on 64 Mb is 333.5 MB instead of 362.0 MB, and on 16 Mb 122.8 instead of 128.1, at the same
  wall time (median ratios 0.997 and 1.001). Outputs equal; 924 tests pass.
- **P7. Step 1: the Ruzzo–Tompa helpers inlined, with tuples on the stack.** Measured only after
  P1 to P3, when what is left of their calls is known; it makes the loop longer and harder to
  read, so it needs a 20% gain.

  *Not run: P2 left nothing for it.* After P2, `_make_open_subsequence` is called 63,126 times
  on 4 Mb (4.6 million before), and it no longer appears among the 14 functions with the most
  own time in cProfile.

**Not tried.** A list comprehension in place of a generator in the loop of step 1 (under the
noise on 1 Mb); a cache of `_is_primitive_motif` (about 2% of 4 Mb); running the six periods of
step 1 on threads (a barrier every 4,096 bases and a more complicated loop, for at most the time of
one contig); cutting a contig at its runs of N (no local rule proves a cut exact, and pvc does
not cut either).

**The owner's.** Brought with their numbers once the candidates above are closed:
- porting pvc's single sweep for absorption (`absorb_by_sweep`) in place of the fixed point;
- bundles without their members, which `doc/specs/str_finder.md` asks the port to keep (13% of
  the peak on 4 Mb);
- `str-regions` on threads, one contig per thread, which needs a `--threads` option and so a
  change to the spec's interface;
- numpy.
