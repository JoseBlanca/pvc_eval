# Plan: the STR finder

*19 October 2026. Approved by the owner on 19 October 2026; under way. It builds `doc/specs/str_finder.md`, approved by
the owner on 19 October 2026, on the branch `str-finder`. Its report will be
`doc/reports/str-finder.md`.*

The plan builds pvc_eval's STR finder, stage 1 of `doc/spec.md` §4: the module `pvc_eval.segments`,
its tests ported from pvc, the `pvc_eval str-regions` command, and the check against `pvc dev
type-regions` over GRCh38, in five milestones that follow the order of the spec's steps.

## In and out

In: everything the spec names, `doc/specs/str_finder.md` and `doc/specs/str_finder_tests.md`.

Out:

- The measurement of time and memory on GRCh38, and any faster step 1, numpy among them. They come
  after this plan, with the `performance` skill, as the brief for this port says.
- The comparison with pvc over GRCh38. It comes after the performance work (the owner, 19 October
  2026); the spec, "The same check over GRCh38", holds what it needs. This plan holds the
  finder to one bound of its own, in milestone 1.
- Threads in `str-regions` (the spec's "Not in this spec").
- `compare` and the loci, which read the segments. Another session builds them against
  `doc/spec.md` §4.


## Before the first milestone

- `uv sync --locked`, then the five checks of the `coding` skill pass: on the starting commit,
  a62de70, `uv run pytest -q` gives "213 passed".
- pvc's source at bbcc9f4 in `/Users/jose/devel/pop_var_caller`, for the tests to be copied from. Its
  HEAD is that commit; `git -C /Users/jose/devel/pop_var_caller rev-parse HEAD` says so.

## Milestones

- [x] 1. The segment contract and the candidates

**What it gives.** `pvc_eval.segments` with `RegionType`, `GenomeSegment` and `StrFinderParameters`,
and step 1, the candidates of a sequence. The loci session can import the module.

**The spec it builds.** `doc/spec.md` §4, "What stage 1 yields"; `doc/specs/str_finder.md`, step 1,
"Parameters" and `StrFinderParameters`.

**Deliverables.**

- `RegionType` and `GenomeSegment` exactly as `doc/spec.md` §4 gives them, in a commit of their own.
  The session sends that commit to the loci session, `pvc-eval-e2` in `ListAgents`, before the
  milestone goes on, as the brief for this port asks.
  Check: the commit's diff holds only `src/pvc_eval/segments.py` and its test, which imports both
  names, builds a segment and fails on `end` ≤ `start`.
- `StrFinderParameters`, with the limits of the spec's table. Check: one parametrised test per field
  at its limit and one past it, each value past the limit raising `ValueError` naming the field.
- Step 1 as a private function, with the 17 ported tests of `tandem_repeat.rs`
  (`doc/specs/str_finder_tests.md`). Check: `uv run pytest tests/test_segments.py -q` runs them all, by
  their names, and they pass.
- Step 1 does not hold every candidate it finds (the spec, "What it holds"). Check: a test runs step 1
  with the default copy floors on 4 million random A, C, G and T bases under `tracemalloc` and its
  peak stays under 30 MB. Holding every candidate would take about 75 MB, by estimate: 750,000
  candidates of period 1, 3 in 16 positions, at about 100 bytes each. The implementer measures that
  peak once, on CPython 3.15t, and writes it into the test's comment.

**What could go wrong.** The ported tests check step 1's whole output, with no copy floor. The bound
needs the copy floor applied as candidates are found. So step 1 takes the floors, and drops a
candidate under them as soon as nothing needs it. `segment_contig` passes the parameters' floors.
With floors of 0 it gives pvc's `find_tandem_repeats` exactly, and the ported tests run it that way.
Step 2, the contig-end flank, only drops candidates, so it runs after step 1, in milestone 2
(changed during milestone 1: it needs the cut of step 6's `minimal_trim`). A candidate under the floor still drops a longer-period one (the spec, step
1, and "What it holds"), so dropping it too early gives more candidates than pvc, silently. The
ported test `a_period_multiple_re_detection_is_dropped_by_geometry` and the brute-force tests guard
it in part; on GRCh38, milestone 5 would show it as STR regions pvc does not have. pvc's brute-force Ruzzo–Tompa oracle may be slow in
Python on 2000 cases of up to 160 scores: the inputs stay, and the test is marked `slow` if it takes
more than 10 s.

- [x] 2. From candidates to STR segments and bundle members

**What it gives.** Steps 2 to 6: which candidates are dropped, which bundle, and which become STR
segments, with the reason of each rejection.

**The spec it builds.** `doc/specs/str_finder.md`, steps 2 to 6.

**Deliverables.**

- The 42 ported tests of `segment_criteria.rs` and the two of `repeat_catalog/row.rs`. Check: `uv run
  pytest tests/test_segments.py -q` runs 44 more tests than after milestone 1, by their names, and
  all pass.
- pvc's fixtures copied into `tests/data/tandem_repeat/`, byte for byte, with a line each in
  `tests/data/README.md`. Check: `cmp` against pvc's copies; the ported
  `differential_vs_production_build_loci_on_scanner_output` reads `synthetic_ref.fa`.

**What could go wrong.** A test of step 4 moved after bundling, which would let a compound or
out-of-range candidate bundle its neighbour; and the copy floor, counted twice, on the detected span
(step 3) and on the cut (step 6), counted once. `the_copy_floor_counts_the_trimmed_tract_not_the_detected_one`
and `each_rejection_reason_names_the_gate_that_fired` guard it. `minimal_trim` fails before it
compares; a test of the port's own guards it (pvc's `minimal_trim_returns_none_without_two_clean_copies`
has no cut either way, found during milestone 2).

- [ ] 3. The segmented genome of a contig

**What it gives.** `segment_contig`: the segments of one contig of an `IndexedFasta`, steps 7 and 8
on the result of milestones 1 and 2.

**The spec it builds.** `doc/specs/str_finder.md`, steps 7 and 8, the worked cases, and
`segment_contig` under "The Python interface".

**Deliverables.**

- The 31 ported tests of `region_typing/mod.rs`, the golden-catalogue test among them. Check: they
  run by name and pass.
- The spec's own cases at `segment_contig`: the four worked cases with their segments as literals,
  the copy floors at their boundary, an empty contig, a contig shorter than the period, and the
  segments of a contig consecutive from 0 to its length, on the fixture and with hypothesis. Check:
  they run by name and pass; the worked cases fail if step 2 is turned off (`near_end` gives a
  bundle) or if the copy floor runs after bundling (`four_a` gives a bundle).
- `segment_contig` called from four threads on the two contigs of the fixture gives what one thread
  gives (`python.md`, "Threads").

**What could go wrong.** Absorption is a fixed point, and a pass that stops one round early gives
satellites too short with no error; `absorption_tracks_the_bundle_radius` and
`a_microsatellite_beside_a_satellite_is_absorbed_into_it` guard it. The periods scanned and the
periods admitted are two ranges; the ported tests set them equal, so only the spec's 120-A case
would show them confused, and it is added as a test here (`scan_min_period` 1, `min_period` 2: a
satellite).

- [ ] 4. `read_str_regions` and the `str-regions` command

**What it gives.** A user can run `uv run pvc_eval str-regions --reference ref.fa --output str.bed`
and read the BED back with `read_str_regions`, for `compare --str-regions` later.

**The spec it builds.** `doc/specs/str_finder.md`, `read_str_regions` and "The `str-regions`
command".

**Deliverables.**

- `src/pvc_eval/cli.py` with the `str-regions` subcommand, its options and its `--help`, and the
  `pvc_eval` entry point in `pyproject.toml`. Check: `uv run pvc_eval str-regions --help` lists every
  option with its default; a test runs the command on the worked cases' FASTA in `tmp_path` and
  compares the BED with the spec's segments, as literals.
- `--regions`: the two cases of the spec, and each error of the command with its message. Check: the
  tests run by name and pass.
- `read_str_regions`: the worked cases' BED read back to the same STR regions, a three-column line to
  `None`, and each refusal of the spec, with its line number in the message. Check: likewise.
- The `--help` and the messages go to the `first-reader` subagent, read as a user of pvc_eval.

**What could go wrong.** The spec's own reader and writer agreeing with each other while both are
wrong; the literals of the worked cases, which came from pvc, guard it.

- [ ] 5. The whole-port check against pvc

**What it gives.** The check of `doc/specs/str_finder.md`, "The whole-port check", on generated
contigs and slices of GRCh38, and its first run.

**The spec it builds.** That section.

**Deliverables.**

- pvc built at bbcc9f4 from `git archive` into a scratch directory, as the spec says. Check: `pvc dev
  type-regions --help` runs, and `git -C /Users/jose/devel/pop_var_caller status --short` shows nothing
  new.
- The test, marked `slow`, over the 600 generated contigs at both settings and over the slices of
  GRCh38. Check: with neither `PVC_EVAL_PVC` nor `PVC_EVAL_GRCH38` set, `uv run pytest -m slow -q -rs`
  reports each part skipped and names what is missing; with both set, it compares, and the slices
  hold each of the three types.
- The first run: identical, or the differences, each traced to the step that makes it and fixed, with
  a test of its own, before the milestone is done. Where pvc's code and the spec disagree, the spec is
  corrected first, in a commit of its own, as the `spec` skill says, and the owner is asked when the
  correction changes a segment. Check: the test passes, and the report gives how long it took and how
  many STR regions it compared.

**What could go wrong.** Slices too small to hold a satellite or a bundle, so that a step is never
compared; the assertion that each type appears guards it. Generated contigs whose planted repeats
never come near each other or a contig end; the generator plants them so, and the test counts the
bundles and the candidates dropped at step 2 and asserts there are some.

## At the end

The five checks of the `coding` skill pass, the whole-port check passes, and the report
`doc/reports/str-finder.md` gives what was built, the counts of tests, and the run of milestone 5.
Then the `performance` skill measures the finder on GRCh38, for the owner's decision on numpy.
