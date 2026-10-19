# Report: the STR finder

*19 October 2026. Under way. The work of `doc/plans/str-finder.md`, which builds
`doc/specs/str_finder.md`, on the branch `str-finder`.*

The plan started from commit bec2800, on which `uv sync --locked`, `ruff format --check`,
`ruff check` and `ty check` pass and `pytest -q` gives "213 passed in 1.38s".

## Milestone 1: the segment contract and the candidates

Commits 90a3c39 to 932d773. `src/pvc_eval/segments.py` and `tests/test_segments.py`.

- **The contract**, `RegionType` and `GenomeSegment` as `doc/spec.md` §4 gives them, in commit
  90a3c39 alone, sent to the loci session on 19 October 2026.
- **`StrFinderParameters`**: each field tested at its limit and one past it; each raises
  `ValueError` naming the field, its value and the limit.
- **Step 1**, with the 17 ported tests of pvc's `tandem_repeat.rs`, all passing.
- **The memory bound.** On 4 million random bases, at the default copy floors, step 1 peaks at
  8.0 MB under `tracemalloc`; holding every candidate, in one block, it peaks at 965.6 MB (Apple
  M5 Pro, CPython 3.15t). The test takes about 34 s and is marked `slow`.
- Checks on 932d773: ruff and ty clean, `pytest -q` "301 passed in 37.56s" (213 before).

**Changed in the plan.** Step 2, the contig-end flank, moved to milestone 2: it needs the cut of
step 6, and it only drops candidates, so it can run after step 1.

**What the review found and is fixed.** Five reviewers, one per category.
- *Tests, a blocker.* Four mutants of step 1's drop rule, and one that tests a candidate before its
  divisors are final, changed the output and passed every test, because two tests took their
  expected values from step 1 itself. Six sequences now carry pvc's candidates as literals, from a
  reviewer's transcription of pvc's code, over block sizes 1 to 60, and each mutant fails them.
- *Design.* The Ruzzo–Tompa "segments" were renamed subsequences, so that the word keeps the
  spec's meaning; two functions that nested four levels were split; internal ends made half-open.
- *Errors.* Three messages of the parameter checks corrected.
- The spec and correctness reviewers found nothing. Each compared step 1 with its own
  transcription of pvc's `find_tandem_repeats`: 7,800 random and adversarial sequences, every block
  size from 1 to 4,096, no difference.

**Not taken.** More checks on `GenomeSegment`: §4 states only that the end is after the start, and
the contract is shared. The duplicate default tests: they are pvc's, ported by name.

**What the owner should know.**
- pvc's test `a_period_multiple_re_detection_is_dropped_by_geometry` never reaches the rule it is
  named for, at bbcc9f4: its period-6 stretches are one copy long or have a motif that is not
  primitive. A test that does reach it was added.
- On input made of long runs of one base, the queue of candidates waiting for their test grows by
  about one per 500 bases, because period 1 keeps one subsequence open. The output stays right, and
  pvc has the same bound. Whether real sequence does this is for the performance work.

## Milestone 2: from candidates to STR segments and bundle members

Commits 098c345 to 221d69c. Steps 2 to 6 in `src/pvc_eval/segments.py`, private; pvc's
fixtures in `tests/data/tandem_repeat/`, byte for byte (`cmp`), with `tests/data/README.md`.

- **The 44 ported tests** (42 of `segment_criteria.rs`, 2 of `repeat_catalog/row.rs`), by pvc's
  names, all passing, plus tests of the port's own for rules no ported test reached.
- **pvc's fixture**, `synthetic_ref.fa`, now pins pvc's counts with `params()`: 15 STR segments on
  ctg1 and 2 on ctg2, no bundle members.
- Every expected value added beyond pvc's own tests equals the output of a line-by-line Python
  transcription of pvc's functions at bbcc9f4, which the tests' comments name.
- Checks: ruff and ty clean, `pytest -q` "370 passed".

**Changed in the plan.** None beyond step 2, moved here from milestone 1. The plan said
`minimal_trim_returns_none_without_two_clean_copies` guards the cut's early failure; its input has no
cut either way, and a test of the port's own guards it instead.

**What the review found and is fixed.**
- *Tests, a blocker.* A candidate under `min_score`, or with a compound motif, could be dropped after
  bundling instead of before, and every test passed, because each case of the test the plan named
  held one candidate. Two tests now put such a candidate beside a clean repeat. The review also
  found three bounds of the cut and the order of the rejection reasons unpinned; each now has a test
  that a mutant fails. One of these needed a second round: a restructuring changed how the motif is
  read, and the first version of the test could no longer fail.
- *Spec.* The rule that a contig-end flank of 0 turns step 2 off had no test; a ported assertion
  had been dropped that the test list keeps.
- *Correctness and spec.* An empty candidate passed the copy floor at a floor of 0, where pvc drops
  it; steps 3 and 4 now test that the end is after the start, as pvc's do.
- *Design.* One rule, "the cut, or else the detected span", was written twice and is now one
  function; "tract", a word the spec forbids, was renamed after what each value is.
- The spec and correctness reviewers compared steps 2 to 6 with two independent transcriptions of
  pvc on 125,000 random and adversarial cases, and steps 1 to 3 with pvc's own binary on 400
  contigs: no difference.

**Not taken.** Purity is kept in 64 bits where pvc keeps 32: the decision can differ only for a cut
of more than about 10 million bases, and pvc_eval never writes purity out. Step 2 keeps a candidate
shorter than its period, where pvc refuses it: step 1 never makes one.

**What the owner should know.** Nothing that changes a result. The comparison with pvc's binary on
400 random contigs agrees for steps 1 to 3; milestone 5 compares the whole finder.

## Milestone 3: the segmented genome of a contig

Commits 6d179b0 to e6961f6. `segment_contig`, steps 7 and 8, and a private filter of STR regions by
requested spans for `str-regions --regions`.

- **The 31 ported tests** of pvc's `region_typing/mod.rs`, the golden-catalogue test among them, run
  with a contig-end flank of 0 and the periods scanned equal to the periods admitted, as pvc's
  scanning route. All pass.
- **The spec's own cases**: the four worked cases, the 120-A case (periods scanned 1 to 6, admitted
  2 to 6: a satellite), the copy floors at their boundary, an empty contig, a contig shorter than the
  period, a property test that the segments of a contig run consecutively from 0 to its length, and
  four threads giving what one gives. Their expected values come from pvc's binary at bbcc9f4.
- `segment_contig` reads and resolves the contig when it is called, so its errors come at the call.
- Checks: ruff and ty clean, `pytest -q` "450 passed".

**Changed in the plan.** The plan left open how pvc's four tests with requested spans are ported;
they run against the private filter that `str-regions --regions` will use, without pvc's assertions on
generic regions.

**What the review found and is fixed.**
- *Spec and tests, a blocker found from both sides.* No test checked that a bundle's hull ends at the
  furthest end of its members rather than its last member's; changed that way, the code cut 7 of 600
  contigs' bundles short and every test passed. A contig whose bundle has a member inside another now
  pins it, with pvc's segments.
- *Errors and tests, a blocker found from both sides.* `segment_contig` passing on a bad FASTA byte or
  a closed FASTA had no test; a version that returned no segments passed.
- *Tests.* Two more rules of step 7 had no test that could fail: a repeat run with one span inside
  another, and the distance at which bundles are regrouped after absorption.
- *Design.* The `--regions` filter took the private segment type and compared every segment with
  every span; it now takes `GenomeSegment` and walks sorted spans once.
- The code matched pvc in every comparison: 4,500 contig runs against pvc's binary at two settings,
  and 51,700 layouts and contigs against two transcriptions of pvc.

**Not taken.** Three comments keep pvc's names `RepeatInterval` and `contig_with_tract`, so that a
reader can find pvc's code.

**What the owner should know.**
- Absorption is the plain fixed point the brief asked for. On each pass it compares every feature
  with every satellite; pvc counted 9.3 billion such comparisons on GRCh38 before it moved to a
  sweep. It is fine at the size of milestone 5's sequences, and is the first thing the performance
  work will meet on a whole chromosome.
- pvc's two tests of absorption that the plan named do not fail when absorption stops after one pass;
  a test of the port's own does.

## Milestone 4: `read_str_regions` and the `str-regions` command

Commits 2835af6 to d16ba92, with main merged in at d06e02c. `src/pvc_eval/cli.py`, pvc_eval's first
command line, with the `pvc_eval` entry point; `read_str_regions` and `write_str_regions` in
`src/pvc_eval/segments.py`; `tests/test_str_regions.py`.

- `uv run pvc_eval str-regions --reference ref.fa --output str.bed [--regions r.bed]` writes the BED
  of the spec; its rows on the worked cases equal pvc's at bbcc9f4. `read_str_regions` reads it back
  to the same STR regions.
- The `--help` and the error messages were read by the `first-reader` subagent as a user would; six
  of seven messages were clear at once, and the help now says what the BED's columns are and what an
  STR segment, a bundle and a satellite are.
- Checks on d16ba92, with main merged: ruff and ty clean, `pytest -q` "679 passed". The loci
  session's test of `build_loci` on `segment_contig`, skipped on main, runs and passes here.

**Changed in the spec during the milestone**, each in a commit of its own before the code:
- `write_str_regions` is public: the command needs one entry point into `pvc_eval.segments`.
- The errors of the command: an `--output` that cannot be written is checked before any contig is
  segmented; a command line that cannot be parsed keeps argparse's message and exit status 2.
- How a BED is read: a line is a header only when its first field is `track` or `browser`, so a
  contig may be named `trackA`; whitespace-only lines, empty trailing columns and a byte-order mark
  read as nothing.
- `--regions` is read by `pvc_eval.bed.parse_bed`, the reader the loci brought into main, so a BED
  means the same to `str-regions` and `compare`, gzipped or not.

**What the review found and is fixed.** No blocker.
- *Errors.* Most faults of `--output` ended in a traceback, one of them after the whole genome was
  segmented; a contig whose name starts with `track` or `browser` was skipped without a word.
- *Tests.* Rows under `--regions` were not tested to come in `.fai` order (a writer following the
  regions file passed, and wrote a BED `read_str_regions` refuses); an empty motif, three error
  messages and Windows line ends had no test.
- *Design.* Each BED line looked its contig up in a tuple, 16 times slower than a set on a
  3,366-contig reference.
- The spec reviewer compared the command's BED with pvc's binary on 2,300 generated contigs, with and
  without `--regions`: identical.

**Not taken.** Moving the writer out of `pvc_eval.segments` (`python.md` keeps file writing in the
table and the command line): the spec places it there, and the move can come with `compare`'s
command line. The pattern that turns field names into option names in the command's messages: a test
pins every one.

**What the owner should know.** An empty `--regions` BED stops the run. For an hour of this
milestone the spec said it wrote a BED of header lines only; adopting the loci's shared reader, which
refuses one, made the two commands agree.
- A `--regions` BED that starts with a byte-order mark is refused, and the remedy the message gives
  would drop its first region. The fault is in `pvc_eval.bed.parse_bed`, the loci's shared reader;
  it was reported to the loci session on 19 October 2026 and not changed here.
