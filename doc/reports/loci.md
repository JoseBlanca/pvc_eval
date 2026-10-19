# Report: the loci

*Done, 19 October 2026. The work of `doc/plans/loci.md`, which builds `doc/specs/loci.md`, on the
branch `loci`.*

The plan is done. pvc_eval now builds its loci: `build_loci` turns the VCF records of file A and
file B and the STR regions into loci as it reads them, and the kind, the stratum and the
`--regions` selection of a locus are read from a locus. The checks give 360 passed and 1 skipped.
The skipped test runs `build_loci` on the STR finder's `segment_contig`, which is not on this
branch; it passed against a copy of it from the `str-finder` branch. `build_loci` gave the same
loci as a brute force of `doc/spec.md` §5.3, written from its four steps alone, on 60,000 random
inputs of both files, two contigs, STR regions and chain distances from 0 to 10. Every time and
memory figure below was measured on the Apple M5 Pro, with CPython 3.15t, on 19 October 2026.

**What the merge brings into main.** The branch holds, besides this plan's code and tests:

- three commits to `doc/spec.md`: `baf1be3` (`build_loci` takes `contigs`), `84d6357` and
  `c0fcb6a` (loci within the chain distance of each other are one, §5.3 step 4). The other
  sessions read the contracts of §4 from main, so until the merge they know these only from
  messages;
- the shared modules, cherry-picked: `pvc_eval.segments` from `str-finder`, `pvc_eval.locus` and
  `Kind` from `reconciliation`. When those branches are merged, git sees the same changes on both
  sides; a conflict, if any, is in those files alone.

Nothing on this branch waits on the other branches: the test that needs `segment_contig` is
skipped until it exists.

**The owner's answers.** On 19 October 2026 the owner agreed to the three points the work raised:
`--regions` and `--strata` read gzipped BEDs and refuse an empty one (milestone 4, below), and
`doc/spec.md` §5.2 now says its cases join "by overlap alone" through the span (`300e085`).

**One point for the owner.** `parse_vcf` lets `zlib.error` through, so a gzipped VCF whose
compressed data is damaged ends in a traceback; the BED reader of milestone 4 catches it. The two
readers hold one copy each of the gzip handling. Making it one private function, used by both,
fixes the VCF reader too, and changes `vcf.py`, which the handoff keeps for the owner. The same
function would drop a leading UTF-8 byte-order mark from a VCF, which `parse_vcf` keeps on its
first header line today, as the BED reader did until milestone 5.

**What went wrong.** While removing its reviewers' worktrees, this session removed 14 worktrees
of the other sessions' reviewers and deleted 15 `worktree-agent-*` branches, by a glob. No commit
was lost: every commit they held is on `reconciliation`, `reconciliation-samples` or `str-finder`,
and the branches had none of their own. The STR finder session reported that its reviewers there
had finished. The rule is now in this session's memory: a worktree is removed by its exact path.

What is left open after the merge: one trim of an ALT against its REF is written twice, in the
kind and in reconciliation, to be made one when the branches meet; and the test on
`segment_contig`, which runs once the STR finder is merged.

## Before the first milestone

- `uv run python -c "import sys; print(sys._is_gil_enabled())"` printed `False`.
- On `cc78799` the four checks of `python.md` gave "8 files already formatted", "All checks
  passed!" twice and "213 passed".
- The shared modules were cherry-picked: `pvc_eval.segments` from the STR finder's `90a3c39` (here
  `89dcb12`) and `pvc_eval.locus` from the reconciliation session's `f603bec` (here `f543ad9`).
  The checks then gave "11 files already formatted", "All checks passed!" twice and "217 passed",
  the 4 new tests being `tests/test_segments.py`'s.

## Milestone 1: `build_loci` on hand-built records

Done, from `2be9d5f` to `ea49f42`. `build_loci` is in `src/pvc_eval/loci.py`, with its tests in
`tests/test_loci.py`.

**Deliverables.**

- `pvc_eval.locus` was not written here: the reconciliation session's `f603bec` and its checks,
  `18ce521`, were cherry-picked (here `f543ad9` and `1f0a785`).
- `uv run pytest tests/test_loci.py -q` gave "60 passed": the spec's 50 cases, the 5 the review
  added, and the default constant, the property test, the two memory tests and the linear-time test.
- The property test fails when the reach is kept from the last widening instead of computed for
  each next record: "3 failed, 49 passed" on the first version, with the shrunk case two SNPs of
  file B at (0, 1) and (4, 5), STR regions at (0, 2) and (3, 4) with a GENERIC segment between
  them, and a chain distance of 1, where `build_loci` gave (0, 2) and the rules give (0, 5).
- The four checks on `ea49f42`: "14 files already formatted", "All checks passed!" twice, "283
  passed".
- Beyond the plan: the correctness reviewer's comparison of `build_loci` with a brute force of
  §5.3 steps 1 to 4 gave 0 differences on 60,000 random cases after the fixes, as before them; the
  spec reviewer's own brute force gave 0 on 140,000 before them. Branch coverage of `loci.py` was
  100%.

**What the review found and is fixed.**

- A GENERIC segment out of order narrowed a locus with no error, because the reading of segments
  stops at a segment of any type. Every segment's order is now checked; the spec said GENERIC
  segments were not checked and was corrected first (`5031cdc`, `974b3be`). Three reviewers found
  it.
- Records with no ALT held over a satellite took time growing with the square of its length:
  5,000, 10,000 and 20,000 of them took 2.56, 10.11 and 40.33 s; after the fix 0.01, 0.01 and
  0.02 s. A test marked `slow` holds 50,000 under 5 s.
- A blocker: no test saw STR regions dropped; with them held forever all 52 tests of that version passed while the
  peak grew from 4,296 bytes to 2,508,752 at 10,000 loci. A memory test now catches it, in 6.4 s:
  its first version took 627 s to fail, which the second look found.
- Three rules guarded only by the property test now have fixed cases: a satellite exactly at the
  chain distance, one on the next contig, and a record with no ALT kept over an STR region.
- A contig named twice in `contigs` now raises; names, docstrings and a test's regex were fixed.

**Not taken.** The property test's brute force uses the same touching and overlap formulas as the
code: they are the spec's formulas, and the brute force imports nothing from `loci.py`.

**What the owner should know.** While cleaning up, this session removed 14 agent worktrees of the
other sessions and deleted 15 `worktree-agent-*` branches by a glob. No commit was lost, since
every commit they held is on `reconciliation`, `reconciliation-samples` or `str-finder`, and the
STR finder session reported no harm; the other sessions were told.

**Tokens**, from the Agent tool's results: implementation 134,000; review 347,000 (spec 77,000,
tests 107,000, correctness 62,000, errors 36,000, design 65,000); fixes, two rounds, 58,000;
second look 29,000.

## Milestone 2: `build_loci` on the spans `parse_vcf` computes

Done, from `ba5b9f9` to `a0f5162`, in `tests/test_loci_from_vcf.py` alone.

**Deliverables.** `uv run pytest tests/test_loci_from_vcf.py -q -rs` gave "5 passed, 1 skipped",
the skip being "segment_contig is not in pvc_eval.segments yet". Every span `parse_vcf` computed
was the one the spec's case assumes. Run against the STR finder's `segments.py` of `6d179b0` in a
scratch copy, the skipped test passed too, with one locus (24, 46) holding the AT repeat. The four
checks on `a0f5162`: "15 files already formatted", "All checks passed!" twice, "288 passed, 1
skipped".

**Changed in the plan.** The check was 4 passed; it is 5, with the case the review added.

**What the review found and is fixed.** Two reviewers found that the two cases of `doc/spec.md`
§5.2 do not show the loci relying on the span: their REFs already touch, and touching spans chain
at any chain distance, so a `build_loci` reading REFs instead of spans passed every active test. A
case the span alone joins was added to the spec (`ab3c064`) and to the tests: over `GCAGAGAGTC`,
`21 CAG C` in file A and `28 T C` in file B are one locus at chain distance 0. It fails when
`build_loci` reads REFs. Two nits, the comments' coordinates and the reference helpers, are fixed.

**What the owner should know.** `doc/spec.md` §5.2 says its second case would fall into two loci
without the span. With the chain distance, touching spans always chain, so the sentence no longer
holds, though it still does "by overlap alone", the phrase §5.3 uses for the traded events. That
correction is the owner's: this session changes doc/spec.md only with the owner's agreement.

**Tokens**: implementation 88,000; review 92,000 (spec and tests 40,000, correctness and design
52,000); fixes 22,000.

## Milestone 3: kind, stratum and regions

Done: the work of `ba5b9f9..dcb6286` on a branch of its own, `loci-m3`, run beside milestone 2, then
cherry-picked onto `loci` (`3d6c366..c26bec1`) and fixed there (`3abdd98..858b7a7`).

**Deliverables.** `locus_kind` beside the shared `Kind` (cherry-picked from the reconciliation
session's `4df38de`), `pvc_eval.regions` and `pvc_eval.strata`, and `pvc_eval.bed`, the BED reader
both share. `uv run pytest` gave "12 passed" for test_kind, "20 passed" for test_strata, "9 passed"
for test_regions and "17 passed" for test_bed. The four checks on `858b7a7`: "23 files already
formatted", "All checks passed!" twice, "346 passed, 1 skipped".

**Changed in the plan.** `tests/test_bed.py` was added, holding once the refusals both readers
share, which had been tested twice.

**What the review found and is fixed.**

- A blocker: no test saw the classes of a stratum sorted and given once across intervals; a
  stratum written `zeta,alpha,zeta` passed all 17 tests. A test now catches every such mutant,
  confirmed by a second look.
- A `--strata` BED of 1,000,000 lines (32 MB) retained 451 MB, one set of classes per piece.
  `Strata` now holds one merged `Regions` per class, as a reviewer showed gives the same answers on
  60,000 queries: 90 MB retained, peak 184 MB, 1.0 s to read instead of 1.6 s.
- What the BED readers refuse, now in the spec (`15d0b39`): a directory or an unreadable file,
  which gave a traceback; a space-separated BED, which was told it had 1 column; a class with
  spaces at its ends, which counted as another class, an empty class, and one holding `,`.
- The message for a contig the reference lacks gave an `awk` command that, for a BED naming `1`
  against `chr1`, kept no line, after which every locus would be skipped. It now says to rename.
- Six more missing tests (a nested region, a contig with no BED line, file B changing the kind, a
  record with no ALT reaching past the locus, and two orders of `pvc_stratum`), and nits.

**Not taken.** The spec's names `locus_kind`, `pvc_stratum` and `overlaps`, which a reviewer found
to be nouns, are the approved interface's. A shared trim of an ALT for the kind and for
reconciliation, which hold one copy each, crosses sessions and waits for the branches to meet. The
non-UTF-8 message names the line its decoder's chunk began at, which is true as it stands.

**Tokens**: implementation 115,000; review 256,000 (spec 57,000, tests 58,000, correctness 36,000,
errors 40,000, design 65,000); fixes 55,000; second look 23,000.

## Milestone 4: gzipped and empty BEDs

Done, from `8b8fb6e` to `26bfdf0`, in `src/pvc_eval/bed.py` and `tests/test_bed.py`: a gzipped or
bgzipped BED is read by its first two bytes, a damaged one is refused with no traceback, and a BED
with no line is refused. `uv run pytest tests/test_bed.py -q` gave "31 passed", 15 of them failing
on `300e085`; the four checks: "23 files already formatted", "All checks passed!" twice, "360
passed, 1 skipped". The review found no defect: every damaged or truncated gzip case it tried gave
the user an error, files close on every path, and each new test fails under its mutant. Not taken:
a one-byte file holding the first byte of a gzip header is told it holds no line, which is unlikely
and still stops the run. The gzip handling copied from `vcf.py` is the point for the owner above.

**Tokens**: implementation 75,000; review, one reviewer for four categories, not given.

## Milestone 5: a BED that starts with a byte-order mark

Done, in `16eff9f`, after the merge into main. The STR finder session found that `parse_bed` kept
a leading UTF-8 byte-order mark on the first contig, so the first line was refused as an unknown
contig. The BED is now decoded as `utf-8-sig`, which drops the mark at the start of the text only.
Three new tests in `tests/test_bed.py` failed on the commit before; `tests/test_bed.py` gave "34
passed", the four checks "27 files already formatted", "All checks passed!" twice and "371 passed,
1 skipped". The review found no defect: 45 error cases gave the same messages before and after,
and a second mark on a block edge of the gzip stream is kept. This commit is on `loci`, not yet on
main.

## The check of the whole plan

`build_loci` over `chr20_5mb_200s.vcf.gz` (18,310 records, 200 samples) against itself, with no STR
regions, read through `parse_vcf` over GRCh38, by a script in the session's scratchpad: 17,641
loci; file A's and file B's records with an ALT summed over the loci, 18,310 each, equal to the
file's; 0 loci with different records in A and B; 47.1 s with tracemalloc on, both files parsed;
traced peak 5.8 MB.

## For whoever next revises the skills

You can stop here if you are the owner deciding the merge.

- A reviewer that copies the worktree to mutate it, `.venv` included, tests the worktree and not the
  copy: the copied environment imports the worktree's `src`. The review of milestone 5 found its
  first mutation runs all passing for that reason. The coding skill's instructions to reviewers
  could say to build a fresh environment in the copy.

- The coding skill's worktree for side-by-side milestones is not said: milestone 3 ran on a branch
  `loci-m3` in a worktree of its own, cherry-picked afterwards, so that two subagents never
  committed to one branch. The skill could say so.
- Cleaning up reviewers' worktrees by a glob removed other sessions' worktrees (see milestone 1).
  The coding skill's "When the plan is done" says to remove "the worktrees of the reviewers"; it
  could say by the exact path the Agent tool returns, never by a pattern, since sessions share
  `.claude/worktrees/`.
- In all three milestones the tests reviewer found rules no test guarded, by mutating the code: 5
  in milestone 1, 1 in milestone 2, 5 in milestone 3. It found the only blockers. Its cost, 58,000
  to 107,000 tokens, was the best spent of the review.
- Milestone 1, at about two days of work, cost 134,000 tokens to implement and 347,000 to review,
  and fitted one subagent; milestone 2, at half a day, cost 88,000 and 92,000.
