# Plan: the loci

*Status: approved by the owner on 19 October 2026; under way. It builds `doc/specs/loci.md`,
approved by the owner on 19 October 2026, on the branch `loci` in `.claude/worktrees/loci`; the
owner's handoff named the branch `loci`, so it is not `plan/loci`. The work report will be
`doc/reports/loci.md`.*

The plan builds stage 3 of the pipeline, `build_loci`, and then the three things read from a locus
once it is built: whether `--regions` keeps it, its kind and its stratum. Three milestones: the
rules of `build_loci` on hand-built inputs, then `build_loci` on the spans `parse_vcf` computes,
then kind, stratum and regions, which can run side by side with the second.

## In and out

In: everything `doc/specs/loci.md` gives, in `pvc_eval.loci`, `pvc_eval.kind`, `pvc_eval.strata`
and `pvc_eval.regions`, with its tests. The spec has no open point.

Out, and where it goes:

- `compare`, the command line, and how `compare` chains `segment_contig` over the contigs that hold
  records: the spec of `compare`, not written yet. So there is no end-to-end run in this plan.
- `segment_contig` and `read_str_regions`: the STR finder's plan, on branch `str-finder`.
- Reconciliation, which reads the loci: the reconciliation session's plan.

**The shared modules.** Three sessions build on the contracts of `doc/spec.md` §4, and each module
is written once, in a commit of its own, which the others bring in with `git cherry-pick`:

- `pvc_eval.segments` (`RegionType`, `GenomeSegment`) exists: commit `90a3c39` on `str-finder`.
- `pvc_eval.locus` (`Locus`) is in the reconciliation session's milestone 3, and the `Kind` enum
  of `pvc_eval.kind` in its milestone 4. Whichever session needs a module first writes it, exactly
  as `doc/spec.md` §4 or `doc/specs/loci.md` gives it, and sends the commit to the other session:
  the reconciliation session is named `reconciliation`, found with `ListAgents` and written to with
  `SendMessage`. Before writing one, the session checks the other branches for it:
  `git log --all --oneline -- src/pvc_eval/locus.py src/pvc_eval/kind.py`.

## Before the first milestone

- The branch `loci` holds `doc/specs/loci.md` approved, at `cc78799` or after.
- `uv run python -c "import sys; print(sys._is_gil_enabled())"` prints `False`.
- The four checks of `python.md` pass on the starting commit: `uv run ruff format --check`,
  `uv run ruff check`, `uv run ty check` and `uv run pytest -q`. On `cc78799` they gave "8 files
  already formatted", "All checks passed!" twice, and "213 passed".
- `git cherry-pick 90a3c39`, then the checks again: `tests/test_segments.py` runs and passes. A
  conflict there is a stop, since the commit touches only files the loci branch lacks.

## Milestones

- [x] 1. `build_loci` on hand-built records

**What it gives.** `build_loci` in `pvc_eval.loci`, which reconciliation and the wrappers of §8 can
call: the loci of two streams of records and the STR regions, yielded as it reads.

**The spec it builds.** `doc/specs/loci.md`, "Building the loci" (the rules, the order of the
loci, what `build_loci` holds and what it checks) and "The Python interface" for `build_loci` and
`DEFAULT_CHAIN_DISTANCE`; `doc/spec.md` §5.3, steps 1 to 4.

**Deliverables.**

- `pvc_eval.locus`. If the `git log --all` above finds it on another branch, that commit is
  cherry-picked. Otherwise `Locus` is written exactly as §4 gives it, in a commit of its own, sent
  to `reconciliation`. Check: `git show --stat` of that commit touches only `src/pvc_eval/locus.py`
  and its test.
- `tests/test_loci.py`, with records built as `VcfRecord` values whose spans are literals and STR
  regions as lists of `GenomeSegment`, holding the cases listed in the spec's "How it is verified"
  under "At `build_loci`, the cases of this spec", "the cases `doc/spec.md` §7.2 lists for loci"
  (by hand; their runs through `parse_vcf` are milestone 2), "What it reads", and "A property".
  Check: `uv run pytest tests/test_loci.py -q` ends with at least 50 passed and nothing skipped
  or failed; the spec's lists name 50 cases, counting each row of a table and each `ValueError`
  as one. On the starting commit it exits with 4, as the file does not exist.
- The property test compares `build_loci` with a slow version of §5.3 written in the test from the
  rules alone, steps 1 to 4, on hypothesis's cases, and checks the invariants the spec lists. It
  is the check that guards the rules: a wrong rule gives plausible loci and no error. Check: it is
  among the tests above, and the implementer shows in the report that it fails when the reach is
  kept from the last widening instead of computed for each next record (the review's case,
  `doc/specs/loci.md`, last row of the table of a locus that reaches back).
- The four checks of `python.md` pass.

**What it stands on.** The cherry-pick of `90a3c39`; `VcfRecord` in `src/pvc_eval/vcf.py`, which
is not changed.

**What could go wrong.** The look-ahead is where the spec's two reviews found faults: the reach,
the furthest position a later record can start at and still join the open locus (the spec's "What
`build_loci` holds"), which must be computed after the next record's segments are read; positions compared without
their contig; records with no ALT and STR regions in a gap of the open locus. The reviewers of
the `spec` and `correctness` categories get these three. The spec's model of these rules, which
agreed with a brute force of §5.3 on 20,000 random cases, is
`/private/tmp/claude-501/-Users-jose-devel-pvc-eval/28493737-955f-4f83-acc5-62de2bc6b142/scratchpad/model2.py`
while that directory lasts; it is not part of the code.

- [x] 2. `build_loci` on the spans `parse_vcf` computes

**What it gives.** The loci tested on the spans `parse_vcf` gives, not on spans written by hand,
so that a disagreement between the two stages shows here.

**The spec it builds.** `doc/specs/loci.md`, "How it is verified": the traded events, the second
and third cases of `doc/spec.md` §5.2 run through `parse_vcf`, and "With the STR finder".

**Deliverables.**

- `tests/test_loci_from_vcf.py`, writing VCF text and a FASTA with its `.fai` in `tmp_path`: the
  traded events of `doc/spec.md` §5.3 at chain distance 1 (one locus) and 0 (two), and the second
  and third cases of §5.2 (one locus each), each locus as a literal. Check: `uv run pytest
  tests/test_loci_from_vcf.py -q` ends with 5 passed and 1 skipped: the review of milestone 2 added a
  case only the span joins, since the two of §5.2 share a locus by their REFs alone.
- The test with the STR finder: `segment_contig` over a FASTA holding an AT repeat and a VCF
  inserting AT before it, one locus holding the repeat. It is skipped, with a message naming
  `segment_contig`, while `pvc_eval.segments` has no `segment_contig`. Check: the skipped test of
  the line above, and `uv run pytest tests/test_loci_from_vcf.py -q -rs` names it.
- The four checks of `python.md` pass.

**What it stands on.** Milestone 1; `parse_vcf`, unchanged.

**What could go wrong.** A span `parse_vcf` computes that differs from the one the spec's case
assumes. That is a finding against `vcf.py` or the spec, not against the loci, and a stop: the
owner's handoff forbids changing `vcf.py` without them. The session commits the tests that pass,
leaves the failing case out of the commit with its output in the report, and asks the owner as the
coding skill's "When to stop and ask" says; milestone 3 can go on meanwhile.

- [ ] 3. Kind, stratum and regions

**What it gives.** The kind and the stratum of a locus, and the selection of `--regions`, which
the table of §6 and `compare` will read.

**The spec it builds.** `doc/specs/loci.md`, "What is read from a locus": Regions, Kind, Stratum,
with their interfaces; `doc/spec.md` §5.6 and §5.7.

**Deliverables.**

- `pvc_eval.kind` with `Kind` and `locus_kind`. If the reconciliation session has written `Kind`,
  its commit is cherry-picked and `locus_kind` added; otherwise `Kind` is written first, in a
  commit of its own, and sent to `reconciliation`. `tests/test_kind.py` holds the four rows of the
  spec's table, an ALT equal to its REF left out, and one locus of each kind, one of them holding 14
  unphased heterozygous records, *too complex* for reconciliation, which still has a kind (§7.2).
  Check: `uv run pytest tests/test_kind.py -q` ends with at least 9 passed, the cases this
  deliverable lists.
- `pvc_eval.strata` with `pvc_stratum`, `read_strata` and `Strata`. `tests/test_strata.py` holds
  the strata cases of `doc/spec.md` §7.2 at `pvc_stratum`, the untyped region that raises, the
  spec's `homopolymer,segdup` case and a locus under no interval (`.`), a BED read unsorted and
  with overlapping lines, and each refusal of `read_strata`, the contig the reference lacks among
  them, with its message. Check: `uv run pytest tests/test_strata.py -q` ends with at least 14
  passed, the cases this deliverable lists, each refusal counted once.
- `pvc_eval.regions` with `read_regions`, `Regions` and `select_loci`. `tests/test_regions.py`
  holds the spec's locus kept by a region that shares positions and skipped by one that lies next
  to it, every locus kept when `regions` is None, a BED read unsorted with overlapping and adjacent
  lines, and each refusal of `read_regions`, the contig the reference lacks among them, with its
  message. Check: `uv run pytest tests/test_regions.py -q` ends with at least 9 passed, the cases this
  deliverable lists, each refusal counted once.
- The four checks of `python.md` pass.

**What it stands on.** `pvc_eval.locus` from milestone 1, or from the reconciliation session;
not `build_loci`, since the tests build `Locus` values by hand. So it can run side by side with
milestone 2, in different files.

**What could go wrong.** Messages a user reads: the two refusals the owner decided on 19 October
2026 (a BED contig the reference lacks, an untyped STR region under `--strata pvc`) say what to do,
as the writing skill's section on errors asks. The `errors` reviewer gets them.

## The check of the whole plan

The milestones' tests use small inputs. At the end, `build_loci` is run over the benchmark VCF of
chr20 against itself, with no STR regions: `chr20_5mb_200s.vcf.gz` (18,310 records, 200 samples),
in `.claude/worktrees/perf-vcf-parser/benchmarks/data/`, made by `benchmarks/make_vcf.py` on the
branch `perf/vcf-parser` as its README says, over GRCh38 at the path the README gives. A script in
the session's scratchpad, not committed, reads both streams with `parse_vcf` and prints the number
of loci, the records of each file the loci hold, and the peak memory tracemalloc traced. Check:
each file's records with an ALT, summed over the loci, equal the file's records with an ALT, and
every locus holds the same records in A and in B. The numbers go into the report. If the file is
gone, the report says so and the check is not run.
