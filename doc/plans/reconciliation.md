# Plan: reconciliation and evaluation

*Status: approved by the owner on 19 October 2026; under way. It builds `doc/specs/reconciliation.md`,
approved by the owner on 19 October 2026, on the branch `reconciliation` in the worktree
`.claude/worktrees/reconciliation`. The work report will be `doc/reports/reconciliation.md`.*

The plan builds stages 4 and 5 of `doc/spec.md` §4: from a `Locus`, the locus alleles and each
paired sample's possible genotypes, then each sample pair's outcome and its row, `OutcomeRow`. It
builds what they need first, the ploidy of the VCF reader and the sample pairs, and ends with the
property that a file against a re-spelling of itself gives only the outcomes a file against
itself allows.

## In and out

In: the five milestones below, which build §2 to §7 and §9 of the spec.

Out, with where it goes:

- Writing `OutcomeRow` as a TSV line and the summary, `compare()` and the command line, with its
  `--ploidy-a`, `--ploidy-b` and `--sample-map`: no spec has them yet (`doc/specs/reconciliation.md`
  §12).
- The kind, the stratum and `build_loci`: the loci's plan. Milestones 3 and 4 write only the
  shared types the loci spec fixes, `Locus`, `GenomeSegment` and `Kind`, and the tests build
  `Locus` values by hand.
- The re-spellings of `doc/spec.md` §7.4 that move records, and the file against itself over a
  real region: they need `build_loci`, and come with the plan that joins the stages.

The spec has no open point.

## Before the first milestone

On commit 76d14dc, on 19 October 2026, with CPython 3.15 free-threaded (`.python-version` holds
`3.15+freethreaded`): `uv sync --locked` is clean, `ruff format --check` and `ruff check` pass,
`ty check` passes, and `uv run pytest -q` gives "213 passed". bcftools 1.24 is at
`/opt/homebrew/bin/bcftools`, needed by no test; the outputs it gave are literals in the spec.

## Milestones

Milestones 1 and 2 touch different files and can run side by side. The others run in order.

- [x] **1. Ploidy in the VCF reader**

  What it gives: `parse_vcf` reads a file with a ploidy, 1 or more, 2 by default; a GT with
  another number of alleles is `OTHER_PLOIDY_GT`, and a GT of only `.` is `MISSING` on every
  haplotype at any ploidy. Reconciliation's tests read their records with it.

  Spec: `doc/specs/reconciliation.md` §2, and the `parse_vcf` item of §9.

  Deliverables:
  - `uv run pytest tests/test_vcf.py -v` runs the seven cases of the `parse_vcf` item of §9, each a
    test named for its case, and all pass. On 76d14dc they fail: `parse_vcf` takes no `ploidy`.
  - The 213 tests of 76d14dc still pass. A test may change only where it names `PLOIDY` or
    `HAPLOID_GT`, which this milestone renames, and the commit message lists each test changed.

  Stands on: nothing.

  What could go wrong: the branch `perf/vcf-parser` is measuring `parse_vcf` and may change
  `vcf.py`. The milestone touches only the GT parsing and the constants, to keep a merge small,
  and the report says which lines changed.

- [x] **2. Sample pairs**

  What it gives: `make_sample_pairs` pairs the samples of two files by name or by a
  `--sample-map` file, with the errors and the warning of the spec.

  Spec: §3 and the `SamplePair` and `make_sample_pairs` block of §7.

  Deliverables:
  - `uv run pytest tests/test_samples.py -v`, over headers written to `tmp_path`, runs the cases of
    the `make_sample_pairs` item of §9: pairing by name in A's order; a map's order, a sample in
    two pairs, a blank line; each of the five errors of §3 with its message; the warning, with the
    first 10 names and the count after them. At least 10 tests, all pass. On 76d14dc the file
    does not exist (pytest exits with 4).

  Stands on: nothing.

- [ ] **3. Reconciliation**

  What it gives: `reconcile` and `reconcile_loci`: for a `Locus`, its locus alleles and the sample
  result of each paired sample, possible genotypes, *too complex* or *invalid* with its reason.

  Spec: §1, §4 and the `pvc_eval.reconcile` block of §7, with `InvalidReason` and `LocusGenotype`
  of `pvc_eval.outcome`.

  Deliverables:
  - First, in a commit of its own, `pvc_eval.locus` with `Locus`, exactly as `doc/spec.md` §4
    gives it; the session sends the commit, with `SendMessage`, to the session building the loci
    (branch `loci`), found with `ListAgents`, which writes `build_loci` beside it. `uv run ty
    check` passes. `pvc_eval.segments`, which `Locus` needs, is the STR finder's commit 90a3c39,
    picked into this branch as 6bbc7bc (changed on 19 October 2026: the STR finder wrote it
    first).
  - `uv run pytest tests/test_reconcile.py -v` runs a test for each row and each case of §4.4,
    R1 to R9, each with its locus alleles and possible genotypes as literals and the case's label
    in its id, and the other cases of the `reconcile` item of §9: the order of the locus alleles
    across two samples and two files; a sample in no pair; a file with no record in the locus; a
    record with no ALT whose REF reaches past the locus, `0/0` and `./.`; ploidy 1 with `1` and
    `0`; and the `ValueError`s of §7. All pass. On 76d14dc the file does not exist.
  - `reconcile_loci` over two loci of one contig gives the same results as `reconcile` on each,
    with the bases read from a FASTA written to `tmp_path`.
  - The verb `reconcile` in the verb table of `.claude/skills/coding/python.md`, with the meaning
    §7 gives it, in the commit of the code that first uses it.

  Stands on: milestone 1, whose `OTHER_PLOIDY_GT` the tests of the invalid reason `ploidy mismatch` read,
  and milestone 2, whose `SamplePair` `reconcile` takes (corrected on 19 October 2026: the plan
  first named milestone 1 alone).

  What could go wrong, and what guards it: a wrong rule here gives a plausible wrong genotype, and
  no later check sees it. The rules most likely to be built wrong are the `*` coverage of §4.2,
  with its placements bounded by a cut span (R2), its missing deletion (R1, last row) and its MNP
  (R1); the collision of an insertion between replaced bases, against one at their end (R4); and
  the order of the assignments, which sets the order of the alleles (R7, the interleaved phase
  sets). Each is a case of §4.4, written as a literal. The reviewers are pointed at these.

- [ ] **4. Evaluation and the row**

  What it gives: `evaluate_pair` and `make_outcome_rows`: each sample pair's outcome, its reported
  genotypes and its `OutcomeRow`, every field filled.

  Spec: §5, §6 and the `pvc_eval.outcome` and `pvc_eval.evaluate` blocks of §7.

  Deliverables:
  - First, in a commit of its own, `pvc_eval.kind` with `Kind` only, exactly as the loci spec
    gives it (`git show loci:doc/specs/loci.md`, Kind); the session sends the commit to the loci
    session, as in milestone 3, which writes `locus_kind` there.
  - `uv run pytest tests/test_evaluate.py -v` runs E1 to E13 of §5 at `evaluate_pair`, one
    parametrized test with the labels as ids, and the cases of the `make_outcome_rows` item of §9:
    every field on one locus with two pairs that share a sample; the partner's genotype in a row
    *invalid in A*; the reason in a row *too complex*; `quality` and `gq` with a QUAL absent and
    with a `*`. All pass. On 76d14dc the file does not exist.
  - The verb `evaluate` in the verb table, as in milestone 3.

  Stands on: milestone 3.

  What could go wrong: the reported pair in rows with several possible genotypes (E6, E7, E10) and
  the equality of genotypes all REF across ploidies (E11); a fall-through order that tests E12 and
  E13 only by accident. The reviewers are pointed at these.

- [ ] **5. A file against its decomposition**

  What it gives: the property of §9: for generated loci, file B made from file A by decomposing
  each record into records of one replaced base or one indel, phased in a phase set of their own,
  gives only *same genotype*, *reference in both*, *missing in both*, *invalid in both* or *too
  complex*. Each of the last three comes only from its own cause in A: *missing in both* where A
  has a missing haplotype, *invalid in both* where A's sample is invalid, *too complex* where it is
  over the cap.

  Spec: the property item of §9.

  Deliverables:
  - The decomposition, a function of the tests, with its own cases: a complex record, a
    multi-ALT record with `1/2`, a `./1`.
  - `uv run pytest tests/test_reconcile_property.py -v` runs the property with hypothesis, at
    least 200 examples, through `parse_vcf`, `reconcile` and `make_outcome_rows`, and passes. On
    76d14dc the file does not exist. A failing example hypothesis finds is added as a fixed case
    to `tests/test_reconcile.py`.

  Stands on: milestone 4.

  What could go wrong: the decomposition can itself make records that collide or a `*` that no
  deletion covers, and then the property fails on the test's input, not on the code. Its cases
  are checked by hand first. When a failure shows that a rule of the spec gives a file and its
  decomposition different outcomes, the session stops and asks the owner, as the coding skill's
  "When to stop and ask" says, since the fix changes an outcome; a failure of the code against
  the spec is fixed like any finding.

## At the end

The five checks of `.claude/skills/coding/python.md` pass, with the test count above 213 by the
tests of the five milestones, and the report gives the count.
