# Plan: pvc's hard cases

*Status: under way, 21 October 2026. It builds `doc/specs/hard_cases.md`, approved by the owner on
21 October 2026 with its 50 expected outcomes, and the owner's answer to that spec's Open 1, a `*`
measured against its record's REF as written, which changed `doc/spec.md` §5.4,
`doc/specs/reconciliation.md` §4.2 and §7, and `doc/specs/respell.md` the same day. The owner
asked for both on the branch `hard-cases`, in the worktree `.claude/worktrees/hard-cases`. The work
report will be `doc/reports/hard-cases.md`.*

The plan makes pvc_eval measure a `*` against its record's REF as written, and then turns the 50
cases of the spec into tests through `compare()`. The rule holds for every record with a `*`; for
one whose alleles share no base at their ends, the REF as written is the normalised REF, so only
padded records, as pvc writes them, read differently. The rule comes first because two of the
cases expect it: S11 and E5, the cases of the spec's §4.5 where pvc pads a record with a `*`. The
spec names its cases by kind, S for GIAB's split records against pvc's repeat record, E for
deletions at a repeat's edge, H for HipSTR and O for the others; R1 to R10 are the cases of
`doc/specs/reconciliation.md` §4.4, R10 the padded `*` record that milestone 1 adds.

## In and out

In: the rule of `doc/spec.md` §5.4 for a `*` and what it changes in `pvc_eval.vcf`,
`pvc_eval.reconcile`, `tests/haplotype_oracle.py` and `tests/respell.py`; §4.5 and §5 of
`doc/specs/hard_cases.md`.

Out: pvc's own errors that the cases show, and the reading of HipSTR's `|` without PS
(`doc/specs/hard_cases.md` §8).

## Before the first milestone

On commit a050b57, on 21 October 2026, CPython 3.15 free-threaded: `uv sync --locked` is clean,
`ruff format --check`, `ruff check` and `ty check` pass, and `uv run pytest -q` gives "1886
passed, 9 skipped in 118.90s" on the Apple M5 Pro (the skipped ones need pvc, GRCh38 or the HG002 VCF). GRCh38 is at
`/Users/jose/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`
and SL4.0 at `/Users/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`, each
with its `.fai`.

## Milestones

They run in order: milestone 2's cases S11 and E5 stand on milestone 1.

- [ ] **1. A `*` measured against its record's REF as written**

  What it gives: a padded record with a `*`, as pvc writes, is read as its caller meant, in
  pvc_eval, in the haplotype oracle and in the re-spelling program.

  Spec: `doc/spec.md` §5.4, the `*` bullet; `doc/specs/reconciliation.md` §4.2, "Impossible", R10
  of §4.4, and `VcfRecord.written_ref` in §7; `doc/specs/respell.md`, the kept record "whose reach
  … touches the REF as written of a record with a `*` ALT", and in §8 a haplotype that carries a
  `*`.

  Deliverables:
  - `VcfRecord.written_ref`: `uv run pytest tests/test_vcf.py` has tests that `parse_vcf` sets it
    to the REF as the file wrote it, for a padded record `10 TATACAC TATAC,*`, whose normalised
    REF is 12–14 and REF as written 10–16 (1-based), and for an unpadded one; that `None` stands
    for the normalised REF; and that a `written_ref` not holding the normalised REF raises
    `ValueError`.
  - R10 at `reconcile`, in `tests/test_reconcile.py`: alleles `TACACATA,TACATA,TACA`, possible
    genotypes `(1, 2)`. It fails on a050b57, which gives invalid, `contradicting records`.
  - The oracle on R10: `check_same_haplotypes` with file A R10's records and file B the same
    haplotypes written plainly, `12 TAC T 0|1` and `15 ACATA A 1|0`, passes; on a050b57 it raises
    `AssertionError`.
  - The re-spelling program: a record that touches the REF as written of a padded `*` record and
    not its normalised REF is kept whole, a test in `tests/test_respell.py` that fails on a050b57.
  - Every test that passed on a050b57 passes unchanged: R1 and R2, whose `*` records are not
    padded, and the property and re-spelling checks among them.

  Stands on: nothing before it.

  What could go wrong: `VcfRecord` is built by hand in many tests, so the field has a default; a
  test that builds a padded record by hand without it gets the old rule silently. The oracle and
  `reconcile` must read the same field, or the property tests of reconciliation compare two rules.

- [ ] **2. The 50 hard cases through `compare()`**

  What it gives: §7.5 of `doc/spec.md` checked, case by case.

  Spec: `doc/specs/hard_cases.md` §1, §4.5 and §5.

  Deliverables:
  - `tests/hard_cases.py`: one literal per case, made from §4.5 of the spec, which is the source of
    every value: the records, written in excerpt positions with the line of §4.5 beside each as a
    comment; the sample; the locus, kind and STR regions; the locus alleles; the possible
    genotypes of each file, as sequences made from §4.5's indices and alleles; and the reported
    genotypes, counts and reasons. The excerpt's 301 bases are read once from GRCh38 or SL4.0 at
    the stretch §4.5 gives, by a script in the session's scratchpad that is not committed; the
    module's docstring says how they were read, and the two excerpt tests below check them against
    the references.
  - `find_possible_haplotypes` in `tests/haplotype_oracle.py`, as §5 of the spec gives it, with a
    test on R10 and on S4's file A, whose two phased records lie in two windows of one phase set:
    taken together they give one possible genotype, where windows taken one by one and joined
    would give a second.
  - `tests/test_hard_cases.py`, as §5 of the spec gives it: `uv run pytest tests/test_hard_cases.py -q`, without `PVC_EVAL_GRCH38` and `PVC_EVAL_TOMATO`,
    gives "50 passed, 2 skipped"; with both set to the paths above, "52 passed". Each of the 50
    is named by its case in the report.
  - On a050b57 `tests/test_hard_cases.py` does not exist, and pytest exits with 4. Once it does,
    48 of the 50 cases would pass on the code of a050b57, which already gives their rows; so each
    of the test's two checks is shown failing on a working copy that is not committed, as §5 of
    the spec says. The oracle's check: an ALT of file B changed by one base in S1. The check of the
    rows `compare()` gives: with the oracle's literal unchanged, the expected outcome of S3
    changed, and one expected locus allele of H1. The report gives the failure line of each.

  Stands on: milestone 1.

  What could go wrong: a literal copied from what pvc_eval prints rather than from §4.5 checks
  nothing. After milestone 1 the two should be equal in all 50 cases. A value of `compare()` that
  differs from §4.5 is a stop for the owner, recorded in the report with the case, the expected
  row and the row given, and the literal is not changed to match.

## When the plan is done

`uv run pytest -q` passes with the 52 new tests among the rest, and with `PVC_EVAL_GRCH38` and
`PVC_EVAL_TOMATO` set, the two excerpt tests pass too. `doc/spec.md` §7.5 already points to the
spec.
