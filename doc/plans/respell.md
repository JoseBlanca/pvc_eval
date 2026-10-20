# Plan: the re-spelling program

*Status: approved by the owner on 20 October 2026; under way. It builds `doc/specs/respell.md`,
approved by the owner on 20 October 2026, on the branch `respell` in the worktree
`.claude/worktrees/respell` (the name the owner gave the branch, in place of `plan/respell`). The
work report will be `doc/reports/respell.md`.*

The plan builds the four re-spellings of `doc/spec.md` §7.4, the oracle that checks them, and the
check of §7.4 itself: a file against each of its re-spellings, through the stages of pvc_eval.
Decomposition comes first, with everything the four share, since it is the simplest and the one
the reconciliation tests already use. Each re-spelling is checked by the oracle before it is
checked through the stages, so that a failure of §7.4 is known to be pvc_eval's and not the
program's.

## In and out

In: §2 to §10 of the spec, in the six milestones below. All of it under `tests/`.

Out, with where it goes:

- The checks over a GIAB HG002 chromosome (spec §9.4) are written in milestone 3 and run in each
  later one, but skipped until the owner says where the VCF and its reference are; the session
  building `compare()` asks the same.
- `compare()` is not on `main`; the check of §7.4 joins the stages itself (spec §9.3). Moving its
  rows to `compare()` is a change for when it arrives, and with it the check takes the reference
  as a path, as `compare()` does, in place of an opened `IndexedFasta` (review of milestone 3).
- Chaining re-spellings, and anything in `src/` or on the command line (spec §12).

The spec has no open point: the owner answered Open 1 on 20 October 2026, keep the rule, which
milestone 5 builds.

## Before the first milestone

On commit 9dbec00, on 20 October 2026, with CPython 3.15 free-threaded: `uv sync --locked` is
clean, `ruff format --check`, `ruff check` and `ty check` pass, and `uv run pytest -q` gives
"914 passed, 4 skipped" (the skipped ones need pvc and GRCh38). bcftools 1.24 is at
`/opt/homebrew/bin/bcftools`; the tests that run it are skipped without it. The cases of the
spec's §3 to §6 were run through `reconcile` on 20 October 2026, input against output, and gave
equal sample results.

## Milestones

The milestones run in order: 4, 5 and 6 each add to `tests/respell.py`, and 2 to 6 stand on
what 1 builds.

- [x] **1. Decomposition, and what the four re-spellings share**

  What it gives: `respell_by_decomposing`, and under it what every re-spelling uses: reading the
  input as lines and as records, the records kept, the touched test with the phase, the reach and
  its slide, the writing of rewritten records, new phase sets, the header, the order and the
  BGZF writer. The reconciliation property uses it, and `_make_decomposed_lines` is removed.

  Spec: §2, §3, the first block and the last paragraph of §7, §10.

  Deliverables:
  - `uv run pytest tests/test_respell.py -v` runs D1 to D12, each a test named for its case with
    its output lines as literals; for each row of §2.2, a test that the record is written as read
    beside a decomposed neighbour, and for a REF mismatch `PvcEvalInputError` with no output file;
    the PS header line added once, and not when present or when an input PS is not a whole
    number; a bgzipped input giving a bgzipped output that `parse_vcf` reads; and an
    output named as its input raising `ValueError`. At least 25 tests, all pass. On 9dbec00 the
    file does not exist (pytest exits with 4).
  - `uv run pytest tests/test_reconcile_property.py -v` passes, with file B made by
    `respell_by_decomposing`, and `grep -rn _make_decomposed_lines tests` finds nothing.
  - `respell` is in the verb table of `.claude/skills/coding/python.md`.
  - The touched test sweeps records by REF: on a file of 20,000 single-sample SNPs, ten bases
    apart, `respell_by_decomposing` takes under 10 s on the owner's Mac; the test that times it is
    marked `slow`. Spec §10 puts the square of 20,000 records at some 96 s.

  Stands on: nothing.

  What could go wrong: the touched test is where a wrong rule is silent, since a record wrongly
  rewritten gives a file B that only the oracle and §7.4 can catch, in milestones 2 and 3. The
  cases D4, D5, D6 and D12 guard it here, and the reviewers should look for a case the phase rule
  of "may share a haplotype" gets wrong. The reach (§2.2, third row) needs a slide of the tests'
  own; it must agree with `parse_vcf`'s spans where those are not cut, which a test checks on
  S1, S2 and the cases of `doc/spec.md` §5.2.

- [x] **2. The oracle and the generated files**

  What it gives: `check_same_haplotypes`, which shows a re-spelling right without pvc_eval's
  reconciliation, and the hypothesis generator of small contigs with repeats; both run on the
  decomposition.

  Spec: §8, the second block of §7, §9.2, and the `check_same_haplotypes` and `bcftools` items
  of §9.1.

  Deliverables:
  - `uv run pytest tests/test_haplotype_oracle.py -v` runs O1 to O6, each named for its case, with
    the window and both sets of a raising case asserted in its message, and the counts of O5; a
    test that `tests/haplotype_oracle.py` imports none of `pvc_eval.reconcile`, `pvc_eval.loci`
    and `pvc_eval.evaluate`; a window over `ORACLE_CAP`, counted and not compared. All pass. On
    9dbec00 the file does not exist (pytest exits with 4).
  - In `tests/test_respell.py`: every D case passes the oracle, input against output; bcftools
    reads every D output as §2.3 says; the oracle property over 300 decomposed files and the
    bcftools property over 50 pass; at least one file in ten has a rewritten record.

  Stands on: 1.

  What could go wrong: an oracle that agrees with a wrong re-spelling. O1, O2, O3 and O6 are the
  wrong re-spellings it must reject; the reviewers should try more, such as the pieces of a
  phased record put in a new phase set, and check that the oracle raises.

- [x] **3. The check of §7.4, for decomposition**

  What it gives: `check_respelling_outcomes`, a file against its re-spelling through `parse_vcf`,
  `build_loci`, `reconcile_loci` and `make_outcome_rows`, and the check of §7.4 for the
  decomposition.

  Spec: §9.3, the third block of §7, the `check_respelling_outcomes` item of §9.1, §9.4.

  Deliverables:
  - `uv run pytest tests/test_respell_check.py -v` runs: every D case through the check; the AT
    repeat case through `segment_contig`, with its one locus holding the STR region and *same
    genotype*; a hand-made B that loses a phase assignment and one with a wrong sequence, each
    raising; the §7.4 property over 300 decomposed files, with the planted repeats as a
    `--str-regions` BED; and the GIAB test, skipped with its message. All pass but the skipped
    one. On 9dbec00 the file does not exist (pytest exits with 4).

  Stands on: 1 and 2.

  What could go wrong: a property that passes because the generator makes few loci where records
  interact. The report gives, from a run outside the suite of 5,000 examples, how many rows had
  each outcome.

- [x] **4. Shift**

  What it gives: `respell_by_shifting`, checked by the oracle and then through the stages.

  Spec: §4, and its parts of §7 and §9.

  Deliverables:
  - In `tests/test_respell.py`: S1 to S9, each named for its case, with literals; the same seed
    giving the same bytes; the §2.2 rows for the shift; every S case through the oracle; the
    oracle property over 300 files and the bcftools property over 50; at least one file in ten
    with a shifted record. On 9dbec00 none exists.
  - In `tests/test_respell_check.py`: every S case, the §7.4 property over 300 shifted files, and
    the AT repeat case of milestone 3 with its file B made by `respell_by_shifting` in place of
    the one milestone 3 writes by hand.

  Stands on: 3.

  What could go wrong: the slide past the contig's end and the span cut by the span limit, which
  S7, S8 and the span condition of §4 guard.

- [x] **5. Spanning deletions**

  What it gives: `respell_with_spanning_deletions`, checked the same way.

  Spec: §5, and its parts of §7 and §9.

  Deliverables:
  - In `tests/test_respell.py`: T1 to T8 with literals; the case of the trimmed REF (`12 TC TA`
    under a deletion of 12, which gets no `*`); the §2.2 rows; the oracle and bcftools properties;
    at least one file in ten with a `*` written, for which the generator draws deletions over
    other records.
  - In `tests/test_respell_check.py`: every T case, and the §7.4 property.

  Stands on: 3.

  What could go wrong: condition 2 of §5, which rewrites the deletion's phase too: T3 guards it,
  and O3 shows the oracle catches a `*` on the wrong haplotype.

- [ ] **6. Merge**

  What it gives: `respell_by_merging`, checked the same way. It ends the program.

  Spec: §6, and its parts of §7 and §9.

  Deliverables:
  - In `tests/test_respell.py`: M1 to M12 with literals; a run all REF, not merged; the same seed
    giving the same bytes; the §2.2 rows; the oracle and bcftools properties; at least one file
    in ten with a merged record.
  - In `tests/test_respell_check.py`: every M case, and the §7.4 property.

  Stands on: 3.

  What could go wrong: a merge that invents a phase, which M3 and O1 guard, and one whose merged
  stretch takes in another record's edit, which M12 guards.

## At the end

The five checks of the `coding` skill pass, and the eight properties, the oracle's and §7.4's for
each of the four re-spellings, pass at 300 examples in the suite; the report gives a run of 5,000 examples of each
outside it, with any failure it found written into the spec as a case. The checks over GIAB run
when the owner says where its VCF is.
