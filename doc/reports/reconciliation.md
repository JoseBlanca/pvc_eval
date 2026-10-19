# Work report: reconciliation and evaluation

*Under way, started 19 October 2026. The plan is `doc/plans/reconciliation.md`, built from
`doc/specs/reconciliation.md`, on the branch `reconciliation`.*

The checks on the commit the work starts from, 76d14dc: `uv sync --locked` clean, `ruff format
--check` and `ruff check` pass, `ty check` passes, `uv run pytest -q` gives "213 passed".

## Milestone 1: ploidy in the VCF reader

Done in 35654af and fixed in 2f0cbc3 to 8f60a5a. `parse_vcf` takes `ploidy`, 2 by default; a GT
of another number of alleles is `OTHER_PLOIDY_GT`, a GT of only `.` and a sample with no GT are
`MISSING` ploidy times.

- The seven cases of the spec's `parse_vcf` item, and two more from the review (a sample with no
  GT at ploidy 1, a GT at ploidy 3): `uv run pytest tests/test_vcf.py -k ploidy` gives "11
  passed". The whole suite gives "225 passed", 213 before the milestone.
- The five checks pass on 8f60a5a.

Changed in the plan: the plan let a test change only where it names `PLOIDY` or `HAPLOID_GT`.
Four tests that pass `span_limit` through a dict also changed, in their type only, because ty
read the dict as possibly filling the new `ploidy`; the session allowed it, since nothing they
check changed.

The review found no blocker. Fixed: two missing tests that a reviewer showed could not fail
(a sample with no GT read with ploidy 1, found by two reviewers; a ploidy of 3), a test that
could not tell an error at the call from one when the stream starts, and three nits of
docstrings and of a branch that repeated the GT parser. Not taken: two styles of whole-number
check in `parse_vcf`, a nit not worth a change in a file another branch is editing; and a note
that `doc/spec.md` §5.4 still describes a haploid GT, which is not so on this branch (76d14dc).

For the owner: the reader gives `MALFORMED_GT` for a GT with a phase mark before its first
allele, `|1`, which VCF 4.4 allows for a phased haploid call. pvc_eval reads VCF 4.2 and 4.3,
where it is not allowed; it is later work if 4.4 files are to be compared. A sample with no GT
in `vcf.py` reached the GT parser through its own branch, now through the GT parser itself.

## Changes to the plan before milestone 3

`pvc_eval.segments` was written first by the STR finder session, in 90a3c39 on its branch, as
`doc/spec.md` §4 gives it. It is picked into this branch as 6bbc7bc, and milestone 3 writes only
`pvc_eval.locus`.

## Milestone 2: sample pairs

Done in adc88a9 and fixed in f5c0cff to 085d259, built on the branch `reconciliation-samples` so
that it could run beside milestone 1, and merged here. `make_sample_pairs` pairs the samples of
two files by name or by a `--sample-map`, with the errors and the warning of the spec's §3.

- `uv run pytest tests/test_samples.py` gives "32 passed"; the plan asked for at least 10. The
  whole suite gave "257 passed" on 085d259.
- The five checks pass on 085d259.

Changed in the spec: the implementer found six small cases §3 left open, and the review six
more; the session wrote them into §3 before the code (90f493a, effc30d). Among them: a map that
is missing or not UTF-8 stops the run with a message, not a traceback; a byte-order mark is
ignored; a line of only spaces and tabs is blank; a name missing from its file is quoted, so a
space in it shows.

The review found no blocker. Fixed: the unreadable map; the byte-order mark, which the fix strips
itself, since decoding with `utf-8-sig` put a bad byte on the wrong line; three behaviours no test
could fail on (the cap of 10 names at exactly 10, an empty first field, a last line with no
newline); pairing by map took time quadratic in the samples, 5.28 s for a map of 30,000 lines,
now 0.035 s; and the count of a noun in a message, `format_count`, is now one function in
`errors.py` that the VCF reader and the sample pairs share. Not taken: a lone carriage return
inside a map line splits the line; Windows line ends read correctly, and nothing else writes it.

For the owner: a VCF that does not exist stops the run with a traceback in `vcf.py`, which the
sample pairs meet first since they read the headers; and a `#CHROM` line that ends in a tab gives
a sample named "". Both are older than this plan and belong to the VCF reader.

Changed in the plan: milestone 3 stands on milestone 2 too, since `reconcile` takes the `SamplePair`
of milestone 2; the plan first named milestone 1 alone, and milestone 3's implementer found it.

## Milestone 3: reconciliation

Done in f603bec (`Locus`), e2a8ac7 (`InvalidReason`, `LocusGenotype`) and f6822e4 (`reconcile`,
`reconcile_loci`), and fixed in 50d1293 to af717e6. For a `Locus`, the locus alleles and each
paired sample's possible genotypes, *too complex*, or *invalid* with its reason.

- `uv run pytest tests/test_reconcile.py` gives "68 passed": every case R1 to R9 of the spec's
  §4.4 with its label in the test's id, the other cases of the `reconcile` item of §9, and the
  cases the review added. `tests/test_locus.py` gives "6 passed". The whole suite gives "337
  passed" on af717e6.
- The five checks pass on af717e6.
- `Locus` (f603bec, with its checks in 18ce521) was sent to the loci session through the
  coordinating session.

Changed in the plan: milestone 3 stands on milestone 2 as well as milestone 1 (above).

The review found one wrong outcome, a blocker found by two reviewers apart: an insertion whose
point lies inside a `*` record's REF counted as covering the `*`, so a sample the spec makes
invalid got a possible genotype. It is fixed, with two tests that fail without the fix. A
reviewer compared the code with a prototype of the rules, corrected to the spec, over 20,000
random loci of up to four records, phased and not, with cut spans among them: with the fix,
every outcome, set of locus alleles and order agreed. The same code gives equal results run
serially and on 8 threads.

The review also found nine plausible wrong rules that no test could fail on, among them `./1`
making a block, the order of the assignments with three blocks, the insertion collision in the
other record order, and two orders of the invalid reasons; each now has a test that fails
under the wrong rule, confirmed by the reviewers' second looks. Further: `Locus` checks its
own invariants and `reconcile` a sample index, both added to the spec's §7 first (10ee2e2);
each edit is made once per locus, 13 times at 13 blocks where it was made 106,496 times; and
the helpers take the verbs of `python.md`. Not taken: a third copy of the whole-number check
(in `vcf.py`, which another branch edits); `kw_only` on `Locus`, which `doc/spec.md` §4 fixes
without it; a read-only mapping for the sample results, which nothing writes; lower-case bases
passed to `reconcile`, which `reconcile_loci` never passes.

For the owner: two cases pass unchecked, both only from a `Locus` built by hand. A sample index
past the genotypes goes unseen at a locus where its file has no record. And a record with no
ALT that does not overlap the locus at all is accepted; with a REF mismatch it would make its
file's samples invalid at a locus it does not touch.

## Milestone 4: evaluation and the row

Done in 4df38de (`Kind`), 3c17c07 and 1c27a9b (`Outcome`, `OutcomeRow`, `evaluate_pair`,
`make_outcome_rows`), and fixed in 24c6558 to f6b6903. Each sample pair's outcome, its reported
genotypes and its row, every field filled.

- `uv run pytest tests/test_evaluate.py` gives "37 passed": E1 to E14 of the spec's §5 with their
  labels as ids, the other side of each rule, the reported pair where A's order and B's order
  differ, and every field of the row. The whole suite gives "374 passed" on f6b6903.
- The five checks pass on f6b6903.
- `Kind` (4df38de), exactly as the loci spec gives it, was sent to the loci session through the
  coordinating session.

Changed in the spec and the plan: E14 joined §5 (149a54b), a case where A's genotypes taken first
and B's taken first report different pairs; the plan checks E1 to E14 (f4eccdb).

The review found no wrong outcome. A reviewer compared `evaluate_pair` with its own reading of §5
over 4,225 combinations of sample results, both ploidies, every outcome reached, with no
difference; another ran seven kinds of bad input through reconciliation and the rows and got the
row §6 gives for each. It found one blocker: no case told the order of "the first such pair"
from its reverse, so swapping two loops passed every test; E14 and two more cases now fail under
it, confirmed by the reviewer's second look. It also found three rules no test could fail on
(the first genotype reported for *different*, *only in B* and *missing in A*; the "holds one
genotype" condition of *only in*; the QUAL and GQ of an invalid sample), each now tested, and a
loop four levels deep, now one helper called once per rule. Not taken: checks in
`OutcomeRow.__post_init__`, which depend on what the other evaluators' wrappers will fill; a
check of its allele indices in `ReconciledLocus`, which only `reconcile` builds; and making
`Kind` a `StrEnum`, which the loci spec owns. That last one was sent to the loci session as a
suggestion: as a plain `Enum`, `Kind` formats as `Kind.SNP` where the other closed sets of a row
format as their words, so the table writer must use `.value` or the loci spec change it.

For the owner: across ploidies, only genotypes all REF count as equal, as you decided. So a
haploid `1` against a diploid `0/1`, like a haploid `1` against `1/1`, is *same alleles, other
genotype*, by the ALT sets alone.

## Milestone 5: stopped for the owner

The decomposition the property needs is built and checked on four hand-worked cases (7ca0a93).
The property itself, over 300 generated cases, found that the collision rules of the spec's §4.2
give a file and its decomposition different outcomes when two records of one haplotype touch:

- reference `TAAC`, `3 AC G,AA 0/1` and `3 A G,AA 1/2`: the insertion of an A after 3 falls
  between the two positions `AC → G` replaces, so A's sample is invalid; decomposed, `AC → G`
  becomes a SNP at 3 and a deletion of 4, the insertion lies between them, and B's is valid;
- reference `GTTT`, `2 TT GCTA 1/1` and `3 T G,TC 2/2`: the bases `TT → GCTA` adds become, decomposed,
  an insertion at the end of the stretch, at the same point as the other record's insertion, so
  B's sample is invalid and A's valid.

Where the records of A do not touch, no row failed in 6,290 generated cases. The property is not
committed; it waits on the owner's answer.
