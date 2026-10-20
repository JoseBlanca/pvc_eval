# Plan: one trim of an ALT against its REF

*Status: under way, 20 October 2026. A change with no spec, planned as one milestone, on the
branch `shared-trim`; the owner asked for it on 20 October 2026.*

Three places trim an ALT against its REF, removing the bases they share first at the left end and
then at the right, as `doc/spec.md` §5.4 defines an edit: `_get_edit_kind` in
`src/pvc_eval/kind.py`, which gives the kind of §5.6; `_make_edit` in `src/pvc_eval/reconcile.py`,
which builds the edit; and `_find_indel_span` in `src/pvc_eval/vcf.py`, which finds an indel's span
(§5.2). The kind is defined "as an edit is" (§5.6), so a change made to one copy and not to the
others would give kinds or spans that disagree with the edits, with no error. This change makes the
three one function. It changes no result.

## Milestone

- [ ] 1. One function that trims an ALT against its REF

**Deliverables.**

- One function, in `src/pvc_eval/vcf.py` beside the VCF record whose alleles it trims, that gives
  how many bases REF and ALT share at the left, and what is left of each once the shared bases are
  removed at the left and then at the right; `kind.py`, `reconcile.py` and `vcf.py` use it, and
  their hand-written trims and `reconcile.py`'s copy of `_count_shared` are gone. Check: `grep -n
  "_count_shared\|while prefix" src/pvc_eval/*.py` finds only the function and its helpers.
- Tests of the function on its own: an ALT that spells REF, a SNP, an insertion and a deletion at
  either end, a complex ALT, and `52 AG → AAG`, where trimming left first puts the insertion after
  the anchor (§5.4). Check: they fail on `main` at `b71c83c`, where the function does not exist.
- Every existing test passes unchanged, which is the check that no result changed.
- `parse_vcf` no slower: the small input of `benchmarks/README.md`, measured with hyperfine before
  and after, within the noise of the two measurements.
- The four checks of `.claude/skills/coding/python.md` pass.

**What could go wrong.** `_find_indel_span` runs once per ALT of every record, so the function call
must not add measurable time to `parse_vcf`.
