# Plan: a REF letter outside ASCII

*Status: under way, 21 October 2026, on the branch `hard-cases`, which the owner asked for. A
small change with no spec of its own; it builds the sentence of `doc/spec.md` §5.2 that says a REF
letter outside ASCII is no base. The work report is `doc/reports/hard-cases.md`.*

`parse_vcf` compares a VCF record's REF with the reference after upper-casing it with Python's
`str.upper`, which upper-cases every letter of Unicode. A REF `ß` at a reference `G` became `SS`,
one base longer than the reference's bases, and the run stopped with a `ValueError` traceback from
`zip(..., strict=True)` in `_check_ref` of `pvc_eval.vcf`, with no file or line. A REF `ſ` became
`S`, the IUPAC code for C or G, and was accepted at a `G`. Both came from the errors review of
milestone 1 of `doc/plans/hard-cases.md`; both are on `main` too. The fix upper-cases only ASCII
letters, so that such a REF differs from the reference and stops the run with the error §5.2 gives
for a REF that differs, naming the file, the line and the position.

Before it: at 95f100b the checks pass and `uv run pytest -q` gives "1968 passed, 11 skipped".

- [ ] **1. A REF letter outside ASCII differs from the reference**

  Spec: `doc/spec.md` §5.2, the sentence on letters outside ASCII.

  Deliverables:
  - In `tests/test_vcf.py`, `parse_vcf` raises `PvcEvalInputError` naming the file, the line and
    the position for a REF `ß` at a reference `G` (a `ValueError` at 95f100b) and for a REF `ſ` at
    a `G` (accepted at 95f100b); and for one in a longer REF, `AßT`. Each fails at 95f100b.
  - A lower-case ASCII REF still agrees with its upper-case reference, and an IUPAC code still
    agrees as before: the tests of §5.2 that pass at 95f100b pass unchanged.
  - What an ALT holding a letter outside ASCII does is found and said in the report, with no change
    to it unless it also gives a traceback.

  What could go wrong: the alleles are upper-cased in another place too, when the record is built;
  a REF must be upper-cased the same way in both, or a REF that passes the check is stored
  otherwise.
