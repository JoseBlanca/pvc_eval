# Plan: a damaged gzip reported as damaged, and an ALT `.` among others refused

*Status: done, 21 October 2026, on the branch `vcf-input-errors` in the worktree
`.claude/worktrees/vcf-input-errors`, from `main` at bd33d55. The owner ordered both changes on
21 October 2026, after the work of the re-spelling program found them; the rules are in
`doc/spec.md` §3 and §5.2 (4a9acae). The work report will be `doc/reports/vcf-input-errors.md`.*

Two faults of a VCF that pvc_eval reports wrongly:

- A gzipped VCF with a damaged byte is reported as a fault of whichever line its garbled text
  makes first, "not sorted" or "the line has 5 columns", since gzip checks its data only at the
  end of each compressed block: of 5,087 copies of one gzipped VCF each with one byte changed,
  4,202 were reported so, 884 as damaged, and 1 not at all. The fix: when reading a gzipped VCF
  or BED stops at a fault of its content, the file is checked whole, and if gzip finds it
  damaged, the message says so.
- An ALT `.` among other ALTs (`G A,.`) is read as an allele. It must stop the run with an error
  naming the file and the line.

## The milestone

- [x] **1. Both faults reported as `doc/spec.md` §5.2 says**

  Spec: `doc/spec.md` §3, "Errors, which stop the run", and §5.2, "Reading" and "The ALT `.`".

  Deliverables:
  - `uv run pytest tests/test_vcf.py -k "alt_dot or damaged" -v` runs at least these, and all
    pass: `G A,.` and `G .,A` each raise `PvcEvalInputError` naming the file and the line;
    `G .` stays a record with no ALT; a gzipped VCF damaged after its header, so that it reads
    as lines out of order or of the wrong number of columns, is reported as damaged, through
    `parse_vcf` and through `parse_vcf_samples` where the damage is in the header; a plain VCF
    with the same bad line keeps its own message. On bd33d55 they fail or do not exist.
  - The same damaged-file test for a gzipped BED read by `pvc_eval.bed`.
  - A property, or a test over every byte of one small gzipped VCF changed in turn: every
    copy that stops the read is reported as damaged, none as a fault of a line. The one copy
    in 5,087 that read without an error on bd33d55 is explained in the report.
  - The five checks of `.claude/skills/coding/python.md` pass.

  Stands on: nothing.

  What could go wrong: the whole check must not run on a good file, nor change a message for a
  plain file; and a damaged file found while the reader's threads parse sample columns must
  still give the damaged message (`doc/spec.md` §5.2, "Threads").

## At the end

The report gives, for the same 5,087 damaged copies, how each is reported after the change.
