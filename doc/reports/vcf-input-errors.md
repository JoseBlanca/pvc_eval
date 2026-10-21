# Work report: a damaged gzip reported as damaged, and an ALT `.` among others refused

*Done, 21 October 2026. The plan is `doc/plans/vcf-input-errors.md`, ordered by the owner on 21
October 2026; the rules are in `doc/spec.md` §3 and §5.2. Branch `vcf-input-errors`, from `main`
at bd33d55, last commit 8e5665a. Nothing is merged into `main` or pushed.*

## For the owner

Both changes you ordered are built, reviewed and fixed, and the branch merges into `main` with
no conflict, since `main` has not moved. The five checks pass on 8e5665a: `uv sync --locked`,
`ruff format --check`, `ruff check`, `ty check`, and `pytest` with "1415 passed, 5 skipped"
(1372 on bd33d55). What is asked of you is the merge.

What a user of pvc_eval now sees:

- **A damaged gzipped VCF or BED is reported as damaged**, where its garbled text used to be
  blamed on a line. Of 5,087 copies of one gzipped VCF, each with one byte changed, 4,202 were
  reported as a fault of a line ("not sorted", "the line has 5 columns"); all 4,202 now read:
  "b.vcf.gz is damaged: gzip's own check of the whole file fails (`gzip -t b.vcf.gz` shows it),
  so the fault pvc_eval found at line 1998 may come from the damaged bytes rather than from the
  file's content. Copy or compress the file again from its source, and run again." 884 keep
  gzip's own "cut short or damaged at or after line N", and 1 still reads, since its changed bit
  is padding after the end of the compressed data that decompresses to the same bytes. This
  holds also where the fault is found after the line is read: sample pairing from a garbled
  `#CHROM` line, a `--strata` class, an end past the contig in `str-regions --regions`.
- **A BGZF file cut between two of its blocks now stops the run.** It read as a valid shorter
  file: on one test, 3,418 of 9,990 records, exit status 0, while bcftools refuses it ("no BGZF
  EOF marker"). A BGZF file, as `bgzip` and `bcftools` write, must now end with its empty end
  block, which a killed `bgzip` leaves out. All 56 cuts of the review's test file now stop the
  run. This is the "cut short" rule §5.2 already had, which the reader did not keep.
- **An ALT `.` among others stops the run**: "b.vcf.gz, line 2204: ALT 'T,.' holds `.` beside
  other alleles, which VCF does not allow: `.` stands alone, for a record with no ALT. Write the
  ALT without the `.`, here 'T', and run again." Also for a symbolic record, `<DEL>,.`, which was
  skipped; and an ALT of dots alone, `..`, which was skipped as a breakend. `.` alone is still a
  record with no ALT.

On good files the reader gives what `main` gives: a reviewer compared every `parse_vcf` and
`parse_bed` call of the test suite under both, 1,765 VCF reads and 156 BED reads, 260,248
records, and generated files with 1 and 4 threads. A good file is decompressed once; a whole
check runs only after a read has already failed. On a 114 MB gzip with a fault at line 3, the
error comes after 0.28 s instead of at once.

Two things you may want to know:

- A reader thread of `compare()` that meets a fault after the run was told to stop finishes one
  whole decompression before it ends, in a background thread. It changes no result; skipping
  it would add a parameter to `parse_vcf`, whose signature `doc/specs/regions_read.md` fixes.
  Left as it is.
- `tests/test_cli_compare.py::test_warnings_and_info_lines_go_to_stderr`, older than this
  branch, failed once in about 20 runs during the review, and passed on rerun.

## What was done

Done in 630eec2, fixed in 2dd5614 to 8e5665a. The check lives in `pvc_eval._text_files`: a
context manager, `open_numbered_lines`, which the VCF and BED readers open their files through,
checks a gzipped file whole when a fault is raised before the read reached its end, and
`check_gzip_whole` does the same for the callers that find a fault after a line was handed over
(`read_strata`, the regions of `str-regions`, `make_sample_pairs`). The ALT check is in
`_parse_alts`.

Changed in the spec (4a9acae, 233e34f): the two rules, then, from the review, that they hold
whatever finds the fault, that a BGZF file without its end block is cut short, and that the `.`
rule holds for a symbolic record too.

The review found three blockers: a fault found by a caller after the line was read kept the
line's message (three reviewers); a BGZF file cut between blocks was read without error; and no
test held a good file to one read. Each is fixed, and the second looks confirmed the first two;
the third is a test of its own. Fixed besides: the two messages, which the first reader, reading
as a user, found not to say what to do or that the fault may come from the damage; a damaged
file larger than the check's 1 MiB chunk; a file removed between the two reads; and docstrings.
Not taken: an error of the reference FASTA met while a damaged VCF is read gives the damaged
message, which is right, since the damage may be its cause.
