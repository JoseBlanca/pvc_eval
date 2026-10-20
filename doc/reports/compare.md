# Report: compare

*Done, 20 October 2026. The work of `doc/plans/compare.md`, which builds `doc/specs/compare.md`,
on the branch `compare`; the last commit of code is 302f5c9, and `main` is merged in at 84aacd9.
Nothing is merged into `main` or pushed.*

## For the owner

**pvc_eval compares two VCF files end to end.** `pvc_eval compare` and the Python function
`compare()` take two VCFs, the reference and the options of `doc/spec.md` §3, and write the table of
§6 and the wide summary you chose. GIAB's HG002 v4.2.1 against itself over chr20 of GRCh38 gives
81,257 rows, all *same genotype*, in 39 s at a peak of 400 MB (Apple M5 Pro, CPython 3.15t), which
is the check of `doc/spec.md` §7.3 on real data, now a test.

What exists now that did not:

- `compare()` and `CompareConfig` (`src/pvc_eval/compare.py`): the stages joined, each VCF read in a
  thread of its own and the STR finder on two threads, scanning only the contigs that have records
  and that `--regions` names, planned from the `.tbi` or `.csi` when there is one;
- the table and the summary (`src/pvc_eval/table.py`), written through `.tmp` files, gzipped for a
  `.gz` path;
- `pvc_eval compare`, with the messages and exit statuses you agreed to;
- the reader of the contigs an index lists (`src/pvc_eval/_vcf_index.py`), with the standard library;
- the test of HG002 against itself (`tests/test_compare_real.py`), skipped without
  `PVC_EVAL_GRCH38` and `PVC_EVAL_HG002_VCF`.

The checks on 84aacd9: `uv sync --locked` clean, ruff and ty pass, `pytest -q` "1306 passed, 5
skipped" (914 passed when the work started), three runs in a row.

**What is asked of you.**

1. **The merge of `compare` into `main`.** `main` was merged into the branch at 84aacd9, after
   `str-regions` on threads, so the merge has no conflict as of now.
2. **Whether `--regions` should also limit the reading.** On HG002, chr20's rows were done at 38 s
   of the 39 s: pvc_eval reads every record of the 22 chromosomes and builds their loci, to raise
   any error the file holds, and skips them. Reading only the contigs `--regions` names, through
   the VCF's index, would bring a chr20 benchmark of a whole-genome VCF to about the 22 s of the
   STR finder's scan of chr20, or a few seconds with `--str-regions`, and an error on another
   chromosome would then go unseen. It needs a reader of a VCF's BGZF blocks from the offsets
   its index gives, in the standard library; about 100 lines, an estimate, not tried. Recommended: yes, as a spec of its own, since it changes
   what the run checks. Nothing waits on it.

**What the spec gained during the work**, each in a commit of its own, none changing a row:

- The rows before an error of the STR finder are those of `compare()` with no thread, which can
  differ by one locus from a chain of `segment_contig` over contigs known in advance (0ad2d19).
- `compare()` no longer repeats the stages' checks of the settings; the command still checks them
  first (eef6533).
- The summary is renamed before the table; any error of the system while writing is a message, not
  a traceback; `--output` and `--summary` that are one file under another spelling are refused
  (b5ac64a).
- The warning for an index that cannot be read says how to make it again (1927bcb).
- `compare()` is annotated as returning a `Generator` (cbcc6e6).

**Left open, none of them changing a row:**

- In one of 15 full runs of the suite, a reviewer saw one test fail without keeping its name; the
  session's 23 runs since have had no failure.
- `pvc_eval str-regions` still ends with a traceback on a full disk; since milestone 3 its Ctrl-C
  ends with "Stopped; nothing was written." and status 130, like `compare`'s.
- Two tests that could be stronger: that the command puts its logger's level back, and that an
  error of a scan not needed is never raised when the loci side passes its contig.

## Before the first milestone

On 1d4a1f8, with `Kind` and `RegionType` already `StrEnum` (ca0463f):
`uv run python -c "import sys; print(sys._is_gil_enabled())"` printed `False`; `uv sync --locked`
was clean; `ruff format --check` gave "40 files already formatted"; `ruff check` and `ty check`
gave "All checks passed!"; `pytest -q` gave "914 passed, 4 skipped".

## Milestone 1: `compare()` with no thread

Built in 37f033b and fixed in 02cb2b9 to da3f025. `compare()`, `CompareConfig` and `PVC_STRATA`
are in `src/pvc_eval/compare.py`; the four-contig case is in `tests/compare_case.py`, its STR
regions checked against `segment_contig` (chr1 100–116 AC and 200–218 AGT, chr2 120–140 ATTC, chr3
100–120 CA, 0-based; none on chr4). On da3f025 `uv run pytest tests/test_compare.py -q` gives "50
passed", the plan asked for at least 15; the five checks pass, with "964 passed, 4 skipped".

- The rows equal those of the nine steps, written out in the test with `segment_contig` over every
  contig with records, with no option and with each option and setting, `--regions` cut inside a
  contig among them. Breaking the segment iterator seven ways fails tests each time (the
  implementer's five and two of the reviewers').
- Only the contigs with records, and with `--regions` only those it names, are scanned, by the
  INFO lines.
- The errors: a VCF out of order and an untyped `--str-regions` region after the nine steps' rows;
  a FASTA byte that is not a letter, raised by the finder, after literal rows.

**The spec changed twice, by the session, after the review.** Both are in commits of their own,
0ad2d19 and eef6533, and neither changes a row of the table.

- Before an error of the STR finder, `compare()` can give one locus more or fewer than a chain of
  `segment_contig` over a list of contigs known before the reading. That chain reads a contig's
  segments at another record than the spec's iterator, which decides a contig only from the
  records `build_loci` has taken. Three reviewers found it, with cases. The spec now names
  `compare()` with `finder_threads=0` as the stages run one after another, which the threads of
  milestone 5 must match. The command line writes no table on an error, so its user sees the same
  message either way.
- `compare()` no longer repeats the stages' checks of the ploidies, the phasing cap, the span limit
  and the chain distance, which were copies with the same messages. The stages raise the same
  `ValueError` before any row, but after an error of a header or a BED.

Also added to the spec: a contig whose records all have no ALT gets no progress line, and the
loci counted as skipped on a contig `--regions` does not name are built without STR regions, so
there can be more of them (be57c54).

**What the review found that mattered, all fixed:** `--regions` cut inside a contig had no test,
so a copy of `select_loci`'s rule could drift (blocker; a second look confirmed the fix); six
settings that no test set; no test that a value at the edge of its range is accepted; no record
on the last contig. Not taken: none.

The implementer added the verb "note", for a pass-through that records something about the items
it hands on, to the table of `python.md`, as that file asks.

## Milestone 4: the index reader

Built in 680c2ce and fixed in 7587c99 to 3a5ba9a. `parse_index_contigs`, in
`src/pvc_eval/_vcf_index.py`, gives the contigs a VCF's `.tbi` and `.csi` list, joined, as a
`frozenset`, reading only their headers with `gzip` and `struct`; the test data in
`tests/data/vcf_index/` were made with bcftools 1.24, with the commands in `tests/data/README.md`.
On 3a5ba9a `uv run pytest tests/test_vcf_index.py -q` gives "32 passed", the plan asked for at
least 6; the five checks pass, with "1109 passed, 4 skipped" three runs in a row (at d756456).

- The reviewer checked the layout on indexes made by `bcftools index -t`, `bcftools index -c -m
  10`, `tabix -p vcf` and `tabix -C -m 20`, over 5,028 header contigs with names such as
  `HLA-A*01:01:01:01`, and over 20,000 contigs spread across 11 BGZF blocks: each gave exactly the
  contigs with records. On GIAB HG002's `.tbi` the reader gives chr1 to chr22 in 0.7 ms, reading
  131,090 of its 1,657,747 bytes.
- **Fixed, a blocker:** a damaged `.tbi` whose name length was below -1 stopped the run with a
  traceback, where the spec ignores it with a warning. A fuzz test of 8,000 damaged indexes found 21
  such files before the fix and none after (second look by the reviewer).
- **Changed in the spec (1927bcb):** the warning for an index that cannot be read now holds also
  when the VCF's other index is read, and says how to make the index again.
- **Changed in the plan (9bd0b91):** the reader does not take the reference, so dropping the
  contigs an index lists that the reference lacks moved to milestone 5, with a test.
- Also fixed: five checks of the header that no test guarded, the contigs returned as a set since
  no caller uses their order, and two helpers renamed for what they do. Not taken: none.

## Milestone 2: the table and the summary

Built in 848a543 and fixed in 3793fd2 to d756456. `TABLE_COLUMNS`, `SUMMARY_COLUMNS`,
`format_table_row` and `write_table` are in `src/pvc_eval/table.py`; the opener of an output's
temporary file is now one function, `open_temporary_output` in `src/pvc_eval/_text_files.py`, used
by the table and by `str-regions`. On d756456 `uv run pytest tests/test_table.py -q` gives "97
passed", the plan asked for at least 12; the five checks pass.

- The worked row of the spec's "The table" comes out as its exact line. QUAL is written 50.0 as
  `50`, 29.53 as `29.53`, 1e-05 as `1e-05`, 1e+20 as `1e+20`, 1234567.0 as `1234567`.
- The wide summary over all sample pairs, as the owner decided, with its order, its zeros and its
  totals as literals; 200,000 rows written with a peak of 2.4 MB under tracemalloc.

**Fixed, two blockers, each confirmed by a second look:**
- `--output` and `--summary` that are one file under another spelling were not refused. On macOS's
  file system, which ignores letter case, `loci.tsv` and `Loci.tsv` overwrote each other and lost
  an earlier table. A path that is the other's `.tmp` gave a wrong file with no error. Both are now
  refused, by comparing the two open files.
- A write error, such as a full disk, left both `.tmp` files, because closing a file whose write
  failed raises again and stopped the cleanup.

**Changed in the spec (b5ac64a):** the summary is renamed before the table. A failure between the
two renames, which nothing can make atomic, then leaves this run's summary beside the earlier
table. Any error of the system while writing is a message naming the option, with no traceback,
since it is a fault of where the output goes. The gzip header holds no file name, and lines end in
`\n` everywhere.

Also fixed: tests that could not fail on the QUAL's form, the outcomes left out, file B's invalid
reason, the strata's order, UTF-8, the line ends and the gzip level. Not taken: none.

**Left open, outside this plan:** in `pvc_eval str-regions`, a full disk after the BED is opened
still ends in a traceback, as before this plan. Its opener now gives messages for every error of
the system, but its writing loop does not.

## Milestone 3: the `pvc_eval compare` command

Built in b800ce0 and fixed in 171a415 to 14d5fed. pvc_eval now runs end to end from the command
line, still with no thread. On 14d5fed `tests/test_cli_compare.py` runs 48 tests, the plan asked
for at least 10; the five checks gave "1157 passed, 4 skipped", run by the implementer on a copy
of the commit, since milestone 5's work in progress was in the shared tree at the time.

- A real Ctrl-C, sent to the installed command at 60 moments of a run, gave status 130, "Stopped;
  nothing was written." and no `.tmp` every time but one. That one came in the first milliseconds,
  while the options were parsed, and gave a traceback; it is fixed.
- `--help` was read by a first reader standing for a bioinformatician new to pvc_eval. It could
  run a comparison but not read its result. The help now lists the outcomes and the kinds, says
  what the STR finder does, which contigs it scans and how long it takes, and how to write the STR
  regions once with `pvc_eval str-regions` for later runs.

**Fixed, a blocker, confirmed by a second look:** no test checked that `--sample-map` or
`--str-regions` reached `compare()`. Dropping either gave a wrong table with every test passing.

Also fixed: every line printed twice when a Python caller had set up logging; no test that a
defect of pvc_eval keeps its traceback; the ploidy options tested at one bound only; `--summary`
and `--all-rows` without their default in `--help`.

**Moved to milestone 5:** `compare()` annotated as returning a generator, which the command closes,
so that its `cast` goes; and a test that the command runs `compare()` with two finder threads,
which cannot be told apart from none before the threads exist.

**Not taken:** a Ctrl-C between the renames of the summary and of the table, a window of
microseconds, prints "nothing was written" after the summary was replaced. The spec admits the
window, and the message is not worth a third state. **Left open, a nit:** no test checks that the
command puts the level of pvc_eval's logger back after a run; the second look found that removing
it still passes.

## Milestone 5: the threads

Built in 7419f7e, cbcc6e6 and 0c43267, and fixed in 5cc7549 to e97e5ec. With `finder_threads` of 1
or more, each VCF is read in a thread of its own into a `queue.Queue` of at most 156 batches of 64
records, and the STR finder runs on a pool of daemon threads, a few contigs ahead of the loci, each
thread with its own handle on the FASTA. On e97e5ec `tests/test_compare_threads.py` runs 113
tests, the plan asked for at least 12, and passed in each of ten runs in a row.

- **Stress, by the correctness reviewer**, on 480 random cases of 50 to 400 contigs, with
  indexes, `--regions`, errors in the FASTA and in the VCFs, 1 to 8 finder threads, run to the
  end and closed at a random row: 7,592 threaded runs before the fixes and 7,592 after, plus
  2,368 and 1,592 with random delays. No hang, no row or first error different from
  `finder_threads=0`, no exception reaching `threading.excepthook`, no reader alive after a close.
  Memory grew with the number of contigs only by the records in flight, not by anything `compare`
  holds.
- **Fixed, a blocker, confirmed by a second look:** three rules had code and no test (planning
  from the indexes; the limit of contigs ahead; a contig the loci side passes leaving the line),
  nor had the finder threads being daemon threads. Each now has a test, and each break fails one.
- Also fixed: a lost wakeup in the reader's buffer that no test could catch (the buffer is now
  `queue.Queue`, whose `shutdown` wakes a blocked reader, 50 lines of locking gone); Ctrl-C while
  the indexes are read became a `RuntimeError`; a stream never closed in a script printed an error
  as Python exited; tests that turned one hang into hours; two paths no input reached, removed.
- From milestone 3's review: `compare()` returns a `Generator`, which the command closes, and a
  test shows the command runs it with two finder threads (spec changed in cbcc6e6).
- **Not taken:** two casts of `parse_vcf`'s result in `compare.py`, which would go only if
  `pvc_eval.vcf` annotated `parse_vcf` as a generator; that module is not this plan's to change.
- **Left open:** no test catches an error of a scan raised when the loci side passes that contig
  rather than never; and in one of 15 full runs of the suite at e97e5ec the second-look reviewer
  saw one test fail, without keeping its name. The session ran the suite 20 times in a row at
  302f5c9 and 3 times at 84aacd9 with no failure.

## Milestone 6: a file against itself over real data

Built in 5bafe18 and fixed in 60abdf3 to 302f5c9. `tests/test_compare_real.py` compares GIAB's
HG002 v4.2.1 benchmark VCF with itself over chr20 of GRCh38, through `PVC_EVAL_GRCH38` and
`PVC_EVAL_HG002_VCF`, and is skipped with a message naming each that is missing. It asserts
81,257 rows, each of the five outcomes doc/spec.md §7.3 allows, and the two sides of each row
equal field by field; fast tests check the same rule on literal rows and on a small file with a
`./.`, a malformed GT and 14 unphased heterozygous records, which give *missing in both*, *invalid
in both* and *too complex*. With the variables it gives "19 passed" in 39.1 s.

- **It can fail:** file B with one GT of HG002 changed (chr20:88108, 0/1 to 1/1) breaks exactly one
  row, *same alleles, other genotype*; a `./.` there gives *missing in B*; the record removed gives
  *only in A*.
- Fixed after review: the rule on the three outcomes above never ran on HG002, which has none of
  them; losing rows did not fail the test; three of the six paired fields were not compared.

## The check of the whole plan

On 84aacd9, the session ran the command on HG002 against itself, with `--regions` naming chr20 and
`--summary`, the table gzipped:

```
uv run pvc_eval compare --reference GRCh38.fna --a HG002.vcf.gz --b HG002.vcf.gz \
    --regions chr20.bed --output self.tsv.gz --summary self_summary.tsv
```

It exited with 0 in 39.1 s, at a peak of 400 MB resident (`/usr/bin/time -l`, Apple M5 Pro,
CPython 3.15.0t with the GIL off, 20 October 2026), at 302f5c9 before the last merge of `main`.
The finder scanned chr20 alone. The summary:

```
outcome         stratum  SNP    indel  STR    complex  total
same genotype   .        65112  3663   12171  311      81257
```

`main` was merged into the branch twice: at d89ba74, which brought the faster STR finder with no
conflict, and at 84aacd9, which brought `str-regions` on threads, with conflicts in the imports of
`cli.py` and in `segments.py`, where the branch's shared opener of an output replaces `main`'s copy.
After it, ruff, ty and "1306 passed, 5 skipped" three runs in a row.

## For whoever next revises the skills

The owner can stop reading here.

- **Two subagents in one worktree work, if each commits by path.** Milestones 2 and 4, and 3 and
  5, and 5's fixes and 6, ran side by side. Their checks, run in the middle, saw the other's
  uncommitted files: two subagents reported a failure that was the other's work in progress. The
  coding skill could say that a subagent sharing a tree runs its checks on a copy of its commit.
- **"A rule with code and no test" was the finding of every milestone.** In five of six
  milestones a reviewer broke a rule the spec gives and every test passed: `--regions` cut inside
  a contig, the options that reach `compare()`, the index's magic, the limit of contigs ahead. The
  implementers tested what they built rather than each sentence of the spec. A line in the
  implementer's prompt, "for each rule of the spec's section, name the test that fails without
  it", might catch these before the review.
- **The spec's review missed what the code review found.** The rows before an error of the
  finder, which the spec promised equal to a chain of `segment_contig`, differ by a locus; three
  code reviewers found it with cases. A spec that promises equality between two runs is worth a
  reviewer who builds both.
- **Tokens, as the `Agent` tool reported them** (implementation including its fixes; the reviews,
  second looks included):

  | milestone | implementation | reviews |
  |---|---|---|
  | 1 | 227 k | 357 k |
  | 2 | 166 k | 237 k |
  | 3 | 166 k | 206 k, a first reader of `--help` included |
  | 4 | 107 k | 161 k |
  | 5 | 259 k | 434 k |
  | 6 | 96 k | 43 k |
