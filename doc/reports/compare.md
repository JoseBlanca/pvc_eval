# Report: compare

*Under way, 20 October 2026. The work of `doc/plans/compare.md`, which builds
`doc/specs/compare.md`, on the branch `compare`; nothing is merged into `main` or pushed.*

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
