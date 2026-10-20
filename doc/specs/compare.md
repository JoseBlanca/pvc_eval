# Compare

*Status: approved by the owner on 20 October 2026, with the answers to its seven open points
listed at the end. No code of this part exists yet. It
joins stages that exist on `main` and are tested: the readers (`doc/spec.md` §5.2), the STR finder
(`doc/specs/str_finder.md`), the loci, kind, strata and regions (`doc/specs/loci.md`), and the
sample pairs, reconciliation and evaluation (`doc/specs/reconciliation.md`). It builds what
`doc/spec.md` §3, §4 and §6 promise above them.*

This part makes pvc_eval usable: `compare()` runs the stages from two VCF files to the rows of the
table, and `pvc_eval compare` writes those rows as the table and, when asked, a summary of counts.
It also decides where the time of a run goes when no `--str-regions` BED is given. The STR finder
then scans the reference, and on GRCh38 that scan is the longest part of a run. One contig takes
from 13 s to 89 s on one thread with the faster finder of the branch `perf/str-finder`, and
five to nine times longer with the finder on `main` (below).
So `compare()` scans only the contigs whose loci need STR regions. It scans them in threads, a few
contigs ahead of the loci, while the two VCF files are read in threads of their own, as the owner
decided on 20 October 2026. The owner's answers to the spec's seven open points, the same day,
are listed at the end, under "Decided by the owner".

Positions in the cases are 1-based and inclusive, as in a VCF and in the table. In the Python
interface they are 0-based and half-open, as everywhere inside pvc_eval.

## Vocabulary

The terms of `doc/spec.md` §2 and of the specs above keep their meaning. These are added:

| term | meaning | not |
|---|---|---|
| **the loci side** | the work done in the caller's thread: `build_loci`, `select_loci`, `reconcile_loci`, `locus_kind`, the stratum and `make_outcome_rows`, one locus after another | |
| **a reader** | `parse_vcf` over one of the two files, run in a thread of its own, handing its records to the loci side through a bounded buffer | "parser" |
| **the finder** | `segment_contig`, run on a pool of threads, one contig per thread at a time | |
| **a contig reached** | a contig on which a reader has read at least one VCF record | |
| **a contig planned** | a contig the finder may scan before any reader reaches it, because the indexes of the VCFs list it (below) | "prefetched" |
| **index** | the `.tbi` or `.csi` file that bcftools or tabix writes beside a bgzipped VCF, `a.vcf.gz.tbi` or `a.vcf.gz.csi` | |

## What `compare()` does

### The stages, joined

`compare()` does what the stages of `doc/spec.md` §4 do one after another. Run alone, those stages
would be:

1. open the reference, an `IndexedFasta`;
2. make the sample pairs, `make_sample_pairs(a, b, sample_map=...)`, which reads the two headers
   and stops the run for the errors of `doc/specs/reconciliation.md` §3;
3. read the `--regions` BED with `read_regions` and the `--strata` BED with `read_strata`, each
   whole, when given;
4. read each file with `parse_vcf(path, reference, span_limit=..., ploidy=...)`;
5. take the segments: from `read_str_regions` when `--str-regions` is given, and otherwise from
   `segment_contig` over the contigs whose loci need them, each decided when `build_loci` reaches
   it, by the iterator of "How the segments reach `build_loci`";
6. `build_loci(records_a, records_b, segments, reference.contigs, chain_distance=...)`;
7. `select_loci(loci, regions)`;
8. `reconcile_loci(loci, reference, sample_pairs, ploidy_a=..., ploidy_b=..., phasing_cap=...)`;
9. for each reconciled locus, `make_outcome_rows(reconciled, sample_pairs, kind=locus_kind(locus),
   stratum=...)`, where the stratum is `pvc_stratum(locus)` with `--strata pvc`,
   `strata.stratum(locus)` with a BED, and `None` without `--strata`.

`compare()` yields the rows in the order these give them: locus by locus, and in each locus in the
order of the sample pairs. Steps 1 to 3 run when the stream starts, at the first `next()`, so that
an error in a header or a BED stops the run before any row. Every other error comes when the
stream reaches it, as in the stages. The stage functions are not changed: the threads belong to
`compare()`.

**Its rows are the rows of the stages run one after another, whatever its threads.** The same rows,
in the same order. The same first error, with the same message, after the same rows. The stages
run one after another are `compare()` with `finder_threads=0`: the nine steps in the caller's
thread, with step 5 decided by that iterator. A chain of `segment_contig` over a list of contigs
known before the reading gives the same rows, but before an error of the finder it can give one
locus more or fewer, since it reads a contig's segments at another record (found by the review
of milestone 1, 20 October 2026). This is
what the threads below must keep, and the first check of "How it is verified".

### Which thread does what

With the default of two finder threads, a run has these threads:

- **The caller's thread** runs steps 1 to 3, then the loci side, steps 6 to 9, and yields the rows.
- **Two reader threads**, one per file, each running `parse_vcf` (step 4) into a buffer of at most
  10,000 VCF records. A reader hands its records over in batches, lists of 64 records, the last of
  a file perhaps shorter, so that the buffer is locked once per batch and not per record. A reader
  that finds its buffer full waits until the loci side takes from it.
  `parse_vcf` keeps its own pool for the sample columns, 4 threads by default (`doc/spec.md` §5.2).
- **The finder's pool**, two threads, each scanning one contig with `segment_contig`. A contig's
  segments are kept as its STR regions alone: `build_loci` drops GENERIC segments, and they are
  half of all segments (on chr20:30,000,001–34,000,000 of GRCh38, 23,026 of 46,051).

With `--str-regions` there is no finder pool. `read_str_regions` reads the BED in the caller's
thread, as `build_loci` asks for its lines.

`finder_threads` sets the size of the pool, 2 by default. With 0, `compare()` makes no thread
of its own. The readers are read as `build_loci` asks for their records, and each contig is
scanned in the caller's thread when the loci side first needs its segments. That is the
stages run one after another, which the threaded run is tested against. 1 is a real choice
too: one finder thread still works while the readers and the loci side do, and that overlap
is what this design is for. So `finder_threads` counts the finder's threads alone, and 0, not
1, means no thread (`.claude/skills/coding/python.md` asks for 1 to mean no pool; it is
counted differently here for that reason).

### Which contigs are scanned

The finder resolves a whole contig before it yields its first segment, and keeps no state between
calls (the docstring of `segment_contig`). So the work can overlap between contigs, and between
the finder and the rest, but not within one contig. What saves most time is scanning only the
contigs whose loci need STR regions: a benchmark on chr20 must not wait for chr1 to chr19.

**A contig needs scanning** when the loci side builds a locus on it that `--regions` does not
skip. So it is scanned only when both hold:

- either VCF has a record on it;
- `--regions` names it, when `--regions` is given. Every locus on another contig is skipped, and
  its loci are built with no STR regions, which changes no row, since no row comes from them.

**When a contig's scan starts.** A contig is put in line for the finder when it is planned, or as
soon as a reader reaches it, whichever comes first, and only if `--regions`, when given, names it.
The finder takes the contig of the lowest rank in the `.fai` among those in line.

- **Planned contigs.** A VCF that is bgzipped with an index beside it plans the contigs its index
  lists, before any record is read. An index lists every contig with a record in its file, in a
  header that the standard library reads (`gzip` and `struct`). The reading was tried on 20
  October 2026 on a `.tbi` and a `.csi` written by bcftools 1.24, and on the `.tbi` of GIAB's
  HG002 v4.2.1 for GRCh38, which lists chr1 to chr22. A contig with a record in one file needs
  scanning whatever the other file holds, so one file's index is enough to plan its contigs.
- **Reached contigs.** For a VCF without an index, a contig is scanned only once a reader reaches
  it, with no guess about the contigs to come.

An index decides only when a scan starts, never whether a contig's loci get their STR regions. A
contig that an out-of-date index leaves out is still scanned when a reader reaches it. A contig it
lists and no file has a record on is scanned for nothing, and its segments are dropped. So the
rows never depend on the index. An index that cannot be read is ignored, with a warning:
"a.vcf.gz.tbi could not be read as the index of a VCF; the contigs of a.vcf.gz are scanned for
STR regions as the reading reaches them, which can take longer." When both a `.tbi` and a `.csi`
are beside a VCF, both are read and their contigs joined, since either may be the one out of date.
Their contigs that the reference lacks are ignored, since a record on one stops the run when the
reader reaches it.

`--regions` never makes a contig planned. A contig it names is scanned only when a record is on
it, known from an index or from a reader. GIAB's benchmark regions name chr1 to chr22, and
scanning them all for a VCF of chr20 would cost the hour and a half this section is there to save
(see "How it runs"). This scans fewer contigs than "scan only the contigs `--regions` names" read
literally, and never more.

**How far ahead.** A contig is **ahead** from the moment its scan starts until the loci side
starts reading its STR regions. The pool starts the next contig in line only while a thread is
free and fewer than `finder_threads + 1` contigs are ahead. So with two threads at most two
contigs are being scanned and at most three are ahead, scanned or not. One exception keeps the
loci side from waiting on the limit: when the loci side waits for a contig that has not started,
that contig starts as soon as a thread is free, whatever the number ahead. It still waits for a
free thread, so behind two long scans of later contigs it can wait for both.

A contig the loci side passes without reading its STR regions, because no record turned out to be
on it, leaves the line if its scan has not started. A scan already started is not interrupted:
its segments are dropped when it ends, and it stops counting as ahead when the loci side passes
it.

### How the segments reach `build_loci`

`build_loci` takes one iterator of segments and reads it only as far as its records need: before
it handles a record, it reads the segments up to the first that comes after the end of the
record's span, by contig rank and then position (`doc/specs/loci.md`, "What `build_loci` holds").
`compare()` hands it an iterator that goes through the contigs of the `.fai` in order, and that
decides for each contig, when it gets there, whether to give its STR regions. It cannot ask the
readers which contigs are still to come without waiting on them, and a reader may be waiting on a
full buffer that only the loci side empties. So the iterator decides from the records the loci
side has taken so far, and from nothing else. A record is **taken** when it is handed to
`build_loci`, one record at a time, and not when its batch leaves a reader's buffer. Counted by
the batch, a record later in the batch could decide a contig while `build_loci` still handles an
earlier one, and an error of that contig's scan would come before rows that the stages run one
after another give before it.

For each contig *X*, in the order of the `.fai`:

1. If a record on *X* has been taken from either file, and `--regions`, when given, names *X*,
   it yields the STR regions of *X*, waiting for the finder if the scan is not done.
2. Otherwise it yields one GENERIC segment at the first base of *X*, `GenomeSegment(X, 0, 1,
   RegionType.GENERIC, None)`. This segment is not part of the segmented genome. It tells
   `build_loci` that no segment before position 0 of *X* is to come. `build_loci` asks for a
   segment past the last one it has read whenever that one starts at or before the end of the
   record it handles, so every record past the last STR region of the contig before *X* gets this
   segment. `build_loci` drops GENERIC segments, so this one changes no locus.
3. When asked for the next segment after that, it decides *X* again by step 1; if *X* still has
   no record taken, it goes on to the next contig.

Step 3 is right because `build_loci` asks for a segment past the one at (*X*, 0) only for a record
whose span ends at or after it, so a record on *X* or a later contig. That record has been taken.
`build_loci` handles records in order, so no record on a contig before it is still to come. If no
record on *X* has been taken by then, *X* has none.

So a contig in the middle of the `.fai` with no record costs one GENERIC segment and no scan. The
iterator ends after the last contig of the `.fai`. With `finder_threads=0`, step 1 runs
`segment_contig` in the caller's thread.

For example, a reference of chr1 to chr4, file A with records on chr1 and chr3, file B on chr2,
neither indexed, no `--regions`. The readers reach chr1 and chr2 at once, and later chr3, and each
contig goes in line as it is reached. When `build_loci` handles its first record, on chr1, the
iterator yields chr1's STR regions, after the scan of chr1 ends. Past them, it yields the GENERIC
segment at (chr2, 0). The first record on chr2 asks again: a record on chr2 has been taken, so it
yields chr2's STR regions, then (chr3, 0), and chr3's STR regions in the same way. A record on
chr3 past chr3's last STR region gets (chr4, 0). No record on chr4 is ever taken, so chr4 is never
scanned, and the iterator is not asked again.

### Errors, and stopping early

**An error in a thread reaches the caller where it would have without threads.** An exception a
reader raises, a VCF out of order at line 2204 say, goes into its buffer after the records read
before it. The loci side raises it when it takes it, the same exception object with its message
and traceback. An exception the finder raises for a contig, a byte of the FASTA that is not a
letter, is raised when the iterator above asks for that contig's segments. The rows before it are
those `compare()` gives with `finder_threads=0` before the same error. An error in a contig that turns out not
to be needed is never raised, as it is not without threads.

Warnings are logged from the thread that finds them, through the module's logger, as the stages
log them. `parse_vcf` logs its warnings when a file's reading ends, which with a reader thread can
come before the rows of the last loci rather than after them. The text is the same.

**The stop.** `compare()` holds one `threading.Event`, the stop, which every thread checks. It is
set when the rows end, so that a scan planned from an index and not needed does not start after
the last row, and when the caller stops reading the rows early, by closing the iterator, by an
exception in the code that reads it, or by Ctrl-C. Once it is set:

- a reader blocked on its full buffer wakes, closes its `parse_vcf` stream, which closes its file
  and its pool, and ends. `compare()` waits for both readers to end, which takes at most one batch
  of `parse_vcf`'s work, before it closes the reference;
- the finder starts no new contig. A contig already being scanned is not waited for. Its thread
  ends when the scan does, up to 89 s later on GRCh38 with the faster finder (below),
  and its segments are dropped. A scan that calls `segment_contig` after the reference is closed
  gets a `ValueError`, which is dropped;
- every thread catches every exception it meets. Before the stop it passes it to the loci side,
  as above; after the stop it drops it. So nothing reaches Python's handler of exceptions in
  threads, which would print a traceback after "Stopped; nothing was written."

The finder's threads are daemon threads, `threading.Thread(daemon=True)`, so that a scan in
progress does not keep Python from exiting for up to 89 s after Ctrl-C, and several minutes with
the finder on `main`. A
`ThreadPoolExecutor`, which `.claude/skills/coding/python.md` asks for, makes threads that Python
waits for when it exits: with one worker sleeping 3 s, the program ended 3.03 s after its last
line (measured by the spec's reviewer on 20 October 2026). The finder departs from `python.md`
for that reason, and keeps its rule that a thread's exception reaches the caller.

## The Python interface

`compare` and `CompareConfig` live in `pvc_eval.compare`; the table and the summary in
`pvc_eval.table`.

`CompareConfig` holds the options of `doc/spec.md` §3 that change the rows, with the defaults of
the stages. `strata` is a path to a BED, or `PVC_STRATA` for `--strata pvc`. A sentinel is used
because a path named `pvc` is a BED (`doc/spec.md` §3). The command line tells them apart by the
text of the option, before it makes a path of it: `Path("./pvc")` equals `Path("pvc")`, so the
`./` that `doc/spec.md` §3 asks of a file named `pvc` is lost once it is a path. `compare()` checks the ploidies, the
phasing cap, the span limit and the chain distance when the stream starts, and raises
`ValueError` naming the field. The command line checks the same values first, and turns a wrong
one into its own message (see "The `compare` command").

```python
PVC_STRATA = sentinel("PVC_STRATA")

@dataclass(frozen=True, slots=True, kw_only=True)
class CompareConfig:
    regions: Path | None = None
    strata: Path | PVC_STRATA | None = None
    str_regions: Path | None = None
    sample_map: Path | None = None
    ploidy_a: int = DEFAULT_PLOIDY              # 1 or 2
    ploidy_b: int = DEFAULT_PLOIDY
    phasing_cap: int = PHASING_CAP              # 1 or more
    span_limit: int | None = DEFAULT_SPAN_LIMIT  # 1 or more, or None for no limit
    chain_distance: int = DEFAULT_CHAIN_DISTANCE  # 0 or more

DEFAULT_FINDER_THREADS = 2

def compare(a: Path, b: Path, reference: Path, config: CompareConfig = CompareConfig(), *,
            finder_threads: int = DEFAULT_FINDER_THREADS) -> Iterator[OutcomeRow]: ...
```

`compare()` yields every row, *reference in both* included, as `doc/spec.md` §3 says. It raises
what the stages raise: `PvcEvalInputError` for a fault of the input, and `ValueError` for a
setting out of range or a `finder_threads` that is not a whole number of 0 or more. It has no
`evaluator` parameter until a wrapper of `doc/spec.md` §8 exists (decided by the owner on 20
October 2026, point 1).

The table and the summary are written by one function, since the command line reads the rows
once and writes both from them. `format_table_row` gives one row's fields as text, for a user who
writes the table another way.

```python
TABLE_COLUMNS: tuple[str, ...]      # the column names of doc/spec.md §6, in its order
SUMMARY_COLUMNS: tuple[str, ...]    # ("outcome", "stratum", "SNP", "indel", "STR", "complex", "total")

def format_table_row(row: OutcomeRow) -> tuple[str, ...]: ...

def write_table(rows: Iterable[OutcomeRow], output: Path, *, summary: Path | None = None,
                all_rows: bool = False) -> None: ...
```

`write_table` writes every row to `output`, except, without `all_rows`, those whose outcome is
*reference in both*; it counts every row in the summary. It raises `PvcEvalInputError` for an
output it cannot write (see "Writing the files"), and lets through what `rows` raises, after
removing what it had written.

## The table

`doc/spec.md` §6 gives the columns and what each holds. This section gives how each field of
`OutcomeRow` (`doc/specs/reconciliation.md` §6) becomes text.

The file is UTF-8 text, one line per row, its fields separated by tabs, each line ending in `\n`.
The first line is the column names, with no `#`, which pandas, R and `csvtk` read as names: `chrom`, `start`, `end`, `kind`,
`stratum`, `str_regions`, `sample_a`, `sample_b`, `outcome`, `alleles`, `genotype_a`,
`genotype_b`, `possible_a`, `possible_b`, `quality_a`, `quality_b`, `gq_a`, `gq_b`, `records_a`,
`records_b`, `invalid_reason_a`, `invalid_reason_b`. A field that is `None` is written `.`.

| column | from the field | written |
|---|---|---|
| `chrom` | `contig` | as it is |
| `start` | `start` | `start + 1`: 1-based |
| `end` | `end` | `end` as it is: one past the last 0-based position is the last 1-based one |
| `kind` | `kind` | the enum's value: `SNP`, `indel`, `STR`, `complex` |
| `stratum` | `stratum` | as it is; `.` for `None` (no `--strata`), and `.` from a BED no interval of which the locus overlaps |
| `str_regions` | `str_regions` | each region as `start+1-end:type`, with `:motif` after an STR segment's, and `start+1-end` alone for a region with no type; joined with `,`; `.` for none |
| `sample_a`, `sample_b` | `sample_a`, `sample_b` | as they are |
| `outcome` | `outcome` | the enum's value, `same alleles, other genotype` for one |
| `alleles` | `alleles` | the locus alleles joined with `,`; a locus allele that is empty, a deletion of the whole locus, is written `-` |
| `genotype_a`, `genotype_b` | `genotype_a`, `genotype_b` | the indices joined with `/`, `MISSING` written `.` |
| `possible_a`, `possible_b` | `possible_a`, `possible_b` | a whole number |
| `quality_a`, `quality_b` | `quality_a`, `quality_b` | the shortest decimal that reads back as the same number, with no `.0` at the end of a whole number |
| `gq_a`, `gq_b` | `gq_a`, `gq_b` | a whole number |
| `records_a`, `records_b` | `records_a`, `records_b` | a whole number |
| `invalid_reason_a`, `invalid_reason_b` | `invalid_reason_a`, `invalid_reason_b` | the enum's value |

One row, worked. The reference at 1-based 100–111 is `CATATATATATG`, and an AT repeat at
101–110 is an STR segment. File A has `100 CAT C 0/1`, QUAL 29.5 and GQ 45, for HG002, and file B
no record there. The deletion of AT slides right as far as 109–110, where the G at 111 stops it, so its span is
100–110; the record touches the repeat, and the locus is 100–110, 0-based (99, 110). With `--strata pvc`:

```
chrom  start  end  kind  stratum  str_regions     sample_a  sample_b  outcome    alleles                  genotype_a  genotype_b  possible_a  possible_b  quality_a  quality_b  gq_a  gq_b  records_a  records_b  invalid_reason_a  invalid_reason_b
chr20  100    110  STR   str      101-110:str:AT  HG002     HG002     only in A  CATATATATAT,CATATATAT  0/1         0/0         1           1           29.5       .          45    .     1          0          .                 .
```

(The columns are aligned here for reading; in the file they are separated by one tab.)

Three fields can look alike and differ. A haploid genotype missing on its one haplotype,
`(MISSING,)`, is written `.`, as is a sample with no genotype in the row (`None`): `possible_a`
tells them apart, 1 or more against `.` or 0. A diploid `(MISSING, MISSING)` is `./.`. A stratum of
`.` means no `--strata` or no class from the BED, and the user knows which they gave.

**The QUAL.** `VcfRecord.qual` is a float, so the text the VCF wrote is not kept. `50` and `50.00`
are both 50.0 and both written `50`; `29.53` is written `29.53`; `1e-05` is written `1e-05`,
as Python's `repr` gives it. Writing QUAL as the VCF wrote it needs `VcfRecord` to keep the text,
which changes `pvc_eval.vcf`; the owner kept the shortest form on 20 October 2026 (point 3).

## The summary

With `--summary`, a second file counts the rows in a wide table (decided by the owner on 20
October 2026, point 2, to be tried first): one line per outcome and stratum, one column per kind,
and a total. Its columns are `outcome`, `stratum`, `SNP`, `indel`, `STR`, `complex` and `total`,
and its first line is these names.

- The counts are summed over every sample pair. A cohort's summary therefore does not say which
  pair disagrees; the table does.
- Every row `compare()` yields is counted, *reference in both* included, with or without
  `--all-rows`. Loci that `--regions` skips give no row and are not counted.
- A line is written for each outcome and stratum that have at least one row; in it, a kind with
  no row is written `0`. Without `--strata` the stratum is `.` on every line, so the table keeps
  one layout.
- The lines come by outcome, in the order of `doc/spec.md` §5.5's table, and then by stratum,
  sorted by code point. `total` is the sum of the four kinds.

```
outcome            stratum  SNP    indel  STR   complex  total
reference in both  generic  1204   88     0     3        1295
same genotype      generic  80311  5102   0     240      85653
same genotype      str      0      0      2210  0        2210
different          str      0      0      61    0        61
```

(Aligned for reading; separated by tabs in the file. The counts are made up for the example.)

## Writing the files

The table and the summary are each written to a temporary file beside it, `<path>.tmp`, and
renamed to their path when every row is written, as `pvc_eval str-regions` writes its BED. A run
stopped by an error or by Ctrl-C removes the temporary files and leaves the table and the summary
as they were before it. A run killed by a signal that Python cannot catch leaves the `.tmp` files
behind.

A path ending in `.gz` is written compressed with gzip, at level 6, which `zcat`, pandas and R
read. It is not bgzipped, so tabix cannot index it; the owner decided on 20 October 2026 that gzip
is enough (point 3).

Before the first row is read, `write_table` opens both temporary files. It raises
`PvcEvalInputError` naming the option for a path that is a directory, whose directory does not
exist, or that cannot be written there, with the messages `str-regions` gives for `--output`, and
for a `--summary` equal to `--output`.

## The `compare` command

```
pvc_eval compare --reference ref.fa --a a.vcf.gz --b b.vcf.gz --output loci.tsv
                 [--regions regions.bed] [--strata strata.bed|pvc] [--str-regions str.bed]
                 [--sample-map pairs.tsv] [--ploidy-a N] [--ploidy-b N] [--phasing-cap N]
                 [--span-limit N|none] [--chain-distance N]
                 [--summary summary.tsv] [--all-rows]
```

These are the options of `doc/spec.md` §3, with the same meanings and defaults, less
`--evaluator` (point 1). The command runs `compare()` with two finder threads and writes its
rows with `write_table`. The STR finder runs with its defaults. A user who wants other parameters
runs `pvc_eval str-regions` with them and passes its BED with `--str-regions`. Nothing else is
set from the command line.

**Values out of range** stop the run before anything is read, with status 1 and a message naming
the option, as `str-regions` does for its own: "--ploidy-a 3: a ploidy is 1 or 2. Correct the
option, and run again." `--span-limit` takes a whole number of 1 or more or the word `none`;
`--chain-distance` a whole number of 0 or more; `--phasing-cap` one of 1 or more. A value that is
not a number at all is argparse's error, with its usage line and status 2, as in `str-regions`.

**Messages**, all on stderr (point 5):

- the warnings of the stages, as they log them: the samples in no pair, the VCF records skipped or
  cut, an index that could not be read;
- one line when the finder starts a contig, since on GRCh38 a contig takes up to 89 s,
  and several minutes with the finder on `main`, and nothing else would show: "chr1: finding STR regions in 248,956,422 bases.";
- one line when the loci of a contig are done: "chr20: 31,204 loci compared, 1,880 skipped by
  --regions; 4 min 12 s since the start." The clause about `--regions` is there only when it is
  given. It counts loci, not rows, since `compare()` does not know which rows the command writes.
  On a contig `--regions` does not name, the loci skipped are built without STR regions, so they
  can be more than with the STR regions: two records in one repeat count as two.
  It comes when the first locus of a later contig is built, or when the rows end, since only then
  does `build_loci` give a contig's last locus; so it can come after the next contig's scan. A
  contig with records and no locus, whose records all have no ALT, gets no line: nothing was
  compared there.

These are logged by `compare()` at INFO through its logger. The command line sends INFO and above
to stderr, the message alone on each line. A Python user sees them only by setting up logging.

**Exit statuses:**

| status | when |
|---|---|
| 0 | the table, and the summary if asked, are written |
| 1 | a fault of the input (`PvcEvalInputError`): its message, with no traceback, as `str-regions` does; or a defect of pvc_eval, with its traceback, as Python ends on an exception it does not catch |
| 2 | an option argparse refuses: its usage line and its message |
| 130 | Ctrl-C: "Stopped; nothing was written." with no traceback, the temporary files removed |

## How it runs

**What is held, and how it grows.** Besides what each stage holds by its own spec, `compare()`
holds:

| what | at most | measured |
|---|---|---|
| each reader's buffer | 10,000 VCF records | 665 bytes per record with one GATK-like sample, 6,350 bytes with 200: so 7 MB and 64 MB per file |
| the finder's scans in progress | `finder_threads` contigs at once | 1.15 GB for chr1, 339 MB for chr20 (below) |
| contigs scanned and waiting for the loci side | `finder_threads + 1` contigs ahead in all, scanned or being scanned, and the contig the loci side reads | STR regions take about 0.85 MiB per Mb of contig: chr1, 249 Mb, about 210 MiB |

None of these grows with the number of contigs or the length of the VCFs. The bytes per record
were measured with tracemalloc on records `parse_vcf` held, from VCFs that
`benchmarks/make_vcf.py` wrote over chr20 of GRCh38 (chr20:1–5,000,000 with 200 samples
`--column gatk`, 18,127 records, and chr20:1–20,000,000 with one sample `--column gatk --mean-gap
650`, 30,401 records). The STR regions per Mb were measured on chr20:30,000,001–34,000,000: 46,051
segments held 6.8 MiB, half of them STR regions. All of it on 20 October 2026, at aaaf31b, on the
Apple M5 Pro, CPython 3.15.0t with the GIL off.

**What the threads gain.** The STR finder's time was measured by the performance review on the
branch `perf/str-finder`, which sent its figures on 20 October 2026. That branch makes the finder
faster and is not merged into `main` yet; the owner decides. Its figures, at its commit c31b3b0,
on the same machine, one run of `segment_contig` per chromosome of GRCh38 one after another in
one process, default parameters:

| chromosome | time | peak of one chromosome alone |
|---|---|---|
| chr1 | 77.7 s | 1.15 GB |
| chr4, the slowest | 89.1 s | 1.01 GB |
| chr20 | 22.1 s | 339 MB |
| chr21, chr22 | 12.9 s, 12.5 s | |
| the 2,556 contigs that are not chromosomes | 3.2 s together | |
| GRCh38 in all | 1,104 s | 1.86 GB for the process |

The peak while a contig is resolved is 4.6 to 5.3 MB per Mb of contig. A process that scans
several contigs peaks above the largest one alone, since CPython's allocator keeps the pages the
earlier ones freed: chr1 then chr2 peaked at 1.65 GB. On four threads, one `segment_contig` per
thread on a 4 Mb contig took 1.38 to 1.42 s each, against 1.18 s on one thread.

The finder on `main` today, at aaaf31b, is slower: 64 Mb of chr1 took about 106 s, and chr1 whole
is forecast at 10 to 12 minutes by the same review. The figures below are for the faster finder;
with `main`'s, a chromosome takes five to nine times longer (64 Mb, 106 s against chr20's 22 s; chr1, 10 to 12 minutes against 78 s), and the conclusions hold.

Reading chr20 with `parse_vcf` takes 0.7 s with one sample and 13 s with 200
(`doc/reports/perf-vcf-parser-2026-10-19.md`). For a run, these mean:

- **One contig, the usual benchmark on chr20.** The finder's scan of chr20, 22 s, is the wait.
  The readers fill their buffers meanwhile. With one sample, 10,000 records are less than a second
  of reading, so the threads gain under a second. What saves time on such a run is that chr1 to
  chr19 are not scanned: 984 s, about 16 minutes, by the table's run. A buffer large enough to
  hold all of chr20 with 200 samples would hide 13 s and cost about 1.5 GB per file, which this
  spec does not take.
- **Many contigs, with indexed VCFs.** Two finder threads scan two contigs at once while the loci
  side works on an earlier one. The whole of GRCh38, 1,104 s on one thread, would take a little
  more than half that on two, since each thread runs a little slower beside another and the last
  contigs leave a thread idle; this is not measured. The cost is a second scan's memory: two of
  the largest chromosomes at once, about 2.2 GB, more than that once the allocator keeps freed
  pages. Without an index the scans run nearly one after another: a reader waits on its full
  buffer, 10,000 records into a contig being scanned, and reaches the next contig only when the
  loci side has used those records (point 7).
- **With `--str-regions`**, nothing is scanned, and only the readers run beside the loci side.

## How it is verified

Each check is made at the highest function at which it can be seen, with its expected result as
a literal. The tests run on small files written to `tmp_path`: a FASTA of four contigs, of a few
hundred bases each, with STR regions the finder admits on three of them, and VCFs written as text,
some of them gzipped by the test with `gzip`, and indexes for some of them.

**The rows are the stages' rows, at `compare()`.** On the four-contig case, the rows of
`compare()` with `finder_threads` 0, 1, 2 and 4 are equal to each other and to those of the nine
steps of "The stages, joined", written out in the test, one after another. Each is run with:
file A with records on contigs 1 and 3 and file B on 2 and 3; both indexed, one indexed, neither;
with `--regions` naming contig 3 alone; and with `--str-regions`. The expected rows of one case,
a locus on contig 3 that an STR region widens, are written as literal `OutcomeRow` values.

**Which contigs are scanned, at `compare()`.** From the INFO lines "finding STR regions", which the
test reads from the log:

- with neither file indexed, the contigs with records, and not contig 4, which has none;
- with `--regions` naming contig 3, contig 3 alone, whatever the indexes;
- with both indexed and an index that lists contig 4 but no record on it, the rows are the same
  as without the index (whether contig 4's scan starts before the rows end depends on the
  threads, so it is not checked);
- with an index that leaves out contig 2, contig 2 is still scanned and the rows are the same;
- with a record past the last STR region of contig 1, one whose span ends on contig 1's last
  base, and the next record on contig 3, the rows are the stages' rows. These are the cases of
  the GENERIC segment at the first base of the next contig.

**The index reader**, on a `.tbi` and a `.csi` written by bcftools 1.24 (`bcftools index -t` and
`bcftools index -c` on a small VCF in `tests/data/`, with the commands in its `README.md`): the
contigs listed are a literal list. An index cut short gives the warning of "Which contigs are
scanned", with its text.

**Errors through the threads, at `compare()`:** a VCF out of order on contig 3; a byte of the
FASTA that is not a letter on contig 2; and a `--strata pvc` with an untyped `--str-regions`
region. With threads and with `finder_threads=0`, the same exception type and message, after the
same rows. For the FASTA error to come from the finder and not from the reader, contig 2 is
longer than 65,536 bases and the bad byte more than 65,536 bases from any record of it, since
the reader's `FastaWindow` reads at least that many bases at a time (`BLOCK_BASES`). A FASTA error
on a contig no record is on is not raised, with an index that lists that contig, so that its scan
is planned and may run.

**The bound and the stop, at `compare()`.** With the finder held on a `threading.Event` the test
sets, through a stand-in for `segment_contig` that the test puts in `pvc_eval.compare`, a reader
of a file of 50,000 records has taken from `parse_vcf` at least 10,000 of them and at most 10,000
plus two batches: the buffer full, one batch the loci side took before it waited on the finder,
since `build_loci` takes the first record of each file before it reads a segment, and one the
reader holds while it waits. After the first row is taken and the iterator closed, the test sets
the event: within 5 seconds no reader thread is alive, and the stand-in has been called no more
times than before the close.

**The table, at `format_table_row`**, on `OutcomeRow` values built in the test, each field as a
literal: 1-based coordinates; `None` as `.`; a `MISSING` haplotype; a haploid genotype; an empty
locus allele as `-`; an STR segment with its motif, an STR bundle, a satellite and a region with no
type; a QUAL of 50.0, 29.53 and 1e-05; each enum's value.

**The files, at `write_table`:** the header and the rows of a short list, with and without
`all_rows`; a `.gz` path, read back with `gzip`; rows whose iterator raises after two rows, which
leave no output, no `.tmp`, and an older file at the path unchanged; each refusal of "Writing the
files" with its message; and the summary of a list of rows over two sample pairs and two
strata, with its lines, their order and their counts as literals, a 0 among them, and
*reference in both* counted without `all_rows`.

**The command, at `main()`:** a run on the four-contig case writes the table that `write_table`
writes for `compare()`'s rows; each option's value out of range, with its message and status 1; a
`PvcEvalInputError` from a stage, its message on stderr, no traceback, status 1; the INFO lines on
stderr, matched by pattern since they hold a time.

**A file against itself over real data**, the check of `doc/spec.md` §7.3: GIAB's HG002 v4.2.1
benchmark VCF for GRCh38 against itself, over chr20, every row *same genotype*, *reference in
both*, *missing in both*, *invalid in both* or *too complex*, the last three only where the file
itself has a missing haplotype, an invalid sample or a sample over the cap. It reads GRCh38 from
the path in `PVC_EVAL_GRCH38` and the VCF from `PVC_EVAL_HG002_VCF`, and is skipped, with a
message naming the variable, without either (point 6). It is marked slow, since the finder's
scan of chr20 takes 22 s, and about 106 s with the finder on `main`.

## Changes to `doc/spec.md` and `doc/specs/loci.md`

After the owner approves this spec, each in a commit of its own:

- `doc/spec.md` §3: a pointer to this spec for `compare()`, the command, the table writer and the
  summary; `--evaluator` and the `evaluator` parameter out until a wrapper of §8 exists;
  `finder_threads` in `compare()`'s signature.
- `doc/spec.md` §6: a pointer to "The table" and "The summary" of this spec, and the summary's
  wide form.
- `doc/spec.md` §4 and `doc/specs/loci.md`: `RegionType` and `Kind` are `StrEnum`, with the code of
  `src/pvc_eval/segments.py` and `src/pvc_eval/kind.py` in the same commit. This changes no output.

## Decided by the owner

The owner answered the seven open points of the draft on 20 October 2026:

1. **`--evaluator`** and the `evaluator` parameter are left out until a wrapper of `doc/spec.md`
   §8 exists. Adding them then breaks nobody, since their default is `pvc`.
2. **The summary** is a wide table, one line per outcome and stratum and one column per kind,
   counted over all sample pairs ("The summary"), to be tried first. The long table, one line per
   sample pair, outcome, kind and stratum, was the other option.
3. **The table**: a header line with no `#`; QUAL as the shortest decimal that reads back as the
   same number, not as the VCF wrote it, which would change `pvc_eval.vcf`; plain gzip for a path
   ending in `.gz`, not bgzip; and the files written to `.tmp` and renamed.
4. **`Kind` and `RegionType` become `StrEnum`**, like `Outcome` and `InvalidReason`, so that
   formatting one gives its word. Whichever of this branch and `perf/str-finder` is merged second
   resolves a one-line conflict in `src/pvc_eval/segments.py`.
5. **The messages** are those of "The `compare` command".
6. **The real data** of the end-to-end check is GIAB's HG002 v4.2.1 benchmark VCF for GRCh38,
   at `/Users/jose/devel/pop_var_caller/benchmarks/giab/all_bench_regions/vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`,
   with GRCh38's analysis set without alts at
   `/Users/jose/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`,
   read through `PVC_EVAL_HG002_VCF` and `PVC_EVAL_GRCH38`, over chr20.
7. **The indexes** are read, with `gzip` and `struct`, so that the contigs they list are planned.

## Not in this spec

- The other evaluators of `doc/spec.md` §8, and the `evaluator` parameter (point 1).
- The re-spellings of `doc/spec.md` §7.4: the next part, once `compare` runs.
- pvc's hard cases of `doc/spec.md` §7.5 (its Open 3).
- The STR finder's speed: the performance session on `perf/str-finder`.
- Running loci on threads, and how many sample-parsing threads each reader uses: the stages keep
  their defaults; both are for the performance work, once a whole run is measured.
- The STR finder's parameters on `compare`: `str-regions` takes them.
