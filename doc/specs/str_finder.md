# The STR finder

*Status: approved by the owner on 19 October 2026, with Open 1 answered as recommended. No code
yet. It depends on `doc/spec.md` §2 (vocabulary), §4 (the segment contract, "What stage 1 yields")
and §5.1 (what the finder does), and ports pvc at the commit of `doc/port_provenance.md`.*

The STR finder is stage 1 of pvc_eval's pipeline. For one contig of the reference, it yields the
segmented genome: consecutive, non-overlapping segments from the first base to the last, each a
generic region, an STR segment, an STR bundle or a satellite. The loci (§5.3) are widened to its STR
regions and `--strata pvc` (§5.7) reads their types. It is a port of pvc's finder, rule for rule.
This spec gives
the rules in pvc's order, with the file and function of pvc each comes from, and the parts of the
interface that `doc/spec.md` §4 leaves open:

- `StrFinderParameters`;
- what `segment_contig` reads and holds;
- `read_str_regions`;
- the `pvc_eval str-regions` command and its BED.

It also says how the port is checked against pvc. Where pvc's documents and its code disagree, the
code is followed. `doc/specs/str_finder_tests.md` lists each of pvc's tests and what the port does
with it, and `doc/port_provenance.md` records the commit of pvc and the function each step comes
from.

pvc finds its regions by two routes. Its **scanning route**, `region_typing::partition_resident_in`,
scans a contig's bases. Its **catalogue route** reads a repeat catalogue instead, a Parquet file of
candidates that `pvc repeat-catalog` builds once per reference. The catalogue route is what `pvc dev
type-regions` runs, which is pvc's command that writes the regions as a BED. The two routes differ in
two places: the catalogue route drops the candidates near a contig end (step 2 below), and it always
scans periods 1 to 6, whatever periods it then admits (step 1). pvc_eval ports the scanning route
with both, so it gives what `dev type-regions` gives, under the conditions the **whole-port check**
states at the end of this spec, where the two are compared.

The owner approved the spec on 19 October 2026, and decided that 52 of pvc's tests are left out
(the decided point at the end). Four sentences of `doc/spec.md` §5.1 that differed from pvc's code were
amended by the owner's order the same day (the section before the end).

## Vocabulary

The terms of `doc/spec.md` §2 keep their meaning. These are added:

| term | meaning | not |
|---|---|---|
| **candidate** | a repeat the scanner reports: a stretch of the contig, its period (1 to 6) and its score. pvc's `RepeatInterval`. | "interval", "tract", "hit" |
| **detected span** | the stretch of a candidate as the scanner reports it | |
| **cut** | a candidate's span after trimming to whole motifs at both ends (pvc `minimal_trim`) | "trim span" |
| **motif** | a candidate's first *period* bases, in upper case, as they appear in the reference (not rotated) | "unit" |
| **copies** | a stretch's length divided by the period, rounded down | "copy number" |
| **copy floor** | the fewest copies a candidate of a given period needs, one value per period | |
| **purity** | the fraction of the bases of the cut that equal the motif repeated from the cut's first base | |
| **bundle distance** | candidates with fewer than this many bases between them bundle; 15 | "flank" |
| **contig-end flank** | the fewest bases a candidate needs between it and each end of the contig; 15 | |
| **repeat run** | a stretch covered by candidates that survive the copy floor, those that overlap or abut merged into one | "coverage run" |

Every position below is 0-based, and every stretch is half-open, `[start, end)`, as in `doc/spec.md`
§4. pvc gives its segments 1-based and inclusive. A pvc segment `[s, e]` is `[s - 1, e)` here.
pvc's BED rows are already 0-based and half-open.

## What the finder does

`segment_contig` reads the contig's bases once, in upper case, and goes through these steps in
order. Two steps differ from what pvc's documents lead a reader to expect. Step 2 is not in pvc's
scanning route. In step 5, pvc's `typed_regions.md` §2.4 says that bundle members skip the copy
floor, but its code applies the floor (step 3) before bundling.

**1. Find the candidates** (pvc `tandem_repeat.rs`, `find_tandem_repeats`,
`maximal_scoring_subsequences`, `is_primitive_motif`). For each period *p* from `scan_min_period`
to `scan_max_period`, 1 to 6, and only when *p* is less than the contig's length:

- Position *j* ≥ *p* scores the match reward (2) when the bases at *j* and *j − p* are both A, C, G
  or T and equal. Otherwise it scores minus the mismatch penalty (7). Any other letter, N included,
  matches nothing, not even itself.
- Every maximal scoring subsequence of these scores (Ruzzo and Tompa, 1999) running from *j0* to
  *j1* inclusive is a candidate `[j0 − p, j1 + 1)` with that score.
- The candidate is kept when it has at least the copies a candidate needs (2) and its motif is
  primitive. A motif is primitive when it is not a shorter motif repeated. For that test, any two
  letters that are not A, C, G or T count as equal.
- The kept candidates are sorted by period, then start. A candidate is then dropped when a kept
  candidate overlaps at least half its length (2 × overlap ≥ length) and has a shorter period that
  divides its period, such as 1, 2 or 3 for a period of 6. Only kept candidates drop others: a
  candidate already dropped drops nothing.

The candidates come out in that order, by period and then start. pvc's code describes the scoring
from its own index, so that its score *k* is position *k + p*; the result is the same.

The periods scanned are not the periods admitted, `min_period` to `max_period` of step 4. pvc's
catalogue is built by scanning periods 1 to 6 (`repeat_catalog/builder.rs`, `scan_contig`), and `dev
type-regions --min-period 2` narrows only what is admitted. The candidates of the other periods
still drop longer ones here, and still make repeat runs at step 7: with `--min-period 2`, a run of
120 A's is a satellite in `dev type-regions` and generic in pvc's scanning route, which scans the
admitted periods only (`region_typing/mod.rs`, `partition_resident_in`). pvc's partition tests that
run at periods 2 to 6 were written for the scanning route, and set the scanned periods to 2 to 6
too.

**2. Drop the candidates near a contig end** (pvc `repeat_catalog/row.rs`, `row_for_interval`).
This is the catalogue route's rule, which pvc's `dev type-regions` applies and its scanning route,
`region_typing::partition_resident_in`, does not. A candidate is dropped when fewer than the
contig-end flank (15) bases lie between the contig's start and its first base, or between its last
base and the contig's end. These are measured on its cut, or on its detected span when it has no
cut (step 6). A candidate dropped here is gone from every step below: it bundles nothing and counts
in no repeat run. A contig-end flank of 0 turns the step off, and gives pvc's scanning route
exactly. The ported tests of pvc's partition run that way.

**3. The copy floor on the detected span** (pvc `segment_criteria.rs`, `prefilter`). A candidate
with fewer copies in its detected span than the copy floor of its period is dropped. The candidates
left are the ones the repeat runs of step 7 are made of.

**4. The scope tests** (pvc `classify`, its loop before `split_bundles`). A candidate is dropped,
in this order, when:

- its period is outside `min_period` to `max_period`;
- its score is under `min_score` (0, which drops nothing the scanner makes);
- its detected span is shorter than its period;
- its motif is compound (`is_compound`), that is, equal to itself shifted by a length that divides
  it.

On the upper-case bases `segment_contig` reads, step 1 already drops every motif `is_compound`
would flag. The test is ported anyway, for pvc's tests of it.

**5. Bundles** (pvc `split_bundles`, `joins_cluster`). The candidates left are sorted by start, then
end, and walked in that order. A cluster's reach is the furthest end of any candidate in it. The
next candidate joins the cluster when fewer than the bundle distance bases lie between the reach
and its start; overlapping counts as no bases. Every candidate of a cluster of two or more is an
STR bundle member. A candidate alone in its cluster goes on to step 6. Bundling sees detected spans,
before any cut, and none of step 6's tests: an impure candidate, or one that has no cut, bundles its
neighbour all the same.

**6. STR segments** (pvc `finish_locus`, `minimal_trim`, `recompute_purity`). A candidate alone in
its cluster becomes an STR segment unless one of these rejects it:

- **No cut.** Let *m* be the motif's length, *L* the candidate's length, and *t* = min(3*m*,
  ⌊*L* / 2⌋). Offsets count from the candidate's first base. The cut starts at the smallest offset
  *s*, from 0 to *t*, at which the 2*m* bases from *s* are two copies of the motif. The search for
  *s* fails, before it compares, as soon as *s* + 2*m* ≥ *L*. The cut ends at the largest offset
  *e*, from *L* down to *L* − *t*, at which the 2*m* bases before *e* are two copies of the motif.
  The search for *e* fails as soon as *e* < 2*m*. There is no cut when either search fails, or when
  *s* ≥ *e*. So a candidate exactly two copies long, *L* = 2*m*, has no cut.
- **Copy floor on the cut.** The cut has fewer copies than the copy floor.
- **Purity.** The cut's purity is under the minimum purity (0.8).
- **Contig ends.** The cut starts on the contig's first base or ends on its last. After step 2,
  with a contig-end flank of 1 or more, this never fires.

The STR segment is the cut, with the candidate's motif. A rejected candidate leaves no segment of its
own. Its bases end up in a generic region, unless step 7 covers them with something else.

**7. Satellites and bundles** (pvc `region_typing/mod.rs`, `resolve_features`, `absorb_into`,
`merge_runs`, `coverage_runs`, `bundle_clusters`).

- A repeat run longer than `max_str_len` (100 bases) is a satellite. The repeat runs are
  made from the candidates of step 3, including those that step 4 drops and step 6 rejects.
- **Absorption**, repeated until a pass changes nothing. Each bundle member, by its detected span,
  and each STR segment, by its cut, is tested against the satellites. A feature that overlaps a
  satellite, or has fewer than the bundle distance bases between them, is absorbed. It is removed,
  and it and every satellite it reaches become one satellite covering all of them, gaps included.
  Satellites that then overlap or abut are merged. This is pvc's plain fixed point over
  `absorb_into`, the form `absorb_by_iteration` runs in pvc's tests. pvc runs a sweep instead,
  `absorb_by_sweep`, which its property test holds equal to the fixed point.
- The members left are grouped again by step 5's rule. A bundle's hull runs from its first
  member's start to the furthest end of any member, both detected spans.
- A hull longer than `max_bundle_len` (100 bases) is a satellite. Absorption is not run again
  for it.
- Otherwise, a hull that starts on the contig's first base or ends on its last is dropped, and its
  bases go to generic regions.
- Otherwise the hull is an STR bundle.

**8. The segmented genome** (pvc `fill_generic_gaps`). The satellites, STR segments and STR bundles,
sorted by start, never overlap. Each stretch between them, and before the first and after the last,
is one generic region. A contig of length 0 has no segments.

The whole contig is resolved before the first segment is yielded, because step 7 can grow a
satellite to either side by any length.

### Worked cases

pvc gives these, and the port must too. Each was run on 19 October 2026 with pvc built at bbcc9f4:
`pvc repeat-catalog --reference cases.fa`, then `pvc dev type-regions --reference cases.fa --min-copies
5,5,4,4,4,3`, which are pvc_eval's copy floors. Each contig is a clean (CA) repeat with a homopolymer
6 bases after it.

```
four_a  GCTAGACATACATACATACAGTCAGCACGACTGTCAGTGATCGCTAGTAGTAGTGATGCACACACACACACACACACACAGTCTGT
        AAAACTACGCTACTGCTGTCACTCGACTGCATACATCAGACAGACTCGTATGACAGTCACGCAGA
five_a  the same with AAAAA in place of AAAA
```

| contig | segments |
|---|---|
| four_a (151 bases) | generic 0–58, str 58–80 `CA`, generic 80–151 |
| five_a (152 bases) | generic 0–58, str-bundle 58–91, generic 91–152 |

The (CA) repeat starts at 58, not 60, because the filler before it ends in `CA`. Four A's are under
the period-1 copy floor of 5 and are dropped at step 3, so the (CA) repeat stands alone. Five A's
survive step 3 and bundle with it at step 5.

Step 2 at work: the same (CA) repeat with six A's 8 bases before it, the A's starting at base 10 of
the contig (`near_end`) or at base 17 (`far_end`).

```
near_end  GCTAGACATGAAAAAACTGTCAGTCACACACACACACACACACAGTCTGTAGCTACTGCTGTCACTCGACTGCATACATCAGACAG
          ACTCGTATGACAGTCACG
far_end   GCTAGACATGCATCAGTAAAAAACTGTCAGTCACACACACACACACACACAGTCTGTAGCTACTGCTGTCACTCGACTGCATACATC
          AGACAGACTCGTATGACAGTCACG
```

| contig | segments |
|---|---|
| near_end | generic 0–24, str 24–44 `CA`, generic 44–104 |
| far_end | generic 0–17, str-bundle 17–51, generic 51–111 |

In `near_end`, step 2 drops the A's, 10 bases from the start, so they bundle nothing. pvc's scanning
route alone would make a bundle of both there.

## Parameters

The defaults are `doc/spec.md` §5.1's, with one difference: §5.1's single "bundle and flank
distance" is two parameters, as it is in pvc, where the bundle distance is
`SsrSegmentCriteria::bundle_threshold` (`DEFAULT_BUNDLE_THRESHOLD`) and the contig-end flank is
`repeat_catalog/criteria.rs` `CATALOG_MIN_FLANK_BP`. pvc's `dev type-regions --flank-bp` moves the
first and not the second. Setting the contig-end flank to 0 is also how pvc's scanning route is
reproduced for its tests.

| field | default | must be |
|---|---|---|
| `scan_min_period`, `scan_max_period` | 1, 6 | 1 ≤ `scan_min_period` ≤ `scan_max_period` ≤ 6 |
| `min_period`, `max_period` | 1, 6 | 1 ≤ `min_period` ≤ `max_period` ≤ 6 |
| `match_reward`, `mismatch_penalty` | 2, 7 | 1 to 2,147,483,647 (2^31 − 1) |
| `scan_min_copies` | 2 | 1 or more |
| `copy_floors` | (5, 5, 4, 4, 4, 3), for periods 1 to 6 | six whole numbers, 0 or more |
| `min_purity` | 0.8 | a finite number from 0 to 1 |
| `min_score` | 0 | a whole number |
| `bundle_distance` | 15 | 1 or more |
| `contig_end_flank` | 15 | 0 or more |
| `max_str_len` | 100 | `bundle_distance` or more |
| `max_bundle_len` | 100 | 0 or more |

The limits on the periods, `min_purity`, `bundle_distance` and `max_str_len` are pvc's, from
`classify`'s asserts and `cli/typed_regions.rs` `walk_config`. A `bundle_distance` of 0 makes every
candidate fail step 6 and bundles nothing, so pvc refuses it. The upper limit of the reward and the
penalty is pvc's too: pvc holds them as `i32` (`tandem_repeat.rs`, `ScanParams`) and its running
totals as `i64`, and pvc_eval's step 1 holds its running totals as 64-bit integers, which at these
limits cannot overflow on a contig shorter than 2^31 bases (decided by the owner on 21 October
2026). The others are pvc_eval's own: pvc allows `scan_min_copies` of 0, a value under which the
scores no longer find repeats.

## The Python interface

`doc/spec.md` §4 fixes `RegionType`, `GenomeSegment`, `segment_contig` and `read_str_regions` in
`pvc_eval.segments`. This adds what it leaves open.

The parameters, a record checked when it is made. A value outside the table above raises
`ValueError` naming the field, the value and its limit.

```python
@dataclass(frozen=True, slots=True, kw_only=True)
class StrFinderParameters:
    scan_min_period: int = 1
    scan_max_period: int = 6
    min_period: int = 1
    max_period: int = 6
    match_reward: int = 2
    mismatch_penalty: int = 7
    scan_min_copies: int = 2
    copy_floors: tuple[int, int, int, int, int, int] = (5, 5, 4, 4, 4, 3)
    min_purity: float = 0.8
    min_score: int = 0
    bundle_distance: int = 15
    contig_end_flank: int = 15
    max_str_len: int = 100
    max_bundle_len: int = 100
```

**`segment_contig(reference, contig, parameters)`** reads the contig once, whole, with
`reference.get_bases(contig, 0, reference.get_length(contig))`. It needs the whole contig because
step 1 runs over all of it and step 7 can reach any distance. A `FastaWindow` would read the same
bases in blocks, and would buy nothing. It raises what `IndexedFasta` raises: `ValueError` for a
contig the `.fai` lacks, and `PvcEvalInputError` for a byte of the FASTA that is not a letter. It
keeps no state between calls and changes nothing it shares, so threads can segment different
contigs of one `IndexedFasta` at once.

What it holds, for the largest human contig (chr1 of GRCh38, 248,956,422 bases):

- The bases, as one `str` of one byte per base, about 250 MB. While `get_bases` reads them it holds
  two copies (`fasta.py`, `_read_upper_bytes`), so about 500 MB at that moment.
- The candidates that survive step 3. pvc's catalogue of GRCh38 holds about 1.8 million of them on
  chr1, at the same copy floors, and 23.6 million on the whole reference. This was counted on the
  catalogue built on 2 September 2026, beside the reference in pvc's `benchmarks/giab`.
- **Not all the candidates of step 1.** Every stretch of two or more equal bases, A, C, G or T, is
  one candidate of period 1. chr1 has 45.4 million of them, counted on its bases on 19 October 2026,
  and the candidates of periods 2 to 6 add more. So step 1 must not hold them all, and must still drop exactly what pvc
  drops. Step 3 drops most of them, but a candidate under the copy floor is used at step 1 first: it
  can drop a candidate of a longer period that it covers by half. So step 1 keeps each candidate
  until every longer-period candidate it could drop has been found, and no longer. How to arrange
  that is the plan's.

How long a contig takes, and how much memory at its peak, are measured after the code exists, with
the `performance` skill. That measurement decides whether the finder needs numpy, which is a runtime
dependency and the owner's decision (`doc/spec.md` §5.1).

**`read_str_regions(path, reference)`** reads a BED of STR regions: the one `str-regions` writes,
or a user's. It reads the file as UTF-8, a byte-order mark at its start ignored. It skips lines
that hold only spaces or tabs, lines that start with `#`, and lines whose first field is `track` or
`browser` (a contig may be named `trackA`). Empty columns at the end of a line are no columns. Each
other line has the contig, start and end, 0-based and half-open, and may have two more columns. The
fourth is the region type, `str`, `str-bundle` or `satellite`. The fifth is the motif, which is
upper-case A, C, G and T for `str` and `.` for the others. It yields one `GenomeSegment` per line.
A line of three columns gives `region_type` `None` and `motif` `None`, as `doc/spec.md` §4 says. A
line of four columns with the type `str` has no motif, and is refused. A motif is 1 to 6 bases.
Columns after the fifth are ignored. pvc's own BED has nine columns and other words for the types,
`ssr_locus`, `ssr_bundle` and `generic`, so it is refused at its first line, as an unknown type.

The BED must be sorted as the reference is: contigs in the order of the `.fai`, and by start within
each. `read_str_regions` reads it as a stream and does not sort it, because a BED of every STR
region of GRCh38 holds millions of lines.

It raises `PvcEvalInputError`, naming the file and the line and saying how to put it right, for:

- a line with fewer than three columns;
- a start or end that is not a whole number, or an end not after its start;
- a contig the reference lacks, or an end past the contig's length;
- an unknown region type, or a motif that does not fit its type;
- a contig whose lines are not all together, or that comes before a contig preceding it in the
  `.fai`;
- a line that starts before the end of the line before it on its contig.

These follow `doc/spec.md` §4, which refuses overlapping regions. `str-regions` never writes
overlapping regions. A user's BED of repeats that overlap can be merged first with `bedtools merge`,
which keeps three columns.

**`write_str_regions(reference, output, parameters, regions, num_threads)`** does what the
`str-regions` command does, below: it segments the contigs, keeps the STR regions that the requested regions overlap, or all
of them when `regions` is `None` (no `--regions`), and writes the BED. The command calls it, and it is public because the
filter of STR regions by requested regions is private to `pvc_eval.segments` (added on 19 October
2026, during the build).

## The `str-regions` command

```
pvc_eval str-regions --reference ref.fa --output str.bed [--regions regions.bed]
                     [--scan-min-period N] [--scan-max-period N]
                     [--min-period N] [--max-period N] [--min-copies N,N,N,N,N,N]
                     [--min-purity X] [--min-score N] [--bundle-distance N]
                     [--contig-end-flank N] [--max-str-len N] [--max-str-bundle-len N]
                     [--scan-match-reward N] [--scan-mismatch-penalty N] [--scan-min-copies N]
                     [--threads N]
```

Each option sets the field of `StrFinderParameters` with the same meaning, with the same default.
Where pvc's `dev type-regions` has the option, the name is pvc's, so that the two commands of the
whole-port check read alike. The exception is `--bundle-distance`, which pvc calls `--flank-bp`;
in pvc that option moves the bundle distance only. pvc has no option for the contig-end flank or for
the periods scanned. `dev type-regions` accepts `--scan-match-reward`, `--scan-mismatch-penalty` and
`--scan-min-copies` and ignores them, because its catalogue was scanned when it was built
(`repeat_catalog/criteria.rs`, `StrRepeatCriteria::from`); in pvc_eval they act.

It segments every contig of the `.fai`, in its order. With `--regions`, it segments only the contigs
the regions name, each still whole, so a region list changes what is written and never what a
segment is. It writes the STR regions, not the generic ones. With `--regions`,
it writes only those that overlap a region by at least one base, each whole, as pvc's `emit_into`
does.

**Threads** (added on 20 October 2026, by the owner's order, after the performance review,
`doc/reports/perf-str-finder.md`). `--threads N` segments up to N contigs at once, each on a
thread of its own; N is a whole number from 1 to 256, 4 when not given, the same default and
limit as the VCF reader's. With 1 there is no pool: the contigs are segmented one after another
in the calling thread, as before. The BED is the same at any N, byte for byte, and its header
does not name N: a contig's rows depend on that contig and its own requested regions only, and
they are written in the order of the `.fai`, each contig's once the contigs before it are
written. The contigs go to a `ThreadPoolExecutor` through `Executor.map` with
`buffersize=2 * N` (Python 3.14 and later), which hands back the results in the order of the
contigs and keeps at most that many submitted ahead of the one being written; so at most N
contigs are being segmented at a time, and the others submitted wait finished. A thread hands back
the formatted rows of its contig, not its segments, so a finished contig that waits holds its
rows only. A long contig early in the `.fai` holds back the writing of those after it, and with
them the submission of more, so the pool can sit partly idle until it is written; on GRCh38, whose
`.fai` lists the chromosomes first, longest first, this costs little.

Each contig being segmented holds 4.6 to 5.3 MB per Mb of its length at its peak (measured on
GRCh38 chromosomes, Apple M5 Pro, CPython 3.15t; chr1: 1.15 GB), so the peak memory grows with
N; the command's `--help` says so. N is checked before the `--regions` BED is read or the `.tmp`
made. An error of the input in a contig stops the run with its message and removes the `.tmp`, as
with one thread, and it is the same error one thread would give, since earlier contigs' errors
come first in order; but it is raised only once the contigs before it are written, and after the
contigs already started have finished. Ctrl-C stops the run once the contigs being segmented at
that moment have finished, which on human chromosomes can take minutes (chr1 alone takes about
75 s on one thread, longer when threads compete), and removes the `.tmp`;
a second Ctrl-C removes it at once, and the process then waits for those threads before it exits.

The BED:

```
##pvc_eval str-regions 0.1.0
##scan_min_period: 1
…                                  one ##field: value line for each field of StrFinderParameters,
                                   in the order of its fields
#chrom	start	end	type	motif
chr1	10432	10470	str	CA
chr1	10519	10601	str-bundle	.
chr1	12011	14002	satellite	.
```

The rows are in the order of the `.fai` and then by start, 0-based and half-open. The `##` lines
let a user see the parameters a cached BED was made with. It is written under `<output>.tmp` and
renamed when complete, so that `compare --str-regions` never reads half a file.

Errors, as `doc/spec.md` §3 and the `writing` skill say, through `PvcEvalInputError`, with exit
status 1:

- the reference or its `.fai` missing, which `IndexedFasta` reports;
- an option value outside the table of the parameters, or a `--threads` outside 1 to 256;
- what `pvc_eval.bed.parse_bed` refuses in a `--regions` BED: a file that is missing or cannot be
  read, a BED with no region, a line with a contig the `.fai` lacks, a start or end that is not a
  whole number, an end not after its start; and a region whose end is past the contig's length;
- an `--output` that cannot be written: a directory, a path through a file, a directory without
  permission. These are checked before any contig is segmented, so that the user does not wait
  through a run to learn it.

A command line that cannot be parsed, an option missing or a value that is not a number, gets
argparse's usage message and exit status 2, as for any program built on it.

A `--regions` BED that holds no region stops the run, as it does for `compare --regions`.

The `--regions` BED is read by `pvc_eval.bed.parse_bed`, the reader `compare --regions` uses, so a
BED means the same to both commands, gzipped or not; `str-regions` adds the test of each end against
its contig's length (changed on 19 October 2026, when the loci brought that reader into main).

## How it is verified

**pvc's tests, rewritten.** pvc has 142 tests of these functions at bbcc9f4: 35 in
`tandem_repeat.rs`, 58 in `segment_criteria.rs` and 49 in `region_typing/mod.rs`, two of those 49
property tests. 90 of them test what is ported, and two more come from `repeat_catalog/row.rs`
for step 2. The other 52 test what is not ported:

- the windowed scan and the region scanner, pvc's other ways of running step 1, 25 tests;
- pvc's BED reading and its region sets (`RegionSet`, `GenomeRegions`), one of them a property
  test, 9;
- Rust types: the motif type and the segment record's own checks, 6;
- the counts pvc keeps of the bases each test turned down, 4;
- the comparison with pvc's old production code, 3;
- the function that turns satellites back into generic regions for pvc's `--call-satellites`,
  `hand_satellites_to_the_generic_path`, 3;
- pvc's own default settings, 1;
- the property test of the absorption sweep against the fixed point, 1, which compares a function
  with itself once only the fixed point is ported.

(**Open 1**.) `doc/specs/str_finder_tests.md` lists every test, with what the Python test changes
and why. Each rewritten test keeps pvc's inputs, copies its inline sequences as they are, and
keeps its expected outputs, moved to 0-based half-open where pvc's are 1-based. It passes explicitly
every parameter pvc took from a `Default`. These tests call the steps above, which are private
functions of `pvc_eval.segments`, as `doc/spec.md` allows a test to do where the rule is clearer
there. A test of pvc's partition uses `contig_end_flank=0`. The steps that pvc tests through their
own output keep that output in the port: the reason step 4 or step 6 rejected a candidate, and the
members of a bundle. `GenomeSegment` does not carry them.

**pvc's fixtures.** `tests/data/tandem_repeat/synthetic_ref.fa`, with its `.fai`, and
`golden.ssr_catalog.bed.gz` are copied into `tests/data/tandem_repeat/`, with a line each in
`tests/data/README.md`. pvc's `the_resident_partition_reproduces_the_golden_catalog` compares the
partition of `synthetic_ref.fa` with the 16 STRs of the golden catalogue, by overlap. It runs with
periods 2 to 6, copy floors (10, 5, 4, 3, 3, 3), bundle distance 50, `max_str_len` 1000 and
`max_bundle_len` 100.

**pvc_eval's own cases**, at `segment_contig`:

- the four worked cases above, with their segments as literals;
- the segments of one contig are consecutive and non-overlapping from 0 to its length. This holds
  on the fixture and, with hypothesis, on generated contigs;
- the copy floors of pvc_eval's defaults at their boundary: a homopolymer of 4 and of 5 bases, and a
  (CA) repeat of 4 and of 5 copies, each 30 bases from the ends;
- an empty contig;
- a contig shorter than the period.

**The `str-regions` command and `read_str_regions`**, on a small FASTA in `tmp_path`:

- the BED of the worked cases, read back to the same STR regions;
- `--regions`: a region that cuts an STR region still writes it whole, and a contig the regions do
  not name is not segmented;
- each error above, with its message;
- each refusal of `read_str_regions`;
- `--threads`: a FASTA of more contigs than threads, with STR regions in each, whose first contig
  is by far the longest, so that the others finish before it, gives the same BED with 1 thread and
  with 3; `--threads 0` and `257` stop with exit status 1 and a message naming `--threads`, and a
  word gets argparse's usage message and status 2, and 1 and 256 are accepted; an error of the
  input in a contig that is not the first, under several threads, stops the run with its message
  and leaves no `.tmp`, the message of the first faulty contig in the `.fai` even when a later
  one fails first; with `--threads 1` every contig is segmented in the calling thread, and with
  N > 1 none is, on at most N threads; and with N threads at most 1 + 2N contigs have started
  when the first is written (`Executor.map` submits one more before it waits on the first).
  These are watched with `sys.monitoring` on `segment_contig`, which records its calls and
  replaces nothing.

**The whole-port check.** pvc_eval's STR regions and `pvc dev type-regions`'s are identical on the
same sequences: same contig, start, end and type, and the same motif for each STR segment. During
development it runs on small sequences, so that it takes minutes and can be run after every change
to a step. It is a module of tests marked `slow`, over two inputs:

- **Generated contigs.** 600 contigs, from a fixed seed, of random bases with planted repeats of
  periods 1 to 6, some impure, some near each other, some near a contig end, and some runs longer
  than 100 bases, and repeats near both of a contig's ends, so that step 2 is compared at each.
  Both inputs are compared twice: at the defaults below, and with `--min-period 2 --max-period 5`,
  which tells periods scanned from periods admitted. The review of this spec ran
  the spec's steps, written in Python, against pvc on 4,800 such contigs with no difference.
- **Slices of GRCh38.** A few stretches of 1 to 2 Mb, each written as a contig of its own into a
  FASTA in `tmp_path`, chosen so that together they hold STR segments, STR bundles and satellites;
  the test asserts that each type appears. Their ends are the ends of the new contigs for both
  programs alike.

For each input, the test writes the FASTA and its `.fai`, builds pvc's catalogue beside it with `pvc
repeat-catalog`, which takes seconds at this size, and runs, with no `--regions`:

```
pvc dev type-regions --reference input.fa --output pvc.bed [--min-period 2 --max-period 5]
pvc_eval str-regions --reference input.fa --min-copies 8,6,6,6,5,4 --output pvc_eval.bed \
                     [--min-period 2 --max-period 5]
```

It reads two paths from the environment: `PVC_EVAL_PVC`, a pvc binary built at the provenance
commit, and `PVC_EVAL_GRCH38`, the GRCh38 FASTA with its `.fai`. Where one is unset or its file is
missing, the part that needs it is skipped, with a message naming what is missing. The pvc binary in
pvc's `target/release` lacks `pvc dev`, so pvc is built from `git archive bbcc9f4` into a scratch
directory, and pvc's repository is not touched.

`dev type-regions` at its defaults has pvc's calling floors, 8, 6, 6, 6, 5, 4, the copy floors pvc
uses when it calls genotypes (`doc/spec.md` §5.1), and pvc_eval's other defaults. It writes
`ssr_locus`, `ssr_bundle` and `satellite`, which are `str`, `str-bundle` and `satellite` here, and
`generic` rows, which are not compared.

pvc_eval equals `dev type-regions` only where both run the same rules: periods scanned 1 to 6, copy
floors at or above the catalogue's 5, 5, 4, 4, 4, 3 (which pvc refuses to go under, `StrRepeatCriteria::
serves`), the scan weights 2, 7, 2 the catalogue was built with, and a contig-end flank of 15. The
check runs there.

**The same check over GRCh38**, over the GIAB benchmark regions, is run once the code exists and its
performance work is done, and is not part of the plan made from this spec (the owner, 19 October
2026). Two facts for it, found by the review of this spec:

- `dev type-regions` is run over the whole genome and its rows that overlap a benchmark region are
  kept. With `--regions`, its catalogue reader loads only the rows within a reach of each region
  (`repeat_catalog/reader.rs`, `rows_for`), and a satellite longer than that reach comes back
  shorter: a 780-base chain of (CA) repeats was a satellite 2000–2770 over the whole contig and
  2567–2770 with a region at 2700–2750. Over HG002's regions the two runs agreed on all 1,530,189
  STR regions, by luck of the data.
- The catalogue in pvc's `benchmarks/giab` was built on 2 September 2026, before bbcc9f4, so it is
  built again beside links to the FASTA and its `.fai`. The review built it with `--threads 4` in
  43 s, at a peak of 22 GB, on this machine on 19 October 2026.

## Decided point

**Open 1, decided. pvc's tests that test what is not ported.** The owner chose to leave them out on
19 October 2026. The brief for this port asks for every one of pvc's 142
tests to be rewritten. 52 test functions the port does not have, and the list is in
`doc/specs/str_finder_tests.md`.

- *Leave them out*, and say so there with the reason for each. The port keeps only what pvc_eval
  runs.
- *Port the functions they test, to port the tests.* That means the windowed scan, the region
  scanner, the counts, pvc's BED normaliser and the absorption sweep, several hundred lines of code
  that pvc_eval would never call.

I recommend leaving them out. Each of the 52 tests either has a subject that pvc_eval does not run,
or, like the sweep property, compares a function with itself once only the fixed point is ported.

## Sentences of `doc/spec.md` amended

Four sentences of `doc/spec.md` §5.1, and the definition of an STR bundle in §2, differed from pvc's
code: the contig-end rule placed among the tests of an STR segment, where pvc drops the candidate
before the copy floor and bundling; "within 15 bases" where pvc's code means "with fewer than 15
bases between them"; and pvc's two routes said to give the same partition, with or without a region
list. The owner had them amended to follow the code on 19 October 2026.

## Not in this spec

- **Threads inside one contig**, and any faster form of the finder beyond
  `doc/reports/perf-str-finder.md`. `str-regions` segments contigs on threads (above); one contig
  is segmented by one thread.
- **How `compare` uses the segments** without `--str-regions`, and which contigs it segments. That
  is `doc/spec.md` §4 and the loci's spec.
- **pvc's counts** of the bases each test turned down. pvc_eval reports segments, not counts.
- **Bundle members** outside the private steps. The BED and `GenomeSegment` do not carry them.
