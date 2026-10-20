# Reading only the contigs `--regions` names

*Status: written on 20 October 2026, after the owner agreed that `--regions` should limit what is
read, and asked for the spec to be reviewed and built without waiting for their approval. No code
of this change exists yet. It changes `parse_vcf` (`doc/spec.md` §5.2) and `compare()`
(`doc/specs/compare.md`).*

With `--regions`, `pvc_eval compare` reads the VCF records of the contigs the regions name, and
of no other contig. Today it reads every record of both files, builds every locus and then skips
those `--regions` leaves out. On GIAB's HG002 v4.2.1 against itself over chr20 of GRCh38, with
`--regions` naming chr20, the rows of chr20 were done at 38 s of a 39 s run: the other 21
chromosomes took 35 s and gave no row (`doc/reports/compare.md`). The table does not change, since
a locus never spans two contigs and every locus of a contig `--regions` does not name is skipped
already.

## What is read, and what is still checked

`parse_vcf` takes the contigs to give, `contigs`. A line of a VCF record on another contig is
skipped after three checks, made on its first two columns alone, without splitting the rest of the
line:

- its contig is in the reference;
- its POS is a whole number;
- it does not come before the line before it, in the `.fai`'s contig order and by position.

So an unsorted file, or one with a contig the reference lacks, still stops the run with the same
message and line number as today. What is no longer checked on a contig `--regions` does not name:
the REF against the reference, the REF within the contig's length, the number of sample columns,
the ALT and the GT. A file made on another build of the reference, whose contigs have the same
names and lengths, is still stopped by its first record on a named contig, since that record's REF
is checked as today.

The warnings at the end of a file count only the records given: a symbolic allele on a contig
`--regions` does not name is not counted, and nor is a QUAL that is not a number there.

Line numbers in messages stay exact, since every line is still read in order. Skipping by the
first two columns costs little: Python reads the 4,048,572 lines of HG002's VCF, testing whether
each starts with `chr20` and a tab, in 2.2 s (measured on 20 October 2026, Apple M5 Pro, CPython
3.15t, `gzip` text reading). A run over chr20 then waits on the STR finder's scan of chr20, 22 s
(`doc/specs/compare.md`, "How it runs"), or on the reading and the loci of chr20, a few seconds,
with `--str-regions`.

This works the same for a VCF with or without an index, plain or bgzipped. The owner's two answers
of 20 October 2026 assumed a reader that seeks through the index to each named contig: whole contigs
read, and a VCF with no index read whole, with a message. Skipping by the first two columns gives
the same rows with no reader of compressed blocks, keeps the line numbers that seeking would lose,
and gives a VCF with no index the same gain, so neither a message nor a refusal is needed. The
index reader of `doc/specs/compare.md` still plans the STR finder's scans.

## The Python interface

`parse_vcf` gains one keyword argument. `None`, the default, gives every record, as today.

```python
def parse_vcf(path: Path, reference: IndexedFasta, *,
              span_limit: int | None = DEFAULT_SPAN_LIMIT,
              num_threads: int = DEFAULT_NUM_THREADS,
              ploidy: int = DEFAULT_PLOIDY,
              contigs: Collection[str] | None = None) -> Iterator[VcfRecord]: ...
```

`contigs` holds contig names; a name the reference lacks gives no record and no error, since the
regions are checked against the reference when they are read. `Regions` gains a property,
`contigs`, the `frozenset` of the contigs it holds a region on, and `compare()` passes it to both
readers when `CompareConfig.regions` is given.

## What `compare()` and the command show

- A contig `--regions` does not name gets no locus, so no progress line ("chr1: 0 loci compared,
  304,084 skipped by --regions" goes), and the count of loci skipped by `--regions` on a named
  contig is unchanged.
- The rows, the summary and the first error on a named contig are as today.
- `doc/specs/compare.md` says, in "The `compare` command", that loci on a contig `--regions` does
  not name are built without STR regions and counted as skipped; that sentence goes, since they are
  not built.

## How it is verified

- **At `parse_vcf`:** with `contigs` naming chr2 of a three-contig file, only chr2's records come,
  in order. A line on chr1 or chr3 with a REF that differs from the reference, too few sample
  columns, a malformed GT or a symbolic ALT gives no error and no warning. A line on chr3 that comes
  before the line before it, a line on a contig the reference lacks, and a POS that is not a number
  each give today's error, with its line number. A name in `contigs` the reference lacks gives no
  error. `contigs=None` gives today's records.
- **At `compare()`:** the rows with `--regions` equal those of today's code on the four-contig case
  of `tests/compare_case.py`, with `--regions` cut inside contig 3 and naming contig 3 alone, at 0
  and 2 finder threads; a REF that differs from the reference on contig 1, with `--regions` naming
  contig 3, does not stop the run, and without `--regions` it does.
- **At `Regions.contigs`:** the contigs of a BED with lines on two contigs.
- **On real data:** the run of HG002 against itself over chr20 gives the same 81,257 rows and
  summary as before, and its time and peak memory go into the report.

## Not in this spec

- Reading only the stretches `--regions` covers within a contig: a locus at a region's edge can
  need records further than any fixed margin (`doc/specs/compare.md`, the owner's answer of 20
  October 2026).
- `pvc_eval str-regions`, which reads no VCF.
