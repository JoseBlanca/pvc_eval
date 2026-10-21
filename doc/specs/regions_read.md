# Reading only the contigs `--regions` names

*Status: written, reviewed and built on 20 October 2026, on the branch `regions-read`; not approved by the owner, who asked for it to be built without waiting. The owner agreed that day that `--regions`
should limit what is read, and asked for the spec to be reviewed and built without waiting for
their approval. No code of this change exists yet. It changes `parse_vcf` (`doc/spec.md` §5.2)
and `compare()` (`doc/specs/compare.md`).*

With `--regions`, `pvc_eval compare` reads the VCF records of the contigs the regions name, and
skips the lines of every other contig after a few checks. Today it reads every record of both
files, builds every locus and then skips the loci `--regions` leaves out. On GIAB's HG002 v4.2.1
against itself over chr20 of GRCh38, with `--regions` naming chr20, the run took 39 s, and the
reading and loci of the 21 other chromosomes, which give no row, took 35 s of them; chr20's last
rows came at 38 s (`doc/reports/compare.md`, Apple M5 Pro, CPython 3.15t). The table does not
change: a locus never spans two contigs, and every locus of a contig `--regions` does not name is
skipped already.

## The contigs read

A contig is read when `--regions` names it: when a region of the BED shares a position with the
contig, which is how `compare()` already decides which contigs the STR finder scans
(`doc/specs/compare.md`, "Which contigs are scanned"). A BED line on a contig that lies wholly past
the contig's end, which the BED reader does not refuse, names no contig.

## What a skipped line is still checked for

A line of a VCF record on a contig that is not read is checked for these, in this order, which is
today's, and with today's messages and line numbers:

1. it has at least the 8 fixed columns of a VCF record;
2. its POS is a whole number written in digits, and not too long to read; 0 is accepted, as it is
   today for a record with a symbolic allele;
3. its contig is in the reference;
4. it does not come before the record line before it in the file, of any contig, in the `.fai`'s
   contig order and then by POS.

So a file out of order, or with a contig the reference lacks, stops the run as today. The rest of
the line is not split. These are no longer checked on a contig that is not read: the number of
sample columns, a POS of 0, the REF against the reference, the REF within the contig's length, the
ALT, the QUAL and the GT. A record with a symbolic allele there is not counted in the warning that
counts them, and nor is a QUAL or GQ that is not a number.

**The REF error on a contig read.** When a REF differs from the reference, the message says how
many records before it agreed with the reference, which tells a VCF made on another reference
(none) from a damaged line (many) (`doc/spec.md` §5.2). The count now covers the records of the
contigs read, since the others are not compared with the reference. For the first record read,
the message says "It is the first record read with no symbolic allele" in place of "the file's
first record". So a damaged first line of chr20, in a file whose other chromosomes agree with the
reference, is reported as from another reference: the file's other contigs, not read, cannot say
otherwise. `doc/spec.md` §5.2 changes to match, in a commit of its own.

**What it costs.** Every line is still read in order, which keeps the line numbers of the messages
exact. Reading HG002's 4,048,572 lines through pvc_eval's line reader and making the four checks on
each takes about 4 s (4.05 s, measured by the spec's reviewer on 20 October 2026 on the same
machine, with the message's location written only when an error is raised). A run over chr20 of a
whole-genome VCF then waits on the STR finder's scan of chr20, 22 s with an index and about 25 s
without, since the scan then starts when the reading reaches chr20, about 3.5 s in; or, with
`--str-regions`, on the reading and the loci of chr20, a few seconds. That is an estimate: the
change is measured on HG002 once it is built.

**Why not seek through the index.** The owner was asked two things: whether to read whole contigs
or only the stretches `--regions` covers, and what to do with a VCF that has no index. The owner
chose whole contigs, read by seeking through the VCF's index, and a VCF with no index read whole
with a message. Skipping by the first columns reads whole contigs too, gives the same rows, needs
no reader of compressed blocks, keeps the line numbers that seeking would lose, and gives a VCF with
no index, plain or bgzipped, the same gain, so neither a message nor a refusal is needed. It costs
the 4 s of reading the skipped lines, which seeking would not. The index reader of
`doc/specs/compare.md` still plans the STR finder's scans.

## The Python interface

`parse_vcf` gains one keyword argument. `None`, the default, reads every record, as today.

```python
def parse_vcf(path: Path, reference: IndexedFasta, *,
              span_limit: int | None = DEFAULT_SPAN_LIMIT,
              num_threads: int = DEFAULT_NUM_THREADS,
              ploidy: int = DEFAULT_PLOIDY,
              contigs: Collection[str] | None = None) -> Iterator[VcfRecord]: ...
```

`contigs` holds the names of the contigs read; a name the reference lacks reads nothing and is no
error. `compare()` passes the contigs `--regions` names, as above, to both readers when
`CompareConfig.regions` is given, and `None` otherwise.

## What `compare()` and the command show

- A contig `--regions` does not name gets no locus, so no progress line: "chr1: 0 loci compared,
  304,084 skipped by --regions", as the HG002 run printed for each chromosome but chr20, goes. On
  a named contig the count of loci skipped by `--regions` is unchanged.
- The rows and the summary are as today. The first error is as today, except the REF error's count
  and wording above.
- In `doc/specs/compare.md`, "The `compare` command", the sentence that loci on a contig `--regions`
  does not name are built without STR regions and counted as skipped goes, since they are not built.

## How it is verified

- **At `parse_vcf`**, on a file with records on chr1, chr2, chr10 and chr3 of a reference in that
  order, with `contigs={"chr1"}`: only chr1's records come, in order, and none of chr10's, whose
  name begins with chr1. A line on another contig with a REF that differs from the reference, too
  few sample columns, a POS of 0, a symbolic ALT or a QUAL that is not a number, each of which
  stops the run or gives a warning on a contig read, gives neither. (A malformed GT gives neither
  on a contig read either: it makes the sample invalid, so it shows nothing here.) Each of the four checks gives today's message and line number on a
  skipped line: fewer than 8 columns; a POS that is not a number; a contig the reference lacks; a
  skipped line before the one before it, and a read line before a skipped one before it. A name in
  `contigs` the reference lacks gives no error. `contigs=None` gives today's records.
- **The REF error**, at `parse_vcf`: a REF that differs on the first record read, after agreeing
  records on a contig not read, says it is the first record read; one after two agreeing records
  read counts 2.
- **At `compare()`**: on the four-contig case of `tests/compare_case.py`, the rows with `--regions`
  equal those of the nine steps of `tests/test_compare.py`, which read every record and select the
  loci with `select_loci`, with `--regions` naming contig 3 alone and cut inside it, at 0 and 2
  finder threads. A REF that differs from the reference on contig 1 does not stop a run with
  `--regions` naming contig 3, and does without `--regions`. A BED line on contig 4 past its end
  reads no record of contig 4.
- **On real data**: HG002 against itself over chr20 gives the same 81,257 rows and summary as
  before, and its time and peak memory go into the report.

## Not in this spec

- Reading only the stretches `--regions` covers within a contig: a locus at a region's edge can
  need records further than any fixed margin, and the owner chose whole contigs.
- `pvc_eval str-regions`, which reads no VCF.
