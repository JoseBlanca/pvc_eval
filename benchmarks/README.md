# Benchmarks

These are the workloads of the performance reviews of pvc_eval, kept here so that each
review runs the same commands on the same inputs as the one before it. The method is
the `performance` skill's (`.claude/skills/performance/SKILL.md`), and the tools that
work on this machine are listed in its `profiling_environment.md`.

## The VCF reader: parse_vcf

The workload is `parse_vcf` (`src/pvc_eval/vcf.py`). It reads the VCF records of one
file, checks and trims them, gives each its span (`doc/spec.md` §5.2), reading the
reference through a `FastaWindow`, and gives them in order of span start.
`run_parse_vcf.py` reads every record and prints how many there were:

    uv run --no-sync python benchmarks/run_parse_vcf.py VCF FASTA [SPAN_LIMIT]

SPAN_LIMIT is a whole number of 1 or more, or `none` for no limit. Without it the
script uses parse_vcf's default, 10,000 bases. Each input is timed with both.

### The inputs

No real VCF on this machine matches a reference that is on it, so `make_vcf.py` writes
one. The file is a bgzipped VCF of one region of the reference, and every REF is copied
from the reference. It holds:

- SNPs, about one every 300 bases;
- one indel of 1 to 20 bases for every nine SNPs; half of them change whole units of a
  repeat run found in the sequence (a unit of 1 to 4 bases repeated at least three
  times), and a third of those are written at the right end of the run, which `bcftools
  norm` moves to the left;
- two ALTs in about 5% of the records, and a `*` ALT in a SNP inside one deletion in 50;
- a GT for each sample. The GTs follow Hardy-Weinberg proportions, with allele
  frequencies drawn for each record, so they are mostly 0/0 and 0/1 and some 1/1. 3% of
  them are `./.`, and a quarter of the rest are phased, with a PS.

The sample columns come in two kinds:

- **GT:GQ:PS**, the default, about 10 bytes per sample. The GQ is drawn apart from the
  GT, and every phased GT of a record has the same PS, the first position of its 100 kb
  block. So most samples of a record share their genotype whole, which real callers'
  output does not.
- **GATK-like**, with `--column gatk`: GT:AD:DP:GQ:PL:PS, about 28 bytes per sample.
  Each sample has a depth near 30, allele depths that follow its GT, PLs worked out from
  those depths with 0 for the GT, and a GQ taken from the PLs. Each sample has phase
  blocks of its own, a new one about every 50 records. The INFO holds AC, AF and AN
  counted from the GTs, DP, MQ, QD and the other annotations GATK writes, about 250
  bytes.

`--mean-gap` sets the mean number of bases between records (default 270, which puts the
SNPs about 300 bases apart).

No record is placed where the reference has an N. The same seed (default 1) and options
give the same file byte for byte. The reference is GRCh38, the analysis set without
alts:

    REF=/Users/jose/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

There are four inputs. Each is made with `uv run python benchmarks/make_vcf.py $REF`
followed by the region and options below:

| input | region and options | records | samples | bytes |
|---|---|---|---|---|
| small | `chr20:1-5,000,000 --samples 200 --output benchmarks/data/chr20_5mb_200s.vcf.gz` | 18,310 | 200 | 7,355,613 |
| small, GATK-like | `chr20:1-5,000,000 --samples 200 --column gatk --output benchmarks/data/chr20_5mb_200s_gatk.vcf.gz` | 18,127 | 200 | 20,250,013 |
| one sample | `chr20 --samples 1 --mean-gap 650 --column gatk --output benchmarks/data/chr20_1s_gatk.vcf.gz` | 98,149 | 1 | 5,654,025 |
| spec size | `chr20 --samples 200 --output benchmarks/data/chr20_200s.vcf.gz` | 234,998 | 200 | 95,055,465 |

The small input runs in a few seconds, which is what a profile needs. The GATK-like
small input has the sample columns and the phase sets of real callers, so it tests
whether a gain that relies on samples sharing their genotype holds on real data. The
one-sample input covers chr20 at one record per 650 bases. It is the common case of
comparing one sample against a truth set, where the work per record weighs more than
the work per sample. The spec-size input is all of chr20 (64,444,167 bases) with 200
samples. `doc/spec.md` §9 says pvc_eval is meant for benchmarks of megabases and up to a
few hundred samples, so this input is at that size.

Each file was checked on 19 October 2026: `bcftools view -H FILE | wc -l` gave the
records above, and `bcftools norm -c e -f $REF FILE` found no REF mismatch. The norm
moved the right-anchored indels to the left: 540 records of the small input, 534 of the
GATK-like one, 2,798 of the one-sample one and 6,565 of the spec-size one. The inputs
are in `benchmarks/data/`, which git ignores. Making them again takes 3 s, 9 s, 1 s and
32 s.

### The output to compare

The timed script prints only the number of records, which a change to a span or a
genotype would leave as it was. `digest_parse_vcf.py`, which is not timed, prints the
number and a sha256 of the repr of every record, in the order parse_vcf gives them, and
every experiment compares its digest with the baseline's:

    uv run --no-sync python benchmarks/digest_parse_vcf.py VCF FASTA [SPAN_LIMIT]

At the default span limit, with VcfRecord without its field `ref_mismatch`:

| input | records | sha256 |
|---|---|---|
| small | 18,310 | 798c27cc82733aa70da7b61bb6a364783cd7427144522a15ff600f36990ac27f |
| small, GATK-like | 18,127 | 2a494a745203118da1588498dd270025d4b48a6a8c31eb85892e13972684b98d |
| one sample | 98,149 | 6b7cc70543d2f6e043b894ef1be72eee1296836b151f25f8c94cf705ce5d6b79 |
| spec size | 234,998 | afefe073faad3daa92e0359e8de78460d49b9a97851dfd32af70bfef772703f8 |

These replaced, on 19 October 2026, the digests of commit 6f8865c (the code of ac58b60),
because the repr of every record lost `ref_mismatch=False` when a REF that differs from
the reference became an error that stops the run (`doc/spec.md` §5.2). The inputs were
made again that day with the commands above, and their bytes and records equal the
table's. The records are the same as before: the repr of each, with
`ref_mismatch=False, ` put back after its line number, gave the old digest on the small
input, 79a790572bb2..., and on the one-sample one, 02860ac4303f....

### Measuring

The wall time is measured with hyperfine, 5 runs after 1 warm-up:

    hyperfine -N --warmup 1 --runs 5 \
      'uv run --no-sync python benchmarks/run_parse_vcf.py benchmarks/data/chr20_5mb_200s.vcf.gz $REF'

The peak memory is the "maximum resident set size" line of `/usr/bin/time -l` on the same
command, the median of 3 runs. The noise is the small input at the default limit
measured twice, apart; the difference between the two is the smallest change an
experiment can claim.

The evidence, on the small input at the default limit:

    uv run python -m cProfile -o prof.out benchmarks/run_parse_vcf.py benchmarks/data/chr20_5mb_200s.vcf.gz $REF
    uv run python -c "import pstats; pstats.Stats('prof.out').sort_stats('tottime').print_stats(30)"
    uv run python -c "import pstats; pstats.Stats('prof.out').sort_stats('cumtime').print_stats(30)"
    uv run python benchmarks/peak_parse_vcf.py benchmarks/data/chr20_5mb_200s.vcf.gz $REF

`peak_parse_vcf.py` prints the peak bytes tracemalloc traced. It also prints the 15
lines that held the most memory at the point between two records when the most was
held. The Tachyon captures are taken by the owner, as `tachyon.md` says, into
`.profiles/`.

### The numbers of each review

Each review writes its baseline, its noise, its profiles and its experiments in its
report, `doc/reports/perf-<scope>-<date>.md`. The numbers of the parser as the last
review left it are below, so that the next review can compare against them, and after
them the baseline of the first review.

#### After the first review: commit 27116a8, 19 October 2026

The first review (`doc/reports/perf-vcf-parser-2026-10-19.md`) kept three changes, P1, P2
and P3. Its code was measured on 19 October 2026 between 16:13 and 16:35, at commit
27116a8 on branch `perf/vcf-parser`, against the code before the review (commit 2911203),
on the same machine and Python. Each input and limit was one hyperfine call that ran both,
10 runs after 1 warm-up (5 for the spec-size input at 10,000, 3 with no limit). The wall
time is the mean ± sd, and the peak RSS the median of 3 (1 for the spec-size input with no
limit). The digests at the default limit equal those of commit 6f8865c, which the ones above replaced.

| input | span limit | wall time | peak RSS | before the review, measured beside it |
|---|---|---|---|---|
| small | 10,000 | 1.054 s ± 0.008 s | 36.3 MB | 3.532 s ± 0.023 s, 35.9 MB |
| small | none | 1.089 s ± 0.008 s | 87.3 MB | 4.070 s ± 0.024 s, 614.9 MB |
| small, GATK-like | 10,000 | 1.463 s ± 0.014 s | 36.8 MB | 4.015 s ± 0.050 s, 35.9 MB |
| small, GATK-like | none | 1.619 s ± 0.024 s | 140.4 MB | 4.570 s ± 0.068 s, 610.1 MB |
| one sample | 10,000 | 735.4 ms ± 3.6 ms | 36.5 MB | 780.6 ms ± 13.4 ms, 34.0 MB |
| one sample | none | 828.2 ms ± 5.6 ms | 118.7 MB | 868.0 ms ± 10.7 ms, 128.8 MB |
| spec size | 10,000 | 13.104 s ± 0.180 s | 37.4 MB | 45.075 s ± 0.297 s, 36.7 MB |
| spec size | none | 14.322 s ± 0.295 s | 706.8 MB | 53.098 s ± 0.326 s, 7,611.9 MB |

The one-minute load average was between 1.34 and 2.89. Test runs and a script of other
sessions held one core at times; each call shared that load between its two sides. The
noise was not measured again: the small inputs' 3% to 5% below is where to start.

#### The baseline of the first review: commit df41b16, 19 October 2026

The baseline was measured on 19 October 2026 between 15:02 and 15:26, at commit df41b16,
whose code of pvc_eval is that of ac58b60. The machine was an Apple M5 Pro (18 cores,
64 GB) running CPython 3.15.0 free-threaded with the GIL off, and hyperfine 2.0.0. The
wall time is the mean ± sd of 5 runs, and the peak RSS the median of 3. The spec-size
input with no limit was run 3 times for the wall time and once for the peak.

| input | span limit | wall time | peak RSS |
|---|---|---|---|
| small | 10,000 | 3.540 s ± 0.141 s | 36.0 MB |
| small | none | 4.013 s ± 0.013 s | 614.9 MB |
| small, GATK-like | 10,000 | 3.947 s ± 0.017 s | 35.8 MB |
| small, GATK-like | none | 4.517 s ± 0.024 s | 610.1 MB |
| one sample | 10,000 | 767.1 ms ± 2.4 ms | 34.0 MB |
| one sample | none | 861.5 ms ± 15.2 ms | 128.8 MB |
| spec size | 10,000 | 44.280 s ± 0.479 s | 36.7 MB |
| spec size | none | 53.102 s ± 1.283 s | 7,611.9 MB |

The noise is each input at the default limit timed again, 10 to 12 minutes after the
first time:

| input | first | second | difference |
|---|---|---|---|
| small | 3.540 s ± 0.141 s | 3.372 s ± 0.023 s | −0.168 s, −4.7% |
| small, GATK-like | 3.947 s ± 0.017 s | 3.833 s ± 0.034 s | −0.114 s, −2.9% |
| one sample | 767.1 ms ± 2.4 ms | 734.0 ms ± 3.9 ms | −33.1 ms, −4.3% |
| spec size | 44.280 s ± 0.479 s | 44.608 s ± 0.312 s | +0.328 s, +0.7% |

All three small inputs ran faster the second time, by 3% to 5%. The first run of the
small input includes one run of 3.79 s among four between 3.45 s and 3.50 s. So an
experiment on the small inputs claims a gain only above about 5%, and only when it is
timed in one hyperfine call that alternates it with the baseline.

The one-minute load average was 1.69 before the first run, between 0.75 and 2.42
during the measurements, and 0.78 after the last. Before and after each run, no
process used more than 50% of a core.

An earlier baseline, at ac58b60 between 14:13 and 14:30, ran while three other
processes each kept one core busy. It gave 3.777 s ± 0.024 s for the small input and
47.666 s ± 0.313 s for the spec-size one at the default limit, 7% to 8% slower than the one
above. Its peak RSS was the same.

Python with `pvc_eval.vcf` imported and nothing read takes 27.4 MB of RSS. tracemalloc
traced a peak of 3.6 MiB on the small input at the default limit.
