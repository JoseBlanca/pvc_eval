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

At the default span limit, on commit 6f8865c (the code of ac58b60):

| input | records | sha256 |
|---|---|---|
| small | 18,310 | 79a790572bb28048339e585db1e51a39a4ec8a32e73db1c96b66ae6bf73a542e |
| small, GATK-like | 18,127 | 74cf45d11f11404781cd14198825e06ce0e0f8a70f8e1a340287f8bc72f02d57 |
| one sample | 98,149 | 02860ac4303f8f8994705370279e8a0d6a06141a2bf1e0db7e226e8954dca361 |
| spec size | 234,998 | 2c627ea747c255e0a12cd8c6c38503e127f174f0cd62028abc1eb7c4d86c4297 |

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
report, `doc/reports/perf-<scope>-<date>.md`. The baseline of the first review is below,
so that the next one can compare against it.

The first baseline was taken on 19 October 2026, on commit ac58b60. The code of
pvc_eval is the same at 6f8865c, which adds only these scripts. The machine was an Apple
M5 Pro (18 cores, 64 GB) running CPython 3.15.0 free-threaded with the GIL off, and
hyperfine 2.0.0. Three other processes each kept one core busy during the whole
measurement: two pytest runs of another worktree and a macOS media service.

| input | span limit | wall time, mean ± sd of 5 | peak RSS, median of 3 |
|---|---|---|---|
| small | 10,000 | 3.777 s ± 0.024 s | 36.0 MB |
| small | none | 4.312 s ± 0.037 s | 614.9 MB |
| spec size | 10,000 | 47.666 s ± 0.313 s | 36.6 MB |
| spec size | none | 59.388 s ± 3.205 s | 7,612.0 MB |

The GATK-like and one-sample inputs were made after this baseline and have not been
timed yet.

The noise: the small input at the default limit, timed again after the other runs, gave
3.675 s ± 0.031 s. That is 0.102 s, or 2.7%, less than the first time.

Python with `pvc_eval.vcf` imported and nothing read takes 27.4 MB of RSS. tracemalloc
traced a peak of 3.6 MiB on the small input at the default limit.
