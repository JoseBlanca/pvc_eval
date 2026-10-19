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
- for each sample, GT:GQ:PS. The GTs follow Hardy-Weinberg proportions, with allele
  frequencies drawn for each record, so they are mostly 0/0 and 0/1 and some 1/1. 3% of
  them are `./.`, and a quarter of the rest are phased, with a PS.

No record is placed where the reference has an N. The same seed (default 1) gives the
same file byte for byte. The reference is GRCh38, the analysis set without alts:

    REF=/Users/jose/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

| input | made with | records | samples | bytes |
|---|---|---|---|---|
| small | `uv run python benchmarks/make_vcf.py $REF chr20:1-5,000,000 --samples 200 --output benchmarks/data/chr20_5mb_200s.vcf.gz` | 18,310 | 200 | 7,355,613 |
| spec size | `uv run python benchmarks/make_vcf.py $REF chr20 --samples 200 --output benchmarks/data/chr20_200s.vcf.gz` | 234,998 | 200 | 95,055,465 |

The small input runs in a few seconds, which is what a profile needs. The spec-size
input is all of chr20 (64,444,167 bases) with 200 samples. `doc/spec.md` §9 says
pvc_eval is meant for benchmarks of megabases and up to a few hundred samples, so this
input is at that size.

Each file was checked on 19 October 2026: `bcftools view -H FILE | wc -l` gave the
records above, and `bcftools norm -c e -f $REF FILE` found no REF mismatch. The norm
moved 540 records of the small input and 6,565 of the spec-size one to the left, the
right-anchored indels. The inputs are in `benchmarks/data/`, which git ignores; they
take 3 s and 32 s to make again.

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

The noise: the small input at the default limit, timed again after the other runs, gave
3.675 s ± 0.031 s. That is 0.102 s, or 2.7%, less than the first time.

Python with `pvc_eval.vcf` imported and nothing read takes 27.4 MB of RSS. tracemalloc
traced a peak of 3.6 MiB on the small input at the default limit.
