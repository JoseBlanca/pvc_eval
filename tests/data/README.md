# Test data

- `tandem_repeat/synthetic_ref.fa`: two contigs, ctg1 of 2,369 bases and ctg2 of 346, of
  aperiodic filler with tandem repeats of periods 2 to 6 planted in it, perfect and impure,
  with a homopolymer, a compound motif and repeats near the contig ends. Copied byte for
  byte from pvc at bbcc9f4 (`tests/data/tandem_repeat/`), where its script
  `generate_synthetic_ref.py` made it.
- `tandem_repeat/synthetic_ref.fa.fai`: its `samtools faidx` index, copied from pvc at
  bbcc9f4.
- `tandem_repeat/golden.ssr_catalog.bed.gz`: the 16 STRs that pvc's former production
  catalogue, built with trf-mod, found in `synthetic_ref.fa`, a gzipped BED with `##`
  header lines. Copied from pvc at bbcc9f4.
