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

## `vcf_index/`: the indexes of a VCF, for `tests/test_vcf_index.py`

Two VCFs of one sample, written by hand on 20 October 2026, whose headers name four contigs of
1,000 bases, ctgA to ctgD, and their bgzipped copies with bcftools 1.24's indexes:

- `records_on_b_d.vcf`: three records, two on ctgB and one on ctgD. No index is beside it, so
  it is also the plain-text VCF of the tests.
- `records_on_b_d.vcf.gz`, `records_on_b_d.vcf.gz.tbi`, `records_on_b_d.vcf.gz.csi`: its
  bgzipped copy and its two indexes, which list ctgB and ctgD, the contigs with a record, and
  not ctgA or ctgC, which only the header names.
- `records_on_a_b.vcf`: two records, one on ctgA and one on ctgB.
- `records_on_a_b.vcf.gz`, `records_on_a_b.vcf.gz.csi`: its bgzipped copy and its `.csi`, which
  lists ctgA and ctgB. The tests put this `.csi` beside `records_on_b_d.vcf.gz`, as an index out
  of date would be.

Made in `tests/data/vcf_index/` with:

```
bcftools view -Oz -o records_on_b_d.vcf.gz records_on_b_d.vcf
bcftools index -t records_on_b_d.vcf.gz
bcftools index -c records_on_b_d.vcf.gz
bcftools view -Oz -o records_on_a_b.vcf.gz records_on_a_b.vcf
bcftools index -c records_on_a_b.vcf.gz
```

The decompressed `.tbi` of `records_on_b_d` is 214 bytes and its `.csi` 234; the tests cut and
change their headers at the byte offsets their comments give.
