# pvc_eval

Compares two VCF files over a reference genome by rebuilding, for every sample, the sequence each
haplotype carries over each locus, so that two spellings of the same genotype agree. Independent of
pvc; it ports only pvc's STR finder.

Status: specification under review, `doc/spec.md`. The reader of VCF records and their normalising (§5.2) exist, in `src/pvc_eval/vcf.py`; the other stages do not yet.
