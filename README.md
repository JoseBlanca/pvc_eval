# pvc_eval

Compares two VCF files over a reference genome by rebuilding, for every sample, the sequence each
haplotype carries over each locus, so that two spellings of the same genotype agree. Independent of
pvc; it ports only pvc's STR finder.

Status: specification under review, `doc/spec.md`. The reader of VCF records (§5.2) exists, in `src/pvc_eval/vcf.py`: it trims each record without moving it, gives it its span and yields the records in order of span start. The other stages do not exist yet.
