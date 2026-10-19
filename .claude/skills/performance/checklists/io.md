# Reading and writing

**Purpose.** pvc_eval reads two VCF files, often bgzipped, a FASTA by the
offsets of its `.fai`, and BED files, and writes a table. When the profile
shows the time in reading, decompressing, splitting lines or writing,
this category finds the cost.

**Sent.** When the profile shows the time in I/O, gzip, or the parsing of
lines.

## Rules

- **Measure the floor.** The time to only read and decompress the file,
  `uv run python -c` over `gzip.open(path, "rt")` with nothing done to the
  lines, is the least the reader can take. A reader far above it spends its
  time in parsing; one close to it can only gain by decompressing less.
- **Lines, not characters.** Iterate over the file object, which reads in
  blocks; never `read()` the whole file, and never read one character at a
  time. `gzip.open` in text mode with `encoding="utf-8"`.
- **Split once, and only what is needed.** A VCF line split on tabs once,
  with `maxsplit` when the sample columns are not needed for this record;
  the INFO field left as a string unless a stage reads it. Parsing the
  genotypes of every sample of a record that no locus needs is work for
  nothing.
- **Text or bytes.** Splitting `bytes` and decoding only the fields that are
  kept can be faster than decoding the whole line. It is more code: only
  with a measured gain.
- **The reference by offsets.** A FASTA is read through its `.fai`: the
  offset of a contig, the bases per line and the bytes per line give the
  byte of any position, so a window is read with one `seek` and one `read`.
  Reading a contig line by line to find a position is a finding.
- **Writing.** Lines collected and written in blocks, or through a file
  opened with the default buffer; never a `write` and a `flush` per row.
  The `csv` module, or `"\t".join`, for the table.
