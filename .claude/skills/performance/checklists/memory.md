# Memory

**Purpose.** The peak memory of a run is one of the two targets.
`doc/spec.md` §4 says that nothing holds more than the current locus and a
short look-ahead, so the peak should not grow with the genome. This
category finds what does grow, and what is held larger than it needs to
be.

**Sent.** Always.

## Rules

- **What grows with the input.** Run the workload at two sizes and compare
  the peaks, from `/usr/bin/time -l`. A peak that grows with the number of
  VCF records, and not only with the length of the longest contig, is a
  finding of the highest severity, and tracemalloc's snapshot at the peak
  names the lines that hold it (`profiling_environment.md`). The usual causes: a stream turned into a list, `list()`,
  `sorted()` or `readlines()` over a file; a dict that collects every locus
  or every row and is never emptied; a cache with no bound; a warning that
  collects every record it is about.
- **The reference.** A contig held whole as a `str` costs one byte per
  base, 250 MB for the human chromosome 1, and a list of its bases or a
  `str` per base costs tens of times that. When the peak is the reference,
  the finding says what the stages read of it: a window around the current
  locus, read through the offsets of the `.fai`, is the alternative, and
  its cost in code goes under the effect on the code.
- **Records.** There can be millions of VCF records alive at once only by
  mistake; there are many short-lived ones always. A dataclass without
  `slots=True` carries a `__dict__` per instance; `python.md` asks for
  `slots`, and a record class without it is a finding. A field that keeps
  the whole line of a VCF record, or all its INFO, when the stages read
  four fields of it, holds what nobody reads.
- **Copies.** A slice of a list or a `str` copies it; `copy.deepcopy`,
  `list(d.items())` to iterate, and a list concatenated with `+` in a loop
  copy too. On a short object this is time; on a contig or a whole table it
  is a second peak.
- **Repeated strings.** The same sample name, contig name or allele held
  once per record. `sys.intern`, or one shared object per value made by the
  reader, holds it once. Only when tracemalloc shows these strings at the
  peak.
- **Threads hold their own.** Each worker thread holds what its task
  holds, so the peak of a stage run on eight threads is about eight times
  the peak of one task. A finding about threads states the peak at one
  thread and at the number the run uses.
- **Garbage.** Much time in Tachyon's `<GC>` frames means many objects made
  and dropped: per base, per sample, per record. Making fewer is a finding
  of this category even when the peak is low.
- **What is freed late.** A large value held by a variable that lives
  until the end of a function, after its last use, adds to the peak of
  whatever comes next. Narrow its scope.

A change that saves memory by adding time is weighed with the table of
`SKILL.md`. Say both numbers.
