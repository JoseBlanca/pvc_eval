# Algorithms, data structures and repeated work

**Purpose.** The largest gains in pure Python are in what is computed, not
in how a line is written: a step whose cost grows as the square of its
input, a lookup in a list, the same work done twice. These changes usually
make the code clearer as well, which is why this category is always sent.

**Sent.** Always.

## Rules

- **Find the complexity of each stage.** For each stage of `doc/spec.md`
  §4, what its cost grows with: VCF records, samples, the length of a
  locus, the number of possible genotypes. A cost that grows as the product
  of two of them, or as a power of one, is where to look first. The phase
  enumeration of §5.4 tries 2^(h−1) assignments for a diploid sample with
  h heterozygous unphased VCF records in a locus: the record cap bounds h,
  and a finding says how large h gets on the workload.
- **Membership and lookup.** `x in some_list` inside a loop is a scan per
  iteration; a `set` or a `dict` built once makes it one lookup. The same
  for `list.index`, and for a search for a record by position that could
  be a `dict` or `bisect` over a sorted list.
- **Queues and windows.** `list.pop(0)` and `list.insert(0, x)` move the
  whole list; `collections.deque` does not. A look-ahead over a stream is a
  `deque`.
- **Sorting.** Sort once, with `key=`, not inside a loop; `min` or `max`
  for one element; `heapq` for the few smallest; `bisect` to keep a list
  sorted. The VCF records are already in order, as §3 requires, so a sort
  of a stream that is sorted is work done for nothing.
- **Work done twice.** The same sequence built, the same reference read,
  the same record parsed, the same locus allele compared, in two places.
  Compute once and pass it on. A pure function called many times on the
  same arguments can take `functools.cache`, but only with a bound,
  `lru_cache(maxsize=...)`, when its arguments come from the input, since
  an unbounded cache is memory that grows with the genome.
- **Strings of sequence.** Building a sequence with `+=` in a loop copies
  it each time; collect the parts and `"".join` them. A slice of a long
  string copies it; compare with `startswith`, `find` or `==` on the
  slices that are needed.
- **Counting.** `collections.Counter` and `defaultdict` over a hand-written
  `if key in d` and increment.

The fix of a finding here is usually shorter than the code it replaces.
When it is longer, say so under the effect on the code.
