# Interpreter cost in hot loops

**Purpose.** When the profile shows the time in a few small Python
functions called millions of times, once per base or per sample of every
record, the cost is the interpreter's: a function call, an attribute
lookup, an object made per item. This category is for those loops and
only those. Most of its patterns make the code a little harder to read,
so a finding here needs a measured gain to be accepted, by the table of
`SKILL.md`.

**Sent.** When the profile shows a few functions with a high call count
and most of the self time.

## Rules

- **Only the loop the profile names.** The call count of cProfile and the
  self time of Tachyon both point at it. A pattern elsewhere is not a
  finding. Tachyon's `--heatmap` shows the lines inside it, and
  `--opcodes` whether the interpreter specialized them (`tachyon.md`).
- **Per-item work moved into a built-in.** A loop in Python over the bases
  of a sequence is slower than one call of a method of `str` or `bytes`
  that does the same in C: `str.count`, `str.find`, `str.translate`,
  `bytes.split`, `==` on two slices, `zip` and `sum` over a generator. The
  same for `any`, `all`, `min`, `max`, `sorted`, and `itertools`. These are
  usually as clear as the loop, or clearer.
- **A comprehension over a loop that appends.** Faster and usually
  shorter.
- **Calls in the innermost loop.** A small helper called per base costs a
  call per base. Inlining it is a finding only with a measured gain larger
  than the noise on the workload, since it makes the code longer.
- **Lookups bound to a local.** `append = out.append` before a hot loop
  saves an attribute lookup per iteration; since Python 3.11 the saving is
  small. Only on the hottest loop, only measured, and the table of
  `SKILL.md` usually rejects it.
- **Exceptions.** A `try` costs nothing when nothing is raised, since 3.11.
  An exception raised per item, as a test for something common, costs
  much; a test before the call does not.
- **Old advice.** Much of what was faster before Python 3.11 is not now.
  A finding that rests on a rule of thumb is measured on the project's
  3.15t before it is believed, since advice measured on the standard build
  need not hold on the free-threaded one.

The fix that leaves pure Python, numpy, Cython, mypyc, a Rust extension, is
not this category's to propose as a change: it goes in one line for the
session, which takes it to the owner.
