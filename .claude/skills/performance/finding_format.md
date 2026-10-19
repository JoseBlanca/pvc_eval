# The format of a performance finding

Every reviewer of a performance review reads this before writing, and the
session reads it before writing the report. It comes from popnei's.

## Severity

A performance finding is not a defect. Its severity says how sure we are
that the change gives a gain that can be measured and is worth its
complexity.

- **Hot**: the place is named in the profile or the benchmark given to
  the reviewer, with a meaningful share of the time or the memory, and the
  change has a clear mechanism. The output that names it is quoted.
- **Likely**: matched by pattern, on a path the profile shows to be hot,
  with no measurement of the place itself yet. Most findings start here.
- **Speculative**: the pattern is there, and it is not clear that the place
  is hot or the gain real. Filed so that it is known, and not tried unless
  something says it matters.
- **Note**: a pointer for future code, or a cold path. Not tried.

A finding from a pattern alone is at most Likely. Code that runs once, at
start up or in an error, is at most a Note.

When the same pattern appears more than about five times, report it once,
with the places.

## Each finding

```
- `src/pvc_eval/file.py:LINE` — [Severity] Title
- Evidence: the lines of the profile or the benchmark, quoted, or "pattern only"
- Mechanism: one to three sentences: what costs, work repeated, a list where a set
  belongs, a stream made into a list, a copy, a Python call per base, and why the
  change removes it
- Axis: time, memory, or both, and what it may cost on the other
- Measurement: the workload, the tool and the number that would confirm it, a count
  where the mechanism can be counted (calls in cProfile, bytes in tracemalloc), and the
  wall time or the peak as the second check
- Effect on the output: none, or how it could change and how that is checked
- Effect on the code: simpler, the same, or more complicated, with what is added or
  removed: lines, a class, a cache, a second path, an invariant to keep
- The change: a diff, or steps complete enough to run
```

Group the findings by severity, Hot first. With nothing to report, write
"No findings." and the list of what was checked.

## Never

- A file or a line that was not read is not cited.
- Output is quoted as it was printed. What was not run is "pattern only".
- No gain is given as a percentage or a multiple unless it was measured.
  Say what is removed: "one dict lookup per base", "the list of all VCF
  records of a contig".
- A finding that belongs to another category is one line at the end, under
  "For another category". A defect of correctness found on the way goes
  there too, marked for the code review, and first when it gives a wrong
  outcome.
