# How Python is written in pvc_eval

This file is read by every subagent that writes, reviews or fixes code,
before it starts. It holds the rules that no tool enforces, each with its
reason, and the configuration of the tools that enforce the rest. It was
written on 19 October 2026 from research on current practice and from the
coding skills of popnei, before any code of pvc_eval existed, and the
configuration below was run on a scratch project on this machine that day.
It is to be revised after the first milestones, with what they found.

The first goal is that every outcome is right, as `doc/spec.md` says; the
second, that the code stays simple enough to check that it is. Speed and
memory come after, and are the `performance` skill's.

## uv, and nothing else

The owner decided that uv runs everything, and that is not reopened.

| what | command |
|---|---|
| make the environment match the lockfile | `uv sync` (`uv sync --locked` in a check: it fails when the lockfile is stale) |
| add a dependency | `uv add <pkg>`, only with the owner's decision, see below |
| add a development tool | `uv add --dev <pkg>` |
| run the tests | `uv run pytest`, `uv run pytest tests/test_loci.py::test_name -x -q` |
| run pvc_eval | `uv run pvc_eval compare ...` |
| run a script | `uv run python scripts/x.py` |
| a tool once, not a dependency | `uvx <tool>` |

- Never `pip`, `python -m venv`, `requirements.txt`, an activated
  environment or `.venv/bin/python`, with one exception: the block of
  commands for a Tachyon capture, where the program's own process id is
  needed (`.claude/skills/performance/tachyon.md`). Never edit `uv.lock`
  by hand.
- The Python is CPython 3.15, free-threaded: `.python-version` holds
  `3.15+freethreaded`, which `uv python pin 3.15t` writes, and
  `pyproject.toml` has `requires-python = ">=3.15"`. The owner chose it on
  19 October 2026 so that pvc_eval's threads run Python code at the same
  time; see "Threads" below. It was checked that day with uv 0.13.1 on
  CPython 3.15.0: ruff 0.17.1, ty 0.0.89, pytest 9 and hypothesis run on
  it with the GIL off.
- `uvx` does not read the pin, and runs a tool in an environment without
  pvc_eval. A tool that needs either is run with `uv run --with <tool>`.
- The layout is the one `uv init --package` makes: `src/pvc_eval/`, with
  the tests in `tests/` run against the installed package, so a file that
  is missing from the package fails a test and not a user.

## Dependencies

pvc_eval is pure Python. A runtime dependency is a decision of the owner,
asked with what it gives, what it costs and what the standard library
would need instead. The standard library reads what pvc_eval needs: a
bgzipped VCF is a series of gzip members, which `gzip` reads as one
stream; a FASTA with its `.fai` is read by the offsets the `.fai` gives.
Development tools, pytest, ruff, ty, hypothesis, a profiler, are added
with `uv add --dev` without asking, and said in the report.

On 3.15t, a package with a C extension that has not declared that it runs
without the GIL turns the GIL back on for the whole process when it is
imported, with only a `RuntimeWarning` to say so; pyinstrument, scalene
and line_profiler did on 19 October 2026. A development tool that does it
is not imported by the tests or by pvc_eval, and the question of a runtime
dependency includes whether it keeps the GIL off.

## The checks

```
uv sync --locked
uv run ruff format --check
uv run ruff check
uv run ty check
uv run pytest -q
```

All five run before a commit is called done, and the last line of each is
what is reported. pytest's last line has the counts, "214 passed in
3.12s": a run that selected nothing says "no tests ran" and exits with 5,
and a count lower than the one before means a test was lost.

The configuration, in `pyproject.toml`:

```toml
[tool.ruff]
include = ["*.py", "*.pyi"]
line-length = 100
src = ["src", "tests"]

[tool.ruff.lint]
select = ["E", "W", "F", "I", "UP", "B", "SIM", "C4", "PTH", "RET", "RSE", "BLE",
          "TRY", "PIE", "PERF", "T20", "PT", "N", "A", "ARG", "ERA", "RUF"]
ignore = ["E501", "TRY003"]

[tool.ruff.lint.per-file-ignores]
"src/pvc_eval/cli.py" = ["T20"]

[tool.ty.terminal]
error-on-warning = true

[tool.pytest]
testpaths = ["tests"]
addopts = ["-ra", "--import-mode=importlib"]
strict = true
markers = ["slow: takes more than a few seconds; deselect with -m 'not slow'"]

[tool.coverage.run]
source = ["pvc_eval"]
branch = true
```

`include` keeps ruff to Python files. ruff 0.17 also formats the Python
blocks inside Markdown, and on 19 October 2026 it wanted to rewrap a
signature in `doc/spec.md`: the code of a spec or a skill is written for
its reader, and the formatter does not decide its lines.

The rule set is listed and not `ALL`, so that a new version of ruff does
not turn on rules nobody chose. `B` catches mutable defaults and closures
over a loop variable, `BLE` and `TRY` broad and swallowed exceptions, `ARG`
unused arguments, `ERA` commented-out code, `T20` a `print` outside the
command line. A rule that proves noisy is removed here with the count that
showed it.

The type checker is ty, Astral's, a development dependency run with `uv
run ty check`. The owner chose it over mypy on 19 October 2026 for its
speed: on a scratch package that day ty 0.0.89 took 25 ms against mypy
2.4.0's 265 ms, or 65 ms with mypy's cache. Both found the same two planted errors, a value returned without
its sentinel check and a wrong return type. ty's version is pinned by the
lockfile, since its releases still change what it reports; a new version
that finds new errors is taken in a commit of its own that fixes them.
`error-on-warning` makes a warning fail the check, so that nothing it
reports is left standing. ty does not do mypy's `--strict`: the rule that
every function has its annotations is this file's, and the `design`
reviewer checks it.

`[tool.pytest]` and `strict = true` need pytest 9, so the development
group pins `pytest>=9,<10`.

A rule is silenced on the line that needs it, with its code and a reason:
`# noqa: ARG001 -- the signature every stage of §4 shares`,
`# ty: ignore[invalid-argument-type]  # reason`. An ignore that names the
wrong rule silences nothing, as tried on 19 October 2026. Never a bare
`noqa` or `ty: ignore`, never a file-wide setting to get a milestone
green.

## Before the code

Code is written from the spec and from a milestone of a plan. Read the
part of the spec, the parts of `doc/spec.md` it stands on, and the code it
calls.

When the spec can be read in two ways that give different outcomes, work
the case where they differ by hand, or run pvc, bcftools or the program
the spec names, before choosing. The choice made without it is a guess
the tests will not catch, because they were written from the same reading.

When the spec does not say what should happen in a case, the choice is
not made in silence. A choice that changes an outcome, the table or the
interface goes to the owner. A smaller one is made, written into the spec
in a commit of its own before the code, and named in the commit message.

## The order of the work

1. The test first. It asserts what the spec gives, as literals, and it
   fails before the change: the new function exists with a body that
   returns a wrong value, `[]`, `None`, and the test fails on its
   assertion, not on an import. A test that reads a stream until the code
   under test says stop has a bound of its own, so that a wrong stub ends
   in a failure and not in a test that never ends.
2. The code, as small as the milestone asks for.
3. The checks above, all of them, with their real output. A check that
   could not be run is reported as not run.
4. The commit, with the message the `writing` skill describes.

## Simple code

What the owner asked for on 19 October 2026: a change that makes the code
clearer, simpler or better structured is worth taking even for a small
cost in speed, and a change that makes it more complicated needs a gain
that pays for it.

- Functions over classes. A module of functions that transform data, one
  stage of `doc/spec.md` §4 per module, named after it: `str_regions.py`,
  `vcf.py`, `loci.py`, `reconcile.py`, `evaluate.py`, `table.py`,
  `cli.py`. No `utils.py`, no `Manager`, `Handler` or `Base`, no class
  that holds one method.
- The concrete thing first. No abstract base class, protocol, factory,
  registry or option for a second implementation that does not exist. The
  third copy of a piece of logic becomes a function; the second may stay.
- No code for a case the spec does not have: no parameter nobody passes,
  no fallback, no backwards compatibility, since there are no users of an
  older version.
- Short functions with early returns. A function longer than a screen, or
  with more than three levels of nesting, is split along a line the spec
  has.
- Input and output at the edges. The stages are pure functions over
  iterators; only the readers open files and only the table writes them,
  so every rule can be tested on values built in the test.

## Records, types and names

- A record is `@dataclass(frozen=True, slots=True)`, with `kw_only=True`
  above three fields: a VCF record, a locus, an STR region. Frozen,
  because a record that a later stage changes is a record whose history
  nobody can follow; slots, because there are millions of them and
  `slots` takes the `__dict__` of each away. Its invariants are checked in
  `__post_init__`, which raises `ValueError`.
- A closed set is an `Enum` or a `StrEnum`: the outcomes, the kinds, the
  reasons a locus is invalid. Never a bare string compared in several
  places.
- Type hints on everything, in the syntax of 3.12 and later: `list[str]`,
  `X | None`, `def first[T](xs: Sequence[T]) -> T`, `type Span = tuple[int,
  int]`. Parameters take the abstract type, `Iterable`, `Sequence`,
  `Mapping` from `collections.abc`; returns give the concrete one, or an
  `Iterator` for a stream. No `Any`, and a `cast` carries a comment that
  says why it holds.
- A name says what the value is, in the words of `doc/spec.md` §2:
  `vcf_record`, `locus_alleles`, `possible_genotypes`, `record_cap`. Never
  `data`, `tmp`, `val`, `rec`, `n`. A word the spec says not to use, a
  variant for a VCF record, truth for file A, is not used in a name
  either.
- A function or a method is named for the action it does, with a verb:
  `get_bases`, `parse_vcf`, `find_span`, never `bases`, `span` or
  `length`. The owner asked for this on 19 October 2026: a call then reads
  as what it does, and a noun is left for values. Private helpers follow
  it too. A property is a value and keeps a noun: `held_region`,
  `contigs`, `end`.
- The verbs come from this table, so that one verb means one thing in the
  whole project, as `doc/spec.md` §2 does for the nouns. The owner asked
  for it on 19 October 2026, after `read` and `parse` had been used for the
  same work in the VCF reader. When no verb fits, a row is added here, with
  its meaning, in the same commit as the code that first uses it.

  | verb | means |
  |---|---|
  | read | take bytes or lines from a file, without interpreting them |
  | parse | turn text into values or objects of pvc_eval |
  | write | put text or bytes into a file |
  | get | return something held, or cheap to work out from it |
  | find | search for something, or work out where it is |
  | check | test a condition and raise when it fails |
  | make | build a new value, object or message |
  | count | return how many of something there are |
  | trim, drop, evict | remove part of a value or of what is held |
  | sort | put in order |
  | format | turn a value into text for a person |
  | is | answer a yes-or-no question about a value, returning a bool |
  | open | open a file and return it, to be read or written |
  | batch | group the items of a stream into lists of a fixed length, the last perhaps shorter |
  | map | apply a function to each item of a stream, giving the results in the order of the items |
  | note | pass the items of a stream on unchanged, adding a fact about each to a value the caller holds |
  | log | send a message to the logger, for the user |
  | put | hand an item to another thread through a buffer or a line they share, waiting while it is full |
  | take | receive an item another thread put, waiting while there is none |
  | start | start threads that run until they end or are stopped |
  | stop | set the stop of a run's threads, or wake those that wait, so that they end |
  | reconcile | build a locus's locus alleles and its samples' sample results (`doc/specs/reconciliation.md` §4) |
  | evaluate | give the outcome of a sample pair (`doc/specs/reconciliation.md` §5) |

- Public is only what another module or the user calls: everything else
  starts with an underscore, functions, classes and constants alike. The
  owner asked for it on 19 October 2026: a small public surface is one a
  user can learn, and one that can change inside without breaking anyone.
  A module's docstring names its public names.
- A default that changes an outcome is a named constant with a comment
  that says where it comes from: `RECORD_CAP = 30  # doc/spec.md §5.3`.
- "No value" that `None` cannot say is a sentinel, below.
- `pathlib.Path` for every path, `encoding="utf-8"` for every text file,
  and every file opened in a `with`.

## Sentinels

Python 3.15 has a builtin, `sentinel(name)`, that makes a unique object
whose repr is its name, which survives `copy` and `pickle` as itself, and
which can stand in a type hint: `int | NOT_GIVEN`. ty narrows
`x is NOT_GIVEN` as it narrows `x is None`, so a value used without the
check is a type error; both were tried on 19 October 2026.

It is for a value that means "nothing here" where `None` or a falsy value
is ambiguous:

- **An argument that was not given**, where `None` is a value the caller
  can pass: `def window(..., end: int | None | NOT_GIVEN = NOT_GIVEN)`.
- **A value that the data can lack, where the empty one is real.** A locus
  allele can be the empty sequence, a deletion of the whole locus
  (`doc/spec.md` §6), and `""` is falsy; a haplotype can be missing (§2).
  So a missing locus allele is not `None` beside `""`, nor `""` itself, but
  a sentinel, and the code that reads an allele says which of the two it
  has with `is`.

The rules:

- Made once, at the top of a module, assigned to a name that matches its
  own: `NOT_GIVEN = sentinel("NOT_GIVEN")`. Pickling finds it by that name,
  and a sentinel made inside a function cannot be pickled.
- Compared with `is`, never `==`, and never tested by truth: a sentinel is
  always truthy, so `if allele:` is wrong for both the sentinel and `""`.
- Its name follows `doc/spec.md` §2. The spec gives *missing* one meaning,
  a haplotype the VCF genotype says nothing about, so `MISSING` is the
  sentinel for that and nothing else; an argument not given is
  `NOT_GIVEN`.
- Not `object()`, and not a one-member `Enum`, which were the ways before
  3.15: the first has no useful repr and no type, the second is more code
  for the same thing.

## Threads

pvc_eval runs on the free-threaded 3.15t, where threads run Python code
at the same time on all the cores. Threads are the way pvc_eval does work
in parallel; processes are the owner's to decide. Code that is not yet
parallel is still written so that it can be:

- **What threads share is read only.** The frozen records, the reference,
  the settings. A thread writes only to objects it made, and hands them
  back as its result. A list, dict or set that two threads change is a
  defect, even where CPython would keep it whole: the order of the changes
  is not fixed, and the outcome then is not either.
- **No module-level state that changes**: no global counter, cache or
  list filled during a run. A cache shared by threads, `functools.cache`
  included, is a decision that is written down, with why it is safe.
- **`concurrent.futures.ThreadPoolExecutor`**, with `map`, which gives the
  results in the order of the input. The output of pvc_eval never depends
  on which thread finished first. An exception raised in a thread reaches
  the caller through the result, and is not lost.
- **The number of threads is a parameter**, with 1 meaning no pool at all,
  so that a test and a profile can run the same code serially.
- **A stage that runs on threads is tested with one thread and with
  several**, over input that gives each thread work, and the outputs are
  compared as equal.
- **The GIL stays off.** A test imports pvc_eval and asserts
  `sys._is_gil_enabled()` is false, so that an import that turns it back on
  fails the checks instead of making the threads take turns in silence.
- `threading.Lock` only where the design says so, never per record: a lock
  in the hot path makes the threads take turns.

## Positions

pvc_eval reads VCF, 1-based, and BED and the `.fai`, 0-based half-open, and
writes a table that is 1-based and inclusive (`doc/spec.md` §6). An error
of one base is the most likely defect it can have, and it gives plausible
wrong loci.

- Inside pvc_eval every position is 0-based and every span half-open,
  `[start, end)`, as Python slices are, so `reference[start:end]` is the
  span. A position is converted where it is read and where it is written,
  in the reader and in the table, and nowhere else.
- The docstring of every function that takes or returns a position says
  so: "0-based, half-open". A name that holds a 1-based position says it,
  `pos_1based`, in the one reader that has it.
- Every rule with a boundary is tested at the boundary: records that touch
  and do not overlap, a record on the first and on the last base of a
  contig, a region that ends where a locus starts.

## Streams

`doc/spec.md` §4 says that nothing holds more than the current locus and a
short look-ahead, so memory does not grow with the genome. So:

- A stage takes an iterator and returns one, written as a generator.
  `list()` over a stream is a defect unless the list is bounded by a locus,
  and its comment says by what.
- An iterator is read once. A function that needs two passes over its
  input takes a sequence, and its type says so.
- A generator that opens a file opens it inside the generator, in a
  `with`, so that the file is open while it is read and closed when the
  stream ends.
- `zip(..., strict=True)`, so that two streams of different lengths raise
  and do not lose the tail of the longer.

## Errors

`doc/spec.md` §3 splits what goes wrong in two. An error of the input that
stops the run, an unsorted VCF, a contig missing from the reference, raises
an exception of pvc_eval's own, `PvcEvalInputError(ValueError)`, whose
message is the one the user reads, written as the `writing` skill says
for errors: the file, the line or the position, and how to put it right.
The command line catches it alone, prints its message to stderr and exits
with 1. Anything else is a defect of pvc_eval and is left to rise with its
traceback.

What the spec calls invalid is not an error: a malformed or haploid GT,
contradicting VCF records. It becomes an outcome of the row, and the run
goes on. A REF that disagrees with the reference is an error that stops the
run (decided by the owner on 19 October 2026, `doc/spec.md` §5.2).

- No silent fallback: no `except` that returns `None`, `[]` or a default,
  no `dict.get(key, default)` where a missing key is a defect, no record
  skipped without a count. The spec says what is skipped and counted, the
  symbolic alleles of §9, and nothing else is.
- Catch the narrowest exception, where something useful can be done with
  it, and `raise ... from err` when one is turned into another.
- `assert` is for what the code itself guarantees and a test would catch,
  never for an input, since `python -O` removes it.
- Logging in the library, `_logger = logging.getLogger(__name__)`, set up
  only in `cli.py`. `print` only in `cli.py`.

## Tests

- The expected values are literals, from the spec's worked cases or from
  a program that was run, with the command in a comment. A test never
  computes its expected value with the code under test or a copy of its
  rule.
- A case is built in the test where it fits: a few bases of reference and
  three VCF lines as strings, written to `tmp_path` when a reader needs a
  file. Larger inputs go in `tests/data/`, small enough to check by hand,
  each with a line in `tests/data/README.md` that says what it holds and
  where it came from. Nothing is written outside `tmp_path`.
- A test has to be able to fail. Check the case against the usual ways it
  cannot: a condition no case reaches; a case where several wrong
  implementations give the same outcome, every sample homozygous, every
  record phased, nothing missing; a parameter left at its default in every
  test, so the code could ignore it. When in doubt, break the code on
  purpose and see the test fail.
- A test calls the public function where the behaviour can be seen
  there. It may call a private function when the rule is much clearer
  there than through the public one, as the spans of `doc/spec.md` §5.2
  are on `_find_span` rather than through a whole VCF file; a private
  function is not made public for its tests.
- One test per behaviour, named for the behaviour and its outcome:
  `test_insertion_before_repeat_pulls_repeat_into_locus`. Cases of one
  rule go in one `pytest.mark.parametrize` with `ids`.
- The properties of `doc/spec.md` §7.3 and §7.4, a file against itself
  and against its re-spellings, are tested with hypothesis over generated
  references and VCF records, as well as on the fixed files: they are what
  finds the case nobody thought of. A failing example hypothesis finds is
  added as a fixed case.
- A table compared whole is a golden file in `tests/data/expected/`, written
  again only when the test is run with `--update-golden`, and the diff is
  read before it is committed.
- No mocks of pvc_eval's own functions, of files or of the code under test.
  The real function on a small input is the test.
- Floats, where there are any, with `pytest.approx` and a tolerance that
  has a reason.
- The suite runs in seconds. A test over a real region of a genome is
  marked `slow`, and the checks above still run it.
