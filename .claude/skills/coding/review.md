# The categories of a code review

A reviewer reads the top of this file and the section of its category,
and `python.md` beside it, where the rules are with their reasons; they
are not repeated here. What is here is what to look for and how to look.

## Every finding

A finding is something the reviewer can show. It has:

- the place, `file:line`, or the part of the spec it is about;
- what is wrong, in two to four sentences;
- the evidence: the command and its output, pasted; the failing test the
  reviewer wrote in the scratchpad; the sentence of the spec it
  contradicts. A defect of correctness that the reviewer could not make
  happen is not a finding: it goes under "Questions" with what would
  settle it;
- what it causes: a wrong outcome in the table, a run that stops when it
  should go on or goes on when it should stop, a test that cannot fail,
  code that the next milestone will have to work around;
- the severity:
  - **blocker**: a wrong outcome on an input a user can give, a departure
    from the spec, a swallowed error, memory that grows with the genome, a
    failing check, a rule of the spec with no test;
  - **should-fix**: a defect on a plausible edge, a test that would not
    catch a likely mistake, complexity that the next milestone will pay
    for;
  - **nit**: what the tools do not catch and a reader would trip on. At
    most five, and they never start another round of review;
- a fix in a sentence or a few lines. It is a suggestion.

A reviewer does not report what ruff, ty or pytest print; it was given
their output. Nothing hypothetical: "this could fail if" without an input
that reaches it is not a finding. Few findings, each shown, are the aim,
and "No findings" with the list of what was checked is a valid report.

## spec

Does the code do what the spec says?

- Go through the part of the spec the milestone builds, rule by rule and
  case by case. For each one find the code that makes it true and the test
  that would fail if it stopped being true. A rule with no code, or with
  code and no test, is a finding.
- Run the cases. The worked cases of the spec, and two or three of your
  own that it names as edges or that sit on a boundary: records that touch,
  a record at the first base of a contig, a `./1`, a multi-ALT record, an
  empty file. Use `uv run python -c` or a scratch test, and paste what it
  printed.
- Look for what the code does that the spec does not say: a default, a
  clamp, a record skipped, an early return, an extra option. Each is a gap
  of the spec or a defect of the code, and the finding says which you
  think it is.
- When the code and the spec disagree and the code looks right, the
  finding is about the spec.

## tests

Can each test fail, and are the numbers the change claims true?

- For each test of the change, break the code it guards and run it: flip
  a comparison, move a boundary by one, drop a branch, return early. A test
  that still passes guards nothing. Before reporting it, show that your
  change did alter the behaviour on some input, because a change that
  alters nothing proves nothing about the test.
- The ways a case hides a defect: a condition no case reaches; values
  under which several wrong implementations give the same outcome, every
  sample homozygous, every record phased, nothing missing, one sample; a
  parameter that is at its default in every test.
- Expected values are literals. A test that computes its expectation with
  the code under test, or with a copy of its rule, is a finding. So is a
  mock of pvc_eval's own code.
- Every number that a docstring, a test, the commit message or the report
  gives about this change, how many cases, how many tests, how long, is
  computed again, not read again, and reported as right or as wrong with
  the right value.
- Run `uv run pytest --cov --cov-report=term-missing` and read the lines of
  the change that no test runs, error branches first.
- Run every experiment that changes the code with
  `PYTHONDONTWRITEBYTECODE=1` and with the `__pycache__` of the changed
  module removed: on 19 October 2026 a stale `.pyc` made a broken line look
  covered, and fooled two reviewers of the VCF reader for one run each.
- Put the code back after each experiment and end with `git status` clean.

## correctness

What `python.md` asks under "Positions", "Streams" and "Errors", read as
a hunt for wrong outcomes.

- Positions: every conversion between 1-based and 0-based, inclusive and
  half-open, and every `<` against `<=` at a boundary. Find where each
  position comes from and check its convention against the docstring.
- Streams: an iterator read twice, the second time empty; a `list()` over a
  stream; a file closed before its generator is read; `zip` without
  `strict=True`; a look-ahead that drops or repeats the record at its
  edge; the last locus of a contig and of the file.
- Python: a mutable default or a mutable class attribute; a closure over a
  loop variable; `if x:` where 0 or an empty string is a valid value, a
  position 0 among them; `is` where `==` is meant; a collection changed
  while it is iterated; a frozen record whose list field is changed in
  place; a builtin shadowed.
- Sets and dicts: an outcome that depends on the order of a set, or of a
  dict filled in an order the input does not fix.
- Sentinels: `None` or `""` used for "no value" where the other is a real
  value, a missing locus allele beside a deleted one; a sentinel tested by
  truth or with `==`; one made inside a function; a sentinel named against
  `doc/spec.md` §2, `MISSING` for anything but a missing haplotype.
- Threads, as "Threads" in `python.md` asks: an object two threads change,
  a module-level cache or counter written during a run, output whose order
  depends on which thread finished first, an exception in a thread that
  never reaches the caller. Run the stage with one thread and with
  several on the same input and compare the outputs; a race shows as a
  difference that comes and goes, so run it several times.

## errors

The section "Errors" of `python.md` and `doc/spec.md` §3.

- For each way the input can be wrong that §3 lists, what does the user
  see? Make it happen: an unsorted VCF, a contig missing from the
  reference, no sample pairs. Paste the message, and check that it names
  the file and the place and says how to put it right, as the `writing`
  skill asks.
- For each content fault that §5.4 makes invalid, a REF mismatch, a
  malformed GT, contradicting records: does the run go on, and does the
  row say invalid with the right reason?
- An `except` that is broader than what it handles, one that returns a
  default, a record skipped without a count, a defect turned into an input
  error, or an input error that shows a traceback.

## design

The section "Simple code" and "Records, types and names" of `python.md`.

- Code that is more complicated than its job: a class for one function, an
  abstraction with one implementation, a parameter nobody passes, a
  fallback nobody asked for, a module that does two stages. Say what the
  simpler version is.
- The same logic in two places that will have to change together. Dead
  code, commented-out code, a `TODO` with no issue.
- Names against `doc/spec.md` §2: a name the spec says not to use, one
  thing with two names, a name that does not say what the value is, a
  function or method not named with a verb for what it does, a verb not
  in the table of `python.md` or used against its meaning there.
- The public surface: a name without an underscore that no other module
  and no user calls; a module docstring that does not list its public
  names.
- Types: a function without its annotations, which ty does not demand
  and `python.md` does; `Any`, a `cast` or an ignore without its reason, a bare string
  for a closed set, a dict carrying what a dataclass should, an `object()`
  or a one-member `Enum` where a `sentinel` says the same with less.
- Docstrings: the coordinates of every position, what is raised; a
  docstring that only repeats the signature, a comment that says what the
  next line does.
- Clear costs in time or memory that need no profile to see: `x in list`
  inside a loop over records, a stream made into a list, a file read whole
  where the spec streams it. Anything subtler belongs to the `performance` skill, and goes in
  one line for it.
