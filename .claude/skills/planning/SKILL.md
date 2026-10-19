---
name: planning
description: How an implementation plan is written in pvc_eval. Use it when writing or revising a plan under doc/plans/, which turns an approved spec, or a part of one, into milestones with deliverables that can be checked, each one built, reviewed and fixed by subagents as the coding skill says. Also for a small change with no spec, which gets a short plan of one milestone. The prose follows the writing skill.
---

# Writing an implementation plan

A plan turns an approved spec into the order of the work. It says what is
built first and what after, in which milestones, and how we will know
that each milestone is done. It is carried out by the `coding` skill: a
session that sends each milestone to a subagent, has it reviewed by other
subagents, has the findings fixed by another, and goes on to the next
one. So a plan has to be something that can be run that way.

A plan decides nothing about the design. When writing it shows that the
spec left something open, the question goes back to the spec, as an open
point for the owner, and the plan waits or is written around it.

The prose follows the `writing` skill. A plan lives in
`doc/plans/<name>.md`. The name says what it builds, in lower case with
hyphens, `str-finder`, `loci`, and it is also the name of the branch,
`plan/<name>`, and of the report, `doc/reports/<name>.md`.

This skill comes from the writing-plans skills of popnei and popnei_web,
ported on 19 October 2026, before any code of pvc_eval. Where they said
work package and task, this one says milestone, and a milestone is the
unit one subagent builds.

## No filler

The number of lines is not the measure of a plan. Filler is. A plan has
as many lines as the session and the subagents need, and none besides.

The test for a sentence: would the session, or the subagent that gets the
milestone, do something different without it? What fails that test is
filler. The usual kinds: what repeats the spec, whose rules and cases the
plan points at, because a copy is wrong as soon as the spec changes; a
part that is there because the shape below has a place for it.

## Before writing

- The spec the plan builds from is approved by the owner, as the `spec`
  skill ends, and its open points are answered or each has a "meanwhile".
  For each one not answered the plan says which milestone the answer
  would change. A small change with no spec, a fix, an option, a
  refactor, gets a plan of one milestone whose first part says what the
  change is and why, in a paragraph.
- Read the spec, the code that exists and its tests. A plan written
  without looking at the code plans work that is done, or builds on
  something that is not there.
- Run what the checks will run. A number a check compares with, the STR
  regions pvc gives over a region, the outcome of a hand-worked case, goes
  into the spec after its command has been run once on this machine, and
  the plan names the part of the spec that holds it.
- Make the sketch the `writing` skill asks for. For a plan, its points are
  the breakdown: the milestones in order, each with what it gives and what
  it stands on, the reason for the order, and what is left out.
- Show the owner the breakdown before the plan is written, as a reply in
  chat. The question is whether the milestones are the right size and in
  the right order. The owner corrects ten lines faster than a finished
  plan.

## Milestones

A milestone ends with something that exists and can be checked, and
leaves pvc_eval in a state that works, with every check of the `coding`
skill passing, whatever comes after. It is the unit given to one
subagent, so it is between half a day and two days of a person's work: a
smaller one is part of its neighbour, because each milestone costs a
prompt, a review by several subagents and a round of fixes; a larger one
is split, because a subagent that runs out of context in the middle
leaves half a milestone that nobody can review.

The work is split the way `doc/spec.md` §4 splits the program, so a
milestone is usually one stage of the pipeline, or a part of one, with
its tests: the STR finder and its ported tests; reading, checking and
normalising VCF records; grouping them into loci; reconciliation;
evaluation; the table and the summary; the command line. A stage is a
function from one stream to the next, so each can be built and tested on
hand-written inputs before the stage before it exists.

Put first the milestone that would change the plan if it failed. For
pvc_eval that is usually the hardest rule on the hardest data: the loci
and the reconciliation over the STR cases of `doc/spec.md` §7.5, before
the table and the command line that only carry their results.

A whole run, end to end, comes early, thin, and is filled in after: a
command that reads two small VCFs and writes a table, with stages that do
the least they can. It finds where the stages do not fit together while
that is cheap to change.

When the code that exists makes the new work hard, the change to that
code is a milestone of its own, before the ones that need it. It changes
no result, and its check is that the tests that were there pass
untouched.

Each milestone has:

- **What it gives**, in a sentence or two: what a user, or the next
  stage, can do when it is done.
- **The part of the spec it builds**, by section, and nothing copied from
  it.
- **Its deliverables, each with the way to check it.** The check can be
  run and can fail: the tests of a part of the spec, named by that part,
  with how many there are; a command that runs pvc_eval on a file of
  `tests/data/` and what its output must hold. "Reconciliation works" is
  not a deliverable. "`uv run pytest tests/test_reconcile.py` runs the
  nine cases of §5.4 and §7.2, the invalid example and its valid variant
  among them, and all pass" is. A check fails on the commit the work
  starts from, because the thing is not there yet. pytest with a filter
  that selects nothing exits with 5 and a test file that is not there with
  4, but a check that only counts passes reads 0 as success, so a check
  names its tests or the number that must run, read from the summary
  line. A check that is true before the work checks nothing.
- **What it stands on**: the milestones before it, and anything outside
  the plan, a release of pvc, a tool that must be installed.
- **What could go wrong**, when something is known: the part of the spec
  that is thinnest, the case that pvc and the spec may disagree on. The
  session reads this to know where to look, and passes it to the
  reviewers.

A milestone whose failure would be silent, a rule that gives a plausible
wrong outcome, says which check guards it.

Milestones that touch different files and do not need each other are
marked as able to run side by side. The rest run in order.

## The shape of the document

1. The opening: what the plan builds and from which spec, the date, and
   the state, draft, approved by the owner, under way or done.
2. In and out: what is built, and what a reader could expect and is not,
   with where it goes. The open points of the spec that are not answered
   go here, each with the milestone its answer would change.
3. What has to be in place before the first milestone, in a form that can
   be checked: the Python version, the tools of the `coding` skill, the
   test data.
4. The milestones in order, numbered, each with the parts above and a box
   the session ticks: `- [ ] 3. Loci from two streams of VCF records`.
5. How the whole plan is checked at the end, when that is more than the
   sum of its milestones: a file against itself and against its
   re-spellings over a real region, as `doc/spec.md` §7 asks.

A plan says what and in what order, not how. It holds no code and no
signatures, which are in the spec, and no reasons of the spec, which it
points to.

## Before handing it over

- Every milestone names its part of the spec. A milestone with nothing of
  the spec behind it is a gap of the spec or work nobody asked for.
- Every deliverable has a check that can fail, and the check was run on
  the commit the work starts from and failed there, or the plan says why
  it could not be run.
- Read it for filler, sentence by sentence, with the test above.
- Could a subagent that reads only one milestone, the spec and the skills
  build it? Read two milestones that way.
- Nothing in the plan decides what the spec left to the owner.
- The plan goes to the `first-reader` subagent, as any document, read as
  the session that will run it, with questions of this kind: what is the
  first thing to do, how do I know milestone 2 is done, what do I do if
  pvc and the spec disagree.
- The owner approves the plan before it is run.
