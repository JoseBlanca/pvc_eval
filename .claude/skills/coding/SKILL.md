---
name: coding
description: How code is built in pvc_eval, the implementation loop. Use it whenever the owner asks to implement, run, continue or fix a plan or a milestone, or to make any change to the code. The session coordinates and does not write the code: for each milestone a subagent implements it test first, reviewer subagents review it in parallel, one per category, the session judges the findings, and a subagent fixes the ones that hold. The rules of the code itself are in python.md beside this file, and the review categories in review.md.
---

# Building pvc_eval: the implementation loop

The code is built milestone by milestone, from a plan written as the
`planning` skill says. For each milestone:

1. a subagent implements it, test first;
2. the session checks what came back, by running the checks, not by
   trusting the report;
3. reviewer subagents review it, one per category, in parallel;
4. the session judges each finding;
5. a subagent fixes the findings that hold, each with a test that failed
   before the fix;
6. the session runs the checks again, has the fix of each blocker, a
   finding of the worst severity in `review.md`, looked at once more by
   its reviewer, and goes on to the next milestone.

The session goes on from one milestone to the next without stopping, and
stops for the owner only in the cases of "When to stop and ask".

## Files beside this one

- `python.md`: how the code is written, the uv commands, the checks, the
  tests. Every subagent of the loop reads it.
- `review.md`: the categories of the review, what each looks for, and the
  form of a finding.

## The session coordinates, and does not write the code

The session sends the implementation, the reviews and the fixes to
subagents, and keeps the plan, the state of the work, the decisions and
the report. The reason is its context: a plan takes many hours, the
details of every file would fill the session's context before the end,
and a coordinator that has forgotten the start of the plan decides badly
at its end. The second reason is the review: the writer of a piece of code
cannot review it, because they know what it was meant to do and read that
into it, and a model rates its own output above other output of the same
quality.

## Before the first milestone

1. The plan is approved by the owner, and its spec too. If not, say so and
   stop. A small change with no spec has the plan of one milestone the
   `planning` skill describes, which the session writes and does not wait
   on the owner for, unless it is one of the stops below.
2. A worktree and a branch, named after the plan, from the local `main`:

       git worktree add .claude/worktrees/<plan> -b plan/<plan> main

   then `EnterWorktree` with that path. All the work happens there and
   `main` is not touched. The reason is the one found in popnei, an
   earlier project of the owner's whose skills these come from, on 21
   September 2026: two sessions that edit one checkout at once leave their
   changes mixed in the same files.
3. `uv sync`, then the checks of `python.md` on the commit the work starts
   from, so that a failure found later belongs to the plan. A check that
   fails already is said in the report and not fixed inside the plan
   unless the plan needs it.
4. The work report, `doc/reports/<plan>.md`, started and committed with the
   plan's state changed to under way.

## 1. The implementation

The session notes the commit the milestone starts from,
`git rev-parse HEAD`, which the review and the report need. The milestone
goes to one subagent, `general-purpose`, on the session's
own model unless the owner has named another, with a prompt that stands
on its own, since the subagent knows nothing of the session:

- the path of the worktree, and that all its work happens there on the
  branch checked out. It is not given `isolation: "worktree"`, which would
  start it from `main` in a tree of its own;
- the path of the plan and the number of the milestone, which it reads
  there with the parts of the spec it names. The milestone is not copied
  into the prompt;
- the skills to follow: `.claude/skills/coding/python.md`, and
  `.claude/skills/writing/SKILL.md` for docstrings, messages and commits;
- today's date, and that every date it writes is that one. A subagent has
  no clock and works the date out from the files: on 24 September 2026, in
  popnei, one subagent working through a list wrote the same wrong date 14
  times;
- what the session knows and the code does not show: a decision of the
  owner, what an earlier milestone left for this one;
- to work test first, as `python.md` says, and to commit when the checks
  pass, in commits that each pass them; not to edit the plan or the
  report; a case it had to add to the spec is a commit of its own before
  the code;
- to stop and report, without working around it, when the milestone finds
  that the spec or the plan does not hold;
- what to send back, in under 300 words: what it built, its commits, the
  last line of each check, what it did differently from the plan and why,
  and any question it could not settle.

## 2. Checking what comes back

What a subagent says is a claim. A claim the next step rests on is
checked by its results, which costs little context:

- `git log` and `git status` in the worktree: the commits are there and
  the tree is clean.
- The checks of `python.md`, run by the session, reading their last lines
  and the counts in them.
- The deliverables of the milestone, each with the check the plan gives.
  A deliverable that does not pass is not done, whatever the report says.
- What the subagent did differently from the plan: is it small, or does it
  change what a later milestone builds on?

When the work is not right, the same subagent gets it back with what is
wrong, through `SendMessage`, since it has the context. After a second
failure a fresh subagent gets the milestone with what was learned. A
third failure is a stop.

## 3. The review

The review is of the milestone's commits, `git diff <start>..HEAD`, where
`<start>` is the commit noted before the implementation, and of what that
code calls and is called by.

Before sending the reviewers, the session gathers what each one is given:
the range of commits and the files in scope; the path of the plan and the
milestone, and the parts of the spec they build; the output of the checks,
as printed, so that no reviewer spends its time finding what ruff or ty
already print; and what the code does not show, an answer of the owner, a
later milestone that will add what looks missing. Most wrong findings come
from missing context, and a line in the prompt avoids them.

The categories are in `review.md`:

| category | sent when |
|---|---|
| `spec` | always: the code against the spec, case by case |
| `tests` | always: whether each test can fail, and the numbers the change claims |
| `correctness` | always: positions, streams, the traps of Python |
| `errors` | the change reads input, raises, catches, or decides what is invalid |
| `design` | always: whether the code is as simple as its job, its names and types |

When in doubt, send it. A reviewer with nothing to report costs little.
The subagent is `code-reviewer`, one per category, all in one message so
that they run in parallel. `spec` and `tests` run and change code, so they
are sent with `isolation: "worktree"`, each in a tree of its own, and told
to check out the last commit of the milestone and confirm it with
`git rev-parse HEAD`, because such a worktree starts on `main`; then
`uv sync`. The other categories read the session's worktree and do not
change it.

## 4. Judging the findings

The session judges each finding itself, and does not hand the owner a
list. A finding is taken seriously: the reviewer read the code without
knowing what it was meant to do, as its next reader will, and usually ran
the case. It is not set aside because it is inconvenient or because the
milestone was nearly done. And it is not the law: a reviewer lacks
context, misreads a line, or assumes an input that cannot occur.

1. Check it. Open the line, run the case, read the part of the spec it
   cites.
2. Decide what it is:
   - it holds: it goes to the fixer, with the session's verdict;
   - it holds and the fix is the owner's, because it changes an outcome,
     the table, the interface or the spec: it is a stop;
   - it holds and belongs to later work: a GitHub issue, or a line in the
     plan of the milestone it belongs to, so that it is not lost;
   - it does not hold: say why, in a line, with what was checked. A reason
     is one the reviewer would accept if it saw it.
3. A finding about a wrong outcome is never closed on a doubt: when it can
   be neither confirmed nor refuted, the reviewer gets the missing context
   and is asked again.

Two reviewers that report the same thing from two sides are evidence, not
repetition: handle it once and count it as more certain.

## 5. The fixes

The findings that hold go to one subagent, the one that implemented the
milestone when it is still there, since it has the context, with:

- each finding as the reviewer wrote it, its number, and the session's
  verdict;
- to fix only those, with the smallest change that settles each one at its
  cause, not at the input that showed it;
- for a defect, a test first that fails on the code as it is, then the
  fix;
- not to weaken a test, mark one skipped or xfail, or add a `noqa` or a
  `ty: ignore` to make a finding go away; a test that has to change says
  why in the commit message;
- not to touch code that no finding names;
- that it may answer that a finding does not hold, with evidence, instead
  of changing the code;
- one commit per finding, or per group of findings with one cause;
- to send back, for each finding, one of: fixed, with the cause and the
  name of its test; disputed, with the evidence; deferred, with the
  reason; and the last line of each check.

A disputed finding is weighed by the session as it weighs a reviewer.

## 6. After the fixes

The session runs the checks again. A fix of a blocker goes back to the
reviewer of its category, scoped to the diff of the fixes: does it settle
the finding, does its test fail without it, and does the fix break
anything else. This second look covers nothing else, so the loop ends.
When it finds the fix wrong, the fixer gets the finding again: at most
two rounds of fixes per milestone, and a finding still open after the
second round of fixes, or one that comes back after it was fixed, is a
stop.

When a finding shows that `python.md`, `review.md`, the spec or the plan
allowed the defect, that document is corrected too, so that the next
writer knows the rule before the next reviewer finds it.

Then the session ticks the milestone in the plan, writes it into the
report, commits both, and goes on.

## Changing the plan

A plan is written before the work and the work knows more. The session
changes the plan without stopping when the change is small: a milestone
split in two, two merged, a file other than the one named, a check
replaced by one at least as strong. Every change goes into the plan and
into the report with its reason.

## When to stop and ask

- A change to what pvc_eval gives: an outcome, the table, the summary, the
  command line, the Python interface, or a rule of the spec. The session
  does not answer an open point of the spec; it follows its "meanwhile".
- A runtime dependency, a new Python version, or anything that is not pure
  Python.
- A deliverable dropped or its check made weaker, a milestone added or
  removed.
- pvc, the VCF specification and the spec contradict each other in a way
  that makes a milestone's approach wrong.
- A milestone that failed three times, or a finding still open after two
  rounds of fixes.
- Anything outside the worktree: a push, a merge, an issue in another
  repository.

Before stopping, do the work that does not depend on the answer. Leave the
tree with the checks passing and everything committed, bring the report up
to date, and ask as the `writing` skill says a decision is asked for: the
options, what each gives and takes, and a recommendation.

## The work report

`doc/reports/<plan>.md`, written as the work goes and not at the end, for
the owner, who was not there. After each milestone, a short section:

- the deliverables, each with the command that checked it and what it
  gave;
- what was changed in the plan, and why;
- what the review found that mattered and is fixed, and each finding not
  taken, with its reason in a line;
- what the owner should know: a case the spec did not foresee, a
  disagreement with pvc or bcftools, a risk for what comes next.

At the end, a section for whoever next revises these skills, with a
heading that says so: what would change a skill, the shape of a plan or
the size of a milestone. A milestone that had to be sent twice and why; a
rule that misled, with its sentence; the tokens the implementation, the
review and the fixes of a milestone each cost, from the results of the
`Agent` tool, side by side, since that is what sizes the next plan's
milestones. Not a log of what happened.

When the plan is done, the top of the report gets what the owner reads
first: whether it is done, what exists now that did not, what is left
open, and what is asked of them, at least the merge. The report goes to
the `first-reader` subagent, as the `writing` skill says, before the
owner is told.

## When the plan is done

1. Every milestone ticked, every deliverable checked, the check of the
   whole plan run, the last part of a plan as the `planning` skill shapes
   it, the state of the plan changed to done, the report
   finished, all committed on the branch.
2. Tell the owner, in a reply as CLAUDE.md describes them: the plan is
   done, on branch `plan/<plan>`, with the report's path and what it asks
   of them.
3. Nothing is merged into `main` and nothing is pushed without the owner's
   order. After the merge, `ExitWorktree` with `keep`, then
   `git worktree remove .claude/worktrees/<plan>` and
   `git branch -d plan/<plan>`, and the worktrees of the reviewers.

## When a session ends before the plan

The state is in the repository: the ticks in the plan, the report, the
commits on the branch. A new session enters the worktree with
`EnterWorktree`, reads the plan and the report, runs the checks, and goes
on from the first milestone not ticked.
