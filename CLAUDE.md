# pvc_eval: how the assistant works here

pvc_eval compares two VCF files over a reference genome by rebuilding, for
every sample, the sequence each haplotype carries over each locus, so that
two spellings of the same genotype agree. What it does, its vocabulary and
its rules are in `doc/spec.md`, which is the source of the names used
everywhere in the project: a VCF record, not a variant; file A and file B,
not truth and query.

It is pure Python, run with uv and nothing else, on the free-threaded
CPython 3.15t, which the owner chose on 19 October 2026 so that its
threads run in parallel. A runtime dependency is the owner's decision.

The skills are under `.claude/skills/` and the subagents under
`.claude/agents/`:

- `writing`, read before anything a person will read is written: a
  document, a spec, a plan, a report, a docstring, a commit message, an
  issue, and the text pvc_eval prints. A document goes to the
  `first-reader` subagent before it is handed over.
- `spec`, for a spec, reviewed by `first-reader` and `spec-reviewer`,
  and approved by the owner before any code is made from it.
- `planning`, for the plan that turns an approved spec into milestones.
- `coding`, the implementation loop: the session coordinates, a subagent
  implements each milestone test first, `code-reviewer` subagents review
  it, one per category, and a subagent fixes what holds. Its
  `python.md` is how the code is written, and every subagent that touches
  code reads it.
- `performance`, for wall time and peak memory: measure, profile with
  Tachyon (its `tachyon.md`), send `perf-reviewer` subagents, then try each
  candidate and keep what the numbers and the code justify.

A session that writes anything in the repository, a spec, a plan, a skill
or code, works in a git worktree and a branch of its own, under
`.claude/worktrees/`, made before its first edit: two sessions that edit
one checkout at once leave their changes mixed in the same files, as
happened in popnei on 21 September 2026. Nothing is merged into `main`,
and nothing is pushed, without the owner's order.

## Replies in chat

The owner reads a reply to decide something or to learn something they
need for their own work. They did not see the session. The principles of
the writing skill hold in chat, and these are the ones that fail most
often there:

- A reply opens with a sentence that says what it is about, also when it
  follows a question: what was asked or what was being worked on, with the
  options or the things named by what they are. The owner may come back to
  it after hours of other work. The answer comes next, in a full sentence.
  When several things are true, the one that changes what the owner does
  goes first.
- A longer reply is sketched before it is written, as the writing skill
  says, so that no term is used before the sentence that explains it.
- A reply holds the decisions that are needed from the owner and what they
  need to know. How the work went stays out, unless it changes what they
  do next.
- A reply carries the numbers the answer turns on and leaves the rest for
  a document, with a line that says where it is or that it can be
  written. A reply that needs headings is usually a document.
- A request for a decision has the options, what each one costs, and a
  recommendation. What can be decided without the owner is decided, and is
  mentioned only when they need to know of it.
- No name before its explanation. Labels from the session, and the names
  of files, functions and tools, are not words the owner has, unless the
  owner used them first.
- A number where an adjective would go, with what it was measured on. A
  comparison has both sides in the same units and says which is better.
- Nothing about the reply itself, and nothing about how interesting or
  important a thing is.
- What failed, or was not done, or was not run, is said as plainly as what
  worked.
