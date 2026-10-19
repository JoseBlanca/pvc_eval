---
name: writing
description: How prose is written in pvc_eval, and how the session talks to the owner. Use it before writing or revising anything a person will read, a document under doc/, a spec, a plan, a report, a skill, a docstring, a commit message, a GitHub issue or a pull request, and the text pvc_eval shows its users, its --help, its warnings and its errors. Chat replies follow the same principles, and the ones that fail most often there are in CLAUDE.md.
---

# Writing

This skill is the writing skill of popnei_web,
`/Users/jose/devel/popnei_web/.claude/skills/writing/SKILL.md`, which
came from popnei's, with the readers, the examples and the forms of
pvc_eval. It was ported on 19 October 2026, when pvc_eval had a spec under
review and no code.

## The reader

What is written in pvc_eval is read by one of three people.

The first is the owner of the project, a population geneticist who
programs in Python and Rust and wrote pvc, the variant caller whose STR
finder pvc_eval ports. They know VCF, the callers, the benchmarks of GIAB
and the decisions recorded in `doc/`. They do not know what happened in
the session: the files that were opened, the commands that were run, the
names that came up along the way. So they need a little context before the
content: what was being worked on, what the question was, and where it
sits in the project. A sentence or two is usually enough.

The second arrives later, a contributor or another session of the
assistant. They may know Python well and VCF less, and they have only the
page.

The third is a user of pvc_eval, a bioinformatician who reads its
`--help`, its warnings, its errors and its summary in the middle of a
benchmark. Their own section is at the end.

The first two read in order to do something: to decide, to build, or to
check. A text is good when its reader gets through it once, without
asking a question and without opening anything else, and can then do what
they came to do.

`doc/spec.md` is the voice to match: facts and numbers, each rule with
its reason, terms defined once in a table and then used in that sense
only, headings that name their subject, and no sentence about the text
itself.

## Before writing: the sketch

The order of a text is decided before its sentences. A text whose order is
found while writing uses its terms before it defines them and gives its
reasons before its subject. So write a sketch first, a few lines that are
not part of the text:

1. The subject in one sentence, with each option or thing named by what it
   is: "whether a locus that holds more VCF records than the record cap is
   split, or reported whole as not compared".
2. What the reader should know or be able to do when they finish, in the
   words you would use with a colleague across the table. When that cannot
   be said yet, the problem is in the thinking and no sentence will fix
   it.
3. The terms the text will need, each with its definition, in the order
   they will be used: "a locus allele: the sequence one haplotype carries
   over the locus".
4. The points, one line each, in the order the reader needs them. A point
   that leans on another comes after it.
5. What is left out, and where it goes.

Then write from the sketch. The opening of the text is the first two lines
of the sketch, in full sentences. When the text is written, find the first
use of each term in line 3 and check that its definition comes before it,
as a sentence of the text. A meaning that the reader could work out from
what came before does not count as a definition.

## The principles

### Write the contents, not the name of the category

"The reconciler handles the genotype logic" gives the reader nothing to
use. "The reconciler builds the locus alleles of a locus and, for each
sample, the locus genotypes its VCF genotypes allow" does. Logic,
handling, machinery, infrastructure, plumbing and layer are names of bags.
Write what is in the bag. When the list cannot be written, the thing is
not yet understood.

A problem is written as what can go wrong, not as the rule it breaks: "two
VCF records that overlap by one base go to different loci, so a deletion
and the SNP inside it are compared apart and both come out as different
genotypes". That says what the user would see. "The chaining is not
transitive" says neither.

An adjective in the place of a number is the same fault: fast, slow,
large, heavy, most. The repair is to put the fact in: "the GIAB benchmark
regions of chromosome 20 take 41 s and 380 MB". Taking the adjective out
and leaving the rest leaves a sentence that still says nothing.

### A number comes with what it was measured on

The dataset, the machine, the Python version, the commit, and the date
when the measurement is not from the work being reported. The date is the
one the session was told, not one worked out from the dates already in the
files: in popnei, on 24 September 2026, sixteen dates one day ahead were
found, written by sessions that had each inferred the same wrong day. A
writer who does not know the date asks.

A comparison has both of its sides in the same units, and says which side
is the better one when the reader could doubt it. A number carried over
from another document keeps all of this.

A trade-off gives up one kind of thing to get another: memory for time,
lines of code for time. Each kind is named, and each option gets its value
in both. A word keeps one meaning through a sentence: "costs" that means
seconds in one half and lines of code in the other has two meanings.

### Context before the name

A term that does work in a sentence is explained before that sentence or
inside it. The reader lacks five kinds of names:

- Labels made up during the session: step numbers, codes of findings, the
  nickname of a function. The reader needs what the thing is.
- Names from the code. That a function is called `chain_loci` does not
  make it a word the reader has. Say what it does, and then the name can
  be used.
- Words of a tool the reader may not use: a fixture, a golden file, a
  sampling profiler, a lockfile.
- Ordinary words that mean something narrower here. A locus, a haplotype,
  missing and invalid mean what the vocabulary of `doc/spec.md` §2 says,
  and a text that leans on them uses them in that sense and no other.
- Notation, and a label in a table, which is a first appearance too: the
  words go in the text before the table.

One name for each thing. Before a new name is added, count the ones the
document already uses for it. The things of pvc_eval keep the names of
`doc/spec.md` §2, and the words its right-hand column lists are not used:
a VCF record, not a variant or a site; file A and file B, not truth and
query.

### What the reader came for goes first

In a document, what it is about and what it decides. In a section, what
the thing does, in words, before the type or the signature. In an issue,
the finding. In a reply, the answer to the question that was asked. The
work usually happened in the opposite order, the answer was found last,
and the temptation is to tell it in the order it happened.

Order carries importance. What matters most is first or has a paragraph of
its own, and then no sentence needs to say that it matters.

### Only what the reader can use

For each sentence, what can the reader do with it? These go out:

- The story of how the work went.
- Sentences about the text: "in one paragraph", "this section describes".
- Sentences about the reader's reaction: "surprisingly", "notably", "the
  key point is".
- Answers to an objection that nobody raised.
- Sentences that only give a verdict, "there is a catch", and leave what
  is good or bad to the next sentence. The sentence that holds the fact
  carries the verdict too.

The test is to delete the clause. When no fact is lost, it stays deleted.

Two things stay although they look like the story of the work: a trap the
next person would fall into, and a measurement that closed an option.
They stay where the next person will look for them, the document or the
issue about that subject.

A few texts have a second reader as well as their own, as the work report
of a plan does, read by the owner to decide whether to merge it and by
whoever next revises a skill. What the second reader needs goes in one
place at the end, under a heading that says who it is for, and the first
reader is told there that they can stop.

### Everything at its true strength

A choice is written as a choice, with the goals it serves and what would
have to be true for another option to win. It is not "forced" or "the only
way". A design is explained by its reason, and what is wrong with the
alternative is not yet a reason.

What was seen is kept apart from what was assumed or read somewhere. "It
works" means nothing until it says what was run: "the 214 tests pass, and
file A against itself over chromosome 20 gives only same genotype and
reference in both; not run on a whole genome". A claim that holds only in
part is given with the condition under which it holds, "for diploid
samples", "outside satellites". Somewhat, relatively and in general tell
the reader that there is a limit and hide where it is.

### Plain sentences that stand on their own

Short sentences with active verbs and ordinary words. A sentence is not
shaped to sound balanced or sharp. The plain version is usually longer,
and it is the one to write.

Every count has its noun and its set: "three of the seven outcomes",
never "three of the seven". Every pronoun has an antecedent the reader can
see, and a "this" that could point at two things gets its noun. After an
edit, read the sentence before and the sentence after, because a deletion
can take an antecedent with it.

A list for parallel items, a table for measurements and for cases, prose
for reasoning. A heading names the subject of its section. Bold marks the
terms that a list or a table defines and nothing else.

## The forms

- **Documents under `doc/`.** The opening paragraph says what the document
  is, its date, what it decides and where the related documents are.
- **Specs and plans** follow the `spec` and `planning` skills, which add
  what goes into them. Performance reports follow the `performance` skill.
- **Skills and agents.** Each rule with its reason, because a rule without
  one is followed where it does not fit and dropped where it does. A
  decision is dated, "decided on 19 October 2026", so that it can be found
  and revised.
- **Docstrings.** What the function or class is or does, in the words of
  `doc/spec.md`; the coordinates of every position, 0-based or 1-based,
  half-open or inclusive; what it raises and when. How it is implemented
  stays out unless the caller sees it. A docstring that repeats the name
  and the type hints is left out.
- **Comments.** Why, never what. A comment that names the rule of the spec
  it follows, "§5.3: touches, so the insertion before the repeat pulls it
  in", is the most useful kind.
- **Commit messages.** A lower case subject that says what changed, and a
  body with the why and the numbers. Nothing about the session.
- **GitHub issues.** The title states the finding or the task. The body
  says what was seen, on what input and how to see it again, what it means
  for the comparison, and what is proposed or asked.
- **Pull requests.** What changed for a user of pvc_eval first, then what
  changed in the code, how it was checked, and what is not done.
- **A request for a decision**, in any form. What is being decided, with
  the options spelled out; what each one gives and takes; the
  recommendation, in a full sentence; what is not known; what happens next
  in each case. It is ready when the reader can answer without asking
  anything back. What the writer can decide alone is decided and not
  asked.

## The text pvc_eval shows its users

The users are bioinformaticians who run pvc_eval on their own VCF files
and read what it prints to decide what to do next. The principles above
hold, and the principle of `doc/spec.md` §1 sets the tone: pvc_eval reports
what it found and does not judge which file is right.

- **`--help`** names each option in the words of `doc/spec.md`, with its
  default and its unit: "--record-cap N: the most VCF records, of both
  files together, a locus may hold (default 30)".
- **An error** stops the run. It says what was found, where, and how to put
  it right, with the names from the user's own files: "b.vcf.gz is not
  sorted: chr1:10432 comes after chr1:10511 (line 2204). Sort it with
  `bcftools sort` and run again." No traceback for a fault of the input, no
  "invalid", no apology, nothing of the code. A traceback is for a defect
  of pvc_eval, and then the message says it is one.
- **A warning** comes from the data and says what was found, what it does to
  the result and what the user can do: "3 samples of a.vcf.gz are not in
  b.vcf.gz and are not compared: HG002, HG003, HG004. Map them with
  --sample-map to compare them." A warning with nothing to do says so.
- **The summary and the table** use the names of `doc/spec.md` §5.5 and §6
  and no others, since a user reads the spec to understand them.

## Before handing a text over

The writer cannot see what is missing from the page, because they know it.
So a document, a spec, a plan, a report or an issue goes to the
`first-reader` subagent, which gets only the text, one sentence that says
who reads it, what for and which documents they already know, and three
to five questions that the reader should be able to answer once they have
read it. It returns what it understood the text to say and to ask, its
answers to the questions, the names it did not have, and the sentences it
could not follow.

Compare its summary and its answers with what was meant. Where they differ
the text is wrong, not the reader. Fix what it reports, and send the text
again when the fix was large.

A commit message, a docstring or a reply in chat does not go to the
subagent. Read it once more against "Only what the reader can use" and
"Context before the name".

## When a text is sent back

The text is corrected, and so is this skill. Find the principle that
allowed the failure, or that is missing, and revise it where it stands.
The case, the situation, the text that was sent and the owner's words, is
saved under `cases/` beside this file, with the date in its name. The
cases are for trying out changes to the skill and are not read when
writing. popnei's cases, which shaped this skill, are at
`/Users/jose/devel/popnei/.claude/skills/writing/cases/`.

The skill does not grow by one entry for each failure. A correction
becomes a new paragraph only when no principle covers it, and each
principle keeps one example. About 300 lines, its size when it was
ported, is the size to stay near, a figure to get oriented by and not a
limit.
