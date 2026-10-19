---
name: spec
description: How a spec is written in pvc_eval. Use it when writing or revising doc/spec.md or a spec under doc/specs/, which says what a part of pvc_eval has to do, the rules that decide each result, its Python interface and how it is verified. A spec is written before the code of its part, reviewed, and approved by the owner; the plan and the tests are made from it. The prose follows the writing skill, which is read first.
---

# Writing a spec

A spec says what one part of pvc_eval has to do and how we will know that
it does it. It is written before the code of that part, and the plan and
the tests are made from it. It does not give the order of the work, which
is the plan's (the `planning` skill), and it does not repeat another spec,
which it points to by section.

pvc_eval has one spec for the whole program, `doc/spec.md`, and will have
smaller ones under `doc/specs/<name>.md` for parts that need more than a
section of it: the wrappers of other evaluators of `doc/spec.md` §8 are the
first one it names. A smaller spec is written when a part needs more than
about a page in `doc/spec.md`, or when it is built and reviewed apart.

The prose follows the `writing` skill. Read that first. This skill says
what goes into a spec, in which order, and what stays out. It comes from
the writing-specs skills of popnei and popnei_web, ported on 19 October
2026.

## Why a spec here

pvc_eval's whole worth is that its outcomes are right: a user reads
*different genotype* and trusts that the two files disagree. A wrong rule
does not crash; it gives a plausible table with wrong rows, and no user can
see it. So the rules are written and reviewed before the code, every rule
has a test, and the code follows the spec. When the code finds the spec
wrong, the spec changes first.

## The two readers

The owner reads a spec to check that the right thing is going to be built
and to decide the points that are theirs. The implementer, a subagent with
nothing but the spec, the skills and the code, reads it to build the part
without having to work out again what a VCF allows or what a rule meant.
Every paragraph is for one of the two.

## No filler

A spec has as many lines as it has things its two readers need, and none
besides. A spec that the owner does not read to the end has failed,
however correct it is, and what stops them is thoroughness for its own
sake: every edge named, every part filled.

The test for a sentence: would an implementer who has read `doc/spec.md`
and the VCF specification get this wrong without it? Would the owner
decide differently without it? When both answers are no, it goes out. The
parts below are the order of a spec and not a form to fill. A part with
nothing to say is a line or is left out.

## What a spec's claims are checked against

pvc_eval has no older program whose output is its specification. Its
claims are checked against:

- **The VCF specification**, version 4.2 and 4.3 (samtools/hts-specs), for
  what a VCF record, a GT, a phase set or a `*` allele means. A claim of
  that kind names the section.
- **pvc's code**, for the STR finder that pvc_eval ports, at the commit of
  `doc/port_provenance.md`. Where pvc's documents and its code disagree, the
  code is followed, as `doc/spec.md` §5.1 says. A claim about it names the
  file and the function, and is checked by running pvc when that is
  cheap.
- **What bcftools and RTG vcfeval do**, for normalising and for the cases
  where pvc_eval is meant to agree or to differ with them. Run them; do not
  recall them. bcftools is at `/opt/homebrew/bin/bcftools`.
- **Worked cases by hand.** A rule is shown on a small case, a few bases of
  reference and two or three VCF records, with the result written out. The
  same case becomes a test with its result as a literal.
- **The properties the design rests on**, which no outside program checks:
  a file against itself gives only *same genotype* or *reference in both*,
  and so does a file against a re-spelling of itself (`doc/spec.md` §7).
  These are written as checks of the spec.

## Before writing

1. Read `doc/spec.md` whole, and the specs the part depends on.
2. Read what the part stands on: pvc's functions it ports, with their
   tests; the sections of the VCF specification it reads; the code of
   pvc_eval that exists.
3. Write down what a reader would not expect: a case the VCF specification
   allows and callers write differently, a default of pvc, a case bcftools
   handles in a way pvc_eval must not copy. These are what the spec is most
   needed for.
4. Make the sketch the writing skill asks for.
5. Ask the owner, in a reply in chat, what would rewrite the spec if it
   were answered after the writing: what the part offers that other parts
   or the user will call, a line each, and any open point whose other
   answer would change more than a sentence. Each is asked as the writing
   skill asks for a decision. An open point that moves one value stays in
   the spec with its "meanwhile". When there is nothing to ask, say so in a
   line and go on.

A claim about what another program does comes from its code or from
running it, and names the program, its version and the command. A claim
that was reasoned out is checked or left out.

## The parts of a spec

**The opening.** What the part gives to the user or to the rest of
pvc_eval, in a paragraph; the date; the status (draft, under review,
approved by the owner on a date); that there is no code yet or what there
is; and the specs it depends on.

**The vocabulary**, when the part adds terms: a table of each term, its
meaning, and the words not to use for it, as `doc/spec.md` §2 does. A term
already in `doc/spec.md` §2 is not defined again.

**What it does.** In words first, what a user would see go right or wrong
because of it, and then the rules, each with a worked case where a reader
could get it wrong. A rule says what happens at each boundary: a position
at the edge of a region, a record that touches and does not overlap, the
first and the last base of a contig, an empty file.

**The Python interface.** The functions and types another part or a user
calls, as signatures in a code block, with a sentence before each one:
what goes in, what comes out, in which coordinates, and what it raises.
Streams are iterators, as `doc/spec.md` §4 says. Nothing private.

**The cases** where the result is not what the rules would make a reader
expect: a multi-ALT record whose ALTs move apart, a `./1`, a phase set
that spans two loci, a reference with an IUPAC code where the VCF has a
plain base, a VCF with no samples.

**How it runs**, when there is something to say: what it holds in memory
and how that grows, and how long it takes on what input, when the spec
has a target.

**How it is verified.** Each check with the function it is made at, the
highest one at which the thing can be seen, because a test pins the
signature it calls. The worked cases with their results as literals; the
properties; and, where another program gives a number, the input, the
command and what it printed, written into the tests as literals.

**Open points** and **Not in this spec**, as below.

## Decided, inherited and open

Every statement in a spec is one of three kinds, and the reader has to be
able to tell which.

Decided: the spec states it, with the reason when there was a real choice.
What the owner decided says so, with the date and the option not taken.

Inherited: it comes from the VCF specification, from pvc, or from another
spec of pvc_eval, and the spec says from where. "An STR bundle joins
repeats within 15 bases of each other, as pvc's STR finder does", with the
file and the function of pvc that sets the 15, is a complete entry, and
tells the reader which values can move and where. The name of a constant
or a function is copied from the code, never recalled.

Open: the owner decides. Where it comes up, the text gives the fact and
marks it, "(**Open 2**, below)", and says no more. The list at the end has
each point once, as a request for a decision in the form the writing skill
gives, with what the implementer does meanwhile. The points are numbered
through the whole spec.

What the writer can decide alone, a name, the shape of a private type, is
decided and not listed. Few open points, each one worth the owner's time:
what changes an outcome, what a user sees or is told, or the interface a
user calls.

**Not in this spec.** What a reader could expect to find here and is
somewhere else or will not be built, with where it goes. A few lines.

## Before handing it over

The spec goes to the `first-reader` subagent as the writing skill says,
once as the owner and once as the implementer, with questions for each.
For the implementer: which outcome does this case give, at which function
is the first test made, what is raised when the GT is malformed. For the
owner: what is asked of me, and what does each answer lead to.

Then read the list of open points against the text: every **Open** of the
text is in the list, and nothing in the list is decided somewhere in the
text. Read it for filler as the owner would, who has a day of other work.

Last, the review. The code is made from the spec, so a wrong sentence in
it becomes wrong code, and the writer who misread the VCF specification
will not find the misreading by reading the spec again. The spec goes to
the `spec-reviewer` subagent, with its path and the files it relies on.
Evaluate each finding before acting on it, since the reviewer can be
wrong, and when the spec is handed to the owner, say which findings were
not taken and why, a line for each.

A change made after the review is checked against every other place that
speaks of the same thing: search the spec for the name and the number.

The owner approves the spec. Its opening then says that it is approved,
with the date, and only then can a plan be written from it.

## Once there is code

The spec changes first. A changed decision, an answered open point, a case
the code found, goes into the spec in a commit of its own, and the commit
that changes the code comes after it, so the spec is never behind the code.
A change that alters an outcome or the interface is the owner's, and is
asked before it is made.

## When a spec is sent back

As in the writing skill: the spec is corrected, the principle here that
allowed the failure is revised, and the case is saved under `cases/`
beside this file.
