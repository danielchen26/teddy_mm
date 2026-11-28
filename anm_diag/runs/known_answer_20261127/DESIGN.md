# A check with a known answer: does a knockout compose with IFN-γ by addition?

27 November 2026. Run 51. Asked by the owner on 27 November. **Nothing is trained.** This document is committed before
any value of the test is read. The scripts are rehearsed on made-up values and committed before the reading. The
reading is made once and no rule is changed after it.

## 1. The question and the known answer

An event is one knockout in one condition (no signal given, IFN-γ, co-culture). The composition under test is
**knockout × IFN-γ**: is the effect of a knockout under IFN-γ the effect of the same knockout with no signal given,
plus the effect of IFN-γ?

The known answer: knocking out **IFNGR1, IFNGR2, JAK1, JAK2 or STAT1** removes the cell's response to IFN-γ, so for
these five the effects do not add. For the great majority of the other knockouts they should.

What this is and is not. Every cell used here was read before: the 20 proteins of these knockouts in batches 1 to 4
and in the final audit, and the five genes were a declared type of the frozen system from the start. The RNA of these
cells was never read gene by gene. So this is a check of a method against an answer known beforehand. It is not a
finding on data never read, and it does not count as a confirmatory run.

## 2. The units

The knockouts whose cells were read before in both conditions: **95** (the five among them), with 100 to 341 cells per
knockout and condition. Control cells: 3,978 with no signal given, 4,000 under IFN-γ. Their names and counts are in
[`results/units.json`](results/units.json), written before any value is read.

Not opened: the 188 units never read (they are kept for a later test), and the 73 reserve donors. The 65 knockouts
read under IFN-γ only are therefore not in the table.

Secondary, described and with no pass line: knockout × co-culture, on the 56 knockouts read both with no signal given
and in co-culture (5,897 control cells there), by the same rule.

## 3. The readouts

All of them, not PD-L1 alone.

| | per group of cells | which |
|---|---|---|
| RNA | the natural log of (mean counts per 10,000 + 0.1) | every gene whose mean among the control cells is at least 0.5 per 10,000 in either condition |
| protein | the mean of the centred log-ratio over the 24 antibody columns, as in every audit | all 20 |

Logs, because a knockout and a signal that act independently multiply. Left out for a knockout: its own gene, and the
pairs the graph declares as direct (the 14 declared edges, such as CD274 → PD-L1). Removing a gene removes its own
product in any condition, by an amount that depends on how much there was; the graph already holds these pairs with a
size per condition. They are judged by the same rule and reported in a column of their own.

## 4. The rule (question 1: the data alone, no model)

**Halves.** Every cell goes to half A or B by a hash of its name (`sha256("half|" + name)`, even or odd), control
cells too. Nothing is ever computed across the halves except the two products below. So the cells that say where a
knockout departs are never the cells that say whether it does.

**The departure from addition.** In each half, for knockout *g* and readout *r*:

    I = (g − control, under IFN-γ) − (g − control, no signal given)

It is zero if the effects add. For the genes it is centred by its median over the genes of the group: RNA levels
are relative, and a knockout that removes an abundant transcript moves every other gene by one common factor, which is
no departure. It is divided by its error from cell counts (the variance among the control cells of
each condition over the number of the knockout's cells), centred by its median over the knockouts and scaled by its
median absolute deviation over the knockouts: *z*. **The unit is the knockout**: what noise does to a readout is
learned from the 95 knockouts, not from cells.

**Per readout**: the smaller of |*z*(A)| and |*z*(B)|, with the sign of their product. Beside it, over all
readouts together, the sum *T* of *z*(A) × *z*(B): the whole profile.

**The interval comes from the mirror side.** If the effects add, the two halves err independently: a product is as
often negative as positive. A real departure can only make it positive. So the negative side shows what noise alone
does, whatever its shape:

- *w* = the largest value of any readout of any knockout whose halves disagree in sign, and not less than 4 (the four
  errors of every other rule here);
- τ = the largest −*T* of any knockout.

A knockout **carries over** if no readout that agrees in sign exceeds *w*. Otherwise it is **refused** (does not
add), and its **witnesses** are the readouts that agree in sign and exceed *w*: a refusal always leaves a witness.
Whether *T* exceeds τ is reported beside the verdict and decides nothing (below). With the witnesses, the share of the
IFN-γ response that the knockout removes (the departure of one half projected on the IFN-γ effect of the other half,
over the readouts IFN-γ changes two-fold; described).

**Resampling.** The 95 knockouts are resampled 2,000 times (seed 20261127); the medians, the deviations, τ and *w*
are computed again each time. Reported: the range of τ and of *w*, and for every knockout the share of resamples in
which it is refused. The verdict is the one of the full set.

By the argument of the mirror about one knockout that adds is expected to be refused, among 90.

**Changed after the rehearsal on made-up values, before any value was read** (the first version of this document is
commit `f1578f4`). Two things, both found on made-up worlds:

1. *The centring over the genes* was not there. In a made-up world with planted knockouts, the planted ones shifted
   every gene the signal does not touch by a common factor (counts per 10,000 are shares), and the sum *T* then
   refused 22 knockouts that add.
2. *T was a second way to be refused.* On made-up worlds with nothing planted it refused 0 to 6 of 40 knockouts in a
   test (8 in six tests; no readout refused any): its end is a single most negative value and varies widely, and a
   refusal through it names no readout. The owner asked that a refusal leave a witness. So the verdict is by readouts,
   and *T* is reported.

The pass line is unchanged.

## 5. Types, from the witnesses

No type is declared by hand. The refused knockouts are grouped by the rule that proposes types (`events.propose_types`,
unchanged: cosine at least 0.5, at least two members): each refused knockout is a flag, and its direction is its *z*
(mean of the halves) over the readouts that are a witness of any refused knockout. A group of two or more is a
class. Only then are the classes compared with the five.

## 6. The pass line

| | holds if |
|---|---|
| **P1** | all five of IFNGR1, IFNGR2, JAK1, JAK2, STAT1 are refused (recall 5 of 5) |
| **P2** | at least 90% of the other 90 knockouts carry over (at most 9 refused) |
| **P3** | one class holds all five, and at most two other knockouts (precision at least 5 of 7) |

**The check is passed if P1, P2 and P3 hold.** Otherwise it is not, and which line failed is reported with the
knockouts concerned. Reported in any case: recall and precision of the refusals and of the class against the five.

What the author expects, written before any value is read: P1 holds 0.95; P2 0.75 (the centred log-ratio couples the
proteins: a knockout that removes an abundant protein moves the others by arithmetic, more so under IFN-γ; B2M and
TAPBP may be refused for that); P3 0.7; passed 0.55. In the co-culture: the five refused, 0.8.

## 7. Question 2: what TEDDY carries of the part that does not add

From what is cached, and no new run of the encoder: the states of all these cells, the head as trained, the RNA
counts.

For every refused knockout (if more than twelve are refused: the members of the classes, then the strongest
readouts, up to twelve), and for PD-L1, MHC-I and HLA-E, every other protein that is a witness, and its ten strongest
RNA witnesses: the account (`account.account`, unchanged) in each condition, with the same halves, at the places
below. Read at Z beside them, `others`: what a readout taught by the other members of the knockout's class reads.

| place | what is read | number |
|---|---|---|
| the input | V: 512 components of the log counts of the cells of the condition | `carried` |
| Z | the state: the mean of the last layer | `carried` (and `read`: what a readout of control cells reads) |
| the head's output | the output as trained; it exists for PD-L1, MHC-I, HLA-E only | `as_trained` |

What a place holds of a change is its share times the change. A change of less than four errors has no account, as
everywhere here: it and what is held of it are taken as zero. **The share of the part that does not add held at a
place = (held under IFN-γ − held with no signal) / (change under IFN-γ − change with no signal).** Per pair, and
pooled over a class (sums over its pairs; interval by resampling the knockouts of the class 2,000 times).

For an RNA readout the input holds the gene itself; V is a view of the input and is reported as what it is. The head
has no output for a gene. This part describes; no pass line rests on it. Expected: the head's output holds little of
it, the state a part, the input most.

## 8. Order

1. This document and `results/units.json`: committed and pushed.
2. On the owner's machine, the RNA counts of the IFN-γ cells of the 95 knockouts and of the 4,000 control cells are
   exported (`scripts/counts.py`; no value is printed). The other counts, the states and the antibody columns are
   already cached.
3. The scripts, rehearsed on made-up values with planted knockouts that do not add
   ([`results/rehearsal.json`](results/rehearsal.json)): committed and pushed.
4. The reading, once: `summarize.py`, `analyze.py`, `teddy.py`.
5. The write-up, as it comes.

Not done here: training of the pooling or of the last layer, the five proposed rule changes, a second model.
