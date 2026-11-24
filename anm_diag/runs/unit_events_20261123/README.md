# Events without control cells in the round

23 November 2026. Step 1 of the plan in [`GENERAL.md`](../../GENERAL.md). A knockout comes with control cells; a
donor, a batch or a disease does not. The audit tool could read such events, and the round could not: its admission,
its repair and its readings all rested on control cells. This run puts them into the general code
(`anm_diag.cycle`, `anm_diag.outer`, `anm_diag.scope`), checks the result on planted systems where every answer is
known, and then runs one round on eight real donors. The design, the rule by which the planted systems are scored and
the forecasts for the real round were committed before the scored runs and before the round read anything
([`DESIGN.md`](DESIGN.md)).

Everything on real cells is exploratory. The donors are development donors that earlier runs have read. No final
unit and no reserve donor was read, no state was computed and no model was trained.

## What came out

**The round now runs on a unit read against the other units, and on a group of units read against another. On
planted systems it says nothing false. On eight real donors the repair of the head removes a good part of what the
head loses, and for most pairs it needs the donor's own record: it does not reach a donor it has not seen.**

### Planted systems

Two systems with a known answer (`anm_diag/toy.py`): a study of units in which no cell is a control cell, with nine
units that differ in a known place, and a disease, patients against healthy donors. Each also in a version with
nothing wrong. Twelve seeds at three sizes, scored by the rule of DESIGN.md, Section 5:

| size | ordinary units, cells per unit | healthy and patient donors, cells per donor | right | undecided | wrong |
|---|---|---|---|---|---|
| large | 72, 300 | 12 and 8, 400 | 192 | 0 | 0 |
| medium | 40, 200 | 10 and 6, 200 | 191 | 1 | 0 |
| small | 20, 100 | 8 and 5, 80 | 112 | 80 | 0 |

By planted fact (right / undecided / wrong, of twelve seeds):

| planted fact | large | medium | small |
|---|---|---|---|
| h1: a loss in the head, closed by the repair | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| h2: a loss in the head, closed by the repair | 12 / 0 / 0 | 12 / 0 / 0 | 4 / 8 / 0 |
| h3: a loss in the head, closed by the repair | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| solo: a loss in the head, nothing placed before the model | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| early: upstream | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |
| enc: upstream | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |
| modb: unreadable | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |
| post: unreadable | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |
| early: lost in the first step | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |
| enc: lost in the second step | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |
| solo: the rows of the other units do not carry | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| weak: refused at selection | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| no other definite item | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| units, nothing wrong: nothing is open and no repair is read | 12 / 0 / 0 | 11 / 1 / 0 | 3 / 9 / 0 |
| disease: closed by the repair | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 |
| disease, nothing wrong: nothing is open | 12 / 0 / 0 | 12 / 0 / 0 | 3 / 9 / 0 |

**576 scores: 495 right, 81 undecided, none wrong.** At the two larger sizes every planted loss is put where it was planted:
in the head and closed by the repair (h1, h2, h3), in the first or the second step of the model (early, enc), outside the
model's input (modb, post). The unit that only its own row reaches (solo) is flagged in all 36 runs and no other unit
ever is. The unit that barely differs is refused at admission in all 36. The planted disease is read as a pool of
patient donors and closed by the repair in all 36.

What shrinks with the data is what the loop is willing to say. With 20 ordinary units of 100 cells a loss before the
model is placed in 3 of 12 runs. In the other 9 no reading of the family is supported, the positive control fails,
and the loop says `invalid` and places nothing. The systems with nothing wrong are undecided in 9 of 12 runs at that
size: the interval of a single unit does not reach the margin, and the loop leaves an item `unresolved`. The property
that was asked for holds: `wrong` stays at zero as the data shrink.

### Eight donors

PBMC of eight healthy donors (Hao et al. 2021), 112 surface proteins, 26,706 cells read (the selection third is read
no further). The head was fitted on bone marrow and has never seen this study. A donor is read against the other
seven, its cells by cell and the other donors by unit. *u* is the share of the measured change that is left
unexplained: 0 is none of it, 1 is all of it.

| | pairs of a donor and a protein |
|---|---|
| admitted on the selection third | 26, for 26 of the 112 proteins (donor P8: 9, P6: 8, P7: 4, P2: 2, P5: 2, P4: 1) |
| the output as trained is rejected | 22 of 26 (median \|*u*\| 0.58) |
| a reader of Z with the donor's own cells: supported / unresolved / rejected | 0 / 11 / 15 |
| the head has a step to remove (A3) | 16 |

For the 16 pairs in which the head has a step to remove:

| | |
|---|---|
| the output as trained, median \|*u*\| | 0.85 |
| the repaired output (`rows`), median \|*u*\| | 0.40 (the median fall is 0.43) |
| grade: relative / absolute | 7 / 0 |
| grade: not tested (gate 3) / loss no longer resolved | 6 / 3 |
| the donor's own record is needed (`rows` reads a step more than `unseen`) | 10; not needed 5; unresolved 1 |
| the same repair without the donor's cells (`unseen`), median \|*u*\| | 0.66 |
| a readout with no record of any donor (`controls`), median \|*u*\| | 0.67 |
| fitted on every cell (`cells`), median \|*u*\| | 0.34 |

- **The repair works where the head is what loses.** In 16 pairs the output as trained leaves a median of 0.85 and
  the repaired output 0.40. Seven reach the relative grade: the repaired output adds no loss beyond a reader of Z that
  had the donor's own cells (PD-L1, CD47, CD103, CD107a, CD95, integrin β7, CD49a).
- **Nothing is absolute, because Z is not.** For no pair is a reader of Z supported, even with the donor's own cells
  (median 0.31 among the 16). What a donor differs by is only partly in the state the head reads.
- **The repair does not reach a donor it has not seen.** Without the donor's cells and record the repair leaves 0.66,
  the same as a readout that has no record of any donor (0.67). The records of the other seven donors say nothing
  about how the eighth differs. This is not particular to the repair: a reader of Z fitted on the other donors leaves
  a median of 0.57, against 0.31 with the donor's own cells. The outer step raises ten items "pooled, and differs"
  (lever S): the type table treats the donors as one type, and the rows do not carry within it.
- **Nothing is placed before the model.** For 24 of the 26 proteins no reading of any stage and scope is supported,
  so the positive control of the input view fails and the items are `invalid`, not `unreadable` or `upstream`. The
  two proteins for which it holds are CD16 and CD20. The planted systems do the same when units are few (9 of 12 runs
  at the small size, which has 29 units; this study has eight). The permuted view reads nothing for every pair (least \|*u*\|
  0.999), as it must.
- **Donor and batch vary together.** Donors P1 to P4 are batch 1 and P5 to P8 batch 2; 23 of the 26 admitted pairs
  are donors of batch 2. The outer step raises the item "confounded" (lever E) for the study, and nothing here is attributed
  to a donor alone. The other side of every event has seven units, so every interval is approximate by the protocol's own
  number.

## What it means

- **Step 1 is done for the shape "a unit against the other units"**: admission, a repair that needs no control cell,
  the readings, the placing and the items are in the general code, checked on planted units and run once on real
  ones. The shape "a group of units against another" (a disease) is in the code and checked on planted cells only.
  The one disease study at hand is the final study of round 1 and was not read.
- **What does not carry is the same thing as across kinds of cell.** Across kinds of cell the rows of the events
  carried and the readout of control cells did not ([`outer_step_20261122`](../outer_step_20261122/)). Across donors
  the part that takes the place of that readout (the variation among the cells of a unit) carries by construction,
  and what is specific to a unit, its row, does not. In both cases the repair is made for what it has seen.
- **For the plan**: a repair that holds for an unseen donor cannot come from this plug-in. It needs the state Z to
  carry what makes a donor differ, which is a matter of training (step 5), or records of many more units. The 73
  reserve donors of the COVID-19 study are the place where this could be tested once; opening them is the owner's
  decision.
- Nothing here is a test. The donors have been read before, and there are eight of them.

## The forecasts

Intervals: 10 of 16 held the value. Outcomes: for 2 of 2 the class given the highest probability happened
([`results/outcome.json`](results/outcome.json)).

| missed | forecast | observed |
|---|---|---|
| admitted pairs | [40, 110] | 26 |
| proteins with an admitted donor | [35, 80] | 26 |
| the output as trained, median \|*u*\| | [0.60, 0.97] | 0.58 |
| a reader of Z with the donor's own cells is supported | [0.05, 0.35] of the pairs | 0 |
| grade absolute, among the pairs with A3 | [0.02, 0.25] | 0 |
| a rejected reference is read as lost at the input view | [0.40, 0.90] | 0 (the positive control fails, so nothing is placed) |

The forecasts were built from a reading of all 112 proteins at once. Protein by protein, with the other side counted
by unit, fewer pairs pass admission and no reference reaches the margin.

## What was not as written

1. **One test assertion was wrong, and the package was not.** The tests were run after the commit of the code, as
   that commit said. One assertion asked that the planted disease's output as trained have a wider interval when the
   healthy side is counted by donor. That output does not follow a donor's level, so its interval need not widen. The
   property is now tested on an output that does miss the level of its donor (the same value, an interval 5.7 times
   wider). No line of the package changed.
2. **The scored runs of the planted systems were started twice.** The first start split the runs unevenly over two
   processes and was stopped after four minutes, when three runs had finished (all of the version with nothing wrong,
   all scored right). The second start made every run again with the same seeds. The runs are deterministic.
3. **Boosting was not used for the real round**, as the design says: 112 proteins times 36 fits is out of reach. A
   reading of the input view here rests on linear readers.

## Files

| file | what |
|---|---|
| [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json) | the design, the rule for the planted systems, the real round and its forecasts |
| [`results/planted.json`](results/planted.json) | every scored run of the planted systems, and the table above |
| [`results/selection.json`](results/selection.json) | admission of every donor for every protein on the selection third |
| [`results/round.json`](results/round.json) | the round: every pair, every arm, the placing, the targets, the open items, the summary |
| [`results/outcome.json`](results/outcome.json) | the forecasts against what was read |
| `scripts/planted.py` | the scored runs of the planted systems |
| `scripts/build.py`, `round.py`, `score.py` | admission and arrays; the round; the score |

To run again:

    python scripts/planted.py --out PLANTED.json && python scripts/planted.py --merge PLANTED.json
    python scripts/build.py --scope-work SCOPE_WORK --pack PACK.npz --work WORK      # the donors' arrays of scope_20261116, and their record
    python -m anm_diag scope WORK/donors.yaml --out OBSERVE
    python scripts/round.py --work WORK --observe OBSERVE && python scripts/round.py --merge
    python scripts/score.py

## What is next

- A disease on real cells, and more than eight units, need data that are either final or the owner's to open.
- A unit event read with a second reader class (boosting) for a short list of proteins, so that "unreadable" can be
  said or withdrawn for them.
- Whether the type table should split donors by batch cannot be asked of this study: batch and donor do not vary
  apart in it.
