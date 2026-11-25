# The account of a change: what came out

25 November 2026. The owner asked why nothing that is specific to the typed event graph (the TEG) had shown any
weight, asked for what is missing from first principles in each of the six links of the loop to be put back as a
repair of the foundation, and asked for the framework to be looked at again. The design is in
[`DESIGN.md`](DESIGN.md). Its Sections 2 and 3 are an exploration on development cells that earlier runs have read;
**explored, not counted**. One reading was declared before it was made: three melanoma targets by event
([`results/declared.json`](results/declared.json)).

Everything here is on development cells. No final unit and no unread unit was opened. **Nothing was trained.**

A step is 0.10. "Stage" means any array a head could read: a view of the model's input, a step inside the model, the
state its head reads.

## What came out

**Two causes, one on top of the other. The measure in use since the run before this one could not see what a stage
holds: it reads near one whenever an event's shift of the stage can be found again. And the events that need the
graph were never targets, while the graph has no edges. With a measure that can prove what it says, the typed state
does carry more than the mean state (about a fifth of IRF1's change, 10 to 16 points of each MHC-I event's), the loss
between the input and the mean state that run 40 had set aside is there (37 points of IRF1's change), and an edge
that was known beforehand forecasts the events that act on the measured protein itself. The loop is rebuilt on three
readings that are each named by their evidence, a chain of four parts cell by cell, and a law of rows with edges
across events. It is a module beside the loop; the protocol is not changed.**

### 1. The first cause: the measure

Run 40 read a repair by letting a readout follow the mean change of an event on one half of its cells and reading the
other half. On the planted system of this run, the event `m` moves the state along a direction the readout does not
follow and changes the readout directly, so no stage holds its change
([`results/mean_fit_demo.json`](results/mean_fit_demo.json), five systems):

| event | planted | u with no record | u with its own record, run 40's rule | left among cells whose state is alike |
|---|---|---|---|---|
| `c1` | through the state | −0.01 | −0.01 | −0.02 |
| **`m`** | **moves the state; the change is direct** | 1.01 | **0.43** | **0.94** |
| `n1` | does not move the state; direct | 1.00 | 1.00 | 1.00 |

The rule "closes" more than half of a change the state does not hold. A row that is closed exactly reads 98% of it.

**The three readings that replace it** (one event, one readout, one stage; fitted on one half of the cells, read on the
other):

| reading | the evidence | what it establishes |
|---|---|---|
| `read` | control cells | what a head reads with no record of the event |
| `carried` | the cells of both arms, inside each arm: the readout given the stage against the arm given the stage | among cells whose stage is alike, how much less the two arms differ; it never fits the difference between the two means |
| `others` | the rows of other events | what their records say of this one: transfer |
| *`found again`* | this event's own mean | that its shift of the stage is found again. **Not a measure of what the stage holds.** |

`overlap` goes with every reading: if a stage tells the event's cells from control cells completely, it holds the
event's identity and `carried` cannot be read from cells; the reading says so.

### 2. The second cause: the events and the edges

The melanoma library under IFN-γ: 206 knockout events, 4,000 control cells, three protein readouts
([`results/explored_library.json`](results/explored_library.json); explored).

| kind | events | u with no record | u with its own record, run 40's rule |
|---|---|---|---|
| through the state, and control cells teach it | IFNGR1, IFNGR2, JAK1, JAK2, STAT1 on PD-L1 | 0.02–0.19 | 0.01–0.08 |
| **not through the state** | CD274 on PD-L1; B2M, TAPBP, HLA-A, HLA-B, HLA-C on MHC-I | 0.92–1.06 | **0.90–1.02** |

- For the first kind a readout of control cells reaches 84% to 100% of the change. Every target of forty runs was of
  this kind, or was IRF1.
- For the second kind the change is after the RNA. CD274's knockout lowers PD-L1 by 1.49; its own transcript falls
  to 54%; the mean state moves by 1.2 times its noise; a readout of control cells reaches 8% from the state and 9%
  with the count of CD274 beside it. With its own record the row is closed on the fitted cells and 0.97 is left.
- 25 of the 206 events move the mean state beyond twice its noise. Of the 173 events whose own transcript can be
  counted, 122 lower it (below 0.8 of control cells and beyond three of its errors; median: to 57%), and the mean
  state moves for 24 of the 173. One site in two thousand positions does not survive a mean.

**Rule S3 treats the relation that matters as a risk.** The declarations name the gene behind each protein
(`events.FORWARD_READOUT`); the rule flags an event at that gene as possible leakage, and such events were kept out
of the targets.

### 3. The law of the rows, with edges (the library; explored)

`Law` fits one readout per protein from the control cells and the rows of all events, with a direct part for what
the state does not explain: small with heavy tails for a pair with no edge, a mean and a spread of its own for the
class "the event's site is a gene of the measured protein" (seven pairs of 618). Every pair is forecast with its
event left out ([`results/library.json`](results/library.json)).

| | with the edges | without |
|---|---|---|
| measured changes inside the 95% interval, all 618 pairs | 97% | 96% |
| the same, the 17 pairs with a change | 71% | 47% |
| error of the forecast over those 17 (no change 1.18; a readout of control cells 0.97) | 0.74 | 0.96 |
| pairs flagged as lacking an edge | TAPBP, SSR2, STAT3, CD44 on MHC-I | those four, and CD274, B2M, HLA-A, HLA-B, HLA-C |

- With its edge, six of the seven pairs are inside their interval (B2M on MHC-I is three errors beyond: −2.87 where
  the class says −0.79 ± 0.69). CD274 on PD-L1: −1.99 measured, −1.00 ± 1.08 forecast.
- **Without any edge the law finds five of the seven pairs by itself**, and **TAPBP on MHC-I** (16.6 errors). Tapasin
  loads MHC-I; the graph did not have that edge. This is the graph growing from the data.
- The rows of other events add nothing through the state in this library (the widening of the control prior is 1.0,
  3.2 and 1.0). What the state holds, control cells already teach.

### 4. THP-1: the account of IRF1 (explored)

[`results/thp1.json`](results/thp1.json). Shares of the measured change of PD-L1 (−0.28).

| the head reads | `carried` | `found again` | `read` | `others` (the four pathway events) |
|---|---|---|---|---|
| `Z`, the mean state | **52%** [37, 66] | 90% [65, 120] | 10% [−12, 30] | 42% [20, 61] |
| `T`, typed by the TEG | **73%** [59, 90] | 95% | 26% [1, 46] | 49% [25, 70] |
| `R`, genes at random | 57% [42, 70] | 83% | 7% | 35% |
| `P`, programs of control cells | 50% [39, 64] | 77% | 2% | 31% |
| `V`, the input view | **89%** [71, 112] | 97% | 3% | 64% [34, 98] |

- **The chain** for IRF1: not in the input 11% [−12, 29]; **lost before the state 37% [17, 61]**; **held, not read 42%
  [22, 65]**; read 10%. Two parts are open: the model before its state, and the head.
- **The typed state** carries 21 points more than the mean state [8, 39] and 17 more than genes at random [7, 28]:
  above zero, not resolved as more than a step. `found again` differs by 5 [−32, +39] between them.
- **The state beside the input view** (`ZV`, what a head that reads both would see) carries 79% [63, 99]: 27 points
  more than the state alone [11, 50], a step.
- **The four pathway events** (the positive control): the mean state carries 85% to 106%, the input view no more;
  a readout of control cells reads 52% to 70%, so "held, not read" is open for each (32% to 37%); the other three
  forecast each at 95% to 107% (STAT1: 154%). The system's own output as trained reads 2% to 7%.
- **The last layers trained on 24 November**: what the state carries of IRF1 changes by 1 to 4 points. One arm (a
  derived record of IRF1) raises `read` by 28 points [12, 48] and `others` by 27 [13, 44]. Training moved what is
  easy to read; it did not add to what a cell's state holds.

### 5. The declared reading: three melanoma targets, by event

[`results/melanoma.json`](results/melanoma.json), [`results/outcome.json`](results/outcome.json). Fifteen accounts.

| | hypothesis | reading |
|---|---|---|
| HA1 | `found again` exceeds `carried` at `Z` by more than a step, for at least four of five MHC-I events | **confirmed**: five of five, by 36 to 44 points |
| HA2 | the typed state carries a step more than the mean state, for at least three of five | **refuted** as declared: none is resolved as a step. All five gains are above zero (10 to 16 points; median 13; intervals from 3–9 to 17–24) |
| HA3 | the typed state carries a step more than genes at random, for at least three of five | **confirmed**: three of five (13 to 19 points for all five) |
| HA4 | "lost before the state" is open for at least three of five MHC-I events | **refuted**: for none. The input view carries what the state carries (48% to 63% against 51% to 60%) |
| HA5 | at least seven of ten PD-L1 accounts have no open part | **unresolved**: four of ten. In co-culture every account has an open first part (28% to 38%), under IFN-γ one of five. (With the chain as amended: 17% to 27% not in the input in co-culture, open for four of five; none under IFN-γ) |
| HA6 | a permuted stage carries nothing; an event made of control cells has no change to account for | **confirmed**: 15 of 15 and 15 of 15 |
| HA7 | `overlap` at `Z` is at least 0.10 | **confirmed**: 0.49 to 0.68 |
| HA8 | `others` at `Z` is at least 0.75, for at least ten of fifteen | **confirmed**: fifteen of fifteen (89% to 125%) |

What the reading showed that the design had not:

- **A protein is not in a cell's RNA in full, and its mean change still follows a law.** For MHC-I in co-culture the
  input view and the state carry about half of each knockout's change cell by cell (37% to 46% is not in the input,
  open for all five), and the other four events forecast the fifth at 91% to 110%. The chain is about single cells
  (what a head could read from a cell). What the event graph must supply is read across events. The first part of
  the chain was named "direct" with the place "an edge"; it is now **not in the input**, and it needs an edge only if
  the other events do not forecast the change.
- **The state can carry more than a view of the input.** For PD-L1 under IFN-γ the state carries 85% to 102% and the
  view `V` 80% to 96%. The chain now takes the input to hold at least what the state holds.
- **The head.** For MHC-I in co-culture a readout of control cells reads 36% to 40% where the state carries 51% to
  60%; "held, not read" is 15 to 23 points (open for three of five). The system's own output as trained reads 3% to
  4% of these changes and 0% to 10% of the PD-L1 ones.

### 6. The planted suite

[`results/planted.json`](results/planted.json). The planted system has four hidden factors (the readout follows three;
the state loses one that the input holds; one is off in control cells) and events whose part is known: read; held
and not read; lost before the state; not in the input, with and without an edge; one that moves the state while its
change is direct; and thirty that fill the library. Truth is the same system with very many cells and exact least
squares. 73 checks per system: every part against the planted one, the place the account points at, the two
controls, and for the law its calibration, its edges and its candidates.

| | systems | checks | wrong | the planted place resolved as a step |
|---|---|---|---|---|
| as planted | 20 | 1,460 | 3 | 173 of 180 |
| every cell of an event edited | 5 | 365 | 0 | 45 of 45 |
| half of the cells edited | 5 | 365 | 1 | 44 of 45 |
| noisy cells (five times the measurement noise) | 5 | 365 | 2 | 44 of 45 |
| few cells (1,500 control cells, 300 for an event) | 5 | 365 | 5 | 40 of 45 |

Eleven of 2,920 checks are wrong. Nine are a planted part outside its interval by a little; two are an account that
points at a neighbouring part (a small change, with few or noisy cells). No event that fills the library was
flagged. With every cell edited the stage tells the arms apart completely; the reading says "by identity" and does
not call the change direct. With noisy cells the planted truth is itself lower: a noisy stage carries less of a
change cell by cell (for `r1` about half), and the reading follows it.

## The forecasts

**Four of seven intervals held the value; the outcome given the higher probability happened for six of eight**
([`results/outcome.json`](results/outcome.json)). I expected the input view to carry more of the MHC-I events than
the state (a median of 0.05 to 0.30; it was 0.01) and more PD-L1 accounts to have nothing open (five to ten; four).
The forecast of the planted suite is not counted: its first full run was seen before the declaration was committed.

## What was not as written

1. **Three amendments after the declared reading** (DESIGN.md, Section 12): the name and the place of the first part
   of the chain; the input holding at least what the state holds; and the numbers of Section 2.1, which came from a
   first planted system of the exploration. The hypotheses are scored as declared, on
   [`results/melanoma_as_first_read.json`](results/melanoma_as_first_read.json). `results/melanoma.json` is the same
   reading put together by the module as amended; `scripts/score.py` checks that the two hold the same readings.
2. **The planted suite was seen before the commit** and two things were changed after its first full run: an account
   lists the parts that are unresolved beside those that are open; and the check of the place asks that the account
   point at the planted place, with "resolved as a step" counted apart.
3. **Numbers given to the owner during the work came from a prototype** (within-arm slopes, intervals from the cells
   that are read only): the typed state over the mean state +23 points [18, 29], "more than two steps"; and "closes
   75%" for the planted event. The module's intervals also draw the control readout from what its fit leaves
   uncertain and are wider: +21 [8, 39], not resolved as a step; and 57% under run 40's rule on the planted system
   of the module.
4. **The measure was chosen on THP-1.** Everything about IRF1 is explored twice over. The melanoma targets are the
   declared test; they hold no event unlike its type and no event that acts on the readout itself.
5. **One class of edge, seven pairs, one library.** The forecast by an edge is explored. Its test is an event that
   acts on a measured protein in cells no run has read.
6. The hyper-parameters of the law are fitted with every event; the forecasts leave one out. The planted suite
   scores the same procedure.
7. The typed states and the trained states are those of runs 40 and of 24 November, made on the owner's machine;
   the counts of the knocked-out genes were read there from the screen's RNA file (labels and counts only).

## What it changes in earlier write-ups

- **Run 40** ([`../typed_repair_20261125/`](../typed_repair_20261125/)): a note at its head and four marked
  passages. What stands: the old rule did not apply a record; the rule that derives a row; the controls. What is
  withdrawn: that IRF1 "closes" at the plain mean state, that the loss at the mean is smaller than reported, that the
  typed state adds 0.05 and never a step, and that the trained layers add nothing in every sense.
- **Run 39** ([`../typed_state_20261125/`](../typed_state_20261125/)): its reading of the typed state (a reader with
  IRF1's cells, 0.40 to 0.27; at one scale 0.41 to 0.24) agrees with this run's (`carried` 52% to 73%).
- **Rule S3** is not changed. Its flag is right as a warning for an audit of the model alone; for the system with
  its graph the same relation is an edge.

## What it means

- **For the TEG.** Its weight is in two places, and neither had been measured. Inside the model: a state typed by the
  recorded events of a type holds more of an event of that type, cell by cell, than a mean does (IRF1: 21 points;
  five MHC-I events: 10 to 16 each). Outside the model: an edge between an event and a readout forecasts what no
  stage holds. As built, the graph is a table of records with labels; the law of this run gives it edges, and a way
  to propose the ones it lacks.
- **For the encoder.** For IRF1 the mean state loses 37 points [17, 61] that the input view holds. That is a loss
  inside the model, cell by cell, and it is the place a state built on the TEG is for. For the melanoma targets the
  state loses nothing the input view holds.
- **For the head.** A readout of control cells reads less than the state carries in THP-1 (32 to 42 points) and
  for MHC-I in co-culture (15 to 23); for PD-L1 in melanoma the two are within about a step (−5 to 11). The head as
  trained reads almost nothing of any of these changes (0% to 10%). A head taught with the cells of recorded events
  is the repair, credited up to what the state carries.
- **For training.** Training last layers by rows did not add to what the state carries. It did raise what a readout
  of control cells reads. The question that is open again is a state that keeps what the input holds, and a training
  is one way to it; typing the state is another, and needs none.
- **For an event without a record.** Three statements can be made, each with its error: through the state, from the
  rows of other events (fifteen of fifteen in melanoma); by its edge, as a class (wide); and "no change", with the
  rate at which that was wrong (4 of 611 pairs without an edge were flagged).
- **For what can be claimed.** Nothing here is confirmed. What is new is that each kind of claim now has a reading
  that can prove or refute it, and a planted system on which that reading was scored.

## What is next, and whose decision it is

1. **Whether the account replaces the round's reading in the protocol** (the measure of a repair, gate 3, the
   reference reader). The owner's decision; nothing in `cycle.py`, `outer.py` or `scope.py` was changed.
2. **A batch that no run has read, forecast before it is opened.** CD274 in THP-1 is the test of an edge across
   kinds of cell; THP-1 knockouts that act through transcription and are unlike the pathway are the test of the
   chain for IRF1. Writing the forecasts needs the RNA of unread units. The owner's decision.
3. **A head taught with the cells of recorded events**, as the round's repair, with the ceiling. Development cells;
   no training of the model.
4. **A state that keeps what the input holds about an event of a type**: typed by the TEG at more than one event and
   more than one kind of cell, judged by `carried` against genes at random. If a training is wanted, it is the
   owner's decision.
5. Of the nine rule changes that were waiting, four are answered by this run (DESIGN.md, Section 10); five remain.

## Files

- [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json): written first. Section 12 of the
  design was added after the readings.
- `anm_diag/account.py` (`reading`, `gain`, `chain`, `must_fail`, `Law`), `anm_diag/toy_account.py` (the planted
  system and its suite), `tests/test_account.py`.
- [`results/melanoma_as_first_read.json`](results/melanoma_as_first_read.json),
  [`results/melanoma.json`](results/melanoma.json), [`results/outcome.json`](results/outcome.json): the declared
  reading and its hypotheses. [`results/thp1.json`](results/thp1.json),
  [`results/library.json`](results/library.json), [`results/explored_library.json`](results/explored_library.json):
  explored. [`results/planted.json`](results/planted.json), [`results/mean_fit_demo.json`](results/mean_fit_demo.json):
  planted systems.
- `scripts/`: `melanoma.py`, `thp1.py`, `library.py`, `explore_library.py`, `mean_fit_demo.py`, `score.py`, `_read.py`;
  `own_gene.py` and `site_counts.py` read the counts of the knocked-out genes from the screen's RNA file.
