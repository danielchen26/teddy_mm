# The repair, designed again: rows that are applied, a typed state, and a row for an event without a record

Written on 25 November 2026, **after an exploration on the THP-1 cells (Section 2) and before any reading of the
melanoma strata by the rules below**. The owner asked for the two open steps to be designed in depth, for the audit to
be run again, and for the open item to be repaired as far as it can be.

Everything here is exploratory in the sense of the ledger: the cells are development cells that earlier runs have
read. No final unit and no unused group is read. Nothing is trained in parts A to D.

## 1. Where things stood

For IRF1 in THP-1 the loop had one item it could not close. With IRF1's own record the repaired output left 0.68 of
the PD-L1 response; with no record 0.90; a reader fitted with IRF1's cells left 0.40 at the state a head reads and
0.18 at the input. [`typed_state_20261125`](../typed_state_20261125/) kept the state by the types of the TEG and
found a gain at the reader (0.27) and, after a check that was not declared, a step in the round once the columns were
at one scale. Two steps were left open: to repeat that by a declared rule on other events, and a repair in which the
types are learned and the layers are trained.

One statement of that run is corrected here. It said that the typed state reads a step more of IRF1 with no record
(0.90 to 0.73) and called that new. That holds for a readout fitted on control cells alone. The tool's reader that is
fitted with the cells of the four pathway genes and never sees a cell of IRF1 leaves 0.64 at the mean state and 0.64
at the typed state. The best reading without IRF1's record was not improved by the typing.

## 2. What the exploration found (THP-1, development cells; not declared; the reason for this design)

The arrays of the run before were read in new ways, in the cloud, before this design was written. The numbers below
are from that exploration and from a trial of this run's script on the two THP-1 targets. They confirm nothing.

**a. The repair did not apply its rows.** `cycle.LinearRows` updates the readout of the control cells by a gain that
weighs the row against the uncertainty of the control fit. For IRF1 that gain closed 27% of the row on the very cells
it was fitted on (gate 3 never passed), because the control fit is far more certain than it has a right to be about
cells it never saw. With the uncertainty of the control fit widened until the row is no longer surprising (Section
3.2), the output with IRF1's record leaves 0.17 [−0.01, 0.29] at the mean state where it left 0.68, on cells the fit
did not see; gate 3 passes and the grade is "relative". A state with its rows permuted stays at 1.00. The four
pathway genes are where they were (0.09). So the item that "one record cannot close" was largely the weight of the
record, not the state.

**b. The same readings say that the reference reader under-reads.** At the mean state the new rows (0.17) leave less
than the tool's reference reader fitted with IRF1's cells (0.40) and less than the round's arm fitted on those cells
(0.42). A ridge reader chosen for its error over all cells does not follow a shift of the mean along a direction in
which control cells hardly vary. "The encoder loses a step for IRF1 (0.40 against 0.18)" is therefore a statement
about that class of readers. This run does not change the tool's reader classes; it reports the difference.

**c. An event-level rule derives IRF1's row from its RNA.** In the input view at one scale, IRF1's shift projected on
the mean shift of the four pathway genes is 0.33 of theirs, and 0.33 times their mean change is −0.26 where the
measured change is −0.28 (92%). Tested on the records with one left out, the same rule gives 85%, 75%, 71% and 117%
of the measured change of the four pathway genes. At the typed state the rule gives 53% for IRF1, at the mean state
19% to 27%. A repair with the derived row leaves 0.22 [0.06, 0.34] for IRF1 at the mean state, with no protein value
of IRF1 anywhere.

**d. Under the new repair the typed state adds little.** With IRF1's record: 0.17 at the mean, 0.10 at the typed
state, 0.09 with the input view beside the mean; the difference is not a step.

## 3. The repaired system

Nothing in the released model changes. What changes is what the head reads and how a record is applied.

### 3.1 The stages a head may read

| stage | what |
|---|---|
| `Z` | the mean of the last layer's output over the gene positions, at unit length: what a head reads today |
| `ZV` | `Z` beside the 512 components of the input view |
| `T` | `Z` kept in two blocks by the TEG: the 128 genes that the recorded events of the type move most, and the rest. **The genes are chosen on the fitted half only** (tokens and labels; no protein), once for each half |
| `R` | the control of `T`: 128 genes at random |
| `P` | types from the input alone (part D): four gene programs found among the control cells of the fitted half, with no record; four blocks and the rest |

Every stage is brought to one scale inside the repair, by the control cells of the fit.

### 3.2 The repair `Full` ([`scripts/repair.py`](scripts/repair.py))

The readout of the control cells is a ridge fit, as before. The rows of the events are then applied with the
uncertainty of that fit widened by one factor *t²* ≥ 1:

    t² = max(1, (r'r − trace(R)) / trace(P)),     readout = controls + (cv/n) t² M D' (t² P + R)⁻¹ r

*r* are the residuals of the rows under the control readout, *R* their errors, *P* = (cv/n) D M D' what the control
fit allows. Where the control fit already explains the rows, *t²* = 1 and this is the old update. Where a row is far
from it, the share of the row that is closed is about 1 − (error of the row)² / (residual)². Arms: `controls` (no
record), `rows` (this rule), `exact` (the rows closed on the fitted cells), `kalman` (the old weight on the scaled
stage), `cells` (fitted on every cell).

The round, its halves, its resamples, *u*, gate 3 and the grade are those of `anm_diag.cycle`, unchanged. The grade
is computed for the new `rows`.

### 3.3 An event without a protein record

| arm | what | what it uses of the event |
|---|---|---|
| `none` | the readout of control cells | nothing |
| `type` | the readout repaired by the rows of the other recorded events of its type | nothing |
| `derived` | the readout of control cells with a derived row for the event | the RNA of its cells in the fitted half |

**The rule that derives a row** (`repair.derive`). In the input view at one scale, on the fitted half: the event's
mean shift from the control cells, projected on the mean shift of the recorded events of the type; that share times
their mean measured change. The derived row is given the error a measured row of as many cells would have. No protein
value of the event enters. The rule is also applied to every recorded event with that event left out, and what it
gives there is reported beside every derived row: a rule is tested on the records before it is used.

### 3.4 Types from the input alone (part D)

The typing of `T` needs the records of a type. `P` needs none: the gene-presence matrix of the control cells of the
fitted half (genes among the tokens of at least 5% of cells) is factorised into four non-negative programs; the 128
genes with the largest loading of each program form a block. The question is whether a typing that the input gives by
itself does what the typing by the records does.

## 4. Step 1: the declared repeat, with nothing trained

**Part A** (the stages `Z` and `ZV` exist; nothing is computed by a model). Five targets of the earlier rounds:

| | target | events | status |
|---|---|---|---|
| T1 | PD-L1, THP-1 | IRF1 on its own | explored before this design: reported, not counted |
| C1 | PD-L1, THP-1 | four pathway genes | explored: reported, not counted |
| **T2** | **MHC-I, melanoma, co-culture** | five pathway genes | **not read by these rules** |
| **C2** | **PD-L1, melanoma, co-culture** | five pathway genes | **not read by these rules** |
| **C3** | **PD-L1, melanoma, IFN-γ** | five pathway genes | **not read by these rules** |

For each: the round as it was (it must repeat the record of 23 November to 0.005), the round with `Full`, the gain of
the new rows over the old on paired resamples, gate 3 and the grade. For each of the fifteen melanoma events (five
genes, three targets): the event on its own with its record, by the old and the new rule; and the event as one
without a protein record (`none`, `type`, `derived`), the other four genes being its recorded type.

**Part B** (states computed on the owner's machine). The same readings with the head reading `T` and `R`, for THP-1
and for the melanoma strata, with the genes chosen on the fitted half.

## 5. Step 2: what is learned, and what is trained

The second step was stated as "types learned from the input, and the layers trained". After Section 2 it is designed
as follows.

- **Types learned from the input: part D.** `P` against `T`, `R` and `Z`, in THP-1 and in the melanoma strata.
- **The layers: part E is a second audit of trainings already made, not a new training.** The five last layers
  trained on 24 November ([`last_layer_trial_20261124`](../last_layer_trial_20261124/)) are read again with `Full`:
  does a layer trained by the rows give a state from which the new repair reads more? A new training is not run in
  this design, for a reason the exploration gives: with the rows applied, the mean state of the released encoder
  already closes IRF1 with its record to a grade, so there is no loss left that a training of the last layer was
  meant to remove; and what is missing for an event without a record is a rule across events, which five recorded
  events cannot teach a layer. If part E shows a trained layer to add a step, that is a reason to train again; if
  not, the next training worth its cost is one over many recorded events, which the development data do not hold.

## 6. Hypotheses (melanoma; THP-1 is reported beside them)

A step is 0.10. "A step" means the interval of a paired gain lies above 0.10; "no step", below it; "harmed" means the
lower end of the interval of a gain is below −0.10.

| | confirmed | refuted |
|---|---|---|
| **H1 the rows are applied** | for T2 with the head reading `Z`, the new rows gain a step over the old | no step |
| **H2 no harm** | C2 and C3 are not harmed by the new rows at `Z`, and at `Z` no melanoma target has a lower grade than it had | one is harmed or has a lower grade |
| **H3 one event, its record** | of the five MHC-I events read on their own, at least three gain a step by the new rule and none is harmed; of the ten PD-L1 events at most one is harmed | fewer than two of five gain a step, or two or more of fifteen are harmed |
| **H4 it is not a trick** | on a state with its rows permuted the new rows leave at least 0.90 for every target; for an event made of control cells the mean change of the output among its cells is below 0.10 of the output's spread | either fails |
| **H5 the guard** | for every target and stage the correlation of the output with the protein among control cells the fit did not see falls by less than 0.03 from `controls` to `rows` | it falls by more for any |
| **H6 without a protein record** | of the five MHC-I events, `derived` gains a step over `none` for at least four, and of all fifteen at most two are harmed | fewer than two of five gain a step, or more than four of fifteen are harmed |
| **H7 the rule on the records** | over the fifteen events the derived change is between 0.67 and 1.5 of the measured one for at least twelve | for fewer than nine |

Unresolved otherwise. Parts B, D and E get their hypotheses and forecasts in an addition to this file that is
committed before their states are read.

## 7. What could make a result wrong, and what the design does about it

| threat | what is done |
|---|---|
| "Closing with its own record" is true of any state that tells the event's cells from control cells. | It is the definition of a repair by a record, and it is said in those words. H4 shows a state without information does not close. The readings that ask more are those without a protein record (H6, H7) and the guard (H5). |
| The new rule was chosen after seeing IRF1 on these cells. | THP-1 is not counted. The rule is fixed in `repair.py` before the melanoma strata are read by it. |
| Applying a row fully follows the noise of a row. | The weight keeps the row's own error (*R*). An event made of control cells is read (H4). `exact`, which drops *R*, is reported beside `rows`. |
| *u* asks for the mean change and a slope among control cells, not for each cell. | The correlation among control cells the fit did not see is reported for every output (H5). |
| The stage of the rule that derives a row (the input view) was chosen on five THP-1 events. | Fixed here before melanoma is read. The rule's result on every recorded event with that event left out is reported (H7). The one event that is not like its type, IRF1, cannot be repeated in the development data: this is said as a limit. |
| The genes of a typed state were chosen with the tokens of the cells that were then read. | In this run they are chosen on the fitted half only. |
| A derived row uses the event's cells. | It uses their RNA in the fitted half and no protein. The arms are named for what they use (Section 3.3). |
| The made-up event and the permuted state use random numbers. | One generator, seed 20261125, in the script. |
| Gate 3 reads the output on the cells it was fitted on. With rows that are closed, that output overshoots by the slope of the control readout on its own cells (in the THP-1 trial, −0.20 at `ZV`), and gate 3 of the protocol can fail although every row is closed. | Gate 3 and the grade are left as the protocol states them and reported as they come. The share of each row that is closed on the fitted cells is reported beside them. A change of gate 3 is the owner's decision. |

## 8. What this run cannot say

- Whether the rule that derives a row holds for an event that is not like its type, beyond IRF1. That needs events
  that no run has read: a batch of knockouts that act through transcription, which is the owner's decision.
- Anything confirmatory. The cells are development cells.
- That the released encoder holds everything the input holds. Section 2b says the earlier localisation depends on
  the class of readers; which reader class the tool should use is a rule of the tool and the owner's decision.

## 9. The forecasts

[`results/declared.json`](results/declared.json): for part A, written before the melanoma strata were read by these
rules.
