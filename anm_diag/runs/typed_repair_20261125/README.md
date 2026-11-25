# The repair, designed again: what came out

25 November 2026. The owner asked for the two open steps to be designed in depth, for the audit to be run again, and
for the open item to be repaired as far as it can be. The design is in [`DESIGN.md`](DESIGN.md). It was written after
an exploration on the THP-1 cells, which it reports and does not count, and **before the melanoma strata were read by
its rules**. Parts B, D and E were declared after part A was read and before their states were read. Thirty
forecasts in all ([`results/declared.json`](results/declared.json)).

Everything here is exploratory in the sense of the ledger: development cells that earlier runs have read. No final
unit and no unused group was read. **Nothing was trained.**

*u* is the share of the measured change that is left unexplained: 0 is none of it, 1 is all of it. A step is 0.10.

## What came out

**The repair of the round had not been applying its rows. With the rows applied, IRF1 closes to a grade with its own
record, at the plain mean state and with nothing trained; and a row derived from its RNA alone brings it almost as
far. On the melanoma strata, where the rule was tested as declared, it helps every event that is read on its own and
does nothing for the one pooled target that was open; it also trips two of the loop's own checks. The typed state
adds about 0.05 wherever the state is the limit and never a step. Trained last layers and types found from the input
alone add nothing.**

### 1. IRF1 in THP-1 (explored before the design: reported, not counted)

| the head reads | no record | the repair of its type | **a row derived from its RNA** | **its own record** | its own record, the repair as it was |
|---|---|---|---|---|---|
| `Z`, the mean | 0.90 | 0.63 | **0.22 [0.06, 0.34]** | **0.17 [−0.01, 0.29]** | 0.68 [0.57, 0.79] |
| `T`, typed by the TEG | 0.74 | 0.54 | 0.17 [0.01, 0.30] | 0.11 [−0.06, 0.24] | |
| `ZV`, the mean beside the input view | 0.96 | 0.62 | 0.15 [−0.01, 0.27] | 0.09 [−0.09, 0.21] | 0.39 |

- **With its own record the item reaches a grade.** At `Z` the new rows remove 0.51 [0.41, 0.65] that the old rows
  left. Gate 3 passes for the first time (the row is closed to 98% on the cells it was fitted on, where the old
  weight closed 27%), and the grade is "relative": the output is no longer rejected and adds nothing to the reference
  reader of the stage.
- **Without a protein value of IRF1.** The rule derives −0.26 where −0.28 was measured (92%), from the RNA of IRF1's
  cells and the records of the four pathway genes. A repair with that row leaves 0.22; against the repair with the
  measured row the difference is 0.06 [0.04, 0.07], no step. The repair of its type, with no cell of IRF1, leaves
  0.63.
- **It is not absolute.** No output is inside the margin with its whole interval.

### 2. The declared test: three melanoma targets and their fifteen events (part A)

| | what was asked | reading |
|---|---|---|
| H1 | for MHC-I in co-culture, pooled, the new rows gain a step over the old | **refuted**: 0.36 and 0.36 |
| H2 | the two PD-L1 targets are not harmed, and none has a lower grade | **refuted**: neither is harmed (0.11 to 0.09, 0.04 to 0.01), but PD-L1 under IFN-γ fails gate 3 with the new rows and has no grade where it was "relative" |
| H3 | events read on their own with their record | **confirmed**: all five MHC-I events gain a step (0.46–0.49 to 0.28–0.36); none of fifteen is harmed |
| H4 | it is not a trick | **confirmed**: on a state with its rows permuted every target stays at 1.00; an event made of control cells moves the output by 0.001 of its spread |
| H5 | the guard: the correlation among control cells falls by less than 0.03 | **refuted**: for MHC-I at `Z` it falls by 0.036 (0.537 to 0.501) |
| H6 | without a protein record, the derived row | **confirmed**: for all five MHC-I events it gains a step over no record (0.61–0.65 to 0.29–0.41); none of fifteen is harmed |
| H7 | the rule on the records, one left out | **confirmed**: for fifteen of fifteen the derived change is between 0.75 and 1.07 of the measured one (median 0.80) |

- **Why H1 fails, and what it says.** For the pooled MHC-I target the old rule had already closed 67% to 87% of
  every row. Applying the rows fully changes nothing in the output, because the price is paid elsewhere: the slope of
  the readout among control cells falls from 4.5 to 3.0. In that state the knockout cells are not far enough from the
  control cells for one linear readout to move them and keep the control cells where they are. There the state is the
  limit; the input view itself leaves 0.23. For IRF1 in THP-1 the same slope falls from 1.51 to 1.42.
- **Why H2 fails.** With rows that are closed, the output on the fitted cells overshoots by the slope of the control
  readout on its own cells (−0.09 here). Gate 3 reads that as a loss. The design named this before the reading
  (Section 7) and left gate 3 as the protocol states it.
- **The derived row is short by a fifth.** For members of a type the rule gives 75% to 107% of the measured change,
  and less than 90% for twelve of fifteen. For PD-L1 in melanoma, where a readout of control cells already leaves
  only 0.03 to 0.34, the derived row gains nothing and the repair of the type does better (−0.16 to 0.18).

### 3. The typed state under the new repair, and types from the input alone (parts B and D)

The genes of a type were chosen on the fitted half only. `T` is two blocks by the TEG, `R` genes at random, `P` four
programs of the control cells found with no record.

| | the rows leave, at `Z` | at `T` | at `R` | at `P` |
|---|---|---|---|---|
| MHC-I, co-culture, pooled | 0.36 | 0.30 | 0.38 | 0.41 |
| IRF1 with its record (THP-1) | 0.17 | 0.11 | 0.21 | 0.26 |
| PD-L1, co-culture | 0.09 | 0.08 | 0.10 | 0.11 |
| PD-L1, IFN-γ | 0.01 | −0.02 | −0.01 | 0.03 |

- **HB1 refuted.** Where the state is the limit, the typed state removes 0.057 [0.04, 0.07]: not a step. The same
  0.05 to 0.06 for each of the five MHC-I events on its own (HB3 refuted), and 0.05 for IRF1.
- **HB2 unresolved.** Against genes at random the typed state removes 0.084 [0.07, 0.10]. The gain is small and it is
  the type: random genes cost, the genes of the events gain.
- **HB4 and HB5 confirmed.** At `T` the guard holds (the correlation falls by 0.023) and the PD-L1 targets are not
  harmed.
- **HD1 refuted.** Four programs of the control cells leave more than the mean at every target. They share 0 to 4
  genes with the genes of the events; in THP-1 one of them is the cell cycle. A typing that the input gives by itself, in this
  simple form, is not the typing the records give.

### 4. The last layers trained on 24 November, read again (part E)

**HE refuted.** With IRF1's record the output leaves 0.17 at the released state and 0.18 to 0.21 at the five trained
ones; with a derived row 0.22 and 0.23 to 0.26. A trained layer still helps a readout that has no record at all
(0.90 to 0.74 for the arm trained with IRF1's record, a step), which is what the trial of 24 November found. Once the
repair applies a record, the trained layer has nothing to add.

## The forecasts

Part A: **4 of 8 intervals held the value; the class given the highest probability happened in 6 of 9**
([`results/outcome_a.json`](results/outcome_a.json)). I expected the pooled MHC-I target to gain (0.05 to 0.28; it
stayed at 0.36) and the guard to hold. Parts B, D and E: **6 of 6 intervals; 6 of 7 outcomes**
([`results/outcome_b.json`](results/outcome_b.json)).

## What was not as written

1. **The rule was found on the cells of IRF1.** Everything about IRF1 and the THP-1 pool here is exploratory twice
   over: development cells, and a rule chosen after seeing them. The melanoma strata are the declared test, and they
   hold no event that is unlike its type. So the two results that matter most for IRF1, the closing with its record
   and the derived row, have no repeat on an event of its kind.
2. **One hypothesis was reworded before the reading.** H2 first asked that no target lose its grade at any stage; a
   trial on THP-1 had shown gate 3 failing at `ZV` for the reason given above, and H2 was restricted to `Z` and its
   forecast lowered, in the commit of the declaration. It was refuted all the same.
3. **The made-up event.** The design first expected the widening factor to stay at one for an event made of control
   cells. Under no change it exceeds one about half the time by chance (3.99 in the trial); what must stay small is
   the change of the output, and H4 was written for that before the reading.
4. **No grid of the tool for the typed states.** A stage that differs by half cannot be given to the tool as one
   array. The grade at `T`, `R` and `P` is taken against the reference reader of `Z`.
5. **Where it was run.** The typed states were made on the owner's machine (one pass of the encoder over 11,171
   melanoma cells; the last layer over the THP-1 cells). Every reading was made in the cloud from those files and from
   the arrays of 23 November.
6. **Part E is a second reading, not a training.** The design gives the reason; the result agrees with it.
7. The gene names of the melanoma sets are written as vocabulary entries; the file of tokens holds no symbols.

## What it means

- **For the item the loop could not close.** IRF1 was not a case of "one record cannot repair it". The repair
  weighed a record against a control fit that was sure of cells it had never seen. With the record applied, the
  released encoder's plain state is enough for the grade the loop calls relative. This is exploratory.
- **For the encoder.** The reading "the state loses a step that the input holds" (0.40 against 0.18) came from ridge
  readers. An output that is made to follow the mean change leaves 0.17 at that same state. The loss at the mean is
  smaller than it was reported, and how small depends on a class of readers the tool does not have.
- **For an event without a record.** Two things now reach it. A record derived from its RNA by a rule that was tested
  on the records (fifteen of fifteen in melanoma, within a fifth to a quarter), and, for members of a type, the repair
  of the type. Neither needs a training. What they need is the event's RNA, and for the derived row a type whose
  shift the event shares in part.
- **For the typed state.** Its place is smaller than the run before made it look: about 0.05 where the state is the
  limit, at IRF1 and at MHC-I alike, never a step, and only when the genes are those of the recorded events.
- **For training.** A last layer trained by the rows adds nothing once the repair applies a record. Nothing found
  here asks for a training of the encoder.
- **What is not repaired.** MHC-I in co-culture stays at 0.30 to 0.36; the input view itself leaves 0.23, so the
  limit there is before the model. And nothing here is confirmed.

## What is next, and whose decision it is

1. **Four rules of the loop, the owner's decision; none was changed.** (a) The round's reader brings a stage to one
   scale. (b) The repair applies its rows, with the uncertainty of the control fit widened until they are not
   surprising. (c) Gate 3 for a repair that closes its rows: stated on the rows, not on an output that overshoots by
   the slope of its own cells. (d) The tool's reader classes: one that follows the mean change.
2. **A batch of events no run has read**, to test what cannot be tested here: knockouts in THP-1 that act through
   transcription and are not like the pathway. For each: the rows applied, and the derived row against the measured
   one. That uses unused groups and is the owner's decision.
3. **No training is proposed.**

## Files

- [`DESIGN.md`](DESIGN.md), [`results/declared.json`](results/declared.json): written first (Section 10 of the
  design and `parts_B_D_E` after part A, before their readings).
- [`results/round_a.json`](results/round_a.json), [`results/round_b.json`](results/round_b.json),
  [`results/round_e.json`](results/round_e.json): the readings. [`results/outcome_a.json`](results/outcome_a.json),
  [`results/outcome_b.json`](results/outcome_b.json): the forecasts. [`results/states/`](results/states/): the genes
  of every set of the typed states, committed before their reading.
- `scripts/`: `repair.py` (the repair `Full`, the rule that derives a row), `round_a.py`, `states.py` (needs the
  model), `round_b.py`, `round_e.py`, `score_a.py`, `score_b.py`.
