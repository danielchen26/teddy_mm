# The round on batch 1: declaration

Written on 23 November 2026, **before the system's states of these cells were computed and before any of their protein
values was read by this run**. It is step 4 of the plan in [`GENERAL.md`](../../GENERAL.md) (Section 7): the states of
batch 1 and one round on them. The batch was listed by the first outer step
([`outer_step_20261122`](../outer_step_20261122/), `results/census.json`, `batch_1`): 8 units and control cells, 6,759
cells, none of them final.

This is the first time the loop is run on units the system has never been read on. Everything is still development:
the 133 final units stay untouched.

## 1. What the batch is, and what it is not

| stratum | units chosen for an open item | calibration unit (drawn at random) | control cells |
|---|---|---|---|
| THP-1, stimulated (`thp1`) | BRD4 (341 cells), CMTM6 (817) | POU2F2 (662) | the 2,386 that have states already |
| melanoma, co-culture (`mel_co`) | CD44 (289), TAPBP (170) | CDH19 (267) | 1,276 new ones |
| melanoma, control condition (`mel_ctl`) | B2M (150) | KLF4 (142) | 2,645 new ones |

The control condition is a stratum the system has not been read in at all.

What the batch can answer: whether the round's repair closes the model's part for an event it has never met, whether
a repair built from the earlier events of a stratum carries to a new event, whether knockouts that act on the measured
complex itself (TAPBP, B2M) read the same way in two more conditions, and how often a unit drawn at random is a
target at all.

What it cannot answer: PD-L1 in melanoma under IFN-γ, the cleanest result so far (0.037 [−0.034, 0.101]), stays at
the margin. That stratum has no unread unit left.

The protein response of every knockout of both screens was read once on all cells (the census of 13 November). The
units were chosen with it. So admission here is made again on a selection third, and those cells enter no reading.

## 2. Roles, sealed before the states are computed

`results/roles.json`. All 8 units are **development**: read, never final. None is used to fit anything that is read
on another unit. The repair that carries the name `earlier` is fitted on events the earlier rounds read (the five
pathway genes of each stratum); those events are the **fit** units of this round.

The selection third of this run: a cell is in it if SHA-256 of `20261123|<cell name>` is divisible by 3. Units and
control cells alike.

## 3. The order

1. Commit this file, the roles and the forecasts.
2. Compute the states of the 6,759 cells (TEDDY-G 70M, the official tokeniser, as the states of the earlier runs). A
   check on cells that have states already has to reproduce them to 1e-4 or the run stops.
3. Admit: A1 (change at least three half-widths) and A2 (not the gene that encodes the readout) on the selection
   third, for every readout of the claim in the stratum (PD-L1; in melanoma also MHC-I and HLA-E).
4. Audit the admitted events with the tool, under the reader specification of package 2: Z by ridge and boosting; the
   input view by ridge and boosting on 512 components of RNA and by ridge on 4,000 genes; a permuted view as control.
5. Read each admitted event on its own (`anm_diag.cycle`, the general code; no system-specific reading code). A new
   event is a type of its own until a reading says otherwise, as IRF1 was: its reference is a reader of Z fitted with
   its own cells, and resampling is by cell.
6. Place what is lost before the model, derive the targets and the open items (`anm_diag.outer`).
7. Compare with the forecasts, one by one.

One reading per event. A script that finds its result file refuses to run again.

## 4. The arms

Every cell is read by a fit that did not see it (two halves).

| arm | what it is | what it answers |
|---|---|---|
| as trained | the system's own output | what the model does today |
| `controls` | a linear readout of Z fitted on control cells of the stratum | what a readout of this kind of cell reads without any record of the event |
| `rows` | the same, updated once by the event's own record (from the fitted half) | the round's repair |
| `earlier` | the same, updated by the records of the stratum's earlier events only | does the earlier repair carry to an event it never met |
| `rows_other`, `other` | the records, and the whole repair, of another stratum | does anything carry across conditions or cell types |
| `cells` | fitted on every cell of the fitted half | the ceiling of this kind of repair |

`earlier`: THP-1: IFNGR1, IFNGR2, JAK2, STAT1 (type A of the first outer step). Co-culture: IFNGR1, IFNGR2, JAK1,
JAK2, STAT1. The control condition has no earlier event. Other stratum, for all three: melanoma under IFN-γ (for
PD-L1 its five pathway genes; for MHC-I, B2M and TAPBP outside that step's selection third).

## 5. What counts as what

- **The model's part is closed for an event**: A3 holds (the output as trained loses a step more than the reference),
  gate 3 passes, and the `rows` arm adds no loss beyond the reference (grade relative), or is supported (absolute).
- **Nothing for the model to close**: A3 does not hold. Where no reader of Z or of the input view reads the event,
  the target is open before the model; where the readings are only unresolved, undecided before the model.
- **The earlier repair carries**: the `earlier` arm adds no loss beyond the reference. Where it adds a loss, the new
  event is not of the earlier type.
- **The rows matter**: `rows` minus `controls` is a step.

Margin 0.10, step margin 0.10, 1,000 resamples, as before. The claim, the protocol (version 3) and the reader
specification are those of package 2. Nothing is amended here.

## 6. The forecasts

`results/declared.json`, written with this file. Probabilities are for outcomes, intervals for numbers; each is scored
after the reading (inside or outside; the class given the highest probability happened or not). The reasons, briefly:

- TAPBP and B2M act on the measured complex. Under IFN-γ no reader of Z or of RNA read them (0.91 to 1.02). The same
  is expected in co-culture and in the control condition: nothing for the model to close.
- CMTM6 keeps PD-L1 at the membrane after it is made. The message of CD274 need not change. Expected to behave like a
  knockout that acts on the protein, with more doubt than for TAPBP and B2M.
- BRD4 acts on transcription, and its knockout moves PD-L1 the other way (up). Expected to be readable at Z with its
  own cells, lost by the head, and closed by its own record but not by the earlier repair.
- A unit drawn at random is expected not to pass A1 for any readout.

## 7. What would change the plan

- The states of cells that have states already are not reproduced: stop, report.
- No event is admitted for a stratum and readout: it is reported as not a target, and nothing else is read there.
- Anything the scripts cannot do as written here is reported as a deviation, with what was done instead.
