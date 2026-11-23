# The head reads Z together with a view of the input: declaration

Written on 23 November 2026, **before any reading of this run was made**. It is step 3 of the plan in
[`GENERAL.md`](../../GENERAL.md) (Section 7): a repair in front of the head at small scale. The question of the plan:
does a loss that the loop placed before the head (IRF1, [`outer_step_20261122`](../outer_step_20261122/)) close
without retraining the encoder, once the head is given a view of the model's own input beside Z?

Everything here is exploratory. The cells are development cells that earlier runs have read. No final unit is read,
no new state is computed, and no model is trained.

## 1. Why this is asked

For IRF1 in THP-1 a linear reader of 512 components of RNA, fitted with IRF1's own cells, leaves 0.18 [0.04, 0.30] of
the PD-L1 response. A reader of Z leaves 0.40 [0.27, 0.50]. The difference, 0.21 [0.12, 0.32], is a loss between the
model's input and the stage its head reads. A repair of the head cannot remove it: the head reads Z. The plan names
two ways on: retrain the encoder (step 5, several GPUs), or let the head read the input as well. This run tries the
second, because it costs no GPU and says whether the first is needed for this item.

MHC-I in co-culture is the second place where a reader of RNA read more than every reader of Z (0.12 [0.09, 0.15]
more). It was not drawn as a loss before the head, because the interval does not clear the step margin and no reading
of that target is supported. It is read here in the same way.

## 2. The joint stage

For every cell of a stratum: `ZV` = the state the head reads (Z at unit length) beside the 512 components of the
stratum's RNA (the array the loop already reads as the input view). Each of the two blocks is divided by the square
root of its total variance among the control cells of the stratum, so that the two carry the same total variance.
The scaling uses control cells only. No readout and no label of an event enters the stage.

A system whose head may read `ZV` is one more system for the same loop (GENERAL.md, Section 2): only the stage array
differs. The code that reads it is `anm_diag.cycle`, unchanged. The repair is the one used so far (`cycle.LinearRows`),
a linear readout of the stage, fitted on control cells and updated once by the records of the events.

## 3. What is read

| | target | events | read as | why it is here |
|---|---|---|---|---|
| T1 | PD-L1, THP-1 | IRF1 | on its own, with its own record; the reference is a reader of the stage with its own cells | the loss placed before the head |
| T2 | MHC-I, melanoma, co-culture | IFNGR1, IFNGR2, JAK1, JAK2, STAT1 | a pool | a reader of RNA read more than every reader of Z |
| C1 | PD-L1, THP-1 | IFNGR1, IFNGR2, JAK2, STAT1 | a pool | control: nothing was lost before the head |
| C2 | PD-L1, melanoma, co-culture | the five pathway genes | a pool | control: RNA and Z read the same (difference −0.00 [−0.04, 0.03]) |
| C3 | PD-L1, melanoma, IFN-γ | the five pathway genes | a pool | control: the reference at Z is supported |

For every target the audit tool reads the stages `Z` and `ZV` in one run (the same cells and the same resamples), with
ridge and boosting, as package 2 declares for Z. For T1 and T2 the run also reads the 512 components alone. Then the
round is run twice with the general code, once with the head reading `Z` and once reading `ZV`: the same cells, the same
halves (seed 0, as the earlier rounds), the same resamples. Every difference between the two is a difference on paired
resamples.

Margin 0.10, step margin 0.10, 1,000 resamples. Claim, protocol (version 3) and reader specification are those of
package 2. Nothing is amended.

## 4. What counts as what

`gain` of an arm = |u| with the head reading Z minus |u| with the head reading ZV, on paired resamples. A gain is
**a step** if its interval lies above 0.10, **no step** if it lies below 0.10, and unresolved otherwise.

- **H1, the plan's question. IRF1 closes without retraining.** Confirmed if, with the head reading ZV, the grade of the
  round is relative or absolute (A3 holds, gate 3 passes, the arm `rows` adds no loss beyond the reference, which is
  now a reader of ZV) **and** the gain of `rows` is a step. Refuted if the gain of `rows` is no step. Otherwise
  unresolved.
- **H2. The same for MHC-I in co-culture**, with the pool.
- **C1 to C3. No harm.** For each control target the arm `rows` at ZV loses no step more than at Z: the interval of
  −gain lies below 0.10. A control that fails means the wider stage costs something where it has nothing to give.
- **What the stage could give at most**: the gain of the reference reader (a reader of Z against a reader of ZV, the
  scope the round uses). Reported for every target.
- **Check**: the readings with the head reading Z repeat those the first outer step recorded (`bridge.json`). The
  largest difference is reported.

What each outcome of H1 means for the plan:

| H1 | reading |
|---|---|
| confirmed | this loss does not need a retraining: giving the head the input view removes it |
| refuted | the head cannot recover it from the input view by this kind of repair; the item stays with the encoder (step 5) |
| unresolved | reported as it is, with the gain and its interval |

One limit is known before the run. In the first outer step the reader of the input view was the best of the readers
tried, which favours the input view. A gain smaller than 0.21 would not be a surprise.

## 5. The order

1. Commit this file and the forecasts (`results/declared.json`).
2. `scripts/build_zv.py`: the arrays with the joint stage and the SPEC files. No readout enters.
3. `anm_diag scope` for the five SPEC files.
4. `scripts/round.py`: the two rounds and their paired differences, once. A result file that exists is not overwritten.
5. `scripts/score.py`: the forecasts against what was read.

## 6. The forecasts

`results/declared.json`. Briefly, and why:

- **IRF1.** With its own record at Z the repaired output left 0.68; fitted on every cell of the fitted half, 0.42. One
  record updates 512 weights from about 540 cells, and it does that poorly. At ZV the same record updates 1,024
  weights. I expect the reader of ZV to read more than the reader of Z (0.15 to 0.33 against 0.40), the arm fitted on
  every cell to gain, and the arm `rows` to gain less than a step. H1: unresolved is the most likely outcome,
  confirmed the least.
- **MHC-I in co-culture.** What RNA read beyond Z was read by boosting, not by a linear reader. A linear readout of
  ZV is expected to gain little. H2: refuted is the most likely outcome.
- **Controls.** No harm is expected for all three.

## 7. What would change the plan

- The readings at Z do not repeat the recorded ones to 0.005: stop and report.
- Anything the scripts cannot do as written is reported as a deviation, with what was done instead.
