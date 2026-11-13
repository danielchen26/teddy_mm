# The loop: protocol

Written 12 November 2026, after three exploratory rounds on the teddy_mm head
([`runs/loop_20261112`](runs/loop_20261112/), [`runs/loop_guard_20261112`](runs/loop_guard_20261112/),
[`runs/repair_20261112`](runs/repair_20261112/)) and a re-reading of their records
([`runs/loop_reading_20261112`](runs/loop_reading_20261112/)). Those rounds decided pass or fail by comparing a point
value with a threshold. Their audits can show that a system loses a response. They can rarely show that a changed system
no longer loses it. This protocol separates the two questions and puts a cheap check in front of every expensive step.

It is general: a system with a chain of stages, a source, a readout, and a change to the system. The teddy_mm head is the
first instance. The readings are made by [`anm_diag/loop.py`](anm_diag/loop.py) from finished runs; it makes no fits.

## Two loops

An audit locates a step: the response is readable before it and not after it. Two different things can be done next, and
they need different evidence.

| | diagnostic loop | repair loop |
|---|---|---|
| question | Is the located step the one that matters for this response? | Is the changed system free of the loss, and still good at its own task? |
| what is compared | two changes against one baseline: one made at the located step, one made elsewhere | the changed system against a control that differs only by the change |
| quantity | the share of the baseline's loss over a span that each change removes, and the difference, each with an interval | the loss the changed system still adds over the span, and the verdict of the span's end stage |
| kind of claim | one change removes more than the other | the loss is absent |
| what it needs | a baseline loss that is large against its interval | intervals narrow against the margin; data no design decision has seen |
| command | `anm_diag loop diagnose` | `anm_diag loop closure` |

The diagnostic loop tests the audit's statement. The repair loop tests the changed system. Passing the first says where to
work. It does not say the system is repaired.

## Quantities

- *u* at a stage: the unexplained share of the response (as in every run). Stage verdicts in the loop are two-sided:
  support needs the whole interval of some declared readout class inside (−margin, margin).
- **Loss added over a span** a → b: |*u*<sub>b</sub>| − |*u*<sub>a</sub>|, computed on every bootstrap draw. Read three ways
  against the step margin: **loss** (lower bound above it), **no added loss** (upper bound below it), **unresolved**.
- **Removed by a change**: the baseline's added loss minus the changed system's, draw by draw. Runs of different systems on
  the same cells, folds and seed share their resamples; the tool checks this and says when a contrast is not paired.
- **Share removed**: one minus the ratio of the two, with its interval.
- The span of a loop runs from the stage before the located step to the system's output. A repair counts by what reaches
  the output. Single steps are reported beside it.

## Readings

`diagnose`, per readout:

- **located step confirmed**: the baseline adds a resolved loss over the span; the change at the located step removes a
  resolved part of it, not confined below the step margin; and it removes more than the control change (interval of the
  difference above zero). *Selective* if whatever the control change removes lies below the step margin.
- **located step not confirmed**: what the located change removes lies below the step margin, or the control removes more.
- **unresolved**, or **not testable** (the baseline adds no resolved loss over the span).

`closure`, per readout, strongest first:

- **absolute**: the end stage of the changed system is supported and the control's is rejected. This is the paper's
  closing of a loop (a rejected stage becomes supported, and the control's does not).
- **relative**: the changed system adds no loss over the span, by its interval, and the control adds a loss.
- **loss no longer resolved**: the changed system's span is unresolved. The loss is not shown. Its absence is not shown
  either. This is not a closure.
- **not closed**: the changed system still adds a resolved loss. **Not testable**: the control shows no loss to close.

## Gates

Each gate is passed before the next step is paid for.

| gate | what is checked | by | cost |
|---|---|---|---|
| 0 admission | The response is resolved and large against its interval. The reference stage is not rejected. The span adds a resolved loss as it stands. How narrow the intervals are: whether "absolute" or "relative" closure could be reached at all. The number of clusters. | `anm_diag loop admit` on the baseline run | none: one existing run |
| 1 location | Sealed declaration, run once. | `anm_diag seal`, `run` | one audit |
| 2 prescription | One change per round. A control: the same procedure without the change, or the same change made elsewhere. The level of generalization claimed (cells, sources, readouts, cell types never used by the change). The span, the readouts and the rule, written before any run. | a design file, committed first | none |
| 3 admission of the change | It converged (a stated criterion, not a step limit). It reached its own objective on the data it was fitted on. No audit cell entered the fit (checked by identifiers). The guard holds on validation data. | the round's fitting script | low; no audit |
| 4 development audit | `diagnose` and `closure` on development data. Rounds may repeat here. Everything read here is exploratory. | `anm_diag run`, `loop` | audits |
| 5 final audit | The recipe is frozen. At least three training seeds per arm. Sealed declaration, run once, on data no earlier gate has used. The outcome must hold in every seed. | `seal`, `run`, `loop closure` | training and audits |
| 6 use | The guard on the system's own task: the mean, and the worst single item. A comparison with a baseline that does not use the changed part. | the system's own evaluation | low |

A readout that fails gate 0 is not a target and is not counted in any outcome. A change that fails gate 3 was not tested.
It did not fail an audit.

## Rules

1. **Decide by intervals.** No outcome rests on a point value compared with a threshold.
2. **One change per round.** Two changes need a factorial design with every arm audited.
3. **Keep the audit's data for the audit.** Prefer a change that does not use the audited sources. If it must use them,
   cross-fit: every audited cell passes through a system fitted without it, and every control cell is used.
4. **Three sets of data.** Fit, development and final. Final data are touched once, under a sealed declaration. Without
   final data the loop stays exploratory, and says so.
5. **Replication is part of the decision.** The rule names the seeds or splits it will be read on.
6. **Stop.** A round states how many further rounds the development data will carry.
7. **One ledger.** Every change, audit and decision is a typed event (`loop.intervention_event`, `audit_event`,
   `decision_event`), so a run is not made twice and a decision can be traced to its runs.

## Evidence grades

Confirmatory: gate 5, sealed, enough clusters. Everything else is exploratory, including every re-reading of a finished
run under this protocol. A re-reading changes no outcome recorded under the rule a run was made with.

## Chosen numbers

| number | value | where |
|---|---|---|
| margin on *u*, step margin | 0.10 | each declaration |
| clusters for intervals that are not marked approximate | 10 | `loop.DEFAULTS` |
| response against the half-width of its interval | 3 | `loop.DEFAULTS` |
| span half-width for "relative" closure to be attainable | half the step margin | `loop.admit` |
| guard on the system's own task | 0.02 on the mean | [`declarations/teddy_mm_step4_intervention_rule.md`](declarations/teddy_mm_step4_intervention_rule.md) |

None of these is derived. Each is reported with every reading.

## What the first three rounds lacked

| lacking | consequence | gate or rule |
|---|---|---|
| a check that a repaired head could be recognized | on THP-1 the reference stage is unresolved, so no round could close there; on melanoma the intervals are about as wide as half the margin | gate 0 |
| readouts admitted by their response | HLA-E and MHC-I were made targets with a response too small to follow | gate 0 |
| a rule on intervals | round 3 "failed" by 0.009, inside every interval | rule 1 |
| a check of the fit | the fits stopped at a step limit; HLA-E was not learned even where it was trained on | gate 3 |
| one change per round | round 3 added a cell type and three outputs together | rule 2 |
| audit data kept for the audit | half the controls went to the fits, and the reference stage fell from supported to unresolved | rule 3 |
| final data | the same audit cells served every round | rule 4 |
| the sealed design of 9 November | three seeds, a matched control retrain and an interval rule were already written there | gate 5 |
