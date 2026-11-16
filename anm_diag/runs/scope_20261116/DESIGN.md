# Two finished pieces of work re-read with the version 2 tool — design

Written 16 November 2026, before any script in this folder was run on data. Exploratory. No model is trained. Nothing
new is audited: the cells, the head and the measurements are those of earlier runs.

The question is whether version 2 of the tool (`anm_diag scope`, hypotheses and tests; `LOOP_PROTOCOL.md`, version 2)
is general: one code and one file format for a perturbation and for a contrast nobody assigned, giving back what is
already known and agreeing with interventions that were already made.

## What is re-read

| | event | cells | readout | earlier work |
|---|---|---|---|---|
| K1 | knockouts of five IFN-γ pathway genes, melanoma | the sealed audit set of `runs/frangieh2021_headfix` (5,171 cells) | PD-L1 protein | located in the head; repaired in `runs/loop_20261112` and after |
| K2 | the same pathway, THP-1 | the sealed audit set of `runs/papalexi2021_headfix` (7,478 cells) | PD-L1 protein | the same |
| D1 | a donor against the other donors, bone marrow | 90,261 cells, 9 donors, 12 donor × site groups | 112 proteins | `runs/generalization_20261112` |
| D2 | a donor against the other donors, external PBMC | 40,000 cells, 8 donors | 112 proteins | the same |

Stages: `Z` (TEDDY's pooled state, the stage the head reads), `HEAD_LAT` (the output of the head's MLP), `HEAD_OUT`
(the head's 134 outputs, log scale). For donors also `RNA512`, 512 principal components of log-normalized counts, as
the input view. `frozen`: the head's own prediction of the readout, as trained.

Scopes: `frozen`; `train` (bone-marrow cells the head was trained on, 67,405; a donor is left out of its own fit);
`other` (for K1 the cells of K2 and the reverse); `units` (the other genes, or the other donors); `own`.

Readout classes: ridge with the penalties of the sealed declarations; for the pooled knockout events also the boosting
class of those declarations. Donor fits use at most 30,000 cells each (a seeded subsample). 1,000 resamples, seed 0,
margin 0.10, step margin 0.10.

## Three criteria, written in advance

**G1. One code, one format.** All four are run by `python -m anm_diag scope SPEC.yaml`. The scripts of this folder only
write arrays and SPEC files. The package has no branch on the kind of event outside the checks of
`events.check_contrast`.

**G2. What is already known comes back.**

- Knockouts, pooled event, `units` scope (the scope of version 1):
  - melanoma: `Z` is not rejected and `HEAD_OUT` is rejected;
  - both cell types: `HEAD_OUT` is rejected, and the loss added from `Z` to `HEAD_OUT` (|u| at `HEAD_OUT` minus |u| at
    `Z`, paired) has its interval above 0.10.
- Donors, `variation` target (the head's ordinary score), at `Z`, in each of the 8 external donors:
  - the gap from `frozen` to `train` is smaller than the gap from `train` to `units` (the head is not the bottleneck;
    cells of the study are);
  - the gap from `train` to `units` is positive;
  - the gap from `units` to `own` is smaller than the gap from `train` to `units`.

**G3. The grid agrees with interventions already made.** For each cell type a hypothesis is written from the pooled
knockout grid alone and tested against the recorded audits of the repaired heads. It is marked retrospective: the
outcomes were known.

- Claim: **coverage**. What the head lacks for this response is cells of this kind, not another map on bone-marrow
  cells and not cells of another cell type.
- Quantity, for every arm: the fall of |u| at `HEAD_OUT` in the audit of the changed head (version 1 statistic).
- Forecast for an arm that learns from cells of scope *s*: |u| at `HEAD_OUT` as trained minus |u| at `Z` in scope *s*
  (both from the grid, paired). `coverage` arm: *s* = `units`. `other study` arm: *s* = `other`. `function` arm:
  *s* = `train`; no head was ever retrained on bone marrow alone, so this arm has a forecast and no test.
- Rule: confirmed if the `coverage` arm's fall lies above 0.20 and the `other study` arm's fall lies below 0.20;
  refuted if the `coverage` arm's fall lies below 0.20 or the `other study` arm's lies above it.
- Arms as recorded: `coverage` is the MLP repair fitted on knockout cells of the same cell type with the audited genes
  held out (`runs/loop_20261112`, pooled); `other study` is the same repair fitted on the other cell type
  (`transfer_*`). For the second, the baseline is the sealed run at its point value, because its resamples were not
  kept; that interval is approximate.
- Each arm is also scored against its forecast interval.

A hypothesis for the donors is written from the D2 grid and left open (not retrospective): the forecast for the three
arms of a training experiment that has not been run. Its rule is fixed here: own arm `coverage`; confirmed if its gain
lies above half the lower bound of its forecast and the gains of `function` and `amount` lie below that value.

## A correction this re-reading has to make

`runs/generalization_20261112` reported how much of a donor's shift in protein level the head leaves unexplained.
Two things were wrong in how that was computed, and both are fixed by the one statistic of version 2:

1. It was a share of squared differences, set against a margin (0.10) that is defined for *u* itself.
2. For the external donors, the head's prediction (trained on the bone-marrow protein scale) was subtracted from a
   measurement on another scale without a slope.

The re-reading reports the corrected values beside the old ones. The conclusions of that run that rest on correlation
within a donor are not touched by this.

## What would count against the tool

Any of G1 to G3 failing. They are reported as they fall.

## Not done here

No training, no repair, no new audit under a sealed declaration, no final data.
