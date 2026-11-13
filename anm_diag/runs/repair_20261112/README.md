# Round 3 of the loop: one head for both cell types — 12 November 2026

Exploratory, and the third round after [`runs/loop_20261112`](../loop_20261112/) and
[`runs/loop_guard_20261112`](../loop_guard_20261112/). The design and the rule are in [`DESIGN.md`](DESIGN.md), committed
before any head of this round was fitted beyond a 20-step timing check and before any audit of this folder ran (commit
`0905d9e`). This file reports what came out.

## Outcome by the rule: not passed

After round 2 the audit still located two things, both at the head's MLP: HLA-E on melanoma, and a head repaired on one
cell type failing on the other. This round repaired the MLP again, in one fit for both cell types and on every measured
protein the head has an output for (PD-L1, MHC-I and HLA-E on melanoma; PD-L1 on THP-1), with the bone-marrow outputs
held as in round 2. Each of six target genes was held out in turn, in both data sets at once, and the audits use only the
held-out knockouts and controls that no fit saw.

| data set | PD-L1, pooled audit, split 0 | mean fall of the six heads on bone marrow (guard 0.02) | readouts still located in the head | outcome by the rule |
|---|---|---|---|---|
| melanoma | not resolved: the head adds 0.109 to *u* at Z (margin 0.10) | 0.005 | HLA-E | **regressed** |
| THP-1 | no loss in the head: the head adds 0.047 to *u* at Z (margin 0.10) | 0.005 | none | **kept; nothing located in the head** |

**One head serves both cell types: no.** On melanoma the head output for PD-L1 is at 0.165 against 0.056 at Z, so
the head adds 0.109; the rule's margin is 0.10, the reading is "not resolved", and the rule counts that as regressed.
In round 2 the same audit read "no loss in the head" (0.121, the head adding 0.064). The miss
against the margin is small. A direct comparison of the two heads on the same cells, with its interval
([`runs/loop_reading_20261112`](../loop_reading_20261112/)), finds the head of this round adding 0.044 [0.002, 0.090] more than
the head of round 2, so the direction of the reading holds. On THP-1 the round passes. The guard on the
original task holds: the six heads score 0.596 to 0.598 on the bone-marrow test cells, against 0.602 as trained.

**HLA-E on melanoma is still located in the head**, from Z to the head latent, with *u* at the output at 0.96
(1.04 in round 2, 1.11 as trained). Putting HLA-E into the loss did not repair it on held-out genes.

What follows, by the stopping rule of the design:

- The repair of round 2 stands: one head per cell type, trained on PD-L1, with the bone-marrow outputs held. It passed its
  rule on both cell types. One head for both cell types is not established.
- HLA-E on melanoma counts as not repaired by changing the MLP with these cells.
- No further variant of the fit is run on these audit cells.

The pooled audits of split 0, every readout (*u* with its interval and the verdict printed by the tool):

**Melanoma, IFN-γ** (1,171 knockout cells of 5 genes, 2,000 audit controls)

| readout | measured change | Z | head latent | head output | root cause printed by the tool | added by the head | reading |
|---|---|---|---|---|---|---|---|
| MHC-I | -0.14 | -0.43 [-1.03, 0.02] support | -0.10 [-0.59, 0.23] unresolved | -0.12 [-0.52, 0.20] unresolved | `unresolved` | +0.30 | not resolved |
| PD-L1 | -0.87 | 0.06 [-0.00, 0.11] unresolved | 0.13 [0.05, 0.20] unresolved | 0.17 [0.10, 0.23] unresolved | `reference_unresolved` | +0.11 | not resolved |
| HLA-E | -0.08 | -0.36 [-1.84, 0.26] unresolved | 0.71 [0.37, 1.48] reject | 0.96 [0.36, 1.50] reject | `relative_loss_at` Z → head latent | +1.32 | located in the head |

**THP-1** (5,092 knockout cells of 5 genes, 1,193 audit controls)

| readout | measured change | Z | head latent | head output | root cause printed by the tool | added by the head | reading |
|---|---|---|---|---|---|---|---|
| PD-L1 | -0.70 | 0.12 [0.01, 0.29] unresolved | 0.18 [0.07, 0.37] unresolved | 0.17 [0.03, 0.38] unresolved | `reference_unresolved` | +0.05 | no loss in the head |
| PD-L2 | -0.31 | 0.16 [-0.08, 0.47] unresolved | 0.04 [-0.33, 0.48] unresolved | 0.24 [-0.14, 0.67] unresolved | `reference_unresolved` | +0.08 | no loss in the head |

PD-L2 has no output in the head and was in no loss. MHC-I and HLA-E change little in these cells, and their intervals
are wide.

## Replication on two fresh splits

Splits 1 and 2 draw the audited half of the controls and the set-aside cells anew and repeat the six fits and the two
pooled audits, with the weight chosen at split 0. The same rule is read on each.

| split | data set | PD-L1: Z | PD-L1: head output | added by the head | reading | still located in the head | mean fall on bone marrow | outcome |
|---|---|---|---|---|---|---|---|---|
| 0 | melanoma | 0.06 | 0.17 | +0.109 | not resolved | HLA-E | 0.005 | regressed |
| 0 | THP-1 | 0.12 | 0.17 | +0.047 | no loss in the head | none | 0.005 | kept; nothing located in the head |
| 1 | melanoma | 0.03 | 0.14 | +0.111 | not resolved | HLA-E | 0.005 | regressed |
| 1 | THP-1 | 0.13 | 0.15 | +0.020 | no loss in the head | none | 0.005 | kept; nothing located in the head |
| 2 | melanoma | 0.01 | 0.16 | +0.140 | not resolved | none | 0.004 | regressed |
| 2 | THP-1 | 0.14 | 0.15 | +0.009 | no loss in the head | none | 0.004 | kept; nothing located in the head |

On the fresh splits melanoma passes in 0 of 2 and THP-1 in 2 of 2, as at split 0. HLA-E on melanoma reads
"located in the head" (0.96), "located in the head" (0.78), "not resolved" (0.87) on splits 0, 1 and 2.

## The loop, round by round

*u* at the head output and the reading, pooled audits on the same cells (split 0). Bone marrow: the head's own evaluation
on 5,000 test cells, lowest to highest of the round's heads.

| readout | Z | round 0: as trained | round 1: MLP repaired on PD-L1 | round 2: bone-marrow outputs held | round 3: one fit, both cell types |
|---|---|---|---|---|---|
| melanoma, MHC-I | -0.43 | 0.19 not resolved | -0.00 not resolved | -0.01 not resolved | -0.12 not resolved |
| melanoma, PD-L1 | 0.06 | 0.45 located in the head | 0.12 no loss in the head | 0.12 no loss in the head | 0.17 not resolved |
| melanoma, HLA-E | -0.36 | 1.11 located in the head | 0.91 located in the head | 1.04 located in the head | 0.96 located in the head |
| THP-1, PD-L1 | 0.12 | 0.54 located in the head | 0.16 no loss in the head | 0.18 no loss in the head | 0.17 no loss in the head |
| THP-1, PD-L2 | 0.16 | 0.65 located in the head | 0.15 no loss in the head | 0.15 no loss in the head | 0.24 no loss in the head |
| bone-marrow test score | | 0.602 | 0.502 to 0.567 | 0.594 to 0.604 | 0.596 to 0.598 |

Rounds 1 and 2 passed the rules written for them; round 3 did not. Read together: the audit located the PD-L1 loss in the
MLP; repairing the MLP removed it and repairing the decoder did not; the accuracy guard found the cost and round 2
removed that; the audit after round 2 located two more things at the same step, and the repair it suggested did not fix
them.

## Outside the rule

**The pilot that chose the weight** (IFNGR1 held out; 2,000 bone-marrow validation cells; the head as trained scores
0.7431 on them):

| weight | validation score | bone-marrow outputs moved (root mean square) | set-aside loss, first check to best | chosen |
|---|---|---|---|---|
| 1 | 0.7409 | 0.070 | 0.5911 to 0.2583 | yes |

The smallest weight of the grid met the condition (within 0.01 on validation cells), so the larger ones were not fitted.

**Held-out genes one at a time, PD-L1** (*u* at the head output and the reading; rounds 0 to 2 from the earlier runs):

| held-out gene | Z | round 0 | round 1 | round 2 | round 3 | bone marrow, this head |
|---|---|---|---|---|---|---|
| melanoma, IFNGR1 | 0.14 | 0.46 located in the head | 0.03 no loss in the head | 0.02 no loss in the head | 0.21 no loss in the head | 0.596 |
| melanoma, IFNGR2 | -0.02 | 0.52 located in the head | 0.13 not resolved | 0.13 not resolved | 0.20 not resolved | 0.596 |
| melanoma, JAK1 | -0.04 | 0.41 located in the head | -0.04 no loss in the head | -0.04 no loss in the head | 0.10 not resolved | 0.597 |
| melanoma, JAK2 | 0.18 | 0.50 located in the head | 0.11 no loss in the head | 0.16 no loss in the head | 0.24 no loss in the head | 0.598 |
| melanoma, STAT1 | 0.04 | 0.52 located in the head | 0.17 not resolved | 0.17 not resolved | 0.27 not resolved | 0.597 |
| THP-1, IFNGR1 | 0.01 | 0.55 located in the head | 0.14 not resolved | 0.14 not resolved | 0.13 not resolved | 0.596 |
| THP-1, IFNGR2 | 0.14 | 0.60 located in the head | 0.15 no loss in the head | 0.14 no loss in the head | 0.15 no loss in the head | 0.596 |
| THP-1, JAK2 | 0.15 | 0.64 located in the head | 0.21 no loss in the head | 0.23 no loss in the head | 0.22 no loss in the head | 0.598 |
| THP-1, STAT1 | -0.16 | 0.36 located in the head | 0.02 not resolved | 0.01 not resolved | -0.05 not resolved | 0.597 |
| THP-1, IRF1 | 0.33 | 1.12 located in the head | 0.60 not resolved | 0.60 not resolved | 0.58 not resolved | 0.597 |

The tool names a head step for 0 of the 10 held-out genes; 4 read "no loss in the head" and 6 "not resolved"
(round 2: 0, 5 and 5). Not settled, with no head step named: melanoma IFNGR2 (output 0.20, -0.02 at Z); melanoma JAK1 (output 0.10, -0.04 at Z); melanoma STAT1 (output 0.27, 0.04 at Z); THP-1 IFNGR1 (output 0.13, 0.01 at Z); THP-1 STAT1 (output -0.05, -0.16 at Z); THP-1 IRF1 (output 0.58, 0.33 at Z).
Where *u* at Z is itself high (THP-1 IRF1), the problem lies before the head.

**The head's own outputs on knockouts**, next to the measured protein. Shift of the knockout cells against the audit
controls, in control standard deviations (`results/own_shift.json`). The first three columns are the knockouts of the
held-out gene; the last three are the knockouts of the other genes of the data set, which that head was fitted on (mean
over genes).

| held-out gene | readout | held out: measured | head as trained | head of this round | fitted on: measured | head as trained | head of this round |
|---|---|---|---|---|---|---|---|
| melanoma, IFNGR1 | PD-L1 | -1.24 | -0.51 | -1.54 | -1.17 | -0.50 | -1.49 |
| melanoma, IFNGR1 | MHC-I | -0.31 | -0.37 | -0.38 | -0.27 | -0.20 | -0.33 |
| melanoma, IFNGR1 | HLA-E | -0.08 | +0.78 | +0.96 | -0.15 | +0.65 | +0.80 |
| melanoma, IFNGR2 | PD-L1 | -1.14 | -0.49 | -1.44 | -1.19 | -0.51 | -1.51 |
| melanoma, IFNGR2 | MHC-I | -0.25 | -0.18 | -0.35 | -0.28 | -0.24 | -0.39 |
| melanoma, IFNGR2 | HLA-E | -0.19 | +0.66 | +0.82 | -0.12 | +0.68 | +0.86 |
| melanoma, JAK1 | PD-L1 | -0.96 | -0.40 | -1.39 | -1.23 | -0.53 | -1.60 |
| melanoma, JAK1 | MHC-I | -0.11 | -0.27 | -0.37 | -0.32 | -0.22 | -0.41 |
| melanoma, JAK1 | HLA-E | -0.11 | +0.61 | +0.77 | -0.14 | +0.69 | +0.88 |
| melanoma, JAK2 | PD-L1 | -1.25 | -0.54 | -1.55 | -1.16 | -0.50 | -1.52 |
| melanoma, JAK2 | MHC-I | -0.31 | -0.09 | -0.32 | -0.27 | -0.26 | -0.39 |
| melanoma, JAK2 | HLA-E | -0.06 | +0.66 | +0.86 | -0.15 | +0.68 | +0.84 |
| melanoma, STAT1 | PD-L1 | -1.31 | -0.59 | -1.69 | -1.15 | -0.48 | -1.53 |
| melanoma, STAT1 | MHC-I | -0.40 | -0.24 | -0.44 | -0.25 | -0.23 | -0.37 |
| melanoma, STAT1 | HLA-E | -0.23 | +0.67 | +0.84 | -0.11 | +0.68 | +0.88 |
| THP-1, IFNGR1 | PD-L1 | -0.94 | +0.20 | -1.21 | -0.86 | +0.19 | -1.08 |
| THP-1, IFNGR2 | PD-L1 | -1.07 | +0.16 | -1.32 | -0.83 | +0.20 | -1.05 |
| THP-1, JAK2 | PD-L1 | -1.15 | +0.17 | -1.39 | -0.81 | +0.20 | -1.01 |
| THP-1, STAT1 | PD-L1 | -0.86 | +0.33 | -1.28 | -0.88 | +0.16 | -1.05 |
| THP-1, IRF1 | PD-L1 | -0.37 | +0.11 | -0.26 | -1.01 | +0.22 | -1.26 |

For PD-L1 the repaired head's output now falls on knockouts it never saw, by 0.71 to 1.49 times the measured fall. For
HLA-E the head's output rises on the knockouts (+0.61 to +0.78 as trained) while the measured protein falls slightly
(-0.23 to -0.06). After the repair it still rises, on the held-out gene (+0.77 to +0.96) and equally on the genes
the head was fitted on (+0.80 to +0.88). So the fit did not learn the HLA-E response even where it was trained on it, within
its 400 steps. The reason was not examined.

**Untrained outputs.** On the audited cells the head's untrained outputs moved by 0.08 to 0.09 (melanoma) and 0.17 to
0.19 (THP-1) root mean square.

**The three trained proteins on bone marrow** (correlation on the 5,000 test cells, as trained and lowest to highest of the
six held-out heads): PD-L1 0.41, then 0.39 to 0.41; MHC-I 0.58, then 0.58 to 0.58; HLA-E 0.33, then 0.30 to 0.31.

**Floor check.** Audits with *u* at the head output more than 0.10 below *u* at Z: none.

**Fits.** Kept step: 400 of at most 400. Every fit ran to the step limit with the set-aside loss still falling, so the repair is not converged.

## The repair as a teddy_mm script

[`scripts/07_repair_head.py`](../../../scripts/07_repair_head.py) makes the repair with teddy_mm's own model code
(PyTorch) and writes a checkpoint in the format of `scripts/04_train.py`. It was added after the design and is not part
of the rule. Run on the same inputs, it gives the heads that were fitted here in NumPy (`results/torch_check.json`;
decoder identical to the head as trained in every case):

| form | compared head | outputs: NumPy against PyTorch (root mean square) | outputs: change made by the repair |
|---|---|---|---|
| round 2: melanoma only, PD-L1, IFNGR1 held out | `runs/loop_guard_20261112` | 0.00043 | 0.065 |
| round 2: THP-1 only, PD-L1, IFNGR1 held out | `runs/loop_guard_20261112` | 0.00034 | 0.201 |
| round 3: joint, all genes, on melanoma cells | `heads/s0__all` of this run | 0.00007 | 0.089 |
| round 3: joint, all genes, on THP-1 cells | `heads/s0__all` of this run | 0.00016 | 0.203 |

Heads made with the script with no gene held out, in teddy_mm's own evaluation (`scripts/05_eval.py`) on the full
bone-marrow test split of 16,750 cells (`results/full_test_eval.json`):

| head | mean correlation, 134 proteins | PD-L1 | MHC-I | HLA-E | proteins lower by more than 0.02 | lowest single change |
|---|---|---|---|---|---|---|
| as trained | 0.6027 | 0.396 | 0.574 | 0.322 | | |
| round 2 form, melanoma, PD-L1 | 0.6043 | 0.394 | 0.581 | 0.317 | 0 | -0.011 |
| round 2 form, THP-1, PD-L1 | 0.5955 | 0.397 | 0.566 | 0.290 | 8 | -0.051 |
| round 3 form, joint (did not pass its rule) | 0.5976 | 0.392 | 0.573 | 0.295 | 9 | -0.050 |

None of these three heads is audited: each was fitted on the knockout cells an audit would use, and the evidence for a
form is the audits of its held-out fits. The mean stays within 0.02 for all three; single proteins move further for the
THP-1 and joint heads. The checkpoint's flow-matching branch (`v_net`) was trained to reproduce the latent of the MLP as
it was; the script copies it unchanged, so only the direct route (MLP, then decoder) is repaired, and the flow-matching
score does not change (0.5821).

## What this round adds, and what it does not

- The loop now has a round that did not pass. The audit that located the loss and confirmed the first repair also read a
  wider repair as no better on melanoma, by a rule fixed before the run.
- A wider fit is not a better one here: with three more outputs and a second cell type in the same loss, PD-L1 on
  melanoma ends slightly further from Z than with the narrower fit of round 2.
- HLA-E and MHC-I were not targets that could be followed: on these cells their response is small against its own
  interval, and they do not pass the admission check of [`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md), written afterwards.
- HLA-E is not a PD-L1-like case. The measured change is small (-0.08), the head's own output moves the other way
  before and after, and the fit did not learn it even on the genes it was trained on. This round therefore shows that this
  fit does not repair HLA-E; it does not show that the MLP cannot carry an HLA-E response. Whether the embedding itself
  carries it is not settled either (Z is "unresolved" with a wide interval).
- Transfer to a cell type outside the fit is not tested by this round, and one head for the two cell types in the fit is
  not established.

## Limits

- Designed after rounds 1 and 2 were read. The same audit controls serve for the fourth time at split 0; splits 1 and 2
  redraw them.
- Two things changed at once against round 2 (a second cell type in the fit, three more trained outputs), so their effects
  cannot be told apart from this round.
- The rule compares a point value with the margin. Read by intervals in `runs/loop_reading_20261112`, the loss this
  round's head adds on melanoma is unresolved on all three splits, and it is 0.044 [0.002, 0.090] above the round-2
  head's in a paired comparison.
- Not the sealed step 4: a fine-tune of one part, the head as trained as comparator, splits instead of training seeds,
  five genes per data set as resampling units (rule `K1`), no positive-control node (rule `C5`), two cell types.
- The fits stopped at the step limit and are not converged. A longer or differently weighted fit is a new variant; by the
  stopping rule it is not run on these audit cells.

## Files

16 folders `pooled_…_joint_s<split>` and `target_…_joint_s0`: `declaration.yaml`, `anm_record.json`, `anm_report.md`.
`results/repair.json`: the audits, the readings, the rule, the replication and the earlier rounds for the same readouts.
`results/fits_s<split>.json`: the pilot, the fits and the side measures. `results/bone_marrow.json`: the test scores of
every head. `results/own_shift.json`, `results/torch_check.json`, `results/full_test_eval.json`: the side measures above.
`scripts/`: `repair.py` (plain NumPy; `--check` compares the gradient with finite differences; every fit asserts that the
decoder stays bit for bit), `make_declarations.py`, `run_one.py`, `run_all.py`, `score_bone_marrow.py`, `analyze.py`,
`own_shift.py`, `compare_torch.py`. Features, heads and the bone-marrow cells are not kept in git (`ANM_FEATURES`,
`ANM_REPAIR_FEATURES`, `ANM_BM_FIT`, `ANM_BM_TEST`).

```
cd scripts
python make_declarations.py
python repair.py pilot                      # chooses the weight; about 5 minutes on one core
python repair.py fit 0 0 & python repair.py fit 0 1 & wait      # six held-out fits of split 0, two processes
python repair.py final && python repair.py features 0
python run_all.py 0                         # 12 audits
# splits 1 and 2: repair.py fit <split> 0|1, repair.py features <split>, run_all.py <split>
python score_bone_marrow.py && python analyze.py && python own_shift.py
```
