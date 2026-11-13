# Round 3 of the loop: one repaired head for both cell types. Design, written before the runs

Status: exploratory. It continues [`runs/loop_20261112`](../loop_20261112/) and
[`runs/loop_guard_20261112`](../loop_guard_20261112/) and was designed after their results were read. This file, the
scripts and the 16 declarations were committed before any head of this round was fitted beyond a 20-step timing check on
training cells, and before any audit of this folder was run. The commit history shows the order.

## The loop so far, and what the audit still locates

Each round is: audit, read where the tool locates a loss, change that place, audit again on cells the change never saw.

| round | change | pooled audit afterwards (held-out genes) | left over |
|---|---|---|---|
| 0 | none: the head as trained | PD-L1 lost in the head on both cell types, from Z to the head latent (the MLP); also HLA-E (melanoma) and PD-L2 (THP-1) | — |
| 1 | repair the MLP on PD-L1 (`loop_20261112`) | PD-L1: no loss in the head, both cell types; repairing the decoder instead changes nothing | bone-marrow score falls by 0.035 to 0.100, found by the accuracy guard |
| 2 | the same with the head's outputs held on bone-marrow cells (`loop_guard_20261112`) | PD-L1 and PD-L2: no loss in the head; bone-marrow score within 0.02 | **HLA-E on melanoma is still located in the head** (Z → head latent); **a head repaired on one cell type is still located in the head on the other** |

Both left-over findings name the head, and for both the step the tool names is the MLP. So the prescription is the same
step again, with what was missing from the fit: the other cell type, and the other measured proteins.

## The repair of this round

`joint`: the `mlp_guarded` repair of round 2, in one fit for both cell types.

- **Cells.** The knockout cells and training controls of both data sets together. One of the six target genes (IFNGR1,
  IFNGR2, JAK1, JAK2, STAT1, IRF1) is held out of the fit in both data sets at once.
- **Trained outputs.** Every measured protein the head has an output for: on melanoma PD-L1 (`CD274`), MHC-I
  (`HLA-A-B-C`) and HLA-E; on THP-1 PD-L1. PD-L2 has no output and stays untrained. Targets as before: the measured
  protein rescaled to the mean and spread of the head's own output on that data set's training controls. Each (data set,
  protein) pair is one mean squared error and the four are summed.
- **Kept as before.** Only the three MLP layers change (the decoder stays bit for bit, asserted); the head's other outputs
  on these cells are held with weight 1; its outputs on 3,000 bone-marrow training cells are held with a weight chosen as
  in round 2 (grid 1, 10, 100; the smallest whose score on 2,000 bone-marrow validation cells is at most 0.01 below the
  head as trained; first gene held out; no audit cell and no test cell is read). Adam, learning rate 3e-4, at most 400
  steps, the same stopping rule.
- **Splits.** Seed 0 is the split of rounds 1 and 2: the same audit controls (checked against their record). Seeds 1 and
  2 draw the audited half of the controls, the set-aside cells and the set-aside bone-marrow cells anew, with the weight
  chosen at seed 0. The fits have no other randomness, so a "seed" here is a split.

## The audits (16) and the guard

- Pooled, per data set and split (6): the pooled declaration of round 1, all readouts; every knockout cell through the
  head fitted with its gene held out; audit controls assigned to that data set's five heads at random.
- Single head, seed 0 (10): one held-out gene, its head, all audit controls; PD-L1.
- **Guard:** the head's own evaluation on the 5,000 bone-marrow test cells of rounds 1 and 2, for the six heads of a split.

## Reading and rule

Readings as in rounds 1 and 2 ("located in the head", "no loss in the head", "not resolved"). Per data set, from the
pooled audit of seed 0:

| PD-L1 | mean fall of the six heads on bone-marrow test cells | readouts still "located in the head" | outcome |
|---|---|---|---|
| no loss in the head | at most 0.02 | none | **kept; nothing located in the head** |
| no loss in the head | at most 0.02 | some | kept; still located in the head: (which) |
| anything else, or | more than 0.02 | | regressed |

**One head serves both cell types** if neither data set is "regressed". The same reading is made for seeds 1 and 2 and
reported as replicated or not, per data set.

**Stopping.** This round uses what the located step offers (the MLP, both cell types, every measured protein). A readout
still located in the head afterwards is reported as not repaired by changing the MLP with these cells; no further variant
of the fit is run on these audit cells. What the single-gene audits leave unresolved with no head step named is listed
with its value at Z, as a problem outside the head.

## The head to use

After the held-out fits: one fit on all six genes at seed 0, with the same weight. It is the head a user would take. It is
scored on the bone-marrow test cells; it is not audited on these knockout cells, because it was fitted on them. The
evidence for it is the held-out audits above.

## Also reported, outside the rule

- The pilot; the fitting curves; how far the untrained outputs moved; the shift of the head's own outputs next to the
  measured ones for each held-out gene.
- The single-gene audits; the floor check; the correlation on bone marrow for the three trained proteins.
- A script that makes the same repair with teddy_mm's own model code (`scripts/07_repair_head.py`), and how closely its
  head matches the one fitted here. It is added after this design and is not part of the rule.

## Audit of this plan

| problem | what the design does about it |
|---|---|
| Two things change at once (both cell types in one fit, three more trained proteins), so their effects cannot be told apart. | Stated. The aim of the round is one head to use; both changes act on the same located step. |
| The same audit controls are used for the fourth time, and each design was written after the previous audit on them. | Seeds 1 and 2 repeat the pooled audits on control halves that were drawn anew. The rule is read at seed 0 and the replication is reported next to it. |
| A head trained on both cell types does not show transfer to a third. | Stated. No third cell type is audited here. |
| HLA-E and MHC-I change little in these cells (−0.08 and −0.14), so their verdicts may stay unresolved whatever the head does. | The rule counts only "located in the head" as left over and reports the other readings as they are. |
| The weight of the holding term could be tuned to pass. | Chosen on bone-marrow validation cells at seed 0 by the rule of round 2, before any audit; reused for seeds 1 and 2. |
| The head to use is fitted on every gene, so no audit of it is clean. | It is not audited. The held-out heads carry the evidence. |
| Holding the bone-marrow outputs says nothing about the head's outputs on other tissues. | Stated as a limit of the head to use. |
| Not the sealed step 4: a fine-tune, the head as trained as comparator, splits instead of training seeds, five genes per data set as resampling units (rule `K1`), no positive-control node (rule `C5`), two cell types. | Stated with the results. |
