# Closing the loop: does repairing the step the audit names repair the system? Design, written before the runs

Status: exploratory. This file, the scripts, the analysis and the 40 declarations were committed before any head was
fitted beyond a 20-step timing check on training cells, and before any audit was run. The commit history shows the order.

## What the audit located

The two sealed head-fix runs read the PD-L1 response to knockouts of the IFN-γ pathway at three stages (quick reruns in
`runs/split_20261111`):

| cells | Z (TEDDY's embedding) | head latent | head output | root cause printed by the tool |
|---|---|---|---|---|
| melanoma, IFN-γ (Frangieh 2021) | 0.04 support | 0.20 reject | 0.46 reject | `loss_at` Z → head latent |
| THP-1 (Papalexi 2021) | 0.15 unresolved | 0.35 reject | 0.47 reject | `relative_loss_at` Z → head latent |

The response is in TEDDY's embedding and is first lost between Z and the head's latent, that is in the head's first
network (the MLP). The decoder comes after that step.

## What that location prescribes, and what it rules out

- **Prescribed:** retrain the MLP. TEDDY need not change, because Z carries the response.
- **Ruled out:** retraining the decoder alone. It only sees the latent, so it cannot bring back what the latent has lost;
  by the audit's reading its output must stay at least as unreadable as the latent is now.

Looking only at the endpoint (the PD-L1 output is wrong) one would do the opposite and refit the output side. The two
repairs are run side by side on the same cells, with the same loss and the same budget. If the audit's location is right,
the first closes the loop and the second does not. If the loss at the latent were an artefact of the readout models
(information present but not read), the decoder-only repair would succeed too and the location would be wrong.

## Data, split, fits

The exported features of `runs/split_20261111`: melanoma (4,000 controls; IFNGR1, IFNGR2, JAK1, JAK2, STAT1) and THP-1
(2,386 controls; IFNGR1, IFNGR2, JAK2, STAT1, IRF1). TEDDY is frozen.

- **Audit controls.** Half of the control cells (at random within batch, seed 0) are never used by any fit.
- **Held-out genes.** For each data set and each of its five targets a head is fitted on the other four targets and the
  training controls. The held-out target's cells, all of them, are audited with that head: knockouts of a gene the head
  has never seen. One more head per data set is fitted on all five targets for the transfer audits.
- **Two repairs.** `mlp_only`: only the three layers of the latent network change; the decoder stays bit for bit.
  `decoder_only`: only the three decoder layers change; the latent stays bit for bit. The script asserts both.
- **Fit.** As in `runs/exposure_20261111`: squared error of the PD-L1 output (CD274) against measured PD-L1 rescaled to the
  output's mean and spread on training controls, plus the mean squared change of the other 133 outputs with weight 1.
  Adam, learning rate 3e-4, full batch, at most 400 steps; checked every 20 steps on a set-aside fifth of the training
  cells, stopped after three checks without improvement, best weights kept. Nothing is tuned on audit cells.

## The audits (40)

All with the branch's `anm_diag`, stages Z, head latent, head output, the folds, bootstrap, margins (0.10) and seed of the
sealed runs.

1. **Pooled, the primary audit** (2 data sets × 3 heads: as trained, `mlp_only`, `decoder_only`). The pooled declaration
   of `runs/split_20261111` on all knockout cells and the audit controls. Every knockout cell is passed through the head
   that was fitted with its gene held out; audit controls are assigned at random to one of the five heads. All readouts.
2. **Single head, the check on 1** (2 × 5 targets × 3 heads). One held-out gene, one head, all audit controls; PD-L1;
   resampling by guide.
3. **Transfer** (2 × 2 repairs). The head fitted on all targets of one cell type, audited on all cells of the other.

## Reading and rule

For one readout of one audit, from the tool's own output:

- **located in the head**: the root cause is `loss_at` or `relative_loss_at` (every step after Z is a head step);
- **no loss in the head**: the root cause is neither, and u at the head output exceeds u at Z by at most 0.10;
- **not resolved**: anything else.

Rule, per data set, on PD-L1 in the pooled audit:

| head as trained | `mlp_only` | `decoder_only` | outcome |
|---|---|---|---|
| not "located in the head" | | | the question cannot be put on these cells |
| located | no loss | not "no loss" | **loop closed at the located step** |
| located | no loss | no loss | closed, but not specific to the located step |
| located | not "no loss" | no loss | closed by the other part: the location was wrong |
| located | not "no loss" | not "no loss" | not closed |

## Also reported, outside the rule

- The single-head audits: how many of the five held-out genes fall in each reading, per head.
- Readouts the repair was not trained on (melanoma: MHC-I, HLA-E; THP-1: PD-L2), in the pooled audit. If the located
  step is the cause, repairing it may repair them too; refitting the decoder on PD-L1 should not.
- Transfer. Expected not to hold: three earlier results (heads retrained with IFN-γ cells and audited on co-culture
  cells; the fine-tune of `runs/exposure_20261111`) improved only where they were trained.
- A floor check: the head's only input is Z, so u at the head output should not be more than 0.10 below u at Z in any
  audit. Violations are listed; they would point to leakage or to an artefact of pooling heads.
- For every held-out gene: the shift of the head's own PD-L1 output next to the measured one, the change of the other
  133 outputs, and whether the latent is identical to the head as trained.
- Cost on the head's original task (bone-marrow cells): measured if those cells can be exported from the machine that
  holds them, and reported either way.

## Audit of this plan: what was wrong with earlier versions, and what remains

| problem found | what the design does about it |
|---|---|
| The first idea was to hold out a whole cell type. With two cell types that tests transfer, not the location, and transfer has failed three times. | The location is tested within a cell type on held-out genes; transfer is a separate, secondary audit. |
| A pass line of u < 0.10 at the head output cannot be met where Z itself is not supported (THP-1). | The reading is the tool's root cause plus the loss the head adds beyond Z. |
| The earlier fine-tune audited a random half of the cells, from guides the head had seen. | Held-out genes; controls no head was fitted on. |
| A head trained on the audited readout might pass for that reason alone. | `decoder_only` is trained on the same readout, cells, loss and budget and is predicted not to pass; readouts not trained on are reported. |
| The pooled audit mixes features from five heads; a readout model could pick up which head a cell went through. | Controls are assigned to heads at random, so the head is not the treatment; the single-head audits are the check; the floor check would show an artefact. |
| A repair can damage everything else the head does. | The other 133 outputs are held in the loss and their change is reported; the bone-marrow cost is attempted. |
| The designer knew related results: an earlier run in another workspace (4b, not in this repository) found on the IFN-γ cells that retraining the MLP removed the loss and retraining the decoder did not; the THP-1 fine-tune was seen yesterday. | Graded exploratory. What is new here: held-out genes, a second cell type, a rule written first, and code and records in the repository. |
| Remaining: five genes per data set as resampling units (rule `K1`), one split of the controls, no positive-control node in the quick declarations (rule `C5`), two cell types only, PD-L1 as the one trained readout. | Stated with the results. |
