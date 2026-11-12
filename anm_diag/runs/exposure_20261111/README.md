# Does the head lose IRF1 because of what it was trained on? — 11 November 2026

Exploratory. The design and the rule are in [`DESIGN.md`](DESIGN.md), committed before the runs. This file reports what came out.

## Outcome by the rule: not answered

The head fine-tuned with IRF1 knockout cells is still rejected at its output for IRF1 and PD-L1, with *u* above half of the
value for the head as trained. By the rule written in advance that leaves the question unanswered.

| IRF1 → PD-L1, audit half (539 knockout cells, 4 guides) | Z | head latent | head output |
|---|---|---|---|
| head as trained | 0.49 [0.27, 0.81] unresolved | 0.95 [0.81, 1.68] reject | 1.13 [0.99, 1.82] reject |
| fine-tuned without IRF1 cells | same | 0.65 [0.47, 0.98] reject | 0.72 [0.54, 0.97] reject |
| fine-tuned with IRF1 cells | same | 0.59 [0.38, 0.83] reject | 0.66 [0.40, 0.87] reject |

| JAK2 → PD-L1, audit half (513 knockout cells, 4 guides) | Z | head latent | head output |
|---|---|---|---|
| head as trained | 0.22 [0.16, 0.29] reject | 0.53 [0.46, 0.58] reject | 0.66 [0.61, 0.70] reject |
| fine-tuned without IRF1 cells | same | 0.32 [0.24, 0.39] reject | 0.36 [0.30, 0.41] reject |
| fine-tuned with IRF1 cells | same | 0.30 [0.22, 0.38] reject | 0.36 [0.30, 0.41] reject |

Half of *u* for the head as trained is 0.57; the head fine-tuned with IRF1 cells has 0.66.

## What the numbers show, outside the rule

- **Most of the gain comes without IRF1 cells.** Fine-tuning on THP-1 controls and the four upstream knockouts takes IRF1
  from 1.13 to 0.72; adding IRF1 knockout cells takes it to 0.66, well inside the intervals. For JAK2 the two fine-tuned
  heads are equal (0.36). The reading that the head has to see small-footprint knockouts is not what these numbers point to;
  that the head had never seen THP-1 cells with their PD-L1 explains more.
- **What remains is mostly already missing at Z on these cells.** After fine-tuning, the head output is 0.17 to 0.23 above Z
  for IRF1 (0.64 above for the head as trained) and 0.14 above Z for JAK2 (0.43 before). This subtraction was not in the design.
- **The head's own PD-L1 output** (three splits). As trained it moves the wrong way for all five targets. Fine-tuned without
  IRF1 cells, it falls for IRF1 knockout cells by 0.61 to 0.82 of the measured fall; with IRF1 cells, by 0.72 to 1.00. For the
  four upstream knockouts both fine-tuned heads overshoot (1.15 to 1.71 of the measured fall).
- **The group shift over the spread at the head output** barely changes for IRF1 (0.08 to 0.11 as trained, 0.10 to 0.13
  fine-tuned): one output of 134 was trained, and the other 133 moved by 0.07 on average.

## Limits

- On half the cells every *u* is larger than in the full runs of `runs/split_20261111` (Z for JAK2: 0.22 here, 0.13 there;
  Z for IRF1: 0.49 here, 0.33 there): with fewer cells the readout models learn less. Compare heads on the same cells, not
  these values with the full runs.
- One data set, one split for the audit, four guides as resampling units (rule `K1`), training and audit cells from the same
  guides. The fine-tune is small: 200 steps, stopped at the step limit in five of six fits while the loss on the set-aside
  cells was still falling.
- The rule's outcome is "not answered". The observations above were read after the result and are not a test.

## Files

`<head>_<target>/` (6): `declaration.yaml`, `anm_record.json`. `results/exposure.json`: the six runs and the rule.
`results/finetune.json`: fitting curves and the side measures for split seeds 0 to 2. `results/audit_half_seed0.json`:
the cells of the audit half. `scripts/`: `finetune_head.py` (plain NumPy; `--check` compares the gradient with finite
differences), `make_declarations.py`, `run_one.py`, `analyze.py`. Features of the three heads are written to
`ANM_EXPOSURE_FEATURES` and are not kept in git.

```
cd scripts
python make_declarations.py
python finetune_head.py            # about ten minutes on one core
python run_one.py as_trained_IRF1  # and the five others; about two minutes each
python analyze.py
```
