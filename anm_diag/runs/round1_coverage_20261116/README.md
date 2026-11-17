# Round 1: a head for donors it has not seen — development (gates 3, 4 and 6) — 16 November 2026

**Exploratory.** The design is [`DESIGN.md`](DESIGN.md), frozen before any script of the round existed. The forecasts ([`results/forecasts.json`](results/forecasts.json)) and the two scripts that read the heads ([`scripts/read.py`](scripts/read.py), [`scripts/grid_changed.py`](scripts/grid_changed.py)) were committed before any fitted head was scored on a PBMC cell. The eight PBMC donors below are the donors the forecast was read on, so **this is development, not the test**. The final study F has not been opened and gate 5 has not been run.

![development](results/development.svg)

*u* is the share of the variation among a donor's cells that the head's output leaves unexplained: 1 minus the mean over 112 proteins of r² within the donor. Every donor is scored only by heads that were not fitted on it. A gain is the fall of *u* from the head fitted as registered (`as_is`) with the same seed, on the same cells.

## Outcome

- **45 of 45 heads are admitted** (gate 3).
- **Adding PBMC cells to the head's training lowers *u* on donors it was not fitted on by 0.054 [0.046, 0.062]**, from 0.636 to 0.582; in every seed (0.058, 0.053, 0.052). In the head's ordinary score, the mean correlation within a donor, that is 0.521 to 0.583.
- **It does so as much when the donor's batch was not in the fit**: 0.053 [0.044, 0.062]; the difference between the two is +0.001 [−0.005, 0.006].
- **No control does**: 20,000 fewer bone-marrow cells −0.004 [−0.009, −0.001]; a wider head trained twice as long −0.001 [−0.009, 0.005]; the PBMC proteins shuffled within the donor −0.341 [−0.353, −0.327].
- **The head's own task is not paid for it**: on bone marrow, site 4, the score of the heads with PBMC cells changes by +0.0005 on average; every head of the round stays within the guard of 0.02.
- **The go rule of gate 4 holds**: `coverage` is admitted in every seed (yes); its gain with the batch not seen has its interval above 0 (yes); the guard holds (yes). Gate 5, the one test on a study nobody has looked at, may be run once F is prepared and the owner of the repository approves opening it.

With eight donors every interval here is approximate.

## The arms

| arm | PBMC cells in the fit | *u* | fall of *u* from `as_is` | by seed | ordinary score | against the rule's bound 0.034 |
|---|---|---|---|---|---|---|
| `as_is`: the registered recipe | none | 0.636 | | | 0.521 | |
| `coverage`, the scored donor's batch in the fit (M1, M2) | 20,000 from 4 donors | 0.582 | +0.054 [0.046, 0.062] | +0.058, +0.053, +0.052 | 0.583 | above |
| `coverage`, the batch not in the fit (B1, B2) | 20,000 from 4 donors | 0.583 | +0.053 [0.044, 0.062] | +0.058, +0.052, +0.050 | 0.583 | above |
| `four_donors` | 5,000 from 4 donors | 0.600 | +0.036 [0.028, 0.044] | +0.040, +0.038, +0.030 | 0.566 | crosses it |
| `one_donor`, of the scored donor's batch | 5,000 from 1 donor | 0.604 | +0.033 [0.024, 0.040] | +0.038, +0.033, +0.026 | 0.564 | crosses it |
| `one_donor`, of the other batch | 5,000 from 1 donor | 0.626 | +0.010 [−0.004, 0.022] | +0.013, +0.010, +0.007 | 0.543 | below |
| `amount`: 20,000 fewer bone-marrow cells | none | 0.640 | −0.004 [−0.009, −0.001] | −0.004, −0.002, −0.007 | 0.516 | below |
| `function`: hidden width 1,024, 80 epochs | none | 0.638 | −0.001 [−0.009, 0.005] | +0.003, +0.001, −0.008 | 0.520 | below |
| `labels`: as `coverage`, proteins shuffled within the donor | 20,000 from 4 donors | 0.977 | −0.341 [−0.353, −0.327] | −0.335, −0.338, −0.349 | 0.004 | below |

The head in use (the phase-1 head, one seed) has *u* = 0.640 on the same cells; `as_is` refits its recipe with three seeds (0.638, 0.638, 0.632).

The rule of `donors_hao_coverage` asks the interval of `coverage` to lie above the bound and those of its controls `function` and `amount` below it. Read on these donors it would be met. These are the donors the forecast came from, so the rule is not decided here. `labels` is the control of a second hypothesis that gate 5 enters before F is opened.

## What development was read for

**D1. Does the arm deliver on donors not seen, and how much of the ceiling?** Yes: 0.054 [0.046, 0.062], 75% of the ceiling: what a ridge readout of the state, fitted on the same 20,000 PBMC cells, gains over the head in use (0.073 [0.066, 0.080]). In *u* itself the head with PBMC cells is at 0.582 and that readout at 0.568. One linear map fitted on bone-marrow and PBMC cells together, without a study label, was forecast to gain 0.036; the head gains more.

**D2. And when the batch was not seen?** The same: 0.053 [0.044, 0.062], 83% of that ceiling (0.064 [0.059, 0.071]; *u* 0.583 against the readout's 0.576; one linear map: 0.029). Paired by donor, seen minus not seen is +0.001 [−0.005, 0.006]. With four donors and 20,000 cells in the fit, whether the scored donor's batch is among them makes no difference that these data resolve.

**D3. Donors or cells?** Cells, as far as this resolves, with one exception. At 5,000 cells, four donors give 0.036 [0.028, 0.044] and one donor of the scored donor's batch 0.033 [0.024, 0.040]: +0.004 [−0.003, 0.010] between them, so more donors at the same number of cells add nothing resolved. Four times the cells from the same four donors add +0.018 [0.013, 0.023]. The exception: one donor of the other batch gives 0.010 [−0.004, 0.022], 0.022 [0.013, 0.034] less than one donor of the same batch; the loss is uneven (with the fit donor from batch 2, the donors of batch 1 gain −0.021 to 0.006; with it from batch 1, those of batch 2 gain 0.020 to 0.037). So with a single donor the batch matters, and with four donors and 20,000 cells (D2) it does not. No arm has 20,000 cells of one donor, or 5,000 cells of four donors of the other batch, so whether donors or cells remove the batch's part is not separated. The three contrasts of this paragraph are not in `scripts/read.py`; they were computed for this README from the per-donor values with the same function.

**D4. Levels or relation?** The relation. With the PBMC proteins shuffled among the cells of each donor, the head's ordinary score on donors it was not fitted on is 0.004 and *u* is 0.977: 0.395 [0.378, 0.412] worse than `coverage`, and 0.341 worse than `as_is`. A head that is shown PBMC cells whose proteins do not follow their RNA stops following the RNA on PBMC cells of other donors too, and keeps its bone-marrow score (below). Read: the head tells the two tissues apart in the state and maps each in its own way. The gain of `coverage` comes from the relation between RNA and protein in the added cells, not from their protein levels.

**D5. What it costs.** Nothing that the guard sees. On bone marrow, site 4 (never fitted), the mean change of the score is +0.0005 for `coverage`; on the 22 proteins the PBMC study does not measure, −0.0002. `amount` loses 0.013. It has the same bone-marrow cells as `coverage` and 30% fewer cells in all, so 30% fewer updates in its 40 epochs; `labels`, with the full number of cells and nothing to learn from the added ones, loses 0.003. So most of the loss of `amount` is not from the bone-marrow cells taken away. The number of updates is the likely reason; no arm tests it. The wider, longer `function` head gains +0.009 on site 4 (validation donor +0.016) and nothing on PBMC donors: more of the same map helps inside the study the head was trained on and not outside it. The largest fall of a single protein in any `coverage` head is CD69, −0.197, in one head of twelve.

**D6. What is left.** The grid of `donors_hao` on the changed head and on `as_is` of the same seed, on the four donors the changed head was not fitted on (fit sets M1 and M2, three seeds; a run has four donors, so numbers compare inside a pair, not with `runs/scope_20261116`).

| *u* of `variation`, mean over 8 donors and 3 seeds | `as_is` | `coverage` | fall |
|---|---|---|---|
| the head's own output | 0.636 | 0.582 | +0.054 [0.046, 0.062] |
| a readout of the head's training cells, at the state | 0.642 | 0.609 | +0.033 [0.025, 0.041] |
| a readout of the other scored donors, at the state | 0.574 | 0.574 | 0.000 [0.000, 0.000] |
| a readout of the donor's own cells, at the state | 0.557 | 0.557 | 0.000 [0.000, 0.000] |
| a readout of the other scored donors, inside the head | 0.572 | 0.570 | +0.002 [0.000, 0.004] |
| a readout of the other scored donors, at the head's output | 0.575 | 0.569 | +0.006 [0.003, 0.009] |

The training cells of `coverage` include its 20,000 PBMC cells, so the second row changes with the head. Rows three and four are readouts of the state, which is TEDDY's and the same for both heads.

- **The gap between the head and what the other donors can teach is mostly closed.** A readout of the other scored donors at the state reads 0.062 [0.055, 0.069] more than the head's own output before the change and 0.008 [0.003, 0.013] after it (computed for this README from the runs' files). That readout is the same for both heads: the state is TEDDY's and was not changed.
- **The changed head reads more than a linear readout of its own training cells.** Its own output leaves 0.027 [0.022, 0.031] less unexplained than a ridge readout of those cells at the state does; for `as_is` the difference is 0.006 [0.001, 0.011].
- **What only the donor's own cells carry is unchanged**: 0.017 [0.014, 0.021].
- Read at the head's output, a readout of the other scored donors gains +0.006 [0.003, 0.009] from the change, and inside the changed head it reads 0.006 more at the output than at the state.
- With the donor's own cells in the readout *u* is still 0.557. That part is not reached by any ridge readout of the state; how much of it is measurement noise of the proteins, these data do not say.
- At the tool's margin of 0.10 every donor reads "unreadable" in both heads (24 of 24 donor-by-seed readings each). That margin does not fit this target; the reading is in the gaps above. The `shift` target is written to the file for completeness and is not read: with three control donors it has almost no resolution.

For a next round this leaves, on these donors: a small remainder between the head and the other donors, and the part specific to a donor. Round 2 is not designed here; the test comes first.

## Against the forecasts

| | forecast, written before any fit | measured |
|---|---|---|
| `coverage`, hypothesis `donors_hao_coverage` (a ceiling) | 0.068 to 0.082 | 0.054 [0.046, 0.062]: **below the forecast**, above the rule's bound |
| `function`, the same hypothesis | −0.004 to 0.001 | −0.001 [−0.009, 0.005] |
| `amount`, the same hypothesis (a chosen band, not read from the grid) | −0.05 to 0.05 | −0.004 [−0.009, −0.001] |
| ceiling with the batch seen (gate 0c) | 0.073 [0.066, 0.080] | 75% of it delivered |
| ceiling with the batch not seen (gate 0c) | 0.064 [0.059, 0.071] | 83% of it delivered |
| one linear map on both studies, batch seen / not seen (gate 0c) | 0.036 / 0.029 | the head delivers 0.054 / 0.053 |

The forecast for `coverage` was read from a readout of the other seven donors (up to 30,000 cells); the arm has four donors and 20,000 cells. Of the ceiling of its own configuration, about a quarter is not delivered with the batch seen. The ceilings are gains over the head in use (*u* = 0.640); the arms' gains are over `as_is` of the same seed (mean 0.636).

## Gate 3: admission

| arm | heads | admitted | 1. validation score, last 5 epochs minus the 5 before (not below −0.005) | 2. own objective | 4. validation donor against `as_is` (not below −0.02) |
|---|---|---|---|---|---|
| `as_is` | 3 | 3 | +0.0008 to +0.0077 | — | — |
| `coverage` | 12 | 12 | −0.0028 to +0.0068 | score on its own PBMC fit cells above `as_is` by 0.071 to 0.093 (at least 0.035) | −0.0051 to +0.0080 |
| `four_donors` | 6 | 6 | −0.0014 to +0.0066 | — | −0.0133 to 0.0000 |
| `one_donor` | 12 | 12 | −0.0010 to +0.0074 | — | −0.0153 to −0.0011 |
| `amount` | 3 | 3 | +0.0002 to +0.0073 | — | −0.0161 to −0.0107 |
| `function` | 3 | 3 | −0.0048 to +0.0038 | training loss below `as_is` of the seed: 1.778 against 2.254, 1.807 against 2.194, 1.774 against 2.235 | +0.0118 to +0.0197 |
| `labels` | 6 | 6 | +0.0003 to +0.0096 | — | −0.0091 to +0.0021 |

Criterion 3 (no scored cell in a fit, by donor and barcode; bone-marrow fit cells are training cells) holds for every head. No cell of F is in any file a fit reads: F has not been downloaded.

The closest call of criterion 1 is `s0_function`: −0.0048 against a limit of −0.005.

## Gate 6: the guard on bone marrow, site 4

| arm | heads | mean change of the score | range over heads | the 22 proteins the PBMC study lacks | largest fall of one protein in one head |
|---|---|---|---|---|---|
| `coverage` | 12 | +0.0005 | −0.0093 to +0.0089 | −0.0002 | CD69 −0.197 |
| `four_donors` | 6 | −0.0023 | −0.0102 to +0.0050 | 0.0000 | CD279 −0.084 |
| `one_donor` | 12 | −0.0030 | −0.0126 to +0.0070 | −0.0014 | TCRVa7.2 −0.111 |
| `amount` | 3 | −0.0131 | −0.0163 to −0.0087 | −0.0137 | HLA-E −0.196 |
| `function` | 3 | +0.0089 | +0.0013 to +0.0178 | +0.0085 | CD223 −0.076 |
| `labels` | 6 | −0.0029 | −0.0151 to +0.0057 | −0.0002 | HLA-E −0.115 |

The score is the mean over 134 proteins of the correlation between measured and predicted on the site-4 cells; the change is against `as_is` of the same seed. The guard is 0.02; every head passes.

## Each donor

| donor | batch | *u*, `as_is` | *u*, `coverage`, batch in the fit | *u*, `coverage`, batch not in the fit | fall, batch in the fit | fall, batch not in the fit |
|---|---|---|---|---|---|---|
| P1 | 1 | 0.661 | 0.615 | 0.619 | 0.045 | 0.042 |
| P2 | 1 | 0.648 | 0.605 | 0.604 | 0.043 | 0.044 |
| P3 | 1 | 0.626 | 0.567 | 0.571 | 0.059 | 0.055 |
| P4 | 1 | 0.672 | 0.622 | 0.629 | 0.050 | 0.043 |
| P5 | 2 | 0.624 | 0.553 | 0.547 | 0.070 | 0.077 |
| P6 | 2 | 0.589 | 0.536 | 0.529 | 0.053 | 0.060 |
| P7 | 2 | 0.620 | 0.554 | 0.561 | 0.066 | 0.059 |
| P8 | 2 | 0.649 | 0.603 | 0.603 | 0.047 | 0.046 |

Means over the three seeds. The fall is positive in every donor.

## How the fits were run

- 45 heads, 11.0 hours of fitting on 2 CPU threads (torch 2.5.1+cu124), in seven jobs of at most two hours. Four times a job's time ended during a head and twice the machine restarted; each time the head went on from the state kept after its last finished epoch. The trainer's test shows a resumed fit gives the same weights as an uninterrupted one.
- The design asks that all arms of a seed are fitted on one machine. Every head records the same device, thread count and library version. Whether the hardware under it stayed the same across the two restarts is not known.
- Nothing in `scripts/train_arm.py`, `scripts/plan.py` or `scripts/run_all.py` changed after the first fit started.

## What this does not show

- It is not the test. The eight donors are the donors the forecast was read on, and the rule is decided on F only.
- Eight donors of one study, all healthy, two batches. Nothing here is about disease, another laboratory or another assay.
- TEDDY's pretraining data probably held the RNA of this study. "Not seen" is a statement about the head, the same for every arm.
- The gain is in the variation among a donor's cells. A donor's shift in level against other donors is not read here.
- The 112 proteins both panels share; one head recipe; targets of the two studies on their own scales.
- A quarter of the ceiling is not delivered, and the ceiling is itself a ridge readout: what a head could reach is not bounded by it.
- D3 does not separate donors from cells for a batch that is not in the fit.

## Next

Gate 5, once: the final study is prepared (criteria checked without reading a protein value, donors listed, cells embedded), the recipe is frozen by a commit, `coverage` and `labels` are refitted on all eight donors, the second hypothesis and the sealed declaration are written, and then the heads are run on F. Opening F needs the owner's approval.

## Files

- [`DESIGN.md`](DESIGN.md)
- [`scripts/plan.py`](scripts/plan.py), [`scripts/train_arm.py`](scripts/train_arm.py), [`scripts/run_all.py`](scripts/run_all.py): the heads and their fits
- [`scripts/forecasts.py`](scripts/forecasts.py), [`results/forecasts.json`](results/forecasts.json): gate 0c
- [`scripts/read.py`](scripts/read.py) → [`results/development.json`](results/development.json) (every head, donor and seed) and [`results/ledger.jsonl`](results/ledger.jsonl) (one intervention event per head)
- [`scripts/grid_changed.py`](scripts/grid_changed.py) → [`results/what_is_left.json`](results/what_is_left.json) and `results/grid_changed/<run>/`
- [`scripts/write_readme.py`](scripts/write_readme.py): this file, the figure and [`results/readme_contrasts.json`](results/readme_contrasts.json)
- The heads' weights are not in the repository (42 files of 5.6 MB and three of 15 MB).
