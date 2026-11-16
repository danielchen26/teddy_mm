# Round 1: a head for donors it has not seen — design (draft)

Draft of 16 November 2026. **Not frozen.** No script of this round exists and no model has been trained. The draft
becomes the design when the owner of the repository approves it; the commit that freezes it will be named here, and
nothing below is changed after that commit except under "Amendments", dated.

This is the first round of the loop of `LOOP_PROTOCOL.md` (version 2) that trains anything for the aim stated in the
root README: a head that works on donors and diseases it has not seen.

## Where the loop stands

1. **Observe** (done, no training). `runs/generalization_20261112` and `runs/scope_20261116`. On the 8 donors of the
   external PBMC study the head's own output leaves *u* = 0.640 of the variation among a donor's cells unexplained
   (mean over 112 proteins of 1 − r² within the donor; range 0.593 to 0.675). A readout fitted on the head's own
   training cells reads no more than the head (gap −0.001 [−0.004, 0.001]). A readout fitted on cells of the other
   donors of that study reads 0.076 [0.068, 0.086] more. The donor's own cells add 0.010.
2. **Hypothesize** (done). `donors_hao_coverage` in `runs/scope_20261116/results/ledger.jsonl`: the shortfall is
   coverage. Forecast for the fall of *u* at the head's own output: coverage arm 0.068 to 0.082 (a ceiling),
   function arm −0.004 to 0.001. Rule: the coverage arm's fall has its interval above 0.034; each control arm's
   below 0.034.
3. **Intervene, test, update**: this round.

## What the sketch lacked

The hypothesis carried a three-line note for the experiment. Checked against the training code and the data on
16 November, it lacked the following. Each row is settled in the section named.

| | found | consequence | settled in |
|---|---|---|---|
| 1 | All 8 PBMC donors were used to write the forecast. The protocol asks that a forecast and its test use different units. | The PBMC study can be fit and development data only. The test needs a study nobody has looked at. | Data and roles; Gate 0 |
| 2 | The PBMC study is registered as evaluation only (`scripts/prepare_external_cite.py`, addendum E6). | Fitting on it ends that role. An addendum has to say so. E6 stays as recorded for the head of that date. | Gate 0 |
| 3 | Its donors P1 to P4 were run in one batch and P5 to P8 in another. | "A donor not seen" has two levels: batch seen, batch not seen. Holding out two donors at random would mix them. | Folds |
| 4 | `scripts/04_train.py` reads one data folder with one split. | A trainer that takes cells of two studies and a list of donors. | Tools |
| 5 | Its loss averages over all 134 proteins. The PBMC study measured 112. | The loss is averaged over the proteins a cell has. | The change |
| 6 | The two studies' targets are on different scales (per-protein mean, PBMC over bone marrow: 0.29 to 4.8, 5th to 95th percentile), and their stored size factors have different definitions (range 0.42 to 1.44 against 0.16 to 14.4). | One definition, written before any fit. | The change |
| 7 | The amount control as noted ("the same number of bone-marrow cells instead") cannot be built: every bone-marrow training cell is already used. | Equal totals: bone-marrow cells are replaced, not added. | Arms |
| 8 | "Another map on the same cells" was not a recipe. | A wider head trained longer. | Arms |
| 9 | No seeds; the seed-to-seed spread of the measured quantity is unknown. | Three seeds per arm; the spread enters every interval. | What is measured |
| 10 | Nothing told a head that learned the PBMC proteins' levels from one that learned how they vary with RNA. | A control with the proteins shuffled among a donor's cells. | Arms |
| 11 | Nothing said whether donors or cells are what is missing. | Arms with one donor or four at the same number of cells. | Arms |
| 12 | `anm_diag scope` takes one output of the system per run. | Several heads on the same cells, compared pairwise with shared resamples. | Tools |
| 13 | A head takes about 4 minutes in the shell this work is run from, whose calls end after 3. | A trainer that stops and resumes. | Tools |
| 14 | No guard, no admission criteria, no statement of what each outcome leads to. | Written below. | Gate 3; Outcomes |

## Data and roles

| data | role | |
|---|---|---|
| bone marrow, training cells (67,405; 6 donors, sites 1 to 3) | fit | as registered |
| bone marrow, donor 18303 | checkpoint choice and the guard of gate 3 | as registered |
| bone marrow, site 4 (3 donors) | the guard of gate 6 | reported, never fitted |
| external PBMC, Hao et al. 2021 (40,000 cells; 8 donors × 5,000; 2 batches × 4 donors; 3 time points per donor; 112 proteins) | fit and development, cross-fitted by donor | **never final** |
| a second PBMC study with patients, **F** | final: opened once, at gate 5 | chosen and locked at gate 0 |
| knockout sets (melanoma, THP-1) | not used | |

Three sets of data, as rule 4 asks. Without F the round stops at gate 4 and stays exploratory.

## Gate 0 (no training)

**0a. Role of the PBMC study.** An addendum under `registration/addenda/` records that from this round its cells may
be fitted on, that E6 is unchanged, and that the study can no longer confirm anything.

**0b. The final study.** Criteria, checked without reading a protein value: open; raw RNA counts per cell; raw protein
counts; a donor and a disease label per cell; at least 20 donors with at least 8 healthy and 8 patients; at least 80
proteins among the 112 the PBMC study shares with the head, by the name-matching rule already used for that study;
never used in this repository. Candidates in order, from `runs/generalization_20261112/results/public_data_census.json`:
Stephenson et al. 2021 (E-MTAB-10026), COMBAT 2022 (Zenodo 5139561), Liu et al. 2021 (GSE161918). From the first
that meets the criteria: 20 healthy donors and 20 patients (all, where fewer), balanced over centres where the study
has several; 1,000 cells per donor, drawn without labels, seed 0. These donors are listed under `roles.final`. The
study's other donors stay unopened as a reserve for one later test. The 40,000 cells are embedded with
`scripts/03_embed_rna.py` as the PBMC study was (about 2 hours on the owner's machine; the PBMC study took 7,610 s).

**0c. Forecasts still to be read**, from existing arrays: the grid on the PBMC donors with the batch as the block (the
ceiling for a donor whose batch was not seen); one readout fitted on bone-marrow and PBMC cells together, without a
study label (what a single map on two scales can deliver).

**0d. Tools**, each with tests on planted data before it sees a measurement:

- `scripts/train_arm.py` in this folder: the loop of `scripts/04_train.py` with lists of cells per study, the masked
  loss, a time budget with resume, and the head's weights written in the format the audits read. With no second
  study it gives the same numbers as `04_train.py` for the same seed and device (a test).
- `anm_diag scope`: named outputs of the system on the same cells; gains between them with shared resamples.
- The four checks of gate 3 as code; the arms file for `anm_diag loop test`; intervention and audit events.
- `scripts/prepare_external_cite.py` extended to the final study.

## The change

One change: cells of the missing kind are added to the head's training. Everything else is the registered recipe:
the model of `teddy_mm/models.py`, AdamW, learning rate 1e-3, 40 epochs, batch 512, checkpoint by the validation
score on bone-marrow donor 18303.

Fixed before any fit:

- **Targets.** Each study's own normalized values as its pack stores them (`adt`): for bone marrow the released
  values, for PBMC the per-cell centred log ratio of `scripts/prepare_external_cite.py`. No rescaling between studies.
- **Loss.** The registered loss, averaged over the proteins measured in the cell. Every cell has weight 1.
- **Size factor.** The recipe's definition applied inside each study: the cell's total over its measured proteins,
  divided by the median of that total over the study's fit cells. The scores below do not depend on it.
- **Total.** Every arm trains on 67,405 cells.
- **Seeds.** 0, 1, 2 for every arm. Subsamples are drawn with the training seed, without labels.
- **What is audited.** The head's MLP and decoder with the training-median size factor, as in every audit so far.

## Arms

| arm | cells | what it is for |
|---|---|---|
| `as_is` | 67,405 bone marrow | the base of every gain; the seed spread |
| `coverage` (own) | 47,405 bone marrow + 20,000 PBMC from the fit donors | the claim |
| `amount` | the same 47,405 bone marrow, nothing added | what removing 20,000 bone-marrow cells does |
| `function` | 67,405 bone marrow; hidden width 1,024, 80 epochs | the rival: the map, not the cells |
| `labels` | as `coverage`, the PBMC proteins shuffled among the cells of each donor | levels learned, relation not |
| `one_donor` | 62,405 bone marrow + 5,000 PBMC from one fit donor | donors or cells |
| `four_donors` | 62,405 bone marrow + 5,000 PBMC from four fit donors (1,250 each) | donors or cells |

The gain of an arm is the fall of *u* from `as_is` with the same seed, on the same cells.

## Folds (development)

Four fit sets of four PBMC donors. Every donor is scored only by heads that did not see it, once with its batch seen
and once with its batch not seen.

| fit set | fit donors | scored donors | batch of the scored donors |
|---|---|---|---|
| M1 | P1 P2 P5 P6 | P3 P4 P7 P8 | seen |
| M2 | P3 P4 P7 P8 | P1 P2 P5 P6 | seen |
| B1 | P1 P2 P3 P4 | P5 P6 P7 P8 | not seen |
| B2 | P5 P6 P7 P8 | P1 P2 P3 P4 | not seen |

`coverage` on all four; `labels`, `one_donor` (the first fit donor of each batch) and `four_donors` on M1 and M2.
45 heads: 3 + 3 + 3 + 12 + 6 + 12 + 6.

## Gate 3: admission of an arm

1. Converged: over the last 5 epochs the validation score changes by less than 0.005.
2. Reached its own objective: on its own PBMC fit cells `coverage` scores at least 0.035 above `as_is` (mean over
   proteins of the correlation within a donor; half of the smallest gain gate 0 read there). `function` scores no
   lower than `as_is` on donor 18303. The other arms: criterion 1 only.
3. No scored cell in the fit, checked by donor and barcode. No cell of F in any file a fit reads.
4. The guard on donor 18303: the mean over 134 proteins of the correlation between measured and predicted is no more
   than 0.02 below `as_is` of the same seed.

An arm that fails was not tested. It did not fail the test.

## What is measured

- **Primary.** *u* of the `variation` target in the `frozen` scope: for a donor, 1 minus the mean over proteins of
  r² between the head's output and the measurement among that donor's cells. Margin on this scale: the rule's 0.034;
  the tool's 0.10 is not used (protocol, "After the first use", 3).
- **Gain and interval.** Per scored donor, the mean over seeds of *u*(`as_is`) − *u*(arm). Over donors, the mean;
  the interval from 1,000 resamples of donors with seeds resampled inside each arm. With 8 donors it is marked
  approximate.
- **In plain terms.** The head's ordinary score on the same cells: the mean over proteins of the correlation within a
  donor, and the share of the gate-0 ceiling delivered.
- **Guard (gate 6).** Bone marrow, site 4: the mean score and the worst single protein, and separately the 22
  proteins the PBMC study does not measure.

## Gate 4: what development is read for (exploratory)

| | question | read from |
|---|---|---|
| D1 | Does the arm deliver on donors not seen, batch seen, and how much of the ceiling (0.075)? | `coverage`, M1 and M2 |
| D2 | And when the batch was not seen? | `coverage`, B1 and B2 against M1 and M2, paired by donor; forecast 0c |
| D3 | Donors or cells? | `one_donor`, `four_donors`, `coverage` |
| D4 | Levels or relation? | `labels` against `coverage` |
| D5 | What it costs | the guard, every arm |
| D6 | What is left | the grid on the changed head, scored donors only |

**Go to gate 5 only if** `coverage` is admitted in every seed, its gain on B1 and B2 has its interval above 0, and
the guard holds. Otherwise the round ends here with its reading, and F stays unopened.

Development carries this round and one revision of the recipe. After that it is spent.

## Gate 5: the test, once

The recipe is frozen by a commit. Then, in this order:

1. New fits, 3 seeds: `coverage` with 20,000 PBMC cells from all 8 donors (2,500 each); `labels` likewise. `as_is`,
   `amount` and `function` are the heads already fitted.
2. A second hypothesis is entered in the ledger before F is opened: `coverage_not_levels`, own arm `coverage`,
   control `labels`, the same bounds, its forecasts read from development (other units than F).
3. A sealed declaration names the hypotheses, the arms, the donors of F and the proteins (those F measures among the
   112).
4. The heads are run on F. `anm_diag loop test` reads `donors_hao_coverage` as it stands in the ledger: the event is
   "a donor the head was not trained on, external PBMC"; the units are new. The outcome has to hold in every seed.
5. Also read, declared here: the gain among healthy donors, among patients, and their difference.
6. Only then the grid on F (donor and disease as contrast events, checked by `events.check_contrast`), for round 2.

## Outcomes

| outcome | what it means | next |
|---|---|---|
| confirmed, guard holds, patients as healthy | the shortfall was coverage of the tissue, and it carries across study and disease | the head is released with these numbers (a fit on all bone-marrow and PBMC cells, reported beside the tested one); round 2 from the grid on F |
| confirmed among healthy, smaller among patients | coverage of the tissue; disease is a further gap | round 2: patients in the fit, tested on the reserve of F |
| own arm below the bound, development had delivered | what was learned belongs to the study | round 2: a second fit study, or study-specific scales in the loss |
| a control arm above its bound | the reading was wrong | back to the grid |
| unresolved | the intervals do not decide | the reserve of F, once |
| guard broken | the change costs the head its own task | not released; the cost is the finding |

## Limits

- TEDDY's pretraining data probably held the RNA of both PBMC studies (addendum E6). "Not seen" is a statement about
  the head. It is the same for every arm.
- PBMC and bone marrow only; the proteins both panels share; one head recipe.
- Eight development donors: every development interval is approximate.
- The normalized targets of two studies are not on one scale. If the plain recipe cannot hold both, the round says so
  and the remedy is a change of its own.

## Work

| step | where | time | needs the owner |
|---|---|---|---|
| tools and tests (0d) | anywhere; no data | about a day | no |
| forecasts (0c) | existing arrays | 1 to 2 hours | no |
| addendum (0a) | repository | minutes | approval |
| final study: check, download, embed (0b) | the owner's machine | about 3 hours, mostly unattended | one command |
| 45 heads, gate 3, development | the owner's machine | about 4 hours through the restricted shell, about 1 hour run directly | the machine awake |
| gate 5 | the owner's machine | under 1 hour | approval to open F |

## Amendments

None.
