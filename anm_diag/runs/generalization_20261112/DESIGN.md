# Generalization to unseen donors: where is the signal the head does not deliver? — design

Written 12 November 2026, before [`scripts/gate0.py`](scripts/gate0.py) was run on data. Exploratory. No model is
trained and no sealed declaration is run. This is an admission check (gate 0 of
[`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md)) for a loop whose aim is a head that works on donors and diseases it has
not seen.

## Why the head's score is not the target

One look was taken before this design: the frozen head's score (mean per-protein Pearson correlation within a group, 112
proteins measured in both data sets) for every donor.

| cells | score by donor |
|---|---|
| trained on, site 1 (2 donors) | 0.75, 0.76 |
| trained on, site 2 (3 donors) | 0.53, 0.55, 0.57 |
| trained on, site 3 (3 donors) | 0.43, 0.47, 0.55 |
| validation donor, site 1 | 0.72 |
| site 4, not trained on (3 donors) | 0.56, 0.64, 0.66 |
| external PBMC, Hao 2021 (8 donors) | 0.48 to 0.56 |

Cells the head was trained on at site 3 score lower than external donors it never saw. One donor (15078) scores 0.76,
0.57, 0.55 at the three training sites and 0.64 at the site that was held out. The score moves more with the site than
with whether the cells were seen, so a difference in score between seen and unseen donors is not a measure of
generalization. What a group's cells allow has to be measured for each group.

## Question

For a donor the head has not seen, split what is readable of each protein into parts that call for different actions:

1. delivered by the frozen head;
2. delivered by a simple readout fitted on the same cells the head was trained on: a better head, no new data;
3. delivered once other donors of the same study are in the fit: more data of that kind;
4. delivered only by the donor's own cells: specific to the donor;
5. readable from the RNA and not from TEDDY's pooled state: lost before the head;
6. not read by any of these.

## Data (existing files only)

- Bone marrow, NeurIPS 2021 CITE-seq: 90,261 cells, 134 proteins, 9 donors at 4 sites, 12 donor × site groups. The
  head was trained on 8 groups (6 donors, sites 1 to 3). Not trained on: donor 18303 (site 1, used to choose the
  checkpoint), donors 13272 and 19593 (site 4), and donor 15078 at site 4 (the donor is in training at other sites).
- External PBMC, Hao et al. 2021 (GSE164378): 40,000 cells, 8 donors × 5,000, 112 of the 134 proteins. Evaluation only
  for the head. Its RNA may have been in TEDDY's pretraining data.
- The head: the phase-1 head as trained (`scripts/04_train.py`), its MLP and decoder, with the training-median size
  factor, as `scripts/05_eval.py` scores it.

Every one of these donors has been looked at before (site 4 many times, the external donors in E6). Nothing here is
final data.

## What is computed

**Views of a cell.** `RAW`: 512 principal components of log-normalized counts (12,052 genes, components fitted on all
cells, no protein involved). `Z`: TEDDY's pooled state, scaled to unit length as the head reads it (512 numbers). Both
standardized with the head's training cells.

**Readout classes.** Ridge regression on the view, and ridge on the view plus 1,024 random Fourier features of it (a
Gaussian-kernel ridge with a linear term). The penalty is chosen from six values inside each fit set, never on the
scored cells. Of the two classes, the one with the better score inside the fit set is reported.

**Fits**, each for both views:

| | fitted on | stands for |
|---|---|---|
| T | the head's training cells (for a training group: without that donor) | what transfers from the head's own data |
| S | as T, also without the group's site (bone marrow) | T with the site unseen |
| L | all other donors of the same study | what more donors of the same kind would give |
| C | the group's own cells, 5-fold | what the group's cells allow |

No fit ever contains a cell of the donor it is scored on, except C, which is cross-fitted within the group.

**Score.** The repository's score: per-protein Pearson correlation within the group, averaged over the 112 shared
proteins. It ignores a protein's level and scale within a group.

**Terms**, per group, each with an interval from 300 resamples of the group's cells (the same resamples for every
series):

- readout minus head: `Z_T − head`
- same-study donors: `Z_L − Z_T`
- own cells: `Z_C − Z_L`
- RNA minus state: `RAW_C − Z_C` (also for T and L)
- cost of an unseen site: `Z_T − Z_S`

**The donor's shift** (the audit's statistic for a donor contrast). For each group, the shift of every protein's mean
against the other donors of the same study, in units of their spread. The unexplained share *u* = Σ(measured −
predicted)² / Σ measured², over proteins, for the head and for the T and L readouts of both views. The L readout is the
audit's reading: a readout learned on the other donors, asked whether the state explains this donor's shift.

## Readings, written in advance

- A term is **resolved** in a group if its interval lies above 0.02, **absent** if its interval lies below 0.02,
  otherwise **unresolved**. 0.02 is the repository's guard on this score; it is a chosen number.
- A reading is stated for a study only if it holds in every donor of that study the head has not seen:
  - *the head leaves transferable signal undelivered* — readout minus head resolved;
  - *more donors of the same kind are needed* — same-study donors resolved;
  - *part is specific to the donor* — own cells resolved;
  - *part is lost before the head* — RNA minus state resolved.
- A donor's shift is **carried** at a stage if the upper bound of *u* is below 0.10, **not carried** if the lower bound
  is above 0.10 (the audit's margin).

## Limits known in advance

- Donors are few (4 unseen groups in bone marrow, 8 external) and all healthy. The intervals describe each donor's
  cells. Nothing is claimed about donors in general, and nothing about disease.
- The readout classes are not the head's architecture. If a readout scores below the head, that does not show that
  nothing more is available.
- `RAW` is 512 components, not the RNA. "RNA minus state" can miss signal that sits outside those components.
- A donor's shift mixes biology with antibody staining and batch. A shift that is not in the RNA is not a failure of
  the model.
- The two studies normalize protein differently, so levels are compared within a study only.
- The head is in-sample on its training groups; its score there is not comparable with a cross-fitted readout.

## Not done here

No training of teddy_mm, no repair, no audit under a sealed declaration, no choice of final data.
