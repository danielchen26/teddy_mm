# The loop on the teddy_mm head, re-read under the protocol — 12 November 2026

Post hoc and exploratory. No fit and no audit was made for this folder. It reads the records and bootstrap draws of the
three rounds ([`loop_20261112`](../loop_20261112/), [`loop_guard_20261112`](../loop_guard_20261112/),
[`repair_20261112`](../repair_20261112/)) with [`anm_diag/loop.py`](../../anm_diag/loop.py), under
[`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md). The outcomes those rounds recorded under their own rules stand as recorded.

Why it was made: the rounds decided by comparing a point value with 0.10. The audits of all heads of a data set use the
same cells, folds and resamples (the draws of the response are identical, checked), so the change between two stages and
the difference between two heads both have intervals that can be read from what was already computed.

The span is Z → head output: from the stage before the located step to what the head delivers.

## In short

- **The located step is confirmed, on both cell types.** Over Z → output the head as trained adds 0.39 (melanoma) and
  0.42 (THP-1) to *u* for PD-L1. Repairing the MLP, the step the audit names, removes 84% [74%, 98%] and 91% [85%, 99%]
  of that. Repairing the decoder with the same cells, loss and budget removes 5% [-4%, 13%] and 15% [8%, 21%].
- **No head is shown to be repaired in the paper's sense.** No changed head has a supported end stage. On THP-1 the changed
  heads add no loss over the span, by their intervals ("relative" closure). On melanoma the loss is no longer resolved, and
  its absence is not shown either.
- **Round 3 is slightly worse than round 2 on melanoma**, in a direct paired comparison, and no different on THP-1.
- **HLA-E and MHC-I should not have been targets.** Their response is too small against its own interval to follow.

## 1. Admission (gate 0), on the audit cells of the loop

| cells | readout | response over the half-width of its interval | reference Z here | reference Z in the sealed run | head as trained, Z → output | half-width | diagnostic loop | closure attainable | clusters |
|---|---|---|---|---|---|---|---|---|---|
| melanoma | MHC-I | 2.9 | unresolved | unresolved | -0.14 [-0.74, +0.31] unresolved | 0.53 | not ready | not attainable | 5 |
| melanoma | PD-L1 | 10.1 | unresolved | support | +0.39 [+0.34, +0.44] loss | 0.05 | ready | relative | 5 |
| melanoma | HLA-E | 1.9 | unresolved | unresolved | +1.09 [-0.47, +1.43] unresolved | 0.95 | not ready | not attainable | 5 |
| THP-1 | PD-L1 | 3.6 | unresolved | unresolved | +0.42 [+0.32, +0.49] loss | 0.09 | ready | marginal | 5 |
| THP-1 | PD-L2 | 5.3 | unresolved | unresolved | +0.48 [+0.37, +0.55] loss | 0.09 | ready | marginal | 5 |

Reference verdicts are two-sided. The sealed runs used every control cell (5,171 and 7,478 cells); the loop's audits keep half of
the controls for the fits. With 5 treated clusters every interval here is approximate, so the grade stays exploratory
whatever else holds. "Relative": a head that adds almost nothing over the span could be recognized. "Marginal": only just.

## 2. The diagnostic loop: is the MLP the step that matters?

Baseline: the head as trained. Located change: the MLP repaired (round 1). Control change: the decoder repaired with the
same cells, loss and budget. Each fit holds one target gene out, and each knockout cell is audited through the head that
never saw its gene.

| cells | readout | head as trained adds | MLP repaired: share removed | decoder repaired: share removed | MLP minus decoder, in *u* | outcome | selective |
|---|---|---|---|---|---|---|---|
| melanoma | PD-L1 | 0.39 [+0.34, +0.44] | 84% [74%, 98%] | 5% [-4%, 13%] | +0.31 [+0.26, +0.37] | located step confirmed | yes |
| THP-1 | PD-L1 | 0.42 [+0.32, +0.49] | 91% [85%, 99%] | 15% [8%, 21%] | +0.32 [+0.25, +0.38] | located step confirmed | yes |
| THP-1 | PD-L2 | 0.48 [+0.37, +0.55] | 104% [66%, 117%] | 25% [16%, 37%] | +0.38 [+0.12, +0.45] | located step confirmed | no |

"Selective": whatever the decoder repair removes lies below the step margin. PD-L2 has no output in the head and was in no
loss; the MLP repair made for PD-L1 removes its loss as well.

**Where the loss shows and where it can be removed are not the same place** (melanoma, PD-L1). As trained the head adds a
resolved loss at both of its steps: +0.16 [+0.11, +0.21] from Z to the latent and +0.23 [+0.19, +0.27] from the latent to the output. With
the MLP repaired and the decoder untouched, the second step adds +0.02 [-0.04, +0.08]. With the decoder repaired and the MLP
untouched, it still adds +0.21 [+0.17, +0.26]. The loss at the decoder's step goes when the MLP is repaired and stays when
the decoder is.

One held-out gene at a time (PD-L1; the head fitted without that gene, audited on its knockouts):

| cells | gene | clusters | gate 0 | MLP repaired: share removed | decoder repaired: share removed | outcome |
|---|---|---|---|---|---|---|
| melanoma | IFNGR1 | 22 | ready | 133% [61%, 181%] | 7% [-40%, 28%] | located step confirmed |
| melanoma | IFNGR2 | 28 | ready | 78% [68%, 190%] | 2% [-39%, 14%] | located step confirmed |
| melanoma | JAK1 | 23 | ready | 109% [80%, 156%] | 3% [-21%, 27%] | located step confirmed |
| melanoma | JAK2 | 24 | ready | 120% [90%, 182%] | -1% [-31%, 12%] | located step confirmed |
| melanoma | STAT1 | 17 | ready | 73% [62%, 120%] | 14% [-0%, 33%] | located step confirmed |
| THP-1 | IFNGR1 | 4 | ready | 76% [72%, 138%] | 22% [15%, 60%] | located step confirmed |
| THP-1 | IFNGR2 | 4 | ready | 98% [86%, 115%] | 25% [18%, 32%] | located step confirmed |
| THP-1 | JAK2 | 4 | ready | 88% [76%, 98%] | 23% [17%, 29%] | located step confirmed |
| THP-1 | STAT1 | 4 | ready | 95% [87%, 168%] | 21% [-7%, 45%] | located step confirmed |
| THP-1 | IRF1 | 4 | not ready: response over half-width 1.2 | — | — | not admitted |

9 of the 9 admitted genes read "located step confirmed". With about 200 knockout cells a share has a wide interval,
and shares above 100% mean the repaired head's *u* at the output came out below its *u* at Z.

## 3. The repair loop: is a changed head free of the loss?

Control: the head as trained, on the same cells. Grades as in the protocol. Bone marrow: the head's own task, mean
correlation over 134 proteins on 5,000 test cells, lowest to highest of the round's heads; 0.602 as trained, and the
guard allows a mean fall of 0.02.

| cells | readout | head | added over Z → output | end stage, two-sided | closure | by step: Z → latent, latent → output | bone marrow | guard |
|---|---|---|---|---|---|---|---|---|
| melanoma | PD-L1 | control: decoder repaired | +0.371 [+0.311, +0.429] | reject | **not closed** | not closed, not closed | 0.579 to 0.595 | kept |
| melanoma | PD-L1 | round 1: MLP repaired | +0.063 [+0.009, +0.110] | unresolved | **loss no longer resolved** | loss no longer resolved, relative | 0.560 to 0.567 | not kept |
| melanoma | PD-L1 | round 2: MLP repaired, bone-marrow outputs held | +0.064 [+0.013, +0.112] | unresolved | **loss no longer resolved** | relative, relative | 0.603 to 0.604 | kept |
| melanoma | PD-L1 | round 3: one fit for both cell types | +0.109 [+0.057, +0.157] | unresolved | **loss no longer resolved** | loss no longer resolved, loss no longer resolved | 0.596 to 0.598 | kept |
| THP-1 | PD-L1 | control: decoder repaired | +0.355 [+0.273, +0.421] | reject | **not closed** | not closed, not testable | 0.364 to 0.491 | not kept |
| THP-1 | PD-L2 | control: decoder repaired | +0.361 [+0.241, +0.429] | reject | **not closed** | not closed, loss no longer resolved | 0.364 to 0.491 | not kept |
| THP-1 | PD-L1 | round 1: MLP repaired | +0.036 [+0.003, +0.066] | unresolved | **relative** | loss no longer resolved, not testable | 0.502 to 0.552 | not kept |
| THP-1 | PD-L2 | round 1: MLP repaired | -0.019 [-0.078, +0.137] | unresolved | **loss no longer resolved** | loss no longer resolved, loss no longer resolved | 0.502 to 0.552 | not kept |
| THP-1 | PD-L1 | round 2: MLP repaired, bone-marrow outputs held | +0.060 [+0.029, +0.096] | unresolved | **relative** | loss no longer resolved, not testable | 0.594 to 0.595 | kept |
| THP-1 | PD-L2 | round 2: MLP repaired, bone-marrow outputs held | -0.017 [-0.085, +0.187] | unresolved | **loss no longer resolved** | loss no longer resolved, loss no longer resolved | 0.594 to 0.595 | kept |
| THP-1 | PD-L1 | round 3: one fit for both cell types | +0.047 [-0.005, +0.098] | unresolved | **relative** | loss no longer resolved, not testable | 0.596 to 0.598 | kept |
| THP-1 | PD-L2 | round 3: one fit for both cell types | +0.079 [-0.060, +0.220] | unresolved | **loss no longer resolved** | loss no longer resolved, loss no longer resolved | 0.596 to 0.598 | kept |

- **Absolute closure: nowhere.** The end stage of every changed head is unresolved. On these audit cells the reference
  stage Z is itself unresolved, so no later stage can be supported.
- **THP-1, PD-L1: relative closure in all three rounds.** The upper bound of the added loss is below 0.10 each time.
- **Melanoma, PD-L1: the loss is no longer resolved, in all three rounds.** Rounds 1 and 2 miss "relative" by their upper
  bounds (0.110 and 0.112 against 0.10). Read one step at a time, round 2 adds no loss at either step.
- **The decoder repair closes nothing**, on either cell type.
- **Round 1 does not keep the guard.** Rounds 2 and 3 do.

## 4. Round 3 against round 2, on the same cells

Both are repairs of the MLP with the bone-marrow outputs held. Round 2 fits one cell type and PD-L1. Round 3 fits both
cell types and every measured protein with a head output. PD-L1, Z → output, paired:

| cells | round 2 adds | round 3 adds | round 3 minus round 2 | reading |
|---|---|---|---|---|
| melanoma | +0.064 [+0.013, +0.112] | +0.109 [+0.057, +0.157] | +0.044 [+0.002, +0.090] | round 3 adds more |
| THP-1 | +0.060 [+0.029, +0.096] | +0.047 [-0.005, +0.098] | -0.014 [-0.053, +0.022] | no difference resolved |

On melanoma the joint head adds more than the per-cell-type head, by an interval that just excludes zero. The rule of
round 3 called this "regressed" from a point value 0.009 over the margin; the direct comparison supports the direction
of that reading and gives its size.

Round 3 on its three splits of the control cells (PD-L1, Z → output; no head as trained was audited on splits 1 and 2):

| cells | split 0 | split 1 | split 2 |
|---|---|---|---|
| melanoma | +0.109 [+0.057, +0.157] unresolved | +0.111 [+0.051, +0.159] unresolved | +0.140 [+0.061, +0.185] unresolved |
| THP-1 | +0.047 [-0.005, +0.098] no added loss | +0.020 [-0.015, +0.057] no added loss | +0.009 [-0.017, +0.043] no added loss |

## 5. A head repaired on the other cell type (round 1)

Audited on every cell of the data set, since the head never saw any of them. PD-L1, Z → output:

| audited on | head repaired on | added | reading | reference Z | end stage |
|---|---|---|---|---|---|
| melanoma (5,171 cells) | THP-1, MLP | +0.39 [+0.33, +0.43] | loss | support | reject |
| THP-1 (7,478 cells) | melanoma, MLP | +0.26 [+0.17, +0.35] | loss | unresolved | reject |

A repair does not carry to the other cell type. With every control cell in the audit, the reference stage on melanoma is
supported: the loop's own audits lose that by keeping half of the controls for the fits.

## What this changes in how the rounds are read

| the rounds' reports say | under the protocol |
|---|---|
| rounds 1 and 2: "closed" on both cell types | the located step is confirmed on both; the head is shown free of added loss on THP-1 only, and only relative to an unresolved reference |
| round 3: "regressed" on melanoma by 0.009 | on melanoma the joint head adds 0.044 [+0.002, +0.090] more than the round-2 head; on THP-1 no difference |
| round 3: HLA-E "still located in the head" | HLA-E does not pass gate 0 on these cells; it was not a target that could be followed |
| "no loss in the head" as a pass | three different things: no added loss shown (THP-1), loss no longer resolved (melanoma), and never a supported head |

## Limits

- Post hoc: the protocol and this reading were written after every result here was known.
- 5 treated clusters in the pooled audits and 4 in the THP-1 single-gene audits: intervals are approximate.
- The same audit cells served all three rounds. Nothing here is a final audit.
- Changes are fine-tunes of a trained head with one training seed, compared with the head as trained.
- Changes in *u* are taken in |*u*|. For PD-L1 and PD-L2 on these cells that is the same as the signed change.

## Files

`scripts/read_loop.py` reads the three rounds and writes `results/loop_reading.json` (every reading, with intervals) and
`results/ledger.jsonl` (the changes, audits and decisions as typed events). `scripts/write_readme.py` writes this file.
The bootstrap draws (`anm_bootstrap.npz` in each audit folder) are not kept in git; re-running an audit writes them again.
