# ANM diagnosis of TEDDY: where does a perturbation response stop being recoverable?

Code and records: [`anm_diag/`](../anm_diag/) (how-to in its [README](../anm_diag/README.md); every number below is in
`anm_diag/runs/*/anm_record.json`). Method: Active Neural Matter (arXiv 2610.14549), §2.2–§2.3 and Table 2.

> [!IMPORTANT]
> **Status (9 Nov 2026).** The sealed tests are complete. TEDDY's own stages (RAW → TOKENS → Z) are never rejected.
> After the head-input fix, the sealed rerun puts a loss between TEDDY's embedding and the teddy_mm head's latent layer
> (Z → HEAD_LAT): absolute for Frangieh PD-L1, relative for the other readouts. Two exploratory follow-ups, each sealed
> and run once, have also run: splitting the head places the loss in its later layers and decoder (the normalization
> and the first layer lose nothing), and a forward test of typed source events shows how the head responds to the
> knockout direction. A rule written before the head split selects the next intervention (retrain the head with
> IFN-γ-varied cells); that step has not run. Details: [Corrections](#corrections-and-deviations).

## The question

A foundation model is useful downstream only if its representation keeps the information that downstream questions
need. Accuracy on a benchmark does not show that: a model can score well on average and still discard the specific
signal a question depends on. ANM asks a sharper question of each stage *z* on a model's computation path:

> When a controlled perturbation (a CRISPR knockout) changes a **measured** readout (a protein, a gene set), can a
> bounded readout class recover that change from *z* alone?

If the change is recoverable from the input but not from a later stage, the information was merged somewhere in
between, and the test names the step.

## Design

| ANM role | here |
|---|---|
| frozen system | TEDDY-G 70M (12 layers, 8 heads, d = 512), followed by the teddy_mm phase-1 protein head |
| source perturbation *A* | CRISPR knockout / CRISPRi of a target gene vs non-targeting controls; folds and bootstrap by target |
| readout *O* | a **measurement** the model never sees (surface protein by CITE-seq), or a declared gene-set score |
| candidate states (chain) | RAW (normalized counts) → TOKENS (TEDDY's rank tokens) → Z (gene-mean of the last layer, 512-D) → HEAD_LAT (the head's latent layer) → HEAD_OUT (its 134 protein predictions) |
| probes (off the chain) | SET (the same genes without rank order), single-gene mRNA (e.g. CD274 for PD-L1), positive control (the readout itself) |
| readout classes | ridge and gradient boosting, fixed hyperparameters |
| tolerance | *u* margin 0.10, step margin 0.10, α = 0.05 |

**The statistic.** For each stage, cross-fit (by target) a model of the readout from *z* and of the treatment from *z*,
and take the slope of the readout residual on the treatment residual:

θ(z) = Σ (O − m̂(z)) (A − ê(z)) / Σ (A − ê(z))²,  *u* = θ(z) / Δ,

where Δ is the readout's raw response (treated minus control). *Intuition:* θ(z) is the part of the knockout's effect
on the readout that is **left over** after everything *z* knows has been used. If *z* retains the response, nothing is
left over (*u* ≈ 0); if *z* has lost it, the whole response is left over (*u* ≈ 1). *u* is a share, so a 0.10 margin
reads "this stage explains at least 90% of the response". A clearly negative *u* means the stage predicts a larger
effect than the one measured.

**Verdicts (sealed rule, one-sided).** Support: upper bound of *u* < 0.10 for some declared readout class. Reject:
lower bound > 0.10 for every class. Otherwise unresolved. A reject is bounded by the readout classes (paper eq. 5): it
says ridge and gradient boosting cannot recover the response from *z*, not that no function could. A step between two
stages adds a resolved loss when the lower bound of the change in *u* exceeds the step margin.

**Root cause** (first match wins): `no_effect` (no response to explain) · `control_failed` (the positive control is not
recovered, so the test itself is not trusted) · `relative_loss_at` (the reference stage is not within tolerance, but one
step adds a resolved loss) · `reference_unresolved` · `loss_at` (the reference is supported and a later stage is
rejected: the step between them merged the distinction) · `unresolved` · `no_deficit`.

**Guardrails before any statistic** (`python -m anm_diag check`): the readout must be a measurement (R1) that does not
flow into the candidate states (R2); every chain stage must be computed from the previous one by the system's own
operator (C1, C2); no analyst compression such as PCA (C3); no input masking or state replacement (C4). Each rule is
there because an earlier teddy_mm analysis made that mistake. **Sealing** (`seal`) freezes the design with the sha256 of
every input file and of the package source; `run` refuses if any of them changed.

## Results

### KS-2: TEDDY on a CRISPRi screen (exploratory)

Replogle et al. 2022 K562 essential-gene CRISPRi; TEDDY-G 70M gene-mean embeddings. The readouts are gene-set scores
from the same counts that feed TEDDY (shared measurement, declared). The design was fixed after looking at these
development cells, so these results are **exploratory**.

**Mitoribosome knockdowns** (6 targets, 240 treated of 1,956 cells; reread with the union rule, see
[Corrections](#corrections-and-deviations)):

| readout | Δ | RAW | TOKENS | Z | SET (no order) | root cause |
|---|---|---|---|---|---|---|
| MT genes | +0.12 [+0.05, +0.19] | −0.01 [−0.08, +0.08] support | +0.01 [−0.14, +0.09] support | **+0.47 [+0.29, +0.64] reject** | +0.78 [+0.47, +0.90] reject | `loss_at` TOKENS → Z |
| PSM genes | −0.07 [−0.09, −0.06] | +0.19 [+0.08, +0.29] unresolved | +0.31 reject | +0.43 reject | +0.46 reject | `reference_unresolved` |

Forward AUC (does the knockdown reach the stage): RAW 0.789, TOKENS 0.717, Z 0.701, SET 0.670.

*Reading.* The MT response is recoverable from the input but not from the same genes without their order (SET) and
not from the pooled embedding (Z), so averaging token states over genes loses it. Under the paper's two-sided rule
(below) the token stage becomes unresolved (its interval reaches −0.14), so the loss is located between RAW and Z
rather than exactly at pooling. For PSM even the input view does not pass, so nothing is attributed to TEDDY.

**Ribosome knockdowns** (12 targets, 480 treated of 2,196 cells): RP and MRP readouts are supported at every stage
(`no_deficit`); the MT readout is supported at RAW and TOKENS and unresolved at Z (+0.07 [+0.03, +0.11]). Forward AUC
RAW 0.855, TOKENS 0.839, Z 0.817, SET 0.844.

### teddy_mm: TEDDY and its protein head on Perturb-CITE-seq and ECCITE-seq (sealed, run once)

Protein readouts measured by CITE-seq, never seen by TEDDY or the head. Cells are labelled from their detected guides;
genes mapped to TEDDY's vocabulary with the 10x GRCh38-3.0.0 table (the files' own Ensembl ids drop the MHC genes).

- **Frangieh 2021** (melanoma, IFN-γ condition): knockouts IFNGR1, IFNGR2, JAK1, JAK2, STAT1 (1,171 treated cells)
  vs 4,000 non-targeting controls; readouts PD-L1, MHC-I (HLA-ABC antibody), HLA-E, CLR across the 20-protein panel.
- **Papalexi 2021** (THP-1, IFN-γ stimulated): knockouts IFNGR1, IFNGR2, JAK2, STAT1, IRF1 (7,478 cells); readouts
  PD-L1, PD-L2, log1p (4-protein panel).

Records: `anm_diag/runs/frangieh2021_headfix/`, `anm_diag/runs/papalexi2021_headfix/` (sealed declarations
`declarations/*_headfix.sealed.yaml`, code hash `cb8fca78dc47…`, the package in this repository).

**Every stage** (*u*, 95% interval, verdict under the sealed rule):

| data | readout | Δ | RAW | TOKENS | Z | HEAD_LAT | HEAD_OUT | Z → HEAD_LAT (change in *u*) | root cause |
|---|---|---|---|---|---|---|---|---|---|
| Frangieh | PD-L1 | −0.85 [−0.93, −0.76] | −0.01 [−0.06, +0.04] **support** | +0.01 [−0.03, +0.05] **support** | +0.04 [−0.01, +0.09] **support** | +0.20 [+0.13, +0.26] **reject** | +0.46 [+0.41, +0.50] **reject** | [+0.12, +0.20] | `loss_at` Z → HEAD_LAT |
| Frangieh | MHC-I | −0.14 [−0.18, −0.09] | −0.27 unresolved | −0.29 support | −0.37 support | +0.01 [−0.34, +0.27] unresolved | +0.19 [−0.10, +0.38] unresolved | [+0.19, +0.65] | `relative_loss_at` Z → HEAD_LAT |
| Frangieh | HLA-E | −0.08 [−0.11, −0.05] | −0.27 unresolved | −0.16 unresolved | +0.18 unresolved | +0.94 [+0.51, +1.27] reject | +1.36 [+1.09, +1.88] reject | [+0.44, +1.21] | `relative_loss_at` Z → HEAD_LAT |
| Papalexi | PD-L1 | −0.69 [−0.86, −0.47] | +0.10 [+0.00, +0.23] unresolved | +0.12 unresolved | +0.15 unresolved | +0.35 [+0.30, +0.45] reject | +0.47 [+0.43, +0.58] reject | [+0.12, +0.28] | `relative_loss_at` Z → HEAD_LAT |
| Papalexi | PD-L2 | −0.29 [−0.34, −0.23] | +0.08 [−0.23, +0.34] unresolved | +0.08 unresolved | +0.13 unresolved | +0.40 [+0.23, +0.62] reject | +0.55 [+0.39, +0.79] reject | [+0.22, +0.35] | `relative_loss_at` Z → HEAD_LAT |

Forward AUC: Frangieh RAW 0.895, TOKENS 0.889, Z 0.884, **HEAD_LAT 0.843, HEAD_OUT 0.811**, SET 0.886; Papalexi RAW
0.767, TOKENS 0.752, Z 0.757, **HEAD_LAT 0.636, HEAD_OUT 0.577**, SET 0.748.

**Single-gene mRNA probes** (the protein's own transcript as the only input): Frangieh PD-L1 from CD274 mRNA
+0.93 [+0.92, +0.95] reject; Papalexi PD-L1 from CD274 +0.69 [+0.63, +0.78] reject; PD-L2 from PDCD1LG2
+0.80 [+0.77, +0.88] reject.

*Reading.*

1. **No TEDDY stage is rejected for any of the five protein readouts.** Where the data resolve the question
   (Frangieh PD-L1), TEDDY's pooled embedding Z retains at least 90% of the protein response to the IFN-γ-pathway
   knockouts.
2. **The protein response is in the transcriptome but not in the protein's own transcript.** CD274 mRNA explains almost
   none of the PD-L1 protein change (*u* +0.93), the whole RNA profile explains all of it (*u* −0.01). The knockouts act
   through the pathway program, which TEDDY's whole-cell embedding captures and a single-gene readout misses.
3. **Rank order is not needed here** (SET supported for PD-L1), unlike the KS-2 MT readout. Whether pooling loses
   information depends on the readout.
4. **The loss is in the teddy_mm head.** For Frangieh PD-L1 the embedding is supported and the head's latent layer is
   rejected, so the distinction merges between Z and HEAD_LAT (`loss_at`). For the other readouts the head stages are
   rejected (MHC-I excepted) and the same step adds a resolved loss, but the input view itself is not within tolerance
   (five targets per dataset give wide intervals), so the claim is relative. The forward AUC falls at the head too:
   the knockout becomes harder to detect after the head than in TEDDY's embedding.
5. **Scope.** The head was trained on bone-marrow CITE-seq; here it is applied unchanged to melanoma and THP-1 cells.
   The finding is "this head, applied to these cells, loses these responses", not a statement about the head on its
   training distribution. In the paper's terms (Table 2, last row) the contrast the task needs is preserved in the
   input to the head and missing downstream; the next experiment is to change or retrain the head and rerun the same
   declaration.

### Two-sided reread (post hoc, exploratory)

The sealed rule is one-sided: it tests only for loss, so a stage whose *u* is clearly negative (it predicts a larger
effect than the one measured) still counts as support. The paper's bounded support is two-sided (the whole interval
inside ±margin). Rereading the same records with the two-sided rule (a stage is supported if some class's interval
lies within ±0.10, rejected if every class's interval lies outside; a step is a loss if |*u*| grows by more than the
step margin) changes three things:

- **Frangieh MHC-I:** TOKENS, Z and SET go from support to unresolved; the root cause becomes `reference_unresolved`.
  Its Z sits on the over-explaining side (*u* ≈ −0.37), so moving to *u* ≈ 0 at HEAD_LAT is not a loss.
- **Frangieh HLA-E:** still `relative_loss_at`, but the resolved step moves to HEAD_LAT → HEAD_OUT.
- **KS-2 (already exploratory):** the MT readout's token stage becomes unresolved, so `loss_at` widens from TOKENS → Z
  to RAW → Z; the ribosome screen's MRP readout goes from `no_deficit` to `reference_unresolved`.

Frangieh PD-L1, Papalexi PD-L1 and PD-L2 are unchanged. The reread fits nothing new; it is labelled exploratory
because the rule was changed after the sealed results were seen. (The reread uses anm_diag 0.2.0, which arrives with
the head-split follow-up.)

### Where in the head? Head split (exploratory, sealed, run once)

The head computes HEAD_LAT = L3(SiLU(L2(SiLU(L1(z / (‖z‖ + 1e-6)))))). The head-split declarations insert each of
these operations as its own stage (Z_UNIT, L1, H1, L2, H2) and use the paper's two-sided rule; every stage shared with
the head-fix rerun matches it digit for digit. Records: `anm_diag/runs/*_headsplit/`.

Increase in |*u*| between neighbouring stages (95% interval; bold = lower bound above 0.10, a resolved step):

| data | readout | Z → Z_UNIT | Z_UNIT → L1 | L1 → H1 | H1 → L2 | L2 → H2 | H2 → HEAD_LAT | HEAD_LAT → HEAD_OUT | root cause |
|---|---|---|---|---|---|---|---|---|---|
| Frangieh | PD-L1 | [−0.005, +0.004] | [−0.031, +0.025] | [−0.014, +0.019] | [−0.020, +0.026] | [−0.034, +0.054] | **[+0.108, +0.187]** | **[+0.223, +0.294]** | `loss_at` H1 → HEAD_LAT (via unresolved L2, H2) |
| Frangieh | MHC-I | [−0.034, +0.011] | [−0.258, +0.152] | [+0.000, +0.000] | [−0.370, +0.237] | [−0.099, +0.082] | [−0.193, +0.178] | [−0.268, +0.229] | `reference_unresolved` |
| Frangieh | HLA-E | [−0.140, +0.231] | [−0.367, +0.353] | [−0.030, +0.031] | [−0.267, +0.282] | **[+0.162, +0.647]** | [−0.057, +0.275] | **[+0.163, +0.993]** | `relative_loss_at` L2 → H2 |
| Papalexi | PD-L1 | [+0.000, +0.006] | [−0.014, +0.032] | [−0.008, +0.030] | [+0.016, +0.044] | [+0.010, +0.038] | [+0.088, +0.166] | **[+0.111, +0.150]** | `relative_loss_at` HEAD_LAT → HEAD_OUT |
| Papalexi | PD-L2 | [−0.011, +0.008] | [−0.052, +0.055] | [−0.024, +0.028] | [−0.001, +0.069] | [+0.044, +0.097] | [+0.099, +0.157] | **[+0.126, +0.182]** | `relative_loss_at` HEAD_LAT → HEAD_OUT |

*Reading.* The resolved losses sit in the later half of the head and in the decoder; none in the first half. Dividing
by the length changes |*u*| by at most 0.04 for PD-L1, PD-L2 and MHC-I; no readout loses anything at the first layer or
its SiLU. Frangieh PD-L1 keeps bounded support through H1, and its resolved step is the third linear layer; HLA-E's is
the second SiLU; for Papalexi each later step adds a little and the resolved step is the decoder.

*Why, from the operators* (`head_operator_stats.json`, weights and TEDDY states only, no readouts): ‖z‖ varies by only
4–5% across cells (coefficient of variation 0.041 and 0.050), so the normalization has almost nothing to discard; no
input to the first SiLU falls on its folding branch (x < −1.2785), so on these cells it is one-to-one; 16% (Frangieh)
and 28% (Papalexi) of the second SiLU's activations do, the first place where information can actually be destroyed.
The third linear layer is full rank (weights in `runs/teddy_mm_head/head_weight_svd.json`: 481 of 512 singular values
below 10% of the largest), so a loss measured there is a direction squeezed below what the declared readout classes
read, not information gone. The decoder ends in 134 outputs and necessarily drops all other directions.

This is exploratory (designed after the head-fix results) and changes no sealed verdict; it says where the loss sits,
not why.

### Forward: typed source events (exploratory, sealed, run once)

The forward half gives the frozen pipeline typed in silico source events (site = gene, channel = count or gene
embedding, sign and amplitude; plus a channel that moves Z itself) and follows the state displacement layer by layer
and the head's predicted protein (log mean). Records: `anm_diag/runs/forward_typed_events/` (100 confirmatory cells per
dataset; detection floors set on 20 separate development cells). The readout is model output, so this half describes
how the model responds, not whether a stage is sufficient.

- **Checks pass.** A routing control (counts of expressed genes outside TEDDY's vocabulary) leaves tokens, states and
  outputs exactly unchanged; the baseline Z recomputed from the weights matches the stored TEDDY states to within 1e-6
  (9.5e-7 and 4.8e-7).
- **The knockout direction.** Moving Z along the measured knockout direction (treated mean − control mean) changes the
  head's predicted PD-L1 by −0.135 [−0.156, −0.113] per full displacement in Frangieh: the measured sign, and larger
  than all 12 random directions of the same length. In Papalexi it is +0.067 [+0.048, +0.086]: the opposite sign to the
  measurement, at the 58th percentile of random directions. Consistent with the inverse result that the head does not
  keep this contrast, most clearly in Papalexi.
- Gene-level events, layer-by-layer displacements and step-refinement checks are in the report.

### Step 4: closing the loop (rule fixed before the head split)

`anm_diag/declarations/teddy_mm_step4_intervention_rule.md` (sha256 `d47f4660…`, written before the head split was sealed
or run) maps the location of Frangieh PD-L1's loss to an intervention: normalization → add log‖z‖ to the head's input;
first layer → retrain with a declared constraint on that layer; later layers or no single resolved step → retrain with
training data that contains IFN-γ variation. Applied to the head-split record, the root-cause span and its only resolved
step both fall in the last group, so the rule selects **retraining with IFN-γ-varied cells**. The criterion: the change
against a control with the same recipe minus the change, at least 3 seeds each; it counts only if every treated seed
brings HEAD_LAT to bounded support and the control does not, with a guardrail of no more than 0.02 drop in mean
per-protein Pearson on teddy_mm's site4 test. This step has not run.

## Corrections and deviations

| when | what | effect | status |
|---|---|---|---|
| 8 Nov, after the first sealed run | **Head input not normalized.** The adapter fed TEDDY's raw gene-mean embedding (L2 norm ≈ 26) to the phase-1 head; teddy_mm trains and evaluates the head on embeddings divided by their L2 norm (`scripts/04_train.py` line 101, `scripts/05_eval.py` line 55). | The first run's head stages tested an operator teddy_mm never uses; they are void (`runs/frangieh2021`, `runs/papalexi2021`). RAW, TOKENS, Z, the probes and the positive controls do not use the head. | Fixed in `PerturbCiteAdapter.head_input`, checked element-wise against teddy_mm's torch code (within 1e-5), `tests/test_head_operator.py`. Same designs resealed together before any rerun result, each run once (`runs/*_headfix`). Root causes and loss locations did not change; magnitudes did (HEAD_LAT, HEAD_OUT, first run → rerun): Frangieh PD-L1 +0.31 → +0.20, +0.58 → +0.46; MHC-I +0.12 → +0.01, +0.56 → +0.19 (output from reject to unresolved); HLA-E +1.19 → +0.94, +1.15 → +1.36; Papalexi within 0.04. The bad input amplified Frangieh's loss but was not its source. Because the fix came after results were seen, a strict reader may treat the head-stage conclusions as exploratory. |
| 8 Nov, rerun | The rerun ran on a second Mac (Python 3.11.13 vs 3.11.15, identical package versions). | RAW, TOKENS and Z match the first run to three decimals, so only the head input changed between runs. | — |
| 8–9 Nov, seals | The head-fix, head-split and first forward seals record, for `teddy_ckpt`, the sha256 of an empty file (`e3b0c442…`). The machine's Documents folder syncs through iCloud, and at sealing time the checkpoint was an undownloaded placeholder that read as zero bytes; loading uses memory mapping, which downloads it, so runs used the real file. | The weights never changed (sha256 `ba0a98f9…`, 284,887,424 bytes). Head-fix and head-split do not load TEDDY (they use the stored states); the forward run recomputes Z from the weights and matches the stored states to within 1e-6, which ties the states to the real checkpoint. | The sealing code now checks bytes read against the file size and refuses to seal otherwise (anm_diag 0.3.0). |
| 8 Nov, first sealed run | Papalexi was sealed by the pipeline after the Frangieh results existed (the script sealed and ran one dataset at a time). | The design was fixed before any result: the sealed file's design fields match the declaration written before the run, field by field (md5 `514693684f…`). | The pipeline now seals all declarations before running any; the rerun's two declarations were sealed together. |
| 9 Nov | Two-sided reread of all records (above). | MHC-I head claim and KS-2 ribosome `no_deficit` do not survive; HLA-E's step moves. | Exploratory; the sealed verdicts stand as recorded. |
| KS-2 | First pass judged each stage by the class with the best point estimate; the MT positive control then failed by 0.003 (`control_failed`) although ridge recovered it exactly. | Reread with the union rule (support if any declared class supports) gives `loss_at` TOKENS → Z for MT. | The union rule is now the package's rule everywhere; because it was adopted after the first pass, KS-2 stays exploratory. |
| before sealing | The files' `ensembl_id` column maps HLA-A/B/C/E, B2M and IFNGR2 to alternate-haplotype ids outside TEDDY's vocabulary. | Used as is, MHC genes would vanish before TEDDY and look like a RAW → TOKENS loss. | Fixed from metadata only, before sealing: 10x GRCh38-3.0.0 gene table (`tools/gene_map.py`). |

## Limitations

- **Few treated targets.** Five per teddy_mm dataset, six and twelve in KS-2. Folds are leave-one-target-out and the
  bootstrap is over targets, so intervals are coarse; most unresolved verdicts reflect this.
- **Bounded claims.** A reject covers the declared readout classes (ridge, gradient boosting); a support covers the
  tested knockouts, readouts and cell context.
- **Support tests an average.** θ is a weighted average of the same-state differences, so opposite-signed differences
  can cancel; a support verdict checks a necessary condition of the paper's mean sufficiency, while a reject is
  unaffected (given correctly fitted readout models).
- **One-sided sealed rule.** See the two-sided reread.
- **KS-2 is exploratory** and its readouts share counts with TEDDY's input (declared; `count_split` is implemented but
  the embeddings have not been recomputed on split counts).
- **Distribution shift for the head.** The head is applied outside its training distribution (bone marrow → melanoma,
  THP-1).
- **One model, one pooling.** Only TEDDY-G 70M with gene-mean pooling was tested.

## What can and cannot be claimed

**Can claim (sealed teddy_mm tests):** for the IFN-γ-pathway knockouts in Frangieh 2021, TEDDY-G 70M's gene-mean
embedding retains the information that explains the PD-L1 protein response within a 0.10 margin (ridge and
gradient-boosting readouts), while the PD-L1 transcript alone does not; no TEDDY stage is rejected for any of the five
protein readouts; and the teddy_mm phase-1 head, applied to these melanoma cells, loses the PD-L1 response between the
embedding and its latent layer (*u* from +0.04 to +0.20 [+0.13, +0.26]).

**Relative (sealed rule):** the same Z → HEAD_LAT loss for Papalexi PD-L1 and PD-L2 and for Frangieh HLA-E and MHC-I,
with the input view not within tolerance (the MHC-I claim does not survive the two-sided reread; HLA-E's step moves to
HEAD_LAT → HEAD_OUT).

**Exploratory:** the KS-2 pooling result; the two-sided reread; the head split (losses in the later layers and decoder, none at the normalization or the first layer); the forward test.

**Cannot claim:** why the head loses the response (the head split locates, it does not explain; the step-4 intervention has not run); that the head loses it on bone-marrow cells;
that TEDDY is sufficient for perturbation responses in general; that a rejected stage is "broken" (a reject names a
step of the map under declared readout classes).

## Reproduce

See [`anm_diag/README.md`](../anm_diag/README.md#reproducing-the-teddy_mm-tests-perturb-cite-seq-and-eccite-seq):
`bash tools/mac_teddy_mm.sh` from `anm_diag/` downloads the scPerturb files (Zenodo 10044268), checks MD5s, seals,
tokenizes, runs TEDDY and runs each formal test once; `ANM_RUN_TAG=_headfix` reproduces the rerun's file names. The
package here has the code hash recorded in the `_headfix` seals. `python -m pytest -q tests` runs the test suite.
