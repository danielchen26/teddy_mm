# ANM diagnosis of TEDDY: where does a perturbation response stop being recoverable?

Code and records: [`anm_diag/`](../anm_diag/) (how-to in its [README](../anm_diag/README.md); every number below is in
`anm_diag/runs/*/anm_record.json`). Method: Active Neural Matter (arXiv 2610.14549), §2.2–§2.3 and Table 2.

> [!IMPORTANT]
> **Status (8 Nov 2026).** The TEDDY-side results (RAW → TOKENS → Z) stand. **The teddy_mm head-stage results
> (HEAD_LAT, HEAD_OUT) are withdrawn**: the first sealed run fed the head un-normalized embeddings (details in
> [Corrections](#corrections-and-deviations)). The head stages are being rerun with the fixed operator; until
> `anm_diag/runs/*_headfix/` exists there is no valid claim about the head.

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
| frozen system | TEDDY-G 70M (12 layers, 8 heads, d = 512), optionally followed by the teddy_mm phase-1 protein head |
| source perturbation *A* | CRISPR knockout / CRISPRi of a target gene vs non-targeting controls; folds and bootstrap by target |
| readout *O* | a **measurement** the model never sees (surface protein by CITE-seq), or a declared gene-set score |
| candidate states (chain) | RAW (normalized counts) → TOKENS (TEDDY's rank tokens) → Z (gene-mean of the last layer, 512-D) → HEAD_LAT → HEAD_OUT |
| probes (off the chain) | SET (the same genes without rank order), single-gene mRNA (e.g. CD274 for PD-L1), positive control (the readout itself) |
| readout classes | ridge and gradient boosting, fixed hyperparameters |
| tolerance | *u* margin 0.10, α = 0.05 |

**The statistic.** For each stage, cross-fit (by target) a model of the readout from *z* and of the treatment from *z*,
and take the slope of the readout residual on the treatment residual:

θ(z) = Σ (O − m̂(z)) (A − ê(z)) / Σ (A − ê(z))²,  *u* = θ(z) / Δ,

where Δ is the readout's raw response (treated minus control). *Intuition:* θ(z) is the part of the knockout's effect
on the readout that is **left over** after everything *z* knows has been used. If *z* retains the response, nothing is
left over (*u* ≈ 0); if *z* has lost it, the whole response is left over (*u* ≈ 1). *u* is a share, so a 0.10 margin
reads "this stage explains at least 90% of the response".

**Verdicts.** Support: upper bound of *u* < 0.10 for some declared readout class. Reject: lower bound > 0.10 for every
class. Otherwise unresolved. A reject is bounded by the readout classes (paper eq. 5): it says ridge and gradient
boosting cannot recover the response from *z*, not that no function could. A support is bounded by the comparisons
tested.

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

*Reading.* The MT response is fully recoverable from TEDDY's rank tokens but not from the same genes without their
order (SET) and not from the pooled embedding (Z). So the information sits in **which genes rank where**, and averaging
token states over genes loses it — consistent with gene-mean pooling averaging away per-position information. For PSM
even the input view does not pass, so nothing is attributed to TEDDY.

**Ribosome knockdowns** (12 targets, 480 treated of 2,196 cells): RP and MRP readouts are supported at every stage
(`no_deficit`); the MT readout is supported at RAW and TOKENS and unresolved at Z (+0.07 [+0.03, +0.11]). Forward AUC
RAW 0.855, TOKENS 0.839, Z 0.817, SET 0.844.

### teddy_mm: TEDDY on Perturb-CITE-seq and ECCITE-seq (sealed, run once)

Protein readouts measured by CITE-seq, never seen by TEDDY or the head. Cells are labelled from their detected guides;
genes mapped to TEDDY's vocabulary with the 10x GRCh38-3.0.0 table (the files' own Ensembl ids drop the MHC genes).

- **Frangieh 2021** (melanoma, IFN-γ condition): knockouts IFNGR1, IFNGR2, JAK1, JAK2, STAT1 (1,171 treated cells)
  vs 4,000 non-targeting controls; readouts PD-L1, MHC-I (HLA-ABC antibody), HLA-E, CLR across the 20-protein panel.
- **Papalexi 2021** (THP-1, IFN-γ stimulated): knockouts IFNGR1, IFNGR2, JAK2, STAT1, IRF1 (7,478 cells); readouts
  PD-L1, PD-L2, log1p (4-protein panel).

**TEDDY stages** (*u*, 95% interval, verdict):

| data | readout | Δ | RAW | TOKENS | Z | SET |
|---|---|---|---|---|---|---|
| Frangieh | PD-L1 | −0.85 [−0.93, −0.76] | −0.01 [−0.06, +0.04] **support** | +0.01 [−0.03, +0.05] **support** | +0.04 [−0.01, +0.09] **support** | +0.01 support |
| Frangieh | MHC-I | −0.14 [−0.18, −0.09] | −0.27 [−0.71, +0.11] unresolved | −0.29 [−0.70, +0.07] **support** | −0.37 [−0.92, +0.02] **support** | −0.30 support |
| Frangieh | HLA-E | −0.08 [−0.11, −0.05] | −0.27 [−1.28, +0.25] unresolved | −0.16 unresolved | +0.18 [−0.49, +0.57] unresolved | unresolved |
| Papalexi | PD-L1 | −0.69 [−0.86, −0.47] | +0.10 [+0.00, +0.23] unresolved | +0.12 [−0.04, +0.37] unresolved | +0.15 [+0.04, +0.32] unresolved | unresolved |
| Papalexi | PD-L2 | −0.29 [−0.34, −0.23] | +0.08 [−0.23, +0.34] unresolved | +0.08 [−0.19, +0.31] unresolved | +0.13 [−0.11, +0.36] unresolved | unresolved |

Forward AUC: Frangieh RAW 0.895, TOKENS 0.889, Z 0.884, SET 0.886; Papalexi RAW 0.767, TOKENS 0.753, Z 0.757,
SET 0.749.

**Single-gene mRNA probes** (the protein's own transcript as the only input): Frangieh PD-L1 from CD274 mRNA
+0.93 [+0.92, +0.95] reject; Papalexi PD-L1 from CD274 +0.69 [+0.63, +0.78] reject; PD-L2 from PDCD1LG2
+0.80 [+0.77, +0.88] reject.

*Reading.*

1. **No TEDDY stage is rejected for any of the five protein readouts.** Where the data resolve the question
   (Frangieh PD-L1; MHC-I from TOKENS on), TEDDY's pooled embedding Z retains at least 90% of the protein response to
   the IFN-γ-pathway knockouts.
2. **The protein response is in the transcriptome but not in the protein's own transcript.** CD274 mRNA explains almost
   none of the PD-L1 protein change (*u* +0.93), the whole RNA profile explains all of it (*u* −0.01). The knockouts act
   through the pathway program, which TEDDY's whole-cell embedding captures and a single-gene readout misses.
3. **Rank order is not needed here** (SET supported for PD-L1 and MHC-I), unlike the KS-2 MT readout. Whether pooling
   loses information depends on the readout.
4. **Papalexi and HLA-E are unresolved, not negative.** Five treated targets give coarse target-level bootstrap
   intervals (K1 warning); even RAW does not reach the margin, so nothing is attributed to TEDDY.

**Head stages: withdrawn.** The first run reported a resolved loss at Z → HEAD_LAT for every readout. That run fed the
head the wrong input scale (see below), so those verdicts do not describe the head; they are excluded here and await the
rerun.

## Corrections and deviations

| when | what | effect | status |
|---|---|---|---|
| 8 Nov, after the first sealed run | **Head input not normalized.** The adapter fed TEDDY's raw gene-mean embedding (L2 norm ≈ 26) to the phase-1 head; teddy_mm trains and evaluates the head on embeddings divided by their L2 norm (`scripts/04_train.py` line 101, `scripts/05_eval.py` line 55). | HEAD_LAT, HEAD_OUT and the Z → HEAD_LAT step are void. RAW, TOKENS, Z, the probes and the positive controls do not use the head and are unaffected. | Fixed in `PerturbCiteAdapter.head_input`, checked element-wise against teddy_mm's torch code (within 1e-5), `tests/test_head_operator.py`. Same declarations resealed and run once into `runs/*_headfix/`; recorded as a formal test with a deviation, since the error was found after seeing results. |
| 8 Nov, first sealed run | Papalexi was sealed by the pipeline after the Frangieh results existed (the script sealed and ran one dataset at a time). | The design was fixed before any result: the sealed file's design fields match the declaration written before the run, field by field (md5 `514693684f…`). | The pipeline now seals all declarations before running any. |
| KS-2 | First pass judged each stage by the class with the best point estimate; the MT positive control then failed by 0.003 (`control_failed`) although ridge recovered it exactly. | Reread with the union rule (support if any declared class supports) gives `loss_at` TOKENS → Z for MT. | The union rule is now the package's rule everywhere; because it was adopted after the first pass, KS-2 stays exploratory. |
| before sealing | The files' `ensembl_id` column maps HLA-A/B/C/E, B2M and IFNGR2 to alternate-haplotype ids outside TEDDY's vocabulary. | Used as is, MHC genes would vanish before TEDDY and look like a RAW → TOKENS loss. | Fixed from metadata only, before sealing: 10x GRCh38-3.0.0 gene table (`tools/gene_map.py`). |

## Limitations

- **Few treated targets.** Five per teddy_mm dataset, six and twelve in KS-2. Folds are leave-one-target-out and the
  bootstrap is over targets, so intervals are coarse; unresolved verdicts mostly reflect this.
- **Bounded claims.** A reject covers the declared readout classes (ridge, gradient boosting); a support covers the
  tested knockouts, readouts and cell context.
- **KS-2 is exploratory** and its readouts share counts with TEDDY's input (declared; `count_split` is implemented
  but the embeddings have not been recomputed on split counts).
- **Distribution shift for the head.** The head was trained on bone-marrow CITE-seq and is applied unchanged to
  melanoma and THP-1 cells, so any head-stage result will be about this head on new cells.
- **One model, one pooling.** Only TEDDY-G 70M with gene-mean pooling was tested.

## What can and cannot be claimed

**Can claim (sealed teddy_mm tests, TEDDY stages):** for the IFN-γ-pathway knockouts in Frangieh 2021, TEDDY-G 70M's
gene-mean embedding retains the information that explains the PD-L1 protein response (and the MHC-I response from the
token stage on), within a 0.10 margin and for ridge and gradient-boosting readouts; that response is not explained by
the PD-L1 transcript alone; and no TEDDY stage was rejected for any of the five protein readouts.

**Exploratory (KS-2):** for mitoribosome knockdowns, the MT gene-set response is carried by TEDDY's rank order and lost
by gene-mean pooling.

**Cannot claim (yet):** anything about where the teddy_mm head loses information (withdrawn pending the rerun); that
TEDDY is sufficient for perturbation responses in general; that a rejected stage is "broken" (a reject names a step of
the map under declared readout classes).

## Reproduce

See [`anm_diag/README.md`](../anm_diag/README.md#reproducing-the-teddy_mm-tests-perturb-cite-seq-and-eccite-seq):
`bash tools/mac_teddy_mm.sh` from `anm_diag/` downloads the scPerturb files (Zenodo 10044268), checks MD5s, seals,
tokenizes, runs TEDDY and runs each formal test once. `python -m pytest -q tests` runs the test suite.
