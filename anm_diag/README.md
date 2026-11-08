# anm_diag: a declaration-driven runner for ANM response-sufficiency diagnoses

Basis: the ANM paper (Active Neural Matter, arXiv 2610.14549), §2.2–§2.3, Table 2 and §8.3.

**What was tested with it and what came out:** [`docs/ANM_DIAGNOSIS.md`](../docs/ANM_DIAGNOSIS.md).
This page is the how-to.

One diagnosis = **one declaration** (YAML; template in `declarations/TEMPLATE.yaml`) + **one adapter** (turns a dataset
into arrays). The declaration names the frozen system, the measurements, the source perturbation, the readouts, the
chain of candidate retained states, the controls, the readout classes and the tolerance. `check` rejects designs that
cannot answer the ANM question before any statistic is computed; `run` computes the forward and inverse tests, the
three-way verdicts and a root cause, and writes a JSON record and a Markdown report. A new system or dataset needs a
new declaration and adapter; the decision rules stay the same.

## Three steps (plus sealing)

Run everything from this folder (relative paths in `paths/*.yaml` are resolved from here).

```bash
cd anm_diag
python -m anm_diag check declarations/ks2_teddy_mitoribo.yaml
python -m anm_diag run declarations/ks2_teddy_mitoribo.yaml --adapter ks2_teddy --paths paths/ks2.yaml --out runs/ks2_mitoribo
python -m anm_diag.reanalyze runs/ks2_mitoribo     # reread with another tolerance or rule, no refitting (exploratory only)
python -m pytest -q tests                          # legality checks, synthetic ground truth, sealing, head operator, adapter smoke test
```

**A formal test is sealed first:** `python -m anm_diag seal decl.yaml --paths paths.yaml --out decl.sealed.yaml`. The
sealed declaration carries every default, the sha256 of each input file (keyed by role, not by path), the sha256 of the
`anm_diag` package source and the sealing time; `run` refuses a sealed declaration that was edited, or whose data or
code changed. Before sealing only metadata may be looked at: `python -m anm_diag inspect --paths paths.yaml` lists obs
column names, the values and counts of design columns (target, condition), var column names and the protein panel. It
never reads a count matrix.

Dependencies: numpy, scipy, scikit-learn, pyyaml; TEDDY tokenization needs torch (the official `torch.topk`); the
Perturb-CITE adapter needs anndata. teddy_mm's `.venv` (`requirements.txt`) has all of them.

## Roles in a declaration (ANM's vocabulary)

| key | ANM role | notes |
|---|---|---|
| `system` | the system under test, held fixed | `frozen: true`; `inputs` lists the measurements it reads |
| `measurements` | measurements (the roots of provenance) | what was measured per cell: RNA counts, ADT counts, … |
| `source` | source perturbation *u* | must be `controlled: true`, with control and treated groups; `cluster` is the unit for folds and resampling (the target gene) |
| `readouts` | readout *O* | taken from a measurement only, never from a model output or a candidate state |
| `nodes` | candidate retained states *z* = *G*(*s*) and probes | `from` builds the computation graph; `operator`: system / analysis / approximation; `role`: chain / probe / positive_control / technical |
| `chain` | the stages on the real computation path, in order | the first is the reference view of the input; every later stage must be the system's own operator |
| `measurement_sharing` | whether a readout and the input share one measurement | if shared, declare `count_split` (preferred) or `declared_shared` |
| `readout_classes` | bounded readout classes | default ridge (RidgeCV / logistic C = 0.05) and gradient boosting, fixed hyperparameters |
| `decision` | decision rule | `margin_u` (default 0.10), `step_margin`, `alpha` |
| `status` | exploratory / sealed | `seal` writes sealed; a design fixed after seeing results can only be exploratory |

## Legality checks (an error refuses to run)

| code | rule | what it prevents |
|---|---|---|
| S1 | a controlled source perturbation with controls is required | a decision layer with no source (teddy_mm Mode B) |
| R1 | a readout must be a measurement, not a node | scoring the model against its own prediction (Mode A used the head's output as the readout) |
| R2 | a readout's measurement must not flow into a candidate state; if it is a system input, a sharing policy must be declared | phase 2 fed measured protein as input and scored against measured protein |
| C1 | each chain stage must be computed from the previous one | comparing stages that are not on one computation path |
| C2 | stages after the reference must be the system's own operators | an approximate tokenizer on the path (a vocabulary mismatch added 340 genes) |
| C3 | no analyst compression (PCA, k-means, …) | artefacts of an undeclared PCA-50 |
| C4 | no input masking, weight edits or state replacement | masking inputs; mean-preserving replacement (a system intervention, not a sufficiency test) |
| C5 | node provenance and roles are legal; warns when no positive control is declared | |
| E1 | the inverse statistic must be a continuous conditional test | a coarse k-means partition |
| K1/K2 | the tolerance must be declared; warns below 10 treated targets; sealed declarations must carry hashes | |
| seal check | `run` verifies the seal digest, data hashes and code hash | changing the design, data or code after sealing |

## Reading the output

For every candidate state × readout the report gives *u* (the share of the readout's response that this stage cannot
explain), a 95% interval and a three-way verdict:

- **support**: upper bound < `margin_u`; **reject**: lower bound > `margin_u` (rejects only the declared readout
  classes, paper eq. 5); otherwise **unresolved**. A stage is judged over the union of the declared readout classes:
  support if any class supports it, reject only if every class rejects it.

With θ(z) the cross-fitted partial-out slope of the readout on the source given *z*
(Σ(O − m̂(z))(A − ê(z)) / Σ(A − ê(z))², folds and bootstrap by target cluster) and Δ the readout's response,
*u* = θ(z)/Δ: 0 means *z* explains the whole response, 1 means it explains none of it.

Root cause per readout (top to bottom; the first match wins):

| root cause | meaning |
|---|---|
| `no_effect` | the readout has no resolvable response (*u* divides by Δ, so no inverse fit) |
| `control_failed` | the positive control failed: the test cannot find information known to be there; no verdict |
| `relative_loss_at` | the reference stage itself is not within tolerance, but one step adds a resolved loss (a weaker, relative claim) |
| `reference_unresolved` | the reference is not within tolerance and no step adds a resolved loss: nothing is attributed to the system. A **rejected** reference means the input itself cannot reproduce the response under the declared classes (a data or readout-class limit, not an internal loss); an **unresolved** one means the data cannot tell |
| `loss_at` | the reference is supported and a later stage is rejected: the distinction merges at this step (a step of the map, not "which part is broken") |
| `unresolved` | a stage after the reference is unresolved |
| `no_deficit` | every stage is supported |

The report also gives the forward AUC (does the perturbation reach each stage, held-out-cluster AUC) and intervals for
the change at each step; the KS-2 adapter adds an exact same-state check (identical tokens, readout moved toward
control).

## Adapters

| adapter | system | data | readout |
|---|---|---|---|
| `ks2_teddy` | TEDDY-G 70M | Replogle 2022 K562 CRISPRi (KS-2 development cells) | mean log1p CP10k of a readout gene set (shares counts with the input; declared) |
| `perturb_cite` | TEDDY-G 70M + teddy_mm phase-1 protein head | Perturb-CITE-seq / ECCITE-seq (RNA + ADT) | measured protein, CLR or log1p (an independent measurement) |
| `synthetic` | toy chain | simulated | checks that the runner finds a planted loss |

A new adapter subclasses `adapters.base.Adapter`, implements `cells` and `readout`, and registers each node `op` in
`ops`. A declaration may only use operations the adapter implements.

## Reproducing the teddy_mm tests (Perturb-CITE-seq and ECCITE-seq)

Declarations: `declarations/teddy_mm_frangieh2021.yaml` (melanoma, IFN-γ condition, knockouts IFNGR1 / IFNGR2 / JAK1 /
JAK2 / STAT1; readouts MHC-I, PD-L1 and HLA-E protein) and `declarations/teddy_mm_papalexi2021.yaml` (THP-1, knockouts
IFNGR1 / IFNGR2 / JAK2 / STAT1 / IRF1; readouts PD-L1 and PD-L2 protein). The chain is RAW → TOKENS → Z → HEAD_LAT →
HEAD_OUT; the readout is measured protein and never enters TEDDY or the head.

**One command (macOS, MPS):** `bash tools/mac_teddy_mm.sh` from this folder. It downloads the four scPerturb files
(Zenodo record 10044268, about 1.7 GB, resumable) and checks the published MD5s → writes metadata only to
`runs/teddy_mm/inspect_*.txt` → waits for `runs/teddy_mm/GO` (written once the declarations are filled from metadata
and pass `check`) → seals all declarations, exports tokens, runs TEDDY and runs each formal test once. Rerunning the
same command skips finished steps; a sealed declaration is never resealed and a finished test is never rerun.
`ANM_RUN_TAG=_headfix` writes new sealed files and result folders with that suffix and leaves earlier ones as they
were. The script assumes this folder sits inside a teddy_mm checkout; set `TEDDY_MM=...` otherwise.

Step by step:

1. **Download** `FrangiehIzar2021_{RNA,protein}.h5ad` and `PapalexiSatija2021_eccite_{RNA,protein}.h5ad` from
   `https://zenodo.org/records/10044268` into `../data/raw/perturb/`, then
   `python tools/verify_downloads.py ../data/raw/perturb` (scPerturb v1.3).
2. **Metadata only:** `python -m anm_diag inspect --paths paths/frangieh2021.yaml`. Fill `obs_target`, `obs_stratum`,
   `obs_guide` and `ensembl_col` (or `symbol_map`) in the paths file and the non-targeting labels, the IFN-γ label and
   the protein names in the declaration. Do not `--show` columns computed from the data (QC metrics, mixscape KO/NP
   calls): they carry outcome information.
3. **Seal:** `python -m anm_diag seal declarations/teddy_mm_frangieh2021.yaml --paths paths/frangieh2021.yaml --out declarations/teddy_mm_frangieh2021.sealed.yaml`
4. **Export tokens:** `python -m anm_diag export-tokens declarations/teddy_mm_frangieh2021.sealed.yaml --adapter perturb_cite --paths paths/frangieh2021.yaml --out frangieh_tokens.npz`
5. **TEDDY states:** `../.venv/bin/python tools/layer_states.py frangieh_tokens.npz ../data/processed/frangieh2021_teddy_states.npz --z-only`
6. **Run once:** `python -m anm_diag run declarations/teddy_mm_frangieh2021.sealed.yaml --adapter perturb_cite --paths paths/frangieh2021.yaml --out runs/frangieh2021`

For ECCITE-seq replace `frangieh2021` with `papalexi2021`.

Design notes fixed before sealing:

- **Cells are labelled from their detected guides** (`source.assign_from_guides`). Controls are cells whose guides are
  all non-targeting; treated cells are those whose guides all hit one target gene. Frangieh: 6,444 controls (4,000
  sampled with a fixed seed), IFNGR1 198, IFNGR2 278, JAK1 248, JAK2 244, STAT1 203 (the Frangieh library has no IRF1).
  Papalexi: one guide per cell, grouped by gene, IFN-γ-stimulated cells only.
- **Genes are mapped to TEDDY's vocabulary with the 10x GRCh38-3.0.0 table** (`tools/gene_map.py`, the table teddy_mm
  used for its own data). The files' `ensembl_id` column gives HLA-A/B/C/E, B2M and IFNGR2 alternate-haplotype ids that
  are not in TEDDY's vocabulary; used as is, the MHC genes would be dropped before TEDDY and misread as a RAW → TOKENS
  loss.
- **Protein statistic chosen in advance:** CLR across the panel for Frangieh (20 proteins); log1p for Papalexi (4
  proteins, where CLR would make unchanged proteins move).
- **The head is used as is.** It was trained on bone-marrow CITE-seq; melanoma and THP-1 cells are outside its training
  distribution, so the head stages test "this head applied to new cells".
- **Few treated targets** (5 per dataset): K1 warns; folds are leave-one-target-out and the cluster-bootstrap intervals
  are coarse.
- **Memory:** the adapter reads only the selected cells (h5ad opened backed). Set `OMP_NUM_THREADS=1` when running
  several diagnoses at once.

## Layout

```text
anm_diag/
├── anm_diag/        # the package (its source hash is part of every seal)
├── declarations/    # TEMPLATE, KS-2 (exploratory), teddy_mm (draft + sealed)
├── paths/           # paths files (relative to this folder); ks2.example.yaml to fill in
├── runs/            # JSON records and Markdown reports of finished runs (see each folder's README)
├── tests/           # pytest suite
└── tools/           # macOS pipeline, gene map, TEDDY layer states, download check
```

> **中文简介**：anm_diag 是按 ANM 论文（§2.2–§2.3、表 2、§8.3）写的通用诊断运行器：一份声明（YAML）+ 一个适配器。`check` 在算任何统计之前拒绝答不了 ANM 问题的设计；`seal` 把设计、数据哈希和代码哈希封存；`run` 给出每个候选状态 × 读出的 *u*、区间、三分判定（support / reject / unresolved）和根因。用它做过的检验和结果见 [`docs/ANM_DIAGNOSIS.md`](../docs/ANM_DIAGNOSIS.md)。
