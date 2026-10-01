# teddy_mm v3 pre-registration

**Status:** frozen before any site4 (test) evaluation of v3. Registration file `registration/registration_v3.json`, sha256 `{{SHA256}}` (also in `registration_v3.json.sha256`). This document is rendered from that file by `bridge_anm/v3_render_registration.py`, so every number below is the registered value.

**What this registers.** The answer key, the evidence, the panels, the questions and their bars, the compared methods, the metrics and the decision rules for six experiments (E1–E6). Every threshold, panel, question, key rule and analysis choice was fixed on the training cells (sites 1–3) and the validation donor (18303, site1) only. `bridge_anm/v3_leakage_check.py` proves that no site4 label, measured site4 protein or site4 embedding value entered any of them (section 12).

**Who uses it.** The builders of E1–E6 read every registered number through `bridge_anm/lib/v3_key.py` (`load_registration` refuses a file whose sha256 differs). A builder that needs a number not given here may fix it only by a procedure written here, only from train/val, and must commit it as an addendum before opening site4 (section 9).

---

## 0. Choices at a glance

| choice | registered value | why |
|---|---|---|
| classes | B, T, NK, myeloid, OUT | the four lineages the evidence can name, plus an outcome where "no call" is correct (the v1/v2 key had neither NK nor OUT) |
| primary key | cells where the annotation class and the protein-gated class agree | two independent-ish sources; a cell the two disagree on is unscored, never forced into a class |
| gate thresholds | 2-component Gaussian mixture on the non-zero training values, posterior 0.5 | label-free; chosen from a declared 12-variant grid by agreement on the 8 training batches, confirmed on val |
| evidence | head prediction / training 95th percentile of the head's own predictions, clipped to [0, 1] | no measured-protein array enters the evidence path; a constant size factor cancels |
| panel | 3 proteins per class, selected on val by worst-pair AUROC | the panel must separate each class from its hardest neighbour (NK from T), not from easy cells |
| Q1 / Q2 / Q3 | soft 4-class / strict 4-class / CD14-anchored (B, T, NK, classical monocyte) | Q3 replaces the synthetic B/T-priority question with a validated, biologically standard one |
| bars | val no-call rates declared up front: Q1 0.15, Q2 0.30, Q3 0.15 | the bar is a declared operating point, not tuned for accuracy |
| ANM | all events at t = 0, equal weights, ANM's default field | snapshot evidence has no time order; with this setup ANM equals the re-coded rule by construction, and E1 checks that cell by cell |
| comparisons | matched coverage; two-stage donor-then-cell bootstrap; win needs the margin overall and in each primary donor | accuracy at different call rates cannot be compared; there are only 2 primary donors |

---

## 1. Fairness contract

1. **Choices on train/val only.** Training cells (`split == train`, sites 1–3) and the validation donor 18303 (`split == val`, site1) are the only cells whose protein, annotation or embedding was read to fix anything. Site4 rows were dropped right after loading; only their split, site and donor metadata were kept to define the evaluation splits and draw label-free subsets.
2. **Same information for every method.** Every compared method reads the same evidence values of the same cells. No method reads the measured protein or the annotation of an evaluated cell; those form the key only.
3. **Matched coverage.** Accuracy is compared only at equal coverage (section 9).
4. **Honest names.** "TEDDY + ANM" only where ANM's field runs; a trust-weighted average or any closed form is a "declared rule". "Zero labels to change the question" holds for any written rule and is reported as such.
5. **Pre-registered margins and failure criteria** for every experiment (section 10); failures and inconclusive results are reported.
6. **Donor-level uncertainty.** Two-stage bootstrap (donors, then cells within donors) plus per-donor rows. With 2 primary donors the interval describes the observed donors and cannot support population claims.

---

## 2. Data inspected

`data/processed/cite/cite_arrays.npz`: {{data.n_cells|,}} cells, {{data.n_adt}} measured proteins (CLR-normalised values), {{data.n_cell_types_total}} annotated cell types (40 of them occur in site4). Absent from the protein panel (checked): {{data.adt_absent_checked}}. So no gate or panel can use CD34 (stem cells), CD235a (erythroid), CD138 (plasma cells) or a pan-γδ TCR antibody; CD71 stands in for erythroid, and plasma cells and γδ T cells are handled through the annotation map.

Embedding: the official TEDDY-G gene-mean embedding `data/processed/cite_official/z_rna.npy` (counts/total × 10⁴, divided by TEDDY gene medians, top 2,048 tokens, fp16, mean over gene tokens; manifest in `z_rna_manifest.json`). Head: phase 1, `outputs/cite_phase1_official/best.pt` (MLP + negative-binomial decoder mean, no flow matching).

Annotated cell types per split (train / val) are listed with the class map in section 4.1. Cell types absent from val: NK CD158e1+, CD8+ T CD57+ CD45RO+, ILC, Plasma cell IGKC-, both plasmablast types, cDC1, CD8+ T naive CD127+ CD26- CD101-, T prog cycling; val has 929 reticulocytes, so its OUT share (26% of val cells by annotation) is higher than training's.

---

## 3. Splits

{{TABLE:splits}}

- **Primary evaluation: site4 donors not in training** (13272: {{splits.test_primary.per_donor.13272|,}} cells; 19593: {{splits.test_primary.per_donor.19593|,}} cells). Why: these are the only cells unseen in both donor and site.
- **Secondary: site4 donor 15078** ({{splits.test_secondary.n|,}} cells). This donor is in training at sites 1–3, so these cells test a new site but not a new donor. Reported with the same rules; never changes a primary verdict.
- **Val: donor 18303 at site1**, a training site. Used for every selection and bar. It is a single donor at site1, the site where the protein gate agrees best with the annotation (section 4.3); this limits every val-selected choice.
- **Cell id:** the global row index in `cite_arrays.npz` (bridge export id `cite_site4_<index>`). `v3_key.split_indices` recomputes each split from metadata and checks its count and index hash.

---

## 4. Answer key

### 4.1 Annotation class (all {{data.n_cell_types_total}} types)

{{TABLE:annotation_map}}

Borderline types and why:

- **gdT CD158b+ → T.** γδ T cells are T by lineage (TCR, CD3 protein). They also carry NK receptors (CD158b, CD56, CD94), which is why they look NK-like in both predicted and measured protein. The protein gate decides cell by cell: CD3-high cells enter the key as T, CD3-low ones gate NK and are unscored. A sensitivity key drops the type entirely (`primary_no_gdT158`).
- **ILC and ILC1 → OUT.** Innate lymphoid cells have no TCR and are not conventional NK cells; a B/T/NK/myeloid call on them is wrong. In this dataset many annotated ILC are NK-like and most annotated ILC1 are CD3-protein high (section 4.3), so they mostly gate NK or T and become unscored; the ones that pass no lineage gate are keyed OUT.
- **pDC → OUT.** Plasmacytoid DCs carry none of the lineage panels (CD123/CD303 high, CD11c low) and are not monocytes or conventional DCs.
- **cDC2 (and cDC1) → myeloid.** Conventional DCs are myeloid (CD11c, CD33 high). cDC1 has no val or site4 cells.
- **Plasma cells and plasmablasts → OUT.** Antibody-secreting cells are a terminal B state outside the mature-B question (CD20 low, CD38 very high); CD138 is not measured. Those that are CD19-high gate B and are unscored.
- **Progenitors (HSC, Lymph prog, G/M prog, MK/E prog) and the erythroid lineage → OUT.** No lineage answer is correct.
- **dnT, T prog cycling, rare CD4/CD8 subsets → T**; the three rare ones have no site4 cells and only affect training labels.

### 4.2 Protein-gated class

Declared logic on measured CLR protein, applied in this order (a cell takes the first rule it passes; "high" means strictly above the frozen threshold):

1. **OUT (erythroid)** if CD71 is high.
2. **T** if CD3 is high, no B marker is high and CD14 is not high.
3. **B** if CD3 is not high, a B marker is high and CD14 is not high.
4. **NK** if CD3, the B markers, CD14 and CD33 are not high and at least {{gate.nk_min_high}} of the NK markers are high.
5. **myeloid** if CD3 and the B markers are not high and any of CD14, CD11c, CD33 is high.
6. **OUT** otherwise (no lineage gate passed).

Chosen variant: B markers = {{gate.b_markers}}; NK markers = {{gate.nk_markers}} (CD56, CD94 and CD335 as the owner asked, plus CD16 for the CD56-dim CD16+ NK majority; CD33 excludes CD16+ monocytes); threshold method = `{{gate.threshold_method}}` fitted on the training cells (seed {{gate.threshold_seed}}).

{{TABLE:gate_thresholds}}

**How the variant was chosen.** The builder scores a fixed 12-variant grid (explored on train/val during development, section 13a): threshold method (Otsu on non-zero values; 2-component Gaussian mixture on non-zero values, posterior 0.5; mixture background mean + 3 SD) × B markers (CD19 or CD20, or CD19 only) × NK markers (with or without CD16). Thresholds were fitted on the training cells; each variant was scored by the mean Cohen kappa (gate vs annotation, 5 classes) over the 8 training batches (site × donor); the best was then confirmed on val. All 12 rows:

{{TABLE:gate_grid}}

**Why one global threshold per protein, and the per-batch alternative.** The training sites stain very differently: in annotated T cells the median measured CD3 ranges from about 1.0 (site3) to 3.1 (site1), so any single threshold behaves differently per site (the chosen gate's lowest training-batch kappa is {{gate.selection.chosen_min_train_batch_kappa}}; table 4.3 lists every batch). The owner asked for thresholds fixed on train/val, so the primary gate uses one frozen threshold per protein. A per-batch variant (the same estimator re-fitted within each site × donor batch, label-free) is registered as a **sensitivity key** only, because on site4 it would read site4's own protein to set the thresholds (verifier side only, never evidence). Training kappa with per-batch thresholds: {{key_quality.per_batch_gate_sensitivity.train.kappa5}}; val: {{key_quality.per_batch_gate_sensitivity.val.kappa5}}.

### 4.3 Primary key and agreement

Primary key = the shared class where annotation and gate agree; otherwise `unscored` (reported separately: the call table of every method on unscored cells, by annotated type). For OUT cells the correct action is no call.

{{TABLE:key_quality}}

Per batch (train and val):

{{TABLE:key_batches}}

Primary-key counts:

{{TABLE:key_counts}}

Val quality targets (written into the builder after the exploratory train/val runs of 13a; a quality bar, not a selection criterion): kappa ≥ {{gate.val_target.kappa5_min|g}}, NK precision ≥ {{gate.val_target.nk_precision_min|g}}, erythroid cells gated OUT ≥ {{gate.val_target.erythroid_out_min|g}}. Met: kappa {{key_quality.val_target_met.kappa5}}, NK precision {{key_quality.val_target_met.nk_precision}}, erythroid OUT {{key_quality.val_target_met.erythroid_out}} (val erythroid OUT rate {{key_quality.val.erythroid_gated_OUT_rate}}).

Per annotated type (train + val):

{{TABLE:key_per_type}}

**On site4** the E1 builder reports, before any method result, the same agreement table, kappa and per-type table for test_primary and test_secondary. Pre-registered key-validity flag: if site4 kappa < {{gate.val_target.kappa5_min|g}}, every result is reported with equal prominence on the annotation-only key. The primary key is not switched on that basis (switching would use site4 labels for a choice).

### 4.4 Sensitivity keys

- `annotation_only`: annotation class for every cell (OUT included).
- `gate_only`: gated class for every cell.
- `primary_per_batch`: agreement of annotation with the per-batch gate (4.2).
- `primary_no_gdT158`: primary key with gdT CD158b+ cells unscored.
- `q3`: the Q3 key (section 7).

---

## 5. Evidence

Embedding: {{evidence.teddy_head.embedding}}. Head: {{evidence.teddy_head.head}}. Evidence for protein *p* = prediction / q95_p, clipped to [0, 1], where q95_p is the 95th percentile of the head's prediction for *p* over the training cells (134 values in the JSON). The decoder's size factor is the constant 1.0: {{evidence.teddy_head.why_size_factor_cancels}}. This replaces v2's measured-protein p95 and its train-median size factor (0.92663), which made the "RNA-only" prediction depend on measured-protein scale factors.

---

## 6. Panels

**Procedure (val only).** Rule: {{panels.selection.rule}}. Score: {{panels.selection.score}}. Candidates: {{panels.selection.candidates}}. Val cells used: {{panels.selection.n_val_cells_used|,}}.

Why worst-pair AUROC. A 4-class call fails at its hardest confusion (NK vs T, B vs B-committed progenitors, myeloid vs G/M progenitors). One-vs-rest AUROC is dominated by easy negatives and saturates on val: many proteins sit within 0.002 of 1, so its ranking among them is noise. Cohen's d depends on scale and on the clip. Both tables are in the dev log. The first build used one-vs-rest AUROC (section 13a); the switch was made on val, before site4.

**NK candidates the owner asked to consider:**

{{TABLE:nk_candidates}}

Top 10 per class by worst-pair AUROC (val):

{{TABLE:top10}}

**Registered panels:**

{{TABLE:panels}}

Notes. CD62P (P-selectin) in the myeloid panel is a platelet protein. Its predicted value separates myeloid cells on val, probably because monocyte–platelet complexes give monocyte RNA a CD62P signal in training; it is kept because the declared rule chose it, and the canonical rows test the dependence. Stability check: the same rule on 20,000 training cells (in-sample for the head) gives B {{panels.selection.stability_same_rule_on_training_cells.panel.B}}, T {{panels.selection.stability_same_rule_on_training_cells.panel.T}}, NK {{panels.selection.stability_same_rule_on_training_cells.panel.NK}}, myeloid {{panels.selection.stability_same_rule_on_training_cells.panel.myeloid}}. The canonical rows reuse the v1/v2 nine markers, which were chosen on site4 Pearson; they are sensitivity rows only.

---

## 7. Questions

Rule score for class *k*: the equal-weight mean of the evidence over the class panel (Q1, Q2) or the evidence of the class anchor (Q3). Call the argmax if the top score reaches the bar; otherwise no call. Bars are the val quantiles matching the declared val no-call rate (over all val cells).

{{TABLE:questions}}

- **Q1 soft.** {{questions.Q1.text}} The declared no-call rate (0.15) is below val's OUT share, so it is soft by design: it calls many OUT cells (val OUT decline {{questions.Q1.val_outcome_rule.out_decline_rate}}).
- **Q2 strict.** {{questions.Q2.text}} This is the nested readout pair used in C2.
- **Q3 CD14-anchored.** {{questions.Q3.text}} Why: it replaces the synthetic B/T-priority question (whose "myeloid" key was mostly NK cells) with the standard four-population counting question (CD19, CD3, CD56, CD14), fixed by biology rather than selected. It changes the correct answer for non-classical monocytes and DCs, so a classifier trained for Q1 cannot answer it without new labels. That makes Experiment 5's label-cost measurement informative.
  - Q3 key: {{questions.Q3.key.rule}}.
  - Val validation of the CD14 gate within primary-key myeloid cells (targets written before this check was first computed; train+val gate profiles per type had been seen: precision ≥ {{questions.Q3.val_target.cd14_gate_precision_for_classical_min|g}}, recall ≥ {{questions.Q3.val_target.cd14_gate_recall_for_classical_min|g}}): precision {{questions.Q3.val_key_check.cd14_gate_precision_for_classical}}, recall {{questions.Q3.val_key_check.cd14_gate_recall_for_classical}}, specificity for other myeloid {{questions.Q3.val_key_check.cd14_gate_specificity_for_other_myeloid}}. Val Q3 key: {{questions.Q3.val_key_check.n_q3_myeloid_val}} classical-monocyte, {{questions.Q3.val_key_check.n_q3_out_from_myeloid_val}} other-myeloid → OUT, {{questions.Q3.val_key_check.n_q3_unscored_from_myeloid_val}} unscored.
  - Predicted CD14 separates classical from other myeloid cells on val with AUROC {{questions.Q3.anchor_auroc_val.CD14.auroc_classical_vs_other_myeloid_val}}.

---

## 8. Methods (same evidence for all)

- **TEDDY + fixed rule (re-coded question):** the rule of section 7.
- **TEDDY + ANM:** ANM `finite_graph_scalar` with ANM's default field (steps {{anm.field_representation.steps}}, retention {{anm.field_representation.retention|g}}, diffusion {{anm.field_representation.diffusion|g}}, source scale {{anm.field_representation.source_scale|g}}). One action per class; one support event per panel protein with value = its evidence; all events at t = 0; equal weights. Because the field is linear and lineages are uncoupled, the action score is field_gain(n) · n · S exactly (field_gain(3) = {{anm.field_gain.3|.6f}}, field_gain(1) = {{anm.field_gain.1|.6f}}). With the readout threshold field_gain(n) · n · bar (`v3_key.anm_readout_threshold`), ANM's calls equal the rule's on every cell. This is a property of the setup, not a finding. E1 checks it cell by cell, and any mismatch counts as a failure.
- **ANM closure readout (C2 only):** ANM's `readout_coordinates` readout over each action's event sites (closure weight {{anm.closure_readout.closure_weight|g}}, mean {{anm.closure_readout.mean_coordinate_weight|g}}, direct {{anm.closure_readout.direct_action_weight|g}}; ANM defaults); val bars at the Q1/Q2 no-call rates: {{anm.closure_bar_Q1|.6f}} / {{anm.closure_bar_Q2|.6f}}. Its re-coded closed form is reported beside it and must match it on every cell.
- **TEDDY + trained classifier:** {{classifier.primary.model}} on the 12 primary-panel evidence values, trained on the training-cell primary key (5 classes, OUT included; {{classifier.primary.n_train_labelled|,}} labelled cells; unscored excluded). C chosen by val log-loss from {{classifier.primary.C_grid}}: C = {{classifier.primary.C|g}}. Q3 classifier: same model on the primary panel plus the Q3 anchors, Q3 labels, C = {{classifier.q3.C|g}}. Secondary classifier: {{classifier.secondary.model}}, no tuning.

{{TABLE:classifier}}

  Decision: call the most probable lineage if its probability reaches the bar (val quantile at the question's no-call rate: Q1 {{classifier.primary.bars.Q1|.6f}}, Q2 {{classifier.primary.bars.Q2|.6f}}). The Q1 bar is tiny: 15% of val cells have a highest lineage probability below it, i.e. the classifier is nearly certain they are OUT. Matching Q1's no-call rate therefore forces it to call some OUT-probable cells, and the matched-coverage comparisons (section 9) handle this the same way for every method.
- **TEDDY margin:** S_top1 − S_top2, a confidence score used in Experiments 3 and 4.

---

## 9. Metrics, matched coverage and statistics

- **Coverage** = fraction of all cells of the evaluated split that get a call (unscored cells included; label-free).
- **Selective accuracy** = among called, scored cells, the fraction whose call equals the key; a call on a key-OUT cell is wrong. **Decision accuracy** = among scored cells, a correct call or no call on a key-OUT cell. **OUT decline rate** = key-OUT cells with no call. Per-class recall and a call table per annotated type are always reported; calls on unscored cells are reported separately.
- **Matched coverage:** at coverage c each method calls its ⌈c · N⌉ most confident cells (its own scores only; ties by cell index), c ∈ {0.95, 0.90, …, 0.50}. **AURC** = trapezoid area of selective accuracy over that grid ÷ 0.45. Each method's deployed bar is also reported with its realised site4 coverage.
- **Bootstrap:** two-stage (donors with replacement, then cells or pairs within each drawn donor), B = 2000, seed {{seeds.bootstrap}}; percentile 95% interval; per-donor rows for 19593, 13272 and (secondary) 15078.
- **Decision rules (all experiments):** *win* = point difference ≥ margin, bootstrap lower bound > 0, and the difference has the same sign and is ≥ margin/2 in each primary donor; *loss* = the mirror; *equivalent* = interval inside (−margin, +margin); otherwise *inconclusive*.
- **Addenda:** numbers an experiment needs and this file does not give (e.g. E4's trust values, E2's z-probe C) are fixed by the procedure written for them, from train/val only, written to `registration/addenda/<experiment>.json` with its sha256 in `registration/addenda/HASHES.txt`, and committed before site4 is opened.
- **Outputs:** `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/v3/<experiment>/` (`<experiment>_results.json` carrying `registration_sha256`, `REPORT.md`, `progress.log`).

---

## 10. Experiments: endpoints, comparators, margins, failure criteria

The machine-readable version is `registration_v3.json → experiments`. Plain summary:

### E1 Mode B rerun on the v3 key (Experiments 1, 3, 4, 5 and C2)

Cells: test_primary (secondary: test_secondary). Keys: primary (Q1, Q2), q3 (Q3); sensitivity keys and panels as in 4.4 and 6.

- **Experiment 1, change the question.** Arms: rule; ANM; the Q1 classifier answering Q3 without relabelling; the Q3 classifier (all training labels, a label-rich ceiling).
  - E1.1a: number of cells where ANM's call differs from the rule's. Registered value 0; any mismatch is a bridge failure, and interpretation stops until it is fixed.
  - E1.1b: share of calls changing Q1→Q2 and Q1→Q3, per type (descriptive).
  - E1.1c / E1.1d: Q3 selective accuracy at the rule's realised Q3 coverage, rule vs each classifier. Margin 0.01.
  - Zero labels to change the question holds for any written rule; it is not an ANM property.
- **Experiment 3, why this call.** ANM leave-one-out over the called class's events vs the rule's closed form (the panel protein with the largest evidence in the called class).
  - E1.3a: agreement of the two top markers (expected 1.000; below 0.99 is reported as a discrepancy).
  - E1.3b: per-donor top-marker shares, against a within-class marker-identity permutation null (1000 permutations; the v2 slot-preserving null is retired).
  - E1.3c: error rate of flip-sensitive calls vs the same number of lowest-margin calls. Win at margin 0.02; otherwise attribution adds no flagging value beyond TEDDY's margin.
- **Experiment 4, which calls to trust.** Scores: ANM soft score (top action score), TEDDY margin, classifier probability, score entropy; questions Q1 and Q2.
  - E1.4a: AURC. E1.4b: selective accuracy at coverage 0.9 and 0.7.
  - Primary comparison: ANM vs TEDDY margin, margin 0.005. Unless ANM wins, the gate claim stays withdrawn.
- **Experiment 5, label cost.** The Q3 classifier is trained on n random Q3-labelled training cells, n ∈ {0 (the Q1 classifier), 25, 50, 100, 200, 500, 1000, 2000, 5000, all}, with 5 draws each.
  - E1.5: the smallest n whose median selective accuracy at the rule's Q3 coverage is ≥ the rule's − 0.005.
  - This is a measurement, not a contest; the rule needs 0 labels as well.
- **C2, nested readouts.** Outer readout P_f = a call is made; inner readout Q_f = the call is right; Q2 calls are nested in Q1 calls. Arms: ANM closure readout, the mean rule (= ANM action readout), the re-coded closure rule (must equal ANM closure), the classifier.
  - E1.C2: key-OUT decline rate at matched coverage c ∈ {rule's realised Q1 coverage, 0.8, 0.7}, with in-scope selective accuracy beside it.
  - Win: ≥ 0.05 more OUT declines than the mean rule, without losing more than 0.005 in-scope accuracy. Any win belongs to the readout form, which a written rule also implements.
  - Sensitivity: annotation-only OUT.
- **Falsification:** "ANM adds decision value in Mode B" is falsified for v3 if C2 and Experiment 4 are not wins.

### E2 Response decomposition (C3) and linearity (C5)

Cells: the registered label-free subset ({{e2_subset.rule}}). Pipeline: official TEDDY-G preprocessing (as the manifest) → frozen head → registered evidence → Q1 rule and ANM.

- **Perturbations:**
  - RNA binomial thinning to 0.95, 0.90, 0.80, 0.60, 0.20 and 0.05 of the UMIs (ε = 1 − kept), two seeds.
  - Gene scaling of each listed gene by (1 − ε), ε ∈ {1.0 (mask), 0.5, 0.25}, in cells that express it. Genes: the coding genes of the 12 panel proteins (MS4A1, CD22, TNFRSF13C, CD3E/D/G, CD2, CD5, IL2RB, KLRD1, NCAM1, SIRPA, ITGAX, SELP) and 13 NK/T genes (NCAM1, KLRD1, NCR1, FCGR3A, KLRF1, KLRC1, SH2D1B, CD3E, CD3D, CD3G, CD5, CD6, CD28).
  - One unperturbed re-embedding per cell, for the fp16 floor.
- **Measurements:** δz, δ evidence (12 panel values), δ rule score and margin, δ ANM score, δ call.
- **Linearity (C5):** pairs (ε, ε/2) = (0.1, 0.05), (0.2, 0.1), (0.4, 0.2), (0.8, 0.4) for thinning and (1.0, 0.5), (0.5, 0.25) for gene scaling. Slope ratio ρ = ‖δ(ε)‖ / (2‖δ(ε/2)‖) and the cosine between the two responses; a cell is linear at a level if |ρ − 1| ≤ 0.1, cosine ≥ 0.9 and ‖δ(ε/2)‖ > 3 × floor. Reported per stage (z, evidence, score).
- **Decomposition (C3):** each cell that was right and called at baseline and is lost after the perturbation is labelled
  - *representation* if a z-space probe (logistic regression on training z, C on val, addendum) is wrong on the perturbed z;
  - else *head* if the rule's argmax on the perturbed evidence is wrong;
  - else *decision* (the right class falls below the bar).
- **Comparator:** the plain perturbation → accuracy curve.
- **Falsified (adds nothing beyond the curve)** unless either holds:
  - (i) two perturbations with accuracy losses within 0.02 of each other differ by ≥ 0.15 in a localisation share;
  - (ii) the small-ε response predicts which cells change call at ε = 0.8 with AUROC ≥ 0.02 above baseline margin alone (bootstrap lower bound > 0).

### E3 NK–T look-alike pairs (C4) and repair

- **Pairs:** {{e3.pairs_rule}}. Flag: {{e3.flag_rule}}, with flag cosine {{e3.flag_cosine|.6f}} from {{e3.n_val_nkt_pairs}} val NK–T neighbour pairs. Variants: all pairs, and no_gdT158.
- **Readouts on frozen TEDDY, trained on training sites only:**
  - our head;
  - R1, an MLP on the 12 layer means;
  - R2, attention pooling over layer-12 gene-token states, trained on 4,000 training cells;
  - a null: the head architecture retrained with seed 1 (regression to the mean, instability).
- **Gap ratio:** median |predicted difference| / median |measured difference| within pairs for CD56, CD94, CD335 and CD3, averaged over the four.
- **E3.H3a:** D_R = [GR_R − GR_head] on flagged pairs minus the same on unflagged pairs, minus the null's D. Must be ≥ 0.05 (win rules). Pooling hypothesis: D_R2 − D_R1 ≥ 0.05.
- **E3.H3b:** NK-vs-T selective accuracy at coverage 0.9 and 0.8 with each readout's evidence (margin 0.01). Trust comparators: head margin, entropy, kNN label disagreement among training cells.
- **Falsified** ("the NK–T loss is in the readout or pooling and is repairable on frozen TEDDY") if neither R1 nor R2 reaches D ≥ 0.05 over the null.

### E4 Two-channel fusion (C6; Experiment 2 redesigned)

- **Key:** the E4 key from the measured primary-panel proteins (Q1 scoring on measured evidence; OUT below τ_K chosen on val by kappa against annotation; addendum). The v3 primary key is secondary and is partly circular for channel 2, disclosed.
- **Channels:**
  - channel 1 = TEDDY + head prediction of the panel from RNA;
  - channel 2 = ridge regression (alpha on val) from the 122 measured non-panel proteins to the panel, trained on training cells. It never sees a panel protein's measured value.
- **Noise (on val and test):** L0 none; L1 RNA thinned to 0.2; L2 50% of channel-2 inputs set to 0; L3 both.
- **Trust:** per channel and level on val, the mean max(0, Pearson) with the measured panel; ANM source scale = trust / max trust (addendum).
- **Methods:** each channel alone; simple average and trust-weighted average (declared rules); a stacker trained on val; learned fusion (logistic regression on the 24 values, training cells); and TEDDY + ANM fusion. The ANM arm has per-channel source scales and contradiction events (each channel's top class contradicts the other classes at 0.5 × its top score), all at t = 0, with the bar from val.
- **Controls:** no propagation; channel 2 at t = 2 vs t = 0 (time-blind); no contradiction events (= trust-weighted average); the re-coded closed form.
- **Endpoint:** selective accuracy at coverage 0.8 (and 0.9), averaged over L0–L3, on the E4 key. ANM vs the best non-ANM fusion (primary) and ANM vs trust-weighted average (the field's own effect). Margin 0.01.
- **Outcomes:** "field adds value" (win), "field adds nothing" (equivalent; an acceptable outcome) or "field hurts" (loss).
- **Untestable** if the two channels' val trust differs by < 0.01 at every level.

### E5 Mode A layer dynamics (C7)

- **Cells:** test_primary cells in E3's NK–T pairs, at most 200 per donor.
- **NK–T direction per layer:** mean training key-NK minus mean training key-T gene-mean layer output.
- **Perturbation:** the unit token embedding of each E2 gene the cell expresses, at that gene's position. Controls: random unit directions and random other genes.
- **Response:** δh_l = J_l δh_{l−1} by `torch.func.jvp` through each actual encoder layer (fp32, CPU, eval mode, bool padding mask). The projection on the NK–T direction gives a layer gain.
- **Check:** central finite differences at ε = 1e-2 and 1e-3. Required: relative error ≤ 0.05 and cosine ≥ 0.99 at 1e-3 for ≥ 95% of cases; otherwise E5 stops with "JVP not validated".
- **Endpoint:** the layer l* with the most negative median log-gain of the NK–T component, compared with the layer where the linear NK–T protein probe's gap ratio drops most (`outputs/mode_a_official/mode_a_results.json`). Agreement within 1 layer.
- **Falsified (no layer localises the NK–T loss)** if the NK–T component's log-gain at l* is within 0.1 of random directions'.

### E6 External confirmation (D11)

- **Candidate:** Hao et al. 2021 PBMC CITE-seq (GEO GSE164378, 3′ data). It is on disk under `data/raw/external_hao2021_gse164378` (downloaded 2026-09-30 by another process) and was **not opened** for this registration: no protein names, values, RNA or labels were read.
- **Frozen:** TEDDY and the head, the evidence normaliser, panels, bars, classifier, gate logic and annotation-map principles.
- **Gate thresholds:** CLR scales differ between datasets, so the BMMC thresholds do not transfer. The E6 primary key uses the registered per-batch estimator within each external donor (label-free), and the annotation-only key is reported too.
- **Missing proteins:** an absent gate protein drops its clause; an absent panel protein leaves that class scored on its remaining panel proteins. Both are listed in the addendum.
- **Annotation map:** written from the external label names only and committed (with sha256) before any external value is read.
- **Contamination:** state whether GSE164378 is in TEDDY's pretraining corpus.
- **Endpoints:** E1.1a, E1.4a, E1.C2 and Q1 selective accuracy at matched coverage, per external donor.
- **Replication:** a conclusion replicates if its sign holds in a majority of external donors and the bootstrap interval excludes 0 in the same direction.
- **Key not validated:** if annotation vs per-batch-gate kappa < 0.70, results are reported on the annotation-only key and labelled "key not validated".

---

## 11. Seeds

Global {{seeds.global}}; bootstrap {{seeds.bootstrap}}; classifier {{seeds.classifier}}; label-cost draws {{seeds.label_cost_draws}}; E2 subset {{seeds.e2_subset}}; E2 thinning {{seeds.e2_thinning}}; E3 token subset {{seeds.e3_token_subset}}; E4 noise {{seeds.e4_noise}}; E5 cells {{seeds.e5_cells}}; E5 random directions {{seeds.e5_random_directions}}; gate thresholds {{gate.threshold_seed}}.

---

## 12. Leakage check

`bridge_anm/v3_leakage_check.py` runs three checks, and all must pass:

1. **Static.** The builder reads the protein matrix, cell types and embedding only on non-test rows, and no earlier site4 result file.
2. **Invariance.** The registration core is rebuilt from inputs whose site4 protein, cell types and embedding are replaced by random values. It must be byte-identical to the core from the real inputs and to the core whose sha256 this registration records (`provenance.core_sha256` = `{{provenance.core_sha256}}`).
3. **Positive control.** The same poisoning applied to the val rows must change the core.

Result: `registration/leakage_check_report.json` and `outputs/v3/registration/LEAKAGE_CHECK.md`.

---

## 13. Disclosures

- **Site4 is not virgin.** v1/v2 experiments, the 2026-09-28 audit and an exploratory answer-key prototype (scratchpad, evaluated on site4) all read site4. v3 reuses none of their thresholds or panels. Its design is still informed by what they showed (for example, that gdT CD158b+ cells look NK-like and that the v2 "myeloid" key-marker class was mostly NK). E6 exists for this reason.
- **The frozen head carries two inherited facts:**
  - its checkpoint was selected on val (`val_fm_pearson`);
  - its training size factors were divided by the median ADT total over all 90,261 cells, site4 included. That is one global scalar; v3's evidence cancels any constant scale, and its effect on the trained weights is a near-constant rescaling that is not testable without retraining.
- **Val is one donor at site1**, where the gate agrees best with the annotation. Every val-selected choice may favour site1-like data; the training-batch diagnostics (4.3) show how the gate degrades at site3.
- **The primary key is an agreement subset.** It keeps the cells both sources call cleanly, which are easier than average. Unscored cells and the annotation-only key are always reported.
- **Only 2 primary donors.** Intervals describe them, not a population.
- **Pretraining contamination** of TEDDY with GSE194122 (this dataset) is not verified (open item from the audit).
- **E6 data prepared in parallel.** After this registration was built (provenance `built_utc`), a separate builder committed `scripts/prepare_external_cite.py` (commit 8318ef4), which parses the external set, maps its ADT to our names and assigns a class from its `celltype.l2` labels. The E6 rule above (class map from label names only, committed before any external value is read) applies to that work; whoever runs E6 must show it was met or report the deviation.
- **E2 compute:** about 0.118 s per cell embedding (the official run's rate); about 1,250 cells × (13 + 3 per expressed gene) embeddings.

### 13a. Development history (train and val only; nothing here touched site4)

So that a reader can see which choices moved after val numbers were seen:

1. **Gate.** Exploratory scripts profiled the gate proteins by annotated type on train+val, then compared Otsu, Gaussian-mixture and mixture-background thresholds, global and per batch, with and without CD20 and CD16, on the training batches and val. The builder then fixed the 12-variant grid above and picks by mean training-batch kappa; the exploration had favoured the same variant.
2. **Panel criterion.** The first build selected panels by one-vs-rest AUROC (B CD20, CD22, CD72; T CD2, CD3, CD5; NK CD122, CD94, CD16; myeloid CD172a, CD11c, CD62P). Because that score saturates on val, it was replaced by worst-pair AUROC (section 6); Cohen's d was also computed. All three were on val only.
3. **Q3 key check.** The first check (Q3-key precision against the annotation) was 1 by construction. It was replaced by the precision and recall of the measured CD14 gate within primary-key myeloid cells (section 7).
4. **Rounding.** Every registered number is rounded to 6 decimals when it is defined and used rounded downstream, so the key and the bars `v3_key` computes are exactly the registered ones. After this change the val log-loss choice of C moved from 10 to 100; the two differ by less than 0.001 in val log-loss (section 8).

---

## 14. Reproduce

```bash
cd <repo>
PY=/private/tmp/claude-501/-Users-tianchichen-Documents-GitHub-teddy-mm/0828c10d-3512-46ba-b296-fafd23e4a760/scratchpad/venv312/bin/python
ANM_ROOT=<ANM v2 fix checkout> $PY bridge_anm/v3_build_registration.py      # writes registration_v3.json + .sha256 (train/val only)
$PY bridge_anm/v3_render_registration.py                                     # this document
ANM_ROOT=<...> $PY -m pytest -q tests/test_v3_key.py
$PY bridge_anm/v3_leakage_check.py                                           # proof of no site4 use
```

Builder `{{provenance.builder}}`, sha256 `{{provenance.builder_sha256}}`; `v3_key.py` sha256 `{{provenance.v3_key_sha256}}`; inputs: `cite_arrays.npz` sha256 `{{provenance.inputs.cite_arrays_sha256}}`, `z_rna.npy` sha256 `{{provenance.inputs.z_sha256}}`, head `best.pt` sha256 `{{provenance.inputs.ckpt_sha256}}`; git HEAD at build `{{provenance.git_head_at_build}}`.
