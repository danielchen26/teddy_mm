# teddy_mm v3 pre-registration

**Status:** frozen before any site4 (test) evaluation of v3. Registration file `registration/registration_v3.json`, sha256 `e4c8a33e5c7c0b292d3a70a9063730c7ea0fd6d1519a705953bcf73ad8212bd3` (also in `registration_v3.json.sha256`). This document is rendered from that file by `bridge_anm/v3_render_registration.py`, so every number below is the registered value.

**Amended before site4 (A1).** An adversarial review found ten problems in the experiment specifications (tie-breaking at matched coverage, a panel-dependent normaliser in E4 channel 2, the E1.3 leave-one-out closed form and tied top markers, E1.4 naming, the C2 verdict point, classifier information in E1.1c/E1.5, the E2 falsification guard, E3 pair construction, the E4 comparator and the scope of addenda). They are fixed in `registration/amendment_A1.json` (rendered as `registration/AMENDMENT_A1.md`), which has its own sha256 and amends this file without changing it. Where the two differ, A1 wins; sections 9 and 10 below show the text as first registered.

**What this registers.** The answer key, the evidence, the panels, the questions and their bars, the compared methods, the metrics and the decision rules for six experiments (E1–E6). Every threshold, panel, question, key rule and analysis choice was fixed on the training cells (sites 1–3) and the validation donor (18303, site1) only. `bridge_anm/v3_leakage_check.py` proves that no site4 label, measured site4 protein or site4 embedding value entered any of them (section 12).

**Who uses it.** The builders of E1–E6 read every registered number through `bridge_anm/lib/v3_key.py` and load the registration with `bridge_anm/lib/v3_amend.py`'s `load_registration_amended()`, which checks the sha256 of this file and of amendment A1 and applies A1's experiment fields. A builder that needs a number not given here may fix it only by a procedure written here, only from train/val, and must commit it as an addendum before opening site4 (section 9).

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

`data/processed/cite/cite_arrays.npz`: 90,261 cells, 134 measured proteins (CLR-normalised values), 45 annotated cell types (40 of them occur in site4). Absent from the protein panel (checked): CD34, CD235a, CD138, TCRgd. So no gate or panel can use CD34 (stem cells), CD235a (erythroid), CD138 (plasma cells) or a pan-γδ TCR antibody; CD71 stands in for erythroid, and plasma cells and γδ T cells are handled through the annotation map.

Embedding: the official TEDDY-G gene-mean embedding `data/processed/cite_official/z_rna.npy` (counts/total × 10⁴, divided by TEDDY gene medians, top 2,048 tokens, fp16, mean over gene tokens; manifest in `z_rna_manifest.json`). Head: phase 1, `outputs/cite_phase1_official/best.pt` (MLP + negative-binomial decoder mean, no flow matching).

Annotated cell types per split (train / val) are listed with the class map in section 4.1. Cell types absent from val: NK CD158e1+, CD8+ T CD57+ CD45RO+, ILC, Plasma cell IGKC-, both plasmablast types, cDC1, CD8+ T naive CD127+ CD26- CD101-, T prog cycling; val has 929 reticulocytes, so its OUT share (26% of val cells by annotation) is higher than training's.

---

## 3. Splits

| split | rule | cells | index sha256 (first 16) |
|---|---|---:|---|
| train | split == 'train' (sites 1-3, all donors except 18303) | 67,405 | `0da6d39783b7ef21` |
| val | split == 'val' (site1, donor 18303) | 6,106 | `0e2787236c64d689` |
| test_primary | site4 cells of donors not in training | 11,294 | `dd5cf5727b2939a3` |
| test_secondary | site4 cells of donor 15078 (in training at sites 1-3) | 5,456 | `4dfab424f5987983` |

- **Primary evaluation: site4 donors not in training** (13272: 7,365 cells; 19593: 3,929 cells). Why: these are the only cells unseen in both donor and site.
- **Secondary: site4 donor 15078** (5,456 cells). This donor is in training at sites 1–3, so these cells test a new site but not a new donor. Reported with the same rules; never changes a primary verdict.
- **Val: donor 18303 at site1**, a training site. Used for every selection and bar. It is a single donor at site1, the site where the protein gate agrees best with the annotation (section 4.3); this limits every val-selected choice.
- **Cell id:** the global row index in `cite_arrays.npz` (bridge export id `cite_site4_<index>`). `v3_key.split_indices` recomputes each split from metadata and checks its count and index hash.

---

## 4. Answer key

### 4.1 Annotation class (all 45 types)

| annotated cell type | class | train | val | note |
|---|---|---:|---:|---|
| B1 B IGKC+ | B | 563 | 41 |  |
| B1 B IGKC- | B | 398 | 48 |  |
| Naive CD20+ B IGKC+ | B | 2,458 | 265 |  |
| Naive CD20+ B IGKC- | B | 1,167 | 171 |  |
| Transitional B | B | 768 | 250 |  |
| CD4+ T CD314+ CD45RA+ | T | 92 | 1 | T; 92 training cells, 1 val, none in site4. |
| CD4+ T activated | T | 5,032 | 547 |  |
| CD4+ T activated integrinB7+ | T | 684 | 141 |  |
| CD4+ T naive | T | 3,924 | 833 |  |
| CD8+ T CD49f+ | T | 615 | 49 |  |
| CD8+ T CD57+ CD45RA+ | T | 692 | 251 |  |
| CD8+ T CD57+ CD45RO+ | T | 431 | 0 |  |
| CD8+ T CD69+ CD45RA+ | T | 474 | 33 |  |
| CD8+ T CD69+ CD45RO+ | T | 517 | 32 |  |
| CD8+ T TIGIT+ CD45RA+ | T | 843 | 25 |  |
| CD8+ T TIGIT+ CD45RO+ | T | 934 | 37 |  |
| CD8+ T naive | T | 1,454 | 878 |  |
| CD8+ T naive CD127+ CD26- CD101- | T | 42 | 0 | T; 42 training cells, none in val or site4. |
| MAIT | T | 379 | 109 |  |
| T prog cycling | T | 24 | 0 | T: annotated as a T-committed cycling population; 24 training cells, none in val or site4. |
| T reg | T | 392 | 36 |  |
| dnT | T | 23 | 33 | T (CD3+ CD4- CD8- T cells). |
| gdT CD158b+ | T | 180 | 62 | T by lineage (gamma-delta TCR, CD3 protein). They carry NK receptors (CD158b, CD56, CD94); the protein gate decides whether each one enters the key as T (CD3 high) or is unscored. Sensitivity key 'primary_no_gdT158' drops them. |
| gdT TCRVD2+ | T | 161 | 10 | T (Vdelta2 gamma-delta T; CD3 protein high). |
| NK | NK | 3,779 | 271 |  |
| NK CD158e1+ | NK | 1,861 | 0 |  |
| CD14+ Mono | myeloid | 20,171 | 309 |  |
| CD16+ Mono | myeloid | 2,348 | 28 |  |
| cDC1 | myeloid | 18 | 0 | myeloid (conventional DC, like cDC2); 18 training cells, none in val or site4. |
| cDC2 | myeloid | 1,449 | 54 | myeloid (conventional DC; CD11c and CD33 protein high). |
| Erythroblast | OUT | 3,472 | 228 | OUT: erythroid lineage. |
| G/M prog | OUT | 1,472 | 45 | OUT: progenitor. |
| HSC | OUT | 1,246 | 71 | OUT: stem cells; no lineage panel applies. |
| ILC | OUT | 135 | 0 | OUT: innate lymphoid cells lack a TCR and are not conventional NK cells. Many are NK-like in protein (CD56/CD94 high); they then gate NK and are unscored, not keyed NK. |
| ILC1 | OUT | 365 | 40 | OUT: as ILC. In this dataset most annotated ILC1 are CD3-protein high (T-like); those gate T and are unscored. |
| Lymph prog | OUT | 931 | 87 | OUT: progenitor (B-committed in part; many are CD19 protein high and then unscored). |
| MK/E prog | OUT | 486 | 8 | OUT: progenitor. |
| Normoblast | OUT | 1,160 | 129 | OUT: erythroid lineage. |
| Plasma cell IGKC+ | OUT | 237 | 9 | OUT: antibody-secreting cells are a terminal B-lineage state outside the mature-B panel question (CD19 partly, CD20 low, CD38 very high). |
| Plasma cell IGKC- | OUT | 220 | 0 | OUT, as plasma cells. |
| Plasmablast IGKC+ | OUT | 216 | 0 | OUT, as plasma cells. |
| Plasmablast IGKC- | OUT | 109 | 0 | OUT, as plasma cells. |
| Proerythroblast | OUT | 865 | 26 | OUT: erythroid lineage. |
| Reticulocyte | OUT | 2,972 | 929 | OUT: erythroid lineage. |
| pDC | OUT | 1,646 | 20 | OUT: plasmacytoid DCs are not monocytes or conventional DCs and carry none of the four lineage panels (CD123/CD303 high, CD11c low); a lineage call on them is wrong. |

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
4. **NK** if CD3, the B markers, CD14 and CD33 are not high and at least 2 of the NK markers are high.
5. **myeloid** if CD3 and the B markers are not high and any of CD14, CD11c, CD33 is high.
6. **OUT** otherwise (no lineage gate passed).

Chosen variant: B markers = CD19; NK markers = CD56, CD94, CD335, CD16 (CD56, CD94 and CD335 as the owner asked, plus CD16 for the CD56-dim CD16+ NK majority; CD33 excludes CD16+ monocytes); threshold method = `gmm_nonzero_post50` fitted on the training cells (seed 0).

| protein | threshold (CLR) | role in the chosen gate |
|---|---:|---|
| CD3 | 0.680750 | T |
| CD19 | 0.622993 | B marker |
| CD20 | 1.254582 | not used (B-marker variant not chosen) |
| CD56 | 1.063468 | NK marker |
| CD94 | 1.532667 | NK marker |
| CD335 | 0.895803 | NK marker |
| CD16 | 0.648140 | NK marker |
| CD14 | 0.601134 | myeloid; excludes T, B, NK; Q3 key |
| CD33 | 0.593140 | myeloid; excludes NK |
| CD11c | 0.436573 | myeloid |
| CD71 | 2.090560 | erythroid -> OUT |

**How the variant was chosen.** The builder scores a fixed 12-variant grid (explored on train/val during development, section 13a): threshold method (Otsu on non-zero values; 2-component Gaussian mixture on non-zero values, posterior 0.5; mixture background mean + 3 SD) × B markers (CD19 or CD20, or CD19 only) × NK markers (with or without CD16). Thresholds were fitted on the training cells; each variant was scored by the mean Cohen kappa (gate vs annotation, 5 classes) over the 8 training batches (site × donor); the best was then confirmed on val. All 12 rows:

| threshold method | B markers | NK markers | mean train-batch kappa | min train-batch kappa | val kappa |
|---|---|---|---:|---:|---:|
| gmm_nonzero_post50 | CD19 | CD56, CD94, CD335, CD16 | 0.8225 | 0.7091 | 0.8905 |
| gmm_nonzero_post50 | CD19, CD20 | CD56, CD94, CD335, CD16 | 0.8210 | 0.7073 | 0.8891 |
| gmm_background_mu_plus_3sd | CD19 | CD56, CD94, CD335, CD16 | 0.8145 | 0.6919 | 0.8990 |
| gmm_background_mu_plus_3sd | CD19, CD20 | CD56, CD94, CD335, CD16 | 0.8141 | 0.6919 | 0.8980 |
| gmm_nonzero_post50 | CD19 | CD56, CD94, CD335 | 0.8022 | 0.6836 | 0.8900 |
| gmm_nonzero_post50 | CD19, CD20 | CD56, CD94, CD335 | 0.8008 | 0.6817 | 0.8886 |
| gmm_background_mu_plus_3sd | CD19 | CD56, CD94, CD335 | 0.7876 | 0.6489 | 0.8977 |
| gmm_background_mu_plus_3sd | CD19, CD20 | CD56, CD94, CD335 | 0.7872 | 0.6492 | 0.8967 |
| otsu_nonzero | CD19, CD20 | CD56, CD94, CD335, CD16 | 0.7706 | 0.4870 | 0.9369 |
| otsu_nonzero | CD19 | CD56, CD94, CD335, CD16 | 0.7667 | 0.4847 | 0.9381 |
| otsu_nonzero | CD19, CD20 | CD56, CD94, CD335 | 0.7578 | 0.4682 | 0.9369 |
| otsu_nonzero | CD19 | CD56, CD94, CD335 | 0.7538 | 0.4658 | 0.9381 |

**Why one global threshold per protein, and the per-batch alternative.** The training sites stain very differently: in annotated T cells the median measured CD3 ranges from about 1.0 (site3) to 3.1 (site1), so any single threshold behaves differently per site (the chosen gate's lowest training-batch kappa is 0.7091; table 4.3 lists every batch). The owner asked for thresholds fixed on train/val, so the primary gate uses one frozen threshold per protein. A per-batch variant (the same estimator re-fitted within each site × donor batch, label-free) is registered as a **sensitivity key** only, because on site4 it would read site4's own protein to set the thresholds (verifier side only, never evidence). Training kappa with per-batch thresholds: 0.8166; val: 0.9258.

### 4.3 Primary key and agreement

Primary key = the shared class where annotation and gate agree; otherwise `unscored` (reported separately: the call table of every method on unscored cells, by annotated type). For OUT cells the correct action is no call.

| split | cells | agreement | kappa (5 classes) | B recall / precision | T | NK | myeloid | OUT |
|---|---:|---:|---:|---|---|---|---|---|
| train | 67,405 | 0.8702 | 0.8254 | 0.801 / 0.797 | 0.853 / 0.935 | 0.784 / 0.929 | 0.927 / 0.941 | 0.857 / 0.730 |
| val | 6,106 | 0.9284 | 0.8905 | 0.788 / 0.926 | 0.972 / 0.980 | 0.937 / 1.000 | 0.890 / 0.895 | 0.921 / 0.836 |

Per batch (train and val):

| batch (site, donor) | cells | agreement | kappa |
|---|---:|---:|---:|
| site1 / 10886 | 4,978 | 0.8831 | 0.8371 |
| site1 / 15078 | 5,227 | 0.9218 | 0.8945 |
| site1 / 18303 | 6,106 | 0.9284 | 0.8905 |
| site2 / 12710 | 5,584 | 0.8582 | 0.8022 |
| site2 / 15078 | 10,465 | 0.9026 | 0.8700 |
| site2 / 16710 | 9,122 | 0.9248 | 0.8951 |
| site3 / 11466 | 11,473 | 0.8592 | 0.7617 |
| site3 / 15078 | 9,521 | 0.8655 | 0.8100 |
| site3 / 28045 | 11,035 | 0.7858 | 0.7091 |

Primary-key counts:

| split | B | T | NK | myeloid | OUT | unscored |
|---|---:|---:|---:|---:|---:|---:|
| train | 4,289 | 14,411 | 4,419 | 22,234 | 13,304 | 8,748 |
| val | 611 | 2,990 | 254 | 348 | 1,466 | 437 |

Val quality targets (written into the builder after the exploratory train/val runs of 13a; a quality bar, not a selection criterion): kappa ≥ 0.85, NK precision ≥ 0.9, erythroid cells gated OUT ≥ 0.75. Met: kappa yes, NK precision yes, erythroid OUT yes (val erythroid OUT rate 1.0000).

Per annotated type (train + val):

| annotated type | class | train+val cells | gated B / T / NK / myeloid / OUT | in primary key |
|---|---|---:|---|---:|
| B1 B IGKC+ | B | 604 | 387 / 0 / 0 / 1 / 216 | 387 |
| B1 B IGKC- | B | 446 | 281 / 2 / 0 / 0 / 163 | 281 |
| Naive CD20+ B IGKC+ | B | 2,723 | 2282 / 6 / 0 / 12 / 423 | 2,282 |
| Naive CD20+ B IGKC- | B | 1,338 | 1130 / 4 / 0 / 3 / 201 | 1,130 |
| Transitional B | B | 1,018 | 820 / 5 / 0 / 2 / 191 | 820 |
| CD4+ T CD314+ CD45RA+ | T | 93 | 0 / 92 / 0 / 0 / 1 | 92 |
| CD4+ T activated | T | 5,579 | 4 / 4860 / 2 / 33 / 680 | 4,860 |
| CD4+ T activated integrinB7+ | T | 825 | 2 / 721 / 0 / 3 / 99 | 721 |
| CD4+ T naive | T | 4,757 | 1 / 4337 / 0 / 3 / 416 | 4,337 |
| CD8+ T CD49f+ | T | 664 | 0 / 601 / 0 / 0 / 63 | 601 |
| CD8+ T CD57+ CD45RA+ | T | 943 | 2 / 750 / 41 / 21 / 129 | 750 |
| CD8+ T CD57+ CD45RO+ | T | 431 | 2 / 370 / 2 / 3 / 54 | 370 |
| CD8+ T CD69+ CD45RA+ | T | 507 | 2 / 390 / 1 / 11 / 103 | 390 |
| CD8+ T CD69+ CD45RO+ | T | 549 | 4 / 426 / 7 / 10 / 102 | 426 |
| CD8+ T TIGIT+ CD45RA+ | T | 868 | 0 / 769 / 3 / 6 / 90 | 769 |
| CD8+ T TIGIT+ CD45RO+ | T | 971 | 1 / 869 / 0 / 3 / 98 | 869 |
| CD8+ T naive | T | 2,332 | 2 / 2089 / 2 / 11 / 228 | 2,089 |
| CD8+ T naive CD127+ CD26- CD101- | T | 42 | 0 / 40 / 0 / 0 / 2 | 40 |
| MAIT | T | 488 | 1 / 412 / 2 / 6 / 67 | 412 |
| T prog cycling | T | 24 | 0 / 12 / 0 / 0 / 12 | 12 |
| T reg | T | 428 | 0 / 287 / 0 / 1 / 140 | 287 |
| dnT | T | 56 | 0 / 48 / 0 / 0 / 8 | 48 |
| gdT CD158b+ | T | 242 | 0 / 166 / 72 / 0 / 4 | 166 |
| gdT TCRVD2+ | T | 171 | 0 / 162 / 0 / 0 / 9 | 162 |
| NK | NK | 4,050 | 47 / 185 / 3102 / 418 / 298 | 3,102 |
| NK CD158e1+ | NK | 1,861 | 23 / 72 / 1571 / 115 / 80 | 1,571 |
| CD14+ Mono | myeloid | 20,480 | 57 / 148 / 9 / 19253 / 1013 | 19,253 |
| CD16+ Mono | myeloid | 2,376 | 25 / 74 / 40 / 2178 / 59 | 2,178 |
| cDC1 | myeloid | 18 | 2 / 2 / 0 / 13 / 1 | 13 |
| cDC2 | myeloid | 1,503 | 26 / 80 / 8 / 1138 / 251 | 1,138 |
| Erythroblast | OUT | 3,700 | 0 / 0 / 0 / 0 / 3700 | 3,700 |
| G/M prog | OUT | 1,517 | 23 / 29 / 2 / 476 / 987 | 987 |
| HSC | OUT | 1,317 | 34 / 27 / 6 / 88 / 1162 | 1,162 |
| ILC | OUT | 135 | 2 / 6 / 121 / 2 / 4 | 4 |
| ILC1 | OUT | 405 | 0 / 344 / 3 / 3 / 55 | 55 |
| Lymph prog | OUT | 1,018 | 570 / 10 / 2 / 12 / 424 | 424 |
| MK/E prog | OUT | 494 | 1 / 3 / 0 / 8 / 482 | 482 |
| Normoblast | OUT | 1,289 | 0 / 0 / 0 / 0 / 1289 | 1,289 |
| Plasma cell IGKC+ | OUT | 246 | 98 / 4 / 0 / 12 / 132 | 132 |
| Plasma cell IGKC- | OUT | 220 | 85 / 1 / 2 / 8 / 124 | 124 |
| Plasmablast IGKC+ | OUT | 216 | 74 / 2 / 2 / 6 / 132 | 132 |
| Plasmablast IGKC- | OUT | 109 | 32 / 1 / 2 / 4 / 70 | 70 |
| Proerythroblast | OUT | 891 | 0 / 1 / 0 / 15 / 875 | 875 |
| Reticulocyte | OUT | 3,901 | 0 / 0 / 1 / 0 / 3900 | 3,900 |
| pDC | OUT | 1,666 | 19 / 55 / 7 / 151 / 1434 | 1,434 |

**On site4** the E1 builder reports, before any method result, the same agreement table, kappa and per-type table for test_primary and test_secondary. Pre-registered key-validity flag: if site4 kappa < 0.85, every result is reported with equal prominence on the annotation-only key. The primary key is not switched on that basis (switching would use site4 labels for a choice).

### 4.4 Sensitivity keys

- `annotation_only`: annotation class for every cell (OUT included).
- `gate_only`: gated class for every cell.
- `primary_per_batch`: agreement of annotation with the per-batch gate (4.2).
- `primary_no_gdT158`: primary key with gdT CD158b+ cells unscored.
- `q3`: the Q3 key (section 7).

---

## 5. Evidence

Embedding: official TEDDY-G gene-mean z (data/processed/cite_official/z_rna.npy), L2-normalised. Head: phase-1 MLP + NB decoder mean (no flow matching), outputs/cite_phase1_official/best.pt. Evidence for protein *p* = prediction / q95_p, clipped to [0, 1], where q95_p is the 95th percentile of the head's prediction for *p* over the training cells (134 values in the JSON). The decoder's size factor is the constant 1.0: prediction = softplus(.) * s for a constant s, and the training q95 scales by the same s, so the evidence is independent of s; no measured-protein array enters the evidence path. This replaces v2's normaliser (the training 95th percentile of measured protein) and its train-median size factor (0.92663, a constant from measured training-cell protein totals). Neither carried per-cell test information, but both put measured-protein scales into the evidence path; v3 has none.

---

## 6. Panels

**Procedure (val only).** Rule: each protein is assigned to the class where its worst-pair AUROC is highest; a class panel is its top 3 assigned proteins with worst-pair AUROC >= 0.8. Score: worst-pair AUROC: for class k, the minimum over every other primary-key class c (OUT included, classes with >= 10 val cells) of the AUROC of the unclipped normalised head prediction for k versus c; ties by one-vs-rest AUROC, then name. Candidates: all 134 measured-panel proteins. Val cells used: 5,669.

Why worst-pair AUROC. A 4-class call fails at its hardest confusion (NK vs T, B vs B-committed progenitors, myeloid vs G/M progenitors). One-vs-rest AUROC is dominated by easy negatives and saturates on val: many proteins sit within 0.002 of 1, so its ranking among them is noise. Cohen's d depends on scale and on the clip. Both tables are in the dev log. The first build used one-vs-rest AUROC (section 13a); the switch was made on val, before site4.

**NK candidates the owner asked to consider:**

| protein | worst-pair AUROC (NK) | NK vs rest AUROC | NK vs T AUROC | assigned class | rank for NK |
|---|---:|---:|---:|---|---:|
| CD56 | 0.9975 | 0.9986 | 0.9975 | NK | 3 |
| CD94 | 0.9982 | 0.9990 | 0.9982 | NK | 2 |
| CD335 | 0.9930 | 0.9979 | 0.9997 | NK | 4 |
| CD16 | 0.9921 | 0.9987 | 0.9985 | NK | 5 |

Top 10 per class by worst-pair AUROC (val):

**B**: CD20 0.9931; CD22 0.9907; CD268 0.9895; CD72 0.9892; CD19 0.9814; HLA-DR 0.9632; CD79b 0.9544; CD21 0.9533; CD24 0.9356; IgM 0.8420

**T**: CD3 0.9983; CD2 0.9965; CD5 0.9937; TCR 0.9872; CD278 0.9769; CD48 0.9719; CD27 0.9548; CD44 0.9295; CD127 0.9268; CD45 0.8884

**NK**: CD122 0.9991; CD94 0.9982; CD56 0.9975; CD335 0.9930; CD16 0.9921; CD161 0.9858; CD328 0.9833; CD244 0.9719; CD45RA 0.9532; CD11a 0.8922

**myeloid**: CD172a 0.9997; CD11c 0.9991; CD62P 0.9989; CD93 0.9989; CD41 0.9985; CD64 0.9940; CD86 0.9918; CD33 0.9870; CD13 0.9512; CD31 0.9452

**Registered panels:**

| panel | B | T | NK | myeloid | role |
|---|---|---|---|---|---|
| primary | CD20, CD22, CD268 | CD3, CD2, CD5 | CD122, CD94, CD56 | CD172a, CD11c, CD62P | primary |
| canonical9_plus_nk | CD19, CD72, CD22 | CD3, CD2, CD5 | CD122, CD94, CD56 | CD16, CD11c, CD36 | old canonical 9-marker panel (chosen on site4 Pearson in v1/v2) plus the val-selected NK panel minus any protein already in it |
| canonical9_variant_plus_nk | CD19, CD20, CD22 | CD3, CD2, CD5 | CD122, CD94, CD56 | CD16, CD11c, CD14 | canonical panel with CD20 for CD72 and CD14 for CD36, plus the val-selected NK panel |
| val_k2 | CD20, CD22 | CD3, CD2 | CD122, CD94 | CD172a, CD11c | same rule, 2 per class |
| val_k5 | CD20, CD22, CD268, CD72, CD19 | CD3, CD2, CD5, TCR, CD278 | CD122, CD94, CD56, CD335, CD16 | CD172a, CD11c, CD62P, CD93, CD41 | same rule, 5 per class |

Notes. CD62P (P-selectin) in the myeloid panel is a platelet protein. Its predicted value separates myeloid cells on val, probably because monocyte–platelet complexes give monocyte RNA a CD62P signal in training; it is kept because the declared rule chose it, and the canonical rows test the dependence. Stability check: the same rule on 20,000 training cells (in-sample for the head) gives B CD22, CD72, CD20, T CD3, CD5, TCR, NK CD56, CD122, CD94, myeloid CD172a, CD11c, CD41. The canonical rows reuse the v1/v2 nine markers, which were chosen on site4 Pearson; they are sensitivity rows only.

---

## 7. Questions

Rule score for class *k*: the equal-weight mean of the evidence over the class panel (Q1, Q2) or the evidence of the class anchor (Q3). Call the argmax if the top score reaches the bar; otherwise no call. Bars are the val quantiles matching the declared val no-call rate (over all val cells).

| question | evidence per class | bar | val no-call target | val no-call realised | val accuracy of calls | val decision accuracy | val OUT decline |
|---|---|---:|---:|---:|---:|---:|---:|
| Q1 soft 4-class lineage question | primary panel, equal weights | 0.298614 | 0.15 | 0.1500 | 0.8657 | 0.8871 | 0.6166 |
| Q2 strict 4-class lineage question | primary panel, equal weights | 0.661255 | 0.30 | 0.3000 | 0.9784 | 0.9425 | 0.9925 |
| Q3 CD14-anchored question (B / T / NK / classical monocyte) | anchor: B CD19, T CD3, NK CD56, myeloid CD14 | 0.273973 | 0.15 | 0.1500 | 0.8549 | 0.8777 | 0.5867 |

- **Q1 soft.** From the TEDDY + head evidence, is this cell B, T, NK or myeloid? Call the class with the highest panel score if that score reaches the bar; otherwise no call. The declared no-call rate (0.15) is below val's OUT share, so it is soft by design: it calls many OUT cells (val OUT decline 0.6166).
- **Q2 strict.** As Q1 with a higher bar (nested: every Q2 call is a Q1 call with the same class). This is the nested readout pair used in C2.
- **Q3 CD14-anchored.** Is this cell a B cell (CD19), T cell (CD3), NK cell (CD56) or a classical CD14+ monocyte (CD14)? Score each class by its one canonical anchor protein; call the highest if it reaches the bar; otherwise no call. Non-classical monocytes and dendritic cells are not an answer (no call). Why: it replaces the synthetic B/T-priority question (whose CD16-anchored "myeloid" class captured many NK cells; audit finding KEY-1) with the standard four-population counting question (CD19, CD3, CD56, CD14), fixed by biology rather than selected. It changes the correct answer for non-classical monocytes and DCs, so a classifier trained for Q1 cannot answer it without new labels. That makes Experiment 5's label-cost measurement informative.
  - Q3 key: primary-key myeloid: annotated classical monocyte AND measured CD14 high -> myeloid; annotated other myeloid AND measured CD14 not high -> OUT (no call correct); else unscored.
  - Val validation of the CD14 gate within primary-key myeloid cells (targets written before this check was first computed; train+val gate profiles per type had been seen: precision ≥ 0.9, recall ≥ 0.7): precision 0.9795, recall 0.9828, specificity for other myeloid 0.8947. Val Q3 key: 286 classical-monocyte, 51 other-myeloid → OUT, 11 unscored.
  - Predicted CD14 separates classical from other myeloid cells on val with AUROC 0.9986.

---

## 8. Methods (same evidence for all)

- **TEDDY + fixed rule (re-coded question):** the rule of section 7.
- **TEDDY + ANM:** ANM `finite_graph_scalar` with ANM's default field (steps 4, retention 0.82, diffusion 0.16, source scale 1). One action per class; one support event per panel protein with value = its evidence; all events at t = 0; equal weights. Because the field is linear and lineages are uncoupled, the action score is field_gain(n) · n · S exactly (field_gain(3) = 0.370041, field_gain(1) = 0.364476). With the readout threshold field_gain(n) · n · bar (`v3_key.anm_readout_threshold`), ANM's calls equal the rule's on every cell. This is a property of the setup, not a finding. E1 checks it cell by cell, and any mismatch counts as a failure.
- **ANM closure readout (C2 only):** ANM's `readout_coordinates` readout over each action's event sites (closure weight 1, mean 0.25, direct 0.05; ANM defaults); val bars at the Q1/Q2 no-call rates: 0.175321 / 0.342623. Its re-coded closed form is reported beside it and must match it on every cell.
- **TEDDY + trained classifier:** multinomial logistic regression (sklearn LogisticRegression, lbfgs, max_iter 3000) on the 12 primary-panel evidence values, trained on the training-cell primary key (5 classes, OUT included; 58,657 labelled cells; unscored excluded). C chosen by val log-loss from 0.01, 0.1, 1.0, 10.0, 100.0: C = 100. Q3 classifier: same model on the primary panel plus the Q3 anchors, Q3 labels, C = 100. Secondary classifier: sklearn MLPClassifier(hidden_layer_sizes=(64,), alpha=1e-4, early_stopping=True, validation_fraction=0.1, max_iter=500, random_state=seed), no tuning.

| C | val log-loss (primary classifier) | val log-loss (Q3 classifier) |
|---:|---:|---:|
| 0.01 | 0.097658 | 0.073782 |
| 0.1 | 0.060488 | 0.045281 |
| 1.0 | 0.051998 | 0.042980 |
| 10.0 | 0.045805 | 0.037127 |
| 100.0 | 0.045278 | 0.035183 |

  Decision: call the most probable lineage if its probability reaches the bar (val quantile at the question's no-call rate: Q1 0.001137, Q2 0.928654). The Q1 bar is tiny: 15% of val cells have a highest lineage probability below it, i.e. the classifier is nearly certain they are OUT. Matching Q1's no-call rate therefore forces it to call some OUT-probable cells, and the matched-coverage comparisons (section 9) handle this the same way for every method.
- **TEDDY margin:** S_top1 − S_top2, a confidence score used in Experiments 3 and 4.

---

## 9. Metrics, matched coverage and statistics

- **Coverage** = fraction of all cells of the evaluated split that get a call (unscored cells included; label-free).
- **Selective accuracy** = among called, scored cells, the fraction whose call equals the key; a call on a key-OUT cell is wrong. **Decision accuracy** = among scored cells, a correct call or no call on a key-OUT cell. **OUT decline rate** = key-OUT cells with no call. Per-class recall and a call table per annotated type are always reported; calls on unscored cells are reported separately.
- **Matched coverage:** at coverage c each method calls its ⌈c · N⌉ most confident cells (its own scores only; ties by cell index), c ∈ {0.95, 0.90, …, 0.50}. **AURC** = trapezoid area of selective accuracy over that grid ÷ 0.45. Each method's deployed bar is also reported with its realised site4 coverage.
- **Bootstrap:** two-stage (donors with replacement, then cells or pairs within each drawn donor), B = 2000, seed 1; percentile 95% interval; per-donor rows for 19593, 13272 and (secondary) 15078.
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

Cells: the registered label-free subset (seeded (seed 7) uniform sample without labels: 500 cells per primary donor, 250 from donor 15078). Pipeline: official TEDDY-G preprocessing (as the manifest) → frozen head → registered evidence → Q1 rule and ANM.

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

- **Pairs:** k = 10 cosine neighbours on the raw L2-normalised final z within the evaluated split; each neighbour edge kept once; NK-T pair = one primary-key NK and one primary-key T cell. Flag: sufficiency-flagged pair: cosine >= flag_cosine (TEDDY puts the two cells closer than the median val NK-T neighbour pair), with flag cosine 0.949402 from 171 val NK–T neighbour pairs. Variants: all pairs, and no_gdT158.
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

Global 20260930; bootstrap 1; classifier 0, 1, 2, 3, 4; label-cost draws 0, 1, 2, 3, 4; E2 subset 7; E2 thinning 11, 12; E3 token subset 13; E4 noise 17; E5 cells 19; E5 random directions 23; gate thresholds 0.

---

## 12. Leakage check

`bridge_anm/v3_leakage_check.py` runs three checks, and all must pass:

1. **Static.** The builder reads the protein matrix, cell types and embedding only on non-test rows, and no earlier site4 result file.
2. **Invariance.** The registration core is rebuilt from inputs whose site4 protein, cell types and embedding are replaced by random values. It must be byte-identical to the core from the real inputs and to the core whose sha256 this registration records (`provenance.core_sha256` = `e30451e7bb4385fa8e0b9a3e2ed44fb929be74dea2eba4b19c0b3ec92ff233f2`).
3. **Positive control.** The same poisoning applied to the val rows must change the core.

Result: `registration/leakage_check_report.json` and `outputs/v3/registration/LEAKAGE_CHECK.md`.

---

## 13. Disclosures

- **Site4 is not virgin.** v1/v2 experiments, the 2026-09-28 audit and an exploratory answer-key prototype (scratchpad, evaluated on site4) all read site4. v3 reuses none of their thresholds or panels. Its design is still informed by what they showed (for example, that gdT CD158b+ cells look NK-like and that the v2 CD16-anchored "myeloid" class captured many NK cells). E6 exists for this reason.
- **The frozen head carries two inherited facts:**
  - its checkpoint was selected on val (`val_fm_pearson`);
  - its training size factors were divided by the median ADT total over all 90,261 cells, site4 included. That is one global scalar; v3's evidence cancels any constant scale, and its effect on the trained weights is a near-constant rescaling that is not testable without retraining.
- **Val is one donor at site1**, where the gate agrees best with the annotation. Every val-selected choice may favour site1-like data; the training-batch diagnostics (4.3) show how the gate degrades at site3.
- **Site4 cell-type counts were looked at.** As the task asked, the cell types per split were inspected, site4 included (counts only). They informed only descriptive notes such as "none in site4" for rare types, never a threshold, panel, bar or key rule (the leakage check covers every registered number).
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

## 13b. Amendment A1

Written after this registration was frozen and before any site4 evaluation; see `registration/AMENDMENT_A1.md` for each finding, its fix, the train/val numbers A1 adds and the review notes. `bridge_anm/v3_build_amendment.py` builds it from train and val rows only, and `bridge_anm/v3_leakage_check.py` checks its computed core the same way as this registration's (site4-poisoned rebuild byte-identical; val-poisoned rebuilds change it).

---

## 14. Reproduce

```bash
cd <repo>
PY=/private/tmp/claude-501/-Users-tianchichen-Documents-GitHub-teddy-mm/0828c10d-3512-46ba-b296-fafd23e4a760/scratchpad/venv312/bin/python
ANM_ROOT=<ANM v2 fix checkout> $PY bridge_anm/v3_build_registration.py      # writes registration_v3.json + .sha256 (train/val only)
$PY bridge_anm/v3_render_registration.py                                     # this document
ANM_ROOT=<...> $PY -m pytest -q tests/test_v3_key.py
ANM_ROOT=<...> $PY bridge_anm/v3_build_amendment.py                          # amendment A1 (train/val only)
ANM_ROOT=<...> $PY -m pytest -q tests/test_v3_amend.py
$PY bridge_anm/v3_leakage_check.py                                           # proof of no site4 use (registration and A1)
```

Builder `v3_build_registration 1.0`, sha256 `bc808b808ab5a9d765ea8ee42a1866526b8d36bdc00e32acb9a0b679e1a02a53`; `v3_key.py` sha256 `45113585c4c2411a966814154ec41a58c0836c51e1f3f16e0fbe31b2e9424602`; inputs: `cite_arrays.npz` sha256 `be408eff537eeeca24129210592605f809724c3282d7c55e008b87b3c193c2da`, `z_rna.npy` sha256 `fb9aa8c77d544d84f716cf66bca8964c265f62d6dc1618a0c627d55268d3cd0b`, head `best.pt` sha256 `b0c543e2851232051552ae80269b7a1791350f80883b9bf3a17b565585c1de15`; git HEAD at build `6e227b037ac9bca56c62d66b82b03843883a2079`.
