# teddy_mm v3: registered results and claims ledger

**Status (2026-10-01).** E1, E2, E3, E4 and E5 are final. **E5-M and E6 are pending.** Their results are not part of this version, and nothing below reports an E5-M or E6 outcome.

**Sources.** Results: `outputs/v3/{E1,E2,E3,E4,E5}/REPORT.md` and `*_results.json`. Registration: `registration/REGISTRATION_v3.md` (`registration_v3.json`, sha256 `e4c8a33e…`), amendments `AMENDMENT_A1.md`, `AMENDMENT_A2.md`, `AMENDMENT_A3.md`, addenda `registration/addenda/*.json`. Earlier (v1/v2) results: `README.md` and `docs/reports/*_OFFICIAL.md`. Every number used here is in [`v3_numbers.json`](v3_numbers.json), each with its source file and key (section [Numbers](#numbers)).

**Precision.** Numbers are printed at the precision of the experiment's own `REPORT.md`: 4 decimals for key validity, E3 and E4; 3 decimals for E1, E2 and E5. Where a number appears only in a results JSON, I round it to the same precision and say so. Full-precision values are in `v3_numbers.json`.

> **中文摘要：** 本文汇总 v3 预注册实验 E1–E5 的最终结果（E5-M、E6 待定，本版未报告任何结果），列出能说、不能说和未检验的结论；所有数字见 `v3_numbers.json`。

---

## Plain names used here

| name used here | what it is |
|---|---|
| **TEDDY embedding** (the gene-mean state, written z_12 in E5-M) | The official TEDDY-G embedding of a cell's RNA: the mean over gene tokens of TEDDY's last-layer (layer 12) token states, 512 numbers, L2-normalised. Our head reads it. TEDDY is never retrained. |
| **head** | The phase-1 head trained in this repo: an MLP with a negative-binomial decoder that predicts 134 surface proteins from the TEDDY embedding (`outputs/cite_phase1_official/best.pt`), frozen. |
| **evidence** | For each panel protein: the head's prediction divided by the 95th percentile of the head's own predictions over the training cells, clipped to [0, 1]. No measured protein enters it. |
| **declared rule** (TEDDY + fixed rule) | For each class, the mean evidence over its panel; call the highest class if its score reaches the question's bar, otherwise "no call". |
| **TEDDY + ANM** | The same evidence fed to ANM's engine (`finite_graph_scalar`, ANM's default field: 4 steps, retention 0.82, diffusion 0.16), one event per panel protein, **all events at t = 0**. Used only where ANM's field actually runs. |
| **decision layer (Mode B)** | ANM used on top of TEDDY's outputs, as here in E1, E2 and E4. **Looking inside TEDDY (Mode A)** means studying TEDDY's layers (E5, E5-M). |
| **OUT** | Out of scope: cells that are none of B, T, NK or myeloid (erythroid cells, progenitors, plasma cells, pDC, ILC). The correct action is no call. |
| **primary key** | A cell's class when the annotated cell type and the protein gate (measured CD3, CD19, CD56, …) agree; otherwise the cell is **unscored**. **Annotation-only key**: the annotated class for every cell. |
| **coverage** / **selective accuracy** / **decision accuracy** | Share of all cells that get a call / share of called, scored cells whose call is right (a call on a key-OUT cell is wrong) / share of scored cells with a right call or a no call on a key-OUT cell. |
| **OUT decline rate** | Share of key-OUT cells that get no call. |
| **matched coverage**, **AURC** | Every method calls its own most confident ⌈c·N⌉ cells, so accuracy is compared at equal coverage; AURC is the area under selective accuracy over coverage 0.95 … 0.50, divided by 0.45. |
| **workability (P_f)** / **exactness (Q_f)** | P_f: share of cells that get a call. Q_f: share of scored cells whose lineage call equals the key (as the E1 code computes it, a key-OUT cell never counts as correct here). |
| **verdict words** | *win*: the difference reaches the margin, the 95% interval's lower bound is above 0, and in each primary donor the difference has the same sign and is at least half the margin. *loss*: the mirror. *equivalent*: the interval lies inside ±margin. Otherwise *inconclusive*. Intervals: two-stage bootstrap (donors, then cells within donors), B = 2000, seed 1. |

> **中文：** 术语表：TEDDY 嵌入、头部、证据、声明规则、TEDDY + ANM（所有事件在 t = 0）、OUT（范围外）、主答案键与仅注释键、覆盖率与选择准确率、判定词（win/loss/equivalent/inconclusive）。

---

## 0. Bottom line

- **ANM equals a declared rule in this setup, and the registered tests found no decision value from ANM's readouts.** ANM's engine made the same call as the declared rule on every test cell for all three questions, and its closure readout matched its re-coded closed form on every cell (E1.1a: 0 mismatches). E1's falsification clause holds: neither the trust-gate test (E1.4) nor the nested-readout test (C2) is a win.
- **Writing a new question with zero labels did not beat training.** On the changed question Q3, a classifier trained only on the old question's labels was more accurate than the declared rule (rule minus classifier −0.078 [−0.134, −0.028], loss), and the label-cost curve gives n* = 0.
- **E2:** splitting perturbation losses into representation, head and decision shares added nothing beyond the plain accuracy curve (falsified). **E3:** the claim that the NK–T loss sits in the readout or pooling and can be repaired on frozen TEDDY is rejected; both trained readouts lost NK-vs-T accuracy to our head. **E4:** the ANM fusion arm lost to a stacker trained on the validation donor (registered outcome label "field hurts"); against the trust-weighted average, which is ANM without its contradiction events, it is inconclusive. **E5:** the registered derivative check failed in float32 (it passes in float64), so E5 stops with "JVP not validated".
- **The answer key is weak for one of the two test donors** (donor 13272: annotation-vs-gate kappa 0.4177). Annotation-only results are therefore reported with equal prominence. In E1 every annotation-only verdict matches the primary-key verdict. In E4 the alternative key reported with equal prominence is the v3 primary key. It gives the same verdicts for the two registered comparisons, but not for every secondary one (see E4).

> **中文：** 在本设置中 ANM 与声明规则逐细胞相同，预注册检验未发现 ANM 读出带来决策价值；零标签改问题不胜训练分类器；E2 被证伪，E3 被拒绝，E4 不敌验证集堆叠器，E5 的 JVP 检验在 float32 下未通过。

---

## 1. Why v3

The v1/v2 experiments (README.md, `docs/reports/*_OFFICIAL.md`) had problems that made their headline claims untestable or invalid. v3 fixes each one before looking at the test site again.

| problem in v1/v2 | evidence (source) | v3 fix |
|---|---|---|
| **The answer key had only three classes (B, T, myeloid): no NK class and no out-of-scope class.** NK cells, which carry CD16, were scored as "myeloid", and there was no class for which "no call" is the right answer. | 1,342 of 1,690 NK cells were keyed "myeloid" under the soft rule (README.md line 84). Out-of-scope cells still got a call 73–95% of the time (line 53). | Five classes: B, T, NK, myeloid and OUT, where no call is the correct answer for OUT. The primary key keeps only cells where the annotation and a protein gate agree; all others are unscored, and the annotation-only key is always reported too. |
| **The B/T-priority question was synthetic.** Its CD16-anchored "myeloid" class captured many NK cells. | Registration section 7, audit finding KEY-1. The question changed 12.48% of expected answers (README.md line 184). | Q3 is the standard four-population question, anchored on CD19 (B), CD3 (T), CD56 (NK) and CD14 (classical monocyte). Biology fixed it; nothing was selected. |
| **The panel was chosen on the test set.** | The 9 panel proteins were chosen by test-set (site4) Pearson (README.md lines 60, 170). | Panels are selected on the validation donor only, by worst-pair AUROC. The old 9-marker panel survives only as a sensitivity row. |
| **The test set was not donor-independent.** | Donor 15078 is in training at sites 1–3 and is 32.6% of the site4 cells. The "out-of-site check" was the validation donor at a training site (README.md line 60). | Primary evaluation uses only the site4 donors that appear in no other split (13272, 19593). Donor 15078 is a secondary split that never changes a verdict. |
| **Measured protein scales sat in the evidence path.** | The first run multiplied every prediction by a per-cell factor computed from measured protein (the answer key). v2 used one train-median factor, 0.92663, and divided by the training 95th percentile of *measured* protein (README.md line 47; registration section 5). | Evidence = prediction ÷ the training 95th percentile of the head's *own* predictions. A constant size factor cancels, so no measured-protein array enters the evidence. |
| **Event order acted as a hidden weight, and the first run's ANM "edge" came from declining more.** | The first run fed markers to ANM one per step in panel order. On B/T-priority, ANM declined 3,373 cells to the rule's 1,796, scoring 0.9603 vs 0.9428; at ANM's coverage the rule scored 0.9657 (README.md lines 46, 49). | All events enter at t = 0, and every accuracy comparison is at matched coverage. |
| **No registered margins, decision rules or donor-level uncertainty.** Experiments 2 and 6 were withdrawn. "0 labels for a new question" holds for any written rule (README.md lines 54–58). | — | Pre-registered endpoints, margins and win/loss/equivalent rules; a two-stage donor-then-cell bootstrap; failures reported. "Zero labels" is credited to any written rule, never to ANM. |

> **中文：** v1/v2 的问题包括：答案键无 NK 和范围外类别、B/T 优先问题是人为设计、面板在测试集上选择、测试供体不独立、证据路径含测量蛋白尺度、事件顺序成为隐藏权重、缺少预注册判定规则；v3 逐一修正。

---

## 2. How v3 was registered

### 2.1 Splits

| split | rule | cells |
|---|---|---:|
| train | sites 1–3, every donor except 18303 | 67,405 |
| val | donor 18303 at site1 | 6,106 |
| test_primary | site4 donors in no other split: 13272 (7,365) and 19593 (3,929) | 11,294 |
| test_secondary | site4 cells of donor 15078 (in training at sites 1–3) | 5,456 |

Every threshold, panel, bar, key rule and analysis choice was fixed on train and val only. The validation set is a single donor at site1, the site where the gate agrees best with the annotation, which limits every val-selected choice (registration section 13).

### 2.2 Answer key

- **Annotation class.** Each of the 45 annotated cell types is mapped to B, T, NK, myeloid or OUT (registration section 4.1). Borderline types are explained there: gdT CD158b+ → T, ILC/ILC1 → OUT, pDC → OUT, plasma cells → OUT, progenitors and erythroid cells → OUT.
- **Protein-gated class.** A declared gate on measured protein, applied in this order: OUT if CD71 is high; T if CD3 is high, CD19 is not high and CD14 is not high; B if CD19 is high, CD3 is not high and CD14 is not high; NK if at least 2 of CD56, CD94, CD335 and CD16 are high and CD3, CD19, CD14 and CD33 are not; myeloid if CD14, CD11c or CD33 is high and CD3 and CD19 are not; otherwise OUT. Thresholds come from a 2-component Gaussian mixture fitted on the training cells. The variant was chosen from a 12-variant grid by mean training-batch kappa (0.8225, minimum 0.7091) and confirmed on val (kappa 0.8905, registration sections 4.2–4.3).
- **Primary key** = the class where annotation and gate agree; otherwise unscored. **Registered key-validity flag:** if site4 kappa < 0.85, every result is also reported, with equal prominence, on the annotation-only key. The primary key is not switched, because switching would use site4 labels to make a choice.

### 2.3 Panels, questions, methods

- **Panels**, selected on val. Each protein is assigned to the class where its worst-pair AUROC (the class against its hardest other class) is highest, and a class panel is its top 3 assigned proteins: B CD20, CD22, CD268; T CD3, CD2, CD5; NK CD122, CD94, CD56; myeloid CD172a, CD11c, CD62P. CD62P is a platelet protein. It was kept because the declared rule chose it.
- **Questions and bars.** Each bar is the val quantile that matches a declared val no-call rate:
  - **Q1** (soft 4-class): bar 0.298614, no-call rate 0.15.
  - **Q2** (strict 4-class, nested in Q1): bar 0.661255, no-call rate 0.30.
  - **Q3** (CD14-anchored B / T / NK / classical monocyte): bar 0.273973, no-call rate 0.15.
- **Methods**, all reading the same evidence of the same cells:
  - The declared rule.
  - TEDDY + ANM. With equal panel sizes and all events at t = 0, ANM's action score is G(3)·3·S, where S is the rule's class score and G(3) = 0.370041 is the field's gain for a class with 3 events. So ANM's calls equal the rule's by construction, and E1 checks this cell by cell.
  - ANM's closure readout (C2 only), checked against its re-coded closed form.
  - A trained classifier: multinomial logistic regression on the 12 evidence values, trained on 58,657 labelled training cells, C = 100 chosen by val log-loss.
  - TEDDY's margin: top score minus second score.
- **Statistics.** Matched coverage with a seeded tie-break (A1.1). Two-stage bootstrap, B = 2000, seed 1, with every endpoint recomputed inside each replicate. Per-donor rows. With 2 primary donors, the intervals describe those donors, not a population.
- **Leakage check.** The registration core rebuilt from inputs whose site4 protein, cell types and embedding were replaced by random values is byte-identical to the real core. Poisoning the val rows changes it (`registration/leakage_check_report.json`). A1, A2, A3 and the E1, E2, E3, E4 and E5-M addenda report the same kind of check as passed.

### 2.4 Amendments A1–A3

All times are 2026-10-01 UTC.

| amendment | what it fixed | built | committed (local git) | pushed to GitHub |
|---|---|---|---|---|
| **A1** (sha256 `1f44be66…`) | Ten problems found by an adversarial review before any site4 evaluation. **A1.1:** tied scores broken by a fixed random permutation instead of row order (val: 0.3662 of cells have a Q1 top score of exactly 1), with matched coverage recomputed per replicate. **A1.2:** E4 channel 2 reads a panel-free renormalisation (the stored per-cell centred log-ratio (CLR) transform mixed panel counts into the non-panel inputs). **A1.3:** exact leave-one-out closed form. The registered form re-coded zeroing an event, while ANM deletes the event's site; ρ = G(2)/G(3) = 0.992479, and 6 val calls change flip status. Top markers become sets. **A1.4:** honest naming in E1.4: ANM's soft score is G(3)·3 × the rule's top score, so an E1.4 win belongs to the top-score readout. **A1.5:** the C2 verdict point is c* (the rule's realised coverage), and "in-scope accuracy" is defined. **A1.6:** information parity in E1.1c/E1.5 (the classifier gets a superset of the rule's 4 anchors). **A1.7:** E2 criterion (i) needs ≥ 30 lost cells and an interval excluding 0. **A1.8:** E3 pairs built within each donor. **A1.9:** the E4 comparator is fixed on val. **A1.10:** an addendum cannot replace a registered endpoint. | 04:25:39 | first version `7c49b7a` 04:24:51; final `ca1c8ec` 04:37:04 | 07:57:00 |
| **A2** (sha256 `20b64174…`) | **A2.1:** the E3 statistic D_R is not scale-free. A readout that shrinks every NK–T difference by a common factor can pass it while repairing nothing, so a win now also needs (b) recovery on flagged pairs and (c) a positive scale-free Dlog. **A2.2:** R2's normaliser is estimated on a seeded 10,000-cell training sample (0.33 h of compute instead of 2.20 h; same estimand). **A2.3:** the per-donor clause of the common win rule is binding for E2's criteria. All three make verdicts stricter or keep them unchanged. | 06:09:15 | `62ba605` 06:11:37 | 07:57:00 |
| **A3** (sha256 `4bb2f856…`) | **A3.1:** the E3 pooling contrast D_R2 − D_R1 has the same scale artefact. A pooling win now also needs (b') flagged-pair recovery beyond R1 and (c') a positive scale-free contrast. Training-pair illustration: a pure rescaling reaches the margin at 1 − s = 0.379949. Stricter only. | 07:00:56 | `36d3584` 07:03:09 | 07:57:00 |

### 2.5 Fairness facts and timeline (UTC)

Local git commit times come from the committing machine. The GitHub push events (`repos/danielchen26/teddy_mm/activity`, branch `exp/registered-v3`) are the independent timestamps. A site4 evaluation here means the first stage that reads site4 protein or cell types; forward passes of site4 RNA alone are label-free.

| file it depends on | committed | pushed to GitHub | site4 evaluation that uses it | pushed before that evaluation? |
|---|---|---|---|---|
| registration_v3.json | 03:54:53 | 04:24:22 (branch created) | first v3 site4 read: E5 cell selection 05:10:36 | **yes** |
| A1 (final) | 04:37:04 | 07:57:00 | E5 05:10:36; E1 05:36:47; E3 08:42:23; E4 09:04:11; E2 10:41:58 | **no for E5 and E1** (committed before, pushed after); yes for E2, E3, E4 |
| E5 addendum (v3, `0c4e028b`) | 04:50:00 | 07:57:00 | E5 05:10:36 (results written 10:42:46) | **no** (committed before; pushed while E5 was running). Its first version (`820cce0`, 04:19:49) was in the 04:24:22 push, but A1.10 made the registered endpoint and cells primary, which the v3 version records |
| E2 addendum | 05:23:26 | 07:57:00 | E2 report 10:41:58 | yes |
| E3 addendum | 05:25:05 | 07:57:00 | E3 evaluate 08:42:23 | yes |
| E1 addendum | 05:30:59 | 07:57:00 | E1 05:36:47–05:44:44 | **no** (committed about 5 minutes before E1 started; pushed 2 h 12 min after it finished) |
| A2 | 06:11:37 | 07:57:00 | E3 08:42:23; E2 10:41:58 | yes |
| E4 addendum | 06:44:49 | 07:57:00 | E4 evaluate 09:04:11 | yes |
| A3 | 07:03:09 | 07:57:00 | E3 08:42:23 | yes |
| E5-M addendum | 07:39:16 | 07:57:00 | E5-M run started 10:43:09 | yes |
| E6 addendum v1 / v2 | 10:07:39 / 10:43:00 | 10:43:44 | E6 external run: none at the time of writing | yes |

So:

- **Pushed before evaluation.** The registration, A2, A3 and the E2, E3, E4 and E5-M addenda reached GitHub before the site4 evaluation they govern. The E6 addenda reached GitHub before any external outcome. The addendum records `outputs/v3/E6` absent at 2026-10-01T10:40:07Z.
- **Committed, not yet pushed.** For **E1** and **E5**, the governing files (A1, the E1 addendum, the E5 addendum v3) were committed locally before the evaluation, and each result file records their sha256 and `committed: true`. They were pushed to GitHub only afterwards (07:57:00). For these two experiments, ordering rests on local commit times and recorded hashes, not on an independent GitHub timestamp.
- **E1 came before A2 and A3.** E1 finished at 05:44:44, before A2 (built 06:09:15) and A3 (built 07:00:56). Both concern only E2 and E3. A3 states that the E1 site4 results had been seen and that nothing in A3 depends on them. A2 does not say whether they had been seen. The E5-M and E6 addenda record the same disclosure (E1, and for E6 also E3 and E4, had been seen).
- **E3's site4 pairs were built early, by E5.** E5's cell selection at 05:10:36 rebuilt E3's within-donor NK–T pairs on site4 from the primary key, so it read site4 protein and cell types. Its log records the pair count (253) and the cells per donor and class. This happened before the E3 addendum (committed 05:25:05), the E1 addendum, A2 and A3. No gap ratio, D statistic or NK-vs-T accuracy was computed then: the rest of E5's log until its results (10:42:46) holds train/val direction diagnostics and per-cell progress lines. A3 describes E3 evaluate as the only *E3* stage that reads site4 protein or cell types, which is true of E3's own stages, but it does not mention this E5 step.
- **Key validity.** On test_primary the annotation-vs-gate kappa is 0.5781 (agreement 0.6847), below the registered 0.85, so the flag is raised. By donor: 13272 kappa 0.4177 (agreement 0.5617) and 19593 kappa 0.8798. In donor 13272 the gate puts 4,148 of 7,365 cells in OUT, against 1,707 annotated OUT cells. It recovers 0.443 of annotated T cells and 0.395 of annotated B cells (3 decimals from the JSON). A per-batch gate would raise kappa to 0.7019, but it is a sensitivity key only, because it reads site4 protein to set thresholds. 3,561 of 11,294 test_primary cells are unscored on the primary key. **Annotation-only rows are reported beside every primary-key row below.**
- **Other disclosures** (registration section 13):
  - Site4 is not virgin: v1/v2 and an exploratory key prototype read it, and v3's design is informed by what they showed.
  - The head's checkpoint was selected on val.
  - The head's training size factors were divided by one global ADT-median scalar over all cells, site4 included. v3's evidence cancels constant scales; the effect on the head's weights cannot be tested without retraining.
  - Pretraining contamination of TEDDY with this dataset (GSE194122) is not verified.

> **中文：** 预注册、A2、A3 及 E2/E3/E4/E5-M 附录在相应 site4 评估之前已推送到 GitHub；但 E1 和 E5 所依赖的 A1、E1 附录、E5 附录只是在评估前本地提交，推送晚于评估。E1 在 A2/A3 之前完成；E5 在 E3 附录、A2、A3 之前已在 site4 上重建了 E3 的配对（仅配对数，未计算 E3 统计量）。供体 13272 的答案键 kappa 仅 0.4177，故仅注释键结果同等呈现。

---

## 3. Experiments

### E1: decision layer rerun on the v3 key (Experiments 1, 3, 4, 5 and C2)

Cells: test_primary (11,294), with test_secondary reported under the same rules. ANM runs through its engine (`finite_field_runner`, events `<lineage>:<protein>` at t = 0).

**E1.1a: does ANM equal the declared rule?** The registered value is 0 mismatches.
- **Result:** call mismatches against the rule Q1 0, Q2 0, Q3 0; closure readout against its re-coded closed form Q1 0, Q2 0; 0 rejected events. Same on test_secondary. **Pass.** This is a property of the setup, not a finding.

**E1.1 (change the question).**
- **Question:** if the question changes from Q1 to Q3, who answers Q3 better at the same coverage: the declared rule, which needs 0 labels, or a classifier?
- **Arms:**
  - The rule on the 4 anchors.
  - E1.1c: a classifier trained on Q1 labels only, with no relabelling, on 14 features that include the 4 anchors.
  - E1.1d: a classifier trained on Q3 labels, the label-rich ceiling.
  - Secondary rows: both classifiers on the 4 anchors only.
- **Margin:** 0.01, at the rule's realised Q3 coverage of 0.935 [0.926, 0.949].
- **Registered verdicts:** E1.1c **loss**, E1.1d **loss**.

| key | arm | selective accuracy [95% CI] | per donor (13272, 19593) | rule minus arm [95% CI] | per donor | verdict |
|---|---|---|---|---|---|---|
| q3 | rule | 0.692 [0.638, 0.751] | 0.645, 0.743 | | | |
| q3 | Q1-label classifier (E1.1c) | 0.770 [0.668, 0.878] | 0.675, 0.872 | −0.078 [−0.134, −0.028] | −0.031, −0.129 | loss |
| q3 | Q3-label classifier (E1.1d) | 0.770 [0.664, 0.880] | 0.671, 0.875 | −0.078 [−0.137, −0.024] | −0.026, −0.132 | loss |
| q3 annotation-only | rule | 0.723 [0.715, 0.733] | | | | |
| q3 annotation-only | Q1-label classifier | 0.777 [0.737, 0.847] | | −0.054 [−0.121, −0.018] | −0.020, −0.116 | loss |
| q3 annotation-only | Q3-label classifier | 0.776 [0.736, 0.847] | | −0.053 [−0.121, −0.018] | −0.019, −0.116 | loss |

*Secondary:*
- The 4-anchor-only classifiers also win against the rule (−0.078 [−0.135, −0.028] and −0.077 [−0.132, −0.028], both loss for the rule).
- E1.1b: changing Q1 → Q2 changes 0.249 of calls, and Q1 → Q3 changes 0.127.
- test_secondary: E1.1c loss on the q3 key (−0.020); inconclusive on annotation-only.

**E1.3 (why this call).**
- **Question:** does ANM's leave-one-out attribution (delete each event of the called class, re-evolve the field) agree with the exact closed form, and do flip-sensitive calls flag errors better than the same number of lowest-margin calls?
- **Margin:** 0.02. Q1 is primary; Q2 is secondary.
- **Registered verdicts:** E1.3a discrepancy **False**; E1.3c **inconclusive**, with the registered reading "attribution adds no flagging value beyond TEDDY's margin".
- **E1.3a:** on 10,715 Q1 calls, top-marker-set agreement is 1.0000 and flip agreement is 1.0000. 0.331 of calls have a multi-marker top set. Flip-sensitive calls: 3,557 (engine 3,557).
- **E1.3c:** error of flip-sensitive calls 0.625 [0.473, 0.755] vs as many lowest-margin calls 0.621 [0.455, 0.760]. Difference 0.004 [−0.008, 0.022], per donor −0.005 and 0.018: **inconclusive**. Annotation-only: 0.008 [−0.000, 0.023], per donor 0.002 and 0.019: inconclusive.
- *Secondary:*
  - Q2: −0.009 [−0.015, 0.003], **equivalent**.
  - E1.3b deciding-marker shares within the called class, both donors pooled, against a permutation null with mean 1/3 (p = 0.0010 for each): CD20 decides 0.484 of B calls, CD2 0.618 of T, CD94 0.572 of NK and CD11c 0.832 of myeloid. The platelet protein CD62P decides only 0.043 of myeloid calls.
  - The B share does not hold by donor: in donor 19593 CD20 decides 0.339 of B calls (p = 0.6494, at the null) and CD22 decides 0.491.
  - These shares describe the panel evidence. Any written rule reproduces them.

**E1.4 (which calls to trust).**
- **Question:** does ranking calls by ANM's top score (= G(3)·3 × the rule's top score, A1.4) order them better than TEDDY's margin?
- **Primary endpoint:** AURC; margin 0.005. Q1 and Q2 share scores and key, so E1.4 is computed once.
- **Registered verdict:** E1.4a **inconclusive**. The registered failure rule, as A1.4 words it: "if the top score (rule = ANM) is not a win against the margin, the gate claim stays withdrawn".

| key | top score AURC | margin AURC | entropy AURC | classifier AURC | top minus margin [95% CI] | per donor | verdict |
|---|---|---|---|---|---|---|---|
| primary | 0.843 [0.818, 0.867] | 0.844 [0.835, 0.854] | 0.847 [0.838, 0.860] | 0.942 [0.894, 0.983] | −0.002 [−0.024, 0.020] | 0.018, −0.021 | inconclusive |
| annotation-only | 0.836 | 0.834 | 0.839 | 0.923 | 0.002 [−0.009, 0.009] | 0.008, −0.007 | inconclusive |

*Secondary* (never decides; any win belongs to the rule's top score, not to ANM):
- Selective accuracy at coverage 0.90, top vs margin: +0.021 [0.017, 0.025], per donor 0.022 and 0.020 → win. At 0.70: −0.002 [−0.038, 0.030] → inconclusive.
- Top score vs classifier AURC: −0.099 [−0.163, −0.034], per donor −0.036 and −0.157 → **loss**. A logistic classifier on the same 12 evidence values orders calls better.
- Sensitivity panels (the old canonical panel and the val k = 2 and k = 5 panels): all inconclusive.

**E1.5 (label cost).**
- **Question:** how many Q3-labelled training cells does the classifier need to reach the rule's Q3 selective accuracy minus 0.005?
- This is a measurement, not a contest; the rule needs 0 labels.
- **Registered result:** n* = **0** on the q3 key and on the annotation-only key. The classifier trained only on Q1 labels (n = 0) already scores 0.770 against the rule's 0.692. n = 25 scores 0.762 and n = all scores 0.770.
- *Secondary:* the 4-anchor curve also gives n* = 0. The Q1 reference (primary key) gives n* = 25.

**C2 (nested readouts).**
- **Question:** at the same coverage, does ANM's closure readout decline more out-of-scope cells than the mean rule, without losing in-scope accuracy?
- **Margin:** 0.05. Guard: in-scope accuracy difference ≥ −0.005. Verdict point: c* = 0.949 [0.942, 0.953].
- **Registered verdict:** **equivalent (in-scope guard also failed)**.

| key | coverage | OUT decline: closure / mean rule / classifier | closure minus rule [95% CI] | per donor | in-scope accuracy difference | verdict |
|---|---|---|---|---|---|---|
| primary | c* | 0.261 / 0.261 / 0.273 | 0.000 [−0.006, 0.003] | 0.001, −0.002 | −0.056 [−0.085, −0.025] | equivalent (guard failed) |
| annotation-only | c* | 0.237 / 0.238 / 0.242 | −0.001 [−0.007, 0.001] | 0.001, −0.005 | −0.041 | equivalent (guard failed) |
| primary (secondary point) | 0.80 | 0.906 / 0.855 / 0.872 | 0.051 [0.011, 0.062] | 0.021, 0.055 | −0.054 | inconclusive (guard failed) |

*Secondary:* the nested pair at the deployed bars. For the rule, every Q2 call is nested in a Q1 call. Workability (P_f) is 0.949 (Q1) and 0.699 (Q2); exactness (Q_f) is 0.649 (Q1) and 0.582 (Q2).

**Deployed operating points (bars from val, test_primary, primary key; descriptive).**
- Q1 rule: coverage 0.949, selective accuracy 0.697 (annotation-only 0.730, 3 decimals from the JSON), OUT decline 0.261. 3,561 test_primary cells are unscored on the primary key, and the Q1 rule calls 3,512 of them (calls on unscored cells).
- ANM closure Q1: 0.947 / 0.654 / 0.271.
- Classifier Q1: 0.968 / 0.757 / 0.172.
- Q2 rule: 0.699 / 0.867 / 0.970. Classifier Q2: 0.728 / 0.969 / 0.942.
- Q3 rule: 0.935 / 0.692 / 0.287.
- *Secondary, per annotated type (Q1 rule):* 385 of 395 annotated gdT CD158b+ cells are called NK. Many annotated erythroid and progenitor cells are called NK as well, for example 467 of 513 proerythroblasts.

**E1 falsification (registered, A1 wording):** "falsified for v3 (A1 wording): neither E1.4 nor C2 is a win, so the readout forms ANM provides (top action score, closure readout) add no decision value over TEDDY's margin and the mean rule. Every ANM arm equals a re-coded rule cell by cell (E1.1a)."

**test_secondary (donor 15078; never changes a verdict).**
- Key kappa 0.8674; E1.1a passes.
- E1.4a top vs margin −0.004 [−0.006, −0.002]: inconclusive.
- C2 at c*: equivalent (guard failed). E1.5 n* = 0.

**Limitations.**
- Two primary donors, and the key is weak in one of them (13272).
- The primary key keeps cells both sources call cleanly, which are easier than average.
- Bars come from one val donor at site1. Site4 coverage at the Q1 bar (0.949) is above the val target (0.85).
- The classifier is a stand-in. TypeSafe AI's Jev was never called.

> **中文：** E1：ANM 与规则逐细胞一致；改为 Q3 后，仅用 Q1 标签的分类器已优于零标签规则（loss，n* = 0）；归因、信任排序、嵌套读出均未取胜，E1 的证伪条款成立。

### E2: response decomposition (C3) and linearity (C5)

- **Question:** when a cell's RNA is degraded (binomial thinning of UMIs) or one gene is scaled down, does splitting each lost call into "representation", "head" or "decision" tell us more than the plain perturbation → accuracy curve?
- **Design:**
  - Cells: a label-free seeded sample of 500 cells per primary donor (1,000 cells, 23,116 cases) plus 250 cells of donor 15078. Official TEDDY-G preprocessing → frozen head → evidence → Q1 rule and ANM.
  - Gene scaling covers the coding genes of the 12 panel proteins plus the registration's 13 NK/T genes. The two lists overlap (for example CD3E, NCAM1, KLRD1). SELP is absent from the data and was dropped, which leaves 20 distinct genes.
  - A lost cell is one that was called right at baseline and is not called right after the perturbation. Each lost cell gets a label in this order:
    - *representation*, if a z-probe (logistic regression on training embeddings, C = 100 chosen on val) is wrong on the perturbed embedding;
    - else *head*, if the rule's argmax on the perturbed evidence is wrong;
    - else *decision* (the right class falls below the bar).
- **Comparator:** the plain curve. **Falsified unless** (i) two perturbations with losses within 0.02 of each other (each with ≥ 30 lost cells) differ by ≥ 0.15 in a localisation share, or (ii) the small-ε response predicts which cells change call at ε = 0.8 with AUROC ≥ 0.02 above the baseline margin. Both criteria use the common win rule including the per-donor clause (A2.3).
- **Registered verdict:** **"adds nothing beyond the curve (falsified)"**, on test_primary and on test_secondary.
- **Key numbers (test_primary):**
  - Baseline Q1 coverage 0.945, selective accuracy 0.695, decision accuracy 0.716.
  - Checks: ANM engine vs rule, 0 mismatches over all cases. The unperturbed re-embedding changes no call.

| perturbation | decision-accuracy loss [95% CI] | per donor (13272, 19593) | lost cells | representation / head / decision share |
|---|---|---|---:|---|
| 80% of UMIs kept | 0.022 [0.004, 0.042] | | 25 | 0.000 / 1.000 / 0.000 |
| 20% kept | 0.162 [0.121, 0.194] | 0.140, 0.174 | 112 | 0.069 [0.00, 0.11] / 0.931 / 0.000 |
| 5% kept | 0.140 [0.112, 0.172] | | 125 | 0.281 [0.21, 0.41] / 0.719 / 0.000 |
| one gene masked (ε = 1) | 0.012 [0.001, 0.024] | | 20 | 0.000 / 1.000 / 0.000 |

  - Criterion (i): 0 perturbation pairs qualify.
  - Criterion (ii): 2,000 cases, 398 change call at ε = 0.8. AUROC: baseline margin 0.854 vs extrapolated small-ε margin 0.763. Difference −0.092 [−0.125, −0.057], per donor −0.090 and −0.090 → **loss**.
- *Secondary (C5 linearity, a registered measurement with no decision rule):*
  - The share of cells whose response is linear in ε (|ρ − 1| ≤ 0.1, cosine ≥ 0.9, above 3 × floor) is at most 0.098 in any row, at any stage (embedding, evidence, score).
  - At the smallest thinning pair (ε 0.1 vs 0.05), the embedding's median slope ratio ρ is 0.667, with median cosine 0.564 between the two responses.
  - At these sizes the pipeline's response is not linear in the perturbation.
- **Limitations.**
  - The z-probe is near-perfect on val (accuracy 0.9977, 4 decimals from the JSON), so "representation" rarely fires. "Head" therefore collects most losses by construction of the label order.
  - "Decision" is 0.000 at every level: no lost call fell below the bar. Coverage rose under thinning (0.961 at 80% kept and 1.000 at 20% kept, as REPORT.md prints them).

> **中文：** E2：把扰动造成的错误分解为表示/头部/决策，并不比普通的扰动-准确率曲线多提供信息（已证伪）；响应在所测尺度下不是线性的。

### E3: NK–T look-alike pairs (C4) and repair on frozen TEDDY

- **Question:**
  - NK and T cells that are near neighbours in TEDDY's embedding: does our head compress their protein differences more on these look-alike pairs?
  - Can a different readout of frozen TEDDY repair this? R1 is an MLP on the 12 layer means. R2 is attention pooling over layer-12 gene-token states.
- **Design:**
  - Pairs: k = 10 cosine neighbours within each donor. One pair = one primary-key NK cell and one primary-key T cell. A pair is *flagged* when its cosine is ≥ 0.949402 (from 171 val NK–T neighbour pairs).
  - 253 pairs (164 flagged; 278 cells): 13272 has 157 (111 flagged), 19593 has 96 (53).
  - Gap ratio GR = median |predicted difference| / median |measured difference| within pairs, averaged over CD56, CD94, CD335 and CD3.
  - H3a: D_R = [GR_R − GR_head](flagged) − [GR_R − GR_head](unflagged) − the same for a null (the head retrained with seed 1). Win ≥ 0.05, plus the A2.1 conditions (b) and (c). Pooling (R2 over R1) uses the A3.1 conditions.
  - H3b: NK-vs-T selective accuracy at coverage 0.9 and 0.8, margin 0.01.
  - E3 has no ANM arm. Its decision rule is the declared Q1 rule.
- **Registered verdicts:**
  - E3.H3a R1 **inconclusive**; R2 **inconclusive**; pooling **inconclusive**.
  - E3.H3b R1 vs head and R2 vs head at 0.9 and at 0.8: **loss** in all four.
  - Falsification (A2.1): **"rejected: neither R1 nor R2 has an A2 win for E3.H3a … so 'the NK-T loss is in the readout or pooling and is repairable on frozen TEDDY' is rejected for this dataset."**

| readout | GR on flagged pairs [95% CI] | GR on unflagged pairs [95% CI] |
|---|---|---|
| head | 0.1169 [0.0472, 0.4677] | 0.3985 [0.3105, 0.6207] |
| R1 | 0.2635 [0.1737, 0.5369] | 0.4058 [0.3481, 0.5586] |
| R2 | 0.2722 [0.2026, 0.3911] | 0.452 [0.3876, 0.6125] |
| null (head, seed 1) | 0.2104 [0.0591, 0.6995] | 0.5023 [0.3547, 0.7488] |

| statistic | point [95% CI] | per donor (13272, 19593) | verdict |
|---|---|---|---|
| D_R1 | 0.1496 [−0.1358, 0.2948] | 0.1728, −0.034 | inconclusive; A2 (b) True, (c) False |
| D_R2 | 0.1121 [−0.5087, 0.4327] | 0.311, −0.417 | inconclusive; A2 (b) False, (c) False |
| D_R2 − D_R1 (pooling) | −0.0375 [−0.4408, 0.22] | 0.1382, −0.383 | inconclusive; A3 (b') False, (c') False |

| H3b (3,861 cells: 740 NK, 3,121 T) | accuracy at 0.9 [95% CI] | difference vs head at 0.9 [95% CI] | per donor | at 0.8 | verdict |
|---|---|---|---|---|---|
| head | 0.8319 [0.7462, 0.9304] | | | | |
| R1 | 0.7911 | −0.0409 [−0.0554, −0.0121] | −0.0503, −0.0208 | −0.0356 [−0.0497, −0.0009] | loss |
| R2 | 0.6394 | −0.1925 [−0.2899, −0.0836] | −0.089, −0.2789 | −0.1978 [−0.2815, −0.0729] | loss |

*Secondary:*
- The head keeps 0.1169 of the measured NK–T gap on flagged pairs, against 0.3985 on unflagged pairs. The two intervals overlap, and this contrast has no registered test of its own.
- R1's flagged-pair recovery term is 0.1465 [0.0493, 0.1858], positive in both donors (0.1166 and 0.1056). Its scale-free contrast Dlog_R1 is 0.4384 [−0.3074, 1.3336], which does not exclude 0. So R1 keeps more of the gap on flagged pairs, but not specifically more than a rescaling would.
- The point estimates of D_R1 and D_R2 exceed 0.05. The intervals include 0, and donor 19593 is negative.
- Without gdT CD158b+ cells (163 pairs), D_R2 is a loss: −0.3165 [−0.5276, −0.1315].
- The head's flagged-pair compression is not robust:
  - By donor (values from the JSON, which stores 4 decimals), the head's gap ratio is 0.0694 flagged vs 0.5345 unflagged in 13272, but 0.3978 vs 0.4214 in 19593.
  - Without gdT CD158b+ cells (163 pairs, 86 flagged) it is 0.4729 [0.3246, 0.7898] flagged vs 0.5137 [0.3258, 0.7393] unflagged.
  - So the pooled 0.1169 vs 0.3985 comes mostly from gdT CD158b+ pairs in donor 13272, the donor with the weak key.
- The seed-1 null head separates NK from T better than the registered head: accuracy 0.8866 at 0.9; null − head +0.0547 [0.0185, 0.0828], per donor 0.0212 and 0.0762. This is descriptive, with no registered rule. Retraining the head with another seed raised NK-vs-T accuracy by more than R1 lowered it (+0.0547 vs −0.0409), so head-training variance is not small next to these comparisons.
- Ranking the head's calls by kNN label disagreement among training cells gives 0.9134 at 0.9, against 0.8501 for the head's margin.
- On test_secondary, every H3a contrast is inconclusive, and R1 vs head is equivalent at 0.9.

**Limitations.**
- 253 pairs from 2 donors give wide intervals.
- R2 fits worse than R1 even on val (val Pearson CD56 0.5514 vs 0.7091).
- The pairs are defined on the primary key, which is weak in donor 13272.

> **中文：** E3："损失位于读出/池化并可在冻结 TEDDY 上修复"被拒绝；R1、R2 的 NK-T 准确率均输给原头部（注册判定）。头部在近邻配对上压缩差异更甚仅为描述性结果，区间重叠，且主要来自供体 13272 的 gdT CD158b+ 配对。

### E4: two-channel fusion (C6; Experiment 2 redesigned)

- **Question:** given two noisy evidence channels, does ANM's field fuse them better than declared rules or trained fusers? The channels:
  - channel 1: TEDDY + head, predicting the panel from RNA;
  - channel 2: a ridge regressor from the 122 non-panel measured proteins, through A1.2's panel-free renormalisation.
- **Design:**
  - E4 key: Q1 scoring on the *measured* panel. A cell is OUT when its top score is below τ_K = 0.55, chosen on val. Neither channel reads the measured panel.
  - Noise levels: L0 none; L1 RNA thinned to 0.2; L2 half of channel-2 inputs set to 0; L3 both.
  - Trust: per channel and level on val (the mean of max(0, Pearson) with the measured panel). ANM's source scale is trust ÷ the larger trust.
- **Methods:**
  - F0: channel 1 alone. F0b: channel 2 alone.
  - F1: simple average (declared rule). F2: trust-weighted average (declared rule).
  - F3: stacker, a logistic regression on the 8 class scores trained on val per noise level.
  - F4: learned fusion, a logistic regression on the 24 evidence values trained on training cells at L0.
  - **F5: TEDDY + ANM fusion**, with per-channel source scale and contradiction events (each channel's top class contradicts the other classes at 0.5 × its top score), all at t = 0. F5's scores equal the closed form R(n_k)·[3 s1 S1_k + 3 s2 S2_k − 0.5 × each disagreeing channel's top score × its scale], where R(n_k) is the field's response per unit source for a class whose star has n_k event sites (its 6 support events plus the contradiction events attached to it), s1 and s2 are the source scales and S1_k, S2_k are the channels' class scores (code: `v3_e4.anm_fusion_closed`).
  - Controls: C1 no propagation; C2 channel 2 at t = 2; C3 no contradiction events (= F2 up to a constant).
- **Comparator** (A1.9, chosen on val): F3. Val endpoints: F3 0.907165, F2 0.877892, F5 0.877021.
- **Endpoint:** selective accuracy at coverage 0.8, averaged over L0–L3; margin 0.01.
- **Registered verdict:** "primary F5 vs F3 at 0.8: **loss (field hurts)**, point −0.0815 [−0.1380, −0.0535]; field effect F5 vs F2: **inconclusive**, point +0.0034 [−0.0007, +0.0103]."
  - "Field hurts" is the registered name for F5 losing to the best non-ANM fusion. The field's own effect, F5 vs the trust-weighted average, is inconclusive.
- **Key validity:** E4 key kappa vs annotation on test_primary 0.7668 (< 0.85), so the v3 primary key is reported with equal prominence. It gives the same verdicts for the two registered comparisons: F5 − F3 −0.1125 [−0.1610, −0.0687] loss; F5 − F2 inconclusive. E4's REPORT.md notes that this key is partly circular for channel 2.

| method (E4 key, test_primary) | mean selective accuracy at 0.8 [95% CI] |
|---|---|
| F0 channel 1 alone | 0.8361 [0.7645, 0.8753] |
| F0b channel 2 alone | 0.7353 [0.7235, 0.7565] |
| F1 simple average | 0.8360 [0.7808, 0.8656] |
| F2 trust-weighted average | 0.8504 [0.7887, 0.8838] |
| **F3 stacker (val-trained)** | **0.9353 [0.9304, 0.9412]** |
| F4 learned fusion | 0.8693 [0.8601, 0.8852] |
| **F5 TEDDY + ANM fusion** | **0.8538 [0.7974, 0.8841]** |
| C2 (channel 2 at t = 2) | 0.8541 [0.7977, 0.8847] |

| comparison at 0.8 | difference [95% CI] | per donor (13272, 19593) | decision |
|---|---|---|---|
| F5 − F3 (primary) | −0.0815 [−0.1380, −0.0535] | −0.0553, −0.1340 | loss |
| F5 − F2 (field effect) | +0.0034 [−0.0007, +0.0103] | +0.0000, +0.0091 | inconclusive |
| F5 − F0 | +0.0177 [+0.0074, +0.0359] | +0.0095, +0.0329 | win |
| F2 − F0 | +0.0143 [+0.0072, +0.0269] | +0.0095, +0.0238 | win |
| F5 − F4 | −0.0155 [−0.0818, +0.0187] | +0.0165, −0.0775 | inconclusive |
| C2 − F5 (time order) | +0.0003 [−0.0003, +0.0009] | +0.0006, +0.0004 | equivalent |

*Secondary:*
- F5 − F2 at coverage 0.9: +0.0030 [+0.0017, +0.0067] → equivalent.
- On the E4 key, fusing the two channels by the declared trust-weighted average beats channel 1 alone (F2 − F0 win). The result does not carry over: on the v3 primary key it is inconclusive (+0.0073 [−0.0036, +0.0183], per donor −0.0009 and +0.0152), and on test_secondary it is equivalent (−0.0030 [−0.0071, +0.0013]).
- ANM's fusion adds nothing measurable over that rule (F5 − F2 inconclusive at 0.8, equivalent at 0.9).
- Declaring channel 2 at a later step (t = 2) changes nothing (equivalent).
- C3 selects exactly the same cells as F2 at every level and coverage (0 differing cells). The ANM engine equals the closed form on every cell (bridge failure False).
- On test_secondary, F5 − F3 is a loss and F5 − F2 is equivalent.

**Limitations.**
- F3 is trained on the validation donor's labels, so it uses information the declared rules do not.
- The E4 key is protein-based and only moderately consistent with the annotation.
- With every event at t = 0, ANM's retention contributes only a constant gain per star size.

> **中文：** E4：ANM 融合输给在验证集上训练的堆叠器（注册标签"field hurts"）；相对于信任加权平均（即去掉矛盾事件的 ANM），差异不确定；时间顺序不改变结果。

### E5: looking inside TEDDY's layers: layer dynamics (C7)

- **Question:** push one gene's token embedding at its input position, and follow the response through each of TEDDY's 12 layers by Jacobian-vector products. In which layer does the NK–T component of the response shrink most, and does that match the layer where a linear NK–T protein probe's gap ratio drops most?
- **Design:**
  - Cells: the unique test_primary cells of E3's within-donor pairs, plus matched random NK/T cells (100 per donor); 478 cells, 7,862 cases.
  - Endpoint: l* = the layer with the most negative median log-gain log|r_l| − log|r_{l−1}|, where r_l is the response's projection on the training NK-minus-T direction.
  - Registered gate: central finite differences at ε = 1e-3 must agree with the JVP (relative error ≤ 0.05 and cosine ≥ 0.99) in ≥ 95% of cases. Otherwise E5 stops with "JVP not validated".
- **Registered verdict:** **"JVP not validated (registered check failed): E5 stops; registered reading: a linear-response description is not valid at this scale. Note: the float64 version of the same check passes (S5), so the failure is float32 rounding; see the addendum's infeasibilities."**
- **Key numbers:**
  - Float32 check at ε = 1e-3: pass share 0.854 of 2,671 cases, below the required 0.95. At ε = 1e-2: 1.000 of 571 cases.
  - Float64 check (S5): 1.000 of 79 cases at both ε.
  - The registered endpoint was still computed and is reported because it is registered; it does not change the verdict. l* = 11 (median log-gain −0.706; l* is layer 11 in 0.9775 of the bootstrap draws). The probe comparator layer is 9, so the two do not agree within the registered 1 layer. Specificity vs random directions D = 0.115 [−0.007, 0.174] (3 decimals from the JSON), per donor 0.117 and 0.087 → inconclusive.
- *Secondary (addendum; never changes the verdict):*
  - S1: for a push along the training NK − T input direction, the largest single-layer drop of its NK–T component is at layer 2. Median step change −0.217 [−0.229, −0.197] (3 decimals from the JSON), per donor −0.224 and −0.208; holds. By layer 5 the median component is 0.285 of the input's (0.503 at layer 12). The probe's argmin layer is 9.
  - S2: look-alike cells do not damp the push more than random NK/T cells.
  - S3: the head does not pass less of the push than the layer-12 state.
  - S4: one-sided finite differences are linear up to ε = 0.02.
- **Limitations.**
  - The addendum lists deviations from the registered path. Apple-GPU (MPS) float32 replaced the CPU, and explicit attention arithmetic replaced the fused kernels, which do not support forward-mode derivatives. No amendment adopting the float64 check was committed before site4, so the registered verdict stands.
  - Secondary findings are descriptive.

> **中文：** E5：注册的有限差分检验在 float32、ε = 1e-3 下通过率 0.854 < 0.95（float64 通过），E5 按注册停止；次要分析不改变结论。

### E5-M: is TEDDY's gene-mean state sufficient for an independent observer? (pending)

**Status: pending.** No E5-M outcome is reported in this version. The registered design follows (addendum `E5M.json`, sha256 `65ab9efe…`, committed 07:39:16 UTC and pushed 07:57:00 UTC, before the site4 run started at 10:43:09 UTC).

- **Question:** at first order, is TEDDY's gene-mean layer-12 state (what our head reads) sufficient for an independent observer of the layer-12 token states, on the NK–T look-alike cells? (The addendum describes these as the cells "where the head loses the NK–T difference". That is its motivation, not a v3 result: E3's flagged-vs-unflagged contrast is descriptive only.)
- **ANM framing (from the addendum):** the gene-mean state is a *declared candidate* retained state, tested against *declared* observers. The estimator H = Λ·pinv(Γ) is ANM's own clamp estimator, computed from clamps of the candidate state alone and never fitted on the tested responses. The readout-visible quotient is reported only as a diagnostic and is never used as the state.
- **Observers:**
  - O_tok: E3's frozen R2.
  - O_rand: a random token observer with R2's form, the negative control.
  - O_head: the head on the gene-mean state, a positive control whose residual is 0 by construction.
  - Outputs: CD56, CD94, CD335 and CD3.
- **Statistic:** per-cell residual r_O = ‖K_O − H·K_z‖ / ‖K_O‖, the share of the observer's response that flows through what gene-mean pooling discards.
- **Cells:** 60 pair cells and 40 random cells per primary donor, from E5's set.
- **Hypothesis H5M:** on NK–T pair cells the residual of O_tok exceeds that on matched random cells. Δ_tok ≥ 0.05 under the common win rule (lower bound > 0, ≥ 0.025 in each donor), and Δ_tok also exceeds the random-observer null (DD_rand lower bound > 0).
- **Falsified** if Δ_tok is equivalent (interval inside ±0.05) or a loss.
- **Verdict labels:** PASS_SUFFICIENT, REJECT_SUFFICIENCY, H5M_FALSIFIED, INCONCLUSIVE, or NOT_VALIDATED (if the finite-difference gate or the positive control fails).

> **中文：** E5-M（待定，本版未报告结果）：检验 TEDDY 的基因均值状态对一个独立的 token 观察者是否充分；使用 ANM 自身的钳制估计量。

### E6: external confirmation on Hao et al. 2021 PBMC CITE-seq (pending)

**Status: pending.** No external outcome is reported in this version (`outputs/v3/E6` did not exist at the time of writing). The registered design follows (addendum `E6.json` v2, sha256 `b82687a7…`, pushed 10:43:44 UTC before any external run).

- **Data:** GSE164378, 3′ data; 40,000 cells, 5,000 per donor from 8 donors, chosen by donor id only.
- **Frozen:** TEDDY, the head, the evidence normaliser, panels, bars, classifier and gate logic.
- **Gate:** thresholds re-estimated within each external donor, label-free. CD33 and CD94 are absent from the external ADT panel, and their gate conditions are removed. CD5 and CD94 are absent panel proteins, but their evidence is the head's prediction from RNA, so it exists; a reduced-panel reading is a sensitivity.
- **Endpoints:** E1.1a, E1.4a, E1.C2 and Q1 selective accuracy at matched coverage, per donor and pooled.
- **Replication rule:** a v3 conclusion replicates if its sign holds in at least 5 of 8 donors and the pooled two-stage interval excludes 0 on the same side. By the addendum's rule, a v3 point of exactly 0 gives "no v3 direction"; E1's C2 point at c* is 0.000.
- **Key not validated** if the kappa between the annotation and the per-donor gate is below 0.70. Results are then reported on the annotation-only key.
- **Contamination label:** "RNA possibly seen by TEDDY in pretraining". The annotation map (31 labels) was rederived independently with 31 of 31 agreeing.
- **Disclosed deviation:** the external map's first version was revised after the prep step had read external values, so the registered order "committed before any external value is read" was not met in time.

> **中文：** E6（待定，尚无外部结果）：在外部 PBMC 数据上冻结所有规则复现 E1 主要终点；注意该数据的 RNA 可能已被 TEDDY 预训练见过。

---

## 4. Claims ledger

### 4.1 What we can claim

Status: **registered** means a registered verdict or check decides the claim. **Secondary** means a reported comparison that never decides a verdict. **Descriptive** means a measurement with no decision rule.

| # | claim | status | evidence row |
|---|---|---|---|
| C1 | In this setup ANM's engine computes exactly the declared rule. This is a check, not a gain. | registered (E1.1a pass) | E1.1a: 0 call mismatches for Q1, Q2, Q3 (vs the rule) and the closure readout (vs its re-coded form), test_primary and test_secondary. E2: 0 mismatches over all perturbation cases. E4: engine = closed form, bridge failure False. |
| C2 | ANM's leave-one-out attribution equals the exact closed form. | registered (E1.3a discrepancy False) | E1.3a: top-set and flip agreement 1.0000 on 10,715 Q1 calls. |
| C3 | Pooled over both primary donors, one panel protein decides more calls than the 1/3 null in each class. This is a property of the evidence, and any rule reproduces it. For B it does not hold by donor. | secondary, descriptive | E1.3b: CD11c 0.832 of myeloid, CD2 0.618 of T, CD94 0.572 of NK, CD20 0.484 of B calls; null 1/3, p = 0.0010. In donor 19593, CD20 0.339 (p = 0.6494) and CD22 0.491 of B calls. |
| C4 | A logistic classifier on the same 12 evidence values orders calls better than the rule's (= ANM's) top score. | secondary (E1.4) | AURC 0.942 vs 0.843; difference −0.099 [−0.163, −0.034], loss for the top score. |
| C5 | A classifier trained only on the old question's labels answers the new CD14-anchored question better than the zero-label rule. | registered (E1.1c loss; E1.5 n*) | E1.1c: −0.078 [−0.134, −0.028], loss for the rule; annotation-only −0.054, loss. E1.5: n* = 0. |
| C6 | Measured operating point of TEDDY + head + Q1 rule on unseen-donor site4 cells: coverage 0.949, selective accuracy 0.697 (primary key) / 0.730 (annotation-only), OUT decline 0.261. | descriptive | E1 deployed operating points. |
| C7 | TEDDY + head + rule loses decision accuracy as RNA is thinned: 0.022 with 80% of UMIs kept, 0.162 with 20% kept (0.140 with 5% kept). | descriptive (E2's comparator curve) | E2 curve table. |
| C8 | The pipeline's response to these perturbations is not linear in their size. | descriptive (C5, a registered measurement with no decision rule) | E2 C5: at most 0.098 of cells in the linear regime at any stage. |
| C9 | Two trained readouts of frozen TEDDY (a layer-mean MLP and attention pooling) separate NK from T worse than the registered head. | registered (E3.H3b losses) | E3.H3b at 0.9: R1 −0.0409, R2 −0.1925; also losses at 0.8. Descriptive context: a seed-1 retrain of the head scored +0.0547 above it. The head's gap ratio of 0.1169 on flagged pairs vs 0.3985 on unflagged pairs is descriptive only: the intervals overlap, and without gdT CD158b+ pairs the ratios are 0.4729 vs 0.5137. |
| C10 | A val-trained stacker fuses the two channels best, and ANM's fusion loses to it. On the E4 key, the declared trust-weighted average also beats TEDDY + head alone. | registered (F5 − F3 loss); secondary (F2 − F0) | E4: F3 0.9353 vs F5 0.8538. F2 − F0 is +0.0143 [+0.0072, +0.0269], a win on the E4 key, but inconclusive on the v3 primary key (+0.0073) and equivalent on test_secondary (−0.0030). |
| C11 | The registered E5 derivative check failed in float32 at ε = 1e-3. The registered verdict attributes this to float32 rounding, because the float64 version passes. | registered (E5 gate) | E5: 0.854 < 0.95; S5 1.000. |

### 4.2 What we cannot claim (and which test rejected it)

| # | claim we cannot make | rejected by |
|---|---|---|
| N1 | "ANM adds decision value as a decision layer on TEDDY (Mode B)." | E1 falsification: neither E1.4 nor C2 is a win; every ANM arm equals a re-coded rule. |
| N2 | "ANM's confidence is a better trust gate than TEDDY's margin." | E1.4a: −0.002 [−0.024, 0.020], inconclusive; the gate claim stays withdrawn. The 0.90-coverage win is secondary and belongs to the rule's top score. |
| N3 | "ANM's nested (closure) readout declines out-of-scope cells better." | C2 at c*: 0.000 [−0.006, 0.003], equivalent; in-scope guard failed (−0.056). |
| N4 | "Editing the question with zero labels beats training" / "ANM saves labels." | E1.1c/E1.1d: loss (−0.078); E1.5: n* = 0. Zero labels holds for any written rule. |
| N5 | "Attribution flags wrong calls better than TEDDY's margin." | E1.3c: 0.004 [−0.008, 0.022], inconclusive (Q2 equivalent). |
| N6 | "The representation / head / decision decomposition tells more than the accuracy curve." | E2: adds nothing beyond the curve (falsified); criterion (ii) is a loss (−0.092). |
| N7 | "The NK–T loss is in the readout or pooling and is repairable on frozen TEDDY." | E3 falsification (A2.1): rejected; H3a inconclusive for R1 and R2; H3b losses. |
| N8 | "The NK–T loss is in mean pooling." | E3.H3a pooling (A3.1): inconclusive, −0.0375 [−0.4408, 0.22]. |
| N9 | "ANM's field fuses evidence channels better." | E4: F5 − F3 −0.0815 loss ("field hurts"); F5 − F2 inconclusive (+0.0034); C2 (t = 2) equivalent. |
| N10 | "A specific TEDDY layer localises the NK–T loss." | E5: JVP not validated, so E5 stops; the computed l* = 11 does not agree with the probe layer 9 either. |
| N11 | Any population-level statement. | Two primary donors; the intervals describe them. |
| N12 | Any v1/v2 headline (for example the first run's ANM accuracy edge, or "0 vs 50–1,000 labels") as a v3 result. | Superseded by E1 under the v3 key; see section 1. |

### 4.3 What was not tested

| # | not tested | why / status |
|---|---|---|
| T1 | ANM's field dynamics: retention over real processing steps, multi-step propagation, histories, feedback. | Every event enters at t = 0, so the field acts as a constant gain per class (G(3)·3·S). E4's t = 2 control is a declared order on one channel, not a history. |
| T2 | Whether TEDDY's gene-mean state is sufficient for an independent token observer on look-alikes. | E5-M: pending; not reported in this version. |
| T3 | External replication on another dataset. | E6: pending; RNA possibly seen by TEDDY in pretraining. |
| T4 | Pretraining contamination of TEDDY with GSE194122 (this dataset). | Not verified (registration section 13). |
| T5 | Effect of the head's global ADT-median size factor (computed over all cells, site4 included) on its weights. | Not testable without retraining. |
| T6 | Retraining TEDDY or the head; other heads, panels or questions beyond the registered ones. | Out of scope; TEDDY and the head are frozen. |
| T7 | Comparison with TypeSafe AI's Jev. | Never called; the classifier is a stand-in. |
| T8 | Validated layer-wise linear response in TEDDY. | E5 stopped at the registered gate; the float64 path passed but was not registered before site4. |

> **中文：** 账本：能说的是 ANM 精确复现声明规则、分类器优于规则、两个替代读出的 NK-T 准确率不如原头部等（每条标明注册判定、次要或描述性）；不能说 ANM 提升决策、信任门控、嵌套读出、节省标签、融合或层定位；未检验的是 ANM 的场动力学、E5-M、E6 及污染问题。

---

## 5. What ANM brings to TEDDY, honestly

**The paper's framing.** In ANM v2, the retained state is a *declared candidate* tested against a declared readout: "a candidate retained coordinate … Its choice does not establish sufficiency" (ANM v2 paper, `sections/protocol.tex` lines 7–11). The readout-visible quotient is a diagnostic in the SI, not "the state ANM finds". The AI workflows in the paper are open-loop with memory: retention, multi-step propagation and request–return histories, with no feedback ("Sources are thus prescribed and propagation is fixed and open-loop: no feedback mechanism acts inside the graph", `sections/ai_evidence.tex` lines 8–10). ANM's role there is a **response-theory verifier**, not a better decision-maker.

**What that means in teddy_mm v3.**

1. **A check that the declared readouts are what the engine computes.** Each readout was stated as a declared object and checked against ANM's engine cell by cell:
   - E1.1a: 0 mismatches.
   - E1.3a: 1.0000 agreement after A1.3. Reading the engine's leave-one-out, which deletes an event's site, showed that the registered closed form re-coded zeroing an event, not deleting it.
   - E4: engine = closed form.
   These checks passed. They are bookkeeping: they show ANM computes the declared rules, not anything new about TEDDY. The paper's response-theory verification corresponds here to E2, E5 and E5-M (point 3).
2. **No better decisions, by construction here.** All events enter at t = 0, so ANM's decision-layer calls equal a re-coded declared rule. Where its readout forms (top action score, closure readout, contradiction fusion) were tested, they did not beat TEDDY's margin, the mean rule, a trained classifier or a val-trained stacker (E1, E4). ANM's fusion did beat channel 1 alone on the E4 key (F5 − F0 +0.0177, win), but so does the declared trust-weighted average (F2 − F0 +0.0143, win), and F5 − F2 is inconclusive. Calling this "TEDDY + ANM decides better" would be wrong. Where ANM's field does not run, the method is a declared rule.
3. **The response-theory tools have not yet produced a validated finding about TEDDY.**
   - E2's decomposition added nothing beyond the accuracy curve.
   - E5's registered derivative gate failed in float32.
   - E5-M is the one test that uses ANM's own sufficiency machinery on TEDDY's retained state: TEDDY's gene-mean as a declared candidate, ANM's clamp estimator, independent observers. It is pending and not reported in this version. Until it is, there is no ANM-based statement about whether TEDDY's embedding is sufficient for the NK–T distinction.
4. **What TEDDY users can take away now** (not ANM-specific):
   - The pipeline's measured operating points under a five-class key with an explicit validity check. That check flagged donor 13272.
   - A trained classifier on the same evidence beats the declared rule even for a changed question.
   - Swapping in a layer-mean MLP or an attention-pooling readout of frozen TEDDY did not improve NK-vs-T accuracy over the registered head (E3.H3b losses). Whether the head compresses the NK–T gap more on look-alike pairs is not established. The pooled point estimates (0.1169 flagged vs 0.3985 unflagged) have overlapping intervals and come mostly from gdT CD158b+ pairs in donor 13272. Without those pairs the ratios are 0.4729 vs 0.5137.
   - The answer key itself is fragile for one test donor.

> **中文：** 论文中 ANM 的角色是"响应理论验证器"。在 teddy_mm 中，引擎与声明读出逐细胞一致（0 不一致、1.0000 一致），这只是核对，不是关于 TEDDY 的新发现；所有事件在 t = 0 时 ANM 不带来更好的决策；E2、E5 未得出经验证的结论，关于 TEDDY 嵌入是否充分的 ANM 检验（E5-M）待定、本版未报告。

---

## Numbers

`docs/reports/v3_numbers.json` is a flat map `{key: {value, source_file, source_key}}` of every number used above (907 entries).

- For a `.json` source, `source_key` is an RFC 6901 JSON Pointer, and the value was read from the file by a script, not typed.
- For `README.md` and `REPORT.md` it is a line number.
- For git it is the commit and its committer date.
- For GitHub it is the push event on `exp/registered-v3`.
- For `progress.log` it is a line number with its local time (EDT = UTC − 4).
- Every timestamp value is in UTC.

`outputs/…` paths are under `/Users/tianchichen/Documents/GitHub/teddy_mm` and are not committed. `registration/…` and `README.md` are in this repository.

> **中文：** 本文用到的每个数字都以 `{value, source_file, source_key}` 形式收录于 `v3_numbers.json`，供之后更新页面使用。
