# teddy_mm v3: registered results and claims ledger

**Status (2026-10-01).** All eight registered experiments are final: E1, E2, E3, E4, E5, E5-M, E6 and E7. The overall conclusion is in [section 6](#6-overall-conclusion).

**Sources.** Results: `outputs/v3/{E1,E2,E3,E4,E5,E5M,E6,E7}/REPORT.md` and `*_results.json` (E6 also `key_validity.json`; E7 also its development selection record `registration/addenda/E7_selection.json`). Registration: `registration/REGISTRATION_v3.md` (`registration_v3.json`, sha256 `e4c8a33e…`), amendments `AMENDMENT_A1.md`, `AMENDMENT_A2.md`, `AMENDMENT_A3.md`, addenda `registration/addenda/*.json` (E7 also `registration/E7_REGISTRATION.md`). Earlier (v1/v2) results: `README.md` and `docs/reports/*_OFFICIAL.md`. Every number used here is in [`v3_numbers.json`](v3_numbers.json), each with its source file and key (section [Numbers](#numbers)).

**Precision.** Numbers are printed at the precision of the experiment's own `REPORT.md`: 4 decimals for key validity, E3 and E4; 3 decimals for E1, E2, E5, E5-M and E7. E6 prints kappa and pooled differences with 4 decimals, and AURC, accuracy, coverage and per-donor differences with 3. Where a number appears only in a results JSON, I round it to the same precision and say so. Full-precision values are in `v3_numbers.json`.

> **中文摘要：** 本文汇总 v3 全部八项预注册实验（E1–E5、E5-M、E6、E7）的最终结果，列出能说、不能说和未检验的结论，并给出可引用的总体结论（第 6 节）；所有数字见 `v3_numbers.json`。

---

## Plain names used here

| name used here | what it is |
|---|---|
| **TEDDY embedding** (the gene-mean state, written z_12 in E5-M) | The official TEDDY-G embedding of a cell's RNA: the mean over gene tokens of TEDDY's last-layer (layer 12) token states, 512 numbers, L2-normalised. Our head reads it. TEDDY is never retrained. |
| **head** | The phase-1 head trained in this repo: an MLP with a negative-binomial decoder that predicts 134 surface proteins from the TEDDY embedding (`outputs/cite_phase1_official/best.pt`), frozen. |
| **evidence** | For each panel protein: the head's prediction divided by the 95th percentile of the head's own predictions over the training cells, clipped to [0, 1]. No measured protein enters it. |
| **declared rule** (TEDDY + fixed rule) | For each class, the mean evidence over its panel; call the highest class if its score reaches the question's bar, otherwise "no call". |
| **TEDDY + ANM** | The same evidence fed to ANM's engine (`finite_graph_scalar`, ANM's default field: 4 steps, retention 0.82, diffusion 0.16), one event per panel protein, **all events at t = 0**. Used only where ANM's field actually runs. |
| **decision layer (Mode B)** | ANM used on top of TEDDY's outputs, as here in E1, E2, E4 and E6. **Looking inside TEDDY (Mode A)** means studying TEDDY's layers (E5, E5-M). |
| **observer** (E5-M) | A declared readout of TEDDY's layer-12 token states whose response is tested. The trained token readout (code O_tok) is E3's attention-pooling readout R2. The random token readout (O_rand) has R2's form with random weights, a negative control. The head readout (O_head) is our head, which reads only the gene-mean, a positive control. |
| **clamp residual** (E5-M) | The share of an observer's first-order response to a gene push that flows through what gene-mean pooling discards. 0 means the gene-mean carries everything this observer uses at first order. |
| **replication rule** (E6) | A v3 conclusion replicates if its sign holds in at least 5 of the 8 external donors and the pooled 95% interval excludes 0 on the same side. The sign comes from E1's site4 point; a point of exactly 0 gives "no v3 direction". |
| **OUT** | Out of scope: cells that are none of B, T, NK or myeloid (erythroid cells, progenitors, plasma cells, pDC, ILC). The correct action is no call. |
| **primary key** | A cell's class when the annotated cell type and the protein gate (measured CD3, CD19, CD56, …) agree; otherwise the cell is **unscored**. **Annotation-only key**: the annotated class for every cell. |
| **coverage** / **selective accuracy** / **decision accuracy** | Share of all cells that get a call / share of called, scored cells whose call is right (a call on a key-OUT cell is wrong) / share of scored cells with a right call or a no call on a key-OUT cell. |
| **OUT decline rate** | Share of key-OUT cells that get no call. |
| **matched coverage**, **AURC** | Every method calls its own most confident ⌈c·N⌉ cells, so accuracy is compared at equal coverage; AURC is the area under selective accuracy over coverage 0.95 … 0.50, divided by 0.45. |
| **workability (P_f)** / **exactness (Q_f)** | P_f: share of cells that get a call. Q_f: share of scored cells whose lineage call equals the key (as the E1 code computes it, a key-OUT cell never counts as correct here). |
| **verdict words** | *win*: the difference reaches the margin, the 95% interval's lower bound is above 0, and in each primary donor the difference has the same sign and is at least half the margin. *loss*: the mirror. *equivalent*: the interval lies inside ±margin. Otherwise *inconclusive*. Intervals: two-stage bootstrap (donors, then cells within donors), B = 2000, seed 1. |

> **中文：** 术语表：TEDDY 嵌入、头部、证据、声明规则、TEDDY + ANM（所有事件在 t = 0）、OUT（范围外）、主答案键与仅注释键、覆盖率与选择准确率、观察者与钳制残差（E5-M）、复现规则（E6）、判定词（win/loss/equivalent/inconclusive）。

---

## 0. Bottom line

- **ANM equals a declared rule in this setup, and the registered tests found no decision value from ANM's readouts.** ANM's engine made the same call as the declared rule on every test cell for all three questions, and its closure readout matched its re-coded closed form on every cell (E1.1a: 0 mismatches). E1's falsification clause holds: neither the trust-gate test (E1.4) nor the nested-readout test (C2) is a win.
- **Writing a new question with zero labels did not beat training.** On the changed question Q3, a classifier trained only on the old question's labels was more accurate than the declared rule (rule minus classifier −0.078 [−0.134, −0.028], loss), and the label-cost curve gives n* = 0.
- **E2:** splitting perturbation losses into representation, head and decision shares added nothing beyond the plain accuracy curve (falsified). **E3:** the claim that the NK–T loss sits in the readout or pooling and can be repaired on frozen TEDDY is rejected; both trained readouts lost NK-vs-T accuracy to our head. **E4:** the ANM fusion arm lost to a stacker trained on the validation donor (registered outcome label "field hurts"); against the trust-weighted average, which is ANM without its contradiction events, it is inconclusive. **E5:** the registered derivative check failed in float32 (it passes in float64), so E5 stops with "JVP not validated".
- **E5-M (does TEDDY's gene-mean state lose what a token readout uses?):** registered verdict **inconclusive**. For E3's attention-pooling readout R2, at first order, the median share of its response that flows through what gene-mean pooling discards is 0.085 on NK–T look-alike cells and 0.127 on matched random cells: smaller on look-alikes, the opposite of the hypothesis. The difference is −0.043 [−0.062, −0.011]. It is not a registered loss (the point does not reach −0.05) and not equivalence (the interval extends past −0.05). So neither "pooling loses the look-alike difference" nor "the gene-mean state is sufficient" is established, and the measurement holds for this one observer only.
- **E6 (external PBMC data; labelled "RNA possibly seen by TEDDY in pretraining" and "key not validated"):** ANM equals the declared rule in every external donor (E1.1a replicates). The primary E1.4a direction replicates: TEDDY's margin ranks the rule's calls better than the rule's (= ANM's) top score, −0.0087 [−0.0182, −0.0017], in 8 of 8 donors. Under E1's own decision rule this is still inconclusive, and on site4 the same comparison was inconclusive (−0.002 [−0.024, 0.020]). E1's nested-readout difference at c* is exactly 0, so the replication rule gives "no v3 direction"; externally it is 0.0000 as well (equivalent). The site4 finding that a classifier orders calls better does not replicate: every arm's AURC lies between 0.983 and 0.993 on these data.
- **E7 (the inverse loop of the ANM paper's v2 revision, PR #16, on TEDDY at the layer-11 cut, consumer fixed):** registered verdict **REJECTED_NO_REVISION**. For a fixed consumer (TEDDY's layer 12, gene-mean pooling and our frozen head), histories with exactly the same layer-11 gene-mean give readouts that differ beyond the 0.05 tolerance on fresh site4 comparisons: median per-cell difference 0.125 [0.090, 0.161] in units of the measured NK–T gap, 0.094 and 0.153 in the two donors. The declared augmentation (the gene-mean plus the token states of the six marker genes) is rejected as well (0.122 [0.086, 0.153]), and no candidate was selected in development, nor in the declared follow-up. So the first step of the inverse loop (reject with the consumer fixed, on fresh comparisons) is done on TEDDY, the second step returned no revision from the declared library, and **no revised representation is validated**. The comparisons were sensitive (positive control and detectability witness detected); the patches are exact state interventions at the layer-11 cut, not biological variation.
- **The answer key is weak for one of the two test donors** (donor 13272: annotation-vs-gate kappa 0.4177). Annotation-only results are therefore reported with equal prominence. In E1 every annotation-only verdict matches the primary-key verdict. In E4 the alternative key reported with equal prominence is the v3 primary key. It gives the same verdicts for the two registered comparisons, but not for every secondary one (see E4). On the external data the key is not validated either (kappa 0.6562 < 0.70), so the annotation-only key decides E6.

> **中文：** 在本设置中 ANM 与声明规则逐细胞相同（site4 与外部数据皆然），预注册检验未发现 ANM 读出带来决策价值；零标签改问题不胜训练分类器；E2 被证伪，E3 被拒绝，E4 不敌验证集堆叠器，E5 的 JVP 检验在 float32 下未通过；E5-M 判定为不确定（R2 观察者的响应流经池化丢弃部分的份额在相似细胞上反而更小）；E6 在外部数据上复现了"TEDDY 的 margin 对调用排序至少与顶分一样好"（按 E1 自身判定规则仍为不确定），分类器优势未复现。E7（论文逆向闭环在 TEDDY 第 11 层切点上的类比，消费者固定）判定为 REJECTED_NO_REVISION：第 11 层 gene-mean 完全相同的历史，读出差异超过 0.05 容差（site4 中位数 0.125 [0.090, 0.161]），加上六个标记基因 token 状态的候选同样被否决，开发阶段（及事先声明的后续）未选出任何候选；因此完成了逆向闭环的第一步（否决），第二步在声明的候选库中没有得到修订，没有验证任何修订后的表示。

---

## 1. Why v3

The v1/v2 experiments (README.md, `docs/reports/*_OFFICIAL.md`) had problems that made their headline claims untestable or invalid. v3 fixes each one before looking at the test site again.

| problem in v1/v2 | evidence (source) | v3 fix |
|---|---|---|
| **The answer key had only three classes (B, T, myeloid): no NK class and no out-of-scope class.** NK cells, which carry CD16, were scored as "myeloid", and there was no class for which "no call" is the right answer. | 1,342 of 1,690 NK cells were keyed "myeloid" under the soft rule (README.md line 176). Out-of-scope cells still got a call 73–95% of the time (line 145). | Five classes: B, T, NK, myeloid and OUT, where no call is the correct answer for OUT. The primary key keeps only cells where the annotation and a protein gate agree; all others are unscored, and the annotation-only key is always reported too. |
| **The B/T-priority question was synthetic.** Its CD16-anchored "myeloid" class captured many NK cells. | Registration section 7, audit finding KEY-1. The question changed 12.48% of expected answers (README.md line 282). | Q3 is the standard four-population question, anchored on CD19 (B), CD3 (T), CD56 (NK) and CD14 (classical monocyte). Biology fixed it; nothing was selected. |
| **The panel was chosen on the test set.** | The 9 panel proteins were chosen by test-set (site4) Pearson (README.md lines 152, 264). | Panels are selected on the validation donor only, by worst-pair AUROC. The old 9-marker panel survives only as a sensitivity row. |
| **The test set was not donor-independent.** | Donor 15078 is in training at sites 1–3 and is 32.6% of the site4 cells. The "out-of-site check" was the validation donor at a training site (README.md line 152). | Primary evaluation uses only the site4 donors that appear in no other split (13272, 19593). Donor 15078 is a secondary split that never changes a verdict. |
| **Measured protein scales sat in the evidence path.** | The first run multiplied every prediction by a per-cell factor computed from measured protein (the answer key). v2 used one train-median factor, 0.92663, and divided by the training 95th percentile of *measured* protein (README.md line 139; registration section 5). | Evidence = prediction ÷ the training 95th percentile of the head's *own* predictions. A constant size factor cancels, so no measured-protein array enters the evidence. |
| **Event order acted as a hidden weight, and the first run's ANM "edge" came from declining more.** | The first run fed markers to ANM one per step in panel order. On B/T-priority, ANM declined 3,373 cells to the rule's 1,796, scoring 0.9603 vs 0.9428; at ANM's coverage the rule scored 0.9657 (README.md lines 138, 141). | All events enter at t = 0, and every accuracy comparison is at matched coverage. |
| **No registered margins, decision rules or donor-level uncertainty.** Experiments 2 and 6 were withdrawn. "0 labels for a new question" holds for any written rule (README.md lines 146–150). | — | Pre-registered endpoints, margins and win/loss/equivalent rules; a two-stage donor-then-cell bootstrap; failures reported. "Zero labels" is credited to any written rule, never to ANM. |

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
- **Leakage check.** The registration core rebuilt from inputs whose site4 protein, cell types and embedding were replaced by random values is byte-identical to the real core. Poisoning the val rows changes it (`registration/leakage_check_report.json`). A1, A2, A3 and the E1, E2, E3, E4, E5-M and E7 addenda report the same kind of check as passed.

### 2.4 Amendments A1–A3

All times are 2026-10-01 UTC.

| amendment | what it fixed | built | committed (local git) | pushed to GitHub |
|---|---|---|---|---|
| **A1** (sha256 `1f44be66…`) | Ten problems found by an adversarial review before any site4 evaluation. **A1.1:** tied scores broken by a fixed random permutation instead of row order (val: 0.3662 of cells have a Q1 top score of exactly 1), with matched coverage recomputed per replicate. **A1.2:** E4 channel 2 reads a panel-free renormalisation (the stored per-cell centred log-ratio (CLR) transform mixed panel counts into the non-panel inputs). **A1.3:** exact leave-one-out closed form. The registered form re-coded zeroing an event, while ANM deletes the event's site; ρ = G(2)/G(3) = 0.992479, and 6 val calls change flip status. Top markers become sets. **A1.4:** honest naming in E1.4: ANM's soft score is G(3)·3 × the rule's top score, so an E1.4 win belongs to the top-score readout. **A1.5:** the C2 verdict point is c* (the rule's realised coverage), and "in-scope accuracy" is defined. **A1.6:** information parity in E1.1c/E1.5 (the classifier gets a superset of the rule's 4 anchors). **A1.7:** E2 criterion (i) needs ≥ 30 lost cells and an interval excluding 0. **A1.8:** E3 pairs built within each donor. **A1.9:** the E4 comparator is fixed on val. **A1.10:** an addendum cannot replace a registered endpoint. | 04:25:39 | first version `7c49b7a` 04:24:51; final `ca1c8ec` 04:37:04 | 07:57:00 |
| **A2** (sha256 `20b64174…`) | **A2.1:** the E3 statistic D_R is not scale-free. A readout that shrinks every NK–T difference by a common factor can pass it while repairing nothing, so a win now also needs (b) recovery on flagged pairs and (c) a positive scale-free Dlog. **A2.2:** R2's normaliser is estimated on a seeded 10,000-cell training sample (0.33 h of compute instead of 2.20 h; same estimand). **A2.3:** the per-donor clause of the common win rule is binding for E2's criteria. All three make verdicts stricter or keep them unchanged. | 06:09:15 | `62ba605` 06:11:37 | 07:57:00 |
| **A3** (sha256 `4bb2f856…`) | **A3.1:** the E3 pooling contrast D_R2 − D_R1 has the same scale artefact. A pooling win now also needs (b') flagged-pair recovery beyond R1 and (c') a positive scale-free contrast. Training-pair illustration: a pure rescaling reaches the margin at 1 − s = 0.379949. Stricter only. | 07:00:56 | `36d3584` 07:03:09 | 07:57:00 |

### 2.5 Fairness facts and timeline (UTC)

Local git commit times come from the committing machine. The GitHub push events (`repos/danielchen26/teddy_mm/activity`, branch `exp/registered-v3`) are the independent timestamps. All times are 2026-10-01. A site4 evaluation here means the first stage that reads site4 protein or cell types; forward passes of site4 RNA alone are label-free.

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
| E5-M addendum | 07:39:16 | 07:57:00 | E5-M run started 10:43:09 (report written 11:06:15) | yes |
| E6 addendum v1 / v2 | 10:07:39 / 10:43:00 | 10:43:44 | E6 run started 11:54:31; its key-validity stage, the first stage of the run that reads external protein and cell types, ran at 11:54:34 | yes (before any external outcome; see the prep note below) |
| E7 addendum v1 / v2 | 22:01:06 / 22:24:05 | 22:25:53 (both) | E7 development on val 22:26:38; site4 confirmation started 22:40:20 | yes (v2 before any development, follow-up, site4 or external forward) |
| E7 selection records, dev / follow-up | 22:33:13 / 22:40:18 | 22:40:19 (both) | site4 confirmation started 22:40:20 by the local clock: it logged the site4 roster (drawn on the primary key) at 22:40:22, and the first site4 forward followed 22:40:23; external confirmation 22:54:05 | yes, by 1 s on two different clocks (see the E7 note below) |

So:

- **Pushed before evaluation.** The registration, A2, A3 and the E2, E3, E4 and E5-M addenda reached GitHub before the site4 evaluation they govern. The E6 addenda reached GitHub before any external outcome. The addendum records `outputs/v3/E6` absent at 2026-10-01T10:40:07Z, and the E6 run started at 11:54:31. It ran at local commit `cac4f8a` (committed 11:32:45), which differs from the pushed `2a55f81` only in `docs/reports`. Earlier, before either addendum, the data-preparation step had read some external protein summaries while the annotation map was being revised, so E6's registered order for the map was not met in time (disclosed in the addendum; see E6 limitations).
- **E7.** Addendum version 2 (`0e40d3e`, committed 22:24:05) reached GitHub at 22:25:53, before E7's first development forward pass (22:26:38, val donor 18303 only). Version 1 (`6dc5999`, 22:01:06) was first pushed in the same push. The dev selection record was committed at 22:33:13, in the same second the select stage wrote it, and the follow-up development run (val only) started then. It reached GitHub only with the follow-up record (committed 22:40:18) at 22:40:19. The site4 confirmation stage started at 22:40:20 by the local clock, 1 s after GitHub's push time. GitHub's server clock and the local machine's clock are not synchronised to the second, so this 1 s margin alone is not an independent ordering. It is supported by the local commit at 22:40:18, by the confirm stage refusing to run without a committed final record, and by the result file, which records the selection record's sha256 (`5721e548…`), `committed: true` and HEAD `871b709`, the pushed commit.
- **Committed before the evaluation, pushed after it.** For **E1** and **E5**, the governing files (A1, the E1 addendum, the E5 addendum v3) were committed locally before the evaluation, and each result file records their sha256 and `committed: true`. They were pushed to GitHub only afterwards (07:57:00). For these two experiments, ordering rests on local commit times and recorded hashes, not on an independent GitHub timestamp.
- **E1 came before A2 and A3.** E1 finished at 05:44:44, before A2 (built 06:09:15) and A3 (built 07:00:56). Both concern only E2 and E3. A3 states that the E1 site4 results had been seen and that nothing in A3 depends on them. A2 does not say whether they had been seen. The E5-M and E6 addenda record the same disclosure (E1, and for E6 also E3 and E4, had been seen); the E7 addendum records that all published v3 results had been read, including E5-M's site4 per-layer clamp residuals.
- **E3's site4 pairs were built early, by E5.** E5's cell selection at 05:10:36 rebuilt E3's within-donor NK–T pairs on site4 from the primary key, so it read site4 protein and cell types. Its log records the pair count (253) and the cells per donor and class. This happened before the E3 addendum (committed 05:25:05), the E1 addendum, A2 and A3. No gap ratio, D statistic or NK-vs-T accuracy was computed then: the rest of E5's log until its results (10:42:46) holds train/val direction diagnostics and per-cell progress lines. A3 describes E3 evaluate as the only *E3* stage that reads site4 protein or cell types, which is true of E3's own stages, but it does not mention this E5 step.
- **Key validity.** On test_primary the annotation-vs-gate kappa is 0.5781 (agreement 0.6847), below the registered 0.85, so the flag is raised. By donor: 13272 kappa 0.4177 (agreement 0.5617) and 19593 kappa 0.8798. In donor 13272 the gate puts 4,148 of 7,365 cells in OUT, against 1,707 annotated OUT cells. It recovers 0.443 of annotated T cells and 0.395 of annotated B cells (3 decimals from the JSON). A per-batch gate would raise kappa to 0.7019, but it is a sensitivity key only, because it reads site4 protein to set thresholds. 3,561 of 11,294 test_primary cells are unscored on the primary key. **Annotation-only rows are reported beside every primary-key row below.**
- **Other disclosures** (registration section 13):
  - Site4 is not virgin: v1/v2 and an exploratory key prototype read it, and v3's design is informed by what they showed.
  - The head's checkpoint was selected on val.
  - The head's training size factors were divided by one global ADT-median scalar over all cells, site4 included. v3's evidence cancels constant scales; the effect on the head's weights cannot be tested without retraining.
  - Pretraining contamination of TEDDY with this dataset (GSE194122) is not verified.

> **中文：** 预注册、A2、A3 及 E2/E3/E4/E5-M 附录在相应 site4 评估之前已推送到 GitHub，E6 附录在任何外部结果之前已推送（但数据准备阶段在此之前读过部分外部蛋白汇总，见 E6 局限）；E7 附录第 2 版在任何开发运行之前推送，两份开发选择记录在 site4 确认开始前 1 秒推送（GitHub 时钟与本地时钟之差，另有本地提交时间与结果文件中记录的哈希佐证）；但 E1 和 E5 所依赖的 A1、E1 附录、E5 附录只是在评估前本地提交，推送晚于评估。E1 在 A2/A3 之前完成；E5 在 E3 附录、A2、A3 之前已在 site4 上重建了 E3 的配对（仅配对数，未计算 E3 统计量）。供体 13272 的答案键 kappa 仅 0.4177，故仅注释键结果同等呈现。

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

### E5-M: does TEDDY's gene-mean state carry what a token readout uses? (observer sufficiency)

Addendum `E5M.json` (sha256 `65ab9efe…`), committed 07:39:16 UTC and pushed 07:57:00 UTC, before the site4 run started at 10:43:09 UTC. The result file records the addendum hash, `committed: true`, and code hashes that match the addendum.

- **Question:** at first order, is TEDDY's gene-mean layer-12 state (what our head reads) sufficient for an independent observer of the layer-12 token states, on the NK–T look-alike cells? (The addendum describes these as the cells "where the head loses the NK–T difference". That is its motivation, not a v3 result: E3's flagged-vs-unflagged contrast is descriptive only.)
- **ANM framing (from the addendum):** the gene-mean state is a *declared candidate* retained state, tested against *declared* observers. The estimator H = Λ·pinv(Γ) is ANM's own clamp estimator, computed from clamps of the candidate state alone and never fitted on the tested responses. The readout-visible quotient is reported only as a diagnostic and is never used as the state.
- **Design:**
  - Pushes: for each of the 7 declared NK/T genes present in a cell (CD3E, CD3D, NCAM1, KLRD1, NCR1, FCGR3A, KLRF1), three unit pushes of the gene's input token: along its own embedding, along E5's training NK-minus-T direction, and along one fixed random direction. Each push is one column; a cell has a median of 9 columns (range 3–18).
  - Observers, each reading the four outputs CD56, CD94, CD335 and CD3:
    - the trained token readout (code O_tok): E3's frozen attention-pooling readout R2;
    - the random token readout (O_rand): R2's form with a random query and random output weights, matched to R2's attention sharpness on training cells; the negative control;
    - the head readout (O_head): our head on the gene-mean; the positive control, whose residual is 0 by construction.
  - Statistic: the per-cell clamp residual r_O = ‖K_O − H·K_z‖ / ‖K_O‖. K_z is the gene-mean's first-order response to the pushes, K_O the observer's, and H ANM's clamp estimator. r_O is the share of the observer's response that flows through what gene-mean pooling discards.
  - Cells: a GPU-budget subset of E5's cells, fixed in the addendum before site4: per primary donor, the first 60 pair cells and the first 40 random cells of E5's seeded order (200 cells). 5 cells contain none of the 7 genes and so have no push; 195 are analysed: pair cells 60 (13272) and 58 (19593), random cells 39 and 38.
  - Hypothesis H5M: Δ_tok = median r_tok on pair cells − median r_tok on random cells is ≥ 0.05 under the common win rule, and it exceeds the random-observer null (DD_rand = Δ_tok − Δ_rand has a lower bound > 0). H5M is falsified if Δ_tok is equivalent or a loss.
  - Verdict labels: REJECT_SUFFICIENCY (H5M supported), PASS_SUFFICIENT (Δ_tok equivalent *and* the pair-cell level's upper bound < 0.05), H5M_FALSIFIED, INCONCLUSIVE, or NOT_VALIDATED (gate or positive control failed).
- **Gates (both passed):**
  - Finite-difference gate: the Richardson extrapolation of central differences at ε = 0.01 and 0.02 agrees with the JVP in a share of 1.000 of the 357 checks for each target (gene-mean, O_tok, O_head, O_rand). The registered requirement is 0.95, so the derivatives are validated. This is E5-M's own registered gate and does not change E5's verdict.
  - Positive control: the largest head-readout residual is 6.02e-07, below the tolerance 0.001. It holds by construction, so it checks the derivative code, not the pooling decomposition. The decomposition is checked separately: R2's residual response K_O − H·K_z equals the direct JVP of R2 along the pooling-discarded part to within 1.04e-07 (largest absolute difference).
- **Registered verdict:** **"inconclusive (Delta_tok neither win, loss nor equivalent)"**.

| observer | median clamp residual on pair cells | on random cells |
|---|---|---|
| trained token readout (O_tok = R2) | 0.085 | 0.127 |
| random token readout (O_rand) | 0.067 | 0.068 |
| head readout (O_head) | 0.000 | 0.000 (by construction) |

| statistic | point [95% CI] | per donor (13272, 19593) | registered role |
|---|---|---|---|
| Δ_tok (pair − random, R2) | −0.043 [−0.062, −0.011] | −0.048, −0.029 | decides; margin 0.05 → inconclusive |
| DD_rand (Δ_tok − Δ_rand) | −0.041 [−0.061, −0.004] | −0.047, −0.021 | H5M needs a lower bound > 0: not met |
| level: median r_tok on pair cells | 0.085 [0.071, 0.096] | 0.081, 0.087 | PASS_SUFFICIENT needs an upper bound < 0.05: not met |
| median r_tok on random cells | 0.127 [0.098, 0.142] | 0.129, 0.116 | reported |
| Δ_rand (pair − random, random readout) | −0.001 [−0.021, 0.012] | −0.001, −0.008 | reported |
| median r_rand on pair cells | 0.067 [0.056, 0.073] | 0.062, 0.069 | reported |

- **Reading the verdict:**
  - Not a win: H5M predicted a positive Δ_tok, and the measured Δ_tok is negative in both donors.
  - Not a registered loss: a loss needs the point at or below −0.05. The point is −0.043, although both donors pass the half-margin clause.
  - Not equivalent: the interval's lower end, −0.062, lies beyond −0.05.
  - So H5M is neither supported nor falsified by the registered rule. The measured direction is the opposite of H5M: on look-alike pair cells a *smaller* share of R2's response flows through what pooling discards (0.085) than on random cells (0.127).
  - The gene-mean is not shown sufficient either. On pair cells 0.085 [0.071, 0.096] of R2's response flows through the discarded part, above the 0.05 margin. The PASS_SUFFICIENT reading, registered as "z_12 suffices for this observer at first order, so the NK-T loss is in our head, not in pooling", is therefore not reached.
- *Secondary (pre-specified; never change the verdict):*
  - S2, the scale-free version on log r: Δ log r_tok −0.408 [−0.620, −0.115], per donor −0.465 and −0.287, the same direction. For the random readout: −0.021 [−0.298, 0.187].
  - S4, ANM's own `analyze_factorization`, per cell:
    - For R2 at rtol 1e-4 the status is REJECTED in 195 of 195 cells, and the source-kernel test is NOT_INFORMATIVE in 195 of 195.
    - The source-kernel result is what the addendum and the ANM paper predict. With m ≤ 21 push columns against d = 512 dimensions, K_z has full column rank (in a share of 1.000 of cells), its near-null space carries none of K_O (share 0.000 at tolerance 0.05 and at 0.01), and a map fitted on the cell itself reproduces K_O to within 9.9e-15. First-order kernel inclusion is therefore automatic and says nothing; the informative statistic is the clamp residual above.
    - At rtol 0.05, 149 cells are REJECTED and 46 INCONCLUSIVE. The random readout is REJECTED in 195 cells at rtol 1e-4. The head readout is NOT_INFORMATIVE, with its factorization passing, in 195.
    - At rtol 1e-4 any relative residual above 1e-4 (in ANM's operator norm) rejects. So the REJECTED status says the gene-mean is not an *exact* first-order retained state for R2; it does not say how much is lost.
  - S5, the readout-visible quotient (a diagnostic, never the state): R2's response has rank 4 in 128 cells, 3 in 62, 2 in 4 and 1 in 1. Its median participation ratio is 1.221 on pair cells and 1.388 on random cells, so the response is effectively one- or two-dimensional.
  - S7, maps calibrated on training cells: the training-average clamp map gives the same Δ_tok, −0.043 [−0.062, −0.011]; for R2-form observers the clamp map is the same in every cell. A ridge map gives −0.031 [−0.056, −0.004], but the ridge fails the positive control (head-readout residual medians 0.856 on pair cells and 1.355 on random cells), so it is reported only.
  - S1, the gene-mean of earlier layers (16 pair and 13 random cells): R2's clamp residual stays near 1 up to layer 9 (pair-cell median 1.006 at layer 9) and falls to 0.874 at layer 10, 0.529 at layer 11 and 0.107 at layer 12. An earlier layer's gene-mean is far from sufficient for this observer. As the addendum notes, the later layers act on every token state, so not even the head is a function of an earlier layer's gene-mean.
- **Limitations.**
  - One trained observer, and a weak one for this question. R2's attention is nearly uniform: on E3's 278 pair cells its median effective number of tokens is 1,141.957 of a median 1,368 tokens (E3, descriptive). An observer that attends almost uniformly departs little from the gene-mean by construction. R2 also separates NK from T worse than the head (E3.H3b loss). A sharper observer could use more of what pooling discards. **Every E5-M statement holds for this observer only.**
  - First order only (JVPs at each cell's own input), 7 declared genes and 3 push directions. Nothing is said about larger pushes.
  - 195 cells from 2 donors, a GPU-budget subset of E5's 478 cells fixed before site4. The pair cells are defined on the primary key, which is weak in donor 13272.
  - Apple-GPU (MPS) float32 with explicit attention arithmetic, as in E5.

> **中文：** E5-M：注册判定为"不确定"。对 E3 的 R2 观察者（一阶，195 个细胞），其响应流经基因均值池化所丢弃部分的份额在 NK–T 相似配对细胞上为 0.085，在匹配随机细胞上为 0.127，方向与假设相反（Δ_tok −0.043 [−0.062, −0.011]）；未达到 −0.05 的损失边界，也不在等价区间内。配对细胞上的份额 0.085 [0.071, 0.096] 高于 0.05，因此也不能宣称基因均值状态对该观察者充分。ANM 的分解检验在 rtol 1e-4 下拒绝精确分解，源核检验不提供信息（这一点与论文预期一致）。结论只针对这一个注意力近乎均匀的观察者。

### E6: external confirmation on Hao et al. 2021 PBMC CITE-seq

**Every E6 result carries two labels: "RNA possibly seen by TEDDY in pretraining" and "key not validated".**

Addendum `E6.json` version 2 (sha256 `b82687a7…`), pushed 10:43:44 UTC. The run started at 11:54:31 UTC (section 2.5). The result file records the addendum hash, `committed: true`, and committed, clean code.

- **Data:** GSE164378, 3′ data; 40,000 cells, 5,000 per donor from 8 donors, chosen by donor id only.
- **Frozen:** TEDDY, the head, the evidence normaliser, panels, bars, classifier and gate logic.
- **Panels (addendum version 2):** every arm reads the full frozen registered panels, as in E1. So every arm reads the same 12 evidence values, and ANM still equals the mean rule. CD5 and CD94 are not in the external protein panel, but their evidence is the head's prediction from RNA, so it exists. The literal reduced-panel reading (version 1's primary) is a sensitivity that never decides.
- **Gate:** thresholds re-estimated within each external donor, label-free. CD33 and CD94 are not measured externally, so their gate conditions are removed.
- **Endpoints:** E1.1a; E1.4a (primary: top score − margin, AURC); E1.C2 (primary: closure − rule, OUT decline at c*); and Q1 selective accuracy at matched coverage; per donor and pooled. Not run externally: E1.1b–d, E1.3, E1.5, the Q2 and Q3 accuracy endpoints, the MLP and E1's sensitivity panels.
- **Replication rule:** see Plain names. The E1 decision rule (win / loss / equivalent / inconclusive, with each external donor as a "primary donor") is reported beside it and never replaces it.

**Key validity (computed before any method result).**
- The pooled 5-class kappa between the annotation and the per-donor gate is 0.6562, below the registered 0.70, so the label is **"key not validated"** and the annotation-only key decides.
- Per donor (P1–P8): 0.7887, 0.7176, 0.7181, 0.7311, 0.6822, 0.5157, 0.4314, 0.7088.
- The reduced gate misses many B and myeloid cells (recall 0.390 and 0.458 of annotated cells) and puts 10,686 cells in OUT, against 949 annotated OUT cells (OUT precision 0.084).

**Registered results (decision key: annotation-only).**

- **E1.1a:** ANM's engine makes the same call as the mean rule for Q1, Q2 and Q3, and its closure readout equals the re-coded closed form, with 0 mismatches in every donor and 0 rejected events. **Replicates.**

| endpoint | v3 point (E1, site4) | E6 pooled [95% CI] | donors with the v3 sign | replication | E1 decision rule, 8 donors |
|---|---:|---|---|---|---|
| **E1.4a top score (rule = ANM) − margin, AURC (primary)** | −0.0015 | −0.0087 [−0.0182, −0.0017] | 8 of 8 | **replicates** | inconclusive |
| **E1.C2 closure − rule, OUT decline at c* (primary)** | 0.0000 | 0.0000 [0.0000, 0.0000] | — | **no v3 direction** | equivalent (in-scope guard failed, −0.0086) |
| E1.4a top score − classifier, AURC | −0.0990 | 0.0010 [−0.0102, 0.0091] | 2 of 8 | does not replicate | inconclusive |
| E1.4a top score − entropy, AURC | −0.0040 | −0.0095 [−0.0213, −0.0013] | 8 of 8 | replicates | inconclusive |
| Q1 accuracy, top score − margin at coverage 0.90 | 0.0210 | −0.0023 [−0.0046, 0.0008] | 3 of 8 | does not replicate | equivalent |
| Q1 accuracy, top score − margin at 0.70 | −0.0021 | −0.0115 [−0.0234, −0.0025] | 8 of 8 | replicates | inconclusive |
| Q1 accuracy, top score − classifier at 0.90 | −0.0871 | −0.0005 [−0.0104, 0.0070] | 2 of 8 | does not replicate | inconclusive |
| Q1 accuracy, top score − classifier at 0.70 | −0.1149 | 0.0003 [−0.0110, 0.0087] | 2 of 8 | does not replicate | inconclusive |
| E1.C2 closure − rule at 0.80 | 0.0512 | 0.0095 [−0.0314, 0.0407] | 4 of 8 | does not replicate | equivalent |
| E1.C2 closure − rule at 0.70 | 0.0069 | 0.0285 [0.0000, 0.0484] | 6 of 8 | does not replicate | equivalent |

Only the two bold rows are primary; the others are secondary. AURC on the annotation-only key: top score (rule = ANM) 0.984 [0.968, 0.995], margin 0.993 [0.986, 0.997], entropy 0.993 [0.989, 0.996], classifier 0.983 [0.973, 0.991]. The rule's realised Q1 coverage c* is 0.993 [0.990, 0.995]. Per donor (P1–P8), the primary E1.4a difference is −0.003, −0.002, −0.000, −0.001, −0.036, −0.003, −0.004 and −0.022.

- **What the primary rows say:**
  - **E1.4a.** The direction "TEDDY's margin ranks the rule's calls at least as well as the rule's (= ANM's) top score" holds on the external data. The margin is ahead in all 8 donors, and the pooled interval lies below 0.
    - Under E1's own decision rule this is still inconclusive. The pooled point passes the 0.005 margin, but several donors fall short of half of it (for example P3 −0.000 and P4 −0.001).
    - On site4 the same endpoint was inconclusive, −0.002 [−0.024, 0.020]. What replicates is the sign of that point, now with an interval that excludes 0.
    - This supports keeping the trust-gate claim withdrawn. It does not make ANM worse than the rule: top score (ANM) − rule is 0.0000, by construction.
  - **E1.C2.** E1's site4 point at c* is exactly 0, so the replication rule gives "no v3 direction". Externally the difference is 0.0000 as well: equivalent under the E1 rule, with the in-scope guard failing (−0.0086), as on site4.
- *Secondary:*
  - The site4 finding that a classifier on the same 12 evidence values orders calls better than the top score (−0.099) does not replicate: 0.0010 [−0.0102, 0.0091]. On these data every arm's AURC lies between 0.983 and 0.993, which leaves little room for any arm to separate.
  - E1's secondary accuracy gain of the top score over the margin at coverage 0.90 does not replicate: −0.0023 [−0.0046, 0.0008].
  - On the other key (the primary key, which leaves 10,097 cells unscored), the two primary rows have the same status: top score − margin −0.0100 [−0.0215, −0.0016] (replicates) and C2 at c* 0.0000 [0.0000, 0.0047] (no v3 direction). One secondary row differs: C2 at 0.70 replicates on this key (0.0311 [0.0010, 0.0521], 6 of 8 donors), while the E1 rule gives inconclusive with the in-scope guard failing (−0.0069).
  - Reduced panels (sensitivity): ANM's engine equals its own closed form (the bridge check passes), but with unequal panel sizes it no longer equals the mean rule (123 Q1 calls differ). ANM's top score − margin, both on the reduced panels, is −0.0116 [−0.0178, −0.0056] (replicates); C2 at c* is 0.0011 [0.0000, 0.0051] (no v3 direction).
  - Deployed operating points (bars from val, annotation-only key, descriptive): the Q1 rule has coverage 0.993, selective accuracy 0.966 and OUT decline 0.292. ANM's closure readout: 0.992 / 0.959 / 0.351. The classifier: 0.998 / 0.967 / 0.103. The rule calls 442 of 486 platelets myeloid, because CD62P is in the myeloid panel.
- **Limitations.**
  - RNA possibly seen by TEDDY in pretraining. The addendum's answer is "likely yes": the dataset was public in CELLxGENE before TEDDY's training data were downloaded. Every arm reads the same embedding, so exposure does not favour one arm, but absolute accuracies may be optimistic.
  - Key not validated (pooled kappa 0.6562; 0.4314 in donor P7). The annotation-only key decides.
  - The data are easy for this pipeline (c* 0.993; every arm's AURC between 0.983 and 0.993), which limits what any comparison can show.
  - The annotation map's first version was revised after the prep step had read external values, so the registered order "committed before any external value is read" was not met in time. An independent re-derivation from the label names alone reproduced the committed map (31 of 31 labels).
  - The version-2 panel revision was decided by the orchestrating session, which had seen site4 results of E1, E3 and E4 but no external value (addendum disclosure).
  - 8 donors; the intervals describe these donors.

> **中文：** E6：所有结果带两个标签——"RNA 可能已被 TEDDY 预训练见过"和"答案键未验证"（kappa 0.6562 < 0.70，故以仅注释键判定）。E1.1a 在每个外部供体上复现（0 不一致）。主要终点 E1.4a（顶分 − margin 的 AURC）为 −0.0087 [−0.0182, −0.0017]，8 个供体全为负，按复现规则"复现"：TEDDY 的 margin 对调用的排序至少与规则（即 ANM）的顶分一样好；按 E1 自身判定规则仍为不确定。C2 在 c* 处为 0.0000，"无 v3 方向"（等价）。分类器排序优势未复现（外部数据过易，各方法 AURC 在 0.983–0.993 之间）。

### E7: the paper's inverse loop on TEDDY, at the layer-11 cut, with the consumer fixed

Addendum `E7.json` version 2 (sha256 `c35e91c5…`; registration page `registration/E7_REGISTRATION.md`) was committed at 22:24:05 UTC and pushed at 22:25:53 UTC, before E7's first development forward pass (22:26:38). The development selection records (dev, then the declared follow-up; both select no candidate) were committed at 22:33:13 and 22:40:18 and pushed together at 22:40:19. The site4 confirmation started at 22:40:20 by the local clock (section 2.5). The result file records the addendum and selection-record hashes, `committed: true` for both, code and a head checkpoint that match the addendum, and git HEAD `871b709` (the pushed commit).

- **Question, in plain words.** Our head reads the mean of TEDDY's last-layer token states. Go one layer back, to the token states that leave layer 11. Is their mean enough to fix what TEDDY's layer 12, mean pooling and our frozen head then predict for the NK/T panel? E7 asks this the way the ANM paper's v2 revision (PR #16, draft, at `aa9e2c20`) tests a description of a model, in three steps:
  1. **Reject.** Find histories with exactly the same candidate state but readouts that differ beyond the tolerance.
  2. **Select.** Choose a revised description from a declared candidate library, by a declared rule, on development data.
  3. **Validate.** Test the choice on fresh comparisons, with the consumer and interface unchanged.

  PR #16's discussion table states the template: "Exactly equal retained states give a readout difference beyond measurement error" → "Reject sufficiency for these histories under the declared conditions" → "Separate the violating histories in a revised representation; test it with consumer and interface unchanged" (`sections/discussion.tex` lines 56–59). It adds that "Selection operates within a supplied candidate library" (`sections/inverse_model.tex` lines 201–202). In its broader development comparison "No candidate qualifies" (line 39); its results table then gives each candidate its own outcome: "Rejected" for the snapshot descriptions, "Bounded support" for the association description and "Uninformative: no matched-history comparison" for full history (lines 176–184).
- **Design** (addendum version 2):
  - **Consumer, fixed and never trained.** Layer-11 token states → TEDDY's own layer 12 (E5's explicit-attention layer, equal to the module's) → mean over the cell's real tokens → our frozen head (L2 normalisation, MLP, negative-binomial mean at the training-median size factor, ÷ the training 95th percentile). Layer 12 has no positional input (checked on training cells), so the consumer sees only the set of layer-11 token states.
  - **Readout.** The head's predictions of CD56, CD94, CD335 and CD3, E3's gap proteins. The difference between two histories is the largest over the 4 proteins of |change| ÷ u_p. Here u_p is the measured NK–T gap of protein p on the val donor (18303: 254 NK and 2,990 T cells): 1.563, 1.614, 1.510 and 1.089. **Tolerance: 0.05**, that is, 5% of the measured NK–T difference of a protein.
  - **Candidate library**, in declared order:
    1. the layer-11 gene-mean;
    2. the gene-mean plus the layer-11 states of the cell's tokens of G = {NCAM1, KLRD1, NCR1, CD3E, CD3D, CD3G}, the coding genes of the 4 panel proteins, fixed from the readout alone;
    3. all token states. No two different histories share them, so this candidate has no matched pair and is uninformative by construction.
  - **Histories: exact mean-preserving state interventions at the layer-11 cut.** Disjoint pairs of non-G tokens have their two states replaced by their average, which leaves the gene-mean unchanged (to rounding). Sizes: a share f of the cell's non-G tokens. Development uses f = 1/32, 1/8 and 1/2. Confirmation uses the same sizes with fresh pairs, plus the fresh size f = 1/4.
    - Per size there are 8 matched patches and 8 G twins: a twin uses the same pairs and also averages every G token with its own non-G partner.
    - Per cell there are 8 G-only patches (every G token averaged with a random non-G partner).
  - **Matched pairs** (equal candidate state):
    - for the gene-mean: (unpatched, matched), (unpatched, twin), (matched, twin) and (unpatched, G-only);
    - for the gene-mean + G: (unpatched, matched).
    - Per cell, D = the largest difference over the candidate's matched pairs, across all sizes of the stage.
  - **Witness and controls:**
    - G witness W: pairs with equal gene-mean but different G states, (matched, twin) and (unpatched, G-only).
    - Positive control P: two clamps move every layer-11 token by ±dmu₁₁, the training NK − T difference of the layer-11 gene-mean (‖dmu₁₁‖ = 0.332).
    - Detectability witness R: the same gene-mean shift, delivered only through the tokens of a matched patch at f = 1/2.
    - Implementation checks on every cell: a padding position, a mean-preserving patch of layer-12 states, a permutation of the layer-11 rows, mean preservation, module vs explicit layer, the official gene-mean, "patch applied" and "intended shift". Float64 re-check on the first 10 cells of each roster (no difference may move by more than 0.005).
  - **Development and selection** (val donor 18303 only; 100 NK + 100 T eligible primary-key cells). A candidate is retained if the 95th percentile over cells of its D is at most 0.05, and the minimal retained candidate is selected. If none is retained, one declared follow-up runs on 200 fresh val cells. It restricts the readout to the NK–T score (mean of CD56, CD94 and CD335 minus CD3; unit 2.579), with the same library and rule.
  - **Confirmation.**
    - Site4 (registered): donors 13272 and 19593, 100 eligible primary-key NK and 100 T cells per donor, excluding every E3, E5 and E5-M cell.
    - External (secondary): E6's 8 donors, 25 NK + 25 T cells per donor.
  - **Decision rules.**
    - Candidate k is rejected if the lower end of the 95% interval of the median of D_k is above 0.05, pooled and in each primary donor.
    - Support needs the upper end of the interval of the 95th percentile of D_k at or below 0.05, with P and R both detected.
    - A witness or control is detected if the lower end of the interval of its median is above 0.05.
    - Intervals: two-stage bootstrap, B = 2000, seed 1.
    - LOOP_COMPLETE needs gene-mean + G selected and supported, the gene-mean rejected and W detected. REJECTED_NO_REVISION means nothing was selected and both the gene-mean and gene-mean + G are rejected on fresh comparisons.
- **Development (val only; recorded before any site4 read):**
  - Primary (panel readout, 200 cells): the 95th percentile of D is 0.276 for the gene-mean and 0.256 for the gene-mean + G, both above 0.05. All token states are uninformative. **No candidate selected.** Median W 0.010, P 0.275, R 0.276; every gate passed.
  - Declared follow-up (NK–T score, 200 fresh val cells): 0.236 and 0.236. **No candidate selected.**
  - The registration had named this path in advance, from code tests on training cells, as the likely one: "null primary selection, the declared follow-up, and then confirmation of the rejections only".
- **Registered verdict (site4):** **"REJECTED_NO_REVISION: no candidate was selected in development; on fresh comparisons both the gene-mean and gene-mean + G are rejected"**. The external family returns the same code. Every implementation check passes on every cell; the largest float64 gap is 9.1e-07 (site4) and 3.3e-06 (external).

| statistic (per-cell, in gap units) | site4, 400 cells [95% CI] | donor 13272 | donor 19593 | external, 400 cells [95% CI] | registered role |
|---|---|---|---|---|---|
| median D, gene-mean | **0.125** [0.090, 0.161] | 0.094 [0.086, 0.106] | 0.153 [0.140, 0.167] | 0.112 [0.097, 0.136] | **rejected** (every lower end > 0.05) |
| median D, gene-mean + G | **0.122** [0.086, 0.153] | 0.092 [0.084, 0.103] | 0.148 [0.136, 0.158] | 0.107 [0.091, 0.128] | **rejected** |
| 95th percentile of D, gene-mean | 0.320 [0.256, 0.391] | 0.259 | 0.375 | 0.296 [0.263, 0.345] | support would need an upper end ≤ 0.05 |
| 95th percentile of D, gene-mean + G | 0.315 [0.252, 0.388] | 0.259 | 0.374 | 0.285 [0.242, 0.344] | support would need an upper end ≤ 0.05 |
| median W, G witness | 0.013 [0.009, 0.018] | 0.009 | 0.017 | 0.009 [0.007, 0.010] | not detected |
| median P, positive control | 0.261 [0.222, 0.282] | 0.265 | 0.243 | 0.342 [0.282, 0.400] | detected |
| median R, detectability witness | 0.260 [0.217, 0.276] | 0.261 | 0.244 | 0.337 [0.283, 0.393] | detected |

All token states: uninformative (no matched pair). Site4: 100 cells per donor and class. External: 25 per donor and class; it carries the labels "RNA possibly seen by TEDDY in pretraining" and "key not validated". The external median D for the gene-mean by donor (P1–P8) is 0.108, 0.078, 0.098, 0.096, 0.139, 0.116, 0.128 and 0.158. Every donor's interval lies above 0.05; the lowest lower end is 0.065 (P2).

- **Reading the verdict:**
  - **Step 1 is done on TEDDY.** With the consumer fixed, histories with exactly the same layer-11 gene-mean give readouts that differ beyond the tolerance. This holds on fresh site4 cells, fresh pairs and a fresh size, in both donors and in the external family. So the layer-11 gene-mean is rejected as a sufficient retained state for this consumer and this class of interventions.
  - **Step 2 returned no revision.** The declared augmentation (the gene-mean plus the six marker genes' token states) did not pass in development, and neither did the declared follow-up. On fresh comparisons it is rejected too (0.122). This parallels PR #16's broader development comparison, which "selects no candidate" (`sections/introduction.tex` lines 68–69). There, a recipient-A scope chosen after that comparison then led to a supported description; E7 did not attempt that.
  - **Step 3 is not reached.** E7 does not validate a revised representation.
  - **What the consumer uses.** For this consumer and these interventions, the readout is fixed neither by the layer-11 gene-mean nor by the gene-mean together with the six marker genes' token states. Mean-preserving changes to the other tokens' states, with the marker states untouched, move it beyond tolerance. Changing the marker states at a fixed gene-mean moves it by a median of only 0.013 (the G witness is not detected). What layer 12 and the head use is spread over the other tokens' states: it is not held by the gene-mean, nor by the marker tokens.
  - **The comparisons could see a change.** Shifting the gene-mean by the training NK–T difference moves the readout by a median 0.261 (P). Delivering the same shift through the patched tokens alone moves it by 0.260 (R). The largest patch, averaging half of the non-G tokens in pairs, moves it by about half as much (median D 0.125).
- *Secondary (never change the verdict; 3 decimals from the JSON):*
  - **By size**, median D for the gene-mean on site4 (external in parentheses):

    | share of non-G tokens averaged | median D | 95th percentile of D |
    |---|---|---|
    | 1/32 | 0.022 (0.020), below the tolerance | 0.083 (0.076) |
    | 1/8 | 0.046 (0.046) | 0.138 (0.177) |
    | 1/4 (fresh size) | 0.073 (0.073) | 0.219 (0.210) |
    | 1/2 | 0.125 (0.111) | 0.320 (0.294) |

    The difference grows with the size of the patch. The registered statistic takes each cell's largest difference over all sizes, so the larger patches drive the rejection. At 1/32 the median cell stays within tolerance, but the 95th percentile (0.083) does not. For the gene-mean + G the same row reads 0.019, 0.045, 0.069 and 0.121.
  - **By class** (site4): median D for the gene-mean is 0.141 in NK cells and 0.091 in T cells (external: 0.162 and 0.070). The G witness is 0.016 and 0.008 (external 0.016 and 0.004).
  - **The other readout**, the NK–T score: median D for the gene-mean is 0.060 (site4) and 0.054 (external).
  - **Input push** (descriptive; the E5-M-style push the review asked for): doubling one G gene's input token embedding and running the full forward pass moves the panel readout by a median 0.005 (50 cells per family). It acts at layer 0, so it moves the cut's mean and deviations together and does not decide.
- **Limitations.**
  - One consumer (TEDDY's layer 12 + gene-mean pooling + our frozen head), one cut (layer 11), one readout panel and one class of interventions. Nothing transfers to other layers, consumers, readouts or intervention classes.
  - **The patches are exact state interventions at the layer-11 cut.** They are off-manifold: a token state averaged with another is not a state that a real cell's RNA produces. So they are not biological variation, and the readout differences are not differences between cells.
  - The library is small and declared: the gene-mean, one marker-token augmentation and the uninformative full state. Selection works only within it. A richer library, another cut or a data-driven candidate were not tested.
  - Only 2 primary donors; the intervals describe these donors.
  - The external family is secondary and carries two labels. "RNA possibly seen by TEDDY in pretraining" is E6's contamination label. "Key not validated" is inherited from E6, whose protein-gated key failed its validity check (kappa 0.6562 < 0.70). E7 does not use that key: it picks external NK and T cells by their annotated cell type through E6's committed annotation map and reads no external protein. Its units and tolerance come from the val donor.
  - E7 does not explain or locate our head's NK–T loss (E3), and it says nothing about what the RNA contains. The interventions act on TEDDY's internal states, not on its input.
- **Disclosures** (all recorded in `E7.json`, `declared.pre_registration_code_test` and `declared.review_revision_v2`):
  - The design was revised twice after code tests on training cells, before registration:
    1. After a test on 2 annotated training cells, the first draft's one-sided-twin positive control and twin-calibrated size were replaced by the clamp and the fixed three-size class.
    2. After a test on 3 training cells, the clamp size κ was raised from 0.2 to 1, the full training NK–T difference.
  - It was revised once more, as addendum version 2, after an adversarial review. Before that revision, version 1's code smoke had run on 5 val smoke-pool cells outside both rosters, and its confirmation family on those 5 cells returned REJECTED_NO_REVISION. Version 2:
    - added the detectability witness R and the "patch applied" and "intended shift" gates;
    - made a control failure remove only a support claim, never a rejection;
    - ended E7 if the development record says stopped.
    It changed no tolerance, unit, candidate, roster, size, selection rule, rejection criterion or bootstrap.
  - The registering session had read every published v3 result, including E5-M's site4 per-layer clamp residuals (the head's 0.709 and 0.72 at layer 11).
  - Timing: the dev selection record was committed before the follow-up development run, but pushed only with the follow-up record; both development runs used val cells only. The push of the selection records (22:40:19, GitHub's clock) and the start of the site4 confirmation (22:40:20, local clock) are 1 s apart on two different clocks (section 2.5).

> **中文：** E7 是论文 v2 修订版（PR #16）逆向闭环在 TEDDY 上的类比：在第 11 层切点，消费者固定（第 12 层 + gene-mean 池化 + 冻结的 head），读出为 head 对 CD56/CD94/CD335/CD3 的预测，单位是 val 供体上实测的 NK–T 差距，容差 0.05。历史是在第 11 层切点做的精确保均值干预（把成对 token 状态替换为其平均）。候选库事先声明：第 11 层 gene-mean；gene-mean 加六个标记基因（面板蛋白的编码基因）的 token 状态；全部 token 状态（无匹配对，无信息）。开发阶段（仅 val）两个候选的 95 分位数都超过容差（0.276、0.256），事先声明的后续（NK–T 分数）同样未选出（0.236、0.236）。注册判定（site4）为 **REJECTED_NO_REVISION**：gene-mean 的中位差异 0.125 [0.090, 0.161]（供体 0.094、0.153），加 G 的候选 0.122 [0.086, 0.153]，均被否决；外部家族 0.112 [0.097, 0.136]，8 个供体的区间都高于 0.05。G 见证未检出（0.013），阳性对照（0.261）和可检测性见证（0.260）检出，实现检查与 float64 复核全部通过。差异随补丁幅度增大（1/32 时约 0.02，低于容差；1/2 时约 0.12），NK 细胞大于 T 细胞。由此在 TEDDY 上完成了逆向闭环的第一步（消费者固定、在新比较上否决），第二步在声明的候选库中没有得到修订，类似论文中较广的开发筛选未选出候选；**没有验证任何修订后的表示**。对这个固定消费者和这些状态干预而言，第 12 层 + head 所用的信息既不由第 11 层 gene-mean 承载，也不由六个标记基因的 token 状态承载，而分布在其他 token 上。局限：补丁是切点上的精确状态干预（偏离数据流形），不是生物学变异；只有 2 个主要供体；外部家族带"可能被 TEDDY 预训练见过"和"答案键未验证"两个标签（E7 不使用 E6 的蛋白门控答案键，只按注释细胞类型选细胞，不读外部蛋白）；E7 不解释 NK–T 损失，也不说明 RNA 含有什么。透明说明：注册前在训练细胞上的代码测试后两次修订设计，在 val 冒烟池检查后又修订一次（第 2 版），均记录于 E7.json。

---

## 4. Claims ledger

Status words: **registered** means a registered verdict, check or replication rule decides the claim. **Secondary** means a reported comparison that never decides a verdict. **Descriptive** means a measurement with no decision rule.

### 4.1 What we can claim

| # | claim | status | evidence |
|---|---|---|---|
| C1 | In this setup ANM's engine computes exactly the declared rule, on site4 and on the external data. This is a check, not a gain. | registered (E1.1a pass; E6 E1.1a replicates) | E1.1a: 0 call mismatches for Q1, Q2, Q3 (vs the rule) and the closure readout (vs its re-coded form), test_primary and test_secondary. E6: 0 mismatches in every external donor, 0 rejected events. E2: 0 mismatches over all perturbation cases. E4: engine = closed form, bridge failure False. |
| C2 | ANM's leave-one-out attribution equals the exact closed form. | registered (E1.3a discrepancy False) | E1.3a: top-set and flip agreement 1.0000 on 10,715 Q1 calls. |
| C3 | Pooled over both primary donors, one panel protein decides more calls than the 1/3 null in each class. This is a property of the evidence, and any rule reproduces it. For B it does not hold by donor. | secondary, descriptive | E1.3b: CD11c 0.832 of myeloid, CD2 0.618 of T, CD94 0.572 of NK, CD20 0.484 of B calls; null 1/3, p = 0.0010. In donor 19593, CD20 0.339 (p = 0.6494) and CD22 0.491 of B calls. |
| C4 | On site4, a logistic classifier on the same 12 evidence values orders calls better than the rule's (= ANM's) top score. This does not replicate on the external data. | secondary (E1.4); E6: does not replicate | Site4: AURC 0.942 vs 0.843; difference −0.099 [−0.163, −0.034], loss for the top score. E6: 0.0010 [−0.0102, 0.0091], v3 sign in 2 of 8 donors; every arm's AURC lies between 0.983 and 0.993. |
| C5 | A classifier trained only on the old question's labels answers the new CD14-anchored question better than the zero-label rule. | registered (E1.1c loss; E1.5 n*) | E1.1c: −0.078 [−0.134, −0.028], loss for the rule; annotation-only −0.054, loss. E1.5: n* = 0. Not run externally. |
| C6 | Measured operating points of TEDDY + head + Q1 rule. On unseen-donor site4 cells: coverage 0.949, selective accuracy 0.697 (primary key) / 0.730 (annotation-only), OUT decline 0.261. On the external PBMC data (annotation-only key; RNA possibly seen by TEDDY in pretraining): coverage 0.993, selective accuracy 0.966, OUT decline 0.292. | descriptive | E1 and E6 deployed operating points. |
| C7 | TEDDY + head + rule loses decision accuracy as RNA is thinned: 0.022 with 80% of UMIs kept, 0.162 with 20% kept (0.140 with 5% kept). | descriptive (E2's comparator curve) | E2 curve table. |
| C8 | The pipeline's response to these perturbations is not linear in their size. | descriptive (C5, a registered measurement with no decision rule) | E2 C5: at most 0.098 of cells in the linear regime at any stage. |
| C9 | Two trained readouts of frozen TEDDY (a layer-mean MLP and attention pooling) separate NK from T worse than the registered head. | registered (E3.H3b losses) | E3.H3b at 0.9: R1 −0.0409, R2 −0.1925; also losses at 0.8. Descriptive context: a seed-1 retrain of the head scored +0.0547 above it. The head's gap ratio of 0.1169 on flagged pairs vs 0.3985 on unflagged pairs is descriptive only: the intervals overlap, and without gdT CD158b+ pairs the ratios are 0.4729 vs 0.5137. |
| C10 | A val-trained stacker fuses the two channels best, and ANM's fusion loses to it. On the E4 key, the declared trust-weighted average also beats TEDDY + head alone. | registered (F5 − F3 loss); secondary (F2 − F0) | E4: F3 0.9353 vs F5 0.8538. F2 − F0 is +0.0143 [+0.0072, +0.0269], a win on the E4 key, but inconclusive on the v3 primary key (+0.0073) and equivalent on test_secondary (−0.0030). |
| C11 | The registered E5 derivative check failed in float32 at ε = 1e-3. The registered verdict attributes this to float32 rounding, because the float64 version passes. | registered (E5 gate) | E5: 0.854 < 0.95; S5 1.000. |
| C12 | On the external data, TEDDY's margin ranks the declared rule's calls at least as well as the rule's (= ANM's) top score (the replicated E1.4a direction): it is ahead in all 8 donors and the pooled interval lies below 0. On site4 the same comparison was inconclusive. | registered (E6 primary E1.4a: replicates); E1 decision rule over the 8 donors: inconclusive | Site4 E1.4a: −0.002 [−0.024, 0.020], inconclusive. E6: −0.0087 [−0.0182, −0.0017], v3 sign in 8 of 8 donors. Key not validated; RNA possibly seen by TEDDY in pretraining. |
| C13 | For E3's attention-pooling readout R2, at first order on 195 site4 cells, the median share of its response that flows through what gene-mean pooling discards is 0.085 [0.071, 0.096] on NK–T look-alike pair cells and 0.127 on matched random cells: smaller, not larger, on look-alikes. | registered measurement; the registered verdict on H5M is inconclusive (E5-M) | Δ_tok −0.043 [−0.062, −0.011], per donor −0.048 and −0.029; DD_rand −0.041 [−0.061, −0.004]. Random token readout: 0.067 / 0.068. Derivative gate validated (pass share 1.000); positive control passed. For this observer only; its attention is nearly uniform. |
| C14 | ANM's own factorization check rejects the gene-mean as an *exact* first-order retained state for R2, and its source-kernel test is uninformative in this setting, as the ANM paper predicts. | secondary (E5-M S4) | rtol 1e-4: REJECTED in 195 of 195 cells; source kernel NOT_INFORMATIVE in 195 of 195 (K_z has full column rank in every cell; near-null share 0.000). |
| C15 | For a fixed consumer (TEDDY's layer 12 + gene-mean pooling + our frozen head) and exact mean-preserving state interventions at the layer-11 cut, the layer-11 gene-mean is rejected as a sufficient retained state on fresh site4 comparisons. So is the gene-mean plus the token states of the six marker genes. No declared candidate was selected in development. This is the first step of the paper's inverse loop on TEDDY, with no revision. | registered (E7: REJECTED_NO_REVISION) | Site4 median D 0.125 [0.090, 0.161], per donor 0.094 [0.086, 0.106] and 0.153 [0.140, 0.167]; gene-mean + G 0.122 [0.086, 0.153]; tolerance 0.05 of the measured val NK–T gap. Development 95th percentiles 0.276 / 0.256 (panel) and 0.236 / 0.236 (score follow-up), above 0.05. G witness 0.013, not detected; positive control 0.261 and detectability witness 0.260, detected; implementation and float64 checks pass. External (secondary; labels): 0.112 [0.097, 0.136], every donor's interval above 0.05. |
| C16 | In E7 the readout difference grows with the size of the mean-preserving patch, and it is larger in NK than in T cells. | secondary (E7) | Site4 median D for the gene-mean: 0.022 at 1/32 of the non-G tokens (below the tolerance), 0.046 at 1/8, 0.073 at 1/4, 0.125 at 1/2. NK 0.141 vs T 0.091 (external 0.162 vs 0.070). |

### 4.2 What we cannot claim (and which test stops it)

| # | claim we cannot make | status | which test |
|---|---|---|---|
| N1 | "ANM adds decision value as a decision layer on TEDDY (Mode B)." | registered: falsified (E1) | E1 falsification: neither E1.4 nor C2 is a win; every ANM arm equals a re-coded rule. E6: the engine equals the rule in every external donor, and the margin ranks calls better than the top score in 8 of 8 donors. |
| N2 | "ANM's confidence is a better trust gate than TEDDY's margin." | registered: inconclusive (E1.4a); E6: the direction against this claim replicates | E1.4a: −0.002 [−0.024, 0.020]; the gate claim stays withdrawn. E6: −0.0087 [−0.0182, −0.0017]. The site4 0.90-coverage win is secondary, belongs to the rule's top score, and does not replicate externally (−0.0023 [−0.0046, 0.0008]). |
| N3 | "ANM's nested (closure) readout declines out-of-scope cells better." | registered: equivalent (E1 C2); E6: no v3 direction, equivalent | C2 at c*: 0.000 [−0.006, 0.003]; in-scope guard failed (−0.056). E6: 0.0000 [0.0000, 0.0000]; guard failed (−0.0086). |
| N4 | "Editing the question with zero labels beats training" / "ANM saves labels." | registered: loss (E1.1c, E1.1d); not run externally | E1.1c/E1.1d: loss (−0.078); E1.5: n* = 0. Zero labels holds for any written rule. |
| N5 | "Attribution flags wrong calls better than TEDDY's margin." | registered: inconclusive (E1.3c) | E1.3c: 0.004 [−0.008, 0.022], inconclusive (Q2 equivalent). |
| N6 | "The representation / head / decision decomposition tells more than the accuracy curve." | registered: falsified (E2) | E2: adds nothing beyond the curve; criterion (ii) is a loss (−0.092). |
| N7 | "The NK–T loss is in the readout or pooling and is repairable on frozen TEDDY." | registered: rejected (E3, A2.1) | E3 falsification: rejected; H3a inconclusive for R1 and R2; H3b losses. |
| N8 | "The NK–T loss is in mean pooling" / "pooling discards NK–T information that a token readout uses on look-alike cells." | registered: inconclusive (E3 pooling; E5-M) | E3.H3a pooling (A3.1): −0.0375 [−0.4408, 0.22]. E5-M: H5M not supported; for R2 the look-alike excess has the opposite sign (Δ_tok −0.043 [−0.062, −0.011]), and DD_rand's lower bound is below 0 (−0.061). |
| N9 | "ANM's field fuses evidence channels better." | registered: loss ("field hurts") and inconclusive (E4) | E4: F5 − F3 −0.0815, loss; F5 − F2 inconclusive (+0.0034); C2 (t = 2) equivalent. |
| N10 | "A specific TEDDY layer localises the NK–T loss." | registered: JVP not validated (E5) | E5 stops; the computed l* = 11 does not agree with the probe layer 9 either. |
| N11 | Any population-level statement. | design | Two site4 donors and 8 external donors; the intervals describe these donors. |
| N12 | Any v1/v2 headline (for example the first run's ANM accuracy edge, or "0 vs 50–1,000 labels") as a v3 result. | superseded | Superseded by E1 under the v3 key; see section 1. |
| N13 | "TEDDY's gene-mean state is sufficient for an independent observer on look-alikes" / "pooling is not the problem; the NK–T loss is in our head." | registered: PASS_SUFFICIENT not reached (E5-M inconclusive) | The pair-cell level is 0.085 [0.071, 0.096], above 0.05, and Δ_tok is not equivalent. Measured for one near-uniform observer only. |
| N14 | "A trained classifier orders calls better than the rule in general" (C4 beyond site4). | secondary: does not replicate (E6) | E6: 0.0010 [−0.0102, 0.0091]. |
| N15 | "E6 confirms v3 on data TEDDY never saw, with a validated key." | registered labels (E6 contamination and key-validity rules) | Contamination answer "likely yes"; key not validated (0.6562 < 0.70). |
| N16 | "E7 validates a revised TEDDY representation" / "the reject → select → validate loop is complete on TEDDY" / "the six marker genes' token states carry what the consumer needs". | registered: REJECTED_NO_REVISION (E7) | No candidate selected in development or in the follow-up; gene-mean + G rejected on fresh comparisons (0.122 [0.086, 0.153]); the G witness is not detected (0.013 [0.009, 0.018]). LOOP_COMPLETE is not reached. |
| N17 | "E7 explains or locates the NK–T loss", says what the RNA contains, or describes biological variation between cells. | design (E7) | The patches are exact state interventions at the layer-11 cut (off-manifold), for one fixed consumer; E7 never changes the input RNA and has no repair or localisation endpoint. |

### 4.3 What was not tested

| # | not tested | status | why |
|---|---|---|---|
| T1 | ANM's field dynamics: retention over real processing steps, multi-step propagation, histories, feedback. | not tested | Every event enters at t = 0, so the field acts as a constant gain per class (G(3)·3·S). E4's t = 2 control is a declared order on one channel, not a history. |
| T2 | Sufficiency of TEDDY's gene-mean for other observers (sharper attention, nonlinear readouts), for larger pushes, and on more cells. | partly tested: one observer, first order, 195 cells (E5-M); one consumer at the layer-11 cut under finite mean-preserving state interventions, 400 site4 cells (E7: rejected) | R2's attention is nearly uniform (median effective number of tokens 1,141.957 of a median 1,368). E7 tests the layer-11 gene-mean for layer 12 + our head, not the layer-12 gene-mean for another observer. |
| T3 | Replication on data TEDDY has certainly not seen, and external replication of E1.1b–d, E1.3, E1.5 and of E2–E5-M. | partly tested (E6; E7's external family) | E6's data are labelled RNA possibly seen by TEDDY in pretraining; E6 ran E1.1a, E1.4a, E1.C2 and Q1 accuracy only. E7 ran on the same donors as a secondary family. |
| T4 | Pretraining contamination of TEDDY with GSE194122 (this dataset). | not verified | Registration section 13. |
| T5 | Effect of the head's global ADT-median size factor (computed over all cells, site4 included) on its weights. | not testable here | Needs retraining. |
| T6 | Retraining TEDDY or the head; other heads, panels or questions beyond the registered ones. | out of scope | TEDDY and the head are frozen. |
| T7 | Comparison with TypeSafe AI's Jev. | not tested | Never called; the classifier is a stand-in. |
| T8 | Validated layer-wise linear response in TEDDY (E5's question). | not tested at the registered gate | E5 stopped at its registered gate; the float64 path passed but was not registered before site4. E5-M's own gate (ε 0.01 and 0.02, Richardson) passed, but it validates E5-M's derivatives and does not answer E5's layer question. |
| T9 | The validation step of the inverse loop on TEDDY: a revised representation selected and then supported on fresh comparisons. Also other cuts, other consumers, a richer candidate library and on-manifold (real-cell) perturbations. | not tested | E7's declared library had no candidate that passed in development, so nothing reached validation. E7 used one cut, one consumer and one class of state interventions. |

> **中文：** 账本（每条都有状态栏）：能说的是 ANM 在 site4 与外部数据上都精确复现声明规则、在外部数据上 TEDDY 的 margin 对调用的排序至少与顶分一样好（8/8 供体；site4 上不确定）、R2 观察者在相似细胞上流经池化丢弃部分的份额反而更小（E5-M，判定不确定，仅限该观察者）、在第 11 层切点且消费者固定时，第 11 层 gene-mean 以及加上六个标记基因 token 状态的候选都被否决、开发阶段未选出任何候选（E7，REJECTED_NO_REVISION）等；不能说 ANM 提升决策、信任门控、嵌套读出、节省标签、融合、层定位，也不能说"池化丢失信息"或"基因均值状态充分"，不能说 TEDDY 上的逆向闭环已完成、某个修订表示已被验证，或 E7 解释了 NK–T 损失；未检验的是 ANM 的场动力学、其他观察者、逆向闭环的验证步骤、真正未见过的数据及污染问题。

---

## 5. What ANM brings to TEDDY, honestly

**The paper's framing.** Paper line numbers in this section refer to the ANM v2 pre-release draft (commit `bab56c34`, not yet released; `paper/arxiv/open_to_public_v2/`). In ANM v2, the retained state is a *declared candidate* tested against a declared readout: "a candidate retained coordinate … Its choice does not establish sufficiency" (`sections/protocol.tex` lines 7–11). The readout-visible quotient is a diagnostic in the SI, not "the state ANM finds". For the AI studies the paper says: "Fixed graph updates propagate and retain these increments; a separate model or program reads the resulting state. The delayed-tool study instead retains request–return histories." (`sections/ai_evidence.tex` lines 7–9; PR #16, at `aa9e2c20`, no longer has this passage). It adds that "the experiments do not require or identify a state jointly formed by interacting live agents" (`sections/introduction.tex` lines 45–46). *Our reading, not the paper's words:* these AI workflows are open-loop in the sense that readouts do not feed back into the graph state or its sources (the constructed graph itself has reciprocal couplings). "Open-loop" (readouts do not feed back) is our label: neither the v2 pre-release draft nor PR #16 calls these workflows open-loop or feedback-free. Earlier versions of this report also called ANM a "verifier" of declared responses; that was our label too, and it understated what ANM measures (next paragraph). The paper's "verifier" is a component inside some studies, e.g. "verifier-selected values" (`sections/ai_evidence.tex` line 67).

**What ANM measures (PR #16, draft, at `aa9e2c20`; line numbers in this paragraph refer to PR #16).** ANM's output is a response chain, not decisions. Push one source by a known amount and measure separately how the retained state and the readout move: "Separate source-to-state and source-to-readout measurements connect forward analysis of propagation and readout sensitivity with inverse tests of what information a representation must preserve" (`abstract.tex` lines 6–8). This locates where information is lost: not carried by the state, or carried but not used by the readout (forward), and which distinctions the state must keep (inverse). The source side is an intervention, so the measured source → readout response is interventional. The paper does not establish that the readout change is mediated by the state: "the dashed state-to-readout relation is a hypothesis, not an established causal path" (`sections/introduction.tex` line 85); "State measurement and causal intervention answer different questions" (`sections/protocol.tex` line 129); "causal mediation, autonomous reduction, cross-model transfer and jointly formed multi-agent states require separate evidence" (`sections/discussion.tex` lines 96–98). For a declared observer (readout, perturbations, window) the same responses test and select the retained state. A candidate is rejected when identical states give resolved differences in output distributions (`abstract.tex` lines 9–11). Then, within a candidate set declared in advance, "selection retains candidates with no observed same-state violation and chooses the fewest distinct discovery states" (`sections/worked_example.tex` lines 105–107); the selected description is fixed and checked on new comparisons (lines 4–10). In the paper's exact construction this gives the coarsest sufficient representation, with two states (line 71; the abstract: it "identifies a sufficient state retaining sensor–recipient associations"). Limits: the candidate set is declared in advance, so the procedure selects rather than discovers new variables ("within a specified candidate set", line 6); `sections/protocol.tex` lines 7–11 still say that choosing a candidate "does not establish sufficiency"; in the real-model study below the broader screen selected none, and support came only after the readout was restricted to recipient A. On TEDDY, E7 ran exactly this reject → select → validate loop (E7's rule: keep a candidate whose development differences stay within the tolerance and select the smallest) and selected nothing: the gene-mean was rejected (median D 0.125 [0.090, 0.161]), so was the gene-mean plus six marker genes' token states (0.122), and development and the declared follow-up selected none. These PR #16 lines are unchanged at its current head `9d5e4909`.

**The paper's v2 revision (PR #16, draft).** PR #16 (draft, at `aa9e2c20`; not merged) adds one restricted inverse test with a fixed model consumer. On recipient A only, 768 new calls (512 base + 256 control) reject a snapshot description: histories with identical snapshots but opposite routing differ by at least 0.8358 in their largest probability difference (total variation). An association description, chosen from a declared candidate list, gets bounded support on four new histories and two continuation combinations: every association-matched comparison has a total-variation upper bound at most 0.0821, below the fixed margin 0.15. The broader development screen selected no candidate; support came only after the readout was restricted to recipient A, a scope chosen after that screen. Every one of the 768 answers matched its reference answer, and the bound counts answers at the reference, so support there coincides with the model answering correctly (our observation; the paper does not say this). The consumer is unchanged, so this validates the analyst's description: "not a repair of the model or evidence of intervention superiority" (PR #16, `sections/inverse_model.tex` lines 204–205).

> **中文：** 以下论文行号均指 ANM v2 待发布稿（commit `bab56c34`，尚未发布）。论文对 AI 研究的原文是："Fixed graph updates propagate and retain these increments; a separate model or program reads the resulting state. The delayed-tool study instead retains request–return histories."（`sections/ai_evidence.tex` 7–9），并说实验 "do not require or identify a state jointly formed by interacting live agents"（`sections/introduction.tex` 45–46）。"开环，指读出不会反馈到图状态或源（构造的图内部本身有双向耦合）"是我们的解读，不是论文原话；本报告早先称 ANM 为声明响应的"验证者"，那也是我们的标签，低估了 ANM 测的东西。ANM 测的是一条响应链（PR #16）：把一个源推动已知的量，分别测量保留状态和读出怎么变，从而看出信息在哪里丢失（状态没带上，或带上了但读出没用上），以及状态必须区分什么。源一端是干预，但论文没有证明读出的变化经由该状态传递（PR #16 `sections/introduction.tex` 85、`sections/protocol.tex` 129、`sections/discussion.tex` 96–98）。声明观测者（读出、扰动、时间窗）之后，同样的响应还能检验并选择保留状态：同一状态下输出分布有可分辨差异的候选被否决；在事先声明的候选集合里，保留没有违例的候选，选出状态数最少的那个，固定后用新的比较验证（PR #16 `sections/worked_example.tex` 4–10、105–107；精确构造中得到最粗的充分表示，2 个状态，71 行）。候选集合是给定的，所以这是选择而不是发现新变量；在真实模型研究中，较广的筛选什么也没选出，只有把读出限定为接收者 A 后才得到支持；TEDDY 上的 E7 运行了这一闭环，什么也没选出。论文 v2 修订版 PR #16（草稿，位于 `aa9e2c20`）新增一项范围受限、消费者固定的逆向检验：仅针对接收者 A，768 次新调用（512 次基础 + 256 次对照）否决了快照描述（快照相同、路由相反的历史，最大概率差的下界至少 0.8358）；从声明候选中选出的关联描述，在四个新历史和两个后续组合上获得有界支持（每个关联匹配比较的上界都不超过 0.0821，低于事先固定的容差 0.15）。较广的开发筛选未选出任何候选，"仅接收者 A" 的范围是在那之后确定的。768 个回答全部等于参考答案，而该界只计参考答案处的回答，所以那里的"支持"与模型答对重合（这是我们的观察，论文没有这样写）。消费者不变，验证的是分析者的描述，论文明确说这不是对模型的修复。

**What that means in teddy_mm v3, across all eight experiments.**

1. **Checks that the declared readouts are what the engine computes.** Each readout was stated as a declared object and checked against ANM's engine cell by cell:
   - E1.1a: 0 mismatches on site4. E6: 0 mismatches in every external donor.
   - E1.3a: 1.0000 agreement after A1.3. Reading the engine's leave-one-out, which deletes an event's site, showed that the registered closed form re-coded zeroing an event, not deleting it.
   - E4: engine = closed form.
   These checks passed. They are bookkeeping: they show ANM computes the declared rules, not anything new about TEDDY. The paper's response measurements (ANM "organizes response theory around three measurable elements", `sections/introduction.tex` line 20) correspond here to E2, E5 and E5-M (point 3). The inverse sequence that PR #16 adds (reject a description; within the declared candidates, select the one with no same-state violation and the fewest states; test it on new comparisons with the consumer fixed) was run on TEDDY in E7 (point 3): its first step held, its second selected nothing, and nothing reached the third.
2. **The decision-gain claims of earlier versions did not hold, by construction here.** These tests checked claims that earlier versions of these pages made; they are not ANM's purpose. All events enter at t = 0, so ANM's decision-layer calls equal a re-coded declared rule.
   - Where its readout forms (top action score, closure readout, contradiction fusion) were tested, they did not beat TEDDY's margin, the mean rule, a trained classifier or a val-trained stacker at any registered primary endpoint (E1, E4). The one secondary site4 gain over the margin (top score at coverage 0.90, +0.021) belongs to the rule's top score and did not replicate externally.
   - On the external data, TEDDY's margin ranked the rule's calls better than the top score in all 8 donors (E6: −0.0087 [−0.0182, −0.0017]; still inconclusive under E1's own decision rule).
   - ANM's fusion did beat channel 1 alone on the E4 key (F5 − F0 +0.0177, win), but so does the declared trust-weighted average (F2 − F0 +0.0143, win), and F5 − F2 is inconclusive. Calling this "TEDDY + ANM decides better" would be wrong.
   - Where ANM's field does not run, the method is a declared rule.
3. **The response-theory tools gave one validated measurement about TEDDY (E5-M) and one registered rejection (E7). Neither locates the NK–T loss.**
   - E2's decomposition added nothing beyond the accuracy curve.
   - E5's registered derivative gate failed in float32, so its layer analysis stopped.
   - E5-M ran ANM's sufficiency machinery on TEDDY end to end: the gene-mean as a declared candidate, ANM's clamp estimator, declared observers, a validated derivative gate (pass share 1.000 for every target) and a passing positive control. For R2 at first order, 0.085 of its response on look-alike pair cells flows through what pooling discards, against 0.127 on random cells. The registered verdict is inconclusive. The look-alike excess that H5M predicted is absent (Δ_tok −0.043 [−0.062, −0.011]), but the gene-mean is not shown sufficient either, because the pair-cell level is above 0.05.
   - ANM's own factorization check rejects exact sufficiency at rtol 1e-4 in every cell (any relative residual above 1e-4 rejects, so it says nothing about size). As the paper predicts, its source-kernel test is not informative when m ≤ 21 pushes meet d = 512 dimensions. The informative quantity is the clamp residual.
   - E7 ran PR #16's inverse sequence at the layer-11 cut, with the consumer fixed (TEDDY's layer 12 + gene-mean pooling + our frozen head) and exact mean-preserving state interventions as histories.
     - **Step 1 held.** The layer-11 gene-mean is rejected on fresh site4 comparisons: median D 0.125 [0.090, 0.161] against a tolerance of 0.05, with the positive control and the detectability witness detected.
     - **Step 2 returned no revision.** Nothing in the declared library passed in development, nor in the declared follow-up, and the marker-token augmentation is rejected on fresh comparisons too (0.122 [0.086, 0.153]).
     - In PR #16's terms, these parallel its broader development comparison, where "No candidate qualifies" (`sections/inverse_model.tex` line 39), and its snapshot descriptions, which are "Rejected" (lines 178–180). **No revised representation is validated.**
   - So we now know how much of one near-uniform observer's first-order response passes through what pooling discards. We also know that, for one fixed consumer and exact state interventions at the layer-11 cut, neither the layer-11 gene-mean nor the gene-mean plus the six marker genes' token states fixes what layer 12 and the head read out. We do not know whether the gene-mean is sufficient for the NK–T distinction in general, or whether the NK–T loss sits in our head or in pooling. E7's interventions are off-manifold state changes, so they do not answer that.
4. **What TEDDY users can take away** (not ANM-specific):
   - The pipeline's measured operating points under a five-class key with an explicit validity check. That check flagged donor 13272 on site4, and the key on the external data (kappa 0.6562).
   - On site4, a trained classifier on the same evidence beats the declared rule even for a changed question. Its call-ordering advantage did not carry over to the external data, where every arm's AURC lies between 0.983 and 0.993.
   - Ranking the rule's calls by TEDDY's own margin instead of the top score: on site4 the AURC comparison was inconclusive (−0.002 [−0.024, 0.020]). On the external data the margin ranked the calls better in 8 of 8 donors (−0.0087 [−0.0182, −0.0017]), which is still inconclusive under E1's own decision rule.
   - Swapping in a layer-mean MLP or an attention-pooling readout of frozen TEDDY did not improve NK-vs-T accuracy over the registered head (E3.H3b losses). Whether the head compresses the NK–T gap more on look-alike pairs is not established. The pooled point estimates (0.1169 flagged vs 0.3985 unflagged) have overlapping intervals and come mostly from gdT CD158b+ pairs in donor 13272. Without those pairs the ratios are 0.4729 vs 0.5137.
   - For the attention-pooling readout, look-alike cells did not depend more, at first order, on what pooling discards than random cells did (median share 0.085 vs 0.127). This is a measurement on one observer and 195 cells; the registered E5-M verdict is inconclusive.
   - Summarising TEDDY's layer-11 token states by their mean loses what layer 12 and our head use, under exact mean-preserving state interventions (E7: rejected on site4; the secondary external family agrees). This is a property of TEDDY's last layer read by our head, not a statement about real cells.
   - The answer key itself is fragile for one site4 donor and on the external data.

> **中文：** ANM 测的是一条响应链：把一个源推动已知的量，分别测量保留状态和读出怎么变，从而看出信息在哪里丢失；声明观测者之后，同样的响应还能检验并选择候选保留状态（源一端是干预，但没有证明读出的变化经由该状态传递）。论文 v2 修订版 PR #16（草稿，位于 `aa9e2c20`）新增的逆向检验顺序（否决描述 → 在声明候选中选出没有违例、状态数最少的那个 → 在新比较上检验，消费者固定）已由 E7 在 TEDDY 第 11 层切点上运行：第一步成立（消费者固定时，第 11 层 gene-mean 在新的 site4 比较上被否决，中位差异 0.125 [0.090, 0.161]，容差 0.05），第二步在声明的候选库中没有得到修订（加六个标记基因 token 状态的候选在开发阶段未通过，在新比较上也被否决），因此没有验证任何修订后的表示。在 teddy_mm 的八项实验中：引擎与声明读出逐细胞一致（site4 与外部数据皆 0 不一致），这只是核对；决策增益检验针对的是早期页面的说法，不是 ANM 的目的，所有事件在 t = 0 时这些说法不成立，外部数据上 TEDDY 的 margin 在 8/8 供体中排序优于顶分（按 E1 判定规则仍为不确定）；E2、E5 未得出经验证的结论；E5-M 首次用 ANM 的充分性工具得到一项经验证的测量（R2 响应在相似细胞上 0.085、随机细胞上 0.127 流经池化丢弃部分），但注册判定为不确定，既不能说池化丢失信息，也不能说基因均值状态充分，且仅针对这一个观察者；E7 的干预是切点上的精确状态干预（偏离数据流形），不能用来解释 NK–T 损失。

---

## 6. Overall conclusion

Eight pre-registered experiments on frozen TEDDY tested the decision-gain claims that earlier versions of these pages made, and found no decision gain from ANM. E1–E5-M and E7 were evaluated on the two site4 donors that appear in no other split. E6 repeated E1's main endpoints on 8 external donors whose RNA TEDDY has possibly seen in pretraining, with an answer key that did not pass its validity check. With all events at t = 0, ANM's engine reproduces a declared rule cell by cell, on site4 and on the external data. At the registered primary endpoints none of its readouts beat TEDDY's own margin, a trained classifier or a validation-trained stacker. On the external data TEDDY's margin ranked the calls better than the rule's (= ANM's) top score in 8 of 8 donors, although under E1's own decision rule this is still inconclusive. ANM's response-theory tools produced one validated measurement: for one trained token readout at first order, a median 0.085 of its response on NK–T look-alike cells flows through what gene-mean pooling discards, less than the 0.127 on matched random cells. The registered verdict on that test is inconclusive, so neither "pooling loses the look-alike difference" nor "TEDDY's gene-mean state is sufficient" is established, and the measurement holds for that one observer only. E7 then ran the inverse loop of the ANM paper's v2 revision on TEDDY, at the layer-11 cut, with the consumer fixed (TEDDY's layer 12, gene-mean pooling and our frozen head). Histories with exactly the same layer-11 gene-mean gave readouts that differ beyond the tolerance on fresh comparisons: a median 0.125 [0.090, 0.161] of the measured NK–T gap, against 0.05. So did histories that also share the six marker genes' token states, and no candidate was selected in development. The rejection step of that loop therefore holds on TEDDY, the selection step returned no revision from the declared library, and no revised representation is validated. These are exact state interventions, not biological variation, and they do not explain the NK–T loss. ANM's role is to measure the response chain and to test and select candidate retained states for a declared readout, not to make decisions. In this project that gave one first-order measurement for one observer (E5-M, inconclusive) and one reject → select → validate run that rejected the layer-11 gene-mean and selected no revision (E7). The decision-gain verdicts test claims of earlier versions; they are not ANM's purpose.

> **中文：** 在冻结的 TEDDY 上进行的八项预注册实验检验了早期页面的决策增益说法，没有发现 ANM 带来决策增益。E1–E5-M 与 E7 在两个未出现在其他划分中的 site4 供体上评估；E6 在 8 个外部供体上重复了 E1 的主要终点，这些供体的 RNA 可能已被 TEDDY 预训练见过，且答案键未通过有效性检验。所有事件在 t = 0 时，ANM 引擎在 site4 和外部数据上都逐细胞复现声明规则。在注册的主要终点上，其读出均未胜过 TEDDY 自身的 margin、训练分类器或验证集堆叠器；在外部数据上，TEDDY 的 margin 在 8 个供体中全部比规则（即 ANM）的顶分排序更好，但按 E1 自身判定规则仍为不确定。ANM 的响应理论工具给出了一项经验证的测量：对一个训练好的 token 读出，在一阶下，NK–T 相似细胞上其响应的中位份额 0.085 流经基因均值池化所丢弃的部分，少于匹配随机细胞上的 0.127。该检验的注册判定为"不确定"，因此既不能说"池化丢失了相似细胞的差异"，也不能说"TEDDY 的基因均值状态是充分的"，且该测量只针对这一个观察者。随后 E7 在 TEDDY 第 11 层切点上运行了论文 v2 修订版的逆向闭环，消费者固定（TEDDY 第 12 层、gene-mean 池化和冻结的 head）：第 11 层 gene-mean 完全相同的历史，在新的比较上读出差异超过容差（中位数为实测 NK–T 差距的 0.125 [0.090, 0.161]，容差 0.05）；同时共享六个标记基因 token 状态的历史也是如此，开发阶段未选出任何候选。因此该闭环的否决步骤在 TEDDY 上成立，选择步骤在声明的候选库中没有得到修订，没有验证任何修订后的表示。这些是精确的状态干预，不是生物学变异，也不能解释 NK–T 损失。ANM 的作用是测量响应链、为声明的读出检验和选择候选保留状态，而不是做决策。在本项目中，这给出了一项针对一个观察者的一阶测量（E5-M，不确定），以及一次"否决 → 选择 → 验证"运行：否决了第 11 层 gene-mean，没有选出修订（E7）。决策增益判定检验的是早期页面的说法，不是 ANM 的目的。

---

## Numbers

`docs/reports/v3_numbers.json` is a flat map `{key: {value, source_file, source_key}}` of every number used above and on the site pages (2,359 entries).

- For a `.json` source, `source_key` is an RFC 6901 JSON Pointer, and the value was read from the file by a script, not typed.
- For `README.md` and `REPORT.md` it is a line number.
- For git it is the commit and its committer date.
- For GitHub it is the push event on `exp/registered-v3`.
- For `progress.log` it is a line number with its local time (EDT = UTC − 4).
- For E7, the git committer dates are read from the commits and the `progress.log` times from the log lines by the builder, not typed.
- Every timestamp value is in UTC.

`outputs/…` paths are under `/Users/tianchichen/Documents/GitHub/teddy_mm` and are not committed. `registration/…` and `README.md` are in this repository.

> **中文：** 本文用到的每个数字都以 `{value, source_file, source_key}` 形式收录于 `v3_numbers.json`，供之后更新页面使用。
