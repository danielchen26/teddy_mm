# E7 registration: the TEDDY analog of ANM's inverse loop (layer-11 cut, consumer fixed)

**Status.** Registered before any development or site4 forward pass of its script. The authoritative record is `registration/addenda/E7.json` (sha256 `e14a80ad591edf52da33b37ed7dad68a608efbc7d68f676bffa8c709d278038a`, line in `registration/addenda/HASHES.txt`). It was built by `scripts/mode_a_inverse_e7.py --stage register` from training and validation cells only, and it passed a leakage check. This page explains the design in prose. **Where this page and the JSON differ, the JSON wins.**

It is an addendum to `registration/registration_v3.json` (sha256 `e4c8a33e…`) as amended by A1 (`1f44be66…`), A2 (`20b64174…`) and A3 (`4bb2f856…`). Nothing registered for E1–E6 or E5-M changes.

**Not run.** Development, selection and confirmation have not been run. The script refuses each stage until the records it depends on are committed (see the order of operations below).

---

## 1. Why E7

ANM v2 (PR #16) adds an inverse step with a fixed model consumer, in three parts:

1. Reject a retained-state description with a witness pair: two histories with equal retained state but different readouts.
2. Select a revised description from a declared candidate library, by a declared rule, on development data.
3. Validate it on fresh comparisons, with the consumer and interface unchanged.

The paper states the template: "Exactly equal retained states give a readout difference beyond measurement error → reject sufficiency → separate the violating histories in a revised representation; test it with consumer and interface unchanged." Selection works only within the supplied library. It does not discover arbitrary missing variables or a unique encoding.

On TEDDY, E5-M performed only the first type of check: first order, one observer, no library, no selection and no fresh validation. Its registered verdict was inconclusive. The PR #16 review (section 7) drafted E7 to do all three steps on TEDDY. This page registers that draft, adapted to the architecture checks below.

## 2. Architecture facts the design relies on (checked)

Each fact was checked on 3 training cells (31464, 36167, 47423; CPU) and recorded in `architecture_check`. No val or site4 cell was used.

- **Layer 12 has no positional input.** TEDDY-G adds a learned position embedding once, at the input (`teddy_encoder.hidden_states`: `embeddings(ids) + position_embeddings(arange(L))`).
  - Every encoder layer is a plain post-norm `torch.nn.TransformerEncoderLayer`: d 512, 8 heads, d_hid 1024, GELU, eval mode. Its only parameters are the attention projections, two linear layers and two LayerNorms.
  - Layer 12 is therefore permutation-equivariant, and gene-mean pooling is invariant. The consumer is a symmetric function of the multiset of layer-11 token states.
  - Numerically: permuting the layer-11 rows changes the layer-12 states by at most 5.8e-06 (summation order) and the head by at most 6e-08. A swap patch changes the head by at most 9e-08.
  - **Consequences.** Swap patches are invisible, so E7 uses none. Gene identity reaches the consumer only through the states themselves. The full set of token states is sufficient by construction but has no matched pair, so it is uninformative.
- **Padding.** The official preprocessing pads to 2,048 tokens with a bool key-padding mask. teddy_mm passes the bool mask; the official `model.py` float mask would be *added* to the logits on MPS, so padding would be attended there.
  - E7 runs every cell unpadded on its own L real tokens, as E5 and E5-M do.
  - On the real tokens, the padded bool-masked module forward equals the unpadded forward to within 5e-06.
  - Changing a padding state changes nothing: 0.0 exactly.
- **Explicit layer.** E5's validated `layer_fn` (explicit attention) equals the module's layer 12 to within 4.8e-06.
  - The explicit float32 layer-11 gene-mean equals the official stored one (an fp16-autocast run) to within 1.9e-04–2.8e-04.
  - float32 and float64 layer-11 states differ by at most 5.3e-06.
- **Mean preservation.** Replacing two token states by their average changes the gene-mean by 1.5e-08 (rounding only). A mean-preserving patch of layer-12 states leaves the head unchanged to within 6e-08.
- **Scale.** A single average patch moves the head by 4e-05 to 2.5e-04 in p95 units, about 1e-04 of the measured NK–T gap. So intervention sizes must scale with the cell (section 4).

## 3. Consumer, observable, tolerance

- **Consumer (fixed, never trained).** Layer-11 token states → E5's `layer_fn` with layer 12's own weights → mean over the cell's real tokens → E5's `load_head`.
  - `load_head`: z / (‖z‖ + 1e-6) → MLP → NB decoder mean at the training-median size factor, divided by the training p95 of measured ADT.
  - The head is the frozen phase-1 head `outputs/cite_phase1_official/best.pt` (sha256 `b0c543e2…`). TEDDY weights sha256 `ba0a98f9…`.
- **Observable.** The head's predictions of CD56, CD94, CD335 and CD3: E3's gap proteins, the same as E5's and E5-M's.
  - The difference between two histories is the maximum over the 4 proteins of |ΔO_p| / u_p.
  - **Units.** u_p = |median over val primary-key NK cells − median over val primary-key T cells| of the measured ADT_p / p95_p. It uses val donor 18303 only: 254 NK and 2,990 T cells.
  - u = 1.5631 (CD56), 1.6140 (CD94), 1.5103 (CD335), 1.0892 (CD3). All are above the 0.05 floor, so no protein is dropped.
  - For comparison, the head's own predicted val NK–T gaps are 0.998, 1.182, 1.009 and 0.788. These are reported only as a re-expression.
- **Follow-up observable.** The NK–T score s = mean(CD56, CD94, CD335) − CD3 (E5's head-score lens), with unit u_s = 2.5786. It is used only by the follow-up.
- **Tolerance:** tol = **0.05** in these units, i.e. 5% of the measured val NK–T difference of a protein. Fixed before any patch output on a val or site4 cell.

## 4. Candidate library and history class

**Candidates**, in declared (dimension) order:

1. **gene-mean**: φ₁(H) = mean_t h_t (512).
2. **gene-mean + G**: the gene-mean plus the layer-11 states of the cell's tokens in G (512·(1+g)).
   - **G = {NCAM1, KLRD1, NCR1, CD3E, CD3D, CD3G}**: the coding genes of the 4 panel proteins (CD56, CD94, CD335, CD3). Ensembl ids are E5's `TOKEN_PROBE_IDS`; all are in TEDDY's vocabulary.
   - G is fixed from the observable alone, with no data used: these are the tokens most directly tied to the readout, as the sensor–recipient association was in ANM's study.
   - At most 6 tokens, so candidate 2 is a minimal augmentation. E5-M's other NK genes (FCGR3A, KLRF1) code proteins outside the panel and are not in G.
3. **all token states** (512·L): no two distinct histories share it, so it has no matched pair. Uninformative by construction.

**Histories** are exact mean-preserving state interventions at the layer-11 cut. A set of disjoint token pairs has each pair's two states replaced by their average. There are no swaps.

- **Magnitude.** f is the share of the cell's non-G tokens averaged: m = max(1, ⌊f·(n_nonG − g)/2⌋) pairs.
  - Development uses f ∈ {1/32, 1/8, 1/2}, all three. Confirmation uses the same three with fresh pairs, plus the fresh magnitude **1/4**, which is never used in development.
- **Per magnitude:** 8 matched patches (m non-G pairs). Each has a **G twin**: the same pairs, plus every G token averaged with its own non-G partner outside the patch. The twin has the same gene-mean; only the G states differ.
- **Per cell:** 8 **G-only** patches (every G token averaged with a random non-G partner), plus 2 **clamp** patches (positive control, section 7).
- **Matched pairs** (equal candidate state):
  - Candidate 1: (unpatched, matched_i), (unpatched, twin_i), (matched_i, twin_i) and (unpatched, G-only_j).
  - Candidate 2: (unpatched, matched_i).
  - Candidate 3: none.
- **G witness pairs**: (matched_i, twin_i) and (unpatched, G-only_j). They have equal gene-mean and differ only in the G states, so they are matched for candidate 1 and not for candidate 2.
  - They are G-specific at every magnitude, because both sides share the same context patch. This applies ANM's lesson that the development roster must be able to separate the candidates.
- **Seeds.** `default_rng([20260930, 7, phase tag, cell id, 1024·f])` for each magnitude and `[…, cell id, 1]` for the G-only patches. The phase tags are dev 1, followup 2, site4 3 and external 4, so confirmation pairs are fresh.
- **Eligibility.** At least 1 G gene among the cell's tokens, and at least 256 non-G tokens. In val, 249 of 254 NK and 2,986 of 2,990 T cells are eligible.
- **Per-cell statistics** cover every magnitude of the stage:
  - D_gene_mean = the maximum difference over candidate 1's matched pairs.
  - D_gene_mean+G = the maximum over candidate 2's matched pairs.
  - W = the maximum over the G witness pairs (W_context and W_pure are also reported).
  - P = the larger of the two clamp differences.

## 5. Development (val donor 18303 only)

- **Roster.** 100 NK + 100 T eligible primary-key cells, drawn with `default_rng([20260930, 7, 0])`. The cells and processing order are listed in the JSON (`rosters.dev`, index sha256 `421c0fb6…`).
- **Selection rule** (as in ANM's study).
  - A candidate is retained if it has at least 1 matched pair and its matched differences on the development roster are within tol. This is summarised as the 95th percentile, over cells, of the per-cell maximum over the whole class.
  - That summary is the point value of the statistic the confirmation bounds. The literal reading (every difference of every cell ≤ tol) is recorded beside it and never decides.
  - The minimal retained candidate is selected, with ties broken by declared order. If none is retained, E7 stops and the record is kept.
- **Gates.** The implementation checks (section 7) must pass on every development cell, and the float64 re-check on the first 10 cells of the order. Otherwise the record says *stopped* and no follow-up runs.
- **One declared follow-up.** It runs only after a null primary selection, and only once. It is labelled "follow-up after a null primary selection".
  - The observable is restricted to the NK–T score s. The library, history class and rule are the same.
  - It uses a fresh roster fixed here: 100 NK + 100 T val cells disjoint from the development roster (`rosters.followup`, sha256 `5f8b0f8d…`).
  - Its develop stage is refused until the null primary record is committed.
- **Record.** `registration/addenda/E7_selection.json`, with a line in HASHES.txt. It must be committed before any site4 or external forward pass. It pins a digest of every development output.
- **Smoke pool.** The val cells in neither roster (49 NK, 2,786 T; sha256 `fa656a2e…`). At most 5 of them are used for smoke runs, so a smoke run never touches a development cell.

## 6. Confirmation

- **Site4 (registered).** Test_primary donors 13272 and 19593, eligible primary-key NK and T cells.
  - Every cell of E5's `e5_cells.npz` is excluded: all of E3's 278 NK–T pair cells (checked against E5's per-donor pair counts) and E5's 200 random cells. So are E5-M's cells in `e5m_cells.npz`. Both files are pinned by sha256.
  - Up to 100 cells per donor and class are drawn with `default_rng([20260930, 7, 3])`.
  - Magnitudes 1/32, 1/8, 1/2 (fresh pairs) and 1/4 (fresh magnitude).
  - The observable is the one the selection was made on (panel, or the score if the follow-up decided); the other is reported as secondary.
- **External (secondary).** E6's 8 donors (Hao 2021, GSE164378; E6's d5k pack, sha256 pinned through E6's addendum).
  - Annotation-only key with E6's committed map, because E6's key was not validated. Up to 25 NK + 25 T eligible cells per donor, drawn with `default_rng([20260930, 7, 4])`.
  - Labelled **"RNA possibly seen by TEDDY in pretraining"** and **"key not validated"**.
  - Pooled bounds over the 8 donors. Per-donor values are reported. It never changes the registered verdict.

## 7. Statistics, controls, verdict

- **Bootstrap.** As in the registration's common.statistics: two-stage (donors, then cells within donor and class), B = 2000, seed 1, percentile 95% intervals. Per donor: a cell bootstrap within the donor and class.
- **Support for candidate k.** The upper end of the 95% interval of the 95th percentile of D_k is ≤ tol, pooled and in each primary donor.
- **Rejection of candidate k.** The lower end of the 95% interval of the median of D_k is > tol, pooled and in each primary donor. Otherwise unresolved.
- **G witness detected.** The lower end of the 95% interval of the median of W is > tol.
- **Positive control (registered).** Two clamp patches move every layer-11 token by ±1 × dmu₁₁.
  - dmu₁₁ is E5's training NK − T difference of the layer-11 gene-mean (`e5_directions.npz` dmu[11], training cells only, ‖dmu₁₁‖ = 0.332; copied into the JSON).
  - The clamp moves the gene-mean by the full NK–T difference and leaves the deviations unchanged. This is an NK/T-sized change of the retained state, as ANM's queried-sensor +1 was.
  - The control is detected if the lower end of the 95% interval of the median of P is > tol. If it fails, any support claim is **UNINFORMATIVE**.
- **Detectability thresholds.** For the positive control and the G witness, the threshold is tol. Numerical thresholds are the gates below.
- **Secondary control (descriptive).** The review's input push: each G gene's token embedding doubled at the input, full forward, on the first 50 cells of each confirmation order. It acts at layer 0, so it moves the cut's mean and deviations together; it does not decide.
- **Implementation checks** (every cell; any failure gives NOT_VALIDATED):
  - padding-position patch ≤ 1e-6 gap units;
  - layer-12 mean-preserving patch ≤ 1e-5 gap units;
  - permutation of the layer-11 rows ≤ 1e-5 gap units;
  - gene-mean change of every pair patch ≤ 1e-5;
  - explicit vs module layer 12 (unpadded and padded-masked) ≤ 1e-4;
  - explicit vs official layer-11 gene-mean ≤ 1e-3.
- **float64 re-check.** The first 10 cells of every roster's order are recomputed in float64 on the CPU, with the same histories. If any difference moves by more than tol/10 = 0.005, the stage is UNRESOLVED_PRECISION.
- **Verdicts** (`v3_e7.verdict`):

| verdict | when |
|---|---|
| LOOP_COMPLETE | gene-mean + G selected, supported, gene-mean rejected and the G witness detected: reject → select (library) → validate |
| SUPPORTED_REVISION_NOT_SHOWN_NECESSARY | gene-mean + G supported, but the gene-mean not rejected or the witness not detected |
| GENE_MEAN_SUPPORTED | the gene-mean selected and supported (positive control passing) |
| UNINFORMATIVE | the selected candidate's support cannot be declared: positive control failed |
| SELECTED_REJECTED / SELECTED_UNRESOLVED | the development selection does not hold up / is undecided on fresh comparisons |
| REJECTED_NO_REVISION | nothing selected; both the gene-mean and gene-mean + G rejected on fresh comparisons |
| REJECTED_GENE_MEAN_NO_SELECTION / NO_SELECTION_UNRESOLVED | nothing selected; gene-mean rejected / not rejected |
| NOT_VALIDATED / UNRESOLVED_PRECISION | an implementation check failed / the float64 rule failed |

## 8. Leakage check

The registration core contains the units, rosters, eligibility and G. It was rebuilt twice:

- with every site4 row poisoned: ADT → U(0, 8), size factors → U(0.1, 10), cell types random, RNA values → U(0.5, 50) on the same sparsity, official z → N(0, 1). The core is **identical** (`aaaf9ab1…`).
- with the val rows' ADT and z poisoned. The core **differs** (`f47bb0a5…`).

The builder reads site4 only as the sha256 of E5's and E5-M's cell files (bytes hashed) and their published counts. dmu₁₁ is a training-only file constant added after the check.

## 9. Disclosure

- **Seen before registration.** The registering session had read all published v3 results (E1–E6, E5-M), including E5-M's per-layer clamp residuals on site4: the head's 0.709 (pair cells) and 0.72 (random cells) at layer 11. It had also read the PR #16 review and the val/train primary-key NK/T cell counts.
- **Computed before registration.** No patch output was computed on any val or site4 cell. Patch outputs were computed on training cells only:
  - **Architecture check:** 3 cells (section 2).
  - **First code run:** 2 training cells without G tokens, to exercise the code paths only.
  - **First code test:** 2 annotated training cells (1 NK, 1 T), with an earlier draft. That draft used one-sided twins as the positive control (only one token of each pair moved) and calibrated the magnitude on them.
    - Median per-cell maxima at f = 1/32 / 1/8 / 1/2: gene-mean 0.014 / 0.065 / 0.166; gene-mean + G 0.013 / 0.055 / 0.149.
    - G witness 0.007 / 0.011 / 0.027; one-sided twins 0.002 / 0.007 / 0.026.
    - A twin-calibrated magnitude (target 2·tol) could never have been reached, because the twin moves the gene-mean only by a random-walk sum. The draft was revised to the clamp control and the fixed three-magnitude class.
  - **Second code test:** 3 training cells (2 NK, 1 T) with the clamp at κ = 0.2: median P 0.049 (at the tolerance), median W 0.019, q95 of the per-cell maximum 0.227 (gene-mean) and 0.202 (gene-mean + G).
    - κ was then set to 1, the full training NK–T difference and the interpretable NK/T-sized change. κ = 0.2 had no rationale beyond being small.
    - The same 3 cells were rerun on MPS to check the device path: float32 and float64 agree to within 1e-6.
  - Nothing else changed after these code tests. They are recorded in the JSON under `pre_registration_code_test`.
- **What the code tests suggest.** On training cells, both the gene-mean and gene-mean + G fail the 0.05 tolerance at f = 1/2, and the G witness is small. If development behaves the same way, the registered path is: null primary selection, the declared follow-up, and then confirmation of the rejections only (REJECTED_NO_REVISION or similar). That would be the analog of the paper's broader screen that "selected no candidate". It is stated here so that a null result is not read as a surprise or as a failure of the procedure.

## 10. Order of operations (enforced by the script)

1. `register` writes E7.json and its HASHES line. It is refused once any development output or selection record exists. **Commit E7.json + HASHES.txt.**
2. `develop --roster dev` is refused unless the addendum is committed, its declared block equals the script's, its leakage check passed, and the reused E5 functions, `v3_e7.py` and the head are the recorded versions.
3. `select --roster dev` writes E7_selection.json and its HASHES line. **Commit them.** If `follow_up_required` is true: run `develop --roster followup` (refused until the null record is committed), then `select --roster followup`, then commit again.
4. `confirm --family site4`, then `confirm --family external`. Both are refused unless the addendum and a final, committed, non-smoke selection record that names this addendum are in place.
5. `report` writes E7_results.json and REPORT.md (refused under the same conditions).

## 11. Run commands and expected compute

From the repo root (`R` = the main checkout):

```
PY=<venv with torch 2.14>/bin/python; OUT=$R/outputs/v3/E7
PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage develop --roster dev --out-dir $OUT --device auto
$PY scripts/mode_a_inverse_e7.py --stage select --roster dev --out-dir $OUT
git add registration/addenda/E7_selection.json registration/addenda/HASHES.txt && git commit -m "E7 development selection record"
# only if the record says follow_up_required:
PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage develop --roster followup --out-dir $OUT --device auto
$PY scripts/mode_a_inverse_e7.py --stage select --roster followup --out-dir $OUT   # then commit the record again
PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage confirm --family site4 --out-dir $OUT --device auto
PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage confirm --family external --out-dir $OUT --device auto
$PY scripts/mode_a_inverse_e7.py --stage report --out-dir $OUT
```

Every device stage is resumable: rerunning the same command skips finished cells. `--max-minutes N` pauses with exit code 75. Progress goes to `$OUT/progress.log` and `*/e7_progress.json`.

**Smoke after registration** (commit `6dc5999`; 5 val smoke-pool cells outside both rosters; CPU float32 with every cell also re-checked in float64):

- It ran every stage: develop, select, the declared follow-up, select, confirm and report.
- Timings: 4.3 s per cell for develop (58 histories) and 6.1 s per cell for confirm (74 histories plus the input push), at 718–1,601 tokens.
- Every implementation gate passed: padding patch 0.0; layer-12 mean-preserving patch 7e-09; permutation 1.1e-07; float64 within 1.3e-06.
- It is a code check only, never a result, and it decides nothing.

Expected compute, from the code tests and the smoke (MPS float32 uncontended is faster than the CPU figures):

| run | cells | histories per cell | expected time |
|---|---|---|---|
| development | 200 | 58 | CPU about 3 s per cell (about 4 s with float64 on the first 10); about 10–15 min |
| follow-up (if needed) | 200 | 58 | the same |
| site4 | ≤ 400 | 74 | CPU about 4–6 s per cell (push on the first 50); about 30–40 min |
| external | ≤ 400 | 74 | longer cells (median about 1,900 tokens, attention cost about 2–3×); about 60–80 min on CPU |

Total about 2–2.5 h on CPU and about 1 h on an uncontended MPS GPU. Select and report are CPU and take minutes.

## 12. What E7 can and cannot show

**Can show.** For one fixed consumer (TEDDY layer 12 + pooling + our frozen head), one observable panel and one declared class of mean-preserving state interventions at the layer-11 cut, E7 can show two things:

- whether the layer-11 gene-mean is rejected as a sufficient retained state on fresh cells and donors;
- whether a declared marker-token augmentation, selected on development data, gets bounded support there.

If both hold, it would be the first complete reject → select → validate step on TEDDY.

**Cannot show:**

- that RNA lacks NK–T information;
- that the head improves (the consumer is unchanged; E3 rejected repair on frozen TEDDY);
- that any representation is unique or minimal beyond the declared library;
- transfer to other consumers, layers, observables or intervention classes;
- behaviour under real-cell (on-manifold) perturbations: these are state interventions, not input histories. ANM v2: "Propagation imposes the state displacements; the readout alone is measured afresh";
- clinical value.

The external family carries its two labels and never decides.

> **中文摘要：** E7 在 TEDDY 上做论文 v2 逆向闭环的类比，消费者固定（第 12 层 + gene-mean 池化 + L2 归一化 + 冻结的 phase-1 head）。观测量是 head 对 CD56/CD94/CD335/CD3 的预测，单位是 val 供体上实测的 NK–T 差距，容差为 0.05。
>
> - **候选库（事先声明）：** 第 11 层 gene-mean；gene-mean 加 G（面板蛋白的 6 个编码基因）的 token 状态；全部 token 状态（无匹配对，无信息）。
> - **历史：** 在第 11 层切点做精确的保均值配对平均，按细胞非 G token 的比例取 1/32、1/8、1/2，确认时另加 1/4。不用交换补丁，因为第 12 层没有位置输入，已核实。
> - **见证对：** 同一补丁加不加 G 配对。
> - **开发：** 只用 val 供体 18303，按预设规则选择；选不出就停止，并只做一次事先声明的后续（读出改为 NK–T 分数）。
> - **确认：** 用新的 site4 细胞（排除 E3/E5/E5-M 的细胞）、新的补丁对和新幅度；外部供体为次要家族，带"可能被 TEDDY 预训练见过"标签。
> - **对照：** 阳性对照是沿训练 NK–T 方向的 gene-mean 钳位移动；不变性对照（padding、第 12 层保均值、置换）用来检查实现；另有 float64 复核。
> - **透明说明：** 在训练细胞上的代码测试提示，两个候选都可能在 f = 1/2 时超出容差，因此"无修订、只确认否决"是一个预期中的可能结果。
