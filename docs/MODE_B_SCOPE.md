# Mode B scope (claim boundary)

**Mode B (implemented in this repo):** TEDDY + our head as a **typed evidence source** → ANM **finite_field** decisions on CITE holdout (site4/test). Editable observer YAML (O0/O1/O2); LOO/flip attribution on **decisions**. With all markers entering at once, ANM's calls equal the fixed rule's for the same question; the claim is the declared, auditable layer.

**Scope in one sentence.** teddy_mm runs ANM's open-loop decision layer: TEDDY is a prescribed evidence source and does not read the field, so ANM's dynamics (state feedback, operator memory, delays) are not tested here.

This is **not** Mode A. Do not read current results as residual-stream / Jacobian / in-silico gene-perturb claims.

> **Corrections (29 Sep 2026): rerun done for Experiments 1, 3, 4 and 5.** See the notice at the top of the [README](../README.md#corrections-29-sep-2026-rerun-done-for-experiments-1-3-4-and-5) and the [site](https://danielchen26.github.io/teddy_mm/#corrections). Reports: [`HARD_PROOF_V2`](reports/HARD_PROOF_V2.md), [`SCOPE_REFINE_PROOF_V2`](reports/SCOPE_REFINE_PROOF_V2.md).
>
> - All markers now enter the engine at once (the first run entered them one per step in panel order, a hidden weight), and predictions use one train-median size factor, 0.92663, instead of a per-cell factor computed from measured protein. The answer key, the fixed rule and ANM share one scoring function.
> - With simultaneous markers ANM's lineage scores are a fixed multiple of the fixed rule's, so ANM makes the same call as TEDDY + fixed rule (re-coded per question) on every site4 cell for O0/O1/O2, by construction: "no calls" 272 / 1,493 / 2,551, accuracy of calls made 0.9429 / 0.9667 / 0.9545, exactness 0.9336 / 0.9312 / 0.8559. ANM's value in Mode B is a declared, auditable layer (written question, P_f/Q_f checks, LOO reasons), not better calls.
> - The O0 and O2 rule shifts (101 → 272 and 1,796 → 2,551 "no calls") come from the export's size factor and the official TEDDY-G preprocessing, not the scoring. O1 is redefined (weighted mean, key weight 2): 881 → 1,277 from the scoring, 1,277 → 1,493 from the size factor and the preprocessing. The first run's ANM accuracy edge (317 / 834 / 3,373 "no calls"; 0.9496 / 0.9723 / 0.9603) came only from declining more.
> - Experiment 3: the most frequent top marker is now CD2 (17.9%), CD5 close behind (23.6% → 17.3%); under simultaneous timing the top marker is the protein with the largest normalised value in the called lineage, so it is not ANM-specific. Flip rate 0.2999; within-lineage permutation p = 0.0099 for CD2's share (15.7% shuffled; the floor for 100 shuffles).
> - Experiment 4: soft_P never gates better than the fixed rule's margin at matched coverage (every difference ≤ 0; O2 at 40%: 0.9988 vs 0.9988); soft_P reaches 1.11, a score, not a probability. Hard cells are mostly out of scope (NK/ILC 3.1×, erythroid 3.4× enriched).
> - Experiment 5, on a disjoint split (11,725 / 5,025): labels to match ANM are logistic 500, MLP 1,000 (was 200), grid 50.
> - The answer key still has no NK or out-of-scope class; out-of-scope cells still get a call 73–95% of the time. Experiments 2 and 6 stay withdrawn pending redesign (prototypes: with degraded RNA, ANM with a trust per level = the fixed rule with a stricter bar per level; the TEDDY + head predictions do separate NK from T in the right direction but compress the differences to about a quarter to a third of the measured ones (CD335 0.88): 1,053 of 1,506 site4 NK–T pairs involve "gdT CD158b+" cells, which look NK-like in both predicted and measured protein, and on the other 453 pairs the median predicted within-pair gap for CD56 / CD94 / CD335 / CD3 is 0.148 / 0.100 / 0.364 / 0.205 vs measured 0.621 / 0.321 / 0.413 / 0.756; first Mode A probe: linear NK-vs-T probes separate NK from T at every TEDDY layer (held-out AUC 0.991–0.997), yet linear protein probes from every layer keep only about a fifth to two-fifths of these pairs' measured differences (final embedding: CD56 0.28, CD94 0.40, CD335 0.40, CD3 0.28) and our head is no worse, so, as far as linear probes can tell, the gap sits in TEDDY's gene-mean representation rather than in our head; linear probes only, pairs chosen for closeness in the embedding). "0 labels for a new question" holds for any written rule, including the fixed rule.
> - Phase 1 re-evaluated with the train-median factor and the official TEDDY-G preprocessing: mean test Pearson 0.603 (MLP), 0.582 (flow matching); was 0.610 / 0.596. The earlier non-official preprocessing gave the same MLP Pearson (0.603).
>
> Disclosures: donor 15078 is in training at other sites and is 32.6% of the site4 test cells; the "out-of-site check" is held-out donor 18303 at a training site, which was also the validation set; the 9 panel proteins were chosen by test-set Pearson; the numbers above use the official TEDDY-G preprocessing (counts / total × 10⁴, divided by TEDDY gene medians, top 2,048 rank tokens, no CLS, mean over gene tokens); the earlier `z_512` (table below) deviated from it.

---

## z_512 loading factor (required definition)

| Item | Value |
|---|---|
| Symbol | `z_512` |
| Source | Frozen TEDDY-G last-layer hidden states |
| Pool | **Mean-pool over real tokens** (attention mask), not token 0, not a disease token |
| Context length | **512** (CITE embed; see `scripts/03_embed_rna.py --seq-len 512`; the stored embeddings reproduce at 512 with cosine 1.000000, vs median 0.970 at 1024) |
| Not | TEDDY pretrain context **2048** |
| Not | A dedicated disease / CLS token readout |
| Normalisation | **No TEDDY gene-median normalisation** (the official pipeline divides each gene by its TEDDY median before ranking) |
| Dim | `d_model = 512` after mean-pool |
| Deviation from official TEDDY-G | The official model returns token 0 by default (no CLS token is added, so token 0 is the top-ranked gene); the official tutorial mean-pools over all 2,048 positions, padding included, after gene-median normalisation. With the medians applied (512 tokens, same masked mean) the cosine to the stored `z_512` is **0.649** (median over 48 site4 cells; fully official 2,048-token setups: 0.205 token 0, 0.419 tutorial mean). The results above are now rerun with the official preprocessing, which gave the same phase-1 Pearson; this table documents the earlier embedding. |

Bridge sidecar `outputs/anm_cite_bridge/z_rna_512.npy` stores **L2-normalized full** `z_512` (export row order). Compact `z_rna_export.npy` keeps only `z_keep=32` leading dims for typed-event JSONL — **not** the reverse-step1 probe space.

---

## Implemented (Mode B)

- RNA→ADT predictions from our head on frozen TEDDY's embedding (`best.pt`) as typed δu evidence into ANM finite_field
- Editable observer YAML: `bridge_anm/readouts/cite_lineage_O{0,1,2}.yaml`
- Missing-modality abstain demos (honest silence under weak channels; Experiment 2, withdrawn pending redesign)
- Closed-form LOO / top-1 flip attribution on **declared decisions**
- Scope gate (never better than the fixed rule's own margin) / zero-label observer transfer vs TEDDY + trained classifier (see proof reports)
- Mode B reverse **decision** loop on must-separate pairs (veto→complement→verify; Experiment 6, withdrawn pending redesign)

Verified numbers only — cite [`HARD_PROOF_V2`](reports/HARD_PROOF_V2.md) and [`SCOPE_REFINE_PROOF_V2`](reports/SCOPE_REFINE_PROOF_V2.md) (corrected rerun); first run: [`HARD_PROOF`](reports/HARD_PROOF.md), [`SCOPE_REFINE_PROOF`](reports/SCOPE_REFINE_PROOF.md), [`MISSING_MODALITY_ANM_DEMO`](reports/MISSING_MODALITY_ANM_DEMO.md) (withdrawn).

## NOT claimed (Mode A and adjacent)

- Mode A: residual stream as field, layer Jacobian, in-silico gene perturbs (the only Mode A result so far is the linear layer probe on NK–T look-alike pairs, `outputs/mode_a_official/`)
- ANM's dynamics: TEDDY does not read the field, so state feedback, operator memory and delays are not exercised
- Gated ±ε / ±ε/2 G1–G4 protocols
- Mode-A reverse closed loop (representation response to perturbs)
- Fusion audit
- ATAC / chromatin foundation

## What current conclusions *are* about

Current proofs speak to **decision-layer sensitivity to the TEDDY + head evidence** (edit observer → abstain / P_f / Q_f / LOO flips), **not** to TEDDY representation response under gene/protein perturbs.

They **cannot** yet answer whether `z_512` (or the export sidecar) is **sufficient for protein readout** in the strong sense. A reverse step-1 “must-separate” probe (near-identical z, different true lineage/protein labels) is a **sufficiency-for-readout** diagnostic — still Mode B / evidence geometry, **not** Mode A perturb-response.

## Honest weak-model note

**ADT-only** is a weak channel in the missing-modality demo (panel Pearson drops; the fixed rule's abstain can stay 0 while Q is mediocre). The abstain story is intentionally honest: ANM silence under weak evidence beats silent over-answer. Not clinical. Experiment 2 is withdrawn pending redesign. Do not claim Pearson > phase-1 full-ADT **~0.60** (0.603 with the train-median size factor).

## Cheap next probes (Mode B only)

Four cheap probes ranked in [`docs/reports/CHEAP_PROBES_RANKING.md`](reports/CHEAP_PROBES_RANKING.md). Still **not** Mode A. Gated ±ε G1–G4 here are **decision-layer** nudges on typed evidence, not residual-stream / gene-perturb protocols.


## Reverse decision loop (Mode B — implemented; Experiment 6 headline withdrawn pending redesign)

Script: [`bridge_anm/reverse_loop_small.py`](../bridge_anm/reverse_loop_small.py)

**Veto → complement → verify** on reverse_step1 `z_512` must-pairs. Winner = true-ADT typed evidence + O1 (not `z+` alone). Deep writeup + NotebookLM posters: [`docs/reports/ANM_HELPS_TEDDY_LOOP.md`](reports/ANM_HELPS_TEDDY_LOOP.md). Still **not** Mode A perturb-response.

## Reverse step 1 (sufficiency probe only)

Script: [`bridge_anm/reverse_step1_must_separate.py`](../bridge_anm/reverse_step1_must_separate.py)

Reads full `z_rna_512.npy` (not compact `z_keep=32`) + true ADT panel, finds near-identical z pairs with different lineage/protein labels (**must-separate**). Report: [`docs/reports/REVERSE_STEP1.md`](reports/REVERSE_STEP1.md).

Still Mode B / evidence geometry — **not** Mode A perturb-response.
