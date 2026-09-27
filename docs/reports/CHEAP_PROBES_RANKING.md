# Cheap next probes — comparative ranking (Mode B)

**Not Mode A.** Evidence-only ranking of four cheap probes on `teddy_mm` (site4/test z_512).

**Strongest signal:** `1_must_separate_autopsy`

## Ranking (by |primary metric|)

| rank | probe | primary metric | value | note |
|---:|---|---|---:|---|
| 1 | `1_must_separate_autopsy` | `max(adt_lineage_disagree_frac, same_celltype_must_frac)` | **0.5810** | same_celltype_must_frac=0.581 (within-type insufficiency) vs adt_lineage_disagree_frac=0.412 |
| 2 | `3_complement_zplus` | `relative_must_reduction_vs_z512` | **-0.2727** | POSITIVE would mean z+ separates held proteins better; observed NEGATIVE => z+ worsens must-rate |
| 3 | `2_gated_loo` | `any_gated_flip_frac` | **0.0273** | higher => more O0 decisions flip under ±ε/±ε/2 on top LOO protein |
| 4 | `4_loading_scan` | `must_frac_spread_across_dims` | **0.0154** | higher spread => must-rate depends on dim-loading proxy (ctx re-encode blocked: no torch) |

## 1) MUST-SEPARATE autopsy

- n_pairs written: **5000**
- ADT-lineage disagree frac: **0.4118** (n=2059)
- key-marker disagree frac: **0.6574**
- coarse disagree frac: **0.1580**
- **same fine cell_type must frac: 0.5810** (within-type insufficiency / annotation-fine drift)
- buckets: `{'bio_adt_lineage': 2059, 'bio_key_marker': 2467, 'bio_coarse_cell_type': 474}`
- top lineage confusion: `{'myeloid|t_lineage': 1490, 'b_lineage|t_lineage': 375, 'abstain|t_lineage': 213, 'b_lineage|myeloid': 194, 'abstain|myeloid': 168, 'abstain|b_lineage': 24}`
- top proteins by mean |Δ|: CD5=0.702, CD3=0.701, CD16=0.694, CD2=0.665, CD36=0.498
- Note: hierarchical `bio_minus_noise` saturates at 1.0 on this written must-set; refined metrics above.

## 2) Gated LOO ±ε / ±ε/2 + G1–G4

- eps=**0.1**, eps/2=**0.05**, kept=**1500** / attr 16433
- existing LOO top-1 flip rate (full): **0.2773**
- subsample LOO flip: **0.2873**
- any gated flip (1−G1): **0.0273**
- gate_counts: **G1=1459, G2=0, G3=0, G4=41**
- All gated flips are **G4** (sign-asymmetric); none at symmetric ±ε/2 scale (G2/G3=0).
- Reading: O0 decisions are **robust to ±0.1 evidence nudges** on the top LOO protein, despite 27.7% LOO-ablate flips — sensitivity is ablation-scale, not ε-scale.

## 3) Complement z⁺ (held lineage = myeloid)

- primary relative must-reduction vs z_512: **-0.2727** (negative ⇒ z+ **worsens**)
- `z512_held_eval`: frac_must=**0.4697** near=177701 must=83469
- `zplus_adt_true_complement`: frac_must=**0.5978** near=308309 must=184313
- `zplus_adt_pred_complement`: frac_must=**0.5982** near=315214 must=188561
- `zplus_adt_pred_all9_leakcheck`: frac_must=**0.5876** near=318336 must=187041
- OOD: `{'zplus_adt_true_complement': 0.31351478811561007, 'zplus_adt_pred_complement': 0.4309013099525077}`
- phase2 joint: {'phase2_joint_embedding': 'not available as exported z; cite_phase2 has best.pt + 4k raw ADT(134) but no aligned site4 joint z sidecar', 'used_instead': 'ADT true/pred complement features on site4 panel (9) with held-lineage protocol'}

## 4) Loading scan

- ctx 256/512/1024 re-encode **blocked**: `["torch_unavailable: ModuleNotFoundError: No module named 'torch'"]`
- ckpt weights present (['ckpt_weights_bytes=284887424', 'processed_cite_arrays_present']); only missing torch in this env.
- available on disk: ctx=**1024** only (`data/processed/cite/z_rna.npy` 90261×512).
- proxy leading-dim subsets must_frac: `{'leading_dim_128': 0.4294535960286514, 'leading_dim_256': 0.44299313418608544, 'leading_dim_512': 0.444839617332583}`
- must_frac_spread: **0.0154** | rel(128 vs 512): **-0.0346**

## Claim boundary

Mode B only (typed evidence / decision sensitivity / evidence geometry). Not Mode A residual-as-field, Jacobian, gene perturbs, fusion audit, or clinical. Do not claim Pearson > ~0.61.
