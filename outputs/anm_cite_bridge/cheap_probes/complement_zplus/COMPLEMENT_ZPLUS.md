# Complement z⁺ (Mode B)

**Not Mode A / not fusion audit.**

- held lineage: `myeloid` proteins ['CD16', 'CD11c', 'CD36']
- z+ proteins: ['CD19', 'CD72', 'CD22', 'CD3', 'CD2', 'CD5']

## site4 variants (must among near @ min_cos)

| variant | near | must | frac must | cells w/ must |
|---|---:|---:|---:|---:|
| z512_held_eval | 177701 | 83469 | 0.4697159835904131 | 11926 |
| zplus_adt_true_complement | 308309 | 184313 | 0.5978190711266943 | 16675 |
| zplus_adt_pred_complement | 315214 | 188561 | 0.5981999530477707 | 16739 |
| zplus_adt_pred_all9_leakcheck | 318336 | 187041 | 0.5875584288299156 | 16736 |

Primary relative must-reduction (zplus_adt_true_complement): **-0.27272456550677143**

## OOD variants

- `zplus_adt_true_complement`: frac_must=0.31351478811561007 near=29686
- `zplus_adt_pred_complement`: frac_must=0.4309013099525077 near=38322

Phase2 note: {'phase2_joint_embedding': 'not available as exported z; cite_phase2 has best.pt + 4k raw ADT(134) but no aligned site4 joint z sidecar', 'used_instead': 'ADT true/pred complement features on site4 panel (9) with held-lineage protocol'}

Do not claim Pearson > ~0.61. Not clinical.
