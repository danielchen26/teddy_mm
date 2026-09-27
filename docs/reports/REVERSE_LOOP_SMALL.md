# Reverse loop (small) — must-separate pairs (Mode B)
**Not Mode A.** Veto → complement → verify on reverse_step1 `z_512` must-pairs.
No TEDDY retrain. Numbers from files / this run only.

## 1) VETO — z_512 insufficient

z_512 (mean-pool last-layer @ ctx 1024, L2 512-D) is INSUFFICIENT to separate these must-pairs: near-identical z carries disagreeing true ADT lineage / protein labels. Downstream readout must abstain, rewrite observer, or admit a complement evidence channel.

- pairs written: **5000**
- myeloid↔T pairs: **1490** (0.2980)
- lineage confusion: `{'myeloid|t_lineage': 1490, 'b_lineage|t_lineage': 375, 'abstain|t_lineage': 213, 'b_lineage|myeloid': 194, 'abstain|myeloid': 168, 'abstain|b_lineage': 24}`
- flag counts: `{'diff_adt_lineage': 2059, 'diff_protein_profile': 4639, 'diff_key_marker_lineage': 3287, 'diff_coarse_cell_type': 790}`
- myeloid↔T focus protein |Δ|: `{'CD5': {'n': 1490, 'mean_abs_delta': 0.7155161401769459, 'median_abs_delta': 0.40774770081043243, 'q75_abs_delta': 1.0865748226642609}, 'CD3': {'n': 1490, 'mean_abs_delta': 0.6981261888506428, 'median_abs_delta': 0.3616558760404587, 'q75_abs_delta': 0.8664438873529434}, 'CD16': {'n': 1490, 'mean_abs_delta': 1.0650999203649343, 'median_abs_delta': 0.907050609588623, 'q75_abs_delta': 1.776064932346344}}`
- myeloid↔T cos_sim quantiles: `{'min': 0.9800012707710266, 'q25': 0.9810009449720383, 'q50': 0.9822191596031189, 'q75': 0.983788788318634, 'max': 0.990542471408844}`

## 2) COMPLEMENT — cheap arms (no retrain)

| arm | description |
|---|---|
| (a) `O1_abstain` | observer rewrite: threshold 0.12→0.28, key_marker_boost×2 |
| (a) `O2_observer` | key-marker priority observer (rewrites expected) |
| (b) `true_ADT_typed` | admit holdout true ADT panel as typed ANM source events (Mode B complement channel; not production no-leakage path) |
| (c) `zplus` | ADT true/pred complement features as z⁺; prior global relative_must_reduction **worsened** must-rate |

### (c) z⁺ on same myeloid↔T pairs

- prior: Prior global z⁺ probe (held myeloid) relative_must_reduction primary=-0.27272456550677143; NEGATIVE => z⁺ alone worsens must-rate among near.
- cos_z512: `{'mean': 0.9826505275780723, 'median': 0.9822191596031189, 'frac_still_ge_0_98': 1.0, 'frac_dropped_below_0_95': 0.0}`
- cos_zplus_true: `{'mean': 0.8255820943005133, 'median': 0.8865448236465454, 'frac_still_ge_0_98': 0.09395973154362416, 'frac_dropped_below_0_95': 0.7208053691275168}`
- cos_zplus_pred: `{'mean': 0.9877049229289061, 'median': 0.9969751238822937, 'frac_still_ge_0_98': 0.8348993288590604, 'frac_dropped_below_0_95': 0.07583892617449664}`
- On FIXED myeloid↔T must-pairs: true-ADT z⁺ lowers cosine (frac still ≥0.98 = 0.09395973154362416; frac <0.95 = 0.7208053691275168), so geometry can separate this fixed set. Pred z⁺ does NOT (frac still ≥0.98 = 0.8348993288590604). Prior GLOBAL probe still found z⁺ worsens must-rate among its own near neighborhoods (relative_must_reduction primary=-0.27272456550677143). Decision winner remains true-ADT typed evidence + O1, not z⁺ alone.

## 3) VERIFY — same pairs (myeloid↔T focus)

| complement | soft_sep | correct_sep | false_agree (both decided) | abstain_pair_frac | Q_among_decided | n_pairs |
|---|---:|---:|---:|---:|---:|---:|
| `baseline_pred_O0` | 0.1738255033557047 | 0.08389261744966443 | 0.8536754507628294 | 0.032214765100671144 | 0.7409695817490495 | 1490 |
| `complement_a_O1_abstain` | 0.22416107382550335 | 0.032214765100671144 | 0.8310567936736161 | 0.06644295302013423 | 0.7777777777777778 | 1490 |
| `complement_a_O2_observer` | 0.3724832214765101 | 0.08389261744966443 | 0.8431018935978359 | 0.2557046979865772 | 0.8824110671936759 | 1490 |
| `complement_b_true_ADT_O0` | 0.5476510067114094 | 0.36845637583892615 | 0.4569491525423729 | 0.010067114093959731 | 0.9716981132075472 | 1490 |
| `complement_b_true_ADT_O1` | 0.625503355704698 | 0.17248322147651007 | 0.39658848614072495 | 0.05570469798657718 | 0.9868421052631579 | 1490 |
| `complement_c_zplus_geometry` | 0.7208053691275168 | None | 0.09395973154362416 | None | None | 1490 |

### All must-pairs (written set)

| complement | soft_sep | false_agree | Q_among_decided | n_pairs |
|---|---:|---:|---:|---:|
| `baseline_pred_O0` | 0.167 | 0.8926275182168881 | 0.8081765087605451 | 5000 |
| `complement_a_O1_abstain` | 0.1938 | 0.9034065441506051 | 0.8665363426337302 | 5000 |
| `complement_a_O2_observer` | 0.5946 | 0.7850503485670023 | 0.8514386376981797 | 5000 |
| `complement_b_true_ADT_O0` | 0.4074 | 0.6040774719673803 | 0.9567061759201497 | 5000 |
| `complement_b_true_ADT_O1` | 0.4694 | 0.5857805255023184 | 0.9899108794350092 | 5000 |

## Winner

**complement_b_true_ADT_O1**

- row: `{'complement': 'complement_b_true_ADT_O1', 'score': 2.8412603305318287, 'soft_separation_rate': 0.625503355704698, 'Q_among_decided': 0.9868421052631579, 'false_agreement_rate_among_both_decided': 0.39658848614072495, 'correct_separation_rate': 0.17248322147651007, 'abstain_pair_frac': 0.05570469798657718}`
- ranking: `[{'complement': 'complement_b_true_ADT_O1', 'score': 2.8412603305318287, 'soft_separation_rate': 0.625503355704698, 'Q_among_decided': 0.9868421052631579, 'false_agreement_rate_among_both_decided': 0.39658848614072495, 'correct_separation_rate': 0.17248322147651007, 'abstain_pair_frac': 0.05570469798657718}, {'complement': 'complement_b_true_ADT_O0', 'score': 2.610050974087993, 'soft_separation_rate': 0.5476510067114094, 'Q_among_decided': 0.9716981132075472, 'false_agreement_rate_among_both_decided': 0.4569491525423729, 'correct_separation_rate': 0.36845637583892615, 'abstain_pair_frac': 0.010067114093959731}, {'complement': 'complement_c_zplus_geometry', 'score': 2.3476510067114096, 'soft_separation_rate': 0.7208053691275168, 'Q_among_decided': None, 'false_agreement_rate_among_both_decided': 0.09395973154362416, 'correct_separation_rate': None, 'abstain_pair_frac': None}, {'complement': 'complement_a_O2_observer', 'score': 1.78427561654886, 'soft_separation_rate': 0.3724832214765101, 'Q_among_decided': 0.8824110671936759, 'false_agreement_rate_among_both_decided': 0.8431018935978359, 'correct_separation_rate': 0.08389261744966443, 'abstain_pair_frac': 0.2557046979865772}, {'complement': 'complement_a_O1_abstain', 'score': 1.3950431317551684, 'soft_separation_rate': 0.22416107382550335, 'Q_among_decided': 0.7777777777777778, 'false_agreement_rate_among_both_decided': 0.8310567936736161, 'correct_separation_rate': 0.032214765100671144, 'abstain_pair_frac': 0.06644295302013423}, {'complement': 'baseline_pred_O0', 'score': 1.2349451376976295, 'soft_separation_rate': 0.1738255033557047, 'Q_among_decided': 0.7409695817490495, 'false_agreement_rate_among_both_decided': 0.8536754507628294, 'correct_separation_rate': 0.08389261744966443, 'abstain_pair_frac': 0.032214765100671144}]`
- Winner = best Mode B decision complement on myeloid↔T must-pairs (soft_separation + Q_among_decided − false_agreement). z⁺ geometry reported separately; prior showed it worsens global must-rate.

## Claim boundary

Mode B only (sufficiency-for-readout / typed evidence / observer edit). Not Mode A residual-as-field, Jacobian, gene perturbs, fusion audit, or clinical. Do not claim Pearson > ~0.61. True-ADT events are an explicit complement probe, not the default no-leakage ANM path.
