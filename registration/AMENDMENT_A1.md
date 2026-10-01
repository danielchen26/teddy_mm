# Amendment A1 to the v3 pre-registration: corrections from the adversarial review

**Status:** written before any site4 (test) evaluation of v3; every number in 'computed' uses train and val cells only (v3_leakage_check.py). File `registration/amendment_A1.json`, sha256 in `amendment_A1.json.sha256`. It amends `registration/registration_v3.json` sha256 `e4c8a33e5c7c0b292d3a70a9063730c7ea0fd6d1519a705953bcf73ad8212bd3`, which is unchanged. This page is rendered from the JSON by `bridge_anm/v3_build_amendment.py`.

**Scope.** registration_v3.json is unchanged; where A1 and registration_v3.json differ, A1 wins. Builders load both with v3_amend.load_registration_amended().

## Findings and fixes

### A1.1 matched coverage (E1.1c/d, E1.3c, E1.4, E1.5, C2, E3.H3b, E4)

- **Finding.** The evidence is clipped at 1, so many cells share the top score (val: 0.3662 of cells have a Q1 top score of exactly 1, 0.5314 for Q3), and the site4 rows are stored in 3 contiguous donor blocks; the registered tie rule (global cell index) would call one donor's tied cells first. Whether matched-coverage selections are redone in each bootstrap replicate was not stated.
- **Fix.** Ties are broken by the rank of the global cell index in numpy default_rng(29).permutation(90261), the same for every method (v3_amend.coverage_select); every endpoint, including the matched-coverage selection and any realised coverage that defines a matched point, is recomputed in every bootstrap replicate.

### A1.2 E4 channel 2

- **Finding.** The stored ADT is a per-cell CLR (expm1 of a row is an integer count vector over one per-cell factor; 1.0000 of 2000 sampled train/val cells pass that check), and the factor is computed from all proteins, the 12 panel proteins included. So the 122 stored non-panel values channel 2 reads carry a function of the measured panel counts, contrary to 'it never sees a panel protein's measured value'. Val ridge diagnostic (mean R2 over the 12 panel targets): stored inputs 0.6480, panel-free inputs 0.6379, panel-free inputs plus the CLR scalar 0.6416.
- **Fix.** Channel 2 reads w_q = log1p(1e4 * u_q / sum_r u_r), u = expm1(stored value), q and r over the 122 non-panel proteins (v3_amend.channel2_inputs). Because u = x / g_c, w equals log1p(1e4 * x_q / sum_r x_r) and depends on the non-panel counts only (unit test: changing a cell's panel counts leaves w unchanged).

### A1.3 E1.3 why this call

- **Finding.** ANM's leave-one-out deletes the event's site (build_graph makes one site per event; v2's _loo_top_for_instance drops the event), so the called action's gain changes from G(n) to G(n-1). The registered closed form (S_called - v_top/n against the bar) re-codes zeroing, not deletion; rho = G(n-1)/G(n) = 0.992479 for n = 3, and on val 6 of 5190 Q1 calls change flip status. Also, 0.5333 of val Q1 calls have two or more panel markers tied at the called class's maximum (mostly clipped at 1), so a single 'top-1 marker' would be set by event order. The question (Q1 or Q2) for E1.3 was not stated.
- **Fix.** Exact closed form flip iff rho * (S_called - v_top / n) < max(bar, max_{k != called} S_k) (v3_amend.loo_flip_exact); top markers are sets (all markers at the class maximum), E1.3a compares sets, E1.3b gives each marker of a set of m markers 1/m credit (and reports single-marker calls separately); E1.3 runs on Q1 (primary) and Q2 (secondary).

### A1.4 E1.4 honest naming

- **Finding.** The ANM soft score is G(3) * 3 times the rule's top score, so its ranking and every E1.4 number equal the declared rule's; the registered text read an E1.4 win as support for an ANM gate claim. Q1 and Q2 share scores and key, so their E1.4 numbers are identical.
- **Fix.** An E1.4 win is reported as a win of the top-score readout (rule = ANM), never as an ANM-specific property; E1.4 is computed once (Q1 = Q2); the E1 falsification sentence is restated accordingly.

### A1.5 E1.C2 primary point

- **Finding.** C2 listed three coverage points without saying which decides the verdict, and did not define 'in-scope selective accuracy'.
- **Fix.** The verdict uses c* = the mean rule's realised Q1 coverage on test_primary (label-free); 0.8 and 0.7 are secondary. In-scope selective accuracy = among called cells whose key class is a lineage, the fraction correct; guard = point difference >= -0.005 at c*.

### A1.6 E1.1c / E1.1d / E1.5 information parity

- **Finding.** The Q3 rule reads CD19, CD3, CD56 and CD14 evidence; the registered 'Q1 classifier without relabelling' reads the 12 primary-panel values only (no CD19, no CD14), while the Q3 classifier reads 14 values. E1.5 therefore changed features as well as labels between n = 0 and n > 0.
- **Fix.** The Q1-label arm for Q3 (E1.1c and E1.5 n = 0) is a classifier on the classifier.q3 features (14 values, a superset of the rule's 4 anchors) trained on Q1 labels, C = 100.0 on val; features are the same at every n. Secondary matched-information rows use exactly the rule's 4 anchor values: Q1-label C = 10.0, Q3-label C = 100.0 (val).

### A1.7 E2 falsification (i)

- **Finding.** Criterion (i) compares every pair of perturbations without a minimum count of lost cells or an uncertainty requirement, so with many pairs and few lost cells per pair noise alone can pass it.
- **Fix.** Only perturbation levels with at least 30 lost cells among the 1,000 test_primary e2_subset cells (in both perturbations) enter (i), and the share difference needs a two-stage bootstrap 95% interval excluding 0; the number of pairs compared is reported.

### A1.8 E3 pairs

- **Finding.** Neighbours 'within the evaluated split' pool the two primary donors, so cross-donor pairs exist, while the registered bootstrap resamples pairs within donors and the flag cosine came from one donor.
- **Fix.** K = 10 cosine neighbours are computed within each donor of the evaluated split; every pair belongs to one donor.

### A1.9 E4 comparator

- **Finding.** 'F5 vs the best of F1-F4' did not say how 'best' is chosen; choosing it on site4 is a post-hoc maximum.
- **Fix.** The comparator is the F1-F4 method with the highest val endpoint (selective accuracy at coverage 0.8, mean over L0-L3; F3 scored by 5-fold cross-fitting on val with seeds.e4_noise), written to the E4 addendum before site4.

### A1.10 addenda

- **Finding.** The untracked E5 addendum replaces the registered primary endpoint (log-gain at l*) with dG_l and uses pooled site4 annotation-key pairs; the registration allows addenda only to fix numbers it leaves open.
- **Fix.** An addendum cannot replace a registered endpoint, key or cell set; the registered E5 endpoint stays primary and dG_l is secondary unless a further amendment is committed before site4.

## Numbers computed for A1 (train and val only)

| item | value |
|---|---|
| tie-break seed / cells | 29 / 90,261 |
| val share with Q1 top score exactly 1 | 0.3662 |
| val share with Q3 top score exactly 1 | 0.5314 |
| site4 donor blocks in row order (metadata) | 3 |
| rho = G(2)/G(3) | 0.992479 |
| val Q1 calls / with tied top marker | 5,190 / 0.5333 |
| val Q1 calls whose flip status the exact form changes | 6 |
| sampled cells passing the per-cell CLR check | 1.0000 of 2,000 |
| val ridge mean R2: stored / panel-free / panel-free + CLR scalar | 0.6480 / 0.6379 / 0.6416 |
| classifier q1labels_q3features (14 features): C | 100.0 (val log-loss 0.01: 0.073346, 0.1: 0.046275, 1.0: 0.040072, 10.0: 0.035518, 100.0: 0.033674) |
| classifier q1labels_anchors (4 features): C | 10.0 (val log-loss 0.01: 0.180354, 0.1: 0.092885, 1.0: 0.065008, 10.0: 0.062309, 100.0: 0.062622) |
| classifier q3labels_anchors (4 features): C | 100.0 (val log-loss 0.01: 0.127127, 0.1: 0.063196, 1.0: 0.053331, 10.0: 0.050793, 100.0: 0.050537) |

## Experiment fields replaced (path in registration_v3.json → experiments)

- `common/matched_coverage/rule`: at coverage c each method calls its ceil(c * N) most confident cells of the evaluated split (N = all cells of the split), using only its own scores, never labels; ties in a method's score are broken by the A1.1 tie-break (rank of the global cell index in numpy default_rng(29).permutation(90261)), the same permutation for every method (v3_amend.coverage_select)
- `common/matched_coverage/tie_break`: {"seed": 29, "n_cells": 90261, "function": "v3_amend.coverage_select"}
- `common/statistics/bootstrap`: two-stage: resample donors with replacement, then cells (or pairs) with replacement within each drawn donor; B = 2000; seed = registration seeds.bootstrap; percentile 95% interval. Every endpoint is recomputed inside each replicate, including the matched-coverage selection (N = replicate size; a repeated cell keeps the A1.1 rank of its original id) and any method's realised coverage that defines a matched point (A1.1)
- `E1/exp1_change_the_question/arms`: ["rule (Q3 anchors CD19, CD3, CD56, CD14)", "anm", "Q1-label classifier for Q3 (A1.6): classifier.q3 features, Q1 labels, C = 100.0; its 'myeloid' is read as the Q3 myeloid call", "Q3-label classifier (registration classifier.q3; all training-site Q3 labels; label-rich ceiling)", "secondary (A1.6, same information as the rule): the two classifiers on the 4 anchor values only (Q1 labels C = 10.0, Q3 labels C = 100.0)"]
- `E1/exp1_change_the_question/endpoints/E1.1c`: Q3 selective accuracy on the q3 key at the rule's realised Q3 coverage: rule vs the A1.6 Q1-label classifier (classifier.q3 features) without relabelling; secondary row with the 4-anchor Q1-label classifier
- `E1/exp1_change_the_question/endpoints/E1.1d`: same, rule vs the Q3-label classifier (label-rich ceiling); secondary row with the 4-anchor Q3-label classifier
- `E1/exp3_why_this_call/question`: Q1 on test_primary (primary); Q2 with its bar (secondary)
- `E1/exp3_why_this_call/anm`: leave-one-out over the called class's events (events at t = 0): removing an event deletes its event site (ANM build_graph makes one site per event, as v2's _loo_top_for_instance); top-1 marker set = every event whose removal lowers the called action's score most (ties kept as a set); flip = removing a top-1 event changes the call (to another class or to no call)
- `E1/exp3_why_this_call/rule`: exact closed form (A1.3, v3_amend.loo_flip_exact): top-1 marker set = the panel proteins at the called class's maximum evidence; flip iff rho * (S_called - v_top / n) < max(bar, max over other classes of S_k), rho = G(n-1)/G(n) = 0.992479 for n = 3
- `E1/exp3_why_this_call/endpoints`: {"E1.3a": "agreement of ANM's top-1 marker set with the closed-form set (pre-registered expectation 1.000; below 0.99 is reported as a discrepancy), and agreement of the flip status (expectation 1.000); the share of calls whose top-1 set has more than one marker is reported", "E1.3b": "per-donor share of each top-1 marker with fractional credit (a set of m markers gives 1/m to each); null = within-class permutation of marker identities per cell, 1000 permutations, same credit; also reported on calls with a single top-1 marker", "E1.3c": "error rate of flip-sensitive calls (exact flip, A1.3) vs the same number of lowest-margin calls (margin ties by the A1.1 permutation), Q1, test_primary"}
- `E1/exp4_which_calls_to_trust/interpretation`: the ANM soft score is G(3) * 3 times the rule's top score, so its ranking and every E1.4 number equal the declared rule's top score; an E1.4 win is reported as a win of the top-score readout (rule = ANM), never as an ANM-specific property. Q1 and Q2 share scores and key, so E1.4 is computed once and reported once
- `E1/exp4_which_calls_to_trust/failure`: if the top score (rule = ANM) is not a win against the margin, the gate claim stays withdrawn; a win supports the top-score readout, not an ANM-specific mechanism
- `E1/c2_nested_readouts/primary_coverage`: the verdict uses c* = the mean rule's realised Q1 coverage on test_primary (label-free; recomputed in each bootstrap replicate); coverage 0.8 and 0.7 are secondary and never change the verdict
- `E1/c2_nested_readouts/in_scope_accuracy`: among called cells whose key class is B, T, NK or myeloid, the fraction whose call equals the key; the guard is the point difference (ANM closure minus mean rule) >= -0.005 at c*
- `E1/falsification`: Mode B cannot show ANM-specific decision value: every ANM arm equals a re-coded rule cell by cell (E1.1a; the closure rule in C2; the top score in E1.4). E1 tests whether the readout forms ANM provides (top action score, closure readout) add decision value over TEDDY's margin and the mean rule; that is falsified for v3 if neither E1.4 nor C2 is a win, and any win is credited to the readout form, which a written rule implements
- `E1/exp5_label_cost/classifier`: features are the same at every n (A1.6): registration classifier.q3 features; n = 0 is the A1.6 Q1-label classifier on them (C = 100.0); n > 0 uses classifier.q3 (C = registration classifier.q3.C) on the drawn Q3 labels; a draw with a single class predicts that class. Secondary curve on the 4 anchor values only (n = 0: Q1 labels, C = 10.0; n > 0: C = 100.0)
- `E2/falsification`: the decomposition adds information beyond the curve only if (i) two perturbations with accuracy losses within 0.02 of each other, each with at least 30 lost cells among the 1,000 test_primary e2_subset cells, differ by >= 0.15 in a localisation share with a two-stage bootstrap 95% interval of the difference excluding 0 (the number of pairs compared is reported), or (ii) the small-eps response (dS at eps = 0.05 thinning, linearly extrapolated) predicts which cells change call at eps = 0.8 with AUROC >= 0.02 above baseline margin alone (bootstrap lower bound > 0); if neither holds, it adds nothing beyond the curve
- `E3/cells`: test_primary; pairs (A1.8): k = 10 cosine neighbours on the raw L2-normalised final z, computed within each donor of the evaluated split; each edge kept once; NK-T pair = one primary-key NK and one primary-key T cell of the same donor; flagged by registration e3.flag_cosine (from val); variants all and no_gdT158
- `E4/channels/channel2`: ridge regression (alpha on val from [0.1, 1, 10, 100, 1000]) from the panel-free renormalised non-panel ADT (A1.2, v3_amend.channel2_inputs: w_q = log1p(1e4 * u_q / sum_r u_r), u = expm1(stored value), q and r over the 122 non-panel proteins) to the measured panel evidence; trained on training cells; normalised by the training q95 of its own predictions; it reads no function of a panel protein's measured value
- `E4/noise_levels/L2`: 50% of the 122 channel-2 inputs w set to 0 per cell (seeds.e4_noise), after the A1.2 transform
- `E4/best_comparator`: the F1-F4 method with the highest val endpoint (selective accuracy at coverage 0.8, mean over L0-L3; F3 scored by 5-fold cross-fitting on val with seeds.e4_noise), fixed in the E4 addendum before site4 (A1.9)
- `common/addenda_scope`: an addendum may only fix numbers the registration leaves open, by the procedure written for them; it cannot replace a registered endpoint, key or cell set (A1.10). The registered E5 endpoint stays primary unless a further amendment is committed before site4

## Review notes (no change needed)

- Class map: all 45 types (40 in site4) reviewed. Every assignment is defensible. pDC -> OUT is the most debatable (some schemes count pDC as myeloid DC); it only enters the primary key when no lineage gate passes, and annotation_only is reported.
- Myeloid panel member CD62P is a platelet protein; in E6 the external set has a 'Platelet' label keyed OUT, so myeloid calls on platelets are a predictable E6 error mode (report it, do not change the panel).
- canonical9 sensitivity panels were chosen on site4 in v1/v2 (already disclosed); they never decide a verdict.
- Splits: train donors 10886, 11466, 12710, 15078, 16710, 28045 (sites 1-3); val 18303 (site1 only); test_primary 13272, 19593 (site4 only, in no other split); test_secondary 15078. No donor-split error found.
- Leakage: no site4 protein, label or embedding enters the registration core or this amendment core (v3_leakage_check.py). The stored ADT CLR is per cell, so no cross-cell (site4) information enters train values.
- E6 prep (commit 8318ef4): L2_CLASS is written from label names, but it was committed with the code that reads external values (already disclosed); the E6 runner must report this.

## Provenance

Builder `v3_build_amendment A1 1.0` sha256 `90181e24c46e316fc2f032dfb82ad0b7a56ddcc538bd883f6f75576d398e9231`; `v3_amend.py` sha256 `0ae5d516236c916283ea83104389b4ecdfc8504d2624f12f9eb9c3b5763d5736`; computed core sha256 `b069f310859bc9d5a073c25563f5bb6b367e41ca3a18626959c8df937a7272b8`.
