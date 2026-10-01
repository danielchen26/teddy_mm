# Amendment A2 to the v3 pre-registration: E3 H3a scale artefact, E3 R2 normaliser sample, E2 donor rule

**Status:** written before any site4 (test) evaluation of E2 and E3; the only numbers computed here use split=train rows (protein and embedding) and the split column (leakage_check). File `registration/amendment_A2.json`, sha256 in `amendment_A2.json.sha256`. It follows amendment A1 (`registration/amendment_A1.json` sha256 `1f44be66bfc22a5bc428cf31a0bcd14047466005f4308386da195f55b767cdee`) and amends `registration/registration_v3.json` sha256 `e4c8a33e5c7c0b292d3a70a9063730c7ea0fd6d1519a705953bcf73ad8212bd3`; both are unchanged. This page is rendered from the JSON by `bridge_anm/v3_build_amendment_A2.py`.

**Scope.** registration_v3.json, amendment_A1.json and the addenda are unchanged; builders apply A1 then A2 (v3_amend.load_registration_amended()). Where A2 differs from the registration as amended by A1, or from the E2 / E3 addenda (supersedes_in_addenda), A2 wins. A2 only makes verdicts stricter (A2.1, A2.3) or estimates the same quantity with less compute (A2.2).

## Findings and fixes

### A2.1 E3.H3a (R1 and R2): scale artefact of the registered D

- **Found by.** The E3 builder, before any site4 evaluation (registration/addenda/E3.json declared.secondary_H3a).
- **Finding.** The registered D_R = [GR_R - GR_head](flagged) - [GR_R - GR_head](unflagged) - D_null is not scale-free: a readout that shrinks every NK-T evidence difference by one common factor s < 1 gets D_R = (1 - s) * (GR_head(unflagged) - GR_head(flagged)) > 0 whenever the head already keeps less of the gap on flagged pairs, and the seed-1 null does not remove this; such a readout repairs nothing on flagged pairs yet can pass the registered win rule.
- **Fix.** The registered D_R computation, margin 0.05 and win / loss / equivalent rules are kept. A readout's H3a win now additionally requires, under the same registered two-stage donor bootstrap (B = 2000, seed 1, the same replicates as D_R) and per-donor rule: (b) recovery on flagged pairs: term_flagged_R = GR_R(flagged) - GR_head(flagged) > 0, its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor; (c) the scale-free contrast Dlog_R = [ln GR_R(fl) - ln GR_head(fl)] - [ln GR_R(unfl) - ln GR_head(unfl)] - (the same for the null), as the E3 builder computes it (addendum secondary_H3a S1), > 0 with its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor. An undefined value fails. If D_R is a win but (b) or (c) fails, the H3a verdict for that readout is 'not supported (scale artefact)'; a loss, equivalent or inconclusive D_R verdict is unchanged. Falsification: the claim is rejected when neither R1 nor R2 has an A2 win. H3b is unchanged; H3a_pooling (D_R2 - D_R1) is unchanged.
- **Effect.** Stricter: every A2 win is a registered win (v3_amend.a2_h3a_verdict).

### A2.2 E3 R2 normaliser (compute)

- **Found by.** The E3 builder: r2_predict over all 67,405 split=train cells is about 2.2 h of TEDDY forward at the official rate, spent on one percentile per target.
- **Finding.** The E3 addendum (declared.readout_evidence) normalises R2 by the 95th percentile of its predictions over all 67,405 split=train cells; a seeded uniform random sample of training cells estimates the same training-population quantity with a sampling error that is small against the H3a margin (computed.q95_sampling_check).
- **Fix.** R2's normaliser = np.percentile(95) of R2's predictions over the A2.2 sample: 10,000 split=train cells, numpy default_rng(31).choice(split == 'train' rows in index order, 10,000, replace=False), sorted; label-free (the split column only), drawn here before any R2 prediction; cell list and sha256 (ee29a7c4cf2928a9...) in computed.e3_r2_normaliser_sample. r2_predict predicts these training cells (shards trainA2_NNNN) instead of every training cell; val and site4 predictions are unchanged. The 4,000-cell R2-subset normaliser stays the registered sensitivity; the unclipped sensitivity is unchanged.
- **Effect.** Same estimand; training-cell forward 0.33 h instead of 2.2 h. Sampling check (train rows only, computed.q95_sampling_check): over the 4 gap proteins the A2 sample's q95 differs from the all-training q95 by at most 0.0221 (relative) on the head and the ridge proxy of R2, the design SD of a 10,000-cell q95 is at most 0.0273 (relative), and the largest |z| over all 13 targets is 1.14 (the A2 sample's q95 is within sampling error of the all-training q95). Implied shift of D_R2 per unit of max_p |r_p(flagged) - r_p(unflagged)|: 0.0062 observed (ridge proxy of R2; head 0.0221), against the margin 0.05. The registered 4,000-cell subset is stratified, and its q95 differs from the all-training q95 by up to 0.32 (relative) on the gap proteins, which is why it is only a sensitivity.

### A2.3 E2 falsification criteria (i) and (ii): common decision rule

- **Found by.** Review of the E2 builder before any site4 evaluation (registration/addenda/E2.json open_choices_fixed.verdict and interval_for_i report the per-donor clause beside the verdict, as secondary).
- **Finding.** The registration's decision rules apply to all experiments (win = point difference >= margin, bootstrap lower bound > 0, and the difference has the same sign and is >= margin / 2 in each primary donor), but the E2 builder computed the E2 verdict without the per-donor clause and reported that clause as secondary.
- **Fix.** The common decision rule is binding for E2's criteria (i) and (ii). (i): a qualifying pair (A1.7: decision-accuracy losses within 0.02, >= 30 lost test_primary e2_subset cells in each unit, on the point estimates) and label counts only if |share difference| >= 0.15, its two-stage 95% interval excludes 0 on the side of the difference, and in each primary donor (13272, 19593) the donor's share difference has the same sign and |difference| >= 0.075 (an undefined donor share fails). (ii): holds only if dAUROC >= 0.02, its lower bound > 0 and dAUROC >= 0.01 in each primary donor. Verdict: 'adds information beyond the curve' if (i) or (ii) holds under this rule, otherwise 'adds nothing beyond the curve (falsified)'. The A1-only reading without the donor clause is reported beside it and never decides.
- **Effect.** Stricter: (i) and (ii) under A2 imply (i) and (ii) under A1 (v3_amend.a2_e2_criterion_i / _ii).

## A2.2 sample and sampling check (training rows only)

Sample: 10,000 of 67,405 split=train cells (0.1484), seed 31, sha256 `ee29a7c4cf2928a95636f9795af1d5513d206183dba679962170257263486a75` (population sha256 `0da6d39783b7ef2121d44e85d47a1a64867ca0b7dff3550245df77295dd59c1a`, as registration splits.train). The 4,000-cell R2 subset is stratified (1/3 NK, 1/3 T, 1/3 other), so its q95 is a biased estimate of the training-population q95; it stays a sensitivity.

q95 (np.percentile 95) over the A2.2 sample vs over all split=train cells, per E3 target, for two readouts that exist before R2 is trained; design_sd_10k = SD of (q95 of a 10,000-cell draw - q95 of all training cells) over other seeded draws (the A2.2 design); z = (A2 sample - all) / that SD. rel_diff_R2_4000_subset = the registered 4,000-cell R2-subset normaliser (stratified 1/3 NK, 1/3 T, 1/3 other) against all training cells, for context.

| readout | target | q95 all train | q95 A2 sample | rel. diff | design rel. SD (10k) | z | rel. diff of the 4,000 subset |
|---|---|---:|---:|---:|---:|---:|---:|
| head | CD20 | 1.4892 | 1.4892 | +0.0001 | 0.0104 | +0.01 | -0.4468 |
| head | CD22 | 1.1521 | 1.1519 | -0.0002 | 0.0282 | -0.01 | -0.7910 |
| head | CD268 | 0.4738 | 0.4712 | -0.0054 | 0.0184 | -0.30 | -0.5934 |
| head | CD3 | 2.0965 | 2.1021 | +0.0026 | 0.0108 | +0.24 | +0.1155 |
| head | CD2 | 2.3828 | 2.3907 | +0.0033 | 0.0074 | +0.44 | +0.0392 |
| head | CD5 | 2.8648 | 2.8565 | -0.0029 | 0.0039 | -0.74 | +0.0376 |
| head | CD122 | 1.2522 | 1.2278 | -0.0195 | 0.0332 | -0.59 | +0.3307 |
| head | CD94 | 2.0649 | 2.0377 | -0.0132 | 0.0158 | -0.83 | +0.2777 |
| head | CD56 | 2.0272 | 1.9824 | -0.0221 | 0.0273 | -0.81 | +0.3204 |
| head | CD172a | 1.6673 | 1.6586 | -0.0052 | 0.0064 | -0.82 | -0.1132 |
| head | CD11c | 1.6504 | 1.6431 | -0.0045 | 0.0092 | -0.48 | -0.1540 |
| head | CD62P | 0.9869 | 0.9859 | -0.0010 | 0.0059 | -0.18 | -0.1448 |
| head | CD335 | 1.2353 | 1.2327 | -0.0021 | 0.0143 | -0.14 | +0.2744 |
| ridge_proxy_of_R2 | CD20 | 0.8300 | 0.8301 | +0.0001 | 0.0111 | +0.01 | -0.2638 |
| ridge_proxy_of_R2 | CD22 | 0.7800 | 0.7787 | -0.0017 | 0.0160 | -0.11 | -0.6772 |
| ridge_proxy_of_R2 | CD268 | 0.6987 | 0.6971 | -0.0023 | 0.0131 | -0.18 | -0.5128 |
| ridge_proxy_of_R2 | CD3 | 0.7994 | 0.7955 | -0.0049 | 0.0097 | -0.51 | +0.0812 |
| ridge_proxy_of_R2 | CD2 | 0.8185 | 0.8118 | -0.0082 | 0.0144 | -0.57 | +0.0981 |
| ridge_proxy_of_R2 | CD5 | 0.8321 | 0.8284 | -0.0045 | 0.0089 | -0.50 | +0.0661 |
| ridge_proxy_of_R2 | CD122 | 0.7716 | 0.7721 | +0.0006 | 0.0068 | +0.09 | +0.1758 |
| ridge_proxy_of_R2 | CD94 | 0.7476 | 0.7456 | -0.0027 | 0.0050 | -0.53 | +0.1417 |
| ridge_proxy_of_R2 | CD56 | 0.8035 | 0.8085 | +0.0062 | 0.0087 | +0.71 | +0.1591 |
| ridge_proxy_of_R2 | CD172a | 0.8626 | 0.8554 | -0.0083 | 0.0073 | -1.14 | -0.1314 |
| ridge_proxy_of_R2 | CD11c | 0.8409 | 0.8383 | -0.0030 | 0.0066 | -0.46 | -0.1346 |
| ridge_proxy_of_R2 | CD62P | 0.6888 | 0.6867 | -0.0032 | 0.0081 | -0.39 | -0.1748 |
| ridge_proxy_of_R2 | CD335 | 0.7427 | 0.7427 | -0.0000 | 0.0049 | -0.00 | +0.1500 |

Implied effect on E3.H3a: GR_R2 on a pair set is a mean over the 4 gap proteins of a median difference divided by R2's q95 for that protein (before clipping), so a relative error d_p of that q95 moves the protein's ratio r_p by about -d_p r_p and D_R2 by about -mean_p d_p (r_p(flagged) - r_p(unflagged)); |shift of D_R2| <= max_p |d_p| * max_p |r_p(flagged) - r_p(unflagged)|. Values below are max_p |d_p| over the gap proteins, i.e. the shift per unit of max_p |r_p(flagged) - r_p(unflagged)|; the ridge proxy is the readout closest to R2 (R2 starts as a linear readout of the gene-mean), the head is shown for a heavier-tailed prediction distribution. head: observed 0.0221, two design SDs 0.0546; ridge_proxy_of_R2: observed 0.0062, two design SDs 0.0194; margin 0.05.

Compute: training-cell forward 0.33 h instead of 2.20 h at the official rate.

## Leakage check

As A1 (v3_leakage_check.py) and the E2 / E3 addenda: rebuild the core from copies of the inputs whose rows of one split carry random protein, cell types and embedding, and compare bytes.

| poisoned rows | role | core identical to real |
|---|---|---|
| test (seed 101) | site4 rows poisoned: core must be byte-identical | True |
| val (seed 303) | val rows poisoned: A2 reads no val value, so identical as well (information) | True |
| train (seed 202) | positive control, training rows poisoned: core must change | False |

Result: **PASS**; real core sha256 `561842f206e06e64ea5aa7ffd4335b44dc954585512c449cff3c270f25c135bf`.

## Experiment fields replaced (path in registration_v3.json → experiments, after A1)

- `E3/endpoints/E3.H3a`: D_R = [GR_R(flagged) - GR_head(flagged)] - [GR_R(unflagged) - GR_head(unflagged)] for R in {R1, R2}, minus the same quantity for the null; pre-registered: >= 0.05 with win rules (pairs resampled within donors). A2.1: a readout's H3a win additionally requires, under the same two-stage donor bootstrap and per-donor rule, (b) recovery on flagged pairs: term_flagged_R = GR_R(flagged) - GR_head(flagged) > 0, its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor and (c) the scale-free contrast Dlog_R = [ln GR_R(fl) - ln GR_head(fl)] - [ln GR_R(unfl) - ln GR_head(unfl)] - (the same for the null), as the E3 builder computes it (addendum secondary_H3a S1), > 0 with its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor; a D_R win without (b) and (c) is 'not supported (scale artefact)'
- `E3/gap_ratio`: per readout and protein p in {CD56, CD94, CD335, CD3}: median |predicted evidence difference| / median |measured evidence difference| within pairs; GR = mean over the 4 proteins; each readout normalised by the training q95 of its own predictions (A2.2: R2's training q95 is estimated on the registered uniform sample of 10,000 split=train cells, numpy default_rng(31), cell list in amendment_A2.json; the 4,000-cell R2-subset normaliser is a sensitivity)
- `E3/falsification`: if neither R1 nor R2 has an A2 win for E3.H3a (D_R >= 0.05 over the null with the win rules, plus (b) recovery on flagged pairs and (c) a positive scale-free Dlog_R, A2.1), 'the NK-T loss is in the readout or pooling and is repairable on frozen TEDDY' is rejected for this dataset
- `E2/falsification`: the decomposition adds information beyond the curve only if (i) two perturbations with accuracy losses within 0.02 of each other, each with at least 30 lost cells among the 1,000 test_primary e2_subset cells, differ by >= 0.15 in a localisation share with a two-stage bootstrap 95% interval of the difference excluding 0 (the number of pairs compared is reported), or (ii) the small-eps response (dS at eps = 0.05 thinning, linearly extrapolated) predicts which cells change call at eps = 0.8 with AUROC >= 0.02 above baseline margin alone (bootstrap lower bound > 0); if neither holds, it adds nothing beyond the curve; A2.3: the registration's common decision rule is binding for (i) and (ii): in (i) the share difference must also have the same sign and be >= 0.075 in each primary donor (13272, 19593), with the interval excluding 0 on the side of the difference; in (ii) dAUROC must also be >= 0.01 in each primary donor

## Addendum fields superseded (files unchanged)

- `registration/addenda/E3.json` (sha256 `09f932ed166cad02...`): declared.readout_evidence: R2's normaliser over all 67,405 split=train cells -> the A2.2 sample; declared.falsification: 'rejected unless E3.H3a for R1 or for R2 is a win' -> an A2 win; declared.secondary_H3a: S1 (Dlog_R) and the flagged term of S2 also enter the A2.1 win conditions
- `registration/addenda/E2.json` (sha256 `d1f0a5077640ecae...`): spec.open_choices_fixed.verdict and interval_for_i: the per-donor clause of the common win rule is binding (A2.3), no longer secondary; the Bonferroni interval stays secondary

## Change log

| file | change | reason |
|---|---|---|
| `registration/amendment_A2.json (+ .sha256, amendment_A2_core.json, AMENDMENT_A2.md)` | new amendment A2 | A2.1-A2.3 above; written before any site4 evaluation of E2 and E3 |
| `bridge_anm/lib/v3_amend.py` | load_registration_amended() applies A1 then A2, verifying the registration, A1 and A2 hashes and that A2 names the registration and A1 on disk; amendment_A2_status() for site4 guards; A2.1 / A2.2 / A2.3 functions | every v3 builder reads the amended registration through this loader |
| `scripts/v3_e3_nkt_repair.py` | site4 stages (r2_predict on site4, evaluate) refuse unless A2 is committed and matches; r2_predict predicts the A2.2 sample instead of all training cells; evaluate normalises R2 by the A2.2 sample, computes (b), (c) and the A2 verdicts on the same replicates, and states the A2 falsification; DECLARED (the addendum's declared part) is unchanged | A2.1, A2.2; earlier stages (null, r1, r2_states, r2_fit) are unchanged |
| `scripts/e2_response_decomposition.py` | site4 stages refuse unless A2 is committed and matches; the report computes (i) and (ii) under the common decision rule and decides the verdict with it; SPEC (the addendum's spec) is unchanged | A2.3 |
| `registration/HASHES.txt` | A2 lines; v3_amend.py line updated to the A2 version (the A1 version's sha256 is in amendment_A1.json provenance) | hash record of the registration files |

## Provenance

Builder `v3_build_amendment_A2 1.0` sha256 `6801499dd07d6144da17fbf7e8b6bb83a5a92ee41f5d23564d329276bf20dd63`; `v3_amend.py` sha256 `59089ef5f9ec8ecb0a4af1d82f37d4fc2b1bd316cc5aeba671877d7ef44abff1`; computed core sha256 `561842f206e06e64ea5aa7ffd4335b44dc954585512c449cff3c270f25c135bf`; built 2026-10-01T06:09:15Z.
