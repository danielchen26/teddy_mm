# Amendment A3 to the v3 pre-registration: E3 H3a_pooling scale artefact

**Status:** written before any site4 (test) evaluation of E3; the only numbers computed here use split=train rows (protein, cell types, embedding) and the split and donor columns (leakage_check). File `registration/amendment_A3.json`, sha256 in `amendment_A3.json.sha256`. It follows amendment A2 (`registration/amendment_A2.json` sha256 `20b64174b063d39997abf09877ae1ba152d935589b43630e15dd6ebe6d21fd36`) and A1 (`registration/amendment_A1.json` sha256 `1f44be66bfc22a5bc428cf31a0bcd14047466005f4308386da195f55b767cdee`) and amends `registration/registration_v3.json` sha256 `e4c8a33e5c7c0b292d3a70a9063730c7ea0fd6d1519a705953bcf73ad8212bd3`; all three are unchanged. This page is rendered from the JSON by `bridge_anm/v3_build_amendment_A3.py`.

**Scope.** registration_v3.json, amendment_A1.json, amendment_A2.json and the addenda are unchanged; builders apply A1, A2, then A3 (v3_amend.load_registration_amended()). Where A3 differs from the registration as amended by A1 and A2, from A2 (supersedes_in_amendments) or from the E3 addendum (supersedes_in_addenda), A3 wins. A3 only makes one verdict stricter (A3.1) and concerns E3 only.

## Disclosure

- **E1 site4 results seen.** Before A3 was written, the orchestrating session had already seen the E1 site4 (test) results (E1 is the Mode B rerun). A3 concerns only E3 (the E3.H3a_pooling verdict) and follows from the algebra of the registered statistic (A2.1's algebra with R1 in the head's place); it reads no E1 output and nothing in it depends on an E1 result.
- **E3 site4 state.** No E3 site4 evaluation had been run: E3 evaluate, the only E3 stage that reads site4 protein or cell types, needs R2's site4 predictions, which E3 r2_predict (a TEDDY + R2 forward of site4 RNA; no protein or label) was still writing while A3 was written; A3 and its builder read none of them. From A3 on, E3 evaluate on site4 is refused unless A3 is committed and matches.
- **E3 val smoke seen.** The val-only smoke runs of E3 evaluate (val donor 18303, small subsets, labelled smoke) had been seen; A3's rule does not depend on them.
- **numbers computed for A3.** The only numbers computed for A3 use split=train rows (protein, cell types, embedding) and the split and donor columns (leakage_check).

## Findings and fixes

### A3.1 E3.H3a_pooling (D_R2 - D_R1, 'the loss is in mean pooling'): scale artefact of the registered contrast

- **Found by.** Review of E3 before any site4 evaluation of E3: A2 left E3.H3a_pooling on the registered rule (A2.1 'H3a_pooling (D_R2 - D_R1) is unchanged'), and commit b4ecad4 then added the A2.1-style conditions (b') and (c') beside the registered pooling verdict as a secondary that never decides (scripts/v3_e3_nkt_repair.py, run_h3a / pooling_scale_check, decides: False).
- **Finding.** The head and the null cancel in the registered pooling contrast: D_R2 - D_R1 = [GR_R2(fl) - GR_R1(fl)] - [GR_R2(unfl) - GR_R1(unfl)], which is not scale-free: if R2 keeps a common fraction s < 1 of R1's NK-T evidence differences (GR_R2 = s GR_R1 on every pair set), D_R2 - D_R1 = (1 - s) * (GR_R1(unfl) - GR_R1(fl)) > 0 whenever R1 keeps less of the gap on flagged pairs. Such an R2 recovers less of the measured gap than R1 on flagged pairs, yet can pass the registered win rule and support 'the loss is in mean pooling'. This is the algebra of A2.1 with R1 in the head's place.
- **Fix.** The registered D_R2 - D_R1 computation, margin 0.05 and win / loss / equivalent rules are kept. The pooling claim wins only if the registered rule gives a win AND, under the same registered two-stage donor bootstrap (B = 2000, seed 1, the same replicates as D_R2 - D_R1) and per-donor rule: (b') recovery on flagged pairs beyond R1: term_flagged_R2_minus_R1 = GR_R2(flagged) - GR_R1(flagged) > 0, its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor; (c') the scale-free contrast Dlog_R2 - Dlog_R1 = [ln GR_R2(fl) - ln GR_R1(fl)] - [ln GR_R2(unfl) - ln GR_R1(unfl)] (Dlog_R as the E3 builder computes it, addendum secondary_H3a S1), > 0 with its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor. An undefined value fails. If the registered pooling verdict is a win but (b') or (c') fails, the E3.H3a_pooling verdict is 'not supported (scale artefact)'; a loss, equivalent or inconclusive registered verdict is unchanged. (b') and (c') are the statistics term_flagged_R2_minus_R1 and Dlog_R2_minus_Dlog_R1 that the E3 evaluate code computes since commit b4ecad4 with decides = False; under A3 they decide. E3.H3a (as amended by A2.1), E3.H3b and the E3 falsification (A2.1) are unchanged.
- **Effect.** Stricter: every A3 pooling win is a registered pooling win (v3_amend.a3_pooling_verdict). Training-pair illustration (computed; 1967 NK-T pairs of the 6 training donors, 1293 flagged): the head's GR is 0.415683 on flagged and 0.54728 on unflagged pairs, so with a head-like R1 a pure rescaling reaches the margin at 1 - s = 0.379949 (point estimate). R2 = 0.5 x R1 has D_R2 - D_R1 = 0.065798 with (b') -0.207842 and (c') 0.0: positive D, nothing repaired; a flagged-pair repair has (b') 0.584317 and (c') 0.877831. The head's flagged / unflagged difference varies in sign across training donors (computed.head_gap_ratio.per_donor), so the size of the artefact depends on R1; the rule does not.

## Illustration on training pairs (training rows only; point values; decides nothing)

Pairs: As E3 evaluate (addendum declared.pairs): k = 10 cosine neighbours of the raw official z within each training donor (v3_e3.nkt_pairs_within_donors), NK-T pair = one key-NK and one key-T cell (registration primary key) of the same donor, flagged if cosine >= 0.949402. 1967 NK-T pairs (1293 flagged) in the training donors 10886 161 (112), 11466 246 (156), 12710 135 (48), 15078 415 (179), 16710 228 (98), 28045 782 (700); pairs sha256 `6a686eeb8d1e9f13...`.

The registered phase-1 head's gap ratio GR (E3 addendum declared.gap_ratio; evidence = prediction / registered q95_train_pred, clipped; measured = ADT / addendum q95_measured_train, clipped) on the training pairs: flagged 0.415683, unflagged 0.54728, difference 0.131596. If R2 keeps a common fraction s of R1's NK-T evidence differences, D_R2 - D_R1 = (1 - s) (GR_R1(unflagged) - GR_R1(flagged)); with a readout like the head as R1 the registered pooling margin 0.05 is reached by a pure rescaling with 1 - s = 0.05 / (GR_head(unflagged) - GR_head(flagged)) when that difference is positive: 1 - s = 0.379949.

| training donor | GR_head flagged | GR_head unflagged | unflagged - flagged |
|---|---:|---:|---:|
| 10886 | None | 0.860957 | None |
| 11466 | 0.350414 | 0.449671 | 0.099257 |
| 12710 | 0.726923 | 0.925779 | 0.198856 |
| 15078 | 0.639979 | 0.512312 | -0.127668 |
| 16710 | 0.925497 | 0.558168 | -0.36733 |
| 28045 | 0.369218 | 0.557603 | 0.188385 |

Point values (pooled over the training pairs and per training donor; no bootstrap, no verdict) of D_R2 - D_R1, (b') term_flagged_R2_minus_R1 and (c') Dlog_R2_minus_Dlog_R1, computed as E3 evaluate does, for two constructed readout pairs with R1 = the head: (uniform_shrinkage) R2 = the head's NK-T evidence differences times s, which repairs nothing; (flagged_repair) R2 = the measured differences on flagged pairs and the head's on unflagged pairs; head and null set to R1 (they cancel). s = 0.5.

| scenario | statistic | pooled | 10886 | 11466 | 12710 | 15078 | 16710 | 28045 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| uniform_shrinkage | D_R2_minus_D_R1 | 0.065798 | None | 0.049629 | 0.099428 | -0.063834 | -0.183665 | 0.094193 |
| uniform_shrinkage | term_flagged_R2_minus_R1 | -0.207842 | None | -0.175207 | -0.363461 | -0.31999 | -0.462749 | -0.184609 |
| uniform_shrinkage | Dlog_R2_minus_Dlog_R1 | 0.0 | None | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| flagged_repair | D_R2_minus_D_R1 | 0.584317 | None | 0.649586 | 0.273077 | 0.360021 | 0.074503 | 0.630782 |
| flagged_repair | term_flagged_R2_minus_R1 | 0.584317 | None | 0.649586 | 0.273077 | 0.360021 | 0.074503 | 0.630782 |
| flagged_repair | Dlog_R2_minus_Dlog_R1 | 0.877831 | None | 1.048641 | 0.318935 | 0.446319 | 0.077424 | 0.996368 |

Uniform shrinkage: |D_R2 - D_R1 - (1 - s)(GR_head(unfl) - GR_head(fl))| = 0.0 (the identity above).

## Leakage check

As A2 (v3_build_amendment_A2.py) and A1 (v3_leakage_check.py): rebuild the core from copies of the inputs whose rows of one split carry random protein, cell types and embedding, and compare bytes.

| poisoned rows | role | core identical to real |
|---|---|---|
| test (seed 101) | site4 rows poisoned: core must be byte-identical | True |
| val (seed 303) | val rows poisoned: A3 reads no val value, so identical as well (information) | True |
| train (seed 202) | positive control, training rows poisoned: core must change | False |

Result: **PASS**; real core sha256 `39de4a238f260cfd214e07f81c90015c9b34c96ab0d74dd908f75035c6da5476`.

## Experiment fields replaced (path in registration_v3.json → experiments, after A1 and A2)

- `E3/endpoints/E3.H3a_pooling`: D_R2 - D_R1 >= 0.05: the loss is in mean pooling (registered; margin 0.05 with the registration win / loss / equivalent rules). A3.1: the pooling claim wins only if the registered rule wins and, under the same two-stage donor bootstrap and per-donor rule, (b') recovery on flagged pairs beyond R1: term_flagged_R2_minus_R1 = GR_R2(flagged) - GR_R1(flagged) > 0, its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor and (c') the scale-free contrast Dlog_R2 - Dlog_R1 = [ln GR_R2(fl) - ln GR_R1(fl)] - [ln GR_R2(unfl) - ln GR_R1(unfl)] (Dlog_R as the E3 builder computes it, addendum secondary_H3a S1), > 0 with its two-stage 95% interval excluding 0 (lower bound > 0) and > 0 in each primary donor; a registered pooling win without (b') and (c') is 'not supported (scale artefact)'

## Addendum and amendment fields superseded (files unchanged)

- `registration/addenda/E3.json` (sha256 `09f932ed166cad02...`): declared.H3a: 'pooling = D_R2 - D_R1; margin 0.05 with the registration win / loss / equivalent rules' -> the A3.1 pooling win (registered win plus (b') and (c')); declared.secondary_H3a: S1 (Dlog_R) enters the A3.1 pooling condition (c') as Dlog_R2 - Dlog_R1
- `registration/amendment_A2.json` (sha256 `20b64174b063d399...`): findings A2.1 fix: 'H3a_pooling (D_R2 - D_R1) is unchanged' -> A3.1

## Change log

| file | change | reason |
|---|---|---|
| `registration/amendment_A3.json (+ .sha256, amendment_A3_core.json, AMENDMENT_A3.md)` | new amendment A3 | A3.1 above; written before any site4 evaluation of E3 |
| `bridge_anm/lib/v3_amend.py` | load_registration_amended() applies A1, A2, then A3, verifying the registration, A1, A2 and A3 hashes and that A3 names the registration, A1 and A2 on disk; amendment_A3_status() for the E3 site4 evaluate guard; a3_pooling_verdict() | every v3 builder reads the amended registration through this loader |
| `scripts/v3_e3_nkt_repair.py` | evaluate on site4 refuses unless A3 is committed and matches; the E3.H3a_pooling verdict is the A3.1 verdict, computed on unrounded values and the same replicates (the b4ecad4 secondary block becomes the deciding A3 block; the registered pooling verdict is reported beside it); DECLARED (the addendum's declared part) is unchanged | A3.1; earlier stages (null, r1, r2_states, r2_fit, r2_predict) are unchanged |
| `tests/test_v3_amend.py, tests/test_e2_response.py` | A3 loading, tamper / chain / missing-file detection, status, the E3 evaluate guard, A3.1 end to end through run_h3a / verdicts; the E2 guard test copies A3 (the loader needs it) | A3 |
| `registration/HASHES.txt` | A3 lines; v3_amend.py line moved to the A3 version (the A2 version's sha256 is in amendment_A2.json provenance) | hash record of the registration files |

## Provenance

Builder `v3_build_amendment_A3 1.0` sha256 `59d7cf1478109ce1a5545f3b66a240412c7fa48bb114f0b4cedc1473665ef4ec`; `v3_amend.py` sha256 `2655ed06ffb6e182d57c52f73d2a782e833ee17540a11762a4aa194a8b1dda91`; `v3_e3.py` sha256 `a1a3bfa08ca41a664e8134d425c011e4876f27ac8ad04ed54a6ee12297d9bede`; computed core sha256 `39de4a238f260cfd214e07f81c90015c9b34c96ab0d74dd908f75035c6da5476`; built 2026-10-01T07:00:56Z.
