# TEDDY × ANM — GitHub Pages site

Pages source: `/docs` on the default branch. The three current pages are written by `python scripts/build_site.py`, which
reads every number from records committed under `anm_diag/runs/`; do not edit them by hand. Each page carries English and
Chinese text (English first).

- **[index.html](index.html)**: overview. Where the PD-L1 response to a knockout is lost (the sealed tests), the loop so far
  in one table, what the repository holds, and what can and cannot be claimed.
- **[anm-framework.html](anm-framework.html)**: how the audit works. Roles, the statistic *u* and its verdicts, the checks that
  refuse illegal designs, typed source events, the loop, evidence grades, formal definitions.
- **[anm-loop.html](anm-loop.html)**: the loop on the teddy_mm head, round by round: the repair at the located step, the
  control that should not work, the cost on the original task, what is left.
- `_src/`: the two explanatory figures and the definitions used by the framework page (not served).

Archived, as they were in October 2026 (the registered decision-layer tests), with a notice at the top:

- `v3.html` (was `index.html`): “Corrections” notice → setup → registered experiments (v3) E1–E7 → the six first-run experiments → claim boundary.
- `v3-anm-framework.html` (was `anm-framework.html`): the framework page of that time (physics, AI workflows, TEDDY).
- `v3-anm-loop.html` (was `anm-loop.html`): interactive report for Experiment 6, look-alike cells (veto → complement → verify; withdrawn, first-run numbers kept as a record; redesigned and registered as v3 E3, E5 and E5-M).
- `assets/site/guide.js`, `assets/site/guide.css` — shared by the archived `v3.html` and `v3-anm-loop.html`: one plain-language definition per term (dotted underline, hover/focus/tap), auto-marking of abbreviations such as O0 or z_512 and of the method names (TEDDY + fixed rule, TEDDY + trained classifier), and the guide column. The column follows the section you are reading: where it sits in the workflow, “how ANM decides here” (set per section in each page’s `GUIDE_CONFIG.sections[id].how`), and “Inside ANM’s engine”, one real cell run step by step (rebuild and check it with `scripts/guide_engine_trace.py`). On wide screens the column floats on the right, can be hidden to a tab on the right edge, and can be resized by dragging its left edge (arrow keys work too); both choices are remembered per browser.
- `assets/infographics/` — icon-first story PNGs (light+dark), still showing first-run numbers; regenerate via `scripts/make_story_infographics.py`
- `assets/hard_proof/`, `assets/missing_modality/` — proof GIFs/dashboards
- `reports/V3_REGISTERED_RESULTS.md`, `reports/v3_numbers.json` — the registered v3 results and claims ledger, and every v3 number with its source file and key; they supersede the v2 reports below where they conflict
- `reports/` — HARD_PROOF_V2 and SCOPE_REFINE_PROOF_V2 (corrected rerun, copied from `outputs/anm_cite_bridge_v2/`); first run: HARD_PROOF, SCOPE_REFINE, MISSING_MODALITY, STAT_PROOF, BAKEOFF, v0 REPORT

- `MODE_B_SCOPE.md` — Mode B vs Mode A claim boundary + `z_512` definition
- `ANM_DIAGNOSIS.md` — ANM response-sufficiency diagnosis of TEDDY (KS-2 exploratory; Perturb-CITE-seq and ECCITE-seq sealed: TEDDY's stages never rejected, loss in the teddy_mm head after the head-input fix); code and runs in `../anm_diag/`
- `reports/PHASE1_PER_PROTEIN.md` — phase-1 per-protein test Pearson / R² (median 0.675 vs mean 0.603) and known limitations
- `reports/REVERSE_STEP1.md` — must-separate z pairs (sufficiency-for-readout; not Mode A)
- `reports/ANM_HELPS_TEDDY_LOOP.md` — deep proof: how Mode B veto→complement→verify helps TEDDY (figures `06`/`07`/`08`)
- `reports/REVERSE_LOOP_SMALL.md` — veto/complement/verify stats on must-pairs
