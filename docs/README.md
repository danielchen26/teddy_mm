# TEDDY × ANM — GitHub Pages site

Open **[index.html](index.html)** for the public landing: a “Corrections” notice (rerun done for Experiments 1, 3, 4 and 5; Experiments 2 and 6 withdrawn) → setup (how ANM connects to TEDDY, what we compare, Mode A vs B) → **registered experiments (v3)**: E1–E6 and E5-M result cards (all final), fairness timeline and claims ledger, numbers from `reports/v3_numbers.json` → problem → comparison protocol → six experiments (each with method, interactive chart, data table and what it proves) → claim boundary → reproduce.

- `anm-loop.html` — interactive report for Experiment 6, look-alike cells (veto → complement → verify, scoreboard by method × metric; withdrawn pending redesign, first-run numbers kept as a record)
- `assets/site/guide.js`, `assets/site/guide.css` — shared by both pages: one plain-language definition per term (dotted underline, hover/focus/tap), auto-marking of abbreviations such as O0 or z_512 and of the method names (TEDDY + fixed rule, TEDDY + trained classifier), and the guide column. The column follows the section you are reading: where it sits in the workflow, “how ANM decides here” (set per section in each page’s `GUIDE_CONFIG.sections[id].how`), and “Inside ANM’s engine”, one real cell run step by step (rebuild and check it with `scripts/guide_engine_trace.py`). On wide screens the column floats on the right, can be hidden to a tab on the right edge, and can be resized by dragging its left edge (arrow keys work too); both choices are remembered per browser.
- `assets/infographics/` — icon-first story PNGs (light+dark), still showing first-run numbers; regenerate via `scripts/make_story_infographics.py`
- `assets/hard_proof/`, `assets/missing_modality/` — proof GIFs/dashboards
- `reports/V3_REGISTERED_RESULTS.md`, `reports/v3_numbers.json` — the registered v3 results and claims ledger, and every v3 number with its source file and key; they supersede the v2 reports below where they conflict
- `reports/` — HARD_PROOF_V2 and SCOPE_REFINE_PROOF_V2 (corrected rerun, copied from `outputs/anm_cite_bridge_v2/`); first run: HARD_PROOF, SCOPE_REFINE, MISSING_MODALITY, STAT_PROOF, BAKEOFF, v0 REPORT

Pages source: `/docs` on the default branch.
- `MODE_B_SCOPE.md` — Mode B vs Mode A claim boundary + `z_512` definition
- `reports/REVERSE_STEP1.md` — must-separate z pairs (sufficiency-for-readout; not Mode A)
- `reports/ANM_HELPS_TEDDY_LOOP.md` — deep proof: how Mode B veto→complement→verify helps TEDDY (figures `06`/`07`/`08`)
- `reports/REVERSE_LOOP_SMALL.md` — veto/complement/verify stats on must-pairs
