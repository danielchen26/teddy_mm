#!/usr/bin/env python3
"""E7 (Mode A): TEDDY analog of ANM's three-step inverse loop, with the consumer fixed.

ANM v2 (PR #16) closes an inverse loop on a model consumer: (1) reject a retained-state description with a
witness pair (equal state, different readout), (2) select a revised description from a DECLARED candidate library
by a DECLARED rule on development data, (3) validate it on fresh comparisons with the consumer and interface
unchanged. E7 does the same on TEDDY at the layer-11 cut:

    consumer (fixed, never trained)  layer-11 token states H -> TEDDY layer 12 (E5's validated explicit-attention
                                     layer_fn) -> gene-mean -> L2 normalisation -> frozen phase-1 head
    observable                       the head's CD56, CD94, CD335, CD3 predictions, each in units of the measured
                                     val NK-T gap of that protein (val donor 18303 only); difference = max over the
                                     4 proteins of |change| / unit
    candidate library (ordered)      (1) layer-11 gene-mean; (2) gene-mean + the layer-11 states of the tokens of G =
                                     {NCAM1, KLRD1, NCR1, CD3E, CD3D, CD3G} (the coding genes of the 4 panel proteins);
                                     (3) all layer-11 token states (no matched pair: uninformative)
    histories                        exact mean-preserving state interventions at the cut: m disjoint token pairs,
                                     each pair's two states replaced by their average (never swaps: layer 12 has no
                                     positional input, so a swap is invisible to the consumer)
    matched / witness pairs          (unpatched, patched) with equal candidate state; witness = equal gene-mean,
                                     different G states (a patch that averages a G token)
    development                      val donor 18303 only: magnitude from the positive control alone, then the
                                     declared selection rule; none -> stop, keep the record, one declared follow-up
    confirmation                     fresh site4 test_primary NK/T cells (none of E3's pairs, E5's or E5-M's cells),
                                     fresh pairs, a fresh magnitude; secondary family on E6's external donors

Stages (each resumable and logged to <out-dir>/progress.log)
  register  CPU, train/val only: G, units, rosters, architecture facts, leakage check (site4 rows poisoned ->
            identical core; val rows poisoned -> different core); writes registration/addenda/E7.json + HASHES line
  develop   device: the development roster (--roster dev) or the follow-up roster (--roster followup), one file per
            cell; refused unless the addendum is committed (the follow-up also needs the committed null selection)
  select    CPU: the declared selection rule, precision and implementation gates on development cells;
            writes registration/addenda/E7_selection.json + HASHES line (to be committed before any site4 forward)
  confirm   device: --family site4 (registered) or external (secondary); refused unless the addendum AND the final
            selection record are committed and match
  report    CPU: bounds, statuses, verdict; E7_results.json + REPORT.md
Smoke runs (--smoke NOTE) use --cell-pool smoke (at most 5 val cells of the declared smoke pool, outside both
development rosters) or --cell-pool train (training cells, code tests only); they never write the registration
directory and never decide anything.

Full run (from the repo root; R = the main checkout; not run by the registering session):
  PY=<venv>/bin/python; OUT=$R/outputs/v3/E7
  $PY scripts/mode_a_inverse_e7.py --stage register                      # then commit E7.json + HASHES.txt
  PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage develop --roster dev --out-dir $OUT
  $PY scripts/mode_a_inverse_e7.py --stage select --roster dev --out-dir $OUT  # commit E7_selection.json + HASHES
  (only if the record says follow_up_required: develop --roster followup, select --roster followup, commit again)
  PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage confirm --family site4 --out-dir $OUT
  PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_inverse_e7.py --stage confirm --family external --out-dir $OUT
  $PY scripts/mode_a_inverse_e7.py --stage report --out-dir $OUT
"""
from __future__ import annotations

import os
import sys


def _early_threads(argv) -> int:
    n = 4
    for i, a in enumerate(argv):
        if a == "--threads" and i + 1 < len(argv):
            n = int(argv[i + 1])
        elif a.startswith("--threads="):
            n = int(a.split("=", 1)[1])
    return max(1, n)


_THREADS = _early_threads(sys.argv[1:])
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = str(_THREADS)

import argparse  # noqa: E402
import copy  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import inspect  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "bridge_anm"))
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location("mode_a_layer_dynamics", HERE / "mode_a_layer_dynamics.py")
mld = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mld)
lp, vk = mld.lp, mld.vk

from lib import v3_amend as va  # noqa: E402
from lib import v3_e7 as v7  # noqa: E402

SCRIPT_VERSION = "mode_a_inverse_e7 1.0"
MAIN = Path("/Users/tianchichen/Documents/GitHub/teddy_mm")
_BASE = MAIN if MAIN.exists() else ROOT
ADDENDUM, SELECTION = "E7.json", "E7_selection.json"
PROTEINS = ("CD56", "CD94", "CD335", "CD3")  # E3's gap proteins = E5's / E5-M's head proteins: the registered NK/T panel
assert tuple(mld.PROTEINS) == PROTEINS
G_GENES = {"NCAM1": "CD56", "KLRD1": "CD94", "NCR1": "CD335", "CD3E": "CD3", "CD3D": "CD3", "CD3G": "CD3"}
SEED, E7_TAG = 20260930, 7  # registration seeds.global, E7's tag
PHASE_TAG = {"dev": 1, "followup": 2, "site4": 3, "external": 4, "train": 8, "smoke": 9, "train_family": 10,
             "smoke_family": 11}
TOL = 0.05
LADDER = (1 / 32, 1 / 8, 1 / 2)  # development magnitudes: share of the cell's non-G tokens averaged in pairs
FRESH = (1 / 4,)  # the confirmation's fresh magnitude (never used in development)
KAPPA = 1.0  # clamp positive control: +- KAPPA x the training NK - T difference of the layer-11 gene-mean
N_MATCHED = 8  # matched patches per magnitude; each has one G twin
N_G_ONLY = 8  # G-only patches per cell (magnitude-free G witnesses)
N_DEV_PER_CLASS = N_F1_PER_CLASS = 100
N_SITE4_PER_CLASS_PER_DONOR = 100
N_EXT_PER_CLASS_PER_DONOR = 25
MIN_NONG_TOKENS = 256
MIN_G_PRESENT = 1
N_F64 = 10
F64_TOL = TOL / 10
N_PUSH_CELLS = 50
N_PAD = 16
UNIT_MIN = 0.05
N_BOOT = 2000
GATES = {"padding_patch_gap_units": 1e-6, "layer12_mean_preserving_gap_units": 1e-5, "permutation_gap_units": 1e-5,
         "mean_preservation_abs": 1e-5, "module_vs_explicit_abs": 1e-4, "padded_vs_unpadded_abs": 1e-4,
         "official_layer11_gene_mean_abs": 1e-3}
N_SMOKE_MAX = 5
EXIT_INCOMPLETE = 75
_LOG = {"path": None}


def say(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG["path"] is not None:
        with open(_LOG["path"], "a") as f:
            f.write(line + "\n")


DECLARED = {
    "question": "TEDDY analog of ANM v2's inverse loop with the consumer fixed: is TEDDY's layer-11 gene-mean a "
                "sufficient retained state for the fixed consumer (layer 12 + gene-mean pooling + L2 normalisation + "
                "the frozen phase-1 head) on the registered NK/T panel, under a declared class of exact mean-preserving "
                "interventions at the layer-11 cut? If it is rejected, does a declared augmentation chosen by a declared "
                "rule on development data get bounded support on fresh comparisons?",
    "anm_template": "ANM v2 discussion: 'Exactly equal retained states give a readout difference beyond measurement "
                    "error -> reject sufficiency -> separate the violating histories in a revised representation; test it "
                    "with consumer and interface unchanged.' Selection operates within a supplied candidate library; "
                    "neither automatic discovery of a missing variable nor a unique state encoding follows.",
    "architecture_facts": {
        "positional_input": "TEDDY-G adds a learned position embedding once, to the input (teddy_mm/teddy_encoder.py "
                            "hidden_states: embeddings(ids) + position_embeddings(arange(L))); each of the 12 encoder layers "
                            "is a plain post-norm torch.nn.TransformerEncoderLayer (d 512, 8 heads, d_hid 1024, GELU, no "
                            "positional or relative term). Layer 12 is therefore equivariant to permutations of the token "
                            "positions, gene-mean pooling is invariant, and the consumer is a symmetric function of the "
                            "multiset of layer-11 token states (checked numerically at register: architecture_check)",
        "consequences": "(a) a swap patch (exchange two tokens' states) is invisible to the consumer, so E7 uses no swaps; "
                        "(b) gene identity reaches the consumer only through the states themselves, so the G-token states "
                        "of candidate 2 are an analyst's labelling of which states are retained; (c) the full set of "
                        "layer-11 states (candidate 3) is sufficient by construction but has no matched pair among distinct "
                        "histories, so it is uninformative, as the full history in ANM's study",
        "padding": "the official preprocessing pads to 2,048 with <pad> and a bool key-padding mask (teddy_mm passes the "
                   "bool mask; the official model.py's float mask would be added to the logits on MPS, so padding would be "
                   "attended there). E7 runs every cell unpadded (its own L real tokens, no mask), as E5 and E5-M; the "
                   "padded bool-masked module forward equals it on the real tokens (control), and a patched padding "
                   "position changes nothing (control)",
        "dropout": "eval mode: dropout (p 0.02 in the config) inactive",
        "cut": "layer-11 cut = output of encoder.layers[10]; the consumer starts at encoder.layers[11]; the layer-11 "
               "gene-mean is the official z_rna_layer_means[10] (fp16-autocast official run; the explicit float32 forward "
               "matches it to the declared gate)",
    },
    "consumer": "fixed, never trained: E5's validated explicit-attention layer_fn with encoder.layers[11]'s own weights "
                "(forward equal to the module's to ~5e-6), mean over the cell's real tokens, then E5's load_head: "
                "z / (||z|| + 1e-6) -> MLP -> NB decoder mean at the training-median size factor, / training p95 of measured "
                "ADT (the frozen phase-1 head outputs/cite_phase1_official/best.pt, sha256 in inputs)",
    "observable": {
        "panel": "O = the head's predictions of CD56, CD94, CD335, CD3 (E3's gap proteins = E5's and E5-M's head proteins); "
                 "difference between two histories = max over the proteins of |O_p - O'_p| / u_p",
        "units": "u_p = |median over val primary-key NK cells - median over val primary-key T cells| of the measured "
                 "m_p = ADT_p / p95_p (p95 = the head's normaliser: lp.p95_replicate, training cells), val donor 18303 "
                 "only, computed at register and frozen here; a protein with u_p < 0.05 would be dropped (recorded)",
        "score_for_follow_up": "s(O) = mean(O_CD56, O_CD94, O_CD335) - O_CD3 (E5's head-score lens), unit u_s = |median "
                               "val NK - median val T| of s(m); used only by the declared follow-up",
        "secondary_units": "the head's own predicted val NK-T gap per protein (official z), reported as a re-expression; "
                           "never decides",
    },
    "tolerance": {"tol": TOL, "units": "per-protein measured val NK-T gap (panel) or u_s (score)",
                  "why": "a change of 5% of the measured NK-T difference of a protein is negligible for NK-vs-T reading; "
                         "fixed before development and before any patch output on a val or site4 cell"},
    "candidates": {
        "order": list(v7.CANDIDATES),
        "gene_mean": "phi_1(H) = mean_t h_t (dimension 512)",
        "gene_mean_plus_G": "phi_2(H) = (mean_t h_t, (h_t) for the tokens of G present in the cell) (dimension 512 (1 + g), "
                            "g = number of G genes present, 1..6)",
        "all_token_states": "phi_3(H) = (h_t)_t (dimension 512 L); no two distinct histories share it -> no matched pair "
                            "-> uninformative by construction",
        "G": {"genes": dict(G_GENES),
              "why": "fixed from the observable alone (it was the orchestrator's example), before any computation: the "
                     "coding genes of the 4 panel proteins (CD56 NCAM1, CD94 KLRD1, CD335 NCR1, CD3 CD3E/CD3D/CD3G; "
                     "Ensembl ids of E5's TOKEN_PROBE_IDS, all in TEDDY's vocabulary and the 12,052-gene universe). These "
                     "are the tokens most directly tied to the readout, as the sensor-recipient association is in ANM's "
                     "study; at most 6 tokens (~0.4% of a 1,400-token cell), so candidate 2 is a minimal augmentation. "
                     "E5-M's other NK genes (FCGR3A, KLRF1) code proteins outside the panel and are not in G. No "
                     "data-driven choice"},
    },
    "histories": {
        "class": "the declared history class = the unpatched cell and exact mean-preserving pair-averaging patches at the "
                 "layer-11 cut at the declared magnitudes: disjoint pairs of real-token positions, each pair's two states "
                 "replaced by their average (exact mean preservation in exact arithmetic; float32 rounding checked by the "
                 "gate mean_preservation_abs)",
        "magnitude": "f = share of the cell's non-G tokens averaged: m = max(1, floor(f (n_nonG - g) / 2)) pairs "
                     "(v3_e7.pairs_for_fraction). Development: f in {1/32, 1/8, 1/2}, all three (no adaptive magnitude). "
                     "Confirmation: the same three with fresh pairs plus the fresh magnitude 1/4",
        "families_per_magnitude": f"{N_MATCHED} matched patches (m pairs among non-G tokens) and, for each, its G twin (the "
                                  "same pairs plus every G token averaged with its own non-G partner outside the patch: same "
                                  "gene-mean, only the G states differ from the matched patch)",
        "g_only": f"{N_G_ONLY} G-only patches per cell (every G token averaged with its own random non-G partner, nothing "
                  "else): magnitude-free G witnesses, counted at every magnitude",
        "draws": "numpy default_rng([20260930, 7, phase tag, global cell id, 1024 f]) per cell and magnitude (matched_i then "
                 "its twin's partners, i = 1..8: v3_e7.draw_patch_sets) and default_rng([20260930, 7, phase tag, cell, 1]) "
                 "for the G-only patches (v3_e7.draw_g_only); phase tags dev 1, followup 2, site4 3, external 4 (fresh "
                 "pairs at confirmation); uniform positions",
        "matched_pairs": "gene_mean: (unpatched, matched_i), (unpatched, twin_i), (matched_i, twin_i), (unpatched, "
                         "G-only_j); gene_mean_plus_G: (unpatched, matched_i); all_token_states: none",
        "witness_pairs": "(matched_i, twin_i) and (unpatched, G-only_j): equal gene-mean, only the G states differ: matched "
                         "under candidate 1, not under candidate 2; G-specific at every magnitude (the context patch is "
                         "the same on both sides)",
        "no_swaps": "a swap leaves the multiset of states unchanged, so the consumer's output is unchanged (architecture "
                    "fact; checked)",
    },
    "eligibility": f"a cell enters any roster only if at least {MIN_G_PRESENT} G gene is among its tokens and it has >= "
                   f"{MIN_NONG_TOKENS} non-G tokens (at least 3 pairs at the smallest magnitude); computed from the official "
                   "tokenisation (label-free)",
    "per_cell_statistics": "d = observable difference between the two histories of a pair; D_gene_mean = max d over "
                           "candidate 1's matched pairs, D_gene_mean_plus_G = max d over candidate 2's, W = max d over the "
                           "G witness pairs (W_context and W_pure reported), P = max d over (unpatched, clamp +-) "
                           "(v3_e7.cell_summary), over every magnitude of the stage (G-only and clamp patches always count)",
    "development": {
        "roster": f"val donor 18303 primary-key NK and T cells, eligible: {N_DEV_PER_CLASS} NK + {N_DEV_PER_CLASS} T drawn "
                  "with numpy default_rng([20260930, 7, 0]) (NK, then T), listed in rosters.dev with the processing order",
        "magnitudes": "1/32, 1/8, 1/2 (all; the history class is their union)",
        "selection_rule": "as ANM's study: a candidate is retained when it has >= 1 matched pair and its matched-pair "
                          "differences on the development roster are <= tol, summarised as the 95th percentile over cells "
                          "of the per-cell maximum over the whole class (the point value of the statistic the confirmation "
                          "bounds; the literal max over every cell is recorded beside it and never decides); the minimal "
                          "retained candidate (dimension order) is selected, ties by declared order; none retained -> stop, "
                          "keep the record, and run the single declared follow-up",
        "gates": "implementation checks on every development cell and the float64 re-check on the first 10 cells of the "
                 "order must pass for a selection to count; otherwise the record says stopped (precision / implementation) "
                 "and no follow-up runs",
        "follow_up": f"only if the development selection is null (not stopped), and only once, labelled 'follow-up after a "
                     f"null primary selection': the observable restricted to the NK-T score s (unit u_s), the same library, "
                     f"history class and rule, on the fresh roster rosters.followup ({N_F1_PER_CLASS} NK + {N_F1_PER_CLASS} "
                     "T val donor 18303 cells disjoint from rosters.dev, fixed here; fewer NK if fewer remain). Its develop "
                     "stage is refused until the null primary record is committed",
        "record": "registration/addenda/E7_selection.json (sha256 line in addenda/HASHES.txt), committed before any site4 "
                  "or external forward of this script; it pins the development digest (sha256 of every development output)",
    },
    "confirmation": {
        "site4": f"test_primary donors 13272 and 19593, primary-key NK and T cells, eligible, excluding every cell of "
                 f"E5's e5_cells.npz (E3's 278 NK-T pair cells and E5's 200 random cells) and of E5-M's e5m_cells.npz "
                 f"(a subset; both files pinned by sha256); per donor and class up to {N_SITE4_PER_CLASS_PER_DONOR} cells "
                 "drawn with numpy default_rng([20260930, 7, 3]) (donors in registration order, NK then T); all eligible "
                 "cells if fewer",
        "magnitudes": "1/32, 1/8, 1/2 with fresh pairs (phase tag 3 or 4) and the fresh magnitude 1/4 (never used in "
                      "development); per-cell statistics over all four",
        "observable": "the one the selection was made on: panel if the primary selection decided (or nothing was "
                      "selected), the NK-T score if the follow-up decided; the other is reported as secondary",
        "external_secondary": f"E6's 8 external donors (Hao 2021, GSE164378, E6's d5k pack; annotation-only key with E6's "
                              f"committed map, because E6's key was not validated): per donor and class up to "
                              f"{N_EXT_PER_CLASS_PER_DONOR} eligible NK and T cells, numpy default_rng([20260930, 7, 4]); "
                              "labelled 'RNA possibly seen by TEDDY in pretraining' and 'key not validated'; pooled bounds "
                              "over the 8 donors (per-donor values reported); never changes the registered verdict",
    },
    "statistics": {
        "bootstrap": f"registration common.statistics: two-stage (donors with replacement, then cells with replacement "
                     f"within each drawn donor and class), B = {N_BOOT}, seed = seeds.bootstrap (1), percentile 95% "
                     "intervals; per donor: a cell bootstrap within the donor and class, same B and seed",
        "support": "candidate k has bounded support if the upper end of the 95% interval of the 95th percentile over cells "
                   "of D_k is <= tol, pooled and in each primary donor",
        "rejection": "candidate k is rejected if the lower end of the 95% interval of the median over cells of D_k is > "
                     "tol, pooled and in each primary donor; otherwise unresolved. all_token_states: uninformative",
        "witness_separation": "the G witness is detected if the lower end of the 95% interval of the median of W is > tol "
                              "(pooled and each primary donor)",
        "positive_control": "detected if the lower end of the 95% interval of the median of P is > tol (pooled and each "
                            "primary donor): the consumer resolves a change of the retained gene-mean of this size; failure "
                            "makes any support UNINFORMATIVE, never support",
        "detectability_thresholds": "positive control and G witness: tol (the instrument must resolve a difference of "
                                    "the tolerance size); numerical: the implementation gates and the float64 rule below",
        "two_donor_caveat": "2 primary donors: the intervals describe these donors, not a population",
        "verdicts": dict(v7.VERDICTS),
    },
    "controls": {
        "positive": {
            "clamp": f"registered: two clamp patches per cell, every layer-11 token moved by +{KAPPA} and -{KAPPA} x "
                     "dmu_11 = E5's training NK - T difference of the layer-11 gene-mean (outputs/v3/E5/e5_directions.npz "
                     "dmu[11], training cells, primary key; pinned by sha256 and copied here): the gene-mean moves by the "
                     "full NK - T difference (an NK/T-sized change of the retained state, as ANM's queried-sensor +1), the "
                     "deviations do not; P = max of the two differences",
            "G_witness": "registered: W (the revision is shown necessary only if W is detected and the gene-mean rejected)",
            "input_push": f"secondary and descriptive (the review's E5-M push): each G gene's token embedding doubled at the "
                          f"input, full forward, first {N_PUSH_CELLS} cells of each confirmation order; it acts at layer 0, "
                          "so it moves the cut's gene-mean and deviations together and does not decide",
        },
        "invariance_implementation": {
            "padding_patch": f"layer 12 run by the module on the cell padded with {N_PAD} extra positions (bool "
                             "key-padding mask, CPU float32); replacing one padding state must change the head by <= "
                             f"{GATES['padding_patch_gap_units']} gap units (observed exactly 0 at register)",
            "layer12_mean_preserving_patch": "averaging two real tokens' layer-12 states leaves the gene-mean, hence the "
                                             f"head, unchanged: <= {GATES['layer12_mean_preserving_gap_units']} gap units",
            "permutation": f"permuting the layer-11 rows changes the head by <= {GATES['permutation_gap_units']} gap units "
                           "(no positional input in layer 12)",
            "mean_preservation": f"max |gene-mean change| of every pair-averaging patch <= {GATES['mean_preservation_abs']}",
            "module_vs_explicit": f"explicit layer_fn vs the module's layer 12 (unpadded) and the padded masked module on "
                                  f"the real tokens: max |state difference| <= {GATES['module_vs_explicit_abs']}",
            "official_state": f"explicit float32 layer-11 gene-mean vs the official stored one (fp16 autocast): <= "
                              f"{GATES['official_layer11_gene_mean_abs']}",
            "rule": "any failure on any cell of a stage -> NOT_VALIDATED (no verdict)",
        },
    },
    "float64_recheck": f"the first {N_F64} cells of every roster's order are recomputed in float64 on the CPU (layers 1-12 "
                       "and the head, the same histories); if any history's difference differs from the float32 one by "
                       f"more than tol / 10 = {F64_TOL}, the stage is UNRESOLVED_PRECISION",
    "can_show": "for one fixed consumer, one observable panel and one declared class of mean-preserving state "
                "interventions at the layer-11 cut: whether the layer-11 gene-mean is rejected as a sufficient retained "
                "state on fresh cells and donors, and whether a declared marker-token augmentation selected on development "
                "data gets bounded support there (if so, the first complete reject -> select -> validate step on TEDDY)",
    "cannot_show": ["that RNA lacks NK-T information", "that the head improves (the consumer is unchanged; E3 rejected "
                    "repair on frozen TEDDY)", "that any representation is unique or minimal beyond the declared library",
                    "transfer to other consumers, layers, observables or intervention classes", "behaviour under on-manifold "
                    "(real-cell) perturbations: these are state interventions, not input histories (ANM v2: 'Propagation "
                    "imposes the state displacements; the readout alone is measured afresh')", "clinical value"],
    "fixed_before_development": "every choice and number here was fixed on training and validation cells and committed "
                                "before any development forward of this script; no patch output on a val or site4 cell "
                                "informed any choice (patch outputs were computed on training cells only: the 3 cells of "
                                "architecture_check and the 7 of pre_registration_code_test)",
    "pre_registration_code_test": "before this addendum, a first code run used 2 training cells without G tokens (code "
                                  "paths only), then the code ran on 2 annotated training cells (1 NK, 1 T; no val or "
                                  "site4 cell) with an earlier draft that used one-sided twins as the positive control and "
                                  "a magnitude calibrated on them. Median per-cell maxima (panel, gap units): gene-mean "
                                  "0.014 / 0.065 / 0.166 and gene-mean + G 0.013 / 0.055 / 0.149 at f = 1/32 / 1/8 / 1/2; "
                                  "G witness 0.007 / 0.011 / 0.027; one-sided twins 0.002 / 0.007 / 0.026, so a "
                                  "twin-calibrated magnitude (target 2 tol) would never have been reached: the twin moves "
                                  "the gene-mean by a random-walk sum. The draft was revised to the clamp positive control "
                                  "and the fixed three-magnitude class, and rerun on 3 annotated training cells (2 NK, 1 T) "
                                  "with the clamp at kappa 0.2: median P 0.049 (at the tolerance), median W 0.019, "
                                  "q95 of the per-cell maximum 0.227 (gene-mean) and 0.202 (gene-mean + G). kappa was then "
                                  "set to 1 (the full training NK - T difference, the interpretable NK/T-sized change; 0.2 "
                                  "had no rationale beyond being small); the same 3 cells were rerun on MPS to check the "
                                  "device path (float32 vs float64 within 1e-6). Nothing else changed after these code tests",
}


# ============================================================================ args
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("register", "develop", "select", "confirm", "report"), required=True)
    p.add_argument("--roster", choices=("dev", "followup"), default="dev")
    p.add_argument("--family", choices=("site4", "external"), default="site4")
    p.add_argument("--processed", type=Path, default=_BASE / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=_BASE / "data/processed/cite_official")
    p.add_argument("--external-pack", type=Path, default=_BASE / "data/processed/external_hao2021_d5k")
    p.add_argument("--external-embed-dir", type=Path, default=_BASE / "data/processed/external_hao2021_d5k_official")
    p.add_argument("--ckpt", type=Path, default=_BASE.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=_BASE / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=_BASE / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--e5-dir", type=Path, default=_BASE / "outputs/v3/E5")
    p.add_argument("--e5m-dir", type=Path, default=_BASE / "outputs/v3/E5M")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--out-dir", type=Path, default=_BASE / "outputs/v3/E7")
    p.add_argument("--device", default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--batch", type=int, default=4, help="histories per layer-12 call")
    p.add_argument("--cell-pool", choices=("registered", "smoke", "train"), default="registered")
    p.add_argument("--n-smoke", type=int, default=N_SMOKE_MAX)
    p.add_argument("--n-boot", type=int, default=N_BOOT)
    p.add_argument("--max-minutes", type=float, default=0, help="develop/confirm: pause after N minutes (exit 75; rerun resumes)")
    p.add_argument("--smoke", default=None, metavar="NOTE")
    a = p.parse_args(argv)
    if (a.cell_pool != "registered" or a.n_boot != N_BOOT) and not a.smoke:
        p.error("--cell-pool smoke/train and --n-boot are for smoke runs only (give --smoke NOTE)")
    if a.smoke and a.cell_pool == "registered" and a.stage != "register":
        p.error("a smoke run uses --cell-pool smoke (val smoke pool) or train; registered cells are never smoked")
    if a.n_smoke > N_SMOKE_MAX:
        p.error(f"--n-smoke is at most {N_SMOKE_MAX}")
    return a


# ============================================================================ helpers
def sha256_file(path: Path) -> str:
    return vk.sha256_file(path)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def write_json(path: Path, obj) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default))
    os.replace(tmp, path)


def atomic_savez(path: Path, **arrays) -> None:
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def git_committed(path: Path) -> bool | None:
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


def git_head() -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True)
        return r.stdout.decode().strip() if r.returncode == 0 else None
    except OSError:
        return None


def code_hashes() -> dict:
    """sha256 of the source of every E5 function reused here and of the E7 statistics module."""
    fns = {"layer_fn": mld.layer_fn, "tokenise": mld.tokenise, "load_head": mld.load_head,
           "primary_keys": mld.primary_keys, "load_donors": mld.load_donors}
    out = {k: sha256_bytes(inspect.getsource(f).encode()) for k, f in fns.items()}
    out["lib_v3_e7"] = sha256_file(ROOT / "bridge_anm/lib/v3_e7.py")
    return out


def load_reg(a) -> dict:
    return va.load_registration_amended(a.registration_dir / "registration_v3.json",
                                        a.registration_dir / "amendment_A1.json")


def addenda_dir(a) -> Path:
    return a.registration_dir / "addenda"


def declared_json() -> dict:
    return json.loads(json.dumps(DECLARED))


def g_tokens(rna_names, vocab: dict) -> dict:
    """G symbol -> (Ensembl id, vocab id); every G gene must be in the processed genes and the vocabulary."""
    rn = set(map(str, rna_names))
    out, bad = {}, []
    for s_ in G_GENES:
        e = mld.TOKEN_PROBE_IDS[s_]
        if e in rn and e in vocab:
            out[s_] = (e, int(vocab[e]))
        else:
            bad.append(s_)
    if bad:
        raise SystemExit(f"G genes not available: {bad}")
    return out


def token_info(toks: list, gtok: set) -> tuple[np.ndarray, np.ndarray]:
    g = np.array([int(np.isin(t, list(gtok)).sum()) for t in toks], np.int64)
    n = np.array([t.size for t in toks], np.int64)
    return g, n - g


# ============================================================================ head mirror (float64 re-check)
def head_mirror(path: Path, meta: dict, device, dtype):
    """E5's load_head arithmetic in a chosen dtype (float64 on the CPU for the re-check); checked against load_head."""
    import torch

    from teddy_mm.models import MLP, AdtDecoder
    blob = torch.load(path, map_location="cpu", weights_only=False)
    w0 = blob["mlp"]["net.0.weight"]
    hidden, z_dim = int(w0.shape[0]), int(w0.shape[1])
    mlp, dec = MLP(z_dim, z_dim, hidden=hidden), AdtDecoder(z_dim, meta["adt"].shape[1], hidden=hidden)
    mlp.load_state_dict(blob["mlp"])
    dec.load_state_dict(blob["dec"])
    mlp.eval().to(device=device, dtype=dtype)
    dec.eval().to(device=device, dtype=dtype)
    sf = float(np.median(meta["adt_size_factor"][meta["split"] == "train"]))
    p95 = lp.p95_replicate(meta)
    col = {n: j for j, n in enumerate(meta["adt_names"])}
    cols = torch.tensor([col[p] for p in PROTEINS], device=device)
    p95v = torch.tensor([p95[p] for p in PROTEINS], dtype=dtype, device=device)

    def yhat(z):
        zn = z / (torch.linalg.vector_norm(z, dim=-1, keepdim=True) + 1e-6)
        mu = dec(mlp(zn), torch.full((z.shape[0],), sf, device=z.device, dtype=dtype))[0]
        return mu.index_select(-1, cols) / p95v
    return yhat


# ============================================================================ consumer context
class Ctx:
    """TEDDY (float32 on the device, CPU float32 layer 12 for the padding control, CPU float64 copy for the re-check)
    and the head in the same three forms."""

    def __init__(self, a, meta_bmmc: dict, with_f64: bool, clamp_dmu11):
        import torch

        from teddy_mm.device import resolve_device
        from teddy_mm.teddy_encoder import load_pad_id, load_teddy, load_vocab
        torch.set_num_threads(_THREADS)
        self.torch = torch
        self.device = resolve_device(a.device)
        self.vocab = load_vocab(a.ckpt)
        self.pad_id = load_pad_id(a.ckpt, self.vocab)
        self.model = load_teddy(a.ckpt, self.device)
        for p_ in self.model.parameters():
            p_.requires_grad_(False)
        self.layers = list(self.model.encoder.layers)
        self.E, self.Pm = self.model.embeddings.weight, self.model.position_embeddings.weight
        self.yhat, _, self.head_info = mld.load_head(a.head_ckpt, meta_bmmc, self.device)
        cpu = torch.device("cpu")
        self.l12_cpu = copy.deepcopy(self.layers[11]).to(cpu).float().eval()
        self.yhat_cpu = head_mirror(a.head_ckpt, meta_bmmc, cpu, torch.float32)
        cv = KAPPA * np.asarray(clamp_dmu11, np.float64)
        self.clamp = torch.as_tensor(cv, dtype=torch.float32, device=self.device)
        self.f64 = None
        if with_f64:
            m64 = copy.deepcopy(self.model).to(cpu).double().eval()
            self.f64 = {"layers": list(m64.encoder.layers), "E": m64.embeddings.weight, "Pm": m64.position_embeddings.weight,
                        "yhat": head_mirror(a.head_ckpt, meta_bmmc, cpu, torch.float64),
                        "clamp": torch.as_tensor(cv, dtype=torch.float64)}
        self.batch = int(a.batch)


def forward_to_cut(torch, layers, E, Pm, ids):
    """Explicit forward (E5's layer_fn) through layers 1..11: the layer-11 token states [L, d]."""
    L = ids.numel()
    h = (E[ids] + Pm[:L])[None]
    for layer in layers[:11]:
        h = mld.layer_fn(layer, h)
    return h[0]


def consumer(ctx_layers, yhat, Hb):
    """Fixed consumer on a batch of layer-11 states [B, L, d] -> head panel outputs [B, 4]."""
    return yhat(mld.layer_fn(ctx_layers[11], Hb).mean(1))


def build_histories(cid: int, ids_np: np.ndarray, gtok: set, fracs, phase: str) -> tuple[list, list, list]:
    """History list of (kind, magnitude code, patch pairs, base index, clamp sign) for one cell: per magnitude the
    matched patches and their G twins (base = the matched patch); then the G-only patches and the two clamp patches
    (code 0)."""
    g_pos = [q for q, t in enumerate(ids_np) if int(t) in gtok]
    nong = [q for q, t in enumerate(ids_np) if int(t) not in gtok]
    hist = []
    for f in fracs:
        code = v7.frac_code(f)
        m = v7.pairs_for_fraction(f, len(nong), len(g_pos))
        rng = np.random.default_rng([SEED, E7_TAG, PHASE_TAG[phase], int(cid), code])
        M, T = v7.draw_patch_sets(rng, g_pos, nong, m, N_MATCHED)
        i0 = len(hist)
        hist += [(v7.KIND_MATCHED, code, p_, -1, 0) for p_ in M]
        hist += [(v7.KIND_TWIN, code, p_, i0 + i, 0) for i, p_ in enumerate(T)]
    rng = np.random.default_rng([SEED, E7_TAG, PHASE_TAG[phase], int(cid), 1])
    hist += [(v7.KIND_G_ONLY, 0, p_, -1, 0) for p_ in v7.draw_g_only(rng, g_pos, nong, N_G_ONLY)]
    none = np.zeros((0, 2), np.int64)
    hist += [(v7.KIND_CLAMP, 0, none, -1, 1), (v7.KIND_CLAMP, 0, none, -1, -1)]
    return hist, g_pos, nong


def patched(H, hk, clamp_vec):
    """The layer-11 states of one history (pair averages, or the clamp shift +- clamp_vec)."""
    k, _, pairs, _, sign = hk
    if k == v7.KIND_CLAMP:
        return v7.apply_clamp(H, sign * clamp_vec)
    return v7.apply_average(H, pairs)


def run_cell(ctx: Ctx, cid: int, ids_np: np.ndarray, gtok: set, mags, phase: str, official_l11, do_f64: bool,
             do_push: bool) -> dict:
    """Every consumer output E7 needs for one cell (float32 on the device; controls; optional float64 re-check and
    input push)."""
    torch = ctx.torch
    t0 = time.time()
    hist, g_pos, nong = build_histories(cid, ids_np, gtok, mags, phase)
    L = int(ids_np.size)
    n = len(hist)
    ids = torch.from_numpy(ids_np).to(ctx.device)
    out = {}
    with torch.no_grad():
        H11 = forward_to_cut(torch, ctx.layers, ctx.E, ctx.Pm, ids)
        z11 = H11.mean(0)
        y12 = mld.layer_fn(ctx.layers[11], H11[None])
        O0 = ctx.yhat(y12.mean(1))[0]
        O = torch.zeros(n, len(PROTEINS), dtype=torch.float64)
        dmean = torch.zeros(n, dtype=torch.float64)
        for b0 in range(0, n, ctx.batch):
            hb = hist[b0:b0 + ctx.batch]
            Hb = torch.stack([patched(H11, hk, ctx.clamp) for hk in hb])
            O[b0:b0 + len(hb)] = consumer(ctx.layers, ctx.yhat, Hb).cpu().double()
            dmean[b0:b0 + len(hb)] = (Hb.mean(1) - z11[None]).abs().amax(-1).cpu().double()
            del Hb
        # ---- controls
        rc = np.random.default_rng([SEED, E7_TAG, PHASE_TAG[phase], int(cid), 0])
        perm = torch.as_tensor(rc.permutation(L), device=ctx.device)
        O_perm = consumer(ctx.layers, ctx.yhat, H11[perm][None])[0]
        ab = rc.choice(L, 2, replace=False)
        Y2 = v7.apply_average(y12[0], np.array([ab]))
        O_l12 = ctx.yhat(Y2.mean(0)[None])[0]
        Hc = H11.detach().float().cpu()
        y12c = y12[0].detach().float().cpu()
        pad_rows = torch.as_tensor(rc.choice(L, N_PAD, replace=True))
        padded = torch.cat([Hc, Hc[pad_rows].clone()], 0)[None]
        mask = torch.zeros(1, L + N_PAD, dtype=torch.bool)
        mask[0, L:] = True
        yp0 = ctx.l12_cpu(padded, src_key_padding_mask=mask)[0, :L]
        padded2 = padded.clone()
        padded2[0, L + 3] = 3.0 * 0.5 * (padded2[0, L + 3] + Hc[0])
        yp1 = ctx.l12_cpu(padded2, src_key_padding_mask=mask)[0, :L]
        O_pad0 = ctx.yhat_cpu(yp0.mean(0)[None])[0]
        O_pad1 = ctx.yhat_cpu(yp1.mean(0)[None])[0]
        ym = ctx.l12_cpu(Hc[None])[0]
        out["controls"] = {"module_vs_explicit_abs": float((ym - y12c).abs().max()),
                           "padded_vs_unpadded_abs": float((yp0 - y12c).abs().max()),
                           "official_layer11_gene_mean_abs": (float(np.abs(np.asarray(official_l11, np.float32)
                                                                           - z11.float().cpu().numpy()).max())
                                                              if official_l11 is not None else None),
                           "head_cpu_mirror_vs_device_abs": float((ctx.yhat_cpu(y12c.mean(0)[None])[0].double()
                                                                   - O0.cpu().double()).abs().max())}
        arr = {"O0": O0.cpu().double().numpy(), "O": O.numpy(), "dmean": dmean.numpy(),
               "O_perm": O_perm.cpu().double().numpy(), "O_l12": O_l12.cpu().double().numpy(),
               "O_pad0": O_pad0.double().numpy(), "O_pad1": O_pad1.double().numpy()}
        # ---- secondary positive control: input push on each G gene (token embedding doubled at the input)
        if do_push and g_pos:
            x0 = (ctx.E[ids] + ctx.Pm[:L])
            Op = []
            for b0 in range(0, len(g_pos), ctx.batch):
                xs = []
                for q in g_pos[b0:b0 + ctx.batch]:
                    x = x0.clone()
                    x[q] = x[q] + ctx.E[ids[q]]
                    xs.append(x)
                h = torch.stack(xs)
                for layer in ctx.layers:
                    h = mld.layer_fn(layer, h)
                Op.append(ctx.yhat(h.mean(1)).cpu().double().numpy())
                del h
            arr["O_push"] = np.concatenate(Op, 0)
        # ---- float64 re-check on the CPU (same histories)
        if do_f64:
            f = ctx.f64
            ids_c = ids.cpu()
            H64 = forward_to_cut(torch, f["layers"], f["E"], f["Pm"], ids_c)
            O0_64 = consumer(f["layers"], f["yhat"], H64[None])[0]
            O64 = torch.zeros(n, len(PROTEINS), dtype=torch.float64)
            for b0 in range(0, n, ctx.batch):
                hb = hist[b0:b0 + ctx.batch]
                Hb = torch.stack([patched(H64, hk, f["clamp"]) for hk in hb])
                O64[b0:b0 + len(hb)] = consumer(f["layers"], f["yhat"], Hb)
            arr["O0_64"], arr["O64"] = O0_64.numpy(), O64.numpy()
    offs = np.cumsum([0] + [len(h[2]) for h in hist]).astype(np.int64)
    arr.update({"kind": np.array([h[0] for h in hist], np.int8), "code": np.array([h[1] for h in hist], np.int32),
                "base": np.array([h[3] for h in hist], np.int32), "sign": np.array([h[4] for h in hist], np.int8),
                "pair_offsets": offs, "pairs": (np.concatenate([h[2] for h in hist]).astype(np.int32) if hist
                                                else np.zeros((0, 2), np.int32)),
                "g_pos": np.array(g_pos, np.int32), "ntokens": np.array(L)})
    out["arrays"] = arr
    out["sec"] = round(time.time() - t0, 2)
    out["n_hist"] = n
    out["g_present"] = int(len(g_pos))
    out["n_nong"] = int(len(nong))
    return out


# ============================================================================ register
def unit_core(meta: dict, key: np.ndarray, zrows: "lp_rows") -> dict:
    """Units from val rows only: measured (primary) and predicted (secondary), per protein, and for the score."""
    split = np.asarray(meta["split"]).astype(str)
    p95 = lp.p95_replicate(meta)
    col = {n: j for j, n in enumerate(meta["adt_names"])}
    nk = np.where((split == "val") & (key == "NK"))[0]
    t = np.where((split == "val") & (key == "T"))[0]
    M = np.stack([np.asarray(meta["adt"][:, col[p]], np.float64) / p95[p] for p in PROTEINS], 1)
    u = np.abs(np.median(M[nk], 0) - np.median(M[t], 0))
    s = v7.score(M)
    us = float(abs(np.median(s[nk]) - np.median(s[t])))
    kept = [p for p, x in zip(PROTEINS, u) if x >= UNIT_MIN]
    pred = zrows(np.concatenate([nk, t]))
    upred = np.abs(np.median(pred[: nk.size], 0) - np.median(pred[nk.size:], 0))
    return {"panel": [float(x) for x in u], "proteins": list(PROTEINS), "proteins_kept": kept, "score": us,
            "p95": {p: float(p95[p]) for p in PROTEINS}, "n_val_NK": int(nk.size), "n_val_T": int(t.size),
            "predicted_panel_secondary": [float(x) for x in upred],
            "predicted_score_secondary": float(abs(np.median(v7.score(pred[: nk.size])) - np.median(v7.score(pred[nk.size:]))))}


def rosters_core(meta: dict, key: np.ndarray, emb, a, pad_id: int, gtok: set) -> dict:
    """Development, follow-up and smoke rosters from val donor 18303 primary-key NK/T cells (label + tokenisation;
    no consumer output)."""
    split = np.asarray(meta["split"]).astype(str)
    pool = {c: np.where((split == "val") & (key == c))[0] for c in ("NK", "T")}
    elig, info = {}, {}
    for c, rows in pool.items():
        toks, ntok = mld.tokenise(meta, rows, a, emb, pad_id)
        if np.any(ntok != emb.ntokens[rows]):
            raise SystemExit("val token counts differ from the official run")
        g, ng = token_info(toks, gtok)
        ok = (g >= MIN_G_PRESENT) & (ng >= MIN_NONG_TOKENS)
        elig[c] = rows[ok]
        info[c] = {"pool": int(rows.size), "eligible": int(ok.sum()), "no_G": int((g < MIN_G_PRESENT).sum()),
                   "too_few_nonG": int(((g >= MIN_G_PRESENT) & (ng < MIN_NONG_TOKENS)).sum()),
                   "g_present_counts": {int(k): int(v) for k, v in zip(*np.unique(g[ok], return_counts=True))}}
    rng = np.random.default_rng([SEED, E7_TAG, 0])
    dev = {c: np.sort(rng.choice(elig[c], min(N_DEV_PER_CLASS, elig[c].size), replace=False)) for c in ("NK", "T")}
    rest = {c: np.setdiff1d(elig[c], dev[c]) for c in ("NK", "T")}
    f1 = {c: np.sort(rng.choice(rest[c], min(N_F1_PER_CLASS, rest[c].size), replace=False)) for c in ("NK", "T")}
    smoke = {c: np.setdiff1d(rest[c], f1[c]) for c in ("NK", "T")}
    dcells = np.concatenate([dev["NK"], dev["T"]])
    fcells = np.concatenate([f1["NK"], f1["T"]])
    dorder = rng.permutation(dcells.size)
    forder = rng.permutation(fcells.size)

    def ro(cells, order, cls):
        return {"cells": [int(x) for x in cells], "class": cls, "order": [int(x) for x in order],
                "sha256": vk.index_hash(cells)}
    return {"eligibility": info,
            "dev": ro(dcells, dorder, ["NK"] * dev["NK"].size + ["T"] * dev["T"].size),
            "followup": ro(fcells, forder, ["NK"] * f1["NK"].size + ["T"] * f1["T"].size),
            "smoke_pool": {"NK": [int(x) for x in smoke["NK"]], "T": [int(x) for x in smoke["T"]],
                           "sha256": vk.index_hash(np.concatenate([smoke["NK"], smoke["T"]]))}}


def core_sha(core: dict) -> str:
    return sha256_bytes(json.dumps(core, sort_keys=True, default=_json_default).encode())


def poisoned_meta(meta: dict, rows: np.ndarray, seed: int, rna: bool) -> dict:
    """Copy of meta whose chosen rows of adt, adt_size_factor, cell_types (and RNA values) are random."""
    pm = dict(meta)
    rng = np.random.default_rng(seed)
    pm["adt"] = meta["adt"].copy()
    pm["adt"][rows] = rng.uniform(0, 8, (rows.size, meta["adt"].shape[1])).astype(np.float32)
    pm["adt_size_factor"] = np.asarray(meta["adt_size_factor"]).copy()
    pm["adt_size_factor"][rows] = rng.uniform(0.1, 10, rows.size)
    pm["cell_types"] = np.asarray(meta["cell_types"]).copy()
    pm["cell_types"][rows] = rng.choice(np.unique(meta["cell_types"]), rows.size)
    if rna and "rna" in meta:
        pm["rna"] = PoisonedRows(meta["rna"], rows, seed)
    return pm


class PoisonedRows:
    """Row access to a CSR matrix whose chosen rows get seeded U(0.5, 50) values on their own sparsity (lazy, so the
    leakage check does not copy the whole RNA matrix)."""

    def __init__(self, X, rows, seed: int):
        self.X, self.rows, self.seed, self.shape = X, set(int(r) for r in rows), int(seed), X.shape

    def __getitem__(self, idx):
        idx = np.asarray(idx, np.int64)
        sub = self.X[idx].tocsr(copy=True)
        for k, i in enumerate(idx):
            if int(i) in self.rows:
                s_, e_ = sub.indptr[k], sub.indptr[k + 1]
                sub.data[s_:e_] = np.random.default_rng([self.seed, int(i)]).uniform(0.5, 50, e_ - s_)
        return sub


def architecture_check(a, meta: dict, emb) -> dict:
    """The architecture facts the design relies on, measured on 3 training cells (CPU float32 / float64)."""
    import torch

    from teddy_mm.teddy_encoder import load_pad_id, load_teddy, load_vocab
    torch.set_num_threads(_THREADS)
    split = np.asarray(meta["split"]).astype(str)
    cells = np.sort(np.random.default_rng([SEED, E7_TAG, 5]).choice(np.where(split == "train")[0], 3, replace=False))
    vocab = load_vocab(a.ckpt)
    pad = load_pad_id(a.ckpt, vocab)
    toks, _ = mld.tokenise(meta, cells, a, emb, pad)
    cpu = torch.device("cpu")
    model = load_teddy(a.ckpt, cpu)
    layers = list(model.encoder.layers)
    yhat, _, _ = mld.load_head(a.head_ckpt, meta, cpu)
    y64 = head_mirror(a.head_ckpt, meta, cpu, torch.float64)
    ym32 = head_mirror(a.head_ckpt, meta, cpu, torch.float32)
    m64 = copy.deepcopy(model).double()
    l12 = layers[11]
    rows = []
    with torch.no_grad():
        for ci, c in enumerate(cells):
            ids = torch.from_numpy(toks[ci])
            L = ids.numel()
            H = forward_to_cut(torch, layers, model.embeddings.weight, model.position_embeddings.weight, ids)
            y = mld.layer_fn(l12, H[None])
            O0 = yhat(y.mean(1))[0]
            rg = np.random.default_rng([SEED, E7_TAG, 6, int(c)])
            perm = torch.as_tensor(rg.permutation(L))
            yp = mld.layer_fn(l12, H[perm][None])
            ab = rg.choice(L, 2, replace=False)
            abt = torch.as_tensor(ab)
            Hs = H.clone()
            Hs[abt] = H[abt.flip(0)]
            Ha = v7.apply_average(H, np.array([ab]))
            padded = torch.cat([H, H[:N_PAD].clone()], 0)[None]
            mask = torch.zeros(1, L + N_PAD, dtype=torch.bool)
            mask[0, L:] = True
            yp0 = l12(padded, src_key_padding_mask=mask)[0, :L]
            padded[0, L + 3] = 3.0 * padded[0, L + 3]
            yp1 = l12(padded, src_key_padding_mask=mask)[0, :L]
            H64 = forward_to_cut(torch, list(m64.encoder.layers), m64.embeddings.weight, m64.position_embeddings.weight, ids)
            rows.append({
                "cell": int(c), "ntokens": int(L),
                "permutation_state_equivariance_abs": float((yp[0] - y[0, perm]).abs().max()),
                "permutation_head_abs": float((yhat(yp.mean(1))[0] - O0).abs().max()),
                "swap_patch_head_abs": float((yhat(mld.layer_fn(l12, Hs[None]).mean(1))[0] - O0).abs().max()),
                "average_patch_head_abs": float((yhat(mld.layer_fn(l12, Ha[None]).mean(1))[0] - O0).abs().max()),
                "average_patch_gene_mean_abs": float((Ha.mean(0) - H.mean(0)).abs().max()),
                "module_vs_explicit_abs": float((l12(H[None])[0] - y[0]).abs().max()),
                "padded_masked_vs_unpadded_abs": float((yp0 - y[0]).abs().max()),
                "padding_patch_state_abs": float((yp1 - yp0).abs().max()),
                "layer12_average_patch_head_abs": float((yhat(v7.apply_average(y[0], np.array([ab])).mean(0)[None])[0] - O0)
                                                        .abs().max()),
                "explicit_vs_official_layer11_gene_mean_abs": float(np.abs(np.asarray(emb.layer_means[10][c], np.float32)
                                                                           - H.mean(0).numpy()).max()),
                "explicit_vs_official_layer12_gene_mean_abs": float(np.abs(np.asarray(emb.layer_means[11][c], np.float32)
                                                                           - y[0].mean(0).numpy()).max()),
                "float32_vs_float64_layer11_abs": float((H64 - H.double()).abs().max()),
                "head_mirror32_vs_load_head_abs": float((ym32(y.mean(1))[0] - O0).abs().max()),
                "head_mirror64_vs_load_head_abs": float((y64(y.mean(1).double())[0] - O0.double()).abs().max()),
            })
    return {"cells": rows, "note": "3 training cells (default_rng([20260930, 7, 5])); one average and one swap patch per "
                                   "cell at positions default_rng([20260930, 7, 6, cell]); no val or site4 cell",
            "layer12": {"type": type(l12).__name__, "norm_first": bool(l12.norm_first),
                        "activation": getattr(l12.activation, "__name__", str(l12.activation)),
                        "parameters": [n for n, _ in l12.named_parameters()], "training_mode": bool(l12.training)}}


def stage_register(a) -> None:
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    if not a.smoke:  # the addendum can only be written before any development output or selection record exists
        done = [str(x) for x in (a.out_dir / "dev" / "cells", a.out_dir / "followup" / "cells", a.out_dir / "site4",
                                 a.out_dir / "external") if x.exists() and any(x.iterdir())]
        if (addenda_dir(a) / SELECTION).exists():
            done.append(str(addenda_dir(a) / SELECTION))
        if done:
            raise SystemExit("register refused: development or confirmation outputs already exist: " + ", ".join(done))
    from teddy_mm.teddy_encoder import load_pad_id, load_vocab
    reg = load_reg(a)
    meta = lp.load_meta(a.processed, need_rna=True)
    split = np.asarray(meta["split"]).astype(str)
    emb = lp.Embedding(a.embed_dir, None, len(split))
    vocab = load_vocab(a.ckpt)
    pad = load_pad_id(a.ckpt, vocab)
    G = g_tokens(meta["rna_names"], vocab)
    gtok = {v[1] for v in G.values()}
    import torch
    yhat_cpu = head_mirror(a.head_ckpt, meta, torch.device("cpu"), torch.float32)

    def zrows_factory(poison_rows=None, pseed=0):
        rv = mld_rowview(emb.z, poison_rows, pseed)

        def zrows(ix):
            with torch.no_grad():
                return yhat_cpu(torch.from_numpy(rv.rows(ix).astype(np.float32))).double().numpy()
        return zrows

    def build(meta_, zr):
        key_ = mld.primary_keys(meta_, reg)
        return {"G": {s_: {"ensembl": e, "vocab_id": t} for s_, (e, t) in G.items()},
                "units": unit_core(meta_, key_, zr), "rosters": rosters_core(meta_, key_, emb, a, pad, gtok)}

    t0 = time.time()
    core = build(meta, zrows_factory())
    real = core_sha(core)
    say(f"register: core {real[:16]} ({time.time() - t0:.0f}s); units {core['units']['panel']} score {core['units']['score']:.4f}; "
        f"eligibility {core['rosters']['eligibility']}")
    test_rows = np.where(split == "test")[0]
    val_rows = np.where(split == "val")[0]
    p_site4 = build(poisoned_meta(meta, test_rows, 101, rna=True), zrows_factory(test_rows, 101))
    pm_val = dict(meta)
    rng = np.random.default_rng(303)
    pm_val["adt"] = meta["adt"].copy()
    pm_val["adt"][val_rows] = rng.uniform(0, 8, (val_rows.size, meta["adt"].shape[1])).astype(np.float32)
    p_val = build(pm_val, zrows_factory(val_rows, 303))
    leak = {"what": "site4 rows (split test) of adt -> U(0, 8), adt_size_factor -> U(0.1, 10), cell_types -> random type, "
                    "RNA values -> U(0.5, 50) on the same sparsity, official z -> N(0, 1) (seed 101); positive control: val "
                    "rows of adt -> U(0, 8) and of z -> N(0, 1) (seed 303)",
            "n_site4_rows": int(test_rows.size), "real_core_sha256": real, "site4_poisoned_core_sha256": core_sha(p_site4),
            "val_poisoned_core_sha256": core_sha(p_val)}
    leak["site4_poisoned_identical"] = leak["site4_poisoned_core_sha256"] == real
    leak["val_poisoned_differs"] = leak["val_poisoned_core_sha256"] != real
    leak["passed"] = bool(leak["site4_poisoned_identical"] and leak["val_poisoned_differs"])
    leak["static"] = ("the builder tokenises val rows only, reads val ADT / cell types / z for the units, E5's "
                      "training-only dmu[11] (added to computed_train_val after the check, as a file constant), and of "
                      "site4 only the sha256 of E5's and E5-M's cell files (bytes hashed) and their published counts; no "
                      "site4 protein, cell type, embedding or RNA value enters the core")
    say(f"register: leakage check site4-poisoned identical {leak['site4_poisoned_identical']}, val-poisoned differs "
        f"{leak['val_poisoned_differs']}")
    if not leak["passed"]:
        raise SystemExit("leakage check failed; nothing written")
    # clamp direction (positive control): E5's training NK - T difference of the layer-11 gene-mean
    dj = json.loads((a.e5_dir / "e5_directions.json").read_text())
    if dj.get("subsampled"):
        raise SystemExit("E5's e5_directions.npz was built from a smoke subsample")
    with np.load(a.e5_dir / "e5_directions.npz") as d_:
        dmu11 = d_["dmu"][11].astype(np.float64)
    z11_norm = float(np.median(np.linalg.norm(np.asarray(emb.layer_means[10][np.where(split == "val")[0][:500]],
                                                         np.float64), axis=1)))
    clamp = {"dmu_11": dmu11.tolist(), "source": str(a.e5_dir / "e5_directions.npz"),
             "source_sha256": sha256_file(a.e5_dir / "e5_directions.npz"), "kappa": KAPPA,
             "norm_dmu_11": float(np.linalg.norm(dmu11)), "norm_clamp_shift": float(KAPPA * np.linalg.norm(dmu11)),
             "median_norm_layer11_gene_mean_val500": z11_norm,
             "definition": "dmu[11] of E5's directions stage = mean over training primary-key NK cells - mean over training "
                           "primary-key T cells of the layer-11 output gene-mean (z_rna_layer_means[10]); training rows only"}
    core = {**core, "clamp": clamp}
    arch = architecture_check(a, meta, emb)
    say("register: architecture check " + json.dumps(arch["cells"][0])[:300])
    # exclusion sources for site4 (cell index files of E5 / E5-M; counts only)
    e5c, e5j = a.e5_dir / "e5_cells.npz", a.e5_dir / "e5_cells.json"
    e5mc = a.e5m_dir / "e5m_cells.npz"
    e5info = json.loads(e5j.read_text())
    excl = {"e5_cells_file": str(e5c), "e5_cells_sha256": sha256_file(e5c), "e5m_cells_file": str(e5mc),
            "e5m_cells_sha256": sha256_file(e5mc), "e5_n_cells": e5info["n_cells"],
            "e3_pair_cells_per_donor": e5info["n_e3_pair_cells_per_donor"], "e3_pairs_sha256": e5info["e3_pairs_sha256"],
            "rule": "exclude every cell in e5_cells.npz (roles pair = all of E3's NK-T pair cells, checked against "
                    "n_e3_pair_cells_per_donor; random) and in e5m_cells.npz"}
    e6 = json.loads((addenda_dir(a) / "E6.json").read_text())
    ext = {"E6_addendum_sha256": sha256_file(addenda_dir(a) / "E6.json"),
           "pack": str(a.external_pack / "cite_arrays.npz"), "pack_sha256_E6": e6["external_data"]["pack_sha256"],
           "annotation_map_sha256_E6": e6["annotation_map"].get("sha256_canonical_json"),
           "embedding_dir": str(a.external_embed_dir), "labels": [e6["contamination"]["label"], "key not validated"],
           "donors": sorted(e6["external_data"]["donors"])}
    amend = {am.name: sha256_file(am) for am in sorted(a.registration_dir.glob("amendment_A*.json"))
             if not am.name.endswith("_core.json")}
    add = {
        "experiment": "E7",
        "addendum_version": 1,
        "title": "E7: TEDDY analog of ANM v2's inverse loop at the layer-11 cut, consumer fixed (Mode A)",
        "addendum_to": "registration/registration_v3.json as amended by A1, A2 and A3: a new pre-specified Mode A experiment "
                       "registered before any development or site4 forward of its script; nothing registered for E1-E6 "
                       "changes",
        "registration_sha256": sha256_file(a.registration_dir / "registration_v3.json"),
        "amendments_sha256": amend,
        "script": f"scripts/mode_a_inverse_e7.py ({SCRIPT_VERSION}); statistics in bridge_anm/lib/v3_e7.py",
        "declared": DECLARED,
        "constants": {"tol": TOL, "ladder": list(LADDER), "ladder_codes": [v7.frac_code(f) for f in LADDER],
                      "fresh": list(FRESH), "fresh_codes": [v7.frac_code(f) for f in FRESH], "kappa": KAPPA,
                      "n_matched": N_MATCHED, "n_twin": N_MATCHED, "n_g_only": N_G_ONLY, "n_clamp": 2,
                      "n_dev_per_class": N_DEV_PER_CLASS, "n_followup_per_class": N_F1_PER_CLASS,
                      "n_site4_per_class_per_donor": N_SITE4_PER_CLASS_PER_DONOR,
                      "n_external_per_class_per_donor": N_EXT_PER_CLASS_PER_DONOR, "min_nonG_tokens": MIN_NONG_TOKENS,
                      "min_G_present": MIN_G_PRESENT, "n_float64": N_F64, "float64_tol": F64_TOL,
                      "n_push_cells": N_PUSH_CELLS, "n_pad": N_PAD, "unit_min": UNIT_MIN, "n_boot": N_BOOT,
                      "bootstrap_seed": int(reg["seeds"]["bootstrap"]), "gates": GATES, "seed_base": [SEED, E7_TAG],
                      "phase_tags": PHASE_TAG, "proteins": list(PROTEINS)},
        "computed_train_val": core,
        "architecture_check": arch,
        "site4_exclusions": excl,
        "external_family": ext,
        "inputs": {"code_sha256": code_hashes(), "head_sha256": sha256_file(a.head_ckpt), "head": str(a.head_ckpt),
                   "teddy_weights_sha256": sha256_file(a.ckpt / "model.safetensors"),
                   "medians_sha256": sha256_file(a.medians), "official_embedding_manifest": str(a.embed_dir / "z_rna_manifest.json")},
        "leakage_check": leak,
        "disclosure": "The registering session had read the published v3 results (E1-E6, E5-M, including E5-M's per-layer "
                      "clamp residuals on site4: the head's 0.709 / 0.72 at layer 11), the PR #16 review and the val / "
                      "train cell counts. Before writing this addendum it ran the architecture check on 3 training cells "
                      "(one average and one swap patch each, recorded in architecture_check) and counted val primary-key "
                      "NK/T cells; it computed no patch output on any val or site4 cell. Site4 confirmation cells exclude "
                      "every E3 / E5 / E5-M cell.",
        "smoke": bool(a.smoke),
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    d = addenda_dir(a) if not a.smoke else a.out_dir
    d.mkdir(parents=True, exist_ok=True)
    f = d / ADDENDUM
    f.write_bytes((json.dumps(add, indent=1, sort_keys=True, default=_json_default) + "\n").encode())
    sha = sha256_file(f)
    hh = d / "HASHES.txt"
    lines = [ln for ln in (hh.read_text().splitlines() if hh.exists() else []) if not ln.rstrip().endswith(f"  {ADDENDUM}")]
    hh.write_text("\n".join(lines + [f"{sha}  {ADDENDUM}"]) + "\n")
    say(f"register: wrote {f} (sha256 {sha}) and its line in {hh}")


class mld_rowview:
    """Row access to an [N, d] array whose chosen rows are replaced by seeded N(0, 1) values (leakage check)."""

    def __init__(self, arr, poison_rows=None, seed: int = 0):
        self.arr = arr
        self.poison = set(int(x) for x in poison_rows) if poison_rows is not None else set()
        self.seed = seed

    def rows(self, ix):
        ix = np.asarray(ix, np.int64)
        X = np.asarray(self.arr[ix], np.float64)
        hit = [k for k, i in enumerate(ix) if int(i) in self.poison]
        if hit:
            X[hit] = np.random.default_rng(self.seed).standard_normal((len(hit), X.shape[1]))
        return X


# ============================================================================ registration guard
def load_addendum(a) -> tuple[dict, Path]:
    f = (a.out_dir / ADDENDUM) if (a.smoke and (a.out_dir / ADDENDUM).exists()) else addenda_dir(a) / ADDENDUM
    if not f.exists():
        raise SystemExit(f"no addendum at {f}: run --stage register first")
    return json.loads(f.read_text()), f


def selection_path(a) -> Path:
    return (a.out_dir / SELECTION) if a.smoke else addenda_dir(a) / SELECTION


def check_registration(a, need: str) -> dict:
    """need in develop / followup / select / confirm / report. Smoke runs only record the state."""
    add, f = load_addendum(a)
    hh = addenda_dir(a) / "HASHES.txt"
    hl = hh.read_text().splitlines() if hh.exists() else []
    sha = sha256_file(f)
    info = {"addendum_file": str(f), "addendum_sha256": sha,
            "addendum_hash_recorded": f"{sha}  {ADDENDUM}" in hl, "addendum_committed": git_committed(f),
            "declared_equals_script": add.get("declared") == declared_json(),
            "leakage_check_passed": add.get("leakage_check", {}).get("passed") is True,
            "smoke_addendum": bool(add.get("smoke")),
            "code_now": code_hashes(), "script_sha256": sha256_file(Path(__file__).resolve()), "git_head": git_head()}
    info["code_matches_addendum"] = info["code_now"] == add.get("inputs", {}).get("code_sha256")
    info["head_matches_addendum"] = Path(a.head_ckpt).exists() and sha256_file(a.head_ckpt) == \
        add.get("inputs", {}).get("head_sha256")
    cf, ch = a.registration_dir / "registration_v3.json", a.registration_dir / "registration_v3.json.sha256"
    info["registration_ok"] = cf.exists() and ch.exists() and ch.read_text().split()[0] == sha256_file(cf) and \
        git_committed(cf) is True
    info["amendments_match"] = all(sha256_file(a.registration_dir / n) == s_ and git_committed(a.registration_dir / n) is True
                                   for n, s_ in add.get("amendments_sha256", {}).items())
    sp = selection_path(a)
    info["selection_file"] = str(sp)
    sel = json.loads(sp.read_text()) if sp.exists() else None
    if sel is not None:
        ss = sha256_file(sp)
        info.update({"selection_sha256": ss, "selection_hash_recorded": f"{ss}  {SELECTION}" in hl,
                     "selection_committed": git_committed(sp), "selection_addendum_matches": sel.get("addendum_sha256") == sha,
                     "selection_final": bool(sel.get("final")), "selection_smoke": bool(sel.get("smoke"))})
    if a.smoke:
        return info
    base = {"addendum committed, sha256 in addenda/HASHES.txt": info["addendum_committed"] is True and info["addendum_hash_recorded"],
            "addendum declared block equals this script's DECLARED": info["declared_equals_script"],
            "addendum leakage check passed, not a smoke addendum": info["leakage_check_passed"] and not info["smoke_addendum"],
            "reused E5 functions and bridge_anm/lib/v3_e7.py are the versions the addendum records": info["code_matches_addendum"],
            "head checkpoint is the one the addendum records": info["head_matches_addendum"],
            "registration_v3.json and the amendments the addendum records are committed and match": info["registration_ok"] and info["amendments_match"]}
    need_sel = {}
    if need == "followup":
        need_sel = {"the null primary selection record is committed (sha256 in HASHES.txt) and names this addendum":
                    sel is not None and info.get("selection_committed") is True and info.get("selection_hash_recorded")
                    and info.get("selection_addendum_matches") and sel.get("primary", {}).get("selected") is None
                    and sel.get("follow_up_required") is True}
    if need in ("confirm", "report"):
        need_sel = {"the final development selection record is committed (sha256 in HASHES.txt), names this addendum and is not a smoke record":
                    sel is not None and info.get("selection_committed") is True and info.get("selection_hash_recorded")
                    and info.get("selection_addendum_matches") and info.get("selection_final") and not info.get("selection_smoke")}
    bad = [k for k, ok in {**base, **need_sel}.items() if not ok]
    if bad:
        what = "site4 / external run" if need in ("confirm", "report") else f"{need} run"
        raise SystemExit(f"{what} refused: " + "; ".join(bad))
    return info


# ============================================================================ cells per stage
def phase_cells(a, add: dict, meta: dict, phase: str) -> tuple[np.ndarray, list, np.ndarray]:
    """(cells, classes, processing order) of a development phase or a smoke pool."""
    R = add["computed_train_val"]["rosters"]
    if a.cell_pool == "smoke":
        rng = np.random.default_rng([SEED, E7_TAG, 9])
        n_nk = (a.n_smoke + 1) // 2
        nk = np.sort(rng.choice(np.asarray(R["smoke_pool"]["NK"]), n_nk, replace=False))
        t = np.sort(rng.choice(np.asarray(R["smoke_pool"]["T"]), a.n_smoke - n_nk, replace=False))
        cells = np.concatenate([nk, t])
        return cells, ["NK"] * nk.size + ["T"] * t.size, rng.permutation(cells.size)
    if a.cell_pool == "train":
        with np.load(a.processed / "cite_arrays.npz", allow_pickle=True) as d:
            split, ctype = d["split"].astype(str), d["cell_types"].astype(str)
        ntok = np.load(a.embed_dir / "z_rna_ntokens.npy")
        cls = np.array([lp.ANNOTATION_MAP.get(t) or "" for t in ctype])
        rng = np.random.default_rng([SEED, E7_TAG, 8])
        cells, classes = [], []
        for c, n_ in (("NK", (a.n_smoke + 1) // 2), ("T", a.n_smoke // 2)):  # annotated NK / T training cells (code tests)
            pick = np.sort(rng.choice(np.where((split == "train") & (ntok >= 600) & (cls == c))[0], n_, replace=False))
            cells += [int(x) for x in pick]
            classes += [c] * n_
        return np.asarray(cells, np.int64), classes, np.arange(len(cells))
    r = R[phase]
    return np.asarray(r["cells"], np.int64), list(r["class"]), np.asarray(r["order"], np.int64)


def family_roster(a, add: dict, family: str, meta: dict, emb, pad_id: int, gtok: set, reg: dict) -> dict:
    """Confirmation roster (deterministic from the declared rule); written once to <out>/<family>/e7_roster.json."""
    out = a.out_dir / family
    f = out / "e7_roster.json"
    if f.exists():
        return json.loads(f.read_text())
    out.mkdir(parents=True, exist_ok=True)
    if family == "site4":
        split = np.asarray(meta["split"]).astype(str)
        donors = mld.load_donors(a.processed)
        key = mld.primary_keys(meta, reg)
        ex = add["site4_exclusions"]
        if sha256_file(Path(ex["e5_cells_file"])) != ex["e5_cells_sha256"] or sha256_file(Path(ex["e5m_cells_file"])) != ex["e5m_cells_sha256"]:
            raise SystemExit("E5 / E5-M cell files differ from the ones the addendum pins")
        with np.load(ex["e5_cells_file"], allow_pickle=False) as d:
            e5c, e5r = d["cells"].astype(np.int64), d["role"].astype(str)
        with np.load(ex["e5m_cells_file"], allow_pickle=False) as d:
            e5m = d["cells"].astype(np.int64)
        for dn, n_ in ex["e3_pair_cells_per_donor"].items():
            if int(np.sum((e5r == "pair") & (donors[e5c] == dn))) != int(n_):
                raise SystemExit("E5's pair cells are not all of E3's pair cells")
        excl = set(map(int, e5c)) | set(map(int, e5m))
        prim = list(reg["splits"]["test_primary"]["donors"])
        rng = np.random.default_rng([SEED, E7_TAG, 3])
        cap = N_SITE4_PER_CLASS_PER_DONOR
        groups = [(dn, c, np.array([i for i in np.where((split == "test") & (donors == dn) & (key == c))[0] if int(i) not in excl]))
                  for dn in prim for c in ("NK", "T")]
        cls_of = {}
    else:
        pack = ext_meta(a, add)
        donors, cls = pack["donors"], pack["class"]
        rng = np.random.default_rng([SEED, E7_TAG, 4])
        cap = N_EXT_PER_CLASS_PER_DONOR
        groups = [(dn, c, np.where((donors == dn) & (cls == c))[0]) for dn in add["external_family"]["donors"] for c in ("NK", "T")]
        meta, emb = pack, ext_embedding(a)
    cells, classes, dons, info = [], [], [], {}
    for dn, c, rows in groups:
        if rows.size:
            toks, ntok = mld.tokenise(meta, rows, a, emb, pad_id)
            if np.any(ntok != emb.ntokens[rows]):
                raise SystemExit(f"{family}: token counts differ from the official run")
            g, ng = token_info(toks, gtok)
            ok = rows[(g >= MIN_G_PRESENT) & (ng >= MIN_NONG_TOKENS)]
        else:
            ok = rows
        pick = np.sort(rng.choice(ok, min(cap, ok.size), replace=False)) if ok.size else ok
        info[f"{dn}/{c}"] = {"pool": int(rows.size), "eligible": int(ok.size), "drawn": int(pick.size)}
        cells += [int(x) for x in pick]
        classes += [c] * pick.size
        dons += [str(dn)] * pick.size
    order = rng.permutation(len(cells))
    R = {"family": family, "cells": cells, "class": classes, "donor": dons, "order": [int(x) for x in order],
         "counts": info, "sha256": vk.index_hash(cells), "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    write_json(f, R)
    say(f"{family} roster: {info}")
    return R


def ext_meta(a, add: dict) -> dict:
    from scipy import sparse

    from lib import v3_e6 as e6
    pf = a.external_pack / "cite_arrays.npz"
    if sha256_file(pf) != add["external_family"]["pack_sha256_E6"]:
        raise SystemExit("external pack differs from the one E6 pins")
    e6add = json.loads((addenda_dir(a) / "E6.json").read_text())
    if sha256_file(addenda_dir(a) / "E6.json") != add["external_family"]["E6_addendum_sha256"]:
        raise SystemExit("E6 addendum differs from the one E7 pins")
    d = np.load(pf, allow_pickle=False)
    meta = {"rna": sparse.csr_matrix((d["rna_data"], d["rna_indices"], d["rna_indptr"]), shape=tuple(d["rna_shape"])),
            "token_ids": d["token_ids"], "rna_names": d["rna_names"], "donors": d["donors"].astype(str),
            "class": e6.annotation_class(d["cell_types"].astype(str), e6add["annotation_map"]["map"])}
    return meta


def ext_embedding(a):
    n = int(np.load(a.external_embed_dir / "z_rna_ntokens.npy").size)
    return lp.Embedding(a.external_embed_dir, None, n)


# ============================================================================ develop / confirm (device)
def run_phase(a, add: dict, phase: str, cells: np.ndarray, classes: list, order: np.ndarray, donors: list | None,
              meta: dict, emb, meta_bmmc: dict, mags) -> int:
    from teddy_mm.teddy_encoder import load_vocab
    out = a.out_dir / (phase if not a.smoke else f"smoke_{phase}")
    cdir = out / "cells"
    cdir.mkdir(parents=True, exist_ok=True)
    sha = sha256_file(load_addendum(a)[1])
    gtok = {int(v["vocab_id"]) for v in add["computed_train_val"]["G"].values()}
    vocab = load_vocab(a.ckpt)
    tok2sym = {int(v["vocab_id"]): s_ for s_, v in add["computed_train_val"]["G"].items()}
    if any(vocab.get(add["computed_train_val"]["G"][s_]["ensembl"]) != t for t, s_ in tok2sym.items()):
        raise SystemExit("G vocabulary ids differ from the addendum")
    f64_set = {int(i) for i in order[: (len(order) if a.smoke else N_F64)]}
    push_set = {int(i) for i in order[:N_PUSH_CELLS]} if phase in ("site4", "external", "smoke_family") else set()
    todo = []
    for i in order:
        f = cdir / f"cell_{int(cells[i])}.npz"
        if f.exists():
            with np.load(f, allow_pickle=False) as d:
                if json.loads(str(d["meta"])).get("addendum_sha256") != sha:
                    raise SystemExit(f"{f} was written under another addendum; use a fresh --out-dir")
        else:
            todo.append(int(i))
    say(f"{phase}: {len(order)} cells, {len(order) - len(todo)} done, {len(todo)} to go; magnitudes {list(mags)}")
    if not todo:
        return 0
    ctx = Ctx(a, meta_bmmc, with_f64=bool(f64_set & set(todo)), clamp_dmu11=add["computed_train_val"]["clamp"]["dmu_11"])
    say(f"{phase}: device {ctx.device}, float32, batch {ctx.batch}")
    toks, ntok = mld.tokenise(meta, cells[todo], a, emb, ctx.pad_id)
    if np.any(ntok != emb.ntokens[cells[todo]]):
        raise SystemExit("token counts differ from the official run")
    t_start, n_done = time.time(), 0
    for jj, i in enumerate(todo):
        if a.max_minutes and time.time() - t_start > 60 * a.max_minutes:
            say(f"{phase}: --max-minutes {a.max_minutes} reached after {n_done} cells; rerun the same command to resume")
            return EXIT_INCOMPLETE
        cid = int(cells[i])
        off = np.asarray(emb.layer_means[10][cid], np.float32) if emb is not None else None
        tag = phase if not a.smoke else (a.cell_pool if phase in ("dev", "followup") else f"{a.cell_pool}_family")
        r = run_cell(ctx, cid, toks[jj], gtok, mags, tag, off, int(i) in f64_set, int(i) in push_set)
        mc = {"cell": cid, "class": classes[i], "donor": str(donors[i]), "phase": phase, "patch_tag": tag,
              "ntokens": int(toks[jj].size), "g_present": r["g_present"], "n_nonG": r["n_nong"],
              "g_genes": sorted({tok2sym[int(toks[jj][q])] for q in r["arrays"]["g_pos"]}), "n_hist": r["n_hist"],
              "magnitudes": [float(m) for m in mags], "controls": r["controls"], "f64": int(i) in f64_set,
              "push": int(i) in push_set, "addendum_sha256": sha, "script_version": SCRIPT_VERSION, "sec": r["sec"],
              "device": str(ctx.device), "smoke": bool(a.smoke)}
        atomic_savez(cdir / f"cell_{cid}.npz", meta=np.array(json.dumps(mc)), **r["arrays"])
        n_done += 1
        rate = (time.time() - t_start) / n_done
        say(f"  {phase} cell {jj + 1}/{len(todo)} ({classes[i]}, {toks[jj].size} tokens, g {r['g_present']}, {r['n_hist']} "
            f"histories, f64 {int(i) in f64_set}) {r['sec']:.1f}s; mean {rate:.1f}s/cell, ETA {rate * (len(todo) - n_done) / 60:.1f} min")
        write_json(out / "e7_progress.json", {"phase": phase, "done_this_call": n_done, "to_do_this_call": len(todo),
                                              "sec_per_cell": round(rate, 2), "eta_min": round(rate * (len(todo) - n_done) / 60, 1),
                                              "updated": time.strftime("%Y-%m-%d %H:%M:%S")})
        if ctx.device.type == "mps":
            ctx.torch.mps.empty_cache()
    say(f"{phase} done: {n_done} cells in {(time.time() - t_start) / 60:.1f} min")
    return 0


def stage_develop(a, add: dict, reg: dict) -> int:
    phase = a.roster
    meta = lp.load_meta(a.processed, need_rna=True)
    emb = lp.Embedding(a.embed_dir, None, len(meta["split"]))
    cells, classes, order = phase_cells(a, add, meta, phase)
    split = np.asarray(meta["split"]).astype(str)
    if a.cell_pool != "train" and not np.all(split[cells] == "val"):
        raise SystemExit("development cells must be val donor 18303 cells")
    donors = mld.load_donors(a.processed)
    return run_phase(a, add, phase, cells, classes, order, list(donors[cells]), meta, emb, meta, LADDER)


def stage_confirm(a, add: dict, reg: dict, sel: dict) -> int:
    meta_b = lp.load_meta(a.processed, need_rna=True)
    mags = tuple(LADDER) + tuple(FRESH)
    if a.smoke:
        cells, classes, order = phase_cells(a, add, meta_b, "smoke")
        emb = lp.Embedding(a.embed_dir, None, len(meta_b["split"]))
        donors = mld.load_donors(a.processed)
        return run_phase(a, add, "smoke_family", cells, classes, order, list(donors[cells]), meta_b, emb, meta_b, mags)
    from teddy_mm.teddy_encoder import load_pad_id, load_vocab
    vocab = load_vocab(a.ckpt)
    pad = load_pad_id(a.ckpt, vocab)
    gtok = {int(v["vocab_id"]) for v in add["computed_train_val"]["G"].values()}
    if a.family == "site4":
        emb = lp.Embedding(a.embed_dir, None, len(meta_b["split"]))
        R = family_roster(a, add, "site4", meta_b, emb, pad, gtok, reg)
        meta = meta_b
    else:
        R = family_roster(a, add, "external", None, None, pad, gtok, reg)
        meta, emb = ext_meta(a, add), ext_embedding(a)
    return run_phase(a, add, a.family, np.asarray(R["cells"], np.int64), R["class"], np.asarray(R["order"], np.int64),
                     R["donor"], meta, emb, meta_b, mags)


# ============================================================================ select / report (CPU)
def load_cell_files(cdir: Path, want: set, sha: str) -> list:
    rows = []
    for f in sorted(cdir.glob("cell_*.npz")):
        with np.load(f, allow_pickle=False) as d:
            mc = json.loads(str(d["meta"]))
            if mc["cell"] not in want:
                continue
            if mc.get("addendum_sha256") != sha:
                raise SystemExit(f"{f} was written under another addendum")
            rows.append((mc, {k: d[k] for k in d.files if k != "meta"}))
    return rows


def units_of(add: dict) -> dict:
    u = add["computed_train_val"]["units"]
    return {"panel": np.asarray(u["panel"], np.float64), "score": float(u["score"])}


def cell_gates(mc: dict, arr: dict, units: dict) -> dict:
    u = units["panel"]
    O0 = arr["O0"]
    km = arr["kind"] != v7.KIND_CLAMP
    c = mc["controls"]
    g = {"padding_patch_gap_units": float(v7.deltas(arr["O_pad1"], arr["O_pad0"], u)[0]),
         "layer12_mean_preserving_gap_units": float(v7.deltas(arr["O_l12"], O0, u)[0]),
         "permutation_gap_units": float(v7.deltas(arr["O_perm"], O0, u)[0]),
         "mean_preservation_abs": float(arr["dmean"][km].max()) if km.any() else 0.0,
         "module_vs_explicit_abs": c["module_vs_explicit_abs"], "padded_vs_unpadded_abs": c["padded_vs_unpadded_abs"],
         "official_layer11_gene_mean_abs": c["official_layer11_gene_mean_abs"]}
    ok = all(v is None or v <= GATES[k] for k, v in g.items())
    return {"values": g, "ok": bool(ok)}


def f64_gap(arr: dict, observable: str, units: dict) -> float | None:
    if "O64" not in arr:
        return None
    d32 = v7.observable_deltas(arr["O"], arr["O0"], observable, units)
    d64 = v7.observable_deltas(arr["O64"], arr["O0_64"], observable, units)
    return float(np.max(np.abs(d32 - d64))) if d32.size else 0.0


def summarise_cells(rows: list, observable: str, units: dict, fracs) -> list:
    codes = [v7.frac_code(f) for f in fracs]

    def dfun(Oa, Ob):
        return v7.observable_deltas(Oa, Ob, observable, units)
    out = []
    for mc, arr in rows:
        s = v7.cell_summary(arr["O"], arr["O0"], arr["kind"], arr["code"], arr["base"], codes, dfun)
        s.update({"cell": mc["cell"], "class": mc["class"], "donor": mc["donor"], "g_present": mc["g_present"],
                  "ntokens": mc["ntokens"]})
        out.append(s)
    return out


def dev_digest(rows: list) -> str:
    h = hashlib.sha256()
    for mc, arr in sorted(rows, key=lambda r: r[0]["cell"]):
        h.update(str(mc["cell"]).encode())
        for k in sorted(arr):
            h.update(k.encode())
            h.update(np.ascontiguousarray(arr[k]).tobytes())
    return h.hexdigest()


def select_phase(a, add: dict, phase: str, observable: str) -> dict:
    meta = lp.load_meta(a.processed, need_rna=False) if a.cell_pool == "train" else None
    cells, classes, order = phase_cells(a, add, meta, phase)
    sha = sha256_file(load_addendum(a)[1])
    cdir = a.out_dir / (phase if not a.smoke else f"smoke_{phase}") / "cells"
    rows = load_cell_files(cdir, {int(c) for c in cells}, sha)
    if len(rows) != len(cells):
        raise SystemExit(f"{phase}: {len(rows)} of {len(cells)} cells done; finish develop first")
    units = units_of(add)
    gates = [cell_gates(mc, arr, units) for mc, arr in rows]
    impl_ok = all(g["ok"] for g in gates)
    f64 = [x for x in (f64_gap(arr, observable, units) for _, arr in rows) if x is not None]
    f64_ok = bool(f64) and max(f64) <= F64_TOL
    # the declared selection rule over the whole development class (all three magnitudes)
    S1 = summarise_cells(rows, observable, units, LADDER)
    D = {c: np.array([s[f"D_{c}"] for s in S1]) for c in v7.CANDIDATES}
    rule = v7.dev_selection(D, TOL)
    stopped = None
    if not impl_ok:
        stopped = "implementation check failed on a development cell"
    elif not f64_ok:
        stopped = "float64 re-check differs by more than tol / 10 (or was not run)"
    selected = rule["selected"] if stopped is None else None
    # descriptive record of every magnitude (never decides)
    per_rung = {}
    for f in LADDER:
        S = summarise_cells(rows, observable, units, [f])
        per_rung[str(v7.frac_code(f))] = {
            "fraction": float(f),
            "q95": {c: v7.q95([s[f"D_{c}"] for s in S]) for c in v7.CANDIDATES},
            "median": {k: v7.finite_median([s[k] for s in S]) for k in ("D_gene_mean", "D_gene_mean_plus_G", "W",
                                                                        "W_context", "W_pure")},
            "literal_max": {c: max([s[f"D_{c}"] for s in S if np.isfinite(s[f"D_{c}"])], default=float("nan"))
                            for c in v7.CANDIDATES[:2]}}
    return {"phase": phase, "observable": observable, "n_cells": len(rows), "magnitudes": [float(f) for f in LADDER],
            "positive_control_median_P": v7.finite_median([s["P"] for s in S1]),
            "G_witness_median_W": v7.finite_median([s["W"] for s in S1]), "per_magnitude": per_rung, "rule": rule,
            "selected": selected, "stopped": stopped, "implementation_ok": impl_ok,
            "implementation_max": {k: max((g["values"][k] for g in gates if g["values"][k] is not None), default=None)
                                   for k in GATES},
            "float64_max_gap": max(f64) if f64 else None, "float64_ok": f64_ok,
            "dims": {"gene_mean": 512, "gene_mean_plus_G_median": 512 * (1 + int(np.median([s["g_present"] for s in S1]))),
                     "all_token_states_median": 512 * int(np.median([s["ntokens"] for s in S1]))},
            "digest_sha256": dev_digest(rows), "cells_sha256": vk.index_hash(cells)}


def stage_select(a, add: dict) -> None:
    sha = sha256_file(load_addendum(a)[1])
    sp = selection_path(a)
    if a.roster == "dev":
        prim = select_phase(a, add, "dev", "panel")
        stop = prim["stopped"] is not None
        rec = {"experiment": "E7", "record_version": 1, "addendum_sha256": sha, "primary": prim,
               "follow_up_required": prim["selected"] is None and not stop, "final": prim["selected"] is not None or stop,
               "selected": prim["selected"], "selected_in": "primary" if prim["selected"] else None,
               "observable": "panel"}
    else:
        if not sp.exists():
            raise SystemExit("no primary selection record")
        rec = json.loads(sp.read_text())
        if rec["primary"]["selected"] is not None or not rec.get("follow_up_required"):
            raise SystemExit("the follow-up runs only after a null primary selection")
        if rec.get("followup") is not None:
            raise SystemExit("the single declared follow-up was already run")
        fu = select_phase(a, add, "followup", "score")
        rec.update({"record_version": 2, "followup": fu, "final": True, "selected": fu["selected"],
                    "selected_in": "followup" if fu["selected"] else None,
                    "observable": "score" if fu["selected"] else "panel",
                    "label": "follow-up after a null primary selection; readout restricted to the NK-T score",
                    "primary_record_version_1_sha256": sha256_file(sp)})
    rec.update({"smoke": bool(a.smoke), "cell_pool": a.cell_pool, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "script_version": SCRIPT_VERSION, "git_head": git_head()})
    sp.parent.mkdir(parents=True, exist_ok=True)
    sp.write_bytes((json.dumps(rec, indent=1, sort_keys=True, default=_json_default) + "\n").encode())
    ssha = sha256_file(sp)
    if not a.smoke:
        hh = addenda_dir(a) / "HASHES.txt"
        lines = [ln for ln in hh.read_text().splitlines() if not ln.rstrip().endswith(f"  {SELECTION}")]
        hh.write_text("\n".join(lines + [f"{ssha}  {SELECTION}"]) + "\n")
        write_json(a.out_dir / f"E7_selection_{a.roster}_copy.json", rec)
    ph = rec.get("followup") or rec["primary"]
    say(f"select ({a.roster}): selected {rec['selected']!r}, stopped {ph['stopped']!r}, follow-up required "
        f"{rec['follow_up_required'] and rec.get('followup') is None}; wrote {sp} (sha256 {ssha})"
        + ("" if a.smoke else "; commit it (and HASHES.txt) before any site4 / external forward"))


def family_block(a, add: dict, sel: dict, family: str, observable: str, rows: list, primary: bool) -> dict:
    units = units_of(add)
    mags = tuple(LADDER) + tuple(FRESH)
    S = summarise_cells(rows, observable, units, mags)
    n = len(S)
    donor = np.array([s["donor"] for s in S]).astype(str)
    cls = np.array([s["class"] for s in S]).astype(str)
    ix = np.arange(n)
    seed = int(add["constants"]["bootstrap_seed"])
    draws = [ix[d_] for d_ in v7.two_stage_draws(donor, cls, a.n_boot, seed)]
    col = {k: np.array([s[k] for s in S], np.float64) for k in ("D_gene_mean", "D_gene_mean_plus_G", "W", "P",
                                                                  "W_context", "W_pure")}

    def med(x):
        return lambda ii: v7.finite_median(x[ii])

    def q(x):
        return lambda ii: v7.q95(x[ii])
    B = {}
    for k in ("D_gene_mean", "D_gene_mean_plus_G", "W", "P"):
        B[f"median_{k}"] = v7.bounds(med(col[k]), ix, donor, cls, a.n_boot, seed, draws)
    # the per-cell difference in W - D_gene_mean_plus_G is reported (G excess) below
    for k in ("D_gene_mean", "D_gene_mean_plus_G"):
        B[f"q95_{k}"] = v7.bounds(q(col[k]), ix, donor, cls, a.n_boot, seed, draws)
    per_donor = primary
    status = {c: v7.candidate_status(B[f"median_D_{c}"], B[f"q95_D_{c}"], TOL, per_donor) for c in v7.CANDIDATES[:2]}
    status["all_token_states"] = "uninformative"
    pc_ok = v7.detects(B["median_P"], TOL, per_donor)
    sep_ok = v7.detects(B["median_W"], TOL, per_donor)
    gates = [cell_gates(mc, arr, units) for mc, arr in rows]
    impl_ok = all(g["ok"] for g in gates)
    f64 = [x for x in (f64_gap(arr, observable, units) for _, arr in rows) if x is not None]
    f64_ok = bool(f64) and max(f64) <= F64_TOL
    V = v7.verdict(sel.get("selected"), status, pc_ok, sep_ok, impl_ok, f64_ok)
    other = "score" if observable == "panel" else "panel"
    S_o = summarise_cells(rows, other, units, mags)
    sec = {"other_observable": other,
           "other_medians": {k: v7.finite_median([s[k] for s in S_o]) for k in ("D_gene_mean", "D_gene_mean_plus_G", "W", "P")},
           "other_q95": {k: v7.q95([s[k] for s in S_o]) for k in ("D_gene_mean", "D_gene_mean_plus_G")},
           "per_magnitude": {}, "per_class": {}, "g_excess_median_W_minus_D2": v7.finite_median(col["W"] - col["D_gene_mean_plus_G"])}
    for f in mags:
        Sm = summarise_cells(rows, observable, units, [f])
        kf = str(v7.frac_code(f))
        sec["per_magnitude"][kf] = {k: v7.finite_median([s[k] for s in Sm]) for k in ("D_gene_mean", "D_gene_mean_plus_G",
                                                                                      "W", "W_context", "W_pure", "P")}
        sec["per_magnitude"][kf].update({f"q95_{k}": v7.q95([s[k] for s in Sm]) for k in ("D_gene_mean", "D_gene_mean_plus_G")})
    for c in ("NK", "T"):
        mk = cls == c
        sec["per_class"][c] = {k: v7.finite_median(col[k][mk]) for k in col} | {"n": int(mk.sum())}
    pu = [v7.observable_deltas(arr["O_push"], arr["O0"], observable, units) for _, arr in rows if "O_push" in arr]
    sec["input_push_on_G"] = {"n_cells": len(pu), "median_of_cell_median": v7.finite_median([np.median(x) for x in pu]) if pu else None}
    u = units_of(add)
    up = np.asarray(add["computed_train_val"]["units"]["predicted_panel_secondary"], np.float64)
    sec["predicted_unit_reexpression_factor_per_protein"] = (u["panel"] / up).tolist()
    return {"family": family, "observable": observable, "n_cells": n, "magnitudes": list(mags),
            "by_donor_class": {f"{d_}/{c}": int(np.sum((donor == d_) & (cls == c))) for d_ in sorted(set(donor)) for c in ("NK", "T")},
            "bounds": B, "status": status, "positive_control_detected": pc_ok, "G_witness_detected": sep_ok,
            "implementation_ok": impl_ok, "implementation_max": {k: max((g["values"][k] for g in gates if g["values"][k] is not None), default=None) for k in GATES},
            "float64_max_gap": max(f64) if f64 else None, "float64_ok": f64_ok, "verdict": V, "secondary": sec,
            "per_cell": [{k: s[k] for k in ("cell", "class", "donor", "g_present", "ntokens", "D_gene_mean",
                                           "D_gene_mean_plus_G", "W", "P")} for s in S]}


def stage_report(a, add: dict, sel: dict, reg_info: dict) -> None:
    sha = sha256_file(load_addendum(a)[1])
    R = {"experiment": "E7", "script_version": SCRIPT_VERSION, "smoke_test": bool(a.smoke), "smoke_note": a.smoke,
         "created": time.strftime("%Y-%m-%d %H:%M:%S"), "addendum_sha256": sha, "registration_check": reg_info,
         "selection": {k: sel.get(k) for k in ("selected", "selected_in", "observable", "label", "final")},
         "families": {}}
    fams = [("smoke_family", True)] if a.smoke else [("site4", True), ("external", False)]
    for fam, primary in fams:
        cdir = a.out_dir / (fam if not a.smoke else "smoke_smoke_family") / "cells"
        if a.smoke:
            meta = lp.load_meta(a.processed, need_rna=False) if a.cell_pool == "train" else None
            cells = phase_cells(a, add, meta, "smoke")[0]
        else:
            rf = a.out_dir / fam / "e7_roster.json"
            if not rf.exists():
                R["families"][fam] = {"status": "not run"}
                continue
            cells = json.loads(rf.read_text())["cells"]
        rows = load_cell_files(cdir, {int(c) for c in cells}, sha)
        if not rows:
            R["families"][fam] = {"status": "no cells done"}
            continue
        blk = family_block(a, add, sel, fam, sel.get("observable") or "panel", rows, primary)
        blk["complete"] = len(rows) == len(cells)
        blk["n_selected_cells"] = len(cells)
        if fam == "external":
            blk["labels"] = add["external_family"]["labels"]
            blk["note"] = "secondary family; never changes the registered verdict"
        if not blk["complete"]:
            blk["verdict"]["verdict"] = f"INTERIM ({len(rows)} of {len(cells)} cells) " + blk["verdict"]["verdict"]
        if sel.get("selected_in") == "followup":
            blk["verdict"]["label"] = sel.get("label")
        R["families"][fam] = blk
    prim = R["families"].get("smoke_family" if a.smoke else "site4", {})
    R["registered_verdict"] = prim.get("verdict") if isinstance(prim, dict) else None
    R["final"] = bool(not a.smoke and isinstance(prim, dict) and prim.get("complete"))
    a.out_dir.mkdir(parents=True, exist_ok=True)
    write_json(a.out_dir / "E7_results.json", R)
    write_report(a.out_dir, R, add)
    say(f"report: verdict {R['registered_verdict']}")


def _f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


def write_report(out: Path, R: dict, add: dict) -> None:
    L = ["# E7: TEDDY analog of ANM's inverse loop (layer-11 cut, consumer fixed)", ""]
    if R["smoke_test"]:
        L += [f"SMOKE TEST: {R['smoke_note']} (never a result)", ""]
    s = R["selection"]
    L += [f"{R['script_version']}; addendum {R['addendum_sha256'][:16]}; selection {s.get('selected')!r} "
          f"(in {s.get('selected_in')}), observable {s.get('observable')}.", ""]
    for fam, b in R["families"].items():
        L += [f"## {fam}", ""]
        if "bounds" not in b:
            L += [f"{b.get('status')}", ""]
            continue
        if b.get("labels"):
            L += ["Labels: " + "; ".join(b["labels"]), ""]
        L += [f"**Verdict:** {b['verdict']['code']}: {b['verdict']['verdict']}", "",
              f"cells {b['n_cells']} of {b['n_selected_cells']}; magnitudes {b['magnitudes']}; by donor/class {b['by_donor_class']}",
              "", "| statistic | value | 95% interval | per donor (value [interval]) |", "|---|---|---|---|"]
        for k, v in b["bounds"].items():
            ci = v["ci"]
            L.append(f"| {k} | {_f(v['value'])} | " + (f"[{_f(ci[0])}, {_f(ci[1])}]" if ci else "n/a") + " | "
                     + "; ".join(f"{d_} {_f(x['value'])} " + (f"[{_f(x['ci'][0])}, {_f(x['ci'][1])}]" if x["ci"] else "")
                                 for d_, x in v["per_donor"].items()) + " |")
        L += ["", f"status {b['status']}; positive control detected {b['positive_control_detected']}; G witness detected "
                  f"{b['G_witness_detected']}; implementation ok {b['implementation_ok']}; float64 max gap "
                  f"{_f(b['float64_max_gap'], 6)} (ok {b['float64_ok']})", "",
              f"implementation maxima `{json.dumps(b['implementation_max'])}`", "",
              f"secondary `{json.dumps(b['secondary'], default=_json_default)[:3000]}`", ""]
    L += [f"tol {TOL} in units of the measured val NK-T gap per protein ({add['computed_train_val']['units']['panel']}); "
          f"score unit {add['computed_train_val']['units']['score']:.4f}."]
    (out / "REPORT.md").write_text("\n".join(L) + "\n")


# ============================================================================ main
def main(argv=None) -> int:
    a = parse_args(argv)
    if a.smoke:
        a.out_dir.mkdir(parents=True, exist_ok=True)
    if a.stage == "register":
        if a.smoke:
            _LOG["path"] = a.out_dir / "progress.log"
        stage_register(a)
        return 0
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG["path"] = a.out_dir / "progress.log"
    say(f"{SCRIPT_VERSION} stage {a.stage} roster {a.roster} family {a.family} out {a.out_dir}"
        + (f" SMOKE: {a.smoke} ({a.cell_pool})" if a.smoke else ""))
    need = {"develop": "followup" if a.roster == "followup" else "develop", "select": "select",
            "confirm": "confirm", "report": "report"}[a.stage]
    info = check_registration(a, need)
    write_json(a.out_dir / f"e7_registration_used_{a.stage}.json", info)
    add, _ = load_addendum(a)
    reg = load_reg(a)
    if a.stage == "develop":
        return stage_develop(a, add, reg)
    if a.stage == "select":
        stage_select(a, add)
        return 0
    sp = selection_path(a)
    if not sp.exists():
        raise SystemExit("no development selection record: run develop and select first")
    sel = json.loads(sp.read_text())
    if not sel.get("final"):
        raise SystemExit("the selection record is not final (the declared follow-up is required first)")
    if a.stage == "confirm":
        return stage_confirm(a, add, reg, sel)
    stage_report(a, add, sel, info)
    return 0


if __name__ == "__main__":
    sys.exit(main())
