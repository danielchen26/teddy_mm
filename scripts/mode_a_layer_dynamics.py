#!/usr/bin/env python3
"""E5 (C7), Mode A layer dynamics: follow a gene-token perturbation through TEDDY-G's 12 layers and
find the layer at which its NK-T component is attenuated or amplified.

The static probes (mode_a_layer_probes.py, mode_a_nonlinear_probes.py, mode_a_token_probes.py) read the
NK-T difference off each depth. This script treats the layers as steps of a dynamical system and asks how
a small push on one gene token propagates:

    h_0 = token embedding + position embedding  (official TEDDY-G input, no CLS)
    h_l = f_l(h_{l-1}),  l = 1..12               (the checkpoint's own encoder layers, post-norm)
    dh_l = J_l dh_{l-1},  J_l = df_l/dh at h_{l-1}

dh_l is the exact linear response (Jacobian-vector product), computed with torch.func.jvp through each
layer in turn. For a post-norm layer f_l(h) = LN2(u + FF(u)), u = LN1(h + SA(h)), so
J_l = dLN2 (I + dFF) dLN1 (I + dSA): the "(I + df/dh)" residual step sits inside the two LayerNorms.
PyTorch's fused attention kernels (CPU fast path, MPS scaled-dot-product kernel) have no forward-mode
derivative, so every layer is evaluated with its own weights and explicit attention arithmetic
(softmax(QK^T / sqrt(d_h)) V), and that forward is checked against the module's forward for every cell.

Cells: site4 cells of the two primary donors (not in training): look-alike NK-T pair cells and a matched
random set (REGISTRATION["cells"]). Per cell (official preprocessing; the cell's real tokens only, so no
padding), per declared gene present among its tokens, per direction v (unit vector), the input state at
that gene's position is pushed by  eps * ||E[gene]|| * v :
  nkt    v = unit(mu_NK - mu_T) of the input gene-mean, training NK and T cells
  self   v = unit(E[gene]) (more or less of the gene's own embedding)
  rand   v = one fixed Gaussian unit vector (control)
  other  2 random non-declared tokens of the cell, each along its own embedding (control)
For each case the script records the JVP response at every depth l = 0..12 (gene-mean pooled dz_l and the
response at the pushed token itself), the head's first-order response (JVP through L2 normalisation and the
Phase-1 head), and finite differences at +-eps and +-eps/2 for a ladder of eps (REGISTRATION["eps"]) as
the linearity check (one-sided slopes, so a quadratic part of the response is not cancelled).

Read-outs per depth l (report stage), all in the unit A0 = ||E[gene]|| / (n ||dmu_0||) (the NK-T component
an aligned input push of the same size has at the input; n = the cell's token count):
  G_l   NK-T component of dz_l along dmu_l = mu_NK,l - mu_T,l (training cells), as a fraction of
        ||dmu_l||^2, over A0. G_0 = 1 for nkt. Gain of the NK-T component; dG_l = G_l - G_{l-1} is the
        share of the injected push gained (> 0) or lost (< 0) at layer l, the localisation statistic.
  P_l   the same along the layer-l logistic NK-vs-T probe (mode_a_layer_probes settings; C on val).
  H_l   head lens: g_head . dz_l over the training NK-T gap of the head score, over A0; g_head = gradient of
        the head score at the cell's z. H_12 is the head's actual first-order response.
  own / spread: the part of G_l carried by the pushed token itself vs by all other tokens (attention).
  linearity pass rate: share of cases whose one-sided FD slopes at eps match the JVP and the slope at eps/2.
  JVP check (registration_v3 experiments.E5): central FD on the unit tangent at eps_abs = 1e-3 and 1e-2.

Every threshold, cell rule, unit and hypothesis is in REGISTRATION below, an addendum to
registration/registration_v3.json (experiments.E5; the deviations from that entry are listed in it).
--stage register writes registration/addenda/E5.json and its sha256 line in registration/addenda/HASHES.txt;
a site4 run is refused unless both, and registration_v3.json, are committed and match. Smoke runs use
val-donor cells (--cell-pool val), never site4.

Stages
  register    write registration/addenda/E5.json + its line in addenda/HASHES.txt
  select      rebuild and check the site4 look-alike pairs, choose the cells (no forward pass)
  directions  CPU: NK-T directions per depth, probe directions, head-score gaps (training + val cells)
  respond     GPU/MPS (CPU works, slowly): per-cell JVP + finite differences, one file per cell (resumable)
  report      CPU: aggregate, two-stage bootstrap, compare with outputs/mode_a_official; E5_results.json + REPORT.md
  all         select, directions, respond, report

Full run (from the repo root; not run inside the workflow). Paths to the main checkout's data and outputs:
  PY=<venv>/bin/python; R=/Users/tianchichen/Documents/GitHub/teddy_mm
  PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_layer_dynamics.py --stage all \
      --processed $R/data/processed/cite --embed-dir $R/data/processed/cite_official \
      --ckpt /Users/tianchichen/Documents/GitHub/teddy_mwe/ckpt/teddy_g_70M \
      --medians $R/data/reference/teddy_gene_medians.json \
      --head-ckpt $R/outputs/cite_phase1_official/best.pt --probe-dir $R/outputs/mode_a_official \
      --out-dir $R/outputs/v3/E5 --device auto --threads 4
Every stage appends to <out-dir>/progress.log; the respond stage logs seconds per cell and an ETA after every
cell (also in e5_progress.json); rerunning the
same command skips finished cells (each cell file carries the registration hash and is refused under another).
Smoke runs: --cell-pool val --smoke NOTE [--limit N --max-genes N --dir-max-train N] (val donor cells only).
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
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location("mode_a_layer_probes", HERE / "mode_a_layer_probes.py")
lp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lp)

SCRIPT_VERSION = "mode_a_layer_dynamics v2"
PROTEINS, GDT158 = lp.PROTEINS, lp.GDT158  # CD56, CD94, CD335, CD3
rnd = lp.rnd
N_DEPTHS = 13  # 0 = input state, 1..12 = layer outputs
DIRECTIONS = ("nkt", "self", "rand", "other")  # "other" = self push of a random non-declared token (control)
# symbol -> (Ensembl GRCh38, class); the same ids as mode_a_token_probes.py
GENES = {
    "CD3E": ("ENSG00000198851", "T"), "CD3D": ("ENSG00000167286", "T"),
    "NCAM1": ("ENSG00000149294", "NK"), "KLRD1": ("ENSG00000134539", "NK"), "NCR1": ("ENSG00000189430", "NK"),
    "FCGR3A": ("ENSG00000203747", "NK"), "KLRF1": ("ENSG00000150045", "NK"),
}
GENE_SYMBOLS = tuple(GENES)
PRIMARY_DONORS = ("13272", "19593")  # site4 donors not in training (registration_v3 splits.test_primary)
SEEDS = {"cells": 19, "random_directions": 23, "bootstrap": 1}  # registration_v3 seeds e5_cells, e5_random_directions, bootstrap
N_BOOT = 2000
EPS_ABS_CHECK = (1e-3, 1e-2)
N_OTHER = 2
N_FULL_LADDER = 100  # the first cells of the seeded processing order get every eps rung
N_FD64_CELLS = 20  # the first cells of the order get the float64 CPU finite-difference JVP check (one case per direction)
_LOG = {"path": None}


def say(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG["path"] is not None:
        with open(_LOG["path"], "a") as f:
            f.write(line + "\n")


REGISTRATION = {
    "experiment": "E5",
    "title": "Mode A layer dynamics (claim C7): Jacobian-vector products through TEDDY-G's layers",
    "addendum_to": "registration/registration_v3.json, experiments.E5. This addendum fixes what experiments_v3.json leaves "
                   "open and lists where this build deviates from it (deviations_from_experiments_v3); the orchestrator "
                   "reconciles the two before the site4 run",
    "registration_version": 2,
    "script": "scripts/mode_a_layer_dynamics.py (" + SCRIPT_VERSION + ")",
    "fixed_before_site4": "every entry is fixed on training and validation cells only and committed before any site4 "
                          "forward pass of this script; smoke runs use val donor 18303 cells only",
    "model": "official TEDDY-G 70M (teddy_mwe/ckpt/teddy_g_70M), official preprocessing (teddy_mm/teddy_encoder.py: "
             "counts/total*1e4 / gene median, torch.topk 2048, zeros dropped, no CLS), float32, the cell's real tokens only "
             "(no padding); each layer with its own weights and explicit attention arithmetic, checked against the module "
             "forward for every cell",
    "head": "Phase-1 head outputs/cite_phase1_official/best.pt: dec(mlp(z / ||z||), train-median size factor), z = gene-mean "
            "of layer 12 (as 04_train.py and mode_a_layer_probes.py)",
    "cells": {
        "split": "test_primary: site4 cells of the donors not in training (13272, 19593); donor 15078 (a training donor) is "
                 "not used",
        "primary_donors": list(PRIMARY_DONORS),
        "pair_cells": "per primary donor: the donor's cells in the no_gdT158 NK-T look-alike pairs of "
                      "mode_a_layer_probes.py (site4 k = 10 cosine neighbours on raw z, fixed annotation key; rebuilt and "
                      "checked against outputs/mode_a_official/mode_a_pair_gaps.npz); pairs taken in one seeded order "
                      "(default_rng(19)), their cells of that donor added until n_pair_cells_per_donor",
        "n_pair_cells_per_donor": 150,
        "random_cells": "per primary donor and class: as many of the donor's site4 NK / T cells as among its chosen pair "
                        "cells, no gdT CD158b+, in no look-alike pair ('all' variant), stratified to those pair cells' "
                        "token-count quintiles (same generator)",
        "max_cells": 600,
        "seed": SEEDS["cells"],
    },
    "genes": {s: {"ensembl": e, "class": c} for s, (e, c) in GENES.items()},
    "genes_rule": "every declared gene among the cell's tokens gets an nkt and a self case; each cell gets one rand case "
                  "(at one of its declared genes drawn with default_rng([23, cell index]), else at a random token) and "
                  "n_other_genes 'other' cases",
    "perturbation": {
        "site": "input state h_0 = E[id] + P[pos] at the token's position; nothing else changes",
        "tangent": "||E[token]|| * v (v unit); finite differences scale it by eps",
        "directions": {
            "nkt": "unit(mu_NK - mu_T) of the input gene-mean over training NK and T cells (split 'train')",
            "self": "unit(E[gene])",
            "rand": "one Gaussian unit vector, numpy default_rng(23), the same for every case (not orthogonalised)",
            "other": "control: n_other_genes tokens of the cell that are not declared genes, drawn with "
                     "default_rng([23, cell index]), each pushed along its own embedding",
        },
        "n_other_genes": N_OTHER,
    },
    "eps": [0.02, 0.1, 1.0],
    "eps_primary": 0.02,
    "eps_cells": "eps_primary on every cell; 0.1 and 1.0 (and eps_abs 1e-2) on the first n_full_ladder cells of the seeded "
                 "processing order (GPU time)",
    "n_full_ladder": N_FULL_LADDER,
    "eps_choice": "fixed on val donor 18303 cells (smoke calibration, 2 cells, KLRD1/CD3E pushes): the one-sided slope "
                  "error grows linearly with eps (about 0.03 at 0.02, 0.08 at 0.05, 0.16 at 0.1), the central error is near "
                  "the float32 floor (about 1e-3) below 0.02; so 0.02 is the small-signal check and 0.1 and 1.0 measure how "
                  "far the linear description reaches",
    "finite_differences": "forwards through the module's own encoder layers (layer(h), no padding) at h_0 + c T for c in "
                          "{0, +-e, +-e/2} for each eps run on the cell and c = +-eps_abs / ||E[token]|| for the JVP check, "
                          "all in one batch shape; one-sided slopes D+(e) = (z(+e) - z0)/e, D-(e) = (z0 - z(-e))/e, central "
                          "Dc = (D+ + D-)/2. The JVP (explicit arithmetic) is thereby checked against the module itself",
    "response": "JVP: (h_l, dh_l) = torch.func.jvp(f_l, (h_{l-1},), (dh_{l-1},)), l = 1..12; dz_l = mean over the cell's "
                "tokens of dh_l; head response by JVP through L2 normalisation + head",
    "jvp_check": {
        "eps_abs": list(EPS_ABS_CHECK),
        "rule": "experiments_v3 E5 check, run in float64: central FD on the unit tangent at eps_abs = 1e-3 and 1e-2, in "
                "float64 on the CPU through the module's own layers, for one case per direction (the first of each) of "
                "the first n_fd64_cells cells of the seeded order; relative error <= 0.05 and cosine >= 0.99 against the "
                "device float32 JVP at every depth 1..12 for >= 95% of the checked cases at 1e-3 -> validated; otherwise "
                "'JVP not validated' and no E5 endpoint is read",
        "n_fd64_cells": N_FD64_CELLS,
        "float32_diagnostic": "the same central FD in float32 on the device (every case at 1e-3, full-ladder cases at "
                              "1e-2) is reported, not used: at 1e-3 it reaches the float32 rounding floor for weak "
                              "responses at depths 9-12 (val smoke: 4 of 36 cases with relative error 0.05-0.12 and "
                              "cosine >= 0.993 there, every case passing at 1e-2)",
        "cpu64_jvp": "the JVP of the first 2 cases of the first 2 cells is recomputed in float64 on the CPU (reported)",
    },
    "readout_directions": {
        "mu_diff": "dmu_l = mean_NK - mean_T of the depth-l gene-mean over training NK and T cells (depth 0 recomputed with "
                   "mode_a_layer_probes.input_layer_means, depths 1..12 from z_rna_layer_means.npy)",
        "probe": "logistic regression NK vs T on the train-standardised depth (mode_a_layer_probes.nkt_probe settings: C in "
                 "(0.01, 0.1, 1) by val log-loss, class_weight balanced, seed 0); direction = coef / sd",
        "head_gradient": "gradient of s_head at the cell's own float32 z",
    },
    "units": "A0 = ||E[token]|| / (n ||dmu_0||), n = the cell's token count. G_l = (dmu_l . dz_l / ||dmu_l||^2) / A0 (G_0 = 1 "
             "for nkt); dG_l = G_l - G_{l-1}; P_l = (w_l . dz_l / w_l . dmu_l) / A0; H_l = (g_head . dz_l / (mean_trainNK "
             "s_head - mean_trainT s_head)) / A0; self pushes are sign-aligned by gene class (NK +1, T -1)",
    "head_score": "s_head = mean(yhat_CD56, yhat_CD94, yhat_CD335) - yhat_CD3, yhat = predicted ADT / train p95 "
                  "(mode_a_layer_probes.p95_replicate), unclipped",
    "linearity": {
        "tau": 0.1,
        "one_sided_pass": "for both signs: ||D(eps) - JVP|| / ||JVP|| <= tau AND ||D(eps) - D(eps/2)|| / ||D(eps/2)|| <= tau, "
                          "on the pooled dz_l (primary)",
        "central_pass": "the same two conditions for Dc (odd part only; secondary)",
        "scalar_pass": "for both signs |G_D(eps) - G_JVP| <= tau * max(|G_JVP|, 0.1) (NK-T component only; secondary)",
        "curvature": "kappa = ||D+ - D-|| / ||D+ + D-|| (descriptive)",
        "linear_range": "largest eps of the ladder with one-sided pass rate >= 0.90 at every depth 1..12, over the cases "
                        "the eps was run on (descriptive)",
    },
    "statistics": {
        "summary": "median over cases (cell x token x direction)",
        "decision_interval": "two-stage bootstrap: donors with replacement, then cells with replacement within each drawn "
                             "donor (all of a cell's cases kept); B = 2000, seed 1, percentile 95% (experiments_v3 "
                             "common.statistics)",
        "cell_interval": "cell-cluster bootstrap, B = 2000, seed 1 (reported, not used for decisions)",
        "per_donor": "every endpoint is also reported per primary donor; a directional verdict needs the same sign in both",
        "two_donor_caveat": "2 primary donors: the interval reflects the spread of the observed donors only",
        "localisation": "by the per-case step change dG_l (the ratio of medians median G_l / median G_{l-1} is reported but "
                        "not used: it is unstable when a median is near 0, as the val smoke showed)",
    },
    "hypotheses": {
        "H1_linearity": "at eps = 0.02 the one-sided pass rate is >= 0.90 at every depth 1..12 (directions nkt, self, rand, "
                        "other; pair and random cells pooled). If not, H2-H5 are read from the one-sided FD gains D+ at "
                        "eps = 0.02 and marked as such.",
        "H2_localised_attenuation": "nkt, pair cells: L_E5 = argmin_{l=1..12} median dG_l has a two-stage 95% upper bound "
                                    "< 0, is the argmin in >= 50% of the two-stage resamples, and has median dG < 0 in each "
                                    "primary donor",
        "H3_localisation_matches_probes": "|L_E5 - L_probe| <= 1, L_probe = argmin_{l=1..12} of the mean over CD56 / CD94 / "
                                          "CD335 / CD3 of the no_gdT158 gap ratio in outputs/mode_a_official/"
                                          "mode_a_results.json (read at run time; layer 9 in the file at registration)",
        "H4_lookalike_specificity": "D = median G_12(pair) - median G_12(random), nkt: holds if the two-stage upper bound "
                                    "< 0 and D < 0 in each primary donor (look-alike cells damp NK-T pushes more); "
                                    "'opposite' if lower bound > 0 and D > 0 in each donor; otherwise no evidence",
        "H5_head_transmission": "median over nkt pair cases of (H_12 - G_12): two-stage upper bound < 0 and < 0 in each "
                                "donor -> the head passes on less of the NK-T push than z carries; descriptive",
        "E5_v3_endpoint": "experiments_v3 E5 endpoint, computed as a secondary endpoint: per case r_l = <dz_l, unit(dmu_l)>, "
                          "log g_l = log|r_l| - log|r_{l-1}|; self pushes in pair cells: l* = argmin_l median log g_l; "
                          "agreement if within 1 layer of the layer of the largest one-step drop of the mean no_gdT158 "
                          "gap ratio; specificity: median log g_{l*}(self) - median log g_{l*}(rand), two-stage interval "
                          "inside (-0.1, 0.1) -> no layer localises the NK-T loss; <= -0.1 with upper bound < 0 and the "
                          "same sign in both donors -> NK-T-specific attenuation at l*; otherwise inconclusive",
    },
    "falsification": "If the JVP check fails, E5 stops ('JVP not validated'). C7 (the layer dynamics localise where the NK-T "
                     "difference is attenuated) is falsified for this model if H1 holds and neither H2 nor H4 holds: the "
                     "NK-T push then passes every layer without a localised attenuation, so the static probe curves are "
                     "not explained by the dynamics of token perturbations. If H1 fails at eps = 0.02 the "
                     "linear-response reading is inconclusive and the FD reading is reported instead.",
    "comparators": "the static layer probes (outputs/mode_a_official/mode_a_results.json: NK-T probe order rate, no_gdT158 "
                   "gap ratios, centred distance ratio per depth) and the token probes (outputs/mode_a_official/token); "
                   "controls: rand and other pushes. E5 adds information beyond the probes only through H2-H5 and the "
                   "own-token vs spread decomposition, which the static probes cannot give",
    "outputs": "outputs/v3/E5/: E5_results.json, REPORT.md, progress.log, cells/ (one file per cell); every result JSON "
               "carries registration_sha256 (registration_v3.json) and addendum_sha256 (this file)",
    "deviations_from_experiments_v3": [
        "cells: the mode_a no_gdT158 pairs (fixed annotation key of mode_a_layer_probes.py, k = 10 over all site4 cells) "
        "restricted to primary-donor cells (150 per donor, <= 200 as experiments_v3), not the E3 pairs under the v3 key "
        "(the v3 key needs registration_v3.json, not built when this was written); plus a matched random set, as the E5 "
        "build instructions ask",
        "genes: the 7 declared NK/T genes, not every E2 gene (GPU time); 2 'other' tokens per cell as the random-gene control",
        "JVP path: explicit attention arithmetic with the layers' own weights, float32 on the GPU (MPS), not the standard "
        "module path on the CPU: forward-mode AD is not implemented for PyTorch's fused attention kernels "
        "(_transformer_encoder_layer_fwd, _native_multi_head_attention, _scaled_dot_product_flash_attention_for_cpu, "
        "_scaled_dot_product_attention_math_for_mps; torch 2.14); the explicit forward is checked against the module "
        "forward for every cell and the JVP against a float64 CPU JVP",
        "endpoint: the primary localisation statistic is the signed per-case step change dG_l (log|r_l| ignores sign "
        "reversals and is heavy-tailed when r_{l-1} is near 0); the experiments_v3 log-gain endpoint and its margins are "
        "computed as a secondary endpoint",
        "linearity: one-sided slopes on an eps ladder are added, because a central difference cancels the quadratic part "
        "of the response, which the val smoke put at 10-20 % of the linear part at eps = 0.1",
        "JVP check: run in float64 on the CPU for a declared subset (1e-3 in float32 hits the rounding floor for weak "
        "responses, see jvp_check.float32_diagnostic); the float32 version is reported",
    ],
}


def registration_bytes() -> bytes:
    return (json.dumps(REGISTRATION, indent=1, sort_keys=True) + "\n").encode()


def registration_sha() -> str:
    return hashlib.sha256(registration_bytes()).hexdigest()


# ============================================================================ args

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("register", "select", "directions", "respond", "report", "all"), default="all")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=ROOT / "data/processed/cite_official")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=ROOT / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--probe-dir", type=Path, default=ROOT / "outputs/mode_a_official",
                   help="mode_a_layer_probes.py out-dir (pairs check, probe curves)")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration",
                   help="holds registration_v3.json(.sha256) and addenda/E5.json + addenda/HASHES.txt")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--cell-pool", choices=("site4", "val"), default="site4",
                   help="val = val donor 18303 NK/T cells (smoke tests; no site4 forward)")
    p.add_argument("--device", default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--jvp-chunk", type=int, default=4, help="tangents per vmap call")
    p.add_argument("--fd-batch", type=int, default=4, help="perturbed sequences per forward call")
    p.add_argument("--cpu64-check", type=int, default=2,
                   help="for the first N processed cells, recompute the JVP of their first 2 cases in float64 on the CPU")
    p.add_argument("--limit", type=int, default=0, help="smoke only: process just the first N cells of the order")
    p.add_argument("--max-genes", type=int, default=0, help="smoke only: at most N genes per cell")
    p.add_argument("--dir-max-train", type=int, default=0,
                   help="smoke only: at most N training (and N val) cells per class for the directions")
    p.add_argument("--n-boot", type=int, default=N_BOOT, help="registered: 2000; smaller only for smoke runs")
    p.add_argument("--n-full-ladder", type=int, default=None, help="smoke only: override the registered n_full_ladder")
    p.add_argument("--smoke", default=None, metavar="NOTE")
    a = p.parse_args(argv)
    if a.stage != "register" and a.out_dir is None:
        p.error("--out-dir is required")
    if (a.limit or a.max_genes or a.dir_max_train or a.n_boot != N_BOOT or a.n_full_ladder is not None) and not a.smoke:
        p.error("--limit / --max-genes / --dir-max-train / --n-boot / --n-full-ladder are for smoke runs only (give --smoke NOTE)")
    if a.cell_pool == "val" and not a.smoke:
        p.error("--cell-pool val is for smoke runs only")
    return a


# ============================================================================ helpers

def fkey(meta):
    missing = sorted(set(meta["cell_types"]) - set(lp.ANNOTATION_MAP))
    if missing:
        raise SystemExit(f"cell types not in the fixed key: {missing[:10]}")
    return np.array([lp.ANNOTATION_MAP[c] or "" for c in meta["cell_types"]])


def load_donors(processed: Path) -> np.ndarray:
    with np.load(processed / "cite_arrays.npz", allow_pickle=True) as d:
        return d["donors"].astype(str)


def atomic_savez(path: Path, **arrays):
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def write_json(path: Path, obj):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1))
    os.replace(tmp, path)


# ============================================================================ register

def stage_register(a):
    d = a.registration_dir / "addenda"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "E5.json"
    f.write_bytes(registration_bytes())
    sha = registration_sha()
    hf = d / "HASHES.txt"
    lines = [ln for ln in (hf.read_text().splitlines() if hf.exists() else []) if not ln.rstrip().endswith("  E5.json")]
    hf.write_text("\n".join(lines + [f"{sha}  E5.json"]) + "\n")
    say(f"wrote {f} and its line in {hf} (sha256 {sha})")


def _git_committed(path: Path) -> bool | None:
    import subprocess
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


def check_registration(a) -> dict:
    sha = registration_sha()
    f = a.registration_dir / "addenda" / "E5.json"
    hf = a.registration_dir / "addenda" / "HASHES.txt"
    info = {"addendum_sha256": sha, "addendum_file": str(f),
            "addendum_matches_script": f.exists() and hashlib.sha256(f.read_bytes()).hexdigest() == sha,
            "addendum_hash_recorded": hf.exists() and f"{sha}  E5.json" in hf.read_text().splitlines(),
            "addendum_committed": _git_committed(f) if f.exists() else False}
    cf = a.registration_dir / "registration_v3.json"
    ch = a.registration_dir / "registration_v3.json.sha256"
    info["registration_sha256"] = None
    if cf.exists():
        h = hashlib.sha256(cf.read_bytes()).hexdigest()
        info["registration_sha256"] = h
        info["registration_hash_file_matches"] = ch.exists() and ch.read_text().split()[0] == h
        info["registration_committed"] = _git_committed(cf)
    if a.cell_pool == "site4" and a.stage != "select":  # select draws cells only (no forward pass, no E5 outcome)
        need = {"addendum file equals the script's REGISTRATION": info["addendum_matches_script"],
                "addendum sha256 recorded in addenda/HASHES.txt": info["addendum_hash_recorded"],
                "addendum committed": info["addendum_committed"] is True,
                "registration_v3.json present, hash file matching, committed":
                    bool(info["registration_sha256"]) and info.get("registration_hash_file_matches") is True
                    and info.get("registration_committed") is True}
        bad = [k for k, ok in need.items() if not ok]
        if bad:
            raise SystemExit("site4 run refused (register and commit before any site4 forward): " + "; ".join(bad))
    return info


# ============================================================================ cells

def stage_select(a, meta, emb, key) -> dict:
    out = a.out_dir
    split, ctype = meta["split"], meta["cell_types"]
    donors = load_donors(a.processed)
    if emb.partial:
        raise SystemExit("needs the full official embedding")
    ntok = emb.ntokens
    C = REGISTRATION["cells"]
    rng = np.random.default_rng(C["seed"])
    if a.cell_pool == "val":
        sel, role = [], []
        for cls in ("NK", "T"):
            pool = np.where((split == "val") & (key == cls))[0]
            m = min(pool.size, C["max_cells"] // 2)
            sel.append(np.sort(rng.choice(pool, m, replace=False)))
            role += ["val"] * m
        sel = np.concatenate(sel)
        role = np.array(role)
        info = {"pool": "val donor 18303 NK/T cells (smoke)"}
    else:
        LR = json.loads((a.probe_dir / "mode_a_results.json").read_text())
        cfg = LR["config"]
        ev = np.where(split == "test")[0]
        z = np.asarray(emb.z, dtype=np.float32)
        t0 = time.time()
        pairs = lp.build_pairs(z[ev], ev, key, ctype, int(cfg["k"]), int(cfg["n_random_per_seed"]), list(cfg["random_seeds"]))
        nk, t = pairs["neighbour"]["no_gdT158"]
        with np.load(a.probe_dir / "mode_a_pair_gaps.npz") as G:
            if not (np.array_equal(G["nk_no_gdT158"], nk) and np.array_equal(G["t_no_gdT158"], t)):
                raise SystemExit("rebuilt no_gdT158 pairs differ from the linear run's")
        say(f"pairs rebuilt and checked ({time.time() - t0:.0f}s): {nk.size} no_gdT158 pairs, "
            f"{np.unique(np.concatenate([nk, t])).size} cells")
        nk_all, t_all = pairs["neighbour"]["all"]
        in_any = np.zeros(len(split), bool)
        in_any[nk_all] = True
        in_any[t_all] = True
        want = C["n_pair_cells_per_donor"]
        order_p = rng.permutation(nk.size)
        pair_cells, rand, pairs_used = [], [], {}
        for dn in C["primary_donors"]:
            chosen, seen, used = [], set(), []
            for i in order_p:
                new = [c for c in (int(nk[i]), int(t[i])) if donors[c] == dn and c not in seen]
                if not new:
                    continue
                for c in new[: want - len(chosen)]:
                    seen.add(c)
                    chosen.append(c)
                used.append(int(i))
                if len(chosen) == want:
                    break
            if len(chosen) < want:
                raise SystemExit(f"donor {dn}: only {len(chosen)} pair cells < {want}")
            chosen = np.sort(np.array(chosen, dtype=np.int64))
            pair_cells.append(chosen)
            pairs_used[dn] = used
            for cls in ("NK", "T"):
                ref = chosen[key[chosen] == cls]
                if ref.size == 0:
                    continue
                pool = ev[(donors[ev] == dn) & (key[ev] == cls) & (ctype[ev] != GDT158) & ~in_any[ev]]
                edges = np.quantile(ntok[ref], [0.2, 0.4, 0.6, 0.8])
                sp, sr = np.digitize(ntok[pool], edges), np.digitize(ntok[ref], edges)
                for q in range(5):
                    m = int(np.sum(sr == q))
                    cand = pool[sp == q]
                    if m > cand.size:
                        raise SystemExit(f"random pool too small: donor {dn} class {cls} quintile {q}")
                    if m:
                        rand.append(rng.choice(cand, m, replace=False))
        pair_cells = np.concatenate(pair_cells)
        rand = np.sort(np.concatenate(rand))
        sel = np.concatenate([pair_cells, rand])
        role = np.array(["pair"] * pair_cells.size + ["random"] * rand.size)
        info = {"pool": "site4 test_primary donors " + ", ".join(C["primary_donors"]), "n_no_gdT158_pairs": int(nk.size),
                "pairs_used_per_donor": {dn: len(v) for dn, v in pairs_used.items()},
                "pair_index_order_used": {dn: v for dn, v in pairs_used.items()}}
    if sel.size > C["max_cells"] or np.unique(sel).size != sel.size:
        raise SystemExit(f"{sel.size} cells (> max_cells or duplicates)")
    order = rng.permutation(sel.size)  # processing order: roles interleaved, so a partial run is a random subset
    info.update({"n_cells": int(sel.size),
                 "by_role_donor_class": {f"{r}/{dn}/{c}": int(np.sum((role == r) & (donors[sel] == dn) & (key[sel] == c)))
                                         for r in np.unique(role) for dn in np.unique(donors[sel]) for c in ("NK", "T")},
                 "ntokens_median_by_role": {str(r): float(np.median(ntok[sel[role == r]])) for r in np.unique(role)}})
    atomic_savez(out / "e5_cells.npz", cells=sel, role=role, order=order, cls=key[sel], donor=donors[sel])
    write_json(out / "e5_cells.json", info)
    say(f"cells: {info['by_role_donor_class']}; token-count medians {info['ntokens_median_by_role']}")
    return info


def load_cells(out: Path):
    with np.load(out / "e5_cells.npz", allow_pickle=False) as d:
        return d["cells"], d["role"].astype(str), d["order"], d["cls"].astype(str)


# ============================================================================ head

def load_head(path: Path, meta: dict, device):
    import torch

    from teddy_mm.models import MLP, AdtDecoder
    blob = torch.load(path, map_location="cpu", weights_only=False)
    w0 = blob["mlp"]["net.0.weight"]
    hidden, z_dim = int(w0.shape[0]), int(w0.shape[1])
    n_adt = meta["adt"].shape[1]
    mlp, dec = MLP(z_dim, z_dim, hidden=hidden), AdtDecoder(z_dim, n_adt, hidden=hidden)
    mlp.load_state_dict(blob["mlp"])
    dec.load_state_dict(blob["dec"])
    mlp.eval().to(device)
    dec.eval().to(device)
    sf = float(np.median(meta["adt_size_factor"][meta["split"] == "train"]))
    p95 = lp.p95_replicate(meta)
    col = {n: j for j, n in enumerate(meta["adt_names"])}
    cols = torch.tensor([col[p] for p in PROTEINS], device=device)
    p95v = torch.tensor([p95[p] for p in PROTEINS], dtype=torch.float32, device=device)

    def yhat(z):  # z: [B, d] raw gene-mean -> [B, 4] predicted ADT / train p95 (unclipped)
        zn = z / (torch.linalg.vector_norm(z, dim=-1, keepdim=True) + 1e-6)
        mu = dec(mlp(zn), torch.full((z.shape[0],), sf, device=z.device))[0]
        return mu.index_select(-1, cols) / p95v

    def score(y):  # [B, 4] -> [B]
        return y[:, :3].mean(-1) - y[:, 3]

    info = {"path": str(path), "hidden": hidden, "z_dim": z_dim, "train_median_size_factor": sf,
            "p95": {p: float(p95[p]) for p in PROTEINS}}
    return yhat, score, info


# ============================================================================ directions

def stage_directions(a, meta, emb, key) -> dict:
    import torch
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss
    torch.set_num_threads(_THREADS)
    out = a.out_dir
    split = meta["split"]
    rng = np.random.default_rng(0)
    rows = {}
    for sp in ("train", "val"):
        r = []
        for cls in ("NK", "T"):
            pool = np.where((split == sp) & (key == cls))[0]
            if a.dir_max_train and pool.size > a.dir_max_train:
                pool = np.sort(rng.choice(pool, a.dir_max_train, replace=False))
            r.append(pool)
        rows[sp] = np.sort(np.concatenate(r))
    both = np.concatenate([rows["train"], rows["val"]])
    o = np.argsort(both)
    both_sorted = both[o]
    t0 = time.time()
    seq_len = int(emb.manifest.get("seq_len", 2048))
    norm = float(emb.manifest.get("normalize_total", 1e4))
    X0, info0 = lp.input_layer_means(meta, both_sorted, a.ckpt, a.medians, seq_len, norm, emb.ntokens[both_sorted])
    say(f"directions: input gene-means of {both.size} train/val NK-T cells in {time.time() - t0:.0f}s")
    pos = {c: i for i, c in enumerate(both_sorted)}
    D = emb.d
    dmu = np.zeros((N_DEPTHS, D), np.float64)
    wraw = np.zeros((N_DEPTHS, D), np.float64)
    probe = []
    ytr = (key[rows["train"]] == "NK").astype(int)
    yva = (key[rows["val"]] == "NK").astype(int)
    for l in range(N_DEPTHS):
        if l == 0:
            Xtr = X0[[pos[c] for c in rows["train"]]]
            Xva = X0[[pos[c] for c in rows["val"]]]
        else:
            Xtr = np.asarray(emb.layer_means[l - 1][rows["train"]], np.float32)
            Xva = np.asarray(emb.layer_means[l - 1][rows["val"]], np.float32)
        dmu[l] = Xtr[ytr == 1].mean(0, dtype=np.float64) - Xtr[ytr == 0].mean(0, dtype=np.float64)
        mu, sd = Xtr.mean(0, dtype=np.float64), Xtr.std(0, dtype=np.float64)
        sd[sd < 1e-8] = 1.0
        best = None
        for C in lp.LR_C_GRID:
            clf = LogisticRegression(C=C, max_iter=3000, class_weight="balanced", random_state=0)
            clf.fit((Xtr - mu) / sd, ytr)
            ll = log_loss(yva, clf.predict_proba((Xva - mu) / sd)[:, 1], labels=[0, 1])
            if best is None or ll < best[0]:
                best = (ll, C, clf)
        w = best[2].coef_[0] / sd
        wraw[l] = w
        probe.append({"depth": l, "C": best[1], "val_logloss": rnd(best[0]), "logit_gap": float(w @ dmu[l]),
                      "dmu_norm": float(np.linalg.norm(dmu[l])),
                      "cos_probe_vs_mudiff": float(w @ dmu[l] / (np.linalg.norm(w) * np.linalg.norm(dmu[l])))})
        say(f"  depth {l:2d}: ||dmu|| {probe[-1]['dmu_norm']:.4f}  probe C {best[1]}  cos(probe, dmu) "
            f"{probe[-1]['cos_probe_vs_mudiff']:.3f}")
    # head-score gap on training NK / T cells (official z, the head's real input)
    yhat, score, hinfo = load_head(a.head_ckpt, meta, torch.device("cpu"))
    z = np.asarray(emb.z[rows["train"]], np.float32)
    with torch.no_grad():
        Y = torch.cat([yhat(torch.from_numpy(z[s:s + 4096])) for s in range(0, z.shape[0], 4096)]).numpy().astype(np.float64)
    S = Y[:, :3].mean(1) - Y[:, 3]
    head_gap_s = float(S[ytr == 1].mean() - S[ytr == 0].mean())
    head_gap_p = (Y[ytr == 1].mean(0) - Y[ytr == 0].mean(0))
    u0 = dmu[0] / np.linalg.norm(dmu[0])
    r = np.random.default_rng(SEEDS["random_directions"]).standard_normal(D)
    u_rand = r / np.linalg.norm(r)
    info = {"n_train": {"NK": int(ytr.sum()), "T": int((1 - ytr).sum())}, "n_val": {"NK": int(yva.sum()), "T": int((1 - yva).sum())},
            "subsampled": bool(a.dir_max_train), "probe": probe, "head": hinfo, "head_gap_s": head_gap_s,
            "head_gap_per_protein": dict(zip(PROTEINS, map(float, head_gap_p))), "input_layer": info0,
            "runtime_sec": round(time.time() - t0, 1)}
    atomic_savez(out / "e5_directions.npz", dmu=dmu, wraw=wraw, u_nkt=u0, u_rand=u_rand,
                 head_gap_s=np.array(head_gap_s), head_gap_p=head_gap_p)
    write_json(out / "e5_directions.json", info)
    say(f"directions done ({info['runtime_sec']}s); head-score train NK-T gap {head_gap_s:.4f}")
    return info


def load_directions(out: Path) -> dict:
    with np.load(out / "e5_directions.npz") as d:
        return {k: d[k] for k in d.files}


# ============================================================================ respond

def layer_fn(layer, x):
    """One nn.TransformerEncoderLayer (eval: dropout inactive), same weights, explicit attention
    softmax(QK^T / sqrt(d_h)) V so that forward-mode AD runs on every device. x: [B, L, d], no padding."""
    import torch
    import torch.nn.functional as F
    mha = layer.self_attn
    B, L, d = x.shape
    H = mha.num_heads
    dh = d // H

    def sa(y):
        q, k, v = F.linear(y, mha.in_proj_weight, mha.in_proj_bias).split(d, dim=-1)
        q = q.reshape(B, L, H, dh).transpose(1, 2)
        k = k.reshape(B, L, H, dh).transpose(1, 2)
        v = v.reshape(B, L, H, dh).transpose(1, 2)
        w = torch.softmax(torch.matmul(q, k.transpose(-1, -2)) * (dh ** -0.5), dim=-1)
        o = torch.matmul(w, v).transpose(1, 2).reshape(B, L, d)
        return F.linear(o, mha.out_proj.weight, mha.out_proj.bias)

    def ff(y):
        return layer.linear2(layer.activation(layer.linear1(y)))

    if layer.norm_first:
        x = x + sa(layer.norm1(x))
        return x + ff(layer.norm2(x))
    x = layer.norm1(x + sa(x))
    return layer.norm2(x + ff(x))


def tokenise(meta, sel, a, emb, vocab_pad):
    from teddy_mm.teddy_encoder import load_gene_medians, median_factors, official_values, rank_encode_official
    seq_len = int(emb.manifest.get("seq_len", 2048))
    norm = float(emb.manifest.get("normalize_total", 1e4))
    factors = median_factors(meta["rna_names"], load_gene_medians(a.medians))
    toks, ntok = [], []
    for s in range(0, len(sel), 256):
        vals = official_values(meta["rna"][sel[s:s + 256]].toarray(), factors, norm)
        t, m = rank_encode_official(vals, meta["token_ids"], max_len=seq_len, pad_id=vocab_pad, pad_to=seq_len)
        for i in range(t.shape[0]):
            n = int(m[i].sum())
            toks.append(t[i, :n].copy())
            ntok.append(n)
    return toks, np.array(ntok)


def stage_respond(a, meta, emb, dirs):
    import torch
    from torch.func import jvp, vmap

    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import load_pad_id, load_teddy, load_vocab
    torch.set_num_threads(_THREADS)
    out = a.out_dir
    cdir = out / "cells"
    cdir.mkdir(exist_ok=True)
    sel, role, order, cls = load_cells(out)
    n_full = a.n_full_ladder if a.n_full_ladder is not None else N_FULL_LADDER
    full_set = {int(i) for i in order[:n_full]}  # cells with every eps rung (the first of the seeded order)
    if a.limit:
        order = order[: a.limit]
    device = resolve_device(a.device)
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    rna_names = set(map(str, meta["rna_names"]))
    gene_tok = {}
    for s in GENE_SYMBOLS:
        ens = GENES[s][0]
        if ens not in rna_names or ens not in vocab:
            raise SystemExit(f"gene {s} {ens} missing from the processed genes or the vocab")
        gene_tok[s] = vocab[ens]
    model = load_teddy(a.ckpt, device)
    for p_ in model.parameters():
        p_.requires_grad_(False)
    layers = list(model.encoder.layers)
    E = model.embeddings.weight
    Pm = model.position_embeddings.weight
    yhat, score, _ = load_head(a.head_ckpt, meta, device)
    f32 = lambda x: torch.as_tensor(np.asarray(x, np.float32), device=device)  # noqa: E731
    dmu, wraw = f32(dirs["dmu"]), f32(dirs["wraw"])
    u_nkt, u_rand = f32(dirs["u_nkt"]), f32(dirs["u_rand"])
    eps_list = [float(e) for e in REGISTRATION["eps"]]
    sha = registration_sha()
    todo = []
    for i in order:
        f = cdir / f"cell_{int(sel[i])}.npz"
        if f.exists():
            with np.load(f, allow_pickle=False) as d:
                if json.loads(str(d["meta"])).get("registration_sha256") != sha:
                    raise SystemExit(f"{f} was written under another registration; use a fresh --out-dir")
        else:
            todo.append(int(i))
    say(f"respond: {len(order)} cells in the order, {len(order) - len(todo)} already done, {len(todo)} to go; "
        f"device {device}, float32, jvp chunk {a.jvp_chunk}, fd batch {a.fd_batch}, eps {eps_list}")
    if not todo:
        return
    toks, ntok = tokenise(meta, sel[todo], a, emb, pad_id)
    ref_ntok = emb.ntokens[sel[todo]]
    if np.any(ntok != ref_ntok):
        raise SystemExit(f"{int(np.sum(ntok != ref_ntok))} cells get another token count than the official run")
    t_start = time.time()
    n_done = 0
    n_cpu64 = 0
    fd64_set = {int(i) for i in order[:N_FD64_CELLS]}
    layers64 = None

    def pooled_and_own(h, pos_k):  # h: [K, 1, L, d] -> pooled [K, d], own [K, d]
        return h[:, 0].mean(1), h[torch.arange(h.shape[0], device=h.device), 0, pos_k]

    def forward_all(h0):  # explicit arithmetic (the JVP's primal): h0 [B, L, d] -> states per depth
        hs = [h0]
        h = h0
        for layer in layers:
            h = layer_fn(layer, h)
            hs.append(h)
        return hs

    def forward_module(h0):  # the module's own layer call (finite differences): pooled per depth [B, 13, d]
        zs = [h0.mean(1)]
        h = h0
        for layer in layers:
            h = layer(h)  # every position is a real token, so there is no padding mask
            zs.append(h.mean(1))
        return torch.stack(zs, 1)

    for jj, i in enumerate(todo):
        tc = time.time()
        g_id = int(sel[i])
        ids_np = toks[jj]
        L = ids_np.size
        ids = torch.from_numpy(ids_np).to(device)
        present = [(gi, s, int(np.nonzero(ids_np == gene_tok[s])[0][0])) for gi, s in enumerate(GENE_SYMBOLS)
                   if np.any(ids_np == gene_tok[s])]
        if a.max_genes:
            present = present[: a.max_genes]
        declared = set(gene_tok.values())
        cand = np.array([q for q in range(L) if int(ids_np[q]) not in declared])
        rng_c = np.random.default_rng([SEEDS["random_directions"], g_id])
        other_pos = np.sort(rng_c.choice(cand, min(N_OTHER, cand.size), replace=False))
        if present:  # the rand push: at one declared gene of the cell, else at a random token
            gi_r, s_r, t_r = present[int(rng_c.integers(len(present)))]
            rand_case = (gi_r, gene_tok[s_r], t_r, 2)
        else:
            q = int(rng_c.integers(L))
            rand_case = (-1, int(ids_np[q]), q, 2)
        full = int(i) in full_set
        meta_cell = {"cell": g_id, "role": str(role[i]), "class": str(cls[i]), "ntokens": int(L), "registration_sha256": sha,
                     "genes_present": [s for _, s, _ in present], "full_ladder": full}
        with torch.no_grad():
            h0 = (E[ids] + Pm[:L])[None]  # [1, L, d]
            hs = forward_all(h0)
            zb = torch.stack([h.mean(1)[0] for h in hs])  # [13, d]
            # module forward (official code path) on the same tokens
            ref_h, ref_lm = model.hidden_states(ids[None], torch.ones(1, L, dtype=torch.long, device=device), return_layer_means=True)
            chk_module = float((ref_lm[:, 0] - zb[1:]).abs().max())
            off = np.stack([np.asarray(emb.layer_means[l][g_id], np.float32) for l in range(emb.n_layers)])
            chk_official = float(np.abs(off - zb[1:].cpu().numpy()).max())
            chk_official_rel = float(np.abs(off - zb[1:].cpu().numpy()).max() / (np.abs(off).max() + 1e-12))
        z12 = zb[12].detach()
        y0 = yhat(z12[None])
        g_head = torch.func.grad(lambda z: score(yhat(z[None]))[0])(z12).detach()
        # case = (gene index or -1, token id, position, direction index)
        cases = [(gi, gene_tok[s], t, di) for gi, s, t in present for di in (0, 1)] + [rand_case]
        cases += [(-1, int(ids_np[q]), int(q), 3) for q in other_pos]
        K = len(cases)
        sg = torch.stack([E[tid].norm() for _, tid, _, _ in cases])  # [K]
        vdir = []
        for _, tid, _, di in cases:
            v = {"nkt": u_nkt, "rand": u_rand}.get(DIRECTIONS[di])
            vdir.append(v if v is not None else E[tid] / E[tid].norm())
        vdir = torch.stack(vdir)  # [K, d]
        pos_k = torch.tensor([t for _, _, t, _ in cases], device=device)
        T0 = torch.zeros(K, 1, L, model.d_model, device=device)
        T0[torch.arange(K, device=device), 0, pos_k] = sg[:, None] * vdir
        # ---------------- JVP chain
        J_pool = torch.zeros(K, N_DEPTHS, model.d_model, device=device)
        J_own = torch.zeros_like(J_pool)
        with torch.no_grad():
            for c0 in range(0, K, a.jvp_chunk):
                c1 = min(c0 + a.jvp_chunk, K)
                dh = T0[c0:c1]
                pk = pos_k[c0:c1]
                J_pool[c0:c1, 0], J_own[c0:c1, 0] = pooled_and_own(dh, pk)
                for li, layer in enumerate(layers):
                    f = lambda h, _layer=layer: layer_fn(_layer, h)  # noqa: E731
                    dh = vmap(lambda t, _h=hs[li], _f=f: jvp(_f, (_h,), (t,))[1])(dh)
                    J_pool[c0:c1, li + 1], J_own[c0:c1, li + 1] = pooled_and_own(dh, pk)
                del dh
            dY = vmap(lambda t: jvp(yhat, (z12[None],), (t[None],))[1][0])(J_pool[:, 12])  # [K, 4]
            dS = dY[:, :3].mean(-1) - dY[:, 3]
            # ---------------- finite differences through the module's own layers, one batch shape for every
            # sequence (the base state is job 0 of the first batch; the last batch is filled with a repeated job).
            # Every cell: eps_primary (+-e, +-e/2) and eps_abs 1e-3; the full-ladder cells: every eps and eps_abs.
            S = len(eps_list)
            cols = [("base", None, 0.0, None)]
            for si, e in enumerate(eps_list):
                if full or e == REGISTRATION["eps_primary"]:
                    cols += [("lad", si, c, q) for q, c in enumerate((e, -e, e / 2, -e / 2))]
            for ei, e in enumerate(EPS_ABS_CHECK):
                if full or ei == 0:
                    cols += [("abs", ei, c, q) for q, c in enumerate((e, -e))]
            colix = {(kind, j, q): ci for ci, (kind, j, _, q) in enumerate(cols)}
            coef_m = torch.tensor([c for _, _, c, _ in cols], device=device)[None].repeat(K, 1)
            for ci, (kind, _, c, _) in enumerate(cols):
                if kind == "abs":
                    coef_m[:, ci] = c / sg  # c * T0 = +-eps_abs * unit tangent
            jobs = [(0, 0)] + [(k, ci) for k in range(K) for ci in range(1, len(cols))]
            nb = a.fd_batch
            while len(jobs) % nb:
                jobs.append(jobs[-1])
            Zfd = torch.zeros(K, len(cols), N_DEPTHS, model.d_model, device=device)
            Yfd = torch.zeros(K, len(cols), len(PROTEINS), device=device)
            for b0 in range(0, len(jobs), nb):
                bj = jobs[b0:b0 + nb]
                kk = torch.tensor([k for k, _ in bj], device=device)
                cc = coef_m[kk, torch.tensor([ci for _, ci in bj], device=device)]
                zz = forward_module(h0.expand(nb, L, model.d_model) + cc[:, None, None] * T0[kk, 0])  # [nb, 13, d]
                yy = yhat(zz[:, 12])
                for r_, (k, ci) in enumerate(bj):
                    Zfd[k, ci] = zz[r_]
                    Yfd[k, ci] = yy[r_]
            z0, yb = Zfd[0, 0], Yfd[0, 0]
            chk_base = float((z0 - zb).abs().max())  # module path vs explicit path, unperturbed
        # linearity statistics on the pooled vector, per case, eps and depth (NaN where an eps was not run):
        # 0 err_plus  ||D+(e) - J|| / ||J||        D+(e) = (z(+e) - z0) / e
        # 1 err_minus ||D-(e) - J|| / ||J||        D-(e) = (z0 - z(-e)) / e
        # 2 slope_plus  ||D+(e) - D+(e/2)|| / ||D+(e/2)||
        # 3 slope_minus ||D-(e) - D-(e/2)|| / ||D-(e/2)||
        # 4 err_central ||Dc(e) - J|| / ||J||,  Dc = (D+ + D-) / 2 (odd part only)
        # 5 slope_central ||Dc(e) - Dc(e/2)|| / ||Dc(e/2)||
        # 6 kappa ||D+(e) - D-(e)|| / ||D+(e) + D-(e)|| (even / odd part)
        lin = np.full((K, S, 7, N_DEPTHS), np.nan, np.float32)
        fd_proj = np.full((K, S, 4, N_DEPTHS, 3), np.nan)  # slopes D+(e), D-(e), D+(e/2), D-(e/2) on (mu_diff, probe, head lens)
        fd_head = np.full((K, S, 4), np.nan)  # the same four slopes of the head score
        Jn = J_pool.norm(dim=-1).clamp(min=1e-30)  # [K, 13]
        nrm = lambda x: x.norm(dim=-1).clamp(min=1e-30)  # noqa: E731
        sc = lambda y: (y[..., :3].mean(-1) - y[..., 3])  # noqa: E731
        s0 = sc(yb)
        for si, e in enumerate(eps_list):
            if ("lad", si, 0) not in colix:
                continue
            qp, qm, qph, qmh = (colix[("lad", si, q)] for q in range(4))
            Dp, Dm = (Zfd[:, qp] - z0) / e, (z0 - Zfd[:, qm]) / e
            Dph, Dmh = (Zfd[:, qph] - z0) / (e / 2), (z0 - Zfd[:, qmh]) / (e / 2)
            Dc, Dch = (Dp + Dm) / 2, (Dph + Dmh) / 2
            st = [nrm(Dp - J_pool) / Jn, nrm(Dm - J_pool) / Jn, nrm(Dp - Dph) / nrm(Dph), nrm(Dm - Dmh) / nrm(Dmh),
                  nrm(Dc - J_pool) / Jn, nrm(Dc - Dch) / nrm(Dch), nrm(Dp - Dm) / nrm(Dp + Dm)]
            lin[:, si] = torch.stack(st, 1).cpu().numpy()
            for hi, D in enumerate((Dp, Dm, Dph, Dmh)):
                fd_proj[:, si, hi, :, 0] = (D * dmu[None]).sum(-1).cpu().double().numpy()
                fd_proj[:, si, hi, :, 1] = (D * wraw[None]).sum(-1).cpu().double().numpy()
                fd_proj[:, si, hi, :, 2] = (D @ g_head).cpu().double().numpy()
            for hi, (qq, sgn, ee) in enumerate(((qp, 1, e), (qm, -1, e), (qph, 1, e / 2), (qmh, -1, e / 2))):
                fd_head[:, si, hi] = (sgn * (sc(Yfd[:, qq]) - s0) / ee).cpu().double().numpy()
        # experiments_v3 JVP check: central FD on the unit tangent at eps_abs; relative error and cosine per depth
        jchk = np.full((K, len(EPS_ABS_CHECK), 2, N_DEPTHS), np.nan, np.float32)
        for ei in range(len(EPS_ABS_CHECK)):
            if ("abs", ei, 0) not in colix:
                continue
            qp, qm = colix[("abs", ei, 0)], colix[("abs", ei, 1)]
            Dc = (Zfd[:, qp] - Zfd[:, qm]) / (2 * coef_m[:, qp])[:, None, None]
            jchk[:, ei, 0] = (nrm(Dc - J_pool) / Jn).cpu().numpy()
            jchk[:, ei, 1] = ((Dc * J_pool).sum(-1) / (nrm(Dc) * Jn)).cpu().numpy()
        # float64 CPU cross-check of the JVP for the first cells (first 2 cases)
        cpu64 = None
        if n_cpu64 < a.cpu64_check:
            cpu64 = cpu64_jvp_check(layers, h0.detach().cpu().double(), T0[:2].detach().cpu().double(),
                                    J_pool[:2].detach().cpu().double())
            n_cpu64 += 1
        fd64 = None
        if int(i) in fd64_set:
            import copy
            if layers64 is None:
                layers64 = [copy.deepcopy(l_).cpu().double() for l_ in layers]
            pick = [next((k for k, c in enumerate(cases) if c[3] == di), None) for di in range(len(DIRECTIONS))]
            pick = [k for k in pick if k is not None]
            fd64 = fd64_check(layers64, h0.detach().cpu().double(), T0[pick].detach().cpu().double(),
                              sg[pick].detach().cpu().double(), J_pool[pick].detach().cpu().double())
            fd64["cases"] = pick
        meta_cell.update({"check_explicit_vs_module_layer_means_max_abs": chk_module,
                          "check_vs_official_fp16_layer_means_max_abs": chk_official,
                          "check_vs_official_fp16_layer_means_max_rel": chk_official_rel,
                          "check_fd_base_vs_jvp_base_max_abs": chk_base,
                          "cpu64": cpu64, "fd64": fd64, "sec": round(time.time() - tc, 2)})
        atomic_savez(cdir / f"cell_{g_id}.npz",
                     meta=np.array(json.dumps(meta_cell)), n_cases=np.array(K),
                     gene=np.array([gi for gi, _, _, _ in cases], np.int16),
                     token_id=np.array([tid for _, tid, _, _ in cases], np.int32), jvp_check=jchk,
                     direction=np.array([di for _, _, _, di in cases], np.int8),
                     pos=np.array([t for _, _, t, _ in cases], np.int32),
                     s_gene=sg.cpu().numpy().astype(np.float64),
                     J_pool=J_pool.cpu().numpy(), J_own=J_own.cpu().numpy(),
                     dY=dY.cpu().double().numpy(), dS=dS.cpu().double().numpy(),
                     g_head=g_head.cpu().numpy(), y0=y0[0].detach().cpu().numpy(), z_base=zb.cpu().numpy(),
                     lin=lin, fd_proj=fd_proj, fd_head=fd_head)
        n_done += 1
        el = time.time() - t_start
        rate = el / n_done
        say(f"  cell {jj + 1}/{len(todo)} ({role[i]}, {cls[i]}, {L} tokens, {K} cases) {time.time() - tc:.1f}s; "
            f"module check {chk_module:.2e}, official fp16 check {chk_official:.2e}"
            + (f", cpu64 rel {cpu64['max_rel_diff_pooled']:.2e}" if cpu64 else "")
            + f"; mean {rate:.1f}s/cell, ETA {rate * (len(todo) - n_done) / 60:.1f} min")
        write_json(out / "e5_progress.json", {"done_this_call": n_done, "to_do_this_call": len(todo),
                                               "sec_per_cell": round(rate, 2), "eta_min": round(rate * (len(todo) - n_done) / 60, 1),
                                               "updated": time.strftime("%Y-%m-%d %H:%M:%S")})
        if device.type == "mps":
            torch.mps.empty_cache()
    say(f"respond done: {n_done} cells in {(time.time() - t_start) / 60:.1f} min")


def fd64_check(layers64, h0, T, sg, J32):
    """experiments_v3 JVP check in float64 on the CPU: central FD through the module's own layers on the unit tangent
    (c = eps_abs / ||E[token]||, c * T = eps_abs * v) against the device float32 JVP; rel. error and cosine per depth."""
    import torch
    t0 = time.time()

    def pooled(x):
        zs = [x.mean(1)]
        h = x
        for l_ in layers64:
            h = l_(h)
            zs.append(h.mean(1))
        return torch.stack(zs, 1)  # [B, 13, d]
    out = np.zeros((T.shape[0], len(EPS_ABS_CHECK), 2, N_DEPTHS))
    with torch.no_grad():
        for ei, e in enumerate(EPS_ABS_CHECK):
            c = (e / sg)[:, None, None]
            x = torch.cat([h0 + c * T[:, 0], h0 - c * T[:, 0]])
            Z = pooled(x)
            n = T.shape[0]
            Dc = (Z[:n] - Z[n:]) / (2 * c)
            Jn = J32.norm(dim=-1).clamp(min=1e-300)
            out[:, ei, 0] = ((Dc - J32).norm(dim=-1) / Jn).numpy()
            out[:, ei, 1] = ((Dc * J32).sum(-1) / (Dc.norm(dim=-1).clamp(min=1e-300) * Jn)).numpy()
    return {"rel_err_cos": out.tolist(), "sec": round(time.time() - t0, 1)}


def cpu64_jvp_check(layers, h0, T, J32):
    """float64 CPU JVP of the same cases; relative difference to the device float32 JVP (pooled, per depth)."""
    import copy

    import torch
    from torch.func import jvp, vmap
    L64 = [copy.deepcopy(l_).cpu().double() for l_ in layers]
    t0 = time.time()
    with torch.no_grad():
        dh = T
        hs = [h0]
        h = h0
        for l_ in L64:
            h = layer_fn(l_, h)
            hs.append(h)
        pooled = [dh[:, 0].mean(1)]
        for li, l_ in enumerate(L64):
            dh = vmap(lambda t, _h=hs[li], _l=l_: jvp(lambda x: layer_fn(_l, x), (_h,), (t,))[1])(dh)
            pooled.append(dh[:, 0].mean(1))
        P = torch.stack(pooled, 1)  # [2, 13, d]
    rel = ((P - J32).norm(dim=-1) / P.norm(dim=-1).clamp(min=1e-300)).numpy()
    return {"max_rel_diff_pooled": float(rel.max()), "rel_diff_by_depth": [float(x) for x in rel.max(0)],
            "sec": round(time.time() - t0, 1)}


# ============================================================================ report

def load_cases(out: Path, sel, role, cls, donors, dirs):
    files = sorted((out / "cells").glob("cell_*.npz"))
    want = {int(c) for c in sel}
    keys = ("cell", "role", "cls", "donor", "gene", "direction", "ntok", "pos", "G", "P", "H", "Gown", "Nrel", "cosmu",
            "dS", "dY", "lin", "G_fd", "jchk", "own_norm_gain", "logr", "gsign")
    rows = {k: [] for k in keys}
    metas = []
    dmu, wraw = dirs["dmu"], dirs["wraw"]
    dmu2 = (dmu ** 2).sum(1)
    dmun = np.sqrt(dmu2)
    logit_gap = (wraw * dmu).sum(1)
    head_gap = float(dirs["head_gap_s"])
    head_gap_p = dirs["head_gap_p"]
    rl = {int(c): (r, k) for c, r, k in zip(sel, role, cls)}
    sha = registration_sha()
    for f in files:
        with np.load(f, allow_pickle=False) as d:
            m = json.loads(str(d["meta"]))
            if m["cell"] not in want:
                continue
            if m.get("registration_sha256") != sha:
                raise SystemExit(f"{f} was written under another registration")
            metas.append(m)
            K = int(d["n_cases"])
            if K == 0:
                continue
            n = m["ntokens"]
            A0 = d["s_gene"] / (n * dmun[0])  # [K]
            Jp, Jo = d["J_pool"].astype(np.float64), d["J_own"].astype(np.float64)
            proj = np.einsum("kld,ld->kl", Jp, dmu)
            Jn = np.linalg.norm(Jp, axis=-1)
            own_n = np.linalg.norm(Jo, axis=-1)
            gene = d["gene"].astype(int)
            r, k = rl[m["cell"]]
            add = {"G": proj / dmu2[None] / A0[:, None],
                   "Gown": np.einsum("kld,ld->kl", Jo, dmu) / dmu2[None] / n / A0[:, None],
                   "P": np.einsum("kld,ld->kl", Jp, wraw) / logit_gap[None] / A0[:, None],
                   "H": (Jp @ d["g_head"].astype(np.float64)) / head_gap / A0[:, None],
                   "Nrel": Jn / dmun[None] / A0[:, None],
                   "cosmu": proj / (Jn * dmun[None] + 1e-300),
                   "dS": d["dS"] / head_gap / A0,
                   "dY": d["dY"] / head_gap_p[None] / A0[:, None],
                   "lin": d["lin"], "jchk": d["jvp_check"],
                   "G_fd": d["fd_proj"][..., 0] / dmu2[None, None, None] / A0[:, None, None, None],  # [K, S, 4, 13]
                   "own_norm_gain": own_n / own_n[:, :1],
                   "logr": np.log(np.abs(proj / dmun[None]) + 1e-30),
                   "gsign": np.array([(1.0 if GENES[GENE_SYMBOLS[g]][1] == "NK" else -1.0) if g >= 0 else 1.0 for g in gene]),
                   "gene": gene, "direction": d["direction"].astype(int), "pos": d["pos"]}
            for kk, v in add.items():
                rows[kk].append(v)
            rows["cell"] += [m["cell"]] * K
            rows["role"] += [r] * K
            rows["cls"] += [k] * K
            rows["donor"] += [donors[m["cell"]]] * K
            rows["ntok"] += [n] * K
    C = {}
    for kk, v in rows.items():
        C[kk] = (np.concatenate(v) if isinstance(v[0], np.ndarray) else np.array(v)) if v else np.zeros(0)
    return C, metas


class Boot:
    """Resampling of cells (all of a cell's cases kept): 'cell' = cells with replacement within each stratum;
    'donor' = two-stage, donors with replacement, then cells with replacement within each drawn donor and stratum."""

    def __init__(self, cell, donor, stratum, n_boot: int, seed: int):
        o = np.argsort(cell, kind="stable")
        cs, start = np.unique(cell[o], return_index=True)
        self.groups = np.split(o, start[1:])
        first = o[start]
        cd, cst = donor[first], stratum[first]
        self.donors = sorted(set(cd.tolist()))
        rng = np.random.default_rng(seed)
        strata = [np.where(cst == s)[0] for s in np.unique(cst)]
        bydon = {d_: [np.where((cd == d_) & (cst == s))[0] for s in np.unique(cst[cd == d_])] for d_ in self.donors}
        self.samples = {"cell": [], "donor": []}
        for _ in range(n_boot):
            self.samples["cell"].append(np.concatenate([g[rng.integers(0, g.size, g.size)] for g in strata]))
        for _ in range(n_boot):
            ds = rng.integers(0, len(self.donors), len(self.donors))
            self.samples["donor"].append(np.concatenate([g[rng.integers(0, g.size, g.size)]
                                                         for x in ds for g in bydon[self.donors[x]]]))
        self.n = sum(g.size for g in self.groups)

    def cases(self, cell_sample):
        return np.concatenate([self.groups[g] for g in cell_sample])


def pct(v):
    v = np.asarray(v, np.float64)
    return np.percentile(v, [2.5, 97.5], axis=0)


def subset_summary(X: dict, donor: np.ndarray, boot: Boot) -> dict:
    """X: name -> [cases, D] per-case quantities of one subset. Median over cases with two-stage ('donor') and
    cell-cluster ('cell') 95% intervals, per-donor medians; for 'dG' also the bootstrap share of each argmin."""
    names = list(X)
    widths = [X[k].shape[1] for k in names]
    M = np.concatenate([X[k] for k in names], 1)
    full = np.median(M, 0)
    res = {}
    boots = {}
    for lev in ("donor", "cell"):
        B = np.array([np.median(M[boot.cases(s)], 0) for s in boot.samples[lev]])
        boots[lev] = B
    per_d = {d_: np.median(M[donor == d_], 0) for d_ in sorted(set(donor.tolist()))}
    o = 0
    for k, w in zip(names, widths):
        sl = slice(o, o + w)
        res[k] = {"median": full[sl].tolist(),
                  "ci_two_stage": pct(boots["donor"][:, sl]).T.tolist(),
                  "ci_cell": pct(boots["cell"][:, sl]).T.tolist(),
                  "per_donor": {d_: v[sl].tolist() for d_, v in per_d.items()}}
        if k in ("dG", "logg"):
            for lev in ("donor", "cell"):
                A = np.argmin(boots[lev][:, sl], 1) + 1
                res[k][f"argmin_share_{'two_stage' if lev == 'donor' else 'cell'}"] = {int(v): float(np.mean(A == v)) for v in np.unique(A)}
        o += w
    return res


def joint_diff(xa, xb, cell, donor, stratum, boot_seed, n_boot):
    """median(xa over stratum 'a') - median(xb over stratum 'b') with a joint two-stage bootstrap (one donor draw,
    cells resampled within donor and stratum) and per-donor values. xa, xb: per-case values over the union."""
    b = Boot(cell, donor, stratum, n_boot, boot_seed)
    A, Bm = stratum == "a", stratum == "b"
    full = float(np.median(xa[A]) - np.median(xb[Bm]))

    def st(ix):
        ia, ib = ix[A[ix]], ix[Bm[ix]]
        return np.median(xa[ia]) - np.median(xb[ib]) if ia.size and ib.size else np.nan
    out = {"diff": full}
    for lev, nm in (("donor", "ci_two_stage"), ("cell", "ci_cell")):
        v = np.array([st(b.cases(s)) for s in b.samples[lev]])
        v = v[np.isfinite(v)]
        out[nm] = pct(v).tolist() if v.size else None
    out["per_donor"] = {d_: float(np.median(xa[A & (donor == d_)]) - np.median(xb[Bm & (donor == d_)]))
                        for d_ in sorted(set(donor.tolist())) if (A & (donor == d_)).any() and (Bm & (donor == d_)).any()}
    return out


def probe_curves(probe_dir: Path) -> dict:
    R = json.loads((probe_dir / "mode_a_results.json").read_text())
    by = {l_["name"]: l_ for l_ in R["layers"]}
    names = ["input"] + [f"layer_{i}" for i in range(1, 13)]
    cur = {"depth_names": names, "order_rate_no_gdT158": [], "auc": [], "gap_ratio_mean_no_gdT158": [],
           "gap_ratio_mean_all": [], "gap_ratio_no_gdT158": {p: [] for p in PROTEINS}, "centred_distance_ratio_no_gdT158": []}
    for nm in names:
        l_ = by[nm]
        cur["order_rate_no_gdT158"].append(l_["nk_t_probe"]["pairs"]["no_gdT158"]["order_rate"])
        cur["auc"].append(l_["nk_t_probe"]["site4_auc"])
        g = l_["protein_probe"]["gap"]["no_gdT158"]
        rs = [g[p]["ratio"] for p in PROTEINS]
        for p, r_ in zip(PROTEINS, rs):
            cur["gap_ratio_no_gdT158"][p].append(r_)
        cur["gap_ratio_mean_no_gdT158"].append(float(np.mean(rs)))
        cur["gap_ratio_mean_all"].append(float(np.mean([l_["protein_probe"]["gap"]["all"][p]["ratio"] for p in PROTEINS])))
        cur["centred_distance_ratio_no_gdT158"].append(l_["cosine"]["centred"]["no_gdT158"]["distance_ratio_neighbour_over_random"])
    gm = np.array(cur["gap_ratio_mean_no_gdT158"])
    cur["L_probe_argmin_gap_ratio"] = int(np.argmin(gm[1:]) + 1)
    cur["L_probe_largest_drop_gap_ratio"] = int(np.argmin(np.diff(gm)) + 1)
    cur["L_probe_largest_drop_gap_ratio_all"] = int(np.argmin(np.diff(np.array(cur["gap_ratio_mean_all"]))) + 1)
    cur["L_probe_largest_drop_distance_ratio"] = int(np.argmin(np.diff(np.array(cur["centred_distance_ratio_no_gdT158"], float))) + 1)
    cur["source"] = str(probe_dir / "mode_a_results.json")
    tp = probe_dir / "token" / "mode_a_token_results.json"
    if tp.exists():
        T = json.loads(tp.read_text())
        cur["token_probes_ridge_gap_ratio_no_gdT158"] = {dn: {fs: {p: T["feature_sets"][dn][fs]["ridge"]["gap_no_gdT158"][p]["ratio"]
                                                                  for p in PROTEINS} for fs in ("genemean", "tokens")}
                                                         for dn in T["feature_sets"]}
    return cur


def spearman(a, b):
    from scipy.stats import spearmanr
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(spearmanr(a[ok], b[ok]).statistic) if ok.sum() >= 3 else None


def _all_donors(per_donor: dict, fn) -> bool:
    return bool(per_donor) and all(fn(v) for v in per_donor.values())


def stage_report(a, meta, reg_info):
    out = a.out_dir
    sel, role, _, cls = load_cells(out)
    donors = load_donors(a.processed)
    dirs = load_directions(out)
    dinfo = json.loads((out / "e5_directions.json").read_text())
    C, metas = load_cases(out, sel, role, cls, donors, dirs)
    tau = REGISTRATION["linearity"]["tau"]
    eps = REGISTRATION["eps"]
    si0 = eps.index(REGISTRATION["eps_primary"])
    seed = SEEDS["bootstrap"]
    R = {"experiment": "E5", "smoke_test": bool(a.smoke), "smoke_note": a.smoke, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
         "script_version": SCRIPT_VERSION, "registration_sha256": reg_info.get("registration_sha256"),
         "addendum_sha256": reg_info["addendum_sha256"], "registration_check": reg_info, "cell_pool": a.cell_pool,
         "n_boot": a.n_boot, "n_cells_selected": int(sel.size), "n_cells_done": len(metas),
         "complete": len(metas) == sel.size, "n_cases": int(C["G"].shape[0]) if C["G"].size else 0,
         "directions_info": dinfo}
    if not C["G"].size:
        write_json(out / "E5_results.json", R)
        say("no cases yet")
        return
    roles = [r for r in ("pair", "random", "val") if np.any(C["role"] == r)]
    R["cases_by_role_direction"] = {f"{r}/{d}": int(np.sum((C["role"] == r) & (C["direction"] == di)))
                                    for r in roles for di, d in enumerate(DIRECTIONS)}
    R["checks"] = {
        "explicit_vs_module_layer_means_max_abs": max(m["check_explicit_vs_module_layer_means_max_abs"] for m in metas),
        "vs_official_fp16_layer_means_max_rel": max(m["check_vs_official_fp16_layer_means_max_rel"] for m in metas),
        "fd_base_vs_jvp_base_max_abs": max(m["check_fd_base_vs_jvp_base_max_abs"] for m in metas),
        "cpu64_jvp_max_rel_diff": max((m["cpu64"]["max_rel_diff_pooled"] for m in metas if m.get("cpu64")), default=None),
        "head_lens_12_vs_head_jvp_max_abs": float(np.max(np.abs(C["H"][:, 12] - C["dS"]))),
        "nkt_G0_is_1_max_abs": float(np.max(np.abs(C["G"][C["direction"] == 0, 0] - 1))) if np.any(C["direction"] == 0) else None,
        "cells_without_declared_genes": int(sum(1 for m in metas if not m["genes_present"])),
        "sec_per_cell_median": float(np.median([m["sec"] for m in metas])),
    }
    # ---------------- experiments_v3 JVP check (gate)
    jc = C["jchk"]  # [cases, eps_abs, (rel err, cos), 13]
    jh = [np.isfinite(jc[:, ei, 0, 1]) for ei in range(len(EPS_ABS_CHECK))]
    jok = [((jc[jh[ei], ei, 0, 1:].max(1) <= 0.05) & (jc[jh[ei], ei, 1, 1:].min(1) >= 0.99)) for ei in range(len(EPS_ABS_CHECK))]
    R["jvp_check"] = {"eps_abs": list(EPS_ABS_CHECK), "n_cases": {str(e): int(jh[ei].sum()) for ei, e in enumerate(EPS_ABS_CHECK)},
                      "pass_share": {str(e): float(jok[ei].mean()) if jok[ei].size else None for ei, e in enumerate(EPS_ABS_CHECK)},
                      "median_rel_err_by_depth": {str(e): np.median(jc[jh[ei], ei, 0], 0).tolist() if jh[ei].any() else None
                                                  for ei, e in enumerate(EPS_ABS_CHECK)},
                      "p05_cos_by_depth": {str(e): np.percentile(jc[jh[ei], ei, 1], 5, axis=0).tolist() if jh[ei].any() else None
                                           for ei, e in enumerate(EPS_ABS_CHECK)},
                      "float32_device_diagnostic": True}
    f64 = [np.array(m["fd64"]["rel_err_cos"]) for m in metas if m.get("fd64")]
    if f64:
        F = np.concatenate(f64)  # [checked cases, eps_abs, (rel, cos), 13]
        fok = [(F[:, ei, 0, 1:].max(1) <= 0.05) & (F[:, ei, 1, 1:].min(1) >= 0.99) for ei in range(len(EPS_ABS_CHECK))]
        R["jvp_check"]["float64_cpu"] = {"n_cases": int(F.shape[0]), "pass_share": {str(e): float(fok[ei].mean()) for ei, e in enumerate(EPS_ABS_CHECK)},
                                         "max_rel_err_by_depth": {str(e): F[:, ei, 0].max(0).tolist() for ei, e in enumerate(EPS_ABS_CHECK)},
                                         "min_cos_by_depth": {str(e): F[:, ei, 1].min(0).tolist() for ei, e in enumerate(EPS_ABS_CHECK)}}
        R["jvp_check"]["validated"] = bool(fok[0].mean() >= 0.95)
    else:
        R["jvp_check"]["validated"] = False
        R["jvp_check"]["note"] = "no float64 check in the processed cells"
    # ---------------- linearity
    LS = C["lin"]  # [cases, S, 7, 13]
    os_pass = np.all(LS[:, :, :4] <= tau, axis=2)
    c_pass = np.all(LS[:, :, 4:6] <= tau, axis=2)
    Gj = C["G"]
    Gp, Gm = C["G_fd"][:, :, 0], C["G_fd"][:, :, 1]  # [cases, S, 13]
    sca = ((np.abs(Gp - Gj[:, None]) <= tau * np.maximum(np.abs(Gj[:, None]), 0.1)) &
           (np.abs(Gm - Gj[:, None]) <= tau * np.maximum(np.abs(Gj[:, None]), 0.1)))
    names = ("err_plus", "err_minus", "slope_plus", "slope_minus", "err_central", "slope_central", "kappa")
    has = np.isfinite(LS[:, :, 0, 1])  # [cases, S]: the rung was run for this case
    lin = {}
    for si, e in enumerate(eps):
        h = has[:, si]
        if not h.any():
            lin[str(e)] = None
            continue
        lin[str(e)] = {"n_cases": int(h.sum()),
                       "one_sided_pass_rate_by_depth": os_pass[h, si].mean(0).tolist(),
                       "central_pass_rate_by_depth": c_pass[h, si].mean(0).tolist(),
                       "scalar_one_sided_pass_rate_by_depth": sca[h, si].mean(0).tolist(),
                       "median_by_depth": {nm: np.median(LS[h, si, j], 0).tolist() for j, nm in enumerate(names)},
                       "one_sided_pass_rate_by_direction_depth": {d: os_pass[h & (C["direction"] == di), si].mean(0).tolist()
                                                                  for di, d in enumerate(DIRECTIONS) if np.any(h & (C["direction"] == di))}}
    R["linearity"] = lin
    h1_rates = np.array(lin[str(eps[si0])]["one_sided_pass_rate_by_depth"][1:])
    H1 = bool(np.all(h1_rates >= 0.90))
    ok_eps = [e for si, e in enumerate(eps) if lin[str(e)] and np.all(np.array(lin[str(e)]["one_sided_pass_rate_by_depth"][1:]) >= 0.90)]
    R["linear_range"] = {"largest_eps_one_sided_pass_ge_0.90_at_every_depth": max(ok_eps) if ok_eps else None, "eps_ladder": eps}
    Guse = Gj if H1 else Gp[:, si0]
    R["gain_source"] = "JVP" if H1 else f"one-sided FD D+ at eps={eps[si0]} (H1 failed)"
    # ---------------- per subset (role x direction)
    curves, idx = {}, {}
    for r in roles:
        for di, d in enumerate(DIRECTIONS):
            ix = np.where((C["role"] == r) & (C["direction"] == di))[0]
            if ix.size < 2:
                continue
            idx[(r, d)] = ix
            sg = C["gsign"][ix][:, None] if d == "self" else np.ones((ix.size, 1))  # sign-align self pushes only
            G = Guse[ix] * sg
            X = {"G": G, "dG": np.diff(G, axis=1), "absG": np.abs(G), "P": C["P"][ix] * sg, "H": C["H"][ix] * sg,
                 "G_own": C["Gown"][ix] * sg, "G_spread": (Guse[ix] - C["Gown"][ix]) * sg, "norm_rel": C["Nrel"][ix],
                 "cos_with_dmu": C["cosmu"][ix] * sg, "own_token_norm_gain": C["own_norm_gain"][ix],
                 "head_minus_z_12": (C["H"][ix, 12:13] - Guse[ix, 12:13]) * sg, "head_per_protein": C["dY"][ix] * sg,
                 "logg": np.diff(C["logr"][ix], axis=1)}
            b = Boot(C["cell"][ix], C["donor"][ix], np.zeros(ix.size, int), a.n_boot, seed)
            blk = {"n_cases": int(ix.size), "n_cells": int(np.unique(C["cell"][ix]).size), **subset_summary(X, C["donor"][ix], b)}
            blk["ratio_of_medians"] = (np.array(blk["G"]["median"][1:]) / np.array(blk["G"]["median"][:-1])).tolist()
            hx = has[ix]
            blk["G_fd_plus_by_eps"] = {str(e): np.median((Gp[ix, si] * sg)[hx[:, si]], 0).tolist() if hx[:, si].any() else None
                                       for si, e in enumerate(eps)}
            blk["G_fd_minus_by_eps"] = {str(e): np.median((Gm[ix, si] * sg)[hx[:, si]], 0).tolist() if hx[:, si].any() else None
                                        for si, e in enumerate(eps)}
            blk["G_jvp_on_full_ladder_cases"] = np.median(G[hx[:, -1]], 0).tolist() if hx[:, -1].any() else None
            if d == "self":
                blk["sign_consistency_by_depth"] = (G > 0).mean(0).tolist()
            blk["by_cell_class"] = {k: {"n_cases": int((C["cls"][ix] == k).sum()), "G_median": np.median(G[C["cls"][ix] == k], 0).tolist()}
                                    for k in ("NK", "T") if (C["cls"][ix] == k).sum() >= 3}
            blk["by_gene"] = {GENE_SYMBOLS[g]: {"n_cases": int((C["gene"][ix] == g).sum()),
                                                "G_12_median": float(np.median(G[C["gene"][ix] == g, 12])),
                                                "H_12_median": float(np.median((C["H"][ix] * sg)[C["gene"][ix] == g, 12])),
                                                "median_rank_position": float(np.median(C["pos"][ix][C["gene"][ix] == g]))}
                              for g in np.unique(C["gene"][ix]) if g >= 0}
            curves[f"{r}/{d}"] = blk
    R["curves"] = curves
    # ---------------- hypotheses
    PC = probe_curves(a.probe_dir)
    R["probe_curves"] = PC
    main = "pair" if "pair" in roles else roles[0]
    hyp = {"JVP_check": {"validated": R["jvp_check"]["validated"],
                         "float64_pass_share": (R["jvp_check"].get("float64_cpu") or {}).get("pass_share"),
                         "float32_device_pass_share": R["jvp_check"]["pass_share"]},
           "H1_linearity": {"holds": H1, "eps": eps[si0], "min_one_sided_pass_rate": float(h1_rates.min()), "rates": h1_rates.tolist()}}
    cm = curves.get(f"{main}/nkt")
    if cm:
        dG = cm["dG"]
        L_e5 = int(np.argmin(dG["median"]) + 1)
        ub = dG["ci_two_stage"][L_e5 - 1][1]
        share = dG["argmin_share_two_stage"].get(L_e5, 0.0)
        pdn = {d_: v[L_e5 - 1] for d_, v in dG["per_donor"].items()}
        hyp["H2_localised_attenuation"] = {"holds": bool(ub < 0 and share >= 0.5 and _all_donors(pdn, lambda v: v < 0)),
                                           "L_E5": L_e5, "median_dG": dG["median"][L_e5 - 1],
                                           "ci_two_stage": dG["ci_two_stage"][L_e5 - 1], "ci_cell": dG["ci_cell"][L_e5 - 1],
                                           "argmin_share_two_stage": share, "argmin_share_cell": dG["argmin_share_cell"].get(L_e5, 0.0),
                                           "per_donor_dG": pdn, "gain_source": R["gain_source"]}
        hyp["H3_localisation_matches_probes"] = {"L_E5": L_e5, "L_probe": PC["L_probe_argmin_gap_ratio"],
                                                 "holds": abs(L_e5 - PC["L_probe_argmin_gap_ratio"]) <= 1,
                                                 "L_probe_largest_drop_gap_ratio": PC["L_probe_largest_drop_gap_ratio"],
                                                 "L_probe_largest_drop_centred_distance_ratio": PC["L_probe_largest_drop_distance_ratio"],
                                                 "spearman_G_vs_probe_gap_ratio_over_depths": spearman(cm["G"]["median"], PC["gap_ratio_mean_no_gdT158"]),
                                                 "spearman_G_vs_probe_order_rate_over_depths": spearman(cm["G"]["median"], PC["order_rate_no_gdT158"])}
        if ("pair", "nkt") in idx and ("random", "nkt") in idx:
            ia, ib = idx[("pair", "nkt")], idx[("random", "nkt")]
            u = np.concatenate([ia, ib])
            st = np.array(["a"] * ia.size + ["b"] * ib.size)
            g12 = Guse[u, 12]
            dd = joint_diff(g12, g12, C["cell"][u], C["donor"][u], st, seed, a.n_boot)
            c = dd["ci_two_stage"]
            hyp["H4_lookalike_specificity"] = {**dd, "holds": bool(c and c[1] < 0 and _all_donors(dd["per_donor"], lambda v: v < 0)),
                                               "opposite": bool(c and c[0] > 0 and _all_donors(dd["per_donor"], lambda v: v > 0))}
        hm = cm["head_minus_z_12"]
        c = hm["ci_two_stage"][0]
        pdh = {d_: v[0] for d_, v in hm["per_donor"].items()}
        hyp["H5_head_transmission"] = {"median_H12_minus_G12": hm["median"][0], "ci_two_stage": c, "ci_cell": hm["ci_cell"][0],
                                       "per_donor": pdh, "head_passes_less": bool(c[1] < 0 and _all_donors(pdh, lambda v: v < 0))}
    cs = curves.get(f"{main}/self")
    if cs:
        lstar = int(np.argmin(cs["logg"]["median"]) + 1)
        v3 = {"l_star": lstar, "median_log_gain_self": cs["logg"]["median"][lstar - 1],
              "comparator_layer_largest_drop_gap_ratio_no_gdT158": PC["L_probe_largest_drop_gap_ratio"],
              "comparator_layer_largest_drop_gap_ratio_all": PC["L_probe_largest_drop_gap_ratio_all"],
              "agreement_within_1_layer": abs(lstar - PC["L_probe_largest_drop_gap_ratio"]) <= 1}
        for ctl in ("rand", "other"):
            if (main, "self") in idx and (main, ctl) in idx:
                ia, ib = idx[(main, "self")], idx[(main, ctl)]
                u = np.concatenate([ia, ib])
                st = np.array(["a"] * ia.size + ["b"] * ib.size)
                lg = np.diff(C["logr"][u], axis=1)[:, lstar - 1]
                dd = joint_diff(lg, lg, C["cell"][u], C["donor"][u], st, seed, a.n_boot)
                c = dd["ci_two_stage"]
                if c and -0.1 < c[0] and c[1] < 0.1:
                    v_ = "no layer localises the NK-T loss (equivalent within 0.1)"
                elif c and dd["diff"] <= -0.1 and c[1] < 0 and _all_donors(dd["per_donor"], lambda v: v <= -0.05):
                    v_ = f"NK-T-specific attenuation at layer {lstar}"
                else:
                    v_ = "inconclusive"
                v3[f"specificity_vs_{ctl}"] = {**dd, "verdict": v_}
        hyp["E5_v3_endpoint"] = v3
    if not R["jvp_check"]["validated"]:
        verdict = "JVP not validated (experiments_v3 check failed): E5 stops; no endpoint is read"
    elif not H1:
        verdict = f"inconclusive for the linear-response reading (H1 failed at eps = {eps[si0]}); FD gains used for H2-H5"
    elif cm:
        h2 = hyp["H2_localised_attenuation"]["holds"]
        h4 = hyp.get("H4_lookalike_specificity", {}).get("holds", False)
        if not h2 and not h4:
            verdict = "C7 falsified for this model: no localised attenuation of the NK-T push (H2 and H4 do not hold)"
        else:
            verdict = "C7 supported: " + ", ".join(x for x, ok in ((f"a localised attenuation at layer {hyp['H2_localised_attenuation']['L_E5']} (H2)", h2),
                                                                     ("look-alike cells damp the push more (H4)", h4)) if ok)
    else:
        verdict = "no nkt cases"
    hyp["verdict"] = verdict if a.cell_pool == "site4" else f"(smoke on {a.cell_pool} cells; not a registered result) {verdict}"
    R["hypotheses"] = hyp
    write_json(out / "E5_results.json", R)
    write_report(out, R)
    say(f"wrote {out / 'E5_results.json'} and REPORT.md ({R['n_cases']} cases, {len(metas)} cells)")


def _f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


def _ci(c, n=2):
    return "" if c is None else f" [{_f(c[0], n)}, {_f(c[1], n)}]"


def write_report(out: Path, R: dict):
    L = []
    eps = REGISTRATION["eps"]
    PC = R["probe_curves"]
    H = R["hypotheses"]
    flag = f"SMOKE TEST: {R['smoke_note']}\n\n" if R["smoke_test"] else ""
    ck, jv = R["checks"], R["jvp_check"]
    L += ["# E5 Mode A layer dynamics (C7)", "", flag +
          f"{SCRIPT_VERSION}; cells {R['n_cells_done']} of {R['n_cells_selected']} ({R['cell_pool']}); cases {R['n_cases']}; "
          f"registration_sha256 {R['registration_sha256']}; addendum_sha256 {R['addendum_sha256']}; gains from {R['gain_source']}.", "",
          f"**Verdict:** {H['verdict']}", "",
          "JVP check (central FD on the unit tangent, rel. error <= 0.05 and cosine >= 0.99 at every depth): float64 CPU "
          + (", ".join(f"{e}: {_f(v, 3)}" for e, v in jv["float64_cpu"]["pass_share"].items()) + f" of {jv['float64_cpu']['n_cases']} cases"
             if jv.get("float64_cpu") else "not run") + f" -> validated {jv['validated']}; float32 device diagnostic "
          + ", ".join(f"{e}: {_f(v, 3)}" for e, v in jv["pass_share"].items()) + ". "
          f"Explicit layers vs module forward max |diff| {ck['explicit_vs_module_layer_means_max_abs']:.1e}; float64 CPU JVP vs "
          f"device JVP max rel. {_f(ck['cpu64_jvp_max_rel_diff'], 7) if ck['cpu64_jvp_max_rel_diff'] is not None else 'n/a'}; "
          f"vs official fp16 layer means max rel. {ck['vs_official_fp16_layer_means_max_rel']:.1e}; "
          f"median {_f(ck['sec_per_cell_median'], 1)} s per cell. Linear range (one-sided pass >= 0.90 at every depth): "
          f"eps <= {R['linear_range']['largest_eps_one_sided_pass_ge_0.90_at_every_depth']}.", ""]
    main = "pair" if "pair/nkt" in R["curves"] else next(iter(R["curves"])).split("/")[0]
    cm, cr = R["curves"].get(f"{main}/nkt"), R["curves"].get("random/nkt")
    if cm:
        L += [f"## NK-T push through the layers ({main} cells) next to the static probes", "",
              "G = NK-T component of the response (fraction of the class-mean distance, per unit of an aligned input push; "
              "G_0 = 1); dG = median per-case step change; r = ratio of medians (not used); H = head lens (depth 12 = the "
              "head's own response); own = share of G at the pushed token; pass = one-sided linearity pass rate. "
              "Intervals: two-stage donor/cell bootstrap.", "",
              "| depth | G | dG | r | G (random) | H | own | " + " | ".join(f"pass {e}" for e in eps)
              + " | probe gap ratio | probe order rate | centred dist. ratio |",
              "|---|---|---|---|---|---|---|" + "---|" * len(eps) + "---|---|---|"]
        g, dG = cm["G"], cm["dG"]
        for l_ in range(N_DEPTHS):
            own = cm["G_own"]["median"][l_] / g["median"][l_] if g["median"][l_] else None
            L.append(f"| {PC['depth_names'][l_]} | {_f(g['median'][l_])}{_ci(g['ci_two_stage'][l_])} | "
                     + (f"{_f(dG['median'][l_ - 1])}{_ci(dG['ci_two_stage'][l_ - 1])} | {_f(cm['ratio_of_medians'][l_ - 1], 2)}" if l_ else " | ")
                     + f" | {_f(cr['G']['median'][l_]) if cr else ''} | {_f(cm['H']['median'][l_])} | {_f(own, 2)} | "
                     + " | ".join(_f(R["linearity"][str(e)]["one_sided_pass_rate_by_depth"][l_], 2) if R["linearity"][str(e)] else "n/a" for e in eps)
                     + f" | {_f(PC['gap_ratio_mean_no_gdT158'][l_])} | {_f(PC['order_rate_no_gdT158'][l_])} | "
                     f"{_f(PC['centred_distance_ratio_no_gdT158'][l_])} |")
        L += ["", "Per donor (median G_12; median dG at L_E5):", ""]
        L2 = H.get("H2_localised_attenuation", {})
        for d_ in cm["G"]["per_donor"]:
            L.append(f"- donor {d_}: G_12 {_f(cm['G']['per_donor'][d_][12])}; dG at layer {L2.get('L_E5')} "
                     f"{_f(L2.get('per_donor_dG', {}).get(d_))}")
        jf = cm.get("G_jvp_on_full_ladder_cases")
        L += ["", f"## Gain vs push size ({main}/nkt, median G; eps {eps[0]} on every cell, the larger eps on the "
              f"first {N_FULL_LADDER} cells of the seeded order)", "", "| push | G_3 | G_6 | G_9 | G_12 |", "|---|---|---|---|---|",
              "| JVP, all cases | " + " | ".join(_f(g["median"][l_]) for l_ in (3, 6, 9, 12)) + " |",
              "| JVP, full-ladder cases | " + (" | ".join(_f(jf[l_]) for l_ in (3, 6, 9, 12)) if jf else "n/a | | |") + " |"]
        for e in eps:
            for nm, k_ in (("+", "G_fd_plus_by_eps"), ("-", "G_fd_minus_by_eps")):
                v = cm[k_][str(e)]
                L.append(f"| {nm}{e} | " + (" | ".join(_f(v[l_]) for l_ in (3, 6, 9, 12)) if v else "n/a | | |") + " |")
    L += ["", "## Registered hypotheses", ""]
    for k, v in H.items():
        if k != "verdict":
            L.append(f"- **{k}**: `" + json.dumps(v, default=float)[:1500] + "`")
    L += ["", "## All subsets (median G at depths 0, 6, 9, 12; H_12)", "",
          "| subset | cases | cells | G_0 | G_6 | G_9 | G_12 | H_12 | self sign consistency at 12 |", "|---|---|---|---|---|---|---|---|---|"]
    for key_, blk in R["curves"].items():
        gm = blk["G"]["median"]
        sc = blk.get("sign_consistency_by_depth")
        L.append(f"| {key_} | {blk['n_cases']} | {blk['n_cells']} | {_f(gm[0])} | {_f(gm[6])} | {_f(gm[9])} | {_f(gm[12])} | "
                 f"{_f(blk['H']['median'][12])} | {_f(sc[12], 2) if sc else ''} |")
    if cm:
        L += ["", f"## Per gene ({main}/nkt)", "", "| gene | cases | median rank | G_12 | H_12 |", "|---|---|---|---|---|"]
        for s_, v in cm["by_gene"].items():
            L.append(f"| {s_} | {v['n_cases']} | {v['median_rank_position']:.0f} | {_f(v['G_12_median'])} | {_f(v['H_12_median'])} |")
    (out / "REPORT.md").write_text("\n".join(L) + "\n")


# ============================================================================ main

def main(argv=None) -> int:
    a = parse_args(argv)
    if a.stage == "register":
        stage_register(a)
        return 0
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG["path"] = a.out_dir / "progress.log"
    say(f"{SCRIPT_VERSION} stage {a.stage} out {a.out_dir}" + (f" SMOKE: {a.smoke}" if a.smoke else ""))
    reg_info = check_registration(a)
    write_json(a.out_dir / "e5_registration_used.json", {**reg_info, "registration": REGISTRATION})
    need_rna = a.stage in ("directions", "respond", "all")
    meta = lp.load_meta(a.processed, need_rna=need_rna)
    key = fkey(meta)
    if a.stage in ("select", "directions", "respond", "all"):
        emb = lp.Embedding(a.embed_dir, None, len(meta["split"]))
    if a.stage in ("select", "all"):
        if (a.out_dir / "e5_cells.npz").exists() and a.stage == "all":
            say("reusing e5_cells.npz (the selection is deterministic; delete the file to rebuild it)")
        else:
            stage_select(a, meta, emb, key)
    if a.stage in ("directions", "all"):
        dj = a.out_dir / "e5_directions.json"
        if (a.out_dir / "e5_directions.npz").exists() and dj.exists() and a.stage == "all":
            if json.loads(dj.read_text()).get("subsampled") and not a.smoke:
                raise SystemExit("e5_directions.npz in this out-dir was built from a smoke subsample; delete it first")
            say("reusing e5_directions.npz")
        else:
            stage_directions(a, meta, emb, key)
    if a.stage in ("respond", "all"):
        stage_respond(a, meta, emb, load_directions(a.out_dir))
    if a.stage in ("report", "all"):
        stage_report(a, meta, reg_info)
    return 0


if __name__ == "__main__":
    sys.exit(main())
