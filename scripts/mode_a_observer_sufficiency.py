#!/usr/bin/env python3
"""E5-M (Mode A): observer-conditioned sufficiency test of TEDDY's gene-mean state.

Question. TEDDY's users and our phase-1 head read the gene-mean of the layer-12 token states, z_12. Does
z_12 retain what an INDEPENDENT observer of the token states uses, on the NK-T look-alike cells where the head
loses the NK-T difference? ANM v2 framing: z_12 is a DECLARED candidate retained state tested against DECLARED
readouts. With m source columns and d = 512, K_z = dz_12/du has full column rank, so first-order kernel inclusion
(ker K_z within ker K_O) is automatic, i.e. NOT_INFORMATIVE; the informative first-order test is the residual of
a map H fixed WITHOUT the tested responses: r = ||K_O - H K_z||_F / ||K_O||_F (ANM paper SI). The readout-visible
quotient C_T R_T (rank <= readout dimension) is reported as a diagnostic only; using it as the state would make
the test automatic.

    sources u_j   unit pushes of the input state h_0 at the position of a declared NK/T gene token present in
                  the cell, along (self) the gene's own embedding, (nkt) the training NK-T direction of the input
                  gene-mean (E5's u_nkt), (rand) one fixed random unit vector; m = 3 x (declared genes present)
    K_z           d x m: dz_12 per column (JVP through every layer, explicit attention as E5)
    observers     O_tok  = E3's frozen R2 (one learned query over layer-12 token states; training cells only)
                  O_head = the phase-1 head on z_12 (a function of z by construction: positive control)
                  O_rand = R2's form, random query (attention-logit spread matched to R2's on training states) and
                           random output rows (norms matched to R2's): negative control
                  outputs: CD56, CD94, CD335, CD3 (E3's gap proteins = E5's), each divided by its training SD
    H (primary)   ANM's clamp estimator H = Lambda pinv(Gamma), Gamma = I (clamps of z_12: every layer-12 token moved
                  by delta, the pooling-discarded deviations h_t - z_12 held), Lambda = the observer's clamp response;
                  no fitted number, never uses K_O. K_O - H K_z = the observer's response to the part of dH_12 that
                  pooling discards, so r = 0 for any function of z (O_head: positive control)
    H (second.)   frozen training maps (E3's 4,000 R2 training cells; registration/addenda/E5M_H.json): ridge of the
                  outputs on z_12 (lambda by val) and the training-average clamp map
    primary       r_O per cell; Delta_O = median r_O(pair cells) - median r_O(matched random cells)

Stages
  register         CPU, train/val only: the training maps + standardisation (E5M_H.json), leakage check (site4 rows poisoned -> identical
                   core; val rows poisoned -> different core), write registration/addenda/E5M.json + HASHES lines
  commit-addendum  git: commit exactly E5M.json, E5M_H.json and their HASHES.txt lines (temporary index)
  select           the cells (E5's registered set, read from E5's out-dir or rebuilt; no forward pass)
  respond          device: per cell one JVP per source column + FD checks; one file per cell (resumable)
  report           CPU: clamp-H residuals, FD gate, positive control, two-stage bootstrap, H5M verdict, diagnostics
  all              select, respond, report
Any site4 respond / report is refused unless registration_v3.json, every amendment (A1, A2, A3 if present) and
the E5M addendum + H file are committed and their hashes match. Smoke runs use val donor 18303 (--cell-pool val).

Full run (from the repo root; R = the main checkout):
  PY=<venv>/bin/python; R=/Users/tianchichen/Documents/GitHub/teddy_mm
  ANM_ROOT=<anm_v2_fix> OMP_NUM_THREADS=4 PYTORCH_ENABLE_MPS_FALLBACK=1 $PY scripts/mode_a_observer_sufficiency.py \
      --stage all --out-dir $R/outputs/v3/E5M --device auto --threads 4 --jvp-chunk 8 --fd-batch 4
(add --max-minutes N to pause after N minutes with exit code 75; the same command resumes)
Every stage appends to <out-dir>/progress.log; respond logs seconds per cell and an ETA after every cell (also in
e5m_progress.json); rerunning the same command skips finished cells (each cell file carries the addendum hash).
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
import inspect  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from types import SimpleNamespace  # noqa: E402

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
from lib import v3_e5m as e5m  # noqa: E402

SCRIPT_VERSION = "mode_a_observer_sufficiency 1.0"
MAIN = Path("/Users/tianchichen/Documents/GitHub/teddy_mm")
_BASE = MAIN if MAIN.exists() else ROOT
PROTEINS = ("CD56", "CD94", "CD335", "CD3")  # E3's gap proteins (= E5's PROTEINS), the registered NK/T readout panel
assert tuple(mld.PROTEINS) == PROTEINS
NKT_GENES = tuple(mld.NKT_GENES)  # CD3E, CD3D, NCAM1, KLRD1, NCR1, FCGR3A, KLRF1 (E5's declared NK/T push genes)
DIRECTIONS = ("self", "nkt", "rand")
OBSERVERS = ("tok", "head", "rand")
N_DEPTHS = 13
EPS_FD = (0.01, 0.02)  # E5's check rung 1e-2 and its primary ladder rung 0.02 (1e-3 hits the float32 floor, E5 addendum)
FD_REL_TOL, FD_COS_TOL, FD_PASS_SHARE = 0.05, 0.99, 0.95
N_FD_ALL_COLUMNS = 10  # the first cells of the processing order get the FD check on every column
N_FD_CELLS = 100  # the first cells of the processing order get the FD check (3 columns of their first present gene)
LAMBDAS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)
MARGIN = 0.05
N_BOOT = 2000
NULL_TOL = 0.05
NULL_TOL_REPORTED = (0.05, 0.01)
QUOTIENT_RTOL = 0.05
ANM_BUDGETS = ((0.0, 1e-4), (0.0, 0.05))  # (atol, rtol): float32-JVP numerical budget; FD-gate measurement budget
N_PAIR_PER_DONOR = 60  # E5's pair cells used: per primary donor, the first 60 of E5's seeded processing order
N_RAND_PER_DONOR = 40  # E5's random cells used: per primary donor, the first 40 of E5's seeded processing order
EXIT_INCOMPLETE = 75
SEED_TAG = 5  # numpy default_rng([seeds.e5_random_directions, 5]) draws v_rand, then q_rand, then W_rand
N_LAYER_CELLS = 30  # secondary per-layer clamp maps (reverse mode) on the first cells of the processing order
POSITIVE_CONTROL_TOL = 1e-3  # max r_head with the clamp map (exactly 0 up to float32 rounding)
ADDENDUM = "E5M.json"
H_ARRAYS = ("H_ridge", "ridge_intercept", "ridge_lambda", "out_mean", "out_sd", "H_train_clamp", "v_rand", "q_rand",
            "W_rand", "u_nkt")
H_FILE = "E5M_H.json"
_LOG = {"path": None}


def say(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG["path"] is not None:
        with open(_LOG["path"], "a") as f:
            f.write(line + "\n")


DECLARED = {
    "question": "is TEDDY's gene-mean layer-12 state z_12 (what TEDDY's users and our head read) sufficient, at first "
                "order, for an independent observer of the layer-12 token states on the NK-T look-alike cells where the "
                "head loses the NK-T difference?",
    "anm_framing": "z_12 is a DECLARED candidate retained state tested against DECLARED readouts. With m <= 21 source "
                   "columns and d = 512, K_z has full column rank, so first-order kernel inclusion is automatic and "
                   "NOT_INFORMATIVE (any K_O equals (K_O pinv K_z) K_z); the informative test is the residual of an H "
                   "fixed without the tested responses (ANM paper SI; referee audit 2026-09-16). Primary H = ANM's own "
                   "estimator H = Lambda pinv(Gamma) from clamps of the candidate state alone (analysis.py: 'Never use Ko "
                   "or held-out endpoints to fit H'). The readout-visible quotient R_T/(R_T cap ker C_T) = C_T R_T has "
                   "rank <= the readout dimension (4) and is reported only as a diagnostic; it is never used as the state.",
    "cells": {
        "site4": "E5's registered cell set (registration E5.cells, E5 addendum 'cells': unique test_primary cells of E3's "
                 "within-donor NK-T pairs under the v3 primary key, <= 200 per donor; E5 secondary 'random_cells': per "
                 "donor 100 matched test_primary NK/T cells in no pair), read from E5's e5_cells.npz (or rebuilt with "
                 f"E5's own select code). Used (GPU budget): per primary donor the first {N_PAIR_PER_DONOR} pair cells and "
                 f"the first {N_RAND_PER_DONOR} random cells of E5's seeded processing order (a seeded random subset; "
                 f"{2 * (N_PAIR_PER_DONOR + N_RAND_PER_DONOR)} cells). Processing order = E5's order restricted to these "
                 "cells. A cell without any declared gene present has no column and is excluded (counted).",
        "smoke": "val donor 18303 primary-key NK/T cells only (seeded draw; roles alternate pair / random so every code "
                 "path runs); never a verdict",
    },
    "sources": {
        "genes": "the 7 declared NK/T genes of E5's nkt push, in this order: " + ", ".join(NKT_GENES) + " (Ensembl ids "
                 "by E5's gene reference; tokens present in the cell's 2,048-token input)",
        "directions": "per present gene, 3 unit pushes of the input state h_0 at the gene's position: self = E[gene] / "
                      "||E[gene]||; nkt = E5's u_nkt (unit mean_NK - mean_T of the training input gene-mean, "
                      "outputs/v3/E5/e5_directions.npz, training cells only; copied into E5M_H.json); rand = one fixed "
                      "Gaussian unit vector (numpy default_rng([seeds.e5_random_directions, 5]), first draw) for every gene "
                      "and cell",
        "columns": "m = 3 x (genes present), column order gene-major (declared gene order), then self, nkt, rand; unit "
                   "input size for every column (eps cancels in the residual)",
    },
    "candidate_state": {
        "primary": "z_12 = gene-mean of the layer-12 token states (float32 explicit forward of E5)",
        "secondary": "z_l = gene-mean of layer l, l = 0..11 (per-layer clamp maps, S1)",
    },
    "observers": {
        "outputs": "CD56, CD94, CD335, CD3 (E3's gap proteins = E5's head proteins), each divided by its SD over the "
                   "calibration training cells (each protein carries equal weight; fixed in E5M_H.json)",
        "O_tok": "E3's frozen R2 (outputs/v3/E3/r2/r2.pt; one learned query over the layer-12 token states, linear "
                 "output; trained on 4,000 training cells, early-stopped on 1,000 val cells) on the float32 layer-12 "
                 "states of the explicit forward (no fp16 rounding: its derivative is zero almost everywhere)",
        "O_head": "the phase-1 head on z_12 (E5's load_head: L2-normalised z, MLP + NB decoder mean at the training "
                  "median size factor, / training p95): a function of z by construction (positive control)",
        "O_rand": "R2's form (R2's mu / sd standardisation, one query, linear output) with a Gaussian query direction "
                  "scaled so that the median over the 4,000 calibration training cells of the per-cell SD of its attention "
                  "logits equals R2's (matched attention sharpness; a query of R2's plain norm attends almost uniformly: "
                  "val smoke r about 0.03) and Gaussian output rows with the norms of R2's rows for the 4 proteins, bias "
                  "0; default_rng([seeds.e5_random_directions, 5]) after v_rand: q then W (negative control: a token "
                  "observer with no NK-T training)",
    },
    "derivatives": {
        "K": "per cell one forward-mode JVP per column through every layer (E5's layer_fn, explicit attention with "
             "each layer's weights; forward AD does not run through fused attention): K_z (d x m) at every depth, the "
             "layer-12 token-state response dH_12, K_O for O_tok and O_rand by JVP of the observer at H_12 along dH_12, "
             "K_O for O_head by JVP of the head at z_12 along dz_12; float32 on the device",
        "fd_gate": f"central finite differences through the module's own layers at eps in {list(EPS_FD)} (unit input "
                   f"push) on the 3 columns of the first present gene of the first {N_FD_CELLS} cells of the order, and on "
                   f"every column of the first {N_FD_ALL_COLUMNS}; targets z_12 (raw), O_tok, O_head, O_rand (standardised); "
                   "the gate derivative is the Richardson extrapolation (4 D(0.01) - D(0.02)) / 3 of the two central "
                   f"differences; pass = relative error <= {FD_REL_TOL} and cosine >= {FD_COS_TOL}; validated if the pass "
                   f"share is >= {FD_PASS_SHARE} for every target; each eps alone is reported; otherwise 'JVP not "
                   "validated' and no verdict",
        "pipeline_checks": "explicit forward vs module forward layer means (every cell); per-cell clamp maps of O_tok and "
                           "O_rand equal to the frozen training-average clamp maps (residuals agree); clamp residual of O_tok also "
                           "computed directly as the JVP of R2 along the pooling-discarded part dH_12 - 1 dz_12' "
                           "(must equal the H_clamp residual); same-cell fitted map K_O pinv(K_z) (~0: why kernel "
                           "inclusion is automatic)",
    },
    "H": {
        "primary_clamp": "per cell and observer, ANM's clamp estimator H = Lambda pinv(Gamma) with Gamma = I_512 (clamps "
                         "of every coordinate of z_12) and Lambda = the observer's response to each clamp: the clamp of "
                         "z_12 by delta moves every layer-12 token state by delta, i.e. changes z_12 by delta and leaves "
                         "the pooling-discarded deviations h_t - z_12 unchanged (the orthogonal complement of the kernel "
                         "of gene-mean pooling). Lambda = d O(H_12 + 1 delta') / d delta at 0 (reverse mode through the "
                         "observer only). Computed before and independently of the source responses; no fitted number. "
                         "Then K_O - H K_z = the observer's response to the pooling-discarded part dH_12 - 1 dz_12' of "
                         "the token-state change, and r = its share of the observer's response. For O_head (a function "
                         "of z) H is the head's Jacobian and r = 0 by construction: the positive control (gate: max r_head "
                         f"<= {POSITIVE_CONTROL_TOL}). For R2-form observers (O_tok, O_rand) the clamp moves every attention "
                         "logit by the same amount, so the attention weights do not change and Lambda = W diag(1/sd) "
                         "exactly: the same map for every cell, equal to the training-average clamp map frozen in "
                         "E5M_H.json (pipeline check). r_tok is then the share of R2's response that comes from its "
                         "attention departing from the gene-mean",
        "secondary_training": "frozen 4 x 512 maps calibrated on E3's 4,000 R2 training cells (training split; E3 addendum "
                              "r2_cells, hashes checked; states = E3's stored float16 layer-12 states, z = z_rna.npy, "
                              "stored in registration/addenda/E5M_H.json: (a) ridge of the "
                              "standardised outputs on raw z_12, lambda (relative to trace(Xc'Xc)/d) from "
                              f"{list(LAMBDAS)} by MSE on E3's 1,000 R2 val cells (the task's example; on val "
                              "smoke cells its O_head residual was 1.7-2.0, i.e. it fails the positive control: reported "
                              "only); (b) the training-average clamp map (mean over the 4,000 cells of each observer's "
                              "Lambda at depth 12)",
        "why_clamp_primary": "fixed on val before site4: a global linear map cannot represent a nonlinear observer's local "
                             "derivative (val smoke: ridge residual of O_head, a function of z, 1.7-2.0), so the training "
                             "map's residual mixes nonlinearity with pooling loss; the clamp map removes the nonlinearity "
                             "exactly (r_head = 0) and isolates the pooling-discarded component",
    },
    "statistics": {
        "residual": "r_O = ||K_O - H_clamp,O K_z||_F / ||K_O||_F per cell (K_O rows standardised), depth 12",
        "contrasts": "Delta_O = median r_O(pair cells) - median r_O(random cells), O in tok, rand (and head, ~0); "
                     "DD_rand = Delta_tok - Delta_rand; level = median r_tok(pair); reported: median r_rand(pair)",
        "bootstrap": "registration common.statistics: two-stage (donors with replacement, then cells with replacement "
                     f"within each drawn donor and role), B = {N_BOOT}, seed = seeds.bootstrap, percentile 95% intervals, "
                     "the same draws for every statistic; per-donor values on each primary donor's own cells",
        "margin": MARGIN,
    },
    "hypothesis": {
        "H5M": "on NK-T pair cells the clamp-H residual of O_tok exceeds that on matched random cells: Delta_tok >= "
               f"{MARGIN} with the common win rule (two-stage lower bound > 0, >= {MARGIN / 2} in each primary donor), and "
               "it exceeds the random-observer null: DD_rand has a two-stage lower bound > 0",
        "falsification": f"H5M is falsified if Delta_tok is 'equivalent' (two-stage 95% interval inside (-{MARGIN}, "
                         f"{MARGIN})) or a common-rule 'loss'",
        "verdicts": {
            "REJECT_SUFFICIENCY": "H5M supported: pooling to z_12 discards information this independent observer uses "
                                  "on look-alike cells, more than on matched random cells and beyond a random token observer",
            "PASS_SUFFICIENT": "H5M falsified by equivalence AND level two-stage upper bound < margin (on pair cells less "
                               "than 5% of O_tok's response norm flows through the pooling-discarded part; the clamp "
                               "residual is an exact decomposition, validated by the FD gate and the positive control): "
                               "z_12 suffices for this observer at first order, so the NK-T loss is in our head, not in "
                               "pooling",
            "H5M_FALSIFIED": "H5M falsified, but the level condition fails: no look-alike-specific excess, yet a share of "
                             "O_tok's response above the margin flows through what pooling discards on pair and random "
                             "cells alike (z_12 not sufficient for this observer, not specifically on look-alikes)",
            "INCONCLUSIVE": "anything else (including a Delta_tok win that fails the O_rand null)",
            "NOT_VALIDATED": "FD gate or positive control failed: E5-M stops",
        },
        "scope": "first order (JVP) at the cells' own inputs, the declared sources and observers only; two primary "
                 "donors (the interval reflects the observed donors); a PASS says nothing about other observers or "
                 "larger pushes",
    },
    "secondary": {
        "note": "pre-specified, reported separately, never change the verdict",
        "S1_layers": f"clamp residual r_O(z_l) for l = 0..12 on the first {N_LAYER_CELLS} cells of the order (Lambda_l by "
                     "reverse mode from the observer to the layer-l token states, summed over tokens: clamping z_l moves "
                     "every layer-l token, then the later layers act; even O_head is not a function of z_l for l < 12)",
        "S2_log": "Delta on log r (scale-free)",
        "S3_directions": "r_O on the self-, nkt- and rand-column subsets",
        "S4_anm_status": "ANM analysis.analyze_factorization(K_z, K_O, Gamma = I_d, Lambda = H_clamp,O) per cell and "
                         f"observer at (atol, rtol) in {[list(b) for b in ANM_BUDGETS]}: counts of overall, source-kernel "
                         "and factorization statuses (source kernel expected NOT_INFORMATIVE when K_z has full column rank)",
        "S5_quotient": f"rank (singular values > {QUOTIENT_RTOL} x max) and participation ratio of K_O (diagnostic only, "
                       "<= 4)",
        "S6_near_null": f"share of ||K_O||_F^2 in the near-null right-singular space of K_z (sigma <= tol x sigma_max), "
                        f"registered tol {NULL_TOL} (also {NULL_TOL_REPORTED[1]})",
        "S7_training_maps": "r_O and the contrasts with the frozen training maps (ridge, training-average clamp) at depth 12",
        "S8_same_cell_fit": "r_tok with H fitted on the cell itself (automatic, ~0)",
    },
    "compute": {
        "device": "MPS float32 (CPU allowed; same arithmetic)",
        "cells": f"{N_PAIR_PER_DONOR} pair + {N_RAND_PER_DONOR} random E5 cells per primary donor",
        "budget": "declared GPU budget <= 2 h at the contention measured when registering (calibration.compute); the run "
                  "is resumable and --max-minutes only pauses it (the declared cells are never changed by time)",
    },
    "fairness": "every number fixed here comes from training / validation cells only (calibration: E3's R2 train / val "
                "subset; directions: E5's training u_nkt); no site4 label, protein, cell type, embedding or RNA "
                "informed any choice; smoke runs used val donor 18303 only",
}


# ============================================================================ args
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("register", "commit-addendum", "select", "respond", "report", "all"), default="all")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=ROOT / "data/processed/cite_official")
    p.add_argument("--ckpt", type=Path, default=_BASE.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=_BASE / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--e3-dir", type=Path, default=_BASE / "outputs/v3/E3", help="E3 out-dir (r2/r2.pt, r2/states)")
    p.add_argument("--e5-dir", type=Path, default=_BASE / "outputs/v3/E5", help="E5 out-dir (e5_cells.npz, e5_directions.npz)")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--gene-reference", type=Path, default=ROOT / mld.GENE_REFERENCE)
    p.add_argument("--anm-root", type=Path, default=Path(os.environ["ANM_ROOT"]) if os.environ.get("ANM_ROOT") else None)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--cell-pool", choices=("site4", "val"), default="site4")
    p.add_argument("--device", default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--jvp-chunk", type=int, default=8)
    p.add_argument("--fd-batch", type=int, default=4)
    p.add_argument("--limit", type=int, default=0, help="smoke only: first N cells of the order")
    p.add_argument("--n-val-cells", type=int, default=6, help="smoke only: val cells per class drawn for the val pool")
    p.add_argument("--n-boot", type=int, default=N_BOOT)
    p.add_argument("--max-minutes", type=float, default=0, help="respond: pause after this many minutes (exit 75; rerun resumes)")
    p.add_argument("--calib-limit", type=int, default=0, help="smoke only: calibrate on the first N train / val cells")
    p.add_argument("--no-extras", action="store_true",
                   help="smoke / timing only: no all-column FD and no per-layer clamps (a typical cell's cost)")
    p.add_argument("--compute-note", default=None, help="register: the measured compute estimate recorded in the addendum")
    p.add_argument("--smoke", default=None, metavar="NOTE")
    a = p.parse_args(argv)
    if a.stage in ("select", "respond", "report", "all") and a.out_dir is None:
        p.error("--out-dir is required")
    if (a.limit or a.n_boot != N_BOOT or a.calib_limit or a.cell_pool == "val" or a.no_extras) and not a.smoke:
        p.error("--limit / --n-boot / --calib-limit / --cell-pool val / --no-extras are for smoke runs only (give --smoke NOTE)")
    return a


# ============================================================================ small helpers
def sha256_file(path: Path) -> str:
    return vk.sha256_file(path)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def write_json(path: Path, obj) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default))
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def atomic_savez(path: Path, **arrays) -> None:
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def git(*args, env=None, input_bytes: bytes | None = None) -> str:
    r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, env=env, input=input_bytes)
    if r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.decode(errors='replace')}")
    return r.stdout.decode()


def git_committed(path: Path) -> bool | None:
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


def code_hashes() -> dict:
    """sha256 of the source of every E5 function reused here (E5's file may change for unrelated reasons)."""
    fns = {"layer_fn": mld.layer_fn, "tokenise": mld.tokenise, "load_head": mld.load_head, "e2_genes": mld.e2_genes,
           "stage_select": mld.stage_select, "primary_keys": mld.primary_keys}
    return {k: sha256_bytes(inspect.getsource(f).encode()) for k, f in fns.items()}


def load_anm(a):
    root = a.anm_root
    if root is None:
        return None, None
    f = Path(root) / "benchmarks" / "paper1_retained_response" / "analysis.py"
    if not f.exists():
        return None, None
    spec = importlib.util.spec_from_file_location("anm_paper1_analysis", f)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m, {"file": str(f), "sha256": sha256_file(f)}


def load_reg(a) -> dict:
    """The registration as every v3 builder reads it: hash-checked registration, then A1, A2 and A3 (A3 if the
    loader requires it); E5-M reads only seeds, splits and E2's gene list from it (none amended)."""
    return va.load_registration_amended(a.registration_dir / "registration_v3.json",
                                        a.registration_dir / "amendment_A1.json")


def addenda_dir(a) -> Path:
    return a.registration_dir / "addenda"


# ============================================================================ observers
def r2_index(targets) -> list[int]:
    t = list(targets)
    return [t.index(p) for p in PROTEINS]


def load_observers(a, meta, device, Hz: dict | None = None):
    """(R2 module, idx, O_rand module, head yhat, info). O_rand parameters come from Hz (the frozen file) when
    given, otherwise they are drawn (register)."""
    import torch

    from lib.v3_e3_readouts import R2AttnPool
    blob = torch.load(str(a.e3_dir / "r2" / "r2.pt"), map_location="cpu", weights_only=False)
    r2 = R2AttnPool(int(blob["d"]), int(blob["d_out"]))
    r2.load_state_dict(blob["state_dict"])
    idx = r2_index(blob["targets"])
    if Hz is None:
        rp = rand_params(r2, idx, a.reg_seed)
    else:
        rp = {"v_rand": Hz["v_rand"], "q_rand": Hz["q_rand"], "W_rand": Hz["W_rand"]}
    rr = R2AttnPool(int(blob["d"]), len(PROTEINS))
    with torch.no_grad():
        rr.mu.copy_(r2.mu)
        rr.sd.copy_(r2.sd)
        rr.q.copy_(torch.as_tensor(rp["q_rand"], dtype=torch.float32))
        rr.out.weight.copy_(torch.as_tensor(rp["W_rand"], dtype=torch.float32))
        rr.out.bias.zero_()
    for m in (r2, rr):
        m.eval().to(device)
        for p_ in m.parameters():
            p_.requires_grad_(False)
    yhat, _, hinfo = mld.load_head(a.head_ckpt, meta, device)
    info = {"r2": str(a.e3_dir / "r2" / "r2.pt"), "r2_sha256": sha256_file(a.e3_dir / "r2" / "r2.pt"),
            "r2_targets": list(blob["targets"]), "r2_index": idx, "r2_train_sha256": blob.get("train_sha256"),
            "r2_val_sha256": blob.get("val_sha256"), "head": hinfo, "head_sha256": sha256_file(a.head_ckpt)}
    return r2, idx, rr, yhat, rp, info


def rand_params(r2, idx, seed) -> dict:
    """v_rand (source direction), q_rand (unit direction here; register rescales it to R2's attention-logit spread on
    the training states), W_rand (O_rand) from default_rng(seed), in this order."""
    rng = np.random.default_rng(seed)
    d = int(r2.q.numel())
    v = rng.standard_normal(d)
    v /= np.linalg.norm(v)
    q = rng.standard_normal(d)
    q /= np.linalg.norm(q)
    W = rng.standard_normal((len(idx), d))
    rows = r2.out.weight.detach().cpu().numpy()[idx].astype(np.float64)
    W *= (np.linalg.norm(rows, axis=1) / np.linalg.norm(W, axis=1))[:, None]
    return {"v_rand": v, "q_rand": q, "W_rand": W}


# ============================================================================ calibration (train / val only)
class StateStore:
    """E3's stored float16 layer-12 states (outputs/v3/E3/r2/states, written by v3_e3_nkt_repair r2_states)."""

    def __init__(self, sd: Path):
        self.views, self.shard_addendum = {}, set()
        for mf in sorted(sd.glob("shard_*_meta.npz")):
            k = mf.name.split("_")[1]
            with np.load(mf) as m:
                cells, offs = m["cells"], m["offsets"]
                self.shard_addendum.add(str(m["addendum_sha256"]))
            st = np.load(sd / f"shard_{k}_states.npy", mmap_mode="r")
            for i, c in enumerate(cells):
                self.views[int(c)] = (st, int(offs[i]), int(offs[i + 1]))

    def get(self, c: int) -> np.ndarray:
        st, s, e = self.views[int(c)]
        return st[s:e]


class RowView:
    """Row access to a [N, d] (or [n_layers, N, d] via layer) array whose chosen rows are replaced by seeded
    random values: the leakage check passes this in place of the real array."""

    def __init__(self, arr, poison_rows: np.ndarray | None = None, seed: int = 0, layer: int | None = None):
        self.arr, self.layer = arr, layer
        self.poison = set(int(x) for x in poison_rows) if poison_rows is not None else set()
        self.seed = seed

    def rows(self, ix: np.ndarray) -> np.ndarray:
        ix = np.asarray(ix, np.int64)
        X = np.asarray(self.arr[self.layer][ix] if self.layer is not None else self.arr[ix], np.float64)
        hit = [k for k, i in enumerate(ix) if int(i) in self.poison]
        if hit:
            rng = np.random.default_rng([self.seed, 0 if self.layer is None else self.layer + 1])
            X[hit] = rng.standard_normal((len(hit), X.shape[1]))
        return X


def clamp_map_tokens(obs, Hs):
    """Lambda = d O(Hs + 1 delta') / d delta at delta = 0 for a token observer obs(Hs [1, L, d]) -> [p]: the response
    to clamping the gene-mean (every token moved by delta; deviations from the mean held). Returns [p, d]."""
    import torch
    d = Hs.shape[-1]
    return torch.func.jacrev(lambda dl: obs(Hs + dl[None, None, :]))(torch.zeros(d, dtype=Hs.dtype, device=Hs.device))


def token_observer_outputs(a, meta, tr: np.ndarray, va: np.ndarray, seed, zview: "RowView") -> tuple[dict, dict, dict]:
    """O_tok, O_rand on the calibration cells (train then val) from E3's stored states (which hold only E3's R2
    train / val cells), their clamp maps summed over the training cells, and the gene-mean of those states against z."""
    import torch
    torch.set_num_threads(_THREADS)
    a.reg_seed = seed
    r2, idx, rr, _, rp, oinfo = load_observers(a, meta, torch.device("cpu"))
    S = StateStore(a.e3_dir / "r2" / "states")
    cells = np.concatenate([tr, va])
    miss = [int(c) for c in cells if int(c) not in S.views]
    if miss:
        raise SystemExit(f"{len(miss)} calibration cells have no stored E3 states")
    zoff = zview.rows(cells)
    # O_rand's query: random direction, scaled so that the median over the training cells of the per-cell SD of its
    # attention logits <q, h~_t> / sqrt(d) equals R2's (matched attention sharpness)
    with torch.no_grad():
        qh = torch.as_tensor(rp["q_rand"], dtype=torch.float32)
        sd_r2, sd_q = np.zeros(tr.size), np.zeros(tr.size)
        for i, c in enumerate(tr):
            ht = (torch.from_numpy(np.asarray(S.get(c), dtype=np.float32)) - r2.mu) / r2.sd
            sd_r2[i] = float((ht @ r2.q).std() / np.sqrt(ht.shape[1]))
            sd_q[i] = float((ht @ qh).std() / np.sqrt(ht.shape[1]))
        scale = float(np.median(sd_r2) / np.median(sd_q))
        rp["q_rand"] = rp["q_rand"] * scale
        rr.q.copy_(torch.as_tensor(rp["q_rand"], dtype=torch.float32))
    oinfo["rand_query"] = {"median_logit_sd_r2": float(np.median(sd_r2)), "median_logit_sd_unit_random": float(np.median(sd_q)),
                           "scale": scale, "norm": float(np.linalg.norm(rp["q_rand"])), "r2_query_norm": float(r2.q.norm())}
    Oc = {k: np.zeros((cells.size, len(PROTEINS))) for k in ("tok", "rand")}
    Lsum = {k: np.zeros((len(PROTEINS), zoff.shape[1])) for k in ("tok", "rand")}
    cos = np.zeros(cells.size)
    idx_t = torch.tensor(idx)
    t0 = time.time()

    def o_tok(Hs):
        return r2(Hs, torch.ones(Hs.shape[0], Hs.shape[1], dtype=torch.bool))[0][0].index_select(-1, idx_t)

    def o_rand(Hs):
        return rr(Hs, torch.ones(Hs.shape[0], Hs.shape[1], dtype=torch.bool))[0][0]
    with torch.no_grad():
        for i, c in enumerate(cells):
            h = torch.from_numpy(np.asarray(S.get(c), dtype=np.float32))[None]
            Oc["tok"][i] = o_tok(h).double().numpy()
            Oc["rand"][i] = o_rand(h).double().numpy()
            if i < tr.size:
                Lsum["tok"] += clamp_map_tokens(o_tok, h).double().numpy()
                Lsum["rand"] += clamp_map_tokens(o_rand, h).double().numpy()
            zs = h[0].mean(0).double().numpy()
            cos[i] = zs @ zoff[i] / (np.linalg.norm(zs) * np.linalg.norm(zoff[i]))
            if (i + 1) % 1000 == 0:
                say(f"  calibration token observers: {i + 1}/{cells.size} cells ({time.time() - t0:.0f}s)")
    info = {"n_train": int(tr.size), "n_val": int(va.size), "states_genemean_cos_vs_z_rna_min": float(cos.min()),
            "states_shard_addendum_sha256": sorted(S.shard_addendum), "observers": oinfo, "sec": round(time.time() - t0, 1)}
    pf = a.e3_dir / "r2" / "pred_subset.npz"
    if pf.exists():
        with np.load(pf) as d:
            pmap = {int(c): j for j, c in enumerate(d["train_cells"])}
            if all(int(c) in pmap for c in tr):
                P = d["train_pred"][[pmap[int(c)] for c in tr]][:, idx]
                info["r2_vs_e3_pred_subset_max_abs"] = float(np.max(np.abs(P - Oc["tok"][:tr.size])))
    return {"O": Oc, "Lsum": Lsum}, rp, info


def head_outputs(a, meta, cells: np.ndarray, ntr: int, zview: "RowView") -> tuple[np.ndarray, np.ndarray]:
    """O_head on the official z rows of the calibration cells (head normalisers from training rows only) and the sum
    of its Jacobians (= its clamp maps) over the first ntr (training) cells."""
    import torch
    from torch.func import jacrev, vmap
    yhat, _, _ = mld.load_head(a.head_ckpt, meta, torch.device("cpu"))
    Z = torch.from_numpy(zview.rows(cells).astype(np.float32))
    with torch.no_grad():
        Y = yhat(Z).double().numpy()
        J = torch.zeros(len(PROTEINS), Z.shape[1], dtype=torch.float64)
        for s in range(0, ntr, 256):
            zz = Z[s:min(s + 256, ntr)]
            J += vmap(jacrev(lambda z: yhat(z[None])[0]))(zz).double().sum(0)
    return Y, J.numpy()


def calibrate_core(meta: dict, Oc: dict, tr: np.ndarray, va: np.ndarray, zview: RowView) -> tuple[dict, dict]:
    """The frozen training maps and the output standardisation: per observer, the training SD of each output; ridge
    of the standardised training outputs on raw z_12 (train rows), lambda by val MSE; the training-average clamp map.
    Reads only the rows tr, va of z."""
    split = np.asarray(meta["split"]).astype(str)
    if not (np.all(split[tr] == "train") and np.all(split[va] == "val")):
        raise SystemExit("calibration cells must be split=train (fit) and split=val (lambda)")
    ntr = tr.size
    no, p_ = len(OBSERVERS), len(PROTEINS)
    H, b, lam = np.zeros((no, p_, 512)), np.zeros((no, p_)), np.zeros(no)
    out_mean, out_sd = np.zeros((no, p_)), np.zeros((no, p_))
    Htc = np.zeros((no, p_, 512))
    summ = {"lambdas": list(LAMBDAS), "observers": {}}
    Xtr, Xva = zview.rows(tr), zview.rows(va)
    for oi, o in enumerate(OBSERVERS):
        Y = Oc["O"][o]
        mu, sd = Y[:ntr].mean(0), Y[:ntr].std(0)
        sd = np.where(sd > 1e-12, sd, 1.0)
        out_mean[oi], out_sd[oi] = mu, sd
        Htc[oi] = Oc["Lsum"][o] / ntr / sd[:, None]
        f = e5m.ridge_select(Xtr, (Y[:ntr] - mu) / sd, Xva, (Y[ntr:] - mu) / sd, LAMBDAS)
        H[oi], b[oi], lam[oi] = f["H"], f["intercept"], f["lambda"]
        summ["observers"][o] = {"out_sd": sd.tolist(), "train_clamp_map_fro": float(np.linalg.norm(Htc[oi])),
                                "ridge": {"lambda": f["lambda"], "val_mse": {str(k): round(v_, 8) for k, v_ in f["val_mse"].items()},
                                          "train_r2": [round(x, 6) for x in f["train_r2"]],
                                          "val_r2": [round(x, 6) for x in f["val_r2"]]}}
    arrays = {"H_ridge": H, "ridge_intercept": b, "ridge_lambda": lam, "out_mean": out_mean, "out_sd": out_sd,
              "H_train_clamp": Htc}
    return arrays, summ


def core_sha(arrays: dict, summ: dict) -> str:
    h = hashlib.sha256()
    for k in sorted(arrays):
        h.update(k.encode())
        h.update(np.ascontiguousarray(arrays[k], dtype=np.float64).tobytes())
    h.update(json.dumps(summ, sort_keys=True).encode())
    return h.hexdigest()


def stage_register(a) -> None:
    """Calibrate H on train/val, run the leakage check, write E5M_H.json, E5M.json and their HASHES.txt lines."""
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    reg = load_reg(a)
    meta = lp.load_meta(a.processed, need_rna=False)
    split = np.asarray(meta["split"]).astype(str)
    seed = [int(reg["seeds"]["e5_random_directions"]), SEED_TAG]
    e3_add = json.loads((addenda_dir(a) / "E3.json").read_text())
    rc = e3_add["computed_train_val"]["r2_cells"]
    with np.load(a.e3_dir / "r2" / "pred_subset.npz") as d:
        tr, va = d["train_cells"].astype(np.int64), d["val_cells"].astype(np.int64)
    if vk.index_hash(tr) != rc["train_sha256"] or vk.index_hash(va) != rc["val_sha256"]:
        raise SystemExit("pred_subset.npz cells are not E3's registered R2 train / val cells")
    if a.calib_limit:
        tr, va = tr[: a.calib_limit], va[: max(10, a.calib_limit // 4)]
    dj = json.loads((a.e5_dir / "e5_directions.json").read_text())
    if dj.get("subsampled"):
        raise SystemExit("E5's e5_directions.npz was built from a smoke subsample")
    with np.load(a.e5_dir / "e5_directions.npz") as d:
        u_nkt = d["u_nkt"].astype(np.float64)
    say(f"register: calibrating on {tr.size} training / {va.size} val cells (E3 R2 subset)")
    emb = lp.Embedding(a.embed_dir, None, len(split))
    cells = np.concatenate([tr, va])

    def views(poison_rows=None, pseed=0):
        return RowView(emb.z, poison_rows, pseed)

    def outputs(m_, zv_):
        O_ = {"O": dict(Otok["O"]), "Lsum": dict(Otok["Lsum"])}
        O_["O"]["head"], O_["Lsum"]["head"] = head_outputs(a, m_, cells, tr.size, zv_)
        return O_

    zv = views()
    Otok, rp, oinfo = token_observer_outputs(a, meta, tr, va, seed, zv)  # stored states: E3's R2 train / val cells only
    Oc = outputs(meta, zv)
    arrays, summ = calibrate_core(meta, Oc, tr, va, zv)
    real = core_sha(arrays, summ)
    say(f"register: calibration core {real[:16]}")
    # leakage check: site4 rows of every array the calibration can read are poisoned (protein, size factor, cell type,
    # z, layer means); the core must not change. Positive control: val rows of z / layer means poisoned -> changes.
    test_rows = np.where(split == "test")[0]
    pm = dict(meta)
    rng = np.random.default_rng(101)
    pm["adt"] = meta["adt"].copy()
    pm["adt"][test_rows] = rng.uniform(0, 8, (test_rows.size, meta["adt"].shape[1])).astype(np.float32)
    pm["adt_size_factor"] = np.asarray(meta["adt_size_factor"]).copy()
    pm["adt_size_factor"][test_rows] = rng.uniform(0.1, 10, test_rows.size)
    pm["cell_types"] = np.asarray(meta["cell_types"]).copy()
    pm["cell_types"][test_rows] = rng.choice(np.unique(meta["cell_types"]), test_rows.size)
    zp = views(test_rows, 101)
    arr_p, summ_p = calibrate_core(pm, outputs(pm, zp), tr, va, zp)
    poisoned = core_sha(arr_p, summ_p)
    val_rows = np.where(split == "val")[0]
    zq = views(val_rows, 303)
    arr_q, summ_q = calibrate_core(meta, outputs(meta, zq), tr, va, zq)
    val_poisoned = core_sha(arr_q, summ_q)
    leak = {"what": "site4 rows of adt -> U(0, 8), adt_size_factor -> U(0.1, 10), cell_types -> random type, z_rna and "
                    "z_rna_layer_means -> N(0, 1) (seed 101); positive control: val rows of z_rna and z_rna_layer_means "
                    "-> N(0, 1) (seed 303)",
            "n_site4_rows": int(test_rows.size), "real_core_sha256": real, "site4_poisoned_core_sha256": poisoned,
            "site4_poisoned_identical": poisoned == real, "val_poisoned_core_sha256": val_poisoned,
            "val_poisoned_differs": val_poisoned != real,
            "static": "the calibration reads E3's stored states only for E3's R2 train / val cells (asserted), never "
                      "outputs/v3/E3/r2/pred (site4 predictions) or any E1-E5 site4 result"}
    leak["passed"] = bool(leak["site4_poisoned_identical"] and leak["val_poisoned_differs"])
    say(f"register: leakage check site4-poisoned identical {leak['site4_poisoned_identical']}, val-poisoned differs "
        f"{leak['val_poisoned_differs']}")
    if not leak["passed"]:
        raise SystemExit("leakage check failed; nothing written")
    d = addenda_dir(a) if not a.smoke else a.out_dir
    d.mkdir(parents=True, exist_ok=True)
    hf = d / H_FILE
    hobj = {"note": "E5-M frozen numbers (train / val only): output standardisation, depth-12 training maps (secondary), "
                    "O_rand parameters, source directions; float64 values, exact JSON round trip",
            "observers": list(OBSERVERS), "proteins": list(PROTEINS), "calib_train_sha256": vk.index_hash(tr),
            "calib_val_sha256": vk.index_hash(va), "n_train": int(tr.size), "n_val": int(va.size),
            **{k: np.asarray(v, np.float64).tolist() for k, v in arrays.items()},
            **{k: np.asarray(rp[k], np.float64).tolist() for k in ("v_rand", "q_rand", "W_rand")},
            "u_nkt": np.asarray(u_nkt, np.float64).tolist()}
    tmp = hf.with_name(hf.name + ".tmp")
    tmp.write_text(json.dumps(hobj, sort_keys=True) + "\n")
    os.replace(tmp, hf)
    h_sha = sha256_file(hf)
    amend = {}
    for am in sorted(a.registration_dir.glob("amendment_A*.json")):
        if not am.name.endswith("_core.json"):
            amend[am.name] = sha256_file(am)
    anm, anm_info = load_anm(a)
    add = {
        "experiment": "E5M",
        "addendum_version": 1,
        "title": "E5-M: observer-conditioned sufficiency of TEDDY's gene-mean state (Mode A)",
        "addendum_to": "registration/registration_v3.json as amended by A1, A2 and A3: a new pre-specified Mode A experiment "
                       "registered before any site4 forward of its script. It reuses E5's registered cell set and JVP "
                       "machinery and E3's frozen R2; nothing registered for E1-E5 changes.",
        "registration_sha256": sha256_file(a.registration_dir / "registration_v3.json"),
        "amendments_sha256": amend,
        "script": f"scripts/mode_a_observer_sufficiency.py ({SCRIPT_VERSION}); statistics in bridge_anm/lib/v3_e5m.py",
        "fixed_before_site4": "every choice and number here was fixed on training and validation cells only and "
                              "committed before any site4 forward pass of this script",
        "disclosure": "the orchestrating session has already seen E1 (Mode B rerun) site4 results; nothing here depends "
                      "on E1. This builder read E5's e5_cells.json counts (cells per role / donor / class, token-count "
                      "medians; no outcome), already printed by E5's run, and no other site4 value. E3's R2 site4 "
                      "predictions (being written) were not read.",
        "declared": DECLARED,
        "calibration": {"H_file": f"registration/addenda/{H_FILE}", "H_file_sha256": h_sha, "core_sha256": real,
                        "summary": summ, "observers_info": oinfo, "seed_v_q_W": seed,
                        "rand_query_norm": float(np.linalg.norm(rp["q_rand"])),
                        "rand_row_norms": np.linalg.norm(rp["W_rand"], axis=1).tolist(),
                        "u_nkt_source": str(a.e5_dir / "e5_directions.npz"),
                        "u_nkt_sha256": sha256_file(a.e5_dir / "e5_directions.npz"),
                        "e3_addendum_sha256": sha256_file(addenda_dir(a) / "E3.json"),
                        "compute": a.compute_note},
        "inputs": {"e5_function_source_sha256": code_hashes(), "lib_v3_e5m_sha256": sha256_file(ROOT / "bridge_anm/lib/v3_e5m.py"),
                   "anm_analysis": anm_info},
        "leakage_check": leak,
        "smoke": bool(a.smoke),
    }
    f = d / ADDENDUM
    f.write_bytes((json.dumps(add, indent=1, sort_keys=True, default=_json_default) + "\n").encode())
    sha = sha256_file(f)
    hh = d / "HASHES.txt"
    lines = [ln for ln in (hh.read_text().splitlines() if hh.exists() else [])
             if not (ln.rstrip().endswith(f"  {ADDENDUM}") or ln.rstrip().endswith(f"  {H_FILE}"))]
    hh.write_text("\n".join(lines + [f"{sha}  {ADDENDUM}", f"{h_sha}  {H_FILE}"]) + "\n")
    say(f"register: wrote {f} (sha256 {sha}) and {hf} (sha256 {h_sha}) and their lines in {hh}")


def stage_commit_addendum(a) -> None:
    """Commit exactly E5M.json, E5M_H.json and their HASHES.txt lines (temporary index; other builders' changes are
    left untouched). Retries if HEAD moves."""
    d = addenda_dir(a)
    add = json.loads((d / ADDENDUM).read_text())
    if add.get("smoke") or not add["leakage_check"]["passed"]:
        raise SystemExit("refusing to commit a smoke addendum or one without a passed leakage check")
    if add["declared"] != json.loads(json.dumps(DECLARED)):
        raise SystemExit("the addendum's declared block differs from this script's DECLARED")
    if sha256_file(d / H_FILE) != add["calibration"]["H_file_sha256"]:
        raise SystemExit("E5M_H.json is not the file the addendum pins")
    names = (ADDENDUM, H_FILE)
    rel = {n: str((d / n).resolve().relative_to(ROOT.resolve())) for n in names}
    hrel = str((d / "HASHES.txt").resolve().relative_to(ROOT.resolve()))
    for attempt in range(10):
        head = git("rev-parse", "HEAD").strip()
        branch = git("symbolic-ref", "--short", "HEAD").strip()
        try:
            old = git("show", f"{head}:{hrel}")
        except RuntimeError:
            old = ""
        keep = [ln for ln in old.splitlines() if not any(ln.rstrip().endswith(f"  {n}") for n in names)]
        new_hashes = ("\n".join(keep + [f"{sha256_file(d / n)}  {n}" for n in names]) + "\n").encode()
        files = {rel[n]: (d / n).read_bytes() for n in names}
        files[hrel] = new_hashes
        msg = (f"E5M addendum: observer-conditioned sufficiency of the gene-mean state, fixed on train/val before site4 "
               f"(sha256 {sha256_file(d / ADDENDUM)[:12]}, H {add['calibration']['H_file_sha256'][:12]}), leakage check "
               "passed\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n")
        with tempfile.TemporaryDirectory() as td:
            env = {**os.environ, "GIT_INDEX_FILE": str(Path(td) / "index")}
            git("read-tree", head, env=env)
            for path, data in files.items():
                blob = git("hash-object", "-w", "--stdin", input_bytes=data).strip()
                git("update-index", "--add", "--cacheinfo", f"100644,{blob},{path}", env=env)
            tree = git("write-tree", env=env).strip()
        commit = git("commit-tree", tree, "-p", head, "-m", msg).strip()
        try:
            git("update-ref", f"refs/heads/{branch}", commit, head)
        except RuntimeError as e:
            say(f"commit-addendum: HEAD moved ({e}); retry {attempt + 1}")
            time.sleep(5)
            continue
        for t in range(60):  # sync the main index for these paths (waits on another process's index.lock)
            try:
                git("reset", "-q", "--", *files.keys())
                break
            except RuntimeError as e:
                if "index.lock" not in str(e):
                    raise
                time.sleep(5)
        say(f"committed the E5M addendum as {commit[:10]} on {branch}: {', '.join(files)}")
        return
    raise SystemExit("commit-addendum: HEAD kept moving; nothing committed")


def check_registration(a, site4: bool) -> dict:
    d = addenda_dir(a)
    f, hf, hh = d / ADDENDUM, d / H_FILE, d / "HASHES.txt"
    info = {"addendum_file": str(f), "addendum_sha256": sha256_file(f) if f.exists() else None}
    hl = hh.read_text().splitlines() if hh.exists() else []
    add = json.loads(f.read_text()) if f.exists() else {}
    info.update({
        "addendum_hash_recorded": bool(info["addendum_sha256"]) and f"{info['addendum_sha256']}  {ADDENDUM}" in hl,
        "addendum_committed": git_committed(f) if f.exists() else False,
        "declared_equals_script": bool(add) and add.get("declared") == json.loads(json.dumps(DECLARED)),
        "H_file_sha256": sha256_file(hf) if hf.exists() else None,
        "leakage_check_passed": bool(add) and add.get("leakage_check", {}).get("passed") is True,
        "smoke_addendum": bool(add.get("smoke")),
    })
    info["H_file_matches_addendum"] = bool(add) and info["H_file_sha256"] == add["calibration"]["H_file_sha256"]
    info["H_file_hash_recorded"] = bool(info["H_file_sha256"]) and f"{info['H_file_sha256']}  {H_FILE}" in hl
    info["H_file_committed"] = git_committed(hf) if hf.exists() else False
    rec = add.get("calibration", {}).get("observers_info", {}).get("observers", {}) if add else {}
    r2p, hp = getattr(a, "e3_dir", None), getattr(a, "head_ckpt", None)
    info["r2_matches_addendum"] = bool(rec) and r2p is not None and (r2p / "r2" / "r2.pt").exists() and \
        sha256_file(r2p / "r2" / "r2.pt") == rec.get("r2_sha256")
    info["head_matches_addendum"] = bool(rec) and hp is not None and Path(hp).exists() and \
        sha256_file(Path(hp)) == rec.get("head_sha256")
    cf, ch = a.registration_dir / "registration_v3.json", a.registration_dir / "registration_v3.json.sha256"
    info["registration_sha256"] = sha256_file(cf) if cf.exists() else None
    info["registration_ok"] = bool(info["registration_sha256"]) and ch.exists() and ch.read_text().split()[0] == \
        info["registration_sha256"] and git_committed(cf) is True
    info["amendments"] = {}
    for am in sorted(a.registration_dir.glob("amendment_A*.json")):
        if am.name.endswith("_core.json"):
            continue
        h = sha256_file(am)
        hf_ = am.with_name(am.name + ".sha256")
        info["amendments"][am.name] = {"sha256": h, "committed": git_committed(am),
                                       "hash_file_matches": hf_.exists() and hf_.read_text().split()[0] == h}
    if site4:
        need = {"registration_v3.json committed with a matching hash file": info["registration_ok"],
                "every amendment (A1, A2, A3 if present) committed with a matching hash file":
                    all(v["committed"] is True and v["hash_file_matches"] for v in info["amendments"].values())
                    and {"amendment_A1.json", "amendment_A2.json"} <= set(info["amendments"]),
                "E5M addendum present, committed, sha256 in addenda/HASHES.txt":
                    bool(info["addendum_sha256"]) and info["addendum_committed"] is True and info["addendum_hash_recorded"],
                "addendum declared block equals this script's DECLARED": info["declared_equals_script"],
                "E5M_H.json matches the addendum, committed, sha256 in HASHES.txt":
                    info["H_file_matches_addendum"] and info["H_file_committed"] is True and info["H_file_hash_recorded"],
                "addendum leakage check passed, not a smoke addendum": info["leakage_check_passed"] and not info["smoke_addendum"],
                "R2 checkpoint and head checkpoint are the ones the addendum records":
                    info["r2_matches_addendum"] and info["head_matches_addendum"]}
        bad = [k for k, ok in need.items() if not ok]
        if bad:
            raise SystemExit("site4 run refused (register and commit the E5M addendum first): " + "; ".join(bad))
    return info


def load_H(a) -> dict:
    src = (a.out_dir / H_FILE) if (a.smoke and (a.out_dir / H_FILE).exists()) else addenda_dir(a) / H_FILE
    raw = json.loads(src.read_text())
    out = {k: (np.asarray(v, np.float64) if k in H_ARRAYS else v) for k, v in raw.items()}
    out["_path"], out["_sha256"] = str(src), sha256_file(src)
    return out


def addendum_sha(a) -> str:
    f = (a.out_dir / ADDENDUM) if (a.smoke and (a.out_dir / ADDENDUM).exists()) else addenda_dir(a) / ADDENDUM
    return sha256_file(f) if f.exists() else "none"


# ============================================================================ cells
def stage_select(a, meta, reg) -> dict:
    out = a.out_dir
    split = np.asarray(meta["split"]).astype(str)
    donors = mld.load_donors(a.processed)
    if a.cell_pool == "val":
        key = mld.primary_keys(meta, reg)  # the registered key (val rows used only)
        rng = np.random.default_rng(int(reg["seeds"]["e5_cells"]))
        sel, cls = [], []
        for c in ("NK", "T"):
            pool = np.where((split == "val") & (key == c))[0]
            pick = np.sort(rng.choice(pool, min(a.n_val_cells, pool.size), replace=False))
            sel.append(pick)
            cls += [c] * pick.size
        sel = np.concatenate(sel)
        cls = np.array(cls)
        order = rng.permutation(sel.size)
        role = np.empty(sel.size, dtype="<U6")
        role[order] = ["pair" if k % 2 == 0 else "random" for k in range(sel.size)]  # alternate along the order
        info = {"pool": "val donor 18303 primary-key NK/T cells (smoke; roles alternate along the order)"}
    else:
        f = a.e5_dir / "e5_cells.npz"
        src = "E5 out-dir"
        if not f.exists():
            key = mld.primary_keys(meta, reg)  # E5's registered selection needs the key (as E5's select stage)
            ns = SimpleNamespace(out_dir=out / "e5_rebuilt", cell_pool="site4", processed=a.processed)
            ns.out_dir.mkdir(parents=True, exist_ok=True)
            emb = lp.Embedding(a.embed_dir, None, len(split))
            mld.stage_select(ns, meta, emb, key, reg)
            f, src = ns.out_dir / "e5_cells.npz", "rebuilt with E5's stage_select"
        with np.load(f, allow_pickle=False) as d:
            e5c, e5r, e5o, e5k = d["cells"], d["role"].astype(str), d["order"], d["cls"].astype(str)
        if not np.all(np.isin(donors[e5c], list(reg["splits"]["test_primary"]["donors"]))):
            raise SystemExit("E5 cells outside test_primary (is this E5's site4 selection?)")
        keep = np.zeros(e5c.size, bool)
        for dn in reg["splits"]["test_primary"]["donors"]:
            for rname, n_ in (("pair", N_PAIR_PER_DONOR), ("random", N_RAND_PER_DONOR)):
                ro = [i for i in e5o if e5r[i] == rname and donors[e5c[i]] == dn]
                keep[ro[:n_]] = True
        o2 = np.array([i for i in e5o if keep[i]])
        remap = {int(i): k for k, i in enumerate(np.where(keep)[0])}
        sel, role, cls = e5c[keep], e5r[keep], e5k[keep]
        order = np.array([remap[int(i)] for i in o2])
        info = {"pool": f"site4 test_primary: per donor the first {N_PAIR_PER_DONOR} pair and {N_RAND_PER_DONOR} random E5 "
                        "cells of E5's order", "e5_cells_file": str(f), "e5_cells_source": src,
                "e5_cells_sha256": sha256_file(f)}
    info.update({"n_cells": int(sel.size),
                 "by_role_donor": {f"{r}/{dn}": int(np.sum((role == r) & (donors[sel] == dn)))
                                   for r in np.unique(role) for dn in np.unique(donors[sel])}})
    atomic_savez(out / "e5m_cells.npz", cells=sel, role=role, order=order, cls=cls, donor=donors[sel])
    write_json(out / "e5m_cells.json", info)
    say(f"select: {info['by_role_donor']}")
    return info


def load_cells(out: Path):
    with np.load(out / "e5m_cells.npz", allow_pickle=False) as d:
        return d["cells"], d["role"].astype(str), d["order"], d["cls"].astype(str), d["donor"].astype(str)


# ============================================================================ respond
def stage_respond(a, meta, reg, Hz: dict) -> int:
    import torch
    from torch.func import jvp, vmap

    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import load_pad_id, load_teddy, load_vocab
    torch.set_num_threads(_THREADS)
    out = a.out_dir
    cdir = out / "cells"
    cdir.mkdir(exist_ok=True)
    sel, role, order, cls, donor = load_cells(out)
    fd_all = {int(i) for i in order[:N_FD_ALL_COLUMNS]} if not a.no_extras else set()
    fd_some = {int(i) for i in order[:N_FD_CELLS]}
    layer_set = {int(i) for i in order[:N_LAYER_CELLS]} if not a.no_extras else set()
    if a.limit:
        order = order[: a.limit]
    device = resolve_device(a.device)
    emb = lp.Embedding(a.embed_dir, None, len(meta["split"]))
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    genes, _, _ = mld.e2_genes(reg, a.gene_reference, meta["rna_names"], vocab)
    gtok = {g[0]: g[2] for g in genes}
    decl = [(s_, gtok[s_]) for s_ in NKT_GENES]
    sha = addendum_sha(a)
    todo = []
    for i in order:
        f = cdir / f"cell_{int(sel[i])}.npz"
        if f.exists():
            with np.load(f, allow_pickle=False) as d:
                if json.loads(str(d["meta"])).get("addendum_sha256") != sha:
                    raise SystemExit(f"{f} was written under another addendum; use a fresh --out-dir")
        else:
            todo.append(int(i))
    say(f"respond: {len(order)} cells in the order, {len(order) - len(todo)} done, {len(todo)} to go; device {device}, "
        f"float32, jvp chunk {a.jvp_chunk}, fd batch {a.fd_batch}")
    if not todo:
        return 0
    toks, ntok = mld.tokenise(meta, sel[todo], a, emb, pad_id)
    if np.any(ntok != emb.ntokens[sel[todo]]):
        raise SystemExit("token counts differ from the official run")
    model = load_teddy(a.ckpt, device)
    for p_ in model.parameters():
        p_.requires_grad_(False)
    layers = list(model.encoder.layers)
    D = model.d_model
    E, Pm = model.embeddings.weight, model.position_embeddings.weight
    a.reg_seed = None
    r2, idx, rr, yhat, _, _ = load_observers(a, meta, device, Hz)
    idx_t = torch.tensor(idx, device=device)
    f32 = lambda x: torch.as_tensor(np.asarray(x, np.float32), device=device)  # noqa: E731
    u_nkt, v_rand = f32(Hz["u_nkt"]), f32(Hz["v_rand"])
    t_start, n_done = time.time(), 0

    def obs_tokens(Hs):  # Hs [B, L, d] -> (tok [B, 4], rand [B, 4])
        m = torch.ones(Hs.shape[0], Hs.shape[1], dtype=torch.bool, device=Hs.device)
        return r2(Hs, m)[0].index_select(-1, idx_t), rr(Hs, m)[0]

    def o_tok(Hs):
        return obs_tokens(Hs)[0][0]

    def o_rand(Hs):
        return obs_tokens(Hs)[1][0]

    def o_head(z):
        return yhat(z[None])[0]

    for jj, i in enumerate(todo):
        if a.max_minutes and time.time() - t_start > 60 * a.max_minutes:
            say(f"respond: --max-minutes {a.max_minutes} reached after {n_done} cells; rerun the same command to resume")
            return EXIT_INCOMPLETE
        tc = time.time()
        g_id = int(sel[i])
        ids_np = toks[jj]
        L = ids_np.size
        ids = torch.from_numpy(ids_np).to(device)
        posmap = {int(t): q for q, t in enumerate(ids_np)}
        present = [(gi, posmap[t]) for gi, (s_, t) in enumerate(decl) if t in posmap]
        cols = []  # (gene index, direction index, position, unit vector)
        for gi, q in present:
            tid = decl[gi][1]
            for di, v in enumerate((E[tid] / E[tid].norm(), u_nkt, v_rand)):
                cols.append((gi, di, q, v))
        m = len(cols)
        meta_cell = {"cell": g_id, "role": str(role[i]), "class": str(cls[i]), "donor": str(donor[i]), "ntokens": int(L),
                     "addendum_sha256": sha, "genes_present": [decl[gi][0] for gi, _ in present], "n_columns": m,
                     "script_version": SCRIPT_VERSION}
        if m == 0:
            atomic_savez(cdir / f"cell_{g_id}.npz", meta=np.array(json.dumps({**meta_cell, "sec": round(time.time() - tc, 2)})),
                         n_columns=np.array(0))
            n_done += 1
            continue
        with torch.no_grad():
            h0 = (E[ids] + Pm[:L])[None]
            hs = [h0]
            h = h0
            for layer in layers:
                h = mld.layer_fn(layer, h)
                hs.append(h)
            zb = torch.stack([x.mean(1)[0] for x in hs])  # [13, d]
            _, ref_lm = model.hidden_states(ids[None], torch.ones(1, L, dtype=torch.long, device=device), return_layer_means=True)
            chk_module = float((ref_lm[:, 0] - zb[1:]).abs().max())
            H12 = hs[12]
            z12 = zb[12]
            O0_tok, O0_rand = obs_tokens(H12)
            O0_head = yhat(z12[None])
            # ---------------- primary H: clamp maps of z_12 (from clamps alone, before any source response)
            Lam12 = torch.stack([clamp_map_tokens(o_tok, H12), torch.func.jacrev(o_head)(z12),
                                 clamp_map_tokens(o_rand, H12)])  # [3, 4, d] in OBSERVERS order
        # ---------------- secondary: clamp maps of z_l for every depth by reverse mode (first cells of the order)
        Lam_layers = None
        if int(i) in layer_set:
            with torch.enable_grad():
                x0 = h0.detach().clone().requires_grad_(True)
                states = [x0]
                hh = x0
                for layer in layers:
                    hh = mld.layer_fn(layer, hh)
                    states.append(hh)
                outs = torch.cat([o_tok(hh), o_head(hh.mean(1)[0]), o_rand(hh)])  # 12 = 3 observers x 4
                LL = torch.zeros(len(OBSERVERS), N_DEPTHS, len(PROTEINS), D, device=device)
                for k in range(outs.numel()):
                    gs = torch.autograd.grad(outs[k], states, retain_graph=k < outs.numel() - 1)
                    for l_ in range(N_DEPTHS):
                        LL[k // len(PROTEINS), l_, k % len(PROTEINS)] = gs[l_][0].sum(0)
                Lam_layers = LL.detach()
                del states, hh, outs, gs, x0
        T0 = torch.zeros(m, 1, L, D, device=device)
        pos_k = torch.tensor([c[2] for c in cols], device=device)
        T0[torch.arange(m, device=device), 0, pos_k] = torch.stack([c[3] for c in cols])
        J_pool = torch.zeros(m, N_DEPTHS, D, device=device)
        K_tok = torch.zeros(m, len(PROTEINS), device=device)
        K_rand, K_tok_dev, K_rand_dev = torch.zeros_like(K_tok), torch.zeros_like(K_tok), torch.zeros_like(K_tok)
        with torch.no_grad():
            for c0 in range(0, m, a.jvp_chunk):
                c1 = min(c0 + a.jvp_chunk, m)
                dh = T0[c0:c1]
                J_pool[c0:c1, 0] = dh[:, 0].mean(1)
                for li, layer in enumerate(layers):
                    f = lambda x, _layer=layer: mld.layer_fn(_layer, x)  # noqa: E731
                    dh = vmap(lambda t, _h=hs[li], _f=f: jvp(_f, (_h,), (t,))[1])(dh)
                    J_pool[c0:c1, li + 1] = dh[:, 0].mean(1)
                dev = dh - dh.mean(dim=2, keepdim=True)  # the pooling-discarded part of dH_12
                K_tok[c0:c1] = vmap(lambda t: jvp(o_tok, (H12,), (t,))[1])(dh)
                K_rand[c0:c1] = vmap(lambda t: jvp(o_rand, (H12,), (t,))[1])(dh)
                K_tok_dev[c0:c1] = vmap(lambda t: jvp(o_tok, (H12,), (t,))[1])(dev)
                K_rand_dev[c0:c1] = vmap(lambda t: jvp(o_rand, (H12,), (t,))[1])(dev)
                del dh, dev
            K_head = vmap(lambda t: jvp(o_head, (z12,), (t,))[1])(J_pool[:, 12])  # [m, 4]
            # ---------------- finite differences through the module's own layers (z_12 and the three observers)
            fd_cols = (list(range(m)) if int(i) in fd_all else
                       [k for k in range(m) if cols[k][0] == cols[0][0]] if int(i) in fd_some else [])
            jobs = [(k, ei, sg) for k in fd_cols for ei in range(len(EPS_FD)) for sg in (1, -1)]
            nb = a.fd_batch
            res = {}
            for b0 in range(0, len(jobs), nb):
                bj = jobs[b0:b0 + nb]
                kk = torch.tensor([j[0] for j in bj], device=device)
                cc = torch.tensor([j[2] * EPS_FD[j[1]] for j in bj], device=device, dtype=torch.float32)
                x = h0.expand(len(bj), L, D) + cc[:, None, None] * T0[kk, 0]
                for layer in layers:
                    x = layer(x)
                zz = x.mean(1)
                ot, orr = obs_tokens(x)
                oh = yhat(zz)
                for r_, j in enumerate(bj):
                    res[j] = (zz[r_], ot[r_], oh[r_], orr[r_])
                del x
            nf = len(fd_cols)
            fdz = np.zeros((nf, len(EPS_FD), D))
            fdo = np.zeros((nf, len(EPS_FD), len(OBSERVERS), len(PROTEINS)))
            for fi, k in enumerate(fd_cols):
                for ei, e in enumerate(EPS_FD):
                    p, q = res[(k, ei, 1)], res[(k, ei, -1)]
                    fdz[fi, ei] = ((p[0] - q[0]) / (2 * e)).cpu().double().numpy()
                    for oi, (pi, qi) in enumerate(((p[1], q[1]), (p[2], q[2]), (p[3], q[3]))):
                        fdo[fi, ei, oi] = ((pi - qi) / (2 * e)).cpu().double().numpy()
        meta_cell.update({"check_explicit_vs_module_layer_means_max_abs": chk_module, "fd_columns": fd_cols,
                          "layer_clamps": Lam_layers is not None, "sec": round(time.time() - tc, 2)})
        extra = {"Lam_layers": Lam_layers.cpu().numpy()} if Lam_layers is not None else {}
        atomic_savez(cdir / f"cell_{g_id}.npz", meta=np.array(json.dumps(meta_cell)), n_columns=np.array(m),
                     col_gene=np.array([c[0] for c in cols], np.int16), col_dir=np.array([c[1] for c in cols], np.int8),
                     col_pos=np.array([c[2] for c in cols], np.int32), J_pool=J_pool.cpu().numpy(),
                     K_tok=K_tok.cpu().double().numpy(), K_head=K_head.cpu().double().numpy(),
                     K_rand=K_rand.cpu().double().numpy(), K_tok_dev=K_tok_dev.cpu().double().numpy(),
                     K_rand_dev=K_rand_dev.cpu().double().numpy(), Lam12=Lam12.cpu().double().numpy(),
                     O0=torch.cat([O0_tok[0], O0_head[0], O0_rand[0]]).cpu().double().numpy(),
                     z12=z12.cpu().double().numpy(), fd_cols=np.array(fd_cols, np.int32), fd_z=fdz, fd_obs=fdo, **extra)
        n_done += 1
        rate = (time.time() - t_start) / n_done
        say(f"  cell {jj + 1}/{len(todo)} ({role[i]}, {cls[i]}, {L} tokens, {m} columns, fd {len(fd_cols)}, layers "
            f"{Lam_layers is not None}) {time.time() - tc:.1f}s; module check {chk_module:.1e}; mean {rate:.1f}s/cell, "
            f"ETA {rate * (len(todo) - n_done) / 60:.1f} min")
        write_json(out / "e5m_progress.json", {"done_this_call": n_done, "to_do_this_call": len(todo),
                                                "sec_per_cell": round(rate, 2),
                                                "eta_min": round(rate * (len(todo) - n_done) / 60, 1),
                                                "updated": time.strftime("%Y-%m-%d %H:%M:%S")})
        if device.type == "mps":
            torch.mps.empty_cache()
    say(f"respond done: {n_done} cells in {(time.time() - t_start) / 60:.1f} min")
    return 0


# ============================================================================ report
def cell_stats(d, Hz: dict, anm) -> dict:
    """Every per-cell statistic from one cell file (K_O and every H in standardised output units)."""
    m = int(d["n_columns"])
    sd = Hz["out_sd"]  # [3, 4]
    Jp = d["J_pool"].astype(np.float64)  # [m, 13, d]
    Kz = Jp[:, 12].T
    K = {o: d[f"K_{o}"].T / sd[oi][:, None] for oi, o in enumerate(OBSERVERS)}  # [4, m]
    Hc = {o: d["Lam12"][oi] / sd[oi][:, None] for oi, o in enumerate(OBSERVERS)}  # primary clamp maps [4, d]
    s = {"m": m}
    for oi, o in enumerate(OBSERVERS):
        s[f"r_{o}"] = e5m.factorization_residual(K[o], Hc[o], Kz)
        s[f"r_{o}_ridge"] = e5m.factorization_residual(K[o], Hz["H_ridge"][oi], Kz)
        s[f"r_{o}_train_clamp"] = e5m.factorization_residual(K[o], Hz["H_train_clamp"][oi], Kz)
        if "Lam_layers" in d.files:
            s[f"r_{o}_clamp_layers"] = [e5m.factorization_residual(K[o], d["Lam_layers"][oi, li_] / sd[oi][:, None], Jp[:, li_].T)
                                        for li_ in range(N_DEPTHS)]
            s[f"lam12_reverse_vs_primary_{o}"] = float(np.max(np.abs(d["Lam_layers"][oi, 12] - d["Lam12"][oi]))
                                                       / (np.max(np.abs(d["Lam12"][oi])) + 1e-30))
        for di, dn in enumerate(DIRECTIONS):
            cm = d["col_dir"] == di
            s[f"r_{o}_{dn}"] = e5m.factorization_residual(K[o][:, cm], Hc[o], Kz[:, cm]) if cm.any() else float("nan")
        q = e5m.quotient_diagnostic(K[o], QUOTIENT_RTOL)
        s[f"quot_rank_{o}"], s[f"quot_pr_{o}"] = q["rank"], q["participation_ratio"]
        for tol in NULL_TOL_REPORTED:
            nn = e5m.near_null_share(Kz, K[o], tol)
            s[f"null_share_{o}_{tol}"] = nn["share"]
            s[f"null_dim_{tol}"], s["Kz_condition"] = nn["near_null_dim"], nn["condition"]
        if anm is not None:
            for atol, rtol in ANM_BUDGETS:
                st = e5m.anm_factorization_status(anm.analyze_factorization, Kz, K[o], Hc[o], atol=atol, rtol=rtol)
                s[f"anm_{o}_{rtol}"] = {k: st[k] for k in ("status", "source_kernel_status", "source_kernel_reason",
                                                           "factorization_status", "Kz_rank", "Kz_nullity")}
    for o in ("tok", "rand"):  # the clamp residual computed directly from the pooling-discarded part of dH_12
        Kd = d[f"K_{o}_dev"].T / sd[OBSERVERS.index(o)][:, None]
        s[f"r_{o}_direct"] = float(np.linalg.norm(Kd) / np.linalg.norm(K[o])) if np.linalg.norm(K[o]) > 0 else float("nan")
    s["r_tok_same_cell_fit"] = e5m.same_cell_fit_residual(K["tok"], Kz)
    s["Kz_rank_full"] = int(np.linalg.matrix_rank(Kz))
    # FD gate: targets z_12 (raw), O_tok, O_head, O_rand (standardised); each eps and the Richardson derivative
    fd = []
    for fi, k in enumerate(d["fd_cols"]):
        Dz = [d["fd_z"][fi, ei] for ei in range(len(EPS_FD))]
        Do = [[d["fd_obs"][fi, ei, oi] / sd[oi] for oi in range(len(OBSERVERS))] for ei in range(len(EPS_FD))]
        for ei in range(len(EPS_FD) + 1):  # the last = Richardson (4 D(0.01) - D(0.02)) / 3
            if ei < len(EPS_FD):
                dz, do = Dz[ei], Do[ei]
            else:
                dz = e5m.richardson(Dz[0], Dz[1])
                do = [e5m.richardson(Do[0][oi], Do[1][oi]) for oi in range(len(OBSERVERS))]
            row = [e5m.fd_check(Kz[:, k], dz, FD_REL_TOL, FD_COS_TOL)]
            row += [e5m.fd_check(K[o][:, k], do[oi], FD_REL_TOL, FD_COS_TOL) for oi, o in enumerate(OBSERVERS)]
            fd.append((int(k), ei, row))
    s["fd"] = fd
    return s


def stage_report(a, meta, reg, reg_info, Hz) -> None:
    out = a.out_dir
    sel, role, _, cls, donor = load_cells(out)
    anm, anm_info = load_anm(a)
    sha = addendum_sha(a)
    rows, metas = [], []
    want = {int(c) for c in sel}
    for f in sorted((out / "cells").glob("cell_*.npz")):
        with np.load(f, allow_pickle=False) as d:
            mc = json.loads(str(d["meta"]))
            if mc["cell"] not in want:
                continue
            if mc.get("addendum_sha256") != sha:
                raise SystemExit(f"{f} was written under another addendum")
            metas.append(mc)
            if int(d["n_columns"]) == 0:
                continue
            s = cell_stats(d, Hz, anm)
        s.update({"cell": mc["cell"], "role": mc["role"], "cls": mc["class"], "donor": mc["donor"]})
        rows.append(s)
    R = {"experiment": "E5M", "script_version": SCRIPT_VERSION, "smoke_test": bool(a.smoke), "smoke_note": a.smoke,
         "cell_pool": a.cell_pool, "created": time.strftime("%Y-%m-%d %H:%M:%S"), "addendum_sha256": sha,
         "registration_sha256": reg_info.get("registration_sha256"), "registration_check": reg_info,
         "H_file": Hz["_path"], "H_file_sha256": Hz["_sha256"], "anm_analysis": anm_info, "n_boot": a.n_boot,
         "n_cells_selected": int(sel.size), "n_cells_done": len(metas), "complete": len(metas) == sel.size,
         "n_cells_without_columns": int(sum(1 for m in metas if m["n_columns"] == 0)),
         "n_cells_analysed": len(rows), "code": {"e5_function_source_sha256": code_hashes(),
                                                  "lib_v3_e5m_sha256": sha256_file(ROOT / "bridge_anm/lib/v3_e5m.py")}}
    if not rows:
        write_json(out / "E5M_results.json", R)
        say("report: no analysed cells yet")
        return
    rl = np.array([r["role"] for r in rows])
    dl = np.array([r["donor"] for r in rows])
    ix = np.arange(len(rows))
    col = lambda k: np.array([np.nan if r.get(k) is None else r[k] for r in rows], np.float64)  # noqa: E731
    roles = sorted(set(rl.tolist()))
    R["columns_per_cell"] = {"median": float(np.median(col("m"))), "min": int(col("m").min()), "max": int(col("m").max())}
    R["by_role_donor"] = {f"{r}/{dn}": int(np.sum((rl == r) & (dl == dn))) for r in roles for dn in np.unique(dl)}
    # ---------------- FD gate (Richardson derivative decides; each eps reported)
    names = [str(e) for e in EPS_FD] + ["richardson"]
    gate = {"derivatives": names, "decides": "richardson", "targets": ["z_12"] + [f"O_{o}" for o in OBSERVERS],
            "rel_tol": FD_REL_TOL, "cos_tol": FD_COS_TOL, "required_pass_share": FD_PASS_SHARE, "by_derivative_target": {}}
    ok_all = True
    for ei, nm in enumerate(names):
        for ti, t in enumerate(gate["targets"]):
            v = [row[ti] for r in rows for (_, e_, row) in r["fd"] if e_ == ei]
            share = float(np.mean([x[2] for x in v])) if v else float("nan")
            gate["by_derivative_target"][f"{nm}/{t}"] = {
                "n": len(v), "pass_share": share, "median_rel_err": float(np.median([x[0] for x in v])) if v else None,
                "max_rel_err": float(np.max([x[0] for x in v])) if v else None, "min_cos": float(np.min([x[1] for x in v])) if v else None}
            if nm == "richardson":
                ok_all &= bool(v) and share >= FD_PASS_SHARE
    gate["validated"] = bool(ok_all)
    R["fd_gate"] = gate
    r_head_max = float(np.nanmax(col("r_head")))
    R["positive_control"] = {"r_head_max": r_head_max, "tol": POSITIVE_CONTROL_TOL, "passed": r_head_max <= POSITIVE_CONTROL_TOL}
    R["pipeline_checks"] = {
        "explicit_vs_module_layer_means_max_abs": max(m.get("check_explicit_vs_module_layer_means_max_abs", 0.0) for m in metas),
        "clamp_vs_direct_tok_max_abs_diff": float(np.nanmax(np.abs(col("r_tok") - col("r_tok_direct")))),
        "clamp_vs_direct_rand_max_abs_diff": float(np.nanmax(np.abs(col("r_rand") - col("r_rand_direct")))),
        "clamp_vs_frozen_train_clamp_max_abs_diff": {o: float(np.nanmax(np.abs(col(f"r_{o}") - col(f"r_{o}_train_clamp"))))
                                                     for o in ("tok", "rand")},
        "r_tok_same_cell_fit_max": float(np.nanmax(col("r_tok_same_cell_fit"))),
        "Kz_full_column_rank_share": float(np.mean([r["Kz_rank_full"] == r["m"] for r in rows])),
        "lam12_reverse_vs_primary_max_rel": {o: float(np.nanmax(col(f"lam12_reverse_vs_primary_{o}")))
                                             for o in OBSERVERS if any(f"lam12_reverse_vs_primary_{o}" in r for r in rows)},
        "sec_per_cell_median": float(np.median([m["sec"] for m in metas])),
    }
    # ---------------- registered contrasts (one set of two-stage draws for every statistic)
    draws = e5m.two_stage_draws(dl, rl, a.n_boot, int(reg["seeds"]["bootstrap"]))
    rv = {o: col(f"r_{o}") for o in OBSERVERS}
    md = {o: e5m.stratum_median_diff(rv[o], rl) for o in OBSERVERS}

    def pair_median(x):
        def f(ii):
            ii = ii[rl[ii] == "pair"]
            return e5m.finite_median(x[ii]) if ii.size else float("nan")
        return f
    stats = {"delta_tok": md["tok"], "delta_rand": md["rand"], "delta_head": md["head"],
             "dd_tok_minus_rand": lambda ii: md["tok"](ii) - md["rand"](ii),
             "level_tok_pair": pair_median(rv["tok"]), "level_tok_random": lambda ii: e5m.finite_median(rv["tok"][ii[rl[ii] == "random"]]),
             "power_rand_pair": pair_median(rv["rand"])}
    R["medians"] = {f"{o}/{r}": e5m.finite_median(rv[o][rl == r]) for o in OBSERVERS for r in roles}
    if {"pair", "random"} <= set(roles):
        C = e5m.contrast(stats, ix, dl, draws)
        V = e5m.h5m_verdict(C, MARGIN, gate["validated"], R["positive_control"]["passed"])
    else:
        C, V = {}, {"code": "INCONCLUSIVE", "verdict": "only one role present"}
    R["registered"] = {"contrasts": C, "margin": MARGIN, **V}
    if a.cell_pool != "site4":
        R["registered"]["verdict"] = f"(smoke on {a.cell_pool} cells; not a registered result) {R['registered']['verdict']}"
    # ---------------- secondary
    sec = {}
    if {"pair", "random"} <= set(roles):
        lg = {o: e5m.log_or_nan(rv[o]) for o in ("tok", "rand")}
        sec["S2_log"] = e5m.contrast({f"delta_log_{o}": e5m.stratum_median_diff(lg[o], rl) for o in lg}, ix, dl, draws)
        S1 = {}
        for kind, n_ in (("clamp_layers", N_DEPTHS),):
            for o in OBSERVERS:
                have = [k for k, r in enumerate(rows) if f"r_{o}_{kind}" in r]
                if not have:
                    continue
                Lr = np.array([rows[k][f"r_{o}_{kind}"] for k in have])
                rr_ = rl[have]
                S1[f"{o}/{kind}"] = {"n_cells": {r: int(np.sum(rr_ == r)) for r in roles},
                                     "depths": list(range(N_DEPTHS)),
                                     "median_pair": [e5m.finite_median(Lr[rr_ == "pair", li_]) for li_ in range(n_)],
                                     "median_random": [e5m.finite_median(Lr[rr_ == "random", li_]) for li_ in range(n_)]}
        sec["S1_layers"] = S1
        sec["S3_directions"] = {f"{o}_{dn}": {r: e5m.finite_median(col(f"r_{o}_{dn}")[rl == r]) for r in roles}
                                for o in OBSERVERS for dn in DIRECTIONS}
        S7 = {}
        for kind in ("ridge", "train_clamp"):
            xs = {o: col(f"r_{o}_{kind}") for o in OBSERVERS}
            S7[kind] = {"medians": {f"{o}/{r}": e5m.finite_median(xs[o][rl == r]) for o in OBSERVERS for r in roles},
                        "contrasts": e5m.contrast({f"delta_{o}": e5m.stratum_median_diff(xs[o], rl) for o in OBSERVERS},
                                                  ix, dl, draws)}
        sec["S7_training_maps"] = S7
    if anm is not None:
        S4 = {}
        for o in OBSERVERS:
            for _, rtol in ANM_BUDGETS:
                k = f"anm_{o}_{rtol}"
                cnt = {}
                for r in rows:
                    for fld in ("status", "source_kernel_status", "factorization_status"):
                        key_ = f"{fld}={r[k][fld]}"
                        cnt[key_] = cnt.get(key_, 0) + 1
                S4[f"{o}/rtol={rtol}"] = cnt
        sec["S4_anm_status"] = S4
    else:
        sec["S4_anm_status"] = "ANM analysis module not available (set ANM_ROOT)"
    sec["S5_quotient"] = {o: {"rank_counts": {int(k): int(v) for k, v in zip(*np.unique(col(f"quot_rank_{o}"), return_counts=True))},
                              "pr_median_by_role": {r: e5m.finite_median(col(f"quot_pr_{o}")[rl == r]) for r in roles}}
                          for o in OBSERVERS}
    sec["S6_near_null"] = {f"{o}/tol={tol}": {r: e5m.finite_median(col(f"null_share_{o}_{tol}")[rl == r]) for r in roles}
                           for o in OBSERVERS for tol in NULL_TOL_REPORTED}
    sec["S6_near_null"]["near_null_dim_median"] = {str(t): float(np.median(col(f"null_dim_{t}"))) for t in NULL_TOL_REPORTED}
    sec["S6_near_null"]["Kz_condition_median"] = e5m.finite_median(col("Kz_condition"))
    sec["S8_same_cell_fit_max"] = R["pipeline_checks"]["r_tok_same_cell_fit_max"]
    R["secondary"] = sec
    R["per_cell"] = [{k: v for k, v in r.items() if k in ("cell", "role", "cls", "donor", "m", "r_tok", "r_head", "r_rand",
                                                          "r_tok_ridge", "r_head_ridge", "r_rand_ridge")} for r in rows]
    write_json(out / "E5M_results.json", R)
    write_report(out, R)
    say(f"report: {len(rows)} cells analysed; verdict: {R['registered']['verdict']}")


def _f(x, n=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{n}f}"


def write_report(out: Path, R: dict) -> None:
    roles = sorted({k.split("/")[1] for k in R["medians"]})
    L = ["# E5-M: observer-conditioned sufficiency of TEDDY's gene-mean state", ""]
    if R["smoke_test"]:
        L += [f"SMOKE TEST: {R['smoke_note']}", ""]
    L += [f"{R['script_version']}; cells {R['n_cells_done']} of {R['n_cells_selected']} ({R['cell_pool']}); analysed "
          f"{R['n_cells_analysed']}; addendum {R['addendum_sha256'][:16]}; H file {R['H_file_sha256'][:16]}.", "",
          f"**Verdict:** {R['registered']['verdict']}", "",
          "Primary H = ANM's clamp estimator (Gamma = I, Lambda = the observer's response to clamping z_12 with the "
          "pooling-discarded token deviations held), so r = the share of the observer's first-order response that flows "
          "through what gene-mean pooling discards; r_head = 0 by construction (positive control).", "",
          f"Positive control: max r_head {R['positive_control']['r_head_max']:.2e} (gate <= {POSITIVE_CONTROL_TOL}): "
          f"{R['positive_control']['passed']}.", "",
          "## FD gate (Richardson derivative decides)", "", "| derivative / target | n | pass share | median rel. err | max rel. err | min cos |",
          "|---|---|---|---|---|---|"]
    for k, v in R["fd_gate"]["by_derivative_target"].items():
        L.append(f"| {k} | {v['n']} | {_f(v['pass_share'])} | {_f(v['median_rel_err'], 4)} | {_f(v['max_rel_err'], 4)} | "
                 f"{_f(v['min_cos'], 4)} |")
    L += ["", f"Validated: {R['fd_gate']['validated']}. Pipeline checks: `{json.dumps(R['pipeline_checks'])}`", "",
          "## Clamp-H residual r = ||K_O - H K_z||_F / ||K_O||_F (median per role)", "",
          "| observer | " + " | ".join(roles) + " |", "|---|" + "---|" * len(roles)]
    for o in OBSERVERS:
        L.append(f"| O_{o} | " + " | ".join(_f(R["medians"].get(f"{o}/{r}")) for r in roles) + " |")
    C = R["registered"].get("contrasts", {})
    if C:
        L += ["", "## Registered contrasts (two-stage 95% intervals, per donor)", "", "| statistic | value | interval | per donor |",
              "|---|---|---|---|"]
        for k, v in C.items():
            ci = v["ci_two_stage"]
            L.append(f"| {k} | {_f(v['value'])} | " + (f"[{_f(ci[0])}, {_f(ci[1])}]" if ci else "n/a") + " | "
                     + ", ".join(f"{d_} {_f(x)}" for d_, x in v["per_donor"].items()) + " |")
    L += ["", "## Secondary (never change the verdict)", ""]
    for k, v in R["secondary"].items():
        L.append(f"- **{k}**: `" + json.dumps(v, default=_json_default)[:2000] + "`")
    (out / "REPORT.md").write_text("\n".join(L) + "\n")


# ============================================================================ main
def main(argv=None) -> int:
    a = parse_args(argv)
    if a.stage == "commit-addendum":
        stage_commit_addendum(a)
        return 0
    if a.stage == "register":
        if a.smoke:
            if a.out_dir is None:
                raise SystemExit("a smoke register needs --out-dir (it never writes the registration directory)")
            a.out_dir.mkdir(parents=True, exist_ok=True)
            _LOG["path"] = a.out_dir / "progress.log"
        stage_register(a)
        return 0
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG["path"] = a.out_dir / "progress.log"
    say(f"{SCRIPT_VERSION} stage {a.stage} out {a.out_dir}" + (f" SMOKE: {a.smoke}" if a.smoke else ""))
    site4 = a.cell_pool == "site4" and a.stage != "select"
    reg_info = check_registration(a, site4=site4)
    write_json(a.out_dir / "e5m_registration_used.json", reg_info)
    reg = load_reg(a)
    meta = lp.load_meta(a.processed, need_rna=a.stage in ("respond", "all"))
    if a.stage in ("select", "all"):
        if (a.out_dir / "e5m_cells.npz").exists() and a.stage == "all":
            say("reusing e5m_cells.npz (the selection is deterministic; delete it to rebuild)")
        else:
            stage_select(a, meta, reg)
    if a.stage in ("respond", "report", "all"):
        Hz = load_H(a)
        if a.stage in ("respond", "all"):
            rc = stage_respond(a, meta, reg, Hz)
            if rc == EXIT_INCOMPLETE:
                stage_report(a, meta, reg, reg_info, Hz)  # interim report on the finished cells
                return EXIT_INCOMPLETE
        stage_report(a, meta, reg, reg_info, Hz)
    return 0


if __name__ == "__main__":
    sys.exit(main())
