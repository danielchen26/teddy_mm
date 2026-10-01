#!/usr/bin/env python3
"""E3 (C4 + repair; Experiment 6 redesigned): NK-T look-alike pairs under the v3 key, and whether a
token-aware or nonlinear readout on frozen TEDDY recovers more of the NK-T protein gap on
sufficiency-flagged pairs than on unflagged ones, and improves NK-vs-T calls at matched coverage.

Registered design: registration/registration_v3.json -> experiments.E3, as amended by A1
(registration/amendment_A1.json: A1.1 tie-break and per-replicate selection, A1.8 within-donor pairs),
loaded with bridge_anm/lib/v3_amend.load_registration_amended(). The details the registration leaves
open are fixed in registration/addenda/E3.json (written by --stage register from train/val only and
committed before any site4 evaluation). Library code: bridge_anm/lib/v3_e3.py, v3_e3_readouts.py.

  pairs     k = 10 cosine neighbours of the raw official z within each test_primary donor; NK-T pair = one
            primary-key NK and one primary-key T cell; flagged if cosine >= 0.949402 (registration e3, val)
  readouts  head (registered phase-1 head), R1 (MLP on the 12 gene-mean layer outputs), R2 (one learned
            query over layer-12 gene-token states, 4,000 training cells), null (phase-1 head retrained with
            seed 1 by scripts/04_train.py); each normalised by the training q95 of its own predictions
  H3a       D_R = [GR_R - GR_head](flagged) - [GR_R - GR_head](unflagged) - D_null >= 0.05 (win rules)
  H3b       NK-vs-T selective accuracy at matched coverage 0.9 / 0.8: Q1 rule on head vs R1 / R2 evidence
            (top score), margin 0.01; trust comparators: head margin, entropy, kNN label disagreement

Stages (each resumable; every stage appends to <out-dir>/progress.log):
  register    write registration/addenda/E3.json + its line in addenda/HASHES.txt (train/val rows only)
  null        scripts/04_train.py --seed 1 (MPS, ~7 min), unchanged settings; skipped when done
  r1          R1 on CPU (training cells, early stopping on val); predictions on train + val
  r2_states   layer-12 token states of R2's 4,000 training + 1,000 val cells (MPS, ~10 min, ~6.5 GB fp16)
  r2_fit      R2 on CPU from the stored states
  r2_predict  stream val, train and site4 cells through TEDDY + R2 (MPS, ~3 h; shards of 1,000 cells)
  evaluate    pairs, gap ratios, H3a, H3b, bootstrap; E3_results.json + REPORT.md (CPU, minutes)
  all         null, r1, r2_states, r2_fit, r2_predict, evaluate
A stage that touches site4 cells (r2_predict on site4, evaluate) is refused unless registration_v3.json,
amendment A1 and addenda/E3.json are committed with matching hashes. Smoke runs (--smoke NOTE) use val
donor 18303 as the evaluated split and never read a site4 row.
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
import contextlib  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "bridge_anm"))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e3 as e3  # noqa: E402
from lib import v3_key as vk  # noqa: E402

SCRIPT_VERSION = "v3_e3_nkt_repair 1.0"
EXIT_INCOMPLETE = 75
N_BOOT = 2000
COVERAGES = (0.9, 0.8)
READOUTS = ("head", "R1", "R2", "null")
TRUST = ("head_margin", "head_entropy", "head_knn")
STATE_SHARD = 250
PRED_SHARD = 1000
OFFICIAL_SEC_PER_CELL = 10615.457 / 90261  # z_rna_manifest.json: official run, MPS, fp16, batch 32
MAIN = Path("/Users/tianchichen/Documents/GitHub/teddy_mm")

# ============================================================================ addendum (declared part)
# Fixed before any site4 evaluation; registration/addenda/E3.json = this dict + the train/val numbers that
# --stage register computes. Nothing here replaces a registered endpoint, key or cell set (A1.10).
DECLARED = {
    "scope": "fixes only details registration_v3.json experiments.E3 (as amended by A1) leaves open, by the "
             "procedures written there; the registered endpoints (E3.H3a, E3.H3a_pooling, E3.H3b), margins (0.05, "
             "0.01), comparators, pairs rule (A1.8), flag cosine and falsification sentence are unchanged",
    "splits": "verdicts on test_primary (donors 13272, 19593); test_secondary (donor 15078) is reported with the "
              "same rules and never changes a verdict",
    "variants": "'all' = registration primary key (primary, decides the verdicts); 'no_gdT158' = sensitivity key "
                "primary_no_gdT158, reported with the same rules, never changes a verdict",
    "pairs": {
        "function": "bridge_anm/lib/knn_cosine.knn_cosine (the function the registration builder used for the val flag) "
                    "on the raw official z rows of one donor of the evaluated split (it L2-normalises), k = 10",
        "edge_cosine": "the float32 cosine knn_cosine returns; each unordered edge once (value of its last listing, as "
                       "the builder's dict)",
        "nkt_pair": "one cell keyed NK and the other keyed T under the variant's key, same donor (A1.8)",
        "flag": "cosine >= registration e3.flag_cosine (0.949402)",
        "check": "the same function on val (one donor) must reproduce the registration's 171 val NK-T pairs and median "
                 "0.949402 (computed_train_val.val_pairs)",
    },
    "targets": ["CD20", "CD22", "CD268", "CD3", "CD2", "CD5", "CD122", "CD94", "CD56", "CD172a", "CD11c", "CD62P",
                "CD335"],
    "targets_note": "registered 'measured panel proteins + CD335 and CD3': the 12 primary-panel proteins (CD3 is one "
                    "of them) plus CD335",
    "measured_evidence": "stored measured ADT value / its 95th percentile over split=train cells (np.percentile), "
                         "clipped to [0, 1]; the 13 values are in computed_train_val.q95_measured_train (6 decimals, "
                         "used rounded)",
    "readout_evidence": "prediction / the 95th percentile of the readout's own predictions over all 67,405 split=train "
                        "cells (np.percentile), clipped to [0, 1] as the registered evidence. head: the registered "
                        "q95_train_pred (registration evidence.teddy_head); R1, R2, null: computed by their stages from "
                        "training cells only. Sensitivities (reported, never a verdict): unclipped evidence; R2 "
                        "normalised by its 4,000 training-subset cells",
    "head_and_null_predictions": "v3_key.head_predict (L2-normalised official z, phase-1 MLP + NB decoder mean, "
                                 "constant size factor 1.0)",
    "R1": {
        "input": "the 12 gene-mean layer outputs of data/processed/cite_official/z_rna_layer_means.npy (layers 1-12 "
                 "in order, concatenated, 6,144 values), float32, standardised per feature with the split=train "
                 "mean and SD (SD floor 1e-6)",
        "architecture": "Linear(6144, 256) - ReLU - Dropout(0.1) - Linear(256, 256) - ReLU - Dropout(0.1) - "
                        "Linear(256, 13)",
        "training": "all split=train cells; MSE on the 13 measured-evidence targets; Adam lr 1e-3 (default betas, no "
                    "weight decay); batch 256, reshuffled every epoch by one torch generator seeded seeds.global; "
                    "torch.manual_seed(seeds.global) before init; at most 30 epochs; after each epoch the MSE on all "
                    "6,106 val cells; stop after 5 epochs without a new best and keep the best-val weights; CPU",
    },
    "R2": {
        "cells": "rng = numpy default_rng(seeds.e3_token_subset = 13); in this order, without replacement: 1,333 "
                 "training cells with primary key NK, 1,333 with primary key T, 1,334 other training cells (primary "
                 "key B, myeloid, OUT or unscored); then 1,000 val cells uniformly from all 6,106 val cells (early "
                 "stopping); lists and hashes in computed_train_val.r2_cells",
        "states": "the last encoder layer's output (layer 12) at every real token, TEDDY-G official preprocessing as "
                  "z_rna_manifest.json (counts/total x 1e4 / gene medians, torch.topk 2,048, zeros dropped, bool "
                  "padding mask), torch.autocast fp16 on MPS, batch 32 length-bucketed, stored as float16; check: the "
                  "gene-mean of the states against z_rna.npy (cosine >= 0.99 per cell, token counts equal)",
        "architecture": "h~_t = (h_t - mu) / sd with mu, sd = per-dimension mean / SD (floor 1e-6) of the 4,000 "
                        "training cells' gene-means (a fixed affine map, a reparametrisation); s_t = <q, h~_t> / "
                        "sqrt(512) over real tokens; a = softmax(s); y = W (sum_t a_t h~_t) + b, 13 outputs; q "
                        "initialised to 0 (uniform attention = the gene-mean), W, b PyTorch default init after "
                        "torch.manual_seed(seeds.global)",
        "training": "MSE on the 13 measured-evidence targets; Adam lr 1e-3; batch 32 cells reshuffled every epoch by "
                    "one torch generator seeded seeds.global; at most 50 epochs; early stopping on the 1,000 val cells' "
                    "MSE, patience 10, keep the best; CPU",
        "predictions": "the same TEDDY forward and R2 on every val, train and site4 cell; the input states are rounded "
                       "to float16 first, as the stored training states",
    },
    "null": "scripts/04_train.py --processed data/processed/cite_official --seed 1 --device mps, every other setting "
            "its default (the official head's command, rerun_official.sh step_train_official, with seed 1); its "
            "best.pt (04_train's val_fm_pearson rule, as the official head). 04_train.py ends by printing aggregate "
            "site4 Pearson metrics for this checkpoint; E3 never reads them",
    "gap_ratio": "per protein p in CD56, CD94, CD335, CD3 and a set of pairs: numpy median over the pairs of |e_R,p(NK "
                 "cell) - e_R,p(T cell)| / numpy median over the same pairs of |m_p(NK cell) - m_p(T cell)|; GR = mean "
                 "over the 4 proteins; undefined (NaN) without pairs or with a zero measured median",
    "H3a": "D_R as registered for R in R1, R2 (algebraically [GR_R - GR_null](flagged) - [GR_R - GR_null](unflagged)); "
           "pooling = D_R2 - D_R1; margin 0.05 with the registration win / loss / equivalent rules; the per-donor "
           "condition uses each primary donor's own pairs",
    "H3b": "cells = key-NK and key-T cells of the split under the variant's key (N = their number); each readout's "
           "Q1 rule on its evidence (equal-weight mean over each primary-panel class, argmax over B, T, NK, myeloid) "
           "gives the call; confidence = top rule score; at coverage c in {0.9, 0.8} the ceil(c * N) most confident "
           "cells, ties by the A1.1 permutation (v3_amend.tie_break_rank of the global cell id); NK-vs-T selective "
           "accuracy = share of those cells whose call equals the key (a B or myeloid call is wrong). Primary "
           "comparisons R1 - head and R2 - head, margin 0.01, one verdict per coverage point (no aggregation)",
    "trust_comparators": "the head's calls re-ranked by (a) margin S_top1 - S_top2, (b) entropy score 1 - H(q) / log 4, "
                         "q = S / sum(S) (uniform if the sum is 0), (c) kNN agreement 1 - d, d = share of the cell's 10 "
                         "nearest scored training cells (cosine on the L2-normalised official z; reference = split=train "
                         "cells whose primary key is B, T, NK, myeloid or OUT) whose key differs from the call. "
                         "Reported: each comparator's selective accuracy and R1 / R2 (top score) minus each comparator "
                         "with the same rules (secondary). The null with its top score is a descriptive row",
    "bootstrap": "registration common.statistics (A1.1): two-stage, donors with replacement then pairs (H3a) or cells "
                 "(H3b) with replacement within each drawn donor; B = 2000; seed = seeds.bootstrap (1); every quantity "
                 "(medians, matched-coverage selections) recomputed in each replicate; the same draws for every "
                 "readout; a replicate in which a quantity is undefined (e.g. no flagged pair drawn) is dropped for "
                 "that quantity and the number of defined replicates is reported",
    "falsification": "registered: 'if neither R1 nor R2 reaches D_R >= 0.05 over the null, the claim is rejected'. "
                     "Evaluated on test_primary, variant all, as: rejected unless E3.H3a for R1 or for R2 is a win; the "
                     "point-only reading (D_R >= 0.05) is reported beside it",
    "secondary_H3a": "pre-specified descriptives (point, two-stage 95% interval, per donor; no verdict, never change a "
                     "registered verdict), added because the registered D is not scale-free: a readout that shrinks every "
                     "NK-T difference by a common factor s < 1 gets D_R = (1 - s) (GR_head(unflagged) - GR_head(flagged)) "
                     "> 0 whenever the head loses more of the gap on flagged pairs, and the seed-1 null does not remove "
                     "this. (S1) Dlog_R = [ln GR_R(fl) - ln GR_head(fl)] - [ln GR_R(unfl) - ln GR_head(unfl)] - (same for "
                     "the null), zero for any uniform rescaling; (S2) the two terms of D, GR_R(fl) - GR_head(fl) and "
                     "GR_R(unfl) - GR_head(unfl), for R in R1, R2, null: a repair in the registered sense needs the first "
                     "term > 0, i.e. more of the measured gap on flagged pairs than the head",
    "smoke": "smoke runs evaluate val donor 18303 only, with small subsets, and are labelled smoke; they never read a "
             "site4 row and never decide anything",
}


def addendum_bytes(addendum: dict) -> bytes:
    return (json.dumps(addendum, indent=1, sort_keys=True) + "\n").encode()


# ============================================================================ helpers
_LOG: Path | None = None


def say(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _LOG is not None:
        with _LOG.open("a") as f:
            f.write(line + "\n")


def write_json(path: Path, obj) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_json_default))
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def atomic_savez(path: Path, **arrays) -> None:
    tmp = path.with_name(path.name + ".tmp.npz")
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


def atomic_save_npy(path: Path, arr: np.ndarray) -> None:
    tmp = path.with_name(path.name + ".tmp.npy")
    np.save(tmp, arr)
    os.replace(tmp, path)


def fnum(x, n=4) -> float | None:
    if x is None:
        return None
    x = float(x)
    return None if not math.isfinite(x) else round(x, n)


def sha256_file(path: Path) -> str:
    return vk.sha256_file(path)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", choices=("register", "null", "r1", "r2_states", "r2_fit", "r2_predict", "evaluate", "all"),
                   default="all")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--embed-dir", type=Path, default=ROOT / "data/processed/cite_official")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M", help="TEDDY-G checkpoint folder")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--head-ckpt", type=Path, default=MAIN / "outputs/cite_phase1_official/best.pt")
    p.add_argument("--null-ckpt", type=Path, default=None, help="default <out-dir>/null_head_seed1/best.pt")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--device", default="auto", help="TEDDY forward device (r2_states, r2_predict)")
    p.add_argument("--null-device", default="mps", help="04_train.py device for the null (registered: mps)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--batch-size", type=int, default=32, help="TEDDY forward batch (official run: 32)")
    p.add_argument("--time-budget-sec", type=float, default=0.0, help="r2_states / r2_predict stop after a shard (exit 75)")
    p.add_argument("--n-boot", type=int, default=N_BOOT)
    p.add_argument("--smoke", default=None, metavar="NOTE", help="smoke run: val only, small subsets")
    p.add_argument("--smoke-train-cells", type=int, default=3000, help="smoke: R1 training cells")
    p.add_argument("--smoke-epochs", type=int, default=2)
    p.add_argument("--smoke-r2-cells", type=int, default=60, help="smoke: R2 subset size (train) ; val = a third")
    p.add_argument("--smoke-predict-train", type=int, default=64, help="smoke: train cells for R2's q95")
    p.add_argument("--smoke-predict-val-extra", type=int, default=60, help="smoke: val cells beyond the pair cells")
    a = p.parse_args(argv)
    if a.stage != "register" and a.out_dir is None:
        p.error("--out-dir is required")
    if a.n_boot != N_BOOT and not a.smoke:
        p.error("--n-boot other than the registered 2000 is for smoke runs only")
    if a.null_ckpt is None and a.out_dir is not None:
        a.null_ckpt = a.out_dir / "null_head_seed1" / "best.pt"
    return a


# ============================================================================ registration checks
def load_reg(a) -> dict:
    return va.load_registration_amended(a.registration_dir / "registration_v3.json",
                                        a.registration_dir / "amendment_A1.json")


def _git_committed(path: Path) -> bool | None:
    try:
        r1 = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--error-unmatch", str(path)], capture_output=True)
        r2 = subprocess.run(["git", "-C", str(ROOT), "diff", "--quiet", "HEAD", "--", str(path)], capture_output=True)
        return r1.returncode == 0 and r2.returncode == 0
    except OSError:
        return None


def _git_head_text(path: Path) -> str:
    try:
        rel = path.resolve().relative_to(ROOT.resolve())
        r = subprocess.run(["git", "-C", str(ROOT), "show", f"HEAD:{rel.as_posix()}"], capture_output=True, text=True)
        return r.stdout if r.returncode == 0 else ""
    except (OSError, ValueError):
        return ""


def check_registration(a, *, site4: bool) -> dict:
    """Hashes and commit state of the registration, A1 and the E3 addendum; refuse a site4 stage unless all hold."""
    rd = a.registration_dir
    info: dict = {}
    reg_f, amd_f = rd / "registration_v3.json", rd / "amendment_A1.json"
    add_f, hf = rd / "addenda" / "E3.json", rd / "addenda" / "HASHES.txt"
    info["registration_sha256"] = sha256_file(reg_f)
    info["registration_hash_file_matches"] = (rd / "registration_v3.json.sha256").read_text().split()[0] == info["registration_sha256"]
    info["registration_committed"] = _git_committed(reg_f)
    info["amendment_A1_sha256"] = sha256_file(amd_f)
    info["amendment_hash_file_matches"] = (rd / "amendment_A1.json.sha256").read_text().split()[0] == info["amendment_A1_sha256"]
    info["amendment_committed"] = _git_committed(amd_f)
    info["addendum_file"] = str(add_f)
    if add_f.exists():
        add = json.loads(add_f.read_text())
        info["addendum_sha256"] = sha256_file(add_f)
        info["addendum_hash_recorded"] = hf.exists() and f"{info['addendum_sha256']}  E3.json" in hf.read_text().splitlines()
        info["addendum_hash_committed"] = f"{info['addendum_sha256']}  E3.json" in _git_head_text(hf).splitlines()
        info["addendum_committed"] = _git_committed(add_f)
        info["addendum_declared_equals_script"] = add.get("declared") == DECLARED
        info["addendum_amends"] = {"registration": add.get("registration_sha256") == info["registration_sha256"],
                                   "amendment_A1": add.get("amendment_A1_sha256") == info["amendment_A1_sha256"]}
    else:
        info.update({"addendum_sha256": None, "addendum_hash_recorded": False, "addendum_hash_committed": False,
                     "addendum_committed": False,
                     "addendum_declared_equals_script": False, "addendum_amends": {}})
    if site4:
        need = {
            "registration_v3.json hash file matches and committed": info["registration_hash_file_matches"] and info["registration_committed"] is True,
            "amendment_A1.json hash file matches and committed": info["amendment_hash_file_matches"] and info["amendment_committed"] is True,
            "addenda/E3.json present and committed, its sha256 in the committed addenda/HASHES.txt":
                bool(info["addendum_sha256"]) and info["addendum_hash_recorded"] and info["addendum_hash_committed"]
                and info["addendum_committed"] is True,
            "addendum declared part equals this script's DECLARED": info["addendum_declared_equals_script"],
            "addendum amends the registration and A1 on disk": all(info["addendum_amends"].values()) if info["addendum_amends"] else False,
        }
        bad = [k for k, ok in need.items() if not ok]
        if bad:
            raise SystemExit("site4 stage refused (register and commit the E3 addendum first): " + "; ".join(bad))
    return info


def load_addendum(a) -> dict:
    f = a.registration_dir / "addenda" / "E3.json"
    if not f.exists():
        raise SystemExit(f"{f} missing: run --stage register first")
    add = json.loads(f.read_text())
    if add.get("declared") != DECLARED:
        raise SystemExit("addenda/E3.json declared part differs from this script's DECLARED")
    return add


# ============================================================================ data
class Data:
    """Arrays of cite_arrays.npz plus the official embedding. Site4 rows of adt / cell_types are only read by
    ``keys_for`` / ``measured`` when the caller passes site4 rows (evaluate, after the registration check)."""

    def __init__(self, a, need_rna: bool = False):
        npz = np.load(a.processed / "cite_arrays.npz", allow_pickle=False)
        self.split = npz["split"].astype(str)
        self.sites = npz["sites"].astype(str)
        self.donors = npz["donors"].astype(str)
        self.adt_names = [str(x) for x in npz["adt_names"]]
        self._npz = npz
        self._adt = None
        self._ct = None
        self.n = int(self.split.size)
        self.z = np.load(a.embed_dir / "z_rna.npy", mmap_mode="r")
        self.ntokens = np.load(a.embed_dir / "z_rna_ntokens.npy")
        self.manifest = json.loads((a.embed_dir / "z_rna_manifest.json").read_text())
        if self.z.shape[0] != self.n:
            raise SystemExit("z rows != cells")
        self.rna = None
        self.token_ids = npz["token_ids"]
        self.rna_names = npz["rna_names"]
        if need_rna:
            from scipy import sparse
            self.rna = sparse.csr_matrix((npz["rna_data"], npz["rna_indices"], npz["rna_indptr"]),
                                         shape=tuple(npz["rna_shape"]))

    def rows(self, name: str) -> np.ndarray:
        return np.where(self.split == name)[0]

    @property
    def adt(self) -> np.ndarray:
        if self._adt is None:
            self._adt = np.asarray(self._npz["adt"], dtype=np.float32)
        return self._adt

    @property
    def cell_types(self) -> np.ndarray:
        if self._ct is None:
            self._ct = self._npz["cell_types"].astype(str)
        return self._ct

    def keys_for(self, rows: np.ndarray, reg: dict) -> dict[str, np.ndarray]:
        return keys_for(self, rows, reg)


class TV:
    """Plain in-memory arrays (the register stage's input; a copy can be poisoned for the leakage check)."""

    def __init__(self, split, donors, adt, cell_types, z, ntokens, adt_names):
        self.split, self.donors, self.adt, self.cell_types = split, donors, adt, cell_types
        self.z, self.ntokens, self.adt_names = z, ntokens, list(adt_names)
        self.n = int(split.size)

    def rows(self, name: str) -> np.ndarray:
        return np.where(self.split == name)[0]


def keys_for(X, rows: np.ndarray, reg: dict) -> dict[str, np.ndarray]:
    """Registered keys (verifier) for the given rows of X; full-length arrays, 'unscored' elsewhere."""
    rows = np.asarray(rows, dtype=np.int64)
    k = vk.build_keys(X.cell_types[rows], X.adt[rows], X.adt_names, reg)
    out = {}
    for name in ("primary", "primary_no_gdT158"):
        full = np.full(X.n, vk.UNSCORED, dtype=object)
        full[rows] = k[name]
        out[name] = full.astype(str)
    return out


def target_cols(names: list[str]) -> list[int]:
    j = {n: i for i, n in enumerate(names)}
    return [j[p] for p in DECLARED["targets"]]


def measured_evidence(D: Data, rows: np.ndarray, q95m: np.ndarray) -> np.ndarray:
    return e3.normalise(D.adt[np.asarray(rows)][:, target_cols(D.adt_names)], q95m)


def panel_in_targets(reg: dict) -> dict[str, list[str]]:
    pan = vk.question_panel(reg, "Q1")
    for c, ps in pan.items():
        for p_ in ps:
            if p_ not in DECLARED["targets"]:
                raise SystemExit(f"panel protein {p_} not among the E3 targets")
    return pan


# ============================================================================ register
def r2_cells(D: Data, keys: dict[str, np.ndarray], reg: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(int(reg["seeds"]["e3_token_subset"]))
    tr = D.rows("train")
    prim = keys["primary"]
    nk_pool = tr[prim[tr] == "NK"]
    t_pool = tr[prim[tr] == "T"]
    ot_pool = tr[~np.isin(prim[tr], ["NK", "T"])]
    nk = rng.choice(nk_pool, 1333, replace=False)
    t = rng.choice(t_pool, 1333, replace=False)
    ot = rng.choice(ot_pool, 1334, replace=False)
    vrows = rng.choice(D.rows("val"), 1000, replace=False)
    train_sub = np.sort(np.concatenate([nk, t, ot])).astype(np.int64)
    val_sub = np.sort(vrows).astype(np.int64)
    info = {"train": [int(x) for x in train_sub], "val": [int(x) for x in val_sub],
            "train_sha256": vk.index_hash(train_sub), "val_sha256": vk.index_hash(val_sub),
            "train_counts": {"NK": 1333, "T": 1333, "other": 1334},
            "pool_sizes": {"NK": int(nk_pool.size), "T": int(t_pool.size), "other": int(ot_pool.size), "val": int(D.rows("val").size)},
            "ntokens_median": {"train": float(np.median(D.ntokens[train_sub])), "val": float(np.median(D.ntokens[val_sub]))},
            "ntokens_total": int(D.ntokens[train_sub].sum() + D.ntokens[val_sub].sum())}
    return train_sub, val_sub, info


def compute_train_val(a, reg: dict, X) -> dict:
    """Every number the addendum adds, from split in {train, val} rows only (X may have poisoned site4 rows)."""
    nt = np.where(X.split != "test")[0]
    keys = keys_for(X, nt, reg)
    tr, vr = X.rows("train"), X.rows("val")
    P = e3.nkt_pairs_within_donors(X.z, vr, X.donors, {"all": keys["primary"]}, 10, reg["e3"]["flag_cosine"])
    pv = P["variants"]["all"]
    med = float(round(float(np.median(pv["cos"])), 6)) if pv["cos"].size else None
    q95m = np.percentile(np.asarray(X.adt[tr])[:, target_cols(X.adt_names)].astype(np.float64), 95, axis=0)
    pred_tr = vk.head_predict(np.asarray(X.z[tr], dtype=np.float32), a.head_ckpt, device="cpu", size_factor=1.0)
    reg_q = np.asarray(reg["evidence"]["teddy_head"]["q95_train_pred"])
    head_q_diff = float(np.max(np.abs(np.round(np.percentile(pred_tr, 95, axis=0), 6) - reg_q)))
    _, _, rinfo = r2_cells(X, keys, reg)
    scored_train = tr[keys["primary"][tr] != vk.UNSCORED]
    return {
        "val_pairs": {"n_val_nkt_pairs": int(pv["nk"].size), "median_cosine": med,
                      "n_flagged_at_flag_cosine": int(pv["flag"].sum()), "n_val_edges": int(sum(P["n_edges"].values())),
                      "reproduces_registration": bool(pv["nk"].size == int(reg["e3"]["n_val_nkt_pairs"])
                                                      and med == float(reg["e3"]["flag_cosine"]))},
        "q95_measured_train": dict(zip(DECLARED["targets"], [float(round(float(x), 6)) for x in q95m])),
        "head_q95_reproduction_max_abs_diff": round(head_q_diff, 6),
        "r2_cells": rinfo,
        "knn_reference": {"n_scored_training_cells": int(scored_train.size), "sha256": vk.index_hash(scored_train)},
        "n_train": int(tr.size), "n_val": int(vr.size),
    }


def _poisoned(X: TV, which: str, seed: int) -> TV:
    """Copy of X whose rows of one split carry random protein, cell types and embedding (leakage check)."""
    rng = np.random.default_rng(seed)
    m = X.split == which
    adt, ct, z, nt = X.adt.copy(), X.cell_types.copy(), np.array(X.z, dtype=np.float32, copy=True), X.ntokens.copy()
    adt[m] = rng.gamma(1.0, 1.0, size=(int(m.sum()), adt.shape[1])).astype(adt.dtype)
    ct[m] = rng.choice(np.unique(X.cell_types), size=int(m.sum()))
    z[m] = rng.normal(size=(int(m.sum()), z.shape[1])).astype(np.float32)
    nt[m] = rng.integers(50, 2049, size=int(m.sum()))
    return TV(X.split, X.donors, adt, ct, z, nt, X.adt_names)


def stage_register(a) -> None:
    reg = load_reg(a)
    D = Data(a)
    t0 = time.time()
    X = TV(D.split, D.donors, D.adt, D.cell_types, np.asarray(D.z, dtype=np.float32), D.ntokens, D.adt_names)
    comp = compute_train_val(a, reg, X)
    vp = comp["val_pairs"]
    if not vp["reproduces_registration"]:
        raise SystemExit(f"val pairs do not reproduce the registration: {vp}")
    if comp["head_q95_reproduction_max_abs_diff"] > 1e-5:
        raise SystemExit(f"head q95 not reproduced (max diff {comp['head_q95_reproduction_max_abs_diff']})")
    # leakage check (as v3_leakage_check.py for the registration): site4 rows poisoned -> identical numbers;
    # val rows poisoned -> the numbers change (positive control)
    same = json.dumps(compute_train_val(a, reg, _poisoned(X, "test", 0)), sort_keys=True) == json.dumps(comp, sort_keys=True)
    changed = json.dumps(compute_train_val(a, reg, _poisoned(X, "val", 0)), sort_keys=True) != json.dumps(comp, sort_keys=True)
    if not (same and changed):
        raise SystemExit(f"leakage check failed: site4-poisoned identical {same}, val-poisoned changed {changed}")
    add = {
        "addendum_id": "E3", "experiment": "E3", "builder": SCRIPT_VERSION,
        "written_before": "any site4 evaluation of E3; --stage register reads no site4 value (leakage_check)",
        "registration_sha256": sha256_file(a.registration_dir / "registration_v3.json"),
        "amendment_A1_sha256": sha256_file(a.registration_dir / "amendment_A1.json"),
        "declared": DECLARED,
        "computed_train_val": comp,
        "leakage_check": {"site4_rows_poisoned_numbers_identical": same, "val_rows_poisoned_numbers_change": changed,
                          "poisoned": "protein, cell types, embedding and token counts of the rows, seed 0"},
    }
    d = a.registration_dir / "addenda"
    d.mkdir(parents=True, exist_ok=True)
    f = d / "E3.json"
    f.write_bytes(addendum_bytes(add))
    sha = sha256_file(f)
    hf = d / "HASHES.txt"
    lines = [ln for ln in (hf.read_text().splitlines() if hf.exists() else []) if not ln.rstrip().endswith("  E3.json")]
    hf.write_text("\n".join(lines + [f"{sha}  E3.json"]) + "\n")
    print(f"wrote {f} (sha256 {sha}) and its line in {hf}; val pairs {vp['n_val_nkt_pairs']}, median "
          f"{vp['median_cosine']}, leakage check pass ({time.time() - t0:.0f}s)")


# ============================================================================ null
def stage_null(a) -> None:
    out = a.out_dir / "null_head_seed1"
    if (out / "metrics.json").exists() and (out / "best.pt").exists():
        say(f"null: {out}/best.pt exists, skipped")
        return
    if a.smoke:
        say("null: smoke run, 04_train.py is not run (it evaluates site4 at its end); --null-ckpt stands in")
        return
    check_registration(a, site4=True)  # 04_train.py prints aggregate site4 metrics when it finishes
    out.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(ROOT / "scripts" / "04_train.py"), "--processed", str(a.embed_dir), "--out", str(out),
           "--seed", "1", "--device", a.null_device]
    say("null: " + " ".join(cmd))
    t0 = time.time()
    with (out / "train.log").open("a") as lf:
        r = subprocess.run(cmd, cwd=str(ROOT), stdout=lf, stderr=subprocess.STDOUT)
    if r.returncode != 0 or not (out / "best.pt").exists():
        raise SystemExit(f"null training failed (exit {r.returncode}); see {out / 'train.log'}")
    say(f"null: done in {time.time() - t0:.0f}s")


# ============================================================================ R1
def r1_inputs(lm, rows: np.ndarray) -> np.ndarray:
    rows = np.asarray(rows, dtype=np.int64)
    o = np.argsort(rows)
    srt = rows[o]
    X = np.concatenate([np.asarray(lm[l][srt], dtype=np.float32) for l in range(lm.shape[0])], axis=1)
    out = np.empty_like(X)
    out[o] = X
    return out


def load_r1(path: Path):
    import torch

    from lib.v3_e3_readouts import R1MLP
    blob = torch.load(str(path), map_location="cpu", weights_only=False)
    m = R1MLP(int(blob["d_in"]), int(blob["d_out"]))
    m.load_state_dict(blob["state_dict"])
    m.eval()
    return m


def r1_predict(model, lm, rows: np.ndarray, batch: int = 4096) -> np.ndarray:
    import torch
    out = []
    with torch.no_grad():
        for s in range(0, len(rows), batch):
            out.append(model(torch.from_numpy(r1_inputs(lm, rows[s:s + batch]))).numpy())
    return np.concatenate(out) if out else np.zeros((0, len(DECLARED["targets"])), np.float32)


def stage_r1(a, reg, add) -> None:
    import torch

    from lib.v3_e3_readouts import R1MLP, fit_with_early_stopping
    torch.set_num_threads(_THREADS)
    od = a.out_dir / "r1"
    od.mkdir(parents=True, exist_ok=True)
    if (od / "r1_done.json").exists():
        say("r1: done, skipped")
        return
    D = Data(a)
    lm = np.load(a.embed_dir / "z_rna_layer_means.npy", mmap_mode="r")
    q95m = np.asarray([add["computed_train_val"]["q95_measured_train"][p] for p in DECLARED["targets"]])
    tr, vr = D.rows("train"), D.rows("val")
    if a.smoke:
        tr = np.sort(np.random.default_rng(0).choice(tr, min(a.smoke_train_cells, tr.size), replace=False))
    seed = int(reg["seeds"]["global"])
    t0 = time.time()
    Xtr, Xva = r1_inputs(lm, tr), r1_inputs(lm, vr)
    Ytr, Yva = measured_evidence(D, tr, q95m).astype(np.float32), measured_evidence(D, vr, q95m).astype(np.float32)
    say(f"r1: inputs {Xtr.shape} train, {Xva.shape} val ({time.time() - t0:.0f}s)")
    torch.manual_seed(seed)
    model = R1MLP(Xtr.shape[1], Ytr.shape[1])
    model.set_standardisation(Xtr.mean(0), Xtr.std(0))
    Xt, Yt, Xv, Yv = (torch.from_numpy(x) for x in (Xtr, Ytr, Xva, Yva))
    gen = torch.Generator().manual_seed(seed)

    def batches(ep):
        perm = torch.randperm(Xt.shape[0], generator=gen)
        for s in range(0, Xt.shape[0], 256):
            ix = perm[s:s + 256]
            yield Xt[ix], Yt[ix]

    def loss_on(m, b):
        return torch.mean((m(b[0]) - b[1]) ** 2)

    def val_mse(m):
        return float(torch.mean((m(Xv) - Yv) ** 2))

    fit = fit_with_early_stopping(model, batches, loss_on, val_mse, lr=1e-3,
                                  max_epochs=a.smoke_epochs if a.smoke else 30, patience=5, log=say)
    torch.save({"state_dict": model.state_dict(), "d_in": Xtr.shape[1], "d_out": Ytr.shape[1],
                "targets": DECLARED["targets"], "train_rows_sha256": vk.index_hash(tr)}, od / "r1.pt")
    rows = np.concatenate([D.rows("train"), vr])
    P = np.full((D.n, Ytr.shape[1]), np.nan, np.float32)
    P[rows] = r1_predict(model, lm, rows)
    atomic_save_npy(od / "pred_trainval.npy", P)
    q = e3.q95(P[D.rows("train")])
    pv = P[vr]
    pear = {p_: fnum(np.corrcoef(pv[:, j], Yva[:, j])[0, 1]) for j, p_ in enumerate(DECLARED["targets"])}
    write_json(od / "r1_done.json", {"fit": fit, "q95_own_train": dict(zip(DECLARED["targets"], map(float, q))),
                                     "val_pearson_vs_measured_evidence": pear, "n_train": int(tr.size),
                                     "smoke": bool(a.smoke), "sec": time.time() - t0})
    say(f"r1: best epoch {fit['best_epoch']} val MSE {fit['best_val_mse']:.6f}; val Pearson CD56 {pear['CD56']} "
        f"CD3 {pear['CD3']} ({time.time() - t0:.0f}s)")


# ============================================================================ TEDDY forward
def teddy_setup(a, D: Data):
    import torch

    from teddy_mm.device import resolve_device
    from teddy_mm.teddy_encoder import load_gene_medians, load_pad_id, load_teddy, load_vocab, median_factors
    torch.set_num_threads(_THREADS)
    device = resolve_device(a.device)
    vocab = load_vocab(a.ckpt)
    pad_id = load_pad_id(a.ckpt, vocab)
    factors = median_factors(D.rna_names, load_gene_medians(a.medians))
    model = load_teddy(a.ckpt, device)
    if D.manifest.get("preprocessing") != "official" or D.manifest.get("autocast") != "fp16":
        raise SystemExit("embedding manifest is not the official fp16 run")
    return model, device, pad_id, factors


def teddy_hidden_batches(a, D: Data, model, device, pad_id, factors, cells: np.ndarray):
    """Official TEDDY-G forward (as scripts/03_embed_rna.py --preprocessing official --length-buckets --autocast
    fp16): yields (rows into cells, last-layer states [B, W, d], bool mask [B, W], token ids [B, W])."""
    import torch

    from teddy_mm.teddy_encoder import official_values, rank_encode_official
    L = int(D.manifest.get("seq_len", 2048))
    norm = float(D.manifest.get("normalize_total", 1e4))
    n = len(cells)
    tokens = np.full((n, L), pad_id, np.int64)
    attn = np.zeros((n, L), np.int64)
    for s in range(0, n, 256):
        e = min(s + 256, n)
        vals = official_values(D.rna[cells[s:e]].toarray(), factors, norm)
        tokens[s:e], attn[s:e] = rank_encode_official(vals, D.token_ids, max_len=L, pad_id=pad_id, pad_to=L)
    ntok = attn.sum(1).astype(np.int32)
    mism = int(np.sum(ntok != D.ntokens[cells]))
    if mism:
        raise SystemExit(f"{mism} cells have another token count than the official run")
    order = np.argsort(ntok, kind="stable")
    amp = torch.float16 if device.type in ("mps", "cuda") else None
    with torch.no_grad():
        for s in range(0, n, a.batch_size):
            rows = order[s:s + a.batch_size]
            width = int(max(1, ntok[rows].max()))
            ids = torch.from_numpy(tokens[rows, :width]).to(device)
            mask = torch.from_numpy(attn[rows, :width]).to(device)
            with (torch.autocast(device_type=device.type, dtype=amp) if amp is not None else contextlib.nullcontext()):
                h = model.hidden_states(ids, mask)
            yield rows, h, mask.bool(), tokens[rows, :width]


def genemean_cos(h, mask, zref: np.ndarray) -> np.ndarray:
    m = mask.unsqueeze(-1).float()
    gm = ((h.float() * m).sum(1) / m.sum(1).clamp(min=1.0)).cpu().numpy().astype(np.float64)
    zr = np.asarray(zref, dtype=np.float64)
    return (gm * zr).sum(1) / (np.linalg.norm(gm, axis=1) * np.linalg.norm(zr, axis=1) + 1e-12)


# ============================================================================ R2 states / fit / predict
def r2_subset(a, D: Data, reg: dict, add: dict) -> tuple[np.ndarray, np.ndarray]:
    rc = add["computed_train_val"]["r2_cells"]
    tsub, vsub = np.asarray(rc["train"], np.int64), np.asarray(rc["val"], np.int64)
    if vk.index_hash(tsub) != rc["train_sha256"] or vk.index_hash(vsub) != rc["val_sha256"]:
        raise SystemExit("R2 cell lists differ from their recorded hashes")
    if a.smoke:
        rng = np.random.default_rng(0)
        tsub = np.sort(rng.choice(tsub, min(a.smoke_r2_cells, tsub.size), replace=False))
        vsub = np.sort(rng.choice(vsub, max(10, a.smoke_r2_cells // 3), replace=False))
    return tsub, vsub


def budget_hit(a, t_start: float, done_now: int) -> bool:
    return bool(a.time_budget_sec) and done_now > 0 and time.time() - t_start > a.time_budget_sec


def stage_r2_states(a, reg, add) -> int:
    import torch
    sd = a.out_dir / "r2" / "states"
    sd.mkdir(parents=True, exist_ok=True)
    D = Data(a, need_rna=True)
    tsub, vsub = r2_subset(a, D, reg, add)
    cells = np.concatenate([tsub, vsub])
    shards = [cells[s:s + STATE_SHARD] for s in range(0, cells.size, STATE_SHARD)]
    todo = [k for k in range(len(shards)) if not (sd / f"shard_{k:04d}_meta.npz").exists()]
    est = OFFICIAL_SEC_PER_CELL * sum(len(shards[k]) for k in todo)
    say(f"r2_states: {cells.size} cells in {len(shards)} shards, {len(todo)} to do (~{est / 60:.1f} min at the official rate)")
    if not todo:
        return 0
    model, device, pad_id, factors = teddy_setup(a, D)
    t_start, done_now = time.time(), 0
    for k in todo:
        if budget_hit(a, t_start, done_now):
            say(f"r2_states: time budget reached before shard {k}; rerun to resume")
            return EXIT_INCOMPLETE
        c = shards[k]
        t0 = time.time()
        per = [None] * len(c)
        cos = np.zeros(len(c))
        for rows, h, mask, _ in teddy_hidden_batches(a, D, model, device, pad_id, factors, c):
            h16 = h.to(torch.float16).cpu().numpy()
            nt = mask.sum(1).cpu().numpy()
            for i, r in enumerate(rows):
                per[r] = h16[i, :nt[i]]
            cos[rows] = genemean_cos(h, mask, D.z[c[rows]])
        offs = np.concatenate([[0], np.cumsum([p_.shape[0] for p_ in per])]).astype(np.int64)
        atomic_save_npy(sd / f"shard_{k:04d}_states.npy", np.concatenate(per, axis=0))
        atomic_savez(sd / f"shard_{k:04d}_meta.npz", cells=c, offsets=offs, genemean_cos=cos,
                     addendum_sha256=np.array(sha256_file(a.registration_dir / "addenda" / "E3.json")))
        done_now += 1
        if cos.min() < 0.99:
            raise SystemExit(f"shard {k}: gene-mean of the recomputed states differs from z_rna.npy (min cos {cos.min():.4f})")
        say(f"r2_states: shard {k + 1}/{len(shards)} {len(c)} cells in {time.time() - t0:.0f}s, gene-mean cos min {cos.min():.6f}")
    return 0


class StateStore:
    def __init__(self, sd: Path):
        self.views = {}
        for mf in sorted(sd.glob("shard_*_meta.npz")):
            k = mf.name.split("_")[1]
            with np.load(mf) as m:
                cells, offs = m["cells"], m["offsets"]
            st = np.load(sd / f"shard_{k}_states.npy", mmap_mode="r")
            for i, c in enumerate(cells):
                self.views[int(c)] = (st, int(offs[i]), int(offs[i + 1]))

    def get(self, c: int) -> np.ndarray:
        st, s, e = self.views[int(c)]
        return st[s:e]


def stage_r2_fit(a, reg, add) -> None:
    import torch

    from lib.v3_e3_readouts import R2AttnPool, fit_with_early_stopping, pad_states
    torch.set_num_threads(_THREADS)
    od = a.out_dir / "r2"
    if (od / "r2_fit.json").exists():
        say("r2_fit: done, skipped")
        return
    D = Data(a)
    tsub, vsub = r2_subset(a, D, reg, add)
    S = StateStore(od / "states")
    miss = [int(c) for c in np.concatenate([tsub, vsub]) if int(c) not in S.views]
    if miss:
        raise SystemExit(f"r2_fit: {len(miss)} subset cells have no stored states (run r2_states)")
    q95m = np.asarray([add["computed_train_val"]["q95_measured_train"][p] for p in DECLARED["targets"]])
    Ytr = torch.from_numpy(measured_evidence(D, tsub, q95m).astype(np.float32))
    Yva = torch.from_numpy(measured_evidence(D, vsub, q95m).astype(np.float32))
    t0 = time.time()
    gms = np.stack([np.asarray(S.get(c), dtype=np.float32).mean(0) for c in tsub])
    seed = int(reg["seeds"]["global"])
    torch.manual_seed(seed)
    model = R2AttnPool(gms.shape[1], Ytr.shape[1])
    model.set_standardisation(gms.mean(0), gms.std(0))
    gen = torch.Generator().manual_seed(seed)

    def batches(ep):
        perm = torch.randperm(tsub.size, generator=gen).numpy()
        for s in range(0, tsub.size, 32):
            ix = perm[s:s + 32]
            H, M = pad_states([np.asarray(S.get(tsub[i])) for i in ix])
            yield H, M, Ytr[ix]

    def loss_on(m, b):
        return torch.mean((m(b[0], b[1])[0] - b[2]) ** 2)

    def predict(m, cells):
        out = []
        for s in range(0, len(cells), 64):
            H, M = pad_states([np.asarray(S.get(c)) for c in cells[s:s + 64]])
            out.append(m(H, M)[0])
        return torch.cat(out)

    def val_mse(m):
        return float(torch.mean((predict(m, vsub) - Yva) ** 2))

    fit = fit_with_early_stopping(model, batches, loss_on, val_mse, lr=1e-3,
                                  max_epochs=a.smoke_epochs if a.smoke else 50, patience=10, log=say)
    torch.save({"state_dict": model.state_dict(), "d": gms.shape[1], "d_out": Ytr.shape[1], "targets": DECLARED["targets"],
                "train_sha256": vk.index_hash(tsub), "val_sha256": vk.index_hash(vsub)}, od / "r2.pt")
    with torch.no_grad():
        pt = predict(model, tsub).numpy()
        pv = predict(model, vsub).numpy()
    atomic_savez(od / "pred_subset.npz", train_cells=tsub, train_pred=pt, val_cells=vsub, val_pred=pv)
    yv = Yva.numpy()
    pear = {p_: fnum(np.corrcoef(pv[:, j], yv[:, j])[0, 1]) for j, p_ in enumerate(DECLARED["targets"])}
    qn = float(torch.linalg.norm(model.q.detach()))
    write_json(od / "r2_fit.json", {"fit": fit, "val_pearson_vs_measured_evidence": pear, "query_norm": qn,
                                    "q95_subset_train": dict(zip(DECLARED["targets"], map(float, e3.q95(pt)))),
                                    "n_train": int(tsub.size), "n_val": int(vsub.size), "smoke": bool(a.smoke),
                                    "r2_sha256": sha256_file(od / "r2.pt"), "sec": time.time() - t0})
    say(f"r2_fit: best epoch {fit['best_epoch']} val MSE {fit['best_val_mse']:.6f}, |q| {qn:.3f}; val Pearson CD56 "
        f"{pear['CD56']} CD3 {pear['CD3']} ({time.time() - t0:.0f}s)")


def predict_plan(a, D: Data, reg: dict) -> list[tuple[str, np.ndarray]]:
    if a.smoke:
        keys = D.keys_for(D.rows("val"), reg)
        vr = D.rows("val")
        P = e3.nkt_pairs_within_donors(D.z, vr, D.donors, {"all": keys["primary"]}, 10, reg["e3"]["flag_cosine"])
        pc = np.unique(np.concatenate([P["variants"]["all"]["nk"], P["variants"]["all"]["t"]]))
        rng = np.random.default_rng(1)
        other = np.setdiff1d(vr[np.isin(keys["primary"][vr], ["NK", "T"])], pc)
        extra = rng.choice(other, min(a.smoke_predict_val_extra, other.size), replace=False)
        trs = rng.choice(D.rows("train"), a.smoke_predict_train, replace=False)
        return [("val", np.sort(np.concatenate([pc, extra]))), ("train", np.sort(trs))]
    return [("val", D.rows("val")), ("train", D.rows("train")), ("test", D.rows("test"))]


def stage_r2_predict(a, reg, add) -> int:
    import torch

    from lib.v3_e3_readouts import R2AttnPool
    od = a.out_dir / "r2"
    pd_ = od / "pred"
    pd_.mkdir(parents=True, exist_ok=True)
    D = Data(a, need_rna=True)
    plan = predict_plan(a, D, reg)
    if any(name == "test" for name, _ in plan):
        check_registration(a, site4=True)
    r2_sha = sha256_file(od / "r2.pt")
    jobs = []
    for name, cells in plan:
        for k, s in enumerate(range(0, cells.size, PRED_SHARD)):
            f = pd_ / f"{name}_{k:04d}.npz"
            if f.exists():
                with np.load(f) as old:
                    if str(old["r2_sha256"]) != r2_sha:
                        raise SystemExit(f"{f} was written by another R2 model; use a new --out-dir")
                continue
            jobs.append((f, cells[s:s + PRED_SHARD]))
    n_todo = sum(len(c) for _, c in jobs)
    say(f"r2_predict: {len(jobs)} shards ({n_todo} cells) to do; ~{OFFICIAL_SEC_PER_CELL * n_todo / 60:.0f} min at the "
        f"official rate")
    if not jobs:
        return 0
    model, device, pad_id, factors = teddy_setup(a, D)
    blob = torch.load(str(od / "r2.pt"), map_location="cpu", weights_only=False)
    r2 = R2AttnPool(int(blob["d"]), int(blob["d_out"]))
    r2.load_state_dict(blob["state_dict"])
    r2.eval().to(device)
    t_start, done_now, cells_done = time.time(), 0, 0
    for f, c in jobs:
        if budget_hit(a, t_start, done_now):
            say(f"r2_predict: time budget reached before {f.name}; rerun to resume")
            return EXIT_INCOMPLETE
        t0 = time.time()
        n = len(c)
        Y = np.zeros((n, int(blob["d_out"])), np.float32)
        cos = np.zeros(n)
        topi = np.zeros((n, 10), np.int64)
        topw = np.zeros((n, 10), np.float32)
        neff = np.zeros(n, np.float32)
        for rows, h, mask, tok in teddy_hidden_batches(a, D, model, device, pad_id, factors, c):
            with torch.no_grad():
                y, att = r2(h.to(torch.float16).float(), mask)
            Y[rows] = y.float().cpu().numpy()
            att = att.float().cpu().numpy()
            kk = min(10, att.shape[1])
            ti = np.argsort(-att, axis=1)[:, :kk]
            topw[rows, :kk] = np.take_along_axis(att, ti, axis=1)
            topi[rows, :kk] = np.take_along_axis(tok, ti, axis=1)
            neff[rows] = 1.0 / np.maximum((att ** 2).sum(1), 1e-12)
            cos[rows] = genemean_cos(h, mask, D.z[c[rows]])
        if cos.min() < 0.99:
            raise SystemExit(f"{f.name}: gene-mean of the recomputed states differs from z_rna.npy (min cos {cos.min():.4f})")
        atomic_savez(f, cells=c, pred=Y, genemean_cos=cos, top_token_ids=topi, top_weights=topw, n_eff=neff,
                     r2_sha256=np.array(r2_sha))
        done_now += 1
        cells_done += n
        rate = (time.time() - t_start) / cells_done
        say(f"r2_predict: {f.name} {n} cells in {time.time() - t0:.0f}s ({1000 * (time.time() - t0) / n:.0f} ms/cell); "
            f"ETA {rate * (n_todo - cells_done) / 60:.0f} min; gene-mean cos min {cos.min():.6f}")
    return 0


def load_r2_pred(od: Path, n: int, splits: tuple[str, ...]) -> tuple[np.ndarray, dict]:
    P = np.full((n, len(DECLARED["targets"])), np.nan, np.float32)
    extra = {"n_eff": np.full(n, np.nan, np.float32), "top_token_ids": np.full((n, 10), -1, np.int64),
             "genemean_cos_min": 1.0}
    for f in sorted((od / "pred").glob("*.npz")):
        if f.name.split("_")[0] not in splits:
            continue
        with np.load(f) as d:
            c = d["cells"]
            P[c] = d["pred"]
            extra["n_eff"][c] = d["n_eff"]
            extra["top_token_ids"][c] = d["top_token_ids"]
            extra["genemean_cos_min"] = min(extra["genemean_cos_min"], float(d["genemean_cos"].min()))
    return P, extra


# ============================================================================ evaluate
def readout_evidence(a, D: Data, reg: dict, add: dict, eval_rows: np.ndarray) -> tuple[dict, dict, dict]:
    """Evidence [n_cells, 13] per readout (NaN rows where a readout has no prediction), clipped and unclipped,
    plus each readout's own training q95."""
    import torch
    torch.set_num_threads(_THREADS)
    tc = target_cols(D.adt_names)
    tr = D.rows("train")
    ev, ev_raw, info = {}, {}, {}
    need = np.unique(np.concatenate([tr, D.rows("val"), eval_rows]))
    z_need = np.asarray(D.z[need], dtype=np.float32)
    # head (registered normaliser)
    ph = np.full((D.n, 134), np.nan, np.float32)
    ph[need] = vk.head_predict(z_need, a.head_ckpt, device="cpu", size_factor=1.0)
    qh = np.asarray(reg["evidence"]["teddy_head"]["q95_train_pred"], np.float64)
    ev["head"], ev_raw["head"] = e3.normalise(ph[:, tc], qh[tc]), e3.normalise(ph[:, tc], qh[tc], clip=False)
    info["head"] = {"q95": "registered evidence.teddy_head.q95_train_pred", "ckpt": str(a.head_ckpt),
                    "ckpt_sha256": sha256_file(a.head_ckpt)}
    # null (own normaliser)
    if not a.null_ckpt.exists():
        raise SystemExit(f"null checkpoint {a.null_ckpt} missing (run --stage null)")
    pn = np.full((D.n, 134), np.nan, np.float32)
    pn[need] = vk.head_predict(z_need, a.null_ckpt, device="cpu", size_factor=1.0)
    qn = e3.q95(pn[tr])
    ev["null"], ev_raw["null"] = e3.normalise(pn[:, tc], qn[tc]), e3.normalise(pn[:, tc], qn[tc], clip=False)
    info["null"] = {"ckpt": str(a.null_ckpt), "ckpt_sha256": sha256_file(a.null_ckpt),
                    "stand_in_for_smoke": bool(a.smoke) and a.null_ckpt.resolve() == a.head_ckpt.resolve(),
                    "q95_own_train": dict(zip(DECLARED["targets"], map(float, qn[tc])))}
    # R1
    P1 = np.load(a.out_dir / "r1" / "pred_trainval.npy")
    miss = eval_rows[np.isnan(P1[eval_rows, 0])]
    if miss.size:
        lm = np.load(a.embed_dir / "z_rna_layer_means.npy", mmap_mode="r")
        P1[miss] = r1_predict(load_r1(a.out_dir / "r1" / "r1.pt"), lm, miss)
    q1 = e3.q95(P1[tr])
    ev["R1"], ev_raw["R1"] = e3.normalise(P1, q1), e3.normalise(P1, q1, clip=False)
    r1i = json.loads((a.out_dir / "r1" / "r1_done.json").read_text())
    info["R1"] = {"q95_own_train": dict(zip(DECLARED["targets"], map(float, q1))), "fit": {k: r1i["fit"][k] for k in ("best_epoch", "best_val_mse", "epochs_run")},
                  "val_pearson_vs_measured_evidence": r1i["val_pearson_vs_measured_evidence"], "n_train": r1i["n_train"]}
    # R2
    P2, x2 = load_r2_pred(a.out_dir / "r2", D.n, ("val", "train", "test") if not a.smoke else ("val", "train"))
    trp = tr[np.isfinite(P2[tr, 0])]
    if not a.smoke and trp.size != tr.size:
        raise SystemExit(f"R2 predictions missing for {tr.size - trp.size} training cells (run r2_predict)")
    q2 = e3.q95(P2[trp])
    ev["R2"], ev_raw["R2"] = e3.normalise(P2, q2), e3.normalise(P2, q2, clip=False)
    sub = np.load(a.out_dir / "r2" / "pred_subset.npz")
    q2s = e3.q95(sub["train_pred"])
    ev["R2_subsetq95"] = e3.normalise(P2, q2s)
    r2i = json.loads((a.out_dir / "r2" / "r2_fit.json").read_text())
    info["R2"] = {"q95_own_train": dict(zip(DECLARED["targets"], map(float, q2))), "n_train_cells_in_q95": int(trp.size),
                  "q95_subset_train": dict(zip(DECLARED["targets"], map(float, q2s))),
                  "fit": {k: r2i["fit"][k] for k in ("best_epoch", "best_val_mse", "epochs_run")}, "query_norm": r2i["query_norm"],
                  "val_pearson_vs_measured_evidence": r2i["val_pearson_vs_measured_evidence"],
                  "recomputed_genemean_cos_min": x2["genemean_cos_min"]}
    return ev, ev_raw, {"info": info, "r2_extra": x2}


def gap_arrays(ev: dict, meas: np.ndarray, nk: np.ndarray, t: np.ndarray, gidx: list[int]) -> tuple[dict, np.ndarray]:
    dE = {r: np.abs(v[nk][:, gidx] - v[t][:, gidx]) for r, v in ev.items()}
    dM = np.abs(meas[nk][:, gidx] - meas[t][:, gidx])
    return dE, dM


def gr_from(dE_r: np.ndarray, dM: np.ndarray, sel: np.ndarray) -> tuple[float, np.ndarray]:
    if sel.size == 0:
        return float("nan"), np.full(dM.shape[1], np.nan)
    num = np.median(dE_r[sel], axis=0)
    den = np.median(dM[sel], axis=0)
    r = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return (float(np.mean(r)) if np.all(np.isfinite(r)) else float("nan")), r


def h3a_stats(dE: dict, dM: np.ndarray, flag: np.ndarray, ix: np.ndarray, readouts=READOUTS) -> dict:
    out = {}
    gr = {}
    for r in readouts:
        gr[r] = {}
        for nm, fm in (("flagged", flag[ix]), ("unflagged", ~flag[ix])):
            g, _ = gr_from(dE[r], dM, ix[fm])
            gr[r][nm] = g
            out[f"GR_{r}_{nm}"] = g
    for R in ("R1", "R2"):
        out[f"D_{R}"] = e3.d_statistic(gr, R)
    out["D_R2_minus_D_R1"] = out["D_R2"] - out["D_R1"]
    if {"R1", "R2", "null"} <= set(readouts):  # secondary S1 / S2 (addendum secondary_H3a)
        for R in ("R1", "R2"):
            out[f"Dlog_{R}"] = e3.d_log_statistic(gr, R)
        for R in ("R1", "R2", "null"):
            out[f"term_flagged_{R}"] = gr[R]["flagged"] - gr["head"]["flagged"]
            out[f"term_unflagged_{R}"] = gr[R]["unflagged"] - gr["head"]["unflagged"]
    return out


def summarise(point: dict, boot: dict, per_donor: dict, margins: dict) -> dict:
    res = {}
    for k_, v in point.items():
        lo, hi, nb = e3.percentile_ci(boot.get(k_, np.array([])))
        row = {"point": fnum(v), "ci95": [fnum(lo), fnum(hi)], "n_boot_defined": nb,
               "per_donor": {d: fnum(pd_[k_]) for d, pd_ in per_donor.items()}}
        if k_ in margins:
            row["margin"] = margins[k_]
            row["verdict"] = e3.verdict(v, lo, hi, {d: pd_[k_] for d, pd_ in per_donor.items()}, margins[k_])
        res[k_] = row
    return res


def run_h3a(pairs_v: dict, ev: dict, meas: np.ndarray, gidx: list[int], donor_order: list[str], n_boot: int,
            seed: int, readouts=READOUTS) -> dict:
    nk, t, flag, dn = pairs_v["nk"], pairs_v["t"], pairs_v["flag"], pairs_v["donor"]
    dE, dM = gap_arrays({r: ev[r] for r in readouts}, meas, nk, t, gidx)
    allix = np.arange(nk.size)
    point = h3a_stats(dE, dM, flag, allix, readouts)
    per_donor = {d: h3a_stats(dE, dM, flag, allix[dn == d], readouts) for d in donor_order}
    boot = e3.boot_stat(nk.size, dn, donor_order, n_boot, seed, lambda ix: h3a_stats(dE, dM, flag, ix, readouts))
    per_protein = {}
    for r in readouts:
        per_protein[r] = {}
        for nm, fm in (("flagged", flag), ("unflagged", ~flag)):
            _, pp = gr_from(dE[r], dM, allix[fm])
            per_protein[r][nm] = {p_: fnum(x) for p_, x in zip(e3.GAP_PROTEINS, pp)}
    margins = {"D_R1": 0.05, "D_R2": 0.05, "D_R2_minus_D_R1": 0.05}
    return {"stats": summarise(point, boot, per_donor, margins), "per_protein_ratio": per_protein,
            "median_abs_measured_diff": {nm: {p_: fnum(x) for p_, x in zip(e3.GAP_PROTEINS, np.median(dM[allix[fm]], axis=0) if fm.any() else [None] * 4)}
                                         for nm, fm in (("flagged", flag), ("unflagged", ~flag))}}


def h3b_methods(ev: dict, cells: np.ndarray, pan: dict, knn_dis: np.ndarray) -> dict:
    names = list(DECLARED["targets"])
    M = {}
    for r in READOUTS:
        S = e3.class_scores(ev[r][cells], names, pan)
        M[r] = (e3.argmax_calls(S), e3.top_score(S))
        if r == "head":
            calls = M[r][0]
            M["head_margin"] = (calls, e3.margin_score(S))
            M["head_entropy"] = (calls, e3.entropy_score(S))
            M["head_knn"] = (calls, 1.0 - knn_dis)
    return M


def h3b_stats(correct: dict, conf: dict, rank: np.ndarray, ix: np.ndarray) -> dict:
    out = {}
    for m in correct:
        for c in COVERAGES:
            sel = e3.select_top(conf[m][ix], rank[ix], c)
            out[f"acc_{m}@{c}"] = float(np.mean(correct[m][ix][sel])) if sel.any() else float("nan")
        out[f"acc_{m}@1.0"] = float(np.mean(correct[m][ix])) if ix.size else float("nan")
    for c in COVERAGES:
        for R in ("R1", "R2"):
            for base in ("head",) + TRUST:
                out[f"{R}-{base}@{c}"] = out[f"acc_{R}@{c}"] - out[f"acc_{base}@{c}"]
        out[f"null-head@{c}"] = out[f"acc_null@{c}"] - out[f"acc_head@{c}"]
    return out


def run_h3b(cells: np.ndarray, truth: np.ndarray, donors: np.ndarray, M: dict, rank_all: np.ndarray,
            donor_order: list[str], n_boot: int, seed: int) -> dict:
    correct = {m: (calls == truth) for m, (calls, _) in M.items()}
    conf = {m: cf for m, (_, cf) in M.items()}
    rank = rank_all[cells]
    allix = np.arange(cells.size)
    point = h3b_stats(correct, conf, rank, allix)
    per_donor = {d: h3b_stats(correct, conf, rank, allix[donors == d]) for d in donor_order}
    boot = e3.boot_stat(cells.size, donors, donor_order, n_boot, seed, lambda ix: h3b_stats(correct, conf, rank, ix))
    margins = {f"{R}-{b}@{c}": 0.01 for c in COVERAGES for R in ("R1", "R2") for b in ("head",) + TRUST}
    calls_tab = {m: {f"{k_}->{c_}": int(np.sum((truth == k_) & (calls == c_))) for k_ in ("NK", "T") for c_ in e3.LINEAGES}
                 for m, (calls, _) in M.items() if m in READOUTS}
    return {"n_cells": int(cells.size), "n_NK": int(np.sum(truth == "NK")), "n_T": int(np.sum(truth == "T")),
            "per_donor_n": {d: int(np.sum(donors == d)) for d in donor_order},
            "stats": summarise(point, boot, per_donor, margins), "call_table_full_coverage": calls_tab}


def stage_evaluate(a, reg, add) -> None:
    reg_info = check_registration(a, site4=not a.smoke)
    D = Data(a)
    amend = reg["amendment_A1"]
    seed = int(reg["seeds"]["bootstrap"])
    flag_cos = float(reg["e3"]["flag_cosine"])
    t0 = time.time()
    if a.smoke:
        splits = {"val_smoke": D.rows("val")}
        orders = {"val_smoke": sorted(set(D.donors[D.rows("val")]))}
    else:
        splits = {nm: vk.split_indices(D.split, D.sites, D.donors, reg, nm) for nm in ("test_primary", "test_secondary")}
        orders = {"test_primary": list(reg["splits"]["test_primary"]["donors"]),
                  "test_secondary": list(reg["splits"]["test_secondary"]["donors"])}
    eval_rows = np.unique(np.concatenate(list(splits.values())))
    tr = D.rows("train")
    keys = D.keys_for(np.unique(np.concatenate([tr, eval_rows])), reg)
    variants = {"all": keys["primary"], "no_gdT158": keys["primary_no_gdT158"]}
    q95m = np.asarray([add["computed_train_val"]["q95_measured_train"][p] for p in DECLARED["targets"]])
    meas = np.full((D.n, len(DECLARED["targets"])), np.nan)
    meas[eval_rows] = measured_evidence(D, eval_rows, q95m)
    ev, ev_raw, ex = readout_evidence(a, D, reg, add, eval_rows)
    say(f"evaluate: evidence for {len(ev)} readouts ({time.time() - t0:.0f}s)")
    gidx = [DECLARED["targets"].index(p_) for p_ in e3.GAP_PROTEINS]
    pan = panel_in_targets(reg)
    rank_all = va.tie_break_rank(np.arange(D.n), amend)
    ref = tr[keys["primary"][tr] != vk.UNSCORED]
    if not a.smoke and vk.index_hash(ref) != add["computed_train_val"]["knn_reference"]["sha256"]:
        raise SystemExit("kNN reference cells differ from the addendum")
    zref = np.asarray(D.z[ref], dtype=np.float32)
    R: dict = {"experiment": "E3", "script_version": SCRIPT_VERSION, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
               "smoke": bool(a.smoke), "smoke_note": a.smoke, "registration_sha256": reg_info["registration_sha256"],
               "amendment_A1_sha256": reg_info["amendment_A1_sha256"], "addendum_sha256": reg_info["addendum_sha256"],
               "registration_check": reg_info, "n_boot": a.n_boot, "bootstrap_seed": seed, "flag_cosine": flag_cos,
               "readouts": ex["info"], "pairs": {}, "H3a": {}, "H3a_sensitivity": {}, "H3b": {}}
    avail = {r: np.isfinite(ev[r][:, 0]) for r in ev}
    for sname, cells_s in splits.items():
        order = orders[sname]
        P = e3.nkt_pairs_within_donors(D.z, cells_s, D.donors, variants, 10, flag_cos)
        R["pairs"][sname] = {"n_edges_per_donor": P["n_edges"], "variants": {}}
        R["H3a"][sname], R["H3a_sensitivity"][sname], R["H3b"][sname] = {}, {}, {}
        # kNN label disagreement for the head's calls (computed once per split, on the split's NK / T cells)
        for vname, key in variants.items():
            pv = P["variants"][vname]
            ok = avail["R2"][pv["nk"]] & avail["R2"][pv["t"]] & avail["R1"][pv["nk"]] & avail["R1"][pv["t"]]
            if not a.smoke and not ok.all():
                raise SystemExit(f"{sname}/{vname}: readout predictions missing for {int((~ok).sum())} pairs")
            pv = {k_: v[ok] for k_, v in pv.items()}
            R["pairs"][sname]["variants"][vname] = {
                "n_pairs": int(pv["nk"].size), "n_flagged": int(pv["flag"].sum()),
                "per_donor": {d: {"n_pairs": int(np.sum(pv["donor"] == d)), "n_flagged": int(np.sum(pv["flag"] & (pv["donor"] == d)))} for d in order},
                "n_unique_cells": int(np.unique(np.concatenate([pv["nk"], pv["t"]])).size),
                "cosine_quantiles": {q: fnum(np.quantile(pv["cos"], q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)} if pv["cos"].size else None,
                "pairs_sha256_E5_form": e3.pairs_sha256(pv["nk"], pv["t"]),
                "dropped_for_missing_predictions_smoke": int((~ok).sum())}
            R["H3a"][sname][vname] = run_h3a(pv, ev, meas, gidx, order, a.n_boot, seed)
            R["H3a_sensitivity"][sname][vname] = {
                "unclipped": run_h3a(pv, ev_raw, meas, gidx, order, a.n_boot, seed)["stats"],
                "R2_subset_q95": run_h3a(pv, {**ev, "R2": ev["R2_subsetq95"]}, meas, gidx, order, a.n_boot, seed)["stats"]}
            say(f"evaluate: {sname}/{vname} H3a: {pv['nk'].size} pairs ({int(pv['flag'].sum())} flagged); "
                f"D_R1 {R['H3a'][sname][vname]['stats']['D_R1']['point']} D_R2 {R['H3a'][sname][vname]['stats']['D_R2']['point']} "
                f"({time.time() - t0:.0f}s)")
            cells = cells_s[np.isin(key[cells_s], ["NK", "T"])]
            okc = avail["R2"][cells] & avail["R1"][cells]
            if not a.smoke and not okc.all():
                raise SystemExit(f"{sname}/{vname}: readout predictions missing for {int((~okc).sum())} NK/T cells")
            cells = cells[okc]
            S_head = e3.class_scores(ev["head"][cells], list(DECLARED["targets"]), pan)
            knn_dis = e3.knn_label_disagreement(np.asarray(D.z[cells], np.float32), zref, keys["primary"][ref],
                                                e3.argmax_calls(S_head), k=10)
            M = h3b_methods(ev, cells, pan, knn_dis)
            R["H3b"][sname][vname] = run_h3b(cells, key[cells], D.donors[cells], M, rank_all, order, a.n_boot, seed)
            R["H3b"][sname][vname]["dropped_for_missing_predictions_smoke"] = int((~okc).sum())
            st = R["H3b"][sname][vname]["stats"]
            say(f"evaluate: {sname}/{vname} H3b: {cells.size} cells; acc@0.9 head {st['acc_head@0.9']['point']} "
                f"R1 {st['acc_R1@0.9']['point']} R2 {st['acc_R2@0.9']['point']} ({time.time() - t0:.0f}s)")
    # R2 attention summary on the primary split's pair cells (descriptive)
    sname0 = next(iter(splits))
    pv0 = P["variants"]["all"] if len(splits) == 1 else e3.nkt_pairs_within_donors(D.z, splits[sname0], D.donors, {"all": variants["all"]}, 10, flag_cos)["variants"]["all"]
    pc = np.unique(np.concatenate([pv0["nk"], pv0["t"]]))
    ne = ex["r2_extra"]["n_eff"][pc]
    R["R2_attention"] = {"split": sname0, "n_pair_cells": int(pc.size), "median_effective_tokens": fnum(np.nanmedian(ne)) if np.isfinite(ne).any() else None,
                         "median_tokens": fnum(np.median(D.ntokens[pc]))}
    top1 = ex["r2_extra"]["top_token_ids"][pc, 0]
    top1 = top1[top1 >= 0]
    if top1.size and (a.ckpt / "vocab.txt").exists():
        from teddy_mm.teddy_encoder import load_vocab
        inv = {i: g for g, i in load_vocab(a.ckpt).items()}
        ids, cnt = np.unique(top1, return_counts=True)
        o = np.argsort(-cnt)[:15]
        R["R2_attention"]["most_frequent_top1_token"] = [{"gene": inv.get(int(ids[i]), str(int(ids[i]))), "n_cells": int(cnt[i])} for i in o]
    R["verdicts"] = verdicts(R, sname0)
    R["elapsed_sec"] = time.time() - t0
    write_json(a.out_dir / "E3_results.json", R)
    write_report(a.out_dir / "REPORT.md", R, sname0)
    say(f"evaluate: wrote {a.out_dir / 'E3_results.json'} and REPORT.md ({time.time() - t0:.0f}s)")


def verdicts(R: dict, sname: str) -> dict:
    h = R["H3a"][sname]["all"]["stats"]
    b = R["H3b"][sname]["all"]["stats"]
    v = {"split": sname, "variant": "all",
         "E3.H3a_R1": h["D_R1"]["verdict"], "E3.H3a_R2": h["D_R2"]["verdict"],
         "E3.H3a_pooling_R2_minus_R1": h["D_R2_minus_D_R1"]["verdict"],
         "E3.H3b": {f"{R_}-head@{c}": b[f"{R_}-head@{c}"]["verdict"] for c in COVERAGES for R_ in ("R1", "R2")},
         "point_only_D_ge_0.05": {R_: (h[f"D_{R_}"]["point"] is not None and h[f"D_{R_}"]["point"] >= 0.05) for R_ in ("R1", "R2")}}
    win = v["E3.H3a_R1"] == "win" or v["E3.H3a_R2"] == "win"
    v["falsification"] = ("not rejected: R1 or R2 has a win for E3.H3a" if win else
                          "rejected: neither R1 nor R2 reaches D_R >= 0.05 over the null with the win rules, so 'the NK-T "
                          "loss is in the readout or pooling and is repairable on frozen TEDDY' is rejected for this dataset")
    return v


def _ci(row) -> str:
    lo, hi = row["ci95"]
    return f"{row['point']} [{lo}, {hi}]"


def write_report(path: Path, R: dict, sname: str) -> None:
    v = R["verdicts"]
    L = [f"# E3: NK-T look-alike pairs and readout repair on frozen TEDDY ({R['script_version']})", ""]
    if R["smoke"]:
        L += [f"**SMOKE RUN ({R['smoke_note']}): val donor only, small subsets; no number here is a result.**", ""]
    L += [f"Registration `{R['registration_sha256'][:16]}`, amendment A1 `{R['amendment_A1_sha256'][:16]}`, addendum E3 "
          f"`{(R['addendum_sha256'] or 'none')[:16]}`. Bootstrap B = {R['n_boot']}, seed {R['bootstrap_seed']}. "
          f"Flag cosine {R['flag_cosine']}.", "",
          "E3 has no ANM arm: the decision rule is the registered Q1 rule (a declared rule; under the registered setup "
          "ANM's action readout ranks cells identically). E3 tests TEDDY readouts, not ANM.", "",
          "## Verdicts (primary split, variant all)", "",
          "| endpoint | verdict |", "|---|---|",
          f"| E3.H3a R1 (D >= 0.05) | {v['E3.H3a_R1']} |", f"| E3.H3a R2 (D >= 0.05) | {v['E3.H3a_R2']} |",
          f"| E3.H3a pooling (D_R2 - D_R1 >= 0.05) | {v['E3.H3a_pooling_R2_minus_R1']} |"]
    for k_, x in v["E3.H3b"].items():
        L.append(f"| E3.H3b {k_} (margin 0.01) | {x} |")
    L += ["", f"Falsification: {v['falsification']}.", ""]
    for s, blk in R["H3a"].items():
        for vn, h in blk.items():
            pr = R["pairs"][s]["variants"][vn]
            st = h["stats"]
            L += [f"## H3a, {s}, variant {vn}", "",
                  f"{pr['n_pairs']} NK-T pairs ({pr['n_flagged']} flagged), {pr['n_unique_cells']} cells; per donor "
                  + ", ".join(f"{d}: {x['n_pairs']} ({x['n_flagged']} flagged)" for d, x in pr["per_donor"].items()) + ".", "",
                  "| readout | GR flagged | GR unflagged |", "|---|---|---|"]
            for r in READOUTS:
                L.append(f"| {r} | {_ci(st[f'GR_{r}_flagged'])} | {_ci(st[f'GR_{r}_unflagged'])} |")
            L += ["", "| statistic | point [95% CI] | per donor | verdict |", "|---|---|---|---|"]
            for k_ in ("D_R1", "D_R2", "D_R2_minus_D_R1"):
                L.append(f"| {k_} | {_ci(st[k_])} | " + ", ".join(f"{d} {x}" for d, x in st[k_]["per_donor"].items())
                         + f" | {st[k_].get('verdict', '')} |")
            L += ["", "Secondary (addendum secondary_H3a; descriptive, no verdict): the registered D rewards a readout that "
                  "shrinks every NK-T difference uniformly, so the scale-free Dlog and the two terms of D are shown.", "",
                  "| statistic | point [95% CI] | per donor |", "|---|---|---|"]
            for k_ in ("Dlog_R1", "Dlog_R2", "term_flagged_R1", "term_unflagged_R1", "term_flagged_R2", "term_unflagged_R2",
                       "term_flagged_null", "term_unflagged_null"):
                if k_ in st:
                    L.append(f"| {k_} | {_ci(st[k_])} | " + ", ".join(f"{d} {x}" for d, x in st[k_]["per_donor"].items()) + " |")
            L.append("")
    for s, blk in R["H3b"].items():
        for vn, h in blk.items():
            st = h["stats"]
            L += [f"## H3b, {s}, variant {vn}", "", f"{h['n_cells']} key-NK / key-T cells ({h['n_NK']} NK, {h['n_T']} T).", "",
                  "| method | acc@0.9 | acc@0.8 | acc@1.0 |", "|---|---|---|---|"]
            for m in READOUTS + TRUST:
                L.append(f"| {m} | {_ci(st[f'acc_{m}@0.9'])} | {_ci(st[f'acc_{m}@0.8'])} | {st[f'acc_{m}@1.0']['point']} |")
            L += ["", "| difference | point [95% CI] | per donor | verdict |", "|---|---|---|---|"]
            for c in COVERAGES:
                for R_ in ("R1", "R2"):
                    for b in ("head",) + TRUST:
                        k_ = f"{R_}-{b}@{c}"
                        L.append(f"| {k_} | {_ci(st[k_])} | " + ", ".join(f"{d} {x}" for d, x in st[k_]["per_donor"].items())
                                 + f" | {st[k_].get('verdict', '')} |")
            L.append("")
    L += ["## Readouts", ""]
    for r, info in R["readouts"].items():
        vp = info.get("val_pearson_vs_measured_evidence")
        L.append(f"- **{r}**: " + (f"val Pearson CD56 {vp['CD56']}, CD94 {vp['CD94']}, CD335 {vp['CD335']}, CD3 {vp['CD3']}; " if vp else "")
                 + (f"fit {info['fit']}" if "fit" in info else f"ckpt {info.get('ckpt')}"))
    ra = R["R2_attention"]
    L += ["", f"R2 attention on {ra['split']} pair cells (descriptive): median effective number of tokens "
          f"{ra['median_effective_tokens']} of median {ra['median_tokens']} tokens"
          + (("; most frequent top-1 tokens " + ", ".join(f"{x['gene']} ({x['n_cells']})" for x in ra["most_frequent_top1_token"][:8]))
             if ra.get("most_frequent_top1_token") else "") + ".", ""]
    path.write_text("\n".join(L) + "\n")


# ============================================================================ main
def main(argv=None) -> int:
    global _LOG
    a = parse_args(argv)
    from threadpoolctl import threadpool_limits
    threadpool_limits(_THREADS)
    if a.stage == "register":
        if a.smoke:
            raise SystemExit("register is not a smoke stage (it is deterministic and train/val only)")
        stage_register(a)
        return 0
    a.out_dir.mkdir(parents=True, exist_ok=True)
    _LOG = a.out_dir / "progress.log"
    say(f"{SCRIPT_VERSION} stage {a.stage}{' (smoke: ' + a.smoke + ')' if a.smoke else ''}")
    reg = load_reg(a)
    add = load_addendum(a)
    stages = ["null", "r1", "r2_states", "r2_fit", "r2_predict", "evaluate"] if a.stage == "all" else [a.stage]
    for st in stages:
        if st == "null":
            stage_null(a)
        elif st == "r1":
            stage_r1(a, reg, add)
        elif st == "r2_states":
            if stage_r2_states(a, reg, add) == EXIT_INCOMPLETE:
                return EXIT_INCOMPLETE
        elif st == "r2_fit":
            stage_r2_fit(a, reg, add)
        elif st == "r2_predict":
            if stage_r2_predict(a, reg, add) == EXIT_INCOMPLETE:
                return EXIT_INCOMPLETE
        elif st == "evaluate":
            stage_evaluate(a, reg, add)
    return 0


if __name__ == "__main__":
    sys.exit(main())
