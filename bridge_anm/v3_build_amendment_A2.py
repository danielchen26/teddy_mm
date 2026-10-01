#!/usr/bin/env python3
"""Build amendment A2 to the v3 registration (registration/amendment_A2.json) from training rows only.

A2 follows A1 and is written before any site4 evaluation of E2 and E3. registration_v3.json, A1 and
the committed addenda are not changed; A2 overrides the experiment fields in ``declared()`` and adds:

* A2.1 (E3.H3a): two extra conditions for a readout's win (recovery on flagged pairs and the scale-free
  Dlog), because the registered D rewards a uniform shrinkage of every NK-T difference;
* A2.2 (E3, compute): R2's training q95 is estimated on a registered, seeded, label-free uniform sample
  of 10,000 split=train cells, drawn here before any R2 prediction (cell list in the core);
* A2.3 (E2): the registration's common decision rule is binding for falsification criteria (i), (ii).

The only numbers computed here are the A2.2 sample (split column only) and a sampling check of the
10,000-cell q95 against the q95 over all 67,405 training cells, on two readouts that exist before R2 is
trained: the registered phase-1 head and a ridge proxy of R2 (linear readout of the standardised gene-mean
z, which is what R2 computes at its initial uniform attention) fitted to the E3 measured-evidence targets.
Protein and embedding values are read on split=train rows only. ``--leakage-check`` rebuilds the core
from copies whose site4 rows (protein, cell type, embedding) are random: the core must be byte-identical;
positive control: poisoning the training rows must change it.

Usage (from the repo root):
  python bridge_anm/v3_build_amendment_A2.py                 # core + leakage check + amendment_A2.json(.sha256) + .md
  python bridge_anm/v3_build_amendment_A2.py --core-only --out-dir D
  python bridge_anm/v3_build_amendment_A2.py --record-hashes # A2 lines in registration/HASHES.txt
"""
from __future__ import annotations

import os
import sys

_THREADS = "4"
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, _THREADS)

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from lib import v3_amend as va_  # noqa: E402
from lib import v3_key as vk  # noqa: E402

BUILDER_VERSION = "v3_build_amendment_A2 1.0"
R2_SAMPLE_SEED = 31            # declared; not used by any other v3 seed (registration seeds 0-23, A1 29)
R2_SAMPLE_N = 10_000
ALT_DRAWS = {"n": 200, "seed_base": 310_000}   # design-based sampling SD: other seeded 10,000-cell draws
RIDGE_PROXY = {"alpha": 1.0, "input": "official z of split=train cells, standardised per dimension (train mean / SD)"}
GAP_PROTEINS = ("CD56", "CD94", "CD335", "CD3")
H3A_MARGIN = 0.05
OFFICIAL_SEC_PER_CELL = 10615.457 / 90261    # z_rna_manifest.json: official run, MPS fp16, batch 32
POISON_SEEDS = {"test": 101, "train": 202, "val": 303}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def r6(x):
    if isinstance(x, (float, np.floating)):
        return float(round(float(x), 6))
    if isinstance(x, dict):
        return {k: r6(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r6(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def content_bytes(d: dict) -> bytes:
    return (json.dumps(d, indent=1, sort_keys=False, ensure_ascii=False) + "\n").encode()


# ============================================================================ inputs
class Inputs:
    """The arrays the core reads. Only ``split`` is read on every row; adt and z are read on split=train rows."""

    def __init__(self, split, adt, adt_names, cell_types, z):
        self.split, self.adt, self.adt_names, self.cell_types, self.z = split, adt, list(adt_names), cell_types, z

    @classmethod
    def load(cls, processed: Path, z_path: Path) -> "Inputs":
        npz = np.load(processed / "cite_arrays.npz", allow_pickle=False)
        z = np.load(z_path, mmap_mode="r")
        return cls(npz["split"].astype(str), np.asarray(npz["adt"], dtype=np.float32),
                   [str(x) for x in npz["adt_names"]], npz["cell_types"].astype(str), np.asarray(z, dtype=np.float32))

    def poisoned(self, which: str, seed: int) -> "Inputs":
        """Copy whose rows of one split carry random protein, cell types and embedding."""
        rng = np.random.default_rng(seed)
        m = self.split == which
        adt, ct, z = self.adt.copy(), self.cell_types.copy(), self.z.copy()
        adt[m] = rng.gamma(1.0, 1.0, size=(int(m.sum()), adt.shape[1])).astype(adt.dtype)
        ct[m] = rng.choice(np.unique(self.cell_types), size=int(m.sum()))
        z[m] = rng.normal(size=(int(m.sum()), z.shape[1])).astype(np.float32)
        return Inputs(self.split, adt, self.adt_names, ct, z)


# ============================================================================ computed core
def ridge_proxy_predict(Ztr: np.ndarray, Ytr: np.ndarray, alpha: float) -> np.ndarray:
    """Ridge (closed form, intercept by centring) on standardised inputs; predictions on the same rows."""
    X = np.asarray(Ztr, dtype=np.float64)
    mu, sd = X.mean(0), np.maximum(X.std(0), 1e-6)
    X = (X - mu) / sd
    Y = np.asarray(Ytr, dtype=np.float64)
    ym = Y.mean(0)
    W = np.linalg.solve(X.T @ X + alpha * np.eye(X.shape[1]), X.T @ (Y - ym))
    return X @ W + ym


def q95_check(P_train: np.ndarray, tr: np.ndarray, sample: np.ndarray, alt: list[np.ndarray], targets: list[str],
              r2_subset: np.ndarray) -> dict:
    """q95 over the A2.2 sample vs over all training cells, with the design-based SD of a 10,000-cell q95
    (SD over other seeded draws of the same size from the same rows)."""
    pos = np.searchsorted(tr, sample)
    q_all = np.percentile(P_train, 95, axis=0)
    q_s = np.percentile(P_train[pos], 95, axis=0)
    q_alt = np.stack([np.percentile(P_train[np.searchsorted(tr, s)], 95, axis=0) for s in alt])
    sd = (q_alt - q_all[None, :]).std(axis=0, ddof=1)
    q_4k = np.percentile(P_train[np.searchsorted(tr, r2_subset)], 95, axis=0)
    per = {}
    for j, t in enumerate(targets):
        per[t] = {"q95_all_train": q_all[j], "q95_A2_sample": q_s[j],
                  "rel_diff_A2_sample": (q_s[j] - q_all[j]) / q_all[j],
                  "design_sd_10k": sd[j], "design_rel_sd_10k": sd[j] / q_all[j],
                  "z_A2_sample": (q_s[j] - q_all[j]) / sd[j] if sd[j] > 0 else 0.0,
                  "rel_diff_R2_4000_subset": (q_4k[j] - q_all[j]) / q_all[j]}
    gap = [t for t in targets if t in GAP_PROTEINS]
    return {"per_target": per,
            "max_abs_rel_diff_A2_sample": max(abs(per[t]["rel_diff_A2_sample"]) for t in targets),
            "max_abs_rel_diff_A2_sample_gap_proteins": max(abs(per[t]["rel_diff_A2_sample"]) for t in gap),
            "max_design_rel_sd_10k_gap_proteins": max(per[t]["design_rel_sd_10k"] for t in gap),
            "max_abs_z_A2_sample": max(abs(per[t]["z_A2_sample"]) for t in targets),
            "max_abs_rel_diff_R2_4000_subset_gap_proteins": max(abs(per[t]["rel_diff_R2_4000_subset"]) for t in gap)}


def build_core(X: Inputs, reg: dict, e3add: dict, head_ckpt: Path) -> dict:
    t0 = time.time()
    split = X.split
    tr = np.where(split == "train")[0].astype(np.int64)   # the only rows whose protein / embedding are read
    spec = reg["splits"]["train"]
    if tr.size != int(spec["n"]) or vk.index_hash(tr) != spec["sha256_indices"]:
        raise SystemExit("split=train rows differ from the registration")
    sample = va_.a2_draw_r2_q95_sample(split, R2_SAMPLE_SEED, R2_SAMPLE_N, "train")
    alt = [va_.a2_draw_r2_q95_sample(split, ALT_DRAWS["seed_base"] + b, R2_SAMPLE_N, "train") for b in range(ALT_DRAWS["n"])]
    targets = list(e3add["declared"]["targets"])
    j = {n: i for i, n in enumerate(X.adt_names)}
    tc = [j[p] for p in targets]
    r2_subset = np.asarray(e3add["computed_train_val"]["r2_cells"]["train"], dtype=np.int64)
    if vk.index_hash(r2_subset) != e3add["computed_train_val"]["r2_cells"]["train_sha256"]:
        raise SystemExit("E3 addendum R2 training cells differ from their hash")
    adt_tr = X.adt[tr][:, tc].astype(np.float64)
    z_tr = X.z[tr]
    # readout 1: the registered phase-1 head (its predictions of the 13 E3 targets)
    head = vk.head_predict(z_tr, head_ckpt, device="cpu", size_factor=1.0)[:, tc]
    # readout 2: ridge proxy of R2 (linear readout of the standardised gene-mean) on the measured-evidence targets
    q95m = np.percentile(adt_tr, 95, axis=0)
    Y = np.clip(adt_tr / np.maximum(q95m, 1e-12)[None, :], 0.0, 1.0)
    ridge = ridge_proxy_predict(z_tr, Y, RIDGE_PROXY["alpha"])
    log(f"core: proxies on {tr.size} training cells ({time.time() - t0:.0f}s)")
    checks = {"head": q95_check(head, tr, sample, alt, targets, r2_subset),
              "ridge_proxy_of_R2": q95_check(ridge, tr, sample, alt, targets, r2_subset)}
    n_test = int(np.sum(split == "test"))
    n_val = int(np.sum(split == "val"))
    core = {
        "amendment_id": va_.A2_ID,
        "builder": BUILDER_VERSION,
        "e3_r2_normaliser_sample": {
            "method": "numpy default_rng(seed).choice(split == 'train' rows in index order, n_cells, replace=False), sorted",
            "seed": R2_SAMPLE_SEED, "n_cells": R2_SAMPLE_N, "population_split": "train",
            "n_population": int(tr.size), "population_sha256": vk.index_hash(tr),
            "sample_sha256": vk.index_hash(sample), "share_of_population": R2_SAMPLE_N / tr.size,
            "cells": [int(c) for c in sample],
        },
        "q95_sampling_check": {
            "what": "q95 (np.percentile 95) over the A2.2 sample vs over all split=train cells, per E3 target, for two "
                    "readouts that exist before R2 is trained; design_sd_10k = SD of (q95 of a 10,000-cell draw - q95 "
                    "of all training cells) over other seeded draws (the A2.2 design); z = (A2 sample - all) / that SD. "
                    "rel_diff_R2_4000_subset = the registered 4,000-cell R2-subset normaliser (stratified 1/3 NK, 1/3 T, "
                    "1/3 other) against all training cells, for context",
            "alternative_draws": {"n": ALT_DRAWS["n"], "seeds": f"{ALT_DRAWS['seed_base']} .. {ALT_DRAWS['seed_base'] + ALT_DRAWS['n'] - 1}"},
            "ridge_proxy": RIDGE_PROXY | {"targets": "measured evidence (stored ADT / its training q95, clipped to [0, 1]), "
                                                     "the 13 E3 targets; fitted and predicted on all training cells"},
            "readouts": checks,
            "implied_bound_on_H3a_D_R2": {
                "statement": "GR_R2 on a pair set is a mean over the 4 gap proteins of a median difference divided by R2's "
                             "q95 for that protein (before clipping), so a relative error d_p of that q95 moves the protein's "
                             "ratio r_p by about -d_p r_p and D_R2 by about -mean_p d_p (r_p(flagged) - r_p(unflagged)); "
                             "|shift of D_R2| <= max_p |d_p| * max_p |r_p(flagged) - r_p(unflagged)|. Values below are "
                             "max_p |d_p| over the gap proteins, i.e. the shift per unit of max_p |r_p(flagged) - "
                             "r_p(unflagged)|; the ridge proxy is the readout closest to R2 (R2 starts as a linear readout "
                             "of the gene-mean), the head is shown for a heavier-tailed prediction distribution",
                "per_readout": {r: {"observed_A2_sample": c["max_abs_rel_diff_A2_sample_gap_proteins"],
                                    "two_design_sd": 2.0 * c["max_design_rel_sd_10k_gap_proteins"]}
                                for r, c in checks.items()},
                "registered_margin_H3a": H3A_MARGIN,
            },
        },
        "compute": {"sec_per_cell_official": OFFICIAL_SEC_PER_CELL,
                    "train_forward_hours_registered_addendum": tr.size * OFFICIAL_SEC_PER_CELL / 3600,
                    "train_forward_hours_A2": R2_SAMPLE_N * OFFICIAL_SEC_PER_CELL / 3600},
        "rows_read": {"split_column": int(split.size), "protein_and_embedding_rows": int(tr.size),
                      "val_rows_values_read": 0, "site4_rows_values_read": 0, "n_val_rows": n_val, "n_site4_rows": n_test},
    }
    return r6(core)


# ============================================================================ leakage check
def leakage_check(X: Inputs, reg: dict, e3add: dict, head_ckpt: Path, real_bytes: bytes) -> dict:
    out = {"pattern": "as A1 (v3_leakage_check.py) and the E2 / E3 addenda: rebuild the core from copies of the inputs "
                      "whose rows of one split carry random protein, cell types and embedding, and compare bytes",
           "real_core_sha256": sha256_bytes(real_bytes)}
    for which, role in (("test", "site4 rows poisoned: core must be byte-identical"),
                        ("val", "val rows poisoned: A2 reads no val value, so identical as well (information)"),
                        ("train", "positive control, training rows poisoned: core must change")):
        b = content_bytes(build_core(X.poisoned(which, POISON_SEEDS[which]), reg, e3add, head_ckpt))
        out[f"poison_{which}"] = {"seed": POISON_SEEDS[which], "role": role, "core_sha256": sha256_bytes(b),
                                  "identical_to_real": b == real_bytes}
        log(f"leakage: {which} rows poisoned -> identical {b == real_bytes}")
    out["passed"] = bool(out["poison_test"]["identical_to_real"] and not out["poison_train"]["identical_to_real"])
    return out


# ============================================================================ declared content
REG_E3_H3A = ("D_R = [GR_R(flagged) - GR_head(flagged)] - [GR_R(unflagged) - GR_head(unflagged)] for R in {R1, R2}, "
              "minus the same quantity for the null; pre-registered: >= 0.05 with win rules (pairs resampled within donors)")
REG_E3_GAP = ("per readout and protein p in {CD56, CD94, CD335, CD3}: median |predicted evidence difference| / median "
              "|measured evidence difference| within pairs; GR = mean over the 4 proteins; each readout normalised by "
              "the training q95 of its own predictions")
WIN_B = ("(b) recovery on flagged pairs: term_flagged_R = GR_R(flagged) - GR_head(flagged) > 0, its two-stage 95% "
         "interval excluding 0 (lower bound > 0) and > 0 in each primary donor")
WIN_C = ("(c) the scale-free contrast Dlog_R = [ln GR_R(fl) - ln GR_head(fl)] - [ln GR_R(unfl) - ln GR_head(unfl)] - "
         "(the same for the null), as the E3 builder computes it (addendum secondary_H3a S1), > 0 with its two-stage 95% "
         "interval excluding 0 (lower bound > 0) and > 0 in each primary donor")


def declared(core: dict, reg: dict, a1: dict, e2add_sha: str, e3add_sha: str) -> dict:
    s = core["e3_r2_normaliser_sample"]
    qc = core["q95_sampling_check"]
    E2_A1 = a1["experiments_overrides"]["E2/falsification"]
    findings = [
        {"id": "A2.1", "area": "E3.H3a (R1 and R2): scale artefact of the registered D",
         "found_by": "the E3 builder, before any site4 evaluation (registration/addenda/E3.json declared.secondary_H3a)",
         "finding": "the registered D_R = [GR_R - GR_head](flagged) - [GR_R - GR_head](unflagged) - D_null is not scale-free: "
                    "a readout that shrinks every NK-T evidence difference by one common factor s < 1 gets D_R = (1 - s) * "
                    "(GR_head(unflagged) - GR_head(flagged)) > 0 whenever the head already keeps less of the gap on flagged "
                    "pairs, and the seed-1 null does not remove this; such a readout repairs nothing on flagged pairs yet can "
                    "pass the registered win rule",
         "fix": "the registered D_R computation, margin 0.05 and win / loss / equivalent rules are kept. A readout's H3a win "
                "now additionally requires, under the same registered two-stage donor bootstrap (B = 2000, seed 1, the same "
                f"replicates as D_R) and per-donor rule: {WIN_B}; {WIN_C}. An undefined value fails. If D_R is a win but (b) "
                "or (c) fails, the H3a verdict for that readout is 'not supported (scale artefact)'; a loss, equivalent or "
                "inconclusive D_R verdict is unchanged. Falsification: the claim is rejected when neither R1 nor R2 has an "
                "A2 win. H3b is unchanged; H3a_pooling (D_R2 - D_R1) is unchanged",
         "effect": "stricter: every A2 win is a registered win (v3_amend.a2_h3a_verdict)"},
        {"id": "A2.2", "area": "E3 R2 normaliser (compute)",
         "found_by": "the E3 builder: r2_predict over all 67,405 split=train cells is about "
                     f"{core['compute']['train_forward_hours_registered_addendum']:.1f} h of TEDDY forward at the official "
                     "rate, spent on one percentile per target",
         "finding": "the E3 addendum (declared.readout_evidence) normalises R2 by the 95th percentile of its predictions over "
                    "all 67,405 split=train cells; a seeded uniform random sample of training cells estimates the same "
                    "training-population quantity with a sampling error that is small against the H3a margin "
                    "(computed.q95_sampling_check)",
         "fix": "R2's normaliser = np.percentile(95) of R2's predictions over the A2.2 sample: "
                f"{s['n_cells']:,} split=train cells, numpy default_rng({s['seed']}).choice(split == 'train' rows in index "
                f"order, {s['n_cells']:,}, replace=False), sorted; label-free (the split column only), drawn here before any "
                f"R2 prediction; cell list and sha256 ({s['sample_sha256'][:16]}...) in computed.e3_r2_normaliser_sample. "
                "r2_predict predicts these training cells (shards trainA2_NNNN) instead of every training cell; val and "
                "site4 predictions are unchanged. The 4,000-cell R2-subset normaliser stays the registered sensitivity; the "
                "unclipped sensitivity is unchanged",
         "effect": f"same estimand; training-cell forward {core['compute']['train_forward_hours_A2']:.2f} h instead of "
                   f"{core['compute']['train_forward_hours_registered_addendum']:.1f} h. Sampling check (train rows only, "
                   "computed.q95_sampling_check): over the 4 gap proteins the A2 sample's q95 differs from the all-training "
                   f"q95 by at most {max(c['max_abs_rel_diff_A2_sample_gap_proteins'] for c in qc['readouts'].values()):.4f} "
                   "(relative) on the head and the ridge proxy of R2, the design SD of a 10,000-cell q95 is at most "
                   f"{max(c['max_design_rel_sd_10k_gap_proteins'] for c in qc['readouts'].values()):.4f} (relative), and the "
                   f"largest |z| over all 13 targets is {max(c['max_abs_z_A2_sample'] for c in qc['readouts'].values()):.2f} "
                   "(the A2 sample's q95 is within sampling error of the all-training q95). Implied shift of D_R2 per unit of "
                   "max_p |r_p(flagged) - r_p(unflagged)|: "
                   f"{qc['implied_bound_on_H3a_D_R2']['per_readout']['ridge_proxy_of_R2']['observed_A2_sample']:.4f} observed "
                   f"(ridge proxy of R2; head {qc['implied_bound_on_H3a_D_R2']['per_readout']['head']['observed_A2_sample']:.4f}), "
                   f"against the margin {H3A_MARGIN}. The registered 4,000-cell subset is stratified, and its q95 differs from "
                   "the all-training q95 by up to "
                   f"{max(c['max_abs_rel_diff_R2_4000_subset_gap_proteins'] for c in qc['readouts'].values()):.2f} (relative) on "
                   "the gap proteins, which is why it is only a sensitivity"},
        {"id": "A2.3", "area": "E2 falsification criteria (i) and (ii): common decision rule",
         "found_by": "review of the E2 builder before any site4 evaluation (registration/addenda/E2.json "
                     "open_choices_fixed.verdict and interval_for_i report the per-donor clause beside the verdict, as "
                     "secondary)",
         "finding": "the registration's decision rules apply to all experiments (win = point difference >= margin, bootstrap "
                    "lower bound > 0, and the difference has the same sign and is >= margin / 2 in each primary donor), but "
                    "the E2 builder computed the E2 verdict without the per-donor clause and reported that clause as "
                    "secondary",
         "fix": "the common decision rule is binding for E2's criteria (i) and (ii). (i): a qualifying pair (A1.7: "
                "decision-accuracy losses within 0.02, >= 30 lost test_primary e2_subset cells in each unit, on the point "
                "estimates) and label counts only if |share difference| >= 0.15, its two-stage 95% interval excludes 0 on "
                "the side of the difference, and in each primary donor (13272, 19593) the donor's share difference has the "
                "same sign and |difference| >= 0.075 (an undefined donor share fails). (ii): holds only if dAUROC >= 0.02, "
                "its lower bound > 0 and dAUROC >= 0.01 in each primary donor. Verdict: 'adds information beyond the curve' "
                "if (i) or (ii) holds under this rule, otherwise 'adds nothing beyond the curve (falsified)'. The A1-only "
                "reading without the donor clause is reported beside it and never decides",
         "effect": "stricter: (i) and (ii) under A2 imply (i) and (ii) under A1 (v3_amend.a2_e2_criterion_i / _ii)"},
    ]
    overrides = {
        "E3/endpoints/E3.H3a":
            REG_E3_H3A + f". A2.1: a readout's H3a win additionally requires, under the same two-stage donor bootstrap and "
            f"per-donor rule, {WIN_B} and {WIN_C}; a D_R win without (b) and (c) is 'not supported (scale artefact)'",
        "E3/gap_ratio":
            REG_E3_GAP + f" (A2.2: R2's training q95 is estimated on the registered uniform sample of {s['n_cells']:,} "
            f"split=train cells, numpy default_rng({s['seed']}), cell list in amendment_A2.json; the 4,000-cell R2-subset "
            "normaliser is a sensitivity)",
        "E3/falsification":
            "if neither R1 nor R2 has an A2 win for E3.H3a (D_R >= 0.05 over the null with the win rules, plus (b) recovery "
            "on flagged pairs and (c) a positive scale-free Dlog_R, A2.1), 'the NK-T loss is in the readout or pooling and "
            "is repairable on frozen TEDDY' is rejected for this dataset",
        "E2/falsification":
            E2_A1 + "; A2.3: the registration's common decision rule is binding for (i) and (ii): in (i) the share "
            "difference must also have the same sign and be >= 0.075 in each primary donor (13272, 19593), with the "
            "interval excluding 0 on the side of the difference; in (ii) dAUROC must also be >= 0.01 in each primary donor",
    }
    e3 = {
        "h3a_win": {
            "registered_D": "unchanged: D_R as registered, margin 0.05, registration win / loss / equivalent rules (E3 addendum declared.H3a)",
            "conditions": {"b": {"statistic": "term_flagged_R = GR_R(flagged) - GR_head(flagged)",
                                 "rule": "point > 0, two-stage 95% lower bound > 0, > 0 in each primary donor"},
                           "c": {"statistic": "Dlog_R (E3 addendum secondary_H3a S1)",
                                 "rule": "point > 0, two-stage 95% lower bound > 0, > 0 in each primary donor"}},
            "bootstrap": "the registered two-stage bootstrap of E3.H3a (A1.1): donors then pairs within donors, B = 2000, "
                         "seed 1, the same replicates as D_R; percentile interval over the defined replicates",
            "verdict": "win if the registered D_R verdict is win and (b) and (c) hold; 'not supported (scale artefact)' if "
                       "D_R is a win and (b) or (c) fails; otherwise the registered D_R verdict",
            "decides": "test_primary, variant all; other splits, variants and sensitivities are reported with the same rule",
            "functions": "v3_amend.a2_positive_condition, v3_amend.a2_h3a_verdict, v3_amend.a2_e3_falsification",
        },
        "r2_normaliser": {"population_split": "train", "n_population": s["n_population"],
                          "population_sha256": s["population_sha256"], "seed": s["seed"], "n_cells": s["n_cells"],
                          "method": s["method"], "sample_sha256": s["sample_sha256"],
                          "prediction_shards": "<out-dir>/r2/pred/trainA2_NNNN.npz",
                          "sensitivity": "q95 over the 4,000 R2 training cells (registered sensitivity, unchanged)",
                          "function": "v3_amend.a2_r2_q95_sample"},
    }
    e2 = {"decision_rule": {
        "source": "registration experiments.common.statistics.win (decision rules, all experiments)",
        "primary_donors": list(reg["splits"]["test_primary"]["donors"]),
        "criterion_i": "a (pair, label) test passes when |diff| >= 0.15, the two-stage 95% interval excludes 0 on the side of "
                       "diff, and every primary donor's diff has the same sign with |diff| >= 0.075; pairs qualify as A1.7",
        "criterion_ii": "dAUROC >= 0.02, lower bound > 0, and dAUROC >= 0.01 in each primary donor",
        "verdict": "'adds information beyond the curve' if (i) or (ii) holds under this rule; otherwise falsified",
        "functions": "v3_amend.a2_e2_criterion_i, v3_amend.a2_e2_criterion_ii, v3_amend.a2_e2_verdict"}}
    supersedes = [
        {"file": "registration/addenda/E3.json", "sha256": e3add_sha, "unchanged_file": True,
         "fields": ["declared.readout_evidence: R2's normaliser over all 67,405 split=train cells -> the A2.2 sample",
                    "declared.falsification: 'rejected unless E3.H3a for R1 or for R2 is a win' -> an A2 win",
                    "declared.secondary_H3a: S1 (Dlog_R) and the flagged term of S2 also enter the A2.1 win conditions"]},
        {"file": "registration/addenda/E2.json", "sha256": e2add_sha, "unchanged_file": True,
         "fields": ["spec.open_choices_fixed.verdict and interval_for_i: the per-donor clause of the common win rule is "
                    "binding (A2.3), no longer secondary; the Bonferroni interval stays secondary"]},
    ]
    change_log = [
        {"file": "registration/amendment_A2.json (+ .sha256, amendment_A2_core.json, AMENDMENT_A2.md)",
         "change": "new amendment A2", "reason": "A2.1-A2.3 above; written before any site4 evaluation of E2 and E3"},
        {"file": "bridge_anm/lib/v3_amend.py",
         "change": "load_registration_amended() applies A1 then A2, verifying the registration, A1 and A2 hashes and that A2 "
                   "names the registration and A1 on disk; amendment_A2_status() for site4 guards; A2.1 / A2.2 / A2.3 functions",
         "reason": "every v3 builder reads the amended registration through this loader"},
        {"file": "scripts/v3_e3_nkt_repair.py",
         "change": "site4 stages (r2_predict on site4, evaluate) refuse unless A2 is committed and matches; r2_predict predicts "
                   "the A2.2 sample instead of all training cells; evaluate normalises R2 by the A2.2 sample, computes (b), "
                   "(c) and the A2 verdicts on the same replicates, and states the A2 falsification; DECLARED (the addendum's "
                   "declared part) is unchanged",
         "reason": "A2.1, A2.2; earlier stages (null, r1, r2_states, r2_fit) are unchanged"},
        {"file": "scripts/e2_response_decomposition.py",
         "change": "site4 stages refuse unless A2 is committed and matches; the report computes (i) and (ii) under the common "
                   "decision rule and decides the verdict with it; SPEC (the addendum's spec) is unchanged",
         "reason": "A2.3"},
        {"file": "registration/HASHES.txt", "change": "A2 lines; v3_amend.py line updated to the A2 version (the A1 version's "
                                                       "sha256 is in amendment_A1.json provenance)",
         "reason": "hash record of the registration files"},
    ]
    return {"title": "Amendment A2 to the v3 pre-registration: E3 H3a scale artefact, E3 R2 normaliser sample, E2 donor rule",
            "written_before": "any site4 (test) evaluation of E2 and E3; the only numbers computed here use split=train rows "
                              "(protein and embedding) and the split column (leakage_check)",
            "scope": "registration_v3.json, amendment_A1.json and the addenda are unchanged; builders apply A1 then A2 "
                     "(v3_amend.load_registration_amended()). Where A2 differs from the registration as amended by A1, or "
                     "from the E2 / E3 addenda (supersedes_in_addenda), A2 wins. A2 only makes verdicts stricter (A2.1, A2.3) "
                     "or estimates the same quantity with less compute (A2.2)",
            "findings": findings, "experiments_overrides": overrides, "e3": e3, "e2": e2,
            "supersedes_in_addenda": supersedes, "change_log": change_log}


# ============================================================================ render
def render_md(full: dict) -> str:
    cap = lambda t: t[:1].upper() + t[1:]  # noqa: E731
    am = full["amends"]
    L = [f"# {full['title']}", "",
         f"**Status:** written before {full['written_before']}. File `registration/amendment_A2.json`, sha256 in "
         "`amendment_A2.json.sha256`. It follows amendment A1 (`registration/amendment_A1.json` sha256 "
         f"`{am['amendment_A1_sha256']}`) and amends `registration/registration_v3.json` sha256 "
         f"`{am['registration_sha256']}`; both are unchanged. This page is rendered from the JSON by "
         "`bridge_anm/v3_build_amendment_A2.py`.", "",
         f"**Scope.** {full['scope']}.", "", "## Findings and fixes", ""]
    for f in full["findings"]:
        L += [f"### {f['id']} {f['area']}", "", f"- **Found by.** {cap(f['found_by'])}.",
              f"- **Finding.** {cap(f['finding'])}.", f"- **Fix.** {cap(f['fix'])}.", f"- **Effect.** {cap(f['effect'])}.", ""]
    c = full["computed"]
    s = c["e3_r2_normaliser_sample"]
    qc = c["q95_sampling_check"]
    L += ["## A2.2 sample and sampling check (training rows only)", "",
          f"Sample: {s['n_cells']:,} of {s['n_population']:,} split=train cells ({s['share_of_population']:.4f}), seed "
          f"{s['seed']}, sha256 `{s['sample_sha256']}` (population sha256 `{s['population_sha256']}`, as registration "
          "splits.train). The 4,000-cell R2 subset is stratified (1/3 NK, 1/3 T, 1/3 other), so its q95 is a biased "
          "estimate of the training-population q95; it stays a sensitivity.", "",
          f"{qc['what']}.", "",
          "| readout | target | q95 all train | q95 A2 sample | rel. diff | design rel. SD (10k) | z | rel. diff of the 4,000 subset |",
          "|---|---|---:|---:|---:|---:|---:|---:|"]
    for rname, rc in qc["readouts"].items():
        for t, v in rc["per_target"].items():
            L.append(f"| {rname} | {t} | {v['q95_all_train']:.4f} | {v['q95_A2_sample']:.4f} | {v['rel_diff_A2_sample']:+.4f} | "
                     f"{v['design_rel_sd_10k']:.4f} | {v['z_A2_sample']:+.2f} | {v['rel_diff_R2_4000_subset']:+.4f} |")
    b = qc["implied_bound_on_H3a_D_R2"]
    L += ["", f"Implied effect on E3.H3a: {b['statement']}. "
          + "; ".join(f"{r}: observed {v['observed_A2_sample']:.4f}, two design SDs {v['two_design_sd']:.4f}"
                      for r, v in b["per_readout"].items())
          + f"; margin {b['registered_margin_H3a']}.", "",
          f"Compute: training-cell forward {c['compute']['train_forward_hours_A2']:.2f} h instead of "
          f"{c['compute']['train_forward_hours_registered_addendum']:.2f} h at the official rate.", ""]
    lk = full["leakage_check"]
    L += ["## Leakage check", "", f"{cap(lk['pattern'])}.", "", "| poisoned rows | role | core identical to real |", "|---|---|---|"]
    for k in ("poison_test", "poison_val", "poison_train"):
        L.append(f"| {k.split('_')[1]} (seed {lk[k]['seed']}) | {lk[k]['role']} | {lk[k]['identical_to_real']} |")
    L += ["", f"Result: **{'PASS' if lk['passed'] else 'FAIL'}**; real core sha256 `{lk['real_core_sha256']}`.", "",
          "## Experiment fields replaced (path in registration_v3.json → experiments, after A1)", ""]
    for p_, v in full["experiments_overrides"].items():
        L.append(f"- `{p_}`: {v}")
    L += ["", "## Addendum fields superseded (files unchanged)", ""]
    for sp in full["supersedes_in_addenda"]:
        L.append(f"- `{sp['file']}` (sha256 `{sp['sha256'][:16]}...`): " + "; ".join(sp["fields"]))
    L += ["", "## Change log", "", "| file | change | reason |", "|---|---|---|"]
    for ch in full["change_log"]:
        L.append(f"| `{ch['file']}` | {ch['change']} | {ch['reason']} |")
    pv = full["provenance"]
    L += ["", "## Provenance", "",
          f"Builder `{pv['builder']}` sha256 `{pv['builder_sha256']}`; `v3_amend.py` sha256 `{pv['v3_amend_sha256']}`; "
          f"computed core sha256 `{pv['core_sha256']}`; built {pv['built_utc']}.", ""]
    return "\n".join(L)


# ============================================================================ hashes
A2_HASH_FILES = ("amendment_A2.json", "amendment_A2_core.json", "AMENDMENT_A2.md",
                 "../bridge_anm/lib/v3_amend.py", "../bridge_anm/v3_build_amendment_A2.py")


def record_hashes(reg_dir: Path) -> None:
    hf = reg_dir / "HASHES.txt"
    lines = hf.read_text().splitlines()
    names = set(A2_HASH_FILES)
    keep = [ln for ln in lines if not (len(ln.split()) == 2 and ln.split()[1] in names) and not ln.startswith("# A2:")]
    if keep and keep[0].startswith("#"):
        keep[0] = ("# sha256 of the v3 registration files (recorded at commit; registration_v3.json is the frozen "
                   "registration, amendment_A1.json amends it before site4, amendment_A2.json follows A1 before any site4 "
                   "evaluation of E2 and E3)")
    keep.append("# A2: v3_amend.py extended by A2 (its A1-era sha256 0ae5d516... is recorded in amendment_A1.json provenance)")
    keep += [f"{vk.sha256_file(reg_dir / n)}  {n}" for n in A2_HASH_FILES]
    hf.write_text("\n".join(keep) + "\n")
    log(f"recorded {len(A2_HASH_FILES)} A2 lines in {hf}")


# ============================================================================ main
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--registration-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--out-dir", type=Path, default=None, help="default: --registration-dir")
    p.add_argument("--core-only", action="store_true", help="write only amendment_A2_core.json")
    p.add_argument("--skip-leakage-check", action="store_true", help="development only; the amendment needs the check")
    p.add_argument("--record-hashes", action="store_true", help="only rewrite the A2 lines of registration/HASHES.txt")
    args = p.parse_args(argv)
    rd = args.registration_dir
    if args.record_hashes:
        record_hashes(rd)
        return
    out_dir = args.out_dir or rd
    out_dir.mkdir(parents=True, exist_ok=True)
    import torch
    torch.set_num_threads(int(_THREADS))
    reg = vk.load_registration(rd / "registration_v3.json")
    a1 = va_.load_amendment(rd / "amendment_A1.json", registration_path=rd / "registration_v3.json")
    e2f, e3f = rd / "addenda" / "E2.json", rd / "addenda" / "E3.json"
    hashes = (rd / "addenda" / "HASHES.txt").read_text().splitlines()
    for f in (e2f, e3f):
        if f"{vk.sha256_file(f)}  {f.name}" not in hashes:
            raise SystemExit(f"{f} does not match addenda/HASHES.txt")
    e3add = json.loads(e3f.read_text())
    if vk.sha256_file(args.ckpt) != reg["provenance"]["inputs"]["ckpt_sha256"]:
        raise SystemExit("head checkpoint differs from the registered input")
    X = Inputs.load(args.processed, args.z)
    core = build_core(X, reg, e3add, args.ckpt)
    core_b = content_bytes(core)
    (out_dir / "amendment_A2_core.json").write_bytes(core_b)
    log(f"core sha256 {sha256_bytes(core_b)}")
    if args.core_only:
        return
    leak = (leakage_check(X, reg, e3add, args.ckpt, core_b) if not args.skip_leakage_check
            else {"passed": False, "skipped": True})
    if not leak["passed"] and not args.skip_leakage_check:
        raise SystemExit(f"leakage check failed: {leak}")
    full = {"amendment_id": va_.A2_ID,
            "amends": {"file": "registration/registration_v3.json", "registration_sha256": vk.sha256_file(rd / "registration_v3.json"),
                       "after": "registration/amendment_A1.json", "amendment_A1_sha256": vk.sha256_file(rd / "amendment_A1.json")},
            **declared(core, reg, a1, vk.sha256_file(e2f), vk.sha256_file(e3f)),
            "computed": core, "leakage_check": leak,
            "provenance": {"builder": BUILDER_VERSION, "builder_sha256": vk.sha256_file(Path(__file__)),
                           "v3_amend_sha256": vk.sha256_file(Path(va_.__file__)),
                           "v3_key_sha256": vk.sha256_file(Path(vk.__file__)),
                           "core_sha256": sha256_bytes(core_b),
                           "inputs": {"cite_arrays_sha256": vk.sha256_file(args.processed / "cite_arrays.npz"),
                                      "z_sha256": vk.sha256_file(args.z), "ckpt_sha256": vk.sha256_file(args.ckpt),
                                      "addendum_E2_sha256": vk.sha256_file(e2f), "addendum_E3_sha256": vk.sha256_file(e3f)},
                           "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                           "test_rows_read": "the split column only (leakage_check)"}}
    out = out_dir / "amendment_A2.json"
    out.write_bytes(content_bytes(full))
    h = vk.sha256_file(out)
    (out_dir / "amendment_A2.json.sha256").write_text(f"{h}  amendment_A2.json\n")
    (out_dir / "AMENDMENT_A2.md").write_text(render_md(full))
    log(f"wrote {out} sha256 {h}")


if __name__ == "__main__":
    main()
