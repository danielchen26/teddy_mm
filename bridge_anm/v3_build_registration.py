#!/usr/bin/env python3
"""Build the v3 pre-registration (registration/registration_v3.json) from train + val only.

Every number this script writes -- gate thresholds, the evidence normaliser, the
panels, the question bars, the classifier settings, the flag threshold for the
NK-T pairs -- is computed from cells with split in {train, val}. The site4 (test)
rows are dropped right after loading; only their split / site / donor metadata is
kept, to define the registered evaluation splits and to draw the label-free E2/E5
subsets. ``bridge_anm/v3_leakage_check.py`` proves this by re-running the build on
inputs whose site4 protein, annotation and embedding rows are scrambled and checking
that the registration content is byte-identical.

Declared choices (made before any number below was computed) are the constants in
this file; each carries its reason in REGISTRATION_v3.md.

Usage (writes registration_v3.json + registration_v3.json.sha256 + the dev log):
  ANM_ROOT=... python bridge_anm/v3_build_registration.py \
      --processed data/processed/cite --z data/processed/cite_official/z_rna.npy \
      --ckpt outputs/outputs/cite_phase1_official/best.pt --out-dir registration
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
import itertools  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from lib import v3_key as vk  # noqa: E402

BUILDER_VERSION = "v3_build_registration 1.0"

# ============================================================================ declared choices
CLASS_MAP: dict[str, str] = {
    # B: mature / transitional B cells
    "Naive CD20+ B IGKC+": "B", "Naive CD20+ B IGKC-": "B", "Transitional B": "B",
    "B1 B IGKC+": "B", "B1 B IGKC-": "B",
    # T: CD4 / CD8 T, Treg, MAIT, gamma-delta T, double-negative T, cycling T
    "CD4+ T activated": "T", "CD4+ T activated integrinB7+": "T", "CD4+ T naive": "T",
    "CD4+ T CD314+ CD45RA+": "T",
    "CD8+ T naive": "T", "CD8+ T naive CD127+ CD26- CD101-": "T", "CD8+ T CD49f+": "T",
    "CD8+ T CD57+ CD45RA+": "T", "CD8+ T CD57+ CD45RO+": "T", "CD8+ T CD69+ CD45RA+": "T",
    "CD8+ T CD69+ CD45RO+": "T", "CD8+ T TIGIT+ CD45RA+": "T", "CD8+ T TIGIT+ CD45RO+": "T",
    "T reg": "T", "MAIT": "T", "gdT TCRVD2+": "T", "gdT CD158b+": "T", "dnT": "T", "T prog cycling": "T",
    # NK
    "NK": "NK", "NK CD158e1+": "NK",
    # myeloid: monocytes and conventional dendritic cells
    "CD14+ Mono": "myeloid", "CD16+ Mono": "myeloid", "cDC2": "myeloid", "cDC1": "myeloid",
    # OUT: no B / T / NK / myeloid call is correct
    "HSC": "OUT", "Lymph prog": "OUT", "G/M prog": "OUT", "MK/E prog": "OUT",
    "Proerythroblast": "OUT", "Erythroblast": "OUT", "Normoblast": "OUT", "Reticulocyte": "OUT",
    "pDC": "OUT", "Plasma cell IGKC+": "OUT", "Plasma cell IGKC-": "OUT",
    "Plasmablast IGKC+": "OUT", "Plasmablast IGKC-": "OUT", "ILC": "OUT", "ILC1": "OUT",
}
CLASS_MAP_NOTES: dict[str, str] = {
    "gdT CD158b+": "T by lineage (gamma-delta TCR, CD3 protein). They carry NK receptors (CD158b, CD56, CD94); "
                   "the protein gate decides whether each one enters the key as T (CD3 high) or is unscored. "
                   "Sensitivity key 'primary_no_gdT158' drops them.",
    "gdT TCRVD2+": "T (Vdelta2 gamma-delta T; CD3 protein high).",
    "dnT": "T (CD3+ CD4- CD8- T cells).",
    "T prog cycling": "T: annotated as a T-committed cycling population; 24 training cells, none in val or site4.",
    "CD4+ T CD314+ CD45RA+": "T; 92 training cells, 1 val, none in site4.",
    "CD8+ T naive CD127+ CD26- CD101-": "T; 42 training cells, none in val or site4.",
    "cDC1": "myeloid (conventional DC, like cDC2); 18 training cells, none in val or site4.",
    "cDC2": "myeloid (conventional DC; CD11c and CD33 protein high).",
    "pDC": "OUT: plasmacytoid DCs are not monocytes or conventional DCs and carry none of the four lineage "
           "panels (CD123/CD303 high, CD11c low); a lineage call on them is wrong.",
    "ILC": "OUT: innate lymphoid cells lack a TCR and are not conventional NK cells. Many are NK-like in "
           "protein (CD56/CD94 high); they then gate NK and are unscored, not keyed NK.",
    "ILC1": "OUT: as ILC. In this dataset most annotated ILC1 are CD3-protein high (T-like); those gate T and "
            "are unscored.",
    "Plasma cell IGKC+": "OUT: antibody-secreting cells are a terminal B-lineage state outside the mature-B "
                         "panel question (CD19 partly, CD20 low, CD38 very high).",
    "Plasma cell IGKC-": "OUT, as plasma cells.", "Plasmablast IGKC+": "OUT, as plasma cells.",
    "Plasmablast IGKC-": "OUT, as plasma cells.",
    "Lymph prog": "OUT: progenitor (B-committed in part; many are CD19 protein high and then unscored).",
    "HSC": "OUT: stem cells; no lineage panel applies.", "G/M prog": "OUT: progenitor.",
    "MK/E prog": "OUT: progenitor.",
    "Proerythroblast": "OUT: erythroid lineage.", "Erythroblast": "OUT: erythroid lineage.",
    "Normoblast": "OUT: erythroid lineage.", "Reticulocyte": "OUT: erythroid lineage.",
}

# Gate proteins (all measured in this panel; CD34, CD235a, CD138 and TCRgd are absent).
GATE_PROTEINS = ["CD3", "CD19", "CD20", "CD56", "CD94", "CD335", "CD16", "CD14", "CD33", "CD11c", "CD71"]
GATE_RULE_TEXT = [
    "OUT (erythroid) if CD71 high",
    "T if CD3 high and no B marker high and CD14 not high",
    "B if CD3 not high and a B marker high and CD14 not high",
    "NK if CD3, B markers, CD14 and CD33 not high and at least nk_min_high of the NK markers high",
    "myeloid if CD3 and B markers not high and any of CD14 / CD11c / CD33 high",
    "OUT otherwise (no lineage gate passed)",
]
NK_MIN_HIGH = 2
THRESHOLD_SEED = 0
# Pre-declared variant grid for the gate, scored on the 8 training batches (site x donor, split=train),
# confirmed on val. Selection: highest mean per-batch Cohen kappa (gate vs annotation, 5 classes);
# ties -> fewer markers.
GATE_VARIANTS = {
    "threshold_method": ["otsu_nonzero", "gmm_nonzero_post50", "gmm_background_mu_plus_3sd"],
    "b_markers": [["CD19", "CD20"], ["CD19"]],
    "nk_markers": [["CD56", "CD94", "CD335", "CD16"], ["CD56", "CD94", "CD335"]],
}
GATE_VAL_TARGET = {"kappa5_min": 0.85, "nk_precision_min": 0.90, "erythroid_out_min": 0.75}

PANEL_K = 3
PANEL_AUROC_FLOOR = 0.80
MIN_PAIR_CELLS = 10
STABILITY_N, STABILITY_SEED = 20000, 3
NK_MUST_CONSIDER = ["CD56", "CD94", "CD335", "CD16"]
CANONICAL9 = {"B": ["CD19", "CD72", "CD22"], "T": ["CD3", "CD2", "CD5"], "myeloid": ["CD16", "CD11c", "CD36"]}
CANONICAL9_VARIANT = {"B": ["CD19", "CD20", "CD22"], "T": ["CD3", "CD2", "CD5"], "myeloid": ["CD16", "CD11c", "CD14"]}

Q_TARGETS = {"Q1": 0.15, "Q2": 0.30, "Q3": 0.15}       # declared val no-call rates (all val cells)
Q3_ANCHORS = {"B": "CD19", "T": "CD3", "NK": "CD56", "myeloid": "CD14"}
Q3_CLASSICAL = ["CD14+ Mono"]
Q3_OTHER_MYELOID = ["CD16+ Mono", "cDC2", "cDC1"]
Q3_VAL_TARGET = {"cd14_gate_precision_for_classical_min": 0.90, "cd14_gate_recall_for_classical_min": 0.70}

ANM_FIELD = {"kind": "finite_graph_scalar", "steps": 4, "retention": 0.82, "diffusion": 0.16, "source_scale": 1.0}
ANM_CLOSURE = {"closure_weight": 1.0, "mean_coordinate_weight": 0.25, "direct_action_weight": 0.05}

CLF_C_GRID = [0.01, 0.1, 1.0, 10.0, 100.0]
CLF_SEED = 0

E3_K = 10
E2_SUBSET = {"seed": 7, "per_primary_donor": 500, "secondary_15078": 250}
E5_CAP_PER_DONOR = 200
SEEDS = {"global": 20260930, "bootstrap": 1, "classifier": [0, 1, 2, 3, 4], "label_cost_draws": [0, 1, 2, 3, 4],
         "e2_subset": 7, "e2_thinning": [11, 12], "e3_token_subset": 13, "e4_noise": 17, "e5_cells": 19,
         "e5_random_directions": 23}


# ============================================================================ helpers
def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def git_head(root: Path) -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True,
                             check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return out + ("+dirty" if dirty else "")
    except Exception:
        return None


def otsu_nonzero(x: np.ndarray, nbins: int = 512) -> float:
    x = np.asarray(x, dtype=np.float64)
    x = x[x > 0]
    hist, edges = np.histogram(x, bins=nbins)
    p = hist / hist.sum()
    c = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(p)
    w1 = 1 - w0
    m0 = np.cumsum(p * c) / np.maximum(w0, 1e-12)
    m1 = ((p * c).sum() - np.cumsum(p * c)) / np.maximum(w1, 1e-12)
    sb = w0 * w1 * (m0 - m1) ** 2
    return float(edges[int(np.argmax(sb[:-1])) + 1])


def gmm_background_mu3sd(x: np.ndarray, seed: int = 0) -> float:
    from sklearn.mixture import GaussianMixture

    x = np.asarray(x, dtype=np.float64)
    g = GaussianMixture(2, random_state=seed).fit(x[:, None])
    m = g.means_.ravel()
    k = int(np.argmin(m))
    return float(m[k] + 3.0 * np.sqrt(g.covariances_.ravel()[k]))


THRESHOLD_METHODS = {
    "otsu_nonzero": lambda x, seed: otsu_nonzero(x),
    "gmm_nonzero_post50": lambda x, seed: vk.gmm_nonzero_threshold(x, seed=seed),
    "gmm_background_mu_plus_3sd": lambda x, seed: gmm_background_mu3sd(x, seed=seed),
}


def rnd(x) -> float:
    """Every registered number is rounded to 6 decimals when it is defined and used rounded downstream."""
    return float(round(float(x), 6))


def r6(x):
    if isinstance(x, (float, np.floating)):
        return float(round(float(x), 6))
    if isinstance(x, dict):
        return {k: r6(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r6(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    return x


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUROC P(pos > neg) + 0.5 P(tie), via average ranks (equals sklearn roc_auc_score)."""
    from scipy.stats import rankdata

    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    n1, n0 = pos.size, neg.size
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(np.r_[pos, neg])
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================================ build
def build(args) -> tuple[dict, dict]:
    t0 = time.time()
    npz = np.load(args.processed / "cite_arrays.npz", allow_pickle=False)
    split = npz["split"].astype(str)
    sites = npz["sites"].astype(str)
    donors = npz["donors"].astype(str)
    adt_names = [str(x) for x in npz["adt_names"]]
    n_all = int(split.size)

    # ---- the only rows whose protein / annotation / embedding are read: split in {train, val}
    nt = np.where(split != "test")[0]
    adt = np.asarray(npz["adt"], dtype=np.float32)[nt]
    ct = npz["cell_types"].astype(str)[nt]
    zmm = np.load(args.z, mmap_mode="r")
    if zmm.shape[0] != n_all:
        raise SystemExit(f"z rows {zmm.shape[0]} != cells {n_all}")
    z = np.asarray(zmm[nt], dtype=np.float32)
    del zmm
    sp, si, do = split[nt], sites[nt], donors[nt]
    tr, va = sp == "train", sp == "val"
    log(f"loaded: {n_all} cells, using {nt.size} train+val rows; adt {adt.shape}; z {z.shape}")

    missing = sorted(set(ct) - set(CLASS_MAP))
    if missing:
        raise SystemExit(f"cell types without a class: {missing}")
    absent = [p for p in ("CD34", "CD235a", "CD138", "TCRgd") if p not in adt_names]

    reg: dict = {"registration_id": "teddy_mm_v3", "builder": BUILDER_VERSION}
    dev: dict = {"builder": BUILDER_VERSION}

    # ---------------------------------------------------------------- splits (metadata only)
    def split_spec(name, sel, rule, **kw):
        idx = np.where(sel)[0]
        return {"rule": rule, "n": int(idx.size), "sha256_indices": vk.index_hash(idx), **kw}

    test_donors = sorted(set(donors[split == "test"]))
    train_donors = sorted(set(donors[split == "train"]))
    prim_donors = [d for d in test_donors if d not in train_donors]
    sec_donors = [d for d in test_donors if d in train_donors]
    reg["splits"] = {
        "cell_id": "global row index in data/processed/cite/cite_arrays.npz (bridge export id cite_site4_<index>)",
        "train": split_spec("train", split == "train", "split == 'train' (sites 1-3, all donors except 18303)",
                            split=["train"], sites=None, donors=None,
                            site_donor_counts={f"{s}|{d}": int(np.sum((split == 'train') & (sites == s) & (donors == d)))
                                               for s in sorted(set(sites[split == 'train']))
                                               for d in sorted(set(donors[(split == 'train') & (sites == s)]))}),
        "val": split_spec("val", split == "val", "split == 'val' (site1, donor 18303)", split=["val"], sites=None,
                          donors=sorted(set(donors[split == "val"]))),
        "test_primary": split_spec("test_primary", (split == "test") & np.isin(donors, prim_donors),
                                   "site4 cells of donors not in training", split=["test"], sites=["site4"],
                                   donors=prim_donors,
                                   per_donor={d: int(np.sum((split == 'test') & (donors == d))) for d in prim_donors}),
        "test_secondary": split_spec("test_secondary", (split == "test") & np.isin(donors, sec_donors),
                                     "site4 cells of donor 15078 (in training at sites 1-3)", split=["test"],
                                     sites=["site4"], donors=sec_donors),
    }
    log(f"splits: primary {reg['splits']['test_primary']['n']} secondary {reg['splits']['test_secondary']['n']}")

    # ---------------------------------------------------------------- annotation map
    reg["classes"] = list(vk.CLASSES)
    reg["annotation_map"] = dict(CLASS_MAP)
    reg["annotation_map_notes"] = dict(CLASS_MAP_NOTES)
    annot = np.asarray([CLASS_MAP[t] for t in ct])
    reg["data"] = {"n_cells": n_all, "n_adt": len(adt_names), "n_cell_types_total": len(CLASS_MAP),
                   "adt_absent_checked": absent, "adt_values": "CLR-normalised (float, as stored in cite_arrays.npz)",
                   "cell_types_per_split": {}}
    for t in sorted(CLASS_MAP):
        reg["data"]["cell_types_per_split"][t] = {"train": int(np.sum((ct == t) & tr)), "val": int(np.sum((ct == t) & va))}

    # ---------------------------------------------------------------- gate: variant grid on training batches
    j = {n: i for i, n in enumerate(adt_names)}
    batch = vk.batch_labels(si, do)
    thr_by_method = {}
    for mname, f in THRESHOLD_METHODS.items():
        thr_by_method[mname] = {p: rnd(f(adt[tr, j[p]], THRESHOLD_SEED)) for p in GATE_PROTEINS}
    grid = []
    train_batches = sorted(set(batch[tr]))
    for mname, bm, nkm in itertools.product(GATE_VARIANTS["threshold_method"], GATE_VARIANTS["b_markers"],
                                            GATE_VARIANTS["nk_markers"]):
        gspec = {"b_markers": bm, "nk_markers": nkm, "nk_min_high": NK_MIN_HIGH}
        H = {p: adt[:, j[p]] > thr_by_method[mname][p] for p in GATE_PROTEINS}
        g = vk.gate_class_from_high(H, gspec)
        kb = {b: vk.cohen_kappa(annot[(batch == b) & tr], g[(batch == b) & tr]) for b in train_batches}
        vrep = vk.agreement_report(annot[va], g[va])
        grid.append({"threshold_method": mname, "b_markers": bm, "nk_markers": nkm,
                     "train_batch_kappa": kb, "mean_train_batch_kappa": float(np.mean(list(kb.values()))),
                     "min_train_batch_kappa": float(np.min(list(kb.values()))),
                     "val_kappa5": vrep["kappa5"], "n_markers": len(bm) + len(nkm)})
    best = sorted(grid, key=lambda r: (-round(r["mean_train_batch_kappa"], 6), r["n_markers"]))[0]
    gate = {
        "proteins": GATE_PROTEINS, "threshold_method": best["threshold_method"], "threshold_fit_on": "train",
        "threshold_seed": THRESHOLD_SEED, "thresholds": thr_by_method[best["threshold_method"]],
        "b_markers": best["b_markers"], "nk_markers": best["nk_markers"], "nk_min_high": NK_MIN_HIGH,
        "high_definition": "measured CLR value > threshold (strictly greater)",
        "rules_in_order": GATE_RULE_TEXT,
        "selection": {"criterion": "highest mean per-batch Cohen kappa (gate vs annotation, 5 classes) over the 8 "
                                   "training batches (site x donor); ties -> fewer markers; confirmed on val",
                      "grid": GATE_VARIANTS, "chosen": {k: best[k] for k in ("threshold_method", "b_markers", "nk_markers")},
                      "chosen_mean_train_batch_kappa": best["mean_train_batch_kappa"],
                      "chosen_min_train_batch_kappa": best["min_train_batch_kappa"],
                      "grid_results": [{k: r[k] for k in ("threshold_method", "b_markers", "nk_markers",
                                                           "mean_train_batch_kappa", "min_train_batch_kappa",
                                                           "val_kappa5")} for r in grid]},
        "val_target": GATE_VAL_TARGET,
    }
    reg["gate"] = gate
    dev["gate_grid"] = grid
    dev["thresholds_by_method"] = thr_by_method
    log(f"gate chosen: {gate['selection']['chosen']} mean train kappa {best['mean_train_batch_kappa']:.4f}")

    # ---------------------------------------------------------------- keys on train / val
    gated = vk.gate_class(adt, adt_names, reg)
    prim = vk.primary_key(annot, gated)
    reg["sensitivity_keys"] = {
        "annotation_only": "annotation class for every cell (OUT included)",
        "gate_only": "gated class for every cell",
        "primary_per_batch": "agreement of annotation with a gate whose thresholds are re-estimated within each "
                             "site x donor batch by the same declared estimator (label-free; uses the batch's own "
                             "measured protein, so it is a verifier variant, never evidence)",
        "no_gdT158": {"excluded_type": "gdT CD158b+", "rule": "primary key with gdT CD158b+ cells set to unscored"},
    }
    kq = {"train": vk.agreement_report(annot[tr], gated[tr]), "val": vk.agreement_report(annot[va], gated[va]),
          "per_batch": {b: vk.agreement_report(annot[batch == b], gated[batch == b]) for b in sorted(set(batch))}}
    ery = np.isin(ct, ["Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte"])
    kq["val"]["erythroid_gated_OUT_rate"] = float(np.mean(gated[va & ery] == "OUT")) if np.any(va & ery) else None
    kq["val"]["primary_key_counts"] = {c: int(np.sum(prim[va] == c)) for c in list(vk.CLASSES) + [vk.UNSCORED]}
    kq["train"]["primary_key_counts"] = {c: int(np.sum(prim[tr] == c)) for c in list(vk.CLASSES) + [vk.UNSCORED]}
    vt = kq["val"]
    kq["val_target_met"] = {
        "kappa5": vt["kappa5"] >= GATE_VAL_TARGET["kappa5_min"],
        "nk_precision": (vt["per_class"]["NK"]["precision_vs_annotation"] or 0) >= GATE_VAL_TARGET["nk_precision_min"],
        "erythroid_out": (vt["erythroid_gated_OUT_rate"] or 0) >= GATE_VAL_TARGET["erythroid_out_min"],
    }
    per_type = {}
    for t in sorted(set(ct)):
        m = ct == t
        per_type[t] = {"n_train_val": int(m.sum()), "annotation_class": CLASS_MAP[t],
                       "gated": {c: int(np.sum(gated[m] == c)) for c in vk.CLASSES},
                       "in_primary_key": int(np.sum(prim[m] != vk.UNSCORED))}
    kq["per_cell_type_train_val"] = per_type
    # per-batch thresholds (sensitivity key), training + val batches only
    reg["gate"]["per_batch_thresholds_train_val"] = vk.per_batch_thresholds(adt, adt_names, batch, reg)
    gb = vk.gate_class_per_batch(adt, adt_names, batch, reg)
    kq["per_batch_gate_sensitivity"] = {"train": vk.agreement_report(annot[tr], gb[tr]),
                                        "val": vk.agreement_report(annot[va], gb[va])}
    reg["key_quality"] = kq
    log(f"val kappa5 {vt['kappa5']:.4f} agreement {vt['agreement']:.4f}; targets {kq['val_target_met']}")

    # ---------------------------------------------------------------- evidence (head on train + val)
    pred = vk.head_predict(z, args.ckpt, device="cpu", size_factor=1.0)
    q95 = np.percentile(pred[tr], 95, axis=0)
    reg["evidence"] = {"teddy_head": {
        "embedding": "official TEDDY-G gene-mean z (data/processed/cite_official/z_rna.npy), L2-normalised",
        "head": "phase-1 MLP + NB decoder mean (no flow matching), outputs/cite_phase1_official/best.pt",
        "size_factor": 1.0,
        "normaliser": "95th percentile of the head's prediction over the training cells (split=train), per protein",
        "clip": [0.0, 1.0],
        "why_size_factor_cancels": "prediction = softplus(.) * s for a constant s, and the training q95 scales by "
                                   "the same s, so the evidence is independent of s; no measured-protein array "
                                   "enters the evidence path",
        "adt_names": adt_names, "q95_train_pred": [rnd(x) for x in q95]}}
    v_all = vk.evidence(pred, reg)                 # clipped
    v_raw = vk.evidence(pred, reg, clip=False)     # unclipped (selection only)

    # ---------------------------------------------------------------- panel selection on val
    vs = va & (prim != vk.UNSCORED)
    pv = prim[vs]

    def score_tables(mask):
        """worst-pair AUROC (selection score), one-vs-rest AUROC and Cohen's d (reported only)."""
        pvm = prim[mask]
        worst, ovr, cd = ({k: {} for k in vk.LINEAGES} for _ in range(3))
        for k in vk.LINEAGES:
            pos = pvm == k
            others = [c for c in vk.CLASSES if c != k and np.sum(pvm == c) >= MIN_PAIR_CELLS]
            for p in adt_names:
                x = v_raw[mask, j[p]]
                ovr[k][p] = auroc(x[pos], x[~pos])
                worst[k][p] = min(auroc(x[pos], x[pvm == c]) for c in others)
                cd[k][p] = float((x[pos].mean() - x[~pos].mean())
                                 / np.sqrt(0.5 * (x[pos].var() + x[~pos].var()) + 1e-12))
        return worst, ovr, cd

    au, au_ovr, au_d = score_tables(vs)
    best_cls = {p: max(vk.LINEAGES, key=lambda k: (au[k][p], -vk.LINEAGES.index(k))) for p in adt_names}

    def select(k_per, tab=au, bcls=best_cls, floor=PANEL_AUROC_FLOOR, tie=au_ovr):
        pan = {}
        for k in vk.LINEAGES:
            cands = [p for p in adt_names if bcls[p] == k and tab[k][p] >= floor]
            cands.sort(key=lambda p: (-tab[k][p], -tie[k][p], p))
            pan[k] = cands[:k_per]
        return pan

    panel = select(PANEL_K)
    nk_vs_t = {}
    mnt = vs & np.isin(prim, ["NK", "T"])
    for p in NK_MUST_CONSIDER + panel["NK"]:
        x = v_raw[mnt, j[p]]
        nk_vs_t[p] = auroc(x[prim[mnt] == "NK"], x[prim[mnt] == "T"])
    ranks = {k: sorted(adt_names, key=lambda p: (-au[k][p], -au_ovr[k][p], p)) for k in vk.LINEAGES}
    # stability diagnostic (not used for selection): the same rule on training cells (in-sample for the head)
    trs = tr & (prim != vk.UNSCORED)
    _rs = np.random.default_rng(STABILITY_SEED)
    _keep = np.zeros_like(trs)
    _keep[_rs.choice(np.where(trs)[0], size=min(STABILITY_N, int(trs.sum())), replace=False)] = True
    trs = _keep
    au_tr, au_ovr_tr, _ = score_tables(trs)
    best_tr = {p: max(vk.LINEAGES, key=lambda k: (au_tr[k][p], -vk.LINEAGES.index(k))) for p in adt_names}
    panel_tr = select(PANEL_K, tab=au_tr, bcls=best_tr, tie=au_ovr_tr)
    reg["panels"] = {
        "selection": {
            "split": "val (donor 18303), primary-key cells (unscored excluded)",
            "score": "worst-pair AUROC: for class k, the minimum over every other primary-key class c (OUT "
                     f"included, classes with >= {MIN_PAIR_CELLS} val cells) of the AUROC of the unclipped "
                     "normalised head prediction for k versus c; ties by one-vs-rest AUROC, then name",
            "rule": f"each protein is assigned to the class where its worst-pair AUROC is highest; a class panel is "
                    f"its top {PANEL_K} assigned proteins with worst-pair AUROC >= {PANEL_AUROC_FLOOR}",
            "candidates": "all 134 measured-panel proteins",
            "criteria_computed_and_not_used": "one-vs-rest AUROC (saturates: many proteins within 0.002 of 1 on "
                                              "val) and Cohen's d (scale- and clip-dependent); both tables are in "
                                              "the dev log",
            "n_val_cells_used": int(vs.sum()),
            "class_counts": {k: int(np.sum(pv == k)) for k in vk.CLASSES},
            "nk_must_consider": {p: {"worst_pair_auroc_nk": au["NK"][p], "auroc_nk_vs_rest": au_ovr["NK"][p],
                                     "assigned_class": best_cls[p], "worst_pair_auroc_assigned": au[best_cls[p]][p],
                                     "rank_for_nk": ranks["NK"].index(p) + 1, "auroc_nk_vs_t": nk_vs_t[p]}
                                 for p in NK_MUST_CONSIDER},
            "top10_per_class": {k: [{"protein": p, "worst_pair_auroc": au[k][p], "auroc_vs_rest": au_ovr[k][p],
                                     "cohen_d": au_d[k][p], "assigned": best_cls[p]} for p in ranks[k][:10]]
                                for k in vk.LINEAGES},
            "stability_same_rule_on_training_cells": {"panel": panel_tr,
                                                      "overlap_with_val_panel": {k: len(set(panel[k]) & set(panel_tr[k]))
                                                                                 for k in vk.LINEAGES},
                                                      "note": "diagnostic only, on a seeded sample of 20,000 "
                                                              "training primary-key cells (in-sample for the head)"},
        },
        "primary": {"proteins": panel, "weights": "equal within a class"},
        "sensitivity": {},
    }
    dev["auroc_val_worst_pair"] = au
    dev["auroc_val_one_vs_rest"] = au_ovr
    dev["cohen_d_val"] = au_d
    nk3 = panel["NK"]

    def with_nk(base):
        used = {p for v in base.values() for p in v}
        return {"B": base["B"], "T": base["T"], "NK": [p for p in nk3 if p not in used], "myeloid": base["myeloid"]}

    reg["panels"]["canonical9_plus_nk"] = {"proteins": with_nk(CANONICAL9),
                                           "note": "old canonical 9-marker panel (chosen on site4 Pearson in v1/v2) "
                                                   "plus the val-selected NK panel minus any protein already in it"}
    reg["panels"]["canonical9_variant_plus_nk"] = {"proteins": with_nk(CANONICAL9_VARIANT),
                                                   "note": "canonical panel with CD20 for CD72 and CD14 for CD36, plus "
                                                           "the val-selected NK panel"}
    reg["panels"]["val_k2"] = {"proteins": select(2), "note": "same rule, 2 per class"}
    reg["panels"]["val_k5"] = {"proteins": select(5), "note": "same rule, 5 per class"}
    reg["panels"]["sensitivity"] = ["canonical9_plus_nk", "canonical9_variant_plus_nk", "val_k2", "val_k5"]
    log(f"panel: {panel}")

    # ---------------------------------------------------------------- questions: bars on val
    reg["questions"] = {}
    names = adt_names
    for q, spec in (("Q1", {"panel": "primary"}), ("Q2", {"panel": "primary"}), ("Q3", {"panel": "anchors", "anchors": Q3_ANCHORS})):
        reg["questions"][q] = dict(spec)
    S = {q: vk.class_scores(v_all, names, reg, q) for q in ("Q1", "Q2", "Q3")}
    for q in ("Q1", "Q2", "Q3"):
        top = S[q][va].max(axis=1)
        bar = rnd(np.quantile(top, Q_TARGETS[q]))
        reg["questions"][q]["bar"] = bar
        reg["questions"][q]["val_no_call_target"] = Q_TARGETS[q]
        reg["questions"][q]["val_no_call_realised"] = float(np.mean(top < bar))
    reg["questions"]["Q1"].update({
        "name": "soft 4-class lineage question",
        "text": "From the TEDDY + head evidence, is this cell B, T, NK or myeloid? Call the class with the "
                "highest panel score if that score reaches the bar; otherwise no call.",
        "actions": list(vk.LINEAGES), "score": "equal-weight mean of the evidence over the class panel",
        "key": "primary"})
    reg["questions"]["Q2"].update({
        "name": "strict 4-class lineage question",
        "text": "As Q1 with a higher bar (nested: every Q2 call is a Q1 call with the same class).",
        "actions": list(vk.LINEAGES), "score": "as Q1", "key": "primary"})
    reg["questions"]["Q3"].update({
        "name": "CD14-anchored question (B / T / NK / classical monocyte)",
        "text": "Is this cell a B cell (CD19), T cell (CD3), NK cell (CD56) or a classical CD14+ monocyte (CD14)? "
                "Score each class by its one canonical anchor protein; call the highest if it reaches the bar; "
                "otherwise no call. Non-classical monocytes and dendritic cells are not an answer (no call).",
        "actions": list(vk.LINEAGES), "score": "evidence of the class anchor protein",
        "key": {"from": "primary", "anchor_protein": "CD14", "classical_types": Q3_CLASSICAL,
                "other_myeloid_types": Q3_OTHER_MYELOID,
                "rule": "primary-key myeloid: annotated classical monocyte AND measured CD14 high -> myeloid; "
                        "annotated other myeloid AND measured CD14 not high -> OUT (no call correct); else unscored"},
        "val_target": Q3_VAL_TARGET})
    # val outcomes of the questions (development only)
    q3k = vk.q3_key(prim, ct, adt, adt_names, reg)
    keys_val = {"Q1": prim[va], "Q2": prim[va], "Q3": q3k[va]}
    for q in ("Q1", "Q2", "Q3"):
        calls = vk.rule_calls(S[q][va], reg["questions"][q]["bar"])
        reg["questions"][q]["val_outcome_rule"] = vk.selective_accuracy(calls, keys_val[q])
    # Q3 key validation on val: does measured CD14 (gate threshold) separate annotated classical monocytes
    # from the other primary-key myeloid cells?
    m14 = va & (prim == "myeloid")
    cd14_hi = adt[:, j["CD14"]] > float(reg["gate"]["thresholds"]["CD14"])
    cls14 = np.isin(ct, Q3_CLASSICAL)
    oth = np.isin(ct, Q3_OTHER_MYELOID)
    reg["questions"]["Q3"]["val_key_check"] = {
        "n_primary_myeloid_val": int(m14.sum()), "n_q3_myeloid_val": int(np.sum(va & (q3k == "myeloid"))),
        "n_q3_out_from_myeloid_val": int(np.sum(m14 & (q3k == "OUT"))),
        "n_q3_unscored_from_myeloid_val": int(np.sum(m14 & (q3k == vk.UNSCORED))),
        "cd14_gate_precision_for_classical": float(np.sum(m14 & cd14_hi & cls14) / max(1, np.sum(m14 & cd14_hi))),
        "cd14_gate_recall_for_classical": float(np.sum(m14 & cd14_hi & cls14) / max(1, np.sum(m14 & cls14))),
        "cd14_gate_specificity_for_other_myeloid": float(np.sum(m14 & ~cd14_hi & oth) / max(1, np.sum(m14 & oth))),
    }
    vk3 = reg["questions"]["Q3"]["val_key_check"]
    reg["questions"]["Q3"]["val_target_met"] = {
        "precision": vk3["cd14_gate_precision_for_classical"] >= Q3_VAL_TARGET["cd14_gate_precision_for_classical_min"],
        "recall": vk3["cd14_gate_recall_for_classical"] >= Q3_VAL_TARGET["cd14_gate_recall_for_classical_min"]}
    anchors_auroc = {}
    for k, p in Q3_ANCHORS.items():
        anchors_auroc[p] = {"class": k, "worst_pair_auroc_val": au[k][p], "auroc_vs_rest_val": au_ovr[k][p]}
    m_my = va & np.isin(q3k, ["myeloid", "OUT"]) & (prim == "myeloid")
    if m_my.any():
        x = v_raw[m_my, j["CD14"]]
        anchors_auroc["CD14"]["auroc_classical_vs_other_myeloid_val"] = auroc(x[q3k[m_my] == "myeloid"], x[q3k[m_my] == "OUT"])
    reg["questions"]["Q3"]["anchor_auroc_val"] = anchors_auroc
    log("bars: " + ", ".join(f"{q}={reg['questions'][q]['bar']:.4f}" for q in ("Q1", "Q2", "Q3")))

    # ---------------------------------------------------------------- ANM criteria
    sizes = {k: len(v) for k, v in panel.items()}
    equal = len(set(sizes.values())) == 1
    n_ev = sizes["B"]
    G = vk.field_gain(n_ev, ANM_FIELD)
    G1 = vk.field_gain(1, ANM_FIELD)
    reg["anm"] = {
        "engine": "ANM finite_field_runner (ANM_ROOT, public v2 fix branch), kind finite_graph_scalar",
        "field_representation": ANM_FIELD,
        "event_timing": "simultaneous: every evidence event of a cell at t = 0",
        "event_value": "evidence v_p in [0, 1] (no weight: equal weights within a class)",
        "actions": list(vk.LINEAGES), "panel_sizes_equal": equal, "panel_size": n_ev,
        "field_gain": {str(n_ev): G, "1": G1},
        "readout_thresholds": {
            "Q1": G * n_ev * reg["questions"]["Q1"]["bar"], "Q2": G * n_ev * reg["questions"]["Q2"]["bar"],
            "Q3": G1 * 1 * reg["questions"]["Q3"]["bar"]},
        "readout_thresholds_note": "informational (rounded); v3_key.anm_readout_threshold recomputes field_gain(n) * n * bar "
                                   "exactly from the field parameters and the registered bar",
        "equivalence": "with equal panel sizes and simultaneous events the ANM action score is field_gain(n) * n * S "
                       "for every class, so ANM's calls equal the re-coded fixed rule's on every cell (checked cell by "
                       "cell in E1; any mismatch is a failure)",
        "closure_readout": ANM_CLOSURE,
        "closure_readout_note": "readout_coordinates over each action's event sites (ANM readout defaults); used only "
                                "for C2 (nested readouts); its re-coded rule is reported beside it",
    }
    if not equal:
        raise SystemExit(f"unequal panel sizes {sizes}: the ANM threshold mapping above assumes equal sizes")
    cl = vk.anm_action_scores(v_all[va], names, reg, "Q1", readout="closure").max(axis=1)
    reg["anm"]["closure_bar_Q1"] = rnd(np.quantile(cl, Q_TARGETS["Q1"]))
    cl2 = vk.anm_action_scores(v_all[va], names, reg, "Q2", readout="closure").max(axis=1)
    reg["anm"]["closure_bar_Q2"] = rnd(np.quantile(cl2, Q_TARGETS["Q2"]))

    # ---------------------------------------------------------------- trained classifier (C on val)
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss

    pan_idx = [j[p] for k in vk.LINEAGES for p in panel[k]]
    Xtr, Xva = v_all[tr][:, pan_idx], v_all[va][:, pan_idx]
    ytr, yva = prim[tr], prim[va]
    mtr, mva = ytr != vk.UNSCORED, yva != vk.UNSCORED
    clf_res = []
    for C in CLF_C_GRID:
        m = LogisticRegression(C=C, max_iter=3000, random_state=CLF_SEED)
        m.fit(Xtr[mtr], ytr[mtr])
        ll = log_loss(yva[mva], m.predict_proba(Xva[mva]), labels=m.classes_)
        clf_res.append({"C": C, "val_log_loss": float(ll)})
    bestC = min(clf_res, key=lambda r: (r["val_log_loss"], r["C"]))["C"]
    m = LogisticRegression(C=bestC, max_iter=3000, random_state=CLF_SEED).fit(Xtr[mtr], ytr[mtr])
    P = m.predict_proba(Xva)
    lin_cols = [list(m.classes_).index(k) for k in vk.LINEAGES]
    conf = P[:, lin_cols].max(axis=1)
    clf_bars = {q: rnd(np.quantile(conf, Q_TARGETS[q])) for q in ("Q1", "Q2")}
    # Q3 classifier (all training labels) C on val
    q3tr, q3va = q3k[tr], q3k[va]
    m3tr, m3va = q3tr != vk.UNSCORED, q3va != vk.UNSCORED
    q3_idx = [j[p] for k in vk.LINEAGES for p in panel[k]]
    for p in Q3_ANCHORS.values():
        if j[p] not in q3_idx:
            q3_idx.append(j[p])
    clf3 = []
    for C in CLF_C_GRID:
        m3 = LogisticRegression(C=C, max_iter=3000, random_state=CLF_SEED).fit(v_all[tr][:, q3_idx][m3tr], q3tr[m3tr])
        clf3.append({"C": C, "val_log_loss": float(log_loss(q3va[m3va], m3.predict_proba(v_all[va][:, q3_idx][m3va]),
                                                            labels=m3.classes_))})
    bestC3 = min(clf3, key=lambda r: (r["val_log_loss"], r["C"]))["C"]
    reg["classifier"] = {
        "primary": {"model": "multinomial logistic regression (sklearn LogisticRegression, lbfgs, max_iter 3000)",
                    "features": "the primary-panel evidence (12 values, same as the rule and ANM)",
                    "feature_proteins": [adt_names[i] for i in pan_idx],
                    "labels": "primary key on training cells (5 classes incl. OUT; unscored excluded)",
                    "n_train_labelled": int(mtr.sum()), "class_weight": None, "seed": CLF_SEED,
                    "C_grid": CLF_C_GRID, "C_selection": "lowest val log-loss on val primary-key cells",
                    "val_log_loss": clf_res, "C": bestC,
                    "decision": "call the lineage with the highest probability if that probability reaches the bar; "
                                "otherwise no call (an OUT-probable cell has low lineage probabilities)",
                    "bars": clf_bars, "bar_rule": "val quantile matching the question's declared val no-call rate"},
        "secondary": {"model": "sklearn MLPClassifier(hidden_layer_sizes=(64,), alpha=1e-4, early_stopping=True, "
                               "validation_fraction=0.1, max_iter=500, random_state=seed)",
                      "features": "as primary", "hyperparameters": "fixed a priori (no tuning)"},
        "q3": {"features": "primary-panel evidence plus the Q3 anchors not already in it",
               "feature_proteins": [adt_names[i] for i in q3_idx],
               "labels": "Q3 key on training cells", "C_grid": CLF_C_GRID, "val_log_loss": clf3, "C": bestC3},
    }
    log(f"classifier C={bestC} (Q3 C={bestC3})")

    # ---------------------------------------------------------------- E3: flag threshold from val NK-T neighbour pairs
    from lib.knn_cosine import knn_cosine

    zv = z[va]
    nn_idx, nn_sim = knn_cosine(zv, E3_K)
    pv_all = prim[va]
    pairs = {}
    for i in range(nn_idx.shape[0]):
        for jj, s in zip(nn_idx[i], nn_sim[i]):
            a_, b_ = (i, int(jj)) if i < jj else (int(jj), i)
            pairs[(a_, b_)] = float(s)
    nkt = [s for (a_, b_), s in pairs.items() if {pv_all[a_], pv_all[b_]} == {"NK", "T"}]
    reg["e3"] = {"pairs_rule": f"k = {E3_K} cosine neighbours on the raw L2-normalised final z within the evaluated "
                                f"split; each neighbour edge kept once; NK-T pair = one primary-key NK and one "
                                f"primary-key T cell",
                 "flag_rule": "sufficiency-flagged pair: cosine >= flag_cosine (TEDDY puts the two cells closer than "
                              "the median val NK-T neighbour pair)",
                 "flag_cosine": rnd(np.median(nkt)) if nkt else None, "n_val_nkt_pairs": len(nkt)}
    log(f"E3 flag cosine {reg['e3']['flag_cosine']} from {len(nkt)} val NK-T pairs")

    # ---------------------------------------------------------------- E2 / E5 subsets (metadata only)
    rng = np.random.default_rng(E2_SUBSET["seed"])
    sub = []
    for d in prim_donors:
        idx = np.where((split == "test") & (donors == d))[0]
        sub += sorted(rng.choice(idx, size=min(E2_SUBSET["per_primary_donor"], idx.size), replace=False).tolist())
    sec = []
    for d in sec_donors:
        idx = np.where((split == "test") & (donors == d))[0]
        sec += sorted(rng.choice(idx, size=min(E2_SUBSET["secondary_15078"], idx.size), replace=False).tolist())
    reg["e2_subset"] = {"rule": f"seeded (seed {E2_SUBSET['seed']}) uniform sample without labels: "
                                f"{E2_SUBSET['per_primary_donor']} cells per primary donor, "
                                f"{E2_SUBSET['secondary_15078']} from donor 15078",
                        "primary_indices": [int(i) for i in sub], "secondary_indices": [int(i) for i in sec],
                        "sha256_primary": vk.index_hash(sub), "sha256_secondary": vk.index_hash(sec)}

    reg["seeds"] = SEEDS
    dev["elapsed_sec"] = time.time() - t0
    return r6(reg), r6(dev)


def content_bytes(reg: dict) -> bytes:
    return (json.dumps(reg, indent=1, sort_keys=False, ensure_ascii=False) + "\n").encode()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--out-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--core-only", action="store_true",
                   help="write only the computed core (registration_v3_core.json); used by the leakage check")
    p.add_argument("--experiments", type=Path, default=ROOT / "registration" / "experiments_v3.json",
                   help="hand-written E1-E6 endpoints / comparators / margins merged into the registration")
    args = p.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    reg, dev = build(args)
    core = content_bytes(reg)
    (args.out_dir / "registration_v3_core.json").write_bytes(core)
    if args.core_only:
        log(f"wrote core {args.out_dir / 'registration_v3_core.json'} sha256 {sha256_bytes(core)}")
        return
    full = dict(reg)
    full["experiments"] = json.loads(args.experiments.read_text())
    flat = sorted(p for v in reg["panels"]["primary"]["proteins"].values() for p in v)
    e2genes = sorted(full["experiments"]["E2"]["perturbations"]["genes"]["panel_coding_genes"])
    if e2genes != flat:
        raise SystemExit(f"experiments_v3.json E2 panel_coding_genes {e2genes} != primary panel {flat}")
    full["experiments"]["E4"]["channels"]["channel2_input_proteins"] = [n for n in reg["evidence"]["teddy_head"]["adt_names"]
                                                                        if n not in flat]
    if len(full["experiments"]["E4"]["channels"]["channel2_input_proteins"]) != 122:
        raise SystemExit("E4 channel 2 should read the 122 non-panel proteins")
    full["provenance"] = {
        "builder": BUILDER_VERSION, "builder_sha256": vk.sha256_file(Path(__file__)),
        "v3_key_sha256": vk.sha256_file(Path(vk.__file__)),
        "experiments_sha256": vk.sha256_file(args.experiments),
        "core_sha256": sha256_bytes(core),
        "inputs": {"cite_arrays.npz": str(args.processed / "cite_arrays.npz"), "z": str(args.z), "ckpt": str(args.ckpt),
                   "ckpt_sha256": vk.sha256_file(args.ckpt), "z_sha256": vk.sha256_file(args.z),
                   "cite_arrays_sha256": vk.sha256_file(args.processed / "cite_arrays.npz")},
        "git_head_at_build": git_head(ROOT),
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "test_rows_read": "split/site/donor metadata only (see v3_leakage_check.py)",
    }
    out = args.out_dir / "registration_v3.json"
    out.write_bytes(content_bytes(full))
    h = vk.sha256_file(out)
    (args.out_dir / "registration_v3.json.sha256").write_text(f"{h}  registration_v3.json\n")
    (args.out_dir / "registration_v3_devlog.json").write_bytes(content_bytes(dev))
    log(f"wrote {out} sha256 {h}")


if __name__ == "__main__":
    main()
