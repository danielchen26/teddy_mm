#!/usr/bin/env python3
"""Build amendment A1 to the v3 registration (registration/amendment_A1.json) from train + val only.

A1 records the corrections of the adversarial review of registration_v3.json, before any
site4 evaluation. registration_v3.json is not changed; A1 overrides the experiment fields
listed in DECLARED_OVERRIDES and adds a few numbers computed here:

* classifier C values (val log-loss) for the information-matched Q3 arms (A1.6);
* val diagnostics that motivate each fix: score ties at the clip (A1.1), tied top markers
  and the flip cells moved by the exact leave-one-out form (A1.3), the per-cell CLR form of
  the stored ADT and the ridge fit with and without the panel-dependent CLR scalar (A1.2).

As in v3_build_registration.py, the protein matrix, the cell types and the embedding are
read only on non-test rows (``nt``); site4 contributes split / site / donor metadata only.
``bridge_anm/v3_leakage_check.py`` rebuilds this core on site4-poisoned inputs and requires
it to be byte-identical.

Usage:
  ANM_ROOT=... python bridge_anm/v3_build_amendment.py            # writes amendment_A1.json(.sha256) + AMENDMENT_A1.md
  python bridge_anm/v3_build_amendment.py --core-only --out-dir D   # computed core only (leakage check)
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

BUILDER_VERSION = "v3_build_amendment A1 1.0"
TIE_BREAK_SEED = 29          # declared; not used by any other v3 seed
CLF_C_GRID = [0.01, 0.1, 1.0, 10.0, 100.0]   # same grid as the registration
CLF_SEED = 0
RIDGE_DIAG = {"alpha": 10.0, "n_train": 20000, "seed": 0}
CLR_DIAG = {"n_cells": 2000, "seed": 5, "max_k": 5, "rel_tol": 1e-3}
E2_MIN_LOST = 30


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


def runs(a: np.ndarray) -> int:
    a = np.asarray(a)
    return int(1 + np.sum(a[1:] != a[:-1])) if a.size else 0


# ============================================================================ computed core
def build_core(args) -> dict:
    t0 = time.time()
    reg = vk.load_registration(args.registration)
    npz = np.load(args.processed / "cite_arrays.npz", allow_pickle=False)
    split = npz["split"].astype(str)
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
    sp = split[nt]
    tr, va = sp == "train", sp == "val"
    j = {n: i for i, n in enumerate(adt_names)}
    log(f"loaded {nt.size} train+val rows of {n_all}")

    keys = vk.build_keys(ct, adt, adt_names, reg)
    prim, q3k = keys["primary"], keys["q3"]
    pred = vk.head_predict(z, args.ckpt, device="cpu", size_factor=1.0)
    v = vk.evidence(pred, reg)
    core: dict = {"amendment_id": "teddy_mm_v3_A1", "builder": BUILDER_VERSION}

    # ---------------------------------------------------------------- A1.1 ties and row order
    S1 = vk.class_scores(v[va], adt_names, reg, "Q1")
    S3 = vk.class_scores(v[va], adt_names, reg, "Q3")
    tie = {}
    for q, S in (("Q1", S1), ("Q3", S3)):
        top = S.max(axis=1)
        _, cnt = np.unique(top, return_counts=True)
        tie[q] = {"share_top_score_exactly_1": float(np.mean(top >= 1.0)),
                  "share_cells_in_a_tie_group": float(cnt[cnt > 1].sum() / top.size),
                  "largest_tie_group_share": float(cnt.max() / top.size)}
    tie["Q1"]["share_margin_exactly_0"] = float(np.mean(vk.margin(S1) == 0.0))
    te_rows = np.where(split == "test")[0]
    core["matched_coverage"] = {
        "tie_break": {"method": "rank of the global cell index in numpy default_rng(seed).permutation(n_cells)",
                      "seed": TIE_BREAK_SEED, "n_cells": n_all},
        "val_ties": tie,
        "test_row_order_metadata": {"n_test_rows": int(te_rows.size),
                                    "test_rows_contiguous": bool(np.all(np.diff(te_rows) == 1)),
                                    "donor_blocks_in_row_order": runs(donors[te_rows])},
    }
    log(f"ties on val: {tie}")

    # ---------------------------------------------------------------- A1.3 leave-one-out on val (Q1)
    fld = reg["anm"]["field_representation"]
    pan = reg["panels"]["primary"]["proteins"]
    n_ev = len(pan["B"])
    bar1 = float(reg["questions"]["Q1"]["bar"])
    calls = vk.rule_calls(S1, bar1)
    called = np.where(calls != vk.NO_CALL)[0]
    ci = np.asarray([vk.LINEAGES.index(c) for c in calls[called]], dtype=np.int64)
    Vc = v[va][called]
    vcal = np.zeros((called.size, n_ev))
    for k, lin in enumerate(vk.LINEAGES):
        rows = ci == k
        vcal[rows] = Vc[rows][:, [j[p] for p in pan[lin]]]
    sets = va_.top_marker_sets(vcal)
    exact = va_.loo_flip_exact(S1[called], vcal, ci, bar1, fld)
    old = va_.loo_flip_registered_form(S1[called], vcal, ci, bar1)
    core["e1_3"] = {
        "rho_G_n_minus_1_over_G_n": {"n": n_ev, "rho": va_.loo_ratio(n_ev, fld),
                                     "G_n": vk.field_gain(n_ev, fld), "G_n_minus_1": vk.field_gain(n_ev - 1, fld)},
        "val_Q1": {"n_calls": int(called.size),
                   "share_calls_with_tied_top_marker": float(np.mean([s.size > 1 for s in sets])) if sets else None,
                   "n_flip_exact": int(exact.sum()), "n_flip_registered_form": int(old.sum()),
                   "n_calls_flip_status_differs": int(np.sum(exact != old))},
    }
    log(f"E1.3 val: {core['e1_3']}")

    # ---------------------------------------------------------------- A1.2 per-cell CLR and channel 2
    rng = np.random.default_rng(CLR_DIAG["seed"])
    samp = rng.choice(adt.shape[0], size=min(CLR_DIAG["n_cells"], adt.shape[0]), replace=False)
    ok = 0
    for i in samp:
        u = np.expm1(adt[i].astype(np.float64))
        nz = u[u > 0]
        if nz.size == 0:
            ok += 1
            continue
        for k in range(1, CLR_DIAG["max_k"] + 1):
            x = u * (k / nz.min())
            if np.all(np.abs(x - np.round(x)) <= CLR_DIAG["rel_tol"] * np.maximum(1.0, x)):
                ok += 1
                break
    flat = [p for k in vk.LINEAGES for p in pan[k]]
    pj = [j[p] for p in flat]
    nj = [i for i, n in enumerate(adt_names) if n not in set(flat)]
    from sklearn.linear_model import Ridge

    a64 = adt.astype(np.float64)
    q95m = np.percentile(a64[tr][:, pj], 95, axis=0)
    Y = np.clip(a64[:, pj] / np.maximum(q95m, 1e-12), 0.0, 1.0)
    X_stored = a64[:, nj]
    X_free, _ = va_.channel2_inputs(a64, adt_names, flat)
    u_np = np.maximum(np.expm1(X_stored), 0.0)
    clr_scalar = np.log(np.maximum(u_np.sum(axis=1), 1e-12))[:, None]  # log(sum_nonpanel x) - log g_c
    tri = np.where(tr)[0]
    rs = np.random.default_rng(RIDGE_DIAG["seed"]).choice(tri, size=min(RIDGE_DIAG["n_train"], tri.size), replace=False)
    vai = np.where(va)[0]

    def r2(X):
        m = Ridge(alpha=RIDGE_DIAG["alpha"]).fit(X[rs], Y[rs])
        P = m.predict(X[vai])
        return float(np.mean([1 - np.sum((Y[vai, k] - P[:, k]) ** 2) / np.sum((Y[vai, k] - Y[vai, k].mean()) ** 2)
                              for k in range(Y.shape[1])]))

    core["e4_channel2"] = {
        "stored_adt_form": "per-cell CLR: expm1(stored) is a cell's integer count vector divided by one per-cell factor",
        "clr_integrality_check": {**CLR_DIAG, "share_cells_integral": ok / samp.size},
        "ridge_val_mean_r2": {**RIDGE_DIAG, "stored_nonpanel_values": r2(X_stored), "panel_free_inputs": r2(X_free),
                              "panel_free_inputs_plus_clr_scalar": r2(np.c_[X_free, clr_scalar])},
        "n_channel2_inputs": len(nj),
    }
    log(f"E4 channel 2 diagnostics: {core['e4_channel2']}")

    # ---------------------------------------------------------------- A1.6 information-matched Q3 classifiers
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss

    q3_feats = list(reg["classifier"]["q3"]["feature_proteins"])
    anchors = [reg["questions"]["Q3"]["anchors"][k] for k in vk.LINEAGES]

    def pick_c(feats, ytr, yva):
        X = v[:, [j[p] for p in feats]]
        mtr, mva = tr & (ytr != vk.UNSCORED), va & (yva != vk.UNSCORED)
        res = []
        for C in CLF_C_GRID:
            m = LogisticRegression(C=C, max_iter=3000, random_state=CLF_SEED).fit(X[mtr], ytr[mtr])
            res.append({"C": C, "val_log_loss": float(log_loss(yva[mva], m.predict_proba(X[mva]), labels=m.classes_))})
        best = min(res, key=lambda r: (r["val_log_loss"], r["C"]))["C"]
        return {"features": feats, "C_grid": CLF_C_GRID, "val_log_loss": res, "C": best,
                "n_train_labelled": int(mtr.sum()), "n_val_labelled": int(mva.sum())}

    core["classifiers"] = {
        "q1labels_q3features": {**pick_c(q3_feats, prim, prim),
                                "labels": "training-cell primary key (Q1 labels, 5 classes incl. OUT)"},
        "q1labels_anchors": {**pick_c(anchors, prim, prim),
                             "labels": "training-cell primary key (Q1 labels, 5 classes incl. OUT)"},
        "q3labels_anchors": {**pick_c(anchors, q3k, q3k), "labels": "training-cell q3 key"},
        "model": "multinomial logistic regression (sklearn LogisticRegression, lbfgs, max_iter 3000, random_state 0); "
                 "C by lowest val log-loss on labelled val cells, ties -> smaller C (as registration classifier)",
    }
    log("classifiers: " + ", ".join(f"{k} C={v_['C']}" for k, v_ in core["classifiers"].items() if isinstance(v_, dict)))
    log(f"core built in {time.time() - t0:.1f}s")
    return r6(core)


# ============================================================================ declared content
def declared(core: dict) -> dict:
    """Findings and overrides. Numbers quoted in text come from the computed core."""
    t1 = core["matched_coverage"]["val_ties"]["Q1"]["share_top_score_exactly_1"]
    t3 = core["matched_coverage"]["val_ties"]["Q3"]["share_top_score_exactly_1"]
    blocks = core["matched_coverage"]["test_row_order_metadata"]["donor_blocks_in_row_order"]
    e13 = core["e1_3"]
    rho = e13["rho_G_n_minus_1_over_G_n"]["rho"]
    d13 = e13["val_Q1"]["n_calls_flip_status_differs"]
    n13 = e13["val_Q1"]["n_calls"]
    tied13 = e13["val_Q1"]["share_calls_with_tied_top_marker"]
    c2 = core["e4_channel2"]["ridge_val_mean_r2"]
    cl = core["classifiers"]
    seed = core["matched_coverage"]["tie_break"]["seed"]
    ncell = core["matched_coverage"]["tie_break"]["n_cells"]

    findings = [
        {"id": "A1.1", "area": "matched coverage (E1.1c/d, E1.3c, E1.4, E1.5, C2, E3.H3b, E4)",
         "finding": f"the evidence is clipped at 1, so many cells share the top score (val: {t1:.4f} of cells have a Q1 top "
                    f"score of exactly 1, {t3:.4f} for Q3), and the site4 rows are stored in {blocks} contiguous donor "
                    "blocks; the registered tie rule (global cell index) would call one donor's tied cells first. "
                    "Whether matched-coverage selections are redone in each bootstrap replicate was not stated",
         "fix": f"ties are broken by the rank of the global cell index in numpy default_rng({seed}).permutation({ncell}), "
                "the same for every method (v3_amend.coverage_select); every endpoint, including the matched-coverage "
                "selection and any realised coverage that defines a matched point, is recomputed in every bootstrap replicate"},
        {"id": "A1.2", "area": "E4 channel 2",
         "finding": "the stored ADT is a per-cell CLR (expm1 of a row is an integer count vector over one per-cell "
                    f"factor; {core['e4_channel2']['clr_integrality_check']['share_cells_integral']:.4f} of "
                    f"{core['e4_channel2']['clr_integrality_check']['n_cells']} sampled train/val cells pass that check), "
                    "and the factor is computed from all proteins, the 12 panel proteins included. So the 122 stored "
                    "non-panel values channel 2 reads carry a function of the measured panel counts, contrary to "
                    "'it never sees a panel protein's measured value'. Val ridge diagnostic (mean R2 over the 12 panel "
                    f"targets): stored inputs {c2['stored_nonpanel_values']:.4f}, panel-free inputs "
                    f"{c2['panel_free_inputs']:.4f}, panel-free inputs plus the CLR scalar "
                    f"{c2['panel_free_inputs_plus_clr_scalar']:.4f}",
         "fix": "channel 2 reads w_q = log1p(1e4 * u_q / sum_r u_r), u = expm1(stored value), q and r over the 122 "
                "non-panel proteins (v3_amend.channel2_inputs). Because u = x / g_c, w equals log1p(1e4 * x_q / sum_r x_r) "
                "and depends on the non-panel counts only (unit test: changing a cell's panel counts leaves w unchanged)"},
        {"id": "A1.3", "area": "E1.3 why this call",
         "finding": "ANM's leave-one-out deletes the event's site (build_graph makes one site per event; v2's "
                    "_loo_top_for_instance drops the event), so the called action's gain changes from G(n) to G(n-1). "
                    f"The registered closed form (S_called - v_top/n against the bar) re-codes zeroing, not deletion; "
                    f"rho = G(n-1)/G(n) = {rho:.6f} for n = 3, and on val {d13} of {n13} Q1 calls change flip status. "
                    f"Also, {tied13:.4f} of val Q1 calls have two or more panel markers tied at the called class's "
                    "maximum (mostly clipped at 1), so a single 'top-1 marker' would be set by event order. The question "
                    "(Q1 or Q2) for E1.3 was not stated",
         "fix": "exact closed form flip iff rho * (S_called - v_top / n) < max(bar, max_{k != called} S_k) "
                "(v3_amend.loo_flip_exact); top markers are sets (all markers at the class maximum), E1.3a compares sets, "
                "E1.3b gives each marker of a set of m markers 1/m credit (and reports single-marker calls separately); "
                "E1.3 runs on Q1 (primary) and Q2 (secondary)"},
        {"id": "A1.4", "area": "E1.4 honest naming",
         "finding": "the ANM soft score is G(3) * 3 times the rule's top score, so its ranking and every E1.4 number equal "
                    "the declared rule's; the registered text read an E1.4 win as support for an ANM gate claim. Q1 and Q2 "
                    "share scores and key, so their E1.4 numbers are identical",
         "fix": "an E1.4 win is reported as a win of the top-score readout (rule = ANM), never as an ANM-specific property; "
                "E1.4 is computed once (Q1 = Q2); the E1 falsification sentence is restated accordingly"},
        {"id": "A1.5", "area": "E1.C2 primary point",
         "finding": "C2 listed three coverage points without saying which decides the verdict, and did not define "
                    "'in-scope selective accuracy'",
         "fix": "the verdict uses c* = the mean rule's realised Q1 coverage on test_primary (label-free); 0.8 and 0.7 are "
                "secondary. In-scope selective accuracy = among called cells whose key class is a lineage, the fraction "
                "correct; guard = point difference >= -0.005 at c*"},
        {"id": "A1.6", "area": "E1.1c / E1.1d / E1.5 information parity",
         "finding": "the Q3 rule reads CD19, CD3, CD56 and CD14 evidence; the registered 'Q1 classifier without relabelling' "
                    "reads the 12 primary-panel values only (no CD19, no CD14), while the Q3 classifier reads 14 values. "
                    "E1.5 therefore changed features as well as labels between n = 0 and n > 0",
         "fix": "the Q1-label arm for Q3 (E1.1c and E1.5 n = 0) is a classifier on the classifier.q3 features (14 values, "
                f"a superset of the rule's 4 anchors) trained on Q1 labels, C = {cl['q1labels_q3features']['C']} on val; "
                "features are the same at every n. Secondary matched-information rows use exactly the rule's 4 anchor "
                f"values: Q1-label C = {cl['q1labels_anchors']['C']}, Q3-label C = {cl['q3labels_anchors']['C']} (val)"},
        {"id": "A1.7", "area": "E2 falsification (i)",
         "finding": "criterion (i) compares every pair of perturbations without a minimum count of lost cells or an "
                    "uncertainty requirement, so with many pairs and few lost cells per pair noise alone can pass it",
         "fix": f"only perturbation levels with at least {E2_MIN_LOST} lost cells among the 1,000 test_primary e2_subset "
                "cells (in both perturbations) enter (i), and the share difference needs a two-stage bootstrap 95% "
                "interval excluding 0; the number of pairs compared is reported"},
        {"id": "A1.8", "area": "E3 pairs",
         "finding": "neighbours 'within the evaluated split' pool the two primary donors, so cross-donor pairs exist, "
                    "while the registered bootstrap resamples pairs within donors and the flag cosine came from one donor",
         "fix": "k = 10 cosine neighbours are computed within each donor of the evaluated split; every pair belongs to "
                "one donor"},
        {"id": "A1.9", "area": "E4 comparator",
         "finding": "'F5 vs the best of F1-F4' did not say how 'best' is chosen; choosing it on site4 is a post-hoc "
                    "maximum",
         "fix": "the comparator is the F1-F4 method with the highest val endpoint (selective accuracy at coverage 0.8, "
                "mean over L0-L3; F3 scored by 5-fold cross-fitting on val with seeds.e4_noise), written to the E4 "
                "addendum before site4"},
        {"id": "A1.10", "area": "addenda",
         "finding": "the untracked E5 addendum replaces the registered primary endpoint (log-gain at l*) with dG_l and "
                    "uses pooled site4 annotation-key pairs; the registration allows addenda only to fix numbers it "
                    "leaves open",
         "fix": "an addendum cannot replace a registered endpoint, key or cell set; the registered E5 endpoint stays "
                "primary and dG_l is secondary unless a further amendment is committed before site4"},
    ]
    review_notes = [
        "Class map: all 45 types (40 in site4) reviewed. Every assignment is defensible. pDC -> OUT is the most "
        "debatable (some schemes count pDC as myeloid DC); it only enters the primary key when no lineage gate passes, "
        "and annotation_only is reported.",
        "Myeloid panel member CD62P is a platelet protein; in E6 the external set has a 'Platelet' label keyed OUT, so "
        "myeloid calls on platelets are a predictable E6 error mode (report it, do not change the panel).",
        "canonical9 sensitivity panels were chosen on site4 in v1/v2 (already disclosed); they never decide a verdict.",
        "Splits: train donors 10886, 11466, 12710, 15078, 16710, 28045 (sites 1-3); val 18303 (site1 only); "
        "test_primary 13272, 19593 (site4 only, in no other split); test_secondary 15078. No donor-split error found.",
        "Leakage: no site4 protein, label or embedding enters the registration core or this amendment core "
        "(v3_leakage_check.py). The stored ADT CLR is per cell, so no cross-cell (site4) information enters train values.",
        "E6 prep (commit 8318ef4): L2_CLASS is written from label names, but it was committed with the code that reads "
        "external values (already disclosed); the E6 runner must report this.",
    ]
    overrides = {
        "common/matched_coverage/rule":
            "at coverage c each method calls its ceil(c * N) most confident cells of the evaluated split (N = all cells "
            "of the split), using only its own scores, never labels; ties in a method's score are broken by the A1.1 "
            f"tie-break (rank of the global cell index in numpy default_rng({seed}).permutation({ncell})), the same "
            "permutation for every method (v3_amend.coverage_select)",
        "common/matched_coverage/tie_break": {"seed": seed, "n_cells": ncell, "function": "v3_amend.coverage_select"},
        "common/statistics/bootstrap":
            "two-stage: resample donors with replacement, then cells (or pairs) with replacement within each drawn donor; "
            "B = 2000; seed = registration seeds.bootstrap; percentile 95% interval. Every endpoint is recomputed inside "
            "each replicate, including the matched-coverage selection (N = replicate size; a repeated cell keeps the "
            "A1.1 rank of its original id) and any method's realised coverage that defines a matched point (A1.1)",
        "E1/exp1_change_the_question/arms": [
            "rule (Q3 anchors CD19, CD3, CD56, CD14)", "anm",
            f"Q1-label classifier for Q3 (A1.6): classifier.q3 features, Q1 labels, C = {cl['q1labels_q3features']['C']}; "
            "its 'myeloid' is read as the Q3 myeloid call",
            "Q3-label classifier (registration classifier.q3; all training-site Q3 labels; label-rich ceiling)",
            f"secondary (A1.6, same information as the rule): the two classifiers on the 4 anchor values only "
            f"(Q1 labels C = {cl['q1labels_anchors']['C']}, Q3 labels C = {cl['q3labels_anchors']['C']})"],
        "E1/exp1_change_the_question/endpoints/E1.1c":
            "Q3 selective accuracy on the q3 key at the rule's realised Q3 coverage: rule vs the A1.6 Q1-label classifier "
            "(classifier.q3 features) without relabelling; secondary row with the 4-anchor Q1-label classifier",
        "E1/exp1_change_the_question/endpoints/E1.1d":
            "same, rule vs the Q3-label classifier (label-rich ceiling); secondary row with the 4-anchor Q3-label classifier",
        "E1/exp3_why_this_call/question": "Q1 on test_primary (primary); Q2 with its bar (secondary)",
        "E1/exp3_why_this_call/anm":
            "leave-one-out over the called class's events (events at t = 0): removing an event deletes its event site "
            "(ANM build_graph makes one site per event, as v2's _loo_top_for_instance); top-1 marker set = every event "
            "whose removal lowers the called action's score most (ties kept as a set); flip = removing a top-1 event "
            "changes the call (to another class or to no call)",
        "E1/exp3_why_this_call/rule":
            "exact closed form (A1.3, v3_amend.loo_flip_exact): top-1 marker set = the panel proteins at the called "
            "class's maximum evidence; flip iff rho * (S_called - v_top / n) < max(bar, max over other classes of S_k), "
            f"rho = G(n-1)/G(n) = {rho:.6f} for n = 3",
        "E1/exp3_why_this_call/endpoints": {
            "E1.3a": "agreement of ANM's top-1 marker set with the closed-form set (pre-registered expectation 1.000; "
                     "below 0.99 is reported as a discrepancy), and agreement of the flip status (expectation 1.000); "
                     "the share of calls whose top-1 set has more than one marker is reported",
            "E1.3b": "per-donor share of each top-1 marker with fractional credit (a set of m markers gives 1/m to each); "
                     "null = within-class permutation of marker identities per cell, 1000 permutations, same credit; "
                     "also reported on calls with a single top-1 marker",
            "E1.3c": "error rate of flip-sensitive calls (exact flip, A1.3) vs the same number of lowest-margin calls "
                     "(margin ties by the A1.1 permutation), Q1, test_primary"},
        "E1/exp4_which_calls_to_trust/interpretation":
            "the ANM soft score is G(3) * 3 times the rule's top score, so its ranking and every E1.4 number equal the "
            "declared rule's top score; an E1.4 win is reported as a win of the top-score readout (rule = ANM), never "
            "as an ANM-specific property. Q1 and Q2 share scores and key, so E1.4 is computed once and reported once",
        "E1/exp4_which_calls_to_trust/failure":
            "if the top score (rule = ANM) is not a win against the margin, the gate claim stays withdrawn; a win "
            "supports the top-score readout, not an ANM-specific mechanism",
        "E1/c2_nested_readouts/primary_coverage":
            "the verdict uses c* = the mean rule's realised Q1 coverage on test_primary (label-free; recomputed in each "
            "bootstrap replicate); coverage 0.8 and 0.7 are secondary and never change the verdict",
        "E1/c2_nested_readouts/in_scope_accuracy":
            "among called cells whose key class is B, T, NK or myeloid, the fraction whose call equals the key; the "
            "guard is the point difference (ANM closure minus mean rule) >= -0.005 at c*",
        "E1/falsification":
            "Mode B cannot show ANM-specific decision value: every ANM arm equals a re-coded rule cell by cell (E1.1a; "
            "the closure rule in C2; the top score in E1.4). E1 tests whether the readout forms ANM provides (top "
            "action score, closure readout) add decision value over TEDDY's margin and the mean rule; that is falsified "
            "for v3 if neither E1.4 nor C2 is a win, and any win is credited to the readout form, which a written rule "
            "implements",
        "E1/exp5_label_cost/classifier":
            "features are the same at every n (A1.6): registration classifier.q3 features; n = 0 is the A1.6 Q1-label "
            f"classifier on them (C = {cl['q1labels_q3features']['C']}); n > 0 uses classifier.q3 (C = "
            "registration classifier.q3.C) on the drawn Q3 labels; a draw with a single class predicts that class. "
            "Secondary curve on the 4 anchor values only (n = 0: Q1 labels, C = "
            f"{cl['q1labels_anchors']['C']}; n > 0: C = {cl['q3labels_anchors']['C']})",
        "E2/falsification":
            "the decomposition adds information beyond the curve only if (i) two perturbations with accuracy losses "
            f"within 0.02 of each other, each with at least {E2_MIN_LOST} lost cells among the 1,000 test_primary "
            "e2_subset cells, differ by >= 0.15 in a localisation share with a two-stage bootstrap 95% interval of the "
            "difference excluding 0 (the number of pairs compared is reported), or (ii) the small-eps response (dS at "
            "eps = 0.05 thinning, linearly extrapolated) predicts which cells change call at eps = 0.8 with AUROC >= 0.02 "
            "above baseline margin alone (bootstrap lower bound > 0); if neither holds, it adds nothing beyond the curve",
        "E3/cells":
            "test_primary; pairs (A1.8): k = 10 cosine neighbours on the raw L2-normalised final z, computed within each "
            "donor of the evaluated split; each edge kept once; NK-T pair = one primary-key NK and one primary-key T "
            "cell of the same donor; flagged by registration e3.flag_cosine (from val); variants all and no_gdT158",
        "E4/channels/channel2":
            "ridge regression (alpha on val from [0.1, 1, 10, 100, 1000]) from the panel-free renormalised non-panel ADT "
            "(A1.2, v3_amend.channel2_inputs: w_q = log1p(1e4 * u_q / sum_r u_r), u = expm1(stored value), q and r over "
            "the 122 non-panel proteins) to the measured panel evidence; trained on training cells; normalised by the "
            "training q95 of its own predictions; it reads no function of a panel protein's measured value",
        "E4/noise_levels/L2": "50% of the 122 channel-2 inputs w set to 0 per cell (seeds.e4_noise), after the A1.2 transform",
        "E4/best_comparator":
            "the F1-F4 method with the highest val endpoint (selective accuracy at coverage 0.8, mean over L0-L3; F3 "
            "scored by 5-fold cross-fitting on val with seeds.e4_noise), fixed in the E4 addendum before site4 (A1.9)",
        "common/addenda_scope":
            "an addendum may only fix numbers the registration leaves open, by the procedure written for them; it cannot "
            "replace a registered endpoint, key or cell set (A1.10). The registered E5 endpoint stays primary unless a "
            "further amendment is committed before site4",
    }
    return {"title": "Amendment A1 to the v3 pre-registration: corrections from the adversarial review",
            "written_before": "any site4 (test) evaluation of v3; every number in 'computed' uses train and val cells "
                              "only (v3_leakage_check.py)",
            "scope": "registration_v3.json is unchanged; where A1 and registration_v3.json differ, A1 wins. Builders "
                     "load both with v3_amend.load_registration_amended()",
            "findings": findings, "review_notes": review_notes, "experiments_overrides": overrides}


# ============================================================================ render
def render_md(full: dict) -> str:
    L = [f"# {full['title']}", "",
         f"**Status:** written before {full['written_before']}. File `registration/amendment_A1.json`, sha256 in "
         "`amendment_A1.json.sha256`. It amends `registration/registration_v3.json` sha256 "
         f"`{full['amends']['registration_sha256']}`, which is unchanged. This page is rendered from the JSON by "
         "`bridge_anm/v3_build_amendment.py`.", "",
         f"**Scope.** {full['scope']}.", "", "## Findings and fixes", ""]
    for f in full["findings"]:
        cap = lambda t: t[:1].upper() + t[1:]  # noqa: E731
        L += [f"### {f['id']} {f['area']}", "", f"- **Finding.** {cap(f['finding'])}.", f"- **Fix.** {cap(f['fix'])}.", ""]
    c = full["computed"]
    L += ["## Numbers computed for A1 (train and val only)", "",
          "| item | value |", "|---|---|",
          f"| tie-break seed / cells | {c['matched_coverage']['tie_break']['seed']} / {c['matched_coverage']['tie_break']['n_cells']:,} |",
          f"| val share with Q1 top score exactly 1 | {c['matched_coverage']['val_ties']['Q1']['share_top_score_exactly_1']:.4f} |",
          f"| val share with Q3 top score exactly 1 | {c['matched_coverage']['val_ties']['Q3']['share_top_score_exactly_1']:.4f} |",
          f"| site4 donor blocks in row order (metadata) | {c['matched_coverage']['test_row_order_metadata']['donor_blocks_in_row_order']} |",
          f"| rho = G(2)/G(3) | {c['e1_3']['rho_G_n_minus_1_over_G_n']['rho']:.6f} |",
          f"| val Q1 calls / with tied top marker | {c['e1_3']['val_Q1']['n_calls']:,} / {c['e1_3']['val_Q1']['share_calls_with_tied_top_marker']:.4f} |",
          f"| val Q1 calls whose flip status the exact form changes | {c['e1_3']['val_Q1']['n_calls_flip_status_differs']} |",
          f"| sampled cells passing the per-cell CLR check | {c['e4_channel2']['clr_integrality_check']['share_cells_integral']:.4f} of {c['e4_channel2']['clr_integrality_check']['n_cells']:,} |",
          f"| val ridge mean R2: stored / panel-free / panel-free + CLR scalar | {c['e4_channel2']['ridge_val_mean_r2']['stored_nonpanel_values']:.4f} / {c['e4_channel2']['ridge_val_mean_r2']['panel_free_inputs']:.4f} / {c['e4_channel2']['ridge_val_mean_r2']['panel_free_inputs_plus_clr_scalar']:.4f} |",
          ]
    for k in ("q1labels_q3features", "q1labels_anchors", "q3labels_anchors"):
        d = c["classifiers"][k]
        ll = ", ".join(f"{r['C']}: {r['val_log_loss']:.6f}" for r in d["val_log_loss"])
        L.append(f"| classifier {k} ({len(d['features'])} features): C | {d['C']} (val log-loss {ll}) |")
    L += ["", "## Experiment fields replaced (path in registration_v3.json → experiments)", ""]
    for p, v in full["experiments_overrides"].items():
        txt = v if isinstance(v, str) else json.dumps(v, ensure_ascii=False)
        L.append(f"- `{p}`: {txt}")
    L += ["", "## Review notes (no change needed)", ""] + [f"- {n}" for n in full["review_notes"]]
    L += ["", "## Provenance", "",
          f"Builder `{full['provenance']['builder']}` sha256 `{full['provenance']['builder_sha256']}`; "
          f"`v3_amend.py` sha256 `{full['provenance']['v3_amend_sha256']}`; computed core sha256 "
          f"`{full['provenance']['core_sha256']}`.", ""]
    return "\n".join(L)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--z", type=Path, default=ROOT / "data/processed/cite_official/z_rna.npy")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/outputs/cite_phase1_official/best.pt")
    p.add_argument("--registration", type=Path, default=ROOT / "registration" / "registration_v3.json")
    p.add_argument("--out-dir", type=Path, default=ROOT / "registration")
    p.add_argument("--core-only", action="store_true", help="write only amendment_A1_core.json (leakage check)")
    args = p.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    core = build_core(args)
    core_b = content_bytes(core)
    (args.out_dir / "amendment_A1_core.json").write_bytes(core_b)
    if args.core_only:
        log(f"wrote core sha256 {sha256_bytes(core_b)}")
        return
    full = {"amendment_id": core["amendment_id"],
            "amends": {"file": "registration/registration_v3.json",
                       "registration_sha256": vk.sha256_file(args.registration)},
            **declared(core), "matched_coverage": {"tie_break": core["matched_coverage"]["tie_break"]},
            "computed": core,
            "provenance": {"builder": BUILDER_VERSION, "builder_sha256": vk.sha256_file(Path(__file__)),
                           "v3_amend_sha256": vk.sha256_file(Path(va_.__file__)),
                           "v3_key_sha256": vk.sha256_file(Path(vk.__file__)),
                           "core_sha256": sha256_bytes(core_b),
                           "inputs": {"cite_arrays_sha256": vk.sha256_file(args.processed / "cite_arrays.npz"),
                                      "z_sha256": vk.sha256_file(args.z), "ckpt_sha256": vk.sha256_file(args.ckpt)},
                           "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                           "test_rows_read": "split / donor metadata only (see v3_leakage_check.py)"}}
    out = args.out_dir / "amendment_A1.json"
    out.write_bytes(content_bytes(full))
    h = vk.sha256_file(out)
    (args.out_dir / "amendment_A1.json.sha256").write_text(f"{h}  amendment_A1.json\n")
    (args.out_dir / "AMENDMENT_A1.md").write_text(render_md(full))
    log(f"wrote {out} sha256 {h}")


if __name__ == "__main__":
    main()
