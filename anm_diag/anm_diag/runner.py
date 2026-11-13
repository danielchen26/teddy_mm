"""run(): check the declaration, compute the tests, write the record and the report."""
import json
import platform
import time
from pathlib import Path

import numpy as np

from . import __version__
from .checks import CheckError, check, errors
from .report import write_report
from .spec import declaration_hash, load_declaration, node_map
from .stats import Crossfit, ci, cluster_bootstrap, forward_auc, theta, usable_folds
from .verdicts import root_cause, verdict_union


def _standardize(S, ref_rows):
    S = np.asarray(S, dtype=np.float32)
    if S.ndim == 1:
        S = S[:, None]
    mu = S[ref_rows].mean(0)
    sd = S[ref_rows].std(0)
    sd[sd < 1e-8] = 1.0
    return ((S - mu) / sd).astype(np.float32)


def run(declaration, adapter, out_dir=None, log=print):
    decl = load_declaration(declaration)
    problems = check(decl)
    if errors(problems):
        raise CheckError(problems)
    from .seal import verify
    sealed_checked = verify(decl, getattr(adapter, "paths", None))     # raises if data or code changed
    t0 = time.time()
    cells = adapter.cells(decl)
    A = np.asarray(cells["A"]).astype(int)
    cl = np.asarray(cells["cluster"]).astype(str)
    n_treated_clusters = len(np.unique(cl[A == 1]))
    problems = check(decl, n_treated_clusters=n_treated_clusters)
    inf, dec = decl["inference"], decl["decision"]
    alpha, margin = float(dec["alpha"]), float(dec["margin_u"])
    classes = decl["readout_classes"]
    fwd_cls = inf.get("forward_class", "ridge")
    nodes = node_map(decl)
    chain = decl["chain"]
    ref_rows = (A == 0) if inf.get("standardize_on", "control") == "control" else slice(None)
    log(f"[{decl['name']}] cells {len(A)} (treated {A.sum()}, clusters {n_treated_clusters}); chain {chain}")

    fold, fold_seed = usable_folds(cl, A, int(inf["folds"]), int(inf["seed"]))
    if fold_seed != int(inf["seed"]):
        msg = (f"with the declared seed {int(inf['seed'])} the treated clusters fell into one fold, which leaves nothing to "
               f"train on; the folds use seed {fold_seed} (the resampling keeps the declared seed). The choice reads cluster "
               f"labels and the treatment indicator only")
        problems.append({"code": "K1", "severity": "warning", "message": msg})
        log("  " + msg)
    boots = cluster_bootstrap(A, cl, int(inf["bootstrap"]), int(inf["seed"]))
    fboots = boots[: int(inf["forward_bootstrap"])]
    cf = Crossfit(A, fold, classes)

    state_ids = [n for n in nodes if nodes[n]["role"] in ("chain", "probe", "technical")]
    states = {}
    for sid in state_ids:
        states[sid] = _standardize(adapter.state(nodes[sid], cells), ref_rows)
        log(f"  state {sid}: {states[sid].shape[1]} features ({time.time() - t0:.0f}s)")

    # forward: does the source reach each state?
    fwd, fwd_b = {}, {}
    for sid in state_ids:
        fwd[sid], fwd_b[sid] = forward_auc(A, cf.ra(sid, states[sid], fwd_cls), fboots)
    fwd_steps = {}
    for a, b in zip(chain[:-1], chain[1:]):
        fwd_steps[f"{a}->{b}"] = ci(fwd_b[b] - fwd_b[a], alpha)
    for sid in state_ids:
        if nodes[sid]["role"] == "probe" and nodes[sid]["from"] in fwd_b:
            fwd_steps[f"{nodes[sid]['from']}->{sid}"] = ci(fwd_b[sid] - fwd_b[nodes[sid]["from"]], alpha)
    log("  forward AUC: " + "  ".join(f"{k} {v:.3f}" for k, v in fwd.items()))

    results = {"readouts": {}}
    boot_store = {}
    for r in decl["readouts"]:
        rid = r["id"]
        Y = np.asarray(adapter.readout(r, cells), dtype=float)
        y0, a0 = Y - Y.mean(), A - A.mean()
        d0 = theta(y0, a0)
        d0b = np.array([theta(y0, a0, s) for s in boots])
        delta_ci = ci(d0b, alpha)
        boot_store[f"{rid}|delta"] = d0b.astype(np.float32)
        if delta_ci[0] <= 0 <= delta_ci[1]:
            # no resolved response: u is undefined (it divides by Delta), so no inverse fits are made
            rc = root_cause(delta_ci, None, chain, {c: "n/a" for c in chain})
            results["readouts"][rid] = {"delta": d0, "delta_ci": delta_ci, "states": {}, "steps": {},
                                        "root_cause": rc, "exact": None}
            log(f"  {rid}: delta {d0:+.3f} {np.round(delta_ci, 3)} | {rc['cause']}")
            continue
        rstates = dict(states)
        for pid, pn in nodes.items():
            if pn["role"] == "positive_control" and pn.get("readout", rid) == rid:
                rstates[pid] = _standardize(adapter.state(pn, cells), ref_rows)
        U, Ub = {}, {}
        for sid, S in rstates.items():
            for cls in classes:
                ry = cf.ry(S, Y, cls)
                ra = cf.ra(sid, S, cls)
                U[f"{sid}/{cls}"] = theta(ry, ra) / d0
                Ub[f"{sid}/{cls}"] = np.array([theta(ry, ra, s) for s in boots]) / d0b
        for k, v in Ub.items():
            boot_store[f"{rid}|{k}"] = v.astype(np.float32)
        per = {}
        for sid in rstates:
            byc = {c: {"u": U[f"{sid}/{c}"], "ci": ci(Ub[f"{sid}/{c}"], alpha)} for c in classes}
            v, best = verdict_union(byc, margin)
            per[sid] = {"role": nodes[sid]["role"], "best_class": best, "u": byc[best]["u"], "ci": byc[best]["ci"],
                        "verdict": v, "by_class": byc}
        steps = {}
        pairs = list(zip(chain[:-1], chain[1:])) + [(nodes[s]["from"], s) for s in rstates
                                                     if nodes[s]["role"] == "probe" and nodes[s]["from"] in rstates]
        for a, b in pairs:
            ba, bb = per[a]["best_class"], per[b]["best_class"]
            steps[f"{a}->{b}"] = {"best_class": ci(Ub[f"{b}/{bb}"] - Ub[f"{a}/{ba}"], alpha),
                                  "same_class": {c: ci(Ub[f"{b}/{c}"] - Ub[f"{a}/{c}"], alpha) for c in classes}}
        pcs = [p for p in rstates if nodes[p]["role"] == "positive_control"]
        pc_v = per[pcs[0]]["verdict"] if pcs else None
        rc = root_cause(delta_ci, pc_v, chain, {c: per[c]["verdict"] for c in chain},
                        step_cis={k: v["best_class"] for k, v in steps.items()},
                        step_margin=float(dec.get("step_margin", margin)))
        exact = adapter.exact(r, cells, d0)
        results["readouts"][rid] = {"delta": d0, "delta_ci": delta_ci, "states": per, "steps": steps,
                                    "root_cause": rc, "exact": exact}
        log(f"  {rid}: delta {d0:+.3f} {np.round(delta_ci, 3)} | " +
            "  ".join(f"{s} {per[s]['u']:+.2f} {per[s]['verdict']}" for s in per) + f" | {rc['cause']}")

    record = {
        "anm_diag_version": __version__, "declaration": decl, "declaration_sha256": declaration_hash(decl),
        "problems": problems, "n_cells": int(len(A)), "n_treated": int(A.sum()),
        "n_treated_clusters": int(n_treated_clusters), "fold_seed": int(fold_seed), "forward": {"auc": fwd, "steps": fwd_steps, "class": fwd_cls},
        **results, "provenance": {**adapter.provenance(), "python": platform.python_version()},
        "runtime_s": round(time.time() - t0, 1),
        "exploratory": decl.get("status") != "sealed", "seal_verified": sealed_checked,
    }
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "anm_record.json").write_text(json.dumps(record, indent=1, default=float))
        np.savez_compressed(out / "anm_bootstrap.npz", **boot_store)    # for re-reading under other margins
        write_report(record, out / "anm_report.md")
    return record
