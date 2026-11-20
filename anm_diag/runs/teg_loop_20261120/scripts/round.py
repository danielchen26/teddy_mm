#!/usr/bin/env python3
"""One round of the loop of LOOP_PROTOCOL.md in which the change is a repair computed from the events of the table.

  python round.py --observe OBSERVE --readouts PDL1 [--tag development]

OBSERVE holds the output of `anm_diag scope` for every SPEC of build.py (the grid on the head as trained).
Environment: ANM_LOOP_WORK (the folder of build.py), ANM_BM_FIT, ANM_BM_TEST.

The five actions of the protocol, for every group (a cell type in a condition) and readout:

  1 observe      the grid of `anm_diag scope`: the head's own output and readouts of Z learned from its training cells, from
                 another group, from the other genes of the group, and with the gene's own cells.
  2 admit        the source chain of the paper (SI S1.4): propose (one typed event per gene, events.source_events), prepare (only
                 the half of the cells that is fitted; no measurement of the half that is read), admit (rules below), assign
                 (the displacement in Z, the measured change and its standard error: computed in the fit).
  3 hypothesize  what is missing is coverage: cells of this kind. Forecast for the repair by rows: the ceiling of the grid
                 (the head's own output against a readout of the other genes). For the repair of another group applied here:
                 the grid's `other` scope. Rule: the repair by rows comes within the step margin of the ceiling, and the
                 repair of another group does not.
  4 intervene    one change: inside a gate typed by cell type and condition, the output is replaced by a linear readout of Z
                 fitted on the group's control cells and updated once by its admitted events (as runs/teg_repair_test_20261119).
                 Two halves: every cell is read by a system fitted without it.
  5 test, update the head's own output again (the `frozen` scope), its gain, the loss it still adds beyond what a readout of
                 the other genes reads at Z, the closure grade, the guard, the ledger and its checks.

Admission (gate 0 of the protocol, on the grid of the head as trained):
  A1  the response is resolved and at least three times the half-width of its interval;
  A2  the event does not act on the readout's own gene (rule S3; events.no_leakage);
  A3  the head loses something a readout of cells of this kind reads: the gap from the head's own output to the `units`
      scope at Z has its interval above the step margin;
  A4  an event enters the repair of its own group only: the fields a pool merges are checked with events.closure, and
      a pool across conditions is flagged.
  Attainable grade: absolute if the `units` scope at Z is supported; otherwise relative (the output can be brought to what Z
  carries, not further).

Closure of the head's own output (the grades of the protocol, read on the frozen scope):
  absolute   the changed output is supported (two-sided) and the output as trained is rejected
  relative   the changed output adds no loss beyond the `units` scope at Z (interval of the difference below the step
             margin), and the output as trained adds one
  loss no longer resolved / not closed   otherwise, by the interval of that difference
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parents[1]))
sys.path.insert(0, str(RUN.parent / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(HERE))
from anm_diag import events as EV, loop as LP, scope as SC  # noqa: E402
from _common import BM_FIT, BM_TEST, TARGETS, forward_np, halves, interval, own_output_u, unit, weights  # noqa: E402
from build import DECISION, GROUPS, OUTPUT, PROBE  # noqa: E402

RES = RUN / "results"
GRID = [10.0 ** (k / 2) for k in range(-6, 5)]
FOLDS, SEED, GUARD, MARGIN, STEP, MIN_SNR = 5, 0, 0.02, DECISION["margin_u"], DECISION["step_margin"], 3.0
ARMS = ["other", "controls", "rows_other", "rows", "cells"]
READOUT_GENES = {"PDL1": ["CD274"], "MHC1": ["HLA-A", "HLA-B", "HLA-C"], "HLAE": ["HLA-E"]}
verdict = lambda iv: SC._verdict(iv, MARGIN, True)


# --- 1 observe --------------------------------------------------------------------------------------------------------
def observed(obs, group, readout):
    d = json.loads((obs / f"{group}_{readout}" / "scope.json").read_text())
    by = {r["event"]["unit"]: r for r in d["results"]}
    return by, np.load(obs / f"{group}_{readout}" / "scope_draws.npz")


def grid_row(r):
    T = r["targets"]["shift"]
    cell = lambda v: {"u": v["u"], "ci": v["ci"], "verdict": v["verdict"], "class": v["class"]}
    row = {"frozen": cell(T["stages"]["system"]["frozen"]), **{s: cell(v) for s, v in T["stages"]["Z"].items()}}
    half = (r["delta"]["ci"][1] - r["delta"]["ci"][0]) / 2
    return {"n_treated": r["event"]["n_treated"], "n_control": r["event"]["n_control"], "approximate": r["event"]["approximate"], "resampled": r["event"]["resampled"],
            "response": {**r["delta"], "snr": abs(r["delta"]["value"]) / half if half > 0 else float("inf")}, "scopes": row,
            "gaps": T["gaps"], "reading": T["reading"]["missing"], "detail": T["reading"]["detail"], "ceiling": SC.forecast(r)}


# --- 2 admit ----------------------------------------------------------------------------------------------------------
def admit(group, readout, grid):
    g = GROUPS[group]
    genes = TARGETS[g["study"]]
    decl = {"name": f"{group}_{readout}", "source": {"kind": "knockout", "treated": genes, "stratum": g["stratum"], "assign_from_guides": True, "control": "non-targeting guides"},
            "system": {"inputs": ["rna"]}, "readouts": [{"id": readout, "measurement": "protein", "features": READOUT_GENES[readout]}]}
    pool = grid["pool"]
    A1 = bool(pool["response"]["resolved"] and pool["response"]["snr"] >= MIN_SNR)
    gap = pool["gaps"]["from_frozen"]["units"]
    A3 = bool(gap["ci"][0] > STEP)
    ref = pool["scopes"]["units"]["verdict"]
    reasons = []
    if not A1:
        reasons.append(f"A1: the response is {'not resolved' if not pool['response']['resolved'] else 'small against its interval'} (|change| / half-width = {pool['response']['snr']:.1f}, wanted {MIN_SNR:g})")
    if not A3:
        reasons.append(f"A3: the head does not lose what cells of this kind read (gap {gap['value']:.2f} [{gap['ci'][0]:.2f}, {gap['ci'][1]:.2f}], step margin {STEP:g})")
    admitted = A1 and A3
    events = []
    for k, e in enumerate(ev for ev in EV.source_events(decl) if ev["role"] == "treated"):
        r = grid[e["site"]]
        ok1 = bool(r["response"]["resolved"] and r["response"]["snr"] >= MIN_SNR)
        ok2 = e["no_leakage"]["routes"]["site"]["status"] == "pass"
        why = [] if admitted else ["the group is not admitted"]
        why += [] if ok1 else ["A1: its own response is not resolved or small against its interval"]
        why += [] if ok2 else ["A2: " + e["no_leakage"]["routes"]["site"]["reason"]]
        e.update(time=k, sign=int(np.sign(r["response"]["value"])), amplitude=float(r["response"]["value"]), amplitude_unit="control SD", readout=readout,
                 cell_type=g["cell_type"], at_Z={s: r["scopes"][s]["verdict"] for s in ("units", "own")},
                 admission={"status": "admitted" if (admitted and ok1 and ok2) else "not admitted", "reason": "; ".join(why) or "A1 to A4 hold",
                            "n_cells": r["n_treated"]})
        events.append(e)
    merged = EV.closure(events).get(f"{group}_{readout}:A=1", [])
    if admitted and not any(e["admission"]["status"] == "admitted" for e in events):
        admitted = False                                    # ready as a pool, but no single event passes A1: there is no row to repair with
        reasons.append("no event of the group passes A1 on its own: there is no row to repair with")
    return {"group": group, "readout": readout, "cell_type": g["cell_type"], "stratum": g["stratum"], "span": ["Z", "own output"],
            "diagnostic": "ready" if admitted else "not ready", "attainable": ("absolute" if ref == "support" else "relative") if admitted else "none",
            "A1": A1, "A3": A3, "reference": {"scope": "units", "stage": "Z", **pool["scopes"]["units"]}, "own_cells": pool["scopes"]["own"],
            "gap_head_to_units": gap, "reasons": reasons, "fields_merged_by_the_pool": merged,
            "rest": None if ref == "support" else f"the `units` scope at Z is {ref}: what it leaves is before the head (upstream, or not readable from RNA; no input view was given)",
            "events": events}


# --- 4 intervene ------------------------------------------------------------------------------------------------------
def cv_weight(F, t, pen):
    fold = np.random.default_rng(SEED).permutation(len(F)) % FOLDS
    err = np.zeros(len(GRID))
    for k in range(FOLDS):
        tr, te = fold != k, fold == k
        M, v = F[tr].T @ F[tr] / tr.sum(), F[tr].T @ t[tr] / tr.sum()
        for i, lam in enumerate(GRID):
            err[i] += ((F[te] @ np.linalg.solve(M + lam * np.diag(pen), v) - t[te]) ** 2).sum()
    err /= len(F)
    return GRID[int(err.argmin())], float(err.min())


class Readout:
    """The readout of one group: fitted on its control cells C, updated by its events (rows), by the events of another group
    (rows_other), or fitted on every knockout cell (cells). x: unit embeddings, y: measured readout, o: the head's log output."""

    def __init__(self, x, y, o, C, events, other=None):
        self.xc = x[C].mean(0)
        F = self.feat(x)
        scale = o[C].std() / y[C].std()
        t = o[C].mean() + scale * (y - y[C].mean())
        pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
        lam, cv = cv_weight(F[C], t[C], pen)
        Minv = np.linalg.inv(F[C].T @ F[C] / len(C) + lam * np.diag(pen))
        self.th = {"controls": Minv @ (F[C].T @ t[C] / len(C))}
        self.log = {"n_controls": int(len(C)), "ridge_weight": lam, "cv_error": cv, "rows": []}

        def update(D, dy, se2):
            S = len(C) * scale ** 2 * se2 / cv
            K = Minv @ D.T @ np.linalg.inv(D @ Minv @ D.T + np.diag(S))
            before = scale * dy - D @ self.th["controls"]
            th = self.th["controls"] + K @ before
            return th, 1 - np.abs(scale * dy - D @ th) / np.maximum(np.abs(before), 1e-12)

        vc = y[C].var(ddof=1) / len(C)
        if events:
            D = np.array([F[T].mean(0) - F[C].mean(0) for _n, T in events])
            dy = np.array([y[T].mean() - y[C].mean() for _n, T in events])
            se2 = np.array([y[T].var(ddof=1) / len(T) + vc for _n, T in events])
            self.th["rows"], closed = update(D, dy, se2)
            self.log["rows"] = [{"event": n, "n_cells": int(len(T)), "measured_change": float(a), "standard_error": float(np.sqrt(s)), "share_of_the_residual_closed": float(c)}
                                for (n, T), a, s, c in zip(events, dy, se2, closed)]
            tr = np.r_[C, np.concatenate([T for _n, T in events])]
            lam2, _cv = cv_weight(F[tr], t[tr], pen)
            self.th["cells"] = np.linalg.solve(F[tr].T @ F[tr] / len(tr) + lam2 * np.diag(pen), F[tr].T @ t[tr] / len(tr))
            self.log["ridge_weight_cells"] = lam2
        if other is not None:                               # events of another group, in units of the control spread
            xo, yo, Ao, to, sites = other
            Fo, Co = self.feat(xo), np.where(Ao == 0)[0]
            conv = y[C].std() / yo[Co].std()
            ev = [(g, np.where((Ao == 1) & (to == g))[0]) for g in sites]
            D = np.array([Fo[T].mean(0) - Fo[Co].mean(0) for _n, T in ev])
            dy = conv * np.array([yo[T].mean() - yo[Co].mean() for _n, T in ev])
            se2 = conv ** 2 * np.array([yo[T].var(ddof=1) / len(T) + yo[Co].var(ddof=1) / len(Co) for _n, T in ev])
            self.th["rows_other"], _c = update(D, dy, se2)

    def feat(self, x):
        return np.hstack([x - self.xc, np.ones((len(x), 1))])

    def out(self, arm, x):
        return self.feat(x) @ self.th[arm]


class System:
    """The head with its repaired outputs: a gate over bone marrow and every group, and behind the gate of a repaired group
    its readout. Where the gate sees bone marrow, or a group without a repair, the head's own output stays.
    The gate assigns a cell to one class, the most probable one (a type is one of several, not a mixture). With soft=True the
    readouts are mixed by the probabilities instead (reported as a sensitivity)."""

    def __init__(self, xb, ctl, readouts, soft=False):
        self.names, self.soft = list(ctl), soft
        X = np.vstack([xb] + [ctl[g] for g in self.names])
        lab = np.concatenate([np.zeros(len(xb))] + [np.full(len(ctl[g]), i + 1) for i, g in enumerate(self.names)])
        self.clf = LogisticRegression(C=1.0, max_iter=3000).fit(X, lab)
        self.readouts = readouts                            # {group: Readout}, the repaired groups

    def gate(self, x):
        p = self.clf.predict_proba(x)
        if not self.soft:
            p = np.eye(p.shape[1])[p.argmax(1)]
        return {g: p[:, i + 1] for i, g in enumerate(self.names)}

    def out(self, arm, x, o):
        p = self.gate(x)
        res, used = np.zeros(len(x)), np.zeros(len(x))
        for g, R in self.readouts.items():
            if arm in R.th:
                res += p[g] * R.out(arm, x)
                used += p[g]
        return res + (1 - used) * o


def bm_score(out_all, j, y, col):
    pred = np.exp(out_all)
    pred[:, j] = np.exp(col)
    rs = np.array([np.corrcoef(y[:, c], pred[:, c])[0, 1] if y[:, c].std() > 1e-8 else np.nan for c in range(y.shape[1])])
    return rs


def run_readout(readout, obs, work, tag, soft=False):
    P0, names, _j = weights()
    j = names.index(OUTPUT[readout])
    b, t_ = np.load(BM_FIT), np.load(BM_TEST)
    xb, xv, yv, xt, yt = (unit(b["z_train"]).astype(np.float64), unit(b["z_val"]).astype(np.float64), b["adt_val"].astype(np.float64),
                          unit(t_["z"]).astype(np.float64), t_["adt"].astype(np.float64))
    head_all = lambda x: forward_np(P0, x)[1]
    cells = {}
    for g in GROUPS:                                         # every group enters the gate, whether or not it has this readout
        d = np.load(work / f"{g}.npz")
        cells[g] = dict(ids=d["ids"].astype(str), x=unit(d["Z"]).astype(np.float64), A=d["A"].astype(int), tg=d["unit"].astype(str),
                        y=d[f"y_{readout}"].astype(np.float64) if f"y_{readout}" in d.files else None)
        cells[g]["o"] = head_all(cells[g]["x"])[:, j]
        rng = np.random.default_rng(SEED)
        cells[g]["first"] = halves(np.where(cells[g]["A"] == 1, cells[g]["tg"], "control"), rng)
    groups = [g for g in GROUPS if readout in GROUPS[g]["readouts"]]
    # 1, 2: the grid on the head as trained, and admission
    grid, draws, adm = {}, {}, {}
    for g in groups:
        by, draws[g] = observed(obs, g, readout)
        grid[g] = {u: grid_row(r) for u, r in by.items()}
        adm[g] = admit(g, readout, grid[g])
        rna = obs / f"{g}_{readout}_rna" / "scope.json"
        if rna.exists():                                    # the same events read from RNA: is what Z leaves there before Z?
            rr = {r["event"]["unit"]: r for r in json.loads(rna.read_text())["results"]}
            dr = np.load(rna.parent / "scope_draws.npz")
            view = {}
            for u_ in grid[g]:
                S_ = rr[u_]["targets"]["shift"]["stages"]["RNA512"]
                zc = grid[g][u_]["scopes"]["own"]["class"]
                d_ = np.abs(draws[g][f"{u_}|shift|Z|own|{zc}"].astype(np.float64)) - np.abs(dr[f"{u_}|shift|RNA512|own|{S_['own']['class']}"].astype(np.float64))
                iv_ = interval(d_[1:])
                view[u_] = {"units": {k: S_["units"][k] for k in ("u", "ci", "verdict")}, "own": {k: S_["own"][k] for k in ("u", "ci", "verdict")},
                            "added_from_RNA_to_Z_own_scope": {"added": float(d_[0]), "ci": iv_, "reading": "loss" if iv_[0] > STEP else ("no added loss" if iv_[1] < STEP else "unresolved")}}
            adm[g]["input_view"] = view
            if adm[g]["rest"]:
                ow, sp_ = view["pool"]["own"]["verdict"], view["pool"]["added_from_RNA_to_Z_own_scope"]["reading"]
                adm[g]["rest"] = (f"the `units` scope at Z is {grid[g]['pool']['scopes']['units']['verdict']}. From RNA, a readout with the gene's own cells is {ow}; "
                                  f"from RNA to Z: {sp_}. " + ("Read as upstream: RNA carries it and Z does not." if (ow == "support" and sp_ == "loss") else
                                                                "Read as not readable from RNA by these readouts." if ow == "reject" else
                                                                "Not placed: neither view reaches a verdict."))
        assert np.allclose(own_output_u(cells[g]["y"], cells[g]["o"], cells[g]["A"], cells[g]["tg"])[0], draws[g]["pool|shift|system|frozen|system"], atol=2e-4), "the resamples do not pair"
    repaired = [g for g in groups if adm[g]["diagnostic"] == "ready"]
    genes_of = {g: [e["site"] for e in adm[g]["events"] if e["admission"]["status"] == "admitted"] for g in groups}
    # the repair of every other group as it is, fitted on all its cells (the arm "other")
    whole = {}
    for g in repaired:
        c = cells[g]
        whole[g] = Readout(c["x"], c["y"], c["o"], np.where(c["A"] == 0)[0], [(s, np.where((c["A"] == 1) & (c["tg"] == s))[0]) for s in genes_of[g]])
    # 4: two halves
    out = {g: {a: np.full(len(cells[g]["A"]), np.nan) for a in ARMS} for g in groups}
    inside = {g: np.full(len(cells[g]["A"]), np.nan) for g in groups}
    gates, fits, guard_val, guard_test = [], {g: [] for g in groups}, [], []
    for fold in (0, 1):
        fit = {g: (cells[g]["first"] if fold == 0 else ~cells[g]["first"]) for g in GROUPS}
        ctl = {g: cells[g]["x"][fit[g] & (cells[g]["A"] == 0)] for g in GROUPS}
        R = {}
        for g in repaired:
            c, og = cells[g], GROUPS[g]["other"]
            C = np.where(fit[g] & (c["A"] == 0))[0]
            ev = [(s, np.where(fit[g] & (c["A"] == 1) & (c["tg"] == s))[0]) for s in genes_of[g]]
            oth = (cells[og]["x"], cells[og]["y"], cells[og]["A"], cells[og]["tg"], genes_of[og]) if (og in groups and genes_of[og]) else None
            R[g] = Readout(c["x"], c["y"], c["o"], C, ev, oth)
            fits[g].append(R[g].log)
        sysm = System(xb, ctl, R, soft)
        for g in groups:
            c, rd = cells[g], ~fit[g]
            assert not set(c["ids"][rd]) & set(c["ids"][fit[g]])
            for a in ("controls", "rows_other", "rows", "cells"):
                if g in R and a in R[g].th:
                    out[g][a][rd] = sysm.out(a, c["x"][rd], c["o"][rd])
            og = GROUPS[g]["other"]
            if og in whole:
                out[g]["other"][rd] = whole[og].out("rows", c["x"][rd])
            if g in R:
                inside[g][fit[g]] = sysm.out("rows", c["x"][fit[g]], c["o"][fit[g]])
            p = sysm.gate(c["x"][rd])
            ko = c["A"][rd] == 1
            gates.append({"fold": fold, "group": g, **{f"on_{k}": float(v.mean()) for k, v in p.items()}, "own_on_control_cells": float(p[g][~ko].mean()),
                          "own_on_knockout_cells": float(p[g][ko].mean())})
        pb = sysm.gate(xt)
        gates.append({"fold": fold, "group": "bone marrow (test cells)", **{f"on_{k}": float(v.mean()) for k, v in pb.items()},
                      "cells_given_to_a_group": int(sum(v.sum() for v in pb.values())) if not soft else None})
        for store, (xq, yq) in ((guard_val, (xv, yv)), (guard_test, (xt, yt))):
            oq = head_all(xq)
            r0, r1 = bm_score(oq, j, yq, oq[:, j]), bm_score(oq, j, yq, sysm.out("rows", xq, oq[:, j]))
            store.append({"readout_as_trained": float(r0[j]), "readout_now": float(r1[j]), "fall": float(r0[j] - r1[j]), "mean_fall": float(np.nanmean(r0) - np.nanmean(r1)),
                          "worst_item_fall": float(np.nanmax(r0 - r1))})
    # 5: the head's own output again
    res = {"readout": readout, "output_of_the_head": OUTPUT[readout], "tag": tag, "margin": MARGIN, "step_margin": STEP, "groups": {}, "gate": gates}
    ledger = []
    for g in groups:
        c, a_, G = cells[g], adm[g], grid[g]
        run_obs, run_arm = f"scope:{g}_{readout}", f"round:{g}_{readout}:{tag}"
        ledger += [dict(e, type="source") for e in a_["events"]]
        ledger.append(LP.decision_event("admission", {**{k: v for k, v in a_.items() if k != "events"}, "readout": readout}, [run_obs]))
        u0, s0 = own_output_u(c["y"], c["o"], c["A"], c["tg"])
        ref = draws[g][f"pool|shift|Z|units|{G['pool']['scopes']['units']['class']}"].astype(np.float64)
        entry = {"admission": {k: v for k, v in a_.items() if k != "events"}, "events": a_["events"], "grid": G["pool"],
                 "as_trained": {"u": float(u0[0]), "ci": interval(u0[1:]), "verdict": verdict(interval(u0[1:])), "slope": s0}, "arms": {}, "fits": fits[g]}
        d0 = np.abs(u0) - np.abs(ref)
        entry["as_trained"]["beyond_units_at_Z"] = {"added": float(d0[0]), "ci": interval(d0[1:])}
        if g not in repaired:
            entry["closure"] = {"grade": "not a target", "why": "; ".join(a_["reasons"])}
            res["groups"][g] = entry
            continue
        for arm in ARMS:
            p = out[g][arm]
            if not np.isfinite(p).all():
                continue
            u, slope = own_output_u(c["y"], p, c["A"], c["tg"])
            gain, span = np.abs(u0) - np.abs(u), np.abs(u) - np.abs(ref)
            iv, si = interval(u[1:]), interval(span[1:])
            T_, C_ = c["A"] == 1, c["A"] == 0
            entry["arms"][arm] = {"u": float(u[0]), "ci": iv, "verdict": verdict(iv), "slope": slope,
                                  "valid": bool(slope > 0 and np.sign(p[T_].mean() - p[C_].mean()) == np.sign(c["y"][T_].mean() - c["y"][C_].mean())),
                                  "gain": float(gain[0]), "gain_ci": interval(gain[1:]),
                                  "beyond_units_at_Z": {"added": float(span[0]), "ci": si, "reading": "loss" if si[0] > STEP else ("no added loss" if si[1] < STEP else "unresolved")},
                                  "events": {}}
            for s in genes_of[g]:
                m = np.where(((c["A"] == 1) & (c["tg"] == s)) | (c["A"] == 0))[0]
                ue, _s = own_output_u(c["y"][m], p[m], c["A"][m], None)
                re = draws[g][f"{s}|shift|Z|units|{G[s]['scopes']['units']['class']}"].astype(np.float64)
                se = interval((np.abs(ue) - np.abs(re))[1:])
                entry["arms"][arm]["events"][s] = {"u": float(ue[0]), "ci": interval(ue[1:]), "verdict": verdict(interval(ue[1:])),
                                                   "beyond_units_at_Z": {"added": float(abs(ue[0]) - abs(re[0])), "ci": se,
                                                                         "reading": "loss" if se[0] > STEP else ("no added loss" if se[1] < STEP else "unresolved")}}
        # gate 3: the change itself
        ui, _s = own_output_u(c["y"], inside[g], c["A"], c["tg"])
        gv = {k: float(np.mean([q[k] for q in guard_val])) for k in guard_val[0]}
        gt = {k: float(np.mean([q[k] for q in guard_test])) for k in guard_test[0]}
        entry["gate3"] = {"closed_form": True, "on_the_cells_it_was_fitted_on": float(ui[0]), "objective_reached": bool(abs(ui[0]) < MARGIN), "no_audit_cell_in_a_fit": True,
                          "guard_on_validation_cells": {**gv, "kept": bool(max(q["fall"] for q in guard_val) <= GUARD and max(q["mean_fall"] for q in guard_val) <= GUARD)}}
        entry["gate6"] = {"guard_on_test_cells": {**gt, "kept": bool(max(q["fall"] for q in guard_test) <= GUARD and max(q["mean_fall"] for q in guard_test) <= GUARD)}}
        ok3 = entry["gate3"]["objective_reached"] and entry["gate3"]["guard_on_validation_cells"]["kept"]
        # closure
        r_ = entry["arms"]["rows"]
        base_loss = entry["as_trained"]["beyond_units_at_Z"]["ci"][0] > STEP
        if not ok3:
            grade = "not tested (gate 3)"
        elif r_["verdict"] == "support" and r_["valid"] and entry["as_trained"]["verdict"] == "reject":
            grade = "absolute"
        elif r_["beyond_units_at_Z"]["reading"] == "no added loss" and r_["valid"] and base_loss:
            grade = "relative"
        elif r_["beyond_units_at_Z"]["reading"] == "unresolved":
            grade = "loss no longer resolved"
        else:
            grade = "not closed"
        order = ["not closed", "loss no longer resolved", "relative", "absolute"]
        entry["closure"] = {"readout": readout, "span": ["Z", "own output"], "grade": grade, "attainable": a_["attainable"],
                            "reached_what_was_attainable": bool(grade in order and order.index(grade) >= order.index(a_["attainable"])),
                            "guard_kept": bool(entry["gate6"]["guard_on_test_cells"]["kept"])}
        # hypothesis and test
        ceil = G["pool"]["ceiling"]["coverage"]
        fo = draws[g][f"pool|shift|Z|other|{G['pool']['scopes']['other']['class']}"].astype(np.float64) if "other" in G["pool"]["scopes"] else None
        forecast = {"rows": ceil["ci"]}
        if fo is not None and "other" in entry["arms"]:
            forecast["other"] = interval((np.abs(draws[g]["pool|shift|system|frozen|system"]) - np.abs(fo))[1:])
        at_least = round(ceil["removes"] - STEP, 3)
        H = LP.hypothesis_event(f"coverage:{g}:{readout}:{tag}", {"event": f"{g}: the admitted knockouts, pooled", "readout": readout, "target": "shift",
                                                                 "quantity": "fall of |u| of the head's own output (frozen scope)"},
                                "coverage", forecast, {"own": "rows", "at_least": at_least, "controls": {"other": at_least} if "other" in forecast else {}},
                                [run_obs], rivals=["function", "upstream"], retrospective=(tag == "development"),
                                note="the forecast is the grid's ceiling; the rule asks the repair by rows to come within the step margin of it, and the repair of another group not to")
        arms = {k: {"gain": v["gain"], "ci": v["gain_ci"], "runs": [run_arm], "admitted": ok3} for k, v in entry["arms"].items() if k in ("rows", "other")}
        Tst = LP.test(H, arms)
        entry["hypothesis"], entry["test"] = H, Tst
        ledger += [LP.intervention_event(f"head+rows:{g}:{readout}", "head as trained", f"output {OUTPUT[readout]}", "replaced inside the gate of the group by a readout of Z updated by the admitted events",
                                         "located", held_fixed=["TEDDY", "the other 133 outputs", "the head outside the gate"], fitted_on=f"one half of {g}; read on the other"),
                   LP.intervention_event(f"head+other:{g}:{readout}", "head as trained", f"output {OUTPUT[readout]}", f"the repair of {GROUPS[g]['other']} applied to these cells", "control",
                                         fitted_on=GROUPS[g]["other"]),
                   H, Tst, LP.decision_event("closure", entry["closure"], [run_obs, run_arm])]
        res["groups"][g] = entry
    pr = obs / f"probe_{readout}" / "scope.json"
    if pr.exists():                                         # an event on the readout's own gene: typed, checked, and not admitted
        r_ = grid_row(json.loads(pr.read_text())["results"][0])
        nl = EV.no_leakage({"source": {"assign_from_guides": True}, "system": {"inputs": ["rna"]},
                            "readouts": [{"id": readout, "measurement": "protein", "features": READOUT_GENES[readout]}]}, PROBE)
        res["probe"] = {"site": PROBE, "group": "mel_ifn", "grid": r_, "no_leakage": nl,
                        "admission": {"status": "not admitted", "reason": "A2: " + nl["routes"]["site"]["reason"]},
                        "goes_to": "what to measure: no readout of Z or of RNA can be asked to carry a change made at the protein's own gene"}
        ledger.append(dict(type="source", event_id=f"probe:{PROBE}:{readout}", site=PROBE, channel="knockout", stratum=GROUPS["mel_ifn"]["stratum"], readout=readout,
                           amplitude=r_["response"]["value"], no_leakage=nl, admission=res["probe"]["admission"]))
    checks = LP.ledger_check(ledger)
    tests = [e for e in ledger if e.get("type") == "test"]
    res["ledger"] = {"events": len(ledger), "checks": checks, "errors": [c for c in checks if c["severity"] == "error"], "forecast_record": LP.forecast_record(tests)}
    res["cross_condition_pool"] = EV.closure([dict(e, pool="melanoma, both conditions") for g in groups if GROUPS[g]["cell_type"] == "melanoma" for e in adm[g]["events"]])
    tgt = [g for g in groups if g in repaired]
    res["loop"] = {"targets": tgt, "not_targets": [g for g in groups if g not in repaired],
                   "every_target_reached_its_attainable_grade": bool(tgt and all(res["groups"][g]["closure"]["reached_what_was_attainable"] for g in tgt)),
                   "guard_kept": bool(all(res["groups"][g]["closure"]["guard_kept"] for g in tgt)), "ledger_errors": len(res["ledger"]["errors"]),
                   "open": [{"group": g, "what": adm[g]["rest"]} for g in groups if adm[g]["rest"]]}
    return res, ledger


def show(res):
    f = lambda v: f"{v:6.3f}"
    print(f"\n=== {res['readout']} ({res['tag']})")
    for g, e in res["groups"].items():
        a, G = e["admission"], e["grid"]
        print(f"\n{g}: response {G['response']['value']:.2f} SD (snr {G['response']['snr']:.1f}) | grid: " + " ".join(f"{s} {v['u']:.2f} {v['verdict']}" for s, v in G["scopes"].items())
              + f" | reading {G['reading']} | {a['diagnostic']}, attainable {a['attainable']}" + ("" if not a["reasons"] else " | " + "; ".join(a["reasons"])))
        if not e["arms"]:
            continue
        print(f"  as trained {f(e['as_trained']['u'])} [{f(e['as_trained']['ci'][0])}, {f(e['as_trained']['ci'][1])}] {e['as_trained']['verdict']:10s} beyond Z {f(e['as_trained']['beyond_units_at_Z']['added'])}")
        for arm, v in e["arms"].items():
            print(f"  {arm:10s} {f(v['u'])} [{f(v['ci'][0])}, {f(v['ci'][1])}] {v['verdict']:10s} gain {f(v['gain'])} [{f(v['gain_ci'][0])}, {f(v['gain_ci'][1])}] | beyond Z {f(v['beyond_units_at_Z']['added'])} "
                  f"[{f(v['beyond_units_at_Z']['ci'][0])}, {f(v['beyond_units_at_Z']['ci'][1])}] {v['beyond_units_at_Z']['reading']:14s} | " + " ".join(f"{s}:{q['u']:.2f}" for s, q in v["events"].items()))
        print(f"  gate 3: on its own cells {f(e['gate3']['on_the_cells_it_was_fitted_on'])}, guard (validation) fall {e['gate3']['guard_on_validation_cells']['fall']:.4f} | gate 6: test fall {e['gate6']['guard_on_test_cells']['fall']:.4f}, "
              f"mean {e['gate6']['guard_on_test_cells']['mean_fall']:.5f}, worst item {e['gate6']['guard_on_test_cells']['worst_item_fall']:.4f}")
        print(f"  hypothesis: forecast {json.dumps({k: [round(x, 2) for x in v] for k, v in e['hypothesis']['forecast'].items()})}, at least {e['hypothesis']['rule']['at_least']:.2f} -> test {e['test']['outcome']} ({e['test']['why']})")
        print(f"  CLOSURE: {e['closure']['grade']} (attainable: {e['closure']['attainable']}; reached: {e['closure']['reached_what_was_attainable']})")
    print("\ngate:", *[json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in q.items()}) for q in res["gate"]], sep="\n  ")
    for g, e in res["groups"].items():
        if e["admission"].get("input_view"):
            v = e["admission"]["input_view"]["pool"]
            print(f"input view, {g}: from RNA units {v['units']['u']:.2f} {v['units']['verdict']}, own {v['own']['u']:.2f} {v['own']['verdict']}; RNA -> Z (own scope) {v['added_from_RNA_to_Z_own_scope']['added']:.3f} "
                  f"{[round(x, 3) for x in v['added_from_RNA_to_Z_own_scope']['ci']]} {v['added_from_RNA_to_Z_own_scope']['reading']}")
    print("ledger:", res["ledger"]["events"], "events;", [f"{c['code']} {c['severity']}" for c in res["ledger"]["checks"]], "| forecasts inside:", res["ledger"]["forecast_record"]["inside"], "of", res["ledger"]["forecast_record"]["n"])
    print("pool across conditions merges:", res["cross_condition_pool"])
    if res.get("probe"):
        q = res["probe"]["grid"]
        print("probe", res["probe"]["site"], "| response", round(q["response"]["value"], 2), "| " + " ".join(f"{k} {v['u']:.2f} {v['verdict']}" for k, v in q["scopes"].items()), "|", res["probe"]["admission"])
    print("LOOP:", json.dumps(res["loop"], ensure_ascii=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--readouts", nargs="+", required=True)
    ap.add_argument("--tag", default="development")
    ap.add_argument("--soft-gate", action="store_true", help="mix the readouts by the gate's probabilities (a sensitivity; writes round_<tag>_soft.json)")
    a = ap.parse_args()
    work = Path(os.environ["ANM_LOOP_WORK"])
    RES.mkdir(exist_ok=True)
    target = RES / (f"round_{a.tag}" + ("_soft" if a.soft_gate else "") + ".json")
    if a.tag != "development" and target.exists():
        raise SystemExit(f"{target.name} exists: this reading was made once and is not made again")
    allres, ledger = {}, []
    for r in a.readouts:
        allres[r], led = run_readout(r, a.observe, work, a.tag, a.soft_gate)
        ledger += led
    target.write_text(json.dumps(allres, indent=1, default=float))
    if not a.soft_gate:
        LP.write_ledger(ledger, RES / f"ledger_{a.tag}.jsonl")
    for r in a.readouts:
        show(allres[r])


if __name__ == "__main__":
    main()
