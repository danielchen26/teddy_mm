"""Read every head of this run, once, by the rule of DESIGN.md. Writes results/reading.json and prints the tables.

  python read.py

For every arm and data set:
  pooled   every knockout cell through the head fitted with its gene held out; the audit controls are split among the five heads as
           in runs/loop_20261112. u of the head's own PD-L1 output (anm_diag.scope, `frozen` scope), treated cells resampled by gene.
  events   one held-out gene, its head, all audit controls; cells resampled.
  fit      for comparison only: the same u on the events the head was fitted on (cells of the fit).
  plain    without the tool: the correlation of the own output with the measured protein on the cells of a held-out event, and the
           shift of the own output against the measured shift, each in units of the spread among the audit controls.
  guard    the head's own evaluation on the bone-marrow test cells.
The audit cells enter here and nowhere else.
"""
import json

import numpy as np

from _common import (ARMS, BM_TEST, CELLS, FALL, GUARD_HEADS, MARGIN, OUT, RES, STANDS_FOR, TARGETS, bm_score, forward_np, interval, load, own_output_u,
                     reading, unit, weights)

ALL = ["as_trained"] + ARMS
REF = "rows"


def main():
    P0, names, j = weights()
    t = np.load(BM_TEST); assert [str(n) for n in t["adt_names"]] == names
    xt, yt = unit(t["z"]), t["adt"].astype(np.float64)
    bm0 = bm_score(P0, xt, yt, names)
    res = {"margin": MARGIN, "guard": FALL, "reference_arm": REF, "stands_for": STANDS_FOR, "bone_marrow_as_trained": bm0, "data_sets": {}, "events": [], "family": {}}
    U_ev = {a: [] for a in ALL}
    for ds, genes in TARGETS.items():
        d = load(ds); A, tg, y, x = d["A"], d["tg"], d["y"], d["x"]
        heads = {a: {g: (P0 if a == "as_trained" else {k: v for k, v in np.load(OUT / "heads" / f"{ds}__{g}__{a}.npz").items()}) for g in genes} for a in ALL}
        out = {a: {g: forward_np(heads[a][g], x)[1][:, j] for g in genes} for a in ALL}
        use = (A == 1) | d["audit_ctl"]
        fitc = d["train_ctl"] & ~d["val"]
        entry = {"cells": CELLS[ds], "n_knockout": int((A == 1).sum()), "n_audit_controls": int(d["audit_ctl"].sum()), "arms": {}}
        draws = {}
        for a in ALL:
            pooled = np.zeros(len(A))
            for k, g in enumerate(genes):
                own = ((A == 1) & (tg == g)) | (d["ctl_group"] == k)
                pooled[own] = out[a][g][own]
            u, slope = own_output_u(y[use], np.exp(pooled[use]), A[use], tg[use])
            draws[a] = u
            iv = interval(u[1:])
            T, C = (A == 1), d["audit_ctl"]
            shift_pred, shift_y = float(np.exp(pooled[T]).mean() - np.exp(pooled[C]).mean()), float(y[T].mean() - y[C].mean())
            sane = bool(slope > 0 and np.sign(shift_pred) == np.sign(shift_y))
            bm = {g: bm_score(heads[a][g], xt, yt, names) for g in genes} if a != "as_trained" else {g: bm0 for g in genes}
            fall = float(np.mean([bm0["mean_pearson"] - bm[g]["mean_pearson"] for g in genes]))
            entry["arms"][a] = {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope_among_audit_controls": slope, "own_shift_has_the_measured_sign": sane,
                                "bone_marrow": {"mean_fall": fall, "lowest": float(min(bm[g]["mean_pearson"] for g in genes)), "highest": float(max(bm[g]["mean_pearson"] for g in genes)),
                                                "pdl1_pearson_mean": float(np.mean([bm[g]["pdl1_pearson"] for g in genes]))},
                                "guard_kept": bool(fall <= FALL), "events": {}}
            for g in genes:
                ug_use = ((A == 1) & (tg == g)) | d["audit_ctl"]
                ug, sg = own_output_u(y[ug_use], np.exp(out[a][g][ug_use]), A[ug_use], None)
                Tg, Cg = (A == 1) & (tg == g), d["audit_ctl"]
                o = out[a][g]
                fit_u = {}
                for g2 in genes:
                    if g2 == g and a != "as_trained":
                        continue
                    fu_use = ((A == 1) & (tg == g2) & ~d["val"]) | fitc
                    fit_u[g2] = float(own_output_u(y[fu_use], np.exp(o[fu_use]), A[fu_use], None)[0][0])
                entry["arms"][a]["events"][g] = {
                    "u": float(ug[0]), "ci": interval(ug[1:]), "reading": reading(interval(ug[1:])), "slope": sg, "n_knockout": int(Tg.sum()),
                    "plain": {"correlation_with_measured": float(np.corrcoef(np.exp(o[ug_use]), y[ug_use])[0, 1]),
                              "own_shift_sd": float((o[Tg].mean() - o[Cg].mean()) / o[Cg].std()), "measured_shift_sd": float((y[Tg].mean() - y[Cg].mean()) / y[Cg].std())},
                    "u_on_the_events_of_the_fit": fit_u}
                U_ev[a].append(abs(float(ug[0])))
        for a in ALL:                                             # every arm against the rows of the TEG, on the same resamples
            diff = np.abs(draws[a]) - np.abs(draws[REF])
            entry["arms"][a]["abs_u_minus_reference"] = {"difference": float(diff[0]), "ci": interval(diff[1:])}
        if GUARD_HEADS:                                           # the refit of arm `data` against the heads of runs/loop_guard_20261112
            from pathlib import Path
            worst = 0.0
            for g in genes:
                Pg = {k: v for k, v in np.load(Path(GUARD_HEADS) / f"{ds}__{g}__mlp_guarded.npz").items()}
                worst = max(worst, float(np.abs(forward_np(Pg, x)[1][:, j] - out["data"][g]).max()))
            entry["arm_data_against_the_heads_of_loop_guard_20261112"] = {"largest_difference_of_the_pdl1_output": worst}
        res["data_sets"][ds] = entry
    res["events"] = [f"{CELLS[ds]}|{g}" for ds, genes in TARGETS.items() for g in genes]
    rng = np.random.default_rng(2)
    idx = rng.integers(0, len(res["events"]), size=(2000, len(res["events"])))
    for a in ALL:                                                 # the ten held-out events as units
        Ua, Ur = np.array(U_ev[a]), np.array(U_ev[REF])
        dm = (Ua - Ur)[idx].mean(axis=1)
        res["family"][a] = {"mean_abs_u": float(Ua.mean()), "per_event_abs_u": [float(v) for v in Ua], "minus_reference": float((Ua - Ur).mean()), "ci": interval(dm),
                            "events_inside_the_margin": int((Ua < MARGIN).sum())}
    # the statements of the design
    well = {a: all(res["data_sets"][ds]["arms"][a]["reading"] == "supported" and res["data_sets"][ds]["arms"][a]["guard_kept"] and res["data_sets"][ds]["arms"][a]["own_shift_has_the_measured_sign"]
                   for ds in TARGETS) for a in ALL}
    res["trained_well"] = well
    res["statements"] = {
        "S1_the_rows_train_the_head_well": well[REF],
        "S1_none_of_class_edge_graph_judge_does": not any(well[a] for a in ("class", "edge1", "graph", "judge")),
        "S2_rows_lower_than_the_repair_already_made": {ds: res["data_sets"][ds]["arms"]["data"]["abs_u_minus_reference"] for ds in TARGETS},
    }
    RES.mkdir(exist_ok=True)
    (RES / "reading.json").write_text(json.dumps(res, indent=1))
    for ds in TARGETS:
        e = res["data_sets"][ds]
        print(f"\n{CELLS[ds]}: {e['n_knockout']} knockout cells of 5 genes, each through the head that never saw its gene; {e['n_audit_controls']} audit controls")
        print(f"{'arm':11s} {'u, own PD-L1 output':28s} {'reading':11s} {'slope':>7s} {'sign':>5s} {'|u| - rows':>24s} {'bone marrow fall':>17s}")
        for a in ALL:
            r = e["arms"][a]; dd = r["abs_u_minus_reference"]
            print(f"{a:11s} {r['u']:6.3f} [{r['ci'][0]:6.3f}, {r['ci'][1]:6.3f}]      {r['reading']:11s} {r['slope_among_audit_controls']:7.2f} {'ok' if r['own_shift_has_the_measured_sign'] else 'NO':>5s} "
                  f"{dd['difference']:7.3f} [{dd['ci'][0]:6.3f}, {dd['ci'][1]:6.3f}] {r['bone_marrow']['mean_fall']:10.4f} {'kept' if r['guard_kept'] else 'BROKEN'}")
    print("\nten held-out events:", {a: (round(v["mean_abs_u"], 3), v["events_inside_the_margin"]) for a, v in res["family"].items()})
    print("trained well, by the rule:", well)
    print("statements:", json.dumps({k: v for k, v in res["statements"].items() if not isinstance(v, dict)}))


if __name__ == "__main__":
    main()
