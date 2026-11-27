#!/usr/bin/env python3
"""A planted twin with the real measurement: version 5 and the candidate for version 6 on it. results/twin.json.

  python twin.py --library LIBRARY.npz --cells batch3_cells.tsv --states B2_STATES_frangieh.npz B3_STATES_frangieh.npz --pack-a batch3_pack_A.npz

Why: batch 3 found faults of the method that the planted system of run 43 could not show, because its readouts were
not measured as the real ones are. Here the cells are real control cells (their states and their 24 antibody columns,
all opened before), so the readouts are centred log-ratios with the real noise and the real background; the events
are made up, and what each one does is planted in the log counts, so the answer is known:

- own-gene units: the unit removes a share k of how far its protein stands above the background (ten proteins that
  stand well above it). The share of a pair differs a little between the two conditions.
- units of one type: their cells' states move along the direction a real pathway knockout moves them, and three
  proteins change in proportion (what a readout of control cells cannot know).
- every pair with no edge has a small direct part of its own, wider in the stratum than in the library, and one pair
  in 25 has a large one: an edge the graph does not hold.

A library (control cells under interferon) and a stratum (control cells of the co-culture) are drawn, in two layouts:
`with teachers of no type` (the stratum has recorded units of every kind, as batch 4 would) and `as in batch 3` (its
only recorded units are three of the type). A method is fit to forecast a batch if it holds here; version 5 should
fail here as it failed on batch 3.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(RUN.parent / "batch3_frozen_20261126" / "scripts"))
import _v6 as V6  # noqa: E402
from _v6 import V  # noqa: E402
import _b3 as B  # noqa: E402
from anm_diag import events as EV  # noqa: E402

OWN = ["CD119", "PDL1", "CD29", "CD44", "CD47", "CD49f", "CD58", "CD59", "MHC1", "CD9"]      # the proteins of the own-gene units
TYPE = "planted type"
BETA = {"MHC1": -0.6, "PDL1": -0.5, "CD119": 0.15}             # what a unit of the type does to three proteins, in log counts, at amplitude one
CELLS, SPREAD = 200, {"library": 0.02, "stratum": 0.05}         # cells per unit; the spread of a direct part with no edge, in spreads of the readout
HIDDEN, CARRIED = 0.04, 0.06                                    # the share of pairs with an edge the graph does not hold; how far a share moves between conditions
LAYOUTS = {"with teachers of no type": {"recorded": [f"O{i:02d}" for i in range(1, 5)] + ["P1", "P2", "P3"] + [f"N{i:02d}" for i in range(1, 9)],
                                        "new": [f"O{i:02d}" for i in range(5, 11)] + ["P4", "P5"] + [f"N{i:02d}" for i in range(9, 21)]},
           "as in batch 3": {"recorded": ["P1", "P2", "P3"], "new": [f"O{i:02d}" for i in range(1, 11)] + ["P4", "P5"] + [f"N{i:02d}" for i in range(1, 21)]}}
GENES = [f"O{i:02d}" for i in range(1, 11)] + [f"P{i}" for i in range(1, 6)] + [f"N{i:02d}" for i in range(1, 26)]
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None


def condition(rng, Z, L, panel, units, n_control, plan, name, direction):
    """The cells of one condition: control cells as they are, and the units with what is planted. Returns the states, the
    readouts, the unit of every cell, the scale of every readout, and what was planted per pair in spreads of the readout."""
    col = {r: panel.index(ab) for r, (ab, _g) in V.MELANOMA.items()}
    iso = [panel.index(k) for k in V.ISOTYPES]
    order = rng.permutation(len(Z))
    ctl, rest = order[:n_control], order[n_control:]
    clr_c = L[ctl] - L[ctl].mean(1, keepdims=True)
    sd = {r: float(clr_c[:, c].std()) for r, c in col.items()}
    excess = {r: float((L[ctl][:, c] - L[ctl][:, iso].mean(1)).mean()) for r, c in col.items()}
    xs, Ls, un, planted = [Z[ctl]], [L[ctl]], ["control"] * n_control, {}
    for u in units:
        idx = rng.choice(rest, CELLS, replace=False)
        x, y = Z[idx].copy(), L[idx].copy()
        p = plan[u]
        for r in V.MELANOMA:
            d = rng.normal(0, SPREAD[name]) * sd[r]                 # a small direct part of its own, in log counts
            kind = "nothing"
            if p.get("own") == r:
                d, kind = -p["share"][name] * excess[r], "its edge"
            elif (u, r) in p["hidden"][name]:
                d, kind = p["hidden"][name][(u, r)] * sd[r], "an edge the graph does not hold"
            if p.get("amp"):
                d += p["amp"] * BETA.get(r, 0.0)
            y[:, col[r]] += d
            planted[(u, r)] = kind
        if p.get("amp"):
            x = V.unit_length(x + p["amp"] * direction)
        xs.append(x); Ls.append(y); un += [u] * CELLS
    x, Lall, un = np.vstack(xs), np.vstack(Ls), np.array(un)
    clr = Lall - Lall.mean(1, keepdims=True)
    bg, C = Lall[:, iso].mean(1), un == "control"
    Y = {r: clr[:, c] for r, c in col.items()}
    scales = {r: V.scale(Lall[C][:, c], bg[C], clr[C][:, c]) for r, c in col.items()}
    return x, Y, un, scales, planted


def world(seed, pools, panel, direction, amp0):
    rng = np.random.default_rng([20261126, seed])
    plan = {}
    for i, u in enumerate(GENES):
        p = {"hidden": {"library": {}, "stratum": {}}}
        if u.startswith("O"):
            k = rng.uniform(0.2, 0.8)
            p.update(own=OWN[i], share={"library": k, "stratum": float(np.clip(k + rng.normal(0, CARRIED), 0.05, 0.95))})
        if u.startswith("P"):
            p["amp"] = float(amp0 * rng.uniform(0.6, 1.4))
        for name in ("library", "stratum"):
            for r in V.MELANOMA:
                if p.get("own") != r and rng.random() < HIDDEN:
                    p["hidden"][name][(u, r)] = float(rng.choice([-1, 1]) * rng.uniform(0.3, 1.5))
        plan[u] = p
    table = {(u, p["own"]): "own" for u, p in plan.items() if "own" in p}
    types = {u: TYPE for u in GENES if u.startswith("P")}
    xl, Yl, ul, scl, _pl = condition(rng, *pools["library"], panel, GENES, 1500, plan, "library", direction)
    out = {}
    for layout, lay in LAYOUTS.items():
        xs, Ys, us, scs, planted = condition(rng, *pools["stratum"], panel, lay["recorded"] + lay["new"], 2000, plan, "stratum", direction)
        rec, new = np.isin(us, ["control"] + lay["recorded"]), np.isin(us, lay["new"])
        Yr = {r: y[rec] for r, y in Ys.items()}
        rows = {}
        # --- version 5
        lib5 = V.library_law(xl, Yl, ul, table, scl, types)
        law, ed, sc = V.stratum_law(xs[rec], Yr, us[rec], table, scs, lib5, types)
        twice = V.measured_twice(lib5, [(law, ed, sc)])
        rows["version 5"] = V.scored(V.forecast(law, ed, sc, xs[new], us[new], lib5, twice, types), law, Ys, us)
        spread5 = law.sigma0
        # --- the candidate for version 6
        lib6 = V6.library_law(xl, Yl, ul, table, scl, types=types)
        laws, _known, _first = V6.stratum_laws({"all": {"x": xs[rec], "Y": Yr, "unit": us[rec], "table": table, "scales": scs}}, lib6, types=types)
        law, ed, sc = laws["all"]
        twice = V.measured_twice(lib6, [(law, ed, sc)])
        first = V6.forecast(law, ed, sc, xs[new], us[new], lib6, twice, types=types)
        known = EV.derived(V6.edge_parts(first), V6.spreads(law), V6.OVER)
        rows["version 6, candidate"] = V.scored(V6.forecast(law, ed, sc, xs[new], us[new], lib6, twice, known=known, types=types), law, Ys, us)
        for k in rows:
            for r in rows[k]:
                r["planted"] = planted[(r["event"], r["readout"])]
                r["unit_kind"] = "own-gene" if r["event"].startswith("O") else "of the type" if r["event"].startswith("P") else "plain"
        out[layout] = {"rows": rows, "spread_with_no_edge": {"version 5": spread5, "version 6, candidate": law.sigma0, "planted": SPREAD["stratum"]},
                       "width_of_version_6": sorted({r["width_from"] for r in rows["version 6, candidate"] if r["edge"] is None})}
    return out


def summary(rows):
    def one(rs):
        if not rs:
            return None
        z2 = np.array([r["z"] ** 2 for r in rs])
        return {"pairs": len(rs), "inside": int(sum(r["inside"] for r in rs)), "share_inside": float(np.mean([r["inside"] for r in rs])), "mean_z2": float(z2.mean()), "median_z2": float(np.median(z2)),
                "error": rms([r["change"] - r["forecast"] for r in rs]), "error_of_the_control_readout": rms([r["change"] - r["read"] for r in rs])}
    sel = {"all pairs": lambda r: True,
           "with its edge": lambda r: r["planted"] == "its edge",
           "no edge, of an own-gene unit": lambda r: r["planted"] == "nothing" and r["unit_kind"] == "own-gene",
           "no edge, of a plain unit": lambda r: r["planted"] == "nothing" and r["unit_kind"] == "plain",
           "no edge, of a unit of the type": lambda r: r["planted"] == "nothing" and r["unit_kind"] == "of the type",
           "an edge the graph does not hold": lambda r: r["planted"] == "an edge the graph does not hold",
           "all but the edges the graph does not hold": lambda r: r["planted"] != "an edge the graph does not hold"}
    return {k: one([r for r in rows if f(r)]) for k, f in sel.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--pack-a", type=Path, required=True)
    ap.add_argument("--worlds", type=int, default=8)
    a = ap.parse_args()
    packs = B.Packs([a.pack_a])
    panel = packs.panel["frangieh2021"]
    logs = lambda ids: np.array([packs.log["frangieh2021"][i] for i in ids])
    d = np.load(a.library)
    c0 = d["A"] == 0
    cells = B.table(a.cells)
    st = B.States(a.states)
    m = (cells["study"] == "frangieh2021") & (cells["condition"] == "Co-culture")
    ids, role, unit = cells["cell"][m], cells["role"][m], cells["target"][m]
    Zs = V.unit_length(st.rows(ids))
    ctl = role == "control"
    pools = {"library": (V.unit_length(d["Z"][c0].astype(np.float64)), logs(d["ids"].astype(str)[c0])), "stratum": (Zs[ctl], logs(ids[ctl]))}
    # how a real knockout of the pathway moves the state: the mean shift of the cells of its recorded units in the co-culture
    shifts = np.array([Zs[unit == g].mean(0) - Zs[ctl].mean(0) for g in V.PATHWAY if (unit == g).sum() >= 30])
    direction, amp0 = shifts.mean(0) / np.linalg.norm(shifts.mean(0)), float(np.median(np.linalg.norm(shifts, axis=1)))
    print(f"pools: library {len(pools['library'][0])} control cells, stratum {len(pools['stratum'][0])}; a pathway knockout moves the state by {amp0:.3f} (median of {len(shifts)})", flush=True)
    res = {"planted": {"cells_per_unit": CELLS, "spread_with_no_edge": SPREAD, "share_of_pairs_with_an_edge_the_graph_does_not_hold": HIDDEN, "how_far_a_share_moves_between_conditions": CARRIED,
                       "the_type_moves_the_state_by": amp0, "the_type_changes": BETA, "own_gene_proteins": OWN, "layouts": LAYOUTS}, "worlds": a.worlds, "layouts": {}}
    pooled = {lay: {"version 5": [], "version 6, candidate": []} for lay in LAYOUTS}
    per_world = {lay: [] for lay in LAYOUTS}
    for seed in range(a.worlds):
        w = world(seed, pools, panel, direction, amp0)
        for lay, v in w.items():
            for k, rows in v["rows"].items():
                pooled[lay][k] += rows
            per_world[lay].append({"world": seed, "spread_with_no_edge": v["spread_with_no_edge"], "width_of_version_6": v["width_of_version_6"],
                                   **{k: {"share_inside": float(np.mean([r["inside"] for r in rows])), "median_z2": float(np.median([r["z"] ** 2 for r in rows]))} for k, rows in v["rows"].items()}})
            print(f"world {seed} | {lay:24s} | " + " | ".join(f"{k}: inside {np.mean([r['inside'] for r in rows]):.3f}, median z2 {np.median([r['z'] ** 2 for r in rows]):.2f}" for k, rows in v["rows"].items())
                  + f" | spread {v['spread_with_no_edge']['version 5']:.3f} / {v['spread_with_no_edge']['version 6, candidate']:.3f}", flush=True)
    for lay in LAYOUTS:
        res["layouts"][lay] = {"by_world": per_world[lay], **{k: summary(rows) for k, rows in pooled[lay].items()}}
    (RUN / "results" / "twin.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    for lay in LAYOUTS:
        print("\n" + lay)
        for k in ("version 5", "version 6, candidate"):
            for name, s in res["layouts"][lay][k].items():
                if s:
                    print(f"  {k:22s} {name:42s} n={s['pairs']:4d} inside {s['share_inside']:.3f} mean z2 {s['mean_z2']:6.2f} median z2 {s['median_z2']:.2f} error {s['error']:.3f} (control readout {s['error_of_the_control_readout']:.3f})")


if __name__ == "__main__":
    main()
