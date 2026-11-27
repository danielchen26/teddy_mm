#!/usr/bin/env python3
"""Batch 3 again, by the candidate for version 6. DEVELOPMENT: the candidate was made from this batch. results/replay.json.

  python replay.py --cells batch3_cells.tsv --states B2_STATES_frangieh.npz B3_STATES_frangieh.npz --states-thp B2_STATES_papalexi.npz
                   --read read_cells.npz --pack2 batch2_pack.npz --library LIBRARY.npz --pack-a batch3_pack_A.npz --pack-b batch3_pack_B.npz

`as it was`: the 698 pairs of runs/batch3_frozen_20261126 forecast again from stage A alone, exactly as forecast.py
did, with the candidate in place of version 5; then scored against the measured changes of reading.json.
`with teachers of no type`: what a batch designed by the candidate has. In every stratum eight units of no type
(seeded) are recorded on the readouts never read as well, as the teachers were, and the other pairs are forecast.
Their values come from stage B, so this layout forecasts fewer pairs.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
B3 = RUN.parent / "batch3_frozen_20261126"
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(B3 / "scripts"))
import _v6 as V6  # noqa: E402
from _v6 import V  # noqa: E402
import _b3 as B  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

PATHWAY, TEACHERS, SEED = "interferon pathway", 8, 20261127
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None
key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"


def stats(rows, f="forecast"):
    if not rows:
        return None
    z2 = np.array([r["z"] ** 2 for r in rows])
    return {"pairs": len(rows), "inside": int(sum(r["inside"] for r in rows)), "share_inside": float(np.mean([r["inside"] for r in rows])), "mean_z2": float(z2.mean()), "median_z2": float(np.median(z2)),
            "error": rms([r["change"] - r[f] for r in rows]), "error_of_the_model_alone": rms([r["change"] - r["read"] for r in rows]), "no_change": rms([r["change"] for r in rows])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--states-thp", type=Path, nargs="+", required=True)
    for k in ("read", "pack2", "library", "pack-a", "pack-b"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    fc = json.loads((B3 / "results" / "forecasts.json").read_text())
    R = json.loads((B3 / "results" / "reading.json").read_text())
    measured, scored5 = {key(r): (r["change"], r["se"]) for r in R["rows"]}, {key(r): r for r in R["rows"]}
    cells = B.table(a.cells)
    st = {"frangieh2021": B.States(a.states), "papalexi2021": B.States(a.states_thp)}
    used = B.used_values(a.read, a.pack2)
    S = {g: B.stratum(g, cells, st[B.STRATA[g][0]]) for g in B.STRATA}
    d = np.load(a.library)
    lib_ids, un = d["ids"].astype(str), np.where(d["A"] == 0, "control", d["target"].astype(str))
    mel_all, mel_new, mel_used = list(B.PANEL["frangieh2021"]), B.NEW["frangieh2021"], B.USED["frangieh2021"]
    scale = fc["scale"]
    sc_of = lambda g: {r: scale[r][g] for r in mel_all}
    sites = sorted((set(un.tolist()) | {u for g in ("mel_co", "mel_ctl") for u in S[g]["unit"].tolist()}) - {"control"})
    _declared, tab = V.declared(sites)
    r43 = json.loads((RUN.parent / "one_loop_20261125" / "results" / "retro.json").read_text())["edges"]
    admitted = [e for e in r43 if e["status"] == "admitted"]
    tab_of = {g: {**tab, **{(e["site"], e["readout"]): e["class"] for e in admitted if e["stratum"] == g}} for g in ("mel_co", "mel_ctl")}

    packs_a = B.Packs([a.pack_a])
    YL, _log, _bg = packs_a.readouts("frangieh2021", lib_ids)
    ZL = V.unit_length(d["Z"].astype(np.float64))
    lib = V6.library_law(ZL, YL, un, tab, sc_of("library"))
    print(f"library: spread with no edge {lib['law'].sigma0:.3f} (version 5: {lib['before'].sigma0:.3f}) | widths {lib['width']} | parts known by declaration {len(lib['known'])}", flush=True)
    flags5, flags6 = lib["before"].flags(), lib["law"].flags()
    res = {"note": "development: the candidate was made from this batch; nothing here is a test",
           "library": {"spread_with_no_edge": {"version 5": lib["before"].sigma0, "candidate": lib["law"].sigma0}, "width": lib["width"], "pairs_with_a_known_part": len(lib["known"]),
                       "flags": {"version 5": [f"{f['event']}->{f['readout']}" for f in flags5], "candidate": [f"{f['event']}->{f['readout']}" for f in flags6]}}, "layouts": {}}

    def melanoma(packs, teachers_of_no_type):
        """The two strata of the melanoma screen: {stratum: rows}. teachers_of_no_type: {stratum: [units]} recorded on every readout."""
        out, first, spec = {}, {}, {}
        for g in ("mel_co", "mel_ctl"):
            c = S[g]
            extra = np.isin(c["unit"], teachers_of_no_type.get(g, []))
            rec_u = c["kind"] != "new"
            Yu = {r: np.array([used[r][i] for i in c["ids"][rec_u]]) for r in mel_used}
            rec_n = np.isin(c["kind"], ["control", "teacher"]) | extra
            test_n = np.isin(c["kind"], ["read", "new"]) & ~extra
            Yn, _l, _b = packs.readouts("frangieh2021", c["ids"][rec_n], mel_new)
            xn, unn, s0 = c["Z"][rec_n], c["unit"][rec_n], None
            if not (set(unn.tolist()) - {"control"}):
                xn, Yn, unn = B.with_placeholder(xn, Yn, unn)
                s0 = lib["law"].sigma0
            spec[g] = {"used": {"x": c["Z"][rec_u], "Y": Yu, "unit": c["unit"][rec_u], "table": V.restrict(tab_of[g], mel_used), "scales": sc_of(g)},
                       "new": {"x": xn, "Y": Yn, "unit": unn, "table": V.restrict(tab, mel_new), "scales": sc_of(g), "sigma0": s0}, "test": test_n}

        def fit(sp, known=None):
            ed, sc = EV.by_scale(sp["table"], sp["scales"])
            law = AC.Law(sp["x"], sp["Y"], sp["unit"], edges=ed, types=V.PATHWAY, scale=sc, known=known)
            return (law.fit(sigma0=sp["sigma0"]) if sp.get("sigma0") is not None else law.fit(fallback=lib["law"].sigma0)), ed, sc

        first = {g: {k: fit(spec[g][k]) for k in ("used", "new")} for g in spec}
        twice = V.measured_twice(lib, [first[g]["used"] for g in first])
        for g in spec:
            c = S[g]
            new, test = c["kind"] == "new", spec[g]["test"]
            # the direct part of every edge of the stratum: as the edge forecasts it, and as measured where the pair is recorded
            direct = {}
            direct.update(V6.edge_parts(V6.forecast(*first[g]["used"], c["Z"][new], c["unit"][new], lib, twice)))
            direct.update(V6.edge_parts(V6.forecast(*first[g]["new"], c["Z"][test], c["unit"][test], lib, twice)))
            for k in ("used", "new"):
                direct.update(V6.edge_parts(first[g][k][0].forecast()))
            known = EV.derived(direct, V6.spreads(first[g]["used"][0], first[g]["new"][0]), V6.OVER)
            laws = {k: fit(spec[g][k], V6.values(known)) for k in ("used", "new")}
            rows = V6.forecast(*laws["used"], c["Z"][new], c["unit"][new], lib, twice, known=known) + V6.forecast(*laws["new"], c["Z"][test], c["unit"][test], lib, twice, known=known)
            out[g] = {"rows": [r for r in rows if r["event"] != B.NOTHING], "spread_with_no_edge": {k: laws[k][0].sigma0 for k in laws}, "spread_from": {k: laws[k][0].spread_from for k in laws},
                      "pairs_with_a_known_part": sum(1 for (e, _r) in known if e in set(c["unit"].tolist()))}
        return out, twice

    def thp(packs, teachers_of_no_type):
        c = S["thp1"]
        y = np.array([used["PDL1"][i] for i in c["ids"]])
        law0 = AC.Law(c["Z"], {"PDL1": y}, c["unit"], types=V.PATHWAY).fit()
        thp_lib = {"law": law0, "classes": {}, "shares": {}, "width": law0.conformal()}
        extra = np.isin(c["unit"], teachers_of_no_type.get("thp1", []))
        rec, test = np.isin(c["kind"], ["control", "teacher"]) | extra, (c["kind"] == "read") & ~extra
        Y, _l, _b = packs.readouts("papalexi2021", c["ids"][rec], B.NEW["papalexi2021"])
        law, ed, scl = V.stratum_law(c["Z"][rec], Y, c["unit"][rec], {}, {}, thp_lib)
        return {"rows": V6.forecast(law, ed, scl, c["Z"][test], c["unit"][test], thp_lib, []), "spread_with_no_edge": {"new": law.sigma0}, "spread_from": {"new": law.spread_from}}

    rng = np.random.default_rng(SEED)
    chosen = {}
    for g in B.STRATA:                                              # eight units of no type per stratum, among those read before: recorded on every readout
        c = S[g]
        plain = sorted({u for u, k in zip(c["unit"].tolist(), c["kind"].tolist()) if k == "read" and u not in V.PATHWAY})
        chosen[g] = sorted(rng.choice(plain, TEACHERS, replace=False).tolist())
    for layout, teachers, packs in (("as it was", {}, packs_a), ("with teachers of no type", chosen, B.Packs([a.pack_a, a.pack_b]))):
        mel, _twice = melanoma(packs, teachers)
        strata = {**mel, "thp1": thp(packs, teachers)}
        rows = []
        for g, v in strata.items():
            for r in v["rows"]:
                r["stratum"] = g
                if key(r) not in measured:
                    continue
                ch, se = measured[key(r)]
                o = scored5[key(r)]
                r.update(change=ch, se=se, z=(ch - r["forecast"]) / r["sd"], inside=bool(r["interval"][0] <= ch <= r["interval"][1]), group=o["group"], kind=o["kind"],
                         version_5={"forecast": o["S2"]["forecast"], "inside": o["S2"]["inside"], "z": o["S2"]["z"]})
                rows.append(r)
        v5 = [{"change": r["change"], "forecast": r["version_5"]["forecast"], "read": r["read"], "z": r["version_5"]["z"], "inside": r["version_5"]["inside"]} for r in rows]
        sets = {"all pairs": lambda r: True,
                "pairs on which the systems can differ": lambda r: key(r) in set(fc["sets"]["the_systems_can_differ"]),
                "own-gene pairs (H3)": lambda r: key(r) in set(fc["sets"]["H3"]),
                "of the declared type, no edge (H4)": lambda r: key(r) in set(fc["sets"]["H4"]),
                "no type, no edge, never read (H5)": lambda r: key(r) in set(fc["sets"]["H5"]),
                "with a change of at least four errors": lambda r: abs(r["change"]) / r["se"] >= 4}
        out = {"teachers_of_no_type": teachers, "pairs": len(rows), "by_set": {}, "by_stratum": {}, "strata": {g: {k: v[k] for k in v if k != "rows"} for g, v in strata.items()},
               "width_from": {}, "outside": []}
        for name, f in sets.items():
            m = [f(r) for r in rows]
            out["by_set"][name] = {"candidate": stats([r for r, k in zip(rows, m) if k]), "version 5": stats([r for r, k in zip(v5, m) if k])}
        for g in B.STRATA:
            m = [r["stratum"] == g for r in rows]
            out["by_stratum"][g] = {"candidate": stats([r for r, k in zip(rows, m) if k]), "version 5": stats([r for r, k in zip(v5, m) if k])}
        for r in rows:
            if r["edge"] is None:
                out["width_from"][r["width_from"]] = out["width_from"].get(r["width_from"], 0) + 1
        out["outside"] = [{"pair": key(r), "edge": r["edge"], "type": r["type"], "change": r["change"], "forecast": r["forecast"], "interval": r["interval"], "z": r["z"], "known_by_declaration": r["known_by_declaration"],
                           "inside_version_5": r["version_5"]["inside"]} for r in sorted(rows, key=lambda r: -abs(r["z"])) if not r["inside"]]
        s = out["by_set"]["all pairs"]
        out["the_gate_as_written_for_batch_3"] = {"H1_share_inside_at_least_0.90": s["candidate"]["share_inside"] >= 0.90, "H1_mean_z2_in_0.5_to_2": 0.5 <= s["candidate"]["mean_z2"] <= 2,
                                                  "median_z2": s["candidate"]["median_z2"], "error_smaller_than_version_5": s["candidate"]["error"] < s["version 5"]["error"]}
        res["layouts"][layout] = out
        print(f"\n{layout}: {len(rows)} pairs", flush=True)
        for name, v in out["by_set"].items():
            c6, c5 = v["candidate"], v["version 5"]
            if c6:
                print(f"  {name:40s} n={c6['pairs']:3d} | inside {c5['inside']:3d} -> {c6['inside']:3d} ({c6['share_inside']:.3f}) | mean z2 {c5['mean_z2']:.2f} -> {c6['mean_z2']:.2f} | median {c5['median_z2']:.2f} -> {c6['median_z2']:.2f} "
                      f"| error {c5['error']:.3f} -> {c6['error']:.3f} (model alone {c6['error_of_the_model_alone']:.3f})")
        for g, v in out["by_stratum"].items():
            if v["candidate"]:
                print(f"  {g:40s} n={v['candidate']['pairs']:3d} | inside {v['version 5']['inside']:3d} -> {v['candidate']['inside']:3d} ({v['candidate']['share_inside']:.3f}) | median z2 {v['version 5']['median_z2']:.2f} -> {v['candidate']['median_z2']:.2f}")
        print("  widths:", out["width_from"])
        print("  spreads:", {g: {k: round(x, 3) for k, x in v["spread_with_no_edge"].items()} for g, v in out["strata"].items()})
        print("  outside:", len(out["outside"]), [(o["pair"], round(o["z"], 1)) for o in out["outside"][:12]])
    (RUN / "results" / "replay.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
