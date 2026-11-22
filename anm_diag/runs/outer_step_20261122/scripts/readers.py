#!/usr/bin/env python3
"""Richer readers of the input view (lever V, step 0): no new data. For a target: every stage and reader of the package's
reader specification, the gap from the input view to Z, the positive control, and the permuted control.

  python readers.py --target "MHC1 mel_co" --z OBS/mel_co_MHC1_z2 --in512 OBS/mel_co_MHC1_in512 --in4k OBS/mel_co_MHC1_in4k \
                    --perm OBS/mel_co_MHC1_perm --out ../results/readers.json

Each folder holds the output of `anm_diag scope` (scope.json, scope_draws.npz) for one view of the same cells; the
resamples of the tool are the same for every view, so gaps are paired. Results are added to the output file by target.

  stage Z        the embedding the head reads: ridge, boosting
  stage RNA512   512 components of RNA: ridge, boosting
  stage RNA4K    the 4,000 genes with the largest variance among the control cells: ridge
  PERM           RNA512 with its rows permuted across cells: ridge (the control: it must read nothing)

For the input view the reader of an event is the best of the readers of RNA512 and RNA4K (smallest |u|, own-cells scope).
dominated        the gap from the input view to Z has its interval below zero
positive control some reader class at some stage supports some event of the family, or their pool
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
from anm_diag import scope as SC  # noqa: E402
from _common import interval  # noqa: E402

MARGIN, STEP = 0.10, 0.10
verdict = lambda iv: SC._verdict(iv, MARGIN, True)


def load(folder):
    d = json.loads((folder / "scope.json").read_text())
    return {r["event"]["unit"]: r["targets"]["shift"]["stages"] for r in d["results"]}, np.load(folder / "scope_draws.npz")


def classes(v):
    """Every reader class of a scope: {class: {u, ci, verdict}}."""
    by = v.get("by_class") or {v["class"]: v}
    return {c: {"u": x["u"], "ci": x["ci"], "verdict": verdict(x["ci"])} for c, x in by.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    for k in ("z", "in512", "in4k", "perm"):
        ap.add_argument(f"--{k}", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    S, D = {}, {}
    for stage, folder in (("Z", a.z), ("RNA512", a.in512), ("RNA4K", a.in4k), ("PERM", a.perm)):
        by, D[stage] = load(folder)
        S[stage] = {u: T[stage] for u, T in by.items()}
    events = list(S["Z"])
    res = {"stages": {"Z": ["ridge", "hgb"], "RNA512": ["ridge", "hgb"], "RNA4K": ["ridge"], "PERM": ["ridge"]}, "events": {}}
    seen = []
    for u in events:
        e = {"by_stage": {st: {sc: classes(S[st][u][sc]) for sc in ("units", "own")} for st in S}}
        seen += [x["verdict"] for st in ("Z", "RNA512", "RNA4K") for sc in ("units", "own") for x in e["by_stage"][st][sc].values()]
        best = {}
        for st in ("Z", "RNA512", "RNA4K"):
            c = min(e["by_stage"][st]["own"], key=lambda k: abs(e["by_stage"][st]["own"][k]["u"]))
            best[st] = (c, D[st][f"{u}|shift|{st}|own|{c}"].astype(np.float64))
        st_in = min(("RNA512", "RNA4K"), key=lambda k: abs(best[k][1][0]))
        zc, zd = best["Z"]
        ic, idr = best[st_in]
        gap = np.abs(zd) - np.abs(idr)                                  # what is added from the input view to Z; below zero: Z reads better
        iv = interval(gap[1:])
        e["Z_own"] = {"reader": zc, **e["by_stage"]["Z"]["own"][zc]}
        e["input_own"] = {"view": st_in, "reader": ic, **e["by_stage"][st_in]["own"][ic]}
        e["from_the_input_view_to_Z"] = {"added": float(gap[0]), "ci": iv, "reading": "loss" if iv[0] > STEP else ("no added loss" if iv[1] < STEP else "unresolved"), "dominated": bool(iv[1] < 0)}
        for st in ("RNA512", "RNA4K"):                                  # each view against Z, for the record
            g = np.abs(zd) - np.abs(best[st][1])
            e.setdefault("each_view_to_Z", {})[st] = {"reader": best[st][0], "added": float(g[0]), "ci": interval(g[1:])}
        if "hgb" in e["by_stage"]["Z"]["own"]:
            for st in ("Z", "RNA512"):
                g = np.abs(D[st][f"{u}|shift|{st}|own|ridge"].astype(np.float64)) - np.abs(D[st][f"{u}|shift|{st}|own|hgb"].astype(np.float64))
                e.setdefault("boosting_reads_better_than_ridge_by", {})[st] = {"value": float(g[0]), "ci": interval(g[1:])}
        e["permuted"] = e["by_stage"]["PERM"]["own"]["ridge"]
        res["events"][u] = e
    p = res["events"]["pool"]
    res["positive_control"] = {"holds": "support" in seen, "readings": len(seen), "supported": seen.count("support"), "unresolved": seen.count("unresolved"), "rejected": seen.count("reject")}
    res["permuted_control"] = {"holds": bool(all(abs(e["permuted"]["u"]) >= 0.9 for e in res["events"].values())), "values": {u: e["permuted"]["u"] for u, e in res["events"].items()}}
    res["input_view"] = {"dominated": p["from_the_input_view_to_Z"]["dominated"], "can_be_read": bool(not p["from_the_input_view_to_Z"]["dominated"] and res["positive_control"]["holds"] and res["permuted_control"]["holds"]),
                         "own": p["input_own"], "gap_to_later": {k: p["from_the_input_view_to_Z"][k] for k in ("added", "ci", "reading")}, "classes": ["ridge", "hgb"], "nonlinear": True}
    out = json.loads(a.out.read_text()) if a.out.exists() else {}
    out[a.target] = res
    a.out.write_text(json.dumps(out, indent=1, default=float) + "\n")
    f = lambda x: f"{x['u']:.2f} [{x['ci'][0]:.2f}, {x['ci'][1]:.2f}] {x['verdict']}"
    print(a.target)
    for u, e in res["events"].items():
        print(f"  {u:7s} " + " | ".join(f"{st}: " + ", ".join(f"{c} {f(x)}" for c, x in e["by_stage"][st]["own"].items()) for st in ("Z", "RNA512", "RNA4K")) + f" | PERM {e['permuted']['u']:.2f}")
        g = e["from_the_input_view_to_Z"]
        print(f"          input ({e['input_own']['view']}, {e['input_own']['reader']}) -> Z ({e['Z_own']['reader']}): {g['added']:.3f} [{g['ci'][0]:.3f}, {g['ci'][1]:.3f}] {g['reading']}{' DOMINATED' if g['dominated'] else ''}"
              + (" | boosting - ridge: " + ", ".join(f"{st} {v['value']:.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}]" for st, v in e.get("boosting_reads_better_than_ridge_by", {}).items())))
    print("  positive control:", res["positive_control"], "| permuted control:", res["permuted_control"]["holds"], "| input view can be read:", res["input_view"]["can_be_read"])


if __name__ == "__main__":
    main()
