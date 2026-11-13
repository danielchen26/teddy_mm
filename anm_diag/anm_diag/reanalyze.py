"""Re-read a finished run under the current verdict rules (or other margins) without recomputing fits.

  python -m anm_diag.reanalyze RUN_DIR [--margin 0.10] [--step-margin 0.10]

Uses anm_record.json and, when present, anm_bootstrap.npz (for step intervals between the chosen classes).
Writes anm_record_reread.json and anm_report_reread.md next to the originals. Re-reading under a margin
chosen after seeing the results is exploratory by definition; the record notes the margins used.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from .report import write_report
from .stats import ci
from .verdicts import root_cause, verdict_union


def reread(record, boots=None, margin=None, step_margin=None):
    rec = copy.deepcopy(record)
    d = rec["declaration"]
    margin = float(margin if margin is not None else d["decision"]["margin_u"])
    step_margin = float(step_margin if step_margin is not None else d["decision"].get("step_margin", margin))
    alpha = float(d["decision"]["alpha"])
    chain = d["chain"]
    nodes = {n["id"]: n for n in d["nodes"]}
    for rid, r in rec["readouts"].items():
        if not r["states"]:          # no resolved response: nothing to re-read
            continue
        for sid, s in r["states"].items():
            v, best = verdict_union(s["by_class"], margin)
            s.update({"verdict": v, "best_class": best, "u": s["by_class"][best]["u"], "ci": s["by_class"][best]["ci"]})
        pairs = list(zip(chain[:-1], chain[1:])) + [(nodes[s]["from"], s) for s in r["states"]
                                                     if nodes[s]["role"] == "probe" and nodes[s]["from"] in r["states"]]
        for a, b in pairs:
            ba, bb = r["states"][a]["best_class"], r["states"][b]["best_class"]
            key = f"{a}->{b}"
            if boots is not None and f"{rid}|{a}/{ba}" in boots:
                r["steps"][key]["best_class"] = ci(boots[f"{rid}|{b}/{bb}"] - boots[f"{rid}|{a}/{ba}"], alpha)
            elif ba == bb:
                r["steps"][key]["best_class"] = r["steps"][key]["same_class"][ba]
            else:
                r["steps"][key]["best_class"] = [float("nan"), float("nan")]
                r["steps"][key]["note"] = "classes differ and no bootstrap file: see same-class steps"
        pcs = [p for p in r["states"] if nodes[p]["role"] == "positive_control"]
        pc_v = r["states"][pcs[0]]["verdict"] if pcs else None
        steps = {k: v["best_class"] for k, v in r["steps"].items() if not np.isnan(v["best_class"][0])}
        r["root_cause"] = root_cause(r["delta_ci"], pc_v, chain, {c: r["states"][c]["verdict"] for c in chain},
                                     step_cis=steps, step_margin=step_margin)
    rec["reread"] = {"margin_u": margin, "step_margin": step_margin, "rule": "union of declared readout classes"}
    rec["declaration"]["decision"]["margin_u"] = margin
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir"); ap.add_argument("--margin", type=float); ap.add_argument("--step-margin", type=float)
    a = ap.parse_args(argv)
    p = Path(a.run_dir)
    rec = json.loads((p / "anm_record.json").read_text())
    bf = p / "anm_bootstrap.npz"
    boots = dict(np.load(bf)) if bf.exists() else None
    out = reread(rec, boots, a.margin, a.step_margin)
    (p / "anm_record_reread.json").write_text(json.dumps(out, indent=1, default=float))
    write_report(out, p / "anm_report_reread.md")
    for rid, r in out["readouts"].items():
        print(f"{rid}: {r['root_cause']['cause']} — {r['root_cause'].get('detail', '')}")


if __name__ == "__main__":
    main()
