#!/usr/bin/env python3
"""Gate 5, before F is read: the second hypothesis and the declaration.

  python declare_final.py --heads HEADS_DIR --final-arrays F_cite_arrays.npz --final-z F_z_rna.npy --pbmc-arrays ext_cite_arrays.npz

Writes ../results/gate5_declaration.json: the two hypotheses, the arms and their heads, the donors and the proteins of
F (from ../results/final_selection.json), how the outcome is stated, and the sha256 of every input of
scripts/read_final.py. Writes ../results/ledger_final.jsonl with the second hypothesis. Reads no protein value of F:
of the pack it reads the donor and status of every cell and which proteins are measured.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import common
import plan
import read_final
from anm_diag import loop

BOUND = 0.03419135502693498          # the rule of donors_hao_coverage


def main():
    p = argparse.ArgumentParser()
    for k in ("heads", "final-arrays", "final-z", "pbmc-arrays"):
        p.add_argument(f"--{k}", type=Path, required=True)
    a = p.parse_args()
    res, here = common.RUN / "results", Path(__file__).resolve().parent
    if (res / "final.json").exists():
        raise SystemExit("F has been read; the declaration is not written again")
    sel = json.loads((res / "final_selection.json").read_text())
    dev = json.loads((res / "development.json").read_text())
    h1 = next(e for e in map(json.loads, open(common.RUN.parent / "scope_20261116" / "results" / "ledger.jsonl")) if e.get("id") == "donors_hao_coverage")
    if h1["rule"]["at_least"] != BOUND:
        raise SystemExit("the bound of donors_hao_coverage is not the one of this round")

    with np.load(a.final_arrays, allow_pickle=False) as z:
        names = [str(x) for x in z["adt_names"]]
        present, donors, status = z["adt_present"].astype(bool), z["donors"].astype(str), z["status"].astype(str)
    with np.load(a.pbmc_arrays, allow_pickle=False) as z:
        in_pbmc = z["adt_present"].astype(bool)
    proteins = [n for n, f, q in zip(names, present, in_pbmc) if f and q]
    if proteins != sel["proteins"]["among_the_112"]:
        raise SystemExit("the proteins of the pack are not those of the selection")
    listed = {d["donor"]: d for d in sel["donors"]}
    if sorted(set(donors.tolist())) != sorted(listed):
        raise SystemExit("the donors of the pack are not those of the selection")
    for d in listed.values():
        me = donors == d["donor"]
        if int(me.sum()) != sel["rule"]["cells_per_donor"] or set(status[me].tolist()) != {d["status"]}:
            raise SystemExit(f"{d['donor']}: the pack does not hold the cells of the selection")

    g = dev["gains_over_as_is"]
    h2 = loop.hypothesis_event(
        "coverage_not_levels",
        {"event": "a donor the head was not trained on, PBMC of a study the head has not seen", "readout": f"{len(proteins)} proteins", "target": "variation",
         "quantity": "fall of u (variation) at the head's own output, mean over the donors of the final study"},
        "coverage", {"coverage": g["coverage_not_seen"]["ci"], "labels": g["labels"]["ci"]},
        {"own": "coverage", "at_least": BOUND, "controls": {"labels": BOUND}}, ["round1_coverage_20261116/development"], rivals=["levels"],
        note="what the added PBMC cells teach is the relation between RNA and protein, not the protein levels of their study. Control: the same cells with the "
             "proteins shuffled among the cells of each donor. The forecasts are the intervals read in development on the 8 PBMC donors (coverage: the donor's "
             "batch not in the fit), other units than those of the final study.")

    heads = {f"s{s}_{read_final.HEAD_OF[k]}": k for s in plan.SEEDS for k in read_final.ARMS}
    reused = json.loads((res / "gate5_heads.json").read_text())["heads"]
    sha = {"final_arrays": read_final.sha256(a.final_arrays), "final_z": read_final.sha256(a.final_z)}
    for n in heads:
        if not (a.heads / n / "DONE").exists():
            raise SystemExit(f"{n} is not finished")
        sha[f"head:{n}"] = read_final.sha256(a.heads / n / "head.npz")
        if n in reused and reused[n] != sha[f"head:{n}"]:
            raise SystemExit(f"{n}: not the weights named in gate5_heads.json")
    for s in read_final.SCRIPTS:
        sha[f"script:{s}"] = read_final.sha256(here / s)
    decl = {
        "what": "round 1, gate 5: what is read on the final study, fixed before any protein value of it is read",
        "study": sel["name"], "run": read_final.RUN, "hypotheses": [h1, h2],
        "arms": {k: [n for n, kk in heads.items() if kk == k] for k in read_final.ARMS},
        "donors": [{k: d[k] for k in ("donor", "status", "centre", "sample")} for d in sel["donors"]],
        "cells": {"per_donor": sel["rule"]["cells_per_donor"], "n": int(len(donors)), "sha256_of_rows": sel["cells"]["sha256_of_rows"]},
        "proteins": proteins,
        "proteins_by_exact_name": [r["ours"] for r in sel["proteins"]["mapping"] if r["rule"] == "exact name" and r["ours"] in proteins],
        "decision": read_final.DECISION,
        "outcome": "the outcome over the three seeds together, stated only if each seed alone reads the same; otherwise unresolved (GATE5.md, part 1)",
        "sha256": sha,
    }
    decl["digest"] = hashlib.sha256(json.dumps(decl, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    (res / "gate5_declaration.json").write_text(json.dumps(decl, indent=1) + "\n")
    loop.write_ledger([h2], res / "ledger_final.jsonl")
    print(f"declared: {len(decl['donors'])} donors ({sum(d['status'] == 'healthy' for d in decl['donors'])} healthy), {len(proteins)} proteins, "
          f"{len(heads)} heads; digest {decl['digest'][:16]}")
    print("coverage_not_levels forecasts:", {k: [round(x, 4) for x in v] for k, v in h2["forecast"].items()})


if __name__ == "__main__":
    main()
