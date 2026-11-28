#!/usr/bin/env python3
"""The rehearsal, in one command: made-up worlds through summarize.py and analyze.py, and what came out against what was
planted. Worlds with planted knockouts and worlds with nothing planted. No real value is read.

  python rehearse.py --work DIR [--planted 11 12 13] [--null 21 22 23] [--teddy DIR_OF_A_PLANTED_WORLD/results]

-> results/rehearsal.json. --teddy: the pooled shares of teddy.py --rehearsal on a planted world, copied in.
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _ka as K  # noqa: E402


def world(work, seed, null):
    d = work / f"{'null' if null else 'planted'}_{seed}"
    if d.exists():
        shutil.rmtree(d)
    (d / "results").mkdir(parents=True)
    run = lambda *a: subprocess.run([sys.executable, *map(str, a)], check=True, capture_output=True, text=True)
    run(HERE / "rehearsal.py", "--out", d, "--seed", seed, *(["--null"] if null else []))
    run(HERE / "summarize.py", "--cells", d / "cells.tsv", "--library", d / "library.npz", "--counts", d / "counts_mel.npz", "--counts-ifn", d / "counts_ifn.npz", "--packs", d / "pack.npz", "--out", d / "summary.npz")
    run(HERE / "analyze.py", "--summary", d / "summary.npz", "--results", d / "results", "--rehearsal")
    R = json.loads((d / "results" / "reading.json").read_text())
    planted = json.loads((d / "planted.json").read_text())
    out = {"seed": seed, "planted": None if null else {"remove_the_signal": planted["remove_the_signal"], "one_gene": planted["one_gene"]}, "tests": {}}
    for t, r in R["tests"].items():
        truth = set() if null else set(planted["remove_the_signal"]) | ({"ONE1"} if t == "ifn" else set())
        ref = set(r["refused"])
        out["tests"][t] = {"knockouts": r["knockouts"], "w": r["w"], "refused": sorted(ref), "planted_and_refused": len(ref & truth), "planted": len(truth), "refused_and_adds": sorted(ref - truth),
                           "classes": [c["members"] for c in r["classes"]], "lines": {k: v["holds"] for k, v in r["lines"].items()}, "passed": r["passed"],
                           "whole_profile_stands_out_and_adds": sorted(set(r["whole_profile_stands_out"]) - truth),
                           "share_of_the_signal_removed": {x["knockout"]: x["share_of_the_signal_removed"]["rna"]["share"] for x in r["rows"] if x["knockout"] in truth}}
    for f in ("counts_mel.npz", "counts_ifn.npz", "library.npz", "states_mel.npz", "pack.npz"):
        (d / f).unlink()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--planted", type=int, nargs="*", default=[11, 12, 13])
    ap.add_argument("--null", type=int, nargs="*", default=[21, 22, 23])
    ap.add_argument("--teddy", type=Path)
    a = ap.parse_args()
    res = {"made_up_values_only": True, "worlds_with_planted_knockouts": [world(a.work, s, False) for s in a.planted], "worlds_with_nothing_planted": [world(a.work, s, True) for s in a.null]}
    if a.teddy:
        T = json.loads((a.teddy / "teddy.json").read_text())
        assert T["rehearsal_on_made_up_values"]
        res["teddy_on_a_planted_world"] = {"knockouts": T["knockouts"], "pairs": len(T["pairs"]), "pooled": T["pooled"]}
    (K.RUN / "results" / "rehearsal.json").write_text(json.dumps(res, indent=1) + "\n")
    for key in ("worlds_with_planted_knockouts", "worlds_with_nothing_planted"):
        for w in res[key]:
            for t, r in w["tests"].items():
                print(key.split("_")[2], w["seed"], t, "| planted and refused", r["planted_and_refused"], "of", r["planted"], "| refused and adds", r["refused_and_adds"], "| classes", r["classes"], "| lines", r["lines"], "| T stands out and adds", len(r["whole_profile_stands_out_and_adds"]))


if __name__ == "__main__":
    main()
