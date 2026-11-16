#!/usr/bin/env python3
"""Fit the heads of the round one after another, in the order of plan.py, skipping those already finished.

  python run_all.py --bm-arrays cite_arrays.npz --bm-z z_rna.npy --pbmc-arrays ext_cite_arrays.npz --pbmc-z ext_z_rna.npy \\
      --out OUT [--seeds 0 1 2] [--final] [--time-budget-sec 6600] [--threads 2]

Exit code 75: the time budget ended; run the same command again and it goes on where it stopped, with the same result.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plan  # noqa: E402
import train_arm  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    for k in ("bm-arrays", "bm-z", "pbmc-arrays", "pbmc-z", "out"):
        p.add_argument(f"--{k}", type=Path, required=True)
    p.add_argument("--seeds", type=int, nargs="+", default=list(plan.SEEDS))
    p.add_argument("--final", action="store_true", help="the heads gate 5 adds, not the development heads")
    p.add_argument("--time-budget-sec", type=float, default=6600)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args()
    t0 = time.time()
    say = lambda *x: print(time.strftime("%H:%M:%S"), *x, flush=True)
    bm, pb = train_arm.load_data(a.bm_arrays, a.bm_z, a.pbmc_arrays, a.pbmc_z)
    todo = plan.heads(seeds=a.seeds, n_train=int((bm["split"] == "train").sum()), final=a.final)
    left = [c for c in todo if not (a.out / c["name"] / "DONE").exists()]
    say(f"{len(todo)} heads, {len(todo) - len(left)} finished, {len(left)} to fit")
    for c in left:
        budget = a.time_budget_sec - (time.time() - t0)
        if budget < 60:
            say("budget ended")
            sys.exit(train_arm.EXIT_BUDGET)
        say(f"{c['name']}: fitting")
        if not train_arm.fit(c, bm, pb, a.out / c["name"], time_budget=budget, threads=a.threads, log=say):
            say(f"{c['name']}: budget ended, state kept")
            sys.exit(train_arm.EXIT_BUDGET)
        say(f"{c['name']}: finished")
    say("all finished")


if __name__ == "__main__":
    main()
