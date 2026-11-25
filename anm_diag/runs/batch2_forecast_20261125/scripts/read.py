#!/usr/bin/env python3
"""The reading of batch 2. results/reading.json. Made once, after forecasts.json was committed.

  python read.py --cells batch2_cells.tsv --states STATES_DIR --read read_cells.npz --pack batch2_pack.npz
                 --counts-frangieh counts_frangieh2021.npz --counts-papalexi counts_papalexi2021.npz

1. Every pair of the batch: the measured change in spreads of the readout among the stratum's control cells, its error,
   and where it falls against the committed forecast.
2. The account (anm_diag.account, as run 41) of every pair of the batch, and of every recorded pair of batch 1, whose
   change is at least four errors: at the input view V, the state Z and, in THP-1 and co-culture, the typed state T and
   genes at random R.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

from _b2 import AC, EDGES, RUN, RUNS, SEED, STRATA, cycle, stratum, table

sys.path.insert(0, str(RUNS / "account_20261125" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
sys.path.insert(0, str(RUNS / "outer_step_20261122" / "scripts"))
from _common import forward_np, unit, weights  # noqa: E402
from _read import account  # noqa: E402
from build import OUTPUT  # noqa: E402
from build2 import lognorm  # noqa: E402

BATCH1 = {"mel_co": ["CD44", "TAPBP", "CDH19"], "mel_ctl": ["B2M", "KLF4"], "thp1": ["BRD4", "CMTM6", "POU2F2"]}
TYPED = {"mel_co": "mel_co", "thp1": "thp1"}                       # the strata with the sets of run 40, and their manifests
PAIRS = (("Z", "T"), ("R", "T"), ("Z", "V"), ("T", "V"))
SEED_SVD = 20261125
csr = lambda c: sp.csr_matrix((c["vals"].astype(np.float32), (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]))


def halves(g, c):
    """The halves: the development cells keep those of run 40 (the first cells of the stratum in the list, in their
    order; checked against the manifest of the typed states), the other cells are split by cycle.halves."""
    first = np.zeros(len(c["ids"]), bool)
    n0 = 0
    if g in TYPED:
        man = json.loads((RUNS / "typed_repair_20261125" / "results" / "states" / f"{TYPED[g]}_manifest.json").read_text())
        n0 = man["cells"]
        first[:n0] = cycle.halves(c["unit"][:n0], SEED)
        assert hashlib.sha256(first[:n0].tobytes()).hexdigest()[:16] == man["sha256_of_the_halves"], "the halves of the development cells are not those of run 40"
    first[n0:] = cycle.halves(c["unit"][n0:], SEED)
    return first


def main():
    ap = argparse.ArgumentParser()
    for k in ("cells", "states", "read", "pack", "counts-frangieh", "counts-papalexi"):
        ap.add_argument("--" + k, type=Path, required=True)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    out_file = RUN / "results" / ("reading.json" if not a.only else f"reading_{a.only}.json")
    if out_file.exists():
        raise SystemExit(f"{out_file.name} exists: the reading was made once and is not made again")
    fc = json.loads((RUN / "results" / "forecasts.json").read_text())
    assert fc["cells_file_sha256"] == hashlib.sha256(a.cells.read_bytes()).hexdigest(), "the forecasts were made for another list of cells"
    cells = table(a.cells)
    P, old = np.load(a.pack), np.load(a.read)
    op = {i: k for k, i in enumerate(old["ids"].astype(str))}
    checks = {}
    for st, rs in (("frangieh", ("PDL1", "MHC1", "HLAE")), ("papalexi", ("PDL1",))):
        kk = np.array([op[i] for i in P[f"check_{st}_ids"].astype(str)])
        for r in rs:
            checks[f"{st} {r}"] = float(np.abs(P[f"check_{st}_{r}"] - old[f"readout__{r}"][kk]).max())
    assert max(checks.values()) < 1e-9, checks
    print("check: protein values of cells read before are identical", checks, flush=True)
    counts = {"frangieh2021": a.counts_frangieh, "papalexi2021": a.counts_papalexi}
    P0, names, _j = weights()
    res = {"checks": checks, "min_snr": AC.MIN_SNR, "strata": {}}
    for g, (study, cond, readouts) in STRATA.items():
        if a.only and g != a.only:
            continue
        c = stratum(g, cells, a.states, read=a.read, new=a.pack)
        C = c["unit"] == "control"
        new = sorted(set(c["unit"][c["role"] == "event"].tolist()))
        recorded = sorted(set(c["unit"][c["role"] == "recorded"].tolist()))
        fcs = {(r["event"], r["readout"]): r for r in fc["strata"][g]["forecasts"]}
        S = {"pairs": [], "accounts": {}}
        # --- 1. every pair against its forecast
        for r in readouts:
            y = c[f"y_{r}"]
            sd = y[C].std()
            for e in new:
                m = c["unit"] == e
                f = fcs[(e.replace("new:", ""), r)]
                change = float((y[m].mean() - y[C].mean()) / sd)
                se = float(np.sqrt(y[m].var(ddof=1) / m.sum() + y[C].var(ddof=1) / C.sum()) / sd)
                S["pairs"].append({"event": f["event"], "readout": r, "role": f["role"], "edge": f["edge"], "cells": int(m.sum()), "change": change, "se": se, "snr": abs(change) / se,
                                   "forecast": f["forecast"], "interval": f["interval"], "sd": f["sd"], "z": (change - f["forecast"]) / f["sd"],
                                   "inside": bool(f["interval"][0] <= change <= f["interval"][1]), "read": f["read"], "from_the_recorded_events": f["from_the_recorded_events"],
                                   "by_its_edge": f["by_its_edge"], "without_its_edge": f["read"] + f["from_the_recorded_events"]})
        for p in S["pairs"]:
            print(f"{g:7s} {p['event']:9s} -> {p['readout']:5s} change {p['change']:+.2f} +- {p['se']:.2f} | forecast {p['forecast']:+.2f} [{p['interval'][0]:+.2f}, {p['interval'][1]:+.2f}] "
                  f"{'inside' if p['inside'] else 'OUTSIDE'} (z {p['z']:+.1f})", flush=True)
        # --- 2. the accounts
        cn = np.load(counts[study])
        cp = {i: k for k, i in enumerate(cn["ids"].astype(str))}
        R = lognorm(csr(cn)[np.array([cp[i] for i in c["ids"]])].tocsr(), np.where(C)[0])
        V = TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED_SVD).fit_transform(R).astype(np.float64)
        first = halves(g, c)
        stages = {"Z": c["Z"], "V": V}
        if g in TYPED:
            nz = np.linalg.norm(c["raw"], axis=1, keepdims=True) + 1e-6
            blk = lambda name, f: c["blocks"][:, c["sets"].index(f"{name}|{f}")].astype(np.float64)
            two = lambda name: tuple(np.hstack([blk(name, f), c["raw"] - blk(name, f)]) / nz for f in (0, 1))
            stages.update(T=two("A"), R=two("R"))
        S["stages"] = {k: int((v[0] if isinstance(v, tuple) else v).shape[1]) for k, v in stages.items()}
        S["cells"] = {"control": int(C.sum()), "recorded": int((c["role"] == "recorded").sum()), "of_the_batch": int((c["role"] == "event").sum()), "genes_of_the_input_view": int(R.shape[1])}
        for r in readouts:
            y = c[f"y_{r}"]
            o = forward_np(P0, unit(c["raw"]))[1][:, names.index(OUTPUT[r])].astype(np.float64)
            cc = {"y": y, "unit": c["unit"], "first": first, "o": o}
            sd = y[C].std()
            for e in new + [h for h in BATCH1[g] if h in recorded]:
                m = c["unit"] == e
                snr = abs(y[m].mean() - y[C].mean()) / np.sqrt(y[m].var(ddof=1) / m.sum() + y[C].var(ddof=1) / C.sum())
                if snr < AC.MIN_SNR:
                    continue
                key = f"{e.replace('new:', '')}|{r}"
                acc = account(stages, cc, e, [h for h in recorded if h != e], pairs=PAIRS, say=lambda s_: print(f"{g:7s} {r:5s} {s_}", flush=True))
                acc.update(of_the_batch=e.startswith("new:"), own_gene=(e.replace("new:", ""), r) in EDGES, change_in_spreads=float((y[m].mean() - y[C].mean()) / sd), snr=float(snr))
                S["accounts"][key] = acc
                res["strata"][g] = S
                out_file.with_suffix(".partial").write_text(json.dumps(res, indent=1, default=float) + "\n")
        res["strata"][g] = S
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    if out_file.with_suffix(".partial").exists():
        out_file.with_suffix(".partial").unlink()


if __name__ == "__main__":
    main()
