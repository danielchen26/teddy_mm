#!/usr/bin/env python3
"""Batch 1: admission on the selection third, and the arrays and SPEC files of the audit. No fit except the components of RNA.

  python build.py --states batch1_states.npz --pack batch1_pack.npz --counts-frangieh COUNTS.npz --counts-papalexi COUNTS.npz --work WORK

Environment as runs/teg_loop_20261120 (ANM_LOOP_WORK, ANM_FEATURES, ANM_TEST_CELLS). Writes results/selection.json and, in WORK,
for every stratum `b1_<stratum>.npz` (the cells outside the selection third: record, Z, readouts, the head's own outputs) and
`b1_<stratum>_in.npz` (the same cells: three views of their RNA), and one SPEC per view for every readout with an admitted event.
The arrays of the earlier rounds (thp1.npz, mel_co.npz, ...) are linked into WORK under their own names and are only read.

Strata: thp1 (the new units of the THP-1 screen with the control cells that have states), mel_co (the new units and the new
control cells of co-culture), mel_ctl (those of the control condition). Admission (A1, A2) is `anm_diag.cycle.select` on the cells
whose SHA-256 of '20261123|<cell name>' is divisible by 3; those cells are written to no array.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
RUN, RUNS = HERE.parent, HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
sys.path.insert(0, str(RUNS / "outer_step_20261122" / "scripts"))
from anm_diag import cycle  # noqa: E402
from _common import FEATURES, forward_np, unit, weights  # noqa: E402
from build import OUTPUT  # noqa: E402
from build2 import BOTH, RIDGE, lognorm, spec, views, write  # noqa: E402

SEED = "20261123"
ROLES = json.loads((RUN / "results" / "roles.json").read_text())
STRATA = {"thp1": ("papalexi2021", "all", "mel_ifn.npz"), "mel_co": ("frangieh2021", "Co-culture", None), "mel_ctl": ("frangieh2021", "Control", None)}
OTHER = {"PDL1": "mel_ifn.npz", "MHC1": "mel_ifn_mhc.npz", "HLAE": None}
OWN = {r: [g.replace("-", "") for g in v] + list(v) for r, v in ROLES["own_gene_of_the_readout"].items()}
csr = lambda c: sp.csr_matrix((c["vals"].astype(np.float32), (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--pack", type=Path, required=True)
    ap.add_argument("--counts-frangieh", type=Path, required=True)
    ap.add_argument("--counts-papalexi", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    a = ap.parse_args()
    out_file = RUN / "results" / "selection.json"
    if out_file.exists():
        raise SystemExit("selection.json exists: admission was made once and is not made again")
    a.work.mkdir(parents=True, exist_ok=True)
    src, src2 = Path(os.environ["ANM_LOOP_WORK"]), Path(os.environ["ANM_STEP_WORK"])
    for f, d in (("bm_train.npz", src), ("mel_ifn.npz", src), ("mel_co.npz", src), ("thp1.npz", src), ("mel_ifn_mhc.npz", src2)):
        if not (a.work / f).exists():
            os.symlink(d / f, a.work / f)
    P, S = np.load(a.pack), np.load(a.states)
    sid = {i: k for k, i in enumerate(S["ids"].astype(str))}
    pid = P["ids"].astype(str)
    assert set(pid) == set(sid), "the pack and the states are not of the same cells"
    Z = S["z_last"][np.array([sid[i] for i in pid])].astype(np.float32)
    # the protein values of cells earlier runs read are the ones those runs used
    old = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "coculture.npz")
    op = {i: k for k, i in enumerate(old["ids"].astype(str))}
    kk = np.array([op[i] for i in P["check_frangieh_ids"].astype(str)])
    checks = {f"frangieh2021 {r}": float(np.abs(P[f"check_frangieh_{r}"] - old[f"readout__{r}"][kk]).max()) for r in ("PDL1", "MHC1", "HLAE")}
    oldp = np.load(FEATURES / "papalexi2021" / "dense.npz")
    op = {i: k for k, i in enumerate(oldp["ids"].astype(str))}
    kk = np.array([op[i] for i in P["check_papalexi_ids"].astype(str)])
    checks["papalexi2021 PDL1"] = float(np.abs(P["check_papalexi_PDL1"] - oldp["readout__PDL1"][kk]).max())
    assert max(checks.values()) < 1e-9, checks
    print("check: protein values of cells already used are identical", checks)
    P0, names, _j = weights()
    head = lambda z, r: forward_np(P0, unit(z))[1][:, names.index(OUTPUT[r])].astype(np.float64)
    study, cond, un = P["study"].astype(str), P["condition"].astype(str), P["unit"].astype(str)
    counts = {"frangieh2021": np.load(a.counts_frangieh), "papalexi2021": np.load(a.counts_papalexi)}
    res = {"rule": f"a cell is in the selection third if sha256('{SEED}|<cell name>') is divisible by 3", "min_snr": cycle.MIN_SNR, "checks": checks, "strata": {}}
    for g, (st, c, _o) in STRATA.items():
        m = (study == st) & (cond == c)
        K = dict(ids=pid[m], unit=un[m], Z=Z[m], **{f"y_{r}": P[f"readout__{r}"][m].astype(np.float64) for r in ROLES["strata"][g]["readouts"]})
        if g == "thp1":                                          # its control cells have states already
            T = np.load(src / "thp1.npz")
            ctl = T["unit"].astype(str) == "control"
            assert ctl.sum() == ROLES["control_cells"]["thp1"]["cells"], "the control cells of the earlier rounds are not all here"
            K = dict(ids=np.r_[K["ids"], T["ids"].astype(str)[ctl]], unit=np.r_[K["unit"], T["unit"].astype(str)[ctl]], Z=np.vstack([K["Z"], T["Z"][ctl]]),
                     y_PDL1=np.r_[K["y_PDL1"], T["y_PDL1"][ctl].astype(np.float64)])
        K["A"] = (K["unit"] != "control").astype(int)
        e = {"cells": {k: int(v) for k, v in zip(*np.unique(K["unit"], return_counts=True))}, "readouts": {}}
        share = None
        for r in ROLES["strata"][g]["readouts"]:
            sel = cycle.select(K["ids"], K["unit"], K[f"y_{r}"], SEED, own=OWN[r])
            share = sel["share"]
            e["readouts"][r] = {"pass_A1": sel["pass_A1"], "refused_A2": [x for x in sel["pass_A1"] if x not in sel["admitted"]], "admitted": sel["admitted"], "table": sel["table"]}
            K[f"frozen_{r}"] = head(K["Z"], r)
        e["cells_in_the_selection_third"] = {k: int(v) for k, v in zip(*np.unique(K["unit"][share], return_counts=True))}
        keep = ~share
        A = {k: v[keep] for k, v in K.items()}
        e["cells_read"] = {k: int(v) for k, v in zip(*np.unique(A["unit"], return_counts=True))}
        assert not (a.work / f"b1_{g}.npz").is_symlink()
        np.savez(a.work / f"b1_{g}.npz", **A)
        # the input view of the cells that are read
        cn = counts[st]
        cp = {i: k for k, i in enumerate(cn["ids"].astype(str))}
        X = csr(cn)
        ctl = np.where(A["A"] == 0)[0]
        if g == "thp1":                                          # control cells: the RNA the earlier rounds read; new units: the screen's RNA file, same genes
            genes = json.loads((FEATURES / "papalexi2021" / "raw_genes.json").read_text())
            R0 = sp.load_npz(FEATURES / "papalexi2021" / "node__RAW.npz").tocsr()
            p0 = {i: k for k, i in enumerate(oldp["ids"].astype(str))}
            cg = {x: k for k, x in enumerate(cn["genes"].astype(str))}
            tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
            R1 = X[:, [cg[x] for x in genes]].tocsr()
            R1.data = np.log1p(R1.data / np.repeat(tot, np.diff(R1.indptr)) * 1e4).astype(np.float32)
            R = sp.vstack([R1[cp[i]] if i in cp else R0[p0[i]] for i in A["ids"]]).tocsr().astype(np.float32)
        else:
            R = lognorm(X[np.array([cp[i] for i in A["ids"]])].tocsr(), ctl)
        V = views(R, ctl)
        np.savez(a.work / f"b1_{g}_in.npz", ids=A["ids"], A=A["A"], unit=A["unit"], **V, **{k: v for k, v in A.items() if k.startswith("y_")})
        for r, x in e["readouts"].items():
            if not x["admitted"]:
                continue
            ev = x["admitted"]
            write(a.work, spec(f"b1_{g}_{r}", st, r, ev, f"b1_{g}.npz", ["Z"], BOTH, other=OTHER[r]))
            write(a.work, spec(f"b1_{g}_{r}_in512", st, r, ev, f"b1_{g}_in.npz", ["RNA512"], BOTH, frozen=False, train=False))
            write(a.work, spec(f"b1_{g}_{r}_in4k", st, r, ev, f"b1_{g}_in.npz", ["RNA4K"], RIDGE, frozen=False, train=False))
            write(a.work, spec(f"b1_{g}_{r}_perm", st, r, ev, f"b1_{g}_in.npz", ["PERM"], RIDGE, frozen=False, train=False))
        res["strata"][g] = e
        print(g, e["cells_read"], {r: x["admitted"] for r, x in e["readouts"].items()})
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
