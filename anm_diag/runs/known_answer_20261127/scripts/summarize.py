#!/usr/bin/env python3
"""The sums the rule stands on: for every condition, every group of cells (a knockout, or the control cells) and each
half, the number of cells and the sums and sums of squares of every readout. No decision is taken here.

  python summarize.py --cells final_cells.tsv --library library.npz --counts MEL_COUNTS... --counts-ifn counts_ifn.npz
                      --packs PACK... --out summary.npz

RNA: counts per 10,000 of every cell (over all genes of the file); kept are the genes whose mean among the control
cells is at least MIN_MEAN in any condition. Protein: the centred log-ratio over the 24 antibody columns, and beside
it the log count against the mean of the isotype antibodies.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ka as K  # noqa: E402


def groups(C, units):
    """Group and half of every cell of a condition: an index into (units x 2), or -1 for a cell of no group."""
    pos = {u: k for k, u in enumerate(units)}
    g = np.array([pos.get(u, -1) for u in C["unit"]])
    return np.where(g >= 0, 2 * g + (~C["A"]).astype(int), -1)


def rna_sums(files, ids, idx, n_groups):
    """Sums and sums of squares of counts per 10,000 by group, over the count files; the first file that holds a cell is used."""
    want = {i: k for k, i in enumerate(ids)}
    seen = np.zeros(len(ids), bool)
    s1 = s2 = genes = None
    for f in files:
        d = np.load(f)
        assert genes is None or (genes == d["genes"]).all(), "the count files list different genes"
        genes = d["genes"].astype(str)
        fid = d["ids"].astype(str)
        X = sp.csr_matrix((d["vals"], (d["rows"], d["cols"])), shape=tuple(int(v) for v in d["shape"]), dtype=np.float32)
        del d
        k = np.array([want.get(i, -1) for i in fid])
        take = np.where((k >= 0) & ~np.r_[seen, False][k])[0] if len(k) else np.array([], int)
        if not len(take):
            continue
        seen[k[take]] = True
        X = X[take]
        tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
        X = sp.diags((1e4 / tot).astype(np.float32)) @ X
        gi = idx[k[take]]
        ok = gi >= 0
        M = sp.csr_matrix((np.ones(int(ok.sum()), np.float64), (gi[ok], np.where(ok)[0])), shape=(2 * n_groups, X.shape[0]))
        a, b = np.asarray((M @ X).todense(), np.float64), np.asarray((M @ X.multiply(X)).todense(), np.float64)
        s1, s2 = (a, b) if s1 is None else (s1 + a, s2 + b)
        del X, M
    assert seen[idx >= 0].all(), "a cell of a group has no RNA counts in the files"
    return s1, s2, genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--counts", type=Path, nargs="+", required=True)
    ap.add_argument("--counts-ifn", type=Path, nargs="+", required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"{a.out.name} exists: the sums are made once")
    C = K.cells(a.cells, a.library)
    ko = {t: K.knockouts(C, t) for t in K.TESTS}
    units = {k: ["control"] + sorted({g for t, (x, y) in K.TESTS.items() if k in (x, y) for g in ko[t]}) for k in C}
    packs = K.B3.Packs(a.packs)
    out = {"proteins": np.array(list(K.B3.PANEL[K.STUDY])), **{f"knockouts_{t}": np.array(v) for t, v in ko.items()}}
    S = {}
    for k, c in C.items():
        idx = groups(c, units[k])
        n = np.bincount(idx[idx >= 0], minlength=2 * len(units[k])).reshape(-1, 2)
        s1, s2, genes = rna_sums(a.counts_ifn if k == "ifn" else a.counts, c["ids"].tolist(), idx, len(units[k]))
        m = idx >= 0
        clr, log, bg = packs.readouts(K.STUDY, c["ids"][m])
        P = np.array([clr[r] for r in out["proteins"]]).T
        Q = np.array([log[r] - bg for r in out["proteins"]]).T
        M = sp.csr_matrix((np.ones(int(m.sum())), (idx[m], np.arange(int(m.sum())))), shape=(2 * len(units[k]), int(m.sum())))
        S[k] = {"units": units[k], "n": n, "rna": (s1, s2), "prot": (M @ P, M @ (P * P)), "iso": (M @ Q, M @ (Q * Q))}
        out["genes_all"] = genes
        print(k, len(units[k]) - 1, "knockouts,", int(n.sum()), "cells, control", n[0].tolist(), flush=True)
    mean_ctl = np.array([S[k]["rna"][0][:2].sum(0) / S[k]["n"][0].sum() for k in S])
    keep = (mean_ctl >= K.MIN_MEAN).any(0)
    out["genes"] = out.pop("genes_all")[keep]
    for k, s in S.items():
        sh = (len(s["units"]), 2, -1)
        out.update({f"{k}_units": np.array(s["units"]), f"{k}_n": s["n"], f"{k}_rna_sum": s["rna"][0][:, keep].reshape(sh), f"{k}_rna_sq": s["rna"][1][:, keep].reshape(sh),
                    f"{k}_prot_sum": s["prot"][0].reshape(sh), f"{k}_prot_sq": s["prot"][1].reshape(sh), f"{k}_iso_sum": s["iso"][0].reshape(sh), f"{k}_iso_sq": s["iso"][1].reshape(sh)})
    np.savez_compressed(a.out, **out)
    print("genes kept", int(keep.sum()), "of", len(keep), "| proteins", len(out["proteins"]), "| sha256", K.sha(a.out)[:16])


if __name__ == "__main__":
    main()
