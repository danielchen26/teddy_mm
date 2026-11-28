#!/usr/bin/env python3
"""A made-up world in the formats of the real files, to find faults of the code before the one reading.

  python rehearsal.py --out DIR [--seed 1] [--null]

Writes DIR/{cells.tsv, library.npz, counts_mel.npz, counts_ifn.npz, pack.npz, states_mel.npz, planted.json}. Then
summarize.py, analyze.py --rehearsal and teddy.py --rehearsal are run on it by hand (README).

The world: 40 knockouts in three conditions. A signal changes 200 genes and three proteins; every knockout lowers its
own gene and moves 15 genes by the same factor in every condition (it adds). Planted, unless --null:
  the five of the known answer   remove the signal's effect entirely (in the co-culture, 60% of it)
  PART1                          removes half of it
  ONE1                           raises one gene five-fold under the signal only
  B2M                            lowers MHC-I in every condition (a declared pair; it adds)
Cells share a hidden factor that moves 300 genes together, so the noise of a group's mean is not independent across genes.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ka as K  # noqa: E402

PANEL = ["CD274", "HLA_A", "HLA_E", "CD117", "CD119", "CD140a", "CD140b", "CD172a", "CD184", "CD202b", "CD29", "CD309", "CD44", "CD47", "CD49f", "CD58", "CD59", "CD61", "CD9", "CD279",
         "Rat_IgG2a", "Mouse_IgG1", "Mouse_IgG2a", "Mouse_IgG2b"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--null", action="store_true")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    ko = list(K.FIVE) + ["PART1", "ONE1", "B2M"] + [f"KO{i:02d}" for i in range(32)]
    genes = np.array(ko + [f"G{i:04d}" for i in range(1500)])
    G = len(genes)
    base = rng.normal(-7.0, 1.2, G)                                                     # log share of a cell's counts
    sig = np.zeros(G); sig[rng.choice(np.arange(len(ko), G), 200, replace=False)] = rng.normal(0.8, 1.2, 200)
    signal = {"ctl": 0.0 * sig, "ifn": sig, "co": 0.6 * sig + np.where(rng.random(G) < 0.05, rng.normal(0, 1, G), 0)}
    load = np.zeros(G); load[rng.choice(G, 300, replace=False)] = rng.normal(0, 0.35, 300)
    eff = {g: np.zeros(G) for g in ko}
    for i, g in enumerate(ko):
        eff[g][rng.choice(np.arange(len(ko), G), 15, replace=False)] = rng.normal(0, 0.5, 15)
        eff[g][i] = np.log(0.2)
    removes = {} if a.null else {**{g: {"ifn": 1.0, "co": 0.6} for g in K.FIVE}, "PART1": {"ifn": 0.5, "co": 0.3}}
    quiet = np.where(sig == 0)[0][len(ko):]
    one_gene = int(quiet[np.argsort(base[quiet])[int(0.9 * len(quiet))]])                # a well-expressed gene the signal does not change
    pk = {r: j for j, r in enumerate(PANEL)}
    pbase = rng.normal(4.0, 0.8, 24)
    psig = np.zeros(24); psig[[pk["CD274"], pk["HLA_A"], pk["HLA_E"]]] = [2.0, 1.5, 1.0]
    Wz = rng.normal(size=(4, 512))
    cells, lib, counts, plog, pid, states = [], {"ids": [], "Z": [], "A": [], "target": []}, {"mel": ([], []), "ifn": ([], [])}, [], [], ([], [])
    n = 0
    for k in ("ctl", "ifn", "co"):
        for g in ["control"] + ko:
            m = 2000 if g == "control" else int(rng.integers(100, 300))
            ids = [f"CELL_{n + i}" for i in range(m)]; n += m
            kept = 1.0 - removes.get(g, {}).get(k, 0.0)
            mu = base + kept * signal[k] + (eff[g] if g != "control" else 0)
            if g == "ONE1" and k == "ifn" and not a.null:
                mu = mu.copy(); mu[one_gene] += np.log(5.0)
            t = rng.normal(size=m)
            lam = np.exp(mu[None, :] + t[:, None] * load[None, :] + rng.normal(0, 0.25, (m, 1)))
            X = sp.csr_matrix(rng.poisson(4000 * lam / lam.sum(1, keepdims=True) * np.exp(rng.normal(0, 0.3, (m, 1)))).astype(np.float32))
            sgn = {"ctl": 0.0, "ifn": 1.0, "co": 0.6}[k] * kept
            P = pbase + sgn * psig + rng.normal(0, 0.4, (m, 24)) + 0.3 * t[:, None] * (np.arange(24) < 6)
            if g == "B2M":
                P[:, pk["HLA_A"]] -= 2.0
            Z = np.c_[sgn + rng.normal(0, 0.5, m), t, rng.normal(size=(m, 2))] @ Wz * 0.6 + rng.normal(size=(m, 512)) + 1.0
            if k == "ifn":
                lib["ids"] += ids; lib["Z"].append(Z.astype(np.float32)); lib["A"] += [int(g != "control")] * m; lib["target"] += [g if g != "control" else "NO_SITE"] * m
            else:
                cells += [f"{K.STUDY}\t{K.CONDITION[k]}\t{g}\t{i}\t{'control' if g == 'control' else 'event'}" for i in ids]
                states[0].extend(ids); states[1].append(Z.astype(np.float32))
            c = counts["ifn" if k == "ifn" else "mel"]
            c[0].extend(ids); c[1].append(X)
            pid += ids; plog.append(P)
    (a.out / "cells.tsv").write_text("\n".join(["study\tcondition\ttarget\tcell\trole"] + cells) + "\n")
    np.savez(a.out / "library.npz", ids=np.array(lib["ids"]), Z=np.vstack(lib["Z"]), A=np.array(lib["A"]), target=np.array(lib["target"]), stratum=np.array([K.CONDITION["ifn"]] * len(lib["ids"])))
    for name, (ids, mats) in counts.items():
        M = sp.vstack(mats).tocoo()
        np.savez_compressed(a.out / f"counts_{name}.npz", rows=M.row.astype(np.int32), cols=M.col.astype(np.int32), vals=M.data.astype(np.float32), ids=np.array(ids), genes=genes, shape=np.array(M.shape))
    np.savez(a.out / "pack.npz", mel_panel=np.array(PANEL), mel_ids=np.array(pid), mel_log=np.vstack(plog), thp_panel=np.array(["CD86", "PDL1", "PDL2", "CD366"]), thp_ids=np.array([], dtype="<U8"), thp_log=np.zeros((0, 4)))
    np.savez(a.out / "states_mel.npz", ids=np.array(states[0]), z_last=np.vstack(states[1]), sets=np.array(["A|0"]))
    (a.out / "planted.json").write_text(json.dumps({"null": bool(a.null), "seed": a.seed, "knockouts": ko, "remove_the_signal": removes, "one_gene": None if a.null else {"ONE1": str(genes[one_gene])},
                                                    "genes_the_signal_changes": int((sig != 0).sum())}, indent=1) + "\n")
    print("cells", n, "| genes", G, "| planted", "nothing" if a.null else list(removes) + ["ONE1"])


if __name__ == "__main__":
    main()
