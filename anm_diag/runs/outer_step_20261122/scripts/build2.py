#!/usr/bin/env python3
"""Arrays and SPEC files of the first outer step, for `anm_diag scope`. No fit is made here except the components of RNA.

  python build2.py --work WORK2 --part bridge --coculture-counts COUNTS.npz
  python build2.py --work WORK2 --part target --selection ../results/selection.json --library-counts COUNTS.npz

Environment as runs/teg_loop_20261120 (ANM_LOOP_WORK: the work folder of its build.py; ANM_FEATURES; ANM_TEST_CELLS).

--part bridge
  thp1a            the THP-1 control cells and the four genes of type A (the provisional split): Z, ridge and boosting
  mel_co_MHC1_z2   MHC-I in co-culture at Z with both reader classes (the round had ridge only)
  mel_co_<r>_in*   the input view of co-culture, one SPEC per view: RNA512 (ridge, boosting), RNA4K (ridge), PERM (RNA512 with its
                   rows permuted across cells; ridge)
--part target
  mel_ifn_mhc      MHC-I in melanoma under IFN-gamma: control cells and admitted events outside the selection third, and the
                   five pathway genes as other units; Z (ridge, boosting) and the same input views
"""
import argparse
import copy
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml
from sklearn.decomposition import TruncatedSVD

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
from _common import FEATURES, forward_np, unit, weights  # noqa: E402
from build import ALPHAS, DECISION, HGB, OUTPUT  # noqa: E402

SEED_SVD, SEED = 20261120, 20261122
RIDGE = {"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}}
BOTH = {**RIDGE, "hgb": HGB}
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
PATHWAY = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
in_third = lambda cell, salt=str(SEED): int(hashlib.sha256(f"{salt}|{cell}".encode()).hexdigest(), 16) % 3 == 0


def lognorm(X, ctl):
    """log1p of counts per 10,000 on the genes detected in at least 1% of the control cells (the RNA view of the round)."""
    tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
    keep = np.asarray((X[ctl] > 0).mean(0)).ravel() >= 0.01
    R = X[:, np.where(keep)[0]].tocsr().astype(np.float32)
    R.data = np.log1p(R.data / np.repeat(tot, np.diff(R.indptr)) * 1e4).astype(np.float32)
    return R


def views(R, ctl):
    """Three views of the RNA of a group. No label and no readout enters."""
    rna512 = TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED_SVD).fit_transform(R).astype(np.float32)
    Rc = R[ctl]
    var = np.asarray(Rc.multiply(Rc).mean(0)).ravel() - np.asarray(Rc.mean(0)).ravel() ** 2
    top = np.sort(np.argsort(-var)[:4000])
    perm = np.random.default_rng(SEED).permutation(R.shape[0])
    return {"RNA512": rna512, "RNA4K": np.asarray(R[:, top].todense(), np.float32), "PERM": rna512[perm]}


def spec(name, study, readout, events, arrays, stages, classes, frozen=True, other=None, train=True):
    s = dict(name=name, kind="perturbation", study=study, readout_name=readout, events=list(events), also_pooled=True, units_exclude=["control"], reads=list(stages)[-1],
             control_folds=5, decision=DECISION, label_basis="design", keep_draws=True, readout_classes=copy.deepcopy(classes),
             keys=dict(readouts=[f"y_{readout}"], stages={k: k for k in stages}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit"), arrays={"audit": arrays})
    if frozen:
        s["keys"]["frozen"] = {"array": f"frozen_{readout}"}
    if train:
        s["arrays"]["train"] = "bm_train.npz"
    if other:
        s["arrays"]["other"] = other
    return s


def write(work, s):
    (work / f"{s['name']}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    print("SPEC", s["name"], "| stages", list(s["keys"]["stages"]), "| classes", list(s["readout_classes"]), "| events", s["events"])


def part_bridge(a, src):
    for f in ("bm_train.npz", "mel_ifn.npz", "mel_co.npz", "thp1.npz"):
        if not (a.work / f).exists():
            os.symlink(src / f, a.work / f)
    K = dict(np.load(src / "thp1.npz"))
    m = np.isin(K["unit"].astype(str), TYPE_A + ["control"])
    np.savez(a.work / "thp1a.npz", **{k: v[m] for k, v in K.items()})
    write(a.work, spec("thp1a_PDL1", "papalexi2021", "PDL1", TYPE_A, "thp1a.npz", ["Z"], BOTH, other="mel_ifn.npz"))
    write(a.work, spec("mel_co_MHC1_z2", "frangieh2021", "MHC1", PATHWAY, "mel_co.npz", ["Z"], BOTH, other="mel_ifn.npz"))
    K = dict(np.load(src / "mel_co.npz"))
    c = np.load(a.coculture_counts)
    assert np.array_equal(c["ids"].astype(str), K["ids"].astype(str))
    X = sp.csr_matrix((c["vals"].astype(np.float32), (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]))
    V = views(lognorm(X, np.where(K["A"] == 0)[0]), np.where(K["A"] == 0)[0])
    old = np.load(src / "mel_co_rna.npz")["RNA512"]
    assert np.abs(V["RNA512"] - old).max() < 1e-3, "the 512 components are not those of the round"
    np.savez(a.work / "mel_co_in.npz", ids=K["ids"], A=K["A"], unit=K["unit"], **V, **{k: v for k, v in K.items() if k.startswith("y_")})
    for r in ("PDL1", "MHC1"):
        write(a.work, spec(f"mel_co_{r}_in512", "frangieh2021", r, PATHWAY, "mel_co_in.npz", ["RNA512"], BOTH, frozen=False, train=False))
        write(a.work, spec(f"mel_co_{r}_in4k", "frangieh2021", r, PATHWAY, "mel_co_in.npz", ["RNA4K"], RIDGE, frozen=False, train=False))
        write(a.work, spec(f"mel_co_{r}_perm", "frangieh2021", r, PATHWAY, "mel_co_in.npz", ["PERM"], RIDGE, frozen=False, train=False))


def part_target(a, src):
    sel = json.loads(a.selection.read_text())
    events = sel["MHC1"]["admitted"]
    assert events, "no event was admitted on the selection third"
    lib = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "library.npz")
    ids, tg = lib["ids"].astype(str), lib["target"].astype(str)
    third = np.array([in_third(i) for i in ids])
    keep = np.isin(tg, events + ["control"] + PATHWAY) & ~third          # the selection third of every candidate is left out
    idx = np.where(keep)[0]
    P0, names, _j = weights()
    K = dict(ids=ids[idx], Z=lib["Z"][idx].astype(np.float32), A=lib["A"][idx].astype(int), unit=tg[idx], y_MHC1=lib["readout__MHC1"][idx].astype(np.float64))
    K["frozen_MHC1"] = forward_np(P0, unit(K["Z"]))[1][:, names.index(OUTPUT["MHC1"])].astype(np.float64)
    np.savez(a.work / "mel_ifn_mhc.npz", **K)
    u, c = np.unique(K["unit"], return_counts=True)
    print("mel_ifn_mhc", dict(zip(u.tolist(), c.tolist())), "| cells of the selection third left out:", int((np.isin(tg, events + ["control"] + PATHWAY) & third).sum()))
    write(a.work, spec("mel_ifn_mhc_MHC1", "frangieh2021", "MHC1", events, "mel_ifn_mhc.npz", ["Z"], BOTH, other="mel_co.npz"))
    # the input view: the round's RNA for the control cells and the pathway genes, the screen's RNA file for the admitted events
    old = np.load(FEATURES / "frangieh2021" / "dense.npz")
    oid = old["ids"].astype(str)
    genes = json.loads((FEATURES / "frangieh2021" / "raw_genes.json").read_text())
    R0 = sp.load_npz(FEATURES / "frangieh2021" / "node__RAW.npz").tocsr()
    c = np.load(a.library_counts)
    cg = {g: k for k, g in enumerate(c["genes"].astype(str))}
    X = sp.csr_matrix((c["vals"].astype(np.float32), (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]))
    tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
    R1 = X[:, [cg[g] for g in genes]].tocsr()
    R1.data = np.log1p(R1.data / np.repeat(tot, np.diff(R1.indptr)) * 1e4).astype(np.float32)
    pos0, pos1 = {i: k for k, i in enumerate(oid)}, {i: k for k, i in enumerate(c["ids"].astype(str))}
    rows = [R0[pos0[i]] if i in pos0 else R1[pos1[i]] for i in K["ids"]]
    R = sp.vstack(rows).tocsr().astype(np.float32)
    V = views(R, np.where(K["A"] == 0)[0])
    np.savez(a.work / "mel_ifn_mhc_in.npz", ids=K["ids"], A=K["A"], unit=K["unit"], y_MHC1=K["y_MHC1"], **V)
    write(a.work, spec("mel_ifn_mhc_MHC1_in512", "frangieh2021", "MHC1", events, "mel_ifn_mhc_in.npz", ["RNA512"], BOTH, frozen=False, train=False))
    write(a.work, spec("mel_ifn_mhc_MHC1_in4k", "frangieh2021", "MHC1", events, "mel_ifn_mhc_in.npz", ["RNA4K"], RIDGE, frozen=False, train=False))
    write(a.work, spec("mel_ifn_mhc_MHC1_perm", "frangieh2021", "MHC1", events, "mel_ifn_mhc_in.npz", ["PERM"], RIDGE, frozen=False, train=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--part", choices=["bridge", "target"], required=True)
    ap.add_argument("--coculture-counts", type=Path)
    ap.add_argument("--library-counts", type=Path)
    ap.add_argument("--selection", type=Path)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    src = Path(os.environ["ANM_LOOP_WORK"])
    (part_bridge if a.part == "bridge" else part_target)(a, src)


if __name__ == "__main__":
    main()
