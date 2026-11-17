#!/usr/bin/env python3
"""Step 6 of gate 5: the grid on the final study F, read after the test. For round 2; it changes no outcome.

  python grid_final.py --heads HEADS_DIR --work WORK --bm-arrays cite_arrays.npz --pbmc-arrays ext_cite_arrays.npz \\
      --final-arrays F_cite_arrays.npz --final-z F_z_rna.npy --tmp TMP [--seeds 0] [--only donors disease] [--summary-only]

Refuses to run before ../results/final.json exists: the grid is read only after the test.

WORK: the arrays runs/scope_20261116/scripts/build.py wrote (bm_train.npz, hao_donors.npz). The RNA view of F is its
log-normalised counts projected on the 512 components of that build: the components are fitted again on the same cells
(bone marrow and the PBMC study) with the same seed, checked against the stored view, and F is projected on them.

For every seed asked for, two runs with the donors of F as events (each donor against the other donors of F): the
`coverage` head of the test and `as_is` of the same seed. The training cells of a run are the cells its head was fitted
on. And one run per head with the disease as the event: the patients against the healthy donors, donors as clusters.
Every run names the donors of F under roles.final and is a final audit of the two hypotheses of the declaration.

Writes ../results/grid_final/<run>/{SPEC.yaml, scope.json, scope_events.jsonl} and ../results/what_is_left_final.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import common
import grid_changed
import plan
import read_final
from anm_diag import scope

HEADS = {"as_is": "as_is", "coverage": "coverage_all"}
STAGES = grid_changed.STAGES
DECISION = grid_changed.DECISION
SEED_SVD = 0                                             # runs/scope_20261116/scripts/build.py


def lognorm(a):
    import scipy.sparse as sp
    X = sp.csr_matrix((a["rna_data"].astype(np.float32), a["rna_indices"], a["rna_indptr"]), shape=tuple(int(v) for v in a["rna_shape"]))
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1.0
    X.data *= np.repeat((1e4 / tot).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def rna_view(bm_arrays, pbmc_arrays, final_arrays, stored, out):
    """The 512 RNA components of F in the basis of the stored views. stored: RNA512 of the PBMC study (hao_donors.npz)."""
    import scipy.sparse as sp
    from sklearn.decomposition import TruncatedSVD
    if Path(out).exists():
        return np.load(out)
    bm, ex, fi = (np.load(p, allow_pickle=False) for p in (bm_arrays, pbmc_arrays, final_arrays))
    if not (np.array_equal(bm["rna_names"], ex["rna_names"]) and np.array_equal(bm["rna_names"], fi["rna_names"])):
        raise SystemExit("the three packs do not list the same genes in the same order")
    R = sp.vstack([lognorm(bm), lognorm(ex)]).tocsr()
    svd = TruncatedSVD(n_components=512, algorithm="randomized", n_iter=2, random_state=SEED_SVD).fit(R)
    again = svd.transform(lognorm(ex)[:2000]).astype(np.float32)
    err = float(np.abs(again - stored[:2000]).max() / np.abs(stored[:2000]).max())
    if err > 1e-3:
        raise SystemExit(f"the components do not reproduce the stored RNA view (relative difference {err:.2e})")
    del R
    V = svd.transform(lognorm(fi)).astype(np.float32)
    np.save(out, V)
    print(f"RNA view of F: {V.shape}; the stored view is reproduced to {err:.1e}", flush=True)
    return V


def arrays_of_run(W, cells, F, rna_f, cols112, wb, wh, bm_split):
    """Audit: the cells of F. Train: the cells the head was fitted on, proteins restricted to those F measures."""
    lat, out, mu = grid_changed.head_stages(W, F["z_raw"])
    audit = dict(RNA512=rna_f, Z=F["z_raw"], HEAD_LAT=lat, HEAD_OUT=out, Y=F["Y"].astype(np.float32), unit=F["donors"], donor=F["donors"],
                 status=F["status"], frozen=mu[:, F["cols"]])
    tr = np.where(np.asarray(bm_split).astype(str) == "train")[0]
    pos = np.searchsorted(tr, cells["bm_index"])
    if pos.max() >= len(tr) or not np.array_equal(tr[pos], cells["bm_index"]):
        raise SystemExit("the head's bone-marrow cells are not training cells of this pack")
    pi = np.asarray(cells["pbmc_index"], np.int64)
    latb, outb, _ = grid_changed.head_stages(W, wb["Z"][pos])
    latp, outp, _ = grid_changed.head_stages(W, wh["Z"][pi]) if len(pi) else (np.zeros((0, lat.shape[1]), np.float32), np.zeros((0, out.shape[1]), np.float32), None)
    train = dict(RNA512=np.vstack([wb["RNA512"][pos], wh["RNA512"][pi]]), Z=np.vstack([wb["Z"][pos], wh["Z"][pi]]), HEAD_LAT=np.vstack([latb, latp]),
                 HEAD_OUT=np.vstack([outb, outp]), Y=np.vstack([wb["Y"][pos][:, cols112], wh["Y"][pi][:, cols112]]),
                 donor=np.concatenate([np.asarray(wb["donor"]).astype(str)[pos], np.asarray(wh["donor"]).astype(str)[pi]]))
    return audit, train


def spec(name, kind, donors, n_proteins):
    base = dict(name=name, study="stephenson2021", arrays={"audit": f"{name}__audit.npz", "train": f"{name}__train.npz"}, reads="Z", input_view="RNA512",
                variation=True, max_fit=30000, decision=DECISION, label_basis="design", readout_name=f"{n_proteins} proteins",
                readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": common.ALPHAS}},
                roles={"final": list(donors)}, final_audit=True, hypotheses=["donors_hao_coverage", "coverage_not_levels"])
    stages = {s: s for s in STAGES}
    if kind == "donor":
        return dict(base, kind="donor", keys=dict(readouts="Y", stages=stages, unit="unit", block="donor", train_block="donor", frozen={"array": "frozen"}),
                    events="each", confounded_with=["sample", "centre", "disease status"])
    return dict(base, kind="disease", keys=dict(readouts="Y", stages=stages, unit="status", block="donor", cluster="donor", train_block="donor", frozen={"array": "frozen"}),
                events=["patient"], confounded_with=["centre", "sample", "age and sex of the donors", "time of sampling"],
                pooling={"site": "the patients of the test, pooled into one event"})


def one_run(name, kind, audit, train, donors, tmp, res):
    import yaml
    np.savez(tmp / f"{name}__audit.npz", **audit)
    np.savez(tmp / f"{name}__train.npz", **train)
    sp = spec(name, kind, donors, audit["Y"].shape[1])
    (tmp / f"{name}.yaml").write_text(yaml.safe_dump(sp, sort_keys=False))
    scope.main([str(tmp / f"{name}.yaml"), "--out", str(res / name)])
    (res / name / "SPEC.yaml").write_text((tmp / f"{name}.yaml").read_text())
    for f in ("__audit.npz", "__train.npz"):
        (tmp / f"{name}{f}").unlink()


def cell(T, scope_, stage):
    return abs(T["stages"]["system" if scope_ == "frozen" else stage][scope_]["u"])


def summarise(res, seeds, status_of):
    """Means over the donors of F, for each head, of u at the head's own output and under every readout; and by group."""
    donors = sorted(status_of)
    cells = [("frozen", "system")] + [(sc, st) for sc in ("train", "units", "own") for st in STAGES]
    key = lambda c: c[0] if c[0] == "frozen" else f"{c[0]}@{c[1]}"
    out = {"donors": len(donors), "seeds": list(seeds), "note": "each donor against the other donors of F; the readouts of the `units` scope are fitted on those"}
    for tg in ("variation", "shift"):
        U = {h: {c: np.full((len(donors), len(seeds)), np.nan) for c in cells} for h in HEADS}
        read = {h: {} for h in HEADS}
        for j, s in enumerate(seeds):
            for h in HEADS:
                r = json.loads((res / f"s{s}_donors_{h}" / "scope.json").read_text())
                for e in r["results"]:
                    T = e["targets"].get(tg)
                    if T is None:
                        continue
                    i = donors.index(e["event"]["unit"])
                    for c in cells:
                        U[h][c][i, j] = cell(T, *c)
                    read[h][T["reading"]["missing"]] = read[h].get(T["reading"]["missing"], 0) + 1
        g = lambda a, b, rows=slice(None): {x: scope.arms_family({"a": a[rows], "b": b[rows]}, base="a")["arms"]["b"][x] for x in ("gain", "ci")}
        hl = np.array([status_of[d] == "healthy" for d in donors])
        fz, un = ("frozen", "system"), ("units", "Z")
        out[tg] = {
            "u": {h: {key(c): float(np.nanmean(v)) for c, v in U[h].items()} for h in HEADS},
            "ceiling_as_is": {"all": g(U["as_is"][fz], U["as_is"][un]), "healthy": g(U["as_is"][fz], U["as_is"][un], hl), "patients": g(U["as_is"][fz], U["as_is"][un], ~hl)},
            "left_after_the_change": {"all": g(U["coverage"][fz], U["coverage"][un]), "healthy": g(U["coverage"][fz], U["coverage"][un], hl),
                                      "patients": g(U["coverage"][fz], U["coverage"][un], ~hl)},
            "fall_from_as_is_to_coverage": {key(c): g(U["as_is"][c], U["coverage"][c]) for c in cells},
            "state_against_rna": {h: g(U[h][("units", "RNA512")], U[h][("units", "Z")]) for h in HEADS},
            "own_cells_after_other_donors": g(U["as_is"][("units", "Z")], U["as_is"][("own", "Z")]),
            "readings": read,
        }
    dis = {}
    for h in HEADS:
        p = res / f"s{seeds[0]}_disease_{h}" / "scope.json"
        if p.exists():
            e = json.loads(p.read_text())["results"][0]
            dis[h] = {tg: {"reading": T["reading"]["missing"], "u": {key(c): cell(T, *c) for c in cells if (c[1] if c[0] != "frozen" else "system") in T["stages"]
                                                                       and c[0] in T["stages"][c[1] if c[0] != "frozen" else "system"]}}
                      for tg, T in e["targets"].items()}
    out["disease"] = dis
    return out


def main():
    p = argparse.ArgumentParser()
    for k in ("heads", "work", "bm-arrays", "pbmc-arrays", "final-arrays", "final-z", "tmp"):
        p.add_argument(f"--{k}", type=Path, required=True)
    p.add_argument("--seeds", type=int, nargs="+", default=[0])
    p.add_argument("--only", nargs="+", default=["donors", "disease"])
    p.add_argument("--heads-only", nargs="+", default=list(HEADS))
    p.add_argument("--summary-only", action="store_true")
    a = p.parse_args()
    res0 = common.RUN / "results"
    if not (res0 / "final.json").exists():
        raise SystemExit("the test has not been read; the grid on F comes after it")
    decl = json.loads((res0 / "gate5_declaration.json").read_text())
    read_final.check_declaration({"sha256": {k: v for k, v in decl["sha256"].items() if k in ("final_arrays", "final_z")}},
                                 {"final_arrays": a.final_arrays, "final_z": a.final_z})
    res = res0 / "grid_final"
    status_of = {d["donor"]: d["status"] for d in decl["donors"]}
    todo = [(f"s{s}_{kind}_{h}", s, kind, h) for s in a.seeds for kind in a.only for h in a.heads_only if kind == "donors" or s == a.seeds[0]]
    left = [t for t in todo if not (res / t[0] / "scope.json").exists()]
    if left and not a.summary_only:
        a.tmp.mkdir(parents=True, exist_ok=True)
        wb, wh = np.load(a.work / "bm_train.npz", allow_pickle=False), np.load(a.work / "hao_donors.npz", allow_pickle=False)
        wb, wh = {k: wb[k] for k in ("RNA512", "Z", "Y", "donor")}, {k: wh[k] for k in ("RNA512", "Z", "Y", "donor")}
        bm, ex, fi = np.load(a.bm_arrays, allow_pickle=False), np.load(a.pbmc_arrays, allow_pickle=False), np.load(a.final_arrays, allow_pickle=False)
        shared = np.where(ex["adt_present"].astype(bool))[0]
        names = [str(x) for x in fi["adt_names"]]
        cols = np.array([names.index(n) for n in decl["proteins"]])
        cols112 = np.array([int(np.where(shared == c)[0][0]) for c in cols])
        Y = fi["adt"][:, cols]
        if np.isnan(Y).any():
            raise SystemExit("a protein of the declaration is not measured in every cell of F")
        F = {"z_raw": np.load(a.final_z).astype(np.float32), "Y": Y, "donors": fi["donors"].astype(str), "status": fi["status"].astype(str), "cols": cols}
        rna_f = rna_view(a.bm_arrays, a.pbmc_arrays, a.final_arrays, wh["RNA512"], a.tmp / "final_rna512.npy")
        for name, s, kind, h in left:
            t0 = time.time()
            print(f"== {name}", flush=True)
            d = a.heads / f"s{s}_{HEADS[h]}"
            audit, train = arrays_of_run(dict(np.load(d / "head.npz")), dict(np.load(d / "fit_cells.npz")), F, rna_f, cols112, wb, wh, bm["split"])
            one_run(name, "donor" if kind == "donors" else "disease", audit, train, sorted(status_of), a.tmp, res)
            print(f"   done in {time.time() - t0:.0f} s", flush=True)
    if any(not (res / f"s{s}_donors_{h}" / "scope.json").exists() for s in a.seeds for h in HEADS):
        raise SystemExit("runs are missing; nothing is summarised before the donor runs of both heads are there")
    out = summarise(res, a.seeds, status_of)
    (res0 / "what_is_left_final.json").write_text(json.dumps(out, indent=1))
    v = out["variation"]
    for h in HEADS:
        print(h, {k: round(x, 4) for k, x in v["u"][h].items() if k in ("frozen", "train@Z", "units@RNA512", "units@Z", "own@Z")}, v["readings"][h])
    print("ceiling (as_is):", {k: round(x["gain"], 4) for k, x in v["ceiling_as_is"].items()}, "left:", {k: round(x["gain"], 4) for k, x in v["left_after_the_change"].items()})


if __name__ == "__main__":
    main()
