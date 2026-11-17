#!/usr/bin/env python3
"""Step 6 of gate 5: the grid on the final study F, read after the test. For a next round; it changes no outcome.

  python grid_final.py --heads HEADS_DIR --work WORK --bm-arrays cite_arrays.npz --pbmc-arrays ext_cite_arrays.npz \\
      --final-arrays F_cite_arrays.npz --final-z F_z_rna.npy --tmp TMP [--seeds 0]

Refuses to run before ../results/final.json exists: the grid is read only after the test.

What it reads (`donors_variation`): for every donor of F, u of the variation among its cells at the head's own output
and under readouts that were not fitted on the donor (the head's training cells; the other donors of F; the donor's
own cells, cross-fitted), on two views of a cell (512 RNA components, TEDDY's state), for `as_is` and the `coverage`
head of the seed. Writes ../results/grid_final/s<seed>_donors_variation.json and ../results/what_is_left_final.json.

WORK: the arrays runs/scope_20261116/scripts/build.py wrote (bm_train.npz, hao_donors.npz). The RNA view is the
log-normalised counts projected on the 512 components of that build: the decomposition is fitted again on the same cells
(bone marrow and the PBMC study) with the same seed, has to return the stored view, and then the cells of F and the
training cells alike are projected on its components (see `rna_view`).

--scope-runs: the runs this script was first written for, scope runs with the donors and with the disease as events.
They were not made: a donor run also reads the shift target, which needs a readout for every pair of donors (about
ten hours a run on the machine of the round), and the disease run would resample the patients' cells one by one and
not by donor. GATE5.md, part 5.
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
SEED_SVD = 20261116                                      # SEED of runs/scope_20261116/scripts/build.py


def lognorm(a):
    import scipy.sparse as sp
    X = sp.csr_matrix((a["rna_data"].astype(np.float32), a["rna_indices"], a["rna_indptr"]), shape=tuple(int(v) for v in a["rna_shape"]))
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1.0
    X.data *= np.repeat((1e4 / tot).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def rna_view(bm_arrays, pbmc_arrays, final_arrays, stored, out):
    """The cells of F, of bone marrow and of the PBMC study projected on the 512 RNA components of the stored views.

    The stored views (runs/scope_20261116/scripts/build.py) are the left singular vectors times the singular values of a
    randomized decomposition, which a cell outside the fit does not have. The decomposition is fitted again on the same
    cells with the same seed and has to return the stored view; then every cell, of the fit and of F alike, is projected
    on its components, so that a readout fitted on training cells and applied to F reads one and the same view.
    stored: RNA512 of the PBMC study (hao_donors.npz). Returns (F, bone marrow, PBMC study)."""
    import scipy.sparse as sp
    from sklearn.decomposition import TruncatedSVD
    if Path(out).exists():
        d = np.load(out)
        return d["final"], d["bm"], d["pbmc"]
    bm, ex, fi = (np.load(p, allow_pickle=False) for p in (bm_arrays, pbmc_arrays, final_arrays))
    if not (np.array_equal(bm["rna_names"], ex["rna_names"]) and np.array_equal(bm["rna_names"], fi["rna_names"])):
        raise SystemExit("the three packs do not list the same genes in the same order")
    Rb, Re = lognorm(bm), lognorm(ex)
    svd = TruncatedSVD(n_components=512, algorithm="randomized", n_iter=2, random_state=SEED_SVD)
    again = svd.fit_transform(sp.vstack([Rb, Re]).tocsr())[Rb.shape[0]:Rb.shape[0] + 2000].astype(np.float32)
    err = float(np.abs(again - stored[:2000]).max() / np.abs(stored[:2000]).max())
    if err > 1e-3:
        raise SystemExit(f"the decomposition does not return the stored RNA view (relative difference {err:.2e})")
    Vb, Ve, Vf = (svd.transform(R).astype(np.float32) for R in (Rb, Re, lognorm(fi)))
    gap = float(np.abs(Ve[:2000] - stored[:2000]).max() / np.abs(stored[:2000]).max())
    np.savez(out, final=Vf, bm=Vb, pbmc=Ve)
    print(f"RNA view of F: {Vf.shape}; the stored view is returned to {err:.1e}; projection against the stored view of the same cells: {gap:.2e}", flush=True)
    return Vf, Vb, Ve


def arrays_of_run(W, cells, F, rna, cols112, wb, wh, bm_split):
    """Audit: the cells of F. Train: the cells the head was fitted on, proteins restricted to those F measures."""
    lat, out, mu = grid_changed.head_stages(W, F["z_raw"])
    rna_f, rna_b, rna_p = rna
    audit = dict(RNA512=rna_f, Z=F["z_raw"], HEAD_LAT=lat, HEAD_OUT=out, Y=F["Y"].astype(np.float32), unit=F["donors"], donor=F["donors"],
                 status=F["status"], frozen=mu[:, F["cols"]])
    tr = np.where(np.asarray(bm_split).astype(str) == "train")[0]
    pos = np.searchsorted(tr, cells["bm_index"])
    if pos.max() >= len(tr) or not np.array_equal(tr[pos], cells["bm_index"]):
        raise SystemExit("the head's bone-marrow cells are not training cells of this pack")
    pi = np.asarray(cells["pbmc_index"], np.int64)
    latb, outb, _ = grid_changed.head_stages(W, wb["Z"][pos])
    latp, outp, _ = grid_changed.head_stages(W, wh["Z"][pi]) if len(pi) else (np.zeros((0, lat.shape[1]), np.float32), np.zeros((0, out.shape[1]), np.float32), None)
    train = dict(RNA512=np.vstack([rna_b[tr[pos]], rna_p[pi]]), Z=np.vstack([wb["Z"][pos], wh["Z"][pi]]), HEAD_LAT=np.vstack([latb, latp]),
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


def prepare(name, kind, audit, train, donors, tmp):
    """Write the arrays and the SPEC of a run."""
    import yaml
    np.savez(tmp / f"{name}__audit.npz", **audit)
    np.savez(tmp / f"{name}__train.npz", **train)
    sp = spec(name, kind, donors, audit["Y"].shape[1])
    (tmp / f"{name}.yaml").write_text(yaml.safe_dump(sp, sort_keys=False))


def run_prepared(name, tmp, res):
    """Read the grid of a run whose arrays and SPEC are written."""
    scope.main([str(tmp / f"{name}.yaml"), "--out", str(res / name)])
    (res / name / "SPEC.yaml").write_text((tmp / f"{name}.yaml").read_text())
    for f in ("__audit.npz", "__train.npz"):
        (tmp / f"{name}{f}").unlink()


VCELLS = (("frozen", "system"), ("train", "RNA512"), ("train", "Z"), ("units", "RNA512"), ("units", "Z"), ("own", "RNA512"), ("own", "Z"))


def vkey(c):
    return c[0] if c[0] == "frozen" else f"{c[0]}@{c[1]}"


def donors_variation(a, decl, seed, out):
    """For every donor of F: u of the `variation` target at the head's own output and under every readout, for
    `as_is` and for `coverage` of a seed. Each donor is read against readouts that were not fitted on it.

    These are the numbers a scope run of kind `donor` gives for that target: the same `scope.Study`, the same readouts
    (`Study.predict`: the head's training cells; the other donors of F; the donor's own cells, cross-fitted) and the
    same resamples of the donor's cells (`scope._u_variation`). The `shift` target of such a run is not read: it needs
    a readout for every pair of donors (780 pairs at every stage, for every head), which this machine does not do in a
    night. The stages inside the head are not read either. The readouts fitted on cells of F do not depend on the head
    and are fitted once."""
    wb, wh = np.load(a.work / "bm_train.npz", allow_pickle=False), np.load(a.work / "hao_donors.npz", allow_pickle=False)
    stored = wh["RNA512"]
    wb, wh = {k: wb[k] for k in ("Z", "Y", "donor")}, {k: wh[k] for k in ("Z", "Y", "donor")}
    bm, ex, fi = np.load(a.bm_arrays, allow_pickle=False), np.load(a.pbmc_arrays, allow_pickle=False), np.load(a.final_arrays, allow_pickle=False)
    shared = np.where(ex["adt_present"].astype(bool))[0]
    names = [str(x) for x in fi["adt_names"]]
    cols = np.array([names.index(n) for n in decl["proteins"]])
    cols112 = np.array([int(np.where(shared == c)[0][0]) for c in cols])
    Y = fi["adt"][:, cols]
    if np.isnan(Y).any():
        raise SystemExit("a protein of the declaration is not measured in every cell of F")
    z_f, donors = np.load(a.final_z).astype(np.float32), fi["donors"].astype(str)
    rna_f, rna_b, rna_p = rna_view(a.bm_arrays, a.pbmc_arrays, a.final_arrays, stored, a.tmp / "rna512.npz")
    del stored
    units = [d["donor"] for d in decl["donors"]]
    sp = spec("donors", "donor", units, len(cols))
    d = {**scope.DEFAULT_DECISION, **sp["decision"]}
    tr = np.where(np.asarray(bm["split"]).astype(str) == "train")[0]
    first, res = None, {"what": "u of the variation among a donor's cells, donor by donor", "seed": seed, "decision": d, "max_fit": sp["max_fit"],
                        "readout_classes": sp["readout_classes"], "heads": {}}
    for h in HEADS:
        t0 = time.time()
        hd = a.heads / f"s{seed}_{HEADS[h]}"
        W, cells = dict(np.load(hd / "head.npz")), dict(np.load(hd / "fit_cells.npz"))
        mu = grid_changed.head_stages(W, z_f)[2]
        pos = np.searchsorted(tr, cells["bm_index"])
        if pos.max() >= len(tr) or not np.array_equal(tr[pos], cells["bm_index"]):
            raise SystemExit("the head's bone-marrow cells are not training cells of this pack")
        pi = np.asarray(cells["pbmc_index"], np.int64)
        train = {"stages": {"RNA512": np.vstack([rna_b[tr[pos]], rna_p[pi]]), "Z": np.vstack([wb["Z"][pos], wh["Z"][pi]])},
                 "Y": np.vstack([wb["Y"][pos][:, cols112], wh["Y"][pi][:, cols112]]).astype(np.float64),
                 "block": np.concatenate([np.asarray(wb["donor"]).astype(str)[pos], np.asarray(wh["donor"]).astype(str)[pi]])}
        study = scope.Study(Y.astype(np.float64), {"RNA512": rna_f, "Z": z_f}, donors, train=train, frozen=mu[:, cols], classes=sp["readout_classes"],
                            folds=5, seed=int(d["seed"]), max_fit=sp["max_fit"])
        if first is not None:                                # the readouts fitted on cells of F are those of the first head
            study._models.update({k: v for k, v in first._models.items() if k[0] == "units"})
            study._own = first._own
        per = {}
        for i, u in enumerate(units):
            ti = np.where(donors == u)[0]
            side = scope._Side(study.Y[ti], None)
            Ct = scope._counts(side.E, d["bootstrap"], np.random.default_rng(d["seed"] + 1))
            per[u] = {}
            for c in VCELLS:
                P = study.frozen[ti] if c[0] == "frozen" else study.predict(c[0], c[1], "ridge", ti, (u,)).astype(np.float64)
                uu, _ = scope._u_variation(side, P, Ct)
                per[u][vkey(c)] = {"u": float(uu[0]), "ci": scope.ci(uu[1:], d["alpha"])}
            if i % 10 == 9:
                print(f"   {h}: {i + 1} donors, {time.time() - t0:.0f} s", flush=True)
        res["heads"][h] = {"head": hd.name, "n_fit": int(len(train["Y"])), "donors": per}
        first = first or study
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    return res


def summarise_variation(res, status_of):
    """Means over the donors of F and the contrasts of the round, from the donor-by-donor values."""
    donors = sorted(status_of)
    U = {h: {vkey(c): np.array([[res["heads"][h]["donors"][u][vkey(c)]["u"]] for u in donors]) for c in VCELLS} for h in HEADS}
    g = lambda a_, b_, rows=slice(None): {x: scope.arms_family({"a": a_[rows], "b": b_[rows]}, base="a")["arms"]["b"][x] for x in ("gain", "ci")}
    hl = np.array([status_of[u] == "healthy" for u in donors])
    groups = lambda a_, b_: {"all": g(a_, b_), "healthy": g(a_, b_, hl), "patients": g(a_, b_, ~hl)}
    return {"donors": len(donors), "seeds": [res["seed"]], "target": "variation",
            "note": "each donor against readouts not fitted on it; `units`: fitted on the other donors of F; one seed, so the intervals resample donors only",
            "variation": {
                "u": {h: {k: float(v.mean()) for k, v in U[h].items()} for h in HEADS},
                "ceiling_as_is": groups(U["as_is"]["frozen"], U["as_is"]["units@Z"]),
                "left_after_the_change": groups(U["coverage"]["frozen"], U["coverage"]["units@Z"]),
                "fall_from_as_is_to_coverage": {k: g(U["as_is"][k], U["coverage"][k]) for k in ("frozen", "train@RNA512", "train@Z")},
                "head_against_readout_of_its_training_cells": {h: g(U[h]["train@Z"], U[h]["frozen"]) for h in HEADS},
                "state_against_rna": {h: g(U[h]["units@RNA512"], U[h]["units@Z"]) for h in HEADS},
                "own_cells_after_other_donors": g(U["as_is"]["units@Z"], U["as_is"]["own@Z"])},
            "not_read": ["the shift target of the donor events (a readout for every pair of donors)", "the stages inside the head",
                         "the disease as an event: the run would resample the patients' cells one by one, not by donor"]}


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
    p.add_argument("--scope-runs", action="store_true",
                   help="the runs this script was first written for: scope runs with the donors and the disease as events. Not used: see donors_variation")
    p.add_argument("--stage", choices=("all", "prepare", "run"), default="all",
                   help="prepare: write the arrays and the SPEC of the runs and stop; run: read the grid of runs prepared before (a process that holds only "
                        "the arrays of its run, so that two runs fit in memory side by side); all: both, run after run")
    a = p.parse_args()
    res0 = common.RUN / "results"
    if not (res0 / "final.json").exists():
        raise SystemExit("the test has not been read; the grid on F comes after it")
    decl = json.loads((res0 / "gate5_declaration.json").read_text())
    read_final.check_declaration({"sha256": {k: v for k, v in decl["sha256"].items() if k in ("final_arrays", "final_z")}},
                                 {"final_arrays": a.final_arrays, "final_z": a.final_z})
    res = res0 / "grid_final"
    status_of = {d["donor"]: d["status"] for d in decl["donors"]}
    if not a.scope_runs:
        a.tmp.mkdir(parents=True, exist_ok=True)
        f = res / f"s{a.seeds[0]}_donors_variation.json"
        r = json.loads(f.read_text()) if f.exists() else donors_variation(a, decl, a.seeds[0], f)
        out = summarise_variation(r, status_of)
        (res0 / "what_is_left_final.json").write_text(json.dumps(out, indent=1))
        v = out["variation"]
        for h in HEADS:
            print(h, {k: round(x, 4) for k, x in v["u"][h].items()})
        print("ceiling (as_is):", {k: round(x["gain"], 4) for k, x in v["ceiling_as_is"].items()}, "left:", {k: round(x["gain"], 4) for k, x in v["left_after_the_change"].items()})
        return
    todo = [(f"s{s}_{kind}_{h}", s, kind, h) for s in a.seeds for kind in a.only for h in a.heads_only if kind == "donors" or s == a.seeds[0]]
    left = [t for t in todo if not (res / t[0] / "scope.json").exists()]
    if left and not a.summary_only and a.stage == "run":
        for name, s, kind, h in left:
            t0 = time.time()
            print(f"== {name}", flush=True)
            run_prepared(name, a.tmp, res)
            print(f"   done in {time.time() - t0:.0f} s", flush=True)
    elif left and not a.summary_only:
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
        rna = rna_view(a.bm_arrays, a.pbmc_arrays, a.final_arrays, wh["RNA512"], a.tmp / "rna512.npz")
        for name, s, kind, h in left:
            t0 = time.time()
            print(f"== {name}", flush=True)
            d = a.heads / f"s{s}_{HEADS[h]}"
            audit, train = arrays_of_run(dict(np.load(d / "head.npz")), dict(np.load(d / "fit_cells.npz")), F, rna, cols112, wb, wh, bm["split"])
            prepare(name, "donor" if kind == "donors" else "disease", audit, train, sorted(status_of), a.tmp)
            if a.stage == "prepare":
                print("   prepared", flush=True)
                continue
            del audit, train
            run_prepared(name, a.tmp, res)
            print(f"   done in {time.time() - t0:.0f} s", flush=True)
    if a.stage == "prepare":
        return
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
