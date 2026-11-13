#!/usr/bin/env python3
"""On donors the head has not seen: where is the protein signal that the head does not deliver?

No model is trained. Ridge readouts are fitted on fixed views of each cell (principal components of the RNA,
TEDDY's pooled state) and scored on cells they were not fitted on. The frozen head is scored on the same cells.
The design, the readings and the limits are in ../DESIGN.md, committed before this script was run on data.

  python gate0.py --bm-arrays cite_arrays.npz --bm-z z_rna.npy \\
      --ext-arrays ext_cite_arrays.npz --ext-z ext_z_rna.npy --head-pred head_pred.npz --out ../results

Inputs: the prepared arrays of scripts/02_prepare_cite.py (RNA counts, measured protein, donor, site, split), the
pooled TEDDY state of scripts/03_teddy_embed.py, and the frozen head's predictions for every cell.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

ALPHAS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)   # ridge penalty as a multiple of the mean eigenvalue of the Gram matrix
N_FOLDS = 5
SEED = 20261112
INNER_CELLS = 3000       # cells per group used to choose the penalty
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ---------------------------------------------------------------- data
def load(arrays: Path, z: Path) -> dict:
    a = np.load(arrays, allow_pickle=True)
    rna = sp.csr_matrix((a["rna_data"].astype(np.float32), a["rna_indices"], a["rna_indptr"]),
                        shape=tuple(int(v) for v in a["rna_shape"]))
    d = dict(rna=rna, adt=a["adt"].astype(np.float64), donors=a["donors"].astype(str), sites=a["sites"].astype(str),
             split=a["split"].astype(str), names=a["adt_names"].astype(str), genes=a["rna_names"].astype(str),
             z=np.load(z).astype(np.float32))
    d["present"] = a["adt_present"].astype(bool) if "adt_present" in a.files else np.ones(d["adt"].shape[1], bool)
    assert d["z"].shape[0] == d["adt"].shape[0] == rna.shape[0]
    return d


def lognorm_inplace(X: sp.csr_matrix) -> sp.csr_matrix:
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1.0
    X.data *= np.repeat((1e4 / tot).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def standardize(X: np.ndarray, ref: np.ndarray) -> np.ndarray:
    mu = X[ref].mean(0, dtype=np.float64)
    sd = X[ref].std(0, dtype=np.float64)
    sd[sd < 1e-8] = 1.0
    return ((X - mu) / sd).astype(np.float32)


def add_rff(X: np.ndarray, ref: np.ndarray, n_rff: int, rng: np.random.Generator) -> np.ndarray:
    """[X, random Fourier features of X]: a ridge on these is a kernel ridge (Gaussian kernel) plus a linear term."""
    sub = X[rng.choice(np.where(ref)[0], size=min(2000, int(ref.sum())), replace=False)].astype(np.float64)
    sq = (sub ** 2).sum(1)
    d2 = sq[:, None] + sq[None, :] - 2 * sub @ sub.T
    med = float(np.median(d2[np.triu_indices(len(sub), 1)]))
    gamma = 1.0 / med
    omega = rng.normal(0.0, np.sqrt(2 * gamma), size=(X.shape[1], n_rff)).astype(np.float32)
    phase = rng.uniform(0, 2 * np.pi, size=n_rff).astype(np.float32)
    out = np.empty((X.shape[0], n_rff), np.float32)
    for i in range(0, X.shape[0], 20000):
        out[i:i + 20000] = np.sqrt(2.0 / n_rff) * np.cos(X[i:i + 20000] @ omega + phase)
    return np.hstack([X, standardize(out, ref)])


# ---------------------------------------------------------------- sufficient statistics and ridge
class Stats:
    __slots__ = ("n", "sx", "sy", "G", "C")

    def __init__(self, n, sx, sy, G, C):
        self.n, self.sx, self.sy, self.G, self.C = n, sx, sy, G, C

    @classmethod
    def of(cls, F: np.ndarray, Y: np.ndarray):
        F = F.astype(np.float64)
        return cls(F.shape[0], F.sum(0), Y.sum(0), F.T @ F, F.T @ Y)

    def __add__(self, o):
        return Stats(self.n + o.n, self.sx + o.sx, self.sy + o.sy, self.G + o.G, self.C + o.C)

    def __sub__(self, o):
        return Stats(self.n - o.n, self.sx - o.sx, self.sy - o.sy, self.G - o.G, self.C - o.C)


def total(items):
    items = list(items)
    out = items[0]
    for s in items[1:]:
        out = out + s
    return out


class Ridge:
    """Ridge with an intercept, from sufficient statistics; one eigendecomposition serves every penalty."""

    def __init__(self, st: Stats, d: int):
        n = st.n
        self.d = d
        self.xbar = st.sx[:d] / n
        self.ybar = st.sy / n
        G = st.G[:d, :d] - np.outer(st.sx[:d], st.sx[:d]) / n
        C = st.C[:d] - np.outer(st.sx[:d], st.sy) / n
        e, V = np.linalg.eigh(G)
        self.e = np.clip(e, 0.0, None)
        self.V = V
        self.VC = V.T @ C
        self.scale = float(self.e.sum() / d)

    def predict(self, F: np.ndarray, alpha: float) -> np.ndarray:
        W = self.V @ (self.VC / (self.e + alpha * self.scale)[:, None])
        return (F[:, :self.d].astype(np.float64) - self.xbar) @ W + self.ybar


def pearson(Y: np.ndarray, P: np.ndarray) -> np.ndarray:
    y = Y - Y.mean(0)
    p = P - P.mean(0)
    den = np.sqrt((y * y).sum(0) * (p * p).sum(0))
    ok = den > 0
    out = np.full(Y.shape[1], np.nan)
    out[ok] = (y * p).sum(0)[ok] / den[ok]
    return out


# ---------------------------------------------------------------- the four fits for one view and one class
class View:
    def __init__(self, name: str, F: np.ndarray, Y: np.ndarray, meta: dict, rng: np.random.Generator):
        self.name, self.F, self.Y, self.m = name, F, Y, meta
        self.groups = meta["groups"]
        self.block = {}
        self.gstat = {}
        for g in self.groups:
            for f in range(N_FOLDS):
                idx = g["folds"][f]
                self.block[(g["id"], f)] = Stats.of(F[idx], Y[idx])
            self.gstat[g["id"]] = total(self.block[(g["id"], f)] for f in range(N_FOLDS))
        self.inner = {g["id"]: (g["idx"] if g["n"] <= INNER_CELLS else np.sort(rng.choice(g["idx"], INNER_CELLS, replace=False)))
                      for g in self.groups}
        log(f"{name}: statistics for {len(self.block)} blocks, {F.shape[1]} features")

    # ---- helpers
    def cols(self, g):
        return self.m["cols"][g["dataset"]]

    def score_groups(self, model: Ridge, groups: list, cols: np.ndarray) -> np.ndarray:
        """mean over groups of the mean per-protein Pearson, for each penalty"""
        out = np.zeros(len(ALPHAS))
        for g in groups:
            idx = self.inner[g["id"]]
            for a, alpha in enumerate(ALPHAS):
                out[a] += np.nanmean(pearson(self.Y[idx][:, cols], model.predict(self.F[idx], alpha)[:, cols]))
        return out / max(len(groups), 1)

    def across(self, d: int, fit_groups: list, targets: list, cache: dict) -> dict:
        """Fit on fit_groups, penalty chosen by leaving one donor of the fit set out in turn; predict targets."""
        cols = self.m["cols"][fit_groups[0]["dataset"]]
        ids = tuple(sorted(g["id"] for g in fit_groups))
        donors = sorted({g["donor"] for g in fit_groups})
        if len(donors) < 3:
            return {}

        def model(key_ids):
            key = (d, key_ids)
            if key not in cache:
                cache[key] = Ridge(total(self.gstat[i] for i in key_ids), d)
            return cache[key]

        sc = np.zeros(len(ALPHAS))
        for dn in donors:
            rest = tuple(sorted(g["id"] for g in fit_groups if g["donor"] != dn))
            held = [g for g in fit_groups if g["donor"] == dn]
            sc += self.score_groups(model(rest), held, cols)
        sc /= len(donors)
        a = int(np.argmax(sc))
        full = model(ids)
        return {t["id"]: dict(pred=full.predict(self.F[t["idx"]], ALPHAS[a]).astype(np.float32), alpha=ALPHAS[a],
                              inner=float(sc[a]), fit_cells=int(sum(g["n"] for g in fit_groups)),
                              fit_donors=len(donors)) for t in targets}

    def within(self, d: int, g: dict) -> dict:
        """Cross-fitted within one group; the penalty is chosen inside each fit set (nested)."""
        gid, cols = g["id"], self.cols(g)
        cache = {}

        def model(excl):
            key = frozenset(excl)
            if key not in cache:
                cache[key] = Ridge(self.gstat[gid] - total(self.block[(gid, f)] for f in excl), d)
            return cache[key]

        pred = np.empty((g["n"], self.Y.shape[1]), np.float32)
        pos = {int(c): i for i, c in enumerate(g["idx"])}
        inner, alphas = [], []
        for k in range(N_FOLDS):
            sc = np.zeros(len(ALPHAS))
            for j in range(N_FOLDS):
                if j == k:
                    continue
                idx = g["folds"][j]
                mdl = model((k, j))
                for a, alpha in enumerate(ALPHAS):
                    sc[a] += np.nanmean(pearson(self.Y[idx][:, cols], mdl.predict(self.F[idx], alpha)[:, cols]))
            sc /= N_FOLDS - 1
            a = int(np.argmax(sc))
            idx = g["folds"][k]
            pred[[pos[int(c)] for c in idx]] = model((k,)).predict(self.F[idx], ALPHAS[a])
            inner.append(float(sc[a]))
            alphas.append(ALPHAS[a])
        return dict(pred=pred, alpha=float(np.median(alphas)), inner=float(np.mean(inner)),
                    fit_cells=int(g["n"] * (N_FOLDS - 1) / N_FOLDS), fit_donors=1)

    def run(self, d: int, label: str) -> dict:
        """returns {fit: {group id: {...}}} for the fits T, S, L, C"""
        G = self.groups
        bm_train = [g for g in G if g["dataset"] == "bm" and g["split"] == "train"]
        out = {k: {} for k in "TSLC"}
        cache = {}
        # T: what the head was trained on. For a training group, without its donor.
        unseen = [g for g in G if not (g["dataset"] == "bm" and g["split"] == "train")]
        out["T"].update(self.across(d, bm_train, unseen, cache))
        for dn in sorted({g["donor"] for g in bm_train}):
            out["T"].update(self.across(d, [g for g in bm_train if g["donor"] != dn],
                                        [g for g in bm_train if g["donor"] == dn], cache))
        del cache
        log(f"{self.name}/{label}: T done")
        # S: as T, and also without the group's site (bone marrow only; equals T where the site is already absent)
        for g in [g for g in G if g["dataset"] == "bm"]:
            fit = [h for h in bm_train if h["donor"] != g["donor"] and h["site"] != g["site"]]
            same = sorted(h["id"] for h in fit) == sorted(h["id"] for h in bm_train
                                                         if not (g["split"] == "train" and h["donor"] == g["donor"]))
            if same:
                out["S"][g["id"]] = out["T"][g["id"]]
            else:
                out["S"].update(self.across(d, fit, [g], {}))
        for g in [g for g in G if g["dataset"] != "bm"]:
            out["S"][g["id"]] = out["T"][g["id"]]
        log(f"{self.name}/{label}: S done")
        # L: the other donors of the same study
        for ds in sorted({g["dataset"] for g in G}):
            cache_l = {}
            for dn in sorted({g["donor"] for g in G if g["dataset"] == ds}):
                out["L"].update(self.across(d, [g for g in G if g["dataset"] == ds and g["donor"] != dn],
                                            [g for g in G if g["dataset"] == ds and g["donor"] == dn], cache_l))
        log(f"{self.name}/{label}: L done")
        # C: the group's own cells
        for g in G:
            out["C"][g["id"]] = self.within(d, g)
        log(f"{self.name}/{label}: C done")
        return out


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bm-arrays", type=Path, required=True)
    p.add_argument("--bm-z", type=Path, required=True)
    p.add_argument("--ext-arrays", type=Path, required=True)
    p.add_argument("--ext-z", type=Path, required=True)
    p.add_argument("--head-pred", type=Path, required=True, help="npz with bm_pred and ext_pred: the frozen head's output")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--svd", type=int, default=512, help="principal components of the RNA (the RAW view)")
    p.add_argument("--rff", type=int, default=1024, help="random Fourier features added in the non-linear class")
    p.add_argument("--boot", type=int, default=300)
    p.add_argument("--min-cells", type=int, default=300)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    bm, ext = load(args.bm_arrays, args.bm_z), load(args.ext_arrays, args.ext_z)
    assert list(bm["names"]) == list(ext["names"]) and list(bm["genes"]) == list(ext["genes"])
    hp = np.load(args.head_pred, allow_pickle=True)
    head = np.vstack([hp["bm_pred"], hp["ext_pred"]]).astype(np.float64)
    nb, ne = bm["adt"].shape[0], ext["adt"].shape[0]
    N, P = nb + ne, bm["adt"].shape[1]
    assert head.shape == (N, P)
    dataset = np.array(["bm"] * nb + ["ext"] * ne)
    donor = np.concatenate([bm["donors"], ext["donors"]])
    site = np.concatenate([bm["sites"], ext["sites"]])
    split = np.concatenate([bm["split"], np.array(["external"] * ne)])
    Y = np.vstack([bm["adt"], np.where(np.isfinite(ext["adt"]), ext["adt"], 0.0)])
    shared = np.where(bm["present"] & ext["present"])[0]
    cols = {"bm": np.where(bm["present"])[0], "ext": np.where(ext["present"])[0]}
    ref = (dataset == "bm") & (split == "train")
    log(f"cells: bone marrow {nb}, external {ne}; proteins {P}, shared {len(shared)}; head-training cells {int(ref.sum())}")

    # groups: bone marrow donor x site (x split), external donor
    groups = []
    keys = sorted({(dataset[i], donor[i], site[i] if dataset[i] == "bm" else "", split[i]) for i in range(N)})
    for ds, dn, st, sp_ in keys:
        idx = np.where((dataset == ds) & (donor == dn) & (split == sp_) & ((site == st) if ds == "bm" else True))[0]
        if len(idx) < args.min_cells:
            log(f"skipped {ds} {dn} {st} {sp_}: {len(idx)} cells")
            continue
        perm = rng.permutation(idx)
        groups.append(dict(id=f"{ds}:{dn}:{st}" if ds == "bm" else f"{ds}:{dn}", dataset=ds, donor=dn, site=st, split=sp_,
                           n=int(len(idx)), idx=idx, folds=[np.sort(perm[f::N_FOLDS]) for f in range(N_FOLDS)],
                           seen_by_head=bool(ds == "bm" and sp_ == "train")))
    assert len({g["id"] for g in groups}) == len(groups)
    meta = dict(groups=groups, cols=cols)
    for g in groups:
        log(f"group {g['id']:24s} {g['split']:8s} n={g['n']}")

    # views
    z = np.vstack([bm["z"], ext["z"]])
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)          # as the head reads it
    Xz = standardize(z, ref)
    del z
    from sklearn.decomposition import TruncatedSVD
    R = sp.vstack([lognorm_inplace(bm.pop("rna")), lognorm_inplace(ext.pop("rna"))]).tocsr()
    svd = TruncatedSVD(n_components=min(args.svd, R.shape[1] - 1), algorithm="randomized", n_iter=2, random_state=SEED)
    Xr = standardize(svd.fit_transform(R).astype(np.float32), ref)
    log(f"RAW view: {Xr.shape[1]} principal components of log-normalized counts, "
        f"{svd.explained_variance_ratio_.sum():.3f} of the variance")
    del R

    preds = {}        # (view, fit, class) -> {group id: dict}
    for vname, X in (("Z", Xz), ("RAW", Xr)):
        F = add_rff(X, ref, args.rff, np.random.default_rng(SEED + (1 if vname == "Z" else 2)))
        V = View(vname, F, Y, meta, np.random.default_rng(SEED + 7))
        for label, d in (("ridge", X.shape[1]), ("rff", F.shape[1])):
            res = V.run(d, label)
            for fit in res:
                preds[(vname, fit, label)] = res[fit]
        del V, F

    # ---------------------------------------------------------------- scores, intervals, terms, donor shifts
    brng = np.random.default_rng(SEED + 99)
    out_groups, per_protein = [], {}
    grp_by_id = {g["id"]: g for g in groups}
    for g in groups:
        gid, idx = g["id"], g["idx"]
        cg = np.intersect1d(cols[g["dataset"]], shared)
        y = Y[idx]
        series = {"A": head[idx]}
        info = {}
        for vname in ("Z", "RAW"):
            for fit in "TSLC":
                cand = {lab: preds[(vname, fit, lab)].get(gid) for lab in ("ridge", "rff")}
                cand = {k: v for k, v in cand.items() if v is not None}
                if not cand:
                    continue
                best = max(cand, key=lambda k: cand[k]["inner"])
                series[f"{vname}_{fit}"] = cand[best]["pred"].astype(np.float64)
                info[f"{vname}_{fit}"] = dict(chosen=best, alpha=cand[best]["alpha"], fit_cells=cand[best]["fit_cells"],
                                              fit_donors=cand[best]["fit_donors"],
                                              by_class={k: float(np.nanmean(pearson(y[:, cg], v["pred"][:, cg].astype(np.float64))))
                                                        for k, v in cand.items()})
        pp = {k: pearson(y[:, cg], v[:, cg]) for k, v in series.items()}
        point = {k: float(np.nanmean(v)) for k, v in pp.items()}
        # the donor shift: this group against the other donors of the same study
        others = np.where((dataset == g["dataset"]) & (donor != g["donor"]))[0]
        ybar_o, sd_o = Y[others][:, cg].mean(0), Y[others][:, cg].std(0)
        sd_o[sd_o < 1e-8] = 1.0
        ref_pred = {"A": head[others][:, cg].mean(0)}
        for vname in ("Z", "RAW"):
            for fit in ("T", "L"):
                k = f"{vname}_{fit}"
                if k not in series:
                    continue
                if fit == "L":
                    ref_pred[k] = ybar_o        # a ridge with an intercept reproduces the mean of its fit set
                else:
                    parts = [preds[(vname, fit, info[k]["chosen"])].get(h["id"]) for h in groups
                             if h["dataset"] == g["dataset"] and h["donor"] != g["donor"]]
                    ref_pred[k] = np.average(np.vstack([q["pred"][:, cg].mean(0) for q in parts if q is not None]), axis=0,
                                             weights=[len(q["pred"]) for q in parts if q is not None])
        d_meas = (y[:, cg].mean(0) - ybar_o) / sd_o

        def shift_u(mean_pred, key):
            d_pred = (mean_pred - ref_pred[key]) / sd_o
            return float(((d_meas - d_pred) ** 2).sum() / max((d_meas ** 2).sum(), 1e-12)), float(np.corrcoef(d_meas, d_pred)[0, 1])

        shift = {k: dict(zip(("u", "r"), shift_u(series[k][:, cg].mean(0), k))) for k in ref_pred}
        # bootstrap over the group's cells, the same resamples for every series
        B = args.boot
        bs = {k: np.empty(B) for k in series}
        bu = {k: np.empty(B) for k in ref_pred}
        ycg = y[:, cg]
        scg = {k: v[:, cg] for k, v in series.items()}
        for b in range(B):
            r = brng.integers(0, g["n"], g["n"])
            yb = ycg[r]
            for k, v in scg.items():
                vb = v[r]
                bs[k][b] = np.nanmean(pearson(yb, vb))
                if k in ref_pred:
                    dm = (yb.mean(0) - ybar_o) / sd_o
                    dp = (vb.mean(0) - ref_pred[k]) / sd_o
                    bu[k][b] = ((dm - dp) ** 2).sum() / max((dm ** 2).sum(), 1e-12)

        def ci(a):
            return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]

        scores = {k: dict(mean=point[k], ci=ci(bs[k]), **info.get(k, {})) for k in series}
        terms = {}

        def term(name, a, b):
            if a in series and b in series:
                terms[name] = dict(value=point[a] - point[b], ci=ci(bs[a] - bs[b]))

        term("readout_minus_head", "Z_T", "A")          # same training cells, a ridge readout instead of the head
        term("unseen_site_cost", "Z_T", "Z_S")          # what having the site in the fit is worth (bone marrow)
        term("same_study_donors", "Z_L", "Z_T")         # what the other donors of the study add
        term("own_cells", "Z_C", "Z_L")                 # what only the donor's own cells give
        term("ceiling_minus_head", "Z_C", "A")
        term("raw_minus_state_own", "RAW_C", "Z_C")     # readable from the RNA components and not from the state
        term("raw_minus_state_T", "RAW_T", "Z_T")
        term("raw_minus_state_L", "RAW_L", "Z_L")
        for k in shift:
            shift[k]["u_ci"] = ci(bu[k])
        call = cols[g["dataset"]]
        out_groups.append(dict(id=gid, dataset=g["dataset"], donor=g["donor"], site=g["site"], split=g["split"], n=g["n"],
                               seen_by_head=g["seen_by_head"], n_proteins=int(len(cg)),
                               head_all_measured_proteins=float(np.nanmean(pearson(y[:, call], head[idx][:, call]))),
                               scores=scores, terms=terms,
                               shift=dict(size_rms_sd=float(np.sqrt((d_meas ** 2).mean())), by=shift)))
        per_protein[gid] = {k: [None if not np.isfinite(x) else round(float(x), 4) for x in v] for k, v in pp.items()}
        log(f"scored {gid}: head {point['A']:.3f}" + "".join(f"  {k} {point[k]:.3f}" for k in ("Z_T", "Z_L", "Z_C", "RAW_C") if k in point))

    res = dict(
        what="per group: the frozen head and ridge readouts of two views, scored by the mean per-protein Pearson "
             "correlation within the group; no model was trained",
        settings=dict(alphas=ALPHAS, folds=N_FOLDS, seed=SEED, svd=int(Xr.shape[1]), rff=args.rff, boot=args.boot,
                      inner_cells=INNER_CELLS, classes=["ridge", "rff"],
                      fits=dict(T="the head's training cells (for a training group: without its donor)",
                                S="as T, also without the group's site",
                                L="the other donors of the same study", C="the group's own cells, 5-fold")),
        proteins_shared=[str(bm["names"][j]) for j in shared], n_cells=dict(bm=nb, ext=ne), groups=out_groups)
    (args.out / "gate0.json").write_text(json.dumps(res, indent=1))
    (args.out / "gate0_per_protein.json").write_text(json.dumps(dict(proteins=[str(bm["names"][j]) for j in shared],
                                                                    note="external groups list the shared proteins; "
                                                                         "bone-marrow groups list the same ones",
                                                                    groups=per_protein)))
    log("wrote", args.out / "gate0.json")


if __name__ == "__main__":
    main()
