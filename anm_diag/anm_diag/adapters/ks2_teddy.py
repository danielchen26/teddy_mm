"""Adapter: TEDDY-G 70M on Replogle 2022 K562 CRISPRi cells from the KS-2 pipeline.

Paths (a dict or YAML file):
  counts        ks2/pipe/out/counts.npz            (CSR counts, rows aligned with prep rows)
  var           ks2/meta/var.csv.gz                (gene symbols `_index`, `ensembl_id`)
  prep          prep.npz                           (P perturbation labels, ctrl, rows)
  cells         rank_need.npy                      (row indices of the cells to use)
  vocab         teddy_g_70M/vocab.txt
  medians       data/reference/teddy_gene_medians.json
  encoder_dir   directory holding teddy_encoder.py (TEDDY's official tokenisation code)
  z_files       [z_V.npz, z_controls.npz]          (stored KS-2 embeddings, arrays rows/z)
  layers        optional layers npz from layer_states.py (layer_means [12,n,512], window_means, rows)
  z_split       optional embeddings recomputed on the state half of split counts (for count_split)

Operations (node `op`):
  log1p_cp10k              raw view: log1p CP10k of genes detected in >=1% of controls
  teddy_official_tokens    official TEDDY tokens -> (2049 - position)/2048 per vocabulary gene, 0 if absent
  drop_order               membership of the official token sequence (order removed)
  teddy_embedding          stored gene-mean embedding Z
  teddy_layer_mean:<l>     gene-mean of encoder layer l (1..12), from `layers`
  teddy_window_mean:<l>    window means (positions 0-63, 64-255, 256+) of layer l (0..12), from `layers`
  readout_genes:<readout>  positive control: log1p CP10k of that readout's genes
  depth                    technical: log total counts, log detected genes
Readout statistic: mean_log1p_cp10k over `features` (gene regex), class targets removed.
"""
import csv
import gzip
import importlib.util
import json
import re

import numpy as np
from scipy import sparse

from .base import Adapter, file_sha256, load_paths

PAD = 43810


def _load_encoder(encoder_dir):
    spec = importlib.util.spec_from_file_location("teddy_encoder", f"{encoder_dir}/teddy_encoder.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class KS2TeddyAdapter(Adapter):
    def __init__(self, paths, log=print):
        paths = load_paths(paths)
        self.paths, self.log = paths, log
        c = np.load(paths["counts"])
        X = sparse.csr_matrix((c["data"], c["indices"], c["indptr"]), shape=tuple(c["shape"]))
        with gzip.open(paths["var"], "rt") as f:
            r = csv.reader(f); h = next(r); rows = list(r)
        self.sym = np.array([row[h.index("_index")] for row in rows])
        self.ens = np.array([row[h.index("ensembl_id")] for row in rows])
        d = np.load(paths["prep"], allow_pickle=True)
        self.need = np.load(paths["cells"])
        self.P_all = d["P"]
        self.rows_all = d["rows"]
        self.X = X[self.need].toarray().astype(np.float64)
        self.te = _load_encoder(paths["encoder_dir"])
        voc = {}
        with open(paths["vocab"]) as f:
            for i, line in enumerate(f):
                voc[line.rstrip("\n")] = i
        meds = self.te.load_gene_medians(paths["medians"])
        tok = np.array([voc.get(e, -1) for e in self.ens])
        self.mapped = np.where(tok >= 0)[0]
        self.token_ids = tok[self.mapped].astype(np.int64)
        self.fac = np.array([1.0 / meds[e] if e in meds else 1.0 for e in self.ens[self.mapped]])
        self._cache = {}
        self.split = None

    # ---- cells and source
    def cells(self, decl):
        src = decl["source"]
        self._rspec_by_id = {r["id"]: r for r in decl["readouts"]}
        P = self.P_all[self.need]
        ctrl = P == src["control"]
        tr = src["treated"]
        isX = np.array([bool(re.match(tr, s)) for s in P]) if isinstance(tr, str) else np.isin(P, tr)
        sel = isX | ctrl
        self.sel = np.where(sel)[0]
        self.Pn = P
        self.targets_of_class = sorted(set(P[isX].tolist()))
        A = isX[sel].astype(int)
        cluster = np.where(A == 1, P[sel], np.array([f"c{i}" for i in range(len(P))])[sel])
        if decl.get("measurement_sharing") == "count_split":
            rng = np.random.default_rng(int(decl["inference"]["seed"]) + 7)
            Xi = self.X.astype(np.int64)
            half = rng.binomial(Xi, 0.5)
            self.split = {"readout": half.astype(np.float64), "state": (Xi - half).astype(np.float64)}
        self._cache = {}
        return {"A": A, "cluster": cluster, "ids": self.rows_all[self.need][sel]}

    def _Xstate(self):
        return self.split["state"] if self.split else self.X

    def _Xread(self):
        return self.split["readout"] if self.split else self.X

    # ---- readouts
    def _readout_mask(self, rspec):
        pat = rspec["features"]
        m = np.array([bool(re.match(pat, s)) for s in self.sym])
        targeted = np.isin(self.sym, self.targets_of_class) if rspec.get("exclude_targets", True) else np.zeros_like(m)
        return m & ~targeted

    def readout(self, rspec, cells):
        if rspec["statistic"] != "mean_log1p_cp10k":
            raise ValueError(f"unsupported readout statistic {rspec['statistic']}")
        X = self._Xread()[self.sel]
        m = self._readout_mask(rspec)
        tot = np.maximum(X.sum(1), 1.0)
        return np.log1p(X[:, m] / tot[:, None] * 1e4).mean(1)

    # ---- states
    def _ctrl_rows(self):
        return self.Pn[self.sel] == "control"

    def _raw(self, node, cells, arg):
        X = self._Xstate()[self.sel]
        tot = np.maximum(X.sum(1), 1.0)
        keep = (X[self._ctrl_rows()] > 0).mean(0) >= 0.01
        return np.log1p(X[:, keep] / tot[:, None] * 1e4).astype(np.float32)

    def tokens(self, X):
        out_tok, out_len = [], []
        for s in range(0, X.shape[0], 256):
            Xm = np.ascontiguousarray(X[s:s + 256][:, self.mapped], dtype=np.float32)
            tok, att = self.te.rank_encode_official(self.te.official_values(Xm, self.fac), self.token_ids, pad_id=PAD)
            out_tok.append(tok); out_len.append(att.sum(1))
        return np.concatenate(out_tok), np.concatenate(out_len)

    def _token_matrix(self):
        if "T" not in self._cache:
            tok, L = self.tokens(self._Xstate()[self.sel])
            tid2m = {int(t): i for i, t in enumerate(self.token_ids)}
            T = np.zeros((len(tok), len(self.mapped)), np.float32)
            for i in range(len(tok)):
                ids = np.array([tid2m[int(t)] for t in tok[i, :L[i]]])
                T[i, ids] = (2049 - np.arange(1, len(ids) + 1)) / 2048
            keep = (T[self._ctrl_rows()] > 0).mean(0) >= 0.01
            self._cache["T"] = T[:, keep]
        return self._cache["T"]

    def _tokens(self, node, cells, arg):
        return self._token_matrix()

    def _set(self, node, cells, arg):
        return (self._token_matrix() > 0).astype(np.float32)

    def _z(self, node, cells, arg):
        if self.split:
            if not self.paths.get("z_split"):
                raise RuntimeError("count_split needs embeddings recomputed on the state half (paths.z_split); "
                                   "run layer_states.py on the split counts on the Mac")
            files = [self.paths["z_split"]]
        else:
            files = self.paths["z_files"]
        zr = {}
        for f in files:
            e = np.load(f)
            zr.update({int(r): v for r, v in zip(e["rows"], e["z"] if "z" in e.files else e["z_last"])})
        rows = self.rows_all[self.need][self.sel]
        return np.array([zr[int(r)] for r in rows], np.float32)

    def _layer(self, node, cells, arg, window=False):
        if not self.paths.get("layers"):
            raise RuntimeError("layer states need paths.layers (output of layer_states.py on the Mac)")
        e = np.load(self.paths["layers"])
        idx = {int(r): i for i, r in enumerate(e["rows"])}
        rows = self.rows_all[self.need][self.sel]
        ii = np.array([idx[int(r)] for r in rows])
        l = int(arg)
        if window:
            W = e["window_means"][l][:, ii]            # [windows, n, d]
            return np.concatenate(list(W), axis=1).astype(np.float32)
        return e["layer_means"][l - 1][ii].astype(np.float32)

    def _layer_mean(self, node, cells, arg):
        return self._layer(node, cells, arg)

    def _window_mean(self, node, cells, arg):
        return self._layer(node, cells, arg, window=True)

    def _readout_genes(self, node, cells, arg):
        X = self._Xstate()[self.sel]
        tot = np.maximum(X.sum(1), 1.0)
        m = self._readout_mask(self._rspec_by_id[arg])
        return np.log1p(X[:, m] / tot[:, None] * 1e4).astype(np.float32)

    def _depth(self, node, cells, arg):
        X = self._Xstate()[self.sel]
        return np.column_stack([np.log(np.maximum(X.sum(1), 1)), np.log(np.maximum((X > 0).sum(1), 1))]).astype(np.float32)

    ops = {"log1p_cp10k": _raw, "teddy_official_tokens": _tokens, "drop_order": _set, "teddy_embedding": _z,
           "teddy_layer_mean": _layer_mean, "teddy_window_mean": _window_mean, "readout_genes": _readout_genes,
           "depth": _depth}

    def bind(self, decl):
        self._rspec_by_id = {r["id"]: r for r in decl["readouts"]}
        return self

    # ---- exact same-state test (minimal family, official tokeniser)
    def exact(self, rspec, cells, delta, max_cells=None):
        X = self.X[self.sel]
        A = np.asarray(cells["A"])
        m = self._readout_mask(rspec)
        ctrl = A == 0
        cp = X / np.maximum(X.sum(1), 1.0)[:, None] * 1e4
        lo_cp, hi_cp = np.percentile(cp[ctrl], 1, axis=0), np.percentile(cp[ctrl], 99, axis=0)
        tid2m = {int(t): i for i, t in enumerate(self.token_ids)}
        direction = -np.sign(delta)
        idx = np.where(A == 1)[0][:max_cells] if max_cells else np.where(A == 1)[0]
        moves, same, resolved = [], [], []
        mm = m[self.mapped]
        for i in idx:
            x = X[i]
            xs, tok0 = self._partner(x, mm, direction, abs(delta), lo_cp, hi_cp, tid2m, m)
            t1, L1 = self.tokens(xs[None])
            same.append(bool(np.array_equal(t1[0, :L1[0]], tok0)))
            Ox = _O(x, m); Os = _O(xs, m)
            moves.append(abs(Ox - Os) / abs(delta))
            resolved.append(abs(Ox - Os) > 2 * _sdO(x, m))
        moves, same, resolved = map(np.array, (moves, same, resolved))
        return {"cells": int(len(idx)), "tokens_identical": round(float(same.mean()), 3),
                "median_move_over_delta": round(float(np.median(moves)), 3),
                "frac_full_delta_hidden": round(float(np.mean((moves >= 1) & same)), 3),
                "frac_resolved_beyond_noise": round(float(np.mean(resolved & same)), 3),
                "note": "minimal family: readout genes move between rank neighbours, vocabulary total fixed"}

    def _partner(self, x, mm, direction, target, lo_cp, hi_cp, tid2m, m, eps=1e-6):
        xm = x[self.mapped].copy(); fac = self.fac; N = x.sum()
        tok, L = self.tokens(x[None]); tok0 = tok[0, :L[0]].copy()
        T = np.array([tid2m[int(t)] for t in tok0])
        inT = np.zeros(len(xm), bool); inT[T] = True
        Rm = np.where((xm > 0) & ~inT)[0]
        mov = mm & (xm > 0)
        g = self.mapped
        lo = np.minimum(lo_cp[g] * N / 1e4, xm); hi = np.maximum(hi_cp[g] * N / 1e4, xm)
        xn = xm.copy()
        if direction < 0:
            for mi in Rm[mov[Rm]]:
                xn[mi] = min(xm[mi], max(1.0, np.ceil(lo[mi])))
            below = float((xn[Rm] * fac[Rm]).max()) if len(Rm) else 0.0
            for pos in range(len(T) - 1, -1, -1):
                mi = T[pos]
                if mov[mi]:
                    c = np.floor(below * (1 + eps) / fac[mi]) + 1.0
                    xn[mi] = min(xm[mi], max(c, np.ceil(lo[mi]), 1.0))
                below = xn[mi] * fac[mi]
            w_last = float((xn[T] * fac[T]).min())
            pool = np.where((xm > 0) & ~inT & ~mm)[0]
            cap = np.maximum(np.ceil(w_last * (1 - eps) / fac[pool]) - 1.0 - xm[pool], 0.0)
        else:
            above = np.inf
            for pos in range(len(T)):
                mi = T[pos]
                if mov[mi]:
                    c = np.floor(hi[mi]) if np.isinf(above) else min(np.ceil(above * (1 - eps) / fac[mi]) - 1.0, np.floor(hi[mi]))
                    xn[mi] = max(xm[mi], c)
                above = xn[mi] * fac[mi]
            for mi in Rm[mov[Rm]]:
                xn[mi] = max(xm[mi], min(np.ceil(above * (1 - eps) / fac[mi]) - 1.0, np.floor(hi[mi])))
            pool = np.where((xm > 0) & ~inT & ~mm)[0]
            cap = np.maximum(xm[pool] - 1.0, 0.0)
        xe = x.copy(); xe[g] = xn
        rnd = np.ceil if direction < 0 else np.floor
        step = xe - x

        def at(lam):
            xs = rnd(x + lam * step)
            D = float(x[g].sum() - xs[g].sum())
            if D == 0:
                return xs
            if cap.sum() < abs(D):
                return None
            left = abs(D); cols = g[pool]
            for j in np.argsort(-cap, kind="stable"):
                if left <= 0:
                    break
                a = min(cap[j], left); xs[cols[j]] += a if D > 0 else -a; left -= a
            return xs

        lam_hi = 1.0
        if at(1.0) is None:
            lo_, hi_ = 0.0, 1.0
            for _ in range(30):
                mid = 0.5 * (lo_ + hi_)
                if at(mid) is None: hi_ = mid
                else: lo_ = mid
            lam_hi = lo_
        O0 = _O(x, m)
        xs = at(lam_hi)
        if abs(_O(xs, m) - O0) > target:
            lo_, hi_ = 0.0, lam_hi
            for _ in range(40):
                mid = 0.5 * (lo_ + hi_)
                if abs(_O(at(mid), m) - O0) < target: lo_ = mid
                else: hi_ = mid
            xs = at(hi_)
        return xs, tok0

    def provenance(self):
        out = {"adapter": "ks2_teddy", "measurement_sharing_split": bool(self.split)}
        for k in ("counts", "prep", "cells", "vocab", "medians"):
            out[f"sha256_{k}"] = file_sha256(self.paths[k])[:16]
        out["sha256_teddy_encoder"] = file_sha256(f"{self.paths['encoder_dir']}/teddy_encoder.py")[:16]
        try:
            import torch
            out["torch"] = torch.__version__
        except Exception:
            pass
        return out


def _O(x, m):
    return float(np.log1p(x[m] / max(x.sum(), 1.0) * 1e4).mean())


def _sdO(x, m):
    N = max(x.sum(), 1.0); c = x[m] / N * 1e4
    return float(np.sqrt(((1e4 / N) ** 2 * (x[m] + 0.5) / (1 + c) ** 2).sum()) / m.sum())
