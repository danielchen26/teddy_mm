"""Adapter that serves features exported by the real perturb_cite adapter for the pooled selection (same controls
for every declaration: the control subsample depends only on the control cells and the seed)."""
import json, os, numpy as np, scipy.sparse as sp
from anm_diag.adapters.base import Adapter

class PrecomputedAdapter(Adapter):
    def __init__(self, feat_dir):
        self.dir = feat_dir; d = np.load(os.path.join(feat_dir, "dense.npz")); self.d = {k: d[k] for k in d.files}; self._n = {}; self.sel = None
    def cells(self, decl):
        src = decl["source"]; tgt = self.d["target"].astype(str); tr = src["treated"]
        isX = np.isin(tgt, tr if isinstance(tr, list) else [tr]); ctrl = tgt == "control"
        keep = np.ones(len(tgt), bool); st = src.get("stratum")
        if st is not None and "stratum" in self.d:      # a declaration may restrict to some strata (e.g. one replicate)
            keep = np.isin(self.d["stratum"].astype(str), [str(v) for v in (st if isinstance(st, list) else [st])])
        sel = np.where(keep & (isX | ctrl))[0]; self.sel = sel; A = isX[sel].astype(int)
        unit = self.d["guide"].astype(str) if src.get("cluster", "target") == "guide" else tgt
        cl = np.where(A == 1, unit[sel], np.array([f"c{g}" for g in self.d["global_pos"]])[sel])
        return {"A": A, "cluster": cl, "ids": self.d["ids"].astype(str)[sel]}
    def readout(self, rspec, cells):
        return self.d["readout__" + rspec["id"]][self.sel]
    def state(self, node, cells):
        k = node["id"]
        if k not in self._n:
            p = os.path.join(self.dir, f"node__{k}.npz"); q = os.path.join(self.dir, f"node__{k}.npy")
            self._n[k] = sp.load_npz(p).toarray().astype(np.float32) if os.path.exists(p) else (np.load(q) if os.path.exists(q) else self.d["node__" + k])
        return self._n[k][self.sel]
    def provenance(self):
        return {"adapter": "precomputed(perturb_cite features exported on the source machine)"}
