"""Adapter: TEDDY-G 70M + the teddy_mm phase-1 protein head on perturbation data with measured protein
(Perturb-CITE-seq, Frangieh et al. 2021; ECCITE-seq, Papalexi et al. 2021; or any AnnData pair).

The readout is measured protein (ADT): a separate measurement that never enters TEDDY or the head, so
the observation is independent of every candidate state (no shared counting noise).

Paths (dict or YAML):
  rna_h5ad       AnnData with RNA counts in .X (or layers[rna_layer]); var has Ensembl ids in `ensembl_col`
  adt_h5ad       AnnData with ADT counts, same cells (obs_names) as the RNA file
  obs_target     obs column with the targeted gene ("NT"/"control" for non-targeting)
  obs_stratum    optional obs column used to restrict cells (e.g. condition), value(s) in declaration.source.stratum
  obs_guide      obs column with the detected guide(s), used when declaration.source.cluster is guide or
                 declaration.source.assign_from_guides is true
  guide_sep      separator between guides in obs_guide (default ";")
  guide_gene_regex  regex whose first group is a guide's gene, e.g. "^(.+)_\\d+$" for IFNGR1_3, "^(.+?)g\\d+$"
                 for IFNGR1g3 (needed with assign_from_guides)
  symbol_map     CSV symbol,ensembl used instead of a var column (tools/gene_map.py)
  vocab, medians, encoder_dir     as in the KS-2 adapter (official tokenisation)
  teddy_states   npz from layer_states.py run on the exported tokens (z_last, layer_means, window_means, rows)
  head_ckpt      teddy_mm outputs/cite_phase1_official/best.pt  (keys mlp, dec)
  head_adt_names optional list or file of the head's 134 protein names (for the report)

Workflow:
  1. adapter.export_tokens(decl, "tokens.npz") here (CPU, official tokeniser)
  2. on the Mac: python layer_states.py tokens.npz teddy_states.npz        (TEDDY forward, MPS)
  3. run(decl, adapter) with paths.teddy_states set

Operations: log1p_cp10k, teddy_official_tokens, drop_order, teddy_embedding, teddy_layer_mean:<l>,
teddy_window_mean:<l>, head_latent, head_output, rna_of_gene:<SYMBOL> (probe: the mRNA of a protein's gene),
adt_readout:<readout id> (positive control: the measured protein itself).
Readout statistic over `features` (protein names), mean over the listed proteins: adt_clr (centred log-ratio
across the panel) or adt_log1p (log1p counts). Only the selected cells' RNA rows are read (backed file), and
source.max_control_cells (with control_seed) declares a seeded subsample of the controls for compute.
"""
import re

import numpy as np

from .base import Adapter, file_sha256, load_paths
from .ks2_teddy import PAD, _load_encoder


class PerturbCiteAdapter(Adapter):
    def __init__(self, paths, log=print):
        paths = load_paths(paths)
        import anndata as ad
        self.paths, self.log = paths, log
        # RNA stays on disk (backed) until cells() knows which rows the declaration selects: the full
        # Perturb-CITE matrix (about 218k cells x 23k genes) does not fit in memory as a dense array.
        self.rna = ad.read_h5ad(paths["rna_h5ad"], backed="r")
        adt = ad.read_h5ad(paths["adt_h5ad"])
        common = self.rna.obs_names.intersection(adt.obs_names, sort=False)
        self.rna_pos = np.asarray(self.rna.obs_names.get_indexer(common))
        adt_pos = np.asarray(adt.obs_names.get_indexer(common))
        self.obs = self.rna.obs.iloc[self.rna_pos]
        self.obs_names = np.asarray(common).astype(str)
        Xa = adt.X[adt_pos]
        self.ADT_all = np.asarray(Xa.todense() if hasattr(Xa, "todense") else Xa, dtype=np.float64)
        self.adt_names = np.asarray(adt.var_names).astype(str)
        self.sym = np.asarray(self.rna.var_names).astype(str)
        col = paths.get("ensembl_col", "gene_ids")
        if col in self.rna.var.columns:
            self.ens = np.asarray(self.rna.var[col]).astype(str)
        elif paths.get("symbol_map"):
            # CSV with columns symbol,ensembl (e.g. exported from Ensembl BioMart); unmapped genes never reach TEDDY
            import csv
            with open(paths["symbol_map"]) as f:
                mp = {row["symbol"]: row["ensembl"] for row in csv.DictReader(f)}
            self.ens = np.array([mp.get(s, "") for s in self.sym])
            log(f"symbol_map: {int((self.ens != '').sum())} of {len(self.sym)} genes mapped to Ensembl")
        else:
            raise KeyError(f"RNA var has no '{col}' column; give paths.symbol_map (CSV symbol,ensembl)")
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
        self.sel = None

    # ---- cells and source
    def _targets(self, src):
        """Per-cell target label: the obs_target column as it is, or, with source.assign_from_guides, a label
        assigned from the detected guides (paths: obs_guide, guide_sep, guide_gene_regex whose first group is
        the guide's gene). "control": every guide is a control guide (source.control_guides); gene G: every
        guide targets G or is a control guide; "" (left out): no guide detected, or guides for two genes."""
        obs = self.obs
        if not src.get("assign_from_guides"):
            return obs[self.paths["obs_target"]].astype(str).to_numpy()
        sep = self.paths.get("guide_sep", ";")
        rx = re.compile(self.paths["guide_gene_regex"])
        ctrl = {str(g) for g in src["control_guides"]}
        out = []
        for g in obs[self.paths["obs_guide"]].to_numpy(dtype=object):
            if not isinstance(g, str) or g in ("nan", "", "None", "NaN"):   # no guide detected (missing value)
                out.append("")
                continue
            genes = set()
            for x in g.split(sep):
                x = x.strip()
                m = rx.match(x)
                genes.add(m.group(1) if m else x)
            tg = genes - ctrl
            out.append("control" if not tg else (next(iter(tg)) if len(tg) == 1 else ""))
        return np.array(out)

    def _select(self, decl):
        src = decl["source"]
        obs = self.obs
        tgt = self._targets(src)
        keep = tgt != ""
        st = src.get("stratum")
        if st is not None and self.paths.get("obs_stratum"):
            vals = [str(v) for v in (st if isinstance(st, list) else [st])]
            keep &= np.isin(obs[self.paths["obs_stratum"]].astype(str).to_numpy(), vals)
        ctrl = np.isin(tgt, src["control"] if isinstance(src["control"], list) else [src["control"]])
        tr = src["treated"]
        isX = np.isin(tgt, tr) if isinstance(tr, list) else np.array([bool(re.match(tr, t)) for t in tgt])
        sel = keep & (isX | ctrl)
        cap = src.get("max_control_cells")
        if cap and (sel & ctrl).sum() > int(cap):
            # a declared, seeded subsample of the controls (compute budget), fixed before any outcome is read
            rng = np.random.default_rng(int(src.get("control_seed", 0)))
            c = np.where(sel & ctrl)[0]
            drop = np.setdiff1d(c, rng.choice(c, int(cap), replace=False))
            sel[drop] = False
        return np.where(sel)[0], isX, tgt

    def _read_rows(self, pos):
        """Rows `pos` (positions among the common cells) of the RNA counts as CSR float32, in the order given."""
        import scipy.sparse as sp
        rows = self.rna_pos[pos]
        order = np.argsort(rows, kind="stable")
        src = self.rna.layers[self.paths["rna_layer"]] if self.paths.get("rna_layer") else self.rna.X
        out = []
        for s in range(0, len(order), 4096):            # h5py reads need increasing indices
            blk = src[rows[order[s:s + 4096]]]
            out.append(sp.csr_matrix(blk, dtype=np.float32))
        M = sp.vstack(out).tocsr() if out else sp.csr_matrix((0, len(self.sym)), dtype=np.float32)
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        return M[inv]

    def cells(self, decl):
        self._rspec_by_id = {r["id"]: r for r in decl["readouts"]}
        sel, isX, tgt = self._select(decl)
        if self.sel is None or len(sel) != len(self.sel) or (sel != self.sel).any():
            self.sel = sel
            self.Xs = self._read_rows(sel)               # counts of the selected cells only (sparse)
            self.ADT = self.ADT_all[sel]
            self._cache = {}
        A = isX[sel].astype(int)
        unit = tgt
        if decl["source"].get("cluster", "target") == "guide":     # resampling and folds by guide
            unit = self.obs[self.paths["obs_guide"]].astype(str).to_numpy()
        cluster = np.where(A == 1, unit[sel], np.array([f"c{i}" for i in range(len(tgt))])[sel])
        self._ctrl = A == 0
        return {"A": A, "cluster": cluster, "ids": self.obs_names[sel]}

    # ---- readout: measured protein
    def _clr(self):
        L = np.log1p(self.ADT)
        return L - L.mean(1, keepdims=True)

    def _adt_matrix(self, statistic):
        if statistic == "adt_clr":
            return self._clr()
        if statistic == "adt_log1p":
            return np.log1p(self.ADT)
        raise ValueError(f"unsupported readout statistic {statistic}")

    def readout(self, rspec, cells):
        """adt_clr: centred log-ratio across the whole panel (removes per-cell ADT depth, but on a small panel a
        drop in some proteins shifts the others: a compositional effect). adt_log1p: log1p counts, no
        cross-protein normalisation (keeps ADT depth). Declare which one before sealing."""
        names = rspec["features"] if isinstance(rspec["features"], list) else [rspec["features"]]
        idx = [int(np.where(self.adt_names == n)[0][0]) for n in names]
        return self._adt_matrix(rspec["statistic"])[:, idx].mean(1)

    # ---- tokens (shared with export)
    def tokens(self, X):
        out_tok, out_len = [], []
        for s in range(0, X.shape[0], 256):
            blk = X[s:s + 256][:, self.mapped]
            Xm = np.ascontiguousarray(blk.toarray() if hasattr(blk, "toarray") else blk, dtype=np.float32)
            tok, att = self.te.rank_encode_official(self.te.official_values(Xm, self.fac), self.token_ids, pad_id=PAD)
            out_tok.append(tok); out_len.append(att.sum(1))
        return np.concatenate(out_tok), np.concatenate(out_len)

    def export_tokens(self, decl, path):
        """Write official tokens of the selected cells for layer_states.py (rows = positions in this adapter)."""
        self.cells(decl)
        tok, L = self.tokens(self.Xs)
        np.savez_compressed(path, tok=tok.astype(np.int32), L=L.astype(np.int32), rows=self.sel, ids=self.obs_names[self.sel])
        self.log(f"exported tokens of {len(L)} cells (median length {int(np.median(L))}) to {path}")
        return path

    # ---- states
    def _raw(self, node, cells, arg):
        X = self.Xs
        tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
        keep = np.asarray((X[np.where(self._ctrl)[0]] > 0).mean(0)).ravel() >= 0.01
        D = X[:, np.where(keep)[0]].toarray()
        return np.log1p(D / tot[:, None] * 1e4).astype(np.float32)

    def _token_matrix(self):
        if "T" not in self._cache:
            import scipy.sparse as sp
            tok, L = self.tokens(self.Xs)
            m = np.full(int(self.token_ids.max()) + 1, -1, np.int64)
            m[self.token_ids] = np.arange(len(self.token_ids))
            r = np.repeat(np.arange(len(L)), L)
            pos = np.concatenate([np.arange(1, l + 1) for l in L]) if len(L) else np.zeros(0, int)
            c = m[np.concatenate([tok[i, :L[i]] for i in range(len(L))])] if len(L) else np.zeros(0, int)
            # exact view of the token sequence: gene at rank k gets (2049 - k) / 2048, absent genes 0
            T = sp.csr_matrix(((2049 - pos) / 2048, (r, c)), shape=(len(L), len(self.mapped)), dtype=np.float32)
            keep = np.asarray((T[np.where(self._ctrl)[0]] > 0).mean(0)).ravel() >= 0.01
            self._cache["T"] = T[:, np.where(keep)[0]].toarray()
        return self._cache["T"]

    def _tokens(self, node, cells, arg):
        return self._token_matrix()

    def _set(self, node, cells, arg):
        return (self._token_matrix() > 0).astype(np.float32)

    def _teddy(self):
        if "teddy" not in self._cache:
            if not self.paths.get("teddy_states"):
                raise RuntimeError("TEDDY states needed: export_tokens(), then run layer_states.py on the Mac and set "
                                   "paths.teddy_states")
            e = np.load(self.paths["teddy_states"])
            if "ids" in e.files:                         # match cells by name when the states carry them
                pos = {str(c): i for i, c in enumerate(e["ids"])}
                ii = np.array([pos[c] for c in self.obs_names[self.sel]])
            else:
                pos = {int(r): i for i, r in enumerate(e["rows"])}
                ii = np.array([pos[int(r)] for r in self.sel])
            self._cache["teddy"] = {k: e[k] for k in e.files}
            self._cache["teddy_idx"] = ii
        return self._cache["teddy"], self._cache["teddy_idx"]

    def _z(self, node, cells, arg):
        e, ii = self._teddy()
        return e["z_last"][ii].astype(np.float32)

    def _layer_mean(self, node, cells, arg):
        e, ii = self._teddy()
        return e["layer_means"][int(arg) - 1][ii].astype(np.float32)

    def _window_mean(self, node, cells, arg):
        e, ii = self._teddy()
        return np.concatenate(list(e["window_means"][int(arg)][:, ii]), axis=1).astype(np.float32)

    def _head(self):
        if "head" not in self._cache:
            import torch
            ck = torch.load(self.paths["head_ckpt"], map_location="cpu", weights_only=False)
            self._cache["head"] = {k: {n: v.numpy() for n, v in ck[k].items()} for k in ("mlp", "dec")}
        return self._cache["head"]

    @staticmethod
    def _mlp(W, x, prefix):
        from scipy.special import expit
        silu = lambda a: a * expit(a)                     # stable for large |a|
        h = silu(x @ W[f"{prefix}net.0.weight"].T + W[f"{prefix}net.0.bias"])
        h = silu(h @ W[f"{prefix}net.3.weight"].T + W[f"{prefix}net.3.bias"])
        return h @ W[f"{prefix}net.6.weight"].T + W[f"{prefix}net.6.bias"]

    @staticmethod
    def head_input(z):
        """The teddy_mm head reads unit-length embeddings: teddy_mm scripts/04_train.py (line 101) and 05_eval.py
        (line 55) divide each z by its L2 norm (+1e-6) before the MLP. Raw gene-mean embeddings have norm ~26."""
        z = np.asarray(z, dtype=np.float64)
        return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)

    def _head_latent(self, node, cells, arg):
        e, ii = self._teddy()
        return self._mlp(self._head()["mlp"], self.head_input(e["z_last"][ii]), "").astype(np.float32)

    def _head_output(self, node, cells, arg):
        lat = self._head_latent(node, cells, arg).astype(np.float64)
        mu = np.logaddexp(0.0, self._mlp(self._head()["dec"], lat, "mu."))     # softplus: NB mean per unit size factor
        return np.log(np.maximum(mu, 1e-6)).astype(np.float32)

    def _rna_of_gene(self, node, cells, arg):
        X = self.Xs
        tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
        genes = arg.split(",")
        idx = [int(np.where(self.sym == g)[0][0]) for g in genes if (self.sym == g).any()]
        if not idx:
            raise KeyError(f"genes {genes} not found in the RNA data")
        return np.log1p(X[:, idx].toarray() / tot[:, None] * 1e4).astype(np.float32)

    def _adt_readout(self, node, cells, arg):
        r = self._rspec_by_id[arg]
        names = r["features"] if isinstance(r["features"], list) else [r["features"]]
        idx = [int(np.where(self.adt_names == n)[0][0]) for n in names]
        return self._adt_matrix(r["statistic"])[:, idx].astype(np.float32)

    ops = {"log1p_cp10k": _raw, "teddy_official_tokens": _tokens, "drop_order": _set, "teddy_embedding": _z,
           "teddy_layer_mean": _layer_mean, "teddy_window_mean": _window_mean, "head_latent": _head_latent,
           "head_output": _head_output, "rna_of_gene": _rna_of_gene, "adt_readout": _adt_readout}

    def provenance(self):
        out = {"adapter": "perturb_cite"}
        for k in ("rna_h5ad", "adt_h5ad", "vocab", "medians", "head_ckpt", "teddy_states"):
            if self.paths.get(k):
                out[f"sha256_{k}"] = file_sha256(self.paths[k])[:16]
        return out
