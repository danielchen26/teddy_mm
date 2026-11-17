#!/usr/bin/env python3
"""Gate 0b for the final study F: the tables of a candidate study, the criteria, the choice of donors and cells, the pack.

  python final_study.py inspect H5AD                      # cell and feature tables only: columns, labels, names
  python final_study.py select  H5AD --train-pack DIR --pbmc-pack DIR --out final_selection.json
  python final_study.py pack    H5AD --selection final_selection.json --train-pack DIR --reference FEATURES.tsv.gz --out DIR

`inspect` and `select` read the cell table (obs) and the feature table (var) of the file and nothing else: no count and
no protein value. `pack` reads the counts of the chosen cells and writes a pack with the keys of data/processed/cite;
it prints sizes and names only.

The rule of DESIGN.md, gate 0b, as it is applied here (written before the file was opened):

  criteria   open; raw RNA counts per cell; raw protein counts; a donor and a disease label per cell; at least 20
             donors with at least 8 healthy and 8 patients; at least 80 of the 112 proteins the PBMC study shares with
             the head, by name; never used in this repository.
  groups     healthy donors, and patients: the donors with the disease the study is about. Other groups are not drawn.
  a donor's sample   the first by the study's own sample name; a donor is eligible if that sample has at least 1,000 cells
             and the donor has one status and one centre over all its cells.
  donors     20 per group (all, where fewer), balanced over centres: the centres in the order of their names, taken
             in turn; inside a centre the eligible donors in one random order (seed 0, the centres numbered in the
             order of their names; healthy and patients drawn apart).
  cells      1,000 per donor from that sample, uniformly at random without labels (seed 0, the donors numbered in the
             order of their names).
  proteins   a protein of the head is measured if the study has a feature of that name once the study's own prefix
             is removed (exact name), or under one of the names of ALIASES below, which were written from the feature
             names only.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

N_PER_GROUP, N_CELLS, MIN_PROTEINS, MIN_GROUP, SEED = 20, 1000, 80, 8, 0


# ------------------------------------------------------------------ the tables of an h5ad, with h5py only
def _text(a):
    a = np.asarray(a)
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a.ravel()]).reshape(a.shape)


def table(f, name):
    """obs or var as {column: array}, with the index under "_index". Both layouts of categoricals are read."""
    g = f[name]
    out = {}
    idx = g.attrs.get("_index", "index")
    idx = idx.decode() if isinstance(idx, bytes) else idx
    cats = g["__categories"] if "__categories" in g else None
    for k in g:
        if k == "__categories":
            continue
        v = g[k]
        if hasattr(v, "keys"):                                   # a categorical stored as a group
            if "categories" in v and "codes" in v:
                c, codes = _text(v["categories"][()]), np.asarray(v["codes"][()])
                out[k] = np.where(codes >= 0, c[np.clip(codes, 0, len(c) - 1)], "nan")
            continue
        a = v[()]
        if cats is not None and k in cats:
            c = _text(cats[k][()])
            out[k] = np.where(np.asarray(a) >= 0, c[np.clip(np.asarray(a), 0, len(c) - 1)], "nan")
        elif a.dtype.kind in "SO":
            out[k] = _text(a)
        else:
            out[k] = np.asarray(a)
    out["_index"] = out.pop(idx)
    return out


def matrix_info(g):
    if hasattr(g, "keys"):
        enc = g.attrs.get("encoding-type", g.attrs.get("h5sparse_format", "?"))
        enc = enc.decode() if isinstance(enc, bytes) else enc
        shape = g.attrs.get("shape", g.attrs.get("h5sparse_shape"))
        return {"encoding": str(enc), "shape": [int(x) for x in shape], "dtype": str(g["data"].dtype), "stored": int(g["data"].shape[0])}
    return {"encoding": "dense", "shape": list(g.shape), "dtype": str(g.dtype)}


def inspect(path):
    import h5py
    with h5py.File(path, "r") as f:
        obs, var = table(f, "obs"), table(f, "var")
        out = {"top": sorted(f.keys()), "layers": sorted(f["layers"].keys()) if "layers" in f else [], "X": matrix_info(f["X"]),
               "layers_info": {k: matrix_info(f["layers"][k]) for k in (f["layers"].keys() if "layers" in f else [])},
               "raw": ({"X": matrix_info(f["raw"]["X"]), "var_columns": sorted(k for k in f["raw"]["var"].keys())} if "raw" in f else None),
               "n_obs": int(len(obs["_index"])), "n_var": int(len(var["_index"])), "obs_columns": {}, "var_columns": {}}
    for k, v in obs.items():
        u = np.unique(v.astype(str))
        out["obs_columns"][k] = {"dtype": str(v.dtype), "n_unique": int(len(u)),
                                 "values": dict(Counter(v.astype(str).tolist()).most_common(40)) if len(u) <= 40 else [str(x) for x in u[:8]]}
    for k, v in var.items():
        u = np.unique(v.astype(str))
        out["var_columns"][k] = {"dtype": str(v.dtype), "n_unique": int(len(u)),
                                 "values": dict(Counter(v.astype(str).tolist()).most_common(12)) if len(u) <= 12 else [str(x) for x in u[:8]]}
    return out, obs, var


# ------------------------------------------------------------------ the studies this script knows, from their tables only
STUDIES = {
    "stephenson2021": dict(
        name="Stephenson et al. 2021, Nat Med 27:904, COVID-19 PBMC CITE-seq (ArrayExpress E-MTAB-10026; haniffa21.processed.h5ad)",
        donor="patient_id", sample="sample_id", status="Status", healthy="Healthy", patient="Covid", centre="Site",
        feature_type="feature_types", protein_type="Antibody Capture", gene_type="Gene Expression", prefix="AB_", counts="layers/raw",
        cell_type="initial_clustering"),
}
# Names of the head's proteins that the study lists under another name. Written on 17 November 2026 from the 192 feature
# names of the study and nothing else. The study names many antibodies by the gene symbol of their target; the head's
# panel uses CD numbers. An entry is made only where the two names are the same target in the standard nomenclature.
ALIASES: dict[str, dict[str, list[str]]] = {"stephenson2021": {
    # among the 112 proteins of the test
    "CD155": ["PVR"], "CD112": ["NECTIN2"], "CD154": ["CD40LG"], "CD11c": ["ITGAX"], "CD194": ["CCR4"], "CD279": ["PD1"], "CD335": ["NCR1"],
    "CD31": ["PECAM1"], "Podoplanin": ["podoplanin"], "CD146": ["MCAM"], "CD195": ["CCR5"], "CD196": ["CCR6"], "CD185": ["CXCR5"], "CD103": ["ITGAE"],
    "CD152": ["CTLA4"], "CD223": ["LAG3"], "CD107a": ["LAMP1"], "CD95": ["FAS"], "CD1c": ["CD1C"], "CD11b": ["ITGAM"], "CD314": ["KLRK1"], "CD35": ["CR1"],
    "CD57": ["B3GAT1"], "CD272": ["BTLA"], "CD278": ["ICOS"], "CD39": ["ENTPD1"], "CD169": ["SIGLEC1"], "integrinB7": ["ITGB7"], "CD268": ["BAFFR"],
    "CD54": ["ICAM1"], "CD62P": ["SELP"], "CD122": ["IL2RB"], "CD41": ["ITGA2B"], "CD137": ["TNFRSF9"], "CD124": ["IL4R"], "CD29": ["ITGB1"],
    "CD49b": ["ITGA2"], "CD18": ["ITGB2"], "CD127": ["IL7R"], "CD26": ["DPP4"], "CD49d": ["ITGA4"], "CD73": ["NT5E"], "TCRVa7.2": ["TCR_Va7.2"],
    "LOX-1": ["OLR1"], "CD158e1": ["KIR3DL1"], "CD319": ["SLAMF7"],
    # of the head's 134, not among the 112: written to the pack, not part of the test
    "HLA-A-B-C": ["HLA-ABC"], "CD49f": ["ITGA6"], "CD11a": ["ITGAL"], "FceRIa": ["FcERIa"], "CD94": ["KLRD1"], "CD85j": ["LILRB1"], "CD23": ["FCER2"],
    "CD328": ["SIGLEC7"], "CD88": ["C5AR1"], "CD224": ["GGT1"], "CD32": ["FCGR2A"],
}}
NOT_MATCHED_NOTE: dict[str, dict[str, str]] = {"stephenson2021": {
    "TCRVd2": "the study has TCR_Vg2 and TCR_Vg9; from the name alone TCR_Vg2 is not TCR V-delta-2, so it is not used",
    "CD134": "the study has OX40L, the ligand, not OX40",
    "CD270": "the study has LIGHT, a ligand of CD270, not CD270",
    "CD48": "", "CD105": "", "CD42b": "", "CD119": "", "CD192": "the study has CCR3 to CCR7, not CCR2", "CD13": "", "CD115": "", "CD63": "",
    "CD172a": "", "CD72": "", "CD93": "", "CD49a": "", "CD9": "", "CD142": "",
}}


def protein_mapping(ours, source, study):
    """ours: the head's protein names. source: the study's protein feature names, prefix removed. One row per protein
    of the head: the study's features it is read from, and by which rule."""
    src, rows = set(source), []
    for p in ours:
        if p in ALIASES[study]:
            feats, rule = ALIASES[study][p], "alias (same target under another name)" if len(ALIASES[study][p]) == 1 else "mean of the raw counts of several features of the same target"
        elif p in src:
            feats, rule = [p], "exact name"
        else:
            feats, rule = [], "not measured" + (": " + NOT_MATCHED_NOTE[study][p] if NOT_MATCHED_NOTE[study].get(p) else "")
        for f in feats:
            if f not in src:
                raise SystemExit(f"mapping for {p} names {f}, which is not a feature of the study")
        rows.append({"ours": p, "source": feats, "rule": rule})
    dup = [f for f, n in Counter(f for r in rows for f in r["source"]).items() if n > 1]
    if dup:
        raise SystemExit(f"features of the study mapped twice: {dup}")
    return rows


def choose(obs, cfg, seed=SEED):
    """Donors, their sample and their cells, by the rule of the module's docstring. obs: {column: array}. Returns
    (per-donor table, chosen donors by group, {donor: rows of obs})."""
    donor, sample, status, centre = (np.asarray(obs[cfg[k]]).astype(str) for k in ("donor", "sample", "status", "centre"))
    rows = []
    for d in sorted(set(donor.tolist())):
        me = donor == d
        st, ce, sm = sorted(set(status[me].tolist())), sorted(set(centre[me].tolist())), sorted(set(sample[me].tolist()))
        n_first = int((me & (sample == sm[0])).sum())
        group = "healthy" if st == [cfg["healthy"]] else ("patient" if st == [cfg["patient"]] else None)
        rows.append({"donor": d, "status_in_the_study": st, "centres": ce, "samples": len(sm), "first_sample": sm[0], "cells_in_first_sample": n_first,
                     "cells": int(me.sum()), "group": group, "eligible": bool(group is not None and len(ce) == 1 and n_first >= N_CELLS)})
    chosen = {}
    for gi, group in enumerate(("healthy", "patient")):
        el = [r for r in rows if r["eligible"] and r["group"] == group]
        queues = []
        for k, c in enumerate(sorted({r["centres"][0] for r in el})):
            ds = sorted(r["donor"] for r in el if r["centres"][0] == c)
            queues.append([ds[i] for i in np.random.default_rng([seed, gi, k]).permutation(len(ds))])
        pick = []
        while len(pick) < N_PER_GROUP and any(queues):
            for q in queues:
                if q and len(pick) < N_PER_GROUP:
                    pick.append(q.pop(0))
        chosen[group] = sorted(pick)
    by = {r["donor"]: r for r in rows}
    cells = {}
    for k, d in enumerate(sorted(chosen["healthy"] + chosen["patient"])):
        idx = np.flatnonzero((donor == d) & (sample == by[d]["first_sample"]))
        cells[d] = np.sort(np.random.default_rng([seed, 9, k]).choice(idx, size=N_CELLS, replace=False))
    return rows, chosen, cells


def rows_digest(cells):
    import hashlib
    h = hashlib.sha256()
    for d in sorted(cells):
        h.update(d.encode())
        h.update(np.asarray(cells[d], np.int64).tobytes())
    return h.hexdigest()


def select(path, study, train_pack, pbmc_pack):
    import h5py
    cfg = STUDIES[study]
    with h5py.File(path, "r") as f:
        obs, var = table(f, "obs"), table(f, "var")
        layout = {"X": matrix_info(f["X"]), "counts": matrix_info(f[cfg["counts"]]) if cfg["counts"] in f else None}
    ft = np.asarray(var[cfg["feature_type"]]).astype(str)
    names = np.asarray(var["_index"]).astype(str)
    prot = [n[len(cfg["prefix"]):] if n.startswith(cfg["prefix"]) else n for n in names[ft == cfg["protein_type"]]]
    with np.load(Path(train_pack) / "cite_arrays.npz", allow_pickle=False) as z:
        ours = [str(x) for x in z["adt_names"]]
    with np.load(Path(pbmc_pack) / "cite_arrays.npz", allow_pickle=False) as z:
        in_pbmc = np.asarray(z["adt_present"]).astype(bool)
    mapping = protein_mapping(ours, prot, study)
    measured = [bool(r["source"]) for r in mapping]
    among = [r["ours"] for r, m, q in zip(mapping, measured, in_pbmc) if m and q]
    rows, chosen, cells = choose(obs, cfg)
    groups = Counter(r["group"] for r in rows)
    crit = {
        "open": {"pass": True, "how": "public file of the study's data portal; no access request"},
        "raw_rna_and_protein_counts_per_cell": {"pass": bool(layout["counts"]), "how": f"{cfg['counts']} of the file: {layout['counts']}; that it holds integers is checked when the pack is written"},
        "donor_and_disease_label_per_cell": {"pass": True, "how": f"obs columns {cfg['donor']} and {cfg['status']}"},
        "at_least_20_donors_8_healthy_8_patients": {"pass": bool(groups["healthy"] >= MIN_GROUP and groups["patient"] >= MIN_GROUP and groups["healthy"] + groups["patient"] >= 20),
                                                     "healthy": groups["healthy"], "patients": groups["patient"], "other_or_mixed": groups[None]},
        "at_least_80_of_the_112_proteins": {"pass": bool(len(among) >= MIN_PROTEINS), "n": len(among)},
        "never_used_in_this_repository": {"pass": True, "how": "no file, script or result of the repository names it before 17 November 2026"},
    }
    by = {r["donor"]: r for r in rows}
    out = {
        "study": study, "name": cfg["name"], "file": Path(path).name, "n_cells_in_file": int(len(names) and len(obs["_index"])), "n_features_in_file": int(len(names)),
        "layout": layout, "criteria": crit, "meets_the_criteria": bool(all(c["pass"] for c in crit.values())),
        "rule": {"per_group": N_PER_GROUP, "cells_per_donor": N_CELLS, "seed": SEED},
        "eligible": {g: sum(r["eligible"] and r["group"] == g for r in rows) for g in ("healthy", "patient")},
        "centres": {g: dict(Counter(by[d]["centres"][0] for d in chosen[g])) for g in ("healthy", "patient")},
        "donors": [{"donor": f"{study}_{d}", "status": g, "centre": by[d]["centres"][0], "sample": by[d]["first_sample"], "cells_in_sample": by[d]["cells_in_first_sample"],
                    "samples_of_the_donor": by[d]["samples"]} for g in ("healthy", "patient") for d in chosen[g]],
        "reserve": sorted(f"{study}_{r['donor']}" for r in rows if r["group"] is not None and r["donor"] not in chosen["healthy"] + chosen["patient"]),
        "cells": {"n": int(sum(len(v) for v in cells.values())), "sha256_of_rows": rows_digest(cells)},
        "proteins": {"features_of_the_study": len(prot), "mapping": mapping, "measured_of_134": int(sum(measured)), "among_the_112": among,
                     "of_the_112_not_measured": [r["ours"] for r, m, q in zip(mapping, measured, in_pbmc) if q and not m]},
    }
    return out, cells, obs, var


# ------------------------------------------------------------------ the pack
def make_unique(names):
    """Names made unique as anndata does for var_names: a repeated name keeps its first use and gets -1, -2 after."""
    dup, seen, out = {n for n, c in Counter(names).items() if c > 1}, Counter(), []
    for n in names:
        out.append(n if n not in dup or seen[n] == 0 else f"{n}-{seen[n]}")
        seen[n] += 1
    return out


def gene_columns(var, cfg, reference, our_genes):
    """For every feature of the study, its column among the head's genes (-1: none). The study names genes by symbol;
    the symbols are those of the 10x GRCh38-3.0.0 reference, whose table gives the Ensembl id of each."""
    import gzip
    names, ft = np.asarray(var["_index"]).astype(str), np.asarray(var[cfg["feature_type"]]).astype(str)
    with gzip.open(reference, "rt") as f:
        ref = [l.rstrip("\n").split("\t") for l in f]
    ref = [r for r in ref if len(r) < 3 or r[2] == "Gene Expression"]
    ens_of = dict(zip(make_unique([r[1] for r in ref]), [r[0] for r in ref]))
    pos = {g: j for j, g in enumerate(our_genes)}
    col = np.full(len(names), -1, np.int64)
    known = 0
    for i, (n, t) in enumerate(zip(names, ft)):
        if t != cfg["gene_type"]:
            continue
        e = ens_of.get(n)
        known += e is not None
        if e in pos:
            col[i] = pos[e]
    used = col[col >= 0]
    if len(set(used.tolist())) != len(used):
        raise SystemExit("two features of the study map to one gene of the head")
    return col, {"genes_in_the_study": int((ft == cfg["gene_type"]).sum()), "with_an_ensembl_id_in_the_reference": int(known),
                 "genes_of_the_head": int(len(our_genes)), "of_them_in_the_study": int(len(used))}


def read_rows(g, ip, rows, n_features):
    """Rows of a CSR matrix stored in an h5ad group, as (data, indices, indptr). ip: the matrix's row pointer;
    rows: sorted row numbers. One slab from the first to the last of the rows is read."""
    lo, hi = int(rows[0]), int(rows[-1]) + 1
    ptr = ip[lo:hi + 1]
    data, idx = np.asarray(g["data"][ptr[0]:ptr[-1]]), np.asarray(g["indices"][ptr[0]:ptr[-1]], np.int64)
    ptr = ptr - ptr[0]
    keep = np.concatenate([np.arange(ptr[r - lo], ptr[r - lo + 1]) for r in rows])
    counts = np.array([ptr[r - lo + 1] - ptr[r - lo] for r in rows], np.int64)
    if idx.size and idx.max() >= n_features:
        raise SystemExit("a column index of the counts is beyond the feature table")
    return data[keep], idx[keep], np.concatenate([[0], np.cumsum(counts)])


def clr(src_counts, counts, present):
    """Per-cell centred log ratio over the whole measured panel, as scripts/prepare_external_cite.py: g = exp(mean over all
    the study's protein features of log1p(count)); value = log1p(count / g)."""
    g = np.exp(np.log1p(src_counts.astype(np.float64)).mean(axis=1))
    g[g <= 0] = 1.0
    out = np.full(counts.shape, np.nan, np.float32)
    out[:, present] = np.log1p(counts[:, present].astype(np.float64) / g[:, None]).astype(np.float32)
    return out, g.astype(np.float32)


def size_factors(tot):
    """teddy_mm.data.size_factors on per-cell totals: the total over its median among cells with a total above 0."""
    tot = np.asarray(tot, np.float64)
    med = np.median(tot[tot > 0]) if np.any(tot > 0) else 1.0
    sf = tot / max(med, 1e-6)
    sf[sf <= 0] = 1.0
    return sf.astype(np.float32)


def pack_blocks(path, study, rows, col, gene_features, n_genes_out, prot_idx, cache, block=20000, seconds=0):
    """The counts of the chosen cells, read in blocks of rows of the file and written to the cache, one file a block.
    rows: all chosen rows, sorted. Returns the number of blocks still to read (0: all there)."""
    import time
    import h5py
    from scipy import sparse
    cfg, t0 = STUDIES[study], time.time()
    starts = sorted({int(r) // block * block for r in rows})
    left = [b for b in starts if not (Path(cache) / f"block_{b:07d}.npz").exists()]
    if not left:
        return 0
    with h5py.File(path, "r") as f:
        g = f[cfg["counts"]]
        enc = g.attrs.get("encoding-type", g.attrs.get("h5sparse_format", ""))
        enc = enc.decode() if isinstance(enc, bytes) else enc
        if "csr" not in str(enc):
            raise SystemExit(f"the counts are stored as {enc!r}; this script reads rows of a CSR matrix")
        ip = np.asarray(g["indptr"][()], np.int64)
        for k, b in enumerate(left):
            if seconds and time.time() - t0 > seconds:
                return len(left) - k
            r = rows[(rows >= b) & (rows < b + block)]
            data, idx, ptr = read_rows(g, ip, r, len(col))
            if not np.all(data == np.round(data)) or (data.size and data.min() < 0):
                raise SystemExit(f"rows {b} to {b + block}: the counts are not non-negative integers")
            X = sparse.csr_matrix((data.astype(np.float32), idx, ptr), shape=(len(r), len(col)))
            genes = np.flatnonzero(col >= 0)
            Xg = X[:, genes].tocoo()
            R = sparse.csr_matrix((Xg.data, (Xg.row, col[genes][Xg.col])), shape=(len(r), n_genes_out), dtype=np.float32)
            R.sort_indices()
            out = Path(cache) / f"block_{b:07d}.npz"
            tmp = out.with_suffix(".tmp.npz")
            np.savez(tmp, rna_data=R.data, rna_indices=R.indices.astype(np.int32), rna_indptr=R.indptr.astype(np.int64),
                     adt_source_counts=np.asarray(X[:, prot_idx].todense(), np.float32),
                     n_count_rna_all_genes=np.asarray(X[:, gene_features].sum(axis=1)).ravel().astype(np.float64), rows=r.astype(np.int64))
            tmp.replace(out)
            print(f"rows {b} to {b + block}: {len(r)} cells ({time.time() - t0:.0f} s)", flush=True)
    return 0


def assemble(study, selection, rows_by_donor, obs, var, train_pack, cache, out_dir, gene_info, compress=True):
    from scipy import sparse
    cfg = STUDIES[study]
    with np.load(Path(train_pack) / "cite_arrays.npz", allow_pickle=False) as z:
        ours = [str(x) for x in z["adt_names"]]
        our_genes = np.array([str(x) for x in z["rna_names"]])
        tokens = np.asarray(z["token_ids"], np.int64)
    names, ft = np.asarray(var["_index"]).astype(str), np.asarray(var[cfg["feature_type"]]).astype(str)
    src_names = [n[len(cfg["prefix"]):] if n.startswith(cfg["prefix"]) else n for n in names[ft == cfg["protein_type"]]]
    mapping = selection["proteins"]["mapping"]
    if [r["ours"] for r in mapping] != ours:
        raise SystemExit("the mapping of the selection is not in the head's protein order")
    pos = {n: j for j, n in enumerate(src_names)}
    files = sorted(Path(cache).glob("block_*.npz"))
    blocks, A, tot_all, row_all = [], [], [], []
    for p in files:
        z = np.load(p)
        blocks.append(sparse.csr_matrix((z["rna_data"], z["rna_indices"], z["rna_indptr"]), shape=(len(z["rows"]), len(our_genes))))
        A.append(z["adt_source_counts"]); tot_all.append(z["n_count_rna_all_genes"]); row_all.append(z["rows"])
    have = np.concatenate(row_all)
    donors, status, centre, sample, want = [], [], [], [], []
    for d in selection["donors"]:
        r = rows_by_donor[d["donor"][len(study) + 1:]]
        want.append(r)
        donors += [d["donor"]] * len(r); status += [d["status"]] * len(r); centre += [f"{study}_{d['centre']}"] * len(r); sample += [d["sample"]] * len(r)
    rows = np.concatenate(want)
    if not (np.all(np.diff(have) > 0) and np.array_equal(np.sort(rows), have)):
        raise SystemExit("the cached cells are not the chosen cells")
    order = np.searchsorted(have, rows)
    X = sparse.vstack(blocks).tocsr()[order]
    A, tot_all = np.vstack(A)[order], [np.concatenate(tot_all)[order]]
    counts = np.full((len(rows), len(ours)), np.nan, np.float32)
    present = np.zeros(len(ours), bool)
    for k, r in enumerate(mapping):
        if r["source"]:
            counts[:, k] = A[:, [pos[f] for f in r["source"]]].astype(np.float64).mean(axis=1)
            present[k] = True
    adt, g = clr(A, counts, present)
    lab = np.asarray(obs[cfg["cell_type"]]).astype(str)[rows] if cfg.get("cell_type") in obs else np.full(len(rows), "", "U1")
    arrays = dict(
        rna_data=X.data.astype(np.float32), rna_indices=X.indices.astype(np.int32), rna_indptr=X.indptr, rna_shape=np.array(X.shape),
        adt=adt, adt_counts=counts, adt_present=present, adt_clr_scale=g, adt_source_counts=A, adt_source_names=np.array(src_names, dtype="U32"),
        token_ids=tokens, split=np.full(len(rows), "test", dtype="U8"), donors=np.array(donors, dtype="U48"), sites=np.array(centre, dtype="U48"),
        status=np.array(status, dtype="U8"), samples=np.array(sample, dtype="U48"), cell_types=lab.astype("U64"),
        barcodes=np.asarray(obs["_index"]).astype(str)[rows].astype("U64"), source_row=rows.astype(np.int64),
        n_count_rna=np.concatenate(tot_all).astype(np.int64), rna_size_factor=size_factors(np.asarray(X.sum(axis=1)).ravel()),
        adt_size_factor=size_factors(np.nan_to_num(counts[:, present]).sum(axis=1)), adt_names=np.array(ours, dtype="U64"), rna_names=our_genes.astype("U32"))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out_dir / "cite_arrays.tmp.npz"
    (np.savez_compressed if compress else np.savez)(tmp, **arrays)
    tmp.replace(out_dir / "cite_arrays.npz")
    meta = {"dataset": cfg["name"], "external": True, "role": "final (round1_coverage_20261116, gate 5): read once", "n_cells": int(len(rows)),
            "n_rna_mapped": int(len(our_genes)), "n_adt": len(ours), "n_adt_present": int(present.sum()), "n_train": 0, "n_val": 0, "n_test": int(len(rows)),
            "test_sites": sorted(set(centre)), "split_definition": "every cell is evaluation-only ('test'); nothing may be fit on it",
            "gene_universe": "train", "cell_selection": f"{N_CELLS} cells per donor from the donor's first sample, without labels, seed {SEED}; "
            f"{len(selection['donors'])} donors chosen by scripts/final_study.py", "rna_values": "raw integer UMI counts (float32)",
            "adt_values": f"per-cell CLR log1p(c/g_i), g_i over all {len(src_names)} measured protein features; NaN = not measured", "genes": gene_info}
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    return {"cells": int(len(rows)), "genes": int(X.shape[1]), "stored_counts": int(X.nnz), "proteins_measured": int(present.sum()),
            "protein_features": int(len(src_names)), "donors": len(selection["donors"])}


def main(argv=None):
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("inspect")
    q.add_argument("h5ad", type=Path)
    q.add_argument("--out", type=Path)
    q = sub.add_parser("select")
    q.add_argument("h5ad", type=Path)
    q.add_argument("--study", default="stephenson2021", choices=sorted(STUDIES))
    q.add_argument("--train-pack", type=Path, required=True)
    q.add_argument("--pbmc-pack", type=Path, required=True)
    q.add_argument("--out", type=Path, required=True)
    q.add_argument("--rows-out", type=Path, help="the chosen rows of the cell table, by donor (npz), for the pack stage")
    q = sub.add_parser("pack")
    q.add_argument("h5ad", type=Path)
    q.add_argument("--study", default="stephenson2021", choices=sorted(STUDIES))
    q.add_argument("--selection", type=Path, required=True)
    q.add_argument("--rows", type=Path, required=True, help="the npz `select --rows-out` wrote")
    q.add_argument("--train-pack", type=Path, required=True)
    q.add_argument("--reference", type=Path, required=True, help="feature table of the 10x GRCh38-3.0.0 reference (tsv.gz: Ensembl id, symbol)")
    q.add_argument("--cache", type=Path, required=True)
    q.add_argument("--out", type=Path, required=True)
    q.add_argument("--seconds", type=float, default=0, help="stop reading after this long; run again to go on")
    q.add_argument("--block", type=int, default=20000, help="rows of the file read at a time")
    q.add_argument("--no-compress", action="store_true")
    a = p.parse_args(argv)
    if a.cmd == "pack":
        import h5py
        cfg = STUDIES[a.study]
        sel = json.loads(a.selection.read_text())
        rows = dict(np.load(a.rows))
        if rows_digest(rows) != sel["cells"]["sha256_of_rows"]:
            raise SystemExit("the rows file is not that of the selection")
        with h5py.File(a.h5ad, "r") as f:
            obs, var = table(f, "obs"), table(f, "var")
        with np.load(a.train_pack / "cite_arrays.npz", allow_pickle=False) as z:
            our_genes = [str(x) for x in z["rna_names"]]
        col, gene_info = gene_columns(var, cfg, a.reference, our_genes)
        ft = np.asarray(var[cfg["feature_type"]]).astype(str)
        prot_idx = np.flatnonzero(ft == cfg["protein_type"])
        a.cache.mkdir(parents=True, exist_ok=True)
        left = pack_blocks(a.h5ad, a.study, np.sort(np.concatenate(list(rows.values()))), col, np.flatnonzero(ft == cfg["gene_type"]), len(our_genes), prot_idx,
                           a.cache, block=a.block, seconds=a.seconds)
        if left:
            print(f"{left} blocks left; run again")
            return 75
        info = assemble(a.study, sel, rows, obs, var, a.train_pack, a.cache, a.out, gene_info, compress=not a.no_compress)
        print(json.dumps({**info, "genes_info": gene_info}, indent=1))
        return 0
    if a.cmd == "inspect":
        out, obs, var = inspect(a.h5ad)
        text = json.dumps(out, indent=1)
        if a.out:
            a.out.write_text(text)
        print(text)
    elif a.cmd == "select":
        out, cells, _, _ = select(a.h5ad, a.study, a.train_pack, a.pbmc_pack)
        a.out.write_text(json.dumps(out, indent=1) + "\n")
        if a.rows_out:
            np.savez(a.rows_out, **{d: v for d, v in cells.items()})
        print(json.dumps({k: out[k] for k in ("criteria", "meets_the_criteria", "eligible", "centres")}, indent=1))
        print("donors", len(out["donors"]), "cells", out["cells"]["n"], "proteins among the 112:", len(out["proteins"]["among_the_112"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
