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
        feature_type="feature_types", protein_type="Antibody Capture", gene_type="Gene Expression", prefix="AB_", counts="layers/raw"),
}
# names of the head's proteins that the study lists under another name; filled from the feature names only (see GATE5.md, part 2)
ALIASES: dict[str, dict[str, list[str]]] = {"stephenson2021": {}}
NOT_MATCHED_NOTE: dict[str, dict[str, str]] = {"stephenson2021": {}}


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
            feats, rule = [], "not measured" + (": " + NOT_MATCHED_NOTE[study][p] if p in NOT_MATCHED_NOTE[study] else "")
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
    a = p.parse_args(argv)
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
