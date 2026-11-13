"""Which knockouts change which measured proteins, by condition: a census of the sources an audit could use.

Reads only the protein file of a Perturb-CITE-seq or ECCITE-seq data set (scPerturb format). Per condition, target gene
with at least --min-cells single-guide cells, and antibody: the shift of the centred-log-ratio protein level against the
control-guide cells, in control standard deviations, and that shift over the half-width of its 95% interval (cells as
units, so the ratio is optimistic). Prints one JSON object.

  python data_census.py --protein-h5ad FILE --guide-col guide_id [--condition-col perturbation_2] [--head-names names.json]
"""
import argparse, json, re
import anndata as ad, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--protein-h5ad", required=True); ap.add_argument("--guide-col", required=True); ap.add_argument("--condition-col")
ap.add_argument("--guide-sep", default=";"); ap.add_argument("--guide-gene-regex", default=r"^(.+?)[_g]?\d+$"); ap.add_argument("--control-regex", default=r"SITE|NT|control|ctrl")
ap.add_argument("--head-names"); ap.add_argument("--head-map", default="CD274=CD274,HLA_A=HLA-A-B-C,HLA_E=HLA-E,PDL1=CD274")
ap.add_argument("--min-cells", type=int, default=100); ap.add_argument("--min-shift", type=float, default=0.3); ap.add_argument("--min-snr", type=float, default=3.0)
a = ap.parse_args()
A = ad.read_h5ad(a.protein_h5ad); X = A.X.toarray() if hasattr(A.X, "toarray") else np.asarray(A.X)
L = np.log1p(X.astype(float)); C = L - L.mean(1, keepdims=True); P = [str(p) for p in A.var_names]
g = A.obs[a.guide_col].astype(str).values
cond = A.obs[a.condition_col].astype(str).values if a.condition_col else np.full(len(g), "all")
def genes(s):
    if s in ("nan", "None", ""): return []
    return sorted({(re.match(a.guide_gene_regex, x).group(1) if re.match(a.guide_gene_regex, x) else x) for x in s.split(a.guide_sep)})
G = [genes(s) for s in g]; single = np.array([len(x) == 1 for x in G]); tgt = np.array([x[0] if len(x) == 1 else "" for x in G])
isctl = np.array([bool(re.search(a.control_regex, t)) for t in tgt])
hmap = dict(x.split("=") for x in a.head_map.split(",")); head = set(json.load(open(a.head_names))) if a.head_names else None
in_head = {p: (hmap.get(p, p) if (head is None or hmap.get(p, p) in head) else None) for p in P}
iso = [p for p in P if "IgG" in p]
out = {"cells": int(len(g)), "antibodies": P, "antibodies_with_a_head_output": {p: v for p, v in in_head.items() if v and head is not None}, "isotype_controls": iso,
       "thresholds": {"min_cells": a.min_cells, "min_shift_sd": a.min_shift, "min_shift_over_half_width": a.min_snr}, "conditions": {}}
for c in sorted(set(cond)):
    m = (cond == c) & single; ctl = m & isctl; n0 = int(ctl.sum())
    if n0 < a.min_cells: continue
    mu0, sd0, v0 = C[ctl].mean(0), C[ctl].std(0), C[ctl].var(0); pairs, targets = [], {}
    for t in sorted(set(tgt[m & ~isctl])):
        T = m & (tgt == t); n1 = int(T.sum()); targets[t] = n1
        if n1 < a.min_cells: continue
        d = (C[T].mean(0) - mu0) / sd0; se = np.sqrt(C[T].var(0) / n1 + v0 / n0) / sd0
        for j, p in enumerate(P):
            r = abs(d[j]) / (1.96 * se[j])
            if p not in iso and abs(d[j]) >= a.min_shift and r >= a.min_snr:
                pairs.append({"target": t, "antibody": p, "cells": n1, "shift_sd": round(float(d[j]), 3), "over_half_width": round(float(r), 1), "head_output": in_head[p] if head is not None else None})
    pairs.sort(key=lambda x: -abs(x["shift_sd"]))
    out["conditions"][c] = {"cells": int((cond == c).sum()), "single_guide_cells": int(m.sum()), "control_guide_cells": n0, "targets": len(targets),
                            "targets_with_min_cells": int(sum(v >= a.min_cells for v in targets.values())), "clear_responses": pairs,
                            "targets_with_a_clear_response": sorted({x["target"] for x in pairs})}
print(json.dumps(out, ensure_ascii=False))
