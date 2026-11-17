#!/usr/bin/env python3
"""Post hoc, after the test was read: is the pack of the final study sound?

  python check_final_pack.py --final-arrays F_cite_arrays.npz --pbmc-arrays ext_cite_arrays.npz

For proteins whose coding gene is plain, the correlation within a donor between the protein (as the pack holds it) and
the log-normalised RNA of its own gene, mean over donors. If the rows of the protein and RNA tables were not the same
cells, or a protein were matched to the wrong feature, these would be near 0. The same pairs are read on the
development study with the same code, as the reference for each pair.

No head enters this. It was written after results/final.json existed and changes nothing in it.
Writes ../results/final_pack_check.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

import common

# protein, as the head's panel names it -> Ensembl id of the gene that codes for it
PAIRS = {"CD3": "ENSG00000198851", "CD4": "ENSG00000010610", "CD8": "ENSG00000153563", "CD19": "ENSG00000177455", "CD14": "ENSG00000170458",
         "CD20": "ENSG00000156738", "CD56": "ENSG00000149294", "CD16": "ENSG00000203747", "CD127": "ENSG00000168685", "HLA-DR": "ENSG00000204287",
         "CD11c": "ENSG00000140678", "CD123": "ENSG00000185291", "CD27": "ENSG00000139193", "CD28": "ENSG00000178562", "CD38": "ENSG00000004468",
         "CD69": "ENSG00000110848", "CD62L": "ENSG00000188404", "CD11b": "ENSG00000169896", "CD1c": "ENSG00000158481", "CD25": "ENSG00000134460",
         "CD86": "ENSG00000114013", "CD33": "ENSG00000105383", "CD7": "ENSG00000173762", "CD52": "ENSG00000169442", "CD40": "ENSG00000101017",
         "CD79b": "ENSG00000007312", "CD161": "ENSG00000111796", "CD5": "ENSG00000110448", "CD2": "ENSG00000116824", "CD36": "ENSG00000135218"}


def pairs_of(path):
    a = np.load(path, allow_pickle=False)
    X = sp.csr_matrix((a["rna_data"].astype(np.float32), a["rna_indices"], a["rna_indptr"]), shape=tuple(int(v) for v in a["rna_shape"]))
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1.0
    names, genes = [str(x) for x in a["adt_names"]], [str(x) for x in a["rna_names"]]
    present, donors = a["adt_present"].astype(bool), a["donors"].astype(str)
    out = {}
    for prot, ens in PAIRS.items():
        if prot not in names or not present[names.index(prot)] or ens not in genes:
            continue
        y = a["adt"][:, names.index(prot)].astype(np.float64)
        g = np.log1p(np.asarray(X[:, genes.index(ens)].todense()).ravel() / tot * 1e4)
        r = [np.corrcoef(y[m], g[m])[0, 1] for m in (donors == d for d in np.unique(donors)) if y[m].std() > 1e-8 and g[m].std() > 1e-8]
        if r:
            out[prot] = float(np.mean(r))
    return {"donors": int(len(np.unique(donors))), "pairs": out}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--final-arrays", type=Path, required=True)
    p.add_argument("--pbmc-arrays", type=Path, required=True)
    a = p.parse_args()
    if not (common.RUN / "results" / "final.json").exists():
        raise SystemExit("the test has not been read; this check comes after it")
    out = {"what": "post hoc, after the test: within a donor, the correlation between a protein and the RNA of its own gene; mean over donors",
           "genes": PAIRS, "final_study": pairs_of(a.final_arrays), "development_study": pairs_of(a.pbmc_arrays)}
    both = [k for k in out["final_study"]["pairs"] if k in out["development_study"]["pairs"]]
    out["mean_over_pairs_in_both"] = {"pairs": len(both), "final_study": float(np.mean([out["final_study"]["pairs"][k] for k in both])),
                                      "development_study": float(np.mean([out["development_study"]["pairs"][k] for k in both]))}
    (common.RUN / "results" / "final_pack_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["mean_over_pairs_in_both"]), {k: round(v, 2) for k, v in out["final_study"]["pairs"].items()})


if __name__ == "__main__":
    main()
