#!/usr/bin/env python3
"""The input view of every group: 512 principal components of the RNA (log1p of counts per 10,000, genes detected in at
least 1% of the group's control cells: the RNA view of the earlier audits). No fit other than the components.

  python build_rna.py --work WORK --coculture-counts coculture_rna_counts.npz

For the two groups audited before, the RNA view is the exported node__RAW.npz under ANM_FEATURES. For the co-culture cells
it is computed here from their counts (rows, cols, vals, ids, genes, shape), taken from the screen's RNA file for the cells
that have embeddings. One SPEC per group and readout, with the single stage RNA512, for `anm_diag scope`: what a readout of
the other genes (`units`) and one with the gene's own cells (`own`) read from RNA.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml
from sklearn.decomposition import TruncatedSVD

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(HERE))
from _common import FEATURES, TARGETS  # noqa: E402
from build import ALPHAS, DECISION, GROUPS  # noqa: E402

SEED = 20261120


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--coculture-counts", type=Path, required=True)
    a = ap.parse_args()
    for name, g in GROUPS.items():
        K = dict(np.load(a.work / f"{name}.npz"))
        if name == "mel_co":
            c = np.load(a.coculture_counts)
            assert np.array_equal(c["ids"].astype(str), K["ids"])
            X = sp.csr_matrix((c["vals"].astype(np.float32), (c["rows"].astype(np.int64), c["cols"].astype(np.int64))), shape=tuple(int(v) for v in c["shape"]))
            tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1.0)
            keep = np.asarray((X[np.where(K["A"] == 0)[0]] > 0).mean(0)).ravel() >= 0.01
            R = X[:, np.where(keep)[0]].tocsr()
            R.data = np.log1p(R.data / np.repeat(tot, np.diff(R.indptr)) * 1e4).astype(np.float32)
        else:
            R = sp.load_npz(FEATURES / g["study"] / "node__RAW.npz").tocsr()
        assert R.shape[0] == len(K["A"])
        rna = TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED).fit_transform(R).astype(np.float32)
        np.savez(a.work / f"{name}_rna.npz", ids=K["ids"], RNA512=rna, A=K["A"], unit=K["unit"], **{k: v for k, v in K.items() if k.startswith("y_")})
        print(name, "RNA view", R.shape, "-> 512 components")
        for r in g["readouts"]:
            s = dict(name=f"{name}_{r}_rna", kind="perturbation", study=g["study"], readout_name=r, events=list(TARGETS[g["study"]]), also_pooled=True,
                     units_exclude=["control"], reads="RNA512", control_folds=5, decision=DECISION, label_basis="design", keep_draws=True,
                     readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": ALPHAS}},
                     keys=dict(readouts=[f"y_{r}"], stages={"RNA512": "RNA512"}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit"),
                     arrays={"audit": f"{name}_rna.npz"})
            (a.work / f"{name}_{r}_rna.yaml").write_text(yaml.safe_dump(s, sort_keys=False))


if __name__ == "__main__":
    main()
