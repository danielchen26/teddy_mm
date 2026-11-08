"""Symbol -> Ensembl map for an RNA h5ad whose genes are named by symbol (paths key symbol_map).

  python tools/gene_map.py RNA.h5ad FEATURES_10x.tsv.gz VOCAB.txt OUT.csv [--var-col ensembl_id]

Why: scPerturb's var column ensembl_id gives alternate-haplotype ids for several MHC genes (HLA-A, -B, -C, -E,
B2M) and for IFNGR2; those ids are not in TEDDY's vocabulary, so those genes would silently never reach TEDDY.
teddy_mm names its own CITE-seq genes with the 10x Cell Ranger GRCh38-3.0.0 feature table, so the same table is
used here. Rule, per gene symbol: the 10x id (the k-th repeat of a symbol is SYMBOL.k or SYMBOL-k); if that id is
not in the TEDDY vocabulary and the file's own id (--var-col) is, the file's own id; otherwise none. Reads only
the var table (h5py), never a count matrix. Prints the coverage and the ids of a few marker genes.
"""
import csv
import gzip
import sys
from collections import Counter

import h5py
import numpy as np


def _col(e):
    if isinstance(e, h5py.Group):                 # categorical
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in e["categories"][:]])
        c = e["codes"][:]
        return np.where(c >= 0, cats[np.maximum(c, 0)], "")
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in e[:]])


def ten_x(features):
    seen, out = Counter(), {}
    with gzip.open(features, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or (len(p) > 2 and p[2] != "Gene Expression"):
                continue
            k = seen[p[1]]
            seen[p[1]] += 1
            for name in ([p[1]] if k == 0 else [f"{p[1]}.{k}", f"{p[1]}-{k}"]):
                out.setdefault(name, p[0])
    return out


def main(argv):
    var_col = "ensembl_id"
    if "--var-col" in argv:
        i = argv.index("--var-col"); var_col = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    rna, features, vocab, out = argv
    voc = {l.rstrip("\n") for l in open(vocab)}
    m10 = ten_x(features)
    with h5py.File(rna, "r") as f:
        v = f["var"]
        idx = v.attrs["_index"]
        idx = idx.decode() if isinstance(idx, bytes) else idx
        sym = _col(v[idx])
        own = _col(v[var_col]) if var_col in v else np.array([""] * len(sym))
    rows, src = [], Counter()
    for s, o in zip(sym, own):
        t = m10.get(s, "")
        if t and (t in voc or o not in voc):
            e, k = t, "10x"
        elif o in voc:
            e, k = o, "file"
        else:
            e, k = "", "none"
        src[k] += 1
        rows.append((s, e))
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["symbol", "ensembl"]); w.writerows(rows)
    inv = sum(e in voc for _, e in rows)
    print(f"{len(sym)} genes: {src['10x']} from the 10x table, {src['file']} from the file's {var_col}, "
          f"{src['none']} unmapped; {inv} in the TEDDY vocabulary -> {out}")
    d = dict(rows)
    print("  " + ", ".join(f"{g} {d.get(g) or '-'}" for g in ("HLA-A", "B2M", "HLA-E", "CD274", "IFNGR2")))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
