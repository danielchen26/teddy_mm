"""Symbol -> Ensembl map (CSV with columns symbol,ensembl) from a 10x Cell Ranger reference feature table.

  python tools/symbol_map_10x.py FEATURES.tsv.gz OUT.csv

For an RNA file that names genes by symbol only (paths key symbol_map). teddy_mm maps its own CITE-seq
genes with the 10x GRCh38-3.0.0 feature table (scripts/prepare_external_cite.py), so the same table keeps the
two consistent. A symbol that the table repeats keeps its first Ensembl id under the bare name; the k-th
repeat is written as SYMBOL.k (Seurat make.unique) and SYMBOL-k (scanpy var_names_make_unique), the two
names a count matrix exported from either tool gives it. Standard library only.
"""
import csv
import gzip
import sys
from collections import Counter


def main(features, out):
    seen = Counter()
    rows = []
    with gzip.open(features, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or (len(p) > 2 and p[2] != "Gene Expression"):
                continue
            ens, sym = p[0], p[1]
            k = seen[sym]
            seen[sym] += 1
            rows += [(sym, ens)] if k == 0 else [(f"{sym}.{k}", ens), (f"{sym}-{k}", ens)]
    names = Counter(s for s, _ in rows)
    clash = [s for s, c in names.items() if c > 1]
    if clash:
        raise SystemExit(f"suffixed names collide with real symbols: {clash[:5]}")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["symbol", "ensembl"])
        w.writerows(rows)
    print(f"{len(seen)} symbols ({sum(c > 1 for c in seen.values())} repeated) -> {len(rows)} rows in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
