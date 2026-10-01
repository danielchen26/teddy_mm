#!/usr/bin/env python3
"""E6 external confirmation data: build a pack compatible with data/processed/cite from a public
CITE-seq dataset that this project has not used.

Dataset (the only one registered here): Hao et al. 2021, Cell 184:3573, "Integrated analysis of
multimodal single-cell data", GEO GSE164378, the 3' CITE-seq PBMC set (161,764 cells, 8 donors
P1-P8, 3 time points, 228 TotalSeq-A antibodies, WNN annotations celltype.l1/l2/l3).

Files (GEO, public; sizes from the series filelist.txt):
  GSM5008737_RNA_3P-{barcodes,features,matrix}  raw UMI counts, Cell Ranger 3.1.0, hg38, symbols only
  GSM5008738_ADT_3P-{barcodes,features,matrix}  raw ADT counts (Alevin)
  GSE164378_sc.meta.data_3P.csv.gz              donor, time, lane, Batch, celltype.l1/l2/l3
  GSE164378_Antibody_HTO_info.xlsx              specificity / clone / Ensembl id of each antibody
and, for symbol -> Ensembl, the feature table of the 10x Cell Ranger GRCh38-3.0.0 reference
(taken from the 10x pbmc_1k_v3 dataset). GEO lists the 33,538 genes in exactly that reference
order; the script checks this row by row (the 24 duplicated symbols carry Seurat's ".1" suffix)
and refuses to run if any row disagrees.

Pack (same keys as teddy_mm.data.prepare_cite, so load_prepared and 03_embed_rna.py read it):
  rna_*        raw integer UMI counts (float32 CSR), so the official TEDDY preprocessing can run
  rna_names / token_ids
               --gene-universe train (default): exactly the 12,052 genes and token ids of
               data/processed/cite, in the same order (genes absent from GRCh38-3.0.0 stay empty);
               --gene-universe vocab: every external gene whose Ensembl id is in the TEDDY vocab
  adt          [n, 134] in our panel order: per-cell CLR log1p(c / g_i), g_i = exp(mean over all
               228 measured antibodies (isotype controls included) of log1p(count)); NaN where the
               protein is not in the external panel (adt_present False)
  adt_counts   [n, 134] harmonised raw counts (NaN where missing); targets measured by two clones
               are the mean of the two clones' raw counts (declared, data-free rule)
  adt_source_counts / adt_source_names   all 228 raw ADT counts as published
  split        "test" for every cell: external cells are evaluation-only, nothing may be fit on them
  donors       "hao2021_P1".. ; sites "hao2021_Batch1|2"; cell_types = celltype.l2
  class_family / key_class   the registered v3 class map principles (B / T / NK / myeloid / OUT, plus unscored)
  extra        cell_types_l1, cell_types_l3, time, lane, barcodes, n_count_rna, n_count_adt

Stages (resumable; each logs progress): download -> parse RNA (cached under <raw>/_cache) ->
parse ADT (cached) -> assemble pack -> report (JSON + markdown). Re-running skips finished stages.

Examples
  python scripts/prepare_external_cite.py --download                       # full pack, all cells
  python scripts/prepare_external_cite.py --max-cells-per-donor 5000 \
      --out data/processed/external_hao2021_d5k                            # label-free subsample
  python scripts/prepare_external_cite.py --smoke 2000 --out /tmp/x --report-dir /tmp/x_report
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.data import size_factors  # noqa: E402

DATASET = "hao2021"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
# name -> (url, expected size in bytes from GEO filelist.txt, or None)
RAW_FILES = {
    "GSM5008737_RNA_3P-barcodes.tsv.gz": (f"{GEO}/samples/GSM5008nnn/GSM5008737/suppl/GSM5008737_RNA_3P-barcodes.tsv.gz", 683095),
    "GSM5008737_RNA_3P-features.tsv.gz": (f"{GEO}/samples/GSM5008nnn/GSM5008737/suppl/GSM5008737_RNA_3P-features.tsv.gz", 233354),
    "GSM5008737_RNA_3P-matrix.mtx.gz": (f"{GEO}/samples/GSM5008nnn/GSM5008737/suppl/GSM5008737_RNA_3P-matrix.mtx.gz", 1125353463),
    "GSM5008738_ADT_3P-barcodes.tsv.gz": (f"{GEO}/samples/GSM5008nnn/GSM5008738/suppl/GSM5008738_ADT_3P-barcodes.tsv.gz", 683095),
    "GSM5008738_ADT_3P-features.tsv.gz": (f"{GEO}/samples/GSM5008nnn/GSM5008738/suppl/GSM5008738_ADT_3P-features.tsv.gz", 1373),
    "GSM5008738_ADT_3P-matrix.mtx.gz": (f"{GEO}/samples/GSM5008nnn/GSM5008738/suppl/GSM5008738_ADT_3P-matrix.mtx.gz", 104745344),
    "GSE164378_sc.meta.data_3P.csv.gz": (f"{GEO}/series/GSE164nnn/GSE164378/suppl/GSE164378_sc.meta.data_3P.csv.gz", 3119866),
    "GSE164378_Antibody_HTO_info.xlsx": (f"{GEO}/series/GSE164nnn/GSE164378/suppl/GSE164378_Antibody_HTO_info.xlsx", None),
}
REF_TAR_URL = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz"
REF_FEATURES = "ref_10x_grch38_3.0.0/features_grch38_3.0.0.tsv.gz"

# ---------------------------------------------------------------------------------------------
# ADT harmonisation: our 134 names (data/processed/cite adt_names) -> Hao 3' feature names.
# Exact name matches are found automatically; the entries below are the non-trivial ones.
# ---------------------------------------------------------------------------------------------
# targets that Hao measured with two clones -> mean of the two raw counts (declared, data-free:
# the clones of our TotalSeq-B panel could not be checked, so neither clone is preferred)
TWO_CLONE = {
    "CD3": ["CD3-1", "CD3-2"],        # UCHT1, SK7
    "CD4": ["CD4-1", "CD4-2"],        # SK3, RPA-T4
    "CD56": ["CD56-1", "CD56-2"],     # 5.1H11, QA17A16
    "CD45": ["CD45-1", "CD45-2"],     # 2D1, HI30
    "CD11b": ["CD11b-1", "CD11b-2"],  # M1/70, ICRF44
    "CD38": ["CD38-1", "CD38-2"],     # HIT2, HB-7
    "CD44": ["CD44-1", "CD44-2"],     # IM7, BJ18
    "CD26": ["CD26-1", "CD26-2"],     # BA5b, BA5b (same clone listed twice)
}
ALIASES = {
    # ours "TCR" is TCR alpha/beta: on our training cells it is high on CD4/CD8 naive T and at
    # background on gdT TCRVD2+ (median CLR 1.29 / 0.93 vs 0.40; NK 0.38) -> Hao TCR-2 (TCR a/b,
    # IP26). Hao TCR-1 is TCR gamma/delta (B1) and is not ours.
    "TCR": ["TCR-2"],
    "integrinB7": ["Integrin-7"],    # Integrin beta7, FIB504
    "TCRVa7.2": ["TCR-V-7.2"],       # TCR Va7.2, 3C10
    "TCRVd2": ["TCR-V-2"],           # TCR Vd2, B6
}
NOT_MAPPED_NOTE = {
    "CD11a": "Hao has only CD11a/CD18 (LFA-1 complex, clone m24, conformation-dependent); not the same target",
    "HLA-A-B-C": "only in Hao's 5' panel (54 antibodies, other cells), not in the 3' panel",
    "CD5": "only in Hao's 5' panel, not in the 3' panel",
    "CD7": "only in Hao's 5' panel, not in the 3' panel",
}

# ---------------------------------------------------------------------------------------------
# Class map: the registered v3 principles (registration/registration_v3.json annotation_map and
# annotation_map_notes, built by bridge_anm/v3_build_registration.py CLASS_MAP) applied to the
# Hao celltype.l2 label names, written from the label names only:
#   B: mature / transitional B cells. T: CD4 / CD8 T, Treg, MAIT, gamma-delta T, double-negative
#   T, cycling T. NK: NK cells. myeloid: monocytes and conventional dendritic cells. OUT (no
#   B / T / NK / myeloid call is correct): progenitors, erythroid lineage, pDC, plasma cells and
#   plasmablasts, ILC. A label that is not a cell type (Doublet) is "unscored".
# class_family keeps the finer groups so another map can be derived without re-running.
# ---------------------------------------------------------------------------------------------
L2_CLASS = {  # Hao celltype.l2 -> (class_family, key_class, rule)
    "B naive": ("B_mature", "B", "mature B (registered: Naive CD20+ B -> B)"),
    "B intermediate": ("B_mature", "B", "mature B"),
    "B memory": ("B_mature", "B", "mature B"),
    "Plasmablast": ("plasma", "OUT", "antibody-secreting (registered: Plasmablast / Plasma cell -> OUT)"),
    "CD4 Naive": ("T", "T", "CD4 T"),
    "CD4 TCM": ("T", "T", "CD4 T"),
    "CD4 TEM": ("T", "T", "CD4 T"),
    "CD4 CTL": ("T", "T", "CD4 T"),
    "CD4 Proliferating": ("T_cycling", "T", "cycling T (registered: T prog cycling -> T)"),
    "Treg": ("T", "T", "Treg (registered: T reg -> T)"),
    "CD8 Naive": ("T", "T", "CD8 T"),
    "CD8 TCM": ("T", "T", "CD8 T"),
    "CD8 TEM": ("T", "T", "CD8 T"),
    "CD8 Proliferating": ("T_cycling", "T", "cycling T (registered: T prog cycling -> T)"),
    "MAIT": ("T", "T", "MAIT (registered: MAIT -> T)"),
    "gdT": ("T_gd", "T", "gamma-delta T (registered: gdT TCRVD2+ / gdT CD158b+ -> T)"),
    "dnT": ("T", "T", "double-negative T (registered: dnT -> T)"),
    "NK": ("NK", "NK", "NK (registered: NK -> NK)"),
    "NK_CD56bright": ("NK", "NK", "NK, CD56-bright subset"),
    "NK Proliferating": ("NK_cycling", "NK", "cycling NK (no BMMC label; as cycling T -> T)"),
    "ILC": ("ILC", "OUT", "innate lymphoid (registered: ILC / ILC1 -> OUT)"),
    "CD14 Mono": ("monocyte", "myeloid", "monocyte (registered: CD14+ Mono -> myeloid)"),
    "CD16 Mono": ("monocyte", "myeloid", "monocyte (registered: CD16+ Mono -> myeloid)"),
    "cDC1": ("cDC", "myeloid", "conventional DC (registered: cDC1 -> myeloid)"),
    "cDC2": ("cDC", "myeloid", "conventional DC (registered: cDC2 -> myeloid)"),
    "pDC": ("pDC", "OUT", "plasmacytoid DC (registered: pDC -> OUT)"),
    "ASDC": ("ASDC", "OUT", "AXL+ SIGLEC6+ DC: neither a monocyte nor a conventional DC (no BMMC label); "
             "OUT as pDC. Hao's celltype.l3 splits them into ASDC_mDC and ASDC_pDC (kept in cell_types_l3)"),
    "Eryth": ("erythroid", "OUT", "erythroid lineage (registered: Erythroblast / Normoblast / Reticulocyte -> OUT)"),
    "HSPC": ("progenitor", "OUT", "progenitor (registered: HSC / G/M prog / Lymph prog / MK/E prog -> OUT)"),
    "Platelet": ("platelet", "OUT", "megakaryocyte lineage, no B / T / NK / myeloid call is correct (no BMMC label)"),
    "Doublet": ("doublet", "unscored", "annotated doublet, not a cell type: unscored"),
}
KEY_CLASSES = ["B", "T", "NK", "myeloid", "OUT"]
KEY_COLUMNS = KEY_CLASSES + ["unscored"]


def annotation_map_addendum() -> tuple[dict[str, str], str]:
    """{celltype.l2: key_class} and the sha256 of its canonical JSON (sorted keys, no spaces)."""
    amap = {k: v[1] for k, v in sorted(L2_CLASS.items())}
    blob = json.dumps(amap, sort_keys=True, separators=(",", ":")).encode()
    return amap, hashlib.sha256(blob).hexdigest()


# proteins named in the task / v3 design (coverage is reported, nothing is selected by it)
CANONICAL_9 = ["CD19", "CD72", "CD22", "CD3", "CD2", "CD5", "CD16", "CD11c", "CD36"]
NK_PANEL = ["CD56", "CD94", "CD335"]
REQUIRED_MIN = [["CD3"], ["CD19", "CD20"], ["CD14"], ["CD16"], ["CD56"], ["CD94", "CD335"]]
# fallback when registration/registration_v3.json is absent (values of 2026-09-30)
REGISTERED_GATE = ["CD3", "CD19", "CD20", "CD56", "CD94", "CD335", "CD16", "CD14", "CD33", "CD11c", "CD71"]
REGISTERED_PRIMARY_PANEL = {"B": ["CD20", "CD22", "CD268"], "T": ["CD3", "CD2", "CD5"],
                            "NK": ["CD122", "CD94", "CD56"], "myeloid": ["CD172a", "CD11c", "CD62P"]}
REGISTRATION = ROOT / "registration" / "registration_v3.json"


def registered_protein_sets() -> dict[str, list[str]]:
    """Protein sets of the v3 registration whose external coverage is reported (nothing is
    selected by it): the gate proteins and every panel. Read from the registration if present."""
    sets: dict[str, list[str]] = {}
    src = "fallback constants"
    if REGISTRATION.exists():
        reg = json.loads(REGISTRATION.read_text())
        src = f"{REGISTRATION.name} sha256 {sha256_file(REGISTRATION)[:16]}"
        sets["gate"] = list(reg.get("gate", {}).get("proteins", []))
        for name, v in reg.get("panels", {}).items():
            if isinstance(v, dict) and isinstance(v.get("proteins"), dict):
                sets[f"panel {name}"] = [p for prots in v["proteins"].values() for p in prots]
    if not sets.get("gate"):
        sets["gate"] = REGISTERED_GATE
    if "panel primary" not in sets:
        sets["panel primary"] = [p for prots in REGISTERED_PRIMARY_PANEL.values() for p in prots]
    sets["_source"] = [src]
    return sets


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


# ------------------------------------------------------------------------------------ download
def _fetch(url: str, dest: Path, expected: int | None) -> None:
    """Resumable download (HTTP Range); refuses a file whose final size differs from GEO's."""
    part = dest.with_name(dest.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(url, headers={"User-Agent": "teddy_mm-prepare-external/1"})
    if have:
        req.add_header("Range", f"bytes={have}-")
    with urllib.request.urlopen(req, timeout=120) as r, part.open("ab" if have else "wb") as f:
        if have and r.status != 206:  # server ignored the range: start again
            f.seek(0)
            f.truncate()
        t0, n = time.time(), 0
        while True:
            block = r.read(1 << 22)
            if not block:
                break
            f.write(block)
            n += len(block)
            if time.time() - t0 > 30:
                log(f"  {dest.name}: {(have + n) / 1e6:.0f} MB")
                t0 = time.time()
    size = part.stat().st_size
    if expected is not None and size != expected:
        raise SystemExit(f"{dest.name}: got {size} bytes, GEO lists {expected}; rerun to resume")
    os.replace(part, dest)


def download(raw: Path) -> None:
    raw.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in RAW_FILES.items():
        dest = raw / name
        if dest.exists() and (expected is None or dest.stat().st_size == expected):
            log(f"have {name}")
            continue
        log(f"download {url}")
        _fetch(url, dest, expected)
    ref = raw / REF_FEATURES
    if not ref.exists():
        ref.parent.mkdir(parents=True, exist_ok=True)
        tar = ref.parent / "pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz"
        log(f"download {REF_TAR_URL} (feature table of the GRCh38-3.0.0 reference)")
        _fetch(REF_TAR_URL, tar, None)
        with tarfile.open(tar) as t:
            member = t.getmember("filtered_feature_bc_matrix/features.tsv.gz")
            with t.extractfile(member) as src:
                data = src.read()
        tmp = ref.with_name(ref.name + ".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, ref)
        tar.unlink()


def check_raw(raw: Path) -> dict:
    missing = [n for n in RAW_FILES if not (raw / n).exists()] + ([REF_FEATURES] if not (raw / REF_FEATURES).exists() else [])
    if missing:
        raise SystemExit(f"missing raw files in {raw}: {missing}\nRun with --download")
    info = {}
    for name, (url, expected) in RAW_FILES.items():
        size = (raw / name).stat().st_size
        if expected is not None and size != expected:
            raise SystemExit(f"{name}: {size} bytes, GEO lists {expected}; rerun --download")
        info[name] = {"url": url, "bytes": size, "geo_listed_bytes": expected}
    info[REF_FEATURES] = {"url": REF_TAR_URL + " (filtered_feature_bc_matrix/features.tsv.gz)",
                          "bytes": (raw / REF_FEATURES).stat().st_size}
    return info


def sha256_all(raw: Path, info: dict, cache: Path) -> None:
    """sha256 of each raw file, cached by (size, mtime)."""
    store = cache / "sha256.json"
    old = json.loads(store.read_text()) if store.exists() else {}
    for name in info:
        p = raw / name
        key = f"{p.stat().st_size}:{int(p.stat().st_mtime)}"
        if old.get(name, {}).get("key") != key:
            old[name] = {"key": key, "sha256": sha256_file(p)}
        info[name]["sha256"] = old[name]["sha256"]
    atomic_write_text(store, json.dumps(old, indent=1))


# ------------------------------------------------------------------------------ small readers
def read_lines_gz(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt") as f:
        return [line.rstrip("\n").split("\t") for line in f]


def read_xlsx_sheet1(path: Path) -> list[list[str]]:
    """First sheet of an .xlsx without openpyxl (shared strings + inline values)."""
    m = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(m + "si"):
            ss.append("".join(x.text or "" for x in si.iter(m + "t")))
    rows = []
    for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(m + "row"):
        vals = []
        for c in r.findall(m + "c"):
            v = c.find(m + "v")
            vals.append((ss[int(v.text)] if c.get("t") == "s" else v.text) if v is not None else "")
        rows.append(vals)
    return rows


def ensembl_for_hao_genes(raw: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Ensembl id of every Hao RNA feature, by position in the GRCh38-3.0.0 reference."""
    hao = read_lines_gz(raw / "GSM5008737_RNA_3P-features.tsv.gz")
    ref = read_lines_gz(raw / REF_FEATURES)
    if len(hao) != len(ref):
        raise SystemExit(f"Hao lists {len(hao)} genes, the GRCh38-3.0.0 reference {len(ref)}")
    seen: Counter = Counter()
    suffixed = 0
    for i, (h, r) in enumerate(zip(hao, ref)):
        sym = r[1]
        seen[sym] += 1
        want = sym if seen[sym] == 1 else f"{sym}.{seen[sym] - 1}"  # Seurat make.unique
        if h[0] != want:
            raise SystemExit(f"row {i}: Hao symbol {h[0]!r} != reference {want!r} ({r[0]})")
        if r[2] != "Gene Expression":
            raise SystemExit(f"row {i}: reference feature type {r[2]!r}")
        suffixed += h[0] != sym
    ens = np.array([r[0] for r in ref])
    sym = np.array([h[0] for h in hao])
    return ens, sym, {"n_genes": len(hao), "position_check": "all rows equal the reference symbol "
                      "(make.unique suffix on repeated symbols)", "n_make_unique_suffixed": int(suffixed)}


def mtx_header(path: Path) -> tuple[int, tuple[int, int, int]]:
    with gzip.open(path, "rt") as f:
        n = 0
        for line in f:
            n += 1
            if not line.startswith("%"):
                a, b, c = (int(x) for x in line.split())
                return n, (a, b, c)
    raise SystemExit(f"{path}: no size line")


def stream_mtx(path: Path, chunk: int, max_col: int | None = None):
    """Yield (row0, col0, value) int32 chunks of a gzip Matrix Market file (1-based -> 0-based).
    Decompression runs in a separate gzip process. With max_col, stop once a chunk passes it
    (10x/Seurat files are written column by column)."""
    n_header, shape = mtx_header(path)
    proc = subprocess.Popen(["gzip", "-dc", str(path)], stdout=subprocess.PIPE, bufsize=1 << 22)
    try:
        reader = pd.read_csv(proc.stdout, sep=" ", header=None, skiprows=n_header, names=["r", "c", "v"],
                             dtype={"r": np.int32, "c": np.int32, "v": np.int32}, chunksize=chunk, engine="c")
        for df in reader:
            r = df["r"].to_numpy() - 1
            c = df["c"].to_numpy() - 1
            v = df["v"].to_numpy()
            yield r, c, v
            if max_col is not None and c.size and c.min() >= max_col:
                break
    finally:
        proc.stdout.close()
        proc.kill()
        proc.wait()


# ------------------------------------------------------------------------------- RNA / ADT parse
def parse_rna(raw: Path, gene_lut: np.ndarray, n_genes_out: int, cache_dir: Path, tag: str,
              chunk: int, max_cells: int | None):
    """Cells x selected genes CSR of raw counts (gene_lut[hao_row] = output column or -1), plus
    each cell's UMI total and number of detected genes over all 33,538 genes (integrity check
    against the published nCount_RNA / nFeature_RNA)."""
    done = cache_dir / f"rna_{tag}.done.json"
    if max_cells is None and done.exists():
        log(f"RNA cache hit ({done.name})")
        meta = json.loads(done.read_text())
        data = np.load(cache_dir / f"rna_{tag}_data.npy", mmap_mode="r")
        indices = np.load(cache_dir / f"rna_{tag}_indices.npy", mmap_mode="r")
        indptr = np.load(cache_dir / f"rna_{tag}_indptr.npy")
        tot_all = np.load(cache_dir / f"rna_{tag}_total_all.npy")
        nfeat_all = np.load(cache_dir / f"rna_{tag}_nfeat_all.npy")
        x = sparse.csr_matrix((np.asarray(data), np.asarray(indices), indptr), shape=tuple(meta["shape"]))
        return x, tot_all, nfeat_all
    path = raw / "GSM5008737_RNA_3P-matrix.mtx.gz"
    _, (n_genes, n_cells, nnz) = mtx_header(path)
    if gene_lut.size != n_genes:
        raise SystemExit(f"gene lookup has {gene_lut.size} rows, matrix {n_genes}")
    n_keep_cells = n_cells if max_cells is None else min(max_cells, n_cells)
    # full run: one allocation of the header's nnz; smoke: start small and grow
    cap = nnz if max_cells is None else min(nnz, 2 * chunk)
    cell = np.empty(cap, dtype=np.int32)
    gene = np.empty(cap, dtype=np.int32)
    val = np.empty(cap, dtype=np.float32)
    tot_all = np.zeros(n_keep_cells, dtype=np.float64)
    nfeat_all = np.zeros(n_keep_cells, dtype=np.int64)
    k = seen = 0
    last_col = -1
    sorted_by_cell = True
    t0 = time.time()
    for r, c, v in stream_mtx(path, chunk, max_col=None if max_cells is None else n_keep_cells):
        seen += r.size
        if c.size:
            sorted_by_cell &= bool(c[0] >= last_col) and bool(np.all(np.diff(c) >= 0))
            last_col = int(c[-1])
        inside = c < n_keep_cells
        tot_all += np.bincount(c[inside], weights=v[inside], minlength=n_keep_cells)[:n_keep_cells]
        nfeat_all += np.bincount(c[inside], minlength=n_keep_cells)[:n_keep_cells]
        g = gene_lut[r]
        keep = (g >= 0) & inside
        m = int(keep.sum())
        if k + m > cell.size:
            new = max(k + m, 2 * cell.size)
            cell, gene, val = (np.concatenate([a[:k], np.empty(new - k, dtype=a.dtype)]) for a in (cell, gene, val))
        cell[k:k + m] = c[keep]
        gene[k:k + m] = g[keep]
        val[k:k + m] = v[keep]
        k += m
        log(f"  RNA {seen / 1e6:.0f}M / {nnz / 1e6:.0f}M entries read, {k / 1e6:.0f}M kept, "
            f"{time.time() - t0:.0f}s")
    if max_cells is None and seen != nnz:
        raise SystemExit(f"read {seen} entries, header says {nnz}")
    cell, gene, val = cell[:k], gene[:k], val[:k]
    if not np.all(val == np.round(val)) or val.min(initial=1) <= 0:
        raise SystemExit("RNA entries are not positive integers; expected raw UMI counts")
    if sorted_by_cell:
        indptr = np.zeros(n_keep_cells + 1, dtype=np.int64)
        np.cumsum(np.bincount(cell, minlength=n_keep_cells), out=indptr[1:])
        x = sparse.csr_matrix((val, gene, indptr), shape=(n_keep_cells, n_genes_out))
    else:
        log("  entries not ordered by cell; building CSR through COO")
        x = sparse.coo_matrix((val, (cell, gene)), shape=(n_keep_cells, n_genes_out)).tocsr()
    del cell, gene, val
    x.sort_indices()
    x.sum_duplicates()
    log(f"RNA parsed: {x.shape[0]} cells x {x.shape[1]} genes, nnz {x.nnz}, {time.time() - t0:.0f}s")
    if max_cells is None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        for name, arr in (("data", x.data), ("indices", x.indices), ("indptr", x.indptr),
                          ("total_all", tot_all), ("nfeat_all", nfeat_all)):
            tmp = cache_dir / f"rna_{tag}_{name}.tmp.npy"
            np.save(tmp, arr)
            os.replace(tmp, cache_dir / f"rna_{tag}_{name}.npy")
        atomic_write_text(done, json.dumps({"shape": list(x.shape), "nnz": int(x.nnz)}))
    return x, tot_all, nfeat_all


def parse_adt(raw: Path, cache_dir: Path, chunk: int, max_cells: int | None) -> tuple[np.ndarray, list[str]]:
    names = [r[0] for r in read_lines_gz(raw / "GSM5008738_ADT_3P-features.tsv.gz")]
    cached = cache_dir / "adt_source_counts.npy"
    if max_cells is None and cached.exists():
        log("ADT cache hit")
        return np.load(cached), names
    path = raw / "GSM5008738_ADT_3P-matrix.mtx.gz"
    _, (n_feat, n_cells, nnz) = mtx_header(path)
    if n_feat != len(names):
        raise SystemExit(f"ADT matrix has {n_feat} rows, features file {len(names)}")
    n_keep = n_cells if max_cells is None else min(max_cells, n_cells)
    out = np.zeros((n_keep, n_feat), dtype=np.float32)
    seen = 0
    for r, c, v in stream_mtx(path, chunk, max_col=None if max_cells is None else n_keep):
        seen += r.size
        keep = c < n_keep
        out[c[keep], r[keep]] = v[keep]
    if max_cells is None and seen != nnz:
        raise SystemExit(f"ADT: read {seen} entries, header says {nnz}")
    log(f"ADT parsed: {out.shape}")
    if max_cells is None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = cache_dir / "adt_source_counts.tmp.npy"
        np.save(tmp, out)
        os.replace(tmp, cached)
    return out, names


# --------------------------------------------------------------------------------- harmonise
def adt_mapping(ours: list[str], source: list[str]) -> list[dict]:
    src = set(source)
    rows = []
    for p in ours:
        if p in TWO_CLONE:
            feats, rule = TWO_CLONE[p], "mean of two clones' raw counts"
        elif p in ALIASES:
            feats, rule = ALIASES[p], "alias (same specificity, renamed)"
        elif p in src:
            feats, rule = [p], "exact name"
        else:
            feats, rule = [], "missing: " + NOT_MAPPED_NOTE.get(p, "not in Hao's 3' panel")
        for f in feats:
            if f not in src:
                raise SystemExit(f"mapping for {p} names {f}, which is not a Hao feature")
        rows.append({"ours": p, "source": feats, "rule": rule})
    used = [f for r in rows for f in r["source"]]
    dup = [f for f, n in Counter(used).items() if n > 1]
    if dup:
        raise SystemExit(f"source features mapped twice: {dup}")
    return rows


def harmonise_adt(src_counts: np.ndarray, src_names: list[str], mapping: list[dict]):
    col = {n: j for j, n in enumerate(src_names)}
    n = src_counts.shape[0]
    counts = np.full((n, len(mapping)), np.nan, dtype=np.float32)
    present = np.zeros(len(mapping), dtype=bool)
    for k, row in enumerate(mapping):
        if row["source"]:
            counts[:, k] = src_counts[:, [col[f] for f in row["source"]]].astype(np.float64).mean(axis=1)
            present[k] = True
    # per-cell CLR over the whole measured panel (isotype controls included), as Seurat/muon
    # CLR across features: g_i = exp(mean_j log1p(c_ij)); value = log1p(c / g_i)
    g = np.exp(np.log1p(src_counts.astype(np.float64)).mean(axis=1))
    g[g <= 0] = 1.0
    clr = np.full_like(counts, np.nan)
    clr[:, present] = np.log1p(counts[:, present].astype(np.float64) / g[:, None]).astype(np.float32)
    return clr, counts, present, g.astype(np.float32)


# ------------------------------------------------------------------------------- cell choice
def choose_cells(donors: np.ndarray, per_donor: int, seed: int) -> np.ndarray:
    """Label-free subsample: per_donor cells per donor, uniformly at random (seeded)."""
    rng = np.random.default_rng(seed)
    pick = []
    for d in sorted(set(donors)):
        idx = np.flatnonzero(donors == d)
        pick.append(idx if idx.size <= per_donor else rng.choice(idx, size=per_donor, replace=False))
    return np.sort(np.concatenate(pick))


# ------------------------------------------------------------------------------------- main
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw", type=Path, default=ROOT / "data/raw/external_hao2021_gse164378")
    p.add_argument("--out", type=Path, default=ROOT / "data/processed/external_hao2021")
    p.add_argument("--report-dir", type=Path, default=ROOT / "outputs/v3/e6_external_prep")
    p.add_argument("--train-pack", type=Path, default=ROOT / "data/processed/cite",
                   help="our pack: adt_names, rna_names and token_ids are read from it (nothing else)")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--medians", type=Path, default=ROOT / "data/reference/teddy_gene_medians.json")
    p.add_argument("--gene-universe", choices=("train", "vocab"), default="train")
    p.add_argument("--download", action="store_true", help="fetch missing raw files (GEO, 10x) first")
    p.add_argument("--max-cells-per-donor", type=int, default=0, help="label-free per-donor subsample; 0 = all")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--smoke", type=int, default=0, help="first N cells only (smoke test; no cache)")
    p.add_argument("--chunk", type=int, default=20_000_000, help="Matrix Market entries per parse chunk")
    p.add_argument("--no-compress", action="store_true", help="np.savez instead of savez_compressed")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args(argv)


def load_vocab(ckpt: Path) -> dict[str, int]:
    with (ckpt / "vocab.txt").open() as f:
        return {line.rstrip("\n"): i for i, line in enumerate(f)}


def main(argv=None) -> int:
    args = parse_args(argv)
    t_start = time.time()
    raw, out = args.raw, args.out
    cache = raw / "_cache"
    cache.mkdir(parents=True, exist_ok=True)
    pack_path = out / "cite_arrays.npz"
    if pack_path.exists() and not args.overwrite:
        log(f"{pack_path} exists; writing the report only (pass --overwrite to rebuild)")
    if args.download:
        download(raw)
    files = check_raw(raw)
    log("raw files present with GEO sizes; sha256 ...")
    sha256_all(raw, files, cache)
    smoke = args.smoke or None
    log("checking gene ids, barcodes, metadata and the antibody panel")

    # ---- our pack: names only
    with np.load(args.train_pack / "cite_arrays.npz", allow_pickle=False) as z:
        our_adt = [str(x) for x in z["adt_names"]]
        our_genes = np.array([str(x) for x in z["rna_names"]])
        our_tokens = np.asarray(z["token_ids"], dtype=np.int64)
    vocab = load_vocab(args.ckpt)
    if not np.array_equal(our_tokens, np.array([vocab[g] for g in our_genes], dtype=np.int64)):
        raise SystemExit("our pack's token_ids do not equal vocab[rna_names]")

    # ---- genes
    ens, sym, gene_check = ensembl_for_hao_genes(raw)
    if len(set(ens)) != len(ens):
        raise SystemExit("duplicate Ensembl ids in the reference")
    pos = {e: i for i, e in enumerate(ens)}
    lut = np.full(len(ens), -1, dtype=np.int32)
    if args.gene_universe == "train":
        rna_names = our_genes
        token_ids = our_tokens
        for j, g in enumerate(our_genes):
            if g in pos:
                lut[pos[g]] = j
        absent = [g for g in our_genes if g not in pos]
    else:
        keep = [i for i, e in enumerate(ens) if e in vocab]
        lut[keep] = np.arange(len(keep), dtype=np.int32)
        rna_names = ens[keep]
        token_ids = np.array([vocab[e] for e in rna_names], dtype=np.int64)
        absent = []
    gene_info = {
        **gene_check, "universe": args.gene_universe, "n_genes_out": int(len(rna_names)),
        "n_external_genes_in_vocab": int(sum(e in vocab for e in ens)),
        "n_external_genes_in_train_universe": int(len(set(ens) & set(our_genes))),
        "n_train_universe_genes": int(len(our_genes)),
        "n_out_genes_absent_from_external_reference": int(len(absent)),
        "absent_examples": absent[:20],
    }
    if args.medians.exists():
        med = json.loads(args.medians.read_text())
        gene_info["n_out_genes_with_teddy_median"] = int(sum(str(g) in med for g in rna_names))

    # ---- barcodes / metadata
    bc_rna = [r[0] for r in read_lines_gz(raw / "GSM5008737_RNA_3P-barcodes.tsv.gz")]
    bc_adt = [r[0] for r in read_lines_gz(raw / "GSM5008738_ADT_3P-barcodes.tsv.gz")]
    if bc_rna != bc_adt:
        raise SystemExit("RNA and ADT barcodes differ")
    meta = pd.read_csv(raw / "GSE164378_sc.meta.data_3P.csv.gz", index_col=0)
    missing_meta = [b for b in bc_rna[:1000] if b not in meta.index]
    if missing_meta or len(meta) != len(bc_rna) or set(meta.index) != set(bc_rna):
        raise SystemExit("metadata rows do not match the matrix barcodes")
    meta = meta.loc[bc_rna]
    l2 = meta["celltype.l2"].astype(str).to_numpy()
    unknown = sorted(set(l2) - set(L2_CLASS))
    if unknown:
        raise SystemExit(f"celltype.l2 labels without a class rule: {unknown}")

    # ---- ADT panel (needs no matrix)
    src_names = [r[0] for r in read_lines_gz(raw / "GSM5008738_ADT_3P-features.tsv.gz")]
    mapping = adt_mapping(our_adt, src_names)
    panel_xlsx = read_xlsx_sheet1(raw / "GSE164378_Antibody_HTO_info.xlsx")
    clone = {r[0].replace("_", "-"): {"specificity": r[4], "clone": r[5]} for r in panel_xlsx[1:] if len(r) > 5}
    for row in mapping:
        row["source_clones"] = [clone.get(f, {}).get("clone") for f in row["source"]]
        row["source_specificity"] = [clone.get(f, {}).get("specificity") for f in row["source"]]

    built_now = False
    if not pack_path.exists() or args.overwrite:
        built_now = True
        x_all, tot_all_genes, nfeat_all_genes = parse_rna(raw, lut, len(rna_names), cache,
                                                          f"{args.gene_universe}", args.chunk, smoke)
        adt_src_all, _ = parse_adt(raw, cache, args.chunk, smoke)
        n_all = x_all.shape[0]
        donors_all = np.array([f"{DATASET}_{d}" for d in meta["donor"].astype(str)])[:n_all]
        if args.max_cells_per_donor and not smoke:
            cells = choose_cells(donors_all, args.max_cells_per_donor, args.seed)
        else:
            cells = np.arange(n_all)
        if cells.size == n_all:  # all cells: no copy of the 2+ GB matrix
            x, adt_src = x_all, adt_src_all
        else:
            x, adt_src = x_all[cells], adt_src_all[cells]
        tot_all_genes, nfeat_all_genes = tot_all_genes[cells], nfeat_all_genes[cells]
        del x_all, adt_src_all
        m = meta.iloc[cells]
        clr, adt_counts, present, g = harmonise_adt(adt_src, src_names, mapping)
        l2c = m["celltype.l2"].astype(str).to_numpy()
        family = np.array([L2_CLASS[t][0] for t in l2c])
        key_class = np.array([L2_CLASS[t][1] for t in l2c])
        pres_counts = np.nan_to_num(adt_counts[:, present])
        arrays = dict(
            rna_data=x.data.astype(np.float32), rna_indices=x.indices.astype(np.int32),
            rna_indptr=x.indptr, rna_shape=np.array(x.shape),
            adt=clr, adt_counts=adt_counts, adt_present=present, adt_clr_scale=g,
            adt_source_counts=adt_src, adt_source_names=np.array(src_names, dtype="U32"),
            token_ids=np.asarray(token_ids, dtype=np.int64),
            split=np.full(len(cells), "test", dtype="U8"),
            donors=np.array([f"{DATASET}_{d}" for d in m["donor"].astype(str)], dtype="U32"),
            sites=np.array([f"{DATASET}_{b}" for b in m["Batch"].astype(str)], dtype="U32"),
            cell_types=l2c.astype("U64"),
            cell_types_l1=m["celltype.l1"].astype(str).to_numpy().astype("U64"),
            cell_types_l3=m["celltype.l3"].astype(str).to_numpy().astype("U64"),
            class_family=family.astype("U16"), key_class=key_class.astype("U16"),
            time=m["time"].to_numpy().astype(np.int16), lane=m["lane"].astype(str).to_numpy().astype("U16"),
            barcodes=np.array(m.index.astype(str), dtype="U32"), source_row=cells.astype(np.int64),
            n_count_rna=m["nCount_RNA"].to_numpy().astype(np.int32),
            n_count_adt=m["nCount_ADT"].to_numpy().astype(np.int32),
            rna_size_factor=size_factors(x), adt_size_factor=size_factors(sparse.csr_matrix(pres_counts)),
            adt_names=np.array(our_adt, dtype="U64"), rna_names=np.array(rna_names, dtype="U32"),
        )
        # integrity: the RNA matrix over all genes must reproduce the published nCount/nFeature_RNA
        tot = np.asarray(x.sum(axis=1)).ravel()
        if not (np.array_equal(tot_all_genes, arrays["n_count_rna"].astype(np.float64))
                and np.array_equal(nfeat_all_genes, m["nFeature_RNA"].to_numpy())):
            raise SystemExit("RNA totals over all genes differ from the published nCount_RNA / nFeature_RNA")
        # the GEO ADT matrix is integer but does not reproduce the published nCount_ADT exactly
        # (computed on another version of the ADT matrix); recorded, not fatal
        adt_tot = adt_src.sum(axis=1).astype(np.float64)
        adt_check = {
            "integer_counts": bool(np.all(adt_src == np.round(adt_src))),
            "share_cells_total_equal_published": float(np.mean(adt_tot == arrays["n_count_adt"])),
            "published_over_matrix_total": {q: float(v) for q, v in zip(
                ("min", "median", "max"), np.quantile(arrays["n_count_adt"] / np.maximum(adt_tot, 1), [0, 0.5, 1]))},
            "published_minus_matrix_nfeature": {q: float(v) for q, v in zip(
                ("min", "median", "max"), np.quantile(m["nFeature_ADT"].to_numpy() - (adt_src > 0).sum(axis=1), [0, 0.5, 1]))},
        }
        out.mkdir(parents=True, exist_ok=True)
        tmp = out / "cite_arrays.tmp.npz"
        log(f"writing {pack_path} ({'uncompressed' if args.no_compress else 'compressed'})")
        (np.savez if args.no_compress else np.savez_compressed)(tmp, **arrays)
        os.replace(tmp, pack_path)
        pack_meta = {
            "dataset": "Hao et al. 2021 PBMC CITE-seq 3' (GEO GSE164378)", "external": True,
            "n_cells": int(len(cells)), "n_rna_mapped": int(len(rna_names)), "n_adt": len(our_adt),
            "n_adt_present": int(present.sum()), "n_train": 0, "n_val": 0, "n_test": int(len(cells)),
            "test_sites": sorted(set(arrays["sites"].tolist())),
            "split_definition": "every external cell is evaluation-only ('test'); nothing may be fit on it",
            "gene_universe": args.gene_universe, "cell_selection":
                ("smoke: first %d cells" % smoke) if smoke else
                (f"label-free {args.max_cells_per_donor} per donor, seed {args.seed}" if args.max_cells_per_donor else "all"),
            "rna_values": "raw integer UMI counts (float32)",
            "adt_values": "per-cell CLR log1p(c/g_i), g_i over all 228 measured antibodies; NaN = not measured",
        }
        atomic_write_text(out / "meta.json", json.dumps(pack_meta, indent=2))
        stats = cell_stats(arrays, tot)
        stats["rna_integrity"] = "per-cell totals and detected genes over all 33,538 genes equal the published nCount_RNA / nFeature_RNA for every cell"
        stats["adt_integrity"] = adt_check
        atomic_write_text(out / "prepare_stats.json", json.dumps(stats, indent=1))
    else:
        pack_meta = json.loads((out / "meta.json").read_text())
        stats = json.loads((out / "prepare_stats.json").read_text())
        if pack_meta.get("gene_universe") != args.gene_universe:
            raise SystemExit(f"{pack_path} was built with --gene-universe {pack_meta.get('gene_universe')}; "
                             f"pass the same flag (or --overwrite)")

    manifest = {
        "script": "scripts/prepare_external_cite.py", "argv": sys.argv[1:], "git_commit": git_commit(),
        "raw_dir": str(raw), "pack": str(pack_path), "files": files, "genes": gene_info,
        "adt_mapping": mapping, "class_map": {k: {"family": v[0], "key_class": v[1], "rule": v[2]}
                                               for k, v in L2_CLASS.items()},
        "annotation_map_addendum": {"map": annotation_map_addendum()[0], "sha256_canonical_json": annotation_map_addendum()[1],
                                    "canonical_json": "json.dumps(map, sort_keys=True, separators=(',', ':'))"},
        "meta": pack_meta, "runtime_sec": round(time.time() - t_start, 1), "built_now": built_now,
    }
    manifest["token_probe"] = token_probe(out, args)
    atomic_write_text(out / "external_manifest.json", json.dumps(manifest, indent=1))
    write_report(args, manifest, stats, our_adt)
    log(f"done in {time.time() - t_start:.0f}s")
    return 0


def token_probe(out: Path, args, n: int = 300) -> dict | None:
    """Official TEDDY tokenisation (no model) of n random cells (seed 0): token counts, compared
    with the official embedding of our pack. Used only for the runtime estimate."""
    try:
        from teddy_mm.data import load_prepared
        from teddy_mm.teddy_encoder import load_gene_medians, median_factors, official_values, rank_encode_official
    except ImportError as e:  # torch missing: skip the probe
        return {"skipped": str(e)}
    if not args.medians.exists():
        return {"skipped": f"no medians at {args.medians}"}
    pack = load_prepared(out)
    rng = np.random.default_rng(0)
    idx = np.sort(rng.choice(pack["rna"].shape[0], size=min(n, pack["rna"].shape[0]), replace=False))
    vals = official_values(pack["rna"][idx].toarray(), median_factors(pack["rna_names"], load_gene_medians(args.medians)))
    _, attn = rank_encode_official(vals, pack["token_ids"], max_len=2048, pad_id=43810)
    nt = attn.sum(axis=1)
    res = {"n_cells": int(idx.size), "seed": 0, "median": float(np.median(nt)), "min": int(nt.min()),
           "share_at_2048": float(np.mean(nt == 2048))}
    ours = ROOT / "data/processed/cite_official/z_rna_ntokens.npy"
    if ours.exists():
        o = np.load(ours)
        res["ours_official"] = {"median": float(np.median(o)), "share_at_2048": float(np.mean(o == 2048))}
    log(f"token probe: median {res['median']:.0f} tokens, {res['share_at_2048']:.3f} at 2048")
    return res


def git_commit() -> str | None:
    try:
        sha = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "scripts/prepare_external_cite.py"],
                               capture_output=True, text=True).stdout.strip()
        return sha + ("+dirty" if dirty else "") if sha else None
    except OSError:
        return None


def cell_stats(a: dict, rna_tot: np.ndarray) -> dict:
    kc, fam, l2, don = a["key_class"], a["class_family"], a["cell_types"], a["donors"]
    per_donor = {d: {c: int(((don == d) & (kc == c)).sum()) for c in KEY_COLUMNS} for d in sorted(set(don))}
    nnz_per_cell = np.diff(a["rna_indptr"])
    return {
        "n_cells": int(len(kc)),
        "key_class": {c: int((kc == c).sum()) for c in KEY_COLUMNS},
        "class_family": dict(sorted(Counter(fam.tolist()).items())),
        "celltype_l2": dict(sorted(Counter(l2.tolist()).items())),
        "per_donor_key_class": per_donor,
        "donors": dict(sorted(Counter(don.tolist()).items())),
        "sites": dict(sorted(Counter(a["sites"].tolist()).items())),
        "time": {str(k): v for k, v in sorted(Counter(a["time"].tolist()).items())},
        "rna_mapped_umis_per_cell": {"median": float(np.median(rna_tot)), "p5": float(np.percentile(rna_tot, 5)),
                                     "p95": float(np.percentile(rna_tot, 95))},
        "rna_mapped_genes_per_cell": {"median": float(np.median(nnz_per_cell)), "min": int(nnz_per_cell.min()),
                                      "max": int(nnz_per_cell.max()),
                                      "share_at_least_2048": float(np.mean(nnz_per_cell >= 2048))},
        "fraction_mapped_umis_of_ncount_rna": float(np.median(rna_tot / np.maximum(a["n_count_rna"], 1))),
    }


def write_report(args, manifest: dict, stats: dict, our_adt: list[str]) -> None:
    rep = args.report_dir
    rep.mkdir(parents=True, exist_ok=True)
    mp = {r["ours"]: r for r in manifest["adt_mapping"]}
    have = lambda p: bool(mp[p]["source"])  # noqa: E731
    coverage = {
        "n_ours": len(our_adt), "n_mapped": sum(have(p) for p in our_adt),
        "n_exact": sum(mp[p]["rule"] == "exact name" for p in our_adt),
        "two_clone_mean": [p for p in our_adt if p in TWO_CLONE],
        "alias": [p for p in our_adt if p in ALIASES],
        "missing": [p for p in our_adt if not have(p)],
        "canonical_9": {p: have(p) for p in CANONICAL_9},
        "nk_panel": {p: have(p) for p in NK_PANEL},
        "required_minimum": {"/".join(alts): any(have(p) for p in alts) for alts in REQUIRED_MIN},
    }
    reg_sets = registered_protein_sets()
    coverage["registered_sets_source"] = reg_sets.pop("_source")[0]
    coverage["registered"] = {name: {p: have(p) for p in prots} for name, prots in reg_sets.items()}
    results = {"dataset": manifest["meta"]["dataset"], "pack": manifest["pack"], "meta": manifest["meta"],
               "genes": manifest["genes"], "adt_coverage": coverage, "cells": stats,
               "files": manifest["files"], "git_commit": manifest["git_commit"]}
    name = Path(manifest["pack"]).parent.name
    atomic_write_text(rep / f"prep_{name}.json", json.dumps(results, indent=1))

    st, gi = stats, manifest["genes"]
    lines = [
        f"# E6 external data: {manifest['meta']['dataset']}", "",
        f"Pack `{manifest['pack']}`; cell selection: {manifest['meta']['cell_selection']}; "
        f"gene universe: {gi['universe']}.", "",
        "## Cells", "",
        f"{st['n_cells']:,} cells, {len(st['donors'])} donors, sites {', '.join(st['sites'])}; "
        f"every cell is evaluation-only (`split = test`).", "",
        "| class | cells |", "|---|---:|",
        *[f"| {c} | {n:,} |" for c, n in st["key_class"].items()], "",
        "Families inside OUT: " + ", ".join(f"{k} {v:,}" for k, v in st["class_family"].items()
                                           if k in {f for f, c, _ in L2_CLASS.values() if c == "OUT"}) + ".", "",
        f"Annotation map (celltype.l2 -> class) sha256 of its canonical JSON: "
        f"`{manifest['annotation_map_addendum']['sha256_canonical_json']}` (map in external_manifest.json).", "",
        "Per donor:", "", "| donor | " + " | ".join(KEY_COLUMNS) + " |",
        "|---|" + "---:|" * len(KEY_COLUMNS),
        *[f"| {d} | " + " | ".join(f"{v[c]:,}" for c in KEY_COLUMNS) + " |"
          for d, v in st["per_donor_key_class"].items()], "",
        "## Genes", "",
        f"GEO gives symbols only; all {gi['n_genes']:,} rows equal the 10x GRCh38-3.0.0 reference in order "
        f"({gi['n_make_unique_suffixed']} repeated symbols carry a `.1` suffix), so each row takes the "
        f"reference's Ensembl id. {gi['n_external_genes_in_vocab']:,} external genes are in the TEDDY vocab; "
        f"{gi['n_external_genes_in_train_universe']:,} of our {gi['n_train_universe_genes']:,} training genes are "
        f"measured here ({gi['n_out_genes_absent_from_external_reference']} absent from the reference)."
        + (f" {gi['n_out_genes_with_teddy_median']:,} pack genes have a TEDDY median." if "n_out_genes_with_teddy_median" in gi else ""),
        "", f"Mapped UMIs per cell: median {st['rna_mapped_umis_per_cell']['median']:,.0f} "
        f"(median share of the published nCount_RNA {st['fraction_mapped_umis_of_ncount_rna']:.3f}); "
        f"mapped genes per cell: median {st['rna_mapped_genes_per_cell']['median']:,.0f}, "
        f"{st['rna_mapped_genes_per_cell']['share_at_least_2048']:.3f} of cells have at least 2,048.", "",
        "## ADT panel", "",
        f"{coverage['n_mapped']} of our {coverage['n_ours']} proteins are measured: {coverage['n_exact']} by exact "
        f"name, {len(coverage['alias'])} by alias ({', '.join(coverage['alias'])}), and "
        f"{len(coverage['two_clone_mean'])} as the mean of two clones' raw counts "
        f"({', '.join(coverage['two_clone_mean'])}).", "",
        f"Missing ({len(coverage['missing'])}): {', '.join(coverage['missing'])}.", "",
        "| set | covered |", "|---|---|",
        f"| canonical 9 | {sum(coverage['canonical_9'].values())}/9; missing "
        f"{', '.join(p for p, v in coverage['canonical_9'].items() if not v) or 'none'} |",
        f"| NK panel | {sum(coverage['nk_panel'].values())}/3; missing "
        f"{', '.join(p for p, v in coverage['nk_panel'].items() if not v) or 'none'} |",
        f"| required minimum | " + ", ".join(f"{k} {'yes' if v else 'NO'}" for k, v in coverage["required_minimum"].items()) + " |",
        *[f"| registered {k} | {sum(v.values())}/{len(v)}; missing {', '.join(q for q, ok in v.items() if not ok) or 'none'} |"
          for k, v in coverage["registered"].items()],
        "", f"Registered sets read from {coverage['registered_sets_source']}.", "",
    ]
    lines += integrity_lines(stats) + caveat_lines(coverage) + embed_lines(args, manifest) + source_lines(manifest)
    atomic_write_text(rep / f"PREP_{name}.md", "\n".join(lines) + "\n")
    log(f"report: {rep / f'PREP_{name}.md'}")


def integrity_lines(st: dict) -> list[str]:
    ai = st.get("adt_integrity", {})
    out = ["## Integrity checks", "", f"- RNA: {st.get('rna_integrity', 'n/a')}."]
    if ai:
        q, f = ai["published_over_matrix_total"], ai["published_minus_matrix_nfeature"]
        out.append(f"- ADT: integer counts {ai['integer_counts']}; the GEO ADT matrix does not reproduce the published "
                   f"nCount_ADT (equal for {ai['share_cells_total_equal_published']:.3f} of cells; published / matrix "
                   f"total min {q['min']:.3f}, median {q['median']:.3f}, max {q['max']:.3f}; nFeature difference "
                   f"{f['min']:+.0f} to {f['max']:+.0f}). The metadata was computed on another version of the ADT matrix; "
                   "the pack uses the GEO matrix.")
    return out + [""]


def caveat_lines(coverage: dict) -> list[str]:
    reg = coverage["registered"]
    miss = lambda d: ", ".join(p for p, ok in d.items() if not ok) or "none"  # noqa: E731
    return [
        "## Caveats and choices the registration must fix", "",
        "- **TEDDY pretraining exposure: label E6 'RNA possibly seen by TEDDY in pretraining'.** TEDDY was "
        "pretrained on CELLxGENE (1,399 datasets at download, arXiv 2503.03485). This dataset has been on CELLxGENE "
        "since 2022-07-15 (collection b0cf0afa-ec40-4d65-b570-ed4ceacc6813, dataset 'nygc multimodal pbmc', 161,764 "
        "cells; checked 2026-09-30 through the CELLxGENE curation API). Our heads, thresholds and rules were never fit "
        "on it. Every compared arm reads the same TEDDY embedding, so the exposure does not favour one arm, but "
        "absolute accuracies may be optimistic.",
        "- **Other chemistry and lab.** TotalSeq-A antibodies and PBMC from NYGC, against our TotalSeq-B bone marrow "
        "from the NeurIPS 2021 sites. Measured-ADT scales differ, and our pack's per-cell CLR scale g_i "
        "was taken over its full panel (about 0.94 of the geometric mean over the 134 stored proteins), here over "
        "Hao's 228. No threshold on measured ADT that was fit on our data transfers automatically. The external answer "
        "key should come from the published annotation (key_class); any protein check must be declared before evaluation.",
        f"- **Missing proteins.** Registered gate proteins not measured here: {miss(reg['gate'])}; registered "
        f"primary-panel proteins not measured here: {miss(reg.get('panel primary', {}))}. TEDDY's predicted values "
        "exist for all 134 proteins, so the rule, ANM and classifier arms run unchanged; only the measured-protein "
        "key loses them. Per E6's registered missing-protein rule, a gate clause that uses an absent protein is "
        "dropped (here: CD33 in the NK and myeloid clauses; CD94 leaves 3 NK markers with nk_min 2).",
        "- **Duplicate clones.** CD3, CD4, CD56, CD45, CD11b, CD38, CD44 and CD26 are the mean of two clones' raw counts "
        "(declared, data-free). Each clone is kept in adt_source_counts for a sensitivity check.",
        "- **Class map.** The registered v3 principles applied to the 31 celltype.l2 names. As registered: plasmablasts, "
        "pDC and ILC are OUT, and gdT and cycling T are T. Choices without a BMMC label: NK Proliferating -> NK, "
        "ASDC -> OUT (neither a monocyte nor a conventional DC), Platelet -> OUT, Doublet -> unscored. class_family "
        "keeps the finer groups, and cell_types_l3 is stored.",
        "- **How the map was made.** It was written from the label names only. Its first version followed the "
        "answer-key prototype (redesign/answer-key/keylib.py: plasma -> B, pDC -> myeloid, ILC -> NK). It was "
        "revised to the registered principles after smoke checks had read external values (per-class medians of 8 "
        "measured proteins on 2,000 cells, ADT totals, token counts). No value informed any mapping, but E6's rule "
        "'committed before any external ADT or RNA value is read' was not met in time order; an independent author "
        "can re-derive the map from the names and compare the sha256 above.",
        "- **Gene universe.** The default pack keeps our 12,052 training genes (the head was fit on embeddings over "
        "that universe); `--gene-universe vocab` builds the alternative. Pick one before embedding.",
        "- **Cells.** The full pack has every published cell; `--max-cells-per-donor 5000 --seed 0` gives a "
        "label-free, donor-balanced subsample for a cheaper embedding. Pick one before evaluation.",
        "- **Time points.** The metadata codes time as 0 / 2 / 7, while GEO's HTO sheet says day 0 / 3 / 7. The pack "
        "stores the metadata value. Each donor contributes 3 time points, so bootstrap over donors (8).",
        "",
    ]


def embed_lines(args, manifest: dict) -> list[str]:
    proc = Path(manifest["pack"]).parent
    out_dir = proc.parent / f"{proc.name}_official"
    n = int(manifest["meta"]["n_cells"])
    est = ""
    off = ROOT / "data/processed/cite_official/z_rna_manifest.json"
    if off.exists():
        m = json.loads(off.read_text())
        if m.get("embed_runtime_sec") and m.get("n_cells"):
            ms = 1000 * float(m["embed_runtime_sec"]) / int(m["n_cells"])
            est = (f" The official run took {ms:.1f} ms per cell on MPS ({int(m['n_cells']):,} cells), so expect at "
                   f"least about {n * ms / 3.6e6:.1f} h for these {n:,} cells.")
            tp = manifest.get("token_probe") or {}
            if "median" in tp and "ours_official" in tp:
                est += (f" Probably more: in a label-free sample of {tp['n_cells']} cells the official tokenisation gives a "
                        f"median of {tp['median']:,.0f} tokens per cell, {tp['share_at_2048']:.2f} of them at the 2,048 cap, "
                        f"against {tp['ours_official']['median']:,.0f} and {tp['ours_official']['share_at_2048']:.2f} on our "
                        "pack.")
    norm = os.path.normpath  # lexical only: resolving would replace the venv python by its base interpreter
    cmd = (f"cd {norm(ROOT)} && PYTHONUNBUFFERED=1 OMP_NUM_THREADS=6 {norm(sys.executable)} scripts/03_embed_rna.py \\\n"
           f"    --processed {norm(proc)} --out-dir {norm(out_dir)} \\\n"
           f"    --ckpt {norm(args.ckpt)} --preprocessing official --medians {norm(args.medians)} \\\n"
           "    --seq-len 2048 --pooling gene-mean --save-poolings first --length-buckets --batch-size 32 \\\n"
           "    --autocast fp16 --save-layer-means --shard-size 5000 --device mps --threads 6")
    return ["## Official embedding (run later, GPU)", "",
            "These are the settings of the official run (scripts/rerun_official.sh step_embed). The run is resumable: "
            "shards go to <out-dir>/shards, a rerun skips finished shards, and `--time-budget-sec N` stops cleanly with "
            "exit 75." + est, "", "```", cmd, "```", ""]


def source_lines(manifest: dict) -> list[str]:
    out = ["## Sources", "", "| file | bytes | sha256 |", "|---|---:|---|"]
    for name, f in manifest["files"].items():
        out.append(f"| [{name}]({f['url'].split(' ')[0]}) | {f['bytes']:,} | `{f.get('sha256', '')[:16]}...` |")
    return out + ["", "Hao Y, Hao S, et al. Integrated analysis of multimodal single-cell data. Cell 184:3573 (2021), "
                  "doi:10.1016/j.cell.2021.04.048; GEO GSE164378 (public).", ""]


if __name__ == "__main__":
    sys.exit(main())
