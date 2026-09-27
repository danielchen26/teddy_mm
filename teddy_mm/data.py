from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def _as_str_series(obs: pd.DataFrame, candidates: list[str]) -> pd.Series:
    for name in candidates:
        if name in obs.columns:
            return obs[name].astype(str)
    raise KeyError(f"none of {candidates} found in obs: {list(obs.columns)[:30]}")


def split_cite_modalities(adata):
    if "feature_types" in adata.var.columns:
        types = adata.var["feature_types"].astype(str)
        gex_mask = types.str.contains("GEX|Gene Expression|RNA", case=False, regex=True)
        adt_mask = types.str.contains("ADT|Antibody|Protein", case=False, regex=True)
        if gex_mask.any() and adt_mask.any():
            return adata[:, gex_mask].copy(), adata[:, adt_mask].copy()
    if "protein" in adata.obsm:
        raise ValueError("unexpected CITE layout: protein in obsm; add a parser if needed")
    raise ValueError("cannot find GEX/ADT feature_types in var")


def gene_ids_from_var(var: pd.DataFrame) -> np.ndarray:
    for col in ("gene_id", "ensembl_id", "ENSEMBL", "gene_ids"):
        if col in var.columns:
            return var[col].astype(str).to_numpy()
    idx = var.index.astype(str).to_numpy()
    if np.any(np.char.startswith(idx, "ENSG")):
        return idx
    return idx


def donor_site_split(
    obs: pd.DataFrame,
    test_sites: list[str],
    val_fraction: float,
    seed: int,
) -> np.ndarray:
    site = _as_str_series(obs, ["Site", "site", "site_id"])
    donor = _as_str_series(obs, ["DonorID", "donor_id", "donor", "Donor"])
    site_norm = site.str.lower().str.replace(" ", "", regex=False)
    test_norm = {s.lower().replace(" ", "") for s in test_sites}
    split = np.array(["train"] * len(obs), dtype=object)
    split[site_norm.isin(test_norm)] = "test"
    train_donors = sorted(donor[split == "train"].unique())
    rng = np.random.default_rng(seed)
    n_val = max(1, int(round(len(train_donors) * val_fraction)))
    val_donors = set(rng.choice(train_donors, size=n_val, replace=False))
    split[(split == "train") & donor.isin(val_donors)] = "val"
    return split


def size_factors(counts) -> np.ndarray:
    if sparse.issparse(counts):
        tot = np.asarray(counts.sum(axis=1)).ravel()
    else:
        tot = np.asarray(counts).sum(axis=1)
    med = np.median(tot[tot > 0]) if np.any(tot > 0) else 1.0
    sf = tot / max(med, 1e-6)
    sf[sf <= 0] = 1.0
    return sf.astype(np.float32)


def to_csr_float32(x):
    if sparse.issparse(x):
        return x.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(x, dtype=np.float32))


def prepare_cite(
    raw_path: Path,
    out_dir: Path,
    vocab: dict[str, int],
    test_sites: list[str],
    val_fraction: float,
    seed: int,
    max_cells: int | None = None,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    import anndata as ad

    adata = ad.read_h5ad(raw_path)
    if max_cells is not None and adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
        adata = adata[idx].copy()

    rna, adt = split_cite_modalities(adata)
    split = donor_site_split(rna.obs, test_sites, val_fraction, seed)
    gene_ids = gene_ids_from_var(rna.var)
    mapped = np.array([gid in vocab for gid in gene_ids])
    if mapped.sum() < 100:
        # last resort: strip version suffix ENSG0001.1 -> ENSG0001
        gene_ids = np.array([g.split(".")[0] for g in gene_ids])
        mapped = np.array([gid in vocab for gid in gene_ids])
    if mapped.sum() < 100:
        raise RuntimeError(
            f"only {mapped.sum()} RNA features map to TEDDY vocab; "
            "check that genes are Ensembl IDs"
        )
    rna = rna[:, mapped].copy()
    token_ids = np.array([vocab[g.split(".")[0] if g.split(".")[0] in vocab else g] for g in gene_ids[mapped]], dtype=np.int64)

    rna_x = to_csr_float32(rna.X)
    adt_x = to_csr_float32(adt.X)
    cell_types = (
        rna.obs["cell_type"].astype(str).to_numpy()
        if "cell_type" in rna.obs
        else np.array(["unknown"] * rna.n_obs)
    )
    donors = _as_str_series(rna.obs, ["DonorID", "donor_id", "donor", "Donor"]).to_numpy()
    sites = _as_str_series(rna.obs, ["Site", "site", "site_id"]).to_numpy()

    np.savez_compressed(
        out_dir / "cite_arrays.npz",
        rna_data=rna_x.data,
        rna_indices=rna_x.indices,
        rna_indptr=rna_x.indptr,
        rna_shape=np.array(rna_x.shape),
        adt=np.asarray(adt_x.todense(), dtype=np.float32),
        token_ids=token_ids,
        split=split.astype("U8"),
        donors=donors.astype("U32"),
        sites=sites.astype("U32"),
        cell_types=cell_types.astype("U64"),
        rna_size_factor=size_factors(rna_x),
        adt_size_factor=size_factors(adt_x),
        adt_names=adt.var_names.astype(str).to_numpy().astype("U64"),
        rna_names=np.array([g.split(".")[0] for g in gene_ids[mapped]], dtype="U32"),
    )
    meta = {
        "n_cells": int(rna.n_obs),
        "n_rna_mapped": int(mapped.sum()),
        "n_adt": int(adt.n_vars),
        "n_train": int((split == "train").sum()),
        "n_val": int((split == "val").sum()),
        "n_test": int((split == "test").sum()),
        "test_sites": test_sites,
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    return meta


def _as_fixed_unicode(arr: np.ndarray, dtype: str = "U64") -> np.ndarray:
    """Coerce object/bytes string arrays to fixed unicode for safe npz reload."""
    if isinstance(arr, np.ndarray) and arr.dtype.kind in ("U", "S"):
        return arr.astype(dtype)
    return np.asarray([str(x) for x in np.asarray(arr).ravel()], dtype=dtype).reshape(np.shape(arr))


def load_prepared(processed_dir: Path) -> dict:
    npz_path = processed_dir / "cite_arrays.npz"
    try:
        z = np.load(npz_path, allow_pickle=False)
        # Touch keys that may be object dtype to fail early into the pickle path.
        _ = z["adt_names"]
    except ValueError:
        z = np.load(npz_path, allow_pickle=True)
    rna = sparse.csr_matrix(
        (z["rna_data"], z["rna_indices"], z["rna_indptr"]),
        shape=tuple(z["rna_shape"]),
    )
    return {
        "rna": rna,
        "adt": z["adt"],
        "token_ids": z["token_ids"],
        "split": z["split"],
        "donors": z["donors"],
        "sites": z["sites"],
        "cell_types": z["cell_types"],
        "rna_size_factor": z["rna_size_factor"],
        "adt_size_factor": z["adt_size_factor"],
        "adt_names": _as_fixed_unicode(z["adt_names"], "U64"),
        "rna_names": z["rna_names"],
        "meta": json.loads((processed_dir / "meta.json").read_text()),
    }
