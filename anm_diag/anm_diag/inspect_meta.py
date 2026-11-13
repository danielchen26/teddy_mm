"""Look at metadata only, before sealing.

  python -m anm_diag inspect --paths PATHS.yaml [--show col1,col2]

Prints, for the RNA and ADT files a paths file names: the shape, every obs column with its type and number
of distinct values (no values), the value counts of the design columns only (obs_target, obs_stratum and
any listed in --show), the var columns, the first var names and which var column holds Ensembl ids, and
the protein panel. It never reads a count matrix, so nothing about a readout or a state is seen before the
declaration is sealed. Do not --show columns derived from the data (QC metrics, response classes such as a
mixscape KO/NP label): their values carry outcome information.
"""
import re

from .adapters.base import load_paths

ENS = re.compile(r"^ENS[A-Z]*G\d{6,}")


def _obs_summary(adata, show, out):
    out.append(f"  cells {adata.n_obs}, features {adata.n_vars}")
    out.append("  obs columns (name, type, distinct values):")
    for c in adata.obs.columns:
        s = adata.obs[c]
        out.append(f"    {c:28s} {str(s.dtype):12s} {s.nunique(dropna=False)}")
    for c in show:
        if c in adata.obs.columns:
            vc = adata.obs[c].astype(str).value_counts()
            out.append(f"  value counts of '{c}' ({len(vc)} values; first 60):")
            for k, v in vc.head(60).items():
                out.append(f"    {k:40s} {v}")
        else:
            out.append(f"  column '{c}' not found")


def _var_summary(adata, out, full=False):
    out.append(f"  var columns: {list(adata.var.columns)}")
    names = [str(x) for x in adata.var_names]
    out.append(f"  var names (first 8): {names[:8]}")
    if full:
        out.append(f"  all var names: {names}")
    for c in adata.var.columns:
        vals = adata.var[c].astype(str).head(200)
        frac = sum(bool(ENS.match(v)) for v in vals) / max(len(vals), 1)
        if frac > 0.8:
            out.append(f"  Ensembl ids in var column '{c}'")
    if sum(bool(ENS.match(n)) for n in names[:200]) > 160:
        out.append("  Ensembl ids are the var names themselves")


def inspect(paths, show=()):
    import anndata as ad
    p = load_paths(paths)
    cols = [c for c in (p.get("obs_target"), p.get("obs_stratum")) if c] + [c for c in show if c]
    out = []
    for key, full in (("rna_h5ad", False), ("adt_h5ad", True)):
        if not p.get(key):
            continue
        a = ad.read_h5ad(p[key], backed="r")
        out.append(f"{key}: {p[key]}")
        _obs_summary(a, list(dict.fromkeys(cols)), out)
        _var_summary(a, out, full=full)
        a.file.close()
    return "\n".join(out)
