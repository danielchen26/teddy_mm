"""Seal a declaration before any outcome is read, and check the seal when the test runs.

  python -m anm_diag seal DECLARATION --paths PATHS.yaml --out SEALED.yaml

The sealed copy is the filled declaration (defaults included, so later changes to the defaults cannot
alter it) with status: sealed, the sha256 of every input file the paths name (data_hashes), the sha256
of the anm_diag source (code_hash) and the time of sealing. run() refuses a sealed declaration whose
files or code no longer match. Only metadata may be looked at before sealing (`python -m anm_diag
inspect`): obs columns and labels, var names, file hashes; never a readout or a state.

Files made after sealing from sealed inputs by the frozen system (the TEDDY states from layer_states.py)
are not in the seal; the TEDDY checkpoint they come from is (paths key teddy_ckpt).
"""
import copy
import datetime
import hashlib
from pathlib import Path

from .adapters.base import file_sha256, load_paths
from .spec import declaration_hash, load_declaration

DATA_KEYS = ("rna_h5ad", "adt_h5ad", "counts", "var", "prep", "cells", "vocab", "medians", "head_ckpt",
             "teddy_ckpt", "symbol_map", "z_files")


class SealError(RuntimeError):
    pass


def code_hash():
    root = Path(__file__).resolve().parent
    h = hashlib.sha256()
    for f in sorted(root.rglob("*.py")):
        h.update(str(f.relative_to(root)).encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def data_hashes(paths):
    p = load_paths(paths)
    out = {}
    for k in DATA_KEYS:
        v = p.get(k)
        files = v if isinstance(v, list) else [v]
        for i, f in enumerate(files):
            if isinstance(f, str) and Path(f).is_file():
                out[k if len(files) == 1 else f"{k}[{i}]"] = file_sha256(f)
    return out


def seal(declaration, paths, out=None):
    import yaml
    decl = copy.deepcopy(load_declaration(declaration))
    dh = data_hashes(paths)
    if not dh:
        raise SealError("no input file found to hash: check the paths file")
    decl.pop("seal_digest", None)
    decl.update({"status": "sealed", "data_hashes": dh, "code_hash": code_hash(),
                 "sealed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")})
    decl["seal_digest"] = _digest(decl)
    if out is not None:
        Path(out).write_text("# SEALED: do not edit. A changed design is a new declaration.\n"
                             + yaml.safe_dump(decl, sort_keys=False, allow_unicode=True))
    return decl


def _digest(decl):
    d = {k: v for k, v in decl.items() if k != "seal_digest"}
    return declaration_hash(d)


def verify(decl, paths):
    """Raise SealError if a sealed declaration was edited, or its data or code changed. Returns the checked keys."""
    if decl.get("status") != "sealed":
        return []
    if decl.get("seal_digest") != _digest(decl):
        raise SealError("sealed declaration was edited after sealing (seal_digest does not match)")
    now = data_hashes(paths) if paths is not None else {}
    bad = [k for k, h in decl["data_hashes"].items() if now.get(k) != h]
    if bad:
        raise SealError(f"sealed declaration: input files changed or missing since sealing: {bad}")
    if decl["code_hash"] != code_hash():
        raise SealError("sealed declaration: anm_diag code changed since sealing")
    return sorted(decl["data_hashes"])
