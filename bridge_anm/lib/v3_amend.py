"""Amendment A1 to the v3 pre-registration (registration/amendment_A1.json).

A1 was written after an adversarial review of ``registration_v3.json`` and before any
site4 (test) evaluation. It leaves the frozen registration file untouched and records,
with its own sha256, the corrections that the review found necessary:

* A1.1 matched coverage: ties are broken by a registered seeded permutation of the cell
  ids (not by row index, which is ordered by donor), and the selection is recomputed in
  every bootstrap replicate;
* A1.2 E4 channel 2 reads a panel-free renormalisation of the stored ADT (the stored
  values are a per-cell CLR whose denominator includes the panel proteins);
* A1.3 E1.3 (why this call): the exact closed form of ANM's leave-one-out (removing an
  event deletes its event site, so the field gain changes from G(n) to G(n-1)), and tied
  top markers are reported as ties;
* A1.4 to A1.10: interpretation and endpoint pins listed in the JSON.

Builders load it with ``load_registration_amended()``, which verifies both hashes and
returns the registration with the amended experiment fields applied. Nothing here
chooses a threshold: the few numbers A1 adds were computed on train/val only by
``bridge_anm/v3_build_amendment.py`` (covered by ``bridge_anm/v3_leakage_check.py``).
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from . import v3_key as vk

REPO = Path(__file__).resolve().parents[2]
DEFAULT_AMENDMENT = REPO / "registration" / "amendment_A1.json"


# --------------------------------------------------------------------------- loading
def load_amendment(path: Path | str | None = None, *, verify_hash: bool = True,
                   registration_path: Path | str | None = None) -> dict[str, Any]:
    """Load amendment A1; check its own sha256 and that it amends the registration on disk."""
    path = Path(path) if path is not None else DEFAULT_AMENDMENT
    if verify_hash:
        hf = path.with_name(path.name + ".sha256")
        if not hf.exists():
            raise FileNotFoundError(f"amendment hash file missing: {hf}")
        recorded = hf.read_text().split()[0].strip()
        actual = vk.sha256_file(path)
        if recorded != actual:
            raise ValueError(f"amendment {path} sha256 {actual} != recorded {recorded}")
    amend = json.loads(path.read_text())
    if verify_hash:
        rp = Path(registration_path) if registration_path is not None else vk.DEFAULT_REGISTRATION
        reg_sha = vk.sha256_file(rp)
        if amend["amends"]["registration_sha256"] != reg_sha:
            raise ValueError(f"amendment amends registration {amend['amends']['registration_sha256']}, "
                             f"but {rp} has sha256 {reg_sha}")
    return amend


def _set_path(d: dict, path: str, value: Any) -> None:
    """Set d[a][b][c] = value for path 'a/b/c' (endpoint names such as 'E1.1c' contain dots)."""
    keys = path.split("/")
    cur = d
    for k in keys[:-1]:
        if k not in cur or not isinstance(cur[k], dict):
            raise KeyError(f"override path {path!r}: {k!r} missing")
        cur = cur[k]
    cur[keys[-1]] = value


def apply_overrides(reg: dict[str, Any], amend: dict[str, Any]) -> dict[str, Any]:
    """Return a deep copy of the registration with the amendment's experiment overrides applied."""
    out = copy.deepcopy(reg)
    for path, value in amend["experiments_overrides"].items():
        _set_path(out["experiments"], path, value)
    out.setdefault("amendments", []).append({k: amend[k] for k in ("amendment_id", "amends", "written_before")})
    out["amendment_A1"] = amend
    return out


def load_registration_amended(registration_path: Path | str | None = None,
                              amendment_path: Path | str | None = None) -> dict[str, Any]:
    """The registration every v3 builder uses: hash-checked registration plus hash-checked A1."""
    reg = vk.load_registration(registration_path)
    amend = load_amendment(amendment_path, registration_path=registration_path)
    return apply_overrides(reg, amend)


def amendment_sha256(path: Path | str | None = None) -> str:
    return vk.sha256_file(Path(path) if path is not None else DEFAULT_AMENDMENT)


# --------------------------------------------------------------------------- A1.1 matched coverage
def tie_break_rank(cell_ids: Iterable[int], amend: dict[str, Any]) -> np.ndarray:
    """Rank of each global cell id in the registered seeded permutation (label-free, same for all methods)."""
    tb = amend["matched_coverage"]["tie_break"]
    perm = np.random.default_rng(int(tb["seed"])).permutation(int(tb["n_cells"]))
    rank = np.empty_like(perm)
    rank[perm] = np.arange(perm.size)
    return rank[np.asarray(list(cell_ids) if not isinstance(cell_ids, np.ndarray) else cell_ids, dtype=np.int64)]


def coverage_select(conf: np.ndarray, cell_ids: np.ndarray, coverage: float,
                    amend: dict[str, Any]) -> np.ndarray:
    """Boolean mask of the ceil(coverage * N) most confident cells; ties by the registered permutation.

    ``cell_ids`` are global row indices (a bootstrap replicate may repeat an id; repeats share the
    rank and are interchangeable). NaN confidence ranks last."""
    conf = np.asarray(conf, dtype=np.float64)
    n = conf.size
    k = int(math.ceil(float(coverage) * n - 1e-9))
    key = np.where(np.isnan(conf), -np.inf, conf)
    rank = tie_break_rank(np.asarray(cell_ids, dtype=np.int64), amend)
    order = np.lexsort((rank, -key))  # last key sorts first: confidence desc, then permutation rank
    mask = np.zeros(n, dtype=bool)
    mask[order[:k]] = True
    return mask


# --------------------------------------------------------------------------- A1.3 leave-one-out
def loo_ratio(n_events: int, field: dict[str, Any]) -> float:
    """G(n-1) / G(n): the called action's field gain after one of its n event sites is deleted."""
    n = int(n_events)
    if n <= 1:
        return 1.0  # nothing is left: the called score becomes 0 whatever the ratio
    return vk.field_gain(n - 1, field) / vk.field_gain(n, field)


def top_marker_sets(v_class: np.ndarray) -> list[np.ndarray]:
    """Per cell, the indices (panel order) of every marker whose evidence equals the class maximum."""
    v = np.asarray(v_class, dtype=np.float64)
    m = v.max(axis=1, keepdims=True)
    return [np.where(row)[0] for row in (v == m)]


def loo_flip_exact(S: np.ndarray, v_called: np.ndarray, called_idx: np.ndarray, bar: float,
                   field: dict[str, Any]) -> np.ndarray:
    """Exact re-coding of ANM's leave-one-out flip for the top marker (events at t = 0).

    Removing the top event deletes its event site, so the called action's score becomes
    G(n-1) * (n * S_called - v_top) while every other action keeps G(n) * n * S_k and the
    readout threshold stays G(n) * n * bar. In rule units: the call flips iff
    rho * (S_called - v_top / n) < max(bar, max_{k != called} S_k), rho = G(n-1) / G(n).

    S: [N, 4] rule scores; v_called: [N, n] evidence of the called class's panel; called_idx: [N]."""
    S = np.asarray(S, dtype=np.float64)
    v = np.asarray(v_called, dtype=np.float64)
    ci = np.asarray(called_idx, dtype=np.int64)
    n = v.shape[1]
    rho = loo_ratio(n, field)
    s_called = S[np.arange(S.shape[0]), ci]
    rem = s_called - v.max(axis=1) / n
    other = S.copy()
    other[np.arange(S.shape[0]), ci] = -np.inf
    return rho * rem < np.maximum(float(bar), other.max(axis=1))


def loo_flip_registered_form(S: np.ndarray, v_called: np.ndarray, called_idx: np.ndarray,
                             bar: float) -> np.ndarray:
    """The form experiments_v3.json wrote (S_called - v_top / n below the bar or another score).

    It is the exact form for *zeroing* the event, not for deleting it; kept only to report how
    many cells the A1.3 correction moves."""
    S = np.asarray(S, dtype=np.float64)
    v = np.asarray(v_called, dtype=np.float64)
    ci = np.asarray(called_idx, dtype=np.int64)
    rem = S[np.arange(S.shape[0]), ci] - v.max(axis=1) / v.shape[1]
    other = S.copy()
    other[np.arange(S.shape[0]), ci] = -np.inf
    return rem < np.maximum(float(bar), other.max(axis=1))


# --------------------------------------------------------------------------- A1.2 E4 channel 2
def channel2_inputs(adt: np.ndarray, adt_names: Iterable[str], panel_proteins: Iterable[str],
                    scale: float = 1e4) -> tuple[np.ndarray, list[str]]:
    """Panel-free channel-2 inputs from the stored per-cell CLR ADT.

    Stored values are v = log1p(x / g_c) with g_c computed per cell from all proteins, the panel
    included. u = expm1(v) = x / g_c, so u_q / sum_{r non-panel} u_r = x_q / sum_r x_r: the
    returned w_q = log1p(scale * x_q / sum_{r non-panel} x_r) depends on the non-panel counts only.
    A cell with no non-panel counts gets w = 0."""
    names = [str(x) for x in adt_names]
    panel = set(map(str, panel_proteins))
    cols = [i for i, n in enumerate(names) if n not in panel]
    u = np.expm1(np.asarray(adt, dtype=np.float64)[:, cols])
    u = np.maximum(u, 0.0)
    tot = u.sum(axis=1, keepdims=True)
    w = np.log1p(scale * np.divide(u, tot, out=np.zeros_like(u), where=tot > 0))
    return w, [names[i] for i in cols]
