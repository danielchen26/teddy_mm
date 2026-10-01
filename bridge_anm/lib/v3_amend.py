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

Amendment A2 (registration/amendment_A2.json, also before any site4 evaluation) follows A1
and only makes verdicts stricter or cheaper to estimate without changing the estimand:

* A2.1 E3 H3a: a readout's H3a win also needs (b) recovery on flagged pairs, GR_R(flagged) -
  GR_head(flagged) > 0, and (c) the scale-free Dlog_R > 0, each with its two-stage 95% interval
  above 0 and positive in each primary donor; a D_R win without them is "not supported (scale
  artefact)" (``a2_positive_condition``, ``a2_h3a_verdict``, ``a2_e3_falsification``);
* A2.2 E3 R2 normaliser: the training q95 of R2's predictions is estimated on a registered,
  seeded, label-free uniform sample of 10,000 split=train cells (``a2_r2_q95_sample``);
* A2.3 E2: the common decision rule (effect beyond the margin overall and in each primary
  donor) is binding for falsification criteria (i) and (ii) (``a2_e2_criterion_i``,
  ``a2_e2_criterion_ii``, ``a2_e2_verdict``).

Builders load both with ``load_registration_amended()``, which verifies the registration, A1
and A2 hashes and returns the registration with A1's and then A2's experiment fields applied.
Nothing here chooses a threshold: the few numbers A1 adds were computed on train/val only by
``bridge_anm/v3_build_amendment.py`` (covered by ``bridge_anm/v3_leakage_check.py``); A2's by
``bridge_anm/v3_build_amendment_A2.py`` (train rows only; its own site4-poisoned rebuild).
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import numpy as np

from . import v3_key as vk

REPO = Path(__file__).resolve().parents[2]
DEFAULT_AMENDMENT = REPO / "registration" / "amendment_A1.json"
DEFAULT_AMENDMENT_A2 = REPO / "registration" / "amendment_A2.json"
A2_ID = "teddy_mm_v3_A2"


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


def _amendment_key(amend: dict[str, Any]) -> str:
    """'amendment_A1' for teddy_mm_v3_A1, 'amendment_A2' for teddy_mm_v3_A2 (the key builders read)."""
    return "amendment_" + str(amend.get("amendment_id", "teddy_mm_v3_A1")).rsplit("_", 1)[-1]


def apply_overrides(reg: dict[str, Any], amend: dict[str, Any]) -> dict[str, Any]:
    """Return a deep copy of the registration with the amendment's experiment overrides applied."""
    out = copy.deepcopy(reg)
    for path, value in amend["experiments_overrides"].items():
        _set_path(out["experiments"], path, value)
    out.setdefault("amendments", []).append({k: amend[k] for k in ("amendment_id", "amends", "written_before")})
    out[_amendment_key(amend)] = amend
    return out


def load_amendment_A2(path: Path | str | None = None, *, verify_hash: bool = True,
                      registration_path: Path | str | None = None,
                      amendment_A1_path: Path | str | None = None) -> dict[str, Any]:
    """Load amendment A2; check its own sha256 and that it amends the registration and A1 on disk."""
    path = Path(path) if path is not None else DEFAULT_AMENDMENT_A2
    if not path.exists():
        raise FileNotFoundError(f"amendment A2 missing: {path} (v3 builders apply A1 then A2)")
    if verify_hash:
        hf = path.with_name(path.name + ".sha256")
        if not hf.exists():
            raise FileNotFoundError(f"amendment hash file missing: {hf}")
        recorded = hf.read_text().split()[0].strip()
        actual = vk.sha256_file(path)
        if recorded != actual:
            raise ValueError(f"amendment {path} sha256 {actual} != recorded {recorded}")
    amend = json.loads(path.read_text())
    if amend.get("amendment_id") != A2_ID:
        raise ValueError(f"{path} is not amendment A2 (amendment_id {amend.get('amendment_id')!r})")
    if verify_hash:
        rp = Path(registration_path) if registration_path is not None else vk.DEFAULT_REGISTRATION
        ap = Path(amendment_A1_path) if amendment_A1_path is not None else DEFAULT_AMENDMENT
        reg_sha, a1_sha = vk.sha256_file(rp), vk.sha256_file(ap)
        if amend["amends"]["registration_sha256"] != reg_sha:
            raise ValueError(f"A2 amends registration {amend['amends']['registration_sha256']}, but {rp} has sha256 {reg_sha}")
        if amend["amends"]["amendment_A1_sha256"] != a1_sha:
            raise ValueError(f"A2 follows A1 {amend['amends']['amendment_A1_sha256']}, but {ap} has sha256 {a1_sha}")
    return amend


def load_registration_amended(registration_path: Path | str | None = None,
                              amendment_path: Path | str | None = None,
                              amendment_A2_path: Path | str | None = None) -> dict[str, Any]:
    """The registration every v3 builder uses: hash-checked registration, then hash-checked A1, then
    hash-checked A2 (default: amendment_A2.json next to the A1 file). A missing or altered A2 raises."""
    reg = vk.load_registration(registration_path)
    amend = load_amendment(amendment_path, registration_path=registration_path)
    out = apply_overrides(reg, amend)
    a1p = Path(amendment_path) if amendment_path is not None else DEFAULT_AMENDMENT
    a2p = Path(amendment_A2_path) if amendment_A2_path is not None else a1p.with_name("amendment_A2.json")
    a2 = load_amendment_A2(a2p, registration_path=registration_path, amendment_A1_path=a1p)
    return apply_overrides(out, a2)


def amendment_A2_status(registration_dir: Path | str, git_committed: Callable[[Path], bool | None]) -> dict[str, Any]:
    """Hash and commit state of amendment A2 for a site4 guard (E2 and E3 refuse site4 stages unless
    ``ok``): A2 and its .sha256 file present and committed, the hash file matching, and A2 naming the
    registration and A1 files on disk."""
    rd = Path(registration_dir)
    f, hf = rd / "amendment_A2.json", rd / "amendment_A2.json.sha256"
    info: dict[str, Any] = {"amendment_A2_file": str(f), "amendment_A2_exists": f.exists(),
                            "amendment_A2_sha256": None, "amendment_A2_hash_file_matches": False,
                            "amendment_A2_committed": False, "amendment_A2_amends_files_on_disk": False}
    if f.exists():
        sha = vk.sha256_file(f)
        info["amendment_A2_sha256"] = sha
        info["amendment_A2_hash_file_matches"] = hf.exists() and hf.read_text().split()[0].strip() == sha
        info["amendment_A2_committed"] = (git_committed(f) is True) and hf.exists() and (git_committed(hf) is True)
        try:
            am = json.loads(f.read_text())["amends"]
            info["amendment_A2_amends_files_on_disk"] = (
                am.get("registration_sha256") == vk.sha256_file(rd / "registration_v3.json")
                and am.get("amendment_A1_sha256") == vk.sha256_file(rd / "amendment_A1.json"))
        except (KeyError, ValueError, FileNotFoundError):
            info["amendment_A2_amends_files_on_disk"] = False
    info["ok"] = bool(info["amendment_A2_exists"] and info["amendment_A2_hash_file_matches"]
                      and info["amendment_A2_committed"] and info["amendment_A2_amends_files_on_disk"])
    return info


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


# --------------------------------------------------------------------------- A2.2 E3 R2 normaliser sample
def a2_draw_r2_q95_sample(split: np.ndarray, seed: int, n_cells: int, population: str = "train") -> np.ndarray:
    """The A2.2 draw: numpy default_rng(seed).choice(rows with split == population, n_cells, replace=False),
    sorted. Reads the split column only (label-free; no protein, cell type or embedding)."""
    rows = np.where(np.asarray(split).astype(str) == population)[0].astype(np.int64)
    rng = np.random.default_rng(int(seed))
    return np.sort(rng.choice(rows, int(n_cells), replace=False)).astype(np.int64)


def a2_r2_q95_sample(split: np.ndarray, a2: dict[str, Any]) -> np.ndarray:
    """The registered A2.2 sample of split=train cells on which R2's training q95 is computed.

    Returns the cell list stored in A2 after checking its hash, that every cell is a split=train row and
    that the registered draw (seed, size) reproduces it."""
    spec = a2["e3"]["r2_normaliser"]
    cells = np.asarray(a2["computed"]["e3_r2_normaliser_sample"]["cells"], dtype=np.int64)
    if vk.index_hash(cells) != spec["sample_sha256"] or cells.size != int(spec["n_cells"]):
        raise ValueError("A2.2 sample differs from its recorded hash or size")
    split = np.asarray(split).astype(str)
    if not np.all(split[cells] == spec["population_split"]):
        raise ValueError("A2.2 sample contains a cell outside split == train")
    redraw = a2_draw_r2_q95_sample(split, int(spec["seed"]), int(spec["n_cells"]), spec["population_split"])
    if not np.array_equal(redraw, cells):
        raise ValueError("the registered A2.2 draw does not reproduce the stored sample on this split column")
    return cells


# --------------------------------------------------------------------------- A2.1 E3 H3a
A2_SCALE_ARTEFACT = "not supported (scale artefact)"


def _finite(x: Any) -> bool:
    return x is not None and isinstance(x, (int, float, np.integer, np.floating)) and math.isfinite(float(x))


def a2_positive_condition(point: float | None, lower: float | None, per_donor: Mapping[str, float | None]) -> dict[str, Any]:
    """A2.1 conditions (b) and (c): positive overall (point > 0), two-stage 95% interval excluding 0 on the
    positive side (lower bound > 0), and the registered per-donor rule with margin 0: positive in each
    primary donor. Undefined values fail (conservative)."""
    vals = dict(per_donor)
    ok = (_finite(point) and float(point) > 0 and _finite(lower) and float(lower) > 0
          and len(vals) > 0 and all(_finite(v) and float(v) > 0 for v in vals.values()))
    return {"point": point, "ci95_lower": lower, "per_donor": vals, "holds": bool(ok)}


def a2_h3a_verdict(registered: str, cond_b: bool, cond_c: bool) -> str:
    """A2.1: a readout's H3a verdict. The registered D_R verdict stands unless it is a win; a win stays a win
    only when (b) and (c) hold, otherwise it is 'not supported (scale artefact)'."""
    if registered != "win":
        return registered
    return "win" if (cond_b and cond_c) else A2_SCALE_ARTEFACT


def a2_e3_falsification(verdict_R1: str, verdict_R2: str) -> dict[str, Any]:
    """A2.1 falsification sentence: rejected when neither R1 nor R2 has an A2 win."""
    win = verdict_R1 == "win" or verdict_R2 == "win"
    return {"rejected": not win,
            "text": ("not rejected: R1 or R2 has an A2 win for E3.H3a" if win else
                     "rejected: neither R1 nor R2 has an A2 win for E3.H3a (D_R >= 0.05 over the null with the win "
                     "rules, plus recovery on flagged pairs and a positive scale-free Dlog_R), so 'the NK-T loss is in "
                     "the readout or pooling and is repairable on frozen TEDDY' is rejected for this dataset")}


# --------------------------------------------------------------------------- A2.3 E2 falsification
def _ci_excludes_zero_on_side(diff: float | None, ci: Iterable[float | None]) -> bool:
    lo, hi = list(ci)[:2]
    if not (_finite(diff) and _finite(lo) and _finite(hi)):
        return False
    return (float(diff) > 0 and float(lo) > 0) or (float(diff) < 0 and float(hi) < 0)


def a2_e2_criterion_i(rows: Iterable[Mapping[str, Any]], share_margin: float = 0.15) -> dict[str, Any]:
    """A2.3 for criterion (i): a (pair, label) test counts only under the common decision rule, i.e. it passes
    A1.7 overall (|diff| >= share_margin, interval excluding 0, now on the side of the difference) and in each
    primary donor the difference has the same sign and |diff| >= share_margin / 2 (``per_donor_diff`` of
    e2_response.criterion_i; a donor with an undefined share fails)."""
    out_rows = []
    for r in rows:
        overall = (bool(r.get("passes")) and _finite(r.get("diff")) and abs(float(r["diff"])) >= share_margin
                   and _ci_excludes_zero_on_side(r.get("diff"), r.get("ci95", [None, None])))
        pdd = dict(r.get("per_donor_diff") or {})
        d = r.get("diff")
        donors_ok = (len(pdd) > 0 and _finite(d) and all(
            _finite(v) and np.sign(float(v)) == np.sign(float(d)) and abs(float(v)) >= share_margin / 2
            for v in pdd.values()))
        out_rows.append({"unit_a": r.get("unit_a"), "unit_b": r.get("unit_b"), "label": r.get("label"),
                         "diff": d, "ci95": r.get("ci95"), "per_donor_diff": pdd,
                         "passes_overall": bool(overall), "holds_in_each_primary_donor": bool(donors_ok),
                         "passes_A2": bool(overall and donors_ok)})
    return {"rule": "A2.3: common decision rule binding (overall and in each primary donor)",
            "n_tests": len(out_rows), "n_passing_overall": sum(r["passes_overall"] for r in out_rows),
            "n_passing_A2": sum(r["passes_A2"] for r in out_rows),
            "rows_passing_A2": [r for r in out_rows if r["passes_A2"]],
            "rows_passing_overall_only": [r for r in out_rows if r["passes_overall"] and not r["passes_A2"]],
            "holds": any(r["passes_A2"] for r in out_rows)}


def a2_e2_criterion_ii(diff_row: Mapping[str, Any], margin: float = 0.02) -> dict[str, Any]:
    """A2.3 for criterion (ii): the registered win rule on dAUROC (point >= margin, lower bound > 0, and
    >= margin / 2 in each primary donor)."""
    point = diff_row.get("point")
    lo = (diff_row.get("ci95") or [None, None])[0]
    pdd = dict(diff_row.get("per_donor") or {})
    overall = _finite(point) and float(point) >= margin and _finite(lo) and float(lo) > 0
    donors_ok = len(pdd) > 0 and all(_finite(v) and float(v) >= margin / 2 for v in pdd.values())
    return {"rule": "A2.3: common decision rule binding (overall and in each primary donor)",
            "point": point, "ci95_lower": lo, "per_donor": pdd, "holds_overall": bool(overall),
            "holds_in_each_primary_donor": bool(donors_ok), "holds": bool(overall and donors_ok)}


def a2_e2_verdict(crit_i_holds: bool, crit_ii_holds: bool) -> str:
    """A2.3 E2 verdict: information beyond the curve only if (i) or (ii) holds under the common decision rule."""
    return ("adds information beyond the curve" if (crit_i_holds or crit_ii_holds)
            else "adds nothing beyond the curve (falsified)")
