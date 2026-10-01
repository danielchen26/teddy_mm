"""E2 (C3 response decomposition, C5 linearity): the pure logic, kept free of torch and file paths.

Everything here is a declared function of arrays; the numbers it needs (bars, panels, seeds, the
z-probe C) come from the registration, amendment A1 and the E2 addendum, read by
``scripts/e2_response_decomposition.py``. The unit tests in ``tests/test_e2_response.py`` cover
every function below on synthetic inputs.

Vocabulary (registration experiments.E2, A1.7):

* case: one perturbed copy of one cell (a thinning draw, or one gene scaled in one cell);
* unit: the set of cases of one perturbation level (thinning to a kept fraction, both seeds; or
  gene scaling at one epsilon, pooled over genes; the finer gene x epsilon units are secondary);
* lost case: the cell was called and right at baseline and is wrong or not called after the
  perturbation; its label is 'representation' (z-probe wrong on z'), else 'head' (rule argmax wrong
  on v'), else 'decision' (the right class falls below the bar).
"""
from __future__ import annotations

import math
from typing import Any, Iterable, Sequence

import numpy as np

LABELS: tuple[str, ...] = ("representation", "head", "decision")
NO_CALL = ""


# --------------------------------------------------------------------------- counts and perturbations
COUNT_TOL_ABS, COUNT_TOL_REL = 1e-3, 4e-7


def recover_counts(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Integer UMI counts from one stored RNA row.

    The stored values are counts times one per-cell factor (checked: every value divided by the row's
    smallest non-zero value is an integer up to float32 rounding). Returns (counts as int64, the unit =
    the smallest non-zero value, the deviation ratio max_g |q_g - rint(q_g)| / (1e-3 + 4e-7 * rint(q_g)),
    q = value / unit). The tolerance covers float32 rounding of the stored value and of the unit (each
    <= 2^-24 relative). A row is valid when the ratio is <= 1; the caller records invalid rows. An
    all-zero row returns zeros, unit 0, ratio 0."""
    x = np.asarray(values, dtype=np.float64)
    nz = x[x > 0]
    if nz.size == 0:
        return np.zeros(x.shape, dtype=np.int64), 0.0, 0.0
    unit = float(nz.min())
    q = x / unit
    r = np.rint(q)
    ratio = np.abs(q - r) / (COUNT_TOL_ABS + COUNT_TOL_REL * r)
    return r.astype(np.int64), unit, float(np.max(ratio))


def thinning_rng(seed: int, cell: int, level_index: int) -> np.random.Generator:
    """The generator of one thinning draw: numpy default_rng([seed, global cell id, level index])."""
    return np.random.default_rng([int(seed), int(cell), int(level_index)])


def thin_counts(counts: np.ndarray, keep: float, rng: np.random.Generator) -> np.ndarray:
    """Binomial thinning: each UMI is kept independently with probability ``keep``."""
    c = np.asarray(counts, dtype=np.int64)
    if keep >= 1.0:
        return c.copy()
    return rng.binomial(c, float(keep)).astype(np.int64)


def scale_gene(values: np.ndarray, col: int, factor: float) -> np.ndarray:
    """Multiply one gene's official value by ``factor`` (= 1 - epsilon) before ranking; factor 0 masks it."""
    v = np.array(values, dtype=np.float32, copy=True)
    v[..., int(col)] = np.float32(v[..., int(col)] * np.float32(factor))
    return v


# --------------------------------------------------------------------------- scores
def l2n(z: np.ndarray) -> np.ndarray:
    """L2 normalisation as the head's input (z / (||z|| + 1e-6))."""
    z = np.asarray(z, dtype=np.float64)
    return z / (np.linalg.norm(z, axis=-1, keepdims=True) + 1e-6)


def top_margin(S: np.ndarray) -> np.ndarray:
    """S_top1 - S_top2 (TEDDY's margin)."""
    s = np.sort(np.asarray(S, dtype=np.float64), axis=-1)
    return s[..., -1] - s[..., -2]


def signed_margin(S: np.ndarray, c0: np.ndarray) -> np.ndarray:
    """S[c0] - max_{k != c0} S[k]: the margin of a fixed class c0 (negative once another class leads)."""
    S = np.asarray(S, dtype=np.float64)
    c0 = np.asarray(c0, dtype=np.int64)
    rows = np.arange(S.shape[0])
    other = S.copy()
    other[rows, c0] = -np.inf
    return S[rows, c0] - other.max(axis=1)


# --------------------------------------------------------------------------- linearity (C5)
def linearity(d_eps: np.ndarray, d_half: np.ndarray, floor: np.ndarray, *, rho_tol: float = 0.1,
              min_cos: float = 0.9, floor_mult: float = 3.0) -> dict[str, np.ndarray]:
    """Slope ratio rho = ||d(eps)|| / (2 ||d(eps/2)||) and cosine(d(eps), d(eps/2)) per row.

    A row is linear when |rho - 1| <= rho_tol, cosine >= min_cos and ||d(eps/2)|| > floor_mult * floor.
    Rows whose half-step response does not clear the floor are flagged ``below_floor`` (never linear);
    rho and cosine are NaN where a norm is 0 or a value is missing."""
    a = np.asarray(d_eps, dtype=np.float64)
    b = np.asarray(d_half, dtype=np.float64)
    fl = np.asarray(floor, dtype=np.float64)
    na = np.linalg.norm(a, axis=-1)
    nb = np.linalg.norm(b, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(nb > 0, na / (2.0 * nb), np.nan)
        cos = np.where((na > 0) & (nb > 0), np.sum(a * b, axis=-1) / (na * nb), np.nan)
    finite = np.isfinite(rho) & np.isfinite(cos)
    above = np.isfinite(nb) & (nb > floor_mult * fl)
    lin = finite & above & (np.abs(rho - 1.0) <= rho_tol) & (cos >= min_cos)
    return {"rho": rho, "cos": cos, "linear": lin, "below_floor": ~above, "norm_eps": na, "norm_half": nb}


# --------------------------------------------------------------------------- decomposition (C3)
def lost_cases(call_base: np.ndarray, call_pert: np.ndarray, key: np.ndarray) -> np.ndarray:
    """Called and right at baseline, wrong or not called after the perturbation."""
    cb = np.asarray(call_base).astype(str)
    cp = np.asarray(call_pert).astype(str)
    k = np.asarray(key).astype(str)
    return (cb != NO_CALL) & (cb == k) & (cp != k)


def label_lost(probe_pert: np.ndarray, argmax_pert: np.ndarray, key: np.ndarray) -> np.ndarray:
    """Label of a lost case: 'representation' if the z-probe's class on z' is wrong; else 'head' if the
    rule's argmax on v' is wrong; else 'decision' (the right class is the argmax but below the bar).
    A case without an embedding (no tokens left) has probe class '' and is 'representation'."""
    pp = np.asarray(probe_pert).astype(str)
    am = np.asarray(argmax_pert).astype(str)
    k = np.asarray(key).astype(str)
    out = np.full(k.shape, "decision", dtype=object)
    out[am != k] = "head"
    out[pp != k] = "representation"
    return out.astype(str)


# --------------------------------------------------------------------------- bootstrap
def two_stage_weights(donor: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    """Multiplicity weights [n_boot, n_cells] of the registered two-stage bootstrap: donors with
    replacement, then cells with replacement within each drawn donor (as many as the donor has).
    One generator numpy default_rng(seed), donors in sorted order."""
    donor = np.asarray(donor).astype(str)
    ds = sorted(set(donor.tolist()))
    idx = {d: np.where(donor == d)[0] for d in ds}
    rng = np.random.default_rng(int(seed))
    W = np.zeros((int(n_boot), donor.size), dtype=np.int32)
    for b in range(int(n_boot)):
        for d in rng.choice(ds, size=len(ds), replace=True):
            cells = rng.choice(idx[str(d)], size=idx[str(d)].size, replace=True)
            np.add.at(W[b], cells, 1)
    return W


def percentile_ci(x: np.ndarray, level: float = 0.95) -> tuple[float | None, float | None, int]:
    """Percentile interval over the finite replicates; returns (lo, hi, number of non-finite replicates)."""
    x = np.asarray(x, dtype=np.float64)
    ok = np.isfinite(x)
    if not ok.any():
        return None, None, int((~ok).sum())
    a = (1.0 - level) / 2.0
    return float(np.quantile(x[ok], a)), float(np.quantile(x[ok], 1.0 - a)), int((~ok).sum())


def weighted_ratio(num_cell: np.ndarray, den_cell: np.ndarray, W: np.ndarray | None = None):
    """sum_c w_c num_c / sum_c w_c den_c (per replicate when W is [B, n]); NaN where the denominator is 0."""
    num = np.asarray(num_cell, dtype=np.float64)
    den = np.asarray(den_cell, dtype=np.float64)
    if W is None:
        d = den.sum()
        return float(num.sum() / d) if d > 0 else float("nan")
    n = W @ num
    d = W @ den
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(d > 0, n / d, np.nan)


def auroc_groups(score: np.ndarray) -> tuple[np.ndarray, int]:
    """Tie groups of a score (inverse index into the sorted unique values), reusable across replicates."""
    _, inv = np.unique(np.asarray(score, dtype=np.float64), return_inverse=True)
    return inv.ravel(), int(inv.max()) + 1 if inv.size else 0


def auroc(score: np.ndarray, y: np.ndarray, weight: np.ndarray | None = None,
          groups: tuple[np.ndarray, int] | None = None) -> float:
    """Weighted Mann-Whitney AUROC (a positive above a negative counts 1, a tie 1/2); NaN without both
    classes. ``groups`` (from auroc_groups) avoids re-sorting when the same score is used per replicate."""
    y = np.asarray(y).astype(bool)
    w = np.ones(y.size) if weight is None else np.asarray(weight, dtype=np.float64)
    inv, n = groups if groups is not None else auroc_groups(score)
    pw = np.bincount(inv, weights=w * y, minlength=n)
    nw = np.bincount(inv, weights=w * ~y, minlength=n)
    tp, tn = pw.sum(), nw.sum()
    if tp <= 0 or tn <= 0:
        return float("nan")
    below = np.cumsum(nw) - nw
    return float(np.sum(pw * (below + 0.5 * nw)) / (tp * tn))


# --------------------------------------------------------------------------- units and falsification
def unit_cell_tables(case_cell: np.ndarray, n_cells: int, *, key: np.ndarray, call_base: np.ndarray,
                     call_pert: np.ndarray, label: np.ndarray) -> dict[str, np.ndarray]:
    """Per-cell sums over one unit's cases; the bootstrap weights act on cells (a cell carries its cases).

    Arrays are per case: ``case_cell`` = the case's cell position (0..n_cells-1), ``key`` = its cell's
    primary key ('unscored' allowed), ``call_*`` = Q1 calls ('' = no call), ``label`` = label_lost of
    the case (ignored where the case is not lost). Selective accuracy = correct calls / calls on scored
    cases; decision accuracy = (correct call or no call on key OUT) / scored cases."""
    cc = np.asarray(case_cell, dtype=np.int64)
    k = np.asarray(key).astype(str)
    cb = np.asarray(call_base).astype(str)
    cp = np.asarray(call_pert).astype(str)
    lab = np.asarray(label).astype(str)
    scored = k != "unscored"
    called_b, called_p = cb != NO_CALL, cp != NO_CALL
    lost = lost_cases(cb, cp, k)

    def s(x):
        return np.bincount(cc, weights=np.asarray(x, dtype=np.float64), minlength=n_cells)

    out = {
        "n_cases": s(np.ones(cc.size)),
        "n_scored": s(scored),
        "called_all_pert": s(called_p),
        "called_all_base": s(called_b),
        "dec_ok_base": s(scored & ((called_b & (cb == k)) | (~called_b & (k == "OUT")))),
        "dec_ok_pert": s(scored & ((called_p & (cp == k)) | (~called_p & (k == "OUT")))),
        "sel_n_base": s(scored & called_b),
        "sel_n_pert": s(scored & called_p),
        "sel_ok_base": s(scored & called_b & (cb == k)),
        "sel_ok_pert": s(scored & called_p & (cp == k)),
        "call_changed": s(cb != cp),
        "lost": s(lost),
    }
    out["lost_any"] = (out["lost"] > 0).astype(np.float64)
    for L in LABELS:
        out[f"lost_{L}"] = s(lost & (lab == L))
    return out


def unit_stats(T: dict[str, np.ndarray], W: np.ndarray | None = None) -> dict[str, Any]:
    """Accuracies, accuracy losses, coverage and localisation shares of one unit: point values (W None)
    or one value per bootstrap replicate (W = [B, n_cells] multiplicity weights)."""
    r = lambda a, b: weighted_ratio(T[a], T[b], W)  # noqa: E731
    out = {"decision_acc_base": r("dec_ok_base", "n_scored"), "decision_acc_pert": r("dec_ok_pert", "n_scored"),
           "selective_acc_base": r("sel_ok_base", "sel_n_base"), "selective_acc_pert": r("sel_ok_pert", "sel_n_pert"),
           "coverage_base": r("called_all_base", "n_cases"), "coverage_pert": r("called_all_pert", "n_cases"),
           "call_change_rate": r("call_changed", "n_cases")}
    out["decision_loss"] = out["decision_acc_base"] - out["decision_acc_pert"]
    out["selective_loss"] = out["selective_acc_base"] - out["selective_acc_pert"]
    for L in LABELS:
        out[f"share_{L}"] = r(f"lost_{L}", "lost")
    if W is None:
        out["n_cases"] = int(T["n_cases"].sum())
        out["n_scored_cases"] = int(T["n_scored"].sum())
        out["n_lost_cases"] = int(T["lost"].sum())
        out["n_lost_cells"] = int(T["lost_any"].sum())
        for L in LABELS:
            out[f"n_lost_{L}"] = int(T[f"lost_{L}"].sum())
    return out


def criterion_i(units: dict[str, dict[str, Any]], boot: dict[str, dict[str, np.ndarray]], *,
                loss_key: str = "decision_loss", loss_tol: float = 0.02, min_lost: int = 30,
                share_margin: float = 0.15, per_donor: dict[str, dict[str, dict[str, Any]]] | None = None,
                level: float = 0.95) -> dict[str, Any]:
    """A1.7 criterion (i).

    ``units[u]`` holds the point stats of unit u (unit_stats without W) and ``boot[u]`` the per-replicate
    arrays (unit_stats with W, same bootstrap weights for every unit). A pair (u, w) qualifies on the point
    estimates when |loss_u - loss_w| <= loss_tol and both units have >= min_lost lost cells; for each
    qualifying pair and label, the share difference passes when |diff| >= share_margin and its percentile
    interval excludes 0. Replicates with an undefined share (no lost case in a unit) are dropped from that
    interval and counted. Per-donor rows (same sign and |diff| >= share_margin / 2 in each donor) are
    reported beside the verdict."""
    names = sorted(units)
    pairs = []
    for i, u in enumerate(names):
        for w in names[i + 1:]:
            pu, pw = units[u], units[w]
            lu, lw = pu[loss_key], pw[loss_key]
            if not (np.isfinite(lu) and np.isfinite(lw)):
                continue
            if abs(lu - lw) > loss_tol or pu["n_lost_cells"] < min_lost or pw["n_lost_cells"] < min_lost:
                continue
            for L in LABELS:
                d = pu[f"share_{L}"] - pw[f"share_{L}"]
                rep = boot[u][f"share_{L}"] - boot[w][f"share_{L}"]
                lo, hi, nbad = percentile_ci(rep, level)
                excl = lo is not None and (lo > 0 or hi < 0)
                row = {"unit_a": u, "unit_b": w, "label": L, "loss_a": lu, "loss_b": lw,
                       "n_lost_cells_a": pu["n_lost_cells"], "n_lost_cells_b": pw["n_lost_cells"],
                       "share_a": pu[f"share_{L}"], "share_b": pw[f"share_{L}"], "diff": d,
                       "ci95": [lo, hi], "n_undefined_replicates": nbad,
                       "passes": bool(np.isfinite(d) and abs(d) >= share_margin and excl)}
                if per_donor is not None:
                    pd = {}
                    for dn, du in per_donor.items():
                        a, b = du.get(u), du.get(w)
                        pd[dn] = (a[f"share_{L}"] - b[f"share_{L}"]) if (a is not None and b is not None) else None
                    vals = [v for v in pd.values()]
                    same = all(v is not None and np.isfinite(v) and np.sign(v) == np.sign(d)
                               and abs(v) >= share_margin / 2 for v in vals) if vals else False
                    row["per_donor_diff"] = pd
                    row["per_donor_same_sign_half_margin"] = bool(same)
                pairs.append(row)
    n_pairs = len({(r["unit_a"], r["unit_b"]) for r in pairs})
    m = max(1, len(pairs))
    # secondary: Bonferroni over every (pair, label) test (reported, never the verdict)
    for r in pairs:
        u, w, L = r["unit_a"], r["unit_b"], r["label"]
        rep = boot[u][f"share_{L}"] - boot[w][f"share_{L}"]
        lo, hi, _ = percentile_ci(rep, 1.0 - (1.0 - level) / m)
        r["ci_bonferroni"] = [lo, hi]
        r["passes_bonferroni"] = bool(r["passes"] and lo is not None and (lo > 0 or hi < 0))
    return {"n_pairs_compared": n_pairs, "n_tests": len(pairs), "rows": pairs,
            "holds": any(r["passes"] for r in pairs),
            "holds_with_per_donor_clause": any(r["passes"] and r.get("per_donor_same_sign_half_margin") for r in pairs),
            "holds_bonferroni": any(r["passes_bonferroni"] for r in pairs)}


def extrapolate(base: np.ndarray, small: np.ndarray, eps_small: float, eps_target: float) -> np.ndarray:
    """Linear extrapolation base + (eps_target / eps_small) * (small - base)."""
    b = np.asarray(base, dtype=np.float64)
    return b + (float(eps_target) / float(eps_small)) * (np.asarray(small, dtype=np.float64) - b)


def win_rule(point: float, ci: Sequence[float | None], per_donor: dict[str, float | None], margin: float) -> str:
    """Registered decision rule (common.statistics): win / loss / equivalent / inconclusive."""
    lo, hi = ci
    if point is None or not np.isfinite(point) or lo is None:
        return "inconclusive"
    vals = list(per_donor.values())
    if point >= margin and lo > 0 and all(v is not None and v >= margin / 2 for v in vals):
        return "win"
    if point <= -margin and hi < 0 and all(v is not None and v <= -margin / 2 for v in vals):
        return "loss"
    if lo > -margin and hi < margin:
        return "equivalent"
    return "inconclusive"


def anm_schema(base_schema: dict, field: dict, readout_threshold: float) -> dict:
    """The bridge schema with the registered field and the Q1 readout threshold (field_gain(3) * 3 * bar)."""
    import copy

    schema = copy.deepcopy(base_schema)
    schema["field_representation"].update({k: v for k, v in field.items() if k != "kind"})
    schema["field_representation"]["readout_threshold"] = float(readout_threshold)
    schema["source_event_schema"]["allowed_modalities"] = ["evidence"]
    return schema


def anm_engine(v12: np.ndarray, ok: np.ndarray, ffr, schema: dict,
               lineages: Sequence[str] = ("B", "T", "NK", "myeloid")) -> tuple[np.ndarray, np.ndarray]:
    """Run ANM's finite_field_runner on each row: one action per lineage, one support event per panel
    protein (value = evidence; ``v12`` is class-major, 3 per class), every event at t = 0.
    Returns (action scores [n, 4], recommended action or '' per row); rows with ok False are skipped."""
    v = np.asarray(v12, dtype=np.float64)
    n = v.shape[0]
    per = v.shape[1] // len(lineages)
    A = np.full((n, len(lineages)), np.nan)
    calls = np.full(n, NO_CALL, dtype=object)
    actions = [{"id": k, "label": k} for k in lineages]
    for i in range(n):
        if not ok[i]:
            continue
        ev = [{"event_id": f"{k}{jj}", "time": 0, "action": k, "modality": "evidence", "polarity": "support",
               "value": float(v[i, per * ci + jj]), "provenance": "e2"}
              for ci, k in enumerate(lineages) for jj in range(per)]
        inst = {"instance_id": str(i), "question": "Q1", "actions": actions, "proposed_source_events": ev}
        val = ffr.validate_source_events(schema, inst)
        graph = ffr.build_graph(inst, val["field_events"])
        field = ffr.evolve_field(schema, graph, val["field_events"])
        rd = ffr.readout(schema, inst, field["state"])
        A[i] = [rd["action_scores"][k] for k in lineages]
        calls[i] = rd["recommended_action"] or NO_CALL
    return A, calls.astype(str)


def ceil_div(a: int, b: int) -> int:
    return int(math.ceil(a / b))


def flatten_genes(panel_coding_genes: dict[str, Any], nk_t_genes: Iterable[str]) -> list[str]:
    """E2 gene symbols in registration order: panel coding genes (lists flattened), then NK/T genes;
    first occurrence kept."""
    syms: list[str] = []
    for v in panel_coding_genes.values():
        syms += list(v) if isinstance(v, (list, tuple)) else [v]
    syms += list(nk_t_genes)
    return list(dict.fromkeys(map(str, syms)))
