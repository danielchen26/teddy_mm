"""E5-M: observer-conditioned sufficiency of TEDDY's gene-mean state (pure functions).

Used by ``scripts/mode_a_observer_sufficiency.py``; nothing here reads data or chooses a number. The
design is declared in ``registration/addenda/E5M.json`` (fixed on train/val before any site4 forward).

Per cell the script measures, for m declared input pushes u_1..u_m (columns),

    K_z  = [dz_12 / du_j]_j          (d x m; z_12 = gene-mean of TEDDY's layer-12 token states)
    K_O  = [dO / du_j]_j             (p x m; O = an observer of the token states or of z)

and this module turns them into the declared statistics:

* ``factorization_residual``: r = ||K_O - H K_z||_F / ||K_O||_F with H fixed WITHOUT the tested responses
  (from clamps of the candidate state, or calibrated on training cells). This is the informative test the
  ANM SI allows when first-order kernel inclusion is automatic: with K_z of full column rank every K_O
  equals H' K_z for some H' (e.g. K_O pinv(K_z)), so only a map fixed independently of K_O can fail.
* ``same_cell_fit_residual``: that automatic H' (fitted on the cell itself), reported to show the point.
* ``quotient_diagnostic``: rank and participation ratio of K_O (the readout-visible quotient
  C_T R_T; its rank is <= the readout dimension p, so it is a diagnostic, never the state).
* ``near_null_share``: share of ||K_O||_F^2 carried by source combinations that move z_12 by at most
  tol * sigma_max(K_z) (the near-null right-singular space of K_z).
* ``ridge_select``: the frozen H (ridge of an observer's training outputs on training z, lambda by val).
* ``two_stage_draws`` / ``contrast``: the registered two-stage donor bootstrap of median differences
  between cell strata, with the common win / loss / equivalent rule.
* ``h5m_verdict``: the pre-registered decision (see the addendum for the wording).
* ``richardson``: the FD-gate derivative from central differences at two step sizes.
"""
from __future__ import annotations

import math
from typing import Any, Callable, Iterable, Sequence

import numpy as np


# --------------------------------------------------------------------------- per-cell linear algebra
def factorization_residual(K_O: np.ndarray, H: np.ndarray, K_z: np.ndarray) -> float:
    """||K_O - H K_z||_F / ||K_O||_F (NaN if K_O is zero). K_O: [p, m], H: [p, d], K_z: [d, m]."""
    K_O = np.asarray(K_O, np.float64)
    R = K_O - np.asarray(H, np.float64) @ np.asarray(K_z, np.float64)
    den = float(np.linalg.norm(K_O))
    return float(np.linalg.norm(R) / den) if den > 0 else float("nan")


def same_cell_fit_residual(K_O: np.ndarray, K_z: np.ndarray, rcond: float = 1e-10) -> float:
    """Residual of the map fitted on the cell itself, H' = K_O pinv(K_z): ~0 whenever K_z has full column
    rank (the automatic case that makes first-order kernel inclusion NOT_INFORMATIVE)."""
    K_O = np.asarray(K_O, np.float64)
    K_z = np.asarray(K_z, np.float64)
    H = K_O @ np.linalg.pinv(K_z, rcond=rcond)
    return factorization_residual(K_O, H, K_z)


def singular_values(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, np.float64)
    return np.linalg.svd(M, compute_uv=False) if M.size else np.zeros(0)


def quotient_diagnostic(K_O: np.ndarray, rtol: float) -> dict[str, Any]:
    """Readout-visible quotient diagnostic: numerical rank of K_O (singular values > rtol * max) and the
    participation ratio (sum s^2)^2 / sum s^4; both are <= min(p, m)."""
    s = singular_values(K_O)
    if not s.size or s[0] <= 0:
        return {"rank": 0, "participation_ratio": 0.0, "singular_values": s.tolist()}
    s2 = s ** 2
    return {"rank": int(np.sum(s > rtol * s[0])), "participation_ratio": float(s2.sum() ** 2 / np.sum(s2 ** 2)),
            "singular_values": s.tolist()}


def near_null_share(K_z: np.ndarray, K_O: np.ndarray, tol: float) -> dict[str, Any]:
    """Share of ||K_O||_F^2 in the near-null right-singular space of K_z: right singular vectors with
    singular value <= tol * sigma_max (plus the exact null space when m > d). Returns the share, the
    near-null dimension and K_z's numerical rank at tol."""
    K_z = np.asarray(K_z, np.float64)
    K_O = np.asarray(K_O, np.float64)
    m = K_z.shape[1]
    if m == 0:
        return {"share": float("nan"), "near_null_dim": 0, "rank": 0, "condition": None}
    _, s, Vt = np.linalg.svd(K_z, full_matrices=True)
    smax = float(s[0]) if s.size else 0.0
    keep = int(np.sum(s > tol * smax)) if smax > 0 else 0
    N = Vt[keep:].T  # [m, m - keep]
    tot = float(np.sum(K_O ** 2))
    share = float(np.sum((K_O @ N) ** 2) / tot) if tot > 0 else float("nan")
    cond = float(smax / s[min(m, s.size) - 1]) if s.size and s[min(m, s.size) - 1] > 0 else None
    return {"share": share, "near_null_dim": int(m - keep), "rank": keep, "condition": cond}


def anm_factorization_status(analyze: Callable, K_z: np.ndarray, K_O: np.ndarray, H: np.ndarray, *,
                             atol: float, rtol: float) -> dict[str, Any]:
    """ANM analysis.analyze_factorization on (K_z, K_O) with the frozen H entered as clamp data:
    Gamma = I_d (clamps of every z coordinate) and Lambda = H (the observer's calibrated response to them),
    so the ANM estimator Lambda pinv(Gamma) is exactly the frozen H and never uses K_O. Returns the
    statuses only (the full report carries d x d bases)."""
    K_z = np.asarray(K_z, np.float64)
    d = K_z.shape[0]
    rep = analyze(K_z, np.asarray(K_O, np.float64), np.eye(d), np.asarray(H, np.float64), atol=atol, rtol=rtol)
    sk = rep["source_kernel_inclusion"]
    kz = rep["conditioning"]["Kz"]
    return {"status": rep["status"], "source_kernel_status": sk["status"], "source_kernel_reason": sk.get("reason"),
            "Kz_rank": kz["rank"], "Kz_nullity": kz["nullity"], "Kz_reliable": kz["reliable"],
            "Kz_rank_ambiguous": kz["rank_ambiguous"], "factorization_status": rep["factorization"]["status"],
            "factorization_residual_spectral": rep["factorization"]["residual_norm"],
            "factorization_tolerance": rep["factorization"]["tolerance"],
            "clamp_consistency_status": rep["clamp_kernel_inclusion"]["status"],
            "coverage_status": rep["source_coverage"]["status"]}


# --------------------------------------------------------------------------- frozen H (ridge, train only)
def ridge_select(Xtr: np.ndarray, Ytr: np.ndarray, Xva: np.ndarray, Yva: np.ndarray,
                 lambdas: Sequence[float]) -> dict[str, Any]:
    """Ridge of Y on X fitted on the training rows only; lambda (relative to trace(Xc'Xc) / d) chosen by the
    validation MSE. Returns H [p, d] (Y units per X unit), the intercept, the chosen lambda, the absolute
    penalty, the validation MSE per lambda and the train / val R^2 per output at the chosen lambda."""
    Xtr = np.asarray(Xtr, np.float64)
    Ytr = np.asarray(Ytr, np.float64)
    Xva = np.asarray(Xva, np.float64)
    Yva = np.asarray(Yva, np.float64)
    xm, ym = Xtr.mean(0), Ytr.mean(0)
    Xc, Yc = Xtr - xm, Ytr - ym
    G = Xc.T @ Xc
    XtY = Xc.T @ Yc
    d = G.shape[0]
    scale = float(np.trace(G) / d)
    mses, fits = [], []
    for lam in lambdas:
        B = np.linalg.solve(G + lam * scale * np.eye(d), XtY)  # [d, p]
        pred = (Xva - xm) @ B + ym
        mses.append(float(np.mean((pred - Yva) ** 2)))
        fits.append(B)
    k = int(np.argmin(mses))
    B = fits[k]

    def r2(X, Y):
        P = (X - xm) @ B + ym
        ss = np.sum((Y - Y.mean(0)) ** 2, 0)
        return (1 - np.sum((Y - P) ** 2, 0) / np.where(ss > 0, ss, np.nan)).tolist()
    return {"H": B.T.copy(), "intercept": (ym - xm @ B).copy(), "lambda": float(lambdas[k]),
            "alpha_abs": float(lambdas[k] * scale), "val_mse": dict(zip(map(float, lambdas), mses)),
            "train_r2": r2(Xtr, Ytr), "val_r2": r2(Xva, Yva)}


# --------------------------------------------------------------------------- bootstrap and decision rules
def two_stage_draws(donor: np.ndarray, stratum: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    """Registered two-stage bootstrap (common.statistics): per replicate, donors with replacement, then cells
    with replacement within each drawn donor and stratum. Returns row-index arrays."""
    donor = np.asarray(donor).astype(str)
    stratum = np.asarray(stratum).astype(str)
    donors = sorted(set(donor.tolist()))
    strata = sorted(set(stratum.tolist()))
    groups = {d_: [np.where((donor == d_) & (stratum == s_))[0] for s_ in strata] for d_ in donors}
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(int(n_boot)):
        ds = rng.integers(0, len(donors), len(donors))
        parts = [g[rng.integers(0, g.size, g.size)] for x in ds for g in groups[donors[x]] if g.size]
        out.append(np.concatenate(parts) if parts else np.zeros(0, np.int64))
    return out


def _median_diff(x: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    xa, xb = x[a], x[b]
    xa, xb = xa[np.isfinite(xa)], xb[np.isfinite(xb)]
    return float(np.median(xa) - np.median(xb)) if xa.size and xb.size else float("nan")


def pct95(v: Iterable[float]) -> list[float] | None:
    v = np.asarray(list(v), np.float64)
    v = v[np.isfinite(v)]
    return np.percentile(v, [2.5, 97.5]).tolist() if v.size else None


def contrast(stats: dict[str, Callable[[np.ndarray], float]], rows: np.ndarray, donor: np.ndarray,
             draws: list[np.ndarray]) -> dict[str, Any]:
    """Point value, two-stage 95% interval (the same draws for every statistic) and per-donor values of each
    statistic. stats: name -> f(row indices) -> float; rows: all row indices; donor: per-row donor."""
    donor = np.asarray(donor).astype(str)
    out = {}
    for nm, f in stats.items():
        bs = [f(ix) for ix in draws]
        ok = np.isfinite(bs)
        out[nm] = {"value": f(rows), "ci_two_stage": pct95(np.asarray(bs)[ok]), "n_defined_replicates": int(ok.sum()),
                   "per_donor": {d_: f(rows[donor[rows] == d_]) for d_ in sorted(set(donor[rows].tolist()))}}
    return out


def stratum_median_diff(x: np.ndarray, stratum: np.ndarray, a: str = "pair", b: str = "random"
                        ) -> Callable[[np.ndarray], float]:
    """f(ix) = median of x over the ix rows of stratum a minus that over the ix rows of stratum b."""
    x = np.asarray(x, np.float64)
    stratum = np.asarray(stratum).astype(str)

    def f(ix: np.ndarray) -> float:
        ix = np.asarray(ix, np.int64)
        return _median_diff(x, ix[stratum[ix] == a], ix[stratum[ix] == b])
    return f


def common_rule(value: float, ci: Sequence[float] | None, per_donor: dict[str, float], margin: float) -> str:
    """registration common.statistics: 'win' = value >= margin, lower bound > 0, and >= margin / 2 in each
    primary donor; 'loss' = the mirror; 'equivalent' = interval inside (-margin, margin); else 'inconclusive'."""
    if ci is None or not np.isfinite(value):
        return "inconclusive"
    pd_ = [v for v in per_donor.values()]
    if value >= margin and ci[0] > 0 and pd_ and all(np.isfinite(v) and v >= margin / 2 for v in pd_):
        return "win"
    if value <= -margin and ci[1] < 0 and pd_ and all(np.isfinite(v) and v <= -margin / 2 for v in pd_):
        return "loss"
    if -margin < ci[0] and ci[1] < margin:
        return "equivalent"
    return "inconclusive"


def h5m_verdict(c: dict[str, Any], margin: float, jvp_validated: bool, control_ok: bool) -> dict[str, Any]:
    """Pre-registered E5-M decision from the contrast block c (keys as built by the script): delta_tok (pair -
    random), dd_tok_minus_rand (delta_tok - delta_rand), level_tok_pair (median r_tok on pair cells). jvp_validated:
    the FD gate; control_ok: the positive control (O_head, a function of z, has r ~ 0 with the clamp map)."""
    if not jvp_validated:
        return {"verdict": "JVP not validated (registered finite-difference gate failed): E5-M stops", "code": "NOT_VALIDATED"}
    if not control_ok:
        return {"verdict": "positive control failed (O_head residual not ~0): E5-M stops", "code": "NOT_VALIDATED"}
    d = c["delta_tok"]
    rule = common_rule(d["value"], d["ci_two_stage"], d["per_donor"], margin)
    ci = lambda k: c[k]["ci_two_stage"] or [float("nan"), float("nan")]  # noqa: E731
    over_rand = bool(ci("dd_tok_minus_rand")[0] > 0)
    level_ok = bool(ci("level_tok_pair")[1] < margin)
    info = {"delta_tok_rule": rule, "exceeds_rand_null": over_rand, "level_condition": level_ok}
    if rule == "win" and over_rand:
        return {**info, "code": "REJECT_SUFFICIENCY",
                "verdict": "H5M supported: on NK-T look-alike cells a larger share of the token readout's response flows "
                           "through what gene-mean pooling discards than on matched random cells, beyond a random token "
                           "observer (pooling discards information this independent observer uses)"}
    if rule == "equivalent" and level_ok:
        return {**info, "code": "PASS_SUFFICIENT",
                "verdict": "H5M falsified and z_12 suffices for this observer at first order: on look-alike cells less "
                           "than the margin of the token readout's response flows through what pooling discards (the "
                           "NK-T loss is in the head, not in pooling, for this observer)"}
    if rule in ("equivalent", "loss"):
        return {**info, "code": "H5M_FALSIFIED",
                "verdict": f"H5M falsified (Delta_tok {rule}): no look-alike-specific excess; but the level condition "
                           "fails, so z_12 is not shown sufficient for this observer (a share above the margin of its "
                           "response flows through what pooling discards, on pair and random cells alike)"}
    if rule == "win":
        return {**info, "code": "INCONCLUSIVE",
                "verdict": "inconclusive: Delta_tok wins but does not exceed the O_rand null"}
    return {**info, "code": "INCONCLUSIVE", "verdict": "inconclusive (Delta_tok neither win, loss nor equivalent)"}


def richardson(D_h: np.ndarray, D_2h: np.ndarray) -> np.ndarray:
    """Central differences at h and 2h -> (4 D(h) - D(2h)) / 3 (removes the h^2 truncation term)."""
    return (4.0 * np.asarray(D_h, np.float64) - np.asarray(D_2h, np.float64)) / 3.0


def fd_check(J: np.ndarray, D: np.ndarray, rel_tol: float, cos_tol: float) -> tuple[float, float, bool]:
    """Relative error ||D - J|| / ||J|| and cosine of a finite-difference derivative D against the JVP J."""
    J = np.asarray(J, np.float64).ravel()
    D = np.asarray(D, np.float64).ravel()
    jn, dn = float(np.linalg.norm(J)), float(np.linalg.norm(D))
    if jn == 0:
        return (0.0 if dn == 0 else float("inf")), (1.0 if dn == 0 else 0.0), dn == 0
    rel = float(np.linalg.norm(D - J) / jn)
    cos = float(J @ D / (jn * dn)) if dn > 0 else 0.0
    return rel, cos, bool(rel <= rel_tol and cos >= cos_tol)


def log_or_nan(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(x > 0, np.log(x), np.nan)


def finite_median(x) -> float:
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def isfinite_or_none(x: float | None) -> float | None:
    return None if x is None or not math.isfinite(x) else float(x)
