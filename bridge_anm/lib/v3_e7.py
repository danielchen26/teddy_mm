"""E7: TEDDY analog of ANM's three-step inverse loop (pure functions).

Used by ``scripts/mode_a_inverse_e7.py``; nothing here reads data or chooses a number. The design is declared in
``registration/addenda/E7.json`` (fixed on train/val before development and before any site4 forward).

The consumer is fixed: TEDDY layer 12 -> gene-mean pooling -> L2 normalisation -> the frozen phase-1 head, read on
the registered NK/T panel (CD56, CD94, CD335, CD3). Histories are exact mean-preserving interventions at the
layer-11 cut: a set of m disjoint token pairs, each pair's two states replaced by their average. For a declared
candidate representation phi of the layer-11 state, a *matched pair* is two histories with equal phi; a *witness
pair* has equal gene-mean but different G-token states (a matched patch and its G twin, or the unpatched cell and a
G-only patch; the G tokens' non-G partners differ too). Two controls make support informative: the clamp (every token
moved by the training NK-T gene-mean difference) and the patch-region detectability witness (only the tokens of a
matched patch moved, so that the gene-mean moves by the same vector: a consumer blind to the patched, non-G tokens
cannot pass it). This module

* sizes and draws the patch sets (``pairs_for_fraction``, ``draw_patch_sets``, ``draw_g_only``) and applies them
  (``apply_average``, ``apply_clamp``, ``apply_region_shift``; numpy or torch);
* turns head outputs into differences in declared units (``deltas``, ``score``);
* summarises each cell (``cell_summary``: D_gene_mean, D_gene_mean_plus_G, W, P, R);
* applies the declared selection rule (``dev_selection``);
* computes the registered two-stage donor bootstrap bounds and the verdict (``two_stage_draws``,
  ``within_donor_draws``, ``bounds``, ``endpoint_flags``, ``verdict``).
"""
from __future__ import annotations

import math
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

CANDIDATES = ("gene_mean", "gene_mean_plus_G", "all_token_states")
KIND_MATCHED, KIND_TWIN, KIND_G_ONLY, KIND_CLAMP, KIND_REGION = 0, 1, 2, 3, 4
KIND_NAMES = {KIND_MATCHED: "matched", KIND_TWIN: "G_twin", KIND_G_ONLY: "G_only", KIND_CLAMP: "clamp",
              KIND_REGION: "patch_region"}
SHIFT_KINDS = (KIND_CLAMP, KIND_REGION)  # histories that move the gene-mean by design (controls, never matched pairs)


# --------------------------------------------------------------------------- patches
def pairs_for_fraction(frac: float, n_nong: int, g: int) -> int:
    """Number of disjoint non-G pairs of a matched patch at magnitude frac (share of the cell's non-G tokens that are
    averaged): m = max(1, floor(frac * (n_nonG - g) / 2)), leaving >= g non-G tokens free as partners of the G tokens."""
    return max(1, int(np.floor(float(frac) * (int(n_nong) - int(g)) / 2.0)))


def frac_code(frac: float) -> int:
    """Integer code of a magnitude (1024 x frac; exact for the declared fractions); 0 = the G-only patches."""
    return int(round(1024 * float(frac)))


def draw_patch_sets(rng: np.random.Generator, g_pos: Sequence[int], nong_pos: Sequence[int], m: int,
                    n_matched: int) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Matched patches at m pairs and their G twins, drawn in this order with one generator, for each i:
    matched_i: 2m distinct non-G positions, pairs (first m, next m);
    twin_i   : matched_i plus every G token averaged with its own non-G partner, the partners drawn among the non-G
               tokens not in matched_i (same gene-mean as matched_i and as the unpatched cell; it differs from
               matched_i only at the G tokens and their g non-G partners, so the G states differ and the candidate-2
               state is not matched). No twin (empty list) if the cell has no G token.
    Each patch is an int64 array [pairs, 2] of real-token positions."""
    g = np.asarray(sorted(int(x) for x in g_pos), np.int64)
    ng = np.asarray(sorted(int(x) for x in nong_pos), np.int64)
    if 2 * m + g.size > ng.size:
        raise ValueError(f"{m} pairs and {g.size} G partners need {2 * m + g.size} non-G tokens, the cell has {ng.size}")
    matched, twins = [], []
    for _ in range(int(n_matched)):
        pick = rng.choice(ng, 2 * m, replace=False)
        pm = np.stack([pick[:m], pick[m:]], 1).astype(np.int64)
        matched.append(pm)
        if g.size:
            free = np.setdiff1d(ng, pick)
            part = rng.choice(free, g.size, replace=False)
            twins.append(np.concatenate([pm, np.stack([g, part], 1)], 0).astype(np.int64))
    return matched, twins


def draw_g_only(rng: np.random.Generator, g_pos: Sequence[int], nong_pos: Sequence[int], n: int) -> list[np.ndarray]:
    """G-only patches (magnitude-free witnesses): every G token averaged with its own distinct non-G partner."""
    g = np.asarray(sorted(int(x) for x in g_pos), np.int64)
    ng = np.asarray(sorted(int(x) for x in nong_pos), np.int64)
    if not g.size:
        return []
    return [np.stack([g, rng.choice(ng, g.size, replace=False)], 1).astype(np.int64) for _ in range(int(n))]


def apply_average(H, pairs: np.ndarray, one_sided: bool = False):
    """Copy of H ([L, d], numpy or torch) with each pair (a, b) replaced by its average (a <- b <- (h_a + h_b) / 2).
    one_sided (tests only): a <- (h_a + h_b) / 2 with b unchanged, so the gene-mean moves by sum (h_b - h_a) / (2 L)."""
    pairs = np.asarray(pairs, np.int64).reshape(-1, 2)
    a, b = pairs[:, 0], pairs[:, 1]
    if hasattr(H, "clone"):
        import torch
        out = H.clone()
        ai = torch.as_tensor(a, device=H.device)
        bi = torch.as_tensor(b, device=H.device)
    else:
        out = np.array(H, copy=True)
        ai, bi = a, b
    avg = 0.5 * (H[ai] + H[bi])
    out[ai] = avg
    if not one_sided:
        out[bi] = avg
    return out


def apply_clamp(H, shift):
    """Clamp patch (positive control): every token state moved by the same vector, so the gene-mean moves by exactly
    that vector and the deviations h_t - mean are unchanged. shift: [d] (numpy or torch, as H)."""
    return H + shift[None, :]


def region_positions(pairs: np.ndarray) -> np.ndarray:
    """The distinct token positions a patch touches (its region)."""
    return np.unique(np.asarray(pairs, np.int64).ravel())


def apply_region_shift(H, pairs: np.ndarray, shift):
    """Patch-region detectability witness: only the tokens of the patch's region (2m positions) are moved, each by
    shift * L / |region|, so the gene-mean moves by exactly shift (as the clamp) while every token outside the region
    keeps its state. shift: [d] (numpy or torch, as H)."""
    pos = region_positions(pairs)
    L = H.shape[0]
    if hasattr(H, "clone"):
        import torch
        out = H.clone()
        pi = torch.as_tensor(pos, device=H.device)
    else:
        out = np.array(H, copy=True)
        pi = pos
    out[pi] = H[pi] + (float(L) / float(pos.size)) * shift[None, :]
    return out


def touches(pairs: np.ndarray, positions: Iterable[int]) -> bool:
    s = set(int(x) for x in positions)
    return any(int(x) in s for x in np.asarray(pairs).ravel())


# --------------------------------------------------------------------------- differences
def deltas(O: np.ndarray, Oref: np.ndarray, units: np.ndarray) -> np.ndarray:
    """Panel difference: max over proteins of |O_p - Oref_p| / u_p. O: [n, p]; Oref: [p] (broadcast) or [n, p]."""
    O = np.atleast_2d(np.asarray(O, np.float64))
    R = np.asarray(Oref, np.float64)
    R = R[None, :] if R.ndim == 1 else R
    return np.max(np.abs(O - R) / np.asarray(units, np.float64)[None, :], axis=1)


def score(O: np.ndarray) -> np.ndarray:
    """NK-T score of the panel outputs in declared order (CD56, CD94, CD335, CD3): mean of the 3 NK proteins - CD3
    (E5's head-score lens)."""
    O = np.asarray(O, np.float64)
    return O[..., :3].mean(-1) - O[..., 3]


def score_deltas(O: np.ndarray, Oref: np.ndarray, unit: float) -> np.ndarray:
    return np.abs(score(np.atleast_2d(O)) - score(np.asarray(Oref, np.float64))) / float(unit)


def observable_deltas(O: np.ndarray, Oref: np.ndarray, observable: str, units: Mapping[str, Any]) -> np.ndarray:
    if observable == "panel":
        return deltas(O, Oref, np.asarray(units["panel"], np.float64))
    if observable == "score":
        return score_deltas(O, Oref, float(units["score"]))
    raise ValueError(f"unknown observable {observable!r}")


# --------------------------------------------------------------------------- per cell
def _max(x: np.ndarray) -> float:
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return float(x.max()) if x.size else float("nan")


def _median(x: np.ndarray) -> float:
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


def cell_summary(O: np.ndarray, O0: np.ndarray, kind: np.ndarray, code: np.ndarray, base: np.ndarray,
                 codes: Iterable[int], dfun: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> dict[str, float]:
    """Per-cell statistics over the histories whose magnitude code is in codes (the magnitude-free G-only, clamp and
    patch-region patches, code 0, always count). dfun(O_a, O_b) = observable difference row by row. Pairs:
      (unpatched, matched_i), (unpatched, twin_i), (matched_i, twin_i), (unpatched, G-only_j): equal gene-mean
      (matched pairs of candidate 1); (unpatched, matched_i): equal gene-mean and G states (matched pairs of
      candidate 2); (matched_i, twin_i) and (unpatched, G-only_j): equal gene-mean, different G states (the G tokens
      and their non-G partners differ: the G witness pairs); (unpatched, clamp_s): the positive control (gene-mean
      moved, deviations unchanged); (unpatched, region_s): the patch-region detectability witness (gene-mean moved by
      the clamp's vector through the tokens of a matched patch only).
    D_gene_mean = max over candidate 1's matched pairs; D_gene_mean_plus_G = max over candidate 2's; W = max over the
    G witness pairs (W_context: (matched, twin) only; W_pure: G-only only); P = max over the clamp patches; R = max
    over the patch-region patches. NaN where a family is empty (the all-token-states candidate has no matched pair by
    construction)."""
    O = np.atleast_2d(np.asarray(O, np.float64))
    kind, code, base = np.asarray(kind), np.asarray(code), np.asarray(base)
    sel = np.isin(code, list(codes)) | np.isin(kind, (KIND_G_ONLY,) + SHIFT_KINDS)
    idx = np.where(sel)[0]
    k = kind[idx]
    d0 = dfun(O[idx], np.asarray(O0, np.float64)) if idx.size else np.zeros(0)
    tw = idx[k == KIND_TWIN]
    d_ctx = dfun(O[tw], O[base[tw]]) if tw.size else np.zeros(0)
    d_m = d0[k == KIND_MATCHED]
    d_tw = d0[k == KIND_TWIN]
    d_go = d0[k == KIND_G_ONLY]
    d_cl = d0[k == KIND_CLAMP]
    d_rg = d0[k == KIND_REGION]
    return {"D_gene_mean": _max(np.concatenate([d_m, d_tw, d_ctx, d_go])),
            "D_gene_mean_plus_G": _max(d_m),
            "D_all_token_states": float("nan"),
            "W": _max(np.concatenate([d_ctx, d_go])), "W_context": _max(d_ctx), "W_pure": _max(d_go),
            "P": _max(d_cl), "R": _max(d_rg),
            "n_matched": int(d_m.size), "n_twin": int(d_tw.size), "n_g_only": int(d_go.size), "n_clamp": int(d_cl.size),
            "n_region": int(d_rg.size)}


# --------------------------------------------------------------------------- development
def q95(x) -> float:
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, 95)) if x.size else float("nan")


def finite_median(x) -> float:
    return _median(x)


def dev_selection(D: Mapping[str, np.ndarray], tol: float, order: Sequence[str] = CANDIDATES) -> dict[str, Any]:
    """Declared selection rule on the development roster. D[candidate] = per-cell maximum over that candidate's
    matched pairs (NaN where the cell has none). A candidate is retained when it has a matched pair in at least one
    cell and the 95th percentile over cells of its per-cell maximum (the point value of the statistic the
    confirmation bounds) is <= tol; among retained candidates the minimal one is selected, ties by declared order
    (the order is by dimension: gene-mean 512 < gene-mean + G 512 (1 + |G|) < all states 512 L). None retained ->
    selected None (stop; the record is kept). The literal reading (every matched difference of every cell <= tol) is
    recorded beside it and never decides."""
    rec: dict[str, Any] = {"tol": float(tol), "order": list(order), "candidates": {}}
    selected = None
    for c in order:
        x = np.asarray(D.get(c, np.zeros(0)), np.float64)
        fin = x[np.isfinite(x)]
        has = bool(fin.size)
        q = q95(fin) if has else float("nan")
        mx = float(fin.max()) if has else float("nan")
        retained = bool(has and q <= tol)
        rec["candidates"][c] = {"cells_with_matched_pairs": int(fin.size), "q95_of_cell_max": q, "max": mx,
                                "median": float(np.median(fin)) if has else float("nan"), "retained": retained,
                                "literal_all_differences_le_tol": bool(has and mx <= tol),
                                "reason": ("no matched pair (uninformative)" if not has else
                                           "q95 <= tol" if retained else "q95 > tol")}
        if retained and selected is None:
            selected = c
    rec["selected"] = selected
    lit = [c for c in order if rec["candidates"][c]["literal_all_differences_le_tol"]]
    rec["literal_rule_selected"] = lit[0] if lit else None
    return rec


# --------------------------------------------------------------------------- bootstrap
def two_stage_draws(donor: np.ndarray, stratum: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    """Registered two-stage bootstrap (common.statistics): per replicate, donors with replacement, then cells with
    replacement within each drawn donor and stratum (the cell's class). Returns row-index arrays."""
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


def within_donor_draws(rows: np.ndarray, stratum: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    """Cell bootstrap within one donor (cells with replacement within each stratum), for the per-donor bounds."""
    rows = np.asarray(rows, np.int64)
    st = np.asarray(stratum).astype(str)[rows]
    groups = [rows[st == s_] for s_ in sorted(set(st.tolist()))]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(int(n_boot)):
        parts = [g[rng.integers(0, g.size, g.size)] for g in groups if g.size]
        out.append(np.concatenate(parts) if parts else np.zeros(0, np.int64))
    return out


def pct95(v: Iterable[float]) -> list[float] | None:
    v = np.asarray(list(v), np.float64)
    v = v[np.isfinite(v)]
    return np.percentile(v, [2.5, 97.5]).tolist() if v.size else None


def bounds(stat: Callable[[np.ndarray], float], rows: np.ndarray, donor: np.ndarray, stratum: np.ndarray,
           n_boot: int, seed: int, draws: list[np.ndarray] | None = None) -> dict[str, Any]:
    """Point value and percentile 95% interval of stat(row indices) under the two-stage draws (pooled), and per donor
    under a cell bootstrap within the donor (same B and seed)."""
    rows = np.asarray(rows, np.int64)
    donor = np.asarray(donor).astype(str)
    if draws is None:
        draws = two_stage_draws(donor[rows], np.asarray(stratum).astype(str)[rows], n_boot, seed)
        draws = [rows[d_] for d_ in draws]
    bs = np.asarray([stat(ix) for ix in draws], np.float64)
    out = {"value": stat(rows), "ci": pct95(bs), "n_defined_replicates": int(np.isfinite(bs).sum()), "per_donor": {}}
    for d_ in sorted(set(donor[rows].tolist())):
        r = rows[donor[rows] == d_]
        wd = within_donor_draws(r, stratum, n_boot, seed)
        bd = np.asarray([stat(ix) for ix in wd], np.float64)
        out["per_donor"][d_] = {"value": stat(r), "ci": pct95(bd), "n_cells": int(r.size)}
    return out


def _lb(b: Mapping[str, Any]) -> float:
    return b["ci"][0] if b.get("ci") else float("nan")


def _ub(b: Mapping[str, Any]) -> float:
    return b["ci"][1] if b.get("ci") else float("nan")


def rejects(b_median: Mapping[str, Any], tol: float, per_donor: bool = True) -> bool:
    """Rejection: the 95% lower bound of the median of the per-cell maximum is > tol, pooled and (per_donor) in each
    donor's own cell bootstrap."""
    ok = np.isfinite(_lb(b_median)) and _lb(b_median) > tol
    if per_donor:
        ok = ok and all(np.isfinite(_lb(v)) and _lb(v) > tol for v in b_median["per_donor"].values())
    return bool(ok)


def supports(b_q95: Mapping[str, Any], tol: float, per_donor: bool = True) -> bool:
    """Bounded support: the 95% upper bound of the 95th percentile of the per-cell maximum is <= tol, pooled and
    (per_donor) in each donor's own cell bootstrap."""
    ok = np.isfinite(_ub(b_q95)) and _ub(b_q95) <= tol
    if per_donor:
        ok = ok and all(np.isfinite(_ub(v)) and _ub(v) <= tol for v in b_q95["per_donor"].values())
    return bool(ok)


def detects(b_median: Mapping[str, Any], threshold: float, per_donor: bool = True) -> bool:
    """Positive control / detectability witness / G witness: 95% lower bound of the median > threshold (pooled, each
    donor)."""
    return rejects(b_median, threshold, per_donor)


def candidate_status(b_median: Mapping[str, Any], b_q95: Mapping[str, Any], tol: float, per_donor: bool = True) -> str:
    if rejects(b_median, tol, per_donor):
        return "rejected"
    if supports(b_q95, tol, per_donor):
        return "supported"
    return "unresolved"


# --------------------------------------------------------------------------- verdict
VERDICTS = {
    "NOT_VALIDATED": "an implementation check failed (padding, layer-12 pooling, permutation, mean preservation, "
                     "module vs explicit layer, official layer-11 state): E7 stops, no verdict",
    "UNRESOLVED_PRECISION": "the float64 re-check differs from float32 by more than tol / 10 on the re-check subset: "
                            "unresolved",
    "LOOP_COMPLETE": "the gene-mean is rejected on fresh comparisons by the G witness, and the selected gene-mean + G "
                     "representation has bounded support on fresh comparisons with the consumer unchanged "
                     "(reject -> select from the declared library -> validate)",
    "SUPPORTED_REVISION_NOT_SHOWN_NECESSARY": "the selected gene-mean + G representation has bounded support, but the "
                                              "gene-mean is not rejected (or the G witness is not detected) on fresh "
                                              "comparisons",
    "GENE_MEAN_SUPPORTED": "the selected gene-mean has bounded support on fresh comparisons with a passing positive "
                           "control: no revision needed for this consumer, intervention class and magnitudes",
    "UNINFORMATIVE": "the selected candidate meets the support bound, but the positive control (clamp) or the "
                     "patch-region detectability witness is not detected (the instrument does not resolve a change of "
                     "the retained gene-mean of the declared size, through every token or through the patched tokens "
                     "only): support cannot be declared",
    "SELECTED_REJECTED": "the candidate selected in development is rejected on fresh comparisons",
    "SELECTED_UNRESOLVED": "the candidate selected in development is neither supported nor rejected",
    "REJECTED_NO_REVISION": "no candidate was selected in development; on fresh comparisons both the gene-mean and "
                            "gene-mean + G are rejected",
    "REJECTED_GENE_MEAN_NO_SELECTION": "no candidate was selected in development; on fresh comparisons the gene-mean "
                                       "is rejected",
    "NO_SELECTION_UNRESOLVED": "no candidate was selected in development and the gene-mean is not rejected",
}


def verdict(selected: str | None, status: Mapping[str, str], pc_ok: bool, sep_ok: bool, impl_ok: bool,
            f64_ok: bool, det_ok: bool) -> dict[str, Any]:
    """Registered verdict of the primary confirmation family. status[c] in rejected / supported / unresolved for
    gene_mean and gene_mean_plus_G; pc_ok = positive control (clamp) detected; det_ok = patch-region detectability
    witness detected; sep_ok = G witness detected. pc_ok and det_ok gate support only: a rejection rests on its own
    witness pairs and never needs them."""
    if not impl_ok:
        code = "NOT_VALIDATED"
    elif not f64_ok:
        code = "UNRESOLVED_PRECISION"
    elif selected is None:
        if status.get("gene_mean") == "rejected" and status.get("gene_mean_plus_G") == "rejected":
            code = "REJECTED_NO_REVISION"
        elif status.get("gene_mean") == "rejected":
            code = "REJECTED_GENE_MEAN_NO_SELECTION"
        else:
            code = "NO_SELECTION_UNRESOLVED"
    elif status.get(selected) == "supported" and not (pc_ok and det_ok):
        code = "UNINFORMATIVE"
    elif status.get(selected) == "supported":
        if selected == "gene_mean_plus_G":
            code = "LOOP_COMPLETE" if (status.get("gene_mean") == "rejected" and sep_ok) else \
                "SUPPORTED_REVISION_NOT_SHOWN_NECESSARY"
        else:
            code = "GENE_MEAN_SUPPORTED"
    elif status.get(selected) == "rejected":
        code = "SELECTED_REJECTED"
    else:
        code = "SELECTED_UNRESOLVED"
    return {"code": code, "verdict": VERDICTS[code]}


def finite_or_none(x: float | None) -> float | None:
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else x
