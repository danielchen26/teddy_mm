"""E3 (C4 + repair): NK-T look-alike pairs under the v3 key and readout repair on frozen TEDDY.

Pure functions used by ``scripts/v3_e3_nkt_repair.py``; nothing here reads data or chooses a number.
The registered design is ``registration_v3.json -> experiments.E3`` as amended by A1 (A1.1 tie-break,
A1.8 within-donor pairs); the open details are fixed in ``registration/addenda/E3.json``.

* pairs: k = 10 cosine neighbours on the raw official z within each donor of the evaluated split
  (``knn_cosine``, the function the registration used for the val flag); each unordered edge once;
  an NK-T pair is one key-NK and one key-T cell of the same donor; flagged if cosine >= flag_cosine.
* gap ratio of a readout on a set of pairs: mean over CD56, CD94, CD335, CD3 of
  median |e_R(NK) - e_R(T)| / median |m(NK) - m(T)| (e = readout evidence, m = measured evidence).
* E3.H3a: D_R = [GR_R(fl) - GR_head(fl)] - [GR_R(unfl) - GR_head(unfl)] - D_null.
* E3.H3b: NK-vs-T selective accuracy among key-NK / key-T cells at matched coverage.
* verdicts: registration common.statistics (win / loss / equivalent / inconclusive).
"""
from __future__ import annotations

import math
from typing import Any, Callable, Iterable, Iterator

import numpy as np

from .knn_cosine import knn_cosine

GAP_PROTEINS: tuple[str, ...] = ("CD56", "CD94", "CD335", "CD3")
LINEAGES: tuple[str, ...] = ("B", "T", "NK", "myeloid")


# --------------------------------------------------------------------------- pairs
def knn_edges(z_rows: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Undirected k-NN edges (local indices, lo < hi) with their cosine, built exactly as the registration
    builder did for the val flag: knn_cosine, then one entry per unordered edge (the value of its last
    listing; the two listings carry the same float32 cosine)."""
    nn_idx, nn_sim = knn_cosine(np.asarray(z_rows, dtype=np.float32), int(k))
    n, kk = nn_idx.shape
    a = np.repeat(np.arange(n, dtype=np.int64), kk)
    b = nn_idx.ravel().astype(np.int64)
    s = nn_sim.ravel().astype(np.float64)
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    code = lo * max(n, 1) + hi
    # last listing wins (as dict assignment in the builder): unique on the reversed arrays
    rev = slice(None, None, -1)
    _, first_in_rev = np.unique(code[rev], return_index=True)
    keep = (code.size - 1 - first_in_rev)
    keep = keep[np.argsort(code[keep], kind="stable")]
    E = np.stack([lo[keep], hi[keep]], axis=1)
    return E, s[keep]


def nkt_pairs_within_donors(z: np.ndarray, cells: np.ndarray, donors: np.ndarray,
                            keys: dict[str, np.ndarray], k: int, flag_cosine: float) -> dict[str, Any]:
    """NK-T pairs of the evaluated cells, built within each donor (A1.8).

    z: full [N_all, d] embedding (raw; knn_cosine L2-normalises); cells: global indices of the evaluated
    split; donors: [N_all] donor labels; keys: {variant: [N_all] key labels}. Returns per variant the
    arrays nk, t (global ids), cos, donor, flag, plus the edge counts per donor."""
    cells = np.asarray(cells, dtype=np.int64)
    donors = np.asarray(donors).astype(str)
    out: dict[str, Any] = {"n_edges": {}, "variants": {}}
    acc = {v: {"nk": [], "t": [], "cos": [], "donor": []} for v in keys}
    for d in sorted(set(donors[cells])):
        rows = np.sort(cells[donors[cells] == d])
        E, s = knn_edges(z[rows], k)
        out["n_edges"][d] = int(E.shape[0])
        ga, gb = rows[E[:, 0]], rows[E[:, 1]]
        for v, key in keys.items():
            ka, kb = np.asarray(key)[ga], np.asarray(key)[gb]
            m = ((ka == "NK") & (kb == "T")) | ((ka == "T") & (kb == "NK"))
            nk = np.where(ka[m] == "NK", ga[m], gb[m])
            t = np.where(ka[m] == "NK", gb[m], ga[m])
            acc[v]["nk"].append(nk)
            acc[v]["t"].append(t)
            acc[v]["cos"].append(s[m])
            acc[v]["donor"].append(np.full(nk.size, d, dtype=object))
    for v, A in acc.items():
        nk = np.concatenate(A["nk"]).astype(np.int64) if A["nk"] else np.zeros(0, np.int64)
        t = np.concatenate(A["t"]).astype(np.int64) if A["t"] else np.zeros(0, np.int64)
        cos = np.concatenate(A["cos"]) if A["cos"] else np.zeros(0)
        dn = np.concatenate(A["donor"]).astype(str) if A["donor"] else np.zeros(0, dtype=str)
        out["variants"][v] = {"nk": nk, "t": t, "cos": cos, "donor": dn, "flag": cos >= float(flag_cosine)}
    return out


def pairs_sha256(nk: np.ndarray, t: np.ndarray) -> str:
    """Hash of an NK-T pair set in the form E5 records (sorted lo * 1e6 + hi codes, v3_key.index_hash)."""
    from .v3_key import index_hash

    P = np.sort(np.stack([np.asarray(nk, np.int64), np.asarray(t, np.int64)], axis=1), axis=1)
    return index_hash(sorted(int(x) * 1_000_000 + int(y) for x, y in P))


# --------------------------------------------------------------------------- evidence
def q95(pred_train: np.ndarray) -> np.ndarray:
    """Per-column 95th percentile over training rows (np.percentile, as the registration builder)."""
    return np.percentile(np.asarray(pred_train, dtype=np.float64), 95, axis=0)


def normalise(pred: np.ndarray, scale: np.ndarray, clip: bool = True) -> np.ndarray:
    v = np.asarray(pred, dtype=np.float64) / np.maximum(np.asarray(scale, dtype=np.float64), 1e-12)[None, :]
    return np.clip(v, 0.0, 1.0) if clip else v


# --------------------------------------------------------------------------- gap ratio and D
def gap_ratio(e_nk: np.ndarray, e_t: np.ndarray, m_nk: np.ndarray, m_t: np.ndarray) -> tuple[float, np.ndarray]:
    """GR over a set of pairs: per protein (columns) median |e_nk - e_t| / median |m_nk - m_t|; returns
    (mean over proteins, per-protein ratios). NaN when there are no pairs or a measured median is 0."""
    e_nk, e_t, m_nk, m_t = (np.asarray(x, dtype=np.float64) for x in (e_nk, e_t, m_nk, m_t))
    if e_nk.shape[0] == 0:
        r = np.full(e_nk.shape[1] if e_nk.ndim == 2 else len(GAP_PROTEINS), np.nan)
        return float("nan"), r
    num = np.median(np.abs(e_nk - e_t), axis=0)
    den = np.median(np.abs(m_nk - m_t), axis=0)
    r = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return float(np.mean(r)) if np.all(np.isfinite(r)) else float("nan"), r


def d_statistic(gr: dict[str, dict[str, float]], readout: str, head: str = "head", null: str = "null") -> float:
    """D_R = [GR_R(fl) - GR_head(fl)] - [GR_R(unfl) - GR_head(unfl)] - D_null (registered E3.H3a).

    gr[name] = {"flagged": GR, "unflagged": GR}. Algebraically the head cancels:
    D_R = [GR_R(fl) - GR_null(fl)] - [GR_R(unfl) - GR_null(unfl)]."""
    def raw(r):
        return (gr[r]["flagged"] - gr[head]["flagged"]) - (gr[r]["unflagged"] - gr[head]["unflagged"])
    return float(raw(readout) - raw(null))


def d_log_statistic(gr: dict[str, dict[str, float]], readout: str, head: str = "head", null: str = "null") -> float:
    """Scale-free secondary of D (addendum secondary_H3a, S1): the same contrast on ln GR, so a readout that rescales
    every NK-T difference by a common factor gets 0. NaN if any GR is not positive."""
    vals = [gr[r][s] for r in (readout, head, null) for s in ("flagged", "unflagged")]
    if not all(np.isfinite(v) and v > 0 for v in vals):
        return float("nan")
    lg = {r: {s: math.log(gr[r][s]) for s in ("flagged", "unflagged")} for r in (readout, head, null)}
    return d_statistic(lg, readout, head, null)


# --------------------------------------------------------------------------- matched coverage
def n_select(coverage: float, n: int) -> int:
    """ceil(c * N) with the same guard as v3_amend.coverage_select."""
    return int(math.ceil(float(coverage) * int(n) - 1e-9))


def select_top(conf: np.ndarray, rank: np.ndarray, coverage: float) -> np.ndarray:
    """Boolean mask of the ceil(c * N) most confident rows; ties by ``rank`` (A1.1 permutation rank of the
    rows' global cell ids, precomputed). Equals v3_amend.coverage_select for the same ids."""
    conf = np.asarray(conf, dtype=np.float64)
    key = np.where(np.isnan(conf), -np.inf, conf)
    order = np.lexsort((np.asarray(rank), -key))
    mask = np.zeros(conf.size, dtype=bool)
    mask[order[:n_select(coverage, conf.size)]] = True
    return mask


def selective_accuracy_at(calls: np.ndarray, truth: np.ndarray, conf: np.ndarray, rank: np.ndarray,
                          coverage: float) -> float:
    m = select_top(conf, rank, coverage)
    return float(np.mean(np.asarray(calls)[m] == np.asarray(truth)[m])) if m.any() else float("nan")


# --------------------------------------------------------------------------- rule and confidence scores
def class_scores(ev: np.ndarray, names: list[str], panel: dict[str, list[str]]) -> np.ndarray:
    """Q1 rule scores [n, 4] (LINEAGES order): equal-weight mean of the evidence over each class panel."""
    j = {n: i for i, n in enumerate(names)}
    ev = np.asarray(ev, dtype=np.float64)
    return np.stack([ev[:, [j[p] for p in panel[c]]].mean(axis=1) for c in LINEAGES], axis=1)


def argmax_calls(S: np.ndarray) -> np.ndarray:
    return np.asarray(LINEAGES, dtype=object)[np.argmax(np.asarray(S), axis=1)].astype(str)


def top_score(S: np.ndarray) -> np.ndarray:
    return np.asarray(S, dtype=np.float64).max(axis=1)


def margin_score(S: np.ndarray) -> np.ndarray:
    s = np.sort(np.asarray(S, dtype=np.float64), axis=1)
    return s[:, -1] - s[:, -2]


def entropy_score(S: np.ndarray) -> np.ndarray:
    """Registration common.confidence_scores.entropy: 1 - H(q) / log(4), q = S / sum(S) (uniform if sum is 0)."""
    S = np.asarray(S, dtype=np.float64)
    tot = S.sum(axis=1, keepdims=True)
    q = np.where(tot > 0, S / np.where(tot > 0, tot, 1.0), 1.0 / S.shape[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.where(q > 0, q * np.log(q), 0.0).sum(axis=1)
    return 1.0 - h / math.log(S.shape[1])


def knn_label_disagreement(z_eval: np.ndarray, z_ref: np.ndarray, ref_labels: np.ndarray, calls: np.ndarray,
                           k: int = 10, batch: int = 512) -> np.ndarray:
    """Share of each evaluated cell's k nearest reference cells (cosine) whose label differs from its call."""
    def l2(x):
        x = np.asarray(x, dtype=np.float32)
        return x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-8)
    ze, zr = l2(z_eval), l2(z_ref)
    lab = np.asarray(ref_labels).astype(str)
    calls = np.asarray(calls).astype(str)
    out = np.zeros(ze.shape[0], dtype=np.float64)
    for s in range(0, ze.shape[0], batch):
        e = min(s + batch, ze.shape[0])
        S = ze[s:e] @ zr.T
        nb = np.argpartition(-S, kth=k - 1, axis=1)[:, :k]
        out[s:e] = np.mean(lab[nb] != calls[s:e, None], axis=1)
    return out


# --------------------------------------------------------------------------- bootstrap and verdicts
def two_stage_indices(groups: np.ndarray, donor_order: Iterable[str], n_boot: int, seed: int) -> Iterator[np.ndarray]:
    """Two-stage bootstrap (registration common.statistics.bootstrap): draw donors with replacement, then
    units (cells or pairs) with replacement within each drawn donor. Yields unit indices per replicate."""
    groups = np.asarray(groups).astype(str)
    by = [np.where(groups == d)[0] for d in donor_order]
    rng = np.random.default_rng(int(seed))
    for _ in range(int(n_boot)):
        dd = rng.integers(0, len(by), size=len(by))
        parts = [by[j][rng.integers(0, by[j].size, size=by[j].size)] if by[j].size else by[j] for j in dd]
        yield np.concatenate(parts) if parts else np.zeros(0, np.int64)


def percentile_ci(x: np.ndarray) -> tuple[float | None, float | None, int]:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return None, None, 0
    lo, hi = np.percentile(x, [2.5, 97.5])
    return float(lo), float(hi), int(x.size)


def verdict(point: float | None, lo: float | None, hi: float | None, per_donor: dict[str, float | None],
            margin: float) -> str:
    """Registration common.statistics: win = point >= margin AND lower bound > 0 AND, in each primary donor,
    the same sign and >= margin/2; loss = the mirror; equivalent = interval inside (-margin, +margin);
    otherwise inconclusive (also when a quantity is undefined)."""
    vals = list(per_donor.values())
    if point is None or lo is None or hi is None or not np.isfinite(point) or any(v is None or not np.isfinite(v) for v in vals):
        return "inconclusive"
    if point >= margin and lo > 0 and all(v >= margin / 2 for v in vals):
        return "win"
    if point <= -margin and hi < 0 and all(v <= -margin / 2 for v in vals):
        return "loss"
    if lo > -margin and hi < margin:
        return "equivalent"
    return "inconclusive"


def boot_stat(n_units: int, groups: np.ndarray, donor_order: list[str], n_boot: int, seed: int,
              fn: Callable[[np.ndarray], dict[str, float]]) -> dict[str, np.ndarray]:
    """Evaluate fn(unit indices) -> {name: value} on every two-stage replicate (same draws for every name)."""
    acc: dict[str, list[float]] = {}
    for ix in two_stage_indices(groups, donor_order, n_boot, seed):
        for k_, v in fn(ix).items():
            acc.setdefault(k_, []).append(v)
    return {k_: np.asarray(v, dtype=np.float64) for k_, v in acc.items()}
