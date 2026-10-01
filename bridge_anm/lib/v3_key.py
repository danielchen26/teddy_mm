"""v3 answer key, evidence, panels and questions, read from the frozen registration.

Everything a v3 experiment needs to grade or to decide is computed here from
``registration/registration_v3.json`` (the pre-registration). Nothing in this
module chooses a threshold, a panel or a bar: those numbers were fixed on the
training and validation cells by ``bridge_anm/v3_build_registration.py`` and are
only read back here. ``load_registration`` checks the file against its recorded
SHA-256 so a silently edited registration fails loudly.

Vocabulary (see REGISTRATION_v3.md):

* classes: B, T, NK, myeloid, OUT. For OUT the correct action is "no call".
* annotation class: the 45 annotated cell types mapped to the five classes.
* gated class: declared gates on measured protein (CLR), thresholds frozen.
* primary key: cells where the annotation class and the gated class agree;
  disagreeing cells are ``"unscored"`` (reported separately, never dropped from
  coverage denominators).
* evidence: the phase-1 head's predicted protein, divided by the 95th
  percentile of the head's predictions over the training cells, clipped to [0, 1].
  A constant ADT size factor cancels in this ratio, so no array derived from
  measured protein enters the evidence path.
* rule: equal-weight mean of the evidence over each class panel; call the
  argmax if the top score reaches the question's bar, else no call.

The measured protein and the annotation of an evaluated cell are used only by
the key functions (``annotation_class``, ``gate_class``, ``primary_key``,
``q3_key``), i.e. as the verifier.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_REGISTRATION = REPO / "registration" / "registration_v3.json"
DEFAULT_HASH_FILE = REPO / "registration" / "registration_v3.json.sha256"

CLASSES: tuple[str, ...] = ("B", "T", "NK", "myeloid", "OUT")
LINEAGES: tuple[str, ...] = ("B", "T", "NK", "myeloid")
UNSCORED = "unscored"
NO_CALL = ""


# --------------------------------------------------------------------------- registration
def sha256_file(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_registration(path: Path | str | None = None, *, verify_hash: bool = True,
                      hash_file: Path | str | None = None) -> dict[str, Any]:
    """Load the frozen registration; with ``verify_hash`` compare it to its .sha256 file."""
    path = Path(path) if path is not None else DEFAULT_REGISTRATION
    if verify_hash:
        hf = Path(hash_file) if hash_file is not None else path.with_name(path.name + ".sha256")
        if not hf.exists():
            raise FileNotFoundError(f"registration hash file missing: {hf}")
        recorded = hf.read_text().split()[0].strip()
        actual = sha256_file(path)
        if recorded != actual:
            raise ValueError(f"registration {path} sha256 {actual} != recorded {recorded} ({hf})")
    return json.loads(path.read_text())


# --------------------------------------------------------------------------- splits
def split_indices(split: np.ndarray, sites: np.ndarray, donors: np.ndarray,
                  reg: dict[str, Any], name: str) -> np.ndarray:
    """Global row indices of a registered split (metadata only: split, site, donor)."""
    spec = reg["splits"][name]
    split = np.asarray(split).astype(str)
    sites = np.asarray(sites).astype(str)
    donors = np.asarray(donors).astype(str)
    m = np.ones(split.shape[0], dtype=bool)
    if spec.get("split") is not None:
        m &= np.isin(split, list(spec["split"]))
    if spec.get("sites") is not None:
        m &= np.isin(sites, list(spec["sites"]))
    if spec.get("donors") is not None:
        m &= np.isin(donors, list(spec["donors"]))
    idx = np.where(m)[0].astype(np.int64)
    if "n" in spec and int(spec["n"]) != idx.size:
        raise ValueError(f"split {name}: {idx.size} cells, registration says {spec['n']}")
    if "sha256_indices" in spec and index_hash(idx) != spec["sha256_indices"]:
        raise ValueError(f"split {name}: index hash differs from the registration")
    return idx


def index_hash(idx: Iterable[int]) -> str:
    arr = np.asarray(sorted(int(i) for i in idx), dtype=np.int64)
    return hashlib.sha256(arr.tobytes()).hexdigest()


# --------------------------------------------------------------------------- annotation
def annotation_class(cell_types: Iterable[str], reg: dict[str, Any]) -> np.ndarray:
    """Annotated cell type -> class (B / T / NK / myeloid / OUT). Unknown types raise."""
    amap = reg["annotation_map"]
    out = []
    for t in cell_types:
        t = str(t)
        if t not in amap:
            raise KeyError(f"cell type {t!r} has no class in the registration annotation map")
        out.append(amap[t])
    return np.asarray(out, dtype=object).astype(str)


# --------------------------------------------------------------------------- gates
def _high(adt: np.ndarray, adt_names: list[str], thresholds: dict[str, float]) -> dict[str, np.ndarray]:
    name_to_j = {str(n): j for j, n in enumerate(adt_names)}
    return {p: np.asarray(adt[:, name_to_j[p]], dtype=np.float64) > float(t) for p, t in thresholds.items()}


def gate_class_from_high(H: dict[str, np.ndarray], gate: dict[str, Any]) -> np.ndarray:
    """Hierarchical declared gate. ``H[p]`` is the boolean 'protein p is high' array."""
    n = len(next(iter(H.values())))
    b_markers = list(gate["b_markers"])
    nk_markers = list(gate["nk_markers"])
    bpos = np.zeros(n, dtype=bool)
    for p in b_markers:
        bpos |= H[p]
    nk_count = np.zeros(n, dtype=np.int64)
    for p in nk_markers:
        nk_count += H[p].astype(np.int64)
    cd3, cd14, cd33, cd11c, cd71 = H["CD3"], H["CD14"], H["CD33"], H["CD11c"], H["CD71"]
    nk_min = int(gate["nk_min_high"])
    # Applied in this order; a cell takes the first rule it passes. NK is tested before
    # myeloid, so a cell reaching the myeloid rule has already failed the NK rule.
    rules = [
        ("OUT", cd71),                                                    # erythroid: CD71 high
        ("T", cd3 & ~bpos & ~cd14),
        ("B", ~cd3 & bpos & ~cd14),
        ("NK", ~cd3 & ~bpos & ~cd14 & ~cd33 & (nk_count >= nk_min)),
        ("myeloid", ~cd3 & ~bpos & (cd14 | cd11c | cd33)),
    ]
    out = np.full(n, "OUT", dtype=object)  # passes no lineage gate -> OUT
    assigned = np.zeros(n, dtype=bool)
    for lab, m in rules:
        sel = m & ~assigned
        out[sel] = lab
        assigned |= sel
    return out.astype(str)


def gate_class(adt: np.ndarray, adt_names: Iterable[str], reg: dict[str, Any],
               thresholds: dict[str, float] | None = None) -> np.ndarray:
    """Gated class from measured ADT with the frozen (global) thresholds."""
    gate = reg["gate"]
    thr = dict(gate["thresholds"] if thresholds is None else thresholds)
    return gate_class_from_high(_high(np.asarray(adt), [str(x) for x in adt_names], thr), gate)


def gmm_nonzero_threshold(x: np.ndarray, seed: int = 0, max_n: int | None = None) -> float:
    """Declared threshold estimator: 2-component Gaussian mixture on the non-zero values;
    the threshold is the first point between the two means where the upper component's
    posterior reaches 0.5. Label-free."""
    from sklearn.mixture import GaussianMixture

    x = np.asarray(x, dtype=np.float64)
    x = x[x > 0]
    if max_n is not None and x.size > max_n:
        rng = np.random.default_rng(seed)
        x = x[rng.choice(x.size, size=max_n, replace=False)]
    g = GaussianMixture(2, random_state=seed).fit(x[:, None])
    m = g.means_.ravel()
    o = np.argsort(m)
    grid = np.linspace(m[o[0]], m[o[1]], 4001)[:, None]
    post = g.predict_proba(grid)[:, o[1]]
    if (post >= 0.5).any():
        return float(grid[int(np.argmax(post >= 0.5)), 0])
    return float(np.mean(m))


def per_batch_thresholds(adt: np.ndarray, adt_names: Iterable[str], batch: np.ndarray,
                         reg: dict[str, Any]) -> dict[str, dict[str, float]]:
    """Sensitivity key only: the declared estimator applied within each site x donor batch."""
    gate = reg["gate"]
    names = [str(x) for x in adt_names]
    j = {n: i for i, n in enumerate(names)}
    seed = int(gate["threshold_seed"])
    out: dict[str, dict[str, float]] = {}
    for b in sorted(set(map(str, batch))):
        m = np.asarray(batch).astype(str) == b
        out[b] = {p: gmm_nonzero_threshold(np.asarray(adt)[m, j[p]], seed=seed) for p in gate["proteins"]}
    return out


def gate_class_per_batch(adt: np.ndarray, adt_names: Iterable[str], batch: np.ndarray,
                         reg: dict[str, Any]) -> np.ndarray:
    names = [str(x) for x in adt_names]
    thr_b = per_batch_thresholds(adt, names, batch, reg)
    out = np.empty(len(adt), dtype=object)
    batch = np.asarray(batch).astype(str)
    for b, thr in thr_b.items():
        m = batch == b
        out[m] = gate_class(np.asarray(adt)[m], names, reg, thresholds=thr)
    return out.astype(str)


def batch_labels(sites: np.ndarray, donors: np.ndarray) -> np.ndarray:
    return np.asarray([f"{s}|{d}" for s, d in zip(np.asarray(sites).astype(str), np.asarray(donors).astype(str))])


# --------------------------------------------------------------------------- keys
def primary_key(annot: np.ndarray, gated: np.ndarray) -> np.ndarray:
    """Agreement key: the shared class where annotation and gate agree, else 'unscored'."""
    annot = np.asarray(annot).astype(str)
    gated = np.asarray(gated).astype(str)
    return np.where(annot == gated, annot, UNSCORED).astype(object).astype(str)


def q3_key(primary: np.ndarray, cell_types: Iterable[str], adt: np.ndarray, adt_names: Iterable[str],
           reg: dict[str, Any]) -> np.ndarray:
    """Q3 (CD14-anchored question) key, derived from the primary key.

    Primary-key myeloid cells: annotated classical monocyte with measured CD14 high
    -> 'myeloid'; annotated other myeloid (CD16+ monocyte, cDC) with measured CD14 not
    high -> 'OUT' (no call is correct under Q3); any other combination -> 'unscored'.
    B / T / NK / OUT / unscored are unchanged."""
    spec = reg["questions"]["Q3"]["key"]
    names = [str(x) for x in adt_names]
    thr = float(reg["gate"]["thresholds"][spec["anchor_protein"]])
    cd14_high = np.asarray(adt)[:, names.index(spec["anchor_protein"])] > thr
    ct = np.asarray([str(t) for t in cell_types])
    out = np.asarray(primary).astype(object).copy()
    my = out == "myeloid"
    classical = np.isin(ct, spec["classical_types"])
    other = np.isin(ct, spec["other_myeloid_types"])
    out[my] = UNSCORED
    out[my & classical & cd14_high] = "myeloid"
    out[my & other & ~cd14_high] = "OUT"
    return out.astype(str)


def build_keys(cell_types: Iterable[str], adt: np.ndarray, adt_names: Iterable[str],
               reg: dict[str, Any], *, sites: np.ndarray | None = None,
               donors: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """All registered keys for a set of cells (verifier only)."""
    ct = [str(t) for t in cell_types]
    names = [str(x) for x in adt_names]
    annot = annotation_class(ct, reg)
    gated = gate_class(adt, names, reg)
    prim = primary_key(annot, gated)
    keys = {"annotation": annot, "gated": gated, "primary": prim,
            "q3": q3_key(prim, ct, adt, names, reg)}
    gd = np.asarray(ct) == reg["sensitivity_keys"]["no_gdT158"]["excluded_type"]
    keys["primary_no_gdT158"] = np.where(gd, UNSCORED, prim).astype(str)
    if sites is not None and donors is not None:
        gb = gate_class_per_batch(adt, names, batch_labels(sites, donors), reg)
        keys["gated_per_batch"] = gb
        keys["primary_per_batch"] = primary_key(annot, gb)
    return keys


def cohen_kappa(a: np.ndarray, b: np.ndarray, labels: Iterable[str] = CLASSES) -> float:
    """Cohen's kappa of two label arrays over the given label set (both arrays must use it)."""
    a = np.asarray(a).astype(str)
    b = np.asarray(b).astype(str)
    if a.size == 0:
        return float("nan")
    labels = list(labels)
    po = float(np.mean(a == b))
    pe = float(sum(np.mean(a == c) * np.mean(b == c) for c in labels))
    return float((po - pe) / (1.0 - pe)) if pe < 1.0 else float("nan")


def agreement_report(annot: np.ndarray, gated: np.ndarray) -> dict[str, Any]:
    annot = np.asarray(annot).astype(str)
    gated = np.asarray(gated).astype(str)
    rep: dict[str, Any] = {"n": int(annot.size), "agreement": float(np.mean(annot == gated)) if annot.size else None,
                           "kappa5": cohen_kappa(annot, gated), "per_class": {}}
    for c in CLASSES:
        tp = int(np.sum((annot == c) & (gated == c)))
        na, ng = int(np.sum(annot == c)), int(np.sum(gated == c))
        rep["per_class"][c] = {"n_annotation": na, "n_gated": ng, "n_agree": tp,
                               "recall_vs_annotation": (tp / na) if na else None,
                               "precision_vs_annotation": (tp / ng) if ng else None}
    return rep


# --------------------------------------------------------------------------- evidence
def head_predict(z: np.ndarray, ckpt: Path | str, *, device: str = "cpu", batch_size: int = 2048,
                 size_factor: float = 1.0) -> np.ndarray:
    """Phase-1 head (MLP + decoder, no flow matching) on L2-normalised z; constant size factor."""
    import sys

    import torch

    sys.path.insert(0, str(REPO))
    from teddy_mm.models import AdtDecoder, MLP

    z = np.asarray(z, dtype=np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)
    blob = torch.load(str(ckpt), map_location="cpu", weights_only=False)
    n_adt = int(blob["dec"]["log_theta"].shape[0])
    d = z.shape[1]
    mlp = MLP(d, d, hidden=512)
    dec = AdtDecoder(d, n_adt, hidden=512)
    mlp.load_state_dict(blob["mlp"])
    dec.load_state_dict(blob["dec"])
    mlp.eval(); dec.eval()
    mlp.to(device); dec.to(device)
    out = []
    with torch.no_grad():
        for s in range(0, z.shape[0], batch_size):
            zt = torch.from_numpy(z[s:s + batch_size]).to(device)
            sf = torch.full((zt.shape[0],), float(size_factor), device=device)
            out.append(dec(mlp(zt), sf)[0].cpu().numpy())
    return np.concatenate(out, axis=0) if out else np.zeros((0, n_adt), dtype=np.float32)


def evidence(pred: np.ndarray, reg: dict[str, Any], *, clip: bool = True,
             channel: str = "teddy_head") -> np.ndarray:
    """Normalised evidence for all 134 proteins: prediction / training-cell q95 (clipped to [0, 1])."""
    q95 = np.asarray(reg["evidence"][channel]["q95_train_pred"], dtype=np.float64)
    v = np.asarray(pred, dtype=np.float64) / np.maximum(q95, 1e-12)[None, :]
    return np.clip(v, 0.0, 1.0) if clip else v


# --------------------------------------------------------------------------- questions
def question_panel(reg: dict[str, Any], question: str, panel: str | None = None) -> dict[str, list[str]]:
    q = reg["questions"][question]
    name = panel or q["panel"]
    if name == "anchors":
        return {k: [p] for k, p in q["anchors"].items()}
    return {k: list(v) for k, v in reg["panels"][name]["proteins"].items()}


def class_scores(v134: np.ndarray, adt_names: Iterable[str], reg: dict[str, Any], question: str,
                 panel: str | None = None) -> np.ndarray:
    """Rule scores S [n, 4] in LINEAGES order: equal-weight mean of evidence over each class panel."""
    names = [str(x) for x in adt_names]
    j = {n: i for i, n in enumerate(names)}
    pan = question_panel(reg, question, panel)
    v134 = np.asarray(v134, dtype=np.float64)
    return np.stack([v134[:, [j[p] for p in pan[k]]].mean(axis=1) for k in LINEAGES], axis=1)


def rule_calls(S: np.ndarray, bar: float) -> np.ndarray:
    """Call argmax lineage when the top score reaches ``bar``; else NO_CALL ('')."""
    S = np.asarray(S, dtype=np.float64)
    top = S.max(axis=1)
    arg = S.argmax(axis=1)
    calls = np.asarray([LINEAGES[i] for i in arg], dtype=object)
    calls[top < float(bar)] = NO_CALL
    return calls.astype(str)


def question_calls(v134: np.ndarray, adt_names: Iterable[str], reg: dict[str, Any], question: str) -> np.ndarray:
    S = class_scores(v134, adt_names, reg, question)
    return rule_calls(S, reg["questions"][question]["bar"])


def margin(S: np.ndarray) -> np.ndarray:
    """TEDDY's own margin: top-1 minus top-2 rule score."""
    s = np.sort(np.asarray(S, dtype=np.float64), axis=1)
    return s[:, -1] - s[:, -2]


# --------------------------------------------------------------------------- ANM field (closed form)
def field_star(values: np.ndarray, field: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    """Exact replay of ANM ``finite_graph_scalar`` (finite_field_runner.evolve_field) for one
    action site fed by n event sites, every event at t = 0 (simultaneous evidence).

    values: [n_cells, n] event values. Edges as build_graph: event -> action weight 1.0,
    action -> event weight 0.15. Returns (action state [n_cells], event states [n_cells, n])
    after ``steps`` ticks. Lineages are uncoupled, so a 4-action instance is four stars."""
    v = np.asarray(values, dtype=np.float64)
    steps, ret = int(field["steps"]), float(field["retention"])
    dif, src = float(field["diffusion"]), float(field["source_scale"])
    a = np.zeros(v.shape[0])
    e = np.zeros_like(v)
    for tick in range(steps + 1):
        na = ret * a
        ne = ret * e + (src * v if tick == 0 else 0.0)
        a = na + dif * 1.0 * ne.sum(axis=1)
        e = ne + dif * 0.15 * na[:, None]
    return a, e


def field_gain(n_events: int, field: dict[str, Any]) -> float:
    """Action score per unit of summed event value (all events at t = 0, n events)."""
    a, _ = field_star(np.ones((1, int(n_events))), field)
    return float(a[0] / n_events)


def anm_readout_threshold(reg: dict[str, Any], question: str) -> float:
    """ANM readout threshold that calls exactly where the rule does: field_gain(n) * n * bar
    (n = events per class; all classes have the same n), computed from the registered field."""
    n = len(next(iter(question_panel(reg, question).values())))
    return field_gain(n, reg["anm"]["field_representation"]) * n * float(reg["questions"][question]["bar"])


def anm_action_scores(v134: np.ndarray, adt_names: Iterable[str], reg: dict[str, Any], question: str,
                      readout: str = "action") -> np.ndarray:
    """ANM field scores [n, 4] (LINEAGES order) for a question, events at t = 0.

    readout='action': the action-site value (ANM default readout). With equal event weights
    it equals field_gain(n) * n * S (the rule score) exactly.
    readout='closure': the readout_coordinates readout (finite_field_runner.readout) over the
    action's event sites: closure_weight*min + mean_coordinate_weight*mean + direct_weight*action."""
    names = [str(x) for x in adt_names]
    j = {n: i for i, n in enumerate(names)}
    pan = question_panel(reg, question)
    fld = reg["anm"]["field_representation"]
    rc = reg["anm"]["closure_readout"]
    v134 = np.asarray(v134, dtype=np.float64)
    out = []
    for k in LINEAGES:
        vals = v134[:, [j[p] for p in pan[k]]]
        a, e = field_star(vals, fld)
        if readout == "action":
            out.append(a)
        elif readout == "closure":
            out.append(float(rc["closure_weight"]) * e.min(axis=1)
                       + float(rc["mean_coordinate_weight"]) * e.mean(axis=1)
                       + float(rc["direct_action_weight"]) * a)
        else:
            raise ValueError(f"unknown readout {readout!r}")
    return np.stack(out, axis=1)


def selective_accuracy(calls: np.ndarray, key: np.ndarray) -> dict[str, Any]:
    """Accuracy of calls made on scored key cells (a call on an OUT cell is wrong) and
    decision accuracy (no call is correct exactly for OUT)."""
    calls = np.asarray(calls).astype(str)
    key = np.asarray(key).astype(str)
    scored = key != UNSCORED
    called = calls != NO_CALL
    n_called = int(np.sum(scored & called))
    correct_called = int(np.sum(scored & called & (calls == key)))
    decision_ok = (called & (calls == key)) | (~called & (key == "OUT"))
    return {
        "n_cells": int(calls.size), "n_scored": int(scored.sum()), "coverage_all": float(called.mean()) if calls.size else None,
        "n_called_scored": n_called,
        "accuracy_of_calls": (correct_called / n_called) if n_called else None,
        "decision_accuracy": float(decision_ok[scored].mean()) if scored.any() else None,
        "out_decline_rate": float((~called[key == "OUT"]).mean()) if np.any(key == "OUT") else None,
    }
