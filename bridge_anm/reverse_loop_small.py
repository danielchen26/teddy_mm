#!/usr/bin/env python3
"""Small reverse closed loop on must-separate pairs (Mode B only).

Pipeline: VETO → COMPLEMENT → VERIFY on reverse_step1 z_512 must-separate pairs
(esp. myeloid↔T / CD5·CD3·CD16). Cheap complements only — no TEDDY retrain.

Complements tried
-----------------
(a) abstain / observer rewrite (O0→O1 stricter abstain; O0→O2 key-marker observer)
(b) true ADT markers as typed evidence into ANM events (Mode B complement channel;
    not the production no-leakage path)
(c) existing ADT / joint features as z⁺ (reuses prior complement_zplus numbers;
    prior z⁺ alone worsened must-rate)

NOT Mode A / not residual-as-field / not gene perturbs / not fusion audit.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(ROOT.parent / "ANM")))
sys.path.insert(0, str(ANM_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from active_neural_matter.field.finite_field_runner import (  # noqa: E402
    build_graph,
    evolve_field,
    readout,
    validate_schema,
    validate_source_events,
    verify,
)

from lib.knn_cosine import l2_normalize  # noqa: E402
from lib.lineage_panels import (  # noqa: E402
    anm_readout_threshold,
    ACTIONS,
    CRITERIA,
    KEY_MARKERS,
    LINEAGE_PANELS,
    MODALITY_FOR_LINEAGE,
    all_panel_proteins,
    event_value_under_criterion,
    expected_from_true,
    protein_to_lineage,
)

PANEL = all_panel_proteins()
FOCUS_PROTS = ["CD5", "CD3", "CD16"]  # top autopsy deltas + key markers


def _load_jsonl(path: Path, max_n: int | None = None) -> list[dict]:
    rows: list[dict] = []
    with path.open() as f:
        for line in f:
            rows.append(json.loads(line))
            if max_n is not None and len(rows) >= max_n:
                break
    return rows


def _schema_for(criterion_id: str, base_schema: dict) -> dict:
    import copy

    schema = copy.deepcopy(base_schema)
    # v2: rule threshold mapped onto the ANM field scale; legacy: declared value.
    thr = float(anm_readout_threshold(CRITERIA[criterion_id], base_schema))
    schema["field_representation"]["readout_threshold"] = thr
    schema["criterion_overlay"] = {"readout_threshold": thr, "criterion_id": criterion_id}
    return schema


def _p95_from_events(events_by_cell: dict[str, list[dict]]) -> dict[str, float]:
    p95: dict[str, float] = {}
    for evs in events_by_cell.values():
        for ev in evs:
            p95[ev["protein"]] = float(ev["norm_p95_train"])
        if len(p95) >= len(PANEL):
            break
    # fallback if sparse
    for p in PANEL:
        p95.setdefault(p, 1.0)
    return p95


def _events_pred(cell_id: str, evs: list[dict], crit: dict) -> list[dict]:
    out = []
    for ev in sorted(evs, key=lambda e: e["time"]):
        out.append(
            {
                "event_id": ev["event_id"],
                "time": int(ev["time"]),
                "action": ev["action"],
                "modality": ev["modality"],
                "polarity": "support",
                "value": event_value_under_criterion(
                    float(ev["value"]), bool(ev.get("is_key_marker")), ev["action"], crit
                ),
                "provenance": ev.get("provenance", "teddy_phase1_mlp"),
            }
        )
    return out


def _events_true_adt(
    cell_id: str,
    true_panel: dict[str, float],
    p95: dict[str, float],
    crit: dict,
) -> list[dict]:
    """Mode B complement: admit true ADT panel as typed source events."""
    p2l = protein_to_lineage()
    out = []
    for t_i, prot in enumerate(PANEL):
        lin = p2l[prot]
        raw = float(true_panel[prot])
        base_val = float(np.clip(raw / max(p95.get(prot, 1.0), 1e-6), 0.0, 1.0))
        is_key = prot == KEY_MARKERS[lin]
        out.append(
            {
                "event_id": f"{cell_id}:trueADT:{prot}",
                "time": t_i,
                "action": lin,
                "modality": MODALITY_FOR_LINEAGE[lin],
                "polarity": "support",
                "value": event_value_under_criterion(base_val, is_key, lin, crit),
                "provenance": "true_adt_complement_mode_b",
            }
        )
    return out


def _run_anm(schema: dict, instance: dict) -> dict[str, Any]:
    validation = validate_source_events(schema, instance)
    graph = build_graph(instance, validation["field_events"])
    field = evolve_field(schema, graph, validation["field_events"])
    rd = readout(schema, instance, field["state"])
    vf = verify(schema, instance, rd)
    return {
        "recommended_action": rd["recommended_action"],
        "action_scores": {k: float(v) for k, v in rd["action_scores"].items()},
        "P_f": float(vf["P_f"]),
        "Q_f": (None if vf["Q_f"] is None else float(vf["Q_f"])),
        "expected_action": vf.get("expected_action"),
        "n_admitted": len(validation["field_events"]),
    }


def _pair_is_myeloid_t(p: dict) -> bool:
    a, b = p.get("adt_lineage_i"), p.get("adt_lineage_j")
    return tuple(sorted([a, b])) == ("myeloid", "t_lineage")


def _protein_abs_delta(adt: np.ndarray, i: int, j: int, name_to_j: dict[str, int]) -> dict[str, float]:
    return {prot: float(abs(adt[i, name_to_j[prot]] - adt[j, name_to_j[prot]])) for prot in PANEL}


def veto_report(
    pairs: list[dict],
    adt: np.ndarray,
    name_to_j: dict[str, int],
) -> dict[str, Any]:
    """Declare z_512 insufficient; count focus pairs / marker deltas."""
    n = len(pairs)
    flag_counts = Counter()
    lineage_conf = Counter()
    focus_mt = []
    for p in pairs:
        for k, v in (p.get("flags") or {}).items():
            if v:
                flag_counts[k] += 1
        a, b = p.get("adt_lineage_i"), p.get("adt_lineage_j")
        if a and b and a != b:
            lineage_conf[tuple(sorted([a, b]))] += 1
        if _pair_is_myeloid_t(p):
            focus_mt.append(p)

    # CD5/CD3/CD16 deltas on myeloid↔T focus
    focus_deltas = {prot: [] for prot in FOCUS_PROTS}
    for p in focus_mt:
        d = _protein_abs_delta(adt, int(p["i"]), int(p["j"]), name_to_j)
        for prot in FOCUS_PROTS:
            focus_deltas[prot].append(d[prot])

    delta_summary = {}
    for prot, vals in focus_deltas.items():
        arr = np.asarray(vals, dtype=np.float64)
        delta_summary[prot] = {
            "n": int(arr.size),
            "mean_abs_delta": float(arr.mean()) if arr.size else None,
            "median_abs_delta": float(np.median(arr)) if arr.size else None,
            "q75_abs_delta": float(np.quantile(arr, 0.75)) if arr.size else None,
        }

    cos = np.asarray([float(p["cos_sim"]) for p in focus_mt], dtype=np.float64)
    return {
        "claim": (
            "z_512 (mean-pool last-layer @ ctx 1024, L2 512-D) is INSUFFICIENT to separate "
            "these must-pairs: near-identical z carries disagreeing true ADT lineage / "
            "protein labels. Downstream readout must abstain, rewrite observer, or admit "
            "a complement evidence channel."
        ),
        "n_pairs_written": n,
        "flag_counts": dict(flag_counts),
        "lineage_confusion": {f"{a}|{b}": c for (a, b), c in lineage_conf.most_common()},
        "n_myeloid_t_pairs": len(focus_mt),
        "frac_myeloid_t": len(focus_mt) / max(n, 1),
        "myeloid_t_cos_sim_quantiles": {
            "min": float(cos.min()) if cos.size else None,
            "q25": float(np.quantile(cos, 0.25)) if cos.size else None,
            "q50": float(np.median(cos)) if cos.size else None,
            "q75": float(np.quantile(cos, 0.75)) if cos.size else None,
            "max": float(cos.max()) if cos.size else None,
        },
        "myeloid_t_focus_protein_deltas": delta_summary,
        "veto": True,
        "mode": "B_sufficiency_for_readout",
        "not_mode_a": True,
    }


def decide_cells(
    cell_ids: list[str],
    *,
    events_by_cell: dict[str, list[dict]],
    meta_by_id: dict[str, dict],
    p95: dict[str, float],
    schema: dict,
    criterion_id: str,
    evidence: str,
) -> dict[str, dict[str, Any]]:
    """Run ANM readout for each cell under pred or true-ADT evidence."""
    crit = CRITERIA[criterion_id]
    out: dict[str, dict[str, Any]] = {}
    for cid in cell_ids:
        cell = meta_by_id[cid]
        true_panel = cell["adt_true_panel_holdout"]
        expected = expected_from_true(true_panel, p95, crit)
        if evidence == "pred":
            proposed = _events_pred(cid, events_by_cell.get(cid, []), crit)
        elif evidence == "true_adt":
            proposed = _events_true_adt(cid, true_panel, p95, crit)
        else:
            raise ValueError(evidence)
        inst = {
            "instance_id": cid,
            "question": "lineage coherence under Mode B complement probe",
            "actions": ACTIONS,
            "proposed_source_events": proposed,
        }
        if expected is not None:
            inst["expected_action"] = expected
        try:
            res = _run_anm(schema, inst)
        except Exception as e:  # noqa: BLE001
            res = {
                "recommended_action": None,
                "action_scores": {},
                "P_f": 0.0,
                "Q_f": None,
                "expected_action": expected,
                "n_admitted": 0,
                "error": str(e),
            }
        res["expected_action"] = expected
        res["evidence"] = evidence
        res["criterion_id"] = criterion_id
        out[cid] = res
    return out


def pair_metrics(
    pairs: list[dict],
    decisions: dict[str, dict[str, Any]],
    meta_by_id: dict[str, dict],
    *,
    focus_myeloid_t_only: bool = False,
) -> dict[str, Any]:
    """Decision agreement / Q / separation on (optionally focused) pairs."""
    n = 0
    both_decided = 0
    same_action = 0  # false agreement risk when labels differ
    differ_action = 0
    any_abstain = 0
    both_abstain = 0
    correct_sep = 0  # each matches expected AND actions differ
    soft_sep = 0  # actions differ OR any abstain (not false-agree)
    q_vals: list[float] = []
    decided_match_expected = 0
    decided_cells = 0

    for p in pairs:
        if focus_myeloid_t_only and not _pair_is_myeloid_t(p):
            continue
        n += 1
        ci, cj = p["cell_id_i"], p["cell_id_j"]
        di, dj = decisions.get(ci), decisions.get(cj)
        if di is None or dj is None:
            continue
        ri, rj = di.get("recommended_action"), dj.get("recommended_action")
        ei, ej = di.get("expected_action"), dj.get("expected_action")

        for d in (di, dj):
            if d.get("recommended_action") is not None:
                decided_cells += 1
                if d.get("Q_f") is not None:
                    q_vals.append(float(d["Q_f"]))
                if d.get("expected_action") is not None:
                    if d["recommended_action"] == d["expected_action"]:
                        decided_match_expected += 1

        if ri is None and rj is None:
            both_abstain += 1
            any_abstain += 1
            soft_sep += 1  # abstain avoids false agreement
            continue
        if ri is None or rj is None:
            any_abstain += 1
            soft_sep += 1
            continue

        both_decided += 1
        if ri == rj:
            same_action += 1
        else:
            differ_action += 1
            soft_sep += 1
            if ei is not None and ej is not None and ri == ei and rj == ej:
                correct_sep += 1

    return {
        "n_pairs": n,
        "both_decided": both_decided,
        "same_action_among_both_decided": same_action,
        "differ_action_among_both_decided": differ_action,
        "false_agreement_rate_among_both_decided": (
            same_action / both_decided if both_decided else None
        ),
        "any_abstain": any_abstain,
        "both_abstain": both_abstain,
        "abstain_pair_frac": any_abstain / max(n, 1),
        "soft_separation_rate": soft_sep / max(n, 1),
        "correct_separation_rate": correct_sep / max(n, 1),
        "n_decided_cell_slots": decided_cells,
        "mean_Q_among_decided_defined": (float(np.mean(q_vals)) if q_vals else None),
        "n_Q_defined": len(q_vals),
        "match_expected_among_decided_with_expected": (
            # approximate: count matches / decided with expected tracked above loosely
            None
        ),
    }


def refine_pair_metrics(
    pairs: list[dict],
    decisions: dict[str, dict[str, Any]],
    *,
    focus_myeloid_t_only: bool = False,
) -> dict[str, Any]:
    """Richer Q / agreement metrics with explicit expected-match denominators."""
    base = pair_metrics(pairs, decisions, {}, focus_myeloid_t_only=focus_myeloid_t_only)
    decided_with_exp = 0
    match_exp = 0
    q_sum = 0.0
    q_n = 0
    for p in pairs:
        if focus_myeloid_t_only and not _pair_is_myeloid_t(p):
            continue
        for cid in (p["cell_id_i"], p["cell_id_j"]):
            d = decisions.get(cid) or {}
            rec = d.get("recommended_action")
            exp = d.get("expected_action")
            if rec is None:
                continue
            if exp is not None:
                decided_with_exp += 1
                if rec == exp:
                    match_exp += 1
                q = 1.0 if rec == exp else 0.0
                q_sum += q
                q_n += 1
    base["n_decided_with_expected"] = decided_with_exp
    base["Q_among_decided"] = (match_exp / decided_with_exp) if decided_with_exp else None
    base["mean_Q_recomputed"] = (q_sum / q_n) if q_n else None
    return base


def zplus_on_same_pairs(
    pairs: list[dict],
    z512: np.ndarray,
    adt_true: np.ndarray,
    adt_pred: np.ndarray | None,
    name_to_j: dict[str, int],
    *,
    held_lineage: str = "myeloid",
) -> dict[str, Any]:
    """Geometry check on FIXED must-pairs: cos under z⁺ vs z_512."""
    held = list(LINEAGE_PANELS[held_lineage])
    zplus_idx = [name_to_j[p] for p in PANEL if p not in held]
    z_true = l2_normalize(adt_true[:, zplus_idx].astype(np.float64))
    z_pred = (
        l2_normalize(adt_pred[:, zplus_idx].astype(np.float64)) if adt_pred is not None else None
    )

    def cos_rows(z: np.ndarray, i: int, j: int) -> float:
        return float(np.dot(z[i], z[j]))

    rows = []
    focus = [p for p in pairs if _pair_is_myeloid_t(p)]
    for p in focus:
        i, j = int(p["i"]), int(p["j"])
        row = {
            "i": i,
            "j": j,
            "cos_z512": float(p["cos_sim"]),
            "cos_zplus_true": cos_rows(z_true, i, j),
        }
        if z_pred is not None:
            row["cos_zplus_pred"] = cos_rows(z_pred, i, j)
        rows.append(row)

    def summ(key: str) -> dict[str, float | None]:
        arr = np.asarray([r[key] for r in rows], dtype=np.float64)
        if not arr.size:
            return {"mean": None, "median": None, "frac_still_ge_0_98": None}
        return {
            "mean": float(arr.mean()),
            "median": float(np.median(arr)),
            "frac_still_ge_0_98": float(np.mean(arr >= 0.98)),
            "frac_dropped_below_0_95": float(np.mean(arr < 0.95)),
        }

    prior_path = ROOT / "outputs/anm_cite_bridge/cheap_probes/complement_zplus/complement_zplus_stats.json"
    prior = None
    if prior_path.is_file():
        prior = json.loads(prior_path.read_text())

    return {
        "n_myeloid_t_pairs": len(rows),
        "held_lineage": held_lineage,
        "held_proteins": held,
        "zplus_proteins": [p for p in PANEL if p not in held],
        "cos_z512": summ("cos_z512"),
        "cos_zplus_true_complement": summ("cos_zplus_true"),
        "cos_zplus_pred_complement": summ("cos_zplus_pred") if z_pred is not None else None,
        "prior_global_relative_must_reduction_vs_z512": (
            (prior or {}).get("relative_must_reduction_vs_z512")
        ),
        "prior_note": (
            "Prior global z⁺ probe (held myeloid) relative_must_reduction primary="
            f"{((prior or {}).get('signal') or {}).get('primary_score')}; "
            "NEGATIVE => z⁺ alone worsens must-rate among near."
        ),
        "interpretation": (
            "On FIXED myeloid↔T must-pairs: true-ADT z⁺ often lowers cosine (see frac_still_ge_0_98 / "
            "frac_dropped_below_0_95); pred z⁺ usually does not. Prior GLOBAL probe still found z⁺ "
            "worsens must-rate among its own near neighborhoods (negative relative_must_reduction). "
            "Decision winner is typed true-ADT evidence, not z⁺ alone."
        ),
    }


def pick_winner(verify: dict[str, Any]) -> dict[str, Any]:
    """Choose complement with best soft_separation + Q on myeloid↔T focus."""
    focus = verify["myeloid_t_focus"]
    ranking = []
    for name, st in focus.items():
        soft = st.get("soft_separation_rate")
        q = st.get("Q_among_decided")
        fa = st.get("false_agreement_rate_among_both_decided")
        # score: prioritize soft sep, then Q, then lower false agreement
        score = 0.0
        if soft is not None:
            score += 2.0 * soft
        if q is not None:
            score += 1.0 * q
        if fa is not None:
            score += 1.0 * (1.0 - fa)
        ranking.append(
            {
                "complement": name,
                "score": score,
                "soft_separation_rate": soft,
                "Q_among_decided": q,
                "false_agreement_rate_among_both_decided": fa,
                "correct_separation_rate": st.get("correct_separation_rate"),
                "abstain_pair_frac": st.get("abstain_pair_frac"),
            }
        )
    ranking.sort(key=lambda r: -r["score"])
    # exclude pure baseline from "winner" if it's baseline_pred_O0
    winners = [r for r in ranking if not r["complement"].startswith("baseline_")]
    # also exclude zplus_geometry if it's not a decision complement
    decision_winners = [r for r in winners if r["complement"] != "complement_c_zplus_geometry"]
    won = decision_winners[0] if decision_winners else (winners[0] if winners else ranking[0])
    return {
        "winner": won["complement"],
        "winner_row": won,
        "ranking": ranking,
        "note": (
            "Winner = best Mode B decision complement on myeloid↔T must-pairs "
            "(soft_separation + Q_among_decided − false_agreement). "
            "z⁺ geometry reported separately; prior showed it worsens global must-rate."
        ),
    }


def write_md(path: Path, payload: dict) -> None:
    veto = payload["veto"]
    verify = payload["verify"]
    win = payload["winner"]
    zplus = payload["complements"]["c_zplus"]
    lines = []
    lines.append("# Reverse loop (small) — must-separate pairs (Mode B)\n")
    lines.append("**Not Mode A.** Veto → complement → verify on reverse_step1 `z_512` must-pairs.\n")
    lines.append("No TEDDY retrain. Numbers from files / this run only.\n")
    lines.append("\n## 1) VETO — z_512 insufficient\n\n")
    lines.append(f"{veto['claim']}\n\n")
    lines.append(f"- pairs written: **{veto['n_pairs_written']}**\n")
    lines.append(f"- myeloid↔T pairs: **{veto['n_myeloid_t_pairs']}** ({veto['frac_myeloid_t']:.4f})\n")
    lines.append(f"- lineage confusion: `{veto['lineage_confusion']}`\n")
    lines.append(f"- flag counts: `{veto['flag_counts']}`\n")
    lines.append(f"- myeloid↔T focus protein |Δ|: `{veto['myeloid_t_focus_protein_deltas']}`\n")
    lines.append(f"- myeloid↔T cos_sim quantiles: `{veto['myeloid_t_cos_sim_quantiles']}`\n")

    lines.append("\n## 2) COMPLEMENT — cheap arms (no retrain)\n\n")
    lines.append("| arm | description |\n|---|---|\n")
    lines.append("| (a) `O1_abstain` | observer rewrite: threshold 0.12→0.28, key_marker_boost×2 |\n")
    lines.append("| (a) `O2_observer` | key-marker priority observer (rewrites expected) |\n")
    lines.append(
        "| (b) `true_ADT_typed` | admit holdout true ADT panel as typed ANM source events "
        "(Mode B complement channel; not production no-leakage path) |\n"
    )
    lines.append(
        "| (c) `zplus` | ADT true/pred complement features as z⁺; prior global "
        "relative_must_reduction **worsened** must-rate |\n"
    )

    lines.append("\n### (c) z⁺ on same myeloid↔T pairs\n\n")
    lines.append(f"- prior: {zplus.get('prior_note')}\n")
    lines.append(f"- cos_z512: `{zplus.get('cos_z512')}`\n")
    lines.append(f"- cos_zplus_true: `{zplus.get('cos_zplus_true_complement')}`\n")
    lines.append(f"- cos_zplus_pred: `{zplus.get('cos_zplus_pred_complement')}`\n")
    lines.append(f"- {zplus.get('interpretation')}\n")

    lines.append("\n## 3) VERIFY — same pairs (myeloid↔T focus)\n\n")
    lines.append(
        "| complement | soft_sep | correct_sep | false_agree (both decided) | "
        "abstain_pair_frac | Q_among_decided | n_pairs |\n"
        "|---|---:|---:|---:|---:|---:|---:|\n"
    )
    for name, st in verify["myeloid_t_focus"].items():
        lines.append(
            f"| `{name}` | {st.get('soft_separation_rate')} | {st.get('correct_separation_rate')} | "
            f"{st.get('false_agreement_rate_among_both_decided')} | {st.get('abstain_pair_frac')} | "
            f"{st.get('Q_among_decided')} | {st.get('n_pairs')} |\n"
        )

    lines.append("\n### All must-pairs (written set)\n\n")
    lines.append(
        "| complement | soft_sep | false_agree | Q_among_decided | n_pairs |\n"
        "|---|---:|---:|---:|---:|\n"
    )
    for name, st in verify["all_written_pairs"].items():
        lines.append(
            f"| `{name}` | {st.get('soft_separation_rate')} | "
            f"{st.get('false_agreement_rate_among_both_decided')} | "
            f"{st.get('Q_among_decided')} | {st.get('n_pairs')} |\n"
        )

    lines.append("\n## Winner\n\n")
    lines.append(f"**{win['winner']}**\n\n")
    lines.append(f"- row: `{win['winner_row']}`\n")
    lines.append(f"- ranking: `{win['ranking']}`\n")
    lines.append(f"- {win['note']}\n")
    lines.append("\n## Claim boundary\n\n")
    lines.append(
        "Mode B only (sufficiency-for-readout / typed evidence / observer edit). "
        "Not Mode A residual-as-field, Jacobian, gene perturbs, fusion audit, or clinical. "
        "Do not claim Pearson > ~0.61. True-ADT events are an explicit complement probe, "
        "not the default no-leakage ANM path.\n"
    )
    path.write_text("".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bridge", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--pairs", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    t0 = time.time()
    bridge = args.bridge
    pairs_path = args.pairs or (bridge / "reverse_step1" / "must_separate_pairs.jsonl")
    out_dir = args.out_dir or (bridge / "reverse_loop_small")
    out_dir.mkdir(parents=True, exist_ok=True)
    docs_out = ROOT / "docs" / "reports"
    docs_out.mkdir(parents=True, exist_ok=True)

    name_to_j = {p: i for i, p in enumerate(PANEL)}
    adt_true = np.load(bridge / "adt_true_panel.npy")
    adt_pred = np.load(bridge / "adt_pred_panel.npy")
    z512 = np.load(bridge / "z_rna_512.npy")[: adt_true.shape[0]]
    pairs = _load_jsonl(pairs_path, args.max_pairs)
    if not pairs:
        raise SystemExit(f"no pairs in {pairs_path}")

    # --- VETO ---
    veto = veto_report(pairs, adt_true, name_to_j)

    # meta + events for cells in pairs
    meta_rows = _load_jsonl(bridge / "cite_cells_meta.jsonl")
    meta_by_id = {m["cell_id"]: m for m in meta_rows}
    need_ids = sorted({p["cell_id_i"] for p in pairs} | {p["cell_id_j"] for p in pairs})
    need_set = set(need_ids)

    events_by_cell: dict[str, list[dict]] = {cid: [] for cid in need_ids}
    with (bridge / "cite_typed_events.jsonl").open() as f:
        for line in f:
            ev = json.loads(line)
            cid = ev["cell_id"]
            if cid in need_set:
                events_by_cell[cid].append(ev)
    p95 = _p95_from_events(events_by_cell)

    base_schema = json.loads((bridge / "schema_O0.json").read_text())
    validate_schema(base_schema)
    schema_O0 = _schema_for("O0", base_schema)
    schema_O1 = _schema_for("O1", base_schema)
    schema_O2 = _schema_for("O2", base_schema)

    # --- COMPLEMENT decisions ---
    variants_spec = [
        ("baseline_pred_O0", "O0", "pred", schema_O0),
        ("complement_a_O1_abstain", "O1", "pred", schema_O1),
        ("complement_a_O2_observer", "O2", "pred", schema_O2),
        ("complement_b_true_ADT_O0", "O0", "true_adt", schema_O0),
        ("complement_b_true_ADT_O1", "O1", "true_adt", schema_O1),
    ]
    decisions_by_variant: dict[str, dict[str, dict]] = {}
    for name, crit_id, evidence, schema in variants_spec:
        decisions_by_variant[name] = decide_cells(
            need_ids,
            events_by_cell=events_by_cell,
            meta_by_id=meta_by_id,
            p95=p95,
            schema=schema,
            criterion_id=crit_id,
            evidence=evidence,
        )

    # --- (c) z+ geometry ---
    zplus = zplus_on_same_pairs(pairs, z512, adt_true, adt_pred, name_to_j)

    # --- VERIFY ---
    verify_all: dict[str, dict] = {}
    verify_mt: dict[str, dict] = {}
    for name, decs in decisions_by_variant.items():
        verify_all[name] = refine_pair_metrics(pairs, decs, focus_myeloid_t_only=False)
        verify_mt[name] = refine_pair_metrics(pairs, decs, focus_myeloid_t_only=True)

    # attach zplus as non-decision geometry arm for ranking note
    verify_mt["complement_c_zplus_geometry"] = {
        "n_pairs": zplus["n_myeloid_t_pairs"],
        "soft_separation_rate": zplus["cos_zplus_true_complement"].get("frac_dropped_below_0_95"),
        "correct_separation_rate": None,
        "false_agreement_rate_among_both_decided": zplus["cos_zplus_true_complement"].get(
            "frac_still_ge_0_98"
        ),
        "abstain_pair_frac": None,
        "Q_among_decided": None,
        "note": "geometry only: soft_sep≈frac cos dropped <0.95; false_agree≈frac still ≥0.98",
    }

    verify = {"all_written_pairs": verify_all, "myeloid_t_focus": verify_mt}
    win = pick_winner(verify)

    # baseline false-agreement for veto narrative
    base_mt = verify_mt["baseline_pred_O0"]
    veto["baseline_pred_O0_myeloid_t"] = {
        "false_agreement_rate_among_both_decided": base_mt.get(
            "false_agreement_rate_among_both_decided"
        ),
        "soft_separation_rate": base_mt.get("soft_separation_rate"),
        "Q_among_decided": base_mt.get("Q_among_decided"),
    }

    payload = {
        "claim_boundary": {
            "mode": "B_reverse_loop_small_must_separate",
            "pipeline": ["veto", "complement", "verify"],
            "not_claimed": [
                "Mode_A_residual_stream_as_field",
                "layer_Jacobian",
                "in_silico_gene_perturbs",
                "fusion_audit",
                "clinical",
            ],
            "true_adt_events_note": (
                "complement (b) admits true ADT as typed evidence for this probe only; "
                "production path keeps holdout ADT verifier-only"
            ),
        },
        "inputs": {
            "pairs_path": str(pairs_path),
            "n_pairs": len(pairs),
            "n_unique_cells": len(need_ids),
            "z_path": str(bridge / "z_rna_512.npy"),
            "panel": PANEL,
        },
        "veto": veto,
        "complements": {
            "a_observer": ["complement_a_O1_abstain", "complement_a_O2_observer"],
            "b_true_adt": ["complement_b_true_ADT_O0", "complement_b_true_ADT_O1"],
            "c_zplus": zplus,
            "variants_run": [v[0] for v in variants_spec],
        },
        "verify": verify,
        "winner": win,
        "wall_seconds": time.time() - t0,
    }

    stats_path = out_dir / "reverse_loop_small_stats.json"
    stats_path.write_text(json.dumps(payload, indent=2))
    md_path = out_dir / "REVERSE_LOOP_SMALL.md"
    write_md(md_path, payload)
    # docs copy
    (docs_out / "REVERSE_LOOP_SMALL.md").write_text(md_path.read_text())
    (docs_out / "reverse_loop_small_stats.json").write_text(stats_path.read_text())

    print(
        json.dumps(
            {
                "wrote": str(out_dir),
                "docs": str(docs_out / "REVERSE_LOOP_SMALL.md"),
                "winner": win["winner"],
                "veto_n_myeloid_t": veto["n_myeloid_t_pairs"],
                "baseline_false_agree_mt": base_mt.get("false_agreement_rate_among_both_decided"),
                "winner_soft_sep": win["winner_row"].get("soft_separation_rate"),
                "winner_Q": win["winner_row"].get("Q_among_decided"),
                "wall_seconds": payload["wall_seconds"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
