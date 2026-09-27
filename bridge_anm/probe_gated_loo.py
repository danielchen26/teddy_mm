#!/usr/bin/env python3
"""Probe 2 — Gated LOO ±ε / ±ε/2 with G1–G4 counts (Mode B decision sensitivity).

Beside existing closed-form LOO / top-1 flip attribution, nudge the top LOO
protein evidence by ±ε and ±ε/2 and partition cells into G1–G4.

Mode B framing: decision-layer sensitivity to typed TEDDY evidence magnitude.
NOT Mode A residual/Jacobian/gene-perturb response.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANM_CANDIDATES = [
    ROOT.parent / "ANM",
    Path("/Users/tianchichen/Documents/GitHub/ANM"),
    Path("/workspace/ANM-pr11"),
]
for p in ANM_CANDIDATES:
    if (p / "active_neural_matter").is_dir():
        sys.path.insert(0, str(p))
        break
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.lineage_panels import CRITERIA  # noqa: E402
from active_neural_matter.field.finite_field_runner import (  # noqa: E402
    build_graph,
    evolve_field,
    readout,
    validate_source_events,
    verify,
)

# Reuse hard_proof helpers without importing sklearn-heavy module top-level.
import run_hard_proof as hp  # noqa: E402


def classify_gated(flips: dict) -> str:
    """Exclusive G1–G4 partition.

    flips keys: plus_eps, minus_eps, plus_half, minus_half (bool action changed)

    G4: asymmetric at ε or at ε/2 (sign-dependent)
    G3: symmetric-ish high sensitivity — flips under ±ε/2 (any half flip, not asym)
    G2: flips under ±ε but not under ±ε/2
    G1: no flips under any gated nudge
    """
    pe, me = flips["plus_eps"], flips["minus_eps"]
    ph, mh = flips["plus_half"], flips["minus_half"]
    asym_eps = pe != me and (pe or me)
    asym_half = ph != mh and (ph or mh)
    if asym_eps or asym_half:
        return "G4"
    if ph or mh:
        return "G3"
    if pe or me:
        return "G2"
    return "G1"


def nudge_and_run(schema, inst, event_idx: int, value: float):
    events = copy.deepcopy(inst["proposed_source_events"])
    events[event_idx]["value"] = float(np.clip(value, 0.0, 1.0))
    trial = {**inst, "proposed_source_events": events}
    return hp.anm_run_one(schema, trial)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridge", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--schema", type=Path, default=None)
    ap.add_argument("--eps", type=float, default=0.10)
    ap.add_argument("--max-cells", type=int, default=2500, help="cheap subsample of attr cells")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()

    bridge = args.bridge
    out_dir = args.out_dir or (bridge / "cheap_probes" / "gated_loo")
    out_dir.mkdir(parents=True, exist_ok=True)
    schema_path = args.schema or (bridge / "schema_O0.json")
    base_schema = json.loads(schema_path.read_text())

    events_path = bridge / "cite_typed_events.jsonl"
    cells_path = bridge / "cite_cells_meta.jsonl"
    z_export = bridge / "z_rna_export.npy"
    attr_path = bridge / "hard_proof" / "attr_compact.npz"
    if not attr_path.is_file():
        raise SystemExit(f"missing {attr_path} (need existing LOO attribution)")

    attr = np.load(attr_path, allow_pickle=True)
    cell_ids = [str(x) for x in attr["cell_id"]]
    top_proteins = [str(x) for x in attr["top_protein"]]
    loo_flipped = np.asarray(attr["flipped"], dtype=np.float32)
    flip_dist = np.asarray(attr["flip_distance"], dtype=np.float32)

    rng = np.random.default_rng(args.seed)
    n_all = len(cell_ids)
    if args.max_cells and args.max_cells < n_all:
        sel = np.sort(rng.choice(n_all, size=args.max_cells, replace=False))
    else:
        sel = np.arange(n_all)

    want = {cell_ids[i] for i in sel}
    print(f"loading cells/events for {len(want)} attr cells...", flush=True)
    cells = hp.load_cells(cells_path, z_export if z_export.is_file() else None, site4_only=True)
    by_cell = hp.events_by_cell(events_path, cell_ids=want)
    p95 = hp.p95_from_events(events_path, cell_ids=want)
    schema0 = hp.schema_for("O0", base_schema)
    crit0 = CRITERIA["O0"]
    eps = float(args.eps)
    half = eps / 2.0

    gate_counts = Counter()
    per_protein_gate = defaultdict_counter = {}
    from collections import defaultdict

    per_protein_gate = defaultdict(Counter)
    loo_flip_among = []
    rows_out = []
    t0 = time.time()
    kept = 0
    for ti, idx in enumerate(sel):
        cid = cell_ids[idx]
        if cid not in by_cell or cid not in cells:
            continue
        inst = hp.make_instance(cid, cells[cid], by_cell[cid], crit0, p95)
        base = hp.anm_run_one(schema0, inst)
        if base["recommended_action"] is None:
            continue
        # locate top LOO protein event (from prior attribution)
        top_p = top_proteins[idx]
        events = inst["proposed_source_events"]
        ev_i = None
        for i, ev in enumerate(events):
            if str(ev["event_id"]).rsplit(":", 1)[-1] == top_p:
                ev_i = i
                break
        if ev_i is None:
            continue
        v0 = float(events[ev_i]["value"])
        base_a = base["recommended_action"]

        flips = {}
        for name, dv in [
            ("plus_eps", +eps),
            ("minus_eps", -eps),
            ("plus_half", +half),
            ("minus_half", -half),
        ]:
            out = nudge_and_run(schema0, inst, ev_i, v0 + dv)
            flips[name] = out["recommended_action"] != base_a

        gate = classify_gated(flips)
        gate_counts[gate] += 1
        per_protein_gate[top_p][gate] += 1
        loo_flip_among.append(float(loo_flipped[idx]))
        kept += 1
        if kept <= 2000:
            rows_out.append(
                {
                    "cell_id": cid,
                    "top_protein": top_p,
                    "value0": v0,
                    "loo_flipped": bool(loo_flipped[idx] > 0.5),
                    "loo_flip_distance": None
                    if flip_dist[idx] < 0
                    else float(flip_dist[idx]),
                    "gate": gate,
                    **{f"flip_{k}": bool(v) for k, v in flips.items()},
                }
            )
        if (ti + 1) % 200 == 0:
            print(f"  {ti+1}/{len(sel)} kept={kept} gates={dict(gate_counts)}", flush=True)

    elapsed = time.time() - t0
    n = max(kept, 1)
    gate_frac = {g: gate_counts.get(g, 0) / n for g in ("G1", "G2", "G3", "G4")}
    # Sensitivity signal: mass on G2+G3+G4 (any gated flip) and G3+G4 (strong)
    any_flip_frac = 1.0 - gate_frac["G1"]
    strong_frac = gate_frac["G3"] + gate_frac["G4"]
    loo_flip_rate = float(np.mean(loo_flip_among)) if loo_flip_among else None

    results = {
        "claim_boundary": {
            "mode": "B_gated_LOO_decision_sensitivity",
            "not_claimed": [
                "Mode_A_residual_stream_as_field",
                "layer_Jacobian",
                "in_silico_gene_perturbs",
                "representation_response_under_gene_protein_perturb",
            ],
            "gate_definitions": {
                "G1": "no decision flip under ±ε or ±ε/2 on top LOO protein",
                "G2": "flips under ±ε but not ±ε/2 (scale-gated)",
                "G3": "flips under ±ε/2 (high sensitivity)",
                "G4": "sign-asymmetric flip at ε and/or ε/2",
            },
            "note": "Nudges typed evidence values in ANM finite_field; not TEDDY hidden-state perturbs.",
        },
        "params": {
            "eps": eps,
            "eps_half": half,
            "max_cells": args.max_cells,
            "n_attr_total": n_all,
            "n_selected": int(len(sel)),
            "n_kept_non_abstain": kept,
            "seed": args.seed,
            "criterion": "O0",
        },
        "existing_LOO_baseline": {
            "source": str(attr_path),
            "full_n_attr": n_all,
            "full_top1_flip_rate": float(np.mean(loo_flipped)),
            "subsample_loo_flip_rate": loo_flip_rate,
        },
        "gate_counts": {g: int(gate_counts.get(g, 0)) for g in ("G1", "G2", "G3", "G4")},
        "gate_fractions": gate_frac,
        "per_protein_gate_counts": {
            p: dict(c) for p, c in sorted(per_protein_gate.items())
        },
        "signal": {
            "any_gated_flip_frac": any_flip_frac,
            "strong_G3_G4_frac": strong_frac,
            "G4_asym_frac": gate_frac["G4"],
            "primary_metric": "any_gated_flip_frac",
            "delta_vs_loo_flip": (
                None if loo_flip_rate is None else float(any_flip_frac - loo_flip_rate)
            ),
        },
        "runtime_sec": elapsed,
    }

    (out_dir / "gated_loo_stats.json").write_text(json.dumps(results, indent=2))
    with (out_dir / "gated_loo_rows.jsonl").open("w") as f:
        for r in rows_out:
            f.write(json.dumps(r) + "\n")

    md = []
    md.append("# Gated LOO ±ε / ±ε/2 (Mode B decision sensitivity)\n\n")
    md.append("**Not Mode A.** Nudges typed evidence on top LOO protein; counts G1–G4.\n\n")
    md.append(f"- eps={eps}, eps/2={half}, kept={kept}, runtime={elapsed:.1f}s\n")
    md.append(f"- existing LOO top1 flip rate (full): **{results['existing_LOO_baseline']['full_top1_flip_rate']:.4f}**\n")
    md.append(f"- subsample LOO flip rate: **{loo_flip_rate}**\n")
    md.append(f"- any gated flip (1-G1): **{any_flip_frac:.4f}**\n")
    md.append(f"- strong G3+G4: **{strong_frac:.4f}**\n\n")
    md.append("## Gate counts\n")
    for g in ("G1", "G2", "G3", "G4"):
        md.append(f"- **{g}**: {gate_counts.get(g, 0)} ({gate_frac[g]:.4f}) — {results['claim_boundary']['gate_definitions'][g]}\n")
    md.append("\nDo not claim Pearson > ~0.61. Not clinical.\n")
    (out_dir / "GATED_LOO.md").write_text("".join(md))
    print(json.dumps({"wrote": str(out_dir), "signal": results["signal"], "gates": results["gate_counts"]}, indent=2))


if __name__ == "__main__":
    main()
