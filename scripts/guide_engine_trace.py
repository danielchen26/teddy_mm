#!/usr/bin/env python3
"""Worked example for the site guide ("Inside ANM's engine"): one real cell, step by step.

Rebuilds the v0 demo cell (cite_site4_73511) from the event values stored in
outputs/anm_cite_bridge/bakeoff/bakeoff_results.json, runs ANM's finite-field engine
under the soft rule (O0), and records each lineage node's level after every step.
The last step must equal the published action_scores; removing CD16 must flip the call.

Usage (from the teddy_mm repo root, ANM checked out next to it):
    PYTHONPATH=../ANM:. python scripts/guide_engine_trace.py
"""
from __future__ import annotations

import copy
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent / "ANM"))

from active_neural_matter.field.finite_field_runner import build_graph, validate_source_events  # noqa: E402
from bridge_anm.lib.lineage_panels import ACTIONS, LINEAGE_PANELS, MODALITY_FOR_LINEAGE  # noqa: E402
from bridge_anm.run_hard_proof import anm_run_one  # noqa: E402

LINEAGES = ["b_lineage", "t_lineage", "myeloid"]


def trace_levels(schema: dict, instance: dict) -> list[list[float]]:
    """Same loop as finite_field_runner.evolve_field, keeping every step's lineage levels."""
    events = validate_source_events(schema, instance)["field_events"]
    graph = build_graph(instance, events)
    cfg = schema["field_representation"]
    by_time = defaultdict(list)
    for ev in events:
        by_time[int(ev["time"])].append(ev)
    state = {site: 0.0 for site in graph["sites"]}
    out = []
    for tick in range(max(by_time) + int(cfg["steps"]) + 1):
        nxt = {site: cfg["retention"] * v for site, v in state.items()}
        for ev in by_time.get(tick, []):
            nxt[ev["site"]] = nxt.get(ev["site"], 0.0) + cfg["source_scale"] * float(ev["amplitude"])
        prop = dict(nxt)
        for src, neighbours in graph["edges"].items():
            for tgt, w in neighbours:
                prop[tgt] = prop.get(tgt, 0.0) + cfg["diffusion"] * w * nxt.get(src, 0.0)
        state = prop
        out.append([state[f"action:{a}"] for a in LINEAGES])
    return out


def main() -> None:
    res = json.loads((ROOT / "outputs/anm_cite_bridge/bakeoff/bakeoff_results.json").read_text())
    demo = res["anm"]["attribution"]["demo_example"]
    schema = json.loads((ROOT / "outputs/anm_cite_bridge/schema_O0.json").read_text())
    cid = demo["instance_id"]
    value = {e["protein"]: float(e["value"]) for e in demo["top_events"]}
    order = [(lin, p) for lin, prots in LINEAGE_PANELS.items() for p in prots]  # export order = event time
    events = [
        {"event_id": f"{cid}:{p}", "time": t, "action": lin, "modality": MODALITY_FOR_LINEAGE[lin],
         "polarity": "support", "value": value[p], "provenance": "teddy"}
        for t, (lin, p) in enumerate(order)
    ]
    inst = {"instance_id": cid, "question": "lineage coherence from TEDDY δu",
            "actions": copy.deepcopy(ACTIONS), "proposed_source_events": events}

    base = anm_run_one(schema, inst)
    levels = trace_levels(schema, inst)
    for a, v in zip(LINEAGES, levels[-1]):
        assert abs(v - base["action_scores"][a]) < 1e-9, "trace must end at the engine's scores"
        assert abs(v - demo["action_scores"][a]) < 1e-9, "engine must reproduce the published scores"
    assert base["recommended_action"] == demo["recommended_action"]

    without = copy.deepcopy(inst)
    without["proposed_source_events"] = [e for e in events if not e["event_id"].endswith(":CD16")]
    loo = anm_run_one(schema, without)

    print(json.dumps({
        "cell": cid,
        "threshold": schema["field_representation"]["readout_threshold"],
        "field": {k: schema["field_representation"][k] for k in ("steps", "retention", "diffusion")},
        "inputs": [[p, round(value[p], 3)] for _, p in order],
        "levels": [[round(v, 4) for v in row] for row in levels],
        "call": base["recommended_action"],
        "without_CD16": {"scores": {a: round(loo["action_scores"][a], 4) for a in LINEAGES},
                         "call": loo["recommended_action"]},
    }, indent=1))


if __name__ == "__main__":
    main()
