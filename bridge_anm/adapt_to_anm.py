#!/usr/bin/env python3
"""Map exported TEDDY CITE events → ANM finite-field instance bundle under a declared criterion."""
from __future__ import annotations

import argparse
import copy
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.lineage_panels import (  # noqa: E402
    anm_readout_threshold,
    ACTIONS,
    CRITERIA,
    event_value_under_criterion,
    expected_from_true,
)


def _load_jsonl(path: Path):
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def build_instances(
    events_path: Path,
    cells_path: Path,
    criterion_id: str,
    max_cells: int | None,
) -> dict:
    crit = CRITERIA[criterion_id]
    by_cell: dict[str, list[dict]] = defaultdict(list)
    for ev in _load_jsonl(events_path):
        by_cell[ev["cell_id"]].append(ev)

    cells = {c["cell_id"]: c for c in _load_jsonl(cells_path)}
    cell_ids = sorted(by_cell.keys())
    if max_cells is not None:
        cell_ids = cell_ids[:max_cells]

    p95: dict[str, float] = {}
    for cid in cell_ids:
        for ev in by_cell[cid]:
            p95[ev["protein"]] = float(ev["norm_p95_train"])

    instances = []
    for cid in cell_ids:
        cell = cells[cid]
        proposed = []
        for ev in sorted(by_cell[cid], key=lambda e: e["time"]):
            value = event_value_under_criterion(
                float(ev["value"]), bool(ev.get("is_key_marker")), ev["action"], crit
            )
            proposed.append(
                {
                    "event_id": ev["event_id"],
                    "time": int(ev["time"]),
                    "action": ev["action"],
                    "modality": ev["modality"],
                    "polarity": "support",
                    "value": value,
                    "provenance": ev["provenance"],
                }
            )
        expected = expected_from_true(cell["adt_true_panel_holdout"], p95, crit)
        inst = {
            "instance_id": cid,
            "question": (
                "From RNA-only TEDDY phase-1 predicted ADT lineage markers, "
                "which major lineage panel is coherent enough to advance?"
            ),
            "actions": copy.deepcopy(ACTIONS),
            "proposed_source_events": proposed,
            "cell_meta": {
                "site": cell["site"],
                "donor": cell["donor"],
                "cell_type": cell["cell_type"],
                "modality_mask": cell["modality_mask"],
                "criterion_id": criterion_id,
                "slice": cell.get("slice", "site4_test"),
            },
            "raw_records": [
                {
                    "record_id": ev["event_id"],
                    "time": ev["time"],
                    "text": (
                        f"TEDDY MLP predicts {ev['protein']}={ev['adt_pred_raw']:.4f} "
                        f"(norm={ev['value']:.3f}, pearson_mlp={ev['phase1_mlp_pearson']:.3f})"
                    ),
                    "hidden_until_endpoint_verification": ["expected_action", "adt_true"],
                }
                for ev in sorted(by_cell[cid], key=lambda e: e["time"])
            ],
        }
        if expected is not None:
            inst["expected_action"] = expected
        instances.append(inst)

    return {
        "run_id": f"cite_anm_bridge_{criterion_id}",
        "source_artifact": str(events_path),
        "criterion_id": criterion_id,
        "criterion": crit,
        "loading_parameters": {
            "family": "cite_phase1_lineage_coherence",
            "purpose": "TEDDY δu → ANM finite-field local response demo",
            "endpoint_hidden_until_verification": True,
            "claim_boundary": "local_response_diagnosis_only",
        },
        "instances": instances,
    }


def schema_with_threshold(base_schema: dict, threshold: float) -> dict:
    schema = copy.deepcopy(base_schema)
    schema["field_representation"]["readout_threshold"] = float(threshold)
    schema["criterion_overlay"] = {"readout_threshold": float(threshold)}
    return schema


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--events", type=Path, default=ROOT / "outputs/anm_cite_bridge/cite_typed_events.jsonl")
    p.add_argument("--cells", type=Path, default=ROOT / "outputs/anm_cite_bridge/cite_cells_meta.jsonl")
    p.add_argument("--schema", type=Path, default=Path(__file__).resolve().parent / "schemas/cite_lineage_finite_field_v0.json")
    p.add_argument("--criterion", choices=sorted(CRITERIA.keys()), default="O0")
    p.add_argument("--max-cells", type=int, default=None, help="None = all cells in export")
    p.add_argument("--out-dir", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    bundle = build_instances(args.events, args.cells, args.criterion, args.max_cells)
    base = json.loads(args.schema.read_text())
    schema = schema_with_threshold(base, anm_readout_threshold(CRITERIA[args.criterion], base))

    inst_path = args.out_dir / f"anm_instances_{args.criterion}.json"
    schema_path = args.out_dir / f"schema_{args.criterion}.json"
    inst_path.write_text(json.dumps(bundle, indent=2))
    schema_path.write_text(json.dumps(schema, indent=2))
    print(
        json.dumps(
            {
                "instances": str(inst_path),
                "schema": str(schema_path),
                "n_instances": len(bundle["instances"]),
                "with_expected": sum(1 for i in bundle["instances"] if "expected_action" in i),
                "criterion": args.criterion,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
