#!/usr/bin/env python3
"""Best-case lineage scores in Block 2's "RNA missing" condition.

With RNA missing, every evidence value is multiplied by the declared trust factor
(0.35, bridge_anm/run_missing_modality_demo.py) before ANM's engine runs. This script
feeds the engine the best possible evidence (every marker at 1.0) under that factor and
prints the highest score each lineage can reach, per question. If the best case stays
below the question's bar, every cell is declined by construction.

Usage (from the teddy_mm repo root, ANM checked out next to it):
    PYTHONPATH=../ANM:. python scripts/missing_modality_bounds.py
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent / "ANM"))

from bridge_anm.lib.lineage_panels import (  # noqa: E402
    ACTIONS, CRITERIA, KEY_MARKERS, LINEAGE_PANELS, MODALITY_FOR_LINEAGE, event_value_under_criterion,
)
from bridge_anm.run_hard_proof import anm_run_one, schema_for  # noqa: E402
from bridge_anm.run_missing_modality_demo import modality_reliability  # noqa: E402


def best_case(crit_id: str, trust: float) -> dict:
    crit = CRITERIA[crit_id]
    events = []
    for t, (lin, prot) in enumerate((lin, p) for lin, prots in LINEAGE_PANELS.items() for p in prots):
        value = event_value_under_criterion(1.0, prot == KEY_MARKERS[lin], lin, crit)
        events.append({"event_id": f"best:{prot}", "time": t, "action": lin, "modality": MODALITY_FOR_LINEAGE[lin],
                       "polarity": "support", "value": max(0.0, min(1.0, value * trust)), "provenance": "bound"})
    base = json.loads((ROOT / "outputs/anm_cite_bridge/schema_O0.json").read_text())
    out = anm_run_one(schema_for(crit_id, base),
                      {"instance_id": "best", "question": "bound", "actions": copy.deepcopy(ACTIONS), "proposed_source_events": events})
    return {"bar": crit["readout_threshold"], "best_scores": {k: round(v, 4) for k, v in out["action_scores"].items()},
            "call_possible": out["recommended_action"] is not None}


def main() -> None:
    trust = modality_reliability("adt_only", {})
    print(json.dumps({"trust_adt_only": trust, **{c: best_case(c, trust) for c in ("O0", "O2")}}, indent=1))


if __name__ == "__main__":
    main()
