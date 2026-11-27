#!/usr/bin/env python3
"""The units of the final audit, from the census of 22 November. results/roles.json. Names only; nothing is opened.

  python roles.py
"""
import json
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent


def main():
    census = json.loads((RUNS / "outer_step_20261122" / "results" / "census.json").read_text())
    final = sorted(census["roles"]["final"])
    read = set(census["batch_1"]["units"])
    for run in ("batch2_forecast_20261125", "batch3_frozen_20261126", "batch4_frozen_20261127"):
        read |= set(json.loads((RUNS / run / "results" / "roles.json").read_text())["units"])
    assert not set(final) & read, "a final unit was read in a batch"
    by = {}
    for u in final:
        by.setdefault(u.split("@")[1], []).append(u.split("@")[0])
    out = {"rule": {"final_units": census["rule"]["final"], "final_control_cells": census["rule"]["final_control_cells"],
                    "reference": "a final unit is read against the final control cells of its stratum, which no run has read; where a stratum has none (THP-1: every control cell was read before the census), against its control cells"},
           "units": final, "strata": {k: {"units": len(v), "genes": v} for k, v in sorted(by.items())}}
    (RUN / "results" / "roles.json").write_text(json.dumps(out, indent=1) + "\n")
    print(len(final), "final units:", {k: v["units"] for k, v in out["strata"].items()})


if __name__ == "__main__":
    main()
