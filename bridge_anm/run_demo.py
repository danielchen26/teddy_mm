#!/usr/bin/env python3
"""Minimal ANM×TEDDY demo: field run + leave-one-out attribution + flip-distance scan.

Uses ANM's finite_field_runner in-process (no retrain). Writes numeric results under
outputs/anm_cite_bridge/ and a REPORT.md with Chinese user-facing summary.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANM_ROOT = ROOT.parent / "ANM"
sys.path.insert(0, str(ANM_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from active_neural_matter.field.finite_field_runner import (  # noqa: E402
    build_graph,
    evolve_field,
    load_instances,
    load_json,
    readout,
    run_instances,
    validate_schema,
    validate_source_events,
    verify,
)

from lib.lineage_panels import CRITERIA  # noqa: E402


def _run_one(schema: dict, instance: dict) -> dict[str, Any]:
    validation = validate_source_events(schema, instance)
    graph = build_graph(instance, validation["field_events"])
    field = evolve_field(schema, graph, validation["field_events"])
    rd = readout(schema, instance, field["state"])
    vf = verify(schema, instance, rd)
    return {
        "recommended_action": rd["recommended_action"],
        "action_scores": rd["action_scores"],
        "P_f": vf["P_f"],
        "Q_f": vf["Q_f"],
        "r_f": vf["r_f"],
        "expected_action": vf.get("expected_action"),
        "n_admitted": len(validation["field_events"]),
    }


def leave_one_out_attribution(schema: dict, instance: dict) -> list[dict[str, Any]]:
    base = _run_one(schema, instance)
    base_action = base["recommended_action"]
    base_scores = base["action_scores"]
    rows = []
    events = list(instance["proposed_source_events"])
    for i, ev in enumerate(events):
        ablated = copy.deepcopy(instance)
        ablated["proposed_source_events"] = [e for j, e in enumerate(events) if j != i]
        out = _run_one(schema, ablated)
        delta_scores = {
            a: float(out["action_scores"].get(a, 0.0) - base_scores.get(a, 0.0))
            for a in set(base_scores) | set(out["action_scores"])
        }
        rows.append(
            {
                "event_id": ev["event_id"],
                "protein": ev["event_id"].rsplit(":", 1)[-1],
                "action": ev["action"],
                "value": ev["value"],
                "base_recommended": base_action,
                "ablated_recommended": out["recommended_action"],
                "flipped": out["recommended_action"] != base_action,
                "delta_scores": delta_scores,
                "delta_chosen_score": float(
                    (out["action_scores"].get(base_action, 0.0) if base_action else 0.0)
                    - (base_scores.get(base_action, 0.0) if base_action else 0.0)
                ),
            }
        )
    rows.sort(key=lambda r: abs(r["delta_chosen_score"]))
    return rows


def flip_distance_scan(schema: dict, instance: dict, event_index: int, grid: int = 21) -> dict[str, Any]:
    """Scan one event's value in [0,1]; report minimal |Δvalue| that flips readout."""
    base = _run_one(schema, instance)
    base_action = base["recommended_action"]
    ev0 = instance["proposed_source_events"][event_index]
    v0 = float(ev0["value"])
    flips = []
    for v in np.linspace(0.0, 1.0, grid):
        trial = copy.deepcopy(instance)
        trial["proposed_source_events"][event_index]["value"] = float(v)
        out = _run_one(schema, trial)
        if out["recommended_action"] != base_action:
            flips.append({"value": float(v), "recommended": out["recommended_action"], "delta": abs(float(v) - v0)})
    if not flips:
        return {
            "event_id": ev0["event_id"],
            "base_value": v0,
            "base_recommended": base_action,
            "flip_exists": False,
            "min_flip_distance": None,
        }
    best = min(flips, key=lambda x: x["delta"])
    return {
        "event_id": ev0["event_id"],
        "base_value": v0,
        "base_recommended": base_action,
        "flip_exists": True,
        "min_flip_distance": best["delta"],
        "flip_to": best["recommended"],
        "flip_value": best["value"],
    }


def summarize_artifact(artifact: dict) -> dict[str, Any]:
    insts = artifact.get("instances", [])
    n = len(insts)
    p_sum = q_sum = q_n = 0.0
    action_counts: dict[str, int] = {}
    abstain = 0
    for inst in insts:
        # artifact structure from run_instances
        rd = inst.get("readout", {}).get("baseline", inst.get("readout", {}))
        vf = inst.get("verifier", {}).get("baseline", inst.get("verifier", {}))
        # tolerate nested condition dicts
        if isinstance(rd, dict) and "recommended_action" not in rd:
            rd = rd.get("baseline") or next(iter(rd.values()), {})
        if isinstance(vf, dict) and "P_f" not in vf:
            vf = vf.get("baseline") or next(iter(vf.values()), {})
        rec = rd.get("recommended_action") if isinstance(rd, dict) else None
        if rec is None:
            abstain += 1
        else:
            action_counts[rec] = action_counts.get(rec, 0) + 1
        p_sum += float(vf.get("P_f") or 0.0)
        if vf.get("Q_f") is not None:
            q_sum += float(vf["Q_f"])
            q_n += 1
    return {
        "n_instances": n,
        "mean_P_f": p_sum / max(n, 1),
        "mean_Q_f_defined": (q_sum / q_n) if q_n else None,
        "n_Q_f_defined": int(q_n),
        "abstain_count": abstain,
        "action_counts": action_counts,
    }


def write_report(
    path: Path,
    *,
    export_manifest: dict,
    o0_summary: dict,
    o1_summary: dict,
    attrib_example: dict,
    flip_example: dict,
    runner_ok: bool,
) -> None:
    pearson = export_manifest.get("phase1_mlp_pearson_quoted", {})
    mean_r = export_manifest.get("phase1_test_mlp_pearson_mean")
    lines = []
    lines.append("# ANM × TEDDY CITE Bridge — v0 REPORT")
    lines.append("")
    lines.append("## 中文摘要（给 Daniel）")
    lines.append("")
    lines.append(
        f"本桥接把 TEDDY phase-1 CITE（RNA→冻结 TEDDY-G→z→MLP→ADT）的 **类型化证据 δu** "
        f"导出为 ANM 有限场实例，并在声明的 P_f / Q_f 读出下跑通场演化、留一归因与翻转距离。"
        f" **不是** 临床优越性声明；成功标准是：事件可导出、schema 已声明、真实 CITE 事件上响应/归因路径可跑、"
        f"准则 O0→O1 可改声明而不重训 TEDDY。"
    )
    lines.append("")
    lines.append(f"- 导出细胞数：**{export_manifest.get('n_cells')}**；事件数：**{export_manifest.get('n_events')}**（test / site4）。")
    lines.append(
        f"- phase-1 测试集整体 MLP Pearson（metrics.json）：**{mean_r:.6f}**。"
        if mean_r is not None
        else "- phase-1 整体 Pearson：见 metrics.json。"
    )
    lines.append("- 面板蛋白 MLP Pearson（test_per_protein.json，原样引用）：")
    for prot, r in sorted(pearson.items(), key=lambda kv: -kv[1]):
        lines.append(f"  - `{prot}`: **{r:.6f}**")
    lines.append(
        f"- O0 场跑通：`mean_P_f={o0_summary.get('mean_P_f'):.4f}`，"
        f"`mean_Q_f(defined)={o0_summary.get('mean_Q_f_defined')}`，"
        f"abstain={o0_summary.get('abstain_count')} / {o0_summary.get('n_instances')}。"
    )
    lines.append(
        f"- O1（仅改声明：threshold 0.12→0.28，key_marker_boost 1→2）："
        f"`mean_P_f={o1_summary.get('mean_P_f'):.4f}`，"
        f"`mean_Q_f(defined)={o1_summary.get('mean_Q_f_defined')}`，"
        f"abstain={o1_summary.get('abstain_count')} / {o1_summary.get('n_instances')}。"
    )
    lines.append(f"- ANM `finite_field_runner` 实际执行：**{'是' if runner_ok else '否'}**。")
    if flip_example.get("flip_exists"):
        lines.append(
            f"- 翻转距离示例（事件 `{flip_example.get('event_id')}`）："
            f" min |Δvalue| = **{flip_example.get('min_flip_distance'):.4f}** "
            f"（{flip_example.get('base_recommended')} → {flip_example.get('flip_to')}）。"
        )
    else:
        lines.append(f"- 翻转距离示例：在 [0,1] 网格上 **未找到** 翻转（event `{flip_example.get('event_id')}`）。")
    top = attrib_example.get("top_events") or []
    if top:
        lines.append("- 留一归因（示例细胞，按 |Δ chosen score|）：")
        for row in top[:5]:
            lines.append(
                f"  - `{row['protein']}` action={row['action']} "
                f"Δscore={row['delta_chosen_score']:.4f} flipped={row['flipped']}"
            )
    lines.append("")
    lines.append("### 准则编辑（O0→O1）而不重训")
    lines.append("")
    lines.append(
        "只改 `bridge_anm/readouts/cite_lineage_O*.yaml` / `lib/lineage_panels.CRITERIA` 与 "
        "实例上的 `readout_threshold` / `key_marker_boost`，重新 `adapt_to_anm.py` + `run_demo.py`。"
        " **不** 触碰 `best.pt`、不重跑 phase-1 训练。"
    )
    lines.append("")
    lines.append("## English technical notes")
    lines.append("")
    lines.append("- TEDDY role: typed evidence engine (δu).")
    lines.append("- ANM role: declared field + nested P_f / Q_f + attribution / flip distance.")
    lines.append("- Holdout `adt_true` is verifier-only; never admitted as source events.")
    lines.append("- Claim boundary: local response diagnosis only.")
    lines.append("")
    lines.append("## Numeric blobs")
    lines.append("")
    lines.append("```json")
    lines.append(
        json.dumps(
            {
                "O0": o0_summary,
                "O1": o1_summary,
                "attribution_example": attrib_example,
                "flip_example": flip_example,
            },
            indent=2,
        )
    )
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--attrib-limit", type=int, default=30, help="instances for LOO attribution")
    p.add_argument("--flip-limit", type=int, default=30, help="instances for flip-distance scan")
    args = p.parse_args()

    out = args.out_dir
    export_manifest = json.loads((out / "export_manifest.json").read_text())

    summaries = {}
    artifacts = {}
    runner_ok = True
    t0 = time.perf_counter()
    for crit in ("O0", "O1"):
        schema_path = out / f"schema_{crit}.json"
        inst_path = out / f"anm_instances_{crit}.json"
        schema = load_json(schema_path)
        validate_schema(schema)
        bundle = load_instances(inst_path)
        art_path = out / f"field_artifact_{crit}.json"
        try:
            artifact = run_instances(
                schema,
                bundle,
                repo_root=ANM_ROOT,
                output_path=art_path,
                seed=args.seed,
                mode="dry_run",
            )
            art_path.write_text(json.dumps(artifact, indent=2, sort_keys=True))
            artifacts[crit] = artifact
            summaries[crit] = summarize_artifact(artifact)
        except Exception as exc:  # noqa: BLE001
            runner_ok = False
            summaries[crit] = {"error": str(exc)}
            print(f"ERROR running {crit}: {exc}")

    # Attribution + flip on O0 first workable instance
    schema0 = load_json(out / "schema_O0.json")
    bundle0 = load_instances(out / "anm_instances_O0.json")
    attrib_example: dict[str, Any] = {}
    flip_example: dict[str, Any] = {"flip_exists": False}
    flip_rows = []
    for inst in bundle0["instances"][: args.attrib_limit]:
        base = _run_one(schema0, inst)
        if base["recommended_action"] is None:
            continue
        rows = leave_one_out_attribution(schema0, inst)
        attrib_example = {
            "instance_id": inst["instance_id"],
            "recommended_action": base["recommended_action"],
            "P_f": base["P_f"],
            "Q_f": base["Q_f"],
            "action_scores": base["action_scores"],
            "top_events": sorted(rows, key=lambda r: abs(r["delta_chosen_score"]), reverse=True),
        }
        # flip-scan the event with largest |Δ|
        if rows:
            # map event_id back to index
            eid = attrib_example["top_events"][0]["event_id"]
            idx = next(i for i, e in enumerate(inst["proposed_source_events"]) if e["event_id"] == eid)
            flip_example = flip_distance_scan(schema0, inst, idx)
        break

    for inst in bundle0["instances"][: args.flip_limit]:
        base = _run_one(schema0, inst)
        if base["recommended_action"] is None:
            continue
        for idx in range(len(inst["proposed_source_events"])):
            fd = flip_distance_scan(schema0, inst, idx, grid=21)
            if fd["flip_exists"]:
                flip_rows.append({"instance_id": inst["instance_id"], **fd})
        if len(flip_rows) >= 50:
            break

    flip_summary = {
        "n_scanned_with_flip": len(flip_rows),
        "median_min_flip_distance": float(np.median([r["min_flip_distance"] for r in flip_rows]))
        if flip_rows
        else None,
        "example": flip_example,
        "samples": flip_rows[:10],
    }

    results = {
        "wall_seconds": time.perf_counter() - t0,
        "runner_executed": runner_ok,
        "summaries": summaries,
        "attribution_example": attrib_example,
        "flip_summary": flip_summary,
        "criteria": {k: CRITERIA[k] for k in ("O0", "O1")},
    }
    (out / "demo_results.json").write_text(json.dumps(results, indent=2))
    write_report(
        out / "REPORT.md",
        export_manifest=export_manifest,
        o0_summary=summaries.get("O0", {}),
        o1_summary=summaries.get("O1", {}),
        attrib_example=attrib_example,
        flip_example=flip_example,
        runner_ok=runner_ok,
    )
    readme = out / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# ANM × TEDDY CITE bridge outputs",
                "",
                "Generated by `teddy_mm/bridge_anm/`.",
                "",
                "## Re-run",
                "",
                "```bash",
                "cd /Users/tianchichen/Documents/GitHub/teddy_mm",
                "PYTORCH_ENABLE_MPS_FALLBACK=1 .venv/bin/python bridge_anm/export_cite_events.py --n-cells 2000",
                ".venv/bin/python bridge_anm/adapt_to_anm.py --criterion O0 --max-cells 500",
                ".venv/bin/python bridge_anm/adapt_to_anm.py --criterion O1 --max-cells 500",
                "PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_demo.py",
                "```",
                "",
                "See `REPORT.md` for quoted numbers.",
                "",
            ]
        )
    )
    print(json.dumps({"demo_results": str(out / "demo_results.json"), "report": str(out / "REPORT.md"), **{k: summaries[k] for k in summaries}}, indent=2, default=str))


if __name__ == "__main__":
    main()
