#!/usr/bin/env python3
"""Run all cheap Mode B probes and write comparative ranking report."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridge", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--gated-max-cells", type=int, default=2000)
    ap.add_argument("--eps", type=float, default=0.10)
    ap.add_argument("--skip-gated", action="store_true")
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    bridge = args.bridge
    out_root = bridge / "cheap_probes"
    out_root.mkdir(parents=True, exist_ok=True)
    py = args.python
    t0 = time.time()

    run([py, str(HERE / "probe_must_separate_autopsy.py"), "--bridge", str(bridge)])
    if not args.skip_gated:
        run(
            [
                py,
                str(HERE / "probe_gated_loo.py"),
                "--bridge",
                str(bridge),
                "--max-cells",
                str(args.gated_max_cells),
                "--eps",
                str(args.eps),
            ]
        )
    run([py, str(HERE / "probe_complement_zplus.py"), "--bridge", str(bridge)])
    run([py, str(HERE / "probe_loading_scan.py"), "--bridge", str(bridge)])

    # Load signals
    autopsy = json.loads((out_root / "must_separate_autopsy" / "autopsy_stats.json").read_text())
    gated_path = out_root / "gated_loo" / "gated_loo_stats.json"
    gated = json.loads(gated_path.read_text()) if gated_path.is_file() else None
    comp = json.loads((out_root / "complement_zplus" / "complement_zplus_stats.json").read_text())
    load = json.loads((out_root / "loading_scan" / "loading_scan_stats.json").read_text())

    # Normalize to comparable "signal strength" in [0,1]-ish magnitudes for ranking.
    # Each probe has its own primary metric; we rank by |primary| with direction notes.
    rows = []

    a_sig = float(autopsy["signal"]["bio_minus_noise"])
    rows.append(
        {
            "probe": "1_must_separate_autopsy",
            "primary_metric": "bio_minus_noise",
            "primary_value": a_sig,
            "abs_for_rank": abs(a_sig),
            "direction": "higher bio_minus_noise => stronger bio (vs noise-like) in must-pairs",
            "extra": {
                "bio_frac": autopsy["signal"]["bio_frac"],
                "noise_like_frac": autopsy["signal"]["noise_like_frac"],
                "top_protein": autopsy["protein_rank_by_mean_abs_delta"][0]["protein"],
                "top_protein_mean_abs_delta": autopsy["protein_rank_by_mean_abs_delta"][0][
                    "mean_abs_delta"
                ],
            },
        }
    )

    if gated is not None:
        g_sig = float(gated["signal"]["any_gated_flip_frac"])
        rows.append(
            {
                "probe": "2_gated_loo",
                "primary_metric": "any_gated_flip_frac",
                "primary_value": g_sig,
                "abs_for_rank": abs(g_sig),
                "direction": "higher => more decisions flip under ±ε/±ε/2 on top LOO protein",
                "extra": {
                    "strong_G3_G4_frac": gated["signal"]["strong_G3_G4_frac"],
                    "G4_asym_frac": gated["signal"]["G4_asym_frac"],
                    "gate_counts": gated["gate_counts"],
                    "loo_flip_rate_full": gated["existing_LOO_baseline"]["full_top1_flip_rate"],
                    "delta_vs_loo_flip": gated["signal"]["delta_vs_loo_flip"],
                },
            }
        )
    else:
        rows.append(
            {
                "probe": "2_gated_loo",
                "primary_metric": "any_gated_flip_frac",
                "primary_value": None,
                "abs_for_rank": -1.0,
                "direction": "skipped",
                "extra": {"skipped": True},
            }
        )

    c_sig = comp["signal"]["primary_score"]
    c_sig_f = float(c_sig) if c_sig is not None else None
    rows.append(
        {
            "probe": "3_complement_zplus",
            "primary_metric": "relative_must_reduction_vs_z512",
            "primary_value": c_sig_f,
            "abs_for_rank": abs(c_sig_f) if c_sig_f is not None else -1.0,
            "direction": "higher => z+ reduces must-among-near vs z_512 on held proteins",
            "extra": {
                "scores": comp.get("relative_must_reduction_vs_z512"),
                "site4_fracs": {
                    k: v.get("frac_must_among_near")
                    for k, v in comp.get("site4_variants", {}).items()
                },
            },
        }
    )

    l_sig = load["signal"]["must_frac_spread_across_dims"]
    l_sig_f = float(l_sig) if l_sig is not None else None
    rows.append(
        {
            "probe": "4_loading_scan",
            "primary_metric": "must_frac_spread_across_dims",
            "primary_value": l_sig_f,
            "abs_for_rank": abs(l_sig_f) if l_sig_f is not None else -1.0,
            "direction": "higher spread => must-rate depends on loading/capacity proxy",
            "extra": {
                "rel_small_vs_512": load["signal"]["rel_must_frac_change_small_vs_512"],
                "ctx_reencode_done": load["signal"]["ctx_reencode_done"],
                "reencode_feasible": load["reencode_feasibility"]["feasible"],
                "blockers": load["reencode_feasibility"].get("blockers"),
                "variant_fracs": {
                    k: v.get("frac_must_among_near") for k, v in load.get("variants", {}).items()
                },
            },
        }
    )

    ranked = sorted(rows, key=lambda r: r["abs_for_rank"], reverse=True)
    for i, r in enumerate(ranked, 1):
        r["rank"] = i

    summary = {
        "claim_boundary": {
            "mode": "B_cheap_next_probes_comparative",
            "not_claimed": [
                "Mode_A_residual_stream_as_field",
                "layer_Jacobian",
                "in_silico_gene_perturbs",
                "reverse_closed_loop_perturb_response",
                "fusion_audit",
                "clinical",
            ],
        },
        "runtime_sec": time.time() - t0,
        "ranking": ranked,
        "strongest_probe": ranked[0]["probe"] if ranked else None,
        "headline_numbers": {
            r["probe"]: {"metric": r["primary_metric"], "value": r["primary_value"]}
            for r in ranked
        },
    }
    (out_root / "comparative_ranking.json").write_text(json.dumps(summary, indent=2))

    # Markdown report for docs/
    md = []
    md.append("# Cheap next probes — comparative ranking (Mode B)\n\n")
    md.append("**Not Mode A.** Evidence-only ranking of four cheap probes on teddy_mm.\n\n")
    md.append(f"Runtime: {summary['runtime_sec']:.1f}s\n\n")
    md.append("## Ranking (by |primary metric|)\n\n")
    md.append("| rank | probe | primary metric | value | note |\n|---:|---|---|---:|---|\n")
    for r in ranked:
        md.append(
            f"| {r['rank']} | `{r['probe']}` | `{r['primary_metric']}` | "
            f"{r['primary_value']} | {r['direction']} |\n"
        )
    md.append(f"\n**Strongest signal:** `{summary['strongest_probe']}`\n\n")

    md.append("## Per-probe evidence\n\n")
    md.append("### 1) MUST-SEPARATE autopsy\n")
    md.append(
        f"- bio_frac={autopsy['signal']['bio_frac']:.4f}, "
        f"noise_like_frac={autopsy['signal']['noise_like_frac']:.4f}, "
        f"bio_minus_noise={a_sig:.4f}\n"
    )
    top3 = autopsy["protein_rank_by_mean_abs_delta"][:3]
    md.append(
        "- top proteins by mean |Δ|: "
        + ", ".join(f"{t['protein']}={t['mean_abs_delta']:.4f}" for t in top3)
        + "\n"
    )
    md.append(f"- buckets: `{dict(autopsy['bucket_counts'])}`\n\n")

    if gated is not None:
        md.append("### 2) Gated LOO ±ε / ±ε/2 + G1–G4\n")
        md.append(
            f"- eps={gated['params']['eps']}, kept={gated['params']['n_kept_non_abstain']}\n"
        )
        md.append(f"- LOO top1 flip (full): {gated['existing_LOO_baseline']['full_top1_flip_rate']:.4f}\n")
        md.append(f"- any gated flip: {gated['signal']['any_gated_flip_frac']:.4f}\n")
        md.append(f"- G3+G4 strong: {gated['signal']['strong_G3_G4_frac']:.4f}\n")
        md.append(f"- gate_counts: `{gated['gate_counts']}`\n\n")
    else:
        md.append("### 2) Gated LOO — skipped\n\n")

    md.append("### 3) Complement z⁺\n")
    md.append(f"- primary relative must-reduction: **{c_sig_f}**\n")
    for k, v in comp.get("site4_variants", {}).items():
        md.append(f"- `{k}` frac_must={v.get('frac_must_among_near')}\n")
    md.append(f"- phase2 note: {comp.get('phase2_note')}\n\n")

    md.append("### 4) Loading scan\n")
    md.append(f"- ctx re-encode feasible: **{load['reencode_feasibility']['feasible']}**\n")
    for b in load["reencode_feasibility"].get("blockers", []):
        md.append(f"- blocker: `{b}`\n")
    md.append(f"- must_frac_spread across dim proxies: **{l_sig_f}**\n")
    for k, v in load.get("variants", {}).items():
        md.append(f"- `{k}` frac_must={v.get('frac_must_among_near')}\n")

    md.append("\n## Claim boundary\n\n")
    md.append("Mode B only: typed evidence / decision sensitivity / evidence geometry. ")
    md.append("Do not claim Mode A residual-as-field, Jacobian, gene perturbs, or clinical utility. ")
    md.append("Do not claim Pearson > ~0.61.\n")

    report_path = ROOT / "docs" / "reports" / "CHEAP_PROBES_RANKING.md"
    report_path.write_text("".join(md))
    # also mirror under outputs
    (out_root / "CHEAP_PROBES_RANKING.md").write_text("".join(md))
    print(json.dumps({"strongest": summary["strongest_probe"], "ranking": ranked, "report": str(report_path)}, indent=2))


if __name__ == "__main__":
    main()
