#!/usr/bin/env python3
"""Probe 3 — Complement z⁺ (Mode B evidence geometry).

Build z⁺ from ADT features (held-protein protocol) and/or TEDDY ADT preds;
measure must-separate ONLY on held protein panel (and OOD if available).

Compares whether a complement evidence channel separates what z_512 cannot.
NOT Mode A / not fusion audit claim.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.knn_cosine import knn_cosine, l2_normalize  # noqa: E402
from lib.lineage_panels import KEY_MARKERS, LINEAGE_PANELS, all_panel_proteins  # noqa: E402

PANEL = all_panel_proteins()


def _adt_lineage_label(row: np.ndarray, name_to_j: dict[str, int], margin: float) -> str:
    scores = {}
    for lin, prots in LINEAGE_PANELS.items():
        scores[lin] = float(np.mean([row[name_to_j[p]] for p in prots]))
    ordered = sorted(scores.items(), key=lambda kv: -kv[1])
    if ordered[0][1] - ordered[1][1] < margin:
        return "abstain"
    return ordered[0][0]


def _key_label(row: np.ndarray, name_to_j: dict[str, int]) -> str:
    best, bv = None, -np.inf
    for lin, m in KEY_MARKERS.items():
        v = float(row[name_to_j[m]])
        if v > bv:
            best, bv = lin, v
    return str(best)


def must_stats(
    z: np.ndarray,
    adt_eval: np.ndarray,
    eval_protein_idx: list[int],
    *,
    k: int,
    min_cos: float,
    protein_l1_rel: float,
    adt_margin: float,
    name_to_j: dict[str, int],
    batch: int = 256,
) -> dict:
    """Near pairs in z; must-separate judged on held eval proteins / lineage."""
    idx, sims = knn_cosine(z, k=k, batch=batch)
    n = z.shape[0]
    # lineage labels from FULL panel (ground truth) — but protein-profile must uses held only
    adt_lin = np.array(
        [_adt_lineage_label(adt_eval[i], name_to_j, adt_margin) for i in range(n)],
        dtype=object,
    )
    key_lin = np.array([_key_label(adt_eval[i], name_to_j) for i in range(n)], dtype=object)

    held = adt_eval[:, eval_protein_idx].astype(np.float64)
    scale = float(np.mean(np.abs(held)) + 1e-8)

    seen = set()
    n_near = n_must = n_diff_lin = n_diff_key = n_diff_prot = 0
    cells_must = set()
    for i in range(n):
        for t in range(idx.shape[1]):
            j = int(idx[i, t])
            if j < 0 or sims[i, t] < min_cos:
                continue
            a, b = (i, j) if i < j else (j, i)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            n_near += 1
            diff_lin = adt_lin[i] != adt_lin[j] and adt_lin[i] != "abstain" and adt_lin[j] != "abstain"
            diff_key = key_lin[i] != key_lin[j]
            rel = float(np.mean(np.abs(held[i] - held[j])) / scale)
            diff_prot = rel >= protein_l1_rel
            if diff_lin:
                n_diff_lin += 1
            if diff_key:
                n_diff_key += 1
            if diff_prot:
                n_diff_prot += 1
            if diff_lin or diff_key or diff_prot:
                n_must += 1
                cells_must.add(i)
                cells_must.add(j)
    return {
        "n_cells": n,
        "z_dim": int(z.shape[1]),
        "n_undirected_near": n_near,
        "n_must_separate": n_must,
        "frac_must_among_near": (n_must / n_near) if n_near else None,
        "n_cells_with_ge1_must": len(cells_must),
        "frac_cells_with_ge1_must": len(cells_must) / max(n, 1),
        "n_diff_adt_lineage": n_diff_lin,
        "n_diff_key_marker": n_diff_key,
        "n_diff_held_protein_profile": n_diff_prot,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridge", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--min-cos", type=float, default=0.98)
    ap.add_argument("--protein-l1-rel", type=float, default=0.35)
    ap.add_argument("--adt-margin", type=float, default=0.05)
    ap.add_argument("--held-lineage", type=str, default="myeloid",
                    help="hold out this lineage's 3 proteins for eval; z+ uses the other 6")
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    bridge = args.bridge
    out_dir = args.out_dir or (bridge / "cheap_probes" / "complement_zplus")
    out_dir.mkdir(parents=True, exist_ok=True)

    name_to_j = {p: i for i, p in enumerate(PANEL)}
    held_prots = list(LINEAGE_PANELS[args.held_lineage])
    zplus_prots = [p for p in PANEL if p not in held_prots]
    held_idx = [name_to_j[p] for p in held_prots]
    zplus_idx = [name_to_j[p] for p in zplus_prots]

    adt_true = np.load(bridge / "adt_true_panel.npy")
    adt_pred_path = bridge / "adt_pred_panel.npy"
    adt_pred = np.load(adt_pred_path) if adt_pred_path.is_file() else None
    z512_path = bridge / "z_rna_512.npy"
    z512 = np.load(z512_path)[: adt_true.shape[0]] if z512_path.is_file() else None

    common = dict(
        k=args.k,
        min_cos=args.min_cos,
        protein_l1_rel=args.protein_l1_rel,
        adt_margin=args.adt_margin,
        name_to_j=name_to_j,
    )

    variants = {}
    # Baseline: z_512 must-separate judged on held proteins only
    if z512 is not None:
        variants["z512_held_eval"] = must_stats(z512, adt_true, held_idx, **common)
        variants["z512_held_eval"]["z_source"] = "z_rna_512 mean-pool ctx1024"

    # z+ = ADT true of complementary (non-held) proteins
    z_plus_true = l2_normalize(adt_true[:, zplus_idx])
    variants["zplus_adt_true_complement"] = must_stats(z_plus_true, adt_true, held_idx, **common)
    variants["zplus_adt_true_complement"]["z_source"] = (
        f"L2 ADT true on non-held proteins {zplus_prots}; eval held {held_prots}"
    )

    # z+ = TEDDY ADT preds on complementary proteins (if available)
    if adt_pred is not None:
        z_plus_pred = l2_normalize(adt_pred[:, zplus_idx])
        variants["zplus_adt_pred_complement"] = must_stats(z_plus_pred, adt_true, held_idx, **common)
        variants["zplus_adt_pred_complement"]["z_source"] = (
            f"L2 TEDDY ADT preds on non-held {zplus_prots}; eval held {held_prots}"
        )
        # also full pred panel as z+ (eval still held) — upper bound leakage check
        z_plus_pred_all = l2_normalize(adt_pred)
        variants["zplus_adt_pred_all9_leakcheck"] = must_stats(
            z_plus_pred_all, adt_true, held_idx, **common
        )
        variants["zplus_adt_pred_all9_leakcheck"]["z_source"] = (
            "L2 TEDDY ADT preds all 9 (includes held — leakage check, not primary)"
        )

    # OOD if available
    ood_true_p = bridge / "ood_adt_true_panel.npy"
    ood_pred_p = bridge / "ood_adt_pred_panel.npy"
    ood = {}
    if ood_true_p.is_file():
        ood_true = np.load(ood_true_p)
        ood_common = dict(common)
        if z512 is not None and z512.shape[0] >= adt_true.shape[0] + ood_true.shape[0]:
            z_ood = z512[adt_true.shape[0] : adt_true.shape[0] + ood_true.shape[0]]
            ood["z512_held_eval"] = must_stats(z_ood, ood_true, held_idx, **ood_common)
        ood["zplus_adt_true_complement"] = must_stats(
            l2_normalize(ood_true[:, zplus_idx]), ood_true, held_idx, **ood_common
        )
        if ood_pred_p.is_file():
            ood_pred = np.load(ood_pred_p)
            ood["zplus_adt_pred_complement"] = must_stats(
                l2_normalize(ood_pred[:, zplus_idx]), ood_true, held_idx, **ood_common
            )

    # phase-2 joint blocker / note
    phase2_note = {
        "phase2_joint_embedding": "not available as exported z; cite_phase2 has best.pt + 4k raw ADT(134) but no aligned site4 joint z sidecar",
        "used_instead": "ADT true/pred complement features on site4 panel (9) with held-lineage protocol",
    }

    # Signal: reduction in frac_must vs z512 baseline (complement separates better => lower must frac)
    base = variants.get("z512_held_eval", {})
    base_frac = base.get("frac_must_among_near")
    scores = {}
    for name, st in variants.items():
        if name.startswith("z512"):
            continue
        frac = st.get("frac_must_among_near")
        if base_frac is None or frac is None or base_frac <= 0:
            scores[name] = None
        else:
            # positive => complement reduces must-among-near (stronger separation signal)
            scores[name] = float((base_frac - frac) / base_frac)

    primary = "zplus_adt_true_complement"
    primary_score = scores.get(primary)

    results = {
        "claim_boundary": {
            "mode": "B_complement_zplus_held_protein_protocol",
            "not_claimed": [
                "Mode_A",
                "fusion_audit",
                "phase2_joint_as_implemented",
            ],
            "held_protocol": (
                f"z+ built from proteins excluding held lineage={args.held_lineage} "
                f"({held_prots}); must-separate judged on held protein profile + lineage labels"
            ),
        },
        "params": {
            "k": args.k,
            "min_cos": args.min_cos,
            "protein_l1_rel": args.protein_l1_rel,
            "held_lineage": args.held_lineage,
            "held_proteins": held_prots,
            "zplus_proteins": zplus_prots,
        },
        "phase2_note": phase2_note,
        "site4_variants": variants,
        "ood_variants": ood,
        "relative_must_reduction_vs_z512": scores,
        "signal": {
            "primary_variant": primary,
            "primary_metric": "relative_must_reduction_vs_z512",
            "primary_score": primary_score,
            "interpretation": (
                "Positive score => z+ near-neighborhoods contain fewer held-protein/lineage "
                "disagreements than z_512 (complement channel separates what z cannot)."
            ),
        },
    }
    (out_dir / "complement_zplus_stats.json").write_text(json.dumps(results, indent=2))

    md = ["# Complement z⁺ (Mode B)\n\n", "**Not Mode A / not fusion audit.**\n\n"]
    md.append(f"- held lineage: `{args.held_lineage}` proteins {held_prots}\n")
    md.append(f"- z+ proteins: {zplus_prots}\n\n")
    md.append("## site4 variants (must among near @ min_cos)\n\n")
    md.append("| variant | near | must | frac must | cells w/ must |\n|---|---:|---:|---:|---:|\n")
    for name, st in variants.items():
        md.append(
            f"| {name} | {st.get('n_undirected_near')} | {st.get('n_must_separate')} | "
            f"{st.get('frac_must_among_near')} | {st.get('n_cells_with_ge1_must')} |\n"
        )
    md.append(f"\nPrimary relative must-reduction ({primary}): **{primary_score}**\n")
    if ood:
        md.append("\n## OOD variants\n\n")
        for name, st in ood.items():
            md.append(f"- `{name}`: frac_must={st.get('frac_must_among_near')} near={st.get('n_undirected_near')}\n")
    md.append(f"\nPhase2 note: {phase2_note}\n")
    md.append("\nDo not claim Pearson > ~0.61. Not clinical.\n")
    (out_dir / "COMPLEMENT_ZPLUS.md").write_text("".join(md))
    print(json.dumps({"wrote": str(out_dir), "signal": results["signal"], "scores": scores}, indent=2))


if __name__ == "__main__":
    main()
