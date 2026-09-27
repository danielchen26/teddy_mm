#!/usr/bin/env python3
"""Probe 4 — Loading scan (Mode B evidence geometry).

Intended: re-encode TEDDY mean-pool z at ctx 256/512/1024 and compare
must-separate rates. If re-encode blocked (missing torch / GPU budget /
ckpt path), document blocker and scan available proxies:
  - leading-dim subsets of existing ctx-1024 z_512 (128/256/512)
  - (optional) PCA dims

NOT Mode A. Dim subsets are NOT true context-length variants.
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


def _labels(adt: np.ndarray, margin: float):
    name_to_j = {p: i for i, p in enumerate(PANEL)}
    adt_lin = []
    key_lin = []
    for i in range(adt.shape[0]):
        scores = {
            lin: float(np.mean([adt[i, name_to_j[p]] for p in prots]))
            for lin, prots in LINEAGE_PANELS.items()
        }
        ordered = sorted(scores.items(), key=lambda kv: -kv[1])
        if ordered[0][1] - ordered[1][1] < margin:
            adt_lin.append("abstain")
        else:
            adt_lin.append(ordered[0][0])
        best, bv = None, -np.inf
        for lin, m in KEY_MARKERS.items():
            v = float(adt[i, name_to_j[m]])
            if v > bv:
                best, bv = lin, v
        key_lin.append(best)
    return np.array(adt_lin, dtype=object), np.array(key_lin, dtype=object)


def must_for_z(z, adt, adt_lin, key_lin, *, k, min_cos, protein_l1_rel, batch=256):
    z = l2_normalize(z)
    idx, sims = knn_cosine(z, k=k, batch=batch)
    n = z.shape[0]
    scale = float(np.mean(np.abs(adt)) + 1e-8)
    seen = set()
    n_near = n_must = 0
    cells = set()
    for i in range(n):
        for t in range(idx.shape[1]):
            j = int(idx[i, t])
            if sims[i, t] < min_cos:
                continue
            a, b = (i, j) if i < j else (j, i)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            n_near += 1
            diff_lin = (
                adt_lin[i] != adt_lin[j]
                and adt_lin[i] != "abstain"
                and adt_lin[j] != "abstain"
            )
            diff_key = key_lin[i] != key_lin[j]
            rel = float(np.mean(np.abs(adt[i] - adt[j])) / scale)
            diff_prot = rel >= protein_l1_rel
            if diff_lin or diff_key or diff_prot:
                n_must += 1
                cells.add(i)
                cells.add(j)
    return {
        "n_undirected_near": n_near,
        "n_must_separate": n_must,
        "frac_must_among_near": (n_must / n_near) if n_near else None,
        "n_cells_with_ge1_must": len(cells),
        "frac_cells_with_ge1_must": len(cells) / max(n, 1),
        "z_dim": int(z.shape[1]),
    }


def check_reencode_feasibility(ckpt: Path) -> dict:
    blockers = []
    notes = []
    try:
        import torch  # noqa: F401
        notes.append("torch_import_ok")
    except Exception as e:
        blockers.append(f"torch_unavailable: {type(e).__name__}: {e}")

    if not ckpt.exists():
        blockers.append(f"ckpt_path_missing: {ckpt}")
    else:
        weights = list(ckpt.glob("*.safetensors")) + list(ckpt.glob("*.bin")) + list(ckpt.glob("*.pt"))
        if not weights:
            # maybe nested
            weights = list(ckpt.rglob("*.safetensors"))
        total = sum(w.stat().st_size for w in weights) if weights else 0
        if total < 1_000_000:
            blockers.append(f"ckpt_weights_too_small_bytes={total}")
        else:
            notes.append(f"ckpt_weights_bytes={total}")

    # raw RNA needed
    rna = ROOT / "data" / "processed" / "cite" / "cite_arrays.npz"
    if not rna.is_file():
        blockers.append(f"processed_cite_arrays_missing: {rna}")
    else:
        notes.append("processed_cite_arrays_present")

    # only ctx1024 z exists
    z_full = ROOT / "data" / "processed" / "cite" / "z_rna.npy"
    available_ctx = []
    if z_full.is_file():
        available_ctx.append({"ctx": 1024, "path": str(z_full), "note": "scripts/03_embed_rna.py default"})
    return {
        "feasible": len(blockers) == 0,
        "blockers": blockers,
        "notes": notes,
        "available_ctx_variants_on_disk": available_ctx,
        "requested_ctx": [256, 512, 1024],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridge", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--min-cos", type=float, default=0.98)
    ap.add_argument("--protein-l1-rel", type=float, default=0.35)
    ap.add_argument("--adt-margin", type=float, default=0.05)
    ap.add_argument(
        "--ckpt",
        type=Path,
        default=ROOT.parent / "teddy_mwe" / "ckpt" / "teddy_g_70M",
    )
    ap.add_argument("--dims", type=str, default="128,256,512")
    args = ap.parse_args()

    bridge = args.bridge
    out_dir = args.out_dir or (bridge / "cheap_probes" / "loading_scan")
    out_dir.mkdir(parents=True, exist_ok=True)

    feasibility = check_reencode_feasibility(args.ckpt)
    # also check Air-typical ckpt path
    air_ckpt = Path("/Users/tianchichen/Documents/GitHub/teddy_mwe/ckpt/teddy_g_70M")
    if air_ckpt.exists() and args.ckpt != air_ckpt:
        feasibility["alt_ckpt_check"] = check_reencode_feasibility(air_ckpt)

    adt = np.load(bridge / "adt_true_panel.npy")
    z = np.load(bridge / "z_rna_512.npy")[: adt.shape[0]]
    adt_lin, key_lin = _labels(adt, args.adt_margin)

    dims = [int(x) for x in args.dims.split(",") if x.strip()]
    variants = {}
    for d in dims:
        d_eff = min(d, z.shape[1])
        st = must_for_z(
            z[:, :d_eff],
            adt,
            adt_lin,
            key_lin,
            k=args.k,
            min_cos=args.min_cos,
            protein_l1_rel=args.protein_l1_rel,
        )
        st["kind"] = "leading_dim_subset_of_ctx1024_z512"
        st["requested_dim"] = d
        st["NOT_true_ctx_variant"] = True
        variants[f"leading_dim_{d_eff}"] = st

    # Signal: how much must-frac changes across loading dims (sensitivity to capacity)
    fracs = [
        variants[f"leading_dim_{min(d, z.shape[1])}"]["frac_must_among_near"]
        for d in dims
        if variants.get(f"leading_dim_{min(d, z.shape[1])}")
    ]
    fracs = [f for f in fracs if f is not None]
    if len(fracs) >= 2:
        spread = float(max(fracs) - min(fracs))
        # also relative change from full 512 to smallest
        full = variants.get("leading_dim_512", {}).get("frac_must_among_near")
        small = variants.get(f"leading_dim_{min(dims)}", {}).get("frac_must_among_near")
        rel = None if full in (None, 0) or small is None else float((small - full) / full)
    else:
        spread, rel = None, None

    results = {
        "claim_boundary": {
            "mode": "B_loading_scan",
            "not_claimed": ["Mode_A", "true_ctx_reencode_unless_feasible"],
            "note": (
                "True ctx 256/512/1024 re-encode is the intended probe. "
                "Leading-dim subsets are capacity proxies on the EXISTING ctx=1024 z_512 only."
            ),
        },
        "reencode_feasibility": feasibility,
        "params": {
            "k": args.k,
            "min_cos": args.min_cos,
            "protein_l1_rel": args.protein_l1_rel,
            "dims": dims,
            "n_cells": int(adt.shape[0]),
            "source_z": "z_rna_512.npy (ctx1024 mean-pool)",
        },
        "variants": variants,
        "signal": {
            "must_frac_spread_across_dims": spread,
            "rel_must_frac_change_small_vs_512": rel,
            "primary_metric": "must_frac_spread_across_dims",
            "interpretation": (
                "Large spread => must-separate rate depends on loading/capacity proxy "
                "(evidence geometry sensitive to how much of z is kept). "
                "Near-zero => insufficiency is stable across dim subsets of this ctx1024 z."
            ),
            "ctx_reencode_done": False,
        },
    }
    (out_dir / "loading_scan_stats.json").write_text(json.dumps(results, indent=2))

    md = ["# Loading scan (Mode B)\n\n"]
    md.append("**Not Mode A.**\n\n")
    md.append("## Re-encode feasibility (ctx 256/512/1024)\n\n")
    md.append(f"- feasible: **{feasibility['feasible']}**\n")
    for b in feasibility.get("blockers", []):
        md.append(f"- blocker: `{b}`\n")
    for n in feasibility.get("notes", []):
        md.append(f"- note: `{n}`\n")
    md.append(f"- available on disk: `{feasibility.get('available_ctx_variants_on_disk')}`\n")
    md.append("\n## Proxy: leading-dim subsets of ctx1024 z_512\n\n")
    md.append("| variant | dim | near | must | frac must |\n|---|---:|---:|---:|---:|\n")
    for name, st in variants.items():
        md.append(
            f"| {name} | {st['z_dim']} | {st['n_undirected_near']} | "
            f"{st['n_must_separate']} | {st['frac_must_among_near']} |\n"
        )
    md.append(f"\nmust_frac_spread: **{spread}** | rel small vs 512: **{rel}**\n")
    md.append("\nDo not claim Pearson > ~0.61. Not clinical.\n")
    (out_dir / "LOADING_SCAN.md").write_text("".join(md))
    print(json.dumps({"wrote": str(out_dir), "signal": results["signal"], "feasible": feasibility["feasible"]}, indent=2))


if __name__ == "__main__":
    main()
