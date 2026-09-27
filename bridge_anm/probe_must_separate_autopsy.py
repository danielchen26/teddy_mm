#!/usr/bin/env python3
"""Probe 1 — MUST-SEPARATE autopsy (Mode B evidence geometry).

From reverse_step1 near-identical z_512 pairs, decompose which proteins /
lineages disagree and score noise-like vs bio signal.

NOT Mode A / not perturb-response.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.lineage_panels import KEY_MARKERS, LINEAGE_PANELS, all_panel_proteins  # noqa: E402

PANEL = all_panel_proteins()


def _load_pairs(path: Path, max_n: int | None = None) -> list[dict]:
    rows = []
    with path.open() as f:
        for line in f:
            rows.append(json.loads(line))
            if max_n is not None and len(rows) >= max_n:
                break
    return rows


def _protein_deltas(adt: np.ndarray, i: int, j: int) -> np.ndarray:
    return np.abs(adt[i] - adt[j]).astype(np.float64)


def classify_pair(flags: dict, coarse_same: bool, key_same: bool) -> str:
    """Noise-like vs bio buckets (Mode B label geometry, not causal)."""
    if flags.get("diff_adt_lineage"):
        return "bio_adt_lineage"
    if flags.get("diff_coarse_cell_type"):
        return "bio_coarse_cell_type"
    if flags.get("diff_key_marker_lineage"):
        return "bio_key_marker"
    if flags.get("diff_protein_profile") and coarse_same and key_same:
        return "noise_like_same_type_protein_drift"
    if flags.get("diff_protein_profile"):
        return "fine_protein_profile"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridge", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--pairs", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--max-pairs", type=int, default=None)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    bridge = args.bridge
    pairs_path = args.pairs or (bridge / "reverse_step1" / "must_separate_pairs.jsonl")
    out_dir = args.out_dir or (bridge / "cheap_probes" / "must_separate_autopsy")
    out_dir.mkdir(parents=True, exist_ok=True)

    adt = np.load(bridge / "adt_true_panel.npy")
    pairs = _load_pairs(pairs_path, args.max_pairs)
    if not pairs:
        raise SystemExit(f"no pairs in {pairs_path}")

    # Per-protein |Δ| among must pairs
    prot_abs = np.zeros((len(pairs), len(PANEL)), dtype=np.float64)
    bucket = Counter()
    lineage_conf = Counter()  # adt_i|adt_j
    key_conf = Counter()
    coarse_conf = Counter()
    celltype_same = 0
    flag_counts = Counter()

    for r, p in enumerate(pairs):
        i, j = int(p["i"]), int(p["j"])
        d = _protein_deltas(adt, i, j)
        prot_abs[r] = d
        flags = p.get("flags") or {}
        for k, v in flags.items():
            if v:
                flag_counts[k] += 1
        coarse_same = p.get("coarse_i") == p.get("coarse_j")
        key_same = p.get("key_marker_lin_i") == p.get("key_marker_lin_j")
        bucket[classify_pair(flags, coarse_same, key_same)] += 1
        if p.get("cell_type_i") == p.get("cell_type_j"):
            celltype_same += 1
        a, b = p.get("adt_lineage_i"), p.get("adt_lineage_j")
        if a and b and a != b:
            lineage_conf[tuple(sorted([a, b]))] += 1
        ka, kb = p.get("key_marker_lin_i"), p.get("key_marker_lin_j")
        if ka and kb and ka != kb:
            key_conf[tuple(sorted([ka, kb]))] += 1
        ca, cb = p.get("coarse_i"), p.get("coarse_j")
        if ca and cb and ca != cb:
            coarse_conf[tuple(sorted([ca, cb]))] += 1

    mean_d = prot_abs.mean(axis=0)
    med_d = np.median(prot_abs, axis=0)
    q75_d = np.quantile(prot_abs, 0.75, axis=0)
    protein_rank = sorted(
        [
            {
                "protein": PANEL[k],
                "lineage": next(lin for lin, ps in LINEAGE_PANELS.items() if PANEL[k] in ps),
                "is_key_marker": PANEL[k] in KEY_MARKERS.values(),
                "mean_abs_delta": float(mean_d[k]),
                "median_abs_delta": float(med_d[k]),
                "q75_abs_delta": float(q75_d[k]),
            }
            for k in range(len(PANEL))
        ],
        key=lambda o: -o["mean_abs_delta"],
    )

    # Bio vs noise-like enrichment
    n = len(pairs)
    bio_keys = {"bio_adt_lineage", "bio_coarse_cell_type", "bio_key_marker"}
    n_bio = sum(bucket[k] for k in bio_keys)
    n_noise = bucket["noise_like_same_type_protein_drift"]
    n_fine = bucket["fine_protein_profile"]

    # Signal score: fraction of must-pairs that are lineage/coarse bio disagreements
    # (stronger bio signal) vs same-type protein drift (noise-like)
    bio_frac = n_bio / n
    noise_frac = n_noise / n
    signal_score = float(bio_frac - noise_frac)  # can be negative

    # Top protein concentration: how much of panel L1 is carried by top-1 / top-3 proteins
    row_sums = prot_abs.sum(axis=1) + 1e-12
    top1_share = float(np.mean(prot_abs.max(axis=1) / row_sums))
    top3_share = float(
        np.mean(np.sort(prot_abs, axis=1)[:, -3:].sum(axis=1) / row_sums)
    )

    results = {
        "claim_boundary": {
            "mode": "B_must_separate_autopsy",
            "not_claimed": [
                "Mode_A_residual_stream_as_field",
                "layer_Jacobian",
                "in_silico_gene_perturbs",
                "perturb_response",
            ],
            "note": "Decomposes existing must-separate pairs; does not claim causal gene/protein response.",
        },
        "inputs": {
            "pairs_path": str(pairs_path),
            "n_pairs": n,
            "adt_shape": list(adt.shape),
            "panel": PANEL,
        },
        "flag_counts": dict(flag_counts),
        "bucket_counts": dict(bucket),
        "bucket_fractions": {k: v / n for k, v in bucket.items()},
        "cell_type_identical_frac": celltype_same / n,
        "lineage_confusion_pairs": {
            f"{a}|{b}": c for (a, b), c in lineage_conf.most_common()
        },
        "key_marker_confusion_pairs": {
            f"{a}|{b}": c for (a, b), c in key_conf.most_common()
        },
        "coarse_confusion_pairs": {
            f"{a}|{b}": c for (a, b), c in coarse_conf.most_common()
        },
        "protein_rank_by_mean_abs_delta": protein_rank,
        "concentration": {
            "mean_top1_share_of_panel_l1": top1_share,
            "mean_top3_share_of_panel_l1": top3_share,
        },
        "signal": {
            "bio_frac": bio_frac,
            "noise_like_frac": noise_frac,
            "fine_protein_frac": n_fine / n,
            "bio_minus_noise": signal_score,
            "primary_metric": "bio_minus_noise",
            "interpretation": (
                "Positive bio_minus_noise => must-pairs dominated by lineage/coarse disagreements "
                "(bio signal in evidence geometry). Near-zero/negative => dominated by same-type "
                "protein drift (noise-like / fine-grained insufficiency)."
            ),
        },
        "signal_refined": {
            "adt_lineage_disagree_frac": flag_counts.get("diff_adt_lineage", 0) / n,
            "key_marker_disagree_frac": flag_counts.get("diff_key_marker_lineage", 0) / n,
            "coarse_disagree_frac": flag_counts.get("diff_coarse_cell_type", 0) / n,
            "same_celltype_must_frac": celltype_same / n,
            "primary_metric": "max(adt_lineage_disagree_frac, same_celltype_must_frac)",
            "primary_value": max(
                flag_counts.get("diff_adt_lineage", 0) / n,
                celltype_same / n,
            ),
            "note": (
                "Written must-pairs often saturate hierarchical bio_minus_noise; "
                "prefer ADT-lineage frac vs same-celltype must frac."
            ),
        },
    }

    stats_path = out_dir / "autopsy_stats.json"
    stats_path.write_text(json.dumps(results, indent=2))

    md = []
    md.append("# MUST-SEPARATE autopsy (Mode B)\n")
    md.append("**Not Mode A.** Decomposes reverse_step1 near-z pairs by protein/lineage.\n")
    md.append(f"- pairs analyzed: **{n}**\n")
    md.append(f"- bio_frac: **{bio_frac:.4f}** | noise_like_frac: **{noise_frac:.4f}** | "
              f"bio_minus_noise: **{signal_score:.4f}**\n")
    md.append(f"- identical fine cell_type frac: **{celltype_same/n:.4f}**\n")
    md.append("\n## Buckets\n")
    for k, v in bucket.most_common():
        md.append(f"- `{k}`: {v} ({v/n:.3f})\n")
    md.append("\n## Protein rank (mean |Δ| among must pairs)\n")
    md.append("| protein | lineage | key? | mean|Δ| | median|Δ| |\n|---|---|---|---:|---:|\n")
    for r in protein_rank:
        md.append(
            f"| {r['protein']} | {r['lineage']} | {r['is_key_marker']} | "
            f"{r['mean_abs_delta']:.4f} | {r['median_abs_delta']:.4f} |\n"
        )
    md.append("\n## Lineage confusion (ADT labels)\n")
    for k, v in list(results["lineage_confusion_pairs"].items())[:12]:
        md.append(f"- {k}: {v}\n")
    md.append("\nDo not claim Pearson > ~0.61. Not clinical.\n")
    (out_dir / "AUTOPSY.md").write_text("".join(md))
    print(json.dumps({"wrote": str(stats_path), "signal": results["signal"]}, indent=2))


if __name__ == "__main__":
    main()
