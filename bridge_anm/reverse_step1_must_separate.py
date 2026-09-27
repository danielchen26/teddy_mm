#!/usr/bin/env python3
"""Reverse step 1 — must-separate pairs in TEDDY evidence space (Mode B only).

Reads bridge export ``z_rna_export.npy`` + true ADT panel / cell meta, finds
near-identical z neighbors that disagree on lineage or protein labels.

Purpose
-------
Sufficiency-for-readout diagnostic: if two cells sit close in z but carry
different true lineage / protein profiles, a downstream readout *must separate*
them somehow (or abstain). This does **not** measure TEDDY representation
response to gene/protein perturbs (Mode A / reverse closed loop) — do not claim
Mode A, residual-as-field, Jacobian, gated ±ε G1–G4, fusion audit, or ATAC.

Notes
-----
``z_rna_export.npy`` is the compact sidecar (L2-normalized full z_512, then
``z_keep`` leading dims; default 32). Full loading factor remains: mean-pool
last-layer tokens at context length 1024.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.lineage_panels import KEY_MARKERS, LINEAGE_PANELS, all_panel_proteins

PANEL = all_panel_proteins()


def _coarse_from_cell_type(ct: str) -> str:
    s = (ct or "").lower()
    if any(
        k in s
        for k in (
            "b1 b",
            "naive b",
            "transitional b",
            "plasma",
            "igd+",
            "igd-",
            "b igk",
        )
    ) or (s.startswith("b ") or " b " in f" {s} "):
        if any(k in s for k in ("prog", "hsc", "gmp", "lmpp")):
            return "progenitor_other"
        return "b_lineage"
    if any(k in s for k in ("mono", " dc", "dc ", "pdc", "cd14", "cd16+ mono", "myeloid")):
        return "myeloid"
    if "nk" in s or "ilc" in s:
        return "nk_ilc"
    if any(k in s for k in ("cd4", "cd8", "treg", "gdt", "gd t", "t activated", "t naive", "t cd")):
        return "t_lineage"
    if any(k in s for k in ("prog", "hsc", "gmp", "lmpp", "normoblast", "eryth", "mk/e", "platelet")):
        return "progenitor_other"
    return "other"


def _adt_lineage_scores(adt_row: np.ndarray, name_to_j: dict[str, int]) -> dict[str, float]:
    out = {}
    for lin, prots in LINEAGE_PANELS.items():
        vals = [float(adt_row[name_to_j[p]]) for p in prots]
        out[lin] = float(np.mean(vals))
    return out


def _adt_label(scores: dict[str, float], margin: float) -> tuple[str, float]:
    ordered = sorted(scores.items(), key=lambda kv: -kv[1])
    top, second = ordered[0], ordered[1]
    gap = float(top[1] - second[1])
    if gap < margin:
        return "abstain", gap
    return top[0], gap


def _key_marker_argmax(adt_row: np.ndarray, name_to_j: dict[str, int]) -> str:
    best_lin, best_v = None, -np.inf
    for lin, marker in KEY_MARKERS.items():
        v = float(adt_row[name_to_j[marker]])
        if v > best_v:
            best_v, best_lin = v, lin
    return str(best_lin)


def knn_cosine(z: np.ndarray, k: int, batch: int) -> tuple[np.ndarray, np.ndarray]:
    norms = np.linalg.norm(z, axis=1, keepdims=True) + 1e-8
    zn = (z / norms).astype(np.float32)
    n = zn.shape[0]
    idx = np.zeros((n, k), dtype=np.int64)
    sims = np.zeros((n, k), dtype=np.float32)
    k_eff = min(k, n - 1)
    for s in range(0, n, batch):
        e = min(s + batch, n)
        S = zn[s:e] @ zn.T
        for i in range(e - s):
            S[i, s + i] = -np.inf
        part = np.argpartition(-S, kth=k_eff - 1, axis=1)[:, :k_eff]
        row = np.take_along_axis(S, part, axis=1)
        order = np.argsort(-row, axis=1)
        idx[s:e, :k_eff] = np.take_along_axis(part, order, axis=1)
        sims[s:e, :k_eff] = np.take_along_axis(row, order, axis=1)
    return idx, sims


def load_site4_meta(meta_path: Path, n_site4: int) -> list[dict]:
    rows = []
    with meta_path.open() as f:
        for line in f:
            o = json.loads(line)
            if o.get("slice") != "site4_test":
                continue
            rows.append(o)
            if len(rows) >= n_site4:
                break
    if len(rows) != n_site4:
        raise SystemExit(f"expected {n_site4} site4_test meta rows, got {len(rows)}")
    return rows


def _pair_flags(
    i: int,
    j: int,
    adt: np.ndarray,
    adt_lin_a: np.ndarray,
    key_lin_a: np.ndarray,
    coarse_a: np.ndarray,
    panel_scale: float,
    protein_l1_rel: float,
) -> dict:
    diff_adt = (
        adt_lin_a[i] != adt_lin_a[j]
        and adt_lin_a[i] != "abstain"
        and adt_lin_a[j] != "abstain"
    )
    diff_key = key_lin_a[i] != key_lin_a[j]
    diff_coarse = (
        coarse_a[i] != coarse_a[j]
        and coarse_a[i] not in ("other", "progenitor_other")
        and coarse_a[j] not in ("other", "progenitor_other")
    )
    l1 = float(np.mean(np.abs(adt[i] - adt[j])))
    rel_l1 = l1 / panel_scale
    diff_prot = rel_l1 >= protein_l1_rel
    must = bool(diff_adt or (diff_key and diff_prot) or (diff_coarse and diff_prot))
    return {
        "diff_adt": bool(diff_adt),
        "diff_key": bool(diff_key),
        "diff_coarse": bool(diff_coarse),
        "diff_prot": bool(diff_prot),
        "must": must,
        "l1": l1,
        "rel_l1": rel_l1,
    }


def sweep_threshold(
    thr: float,
    nn_idx: np.ndarray,
    nn_sim: np.ndarray,
    adt: np.ndarray,
    adt_lin_a: np.ndarray,
    key_lin_a: np.ndarray,
    coarse_a: np.ndarray,
    panel_scale: float,
    protein_l1_rel: float,
) -> dict:
    n = adt.shape[0]
    k = nn_idx.shape[1]
    seen = set()
    n_near_u = 0
    n_must = 0
    n_diff_adt = n_diff_key = n_diff_coarse = n_diff_prot = 0
    cells_m = np.zeros(n, dtype=bool)
    for i in range(n):
        for t in range(k):
            j = int(nn_idx[i, t])
            sim = float(nn_sim[i, t])
            if sim < thr:
                continue
            a, b = (i, j) if i < j else (j, i)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            n_near_u += 1
            fl = _pair_flags(
                i, j, adt, adt_lin_a, key_lin_a, coarse_a, panel_scale, protein_l1_rel
            )
            if fl["diff_adt"]:
                n_diff_adt += 1
            if fl["diff_key"]:
                n_diff_key += 1
            if fl["diff_coarse"]:
                n_diff_coarse += 1
            if fl["diff_prot"]:
                n_diff_prot += 1
            if fl["must"]:
                n_must += 1
                cells_m[i] = True
                cells_m[j] = True
    return {
        "min_cos": thr,
        "n_undirected_near_pairs": n_near_u,
        "n_must_separate": n_must,
        "frac_must_among_near": (n_must / n_near_u) if n_near_u else None,
        "n_cells_with_ge1_must": int(cells_m.sum()),
        "n_diff_adt_lineage": n_diff_adt,
        "n_diff_key_marker": n_diff_key,
        "n_diff_coarse": n_diff_coarse,
        "n_diff_protein_profile": n_diff_prot,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bridge-dir", type=Path, default=ROOT / "outputs" / "anm_cite_bridge")
    p.add_argument("--k", type=int, default=30)
    p.add_argument("--min-cos", type=float, default=0.98)
    p.add_argument("--adt-margin", type=float, default=0.05)
    p.add_argument("--protein-l1-rel", type=float, default=0.35)
    p.add_argument("--batch", type=int, default=256)
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--max-pairs-write", type=int, default=5000)
    p.add_argument("--cos-sweep", type=str, default="0.995,0.99,0.98,0.97,0.95")
    args = p.parse_args()

    bridge = args.bridge_dir
    out_dir = args.out_dir or (bridge / "reverse_step1")
    out_dir.mkdir(parents=True, exist_ok=True)

    z_all = np.load(bridge / "z_rna_export.npy").astype(np.float32)
    adt = np.load(bridge / "adt_true_panel.npy").astype(np.float32)
    n_site4 = adt.shape[0]
    if z_all.shape[0] < n_site4:
        raise SystemExit(f"z rows {z_all.shape[0]} < adt rows {n_site4}")
    z = z_all[:n_site4]
    meta = load_site4_meta(bridge / "cite_cells_meta.jsonl", n_site4)
    name_to_j = {name: i for i, name in enumerate(PANEL)}
    if adt.shape[1] != len(PANEL):
        raise SystemExit(f"adt panel width {adt.shape[1]} != {len(PANEL)}")

    adt_lin, key_lin, coarse, cell_types, cell_ids = [], [], [], [], []
    for i, m in enumerate(meta):
        scores = _adt_lineage_scores(adt[i], name_to_j)
        lab, _gap = _adt_label(scores, args.adt_margin)
        adt_lin.append(lab)
        key_lin.append(_key_marker_argmax(adt[i], name_to_j))
        ct = str(m.get("cell_type", ""))
        cell_types.append(ct)
        coarse.append(_coarse_from_cell_type(ct))
        cell_ids.append(str(m["cell_id"]))
    adt_lin_a = np.array(adt_lin)
    key_lin_a = np.array(key_lin)
    coarse_a = np.array(coarse)

    print(f"site4 n={n_site4} z_export_dim={z.shape[1]} k={args.k} min_cos={args.min_cos}")
    print("adt_lineage counts:", {k: int((adt_lin_a == k).sum()) for k in sorted(set(adt_lin))})
    print("coarse counts:", {k: int((coarse_a == k).sum()) for k in sorted(set(coarse))})

    nn_idx, nn_sim = knn_cosine(z, k=args.k, batch=args.batch)
    panel_scale = float(np.mean(np.abs(adt)) + 1e-6)

    pairs = []
    seen = set()
    n_near = 0
    n_diff_adt_lin = n_diff_key = n_diff_coarse = n_diff_protein = 0
    n_must = 0
    cells_with_must = np.zeros(n_site4, dtype=bool)

    for i in range(n_site4):
        for t in range(args.k):
            j = int(nn_idx[i, t])
            sim = float(nn_sim[i, t])
            if sim < args.min_cos:
                continue
            n_near += 1
            a, b = (i, j) if i < j else (j, i)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            fl = _pair_flags(
                i, j, adt, adt_lin_a, key_lin_a, coarse_a, panel_scale, args.protein_l1_rel
            )
            if fl["diff_adt"]:
                n_diff_adt_lin += 1
            if fl["diff_key"]:
                n_diff_key += 1
            if fl["diff_coarse"]:
                n_diff_coarse += 1
            if fl["diff_prot"]:
                n_diff_protein += 1
            if not fl["must"]:
                continue
            n_must += 1
            cells_with_must[i] = True
            cells_with_must[j] = True
            if len(pairs) < args.max_pairs_write:
                pairs.append(
                    {
                        "i": i,
                        "j": j,
                        "cell_id_i": cell_ids[i],
                        "cell_id_j": cell_ids[j],
                        "cos_sim": sim,
                        "l2_export": float(np.linalg.norm(z[i] - z[j])),
                        "adt_lineage_i": adt_lin_a[i],
                        "adt_lineage_j": adt_lin_a[j],
                        "key_marker_lin_i": key_lin_a[i],
                        "key_marker_lin_j": key_lin_a[j],
                        "coarse_i": coarse_a[i],
                        "coarse_j": coarse_a[j],
                        "cell_type_i": cell_types[i],
                        "cell_type_j": cell_types[j],
                        "adt_panel_l1_mean": fl["l1"],
                        "adt_panel_rel_l1": fl["rel_l1"],
                        "flags": {
                            "diff_adt_lineage": fl["diff_adt"],
                            "diff_key_marker_lineage": fl["diff_key"],
                            "diff_coarse_cell_type": fl["diff_coarse"],
                            "diff_protein_profile": fl["diff_prot"],
                        },
                    }
                )

    # `seen` already holds undirected unique near pairs at --min-cos
    n_near_unique = len(seen)

    must_sims = [x["cos_sim"] for x in pairs]
    must_l1 = [x["adt_panel_rel_l1"] for x in pairs]

    sweep_rows = [
        sweep_threshold(
            float(x),
            nn_idx,
            nn_sim,
            adt,
            adt_lin_a,
            key_lin_a,
            coarse_a,
            panel_scale,
            args.protein_l1_rel,
        )
        for x in args.cos_sweep.split(",")
        if x.strip()
    ]

    stats = {
        "claim_boundary": {
            "mode": "B_sufficiency_for_readout_probe",
            "not_claimed": [
                "Mode_A_residual_stream_as_field",
                "layer_Jacobian",
                "in_silico_gene_perturbs",
                "gated_pm_eps_G1_G4",
                "reverse_closed_loop_perturb_response",
                "fusion_audit",
                "ATAC_chromatin",
            ],
            "distinction": (
                "sufficiency-for-readout (near-identical z, different true labels) "
                "≠ perturb-response (how z moves under interventions)"
            ),
            "z_512_definition": (
                "mean-pool last-layer tokens at context length 1024 "
                "(not TEDDY pretrain 2048, not disease token); "
                "this script uses export sidecar z_keep dims after L2-normalize"
            ),
        },
        "inputs": {
            "z_path": str(bridge / "z_rna_export.npy"),
            "adt_true_panel": str(bridge / "adt_true_panel.npy"),
            "meta": str(bridge / "cite_cells_meta.jsonl"),
            "n_site4": n_site4,
            "z_export_dim": int(z.shape[1]),
            "z_all_rows": int(z_all.shape[0]),
            "panel": PANEL,
        },
        "params": {
            "k": args.k,
            "min_cos": args.min_cos,
            "adt_margin": args.adt_margin,
            "protein_l1_rel": args.protein_l1_rel,
            "batch": args.batch,
            "seed": args.seed,
            "max_pairs_write": args.max_pairs_write,
            "cos_sweep": args.cos_sweep,
        },
        "label_counts": {
            "adt_lineage": {k: int((adt_lin_a == k).sum()) for k in sorted(set(adt_lin))},
            "key_marker_lineage": {k: int((key_lin_a == k).sum()) for k in sorted(set(key_lin))},
            "coarse_cell_type": {k: int((coarse_a == k).sum()) for k in sorted(set(coarse))},
        },
        "results": {
            "n_directed_near_edges": n_near,
            "n_undirected_near_pairs_approx": n_near_unique,
            "n_must_separate_pairs_written": len(pairs),
            "n_must_separate_pairs_total_undirected": n_must,
            "n_cells_with_ge1_must_neighbor": int(cells_with_must.sum()),
            "frac_cells_with_ge1_must_neighbor": float(cells_with_must.mean()),
            "n_undirected_diff_adt_lineage": n_diff_adt_lin,
            "n_undirected_diff_key_marker": n_diff_key,
            "n_undirected_diff_coarse": n_diff_coarse,
            "n_undirected_diff_protein_profile": n_diff_protein,
            "must_cos_sim_quantiles": (
                {
                    "q25": float(np.quantile(must_sims, 0.25)),
                    "q50": float(np.quantile(must_sims, 0.50)),
                    "q75": float(np.quantile(must_sims, 0.75)),
                    "min": float(np.min(must_sims)),
                    "max": float(np.max(must_sims)),
                }
                if must_sims
                else None
            ),
            "must_rel_l1_quantiles": (
                {
                    "q25": float(np.quantile(must_l1, 0.25)),
                    "q50": float(np.quantile(must_l1, 0.50)),
                    "q75": float(np.quantile(must_l1, 0.75)),
                }
                if must_l1
                else None
            ),
        },
        "cos_sweep": sweep_rows,
    }

    pairs_path = out_dir / "must_separate_pairs.jsonl"
    with pairs_path.open("w") as f:
        for rec in pairs:
            f.write(json.dumps(rec) + "\n")

    stats_path = out_dir / "reverse_step1_stats.json"
    stats_path.write_text(json.dumps(stats, indent=2) + "\n")

    r = stats["results"]
    sweep_lines = [
        "| min_cos | near pairs | must-separate | frac must | cells w/ must | diff ADT lin | diff key | diff coarse | diff protein |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sweep_rows:
        frac = row["frac_must_among_near"]
        frac_s = f"{frac:.4f}" if frac is not None else "nan"
        sweep_lines.append(
            f"| {row['min_cos']} | {row['n_undirected_near_pairs']} | {row['n_must_separate']} | "
            f"{frac_s} | {row['n_cells_with_ge1_must']} | {row['n_diff_adt_lineage']} | "
            f"{row['n_diff_key_marker']} | {row['n_diff_coarse']} | {row['n_diff_protein_profile']} |"
        )
    sweep_table = "\n".join(sweep_lines)

    md = out_dir / "REVERSE_STEP1.md"
    md.write_text(
        f"""# Reverse step 1 — must-separate pairs (Mode B)

**Not Mode A.** This probe asks whether nearby points in the TEDDY evidence
export can carry different true lineage / protein labels
(**sufficiency-for-readout**). It does **not** measure how `z` responds to
gene or protein perturbs (**perturb-response** / reverse closed loop).

## z definition

- Full loading factor: `z_512` = mean-pool last-layer tokens @ context **1024**
  (not pretrain 2048, not disease token).
- This run used export sidecar `{bridge / "z_rna_export.npy"}`:
  site4 n={n_site4}, dim={z.shape[1]} (compact `z_keep`, not full 512).

## Params (primary pair list)

- k={args.k}, min_cos={args.min_cos}, adt_margin={args.adt_margin},
  protein_l1_rel={args.protein_l1_rel}

## Preliminary stats (site4/test n={n_site4}, primary min_cos={args.min_cos})

| Metric | Value |
|---|---:|
| Undirected near pairs (i&lt;j, cos≥min) | {r['n_undirected_near_pairs_approx']} |
| Must-separate pairs (undirected) | {r['n_must_separate_pairs_total_undirected']} |
| Pairs written | {r['n_must_separate_pairs_written']} |
| Cells with ≥1 must neighbor | {r['n_cells_with_ge1_must_neighbor']} ({r['frac_cells_with_ge1_must_neighbor']:.4f}) |
| Diff ADT lineage (among unique near) | {r['n_undirected_diff_adt_lineage']} |
| Diff key-marker lineage | {r['n_undirected_diff_key_marker']} |
| Diff coarse cell_type | {r['n_undirected_diff_coarse']} |
| Diff protein profile (rel L1) | {r['n_undirected_diff_protein_profile']} |

Must-pair cos quantiles: `{r['must_cos_sim_quantiles']}`

## Cosine threshold sweep (same kNN)

{sweep_table}

Outputs: `{pairs_path.name}`, `{stats_path.name}`.

Do not claim Pearson > ~0.61. Not clinical.
"""
    )

    print(json.dumps({"results": r, "cos_sweep": sweep_rows}, indent=2))
    print(f"wrote {pairs_path}")
    print(f"wrote {stats_path}")
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
