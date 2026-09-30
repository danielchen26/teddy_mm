#!/usr/bin/env python3
"""Export TEDDY phase-1 CITE typed events (δu) for the ANM bridge.

Samples test/site4 cells (prefer full test when feasible), runs frozen MLP+decoder
from best.pt (no FM, no retrain), writes JSONL events + compact sidecar arrays.
Optionally exports an OOD slice from non-site4 (val preferred) for transfer check.
Holdout adt_true is stored for post-hoc Q_f only and marked as such.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from teddy_mm.device import resolve_device
from teddy_mm.models import AdtDecoder, MLP

from lib.lineage_panels import (
    CRITERIA,
    KEY_MARKERS,
    LINEAGE_PANELS,
    MODALITY_FOR_LINEAGE,
    all_panel_proteins,
    protein_to_lineage,
)


EVENT_TIMINGS = ("simultaneous", "panel-order")
SIZE_FACTORS = ("train-median", "one", "measured")


def event_time(timing: str, t_i: int) -> int:
    """Event time for the t_i-th panel protein.

    simultaneous (default): every panel event of a cell is observed at t=0 —
    the panel is one simultaneous readout, not a sequence.
    panel-order (legacy): t = position in the panel list.
    """
    if timing == "simultaneous":
        return 0
    if timing == "panel-order":
        return int(t_i)
    raise ValueError(f"unknown event timing {timing!r}")


def resolve_size_factor(mode: str, sf_all, split, pick) -> tuple[np.ndarray, dict]:
    """ADT size factor fed to the decoder for the exported (RNA-only) cells.

    train-median (default): one constant = median adt_size_factor over
        split=='train'; never reads the exported cells' measured ADT depth.
    one: constant 1.0.
    measured (legacy): the exported cells' own measured adt_size_factor
        (leaks ADT library size into an "RNA-only" prediction).
    """
    pick = np.asarray(pick, dtype=np.int64)
    n = int(pick.size)
    if mode == "train-median":
        train_idx = np.where(np.asarray(split) == "train")[0]
        # Exclude exported cells (only possible for a train-split OOD fallback)
        # so their own measured depth never enters the constant.
        n_excl = int(np.intersect1d(train_idx, pick).size)
        if n_excl:
            train_idx = np.setdiff1d(train_idx, pick)
        if train_idx.size == 0:
            raise SystemExit("--size-factor train-median: no split=='train' cells")
        val = float(np.median(np.asarray(sf_all[train_idx], dtype=np.float64)))
        return np.full(n, val, dtype=np.float32), {
            "mode": mode, "value": val, "n_train_cells": int(train_idx.size),
            "n_exported_excluded": n_excl,
        }
    if mode == "one":
        return np.ones(n, dtype=np.float32), {"mode": mode, "value": 1.0}
    if mode == "measured":
        return np.asarray(sf_all[pick], dtype=np.float32), {
            "mode": mode, "value": None,
            "note": "per-cell measured adt_size_factor of exported cells (legacy)",
        }
    raise ValueError(f"unknown size factor mode {mode!r}")


def git_commit(root: Path = ROOT) -> str | None:
    import subprocess

    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        return out + ("+dirty" if dirty else "")
    except Exception:
        return None


def _phase1_mean(path: Path) -> float | None:
    try:
        return float(json.loads(path.read_text())["test"]["test_mlp_pearson"])
    except (OSError, KeyError, ValueError):
        print(f"WARN could not read phase-1 metrics at {path}")
        return None


def _load_pearson(path: Path) -> dict[str, float]:
    rows = json.loads(path.read_text())
    return {r["protein"]: float(r["pearson_mlp"]) for r in rows}


def _verify_panels(pearson: dict[str, float], min_r: float = 0.85) -> dict[str, float]:
    used = {}
    missing = []
    weak = []
    for lin, prots in LINEAGE_PANELS.items():
        for p in prots:
            if p not in pearson:
                missing.append(p)
                continue
            r = pearson[p]
            used[p] = r
            if r < min_r:
                weak.append((p, r))
    if missing:
        raise SystemExit(f"panel proteins missing from test_per_protein.json: {missing}")
    if weak:
        print(f"WARN panel proteins below {min_r}: {weak}")
    return used


def _infer(mlp, dec, z, sf, device, batch_size):
    zt = torch.from_numpy(z).to(device)
    sft = torch.from_numpy(sf.astype(np.float32)).to(device)
    pred_chunks = []
    n = z.shape[0]
    with torch.no_grad():
        for s in range(0, n, batch_size):
            e = min(s + batch_size, n)
            pred_chunks.append(dec(mlp(zt[s:e]), sft[s:e])[0].cpu().numpy())
    return np.concatenate(pred_chunks, axis=0)


def _write_slice(
    *,
    fe,
    fc,
    pick,
    sites,
    donors,
    cell_types,
    split_arr,
    z,
    adt_pred,
    adt_true,
    name_to_j,
    panel,
    p95,
    used_r,
    p2l,
    ckpt_name,
    z_keep,
    slice_name,
    id_prefix,
    compact,
    store_z_rows,
    event_timing="simultaneous",
    size_factor_mode="train-median",
):
    n_events = 0
    for local_i, global_i in enumerate(pick):
        cell_id = f"{id_prefix}_{int(global_i)}"
        pred_panel = {prot: float(adt_pred[local_i, name_to_j[prot]]) for prot in panel}
        true_panel = {prot: float(adt_true[local_i, name_to_j[prot]]) for prot in panel}
        cell_rec = {
            "cell_id": cell_id,
            "global_index": int(global_i),
            "split": str(split_arr[global_i]),
            "site": str(sites[global_i]),
            "donor": str(donors[global_i]),
            "cell_type": str(cell_types[global_i]),
            "modality_mask": "rna_only",
            "slice": slice_name,
            "adt_pred_panel": pred_panel,
            "adt_true_panel_holdout": true_panel,
            "holdout_label_policy": "posthoc_Q_f_only_do_not_fit_ANM",
            "protein_names_panel": panel,
            "phase1_mlp_pearson": used_r,
        }
        if not compact:
            z_comp = z[global_i, :z_keep].astype(np.float32)
            cell_rec["z_rna_compressed"] = [float(x) for x in z_comp]
        else:
            store_z_rows.append(z[global_i, :z_keep].astype(np.float32))
            cell_rec["z_rna_row"] = len(store_z_rows) - 1
        fc.write(json.dumps(cell_rec) + "\n")

        for t_i, prot in enumerate(panel):
            lin = p2l[prot]
            raw = pred_panel[prot]
            value = float(np.clip(raw / max(p95[prot], 1e-6), 0.0, 1.0))
            ev = {
                "event_id": f"{cell_id}:{prot}",
                "cell_id": cell_id,
                "time": event_time(event_timing, t_i),
                "panel_order": t_i,
                "protein": prot,
                "action": lin,
                "modality": MODALITY_FOR_LINEAGE[lin],
                "polarity": "support",
                "value": value,
                "adt_pred_raw": raw,
                "norm_p95_train": p95[prot],
                "provenance": f"teddy_phase1_mlp:{ckpt_name}:rna_only",
                "size_factor_mode": size_factor_mode,
                "modality_mask": "rna_only",
                "split": str(split_arr[global_i]),
                "site": str(sites[global_i]),
                "slice": slice_name,
                "phase1_mlp_pearson": used_r[prot],
                "is_key_marker": prot == KEY_MARKERS[lin],
                "holdout_adt_true": true_panel[prot],
                "holdout_note": "verifier_only",
            }
            fe.write(json.dumps(ev) + "\n")
            n_events += 1
    return n_events


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/cite_phase1/best.pt")
    p.add_argument(
        "--per-protein",
        type=Path,
        default=ROOT / "outputs/cite_phase1/test_per_protein.json",
    )
    p.add_argument(
        "--metrics",
        type=Path,
        default=None,
        help="phase-1 metrics.json (default: next to --per-protein)",
    )
    p.add_argument("--out-dir", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    p.add_argument(
        "--event-timing",
        choices=EVENT_TIMINGS,
        default="simultaneous",
        help="simultaneous: all panel events at t=0 (default); "
        "panel-order: t = panel index (legacy)",
    )
    p.add_argument(
        "--size-factor",
        choices=SIZE_FACTORS,
        default="train-median",
        help="decoder ADT size factor: train-median (default, constant median over "
        "split=='train'), one, or measured (legacy per-cell measured ADT depth)",
    )
    p.add_argument(
        "--n-cells",
        type=int,
        default=0,
        help="site4/test cells to export; 0 = all feasible",
    )
    p.add_argument(
        "--ood-n",
        type=int,
        default=2000,
        help="OOD cells from non-site4 (val preferred); 0 = skip",
    )
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--device", default="auto")
    p.add_argument("--z-keep", type=int, default=32)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument(
        "--compact",
        action="store_true",
        default=True,
        help="omit z vectors from JSONL; store z_rna_export.npy (default True)",
    )
    p.add_argument("--no-compact", action="store_true")
    args = p.parse_args()
    compact = False if args.no_compact else args.compact

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pearson = _load_pearson(args.per_protein)
    used_r = _verify_panels(pearson)

    npz = np.load(args.processed / "cite_arrays.npz", allow_pickle=False)
    z = np.load(args.processed / "z_rna.npy").astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)

    split = npz["split"]
    sites = npz["sites"]
    mask = (split == "test") & (sites == "site4")
    idx_all = np.where(mask)[0]
    if idx_all.size == 0:
        mask = split == "test"
        idx_all = np.where(mask)[0]
        print("WARN no site4 test cells; falling back to all test")

    rng = np.random.default_rng(args.seed)
    n_avail = int(idx_all.size)
    n = n_avail if args.n_cells <= 0 else min(args.n_cells, n_avail)
    if n < n_avail:
        pick = rng.choice(idx_all, size=n, replace=False)
        pick.sort()
        capped = True
    else:
        pick = idx_all.copy()
        capped = False
    print(f"site4/test export: n={n} / available={n_avail} capped={capped}")

    ood_pick = np.array([], dtype=np.int64)
    ood_slice_name = None
    if args.ood_n > 0:
        ood_mask = (split == "val") & (sites != "site4")
        ood_idx = np.where(ood_mask)[0]
        ood_slice_name = "val_non_site4"
        if ood_idx.size == 0:
            ood_mask = (split == "train") & (sites != "site4")
            ood_idx = np.where(ood_mask)[0]
            ood_slice_name = "train_non_site4"
        n_ood = min(args.ood_n, int(ood_idx.size))
        ood_pick = rng.choice(ood_idx, size=n_ood, replace=False)
        ood_pick.sort()
        print(f"OOD export: n={n_ood} slice={ood_slice_name}")

    adt_names = [str(x) for x in npz["adt_names"]]
    name_to_j = {n: j for j, n in enumerate(adt_names)}
    panel = all_panel_proteins()
    for prot in panel:
        if prot not in name_to_j:
            raise SystemExit(f"protein {prot} not in adt_names")

    train_mask = split == "train"
    train_idx = np.where(train_mask)[0]
    train_sub = rng.choice(train_idx, size=min(8000, train_idx.size), replace=False)
    train_adt = npz["adt"][train_sub]
    p95 = {}
    for prot in panel:
        j = name_to_j[prot]
        p95[prot] = float(np.percentile(train_adt[:, j], 95))
        if p95[prot] < 1e-4:
            p95[prot] = 1.0

    device = resolve_device(args.device)
    print(f"device={device}")
    z_dim = z.shape[1]
    n_adt = len(adt_names)
    mlp = MLP(z_dim, z_dim, hidden=512).to(device)
    dec = AdtDecoder(z_dim, n_adt, hidden=512).to(device)
    blob = torch.load(args.ckpt, map_location=device, weights_only=False)
    mlp.load_state_dict(blob["mlp"])
    dec.load_state_dict(blob["dec"])
    mlp.eval()
    dec.eval()

    sf_all = npz["adt_size_factor"]
    sf_pick, sf_info = resolve_size_factor(args.size_factor, sf_all, split, pick)
    print(f"size_factor={sf_info}")
    adt_pred = _infer(mlp, dec, z[pick], sf_pick, device, args.batch_size)
    adt_true = npz["adt"][pick].astype(np.float32)
    if ood_pick.size:
        sf_ood, sf_ood_info = resolve_size_factor(args.size_factor, sf_all, split, ood_pick)
        ood_pred = _infer(mlp, dec, z[ood_pick], sf_ood, device, args.batch_size)
        ood_true = npz["adt"][ood_pick].astype(np.float32)
    else:
        ood_pred = ood_true = None

    p2l = protein_to_lineage()
    events_path = args.out_dir / "cite_typed_events.jsonl"
    cells_path = args.out_dir / "cite_cells_meta.jsonl"
    panel_idx = np.array([name_to_j[p] for p in panel], dtype=np.int64)
    np.save(args.out_dir / "adt_pred_panel.npy", adt_pred[:, panel_idx])
    np.save(args.out_dir / "adt_true_panel.npy", adt_true[:, panel_idx])
    if ood_pred is not None:
        np.save(args.out_dir / "ood_adt_pred_panel.npy", ood_pred[:, panel_idx])
        np.save(args.out_dir / "ood_adt_true_panel.npy", ood_true[:, panel_idx])

    store_z_rows: list = []
    n_events = 0
    with events_path.open("w") as fe, cells_path.open("w") as fc:
        n_events += _write_slice(
            fe=fe,
            fc=fc,
            pick=pick,
            sites=sites,
            donors=npz["donors"],
            cell_types=npz["cell_types"],
            split_arr=split,
            z=z,
            adt_pred=adt_pred,
            adt_true=adt_true,
            name_to_j=name_to_j,
            panel=panel,
            p95=p95,
            used_r=used_r,
            p2l=p2l,
            ckpt_name=args.ckpt.name,
            z_keep=args.z_keep,
            slice_name="site4_test",
            id_prefix="cite_site4",
            compact=compact,
            store_z_rows=store_z_rows,
            event_timing=args.event_timing,
            size_factor_mode=args.size_factor,
        )
        if ood_pick.size:
            n_events += _write_slice(
                fe=fe,
                fc=fc,
                pick=ood_pick,
                sites=sites,
                donors=npz["donors"],
                cell_types=npz["cell_types"],
                split_arr=split,
                z=z,
                adt_pred=ood_pred,
                adt_true=ood_true,
                name_to_j=name_to_j,
                panel=panel,
                p95=p95,
                used_r=used_r,
                p2l=p2l,
                ckpt_name=args.ckpt.name,
                z_keep=args.z_keep,
                slice_name=ood_slice_name,
                id_prefix="cite_ood",
                compact=compact,
                store_z_rows=store_z_rows,
                event_timing=args.event_timing,
                size_factor_mode=args.size_factor,
            )

    z_path = None
    z_512_path = None
    if compact and store_z_rows:
        z_path = args.out_dir / "z_rna_export.npy"
        np.save(z_path, np.stack(store_z_rows, axis=0))
        # Full L2-normalized z_512 in the same export row order (site4 then OOD).
        full_idx = np.concatenate([pick, ood_pick]) if ood_pick.size else pick
        z_512_path = args.out_dir / "z_rna_512.npy"
        np.save(z_512_path, z[full_idx].astype(np.float32))
        print(f"wrote {z_512_path} shape={(len(full_idx), z.shape[1])}")

    n_site4 = int(pick.size)
    n_ood = int(ood_pick.size)
    manifest = {
        "n_cells_site4_test": n_site4,
        "n_cells_ood": n_ood,
        "n_cells_total": n_site4 + n_ood,
        "n_cells": n_site4,
        "n_events": n_events,
        "n_available_site4_test": n_avail,
        "capped": capped,
        "split": "test",
        "sites": ["site4"],
        "ood_slice": ood_slice_name,
        "panel_proteins": panel,
        "lineage_panels": LINEAGE_PANELS,
        "key_markers": KEY_MARKERS,
        "phase1_mlp_pearson_quoted": used_r,
        "phase1_test_mlp_pearson_mean": _phase1_mean(
            args.metrics or args.per_protein.parent / "metrics.json"
        ),
        "event_timing": args.event_timing,
        "size_factor": sf_info,
        "size_factor_ood": sf_ood_info if ood_pick.size else None,
        "flags": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "git_commit": git_commit(),
        "processed": str(args.processed),
        "ckpt": str(args.ckpt),
        "seed": args.seed,
        "z_keep": args.z_keep,
        "compact": compact,
        "z_rna_export": str(z_path) if z_path else None,
        "z_rna_512": str(z_512_path) if z_512_path else None,
        "z_512_definition": (
            "L2-normalized mean-pool last-layer TEDDY-G @ ctx 1024 "
            "(from data/processed/cite/z_rna.npy); full d_model=512"
        ),
        "z_512_shape": (
            [int(n_site4 + n_ood), int(z.shape[1])] if z_512_path else None
        ),
        "device": str(device),
        "criteria_available": list(CRITERIA.keys()),
        "events_path": str(events_path),
        "cells_path": str(cells_path),
        "claim_boundary": "typed_evidence_export_only_not_clinical",
    }
    (args.out_dir / "export_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: manifest[k] for k in (
        "n_cells_site4_test", "n_cells_ood", "n_events", "capped", "device", "compact",
        "event_timing", "size_factor", "git_commit",
    )}, indent=2))


if __name__ == "__main__":
    main()
