#!/usr/bin/env python3
"""Export typed events under modality_mask ∈ {rna_only, adt_only, joint}.

Preference order:
  1) Phase-2 BidirectionalCite checkpoint (outputs/cite_phase2/best.pt) —
     three-mode predict_adt on the SAME site4/test cells.
  2) Fallback: simulate masks on phase-1 MLP predictions + ADT-channel
     identity/prior stand-in (documented in manifest). Phase-2 Pearson is
     expected weak (~0.25–0.28 historically); do NOT claim win over phase-1 0.61.

Size factor: every mask decodes with the train-median constant by default
(--size-factor for rna_only, --size-factor-observed for adt_only/joint);
``--size-factor-observed measured`` restores the legacy per-cell measured depth
for the ADT-observed masks. The phase-2 ADT input transform comes from the
checkpoint (legacy checkpoints: clr). The phase-2 vs hybrid choice is gated on
the checkpoint's recorded validation Pearson by default (--phase2-gate val);
--phase2-gate test is the legacy gate on the exported test cells.

adt_only / joint read the exported cell's measured ADT as input; the same
values are the holdout answer key (manifest: answer_key_is_encoder_input).

Writes under outputs/anm_cite_bridge/missing_modality/.
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

from teddy_mm.bidirectional import load_bidirectional
from teddy_mm.device import resolve_device
from teddy_mm.models import AdtDecoder, MLP

from lib.lineage_panels import (
    KEY_MARKERS,
    LINEAGE_PANELS,
    MODALITY_FOR_LINEAGE,
    all_panel_proteins,
    protein_to_lineage,
)
from export_cite_events import (
    EVENT_TIMINGS,
    SIZE_FACTORS,
    event_time,
    git_commit,
    resolve_size_factor,
)

MASKS = ("rna_only", "adt_only", "joint")


def _pearson_mean(a: np.ndarray, b: np.ndarray) -> float:
    scores = []
    for j in range(a.shape[1]):
        if a[:, j].std() < 1e-8 or b[:, j].std() < 1e-8:
            continue
        scores.append(float(np.corrcoef(a[:, j], b[:, j])[0, 1]))
    return float(np.mean(scores)) if scores else float("nan")


def _infer_phase1(mlp, dec, z, sf, device, batch_size):
    zt = torch.from_numpy(z).to(device)
    sft = torch.from_numpy(sf.astype(np.float32)).to(device)
    chunks = []
    with torch.no_grad():
        for s in range(0, z.shape[0], batch_size):
            e = min(s + batch_size, z.shape[0])
            chunks.append(dec(mlp(zt[s:e]), sft[s:e])[0].cpu().numpy())
    return np.concatenate(chunks, axis=0)


def _infer_phase2(model, z, adt, sf, mode, device, batch_size, fm_steps, use_fm=True):
    zt = torch.from_numpy(z).to(device)
    yt = torch.from_numpy(adt).to(device)
    sft = torch.from_numpy(sf.astype(np.float32)).to(device)
    chunks = []
    with torch.no_grad():
        for s in range(0, z.shape[0], batch_size):
            e = min(s + batch_size, z.shape[0])
            pred = model.predict_adt(
                zt[s:e],
                yt[s:e],
                sft[s:e],
                mode=mode,
                fm_steps=fm_steps,
                use_fm=use_fm,
            )
            chunks.append(pred.cpu().numpy())
    return np.concatenate(chunks, axis=0)


def _simulate_mask_panel(phase1_pred, adt_true, mask, rng):
    """Fallback when phase-2 weights unavailable/weak.

    rna_only : phase-1 RNA→ADT pred
    adt_only : ADT-channel stand-in = observed ADT (identity evidence when ADT
               measured; NOT used as verifier label here — stored separately)
    joint    : mean of rna pred + observed ADT (both channels present)
    """
    if mask == "rna_only":
        return phase1_pred.copy()
    if mask == "adt_only":
        return adt_true.copy()
    # joint
    return 0.5 * phase1_pred + 0.5 * adt_true


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--phase1-ckpt", type=Path, default=ROOT / "outputs/cite_phase1/best.pt")
    p.add_argument("--phase2-ckpt", type=Path, default=ROOT / "outputs/cite_phase2/best.pt")
    p.add_argument("--out-dir", type=Path, default=ROOT / "outputs/anm_cite_bridge/missing_modality")
    p.add_argument("--n-cells", type=int, default=0, help="0 = all site4/test")
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--device", default="auto")
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--fm-steps", type=int, default=12)
    p.add_argument(
        "--force-simulate",
        action="store_true",
        help="Ignore phase-2 ckpt; simulate masks on phase-1 + ADT",
    )
    p.add_argument(
        "--phase2-min-pearson",
        type=float,
        default=0.35,
        help="If phase-2 rna_only pearson below this, fall back to simulate",
    )
    p.add_argument(
        "--event-timing",
        choices=EVENT_TIMINGS,
        default="simultaneous",
        help="simultaneous: all panel events at t=0 (default); panel-order: legacy",
    )
    p.add_argument(
        "--size-factor",
        choices=SIZE_FACTORS,
        default="train-median",
        help="ADT size factor for the RNA-only path (phase-1 and phase-2 rna_only): "
        "train-median (default), one, or measured (legacy).",
    )
    p.add_argument(
        "--size-factor-observed",
        choices=SIZE_FACTORS,
        default="train-median",
        help="ADT size factor for the phase-2 adt_only/joint masks: train-median "
        "(default), one, or measured (legacy: the exported cells' own measured ADT "
        "depth, i.e. the answer key's library size)",
    )
    p.add_argument(
        "--phase2-decode",
        choices=("fm", "direct"),
        default="fm",
        help="how phase-2 predictions are decoded: fm (default, legacy) = one flow-matching sample per "
        "cell (unseeded, so it varies run to run); direct = decode the encoder latent without sampling",
    )
    p.add_argument(
        "--phase2-gate",
        choices=("val", "test"),
        default="val",
        help="Pearson compared with --phase2-min-pearson: val (default) = the "
        "checkpoint's recorded val_rna_only_pearson; test (legacy) = rna_only "
        "Pearson on the exported test cells",
    )
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    npz = np.load(args.processed / "cite_arrays.npz", allow_pickle=False)
    z = np.load(args.processed / "z_rna.npy").astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)
    split, sites = npz["split"], npz["sites"]
    mask = (split == "test") & (sites == "site4")
    idx_all = np.where(mask)[0]
    if idx_all.size == 0:
        idx_all = np.where(split == "test")[0]
        print("WARN no site4 test; using all test")

    rng = np.random.default_rng(args.seed)
    n_avail = int(idx_all.size)
    n = n_avail if args.n_cells <= 0 else min(args.n_cells, n_avail)
    pick = idx_all if n == n_avail else np.sort(rng.choice(idx_all, size=n, replace=False))
    print(f"site4/test n={n} / available={n_avail}")

    adt_names = [str(x) for x in npz["adt_names"]]
    name_to_j = {nm: j for j, nm in enumerate(adt_names)}
    panel = all_panel_proteins()
    for prot in panel:
        if prot not in name_to_j:
            raise SystemExit(f"missing protein {prot}")
    panel_idx = np.array([name_to_j[p] for p in panel], dtype=np.int64)

    train_idx = np.where(split == "train")[0]
    train_sub = rng.choice(train_idx, size=min(8000, train_idx.size), replace=False)
    train_adt = npz["adt"][train_sub]
    p95 = {}
    for prot in panel:
        j = name_to_j[prot]
        p95[prot] = float(np.percentile(train_adt[:, j], 95))
        if p95[prot] < 1e-4:
            p95[prot] = 1.0

    device = resolve_device(args.device)
    z_dim = int(z.shape[1])
    n_adt = len(adt_names)
    adt_pick = npz["adt"][pick].astype(np.float32)
    # Measured ADT depth of the exported cells = the answer key's library size.
    # Only read under the legacy --size-factor(-observed) measured.
    sf_rna, sf_info = resolve_size_factor(args.size_factor, npz["adt_size_factor"], split, pick)
    sf_obs, sf_obs_info = resolve_size_factor(args.size_factor_observed, npz["adt_size_factor"], split, pick)
    print(f"rna_only size_factor={sf_info}")
    print(f"adt_only/joint size_factor={sf_obs_info}")
    z_pick = z[pick]

    # Always compute phase-1 baseline (strong unidirectional ~0.61)
    mlp = MLP(z_dim, z_dim, hidden=512).to(device)
    dec = AdtDecoder(z_dim, n_adt, hidden=512).to(device)
    blob1 = torch.load(args.phase1_ckpt, map_location=device, weights_only=False)
    mlp.load_state_dict(blob1["mlp"])
    dec.load_state_dict(blob1["dec"])
    mlp.eval()
    dec.eval()
    phase1_pred = _infer_phase1(mlp, dec, z_pick, sf_rna, device, args.batch_size)
    phase1_pearson = _pearson_mean(adt_pick, phase1_pred)
    print(f"phase-1 rna_only pearson (full ADT)={phase1_pearson:.4f}")

    source = "simulate_phase1_plus_adt"
    phase2_metrics = {}
    phase2_info: dict = {}
    preds_by_mask: dict[str, np.ndarray] = {}

    use_phase2 = (not args.force_simulate) and args.phase2_ckpt.exists()
    if use_phase2:
        try:
            blob2 = torch.load(args.phase2_ckpt, map_location=device, weights_only=False)
            blob2.setdefault("z_dim", z_dim)
            blob2.setdefault("n_adt", n_adt)
            model, phase2_info = load_bidirectional(blob2, device)
            print(f"phase-2 model: {phase2_info}")
            ck_sf = phase2_info.get("train_median_size_factor")

            def _phase2_sf(mode_name: str, sf_arr: np.ndarray) -> np.ndarray:
                # train-median: prefer the constant over the cells the phase-2 model was trained on
                if mode_name == "train-median" and ck_sf is not None:
                    return np.full(n, float(ck_sf), dtype=np.float32)
                return sf_arr

            sf_p2 = {"rna_only": _phase2_sf(args.size_factor, sf_rna)}
            sf_p2["adt_only"] = sf_p2["joint"] = _phase2_sf(args.size_factor_observed, sf_obs)
            phase2_info["size_factor_by_mask"] = {
                "rna_only": args.size_factor, "adt_only": args.size_factor_observed, "joint": args.size_factor_observed,
            }
            phase2_info["phase2_decode"] = args.phase2_decode
            phase2_info["train_median_source"] = "checkpoint" if ck_sf is not None else "split=='train' median"
            for m in MASKS:
                preds_by_mask[m] = _infer_phase2(
                    model, z_pick, adt_pick, sf_p2[m],
                    m, device, args.batch_size, args.fm_steps,
                    use_fm=(args.phase2_decode == "fm"),
                )
                phase2_metrics[f"{m}_pearson"] = _pearson_mean(adt_pick, preds_by_mask[m])
                print(f"phase-2 {m} pearson={phase2_metrics[f'{m}_pearson']:.4f}")
            # If phase-2 rna_only is weak, prefer simulate for rna/joint honesty
            # but KEEP phase-2 adt_only (real ADT-channel path without GT identity).
            # Gate on validation (recorded in the checkpoint) unless --phase2-gate test.
            if args.phase2_gate == "val":
                rna_r = float((blob2.get("metrics") or {}).get("val_rna_only_pearson", float("nan")))
            else:
                rna_r = phase2_metrics.get("rna_only_pearson", float("nan"))
            phase2_info["gate"] = {"split": args.phase2_gate, "rna_only_pearson": rna_r,
                                   "min_pearson": args.phase2_min_pearson}
            if not (rna_r == rna_r) or rna_r < args.phase2_min_pearson:
                print(
                    f"phase-2 rna_only pearson={rna_r} < {args.phase2_min_pearson}; "
                    "hybrid: phase-1 for rna_only/joint RNA channel, phase-2 for adt_only"
                )
                source = "hybrid_phase1_rna_phase2_adt"
                preds_by_mask["rna_only"] = phase1_pred
                # joint: fuse phase-1 RNA pred with phase-2 adt_only recon
                preds_by_mask["joint"] = 0.5 * phase1_pred + 0.5 * preds_by_mask["adt_only"]
                phase2_metrics["rna_only_pearson_phase1_override"] = phase1_pearson
                phase2_metrics["joint_pearson_hybrid"] = _pearson_mean(
                    adt_pick, preds_by_mask["joint"]
                )
            else:
                source = "phase2_bidirectional"
        except Exception as e:
            print(f"phase-2 load/infer failed ({e}); falling back to simulate")
            use_phase2 = False

    if not use_phase2 or source == "simulate_phase1_plus_adt" and not preds_by_mask:
        source = "simulate_phase1_plus_adt"
        print("using simulate masks on phase-1 pred + observed ADT channel")
        for m in MASKS:
            preds_by_mask[m] = _simulate_mask_panel(phase1_pred, adt_pick, m, rng)
            phase2_metrics[f"sim_{m}_pearson"] = _pearson_mean(adt_pick, preds_by_mask[m])

    p2l = protein_to_lineage()
    donors = npz["donors"]
    cell_types = npz["cell_types"]

    events_path = args.out_dir / "missing_modality_events.jsonl"
    cells_path = args.out_dir / "missing_modality_cells.jsonl"
    n_events = 0
    n_cells_written = 0

    # Save compact panels per mask
    for m in MASKS:
        np.save(args.out_dir / f"adt_pred_panel_{m}.npy", preds_by_mask[m][:, panel_idx])
    np.save(args.out_dir / "adt_true_panel.npy", adt_pick[:, panel_idx])
    np.save(args.out_dir / "global_indices.npy", pick.astype(np.int64))

    with events_path.open("w") as fe, cells_path.open("w") as fc:
        for local_i, global_i in enumerate(pick):
            base_id = f"cite_site4_{int(global_i)}"
            true_panel = {prot: float(adt_pick[local_i, name_to_j[prot]]) for prot in panel}
            for m in MASKS:
                cell_id = f"{base_id}__{m}"
                pred = preds_by_mask[m]
                pred_panel = {prot: float(pred[local_i, name_to_j[prot]]) for prot in panel}
                cell_rec = {
                    "cell_id": cell_id,
                    "base_cell_id": base_id,
                    "global_index": int(global_i),
                    "modality_mask": m,
                    "split": str(split[global_i]),
                    "site": str(sites[global_i]),
                    "donor": str(donors[global_i]),
                    "cell_type": str(cell_types[global_i]),
                    "slice": "site4_test",
                    "adt_pred_panel": pred_panel,
                    "adt_true_panel_holdout": true_panel,
                    "holdout_label_policy": "posthoc_Q_f_only_do_not_fit_ANM",
                    "protein_names_panel": panel,
                    "source": source,
                }
                fc.write(json.dumps(cell_rec) + "\n")
                n_cells_written += 1
                for t_i, prot in enumerate(panel):
                    lin = p2l[prot]
                    raw = pred_panel[prot]
                    value = float(np.clip(raw / max(p95[prot], 1e-6), 0.0, 1.0))
                    # Source ablation semantics: tag channel; admission by mask
                    if m == "rna_only":
                        channel = "rna"
                    elif m == "adt_only":
                        channel = "adt"
                    else:
                        channel = "joint"
                    # Missing-modality penalty on declared workability:
                    # rna_only keeps full RNA evidence; adt_only has no RNA→protein
                    # typed path unless phase-2 ADT encoder; joint keeps both.
                    # Value already from the mask-specific prediction.
                    ev = {
                        "event_id": f"{cell_id}:{prot}",
                        "cell_id": cell_id,
                        "base_cell_id": base_id,
                        "time": event_time(args.event_timing, t_i),
                        "panel_order": t_i,
                        "protein": prot,
                        "action": lin,
                        "modality": MODALITY_FOR_LINEAGE[lin],
                        "polarity": "support",
                        "value": value,
                        "adt_pred_raw": raw,
                        "norm_p95_train": p95[prot],
                        "provenance": f"teddy_missing_modality:{source}:{m}",
                        "modality_mask": m,
                        "source_channel": channel,
                        "split": str(split[global_i]),
                        "site": str(sites[global_i]),
                        "slice": "site4_test",
                        "is_key_marker": prot == KEY_MARKERS[lin],
                        "holdout_adt_true": true_panel[prot],
                        "holdout_note": "verifier_only",
                    }
                    fe.write(json.dumps(ev) + "\n")
                    n_events += 1

    # Per-mask panel pearson (honest)
    panel_pearson = {}
    for m in MASKS:
        panel_pearson[m] = _pearson_mean(
            adt_pick[:, panel_idx], preds_by_mask[m][:, panel_idx]
        )

    manifest = {
        "n_base_cells_site4": int(n),
        "n_cells_written": n_cells_written,
        "n_events": n_events,
        "n_available_site4_test": n_avail,
        "masks": list(MASKS),
        "source": source,
        "phase1_full_adt_pearson": phase1_pearson,
        "phase1_quoted_unidirectional": 0.61,
        "phase2_metrics": phase2_metrics,
        "panel_pearson_by_mask": panel_pearson,
        "claim_boundary": (
            "missing_modality_ANM_demo_only; do_not_claim_win_over_phase1_0.61; "
            "not_clinical"
        ),
        "panel_proteins": panel,
        "lineage_panels": LINEAGE_PANELS,
        "key_markers": KEY_MARKERS,
        "p95_train": p95,
        "phase1_ckpt": str(args.phase1_ckpt),
        "phase2_ckpt": str(args.phase2_ckpt) if args.phase2_ckpt.exists() else None,
        "processed": str(args.processed),
        "seed": args.seed,
        "device": str(device),
        "events_path": str(events_path),
        "cells_path": str(cells_path),
        "same_cells_across_masks": True,
        "event_timing": args.event_timing,
        "size_factor_rna_only": sf_info,
        "size_factor_adt_observed_masks": sf_obs_info,
        "phase2_model": phase2_info or None,
        # adt_only/joint values are built from the exported cells' measured ADT in every
        # source (phase-2 input, hybrid via phase-2 adt_only, simulate by copy); the same
        # values are adt_true_panel_holdout.
        "answer_key_is_encoder_input": {"rna_only": False, "adt_only": True, "joint": True},
        "adt_only_is_answer_key_copy": source == "simulate_phase1_plus_adt",
        "flags": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
        "git_commit": git_commit(),
        "anm_semantics": (
            "missing modality = source ablation / channel-restricted typed evidence; "
            "P_f workability drops when key evidence absent; Q_f vs holdout GT only"
        ),
    }
    (args.out_dir / "export_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({
        "source": source,
        "n_base_cells": n,
        "n_cells_written": n_cells_written,
        "n_events": n_events,
        "panel_pearson_by_mask": panel_pearson,
        "phase1_pearson": phase1_pearson,
    }, indent=2))


if __name__ == "__main__":
    main()
