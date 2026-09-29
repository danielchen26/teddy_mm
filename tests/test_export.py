"""Tests for bridge_anm/export_cite_events.py event timing and size-factor flags."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
BRIDGE = ROOT / "bridge_anm"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(BRIDGE))

import export_cite_events as ex  # noqa: E402
from lib.lineage_panels import all_panel_proteins  # noqa: E402
from teddy_mm.models import AdtDecoder, MLP  # noqa: E402


# ---------------------------------------------------------------- unit tests
def test_event_time_default_is_zero():
    assert [ex.event_time("simultaneous", i) for i in range(9)] == [0] * 9
    assert [ex.event_time("panel-order", i) for i in range(4)] == [0, 1, 2, 3]
    with pytest.raises(ValueError):
        ex.event_time("bogus", 1)


def test_train_median_never_reads_exported_cells():
    split = np.array(["train"] * 5 + ["test"] * 3)
    sf = np.array([1.0, 2.0, 3.0, 4.0, 5.0, np.nan, np.nan, np.nan])
    pick = np.array([5, 6, 7])
    out, info = ex.resolve_size_factor("train-median", sf, split, pick)
    assert np.all(np.isfinite(out)) and np.allclose(out, 3.0)
    assert info["mode"] == "train-median" and info["value"] == 3.0


def test_train_median_excludes_exported_train_cells():
    split = np.array(["train"] * 6)
    sf = np.array([1.0, 2.0, 3.0, np.nan, 1e9, 1e9])
    out, info = ex.resolve_size_factor("train-median", sf, split, np.array([3, 4, 5]))
    assert np.allclose(out, 2.0)
    assert info["n_exported_excluded"] == 3


def test_one_and_measured_modes():
    split = np.array(["train", "test", "test"])
    sf = np.array([7.0, 2.0, 3.0], dtype=np.float32)
    one, _ = ex.resolve_size_factor("one", sf, split, np.array([1, 2]))
    meas, _ = ex.resolve_size_factor("measured", sf, split, np.array([1, 2]))
    assert np.allclose(one, 1.0)
    assert np.allclose(meas, [2.0, 3.0])


# ---------------------------------------------------------- end-to-end tests
def _make_fixture(tmp: Path, poison_exported_sf: bool = True) -> dict[str, Path]:
    rng = np.random.default_rng(0)
    panel = all_panel_proteins()
    adt_names = panel + ["EXTRA1", "EXTRA2"]
    n, z_dim = 60, 8
    split = np.array(["train"] * 30 + ["val"] * 10 + ["test"] * 20)
    sites = np.array(["site1"] * 40 + ["site4"] * 20)
    adt = rng.poisson(5.0, size=(n, len(adt_names))).astype(np.float32)
    sf = rng.uniform(0.5, 2.0, size=n).astype(np.float32)
    if poison_exported_sf:
        # Exported cells = site4/test and val/non-site4 (OOD). Any read of their
        # measured size factor propagates NaN into the predictions.
        sf[30:] = np.nan
    processed = tmp / "processed"
    processed.mkdir()
    np.savez(
        processed / "cite_arrays.npz",
        split=split,
        sites=sites,
        donors=np.array([f"d{i % 3}" for i in range(n)]),
        cell_types=np.array(["T"] * n),
        adt_names=np.array(adt_names),
        adt=adt,
        adt_size_factor=sf,
    )
    np.save(processed / "z_rna.npy", rng.normal(size=(n, z_dim)).astype(np.float32))
    phase1 = tmp / "phase1"
    phase1.mkdir()
    torch.manual_seed(0)
    mlp = MLP(z_dim, z_dim, hidden=512)
    dec = AdtDecoder(z_dim, len(adt_names), hidden=512)
    torch.save({"mlp": mlp.state_dict(), "dec": dec.state_dict()}, phase1 / "best.pt")
    (phase1 / "test_per_protein.json").write_text(
        json.dumps([{"protein": p, "pearson_mlp": 0.9, "pearson_fm": 0.9} for p in panel])
    )
    (phase1 / "metrics.json").write_text(json.dumps({"test": {"test_mlp_pearson": 0.61}}))
    return {"processed": processed, "ckpt": phase1 / "best.pt",
            "per_protein": phase1 / "test_per_protein.json", "out": tmp / "out"}


def _run(fx: dict[str, Path], *extra: str) -> tuple[list[dict], dict]:
    cmd = [
        sys.executable, str(BRIDGE / "export_cite_events.py"),
        "--processed", str(fx["processed"]), "--ckpt", str(fx["ckpt"]),
        "--per-protein", str(fx["per_protein"]), "--out-dir", str(fx["out"]),
        "--ood-n", "5", "--device", "cpu", *extra,
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    events = [json.loads(l) for l in (fx["out"] / "cite_typed_events.jsonl").open()]
    manifest = json.loads((fx["out"] / "export_manifest.json").read_text())
    return events, manifest


def test_default_export_simultaneous_and_train_median(tmp_path):
    fx = _make_fixture(tmp_path, poison_exported_sf=True)
    events, man = _run(fx)
    assert len(events) == (20 + 5) * len(all_panel_proteins())
    assert all(e["time"] == 0 for e in events)
    assert all(np.isfinite(e["adt_pred_raw"]) for e in events)
    assert np.all(np.isfinite(np.load(fx["out"] / "adt_pred_panel.npy")))
    assert np.all(np.isfinite(np.load(fx["out"] / "ood_adt_pred_panel.npy")))
    assert man["event_timing"] == "simultaneous"
    assert man["size_factor"]["mode"] == "train-median"
    sf_train = np.load(fx["processed"] / "cite_arrays.npz")["adt_size_factor"][:30]
    assert man["size_factor"]["value"] == pytest.approx(float(np.median(sf_train)))
    assert man["flags"]["event_timing"] == "simultaneous"
    assert man["flags"]["size_factor"] == "train-median"
    assert "git_commit" in man
    assert man["phase1_test_mlp_pearson_mean"] == pytest.approx(0.61)


def test_legacy_flags_reachable(tmp_path):
    fx = _make_fixture(tmp_path, poison_exported_sf=False)
    events, man = _run(fx, "--event-timing", "panel-order", "--size-factor", "measured")
    n_panel = len(all_panel_proteins())
    first = [e["time"] for e in events[:n_panel]]
    assert first == list(range(n_panel))
    assert man["size_factor"]["mode"] == "measured"
