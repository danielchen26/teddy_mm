"""Phase-2 (BidirectionalCite) fixes: ADT transform applied once, no test-time
input read from the measured test ADT, and legacy flags reproducing old outputs.

Legacy-reproduction tests compare against the pre-fix code itself, extracted with
``git archive`` from LEGACY_COMMIT (the commit this fix branched from); they skip
when that commit is not reachable.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import io
import tarfile
from pathlib import Path

import numpy as np
import pytest
import torch
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
BRIDGE = ROOT / "bridge_anm"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(BRIDGE))

from teddy_mm import bidirectional as bd  # noqa: E402
from teddy_mm.bidirectional import (  # noqa: E402
    AdtEncoder,
    BidirectionalCite,
    load_bidirectional,
    resolve_adt_input_transform,
)
from teddy_mm.data import eval_size_factor  # noqa: E402
from teddy_mm.modality_dropout import clr_log_adt  # noqa: E402

LEGACY_COMMIT = "854cac6"
Z_DIM, N_ADT = 16, 12
ENV = {**os.environ, "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def legacy_tree(tmp_path_factory) -> Path:
    """Pre-fix copy of teddy_mm/, the phase-2 scripts and the missing-modality export.

    From ``git archive LEGACY_COMMIT``; if git is unavailable, from the directory in
    $PHASE2_LEGACY_TREE (an extracted copy of that commit).
    """
    env_tree = os.environ.get("PHASE2_LEGACY_TREE")
    try:
        blob = subprocess.run(
            ["git", "-C", str(ROOT), "archive", "--format=tar", LEGACY_COMMIT,
             "teddy_mm", "scripts/06_train_bidirectional.py", "scripts/phase2_decode_check.py",
             "bridge_anm/export_missing_modality_events.py"],
            capture_output=True, check=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        if env_tree and (Path(env_tree) / "teddy_mm/bidirectional.py").exists():
            return Path(env_tree)
        pytest.skip(f"legacy commit {LEGACY_COMMIT} not reachable (set PHASE2_LEGACY_TREE)")
    out = tmp_path_factory.mktemp("legacy")
    with tarfile.open(fileobj=io.BytesIO(blob)) as tf:
        tf.extractall(out, filter="data")
    return out


def _batch(seed: int = 0, n: int = 24, counts: bool = False):
    g = torch.Generator().manual_seed(seed)
    z = torch.nn.functional.normalize(torch.randn(n, Z_DIM, generator=g), dim=-1)
    if counts:
        adt = torch.poisson(torch.full((n, N_ADT), 3.0), generator=g)
    else:  # already-normalized: non-negative floats with exact zeros, like data/processed/cite
        adt = torch.log1p(torch.poisson(torch.full((n, N_ADT), 3.0), generator=g) * 0.37)
    sf = torch.rand(n, generator=g) + 0.5
    return z, adt, sf


# ------------------------------------------------------------ transform once
def _encoder_input(enc: AdtEncoder, adt: torch.Tensor) -> torch.Tensor:
    seen = []
    h = enc.mlp.register_forward_pre_hook(lambda m, a: seen.append(a[0].detach().clone()))
    enc.eval()
    enc(adt)
    h.remove()
    assert len(seen) == 1
    return seen[0]


def test_default_encoder_does_not_retransform_stored_adt():
    _, adt, _ = _batch()
    enc = AdtEncoder(N_ADT, Z_DIM)
    assert enc.input_transform == "none"
    assert torch.equal(_encoder_input(enc, adt), adt)


def test_each_transform_is_applied_exactly_once():
    _, adt, _ = _batch(counts=True)
    for name, ref in (("clr", clr_log_adt(adt)), ("log1p", torch.log1p(adt)), ("none", adt)):
        x = _encoder_input(AdtEncoder(N_ADT, Z_DIM, input_transform=name), adt)
        assert torch.allclose(x, ref), name
    # A second clr on top of the first changes the values (the legacy double transform).
    assert not torch.allclose(clr_log_adt(clr_log_adt(adt)), clr_log_adt(adt))


def test_transform_called_once_per_encoder_pass(monkeypatch):
    calls = []
    real = bd.transform_adt_input
    monkeypatch.setattr(bd, "transform_adt_input", lambda a, t: calls.append(t) or real(a, t))
    z, adt, sf = _batch()
    torch.manual_seed(0)
    model = BidirectionalCite(Z_DIM, N_ADT, hidden=32, adt_hidden=16)
    model.forward_train(z, adt, sf, fm_steps=4)
    assert calls == ["none"]
    calls.clear()
    for mode in ("rna_only", "adt_only", "joint"):
        model.predict_adt(z, adt, sf, mode=mode, fm_steps=4)
    assert calls == ["none"] * 3


def test_auto_resolves_against_stored_adt():
    _, normed, _ = _batch()
    _, counts, _ = _batch(counts=True)
    assert resolve_adt_input_transform("auto", normed.numpy()) == "none"
    assert resolve_adt_input_transform("auto", counts.numpy()) == "log1p"
    assert resolve_adt_input_transform("clr", normed.numpy()) == "clr"
    with pytest.raises(ValueError):
        resolve_adt_input_transform("bogus", normed.numpy())


def test_train_script_feeds_stored_adt_unchanged(tmp_path):
    tr = _load_module("train06", ROOT / "scripts/06_train_bidirectional.py")
    pack = _write_pack(tmp_path, counts=False)
    from teddy_mm.data import load_prepared
    loaded = load_prepared(pack)
    z = np.load(pack / "z_rna.npy")
    idx = np.array([0, 3, 5])
    sub = tr.subset(loaded, z, idx, np.ones(3))
    assert np.array_equal(sub["adt"].numpy(), loaded["adt"][idx])


# --------------------------------------------- no test-time answer-key input
def test_rna_only_prediction_never_reads_adt():
    z, adt, sf = _batch()
    torch.manual_seed(0)
    model = BidirectionalCite(Z_DIM, N_ADT, hidden=32, adt_hidden=16).eval()
    outs = []
    for a in (adt, torch.full_like(adt, float("nan")), None):
        for use_fm in (True, False):
            torch.manual_seed(1)
            outs.append((use_fm, model.predict_adt(z, a, sf, mode="rna_only", use_fm=use_fm, fm_steps=4)))
    for use_fm, o in outs:
        ref = outs[0][1] if use_fm else outs[1][1]
        assert torch.isfinite(o).all() and torch.equal(o, ref)


def test_eval_size_factor_train_median_ignores_eval_cells():
    sf = np.array([1.0, 2.0, 3.0, 4.0, np.nan, np.nan], dtype=np.float32)
    out, info = eval_size_factor("train-median", sf, np.arange(4), np.array([4, 5]))
    assert np.allclose(out, 2.5) and info["value"] == 2.5 and info["n_train_cells"] == 4
    meas, _ = eval_size_factor("measured", sf, np.arange(4), np.array([4, 5]))
    assert np.isnan(meas).all()
    with pytest.raises(ValueError):
        eval_size_factor("train-median", sf, np.arange(5), np.array([4, 5]))


def _write_pack(tmp: Path, *, counts: bool, poison_eval_sf: bool = False, n: int = 90,
                panel: list[str] | None = None) -> Path:
    rng = np.random.default_rng(0)
    names = (panel or []) + [f"P{i}" for i in range(N_ADT)]
    programs = rng.normal(size=(n, 3))
    lam = np.exp(programs @ rng.normal(scale=0.6, size=(3, len(names)))) * 2.0
    adt = rng.poisson(lam).astype(np.float32)
    if not counts:
        adt = np.log1p(adt / np.exp(np.log1p(adt).mean(axis=1, keepdims=True))).astype(np.float32)
    split = np.array(["train"] * 50 + ["val"] * 15 + ["test"] * (n - 65), dtype="U8")
    sites = np.array(["site1"] * 65 + ["site4"] * (n - 65), dtype="U32")
    tot = adt.sum(axis=1)
    sf = (tot / np.median(tot)).astype(np.float32)
    if poison_eval_sf:
        sf[50:] = np.nan  # any read of a val/test cell's measured depth turns into NaN
    rna = sparse.csr_matrix(rng.poisson(1.0, size=(n, 5)).astype(np.float32))
    z = np.concatenate([programs, rng.normal(scale=0.1, size=(n, Z_DIM - 3))], axis=1).astype(np.float32)
    out = tmp / ("pack_counts" if counts else "pack_norm") / ("poison" if poison_eval_sf else "clean")
    out.mkdir(parents=True)
    np.savez(
        out / "cite_arrays.npz",
        rna_data=rna.data, rna_indices=rna.indices, rna_indptr=rna.indptr, rna_shape=np.array(rna.shape),
        adt=adt, token_ids=np.arange(5), split=split,
        donors=np.array([f"d{i % 3}" for i in range(n)], dtype="U32"), sites=sites,
        cell_types=np.array(["T"] * n, dtype="U64"), rna_size_factor=np.ones(n, dtype=np.float32),
        adt_size_factor=sf, adt_names=np.array(names, dtype="U64"),
        rna_names=np.array([f"ENSG{i:011d}" for i in range(5)], dtype="U32"),
    )
    np.save(out / "z_rna.npy", z)
    (out / "meta.json").write_text(json.dumps({"n_cells": n}))
    return out


def _run06(script: Path, pack: Path, out: Path, *extra: str, new: bool = True) -> dict:
    cmd = [sys.executable, str(script), "--processed", str(pack), "--out", str(out),
           "--epochs", "2", "--batch-size", "16", "--hidden", "32", "--adt-hidden", "16",
           "--fm-steps", "4", "--device", "cpu", "--seed", "0", *extra]
    if new:
        cmd += ["--threads", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads((out / "metrics.json").read_text())


def test_train_script_eval_does_not_read_measured_eval_adt(tmp_path):
    pack = _write_pack(tmp_path, counts=False, poison_eval_sf=True)
    m = _run06(ROOT / "scripts/06_train_bidirectional.py", pack, tmp_path / "fixed")
    assert m["run"]["test_size_factor"]["mode"] == "train-median"
    assert m["run"]["adt_input_transform"] == "none" and m["run"]["stored_adt_count_like"] is False
    assert all(np.isfinite(v) for v in m["test"].values())
    assert all(np.isfinite(h["val_rna_only_pearson"]) for h in m["history"])
    ck = torch.load(tmp_path / "fixed/best.pt", map_location="cpu", weights_only=False)
    assert ck["config"]["adt_input_transform"] == "none"
    sf = np.load(pack / "cite_arrays.npz")["adt_size_factor"]
    assert ck["train_median_size_factor"] == pytest.approx(float(np.median(sf[:50].astype(np.float64))))
    # The legacy flag reads the evaluated cells' measured depth, so the poison shows up.
    legacy = _run06(ROOT / "scripts/06_train_bidirectional.py", pack, tmp_path / "legacy", "--legacy")
    assert legacy["run"]["adt_input_transform"] == "clr"
    assert all(np.isnan(v) for v in legacy["test"].values())


def test_train_script_auto_uses_log1p_on_counts(tmp_path):
    pack = _write_pack(tmp_path, counts=True)
    m = _run06(ROOT / "scripts/06_train_bidirectional.py", pack, tmp_path / "o", "--max-train-cells", "20")
    assert m["run"]["adt_input_transform"] == "log1p" and m["run"]["stored_adt_count_like"] is True
    assert m["run"]["n_train_cells"] == 20


def test_legacy_flag_conflict_is_rejected(tmp_path):
    pack = _write_pack(tmp_path, counts=False)
    r = subprocess.run([sys.executable, str(ROOT / "scripts/06_train_bidirectional.py"), "--processed", str(pack),
                        "--out", str(tmp_path / "o"), "--legacy", "--adt-input-transform", "none"],
                       capture_output=True, text=True, env=ENV)
    assert r.returncode != 0 and "conflicts" in r.stderr


# --------------------------------------------------- legacy reproduces old code
def test_legacy_model_matches_pre_fix_model(legacy_tree):
    old = _load_module("old_bidirectional", legacy_tree / "teddy_mm/bidirectional.py")
    z, adt, sf = _batch(seed=3)
    torch.manual_seed(0)
    m_old = old.BidirectionalCite(Z_DIM, N_ADT, hidden=32, adt_hidden=16)
    torch.manual_seed(0)
    m_new = BidirectionalCite(Z_DIM, N_ADT, hidden=32, adt_hidden=16, adt_input_transform="clr")
    so, sn = m_old.state_dict(), m_new.state_dict()
    assert so.keys() == sn.keys() and all(torch.equal(so[k], sn[k]) for k in so)
    m_old.train(); m_new.train()
    torch.manual_seed(5)
    lo = m_old.forward_train(z, adt, sf, fm_steps=4)
    torch.manual_seed(5)
    ln = m_new.forward_train(z, adt, sf, fm_steps=4)
    for k in lo:
        assert torch.equal(lo[k], ln[k]), k
    m_old.eval(); m_new.eval()
    for mode in ("rna_only", "adt_only", "joint"):
        for use_fm in (True, False):
            torch.manual_seed(7)
            po = m_old.predict_adt(z, adt, sf, mode=mode, fm_steps=4, use_fm=use_fm)
            torch.manual_seed(7)
            pn = m_new.predict_adt(z, adt, sf, mode=mode, fm_steps=4, use_fm=use_fm)
            assert torch.equal(po, pn), (mode, use_fm)
    # An old checkpoint (no config) is rebuilt with the legacy clr encoder.
    torch.manual_seed(0)
    m_def = old.BidirectionalCite(Z_DIM, N_ADT)
    m_load, info = load_bidirectional({"model": m_def.state_dict(), "z_dim": Z_DIM, "n_adt": N_ADT})
    assert info["adt_input_transform"] == "clr" and info["adt_input_transform_source"] == "legacy_default_clr"
    m_def.eval()
    torch.manual_seed(7)
    a = m_def.predict_adt(z, adt, sf, mode="adt_only", fm_steps=4)
    torch.manual_seed(7)
    b = m_load.predict_adt(z, adt, sf, mode="adt_only", fm_steps=4)
    assert torch.equal(a, b)


def test_legacy_train_script_matches_pre_fix_script(tmp_path, legacy_tree):
    pack = _write_pack(tmp_path, counts=False)
    old = _run06(legacy_tree / "scripts/06_train_bidirectional.py", pack, tmp_path / "old", new=False)
    new = _run06(ROOT / "scripts/06_train_bidirectional.py", pack, tmp_path / "new", "--legacy")
    assert len(old["history"]) == len(new["history"])
    for ho, hn in zip(old["history"], new["history"]):
        for k, v in ho.items():
            assert hn[k] == pytest.approx(v, rel=1e-6, abs=1e-7), k
    for k, v in old["test"].items():
        assert new["test"][k] == pytest.approx(v, rel=1e-6, abs=1e-7), k
    fixed = _run06(ROOT / "scripts/06_train_bidirectional.py", pack, tmp_path / "fixed")
    assert fixed["test"]["test_adt_only_pearson"] != pytest.approx(old["test"]["test_adt_only_pearson"], abs=1e-9)


# ----------------------------------------------- missing-modality export path
def _export_fixture(tmp: Path, *, legacy_ckpt: bool, poison: bool = True):
    from lib.lineage_panels import all_panel_proteins
    from teddy_mm.models import AdtDecoder, MLP

    panel = all_panel_proteins()
    pack = _write_pack(tmp, counts=False, poison_eval_sf=poison, panel=panel)
    n_adt = len(panel) + N_ADT
    torch.manual_seed(0)
    p1 = tmp / "p1.pt"
    torch.save({"mlp": MLP(Z_DIM, Z_DIM, hidden=512).state_dict(),
                "dec": AdtDecoder(Z_DIM, n_adt, hidden=512).state_dict()}, p1)
    model = BidirectionalCite(Z_DIM, n_adt, adt_input_transform="clr" if legacy_ckpt else "none")
    blob = {"model": model.state_dict(), "z_dim": Z_DIM, "n_adt": n_adt,
            "metrics": {"val_rna_only_pearson": 0.9}}
    if not legacy_ckpt:
        blob.update(config=model.config, train_median_size_factor=1.25)
    p2 = tmp / "p2.pt"
    torch.save(blob, p2)
    return pack, p1, p2


# Runs a script with torch seeded (the export draws flow-matching noise without seeding
# torch) and MPS hidden (the pre-fix decode check picks MPS whenever it is available).
_SEEDED = (
    "import runpy, sys, torch\n"
    "torch.backends.mps.is_available = lambda: False\n"
    "torch.manual_seed(0)\n"
    "sys.argv = sys.argv[1:]\n"
    "runpy.run_path(sys.argv[0], run_name='__main__')\n"
)


def _run_script(script: Path, *args: str, seeded: bool = False, cwd: Path | None = None,
                pythonpath: Path | None = None) -> subprocess.CompletedProcess:
    env = dict(ENV)
    if pythonpath is not None:
        env["PYTHONPATH"] = str(pythonpath)
    head = [sys.executable, "-c", _SEEDED] if seeded else [sys.executable]
    r = subprocess.run([*head, str(script), *args], capture_output=True, text=True, env=env, cwd=cwd)
    assert r.returncode == 0, r.stderr[-2000:]
    return r


def _run_export(tmp: Path, pack, p1, p2, *extra, script: Path | None = None, tag: str = "",
                seeded: bool = False, pythonpath: Path | None = None):
    name = "out_" + (tag or ("_".join(x.strip("-") for x in extra) if extra else "default"))
    out = tmp / name
    _run_script(script or (BRIDGE / "export_missing_modality_events.py"), "--processed", str(pack),
                "--phase1-ckpt", str(p1), "--phase2-ckpt", str(p2), "--out-dir", str(out),
                "--device", "cpu", "--fm-steps", "4", *extra,
                seeded=seeded, cwd=tmp, pythonpath=pythonpath)
    man = json.loads((out / "export_manifest.json").read_text())
    preds = {m: np.load(out / f"adt_pred_panel_{m}.npy") for m in ("rna_only", "adt_only", "joint")}
    return man, preds, out


def test_missing_modality_export_uses_no_measured_test_depth(tmp_path):
    pack, p1, p2 = _export_fixture(tmp_path, legacy_ckpt=False)
    man, preds, _ = _run_export(tmp_path, pack, p1, p2)
    assert man["source"] == "phase2_bidirectional"
    assert all(np.isfinite(v).all() for v in preds.values())
    assert man["size_factor_adt_observed_masks"]["mode"] == "train-median"
    assert man["phase2_model"]["adt_input_transform"] == "none"
    assert man["phase2_model"]["train_median_source"] == "checkpoint"
    assert man["phase2_model"]["gate"]["split"] == "val"
    assert man["answer_key_is_encoder_input"] == {"rna_only": False, "adt_only": True, "joint": True}
    # Legacy: ADT-observed masks read the exported cells' measured depth -> poison propagates.
    man_l, preds_l, _ = _run_export(tmp_path, pack, p1, p2, "--size-factor-observed", "measured")
    assert np.isfinite(preds_l["rna_only"]).all()
    assert np.isnan(preds_l["adt_only"]).all() and np.isnan(preds_l["joint"]).all()


def test_missing_modality_export_rebuilds_legacy_ckpt_with_clr(tmp_path):
    pack, p1, p2 = _export_fixture(tmp_path, legacy_ckpt=True)
    man, preds, _ = _run_export(tmp_path, pack, p1, p2)
    assert man["phase2_model"]["adt_input_transform"] == "clr"
    assert man["phase2_model"]["adt_input_transform_source"] == "legacy_default_clr"
    assert man["phase2_model"]["train_median_source"] == "split=='train' median"
    assert all(np.isfinite(v).all() for v in preds.values())


def test_missing_modality_export_decodes_with_checkpoint_train_median(tmp_path):
    # Phase-2 masks decode with the checkpoint's train median (1.25), not the pack's
    # split=='train' median: every prediction is exactly 1.25x the size-factor-one run.
    pack, p1, p2 = _export_fixture(tmp_path, legacy_ckpt=False)
    man, preds, _ = _run_export(tmp_path, pack, p1, p2, seeded=True, tag="median")
    man_1, preds_1, _ = _run_export(tmp_path, pack, p1, p2, "--size-factor", "one",
                                    "--size-factor-observed", "one", seeded=True, tag="one")
    assert man["source"] == man_1["source"] == "phase2_bidirectional"
    for m in ("rna_only", "adt_only", "joint"):
        assert np.allclose(preds[m], 1.25 * preds_1[m], rtol=1e-5, atol=0), m


@pytest.mark.parametrize("min_pearson", [None, "-1"])
def test_legacy_export_flags_match_pre_fix_export(tmp_path, legacy_tree, min_pearson):
    # Pre-fix export vs --size-factor-observed measured --phase2-gate test, both seeded,
    # on a pre-fix-style checkpoint. min_pearson=-1 forces the phase-2 source; the default
    # threshold sends the untrained model to the hybrid source.
    old_script = legacy_tree / "bridge_anm/export_missing_modality_events.py"
    if not old_script.exists():
        pytest.skip("legacy tree has no bridge_anm/export_missing_modality_events.py")
    pack, p1, p2 = _export_fixture(tmp_path, legacy_ckpt=True, poison=False)
    extra = () if min_pearson is None else ("--phase2-min-pearson", min_pearson)
    # The pre-fix script imports lib/ and export_cite_events (unchanged by the fix) from its
    # own folder; the legacy tree only holds the export itself, so they come from BRIDGE.
    man_o, preds_o, out_o = _run_export(tmp_path, pack, p1, p2, *extra, script=old_script, tag="old",
                                        seeded=True, pythonpath=BRIDGE)
    man_n, preds_n, out_n = _run_export(tmp_path, pack, p1, p2, *extra, "--size-factor-observed", "measured",
                                        "--phase2-gate", "test", tag="new", seeded=True)
    assert man_n["source"] == man_o["source"]
    assert man_o["source"] == ("hybrid_phase1_rna_phase2_adt" if min_pearson is None else "phase2_bidirectional")
    for m in preds_o:
        assert np.array_equal(preds_o[m], preds_n[m]), m
    for k in ("phase2_metrics", "panel_pearson_by_mask", "phase1_full_adt_pearson", "p95_train"):
        assert man_n[k] == man_o[k], k
    for f in ("missing_modality_events.jsonl", "missing_modality_cells.jsonl", "adt_true_panel.npy",
              "global_indices.npy"):
        assert (out_o / f).read_bytes() == (out_n / f).read_bytes(), f


# ------------------------------------------------------- phase-2 decode check
def _decode_ckpt(tmp: Path, *, legacy_ckpt: bool) -> Path:
    torch.manual_seed(0)
    model = BidirectionalCite(Z_DIM, N_ADT, adt_input_transform="clr" if legacy_ckpt else "none")
    blob = {"model": model.state_dict(), "z_dim": Z_DIM, "n_adt": N_ADT}
    if not legacy_ckpt:
        blob.update(config=model.config, train_median_size_factor=1.25)
    d = tmp / ("ck_legacy" if legacy_ckpt else "ck_new")
    d.mkdir()
    torch.save(blob, d / "best.pt")
    (d / "metrics.json").write_text(json.dumps({"test": {"test_rna_only_pearson": 0.0}}))
    return d / "best.pt"


def _run_decode(pack: Path, ckpt: Path, out: Path, *extra: str) -> dict:
    _run_script(ROOT / "scripts/phase2_decode_check.py", "--processed", str(pack), "--ckpt", str(ckpt),
                "--out", str(out), "--device", "cpu", *extra)
    return json.loads(out.read_text())


def _finite_by_mode(res: dict) -> bool:
    return all(np.isfinite(v) for row in res["by_mode"].values() for v in row.values())


def test_decode_check_uses_no_measured_test_depth(tmp_path):
    pack = _write_pack(tmp_path, counts=False, poison_eval_sf=True)
    new_ck = _decode_ckpt(tmp_path, legacy_ckpt=False)
    res = _run_decode(pack, new_ck, tmp_path / "new.json")
    assert res["size_factor"] == {"mode": "train-median", "value": 1.25, "source": "checkpoint"}
    assert res["model"]["adt_input_transform"] == "none" and _finite_by_mode(res)
    old_ck = _decode_ckpt(tmp_path, legacy_ckpt=True)
    res_l = _run_decode(pack, old_ck, tmp_path / "old.json")
    sf = np.load(pack / "cite_arrays.npz")["adt_size_factor"]
    assert res_l["size_factor"]["value"] == pytest.approx(float(np.median(sf[:50].astype(np.float64))))
    assert res_l["model"]["adt_input_transform"] == "clr" and _finite_by_mode(res_l)
    # Legacy --size-factor measured reads each test cell's measured depth -> poison propagates.
    res_m = _run_decode(pack, new_ck, tmp_path / "measured.json", "--size-factor", "measured")
    assert all(np.isnan(v) for row in res_m["by_mode"].values() for v in row.values())


def test_legacy_decode_check_matches_pre_fix_script(tmp_path, legacy_tree):
    # The pre-fix script hard-codes ROOT/data/processed/cite and ROOT/outputs/cite_phase2,
    # so it runs from a scratch ROOT holding the legacy teddy_mm, the pack and the checkpoint.
    import shutil

    pack = _write_pack(tmp_path, counts=False)
    ck = _decode_ckpt(tmp_path, legacy_ckpt=True)
    lroot = tmp_path / "legacy_root"
    (lroot / "scripts").mkdir(parents=True)
    shutil.copy(legacy_tree / "scripts/phase2_decode_check.py", lroot / "scripts")
    shutil.copytree(legacy_tree / "teddy_mm", lroot / "teddy_mm", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(pack, lroot / "data/processed/cite")
    shutil.copytree(ck.parent, lroot / "outputs/cite_phase2")
    _run_script(lroot / "scripts/phase2_decode_check.py", seeded=True, cwd=tmp_path)
    old = json.loads((lroot / "outputs/cite_phase2/decode_check.json").read_text())
    new = _run_decode(pack, ck, tmp_path / "new.json", "--size-factor", "measured")
    assert new["test_cells"] == old["test_cells"] and new["recorded"] == old["recorded"]
    assert new["by_mode"] == old["by_mode"]
