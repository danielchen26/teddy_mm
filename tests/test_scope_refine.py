"""Tests for bridge_anm/run_scope_refine_proof.py.

Run with ANM_ROOT pointing at the ANM checkout. Data-backed tests skip when the
gitignored bridge export is not present.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BRIDGE_PY = REPO / "bridge_anm"
SCRIPT = BRIDGE_PY / "run_scope_refine_proof.py"
MAIN_REPO = Path("/Users/tianchichen/Documents/GitHub/teddy_mm")
DATA = Path(os.environ.get("SCOPE_REFINE_BRIDGE_DIR", str(MAIN_REPO / "outputs/anm_cite_bridge")))
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))

sys.path.insert(0, str(BRIDGE_PY))
os.environ.setdefault("ANM_ROOT", str(ANM_ROOT))

srp = pytest.importorskip("run_scope_refine_proof")

HAVE_DATA = (DATA / "cite_cells_meta.jsonl").exists() and (DATA / "hard_proof/attr_compact.npz").exists()
needs_data = pytest.mark.skipif(not HAVE_DATA, reason=f"bridge export not found under {DATA}")

N_SUBSET = 240


# ---------------------------------------------------------------- pure helpers
def test_q_at_coverage_ranks_by_score_and_breaks_ties_by_id():
    rows = [
        {"cell_id": "c", "s": 0.9, "pred": "a", "label": "a"},
        {"cell_id": "a", "s": 0.5, "pred": "a", "label": "b"},
        {"cell_id": "b", "s": 0.5, "pred": "a", "label": "a"},
        {"cell_id": "d", "s": 0.1, "pred": "a", "label": None},
    ]
    out = srp.q_at_coverage(rows, "pred", "s", [0.25, 0.5, 1.0])
    assert [o["n_answered"] for o in out] == [1, 2, 4]
    assert out[0]["Q"] == 1.0
    assert out[1]["Q"] == 0.5  # tie at 0.5 -> "a" (wrong) comes before "b"
    assert out[2]["n_labeled_answered"] == 3 and abs(out[2]["Q"] - 2 / 3) < 1e-12


def test_teddy_lineage_margin():
    assert srp.teddy_lineage_margin({"x": 0.7, "y": 0.2, "z": 0.5}) == pytest.approx(0.2)


def test_coarse_cell_type():
    cases = {
        "CD4+ T activated": "T",
        "gdT CD158b+": "T",
        "MAIT": "T",
        "Naive CD20+ B IGKC+": "B / plasma",
        "Transitional B": "B / plasma",
        "Plasmablast IGKC-": "B / plasma",
        "NK CD158e1+": "NK / ILC",
        "ILC1": "NK / ILC",
        "CD14+ Mono": "myeloid (mono/DC)",
        "pDC": "myeloid (mono/DC)",
        "HSC": "progenitor",
        "G/M prog": "progenitor",
        "Proerythroblast": "erythroid",
        "Reticulocyte": "erythroid",
        None: "unknown",
    }
    for ct, want in cases.items():
        assert srp.coarse_cell_type(ct) == want, ct


def test_cell_type_composition_counts_and_enrichment():
    cells = {f"c{i}": {"cell_type": "NK" if i < 2 else "CD14+ Mono"} for i in range(8)}
    comp = srp.cell_type_composition(["c0", "c1", "c2"], cells)
    fine = {r["cell_type"]: r for r in comp["fine"]}
    assert comp["n_flagged"] == 3 and comp["n_all"] == 8
    assert fine["NK"]["n_flagged"] == 2 and fine["NK"]["flag_rate"] == 1.0
    assert fine["NK"]["enrichment"] == pytest.approx((2 / 3) / (2 / 8))
    assert comp["fine"][0]["cell_type"] == "NK"  # sorted by n_flagged


# ---------------------------------------------------------------- data-backed
@pytest.fixture(scope="module")
def subset_bridge(tmp_path_factory):
    """Small bridge dir: first N site4 cells + their events + attr; no missing_modality/."""
    if not HAVE_DATA:
        pytest.skip("no data")
    out = tmp_path_factory.mktemp("bridge_subset")
    rows = []
    with (DATA / "cite_cells_meta.jsonl").open() as f:
        for line in f:
            if '"site4_test"' in line:
                rows.append(json.loads(line))
    rows.sort(key=lambda r: r["cell_id"])
    rows = rows[:N_SUBSET]
    want = {r["cell_id"] for r in rows}
    with (out / "cite_cells_meta.jsonl").open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    with (DATA / "cite_typed_events.jsonl").open() as fi, (out / "cite_typed_events.jsonl").open("w") as fo:
        for line in fi:
            if json.loads(line)["cell_id"] in want:
                fo.write(line)
    (out / "hard_proof").mkdir()
    (out / "hard_proof/attr_compact.npz").write_bytes((DATA / "hard_proof/attr_compact.npz").read_bytes())
    return out


def _strip_timing(obj):
    if isinstance(obj, dict):
        return {
            k: _strip_timing(v)
            for k, v in obj.items()
            if k not in {"elapsed_sec", "runtime_sec", "timestamp_local", "workers"}
        }
    if isinstance(obj, list):
        return [_strip_timing(v) for v in obj]
    return obj


@needs_data
def test_gate_claim_workers_1_vs_4_identical(subset_bridge):
    cells = srp.load_site4_cells(subset_bridge / "cite_cells_meta.jsonl")
    by, p95 = srp.load_events_for(subset_bridge / "cite_typed_events.jsonl", set(cells))
    schema = json.loads((REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    r1 = srp.run_gate_claim(cells, by, p95, schema, ["O0", "O2"], n_tau=8, workers=1)
    r4 = srp.run_gate_claim(cells, by, p95, schema, ["O0", "O2"], n_tau=8, workers=4)
    assert _strip_timing(r1) == _strip_timing(r4)
    blk = r1["criteria"]["O0"]
    assert blk["coverage_Q_curve_teddy_gated_by_teddy_margin"]
    assert blk["matched_coverage_soft_P_vs_teddy_margin"][0]["coverage"] == 1.0
    # at full coverage both gates answer every cell, so Q must agree
    m0 = blk["matched_coverage_soft_P_vs_teddy_margin"][0]
    assert m0["Q_gate_anm_soft_P"] == m0["Q_gate_teddy_margin"] == blk["baseline"]["teddy_always_Q"]


def _run_main(bridge: Path, out: Path, *extra: str) -> dict:
    env = {**os.environ, "ANM_ROOT": str(ANM_ROOT)}
    cmd = [
        sys.executable, str(SCRIPT), "--bridge-dir", str(bridge), "--out-dir", str(out),
        "--n-tau", "8", "--n-random", "4", "--n-seeds", "2", *extra,
    ]
    res = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=600)
    assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-4000:]
    assert (out / "SCOPE_REFINE_PROOF.md").exists()
    return json.loads((out / "scope_refine_results.json").read_text())


@needs_data
def test_main_end_to_end_without_missing_modality(subset_bridge, tmp_path):
    a = _run_main(subset_bridge, tmp_path / "w1", "--skip-adt-only", "--workers", "1")
    b = _run_main(subset_bridge, tmp_path / "w4", "--skip-adt-only", "--workers", "4")
    assert a["run_args"]["adt_only_status"] == "skipped_by_flag"
    assert "claim1_gate_adt_only" not in a
    assert a["paths"]["bridge_dir"] == str(subset_bridge)
    comp = a["claim2_attribution"]["hard_cell_type_composition"]
    assert comp["n_flagged"] == a["claim2_attribution"]["n_flip_sensitive"]
    assert sum(r["n_flagged"] for r in comp["coarse"]) == comp["n_flagged"]
    for doc in (a, b):
        doc.pop("paths")
        doc["run_args"].pop("workers")
    assert _strip_timing(a) == _strip_timing(b)
    # without the flag, a bridge dir lacking missing_modality/ still runs
    c = _run_main(subset_bridge, tmp_path / "nf", "--max-cells", "60")
    assert c["run_args"]["adt_only_status"] == "skipped_missing_files"
