"""Tests for E6 (bridge_anm/lib/v3_e6.py, scripts/v3_e6_external.py, registration/addenda/E6.json).

Synthetic tests always run. The engine test and the val run through the E6 code path skip when ANM_ROOT or the BMMC
data are unavailable. No test opens the external pack values or the external embedding.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e1 as e1  # noqa: E402
from lib import v3_e6 as e6  # noqa: E402
from lib import v3_key as vk  # noqa: E402

ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))
ADDENDUM = REPO / "registration" / "addenda" / "E6.json"
HAVE_BMMC = (REPO / "data/processed/cite/cite_arrays.npz").exists() and (REPO / "data/processed/cite_official/z_rna.npy").exists()


@pytest.fixture(scope="module")
def reg():
    return va.load_registration_amended()


@pytest.fixture(scope="module")
def add():
    return json.loads(ADDENDUM.read_text())


def _runner():
    spec = importlib.util.spec_from_file_location("v3_e6_external", REPO / "scripts" / "v3_e6_external.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _random_H(rng, proteins, n=4000):
    return {p: rng.random(n) < 0.35 for p in proteins}


# ----------------------------------------------------------------------------- gate with absent proteins
def test_reduced_gate_with_every_protein_equals_registered_gate(reg):
    rng = np.random.default_rng(0)
    gate = reg["gate"]
    spec = e6.reduced_gate_spec(gate, gate["proteins"])
    assert not spec["dropped_clauses"] and not spec["removed_conditions"]
    for _ in range(5):
        H = _random_H(rng, gate["proteins"])
        np.testing.assert_array_equal(e6.gate_class_from_spec(H, spec), vk.gate_class_from_high(H, gate))


def test_reduced_gate_external_panel_clauses(reg, add):
    gate = reg["gate"]
    present = [p for p in gate["proteins"] if p not in ("CD94", "CD33")]
    spec = e6.reduced_gate_spec(gate, present)
    assert spec["absent_gate_proteins"] == ["CD33", "CD94"]
    assert spec["dropped_clauses"] == []
    cl = {c["name"]: c for c in spec["clauses"]}
    assert cl["NK"]["count_high"] == ["CD56", "CD335", "CD16"] and cl["NK"]["count_min"] == 2
    assert "CD33" not in cl["NK"]["none_high"] and cl["myeloid"]["any_high"] == ["CD14", "CD11c"]
    assert spec["clauses"] == add["key"]["clauses_remaining_in_order"]
    # behaviour on hand-made cells
    names = ["CD3", "CD19", "CD20", "CD56", "CD335", "CD16", "CD14", "CD11c", "CD71"]
    rows = {"cd56_cd16": {"CD56", "CD16"}, "cd56_only": {"CD56"}, "cd11c_only": {"CD11c"}, "none": set(),
            "cd71": {"CD71", "CD3"}, "t": {"CD3"}, "b": {"CD19"}, "nk_cd14": {"CD56", "CD16", "CD14"}}
    H = {p: np.asarray([p in s for s in rows.values()]) for p in names}
    got = dict(zip(rows, e6.gate_class_from_spec(H, spec)))
    assert got == {"cd56_cd16": "NK", "cd56_only": "OUT", "cd11c_only": "myeloid", "none": "OUT", "cd71": "OUT",
                   "t": "T", "b": "B", "nk_cd14": "myeloid"}


def test_registered_examples_of_the_missing_protein_rule(reg):
    gate = reg["gate"]
    spec = e6.reduced_gate_spec(gate, [p for p in gate["proteins"] if p not in ("CD335", "CD71")])
    assert [d["name"] for d in spec["dropped_clauses"]] == ["erythroid"]            # CD71 absent -> clause dropped
    nk = [c for c in spec["clauses"] if c["name"] == "NK"][0]
    assert nk["count_high"] == ["CD56", "CD94", "CD16"] and nk["count_min"] == 2     # CD335 absent -> count over the rest
    rng = np.random.default_rng(1)
    H = _random_H(rng, gate["proteins"])
    H0 = dict(H)
    H0["CD71"] = np.zeros_like(H["CD71"])
    H0["CD335"] = np.zeros_like(H["CD335"])
    np.testing.assert_array_equal(e6.gate_class_from_spec(H, spec), vk.gate_class_from_high(H0, gate))
    with pytest.raises(KeyError):
        e6.gate_class_from_spec({p: v for p, v in H.items() if p != "CD56"}, spec)


def test_per_donor_gate_matches_registered_estimator(reg):
    rng = np.random.default_rng(2)
    gate = reg["gate"]
    names = list(gate["proteins"])
    n = 1200
    adt = np.abs(rng.normal(0.4, 0.5, (n, len(names)))) * (rng.random((n, len(names))) < 0.8)
    adt[rng.random((n, len(names))) < 0.3] += 2.0
    donors = np.where(np.arange(n) < 500, "dA", "dB")
    spec = e6.reduced_gate_spec(gate, names)
    thr = e6.per_donor_thresholds(adt, names, donors, spec["proteins_used"], seed=0)
    for d in ("dA", "dB"):
        m = donors == d
        for p in spec["proteins_used"]:
            assert thr[d][p] == vk.gmm_nonzero_threshold(adt[m, names.index(p)], seed=0)
    got = e6.gate_per_donor(adt, names, donors, spec, thr)
    for d in ("dA", "dB"):
        m = donors == d
        full_thr = {p: thr[d].get(p, 99.0) for p in gate["proteins"]}
        np.testing.assert_array_equal(got[m], vk.gate_class(adt[m], names, reg, thresholds=full_thr))
    adt_nan = adt.copy()
    adt_nan[:, names.index("CD94")] = np.nan
    H = e6.high_arrays(adt_nan, names, {"CD94": 0.0})
    assert not H["CD94"].any()


# ----------------------------------------------------------------------------- keys and validation
def test_keys_doublet_and_validation_rule():
    annot = np.asarray(["B", "T", "unscored", "NK", "OUT", "myeloid"] * 50)
    gated = np.asarray(["B", "T", "T", "NK", "OUT", "OUT"] * 50)
    donors = np.repeat(["d1", "d2", "d3"], 100)
    k = e6.e6_keys(annot, gated)
    assert np.all(k["primary"][annot == "unscored"] == "unscored")
    assert np.all(k["annotation_only"][annot == "unscored"] == "unscored")
    assert np.all(k["primary"][annot == "myeloid"] == "unscored")
    kv = e6.key_validation(annot, gated, donors)
    assert kv["n_annotation_unscored"] == 50 and kv["validated"] and kv["decision_key"] == "primary"
    bad = np.where(np.arange(annot.size) % 2 == 0, "OUT", gated)
    kv2 = e6.key_validation(annot, bad, donors)
    assert kv2["kappa5"] < 0.70 and not kv2["validated"]
    assert kv2["decision_key"] == "annotation_only" and kv2["label"] == "key not validated"


def test_annotation_map_hash_rederivation_and_addendum(add):
    sys.path.insert(0, str(REPO / "bridge_anm"))
    import v3_e6_addendum as ad

    assert e6.map_sha256(ad.INDEPENDENT_MAP) == ad.COMMITTED_MAP_SHA
    assert add["annotation_map"]["map"] == ad.INDEPENDENT_MAP
    assert add["annotation_map"]["sha256_canonical_json"] == ad.COMMITTED_MAP_SHA
    assert add["annotation_map"]["independent_rederivation"]["disagreements"] == []
    assert add["annotation_map"]["independent_rederivation"]["n_agree"] == 31
    assert {k: v["decision"] for k, v in add["annotation_map"]["independent_rederivation"]["flagged_by_prep_for_owner_confirmation"].items()} \
        == {"Platelet": "OUT", "ASDC": "OUT", "NK Proliferating": "NK", "Doublet": "unscored"}
    assert add["contamination"]["label"] == e6.CONTAMINATION_LABEL
    assert add["panels"]["remaining_E6_primary"] == {"B": ["CD20", "CD22", "CD268"], "T": ["CD3", "CD2"],
                                                     "NK": ["CD122", "CD56"], "myeloid": ["CD172a", "CD11c", "CD62P"]}
    sha = vk.sha256_file(ADDENDUM)
    assert f"{sha}  E6.json" in (ADDENDUM.parent / "HASHES.txt").read_text().splitlines()
    out = e6.annotation_class(["B naive", "Doublet", "Platelet"], add["annotation_map"]["map"])
    assert out.tolist() == ["B", "unscored", "OUT"]
    with pytest.raises(KeyError):
        e6.annotation_class(["not a label"], add["annotation_map"]["map"])


# ----------------------------------------------------------------------------- panels, scores, ANM closed form
def test_reduced_panel_and_recoded_field(reg):
    names = list(reg["evidence"]["teddy_head"]["adt_names"])
    rng = np.random.default_rng(3)
    v = rng.random((300, len(names)))
    pan = vk.question_panel(reg, "Q1")
    fld, rc = reg["anm"]["field_representation"], reg["anm"]["closure_readout"]
    np.testing.assert_allclose(e6.anm_recoded_scores(v, names, pan, fld), vk.anm_action_scores(v, names, reg, "Q1"), rtol=0, atol=1e-15)
    np.testing.assert_allclose(e6.anm_recoded_scores(v, names, pan, fld, rc),
                               vk.anm_action_scores(v, names, reg, "Q1", readout="closure"), rtol=0, atol=1e-15)
    np.testing.assert_allclose(e6.class_scores_panel(v, names, pan), vk.class_scores(v, names, reg, "Q1"), atol=1e-15)
    red = e6.reduced_panel(pan, [p for p in names if p not in ("CD5", "CD94")])
    S = e6.class_scores_panel(v, names, red)
    A = e6.anm_recoded_scores(v, names, red, fld)
    for k, lin in enumerate(vk.LINEAGES):
        n = len(red[lin])
        np.testing.assert_allclose(A[:, k], vk.field_gain(n, fld) * n * S[:, k], rtol=1e-12)
    with pytest.raises(ValueError):
        e6.reduced_panel(pan, ["CD3"])


@pytest.mark.skipif(not (ANM_ROOT / "active_neural_matter").exists(), reason="ANM_ROOT not available")
def test_engine_equals_recoded_closed_form_on_unequal_panels(reg):
    names = list(reg["evidence"]["teddy_head"]["adt_names"])
    j = {n: i for i, n in enumerate(names)}
    bridge = e1.AnmBridge(reg, REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json", ANM_ROOT)
    pan = e6.reduced_panel(vk.question_panel(reg, "Q1"), [p for p in names if p not in ("CD5", "CD94")])
    fld, rc = reg["anm"]["field_representation"], reg["anm"]["closure_readout"]
    rng = np.random.default_rng(4)
    v = np.clip(rng.random((60, len(names))) * 1.2, 0, 1)
    rec = e6.anm_recoded_scores(v, names, pan, fld)
    recc = e6.anm_recoded_scores(v, names, pan, fld, rc)
    for i in range(v.shape[0]):
        o = bridge.run(pan, {k: [v[i, j[p]] for p in pan[k]] for k in vk.LINEAGES}, readouts=("Q1",),
                       closure_readouts=("closure_Q1",))
        np.testing.assert_allclose(o["Q1"][0], rec[i], rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(o["closure_Q1"][0], recc[i], rtol=1e-12, atol=1e-14)
        assert e6.engine_calls([o["Q1"][1]])[0] == e6.calls_at(rec[i:i + 1], bridge.thresholds["Q1"])[0]


# ----------------------------------------------------------------------------- evaluator and replication
def test_evaluator_matches_materialised_replicates(reg):
    amend = reg["amendment_A1"]
    rng = np.random.default_rng(5)
    n = 400
    ids = np.sort(rng.choice(40000, n, replace=False))
    S = np.round(rng.random((n, 4)), 2)
    key = rng.choice(["B", "T", "NK", "myeloid", "OUT", "unscored"], n)
    grid = [0.9, 0.7, 0.5]
    m = e6.Method("rule", S.max(axis=1), S.argmax(axis=1), ids, amend, ["g0.90", "g0.70", "g0.50", "c1"], ["k"], ["c1", "g0.70"])
    called = (S.max(axis=1) >= 0.6).astype(float)
    ev = e6.Evaluator([m], {"k": key}, grid, {"c1": called}, ["k"])
    calls = np.asarray(vk.LINEAGES)[S.argmax(axis=1)]
    for rep in range(15):
        w = np.bincount(rng.integers(0, n, n), minlength=n) if rep else np.ones(n, dtype=np.int64)
        st = ev.stats(w)
        pos = np.repeat(np.arange(n), w)
        c1 = called[pos].mean()
        for name, c in (("g0.90", 0.9), ("g0.70", 0.7), ("c1", c1)):
            sel = va.coverage_select(S.max(axis=1)[pos], ids[pos], c, amend)
            kk, cc = key[pos][sel], calls[pos][sel]
            sc = kk != "unscored"
            want = np.mean(cc[sc] == kk[sc]) if sc.any() else np.nan
            assert np.isclose(st[f"rule|{name}|k|acc"], want, equal_nan=True)
            if name in ("c1", "g0.70"):
                n_out = np.sum(key[pos] == "OUT")
                assert np.isclose(st[f"rule|{name}|k|outdecl"], 1 - np.sum(kk == "OUT") / n_out)
                ins = np.isin(kk, vk.LINEAGES)
                assert np.isclose(st[f"rule|{name}|k|inscope"], np.mean(cc[ins] == kk[ins]))


def test_replication_rule():
    pd = {f"d{i}": v for i, v in enumerate([0.1, 0.2, 0.05, 0.3, 0.01, -0.1, -0.2, -0.05])}
    r = e6.replication(0.02, 0.1, [0.01, 0.2], pd)
    assert r["replicates"] and r["n_donors_same_sign"] == 5
    assert not e6.replication(0.02, 0.1, [-0.01, 0.2], pd)["replicates"]                # interval includes 0
    pd4 = dict(pd, d4=-0.01)
    assert not e6.replication(0.02, 0.1, [0.01, 0.2], pd4)["replicates"]               # 4 of 8 is no majority
    assert e6.replication(-0.02, -0.1, [-0.2, -0.01], {k: -v for k, v in pd.items()})["replicates"]
    assert not e6.replication(-0.02, 0.1, [0.01, 0.2], pd)["replicates"]                # opposite direction
    assert e6.replication(0.0, 0.1, [0.01, 0.2], pd)["status"] == "no v3 direction"
    assert e6.replication(None, 0.1, [0.01, 0.2], pd)["status"] == "no v3 direction"
    pdn = dict(pd, d0=None, d1=None)
    assert not e6.replication(0.02, 0.1, [0.01, 0.2], pdn)["replicates"]               # undefined donors count against


# ----------------------------------------------------------------------------- guards
def test_guard_refuses_incomplete_embedding_and_uncommitted_addendum(tmp_path, add):
    R = _runner()
    st = {"addendum_hash_recorded": True, "addendum_names_these_files": True,
          "registration_and_amendments_committed": True, "addendum_committed": False, "code_committed_clean": True}
    with pytest.raises(SystemExit, match="commit the E6 addendum"):
        R.guard_registration(st, smoke=False)
    R.guard_registration(st, smoke=True)
    st_code = dict(st, addendum_committed=True, code_committed_clean=False)
    with pytest.raises(SystemExit, match="runner and libraries committed"):
        R.guard_registration(st_code, smoke=False)
    R.guard_registration(dict(st_code, code_committed_clean=True), smoke=False)
    ext, emb = tmp_path / "ext", tmp_path / "emb"
    (emb / "shards").mkdir(parents=True)
    ext.mkdir()
    for k in range(4):
        np.savez(emb / "shards" / f"shard_{k:04d}.npz", cfg_hash=np.asarray("x"))
    a = R.parse_args(["--external", str(ext), "--external-embedding", str(emb)])
    with pytest.raises(SystemExit) as exc:
        R.guard_external(a, add)
    msg = str(exc.value)
    assert "pack sha256" in msg and "4 of 8 shards" in msg and "z_rna.npy" in msg


def test_guard_refuses_frozen_inputs_other_than_the_addendums(tmp_path, add):
    R = _runner()
    proc = tmp_path / "proc"
    proc.mkdir()
    for p in (proc / "cite_arrays.npz", tmp_path / "z.npy", tmp_path / "best.pt"):
        p.write_bytes(b"not the registered file")
    a = R.parse_args(["--processed", str(proc), "--z", str(tmp_path / "z.npy"), "--ckpt", str(tmp_path / "best.pt")])
    with pytest.raises(SystemExit, match="frozen inputs differ") as exc:
        R.guard_frozen_inputs(a, add)
    assert all(k in str(exc.value) for k in ("ckpt_sha256", "bmmc_cite_arrays_sha256", "bmmc_z_sha256"))
    if HAVE_BMMC:
        got = R.guard_frozen_inputs(R.parse_args([]), add)
        assert got == {k: add["provenance"]["inputs"][k] for k in R.FROZEN_INPUTS}


def test_n_boot_override_refused_outside_smoke():
    R = _runner()
    with pytest.raises(SystemExit):
        R.parse_args(["--n-boot", "10"])


# ----------------------------------------------------------------------------- val through the E6 code path
@pytest.mark.skipif(not (HAVE_BMMC and (ANM_ROOT / "active_neural_matter").exists()), reason="BMMC data or ANM_ROOT missing")
def test_val_smoke_through_e6_path(tmp_path):
    R = _runner()
    out = tmp_path / "e6_smoke"
    rc = R.main(["--smoke", "--stage", "all", "--out-dir", str(out), "--n-boot", "12", "--chunk", "6",
                 "--anm-root", str(ANM_ROOT), "--external", str(tmp_path / "no_pack"),
                 "--external-embedding", str(tmp_path / "no_embedding")])
    assert rc == 0
    res = json.loads((out / "E6_results.json").read_text())
    kv = json.loads((out / "key_validity.json").read_text())
    assert res["smoke"] and res["contamination_label"] == e6.CONTAMINATION_LABEL
    assert res["state"]["frozen_inputs_sha256"]["ckpt_sha256"] == json.loads(ADDENDUM.read_text())["provenance"]["inputs"]["ckpt_sha256"]
    assert kv["gate_spec"]["absent_gate_proteins"] == ["CD33", "CD94"]
    assert res["E1.1a"]["bridge_passed"] and res["E1.1a"]["full_panel_passed"]
    assert res["consistency_checks"]["rule_Q1_at_c_star_equals_operating_point"]
    assert set(res["cells"]["per_donor"]) == {f"valP{i}" for i in range(1, 9)}
    assert sum(res["cells"]["per_donor"].values()) == 6106
    rep = res["summary"]["replication"]
    assert "E1.4a anm - margin AURC" in rep and "E1.C2 closure - rule OUT decline @c*" in rep
    txt = (out / "REPORT.md").read_text()
    assert "SMOKE TEST" in txt and e6.CONTAMINATION_LABEL in txt and txt.index("Key validity") < txt.index("Replication")
    # resume: a second invocation reuses every cache
    rc = R.main(["--smoke", "--stage", "all", "--out-dir", str(out), "--n-boot", "12", "--chunk", "6",
                 "--anm-root", str(ANM_ROOT)])
    log = (out / "progress.log").read_text()
    assert rc == 0 and "prepare: cached" in log and "anm: cached" in log and "boot: cached" in log
