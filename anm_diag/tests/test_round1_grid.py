"""D6 of runs/round1_coverage_20261116: the grid on the changed head, scored donors only (scripts/grid_changed.py)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

RUNS = Path(__file__).resolve().parents[1] / "runs"
SCRIPTS = RUNS / "round1_coverage_20261116" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import common  # noqa: E402
import plan  # noqa: E402

pytest.importorskip("torch")                          # common.load_studies imports the trainer; nothing here trains
import grid_changed as gc  # noqa: E402
from anm_diag import scope  # noqa: E402

ZD, HID, NA, RD = 12, 8, 6, 10
SHARED = np.array([0, 1, 2, 4, 5])


def weights(seed):
    r = np.random.default_rng(seed)
    W = {}
    for part, (i, o) in (("mlp__net", (ZD, ZD)), ("dec__mu.net", (ZD, NA))):
        for k, (a, b) in zip((0, 3, 6), ((HID, i), (HID, HID), (o, HID))):
            W[f"{part}.{k}.weight"], W[f"{part}.{k}.bias"] = r.normal(size=(a, b)).astype(np.float32), r.normal(size=a).astype(np.float32)
    return W


def planted(seed=0, per_donor=60, n_bm=400):
    """Two tiny studies in the layout of the work arrays: bone-marrow training cells, and eight donors."""
    r = np.random.default_rng(seed)
    split = np.array(["train"] * n_bm + ["val"] * 40 + ["test"] * 40)
    r.shuffle(split)
    tr = np.where(split == "train")[0]
    A = r.normal(size=(ZD, len(SHARED)))
    zb, zh = r.normal(size=(n_bm, ZD)).astype(np.float32), r.normal(size=(8 * per_donor, ZD)).astype(np.float32)
    wb = dict(RNA512=r.normal(size=(n_bm, RD)).astype(np.float32), Z=zb, Y=(zb @ A + r.normal(size=(n_bm, len(SHARED)))).astype(np.float32),
              donor=r.choice(["d1", "d2", "d3"], n_bm))
    wh = dict(RNA512=r.normal(size=(len(zh), RD)).astype(np.float32), Z=zh, Y=(zh @ A + 1.0 + r.normal(size=(len(zh), len(SHARED)))).astype(np.float32),
              donor=np.repeat(plan.DONORS, per_donor))
    return split, tr, wb, wh


def test_the_head_is_the_head_the_reading_script_runs():
    W = weights(0)
    z = np.random.default_rng(1).normal(size=(50, ZD)).astype(np.float32)
    lat, out, mu = gc.head_stages(W, z)
    unit = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)
    np.testing.assert_allclose(mu, common.head_forward(W, unit), rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(np.exp(out) * common.SF, mu, rtol=1e-4, atol=1e-5)
    assert lat.shape == (50, ZD) and out.shape == (50, NA)


def test_the_cells_of_a_run():
    split, tr, wb, wh = planted()
    r = np.random.default_rng(2)
    bm_index = np.sort(r.choice(tr, 250, replace=False))
    fit = np.where(np.isin(wh["donor"], plan.fit_donors("M1")))[0]
    pbmc_index = np.sort(r.choice(fit, 100, replace=False))
    cells = {"bm_index": bm_index, "pbmc_index": pbmc_index}
    scored = plan.scored_donors("M1")
    audit, train = gc.run_arrays(weights(0), cells, scored, wb, wh, split, SHARED)
    assert sorted(set(audit["donor"])) == sorted(scored) and len(audit["Y"]) == 4 * 60
    assert not set(train["donor"]) & set(scored)
    pos = np.searchsorted(tr, bm_index)
    np.testing.assert_array_equal(train["Z"], np.vstack([wb["Z"][pos], wh["Z"][pbmc_index]]))
    np.testing.assert_array_equal(train["Y"], np.vstack([wb["Y"][pos], wh["Y"][pbmc_index]]))
    assert audit["frozen"].shape == audit["Y"].shape and train["HEAD_LAT"].shape == (350, ZD) and train["HEAD_OUT"].shape == (350, NA)
    lat, _, mu = gc.head_stages(weights(0), wh["Z"])
    m = np.isin(wh["donor"], scored)
    np.testing.assert_array_equal(audit["frozen"], mu[m][:, SHARED])
    np.testing.assert_array_equal(train["HEAD_LAT"][250:], lat[pbmc_index])
    with pytest.raises(SystemExit):                    # a scored donor among the fit cells
        gc.run_arrays(weights(0), {"bm_index": bm_index, "pbmc_index": np.where(m)[0][:5]}, scored, wb, wh, split, SHARED)
    with pytest.raises(SystemExit):                    # a bone-marrow cell that is not a training cell
        gc.run_arrays(weights(0), {"bm_index": np.sort(np.r_[bm_index[:-1], np.where(split == "val")[0][:1]]), "pbmc_index": pbmc_index},
                      scored, wb, wh, split, SHARED)


def test_the_spec_is_that_of_donors_hao():
    old = yaml.safe_load((RUNS / "scope_20261116" / "results" / "donors_hao" / "SPEC.yaml").read_text())
    new = gc.spec("x")
    assert new["keys"].pop("train_block") == "donor"
    for s in (old, new):
        for k in ("name", "arrays", "keep_draws", "readout_name"):          # the last only names the readouts in the events
            s.pop(k, None)
    assert new == old


def test_runs_pair_a_changed_head_with_as_is_on_the_donors_it_did_not_see():
    by_name = {c["name"]: c for c in plan.heads()}
    todo = gc.runs(plan.SEEDS, ["M1", "M2"])
    assert len(todo) == 12
    for r in todo:
        c = by_name[r["weights"]]
        assert c["seed"] == r["seed"] and not set(c["pbmc"]) & set(r["scored"]) and len(r["scored"]) == 4
        assert c["arm"] == r["head"]
    for s in plan.SEEDS:                               # every donor is scored once per seed and head
        for h in gc.HEADS:
            assert sorted(d for r in todo if r["seed"] == s and r["head"] == h for d in r["scored"]) == sorted(plan.DONORS)


def test_end_to_end_on_planted_cells(tmp_path):
    split, tr, wb, wh = planted()
    r = np.random.default_rng(3)
    heads = tmp_path / "heads"
    for name, n_pb in (("s0_as_is", 0), ("s0_coverage_M1", 100), ("s0_coverage_M2", 100)):
        (heads / name).mkdir(parents=True)
        fs = name.split("_")[-1]
        fit = np.where(np.isin(wh["donor"], plan.fit_donors(fs)))[0] if n_pb else np.zeros(0, int)
        np.savez(heads / name / "head.npz", **weights(hash(name) % 1000), adt_names=np.array([f"p{i}" for i in range(NA)]))
        np.savez(heads / name / "fit_cells.npz", bm_index=tr if not n_pb else np.sort(r.choice(tr, 300, replace=False)),
                 pbmc_index=np.sort(r.choice(fit, n_pb, replace=False)) if n_pb else fit, pbmc_row_from=fit[:0])
    res, tmp = tmp_path / "res", tmp_path / "tmp"
    tmp.mkdir()
    todo = gc.runs([0], ["M1", "M2"])
    for run in todo:
        gc.one_run(run, heads, wb, wh, split, SHARED, tmp, res, decision={"bootstrap": 40})
    assert not list(tmp.glob("*.npz"))                 # the arrays are removed again
    results = {run["name"]: json.loads((res / run["name"] / "scope.json").read_text()) for run in todo}
    out = gc.summarise(results, [0], ["M1", "M2"])
    v = out["variation"]
    assert out["donors"] == sorted(plan.DONORS) and set(v["u"]) == set(gc.HEADS)
    # the frozen cell is the quantity read.py reads: u of the head's own output on a donor's cells
    for h, wname in (("as_is", lambda fs: "s0_as_is"), ("coverage", lambda fs: f"s0_coverage_{fs}")):
        us = []
        for fs in ("M1", "M2"):
            _, _, mu = gc.head_stages(dict(np.load(heads / wname(fs) / "head.npz")), wh["Z"])
            for d in plan.scored_donors(fs):
                m = wh["donor"] == d
                us.append(scope.arms(wh["Y"][m], {"x": mu[m][:, SHARED]}, base="x", decision={"bootstrap": 40})["arms"]["x"]["u"])
        assert v["u"][h]["frozen"] == pytest.approx(np.mean(np.abs(us)), rel=1e-6)
    fall = v["fall_from_as_is_to_coverage"]["frozen"]
    assert fall["gain"] == pytest.approx(v["u"]["as_is"]["frozen"] - v["u"]["coverage"]["frozen"], rel=1e-9, abs=1e-12)
    # stages and scopes the head does not enter are the same in the two runs of a pair
    for k in ("units@RNA512", "units@Z", "own@RNA512", "own@Z"):
        assert v["u"]["as_is"][k] == pytest.approx(v["u"]["coverage"][k], rel=1e-9)
    assert sum(v["readings"]["coverage"].values()) == 8 and set(v["steps_at_the_stage_the_head_reads"]["coverage"]) == {
        "frozen->train", "train->units", "units->own"}
