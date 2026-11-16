"""The trainer of runs/round1_coverage_20261116: it is the registered recipe when nothing is changed, a stopped fit goes on
with the same result, and the cells of an arm are what the design says. Planted data only."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "anm_diag" / "runs" / "round1_coverage_20261116" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import plan  # noqa: E402
import train_arm  # noqa: E402

N_TRAIN, ZD, P = 1300, 24, 9


def _packs(tmp, seed=0):
    g = np.random.default_rng(seed)
    n = N_TRAIN + 400
    z = g.normal(size=(n, ZD)).astype(np.float32)
    W = g.normal(size=(ZD, P))
    adt = np.log1p(np.exp(z @ W * 0.5 + g.normal(size=(n, P)) * 0.3)).astype(np.float32)
    tot = adt.sum(1)
    split = np.array(["train"] * N_TRAIN + ["val"] * 200 + ["test"] * 200)
    d = tmp / "bm"
    d.mkdir()
    names = np.array([f"p{j}" for j in range(P)]).astype("U64")
    np.savez(d / "cite_arrays.npz", rna_data=np.zeros(1, np.float32), rna_indices=np.zeros(1, np.int32), rna_indptr=np.r_[0, np.ones(n)].astype(np.int32),
             rna_shape=np.array([n, 3]), adt=adt, token_ids=np.arange(3), split=split, donors=np.array(["d"] * n), sites=np.array(["s"] * n),
             cell_types=np.array(["c"] * n), rna_size_factor=np.ones(n, np.float32), adt_size_factor=(tot / np.median(tot)).astype(np.float32),
             adt_names=names, rna_names=np.array(["g0", "g1", "g2"]))
    (d / "meta.json").write_text("{}")
    np.save(d / "z_rna.npy", z)
    # a second study: 8 donors x 60 cells, the last two proteins not measured
    m = 8 * 60
    ze = g.normal(size=(m, ZD)).astype(np.float32)
    ae = np.log1p(np.exp(ze @ W * 0.5 + 1.0 + g.normal(size=(m, P)) * 0.3)).astype(np.float32)
    ae[:, -2:] = np.nan
    e = tmp / "pb"
    e.mkdir()
    np.savez(e / "cite_arrays.npz", adt=ae, adt_names=names, donors=np.repeat(plan.DONORS, 60), barcodes=np.array([f"c{i}" for i in range(m)]))
    np.save(e / "z_rna.npy", ze)
    return d, e


def _data(d, e):
    return train_arm.load_data(d / "cite_arrays.npz", d / "z_rna.npy", e / "cite_arrays.npz", e / "z_rna.npy")


def _cfg(**kw):
    return {"name": "h", "seed": 5, "arm": "as_is", "n_bm": None, "pbmc": {}, "shuffle_labels": False, "hidden": 32, "epochs": 3, "fit_set": None, **kw}


def _same(a, b):
    A, B = torch.load(a, weights_only=False), torch.load(b, weights_only=False)
    return all(torch.equal(A[part][k], B[part][k]) for part in ("mlp", "v_net", "dec") for k in A[part]) and A["metrics"] == B["metrics"]


def test_as_is_is_the_registered_recipe(tmp_path):
    d, e = _packs(tmp_path)
    ref = tmp_path / "ref"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "04_train.py"), "--processed", str(d), "--out", str(ref), "--epochs", "3",
                    "--hidden", "32", "--device", "cpu", "--seed", "5"], check=True, capture_output=True)
    bm, pb = _data(d, e)
    assert train_arm.fit(_cfg(), bm, pb, tmp_path / "mine", log=lambda *a: None)
    assert _same(ref / "best.pt", tmp_path / "mine" / "best.pt")
    hist = json.loads((ref / "metrics.json").read_text())["history"]
    mine = json.loads((tmp_path / "mine" / "metrics.json").read_text())
    assert mine["history"] == hist and mine["n_fit"] == N_TRAIN and not mine["masked_loss"]
    W = np.load(tmp_path / "mine" / "head.npz")
    best = torch.load(ref / "best.pt", weights_only=False)
    assert np.array_equal(W["mlp__net.0.weight"], best["mlp"]["net.0.weight"].numpy()) and "dec__mu.net.6.bias" in W.files


@pytest.mark.parametrize("arm", ["as_is", "coverage"])
def test_a_stopped_fit_goes_on_with_the_same_result(tmp_path, arm):
    d, e = _packs(tmp_path)
    bm, pb = _data(d, e)
    kw = {} if arm == "as_is" else {"arm": "coverage", "n_bm": 1100, "pbmc": {x: 50 for x in plan.fit_donors("M1")}}
    quiet = lambda *a: None
    assert train_arm.fit(_cfg(**kw), bm, pb, tmp_path / "whole", log=quiet)
    assert not train_arm.fit(_cfg(**kw), bm, pb, tmp_path / "parts", max_epochs_this_call=1, log=quiet)
    assert (tmp_path / "parts" / "state.pt").exists()
    assert not train_arm.fit(_cfg(**kw), bm, pb, tmp_path / "parts", max_epochs_this_call=1, log=quiet)
    assert train_arm.fit(_cfg(**kw), bm, pb, tmp_path / "parts", log=quiet)
    assert _same(tmp_path / "whole" / "best.pt", tmp_path / "parts" / "best.pt")
    a, b = (json.loads((tmp_path / k / "metrics.json").read_text())["history"] for k in ("whole", "parts"))
    assert a == b and not (tmp_path / "parts" / "state.pt").exists()
    with pytest.raises(SystemExit):                                   # a state is not continued under another head
        train_arm.fit(_cfg(**kw), bm, pb, tmp_path / "x", max_epochs_this_call=1, log=quiet)
        train_arm.fit(_cfg(**{**kw, "seed": 6}), bm, pb, tmp_path / "x", log=quiet)


def test_the_loss_leaves_out_what_was_not_measured():
    from teddy_mm.models import nb_nll
    g = torch.Generator().manual_seed(0)
    x, mu, th = torch.rand(40, 7, generator=g) * 3, torch.rand(40, 7, generator=g) + 0.1, torch.rand(7, generator=g) + 0.5
    ones = torch.ones(40, 7)
    assert torch.allclose(train_arm.nb_nll_cells(x, mu, th, ones), nb_nll(x, mu, th), atol=1e-6)
    m = ones.clone()
    m[20:, 5:] = 0
    mu2 = mu.clone().requires_grad_(True)
    L = train_arm.nb_nll_cells(x, mu2, th, m)
    L.backward()
    assert torch.all(mu2.grad[20:, 5:] == 0) and torch.all(mu2.grad[:20] != 0)
    x2 = x.clone()
    x2[20:, 5:] = 99.0
    assert torch.equal(train_arm.nb_nll_cells(x2, mu, th, m), train_arm.nb_nll_cells(x, mu, th, m))
    # every cell has weight 1: the mean over cells of each cell's mean over its measured proteins
    per_cell = torch.stack([nb_nll(x[i:i + 1, m[i] > 0], mu[i:i + 1, m[i] > 0], th[m[i] > 0]) for i in range(40)])
    assert torch.allclose(L.detach(), per_cell.mean(), atol=1e-6)


def test_the_cells_of_an_arm(tmp_path):
    d, e = _packs(tmp_path)
    bm, pb = _data(d, e)
    fit = plan.fit_donors("M1")
    (z, y, sf, m), rec = train_arm.build_train(_cfg(n_bm=1100, pbmc={x: 50 for x in fit}), bm, pb)
    assert len(z) == 1100 + 200 and m.shape == y.shape and torch.all(m[:1100] == 1) and torch.all(m[1100:, -2:] == 0) and torch.all(m[1100:, :-2] == 1)
    assert not torch.isnan(y).any() and set(pb["donors"][rec["pbmc_index"]]) == set(fit)
    assert np.all(bm["split"][rec["bm_index"]] == "train") and len(set(rec["bm_index"])) == 1100
    assert float(np.median(sf[1100:].numpy())) == pytest.approx(1.0, abs=1e-6)
    assert np.array_equal(rec["pbmc_index"], rec["pbmc_row_from"])
    # nested: fewer bone-marrow cells are a subset; fewer cells of a donor are a subset
    (_, _, _, _), small = train_arm.build_train(_cfg(n_bm=900, pbmc={fit[0]: 20}), bm, pb)
    assert set(small["bm_index"]) <= set(rec["bm_index"]) and set(small["pbmc_index"]) <= set(rec["pbmc_index"])
    # all training cells: the pack's order, and no mask (the registered path)
    (z0, y0, sf0, m0), full = train_arm.build_train(_cfg(), bm, pb)
    assert m0 is None and np.array_equal(full["bm_index"], np.where(bm["split"] == "train")[0]) and len(full["pbmc_index"]) == 0
    # shuffled proteins: every cell gets the row of another cell of its own donor; a donor's rows are all still there
    (zs, ys, sfs, ms), sh = train_arm.build_train(_cfg(n_bm=1100, pbmc={x: 50 for x in fit}, shuffle_labels=True), bm, pb)
    assert np.array_equal(sh["pbmc_index"], rec["pbmc_index"]) and torch.equal(zs, z)
    assert np.array_equal(pb["donors"][sh["pbmc_index"]], pb["donors"][sh["pbmc_row_from"]])
    assert np.mean(sh["pbmc_index"] != sh["pbmc_row_from"]) > 0.9 and sorted(sh["pbmc_index"]) == sorted(sh["pbmc_row_from"])
    assert not torch.equal(ys, y) and torch.allclose(ys[1100:].sum(0), y[1100:].sum(0), atol=1e-3)
    row = int(np.where(sh["pbmc_index"] == sh["pbmc_row_from"][0])[0][0])
    assert torch.equal(ys[1100], y[1100 + row]) and sfs[1100] == sf[1100 + row]
    # another seed draws other cells
    (_, _, _, _), other = train_arm.build_train(_cfg(seed=6, n_bm=1100, pbmc={x: 50 for x in fit}), bm, pb)
    assert set(other["bm_index"]) != set(rec["bm_index"])


def test_the_plan_is_the_design():
    H = plan.heads()
    assert len(H) == 45 and len({c["name"] for c in H}) == 45
    count = lambda arm: sum(c["arm"] == arm for c in H)
    assert [count(a) for a in ("as_is", "amount", "function", "coverage", "labels", "one_donor", "four_donors")] == [3, 3, 3, 12, 6, 12, 6]
    for c in H:
        n_pb = sum(c["pbmc"].values())
        total = (67405 if c["n_bm"] is None else c["n_bm"]) + n_pb
        assert total == (47405 if c["arm"] == "amount" else 67405)
        assert (c["hidden"], c["epochs"]) == ((1024, 80) if c["arm"] == "function" else (512, 40))
        assert c["shuffle_labels"] == (c["arm"] == "labels")
        if c["fit_set"]:
            assert set(c["pbmc"]) <= set(plan.fit_donors(c["fit_set"])) and not set(c["pbmc"]) & set(plan.scored_donors(c["fit_set"]))
            assert len(c["pbmc"]) == (1 if c["arm"] == "one_donor" else 4)
    for fs, seen in plan.BATCH_SEEN.items():
        fit_batches = {plan.BATCH_OF[x] for x in plan.fit_donors(fs)}
        assert all((plan.BATCH_OF[x] in fit_batches) == seen for x in plan.scored_donors(fs)) and len(plan.scored_donors(fs)) == 4
    assert sorted(x for fs in ("M1", "M2") for x in plan.scored_donors(fs)) == sorted(plan.DONORS)
    assert sorted(x for fs in ("B1", "B2") for x in plan.scored_donors(fs)) == sorted(plan.DONORS)
    F = plan.heads(final=True)
    assert len(F) == 6 and all(sum(c["pbmc"].values()) == 20000 and len(c["pbmc"]) == 8 for c in F)
    assert plan.head_by_name("s1_coverage_B2")["seed"] == 1 and plan.head_by_name("s2_labels_all")["shuffle_labels"]
