"""The reading script of runs/round1_coverage_20261116: which head may score which donor, and the checks of gate 3."""
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "runs" / "round1_coverage_20261116" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import common  # noqa: E402
import plan  # noqa: E402

pytest.importorskip("torch")
import read  # noqa: E402


def test_no_head_scores_a_donor_it_was_fitted_on():
    by_name = {c["name"]: c for c in plan.heads()}
    used = set()
    for s in plan.SEEDS:
        for d in plan.DONORS:
            H = read.heads_for(s, d)
            assert set(H) == set(read.ARMS)
            for label, name in H.items():
                c = by_name[name]
                used.add(name)
                assert c["seed"] == s and d not in c["pbmc"], (label, name, d)
            seen = {plan.BATCH_OF[x] for x in by_name[H["coverage_seen"]]["pbmc"]}
            unseen = {plan.BATCH_OF[x] for x in by_name[H["coverage_not_seen"]]["pbmc"]}
            assert plan.BATCH_OF[d] in seen and plan.BATCH_OF[d] not in unseen
            own, = by_name[H["one_donor_own_batch"]]["pbmc"]
            other, = by_name[H["one_donor_other_batch"]]["pbmc"]
            assert plan.BATCH_OF[own] == plan.BATCH_OF[d] != plan.BATCH_OF[other]
    assert used == set(by_name)                       # every head of the plan scores some donor


def test_ran_without_falling():
    h = lambda v: [{"val_mlp_pearson": x} for x in v]
    assert read.ran_without_falling(h([0.5] * 5 + [0.6] * 5))["pass"]
    assert read.ran_without_falling(h([0.6] * 5 + [0.596] * 5))["pass"]
    assert not read.ran_without_falling(h([0.6] * 5 + [0.59] * 5))["pass"]


def _toy():
    g = np.random.default_rng(0)
    n, P = 400, 6
    donors = np.repeat(plan.DONORS, n // 8)
    Y = g.normal(size=(n, P)).astype(np.float32)
    Y[:, -1] = np.nan
    pb = {"adt": Y, "donors": donors, "barcodes": np.array([f"c{i}" for i in range(n)])}
    nb = 300
    bm = {"adt": g.normal(size=(nb, P)).astype(np.float32), "split": np.array(["train"] * 200 + ["val"] * 50 + ["test"] * 50),
          "names": np.array([f"p{j}" for j in range(P)])}
    return g, pb, bm, np.arange(P - 1)


def test_gate3_checks():
    g, pb, bm, shared = _toy()
    cfg = {"arm": "coverage", "fit_set": "M1", "pbmc": {d: 50 for d in plan.fit_donors("M1")}}
    fit = np.where(np.isin(pb["donors"], plan.fit_donors("M1")))[0]
    cells = {"pbmc_index": fit, "bm_index": np.arange(200)}
    hist = {"history": [{"val_mlp_pearson": 0.6, "train_loss": 2.0}] * 10}
    Yp = pb["adt"][:, shared]
    good = {"pbmc": Yp + 0.5 * g.normal(size=Yp.shape), "bm": bm["adt"] + 0.5 * g.normal(size=bm["adt"].shape), "metrics": hist}
    base = {"pbmc": Yp + 3.0 * g.normal(size=Yp.shape), "bm": good["bm"], "metrics": hist}
    r = read.gate3(cfg, hist, cells, good, base, bm, pb, shared)
    assert r["admitted"] and r["own_objective"]["gain"] > read.OWN_OBJECTIVE and r["guard_validation"]["change"] == 0.0
    assert not read.gate3(cfg, hist, cells, base, base, bm, pb, shared)["own_objective"]["pass"]          # learned nothing of its own cells
    leak = {"pbmc_index": np.r_[fit, np.where(pb["donors"] == plan.scored_donors("M1")[0])[0][:1]], "bm_index": np.arange(200)}
    assert not read.gate3(cfg, hist, leak, good, base, bm, pb, shared)["no_scored_cell_in_fit"]["pass"]      # a scored donor's cell in the fit
    assert not read.gate3(cfg, hist, {"pbmc_index": fit, "bm_index": np.arange(190, 210)}, good, base, bm, pb, shared)["no_scored_cell_in_fit"]["pass"]
    worse = {**good, "bm": bm["adt"] + 5.0 * g.normal(size=bm["adt"].shape)}
    assert not read.gate3(cfg, hist, cells, worse, base, bm, pb, shared)["guard_validation"]["pass"]
    f = {"arm": "function", "fit_set": None, "pbmc": {}}
    none = {"pbmc_index": np.zeros(0, int), "bm_index": np.arange(200)}
    lower = {"history": [{"val_mlp_pearson": 0.6, "train_loss": 1.9}] * 10}
    assert read.gate3(f, lower, none, good, base, bm, pb, shared)["admitted"]
    assert not read.gate3(f, hist, none, good, base, bm, pb, shared)["own_objective"]["pass"]
    g6 = read.guard6(bm, worse, base, shared)
    assert not g6["pass"] and g6["worst_protein"]["change"] < 0 and np.isfinite(g6["proteins_the_second_study_lacks"])


def test_head_forward_is_the_torch_head():
    torch = pytest.importorskip("torch")
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from teddy_mm.models import AdtDecoder, MLP
    torch.manual_seed(0)
    for hidden in (16, 40):
        mlp, dec = MLP(12, 12, hidden=hidden).eval(), AdtDecoder(12, 5, hidden=hidden).eval()
        z = torch.nn.functional.normalize(torch.randn(30, 12), dim=-1)
        W = {f"mlp__{k}": v.numpy() for k, v in mlp.state_dict().items()}
        W.update({f"dec__{k}": v.numpy() for k, v in dec.state_dict().items()})
        with torch.no_grad():
            ref = dec(mlp(z), torch.full((30,), common.SF))[0].numpy()
        assert np.allclose(common.head_forward(W, z.numpy()), ref, atol=1e-5)
