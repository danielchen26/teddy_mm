"""Several outputs of a system on the same cells: scope.arms and scope.arms_family, on planted data."""
import numpy as np
import pytest

from anm_diag import scope


def _unit(n=1500, r=6, seed=0):
    g = np.random.default_rng(seed)
    f = g.normal(size=(n, r))
    Y = f + 0.5 * g.normal(size=(n, r))
    return g, f, Y


def test_arms_orders_outputs_and_pairs_gains():
    g, f, Y = _unit()
    good, poor = f + 0.2 * g.normal(size=f.shape), 0.3 * f + g.normal(size=f.shape)
    out = scope.arms(Y, {"poor": poor, "good": good, "same": poor.copy(), "scaled": 3.0 * poor + 7.0}, base="poor",
                     decision={"bootstrap": 300})
    a = out["arms"]
    assert a["good"]["u"] < a["poor"]["u"] and a["good"]["gain_ci"][0] > 0
    assert a["poor"]["gain"] == 0.0 and a["same"]["gain"] == 0.0 and a["same"]["gain_ci"] == [0.0, 0.0]
    assert abs(a["scaled"]["gain"]) < 1e-12                       # u does not depend on the level or the scale of an output
    r2 = np.mean([np.corrcoef(Y[:, j], good[:, j])[0, 1] ** 2 for j in range(Y.shape[1])])
    assert abs(a["good"]["u"] - (1 - r2)) < 1e-9


def test_arms_agrees_with_the_grid():
    g, f, Y = _unit(n=600, r=4, seed=1)
    P = f + g.normal(size=f.shape)
    n2 = 600
    Yall = np.vstack([Y, g.normal(size=(n2, 4))])
    Pall = np.vstack([P, g.normal(size=(n2, 4))])
    block = np.array(["a"] * 600 + ["b"] * n2)
    st = scope.Study(Yall, {"Z": g.normal(size=(1200, 3))}, block, frozen=Pall)
    dec = {"bootstrap": 200, "seed": 3}
    grid = scope.scope_grid(st, block == "a", block == "b", scopes=["frozen"], variation=True, decision=dec)
    cell = grid["targets"]["variation"]["stages"]["system"]["frozen"]
    one = scope.arms(Y, {"head": P}, base="head", decision=dec)["arms"]["head"]
    assert one["u"] == pytest.approx(cell["u"], abs=1e-12) and one["ci"] == pytest.approx(cell["ci"], abs=1e-12)


def test_arms_refuses_a_missing_base_and_a_wrong_shape():
    _, f, Y = _unit(n=100)
    with pytest.raises(scope.ScopeError):
        scope.arms(Y, {"a": f}, base="b")
    with pytest.raises(scope.ScopeError):
        scope.arms(Y, {"a": f[:, :2]}, base="a")


def test_family_gain_and_interval():
    g = np.random.default_rng(0)
    units, fits = 12, 3
    base = 0.6 + 0.02 * g.normal(size=(units, 1)) + 0.005 * g.normal(size=(units, fits))
    arm = base - 0.05 + 0.005 * g.normal(size=(units, fits))
    nul = base + 0.005 * g.normal(size=(units, fits))
    out = scope.arms_family({"base": base, "arm": arm, "nul": nul, "two": arm[:, :2]}, base="base")
    a = out["arms"]
    assert a["arm"]["gain"] == pytest.approx(float(base.mean() - arm.mean()))
    assert a["arm"]["ci"][0] < 0.05 < a["arm"]["ci"][1] and a["arm"]["ci"][0] > 0.03
    assert a["nul"]["ci"][0] < 0 < a["nul"]["ci"][1]
    assert a["base"]["gain"] == 0.0 and len(a["arm"]["per_fit"]) == 3 and len(a["two"]["per_fit"]) == 2
    assert all(f["ci"][0] > 0 for f in a["arm"]["per_fit"]) and not out["approximate"]
    assert scope.arms_family({"base": base[:8], "arm": arm[:8]}, base="base")["approximate"]
