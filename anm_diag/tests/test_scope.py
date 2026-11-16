"""Ground truth for the scope axis: six planted situations, and the reading must name each one.

Every case has one treated unit g, control cells, three other treated units of the same kind, and a separate set of
training cells. The readout responds to the event by delta. What differs is where a readout could have learned it.
"""
import numpy as np
import pytest

from anm_diag import scope as sc

D, C, S, DELTA = 8, 1.0, 2.0, 2.0          # state dimension, effect of z1 on the readout, shift of z1 in a treated unit


def make(case, seed=0, n_t=400, n_c=2000, n_o=300, n_train=3000, R=1):
    rng = np.random.default_rng(seed)
    units = ["g"] * n_t + ["h1"] * n_o + ["h2"] * n_o + ["h3"] * n_o
    n = len(units) + n_c
    block = np.array(units + list(sc.fold_blocks(n_c, 5, seed)))
    treated_any = np.array([True] * len(units) + [False] * n_c)
    is_g = block == "g"
    z = rng.normal(size=(n, D))
    x = rng.normal(size=(n, D))                               # an input view
    y = z[:, 2] + 0.1 * rng.normal(size=n)
    zt = rng.normal(size=(n_train, D))
    yt = zt[:, 2] + 0.1 * rng.normal(size=n_train)
    frozen = z[:, 2].copy()
    if case in ("none", "function", "coverage"):
        z[treated_any, 1] += S                                # the state carries the event
        y += C * z[:, 1]
        if case != "coverage":
            yt += C * zt[:, 1]                                # the training cells show how z1 moves the readout
        if case == "none":
            frozen = z[:, 2] + C * z[:, 1]
    elif case == "unit":
        z[:, 3:7] *= 0.05
        for j, u_ in enumerate(["g", "h1", "h2", "h3"]):      # each unit has its own mark in the state and its own shift
            z[block == u_, 3 + j] += 1.0
            y[block == u_] += DELTA * (1 if j % 2 == 0 else -1)
    elif case == "upstream":
        x[treated_any, 1] += S                                # the input view carries the event, the state does not
        y += C * x[:, 1]
    elif case == "unreadable":
        y[treated_any] += DELTA
    else:
        raise ValueError(case)
    if R > 1:                                                 # more readouts: the same response at different sizes
        w = np.linspace(1.0, 0.5, R)
        resp = y - z[:, 2]
        y = np.column_stack([z[:, 2] + w[j] * resp + 0.05 * rng.normal(size=n) for j in range(R)])
        yt = np.column_stack([zt[:, 2] + w[j] * (yt - zt[:, 2]) for j in range(R)])
        frozen = np.column_stack([z[:, 2] + w[j] * (frozen - z[:, 2]) for j in range(R)])
    study = sc.Study(y, {"X": x, "Z": z}, block, train={"stages": {"X": rng.normal(size=(n_train, D)), "Z": zt}, "Y": yt, "block": None},
                     frozen=frozen, seed=seed)
    return study, is_g, ~treated_any


DEC = {"bootstrap": 300, "seed": 0}


@pytest.mark.parametrize("case", ["none", "function", "coverage", "unit", "upstream", "unreadable"])
def test_the_reading_names_what_was_planted(case):
    study, t, c = make(case)
    r = sc.scope_grid(study, t, c, reads="Z", input_view="X", decision=DEC)
    assert r["delta"]["resolved"]
    assert r["targets"]["shift"]["reading"]["missing"] == case, r["targets"]["shift"]["reading"]


def test_other_seeds_agree():
    for case in ("function", "coverage", "unit"):
        for seed in (1, 2):
            study, t, c = make(case, seed=seed)
            assert sc.scope_grid(study, t, c, reads="Z", input_view="X", decision=DEC)["targets"]["shift"]["reading"]["missing"] == case


def test_u_is_near_zero_and_one_where_it_should_be():
    study, t, c = make("coverage")
    S_ = sc.scope_grid(study, t, c, reads="Z", decision=DEC)["targets"]["shift"]["stages"]
    assert abs(S_["system"]["frozen"]["u"] - 1.0) < 0.1          # the frozen output ignores z1
    assert abs(S_["Z"]["train"]["u"] - 1.0) < 0.1                # the training cells never showed the dependence
    assert abs(S_["Z"]["units"]["u"]) < 0.1 and abs(S_["Z"]["own"]["u"]) < 0.1


def test_forecast_follows_the_gaps():
    study, t, c = make("coverage")
    f = sc.forecast(sc.scope_grid(study, t, c, reads="Z", decision=DEC))
    assert f["coverage"]["removes"] > 0.8 and abs(f["function"]["removes"]) < 0.15
    study, t, c = make("function")
    f = sc.forecast(sc.scope_grid(study, t, c, reads="Z", decision=DEC))
    assert f["function"]["removes"] > 0.8


def test_prediction_in_other_units_is_put_on_the_measured_scale():
    """A frozen output that is right up to a scale and an offset carries the response."""
    study, t, c = make("none")
    study.frozen = 0.2 * study.frozen - 7.0
    r = sc.scope_grid(study, t, c, reads="Z", decision=DEC)
    assert r["targets"]["shift"]["stages"]["system"]["frozen"]["verdict"] == "support"


def test_several_readouts_and_the_variation_target():
    study, t, c = make("coverage", R=4)
    r = sc.scope_grid(study, t, c, reads="Z", decision=DEC, variation=True)
    sh = r["targets"]["shift"]
    assert sh["reading"]["missing"] == "coverage"
    assert len(sh["stages"]["Z"]["units"]["by_class"]["ridge"]["per_readout"]) == 4
    va = r["targets"]["variation"]["stages"]
    assert va["Z"]["own"]["u"] < 0.2 and va["Z"]["own"]["u"] <= va["Z"]["train"]["u"] + 0.02


def test_no_cell_of_the_scored_unit_enters_a_units_fit():
    study, t, c = make("unit")
    study.predict("units", "Z", "ridge", np.where(t)[0], treated_blocks=("g",))
    assert all("g" in excl for (where, _, _, excl) in study._models if where == "units")


def test_a_response_that_is_not_resolved_is_not_read():
    study, t, c = make("unreadable")
    study.Y[:] = np.random.default_rng(3).normal(size=study.Y.shape)      # no response at all
    r = sc.scope_grid(study, t, c, reads="Z", decision=DEC)
    assert not r["delta"]["resolved"] and r["targets"]["shift"]["reading"]["missing"] == "not testable"


def test_family_counts_events():
    rs = []
    for seed in range(4):
        study, t, c = make("coverage", seed=seed)
        rs.append(sc.scope_grid(study, t, c, reads="Z", decision=DEC))
    fam = sc.family(rs)
    assert fam["n_events"] == 4 and fam["approximate"] and fam["missing"] == {"coverage": 4}
    assert fam["gaps"]["train->units"]["resolved_in"] == 4
