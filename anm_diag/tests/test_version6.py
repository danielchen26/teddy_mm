"""What batch 3 showed to be missing (LOOP_PROTOCOL.md, the candidate for version 6): the arithmetic of a readout that is a
ratio, a part of a direct part that is known by declaration, and the width of an interval from the recorded pairs."""
import numpy as np

from anm_diag import account as AC
from anm_diag import events as EV


def _world(seed, n_events=30, own=4, K=8, p=12, n=150, nc=1500, plain=0.03, far=0.0, big=-1.2):
    """Cells with a stage of p dimensions and K columns of log counts that the stage partly tells. The readouts are the
    columns minus their mean (a ratio over K columns). The first `own` events each lower one column by `big`; every pair
    has a small direct part of its own, and a share `far` of them a large one."""
    rng = np.random.default_rng(seed)
    Bm = rng.normal(0, 0.25, (p, K))
    names = ["control"] + [f"e{i:02d}" for i in range(n_events)]
    xs, Ls, un = [], [], []
    for i, u in enumerate(names):
        m = nc if u == "control" else n
        x = rng.normal(0, 1, (m, p)) + (0 if u == "control" else rng.normal(0, 0.15, p))
        L = 3 + x @ Bm + rng.normal(0, 0.4, (m, K))
        if u != "control":
            L += rng.normal(0, plain, K) * 0.4
            L += np.where(rng.random(K) < far, rng.choice([-1, 1], K) * rng.uniform(0.5, 1.5, K), 0.0) * 0.4
            if i - 1 < own:
                L[:, i - 1] += big
        xs.append(x); Ls.append(L); un += [u] * m
    x, L, un = np.vstack(xs), np.vstack(Ls), np.array(un)
    clr = L - L.mean(1, keepdims=True)
    Y = {f"r{j}": clr[:, j] for j in range(K)}
    edges = {(f"e{i:02d}", f"r{i}"): "own" for i in range(own)}
    return x, Y, un, edges, {r: ("the panel", K) for r in Y}


def _parts(rows):
    return {(r["event"], r["readout"]): (r["direct"], r["direct_sd"]) for r in rows if r["edge"] is not None and r["no_forecast"] is None}


def test_the_arithmetic_of_a_ratio_is_exact():
    rng = np.random.default_rng(0)
    K, L = 6, rng.normal(3, 0.5, (400, 6))
    clr = lambda M: M - M.mean(1, keepdims=True)
    sd = {f"r{j}": float(clr(L)[:, j].std()) for j in range(K)}
    over = {f"r{j}": ("the panel", K) for j in range(K)}
    for lowered in ({0: -1.0}, {0: -1.0, 3: 0.4}):                  # one column changes; two columns change
        M = L.copy()
        for j, d in lowered.items():
            M[:, j] += d
        change = {f"r{j}": float((clr(M)[:, j] - clr(L)[:, j]).mean() / sd[f"r{j}"]) for j in range(K)}
        direct = {("e", f"r{j}"): (change[f"r{j}"], 0.01) for j in lowered}
        got = EV.derived(direct, sd, over)
        assert set(got) == {("e", f"r{j}") for j in range(K) if j not in lowered}
        for (_e, r), (v, s) in got.items():
            assert abs(v - change[r]) < 1e-12 and s > 0             # what the other readouts move by, with nothing fitted
    assert EV.derived({("e", "r0"): (-2.0, 0.0)}, sd, {}) == {}      # a readout that is no ratio moves nothing else


def test_a_part_known_by_declaration_is_taken_out_of_a_row_and_put_back_into_its_forecast():
    x, Y, un, edges, over = _world(1)
    first = AC.Law(x, Y, un, edges=edges, types={}).fit()
    known = EV.derived(_parts(first.forecast()), {r: q["sd"] for r, q in first.r.items()}, over)
    law = AC.Law(x, Y, un, edges=edges, types={}, known={k: v for k, (v, _s) in known.items()}).fit()
    owners = {e for e, _r in edges}
    far = lambda L: np.median([abs(r["z"]) for r in L.forecast() if r["event"] in owners and r["edge"] is None])
    assert far(first) > 2.5 and far(law) < 1.2                      # the other readouts of an event that lowers one column: far out, then where they belong
    row = next(r for r in law.forecast() if r["event"] == "e00" and r["readout"] == "r5")
    assert abs(row["known_by_declaration"] - known[("e00", "r5")][0]) < 1e-12 and row["known_by_declaration"] > 0.1
    assert abs(row["forecast"] - (row["read"] + row["from_the_other_events"] + row["known_by_declaration"])) < 1e-9
    assert len(law.flags()) < len(first.flags())                    # and they are no longer flagged as edges the graph lacks
    m = un == "e00"
    new = {r["readout"]: r for r in law.new(x[m], np.full(int(m.sum()), "n"), known={("n", "r5"): (0.3, 0.04)})}
    assert abs(new["r5"]["known_by_declaration"] - 0.3) < 1e-12 and abs(new["r5"]["forecast"] - new["r5"]["read"] - 0.3) < 1e-9
    assert new["r5"]["sd"] > new["r4"]["sd"] - 1e-9 and new["r4"]["known_by_declaration"] == 0.0


def test_the_width_of_an_interval_comes_from_the_recorded_pairs():
    x, Y, un, _edges, _over = _world(2, n_events=60, own=0, far=0.03)
    events = sorted(set(un.tolist()) - {"control"})
    rec, new = np.isin(un, ["control"] + events[:30]), np.isin(un, events[30:])
    law = AC.Law(x[rec], {r: y[rec] for r, y in Y.items()}, un[rec], types={}).fit()
    w = law.conformal()
    assert w["no type"]["pairs"] == 240 and w["no type"]["events"] == 30 and w["of a type"]["quantile"] is None
    assert 1.96 < w["no type"]["quantile"] < 6                       # a few pairs stand far out: the interval is wider than a normal's
    inside = lambda rows: np.mean([r["interval"][0] <= law.measure(Y[r["readout"]][new], un[new], r["readout"])[r["event"]][0] <= r["interval"][1] for r in rows])
    plain, wide = law.new(x[new], un[new]), law.new(x[new], un[new], width=w)
    assert all(r["width_from"].startswith("the 240 recorded pairs") for r in wide) and all(r["width_from"].startswith("1.96") for r in plain)
    assert inside(wide) >= 0.93 and inside(wide) > inside(plain)
    few = np.isin(un, ["control"] + events[:5])                      # five events: too few to say how wide; the law's own mixture, and the row says so
    small = AC.Law(x[few], {r: y[few] for r, y in Y.items()}, un[few], types={}).fit()
    assert small.conformal()["no type"]["quantile"] is None
    row = small.new(x[new], un[new], width=small.conformal())[0]
    assert row["width_from"].startswith("assumed") and (row["interval"][1] - row["interval"][0]) / 2 > 1.96 * row["sd"]


def test_the_central_interval_of_a_mixture():
    assert abs(AC._mixture_half(0.04, 1.04, 0.0) - 1.96 * 0.2) < 1e-3
    c = AC._mixture_half(0.04, 1.04, 0.03)
    rng = np.random.default_rng(0)
    d = np.where(rng.random(400000) < 0.03, rng.normal(0, np.sqrt(1.04), 400000), rng.normal(0, 0.2, 400000))
    assert c > 1.96 * 0.2 and abs(np.mean(np.abs(d) > c) - 0.05) < 0.003
