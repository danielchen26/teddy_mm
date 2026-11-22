"""A planted system with two modalities, two steps inside the model and two types of cell, where the place of every
loss is known. It is the check that the general cycle (cycle.py, outer.py) needs nothing of the system it was built
on: no name here occurs anywhere else in the package.

  modality a   carries factors 1, 2, 4 and 5            (the model's own input)
  modality b   carries factor 3                         (measured, and not given to the model)
  first step   keeps factors 1, 2 and 4, drops factor 5 (mid: the output of the encoder of modality a)
  second step  keeps factors 1 and 4, drops factor 2    (z: the stage the head reads)
  head         follows an unrelated direction of z, and a tenth of factor 1   (the system's own output, as trained)
  readout      type one: factors 1 + 2 + 3 + 5;  type two: factors 4 + 2 + 3 + 5;  plus noise; measured outside the system

  type  event      acts on             so the loss is                              and the loop should say
  one   p1 p2 p3   factor 1            in the head                                 closed by a repair of the head
  one   early      factor 5            in the first step                           upstream, between the input and mid
  one   enc        factor 2            in the second step                          upstream, between mid and z
  one   modb       factor 3            not in the model's input; b carries it      unreadable; another measured input reads it
  one   post       the readout itself  in no input at all                          unreadable; no input reads it
  one   weak       factor 1, barely    nowhere that can be followed                refused at selection
  two   q1 q2 q3   factor 4            in the head, along another direction        the repair of type one does not carry over
  two   odd        factor 2            in the second step, pooled with q1 q2 q3    one event inside the pool differs; upstream

`sound=True` plants a system with nothing wrong (a head that reads the right factor, and only events it can read): the
loop should leave no open item and claim no repair.

  python -m anm_diag.toy [folder]      (runs the cycle and prints what it says)
"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

from . import cycle, outer
from . import scope as SC

SD = np.array([1.0, 0.4, 0.4, 1.0, 0.4])                     # the spread of the five factors among control cells; the hidden ones are kept small so that a stage that carries a factor reads it within the margin
TYPES = {"one": {"reads": 0, "pool": {"p1": (0, -3.0), "p2": (0, -2.5), "p3": (0, -3.5)}, "alone": {"early": (4, -3.0), "enc": (1, -3.0), "modb": (2, -3.0), "post": (None, -3.0)},
                 "weak": {"weak": (0, -0.03)}},
         "two": {"reads": 3, "pool": {"q1": (3, -3.0), "q2": (3, -2.5), "q3": (3, -3.5), "odd": (1, -3.0)}, "alone": {}, "weak": {}}}
TRUTH = {("one", "early"): "upstream", ("one", "enc"): "upstream", ("one", "modb"): "unreadable", ("one", "post"): "unreadable", ("two", "odd"): "upstream"}
LOST = {("one", "early"): ["the input", "mid"], ("one", "enc"): ["mid", "z"], ("two", "odd"): ["mid", "z"]}        # between which two stages
DECISION = {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 1000, "seed": 0, "min_clusters": 10}
CLASSES = {"ridge": {"regressor": "RidgeCV", "alphas": [float(a) for a in np.logspace(-1, 5, 13)]},
           "hgb": {"regressor": "HistGradientBoostingRegressor", "max_iter": 40, "learning_rate": 0.1, "max_leaf_nodes": 8}}
STAGES = ("view_a", "view_b", "mid", "z", "perm")


def events_of(name, sound=False):
    T = TYPES[name]
    return {k: v for k, v in T["pool"].items() if v[0] == T["reads"]} if sound else {**T["pool"], **T["alone"], **T["weak"]}


def planted(seed=0, n_control=4000, n_event=800, sound=False):
    """The cells of both types, by type. One system reads them all: the loadings and the encoder are the same."""
    rng = np.random.default_rng(seed)
    typ, unit, reads, shift = [], [], [], []
    for name in TYPES:
        for u, (factor, d) in [("control", (None, 0.0)), *events_of(name, sound).items()]:
            k = n_control if u == "control" else n_event
            typ += [name] * k
            unit += [u] * k
            reads += [TYPES[name]["reads"]] * k
            shift += [(factor, d)] * k
    typ, unit, reads, n = np.array(typ), np.array(unit), np.array(reads), len(unit)
    f, own = rng.normal(size=(n, 5)) * SD, np.zeros(n)
    for i, (factor, d) in enumerate(shift):
        if factor is None:
            own[i] = d
        else:
            f[i, factor] += d
    def L(dims):                                             # loadings of a fixed size, so that a view carries its factor at every seed
        v = rng.normal(size=(1, dims))
        return v * np.sqrt(dims) / np.linalg.norm(v)

    noise = lambda dims, s=0.2: s * rng.normal(size=(n, dims))
    a = np.hstack([f[:, [0]] * L(8) + noise(8), f[:, [1]] * L(8) + noise(8), f[:, [3]] * L(8) + noise(8), noise(4, 1.0), f[:, [4]] * L(8) + noise(8)])
    b = np.hstack([f[:, [2]] * L(6) + noise(6), noise(4, 1.0)])
    rot = lambda x: x @ np.linalg.qr(rng.normal(size=(x.shape[1], x.shape[1])))[0]
    mid = rot(a[:, :28])                                     # the first step drops factor 5
    z = rot(np.hstack([a[:, :8], a[:, 16:28]]))              # the second step drops factor 2
    y = f[np.arange(n), reads] + f[:, 1] + f[:, 2] + f[:, 4] + own + 0.3 * rng.normal(size=n)
    o = f[np.arange(n), reads] + 0.1 * rng.normal(size=n) if sound else a[:, 24] + 0.1 * f[:, 0]
    cells = {"ids": np.array([f"cell{i:05d}" for i in range(n)]), "unit": unit, "treated": unit != "control", "view_a": a, "view_b": b, "mid": mid, "z": z, "y": y, "o": o}
    out = {}
    for name in TYPES:
        m = np.where(typ == name)[0]
        out[name] = {k: v[m] for k, v in cells.items()}
        out[name]["perm"] = out[name]["view_a"][rng.permutation(len(m))]
    return out


def spec(name, arrays, events, pooled):
    return dict(name=name, kind="perturbation", study="planted", readout_name="readout", events=list(events), also_pooled=pooled, units_exclude=["control"], reads="z", control_folds=5,
                decision=DECISION, label_basis="design", keep_draws=True, readout_classes=json.loads(json.dumps(CLASSES)),
                keys=dict(readouts=["y"], stages={k: k for k in STAGES}, unit="unit", control={"key": "A", "equals": 0}, cluster="unit", frozen={"array": "o"}), arrays={"audit": arrays})


def run(folder=None, seed=0, sound=False, **size):
    folder = Path(folder or tempfile.mkdtemp())
    folder.mkdir(parents=True, exist_ok=True)
    cells, sel, sites, alone = {}, {}, {}, {}
    for name, d in planted(seed, sound=sound, **size).items():
        sel[name] = cycle.select(d["ids"], d["unit"], d["y"], seed=1)       # admission on a selection share; those cells are read no further
        k = ~sel[name]["share"]
        d = cells[name] = {q: v[k] for q, v in d.items()}
        np.savez(folder / f"{name}.npz", A=d["treated"].astype(int), **{q: d[q] for q in ("ids", "unit", "y", "o", *STAGES)})
        ev = events_of(name, sound)
        sites[name] = [s for s in ev if s in TYPES[name]["pool"] and s in sel[name]["admitted"]]
        alone[name] = [s for s in ev if s not in TYPES[name]["pool"] and s in sel[name]["admitted"]]
    whole = {name: cycle.LinearRows(d["z"], d["y"], d["o"], np.where(~d["treated"])[0], [(s, np.where(d["unit"] == s)[0]) for s in sites[name]]) for name, d in cells.items() if sites[name]}
    res, targets = {"admitted": {k: v["admitted"] for k, v in sel.items()}, "types": {}}, []
    for name, d in cells.items():
        if not sites[name]:                                  # no event of the pool passes the screen: not a target (protocol 3c)
            t = {"readout": "readout", "group": name, "type": {"type": name}, "events_n": 0, "admitted": False, "response": {"snr": max(x["snr"] for x in sel[name]["table"])}}
            res["types"][name] = {"read_as": "not read", "target": t, "status": outer.target_status(t)}
            targets.append(t)
            continue
        oth = next(k for k in cells if k != name)
        do = cells[oth]
        grids = [cycle.grid(SC.run_spec(spec(name, f"{name}.npz", sites[name], True), folder))]
        if alone[name]:
            grids.append(cycle.grid(SC.run_spec(spec(name + "_alone", f"{name}.npz", alone[name], False), folder)))
        across = dict(other=(do["z"], do["y"], do["treated"].astype(int), do["unit"], sites[oth]), fixed={"other": whole[oth]}) if sites[oth] else {}
        T = cycle.one_type("readout", name, {"type": name}, d["z"], d["y"], d["o"], d["treated"], d["unit"], sites[name], *cycle.merge(*grids), "z", inputs=["view_a"], between=["mid"],
                           others=["view_b"], control="perm", alone=alone[name], pair="one and two", evidence="planted", **across)
        targets.append(T["target"])
        res["types"][name] = {**T, "status": outer.target_status(T["target"])}
    res["items"] = outer.open_items(targets, {"final_units": 1})
    return res


def main():
    r = run(sys.argv[1] if len(sys.argv) > 1 else None)
    f = lambda x: f"{x['u']:.2f} [{x['ci'][0]:.2f}, {x['ci'][1]:.2f}] {x['verdict']}"
    print("admitted on the selection share:", r["admitted"])
    for name, T in r["types"].items():
        print(f"\ntype {name}, read as {T['read_as']}: status {T['status']}")
        if T["read_as"] == "a pool":
            e = T["reading"]
            print(f"  pool of {e['events']} | reference {f(e['reference'])} | as trained {f(e['as_trained'])} | repaired {f(e['arms']['rows'])} | gate 3 {e['gate3']['passes']}"
                  + (f" | repair of the other type {f(e['arms']['other'])}" if "other" in e["arms"] else "") + f" | repair: {T['target']['repair']}")
        for u, p in (T.get("placed") or {"events": {}})["events"].items():
            print(f"  {u:5s} input {f(p['input'])} | z {f(p['model_stage'])} | other input " + ", ".join(f"{k} {f(v)}" for k, v in p["other_inputs"].items()) + f" | permuted {p['control']['u']:.2f} | "
                  + ", ".join(f"{x['from']} to {x['to']} {x['added']:.2f} {x['reading']}" for x in p["along_the_model"]))
        if T.get("placed"):
            print("  positive control:", T["placed"]["positive_control"], "| control view:", T["placed"]["control_view"]["holds"])
    print()
    for i in r["items"]:
        print(f"  {i['id']} {i['group']:4s} {i['event']:5s} {i['reading']:20s} {i['lever']} | {i['question']} | {i.get('forecast') or i['next_step_needs']}"
              + (f" | {i['note']}" if i["reading"] == "upstream" and i.get("note") else ""))


if __name__ == "__main__":
    main()
