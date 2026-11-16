#!/usr/bin/env python3
"""D6 of gate 4: what is left. The grid of `donors_hao` (runs/scope_20261116) on the changed head, scored donors only.

  python grid_changed.py --heads HEADS_DIR --work WORK --bm-arrays cite_arrays.npz --pbmc-arrays ext_cite_arrays.npz \\
      --tmp TMP [--seeds 0 1 2] [--fit-sets M1 M2] [--summary-only]

WORK: the arrays runs/scope_20261116/scripts/build.py wrote (bm_train.npz, hao_donors.npz): the RNA view and the
state of every cell, in the order of the packs. Nothing is fitted here but the readouts of `anm_diag scope`.

For every seed and fit set, two runs with the SPEC of `donors_hao`: the `coverage` head, and `as_is` of the same seed,
both on the four donors the `coverage` head was not fitted on. The training cells of a run are the cells its head was
fitted on. Between the two runs of a pair only what the head computes differs: its two stages, its output, and its
training cells. The fit sets read are M1 and M2 (batch seen), the configuration of the fits gate 5 makes.

A run here has four donors, so the other units of a donor are three donors, not seven: numbers compare inside a pair
and with each other, not with runs/scope_20261116. With three control donors the `shift` target has almost no
resolution; it is written for completeness and the summary is about `variation`.

Writes ../results/grid_changed/<run>/{SPEC.yaml, scope.json, scope_events.jsonl} and ../results/what_is_left.json.
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path

import numpy as np

import common
import plan
from anm_diag import scope

HEADS = ("as_is", "coverage")
STAGES = ("RNA512", "Z", "HEAD_LAT", "HEAD_OUT")
DECISION = {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 1000, "seed": 0, "min_clusters": 10}     # as runs/scope_20261116
STEPS = (("frozen", "train"), ("train", "units"), ("units", "own"))


def head_stages(W, z, chunk=20000):
    """The head from its weights, as runs/scope_20261116 ran it: state -> HEAD_LAT -> log mean (HEAD_OUT) -> mean."""
    W = {k: np.asarray(v, np.float64) for k, v in W.items() if k != "adt_names"}
    lat, out, mu = [], [], []
    for i in range(0, len(z), chunk):
        h = np.asarray(z[i:i + chunk], np.float64)
        h = h / (np.linalg.norm(h, axis=1, keepdims=True) + 1e-6)
        got = []
        for part in ("mlp__net", "dec__mu.net"):
            h = common.silu(h @ W[f"{part}.0.weight"].T + W[f"{part}.0.bias"])
            h = common.silu(h @ W[f"{part}.3.weight"].T + W[f"{part}.3.bias"])
            h = h @ W[f"{part}.6.weight"].T + W[f"{part}.6.bias"]
            got.append(h)
        spl = np.logaddexp(0.0, got[1])
        lat.append(got[0].astype(np.float32)); out.append(np.log(np.maximum(spl, 1e-6)).astype(np.float32)); mu.append((common.SF * spl).astype(np.float32))
    return np.vstack(lat), np.vstack(out), np.vstack(mu)


def run_arrays(W, cells, scored, wb, wh, bm_split, shared):
    """The two arrays of one run. wb / wh: the work arrays of the bone-marrow training cells and of the second study.
    cells: the head's record of the cells it was fitted on (rows of the packs). Returns (audit, train)."""
    donors = np.asarray(wh["donor"]).astype(str)
    lat, out, mu = head_stages(W, wh["Z"])
    m = np.isin(donors, list(scored))
    audit = dict(RNA512=wh["RNA512"][m], Z=wh["Z"][m], HEAD_LAT=lat[m], HEAD_OUT=out[m], Y=wh["Y"][m], unit=donors[m], donor=donors[m],
                 frozen=mu[m][:, shared])
    tr = np.where(np.asarray(bm_split).astype(str) == "train")[0]
    pos = np.searchsorted(tr, cells["bm_index"])
    if len(wb["Z"]) != len(tr) or pos.max() >= len(tr) or not np.array_equal(tr[pos], cells["bm_index"]):
        raise SystemExit("the head's bone-marrow cells are not training cells of this pack")
    pi = np.asarray(cells["pbmc_index"], np.int64)
    if set(donors[pi]) & set(scored):
        raise SystemExit("a scored donor is among the cells the head was fitted on")
    latb, outb, _ = head_stages(W, wb["Z"][pos])
    train = dict(RNA512=np.vstack([wb["RNA512"][pos], wh["RNA512"][pi]]), Z=np.vstack([wb["Z"][pos], wh["Z"][pi]]),
                 HEAD_LAT=np.vstack([latb, lat[pi]]), HEAD_OUT=np.vstack([outb, out[pi]]), Y=np.vstack([wb["Y"][pos], wh["Y"][pi]]),
                 donor=np.concatenate([np.asarray(wb["donor"]).astype(str)[pos], donors[pi]]))
    return audit, train


def spec(name):
    return dict(name=name, kind="donor", study="hao2021", arrays={"audit": f"{name}__audit.npz", "train": f"{name}__train.npz"},
                keys=dict(readouts="Y", stages={s: s for s in STAGES}, unit="unit", block="donor", train_block="donor", frozen={"array": "frozen"}),
                events="each", reads="Z", input_view="RNA512", variation=True, max_fit=30000, decision=DECISION, label_basis="design",
                confounded_with=[], readout_name="112 proteins", readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": common.ALPHAS}})


def one_run(r, heads, wb, wh, bm_split, shared, tmp, res, decision=None):
    """Arrays and SPEC of one run into tmp, `anm_diag scope` on them, the arrays removed again."""
    import yaml
    d = Path(heads) / r["weights"]
    audit, train = run_arrays(dict(np.load(d / "head.npz")), dict(np.load(d / "fit_cells.npz")), r["scored"], wb, wh, bm_split, shared)
    np.savez(tmp / f"{r['name']}__audit.npz", **audit)
    np.savez(tmp / f"{r['name']}__train.npz", **train)
    del audit, train
    sp = spec(r["name"])
    if decision:
        sp["decision"] = {**sp["decision"], **decision}
    (tmp / f"{r['name']}.yaml").write_text(yaml.safe_dump(sp, sort_keys=False))
    scope.main([str(tmp / f"{r['name']}.yaml"), "--out", str(res / r["name"])])                  # the same entry as `python -m anm_diag scope`
    (res / r["name"] / "SPEC.yaml").write_text((tmp / f"{r['name']}.yaml").read_text())
    for f in ("__audit.npz", "__train.npz"):
        (tmp / f"{r['name']}{f}").unlink()


def runs(seeds, fit_sets):
    return [{"name": f"s{s}_{fs}_{h}", "seed": s, "fit_set": fs, "head": h, "weights": f"s{s}_as_is" if h == "as_is" else f"s{s}_coverage_{fs}",
             "scored": plan.scored_donors(fs)} for s in seeds for fs in fit_sets for h in HEADS]


def cell(T, scope_, stage):
    return abs(T["stages"]["system" if scope_ == "frozen" else stage][scope_]["u"])


def summarise(results, seeds, fit_sets):
    """results: {run name: the dict `anm_diag scope` wrote}. Every donor is scored by one fit set, so the arrays are
    (donors, seeds); intervals resample the donors, and the seeds inside each head."""
    todo = runs(seeds, fit_sets)
    donors = sorted({d for r in todo for d in r["scored"]})
    out = {"donors": donors, "seeds": list(seeds), "fit_sets": list(fit_sets), "runs": [r["name"] for r in todo],
           "note": "four donors per run: the other units of a donor are three donors; compare inside a pair, not with runs/scope_20261116"}
    for tg in ("variation", "shift"):
        cells = [("frozen", "system")] + [(sc, st) for sc in ("train", "units", "own") for st in STAGES]
        U = {h: {c: np.full((len(donors), len(seeds)), np.nan) for c in cells} for h in HEADS}
        readings = {h: Counter() for h in HEADS}
        for r in todo:
            for e in results[r["name"]]["results"]:
                T = e["targets"].get(tg)
                if T is None:
                    continue
                i, j = donors.index(e["event"]["unit"]), list(seeds).index(r["seed"])
                for c in cells:
                    U[r["head"]][c][i, j] = cell(T, *c)
                readings[r["head"]][T["reading"]["missing"]] += 1
        key = lambda c: c[0] if c[0] == "frozen" else f"{c[0]}@{c[1]}"
        fam = lambda a, b: {x: scope.arms_family({"a": a, "b": b}, base="a")["arms"]["b"][x] for x in ("gain", "ci")}      # mean of a minus mean of b
        at = lambda h, sc: U[h][("frozen", "system") if sc == "frozen" else (sc, "Z")]
        out[tg] = {
            "u": {h: {key(c): float(np.nanmean(v)) for c, v in U[h].items()} for h in HEADS},
            "steps_at_the_stage_the_head_reads": {h: {f"{a}->{b}": fam(at(h, a), at(h, b)) for a, b in STEPS} for h in HEADS},
            "added_inside_the_head_units_scope": {h: fam(U[h][("units", "HEAD_OUT")], U[h][("units", "Z")]) for h in HEADS},
            "fall_from_as_is_to_coverage": {key(c): fam(U["as_is"][c], U["coverage"][c]) for c in cells},
            "readings": {h: dict(readings[h]) for h in HEADS},
            "approximate": True,
        }
    return out


def main():
    p = argparse.ArgumentParser()
    for k in ("heads", "work", "bm-arrays", "pbmc-arrays", "tmp"):
        p.add_argument(f"--{k}", type=Path, required=True)
    p.add_argument("--seeds", type=int, nargs="+", default=list(plan.SEEDS))
    p.add_argument("--fit-sets", nargs="+", default=["M1", "M2"])
    p.add_argument("--summary-only", action="store_true")
    a = p.parse_args()
    res = common.RUN / "results" / "grid_changed"
    todo = runs(a.seeds, a.fit_sets)
    left = [r for r in todo if not (res / r["name"] / "scope.json").exists()]
    if left and not a.summary_only:
        missing = [r["weights"] for r in left if not (a.heads / r["weights"] / "DONE").exists()]
        if missing:
            raise SystemExit(f"heads not finished: {', '.join(sorted(set(missing)))}")
        a.tmp.mkdir(parents=True, exist_ok=True)
        wb, wh = np.load(a.work / "bm_train.npz", allow_pickle=False), np.load(a.work / "hao_donors.npz", allow_pickle=False)
        wb, wh = {k: wb[k] for k in ("RNA512", "Z", "Y", "donor")}, {k: wh[k] for k in ("RNA512", "Z", "Y", "donor")}
        bm, ex = np.load(a.bm_arrays, allow_pickle=False), np.load(a.pbmc_arrays, allow_pickle=False)
        shared = np.where(ex["adt_present"].astype(bool))[0]
        tr = bm["split"].astype(str) == "train"
        if not (np.array_equal(wh["donor"].astype(str), ex["donors"].astype(str)) and np.array_equal(wb["donor"].astype(str), bm["donors"].astype(str)[tr])
                and np.array_equal(wh["Y"], ex["adt"][:, shared].astype(np.float32)) and np.array_equal(wb["Y"], bm["adt"][tr][:, shared].astype(np.float32))):
            raise SystemExit("the work arrays are not those of these packs")
        for r in left:
            t0 = time.time()
            print(f"== {r['name']}", flush=True)
            one_run(r, a.heads, wb, wh, bm["split"], shared, a.tmp, res)
            print(f"   done in {time.time() - t0:.0f} s", flush=True)
    if any(not (res / r["name"] / "scope.json").exists() for r in todo):
        raise SystemExit("runs are missing; nothing is summarised before all are there")
    out = summarise({r["name"]: json.loads((res / r["name"] / "scope.json").read_text()) for r in todo}, a.seeds, a.fit_sets)
    (common.RUN / "results" / "what_is_left.json").write_text(json.dumps(out, indent=1))
    for h in HEADS:
        v = out["variation"]
        print(h, {k: round(x, 4) for k, x in v["u"][h].items() if k in ("frozen", "train@Z", "units@Z", "own@Z", "units@HEAD_OUT")}, v["readings"][h])


if __name__ == "__main__":
    main()
