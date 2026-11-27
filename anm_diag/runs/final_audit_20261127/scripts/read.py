#!/usr/bin/env python3
"""The reading of the final audit, after forecasts.json was committed and the final cells were exported. Each part once.

  python read.py --part pairs    --cells final_cells.tsv --states MEL... --states-thp THP... --packs PACK_A PACK_B PACK_4 FINAL_PACK
  python read.py --part accounts (the same) --counts MEL_COUNTS... --counts-thp THP_COUNTS...

pairs -> results/reading.json
  Every forecast pair: the measured change of the unit against the reference of its stratum (the final control cells;
  in THP-1 its control cells), in spreads of the readout among the control cells the law was fitted on, and its error;
  inside or outside the intervals of S2 and of S3; the errors of S0, S2 and S3. The edges listed in forecasts.json,
  asked of the rows of the audit in each stratum; those rows were forecast without them (`events.test_edges`). Also
  every unit on every readout as a record, and every melanoma pair measured against the isotype antibodies.
accounts -> results/accounts.json
  The account (`account.account`) of every forecast pair whose change is at least four errors: at the input view V and
  the state Z, and in the co-culture also at the typed state T and its control R. Here an event is read against every
  control cell of its stratum. Then `outer.open_parts`. It describes; no hypothesis rests on it.

--rehearsal --results DIR: on the made-up file of rehearsal.py, to find faults of the code before the one reading.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fa as F  # noqa: E402
from _fa import V  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402
from anm_diag import outer as OU  # noqa: E402

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
RUNS = F.RUNS
TYPED = {"mel_co": "mel_co"}
PAIRS = (("Z", "T"), ("R", "T"), ("Z", "V"), ("T", "V"))
OUTPUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}   # melanoma readout -> output of the head as trained
SEED_SVD = 20261125


def inputs(a):
    fc = json.loads((F.RUN / "results" / "forecasts.json").read_text())
    assert fc["cells_file_sha256"] == sha(a.cells), "the forecasts were made for another list of cells"
    assert a.rehearsal or fc["packs_sha256"] == [sha(p) for p in a.packs[:-1]], "the forecasts stand on other antibody files"
    cells = F.table(a.cells)
    st, packs = {"frangieh2021": F.States(a.states), "papalexi2021": F.States(a.states_thp)}, F.Packs(a.packs)
    S = {}
    for g in F.STRATA:
        study = F.STRATA[g][0]
        c = F.stratum(g, cells, st, F.final_units(g))
        assert packs.has(study, c["ids"]).all(), "a cell of the list has no antibody columns"
        c["Y"], log, bg = packs.readouts(study, c["ids"])
        c["against_the_background"] = {r: log[r] - bg for r in log} if bg is not None else None
        assert int(c["ref"].sum()) == fc["strata"][g]["reference_cells"] and int((c["kind"] == "control").sum()) == fc["strata"][g]["control_cells"]
        S[g] = c
    return fc, S


def measure(y, m, ref, sd):
    """The change of the cells m against the reference cells, in the spread sd, and its error."""
    v, w = y[m], y[ref]
    return float((v.mean() - w.mean()) / sd), float(np.sqrt(v.var(ddof=1) / len(v) + w.var() / len(w)) / sd)


def part_pairs(a, out_file):
    fc, S = inputs(a)
    res = {"cells_file_sha256": fc["cells_file_sha256"], "pack_sha256": sha(a.packs[-1]), "rehearsal_on_made_up_values": bool(a.rehearsal), "min_snr": AC.MIN_SNR, "spread": {}, "records": {}, "rows": []}
    for g, c in S.items():
        C = c["kind"] == "control"
        res["spread"][g] = {r: float((y[C] - y[C].mean()).std()) for r, y in c["Y"].items()}
        res["records"][g] = {}
        for e in sorted(set(c["unit"][c["kind"] == "new"].tolist())):
            m = c["unit"] == e
            res["records"][g][e] = {"cells": int(m.sum()), "changes": {r: list(measure(y, m, c["ref"], res["spread"][g][r])) for r, y in c["Y"].items()}}
    for f in fc["rows"]:
        g, e, r = f["stratum"], f["event"], f["readout"]
        c = S[g]
        m = c["unit"] == e
        assert int(m.sum()) == f["cells"] and set(c["kind"][m].tolist()) == {"new"} and int(c["ref"].sum()) == f["reference_cells"], (g, e)
        ch, se = res["records"][g][e]["changes"][r]
        row = {k: f[k] for k in ("stratum", "event", "readout", "cells", "reference_cells", "edge", "type", "scale", "size_of_the_edge_from", "shift_times_noise")}
        row.update(change=ch, se=se, snr=abs(ch) / se, S0={"forecast": f["S0"]["forecast"], "error": ch - f["S0"]["forecast"]})
        for s in ("S2", "S3"):
            p = f[s]
            row[s] = {**{k: p[k] for k in p if k != "parts_of_the_spread"}, "error": ch - p["forecast"], "z": (ch - p["forecast"]) / p["sd"], "inside": bool(p["interval"][0] <= ch <= p["interval"][1])}
        if c["against_the_background"] is not None:
            cb, sb = measure(c["against_the_background"][r], m, c["ref"], res["spread"][g][r])
            row["against_the_background"] = {"change": cb, "se": sb}
        res["rows"].append(row)
    assert len(res["rows"]) == fc["pairs"]
    tested, without = [], []
    for g in F.MEL:
        rows = [{"event": r["event"], "readout": r["readout"], "change": r["change"], "forecast": r["S3"]["forecast"], "sd": r["S3"]["sd"]} for r in res["rows"] if r["stratum"] == g and r["edge"] is None]
        for e in EV.test_edges(fc["edges"]["to_test_on_the_rows_of_the_audit"][g], rows):
            if "z" not in e.get("test", {}):
                without.append(dict(e, edge_id=f"{e['edge_id']}@{g}", stratum=g))
                continue
            e = dict(e, edge_id=f"{e['edge_id']}@{g}", stratum=g)
            e["test"]["scale_of_the_readout_here"] = fc["scale"][e["readout"]][g]
            e["test"]["direct_part"] = e["test"]["change"] - e["test"]["forecast_without_the_edge"]
            row = next(r for r in res["rows"] if r["stratum"] == g and r["event"] == e["site"] and r["readout"] == e["readout"])
            e["test"]["inside_the_interval_of_S3"] = row["S3"]["inside"]
            tested.append(e)
    res["edges"] = {"tested": tested, "with_no_row": without}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    rows = res["rows"]
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
    print(f"pairs {len(rows)} | inside S2 {sum(r['S2']['inside'] for r in rows)} S3 {sum(r['S3']['inside'] for r in rows)} | median z2 S2 {np.median([r['S2']['z'] ** 2 for r in rows]):.2f} "
          f"S3 {np.median([r['S3']['z'] ** 2 for r in rows]):.2f} | error S0 {rms([r['S0']['error'] for r in rows]):.3f} S2 {rms([r['S2']['error'] for r in rows]):.3f} S3 {rms([r['S3']['error'] for r in rows]):.3f} "
          f"| pairs with a change {sum(r['snr'] >= AC.MIN_SNR for r in rows)}")
    print("edges:", [(e["edge_id"], e["status"], round(e["test"]["z"], 1)) for e in tested], "| no row:", [e["edge_id"] for e in without])


def halves(g, c, unit):
    """The halves, as in batches 2 to 4: the development cells keep those of run 40 (the first cells of the stratum in the
    list, in their order; checked against the manifest of the typed states), the other cells are split by cycle.halves."""
    sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
    from _common import SEED
    from anm_diag import cycle
    first, n0 = np.zeros(len(c["ids"]), bool), 0
    if g in TYPED:
        man = json.loads((RUNS / "typed_repair_20261125" / "results" / "states" / f"{TYPED[g]}_manifest.json").read_text())
        n0 = man["cells"]
        first[:n0] = cycle.halves(unit[:n0], SEED)
        assert hashlib.sha256(first[:n0].tobytes()).hexdigest()[:16] == man["sha256_of_the_halves"], "the halves of the development cells are not those of run 40"
    first[n0:] = cycle.halves(unit[n0:], SEED)
    return first, SEED


def input_view(files, ids, control):
    """V: 512 components of the log of counts per 10,000 of the stratum's cells (genes detected in 1% of control cells)."""
    import scipy.sparse as sp
    from sklearn.decomposition import TruncatedSVD
    sys.path.insert(0, str(RUNS / "outer_step_20261122" / "scripts"))
    from build2 import lognorm
    mats, pos, genes, n = [], {}, None, 0
    for f in files:
        d = np.load(f)
        assert genes is None or (genes == d["genes"]).all(), "the count files list different genes"
        genes = d["genes"]
        M = sp.csr_matrix((d["vals"].astype(np.float32), (d["rows"].astype(np.int64), d["cols"].astype(np.int64))), shape=tuple(int(v) for v in d["shape"]))
        for k, i in enumerate(d["ids"].astype(str)):
            pos.setdefault(i, n + k)
        mats.append(M)
        n += M.shape[0]
    X = sp.vstack(mats).tocsr()[np.array([pos[i] for i in ids])].tocsr()
    R = lognorm(X, np.where(control)[0])
    return TruncatedSVD(n_components=512, algorithm="randomized", n_iter=4, random_state=SEED_SVD).fit_transform(R).astype(np.float64), int(R.shape[1])


def part_accounts(a, out_file):
    fc, S = inputs(a)
    R = json.loads((out_file.parent / "reading.json").read_text())
    assert R["pack_sha256"] == sha(a.packs[-1]), "the pairs were read from another file"
    want = [r for r in R["rows"] if r["snr"] >= AC.MIN_SNR]
    partial = out_file.with_suffix(".partial")
    res = {"pack_sha256": R["pack_sha256"], "min_snr": AC.MIN_SNR, "pairs_with_a_change": len(want), "strata": {}, "accounts": {}}
    if partial.exists():
        old = json.loads(partial.read_text())
        if old["pack_sha256"] == res["pack_sha256"]:
            res["accounts"], res["strata"] = old["accounts"], old["strata"]
            print("taken up again after", len(res["accounts"]), "accounts", flush=True)
    print(len(want), "pairs with a change of at least", AC.MIN_SNR, "errors:", {g: sum(r["stratum"] == g for r in want) for g in F.STRATA}, flush=True)
    head = None
    for g in F.STRATA:
        rows = [r for r in want if r["stratum"] == g and f"{g}|{r['event']}|{r['readout']}" not in res["accounts"]]
        if not rows:
            continue
        c = S[g]
        unit = np.where(c["kind"] == "reference", "control", c["unit"])     # here an event is read against every control cell of its stratum
        C = unit == "control"
        first, seed = halves(g, c, unit)
        Vw, genes = input_view(a.counts if g in F.MEL else a.counts_thp, c["ids"], C)
        stages = {"V": Vw, "Z": c["Z"]}
        if g in TYPED:
            blocks, sets = c["states"].rows(c["ids"], "blocks"), c["states"].sets
            nz = np.linalg.norm(c["raw"], axis=1, keepdims=True) + 1e-6
            blk = lambda name, f: blocks[:, sets.index(f"{name}|{f}")]
            two = lambda name: tuple(np.hstack([blk(name, f), c["raw"] - blk(name, f)]) / nz for f in (0, 1))
            stages.update(T=two("A"), R=two("R"))
            del blocks
        res["strata"][g] = {"stages": {k: int((v[0] if isinstance(v, tuple) else v).shape[1]) for k, v in stages.items()}, "genes_of_the_input_view": genes, "cells": {"control": int(C.sum()), "all": int(len(C))}}
        out = {}
        if g in F.MEL and any(r["readout"] in OUTPUT for r in rows):
            if head is None:
                sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
                from _common import forward_np, unit as unit_, weights
                P0, names, _j = weights()
                head = lambda raw: forward_np(P0, unit_(raw))[1]
            o_all = head(c["raw"])
            out = {r: o_all[:, names.index(n)].astype(np.float64) for r, n in OUTPUT.items()}
        for r in sorted(rows, key=lambda r: (r["readout"], r["event"])):
            e, ro = r["event"], r["readout"]
            acc = AC.account(stages, c["Y"][ro], unit, first, e, others=[], inputs=("V",), state="Z", o=out.get(ro), pairs=PAIRS if g in TYPED else (("Z", "V"),), seed=seed)
            acc.update(stratum=g, event=e, readout=ro, edge=r["edge"], type=r["type"], change_in_spreads=r["change"], snr=r["snr"])
            res["accounts"][f"{g}|{e}|{ro}"] = acc
            partial.write_text(json.dumps(res, default=float) + "\n")
            print(f"{g:7s} {e:9s} -> {ro:6s} change {r['change']:+.2f} ({r['snr']:.1f} errors) | {acc['chain']['says']}", flush=True)
    graph = fc["edges"]["declared"] + [e for e in fc["edges"]["on_record"] if e["status"] == "admitted"] + [e for e in R["edges"]["tested"] if e["status"] == "admitted"]
    res["open_parts"] = OU.open_parts([{"readout": v["readout"], "group": v["stratum"], "type": v["type"], "event": v["event"], "account": v} for v in res["accounts"].values()], edges=graph)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    if partial.exists():
        partial.unlink()
    print("accounts", len(res["accounts"]), "| open parts", OU.counts(res["open_parts"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["pairs", "accounts"], required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--states-thp", type=Path, nargs="+", required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--counts", type=Path, nargs="+")
    ap.add_argument("--counts-thp", type=Path, nargs="+")
    ap.add_argument("--results", type=Path, default=F.RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (F.RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out_file = a.results / {"pairs": "reading.json", "accounts": "accounts.json"}[a.part]
    if out_file.exists():
        raise SystemExit(f"{out_file.name} exists: this part of the reading was made once and is not made again")
    {"pairs": part_pairs, "accounts": part_accounts}[a.part](a, out_file)


if __name__ == "__main__":
    main()
