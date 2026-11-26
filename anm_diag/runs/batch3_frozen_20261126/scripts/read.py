#!/usr/bin/env python3
"""The reading of batch 3, after forecasts.json was committed and stage B was exported. Each part is made once.

  python read.py --part pairs    --cells batch3_cells.tsv --states B2_STATES_frangieh.npz B3_STATES_frangieh.npz --states-thp B2_STATES_papalexi.npz
                                 --read read_cells.npz --pack2 batch2_pack.npz --pack-a batch3_pack_A.npz --pack-b batch3_pack_B.npz
  python read.py --part accounts (the same) --counts B2_COUNTS_frangieh.npz B3_COUNTS_frangieh.npz --counts-thp COUNTS_papalexi.npz

pairs -> results/reading.json
  1. Every forecast pair: the measured change in spreads of the readout among the control cells of its stratum and its
     error (as `Law.measure`); inside or outside the intervals of S1 and of S2; the errors of S0, S1 and S2.
  2. The edges listed in forecasts.json before stage B was opened, asked of the rows of stage B of each melanoma stratum.
     Those rows were forecast without them (`events.test_edges`, rule G2).
  Also what the batch leaves as records: the measured change of every unit on every readout, and the spread of every
  readout among control cells. And, for the melanoma screen, the change of every forecast pair measured against the
  isotype antibodies instead of the mean of the 24 columns (score.py says why; it enters no hypothesis).
accounts -> results/accounts.json
  3. The account (`account.account`) of every forecast pair whose change is at least four errors: at the input view V
     (512 components of the RNA) and the state Z, and in the co-culture and THP-1 also at the typed state T and its
     control R. `others`: the teachers of the stratum. The head as trained has an output for the three readouts read
     before only, so `as_trained` is read for those. Then `outer.open_parts`.

--rehearsal --results DIR: a rehearsal on the made-up stage-B file of rehearsal.py, to find faults of the code before
the one reading. It writes outside this run's results and does not stop at the check of the values read before.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b3 as B  # noqa: E402
from _b3 import V  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402
from anm_diag import outer as OU  # noqa: E402

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
RUNS = B.RUN.parent
TYPED = {"mel_co": "mel_co", "thp1": "thp1"}                       # the strata with the sets of run 40, and their manifests
PAIRS = (("Z", "T"), ("R", "T"), ("Z", "V"), ("T", "V"))
OUTPUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}   # readout -> output of the head as trained (melanoma)
SEED_SVD = 20261125
MELANOMA = ("mel_co", "mel_ctl")


def inputs(a):
    fc = json.loads((B.RUN / "results" / "forecasts.json").read_text())
    assert fc["cells_file_sha256"] == sha(a.cells), "the forecasts were made for another list of cells"
    assert fc["pack_a_sha256"] == sha(a.pack_a), "the forecasts were made from another stage-A file"
    cells = B.table(a.cells)
    st = {"frangieh2021": B.States(a.states), "papalexi2021": B.States(a.states_thp)}
    packs = B.Packs([a.pack_a, a.pack_b])
    S = {}
    for g, (study, _cond) in B.STRATA.items():
        c = B.stratum(g, cells, st[study])
        assert packs.has(study, c["ids"]).all(), "a cell of the list is in neither stage"
        c["Y"], log, bg = packs.readouts(study, c["ids"])
        c["against_the_background"] = None if bg is None else {r: log[r] - bg for r in log}
        S[g] = c
    return fc, S


def measure(y, m, C, sd=None):
    """As `Law.measure`: the change of the cells m against the control cells C in spreads of y among C, and its error."""
    centre = y[C].mean()
    sd = float((y[C] - centre).std()) if sd is None else sd
    v = y[m]
    return float((v.mean() - centre) / sd), float(np.sqrt(v.var(ddof=1) / len(v) + y[C].var() / C.sum()) / sd)


def part_pairs(a, out_file):
    fc, S = inputs(a)
    used = B.used_values(a.read, a.pack2)
    res = {"cells_file_sha256": fc["cells_file_sha256"], "pack_a_sha256": fc["pack_a_sha256"], "pack_b_sha256": sha(a.pack_b), "min_snr": AC.MIN_SNR, "checks": {}, "spread": {}, "records": {}, "rows": []}
    # --- the readouts read before, formed here from the antibody columns, against the values read before
    for g, c in S.items():
        study = B.STRATA[g][0]
        before = c["kind"] != "new"
        for r in B.USED[study]:
            old = np.array([used[r][i] for i in c["ids"][before]])
            res["checks"][f"{g} {r}"] = float(np.abs(c["Y"][r][before] - old).max())
    assert a.rehearsal or max(res["checks"].values()) < 1e-4, res["checks"]
    res["rehearsal_on_made_up_values"] = bool(a.rehearsal)
    print("check: the readouts read before, formed again from the antibody columns, differ by at most", max(res["checks"].values()), flush=True)
    # --- what the batch leaves as records: every unit on every readout
    for g, c in S.items():
        C = c["kind"] == "control"
        res["spread"][g] = {r: float((y[C] - y[C].mean()).std()) for r, y in c["Y"].items()}
        res["records"][g] = {}
        for e in sorted(set(c["unit"][~C].tolist())):
            m = c["unit"] == e
            res["records"][g][e] = {"kind": str(c["kind"][m][0]), "cells": int(m.sum()), "changes": {r: list(measure(y, m, C)) for r, y in c["Y"].items()}}
    # --- 1. every forecast pair
    for f in fc["rows"]:
        g, e, r = f["stratum"], f["event"], f["readout"]
        c = S[g]
        C, m = c["kind"] == "control", c["unit"] == e
        assert int(m.sum()) == f["cells"] and set(c["kind"][m].tolist()) == {f["kind"]}, (g, e)
        ch, se = res["records"][g][e]["changes"][r]
        row = {k: f[k] for k in ("stratum", "event", "kind", "readout", "group", "cells", "edge", "type", "taught_by", "scale", "size_of_the_edge_from", "shift_times_noise")}
        row.update(change=ch, se=se, snr=abs(ch) / se, S0={"forecast": f["S0"]["forecast"], "error": ch - f["S0"]["forecast"]})
        for s in ("S1", "S2"):
            p = f[s]
            assert p["no_forecast"] is None
            row[s] = {"forecast": p["forecast"], "sd": p["sd"], "interval": p["interval"], "from_the_recorded_events": p["from_the_recorded_events"], "by_its_edge": p["by_its_edge"],
                      "error": ch - p["forecast"], "z": (ch - p["forecast"]) / p["sd"], "inside": bool(p["interval"][0] <= ch <= p["interval"][1])}
        if c["against_the_background"] is not None:                 # the same pair with the isotype antibodies as the reference, in the same spreads
            cb, sb = measure(c["against_the_background"][r], m, C, sd=res["spread"][g][r])
            row["against_the_background"] = {"change": cb, "se": sb}
        res["rows"].append(row)
    assert len(res["rows"]) == fc["pairs"]
    # --- 2. the edges proposed before stage B, asked of rows that were forecast without them
    to_test = fc["edges"]["to_test_on_the_rows_of_stage_B"]
    tested, have = [], set()
    for g in MELANOMA:
        rows = [{"event": r["event"], "readout": r["readout"], "change": r["change"], "forecast": r["S2"]["forecast"], "sd": r["S2"]["sd"]} for r in res["rows"] if r["stratum"] == g and r["edge"] is None]
        for e in EV.test_edges(to_test, rows):
            if "z" in e.get("test", {}):
                sc = fc["scale"][e["readout"]][g]
                e = dict(e, edge_id=f"{e['edge_id']}@{g}", stratum=g)
                e["test"]["scale_of_the_readout_here"] = sc
                e["test"]["direct_part"] = e["test"]["change"] - e["test"]["forecast_without_the_edge"]
                tested.append(e)
                have.add((e["site"], e["readout"]))
    res["edges"] = {"tested": tested, "with_no_row_in_stage_B": [e for e in to_test if (e["site"], e["readout"]) not in have]}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    rows = res["rows"]
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
    print(f"pairs {len(rows)} | inside S1 {sum(r['S1']['inside'] for r in rows)} S2 {sum(r['S2']['inside'] for r in rows)} | mean z2 S1 {np.mean([r['S1']['z'] ** 2 for r in rows]):.2f} "
          f"S2 {np.mean([r['S2']['z'] ** 2 for r in rows]):.2f} | error S0 {rms([r['S0']['error'] for r in rows]):.3f} S1 {rms([r['S1']['error'] for r in rows]):.3f} S2 {rms([r['S2']['error'] for r in rows]):.3f} "
          f"| pairs with a change {sum(r['snr'] >= AC.MIN_SNR for r in rows)}")
    print("edges:", [(e["edge_id"], e["status"], round(e["test"]["z"], 1)) for e in tested], "| no row:", [e["edge_id"] for e in res["edges"]["with_no_row_in_stage_B"]])


def halves(g, c):
    """The halves, as in batch 2: the development cells keep those of run 40 (the first cells of the stratum in the list,
    in their order; checked against the manifest of the typed states), the other cells are split by cycle.halves."""
    sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
    from _common import SEED
    from anm_diag import cycle
    first, n0 = np.zeros(len(c["ids"]), bool), 0
    if g in TYPED:
        man = json.loads((RUNS / "typed_repair_20261125" / "results" / "states" / f"{TYPED[g]}_manifest.json").read_text())
        n0 = man["cells"]
        first[:n0] = cycle.halves(c["unit"][:n0], SEED)
        assert hashlib.sha256(first[:n0].tobytes()).hexdigest()[:16] == man["sha256_of_the_halves"], "the halves of the development cells are not those of run 40"
    first[n0:] = cycle.halves(c["unit"][n0:], SEED)
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
    assert R["pack_b_sha256"] == sha(a.pack_b), "the pairs were read from another stage-B file"
    want = [r for r in R["rows"] if r["snr"] >= AC.MIN_SNR]
    partial = out_file.with_suffix(".partial")
    res = {"pack_b_sha256": R["pack_b_sha256"], "min_snr": AC.MIN_SNR, "pairs_with_a_change": len(want), "strata": {}, "accounts": {}}
    if partial.exists():
        old = json.loads(partial.read_text())
        if old["pack_b_sha256"] == res["pack_b_sha256"]:
            res["accounts"], res["strata"] = old["accounts"], old["strata"]
            print("taken up again after", len(res["accounts"]), "accounts", flush=True)
    print(len(want), "pairs with a change of at least", AC.MIN_SNR, "errors:", {g: sum(r["stratum"] == g for r in want) for g in B.STRATA}, flush=True)
    head = None
    for g, (study, _cond) in B.STRATA.items():
        rows = [r for r in want if r["stratum"] == g and f"{g}|{r['event']}|{r['readout']}" not in res["accounts"]]
        if not rows:
            continue
        c = S[g]
        C = c["kind"] == "control"
        first, seed = halves(g, c)
        Vw, genes = input_view(a.counts if study == "frangieh2021" else a.counts_thp, c["ids"], C)
        stages = {"V": Vw, "Z": c["Z"]}
        if g in TYPED:
            blocks, sets = c["states"].rows(c["ids"], "blocks"), c["states"].sets
            nz = np.linalg.norm(c["raw"], axis=1, keepdims=True) + 1e-6
            blk = lambda name, f: blocks[:, sets.index(f"{name}|{f}")]
            two = lambda name: tuple(np.hstack([blk(name, f), c["raw"] - blk(name, f)]) / nz for f in (0, 1))
            stages.update(T=two("A"), R=two("R"))
            del blocks
        res["strata"][g] = {"stages": {k: int((v[0] if isinstance(v, tuple) else v).shape[1]) for k, v in stages.items()}, "genes_of_the_input_view": genes,
                            "cells": {"control": int(C.sum()), "all": int(len(C))}, "others": B.teachers(g)}
        out = {}
        if study == "frangieh2021" and any(r["readout"] in OUTPUT for r in rows):
            if head is None:
                sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
                from _common import forward_np, unit, weights
                P0, names, _j = weights()
                head = lambda raw: forward_np(P0, unit(raw))[1]
            o_all = head(c["raw"])
            out = {r: o_all[:, names.index(n)].astype(np.float64) for r, n in OUTPUT.items()}
        for r in sorted(rows, key=lambda r: (r["readout"], r["event"])):       # by readout: the fit of the control cells is kept between events
            e, ro = r["event"], r["readout"]
            acc = AC.account(stages, c["Y"][ro], c["unit"], first, e, others=[h for h in B.teachers(g) if h != e], inputs=("V",), state="Z", o=out.get(ro), pairs=PAIRS, seed=seed)
            acc.update(stratum=g, event=e, readout=ro, kind=r["kind"], group=r["group"], edge=r["edge"], type=r["type"], change_in_spreads=r["change"], snr=r["snr"])
            res["accounts"][f"{g}|{e}|{ro}"] = acc
            partial.write_text(json.dumps(res, default=float) + "\n")
            print(f"{g:7s} {e:7s} -> {ro:6s} change {r['change']:+.2f} ({r['snr']:.1f} errors) | {acc['chain']['says']}", flush=True)
    r43 = json.loads((RUNS / "one_loop_20261125" / "results" / "retro.json").read_text())["edges"]
    graph = fc["edges"]["declared"] + [e for e in r43 if e["status"] == "admitted"] + [e for e in R["edges"]["tested"] if e["status"] == "admitted"]
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
    for k in ("read", "pack2", "pack-a", "pack-b"):
        ap.add_argument("--" + k, type=Path, required=True)
    ap.add_argument("--counts", type=Path, nargs="+")
    ap.add_argument("--counts-thp", type=Path, nargs="+")
    ap.add_argument("--results", type=Path, default=B.RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (B.RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out_file = a.results / {"pairs": "reading.json", "accounts": "accounts.json"}[a.part]
    if out_file.exists():
        raise SystemExit(f"{out_file.name} exists: this part of the reading was made once and is not made again")
    (part_pairs if a.part == "pairs" else part_accounts)(a, out_file)


if __name__ == "__main__":
    main()
