#!/usr/bin/env python3
"""The reading of batch 4, after forecasts.json was committed and the batch was exported. Each part is made once.

  python read.py --part pairs    --cells batch4_cells.tsv --states B2 B3 B4 --packs batch3_pack_A.npz batch3_pack_B.npz batch4_pack.npz
  python read.py --part rule     (the same) --library LIBRARY.npz
  python read.py --part accounts (the same) --counts B2_COUNTS B3_COUNTS B4_COUNTS

pairs -> results/reading.json
  1. Every forecast pair: the measured change in spreads of the readout among the control cells of its stratum and its
     error (as `Law.measure`); inside or outside the intervals of S2 and of S3; the errors of S0, S2 and S3.
  2. The edges listed in forecasts.json before the batch was opened, asked of its rows in each stratum. Those rows were
     forecast without them (`events.test_edges`, rule G2).
  Also the records the batch leaves (every unit on every readout), and every pair measured against the isotype
  antibodies instead of the mean of the 24 columns (it enters no hypothesis).
rule -> results/by_rule.json
  What the system does with the batch by rule and by nothing else: the rows are records; each stratum's law of version 6
  is fitted again with them; what it flags is proposed as an edge and asked of the library's row of the same pair,
  which was read under interferon and took no part in the flag.
accounts -> results/accounts.json
  3. The account (`account.account`) of every forecast pair whose change is at least four errors: at the input view V
     and the state Z, and in the co-culture also at the typed state T and its control R. No event of the batch has a
     type and no stratum is asked for `others`. Then `outer.open_parts`.

--rehearsal --results DIR: a rehearsal on the made-up file of rehearsal.py, to find faults of the code before the one
reading. It writes outside this run's results.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b4 as B  # noqa: E402
from _b4 import V, V6  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402
from anm_diag import outer as OU  # noqa: E402

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
RUNS = B.RUNS
TYPED = {"mel_co": "mel_co"}
PAIRS = (("Z", "T"), ("R", "T"), ("Z", "V"), ("T", "V"))
OUTPUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}   # readout -> output of the head as trained
SEED_SVD = 20261125


def inputs(a):
    fc = json.loads((B.RUN / "results" / "forecasts.json").read_text())
    assert fc["cells_file_sha256"] == sha(a.cells), "the forecasts were made for another list of cells"
    assert a.rehearsal or fc["packs_sha256"] == [sha(p) for p in a.packs[:-1]], "the forecasts stand on other antibody files"
    cells = B.table(a.cells)
    st, packs = B.States(a.states), B.Packs(a.packs)
    S = {}
    for g in B.STRATA:
        c = B.stratum(g, cells, st)
        assert packs.has(B.STUDY, c["ids"]).all(), "a cell of the list has no antibody columns"
        c["Y"], log, bg = packs.readouts(B.STUDY, c["ids"])
        c["against_the_background"] = {r: log[r] - bg for r in log}
        S[g] = c
    return fc, S, packs


def measure(y, m, C, sd=None):
    """As `Law.measure`: the change of the cells m against the control cells C in spreads of y among C, and its error."""
    centre = y[C].mean()
    sd = float((y[C] - centre).std()) if sd is None else sd
    v = y[m]
    return float((v.mean() - centre) / sd), float(np.sqrt(v.var(ddof=1) / len(v) + y[C].var() / C.sum()) / sd)


def part_pairs(a, out_file):
    fc, S, _packs = inputs(a)
    res = {"cells_file_sha256": fc["cells_file_sha256"], "pack_sha256": sha(a.packs[-1]), "rehearsal_on_made_up_values": bool(a.rehearsal), "min_snr": AC.MIN_SNR, "spread": {}, "records": {}, "rows": []}
    for g, c in S.items():
        C = c["kind"] == "control"
        res["spread"][g] = {r: float((y[C] - y[C].mean()).std()) for r, y in c["Y"].items()}
        res["records"][g] = {}
        for e in sorted(set(c["unit"][c["kind"] == "new"].tolist())):
            m = c["unit"] == e
            res["records"][g][e] = {"cells": int(m.sum()), "changes": {r: list(measure(y, m, C)) for r, y in c["Y"].items()}}
    for f in fc["rows"]:
        g, e, r = f["stratum"], f["event"], f["readout"]
        c = S[g]
        C, m = c["kind"] == "control", c["unit"] == e
        assert int(m.sum()) == f["cells"] and set(c["kind"][m].tolist()) == {"new"}, (g, e)
        ch, se = res["records"][g][e]["changes"][r]
        row = {k: f[k] for k in ("stratum", "event", "chosen", "readout", "cells", "edge", "type", "scale", "size_of_the_edge_from", "shift_times_noise")}
        row.update(change=ch, se=se, snr=abs(ch) / se, S0={"forecast": f["S0"]["forecast"], "error": ch - f["S0"]["forecast"]})
        for s in ("S2", "S3"):
            p = f[s]
            row[s] = {**{k: p[k] for k in p if k != "parts_of_the_spread"}, "error": ch - p["forecast"], "z": (ch - p["forecast"]) / p["sd"], "inside": bool(p["interval"][0] <= ch <= p["interval"][1])}
        cb, sb = measure(c["against_the_background"][r], m, C, sd=res["spread"][g][r])
        row["against_the_background"] = {"change": cb, "se": sb}
        res["rows"].append(row)
    assert len(res["rows"]) == fc["pairs"]
    to_test = fc["edges"]["to_test_on_the_rows_of_the_batch"]
    tested, have = [], set()
    for g in S:
        rows = [{"event": r["event"], "readout": r["readout"], "change": r["change"], "forecast": r["S3"]["forecast"], "sd": r["S3"]["sd"]} for r in res["rows"] if r["stratum"] == g and r["edge"] is None]
        for e in EV.test_edges(to_test, rows):
            if "z" in e.get("test", {}):
                e = dict(e, edge_id=f"{e['edge_id']}@{g}", stratum=g)
                e["test"]["scale_of_the_readout_here"] = fc["scale"][e["readout"]][g]
                e["test"]["direct_part"] = e["test"]["change"] - e["test"]["forecast_without_the_edge"]
                row = next(r for r in res["rows"] if r["stratum"] == g and r["event"] == e["site"] and r["readout"] == e["readout"])
                e["test"]["inside_the_interval_of_S3"] = row["S3"]["inside"]   # reported; the rule (G2) asks 1.96 spreads
                tested.append(e)
                have.add((e["site"], e["readout"]))
    res["edges"] = {"tested": tested, "with_no_row_in_the_batch": [e for e in to_test if (e["site"], e["readout"]) not in have]}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    rows = res["rows"]
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
    print(f"pairs {len(rows)} | inside S2 {sum(r['S2']['inside'] for r in rows)} S3 {sum(r['S3']['inside'] for r in rows)} | median z2 S2 {np.median([r['S2']['z'] ** 2 for r in rows]):.2f} "
          f"S3 {np.median([r['S3']['z'] ** 2 for r in rows]):.2f} | error S0 {rms([r['S0']['error'] for r in rows]):.3f} S2 {rms([r['S2']['error'] for r in rows]):.3f} S3 {rms([r['S3']['error'] for r in rows]):.3f} "
          f"| pairs with a change {sum(r['snr'] >= AC.MIN_SNR for r in rows)}")
    print("edges:", [(e["edge_id"], e["status"], round(e["test"]["z"], 1)) for e in tested], "| no row:", [e["edge_id"] for e in res["edges"]["with_no_row_in_the_batch"]])


def part_rule(a, out_file):
    fc, S, packs = inputs(a)
    R = json.loads((out_file.parent / "reading.json").read_text())
    assert R["pack_sha256"] == sha(a.packs[-1])
    d = np.load(a.library)
    lib_ids, un = d["ids"].astype(str), np.where(d["A"] == 0, "control", d["target"].astype(str))
    scale = fc["scale"]
    sc_of = lambda g: {r: scale[r][g] for r in B.READOUTS}
    sites = sorted((set(un.tolist()) | {u for g in S for u in S[g]["unit"].tolist()}) - {"control"})
    declared, tab = V.declared(sites)
    on_record = fc["edges"]["on_record"] + R["edges"]["tested"]
    admitted = [e for e in on_record if e["status"] == "admitted"]
    YL, _log, _bg = packs.readouts(B.STUDY, lib_ids)
    lib = V6.library_law(V.unit_length(d["Z"].astype(np.float64)), YL, un, tab, sc_of("library"))
    library_rows = lib["law"].forecast()                            # every pair of the library with its event left out, forecast without any proposed edge
    res = {"pack_sha256": R["pack_sha256"], "strata": {}, "edges": {"proposed": [], "tested_on_the_library": []}}
    known = declared + on_record + fc["edges"]["to_test_on_the_rows_of_the_batch"]
    for g, c in S.items():
        tab_g = {**tab, **{(e["site"], e["readout"]): e["class"] for e in admitted if e["stratum"] == g}}
        laws, _known, _first = V6.stratum_laws({"all": {"x": c["Z"], "Y": c["Y"], "unit": c["unit"], "table": tab_g, "scales": sc_of(g)}}, lib)
        law = laws["all"][0]
        new = set(c["unit"][c["kind"] == "new"].tolist())
        flags = law.flags()                                         # every flag of the law, of the units of the batch and of the units recorded before
        prop = EV.propose_edges(flags, known=known, stratum=g)      # a pair that has an edge, or was proposed before, is not proposed again
        tested = EV.test_edges(prop, library_rows)
        res["strata"][g] = {"recorded_events": len(law.events), "spread_with_no_edge": law.sigma0, "spread_from": law.spread_from, "width": law.conformal(),
                            "flags": [{**{k: f[k] for k in ("event", "readout", "flag", "change", "forecast", "z")}, "a_unit_of_the_batch": f["event"] in new} for f in flags]}
        res["edges"]["proposed"] += prop
        res["edges"]["tested_on_the_library"] += tested
        print(g, "| recorded events", len(law.events), "| spread", round(law.sigma0, 3), "| flags", [(f["event"], f["readout"], f["flag"], round(f["z"], 1)) for f in flags], flush=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("proposed", len(res["edges"]["proposed"]), "| on the library's rows:", [(e["edge_id"], e["status"]) for e in res["edges"]["tested_on_the_library"]])


def halves(g, c):
    """The halves, as in batches 2 and 3: the development cells keep those of run 40 (the first cells of the stratum in the
    list, in their order; checked against the manifest of the typed states), the other cells are split by cycle.halves."""
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
    fc, S, _packs = inputs(a)
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
    print(len(want), "pairs with a change of at least", AC.MIN_SNR, "errors:", {g: sum(r["stratum"] == g for r in want) for g in B.STRATA}, flush=True)
    head = None
    for g in B.STRATA:
        rows = [r for r in want if r["stratum"] == g and f"{g}|{r['event']}|{r['readout']}" not in res["accounts"]]
        if not rows:
            continue
        c = S[g]
        C = c["kind"] == "control"
        first, seed = halves(g, c)
        Vw, genes = input_view(a.counts, c["ids"], C)
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
        if any(r["readout"] in OUTPUT for r in rows):
            if head is None:
                sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
                from _common import forward_np, unit, weights
                P0, names, _j = weights()
                head = lambda raw: forward_np(P0, unit(raw))[1]
            o_all = head(c["raw"])
            out = {r: o_all[:, names.index(n)].astype(np.float64) for r, n in OUTPUT.items()}
        for r in sorted(rows, key=lambda r: (r["readout"], r["event"])):
            e, ro = r["event"], r["readout"]
            acc = AC.account(stages, c["Y"][ro], c["unit"], first, e, others=[], inputs=("V",), state="Z", o=out.get(ro), pairs=PAIRS, seed=seed)
            acc.update(stratum=g, event=e, readout=ro, chosen=r["chosen"], edge=r["edge"], type=r["type"], change_in_spreads=r["change"], snr=r["snr"])
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
    ap.add_argument("--part", choices=["pairs", "rule", "accounts"], required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--library", type=Path)
    ap.add_argument("--counts", type=Path, nargs="+")
    ap.add_argument("--results", type=Path, default=B.RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (B.RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out_file = a.results / {"pairs": "reading.json", "rule": "by_rule.json", "accounts": "accounts.json"}[a.part]
    if out_file.exists():
        raise SystemExit(f"{out_file.name} exists: this part of the reading was made once and is not made again")
    {"pairs": part_pairs, "rule": part_rule, "accounts": part_accounts}[a.part](a, out_file)


if __name__ == "__main__":
    main()
