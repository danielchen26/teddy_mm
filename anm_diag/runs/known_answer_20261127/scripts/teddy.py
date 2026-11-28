#!/usr/bin/env python3
"""Question 2: what TEDDY carries of the part that does not add, at its input, at its state and at its head's output.
From what is cached: the states, the head as trained, the RNA counts. No run of the encoder. Once, after analyze.py.

  python teddy.py --cells final_cells.tsv --library library.npz --states MEL_STATES... --counts MEL_COUNTS...
                  --counts-ifn counts_ifn.npz --packs PACK... [--results DIR] [--rehearsal]

-> results/teddy.json

For every refused knockout of the composition knockout x IFN-γ (if more than twelve: the members of the classes, then
the strongest readouts, up to twelve) and for PD-L1, MHC-I, HLA-E, every other protein that is a witness and its ten strongest
RNA witnesses: `account.account`, unchanged, in each of the two conditions, with the halves of the run, at V (512
components of the log counts of the cells of the condition that enter this run) and at the state Z; for PD-L1, MHC-I
and HLA-E also the head's output as trained. `others`: what a readout taught by the other members of the knockout's
class reads. A pair the graph declares as direct is left out, as in question 1.

What a place holds of a change is its share times the change. A change of less than four errors has no account: it and
what is held of it are taken as zero. The share of the part that does not add held at a place is
(held under IFN-γ - held with no signal) / (change under IFN-γ - change with no signal).
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ka as K  # noqa: E402
from analyze import PROTEIN  # noqa: E402
from anm_diag import account as AC  # noqa: E402

RUNS = K.RUNS
OUTPUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}   # melanoma readout -> output of the head as trained (as in the final audit)
MOST, GENES = 12, 10
PLACES = {"input": ("V", "carried"), "Z": ("Z", "carried"), "Z_read_by_a_readout_of_control_cells": ("Z", "read"), "Z_read_by_a_readout_taught_by_its_class": ("Z", "others"), "head_output": ("Z", "as_trained")}


def by_path(name, path):
    """A script of another run, loaded by its path (several runs have scripts of the same name)."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def chosen(test):
    """The refused knockouts whose accounts are made, and the class of each."""
    cls = {g: k + 1 for k, c in enumerate(test["classes"]) for g in c["members"]}
    ref = [r for r in test["rows"] if r["verdict"] == "refused"]
    if len(ref) > MOST:
        ref = sorted(ref, key=lambda r: (r["knockout"] not in cls, -r["strongest_readout"]["value"]))[:MOST]
    return {r["knockout"]: r for r in ref}, cls


def readouts_of(row):
    left = {d["readout"] for d in row["declared_pairs_left_out"]}
    prot = [PROTEIN + r for r in OUTPUT] + [w["readout"] for w in row["witnesses"]["proteins"]]
    return [r for r in dict.fromkeys(prot + [w["readout"] for w in row["witnesses"]["genes"][:GENES]]) if r not in left]


def gene_values(files, ids, genes):
    """log1p of counts per 10,000 of the named genes, cell by cell; the first file that holds a cell is used."""
    want = {i: k for k, i in enumerate(ids)}
    out, seen = {g: np.zeros(len(ids)) for g in genes}, np.zeros(len(ids), bool)
    for f in files:
        d = np.load(f)
        names = d["genes"].astype(str).tolist()
        k = np.array([want.get(i, -1) for i in d["ids"].astype(str)])
        take = (k >= 0) & ~np.r_[seen, False][k]
        if not take.any():
            continue
        rows, cols, vals = d["rows"], d["cols"], d["vals"].astype(np.float64)
        tot = np.maximum(np.bincount(rows, weights=vals, minlength=len(k)), 1.0)
        for g in genes:
            m = cols == names.index(g)
            v = np.zeros(len(k))
            v[rows[m]] = np.log1p(vals[m] / tot[rows[m]] * 1e4)
            out[g][k[take]] = v[take]
        seen[k[take]] = True
        del d, rows, cols, vals
    assert seen.all(), "a cell has no RNA counts in the files"
    return out


def condition(k, C, ko, a, genes, packs):
    """The cells of one condition that enter the run (control cells and the knockouts of the table), with V, Z, the
    readouts and the head's output."""
    FA = by_path("final_audit_read", RUNS / "final_audit_20261127" / "scripts" / "read.py")       # the input view, exactly as the final audit made it
    TG = by_path("teg_common", RUNS / "teg_training_20261119" / "scripts" / "_common.py")          # the head as trained
    forward_np, unit_, weights = TG.forward_np, TG.unit, TG.weights
    c = C[k]
    m = np.isin(c["unit"], ["control"] + ko)
    ids, unit, first = c["ids"][m], c["unit"][m], c["A"][m]
    if k == "ifn":
        L = np.load(a.library)
        pos = {i: j for j, i in enumerate(L["ids"].astype(str))}
        raw = L["Z"][np.array([pos[i] for i in ids])].astype(np.float64)
    else:
        raw = K.B3.States(a.states).rows(ids)
    files = a.counts_ifn if k == "ifn" else a.counts
    Vw, n_genes = FA.input_view(files, ids, unit == "control")
    clr, _log, _bg = packs.readouts(K.STUDY, ids)
    P0, names, _j = weights()
    o_all = forward_np(P0, unit_(raw))[1]
    Y = {PROTEIN + r: v for r, v in clr.items()}
    Y.update(gene_values(files, ids.tolist(), genes))
    return {"unit": unit, "first": first, "V": Vw, "Z": K.B3.V.unit_length(raw), "Y": Y, "o": {PROTEIN + r: o_all[:, names.index(n)].astype(np.float64) for r, n in OUTPUT.items()},
            "genes_of_the_input_view": n_genes, "cells": {"control": int((unit == "control").sum()), "all": int(len(unit))}}


def one(c, g, r, mates):
    """The account of knockout g on readout r in one condition, or None where the change is under four errors."""
    m = np.isin(c["unit"], ["control", g] + mates)
    o = c["o"].get(r)
    acc = AC.account({"V": c["V"][m], "Z": c["Z"][m]}, c["Y"][r][m], c["unit"][m], c["first"][m], g, others=mates, inputs=("V",), state="Z", o=None if o is None else o[m], pairs=(("Z", "V"),), seed=K.SEED, controls=False)
    z = acc["stages"]["Z"]
    out = {"change": z["change"], "errors": z["snr"], "cells": z["cells"], "account": z["snr"] >= AC.MIN_SNR}
    if out["account"]:
        out["shares"] = {p: acc["stages"][s]["shares"][key] for p, (s, key) in PLACES.items() if key in acc["stages"][s]["shares"]}
        out["carried_by"] = {s: acc["stages"][s]["carried_by"] for s in ("V", "Z")}
        out["parts"] = {p: v["share"] for p, v in acc["chain"]["parts"].items()}
        out["says"] = acc["chain"]["says"]
    return out


def held(pair):
    """Per place: what is held of the part that does not add, and that part. None where neither condition has an account."""
    a, b = pair["IFNγ"], pair["Control"]
    if not (a["account"] or b["account"]):
        return None
    N = (a["change"] if a["account"] else 0.0) - (b["change"] if b["account"] else 0.0)
    h = {}
    for p in PLACES:
        v = [x["shares"][p]["share"] * x["change"] if x["account"] and p in x["shares"] else (0.0 if not x["account"] else None) for x in (a, b)]
        h[p] = None if v[0] is None or v[1] is None else v[0] - v[1]
    return {"part_that_does_not_add": N, "held": h, "share": {p: (None if v is None else v / N) for p, v in h.items()}}


def pooled(pairs, seed):
    """Sums over the pairs of a group; the interval by resampling its knockouts."""
    out = {}
    kos = sorted({p["knockout"] for p in pairs})
    rng = np.random.default_rng(seed)
    draws = [kos] + [[kos[i] for i in rng.integers(0, len(kos), len(kos))] for _ in range(K.RESAMPLES)]
    for place in PLACES:
        use = [p for p in pairs if p["non_additive"] and p["non_additive"]["held"][place] is not None]
        if not use:
            continue
        by = {}
        for p in use:
            N = p["non_additive"]["part_that_does_not_add"]
            s = by.setdefault(p["knockout"], [0.0, 0.0])
            s[0] += np.sign(N) * p["non_additive"]["held"][place]; s[1] += abs(N)
        val = [sum(by[g][0] for g in d if g in by) / max(sum(by[g][1] for g in d if g in by), 1e-300) for d in draws]
        out[place] = {"share": float(val[0]), "ci": [float(x) for x in np.percentile(val[1:], [2.5, 97.5])] if len(by) > 1 else None, "pairs": len(use), "knockouts": len(by)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--counts", type=Path, nargs="+", required=True)
    ap.add_argument("--counts-ifn", type=Path, nargs="+", required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--results", type=Path, default=K.RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (K.RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out_file = a.results / "teddy.json"
    if out_file.exists():
        raise SystemExit("teddy.json exists: this part of the reading was made once and is not made again")
    R = json.loads((a.results / "reading.json").read_text())
    assert R["rehearsal_on_made_up_values"] == bool(a.rehearsal)
    test = R["tests"]["ifn"]
    rows, cls = chosen(test)
    want = {g: readouts_of(r) for g, r in rows.items()}
    genes = sorted({r for v in want.values() for r in v if not r.startswith(PROTEIN)})
    C = K.cells(a.cells, a.library)
    ko = K.knockouts(C, "ifn")
    res = {"rehearsal_on_made_up_values": bool(a.rehearsal), "knockouts": list(rows), "min_errors": AC.MIN_SNR, "head_outputs": OUTPUT, "conditions": {}, "pairs": []}
    acc, packs = {}, K.B3.Packs(a.packs)
    for k in ("ifn", "ctl"):
        c = condition(k, C, ko, a, genes, packs)
        res["conditions"][K.CONDITION[k]] = {"cells": c["cells"], "genes_of_the_input_view": c["genes_of_the_input_view"]}
        for g, rs in want.items():
            mates = [h for h in rows if h != g and cls.get(h) is not None and cls.get(h) == cls.get(g)]
            for r in rs:
                acc[(g, r, k)] = one(c, g, r, mates)
                x = acc[(g, r, k)]
                print(f"{K.CONDITION[k]:8s} {g:8s} -> {r:16s} change {x['change']:+.3f} ({x['errors']:.1f} errors)" + (f" | {x['says']}" if x["account"] else " | no account"), flush=True)
        del c
    for g, rs in want.items():
        for r in rs:
            p = {"knockout": g, "readout": r, "kind": "protein" if r.startswith(PROTEIN) else "gene", "class": cls.get(g), "IFNγ": acc[(g, r, "ifn")], "Control": acc[(g, r, "ctl")]}
            p["non_additive"] = held(p)
            res["pairs"].append(p)
    groups = {}
    for p in res["pairs"]:
        who = f"class {p['class']}" if p["class"] else "refused, in no class"
        what = {PROTEIN + "PDL1": "PD-L1", PROTEIN + "MHC1": "MHC-I", PROTEIN + "HLAE": "HLA-E"}.get(p["readout"], "other proteins" if p["kind"] == "protein" else "genes")
        for key in {(who, what), (who, "the three proteins the head has an output for" if what in ("PD-L1", "MHC-I", "HLA-E") else what)}:
            groups.setdefault(" | ".join(key), []).append(p)
    res["pooled"] = {k: pooled(v, K.SEED) for k, v in sorted(groups.items())}
    out_file.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n")
    for k, v in res["pooled"].items():
        print(k, {p: (round(x["share"], 2), x["ci"] and [round(y, 2) for y in x["ci"]], x["pairs"]) for p, x in v.items()})


if __name__ == "__main__":
    main()
