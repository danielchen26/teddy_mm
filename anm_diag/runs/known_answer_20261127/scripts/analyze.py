#!/usr/bin/env python3
"""Question 1, the data alone: does every knockout compose with the signal by addition? The rule of DESIGN.md, once.

  python analyze.py --summary summary.npz [--results DIR] [--rehearsal]

-> results/reading.json, results/table_ifn.tsv, results/table_co.tsv, results/outcome.json

The verdict of a knockout is by its readouts: refused if a readout stands outside in both halves with one sign. The
sum T over all readouts is reported beside it and decides nothing (DESIGN.md, changed after the rehearsal).

Refuses to run if reading.json exists. --rehearsal --results DIR: on made-up sums, never into the results of the run.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _ka as K  # noqa: E402
from anm_diag import events as EV  # noqa: E402

LOG2 = float(np.log(2.0))
PROTEIN = "protein:"                                           # a protein readout is named protein:PDL1, a gene by its symbol


def levels(S, k):
    """Per group and half of condition k: the level of every readout (genes, then proteins) and, from the control cells
    of the condition, the variance of that level for one cell."""
    n = S[f"{k}_n"].astype(float)[:, :, None]
    m = S[f"{k}_rna_sum"] / n
    nc = S[f"{k}_n"][0].sum()
    mc, qc = S[f"{k}_rna_sum"][0].sum(0) / nc, S[f"{k}_rna_sq"][0].sum(0) / nc
    v_rna = (qc - mc ** 2) * nc / (nc - 1) / (mc + K.EPS) ** 2
    out = {"mean_ctl": mc}
    for name in ("prot", "iso"):
        pc, pq = S[f"{k}_{name}_sum"][0].sum(0) / nc, S[f"{k}_{name}_sq"][0].sum(0) / nc
        out[name] = (np.concatenate([np.log(m + K.EPS), S[f"{k}_{name}_sum"] / n], 2), np.concatenate([v_rna, (pq - pc ** 2) * nc / (nc - 1)]))
    return out


def departure(S, test, scale="prot"):
    """u[g, half, readout]: the departure from addition over its error from cell counts. Also the departure itself, the
    effect of the signal on the control cells, the names of the readouts and the pairs left out for each knockout."""
    a, b = K.TESTS[test]
    ko = [str(g) for g in S[f"knockouts_{test}"]]
    La, Lb = levels(S, a), levels(S, b)
    genes, prot = S["genes"].astype(str), [str(p) for p in S["proteins"]]
    keep = np.r_[(La["mean_ctl"] >= K.MIN_MEAN) | (Lb["mean_ctl"] >= K.MIN_MEAN), np.ones(len(prot), bool)]
    names = np.array(list(genes) + [PROTEIN + p for p in prot])[keep]
    row = lambda k: np.array([[str(u) for u in S[f"{k}_units"]].index(g) for g in ko])
    (xa, va), (xb, vb) = La[scale], Lb[scale]
    xa, xb, va, vb = xa[:, :, keep], xb[:, :, keep], va[keep], vb[keep]
    f = xa[0] - xb[0]                                                                    # the signal on the control cells, per half
    I = (xa[row(a)] - xa[0]) - (xb[row(b)] - xb[0])
    gene = ~np.char.startswith(names, PROTEIN)
    I[:, :, gene] -= np.median(I[:, :, gene], 2, keepdims=True)                          # RNA levels are relative: the factor common to all genes of a group is no departure
    f[:, gene] -= np.median(f[:, gene], 1, keepdims=True)
    s2 = va[None, None, :] / S[f"{a}_n"][row(a)][:, :, None] + vb[None, None, :] / S[f"{b}_n"][row(b)][:, :, None]
    ok = (s2 > 0).all((0, 1))
    left = np.zeros((len(ko), len(names)), bool)
    dec = K.declared()
    for i, g in enumerate(ko):
        left[i] = np.isin(names, [g] + [PROTEIN + r for r in dec.get(g, [])])
    return {"ko": ko, "names": names[ok], "u": I[:, :, ok] / np.sqrt(s2[:, :, ok]), "I": I[:, :, ok], "f": f[:, ok], "left": left[:, ok],
            "cells": {g: {K.CONDITION[k]: S[f"{k}_n"][r].tolist() for k, r in ((a, ra), (b, rb))} for g, ra, rb in zip(ko, row(a), row(b))}}


def standardize(u, idx=None):
    """z: centred by the median over the knockouts and scaled by their median absolute deviation, per readout and half.
    idx: the knockouts the medians are taken over (a resample); every knockout is standardized."""
    w = u if idx is None else u[idx]
    med = np.median(w, 0)
    mad = 1.4826 * np.median(np.abs(w - med), 0)
    return (u - med) / np.where(mad > 0, mad, np.inf)


def judge(z, left, idx=None):
    """The two numbers of every knockout and the two ends from the mirror side (over the knockouts idx)."""
    zA, zB = np.where(left, 0.0, z[:, 0]), np.where(left, 0.0, z[:, 1])
    P = zA * zB
    T, m = P.sum(1), np.minimum(np.abs(zA), np.abs(zB))
    rows = slice(None) if idx is None else idx
    tau = float(max(np.max(-T[rows]), 0.0))
    w = float(max(K.FOUR, np.max(np.where(P[rows] < 0, m[rows], 0.0))))
    wit = (P > 0) & (m > w)
    return {"T": T, "m": np.where(P > 0, m, -m), "tau": tau, "w": w, "witness": wit, "by_T": T > tau, "by_readout": wit.any(1)}


def removed(D, i):
    """The share of the signal's effect that knockout i removes: its departure in one half projected on the signal's
    effect of the other half, over the readouts the signal changes two-fold in both halves, genes and proteins apart."""
    fA, fB = D["f"]
    out = {}
    for key, sel in (("rna", ~np.char.startswith(D["names"], PROTEIN)), ("protein", np.char.startswith(D["names"], PROTEIN))):
        r = sel & (np.abs(fA) >= LOG2) & (np.abs(fB) >= LOG2) & (fA * fB > 0) & ~D["left"][i]
        den = float((fA[r] * fB[r]).sum())
        out[key] = {"readouts": int(r.sum()), "share": float(-(D["I"][i, 0, r] * fB[r] + D["I"][i, 1, r] * fA[r]).sum() / (2 * den)) if den > 0 else None}
    return out


def classes(D, z, J, refused):
    """The refused knockouts grouped by `events.propose_types`, unchanged: each is a flag; its direction is its z (mean
    of the halves) over the readouts that are a witness of any refused knockout (every readout if there is none)."""
    if not refused:
        return [], []
    W = J["witness"][refused].any(0)
    W = W if W.any() else np.ones(len(W), bool)
    zm = np.where(D["left"], 0.0, z.mean(1))
    flags = [{"flag": "a direction the law does not cover", "event": D["ko"][i], "readout": r} for i in refused for r in (D["names"][J["witness"][i]].tolist() or ["(the whole profile)"])]
    types = EV.propose_types(flags, {D["ko"][i]: zm[i, W] for i in refused})
    return [{"members": t["members"], "cosine_inside": t["evidence"]["cosine_of_the_shifts"], "largest_cosine_with_a_knockout_outside": t["evidence"]["largest_cosine_with_an_event_outside"]} for t in types], D["names"][W].tolist()


def one_test(S, test, draws):
    D = departure(S, test)
    ko, names = D["ko"], D["names"]
    z = standardize(D["u"])
    J = judge(z, D["left"])
    ref = J["by_readout"]                                                               # the verdict: a readout stands outside. The sum T is reported beside it
    refused = [i for i in range(len(ko)) if ref[i]]
    Di = departure(S, test, "iso")
    assert (Di["names"] == names).all()
    zi = standardize(Di["u"])
    rng = np.random.default_rng(K.SEED)
    share, taus, ws = np.zeros(len(ko)), [], []
    for _ in range(draws):
        idx = rng.integers(0, len(ko), len(ko))
        Jb = judge(standardize(D["u"], idx), D["left"], idx)
        share += Jb["by_readout"]
        taus.append(Jb["tau"]); ws.append(Jb["w"])
    cls, used = classes(D, z, J, refused)
    of = {g: k + 1 for k, c in enumerate(cls) for g in c["members"]}
    isp = np.char.startswith(names, PROTEIN)
    rows = []
    for i, g in enumerate(ko):
        wi = np.where(J["witness"][i])[0]
        wi = wi[np.argsort(-np.abs(J["m"][i, wi]))]
        wit = [{"readout": str(names[r]), "z": [float(z[i, 0, r]), float(z[i, 1, r])], "departure": float(D["I"][i, :, r].mean()),
                **({"also_against_the_isotypes": bool(zi[i, 0, r] * zi[i, 1, r] > 0 and min(abs(zi[i, 0, r]), abs(zi[i, 1, r])) >= K.FOUR)} if isp[r] else {})} for r in wi]
        best = int(np.argmax(J["m"][i]))
        rows.append({"knockout": g, "cells": D["cells"][g], "verdict": "refused" if ref[i] else "carries over",
                     "T": float(J["T"][i]), "T_over_tau": float(J["T"][i] / J["tau"]) if J["tau"] > 0 else None, "the_whole_profile_stands_out": bool(J["by_T"][i]),
                     "strongest_readout": {"readout": str(names[best]), "value": float(J["m"][i, best]), "over_w": float(J["m"][i, best] / J["w"])},
                     "share_of_resamples_refused": float(share[i] / draws) if draws else None,
                     "witnesses": {"n": int(len(wi)), "genes_down": int(sum(1 for r in wi if not isp[r] and z[i, 0, r] < 0)), "genes_up": int(sum(1 for r in wi if not isp[r] and z[i, 0, r] > 0)),
                                   "proteins": [w for w in wit if w["readout"].startswith(PROTEIN)], "genes": [w for w in wit if not w["readout"].startswith(PROTEIN)]},
                     "share_of_the_signal_removed": removed(D, i),
                     "declared_pairs_left_out": [{"readout": str(names[r]), "z": [float(z[i, 0, r]), float(z[i, 1, r])], "beyond_w_in_both_halves": bool(z[i, 0, r] * z[i, 1, r] > 0 and min(abs(z[i, 0, r]), abs(z[i, 1, r])) > J["w"])}
                                                 for r in np.where(D["left"][i])[0]],
                     "class": of.get(g), "of_the_five": g in K.FIVE})
    five = [r for r in rows if r["of_the_five"]]
    other = [r for r in rows if not r["of_the_five"]]
    n_ref = sum(r["verdict"] == "refused" for r in rows)
    hit = sum(r["verdict"] == "refused" for r in five)
    best = max(cls, key=lambda c: len(set(c["members"]) & set(K.FIVE)), default=None)
    in_best = len(set(best["members"]) & set(K.FIVE)) if best else 0
    carried = sum(r["verdict"] == "carries over" for r in other)
    lines = {"P1": {"says": "all five are refused", "refused": hit, "of": len(five), "holds": hit == len(five) == len(K.FIVE)},
             "P2": {"says": "at least 90% of the other knockouts carry over", "carry_over": carried, "of": len(other), "share": carried / len(other), "holds": carried / len(other) >= 0.9},
             "P3": {"says": "one class holds all five and at most two others", "class": best["members"] if best else None, "of_the_five_in_it": in_best, "others_in_it": (len(best["members"]) - in_best) if best else None,
                    "holds": bool(best) and in_best == len(K.FIVE) and len(best["members"]) - in_best <= 2}}
    pct = lambda v: [float(x) for x in np.percentile(v, [2.5, 97.5])] if len(v) else None
    return {"conditions": [K.CONDITION[k] for k in K.TESTS[test]], "knockouts": len(ko), "readouts": {"genes": int((~isp).sum()), "proteins": int(isp.sum())},
            "signal_changes_two_fold": {"genes": int(((np.abs(D["f"]) >= LOG2).all(0) & (D["f"][0] * D["f"][1] > 0) & ~isp).sum()), "proteins": int(((np.abs(D["f"]) >= LOG2).all(0) & (D["f"][0] * D["f"][1] > 0) & isp).sum())},
            "tau": J["tau"], "w": J["w"], "resamples": {"n": draws, "tau": pct(taus), "w": pct(ws)},
            "refused": [r["knockout"] for r in rows if r["verdict"] == "refused"], "whole_profile_stands_out": [r["knockout"] for r in rows if r["the_whole_profile_stands_out"]], "classes": cls, "readouts_the_classes_stand_on": len(used),
            "against_the_five": {"refusals": {"recall": hit / len(K.FIVE), "precision": hit / n_ref if n_ref else None, "refused": n_ref},
                                 "class": {"recall": in_best / len(K.FIVE), "precision": in_best / len(best["members"]) if best else None, "members": len(best["members"]) if best else 0}},
            "lines": lines, "passed": all(v["holds"] for v in lines.values()), "rows": rows}


def table(res, path):
    head = ["knockout", "verdict", "strongest readout", "its value/w", "refused in resamples", "witness genes (down/up)", "witness proteins", "strongest witness genes",
            "signal removed (RNA)", "(protein)", "class", "whole profile T/tau", "cells no signal (A,B)", "cells with signal (A,B)"]
    f3 = lambda v: "" if v is None else f"{v:.2f}"
    lines = ["\t".join(head)]
    for r in sorted(res["rows"], key=lambda r: (r["verdict"] != "refused", -r["strongest_readout"]["over_w"])):
        w = r["witnesses"]
        c = list(r["cells"].values())
        lines.append("\t".join([r["knockout"], r["verdict"], r["strongest_readout"]["readout"], f3(r["strongest_readout"]["over_w"]), f3(r["share_of_resamples_refused"]),
                                f"{w['n'] - len(w['proteins'])} ({w['genes_down']}/{w['genes_up']})", ", ".join(p["readout"][len(PROTEIN):] for p in w["proteins"]), ", ".join(x["readout"] for x in w["genes"][:8]),
                                f3(r["share_of_the_signal_removed"]["rna"]["share"]), f3(r["share_of_the_signal_removed"]["protein"]["share"]), str(r["class"] or ""), f3(r["T_over_tau"]), str(c[1]), str(c[0])]))
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=K.RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    ap.add_argument("--draws", type=int, default=K.RESAMPLES)
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (K.RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out = a.results / "reading.json"
    if out.exists():
        raise SystemExit("reading.json exists: the reading was made once and is not made again")
    S = dict(np.load(a.summary))
    res = {"summary_sha256": K.sha(a.summary), "rehearsal_on_made_up_values": bool(a.rehearsal), "rule": {"min_mean": K.MIN_MEAN, "eps": K.EPS, "four": K.FOUR, "seed": K.SEED, "min_cosine": EV.MIN_COSINE, "min_class": EV.MIN_TYPE},
           "the_five": list(K.FIVE), "tests": {t: one_test(S, t, a.draws) for t in K.TESTS}}
    a.results.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    for t, r in res["tests"].items():
        table(r, a.results / f"table_{t}.tsv")
    main_ = res["tests"]["ifn"]
    (a.results / "outcome.json").write_text(json.dumps({"rehearsal_on_made_up_values": bool(a.rehearsal), "the_check": "knockout x IFN-γ", "passed": main_["passed"], "lines": main_["lines"],
                                                        "against_the_five": main_["against_the_five"], "refused": main_["refused"], "classes": [c["members"] for c in main_["classes"]], "whole_profile_stands_out": main_["whole_profile_stands_out"],
                                                        "secondary_co_culture": {k: res["tests"]["co"][k] for k in ("knockouts", "refused", "against_the_five")} | {"classes": [c["members"] for c in res["tests"]["co"]["classes"]]}},
                                                       indent=1, ensure_ascii=False) + "\n")
    for t, r in res["tests"].items():
        print(f"{t}: {r['knockouts']} knockouts, {r['readouts']} | w {r['w']:.2f} {r['resamples']['w']} | refused {len(r['refused'])}: {r['refused']} | whole profile stands out (described; tau {r['tau']:.0f}): {r['whole_profile_stands_out']}")
        print("   classes:", [c["members"] for c in r["classes"]], "|", {k: v["holds"] for k, v in r["lines"].items()}, "| passed" if r["passed"] else "| NOT passed", "|", r["against_the_five"])


if __name__ == "__main__":
    main()
