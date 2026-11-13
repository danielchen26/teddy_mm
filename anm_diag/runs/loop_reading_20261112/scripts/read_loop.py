"""Re-read the three rounds of the loop on the teddy_mm head under LOOP_PROTOCOL.md, with anm_diag/loop.py.

No fit and no audit is made here. The script reads the records and bootstrap draws of runs/loop_20261112,
runs/loop_guard_20261112 and runs/repair_20261112 (anm_record.json and anm_bootstrap.npz in each audit folder; the draws
are not kept in git and come back by re-running the audits) and writes results/loop_reading.json and results/ledger.jsonl.

  python read_loop.py

Post hoc and exploratory: the outcomes those rounds recorded under their own rules stand as recorded.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
from anm_diag import loop  # noqa: E402

DS = {"frangieh2021": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "papalexi2021": ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]}
SPAN = ("Z", "HEAD_OUT")                                        # from the stage before the located step to the output
STEPS = [("Z", "HEAD_LAT"), ("HEAD_LAT", "HEAD_OUT")]
ARMS = {"as_trained": "loop_20261112", "decoder_only": "loop_20261112", "mlp_only": "loop_20261112", "mlp_guarded": "loop_guard_20261112", "joint_s0": "repair_20261112"}
key = lambda s: f"{s[0]}->{s[1]}"


def run(folder, name):
    return loop.Run.load(RUNS / folder / name, name)


def pub(d):
    return {k: v for k, v in d.items() if k != "draws"}


out = {"protocol": "LOOP_PROTOCOL.md", "span": list(SPAN), "grade": "post hoc re-reading, exploratory",
       "admission": {}, "reference_in_the_sealed_runs": {}, "diagnosis": {"pooled": {}, "single_gene": {}}, "closure": {}, "round3_splits": {},
       "round3_against_round2": {}, "transfer": {}, "guard": {}}
ledger = [loop.intervention_event("head/as_trained", None, "head", "none: the phase-1 head as trained", "baseline")]
for arm, change, role, comp, fixed in (
        ("decoder_only", "decoder fine-tuned on measured PD-L1 of knockout cells, one target gene held out", "control", "head.decoder", ["head.mlp"]),
        ("mlp_only", "MLP fine-tuned on measured PD-L1 of knockout cells, one target gene held out", "located", "head.mlp", ["head.decoder"]),
        ("mlp_guarded", "as mlp_only, with the head's outputs held on bone-marrow training cells", "located", "head.mlp", ["head.decoder"]),
        ("joint", "as mlp_guarded, one fit for both cell types and every measured protein with a head output", "located", "head.mlp", ["head.decoder"])):
    ledger.append(loop.intervention_event(f"head/{arm}", "head/as_trained", comp, change, role, held_fixed=fixed, fitted_on="knockout cells of the other target genes and half of the control cells"))

for ds, genes in DS.items():
    P = {a: run(f, f"pooled_{ds}_{a}" if a != "joint_s0" else f"pooled_{ds}_joint_s0") for a, f in ARMS.items()}
    base = P["as_trained"]
    for a, r in P.items():
        ledger.append(loop.audit_event(r, "head/" + ("joint" if a == "joint_s0" else a)))
    # gate 0, on the audit cells of the loop
    out["admission"][ds] = {rid: loop.admit(base, rid, SPAN) for rid in base.record["readouts"]}
    ok = [rid for rid, a in out["admission"][ds].items() if a["diagnostic"] == "ready"]
    for rid, a in out["admission"][ds].items():
        ledger.append(loop.decision_event("admission", a, [base.label]))
    # the same reference stage in the sealed runs, which used every control cell
    sealed = loop.Run(json.loads((RUNS / f"{ds}_headfix" / "anm_record.json").read_text()), None, f"{ds}_headfix")
    out["reference_in_the_sealed_runs"][ds] = {"n_cells": sealed.record["n_cells"], "n_treated": sealed.record["n_treated"],
                                               "stages": {rid: {s: loop.stage(sealed, rid, s) for s in ("Z", "HEAD_LAT", "HEAD_OUT")} for rid in sealed.readouts()}}
    # the diagnostic loop: MLP (located) against decoder (control), pooled and one held-out gene at a time
    out["diagnosis"]["pooled"][ds] = {rid: {key(s): loop.diagnose(base, P["mlp_only"], P["decoder_only"], rid, s) for s in [SPAN] + STEPS} for rid in ok}
    out.setdefault("steps", {})[ds] = {rid: {a: {key(s): pub(loop.added(P[a], rid, s)) for s in [SPAN] + STEPS} for a in P} for rid in ok}
    for rid in ok:
        ledger.append(loop.decision_event("diagnosis", out["diagnosis"]["pooled"][ds][rid][key(SPAN)], [base.label, P["mlp_only"].label, P["decoder_only"].label]))
    out["diagnosis"]["single_gene"][ds] = {}
    for g in genes:
        b, l, c = (run("loop_20261112", f"target_{ds}_{g}_{a}") for a in ("as_trained", "mlp_only", "decoder_only"))
        adm = loop.admit(b, "PDL1", SPAN)                              # a gene that fails gate 0 is not read further
        out["diagnosis"]["single_gene"][ds][g] = {"admission": adm, "diagnosis": loop.diagnose(b, l, c, "PDL1", SPAN) if adm["diagnostic"] == "ready" else {"outcome": "not admitted"}}
    # the repair loop: closure of every changed head against the head as trained
    out["closure"][ds] = {}
    for a in ("decoder_only", "mlp_only", "mlp_guarded", "joint_s0"):
        out["closure"][ds][a] = {rid: {key(s): loop.closure(P[a], base, rid, s) for s in [SPAN] + STEPS} for rid in ok}
        for rid in ok:
            ledger.append(loop.decision_event("closure", out["closure"][ds][a][rid][key(SPAN)], [P[a].label, base.label]))
    # round 3 on its three splits (no head as trained was audited on splits 1 and 2, so the span is read alone)
    out["round3_splits"][ds] = {}
    for s in (0, 1, 2):
        r = run("repair_20261112", f"pooled_{ds}_joint_s{s}")
        out["round3_splits"][ds][str(s)] = {rid: {"end": loop.stage(r, rid, SPAN[1]), "reference": loop.stage(r, rid, SPAN[0]), **{key(x): pub(loop.added(r, rid, x)) for x in [SPAN] + STEPS}} for rid in ok}
    # round 3 against round 2 on the same cells: does the joint head add more than the per-cell-type head?
    out["round3_against_round2"][ds] = {rid: pub(loop.contrast(P["mlp_guarded"], P["joint_s0"], rid, SPAN)) for rid in ok}
    # a head repaired on the other cell type (round 1), audited here on every cell
    out["transfer"][ds] = {}
    for a in ("mlp_only", "decoder_only"):
        r = run("loop_20261112", f"transfer_{ds}_{a}")
        out["transfer"][ds][a] = {"n_cells": r.record["n_cells"], "PDL1": {"end": loop.stage(r, "PDL1", SPAN[1]), "reference": loop.stage(r, "PDL1", SPAN[0]), key(SPAN): pub(loop.added(r, "PDL1", SPAN))}}

# gate 6: the head's own task, from the rounds' result files
c1 = json.load(open(RUNS / "loop_20261112" / "results" / "cost_bone_marrow.json")); c2 = json.load(open(RUNS / "loop_guard_20261112" / "results" / "bone_marrow.json"))
c3 = json.load(open(RUNS / "repair_20261112" / "results" / "bone_marrow.json")); b0 = c1["as_trained"]["mean_pearson"]
groups = {f"{ds}|{arm}": [v["mean_pearson"] for k, v in c1["heads"].items() if k.startswith(ds) and k.endswith(arm)] for ds in DS for arm in ("decoder_only", "mlp_only")}
groups.update({f"{ds}|mlp_guarded": [v["mean_pearson"] for k, v in c2["heads"].items() if k.startswith(ds)] for ds in DS})
groups["both|joint_s0"] = [v["mean_pearson"] for k, v in c3["heads"].items() if k.startswith("s0__") and k != "s0__all"]
out["guard"] = {"as_trained": b0, "allowed_fall": 0.02, "groups": {k: {"lowest": min(v), "highest": max(v), **loop.guard(b0, sum(v) / len(v), 0.02)} for k, v in groups.items()}}

(HERE.parent / "results").mkdir(exist_ok=True)
json.dump(out, open(HERE.parent / "results" / "loop_reading.json", "w"), indent=1, default=float)
n = loop.write_ledger(ledger, HERE.parent / "results" / "ledger.jsonl")
f3 = lambda iv: f"[{iv[0]:+.3f}, {iv[1]:+.3f}]"
for ds in DS:
    print("==", ds)
    for rid, a in out["admission"][ds].items():
        print(f"  admission {rid}: diagnostic {a['diagnostic']}, certification {a['certification']}")
    for rid, d in out["diagnosis"]["pooled"][ds].items():
        x = d[key(SPAN)]
        print(f"  diagnosis {rid}: {x['outcome']}; baseline adds {x['baseline']['added']:.3f} {f3(x['baseline']['ci'])}; MLP removes {x['located']['share_removed']['point']:.2f} {f3(x['located']['share_removed']['ci'])}; "
              f"decoder {x['control']['share_removed']['point']:.2f} {f3(x['control']['share_removed']['ci'])}; difference {f3(x['located_minus_control']['ci'])}; selective {x['selective']}")
    print("  single genes:", {g: v["diagnosis"]["outcome"] for g, v in out["diagnosis"]["single_gene"][ds].items()})
    for a, v in out["closure"][ds].items():
        print(f"  closure {a}:", {rid: {k: c["grade"] for k, c in x.items()} for rid, x in v.items()})
    print("  round 3 splits, PD-L1 span:", {s: (round(v["PDL1"][key(SPAN)]["added"], 3), f3(v["PDL1"][key(SPAN)]["ci"]), v["PDL1"][key(SPAN)]["reading"]) for s, v in out["round3_splits"][ds].items()})
    c = out["round3_against_round2"][ds]["PDL1"]
    print(f"  round 3 against round 2, PD-L1: round 2 removes {c['removed']:+.3f} {f3(c['removed_ci'])} more than round 3 ({c['reading']}); paired {c['paired']}")
print("ledger events:", n)
