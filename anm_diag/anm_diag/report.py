"""Markdown report of a run record."""

CAUSE_TEXT = {
    "no_effect": "No resolved readout response: nothing to explain (data).",
    "control_failed": "Positive control not supported: the test cannot recover known information; no verdict.",
    "reference_unresolved": "The reference stage is not supported. Rejected: the input view itself does not "
                            "reproduce the response for any declared readout class (a limit of the data or the "
                            "classes). Unresolved: the data cannot tell. Either way nothing is charged to the system.",
    "loss_at": "Representation loss: a stage is rejected after a supported stage.",
    "relative_loss_at": "Additional loss: the reference does not fully explain the readout, but one step of the "
                        "chain adds a resolved loss on top (a weaker, relative claim).",
    "unresolved": "Unresolved: later stages neither supported nor rejected.",
    "no_deficit": "No deficit: every stage supported.",
}


def _f(x, d=2):
    return f"{x:+.{d}f}" if isinstance(x, (int, float)) else str(x)


def _iv(iv, d=2):
    return f"[{iv[0]:+.{d}f}, {iv[1]:+.{d}f}]"


def write_report(rec, path):
    d = rec["declaration"]
    L = []
    L.append(f"# ANM diagnosis: {d['name']}")
    L.append("")
    L.append(f"- status: **{d.get('status')}**" + (" (exploratory: margins and design were not sealed before these data)"
                                                   if rec["exploratory"] else ""))
    L.append(f"- system: {d['system']['id']} (frozen) · inputs {d['system'].get('inputs')}")
    L.append(f"- source: {d['source'].get('kind', '')} · treated `{d['source']['treated']}` vs `{d['source']['control']}`"
             f" · {rec['n_treated']} treated cells in {rec['n_treated_clusters']} clusters, {rec['n_cells']} cells in all")
    L.append(f"- chain: {' → '.join(d['chain'])}")
    L.append(f"- margin u = {d['decision']['margin_u']}, alpha = {d['decision']['alpha']}; "
             f"declaration sha256 `{rec['declaration_sha256'][:16]}…`; anm_diag {rec['anm_diag_version']}")
    probs = [p for p in rec["problems"]]
    if probs:
        L.append("")
        L.append("**Checks**")
        for p in probs:
            L.append(f"- [{p['code']}] {p['severity']}: {p['message']}")
    L.append("")
    L.append("## Forward: does the source reach each state? (held-out-cluster AUC)")
    L.append("")
    L.append("| state | AUC |")
    L.append("|---|---|")
    for k, v in rec["forward"]["auc"].items():
        L.append(f"| {k} | {v:.3f} |")
    L.append("")
    L.append("| step | change, 95% interval |")
    L.append("|---|---|")
    for k, v in rec["forward"]["steps"].items():
        L.append(f"| {k} | {_iv(v, 3)} |")
    for rid, r in rec["readouts"].items():
        L.append("")
        L.append(f"## Readout {rid}")
        L.append("")
        L.append(f"Response Delta = {r['delta']:+.3f} {_iv(r['delta_ci'], 3)}")
        if not r["states"]:
            rc = r["root_cause"]
            L.append("")
            L.append(f"**Root cause: {rc['cause']}.** {CAUSE_TEXT[rc['cause']]} No inverse test is made.")
            continue
        L.append("")
        L.append("| state | role | u (best class) | 95% interval | verdict | all classes |")
        L.append("|---|---|---|---|---|---|")
        for s, v in r["states"].items():
            allc = "; ".join(f"{c} {_f(x['u'])} {_iv(x['ci'])}" for c, x in v["by_class"].items())
            L.append(f"| {s} | {v['role']} | {_f(v['u'])} ({v['best_class']}) | {_iv(v['ci'])} | **{v['verdict']}** | {allc} |")
        L.append("")
        L.append("| step | change in u (best class) | same class |")
        L.append("|---|---|---|")
        for k, v in r["steps"].items():
            sc = "; ".join(f"{c} {_iv(x)}" for c, x in v["same_class"].items())
            L.append(f"| {k} | {_iv(v['best_class'])} | {sc} |")
        rc = r["root_cause"]
        L.append("")
        L.append(f"**Root cause: {rc['cause']}.** {CAUSE_TEXT[rc['cause']]} {rc.get('detail', '')}")
        if r.get("exact"):
            e = r["exact"]
            L.append("")
            L.append("**Exact same-state test**: " + "; ".join(f"{k} = {v}" for k, v in e.items()))
    L.append("")
    L.append("## Provenance")
    for k, v in rec["provenance"].items():
        L.append(f"- {k}: `{v}`")
    L.append(f"- runtime: {rec['runtime_s']} s")
    path.write_text("\n".join(L) + "\n")
