"""Typed source events: one table for the sources of every kind of run.

The ANM paper builds its typed event graph through a source chain (SI S1.4):
propose -> prepare -> admit -> assign, and asks when an event representation is
sufficient (SI S1.4.1): a field of the event may be dropped only if no declared
readout tells its values apart. teddy_mm's earlier bridge implemented that chain
with a schema; the sufficiency runs here declared their source as a list of
treated labels, and the forward run declared events of its own. This module
writes all of them as the same record, so that three things can be checked:

  closure    which event fields a pooled source drops (rule S2 in checks.py)
  no leakage whether the readout's measurement could have entered an event (rules R2, S3, S4)
  admission  for how many cells an event could be applied at all
  join       whether the forward and the inverse half agree on the same source

An event has: time, site, channel, sign, amplitude (the fields a verifier can
depend on), plus role, stratum, provenance, the pool it is merged into, and an
admission entry. Nothing here touches cell data: events are built from
declarations and finished records.
"""
import json
import re
import statistics

VERIFIER_FIELDS = ("time", "site", "channel", "sign", "amplitude")
POOLABLE = ("site", "stratum", "channel")       # fields a sufficiency source may merge
FORWARD_READOUT = {"PDL1": "CD274", "MHC1": "HLA-A-B-C", "HLAE": "HLA-E"}   # inverse readout id -> head output


def _as_list(x):
    if x is None:
        return []
    return list(x) if isinstance(x, (list, tuple)) else [x]


def _is_pattern(label):
    return isinstance(label, str) and bool(re.search(r"[\^\$\[\]\\\*\+\?\(\)\|]", label))


LEAK_ROUTES = ("assignment", "measurement", "site")
_LEAK_ORDER = {"pass": 0, "flag": 1, "unknown": 2, "fail": 3}


def _norm_gene(s):
    return str(s).upper().replace("_", "-")


def _site_in_readout(site, feats):
    for f in feats:
        if _is_pattern(f):
            if site == f or (not _is_pattern(site) and re.search(f, str(site))):
                return True
        elif not _is_pattern(site) and _norm_gene(site) == _norm_gene(f):
            return True
    return False


def no_leakage(decl, site=None):
    """No-leakage entry of one source event of a sufficiency declaration.

    The paper keeps the endpoint away from the source: hashes are fixed before the verifier runs (SI S7.1.1),
    and teddy_mm's earlier bridge rejects any event that uses the hidden endpoint. Here the hidden endpoint
    is the measurement a readout is computed from. It can reach a source event by three routes, and each
    route gets a status and a reason:

      assignment   how the cells carrying the event were chosen. pass: from detected guides, or from labels
                   declared as fixed by the design (source.label_basis: design). fail: from classes derived
                   from the measured response (source.label_basis: response). unknown: a dataset label whose
                   origin is not declared.
      measurement  whether a readout's measurement is also an input of the audited system (rule R2).
                   pass: separate measurement, or count-split. flag: declared_shared. fail: shared, undeclared.
      site         whether the event acts on a readout's own gene (rule S3). flag unless the readout leaves
                   the targeted gene out (exclude_targets).

    status is the worst route: fail > unknown > flag > pass. `site=None` gives the entry of a control event.
    """
    src = decl.get("source") or {}
    routes = {}
    if src.get("assign_from_guides"):
        routes["assignment"] = ("pass", "cells are assigned from their detected guides, a measurement that enters no readout")
    else:
        basis = src.get("label_basis")
        if basis == "design":
            routes["assignment"] = ("pass", "labels are declared as fixed by the experimental design")
        elif basis == "response":
            routes["assignment"] = ("fail", "labels are derived from the measured response: the choice of cells depends on the outcome")
        else:
            routes["assignment"] = ("unknown", "a dataset label whose origin is not declared (source.label_basis: design | response)")
    sysinp = set(_as_list((decl.get("system") or {}).get("inputs")))
    pol = decl.get("measurement_sharing", "none")
    worst = ("pass", "no readout uses an input of the audited system")
    for r in decl.get("readouts") or []:
        rm = r.get("measurement")
        if rm in sysinp:
            if pol == "count_split":
                cur = ("pass", f"readout '{r.get('id')}' shares the input '{rm}', count-split")
            elif pol == "declared_shared":
                cur = ("flag", f"readout '{r.get('id')}' shares the system input '{rm}' (declared_shared): counting noise is shared")
            else:
                cur = ("fail", f"readout '{r.get('id')}' is computed from the system input '{rm}' with no sharing policy")
            if _LEAK_ORDER[cur[0]] > _LEAK_ORDER[worst[0]]:
                worst = cur
    routes["measurement"] = worst
    if site is None:
        routes["site"] = ("pass", "a control event has no site")
    else:
        own = [r for r in decl.get("readouts") or [] if _site_in_readout(site, _as_list(r.get("features")))]
        if not own:
            routes["site"] = ("pass", "the event does not act on a readout's own gene")
        elif all(r.get("exclude_targets") for r in own):
            routes["site"] = ("pass", "the site is among a readout's genes, and the readout leaves the targeted gene out (exclude_targets)")
        else:
            routes["site"] = ("flag", "the event acts on the own gene of readout " + ", ".join(str(r.get("id")) for r in own) + " (rule S3)")
    status = max((v[0] for v in routes.values()), key=_LEAK_ORDER.get)
    return {"status": status, "routes": {k: {"status": v[0], "reason": v[1]} for k, v in routes.items()}}


def forward_no_leakage(site, channel, readout_genes=()):
    """No-leakage entry of an in-silico event. Its amplitude and sign are fixed in the declaration and it is
    applied to control cells, so no measured readout can enter it; the readout is the model's own output.
    An event on a readout's own gene is flagged, as in rule S3."""
    how = ("the direction is estimated from measured RNA embeddings (treated against control cells); no readout measurement is used"
           if channel == "z_direction" else "set in silico on control cells; amplitude and sign are fixed in the declaration")
    routes = {"assignment": ("pass", how),
              "measurement": ("pass", "the readout is the model's own output, not a measurement"),
              "site": (("flag", "the event acts on a readout's own gene: a direct path is expected")
                       if any(_norm_gene(site) == _norm_gene(g) for g in readout_genes) else ("pass", "the event does not act on a readout's own gene"))}
    status = max((v[0] for v in routes.values()), key=_LEAK_ORDER.get)
    return {"status": status, "routes": {k: {"status": v[0], "reason": v[1]} for k, v in routes.items()}}


def source_events(decl, record=None, trail=None):
    """Events of a sufficiency (inverse) declaration: one per treated label and stratum, pooled into A=1.

    `trail` is an optional admission trail {site: {"admitted": n, "rejected_with_other_target": m, ...}} written
    when the cells were selected; with it every event says how many cells carrying it were admitted and why
    the others were not.
    """
    src = decl["source"]
    name = decl.get("name", "run")
    strata = _as_list(src.get("stratum")) or [None]
    treated = _as_list(src.get("treated"))
    prov = "guide assignment" if src.get("assign_from_guides") else "dataset label"
    n_treated = (record or {}).get("n_treated")
    out = []
    for site in treated:
        for st in strata:
            tr = (trail or {}).get(site)
            if tr and len(strata) == 1:
                adm = dict(status="admitted", n_admitted=tr["admitted"], n_attempted=tr["admitted"] + tr.get("rejected_with_other_target", 0),
                           reason="cells that also carry a guide for another target are not admitted")
            else:
                adm = dict(status="admitted", n_admitted=None, n_attempted=None,
                           reason=(f"the record is pooled ({n_treated} treated cells); cells per event are not recorded"
                                   if n_treated is not None else "no record given"))
            out.append(dict(
                event_id=f"{name}:{src.get('kind', 'source')}:{site}" + (f":{st}" if st is not None else ""),
                origin=f"inverse:{name}", time=0, site=site, site_is_pattern=_is_pattern(site),
                channel=src.get("kind", "source"), sign=-1, amplitude=1.0, role="treated", stratum=st,
                provenance=prov, pool=f"{name}:A=1", admission=adm, no_leakage=no_leakage(decl, site)))
    n_ctrl = None
    if record and "n_cells" in record and n_treated is not None:
        n_ctrl = record["n_cells"] - n_treated
    for st in strata:
        out.append(dict(
            event_id=f"{name}:control" + (f":{st}" if st is not None else ""), origin=f"inverse:{name}", time=0,
            site=None, site_is_pattern=False, channel=src.get("kind", "source"), sign=0, amplitude=0.0, role="control",
            stratum=st, provenance="control labels " + ", ".join(map(str, _as_list(src.get("control_guides") or src.get("control")))),
            pool=f"{name}:A=0", admission=dict(status="admitted", n_admitted=n_ctrl if len(strata) == 1 else None,
                                               n_attempted=None, reason=""), no_leakage=no_leakage(decl, None)))
    return out


def forward_events(record, dataset):
    """Events of a finished forward run: one per site, channel and amplitude, with its admission count."""
    d = record["datasets"][dataset]
    n = d["n_confirmatory"]
    rg = _as_list((((record.get("declaration") or {}).get("events") or {}).get("sites") or {}).get("readout_genes"))
    out = []
    for key, e in d["sites"].items():
        site, channel = key.split("|")
        for eps, v in e["by_eps"].items():
            k = v["n_cells"]
            out.append(dict(
                event_id=f"forward:{dataset}:{site}:{channel}:{eps}", origin=f"forward:{dataset}", time=0, site=site,
                site_is_pattern=False, channel=channel, sign="+/-", amplitude=float(eps), role="treated",
                stratum="control cells", provenance="in silico", pool=None, no_leakage=forward_no_leakage(site, channel, rg),
                admission=dict(status="admitted" if k else "skipped", n_admitted=k, n_attempted=n,
                               reason="" if k == n else "zero count, not among the cell's tokens, or outside the vocabulary")))
    return out


def closure(events, fields=VERIFIER_FIELDS + ("stratum",)):
    """For every pool, the fields whose values differ among the pooled events (or are an unresolved pattern)."""
    out = {}
    for pool in sorted({e["pool"] for e in events if e.get("pool")}):
        evs = [e for e in events if e["pool"] == pool]
        flagged = [f for f in fields if len({json.dumps(e.get(f), sort_keys=True) for e in evs}) > 1]
        if any(e.get("site_is_pattern") for e in evs) and "site" not in flagged:
            flagged.append("site")
        out[pool] = flagged
    return out


def check_source(decl):
    """Rules S2 and S3 on the source of a sufficiency declaration.

    S2: a source that pools events differing in site, stratum or channel must say so.
    S3: a source that acts on the gene of a declared readout is flagged (warning).
    S4: the cells carrying a source event must not be chosen from the measured response (error); a dataset
        label whose origin is not declared is a warning (declare source.label_basis: design | response).

    Declared under source.pooling, e.g.  pooling: {site: "one pathway; tested per target in run X"}.
    Undeclared pooling is a warning, not an error: the run still answers its question for the pooled
    source, but the claim does not extend to each pooled value (SI S1.4.1).
    """
    src = decl.get("source") or {}
    if not src.get("treated"):
        return []
    declared = src.get("pooling") or {}
    flags = closure(source_events(decl)).get(f"{decl.get('name', 'run')}:A=1", [])
    out = []
    for f in flags:
        if f in POOLABLE and f not in declared:
            vals = sorted({str(e.get(f)) for e in source_events(decl) if e["role"] == "treated"})
            what = "a label pattern" if (f == "site" and any(_is_pattern(v) for v in vals)) else f"{len(vals)} values"
            out.append({"code": "S2", "severity": "warning",
                        "message": f"the source pools events that differ in '{f}' ({what}: {', '.join(vals[:6])}"
                                   f"{', ...' if len(vals) > 6 else ''}) without declaring it under source.pooling: "
                                   f"the verdicts hold for the pooled source, not for each {f}"})
    asg = no_leakage(decl)["routes"]["assignment"]
    if asg["status"] == "fail":
        out.append({"code": "S4", "severity": "error", "message": "the source is assigned from the measured response (source.label_basis: "
                    "response): cells chosen by their outcome make every later verdict depend on that outcome"})
    elif asg["status"] == "unknown":
        out.append({"code": "S4", "severity": "warning", "message": "the source uses a dataset label whose origin is not declared: say "
                    "source.label_basis: design (fixed by the experiment) or response (derived from the measured response)"})
    norm = lambda s: str(s).upper().replace("_", "-")
    sites = {norm(s) for s in _as_list(src.get("treated")) if not _is_pattern(s)}
    for r in decl.get("readouts") or []:
        own = sorted(sites & {norm(f) for f in _as_list(r.get("features"))})
        if own:
            out.append({"code": "S3", "severity": "warning",
                        "message": f"the source acts on the readout's own gene ({', '.join(own)} for readout '{r.get('id')}'): a knockout "
                                   f"can remove the protein while its RNA barely changes, so an RNA-only system may have nothing to "
                                   f"read; a rejection at the input view then describes the input, not the system"})
    return out


def admission(record, dataset, sites=None, channels=("count", "embed")):
    """Share of confirmatory cells in which each site could be perturbed, by channel."""
    d = record["datasets"][dataset]
    n = d["n_confirmatory"]
    out = {}
    for key, e in d["sites"].items():
        site, ch = key.split("|")
        if ch in channels and (sites is None or site in sites):
            out.setdefault(site, {})[ch] = e["n_cells"] / n
    if sites:
        for s in sites:
            out.setdefault(s, {})
            for ch in channels:
                out[s].setdefault(ch, 0.0)
    return out


def _slope(entry, readout):
    r = entry["readout"][readout]
    return r["slope_mean"], list(r["slope_ci"])


def site_responses(record, dataset, sites, readout, eps=None):
    """The head's response to in-silico events at each site (largest declared amplitude unless eps is given)."""
    d = record["datasets"][dataset]
    rows = []
    for site in sites:
        for ch in ("count", "embed"):
            e = d["sites"].get(f"{site}|{ch}")
            if not e or not e["n_cells"]:
                rows.append(dict(site=site, channel=ch, n=0))
                continue
            key = eps if eps is not None else max(e["by_eps"], key=float)
            v = e["by_eps"][key]
            m, ci = _slope(v, readout)
            rows.append(dict(site=site, channel=ch, eps=float(key), n=v["n_cells"], slope=m, ci=ci,
                             resolved=bool(ci[0] > 0 or ci[1] < 0), z_move=v["state"]["Z"]["mean"]))
    return rows


def distinguishable(rows):
    """Pairs of pooled sites whose response intervals do not overlap, per channel."""
    out = []
    ok = [r for r in rows if r.get("n")]
    for i, a in enumerate(ok):
        for b in ok[i + 1:]:
            if a["channel"] == b["channel"] and (a["ci"][1] < b["ci"][0] or b["ci"][1] < a["ci"][0]):
                out.append((a["site"], b["site"], a["channel"]))
    return out


def join(forward_record, inverse_record, dataset, readout):
    """Set the forward response along the measured knockout direction beside the inverse verdicts.

    The forward number is the head's own output (log mean); the measured change is in the readout's
    statistic. The two are not in identical units, so their ratio is an order of magnitude only.
    """
    d = forward_record["datasets"][dataset]
    fro = FORWARD_READOUT.get(readout, readout)
    if fro not in d["readouts"]:
        return None
    ko = d["sites"]["KO_direction|z_direction"]["by_eps"]
    big, small = max(ko, key=float), min(ko, key=float)
    s1, ci = _slope(ko[big], fro)
    s2, _ = _slope(ko[small], fro)
    rnd = [k for k in d["sites"] if k.startswith("random_dir")]
    rnd_slopes = [abs(_slope(d["sites"][k]["by_eps"][big], fro)[0]) for k in rnd]
    stages = [s for s, v in ko[big]["state"].items() if v and v.get("mean") == v.get("mean") and v.get("mean") is not None]
    profile = {s: dict(ko=ko[big]["state"][s]["mean"],
                       random_median=statistics.median(d["sites"][k]["by_eps"][big]["state"][s]["mean"] for k in rnd))
               for s in stages}
    R = inverse_record["readouts"][readout]
    return dict(dataset=dataset, readout=readout, forward_readout=fro,
                measured_delta=R["delta"], measured_ci=list(R["delta_ci"]),
                head_response=s1, head_response_ci=ci, half_step=s2,
                step_difference=abs(s1 - s2) / max(abs(s1), abs(s2), 1e-12),
                same_sign=bool(s1 * R["delta"] > 0), share_of_measured=s1 / R["delta"],
                random_abs_median=statistics.median(rnd_slopes), random_abs_max=max(rnd_slopes), n_random=len(rnd),
                percentile_among_random=d["reference_z_direction"]["KO_direction|z_direction"][fro]["percentile_among_random"],
                ko_direction_norm=d["z_direction_norm"], state_profile=profile,
                inverse={s: dict(u=v["u"], ci=list(v["ci"]), verdict=v["verdict"], best_class=v["best_class"])
                         for s, v in R["states"].items() if v.get("role") == "chain"})


def pooling_evidence(pooled, per_site, readout):
    """Does a pooled source hide differences between the pooled sites? (the measured side of rule S2)

    `pooled` is the record of the pooled run; `per_site` maps each site to the record of the same declaration
    with that site alone as the treated class (same controls). Returns, for one readout, each site's measured
    change, a heterogeneity test across sites, and for every chain stage each site's verdict next to the pool's.

      holds    every site responds in the pool's direction and no two sites have opposite resolved verdicts
      weak     no opposite verdicts, but some site's own response is not resolved from zero (the runner
               gives such a site no stage verdicts; it is listed under no_response and left out of the stages)
      differs  a site responds in the other direction, or two sites have opposite resolved verdicts at a stage

    The heterogeneity statistic is Cochran's Q with standard errors read off the interval widths, so its
    p-value is approximate; the verdict comparison does not depend on it.
    """
    from scipy.stats import chi2
    z = 1.959963984540054
    sites = list(per_site)
    P = pooled["readouts"][readout]

    def het(vals, cis):
        ses = [max((c[1] - c[0]) / (2 * z), 1e-12) for c in cis]
        w = [1 / s ** 2 for s in ses]
        m = sum(wi * v for wi, v in zip(w, vals)) / sum(w)
        q = sum(wi * (v - m) ** 2 for wi, v in zip(w, vals))
        k = len(vals) - 1
        return dict(Q=q, df=k, p=float(chi2.sf(q, k)) if k > 0 else 1.0, weighted_mean=m, n=len(vals))

    delta = {s: dict(v=per_site[s]["readouts"][readout]["delta"], ci=list(per_site[s]["readouts"][readout]["delta_ci"]),
                     n_treated=per_site[s]["n_treated"]) for s in sites}
    responds = {s: bool(d["ci"][0] > 0 or d["ci"][1] < 0) for s, d in delta.items()}
    same_dir = {s: bool(responds[s] and d["v"] * P["delta"] > 0) for s, d in delta.items()}
    opposite = [s for s, d in delta.items() if responds[s] and d["v"] * P["delta"] < 0]
    staged = [s for s in sites if per_site[s]["readouts"][readout].get("states")]
    stages, conflicts = {}, []
    for st, pv in P["states"].items():
        if pv.get("role") != "chain":
            continue
        per = {s: per_site[s]["readouts"][readout]["states"][st] for s in staged}
        if not per:
            continue
        verd = {s: v["verdict"] for s, v in per.items()}
        if "support" in verd.values() and "reject" in verd.values():
            conflicts.append(st)
        stages[st] = dict(pooled=dict(u=pv["u"], ci=list(pv["ci"]), verdict=pv["verdict"]),
                          sites={s: dict(u=v["u"], ci=list(v["ci"]), verdict=v["verdict"]) for s, v in per.items()},
                          heterogeneity=het([v["u"] for v in per.values()], [v["ci"] for v in per.values()]),
                          agree_with_pool=[s for s in staged if verd[s] == pv["verdict"]])
    label = "differs" if (opposite or conflicts) else ("holds" if all(same_dir.values()) else "weak")
    return dict(readout=readout, label=label, pooled_delta=P["delta"], pooled_delta_ci=list(P["delta_ci"]), delta=delta,
                delta_heterogeneity=het([d["v"] for d in delta.values()], [d["ci"] for d in delta.values()]),
                responds=responds, no_response=[s for s in sites if not responds[s]], opposite_direction=opposite,
                verdict_conflicts=conflicts, stages=stages)


def verdict_events(record, run, site="pool", context=None):
    """The verdicts of a finished sufficiency record as typed events: one per readout and chain stage."""
    out = []
    for rid, R in record["readouts"].items():
        for st, v in (R.get("states") or {}).items():
            if v.get("role") == "chain":
                out.append(dict(event_id=f"verdict:{run}:{rid}:{st}", origin=f"verdict:{run}", run=run, readout=rid, site=st,
                                source_site=site, context=context, channel=v["best_class"], sign={"support": 1, "reject": -1}.get(v["verdict"], 0),
                                amplitude=v["u"], ci=list(v["ci"]), verdict=v["verdict"]))
    return out


def ledger_conflicts(events):
    """Pairs of runs with opposite resolved verdicts (support against reject) at the same readout and stage,
    with the fields in which the two runs differ: a difference in source_site or context is a dependence to
    report, not a contradiction; identical fields on both sides would be one."""
    out = []
    key = lambda e: (e["readout"], e["site"])
    by = {}
    for e in events:
        by.setdefault(key(e), []).append(e)
    for (rid, st), evs in sorted(by.items()):
        for i, a in enumerate(evs):
            for b in evs[i + 1:]:
                if a["sign"] * b["sign"] == -1:
                    out.append(dict(readout=rid, stage=st, runs=[a["run"], b["run"]], verdicts=[a["verdict"], b["verdict"]], u=[a["amplitude"], b["amplitude"]],
                                    differ_in=[f for f in ("source_site", "context") if a.get(f) != b.get(f)]))
    return out


def write_jsonl(events, path):
    with open(path, "w") as f:
        for e in events:
            f.write(json.dumps(e, ensure_ascii=False) + "\n")
