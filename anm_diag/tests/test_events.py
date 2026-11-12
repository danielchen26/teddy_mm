"""Typed source events: closure (rule S2), admission, and the forward/inverse join, on synthetic input and on
the finished records kept in runs/."""
import json
import os

import pytest

from anm_diag import events as ev
from anm_diag.checks import check, errors
from anm_diag.spec import load_declaration

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _decl(treated, stratum=None, pooling=None, name="t"):
    src = {"controlled": True, "kind": "crispr_ko", "control": ["NT"], "treated": treated, "stratum": stratum, "assign_from_guides": True}
    if pooling:
        src["pooling"] = pooling
    return {"name": name, "source": src}


def test_single_site_is_not_flagged():
    d = _decl(["GENE1"])
    assert ev.closure(ev.source_events(d)) == {"t:A=0": [], "t:A=1": []}
    assert ev.check_source(d) == []


def test_pooled_sites_are_flagged_until_declared():
    d = _decl(["GENE1", "GENE2"])
    assert ev.closure(ev.source_events(d))["t:A=1"] == ["site"]
    probs = ev.check_source(d)
    assert [p["code"] for p in probs] == ["S2"] and probs[0]["severity"] == "warning"
    assert ev.check_source(_decl(["GENE1", "GENE2"], pooling={"site": "one pathway; tested per target elsewhere"})) == []


def test_pooled_strata_are_flagged_separately():
    d = _decl(["GENE1"], stratum=["rep1", "rep2"])
    assert ev.closure(ev.source_events(d))["t:A=1"] == ["stratum"]
    assert "'stratum'" in ev.check_source(d)[0]["message"]


def test_label_pattern_counts_as_pooling():
    d = _decl("^RP[LS]\\d")
    evs = ev.source_events(d)
    assert evs[0]["site_is_pattern"] and ev.closure(evs)["t:A=1"] == ["site"]
    assert "pattern" in ev.check_source(d)[0]["message"]


def test_s2_is_a_warning_inside_check_and_never_an_error():
    decl = load_declaration(os.path.join(ROOT, "declarations", "teddy_mm_frangieh2021_headfix.sealed.yaml"))
    probs = check(decl)
    assert any(p["code"] == "S2" for p in probs)
    assert not any(p["code"] == "S2" for p in errors(probs))


def _records(ds):
    fwd = json.load(open(os.path.join(ROOT, "runs", "forward_typed_events", "forward_record.json")))
    inv = json.load(open(os.path.join(ROOT, "runs", f"{ds}_headfix", "anm_record.json")))
    decl = load_declaration(os.path.join(ROOT, "declarations", f"teddy_mm_{ds}_headfix.sealed.yaml"))
    return fwd, inv, decl


@pytest.mark.parametrize("ds,flags", [("frangieh2021", ["site"]), ("papalexi2021", ["site", "stratum"])])
def test_sealed_runs_pool_what_the_closure_says(ds, flags):
    fwd, inv, decl = _records(ds)
    evs = ev.source_events(decl, inv)
    assert sorted(ev.closure(evs)[f"{decl['name']}:A=1"]) == flags
    assert all(e["admission"]["n_admitted"] is None for e in evs if e["role"] == "treated")


def test_forward_events_carry_their_admission_counts():
    fwd, inv, decl = _records("frangieh2021")
    evs = ev.forward_events(fwd, "frangieh2021")
    by = {e["event_id"]: e for e in evs}
    assert by["forward:frangieh2021:STAT1:count:1.0"]["admission"]["n_admitted"] == 100
    jak2 = by["forward:frangieh2021:JAK2:embed:0.2"]["admission"]
    assert jak2["n_admitted"] < 10 and jak2["n_attempted"] == 100 and jak2["reason"]
    adm = ev.admission(fwd, "frangieh2021", decl["source"]["treated"])
    assert adm["STAT1"]["embed"] > 0.9 and adm["JAK2"]["embed"] < 0.1


def test_join_reads_the_same_source_from_both_halves():
    fwd, inv, decl = _records("frangieh2021")
    j = ev.join(fwd, inv, "frangieh2021", "PDL1")
    assert j["forward_readout"] == "CD274" and j["same_sign"]
    assert j["measured_delta"] == pytest.approx(-0.849, abs=1e-3)
    assert j["head_response"] == pytest.approx(-0.135, abs=1e-3)
    assert 0.1 < j["share_of_measured"] < 0.25 and j["step_difference"] < 0.25
    assert j["percentile_among_random"] == 100.0
    assert j["inverse"]["Z"]["verdict"] == "support" and j["inverse"]["HEAD_OUT"]["verdict"] == "reject"
    # the knockout direction is carried through the head more strongly than a random direction of the same size
    assert j["state_profile"]["HEAD_LAT"]["ko"] > 2 * j["state_profile"]["HEAD_LAT"]["random_median"]


def test_join_is_none_for_a_readout_the_forward_run_did_not_declare():
    fwd, inv, decl = _records("papalexi2021")
    assert ev.join(fwd, inv, "papalexi2021", "PDL2") is None
    assert ev.join(fwd, inv, "papalexi2021", "PDL1")["forward_readout"] == "CD274"


def _rec(delta, ci, states, n=100):
    return {"n_treated": n, "readouts": {"R": {"delta": delta, "delta_ci": ci, "states": {
        s: {"role": "chain", "u": u, "ci": [u - 0.03, u + 0.03], "verdict": v} for s, (u, v) in states.items()}}}}


def test_pooling_holds_when_sites_agree():
    pool = _rec(-0.8, [-0.9, -0.7], {"Z": (0.02, "support"), "OUT": (0.4, "reject")})
    sites = {g: _rec(-0.8 + d, [-0.95 + d, -0.65 + d], {"Z": (0.03, "support"), "OUT": (0.4, "reject")})
             for g, d in (("A", 0.0), ("B", 0.05))}
    e = ev.pooling_evidence(pool, sites, "R")
    assert e["label"] == "holds" and e["verdict_conflicts"] == [] and e["delta_heterogeneity"]["p"] > 0.05
    assert e["stages"]["OUT"]["agree_with_pool"] == ["A", "B"]


def test_pooling_differs_on_opposite_verdicts_or_direction():
    pool = _rec(-0.8, [-0.9, -0.7], {"Z": (0.02, "support")})
    a = _rec(-0.8, [-0.9, -0.7], {"Z": (0.02, "support")})
    b = _rec(-0.8, [-0.9, -0.7], {"Z": (0.5, "reject")})
    assert ev.pooling_evidence(pool, {"A": a, "B": b}, "R")["verdict_conflicts"] == ["Z"]
    c = _rec(+0.6, [0.4, 0.8], {"Z": (0.02, "support")})
    e = ev.pooling_evidence(pool, {"A": a, "C": c}, "R")
    assert e["label"] == "differs" and e["opposite_direction"] == ["C"] and e["delta_heterogeneity"]["p"] < 1e-6


def test_pooling_is_weak_when_a_site_does_not_respond():
    pool = _rec(-0.8, [-0.9, -0.7], {"Z": (0.02, "support")})
    a = _rec(-0.8, [-0.9, -0.7], {"Z": (0.02, "support")})
    n = _rec(-0.05, [-0.2, 0.1], {})          # the runner gives a site with no resolved response no stage verdicts
    e = ev.pooling_evidence(pool, {"A": a, "N": n}, "R")
    assert e["label"] == "weak" and e["no_response"] == ["N"] and list(e["stages"]["Z"]["sites"]) == ["A"]



def test_admission_trail_fills_the_real_events():
    d = _decl(["GENE1", "GENE2"])
    trail = {"GENE1": {"admitted": 198, "rejected_with_other_target": 268}, "GENE2": {"admitted": 278, "rejected_with_other_target": 321}}
    evs = {e["site"]: e for e in ev.source_events(d, trail=trail) if e["role"] == "treated"}
    assert evs["GENE1"]["admission"]["n_admitted"] == 198 and evs["GENE1"]["admission"]["n_attempted"] == 466
    assert "another target" in evs["GENE2"]["admission"]["reason"]


def test_s3_flags_a_source_on_the_readouts_own_gene():
    d = _decl(["CD274", "IFNGR1"]); d["source"]["pooling"] = {"site": "declared"}
    d["readouts"] = [{"id": "PDL1", "features": ["CD274"]}, {"id": "MHC1", "features": ["HLA_A"]}]
    probs = ev.check_source(d)
    assert [p["code"] for p in probs] == ["S3"] and "CD274" in probs[0]["message"] and probs[0]["severity"] == "warning"
    d["source"]["treated"] = ["IFNGR1", "HLA-A"]
    assert "MHC1" in ev.check_source(d)[0]["message"]


def test_ledger_reports_opposite_verdicts_and_what_differs():
    a = _rec(-0.8, [-0.9, -0.7], {"Z": (0.02, "support"), "OUT": (0.4, "reject")})
    b = _rec(-0.8, [-0.9, -0.7], {"Z": (0.5, "reject"), "OUT": (0.4, "reject")})
    for r in (a, b):
        for v in r["readouts"]["R"]["states"].values():
            v["best_class"] = "ridge"
    evs = ev.verdict_events(a, "run_a", site="pool", context="cells_1") + ev.verdict_events(b, "run_b", site="pool", context="cells_2")
    c = ev.ledger_conflicts(evs)
    assert len(c) == 1 and c[0]["stage"] == "Z" and c[0]["differ_in"] == ["context"] and c[0]["verdicts"] == ["support", "reject"]


def test_no_leakage_passes_for_guides_and_a_separate_measurement():
    d = _decl(["IFNGR1", "JAK1"]); d["system"] = {"inputs": ["rna_counts"]}
    d["readouts"] = [{"id": "PDL1", "measurement": "adt_counts", "features": ["PDL1"]}]
    evs = ev.source_events(d)
    assert all(e["no_leakage"]["status"] == "pass" for e in evs)
    assert set(evs[0]["no_leakage"]["routes"]) == set(ev.LEAK_ROUTES)
    assert not [p for p in ev.check_source({**d, "source": {**d["source"], "pooling": {"site": "declared"}}}) if p["code"] == "S4"]


def test_no_leakage_assignment_depends_on_where_the_labels_come_from():
    d = _decl(["GENE1"]); del d["source"]["assign_from_guides"]
    assert ev.no_leakage(d, "GENE1")["routes"]["assignment"]["status"] == "unknown"
    s4 = [p for p in ev.check_source(d) if p["code"] == "S4"]
    assert len(s4) == 1 and s4[0]["severity"] == "warning"
    d["source"]["label_basis"] = "design"
    assert ev.no_leakage(d, "GENE1")["status"] == "pass" and not [p for p in ev.check_source(d) if p["code"] == "S4"]
    d["source"]["label_basis"] = "response"
    assert ev.no_leakage(d, "GENE1")["status"] == "fail"
    s4 = [p for p in ev.check_source(d) if p["code"] == "S4"]
    assert len(s4) == 1 and s4[0]["severity"] == "error"


def test_no_leakage_flags_a_shared_measurement_and_the_readouts_own_gene():
    d = _decl(["CD274"]); d["system"] = {"inputs": ["rna_counts"]}; d["measurement_sharing"] = "declared_shared"
    d["readouts"] = [{"id": "R", "measurement": "rna_counts", "features": ["CD274"]}]
    n = ev.no_leakage(d, "CD274")
    assert n["routes"]["measurement"]["status"] == "flag" and n["routes"]["site"]["status"] == "flag" and n["status"] == "flag"
    d["measurement_sharing"] = "none"
    assert ev.no_leakage(d, "CD274")["routes"]["measurement"]["status"] == "fail"
    d["measurement_sharing"] = "count_split"; d["readouts"][0]["exclude_targets"] = True
    assert ev.no_leakage(d, "CD274")["status"] == "pass"
    assert ev.no_leakage(d, None)["routes"]["site"]["status"] == "pass"


def test_no_leakage_on_the_sealed_runs():
    pc = load_declaration(os.path.join(ROOT, "declarations", "teddy_mm_frangieh2021_headfix.sealed.yaml"))
    assert {e["no_leakage"]["status"] for e in ev.source_events(pc)} == {"pass"}
    ks = load_declaration(os.path.join(ROOT, "declarations", "ks2_teddy_ribo.yaml"))
    n = ev.source_events(ks)[0]["no_leakage"]["routes"]
    assert n["measurement"]["status"] == "flag" and n["assignment"]["status"] == "unknown" and n["site"]["status"] == "pass"
    fwd = json.load(open(os.path.join(ROOT, "runs", "forward_typed_events", "forward_record.json")))
    fe = {(e["site"], e["channel"]): e["no_leakage"] for e in ev.forward_events(fwd, "frangieh2021")}
    assert fe[("IFNGR1", "count")]["status"] == "pass" and fe[("HLA-E", "count")]["status"] == "flag"
