"""Gate 5 of runs/round1_coverage_20261116 (scripts/read_final.py): how the outcome of the round is stated, the
difference between two groups of donors, and the refusal to read files that are not those of the declaration."""
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "runs" / "round1_coverage_20261116" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import common  # noqa: E402,F401
import plan  # noqa: E402

pytest.importorskip("torch")
import read_final  # noqa: E402
from anm_diag import loop, scope  # noqa: E402

HYP = loop.hypothesis_event("h", {"event": "e", "readout": "r", "target": "variation", "quantity": "q"}, "coverage",
                            {"coverage": [0.04, 0.06], "function": [-0.01, 0.01], "amount": [-0.05, 0.05]},
                            {"own": "coverage", "at_least": 0.034, "controls": {"function": 0.034, "amount": 0.034}}, ["basis"])


def _family(gain_by_seed, n=40, noise=0.004, seed=0):
    g = np.random.default_rng(seed)
    U = {"as_is": 0.6 + g.normal(0, 0.02, (n, 1)) + g.normal(0, noise, (n, 3))}
    for k, by in gain_by_seed.items():
        U[k] = U["as_is"] - np.asarray(by)[None, :] + g.normal(0, noise, (n, 3))
    return scope.arms_family(U, base="as_is"), U


def test_the_outcome_is_stated_only_if_every_seed_reads_the_same():
    ok = {k: [True] * 3 for k in ("coverage", "function", "amount")}
    fam, _ = _family({"coverage": [0.06, 0.06, 0.06], "function": [0.0, 0.0, 0.0], "amount": [0.0, 0.0, 0.0]})
    t = read_final.tests_of(HYP, fam, ok, ("coverage", "function", "amount"))
    assert t["together"]["outcome"] == "confirmed" and [x["outcome"] for x in t["by_seed"]] == ["confirmed"] * 3
    assert t["outcome_of_the_round"] == "confirmed"
    fam, _ = _family({"coverage": [0.08, 0.08, 0.034], "function": [0.0, 0.0, 0.0], "amount": [0.0, 0.0, 0.0]})     # one seed sits on the bound
    t = read_final.tests_of(HYP, fam, ok, ("coverage", "function", "amount"))
    assert t["together"]["outcome"] == "confirmed" and t["by_seed"][2]["outcome"] == "unresolved"
    assert t["outcome_of_the_round"] == "unresolved"
    fam, _ = _family({"coverage": [0.0, 0.0, 0.0], "function": [0.0, 0.0, 0.0], "amount": [0.0, 0.0, 0.0]})
    t = read_final.tests_of(HYP, fam, ok, ("coverage", "function", "amount"))
    assert t["outcome_of_the_round"] == "refuted"


def test_an_arm_that_was_not_admitted_is_not_tested():
    fam, _ = _family({"coverage": [0.06, 0.06, 0.06], "function": [0.0, 0.0, 0.0], "amount": [0.0, 0.0, 0.0]})
    adm = {"coverage": [True, False, True], "function": [True] * 3, "amount": [True] * 3}
    t = read_final.tests_of(HYP, fam, adm, ("coverage", "function", "amount"))
    assert t["together"]["outcome"] == "not testable" and t["by_seed"][1]["outcome"] == "not testable"
    assert t["by_seed"][0]["outcome"] == "confirmed"
    assert t["outcome_of_the_round"] == "unresolved"              # GATE5.md, part 1: the seeds do not all read as the three together


def test_tests_use_a_run_of_their_own():
    fam, _ = _family({"coverage": [0.06] * 3, "function": [0.0] * 3, "amount": [0.0] * 3})
    t = read_final.tests_of(HYP, fam, {k: [True] * 3 for k in ("coverage", "function", "amount")}, ("coverage", "function", "amount"))
    assert t["together"]["runs"] == [read_final.RUN]
    assert not [c for c in loop.ledger_check([HYP, t["together"]] + t["by_seed"]) if c["severity"] == "error"]


def test_difference_between_two_groups_of_donors():
    g = np.random.default_rng(1)
    n = 40
    first, second = np.arange(20), np.arange(20, 40)
    Ub = 0.6 + g.normal(0, 0.02, (n, 1)) + g.normal(0, 0.003, (n, 3))
    same = Ub - 0.05 + g.normal(0, 0.003, (n, 3))
    d = read_final.group_difference(Ub, same, first, second)
    assert d["ci"][0] < 0 < d["ci"][1] and abs(d["gain"]) < 0.005
    less = same.copy()
    less[first] += 0.03                                           # the first group gains 0.02 where the second gains 0.05
    d = read_final.group_difference(Ub, less, first, second)
    assert d["ci"][1] < 0 and abs(d["gain"] + 0.03) < 0.005
    a, b = read_final.group_gain(Ub, less, first), read_final.group_gain(Ub, less, second)
    assert abs(a["gain"] - 0.02) < 0.005 and abs(b["gain"] - 0.05) < 0.005 and a["donors"] == 20
    assert abs((a["gain"] - b["gain"]) - d["gain"]) < 1e-9


def test_files_that_are_not_those_of_the_declaration_are_refused(tmp_path):
    f1, f2 = tmp_path / "a.npz", tmp_path / "b.npy"
    f1.write_bytes(b"pack"); f2.write_bytes(b"states")
    decl = {"sha256": {"final_arrays": read_final.sha256(f1), "final_z": read_final.sha256(f2)}}
    read_final.check_declaration(decl, {"final_arrays": f1, "final_z": f2})
    f2.write_bytes(b"other states")
    with pytest.raises(SystemExit, match="final_z"):
        read_final.check_declaration(decl, {"final_arrays": f1, "final_z": f2})
    with pytest.raises(SystemExit, match="final_arrays"):
        read_final.check_declaration(decl, {"final_z": f2})


def test_the_arms_of_the_test_are_heads_of_the_plan():
    names = {c["name"] for c in plan.heads()} | {c["name"] for c in plan.heads(final=True)}
    assert {f"s{s}_{read_final.HEAD_OF[k]}" for s in plan.SEEDS for k in read_final.ARMS} <= names
    assert {c["arm"] for c in plan.heads(final=True)} == {"coverage", "labels"}
