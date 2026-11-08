"""The ordered root-cause rules on hand-made verdicts."""
from anm_diag.verdicts import root_cause, verdict

CH = ["RAW", "TOKENS", "Z"]


def test_three_way_verdict():
    assert verdict([-0.1, 0.05], 0.10) == "support"
    assert verdict([0.2, 0.5], 0.10) == "reject"
    assert verdict([0.05, 0.3], 0.10) == "unresolved"


def test_union_of_readout_classes():
    from anm_diag.verdicts import verdict_union
    # ridge recovers the readout exactly; boosting does not: the state is supported (some class reproduces it)
    v, c = verdict_union({"ridge": {"u": 0.0, "ci": [-0.001, 0.001]}, "hgb": {"u": -0.01, "ci": [-0.03, 0.15]}}, 0.10)
    assert (v, c) == ("support", "ridge")
    # rejected only when every class fails
    v, _ = verdict_union({"ridge": {"u": 0.47, "ci": [0.29, 0.64]}, "hgb": {"u": 0.78, "ci": [0.64, 0.86]}}, 0.10)
    assert v == "reject"
    v, _ = verdict_union({"ridge": {"u": 0.47, "ci": [0.29, 0.64]}, "hgb": {"u": 0.2, "ci": [0.05, 0.35]}}, 0.10)
    assert v == "unresolved"


def test_order_of_rules():
    sup = {c: "support" for c in CH}
    assert root_cause([-0.01, 0.02], "support", CH, sup)["cause"] == "no_effect"
    assert root_cause([0.1, 0.2], "unresolved", CH, sup)["cause"] == "control_failed"
    assert root_cause([0.1, 0.2], "support", CH, sup)["cause"] == "no_deficit"
    v = {"RAW": "support", "TOKENS": "support", "Z": "reject"}
    rc = root_cause([0.1, 0.2], "support", CH, v)
    assert rc["cause"] == "loss_at" and (rc["from"], rc["to"]) == ("TOKENS", "Z")
    v = {"RAW": "support", "TOKENS": "unresolved", "Z": "reject"}
    rc = root_cause([0.1, 0.2], "support", CH, v)
    assert rc["cause"] == "loss_at" and rc["span"] == ["RAW", "TOKENS", "Z"]
    v = {"RAW": "support", "TOKENS": "unresolved", "Z": "unresolved"}
    assert root_cause([0.1, 0.2], "support", CH, v)["cause"] == "unresolved"


def test_reference_unresolved_and_relative_loss():
    v = {"RAW": "unresolved", "TOKENS": "unresolved", "Z": "reject"}
    steps = {"RAW->TOKENS": [-0.03, 0.30], "TOKENS->Z": [-0.04, 0.22]}
    assert root_cause([-0.1, -0.05], "support", CH, v, steps, 0.10)["cause"] == "reference_unresolved"
    steps = {"RAW->TOKENS": [-0.03, 0.30], "TOKENS->Z": [0.15, 0.40]}
    rc = root_cause([-0.1, -0.05], "support", CH, v, steps, 0.10)
    assert rc["cause"] == "relative_loss_at" and (rc["from"], rc["to"]) == ("TOKENS", "Z")


def test_paths_expand_home(tmp_path, monkeypatch):
    from anm_diag.adapters.base import load_paths
    monkeypatch.setenv("HOME", str(tmp_path))
    p = load_paths({"a": "~/x.npz", "b": ["~/y", 3], "c": 1})
    assert p["a"] == f"{tmp_path}/x.npz" and p["b"][0] == f"{tmp_path}/y" and p["c"] == 1
