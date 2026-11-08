"""Ground truth: the runner must localise a known loss, report no deficit and no effect correctly."""
from anm_diag import run
from anm_diag.adapters.synthetic import SyntheticAdapter, synthetic_declaration


def _cause(adapter):
    rec = run(synthetic_declaration(), adapter, log=lambda *a: None)
    return rec["readouts"]["O"]


def test_loss_is_localised():
    r = _cause(SyntheticAdapter(loss_stage="LOSE"))
    assert r["root_cause"]["cause"] == "loss_at"
    assert (r["root_cause"]["from"], r["root_cause"]["to"]) == ("KEEP", "LOSE")
    assert r["states"]["PROGRAM"]["verdict"] == "support"


def test_no_deficit():
    r = _cause(SyntheticAdapter(loss_stage=None))
    assert r["root_cause"]["cause"] == "no_deficit"


def test_no_effect():
    r = _cause(SyntheticAdapter(effect=0.0, loss_stage=None))
    assert r["root_cause"]["cause"] == "no_effect"


def test_loss_at_first_stage_after_reference():
    r = _cause(SyntheticAdapter(loss_stage="KEEP"))
    assert r["root_cause"]["cause"] == "loss_at"
    assert (r["root_cause"]["from"], r["root_cause"]["to"]) == ("INPUT_VIEW", "KEEP")
