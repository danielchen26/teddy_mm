"""Sealing: a sealed declaration runs only on the data and code it was sealed with, and only unedited."""
import pytest

from anm_diag import run
from anm_diag.adapters.synthetic import SyntheticAdapter, synthetic_declaration
from anm_diag.seal import SealError, seal, verify
from anm_diag.spec import load_declaration


def test_seal_verify_and_refuse(tmp_path):
    f = tmp_path / "vocab.txt"
    f.write_text("a\nb\n")
    paths = {"vocab": str(f)}
    d = seal(synthetic_declaration(), paths, tmp_path / "sealed.yaml")
    assert d["status"] == "sealed" and set(d["data_hashes"]) == {"vocab"}
    d2 = load_declaration(tmp_path / "sealed.yaml")
    assert verify(d2, paths) == ["vocab"]
    ad = SyntheticAdapter(loss_stage="LOSE")
    ad.paths = paths
    rec = run(d2, ad, log=lambda *a: None)
    assert rec["exploratory"] is False and rec["seal_verified"] == ["vocab"]
    assert rec["readouts"]["O"]["root_cause"]["cause"] == "loss_at"
    edited = load_declaration(tmp_path / "sealed.yaml")
    edited["decision"]["margin_u"] = 0.2                 # a design change after sealing
    with pytest.raises(SealError):
        run(edited, ad, log=lambda *a: None)
    f.write_text("a\nc\n")                                # the data changed after sealing
    with pytest.raises(SealError):
        run(d2, ad, log=lambda *a: None)


def test_unsealed_runs_are_exploratory():
    rec = run(synthetic_declaration(), SyntheticAdapter(loss_stage=None), log=lambda *a: None)
    assert rec["exploratory"] is True and rec["seal_verified"] == []
