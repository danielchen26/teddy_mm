"""Replay of the invalid designs met so far: the checker must refuse each, and accept the valid ones."""
import copy
from pathlib import Path

import pytest

from anm_diag import check, load_declaration
from anm_diag.checks import errors

HERE = Path(__file__).resolve().parent.parent
BASE = load_declaration(HERE / "declarations" / "ks2_teddy_mitoribo.yaml")
CITE = load_declaration(HERE / "declarations" / "teddy_mm_frangieh2021.yaml")


def codes(decl):
    return {p["code"] for p in errors(check(decl))}


def mutate(base, fn):
    d = copy.deepcopy(base)
    fn(d)
    return d


def test_valid_declarations_pass():
    assert codes(BASE) == set()
    assert codes(CITE) == set()


def test_readout_from_model_output_rejected():
    # the 65 earlier rounds: observation = the model's own prediction
    d = mutate(BASE, lambda d: d["readouts"].append({"id": "PRED", "measurement": "Z", "statistic": "x", "features": "x"}))
    assert "R1" in codes(d)


def test_approximate_tokens_on_path_rejected():
    def f(d):
        for n in d["nodes"]:
            if n["id"] == "TOKENS":
                n["operator"] = "approximation"
    assert "C2" in codes(mutate(BASE, f))


def test_undeclared_pca_rejected():
    def f(d):
        d["nodes"].append({"id": "Z_PCA50", "from": "Z", "op": "pca:50", "operator": "analysis", "role": "probe"})
    assert "C3" in codes(mutate(BASE, f))


def test_masked_input_rejected():
    def f(d):
        d["nodes"].append({"id": "Z_MASKED", "from": "TOKENS", "op": "mask_input:readout_genes", "operator": "system",
                           "role": "probe"})
    assert "C4" in codes(mutate(BASE, f))


def test_state_patch_rejected():
    # teddy_mm E7: mean-preserving patches at a layer cut are system interventions, not histories
    def f(d):
        d["nodes"].append({"id": "L11_PATCH", "from": "TOKENS", "op": "mean_preserving_patch:11", "operator": "system",
                           "role": "probe"})
    assert "C4" in codes(mutate(BASE, f))


def test_kmeans_partition_rejected():
    assert "E1" in codes(mutate(BASE, lambda d: d["inference"].__setitem__("inverse", "kmeans_partition")))


def test_off_path_stage_rejected():
    # comparing a stage that is not on the computation path (e.g. a re-implemented tokeniser branch)
    def f(d):
        d["nodes"].append({"id": "TOK2", "from": "RAW", "op": "teddy_official_tokens", "operator": "system", "role": "chain"})
        d["chain"] = ["RAW", "TOKENS", "TOK2", "Z"]
    assert "C1" in codes(mutate(BASE, f))


def test_no_controlled_source_rejected():
    # teddy_mm Mode B: a decision layer with no source perturbation
    d = mutate(BASE, lambda d: d["source"].__setitem__("controlled", False))
    assert "S1" in codes(d)


def test_decision_layer_readout_rejected():
    # teddy_mm Mode B: the "readout" is the engine's own call computed from predicted proteins
    def f(d):
        d["nodes"].append({"id": "ANM_CALL", "from": "Z", "op": "decision_engine", "operator": "analysis", "role": "probe"})
        d["readouts"] = [{"id": "CALL", "measurement": "ANM_CALL", "statistic": "label", "features": "x"}]
    assert "R1" in codes(mutate(BASE, f))


def test_shared_measurement_needs_declaration():
    d = mutate(BASE, lambda d: d.__setitem__("measurement_sharing", "none"))
    assert "R2" in codes(d)


def test_protein_input_leak_rejected():
    # teddy_mm phase-2: a model that takes measured protein as input, tested against measured protein
    def f(d):
        d["nodes"].append({"id": "PH2", "from": ["Z", "adt_counts"], "op": "phase2_fusion", "operator": "system",
                           "role": "probe"})
    assert "R2" in codes(mutate(CITE, f))


def test_sealed_needs_hashes():
    assert "K2" in codes(mutate(BASE, lambda d: d.__setitem__("status", "sealed")))
