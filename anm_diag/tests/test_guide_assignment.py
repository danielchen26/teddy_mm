"""Cell selection from detected guides and list-valued strata (Perturb-CITE adapter), on a hand-made obs table."""
import numpy as np
import pandas as pd

from anm_diag.adapters.perturb_cite import PerturbCiteAdapter


def _adapter(obs, **paths):
    a = object.__new__(PerturbCiteAdapter)
    a.obs, a.paths = obs, paths
    return a


def test_frangieh_style_guides():
    obs = pd.DataFrame({
        "guide_id": ["IFNGR1_1", "IFNGR1_3;KCNN4_2", "IFNGR1_2;NO_SITE_664", "NO_SITE_47", "NO_SITE_1;ONE_NON-GENE_SITE_2",
                     np.nan, "HLA-A_3", "ONE_NON-GENE_SITE_531"],
        "perturbation": ["IFNGR1", "IFNGR1", "IFNGR1", "control", "control", "control", "HLA-A", "control"],
        "perturbation_2": ["IFNγ", "IFNγ", "IFNγ", "IFNγ", "Control", "IFNγ", "IFNγ", "IFNγ"]})
    obs["guide_id"] = obs["guide_id"].astype("category")
    a = _adapter(obs, obs_guide="guide_id", guide_gene_regex=r"^(.+)_\d+$", obs_stratum="perturbation_2")
    src = {"assign_from_guides": True, "control_guides": ["NO_SITE", "ONE_NON-GENE_SITE"], "control": ["control"],
           "treated": ["IFNGR1", "HLA-A"], "stratum": "IFNγ"}
    assert a._targets(src).tolist() == ["IFNGR1", "", "IFNGR1", "control", "control", "", "HLA-A", "control"]
    sel, isX, tgt = a._select({"source": src})
    # two-gene cell, no-guide cell and the other condition are left out
    assert sel.tolist() == [0, 2, 3, 6, 7]
    assert isX[sel].tolist() == [True, True, False, True, False]


def test_papalexi_style_guides_and_stratum_list():
    obs = pd.DataFrame({"guide_id": ["IFNGR1g1", "IFNGR1g4", "NTg7", "STAT1g2", "NTg1", "JAK2g3"],
                        "hto": ["rep1-tx", "rep3-tx", "rep4-tx", "rep2-ctrl", "rep2-tx", "rep1-tx"]})
    a = _adapter(obs, obs_guide="guide_id", guide_gene_regex=r"^(.+?)g\d+$", obs_stratum="hto")
    src = {"assign_from_guides": True, "control_guides": ["NT"], "control": ["control"], "treated": ["IFNGR1", "STAT1"],
           "stratum": ["rep1-tx", "rep2-tx", "rep3-tx", "rep4-tx"]}
    sel, isX, tgt = a._select({"source": src})
    assert tgt.tolist() == ["IFNGR1", "IFNGR1", "control", "STAT1", "control", "JAK2"]
    assert sel.tolist() == [0, 1, 2, 4]        # unstimulated cell and the untreated target are left out
    assert isX[sel].tolist() == [True, True, False, False]


def test_obs_target_unchanged_without_guides():
    obs = pd.DataFrame({"perturbation": ["A", "control", "B"], "cond": ["x", "x", "y"]})
    a = _adapter(obs, obs_target="perturbation", obs_stratum="cond")
    sel, isX, _ = a._select({"source": {"control": ["control"], "treated": ["A", "B"], "stratum": "x"}})
    assert sel.tolist() == [0, 1] and isX[sel].tolist() == [True, False]
