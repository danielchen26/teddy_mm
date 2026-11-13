"""Smoke test of the Perturb-CITE adapter on a tiny fabricated AnnData pair (code paths only, no biology).

Skipped unless the TEDDY reference files are present (vocab, medians, encoder, phase-1 head). Defaults follow the
teddy_mm layout (TEDDY weights in a sibling teddy_mwe checkout); override with ANM_TEDDY_VOCAB, ANM_TEDDY_MEDIANS,
ANM_TEDDY_ENCODER and ANM_HEAD_CKPT.
"""
import os
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]          # the teddy_mm checkout that holds anm_diag/
NEEDED = [Path(os.environ.get("ANM_TEDDY_VOCAB", REPO.parent / "teddy_mwe/ckpt/teddy_g_70M/vocab.txt")),
          Path(os.environ.get("ANM_TEDDY_MEDIANS", REPO / "data/reference/teddy_gene_medians.json")),
          Path(os.environ.get("ANM_TEDDY_ENCODER", REPO / "teddy_mm/teddy_encoder.py"))]
HEAD = Path(os.environ.get("ANM_HEAD_CKPT", REPO / "outputs/cite_phase1_official/best.pt"))


@pytest.mark.skipif(not all(p.exists() for p in NEEDED) or not HEAD.exists(), reason="TEDDY reference files not found")
def test_perturb_cite_end_to_end(tmp_path):
    import anndata as ad
    import pandas as pd
    from anm_diag import run, load_declaration
    from anm_diag.adapters.perturb_cite import PerturbCiteAdapter

    rng = np.random.default_rng(0)
    vocab = [l.strip() for l, _ in zip(open(NEEDED[0]), range(3000))][10:610]
    genes = pd.DataFrame({"gene_ids": vocab}, index=[f"G{i}" for i in range(len(vocab))])
    genes.index = list(genes.index[:-3]) + ["HLA-A", "B2M", "CD274"]
    targets = ["NT"] * 160 + sum([[t] * 24 for t in ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1"]], [])
    n = len(targets)
    A = np.array([t != "NT" for t in targets])
    lam = rng.gamma(0.6, 2.0, size=(1, len(vocab)))
    X = rng.poisson(lam * rng.uniform(0.5, 1.5, (n, 1))).astype(np.float32)
    X[:, -3:] = rng.poisson(np.where(A[:, None], 3.0, 8.0), (n, 3))
    obs = pd.DataFrame({"target": targets, "condition": "IFNg"}, index=[f"c{i}" for i in range(n)])
    # HLA-E: no response by construction (every target group and the controls hold the same values)
    hlae = np.empty(n)
    for t in set(targets):
        ii = np.where(np.array(targets) == t)[0]
        hlae[ii] = np.resize([20, 27, 33, 40], len(ii))
    adt = np.column_stack([rng.poisson(np.where(A, 20, 60)), rng.poisson(np.where(A, 10, 30)), hlae]).astype(np.float32)
    import scipy.sparse as sp
    ad.AnnData(sp.csr_matrix(X), obs=obs, var=genes).write_h5ad(tmp_path / "rna.h5ad")
    ad.AnnData(adt, obs=obs, var=pd.DataFrame(index=["HLA-A", "CD274", "HLA-E"])).write_h5ad(tmp_path / "adt.h5ad")
    paths = {"rna_h5ad": str(tmp_path / "rna.h5ad"), "adt_h5ad": str(tmp_path / "adt.h5ad"), "obs_target": "target",
             "obs_stratum": "condition", "vocab": str(NEEDED[0]), "medians": str(NEEDED[1]),
             "encoder_dir": str(NEEDED[2].parent), "head_ckpt": str(HEAD), "ensembl_col": "gene_ids"}
    decl = load_declaration(Path(__file__).resolve().parent.parent / "declarations" / "teddy_mm_frangieh2021.yaml")
    decl["inference"].update({"bootstrap": 60, "forward_bootstrap": 30, "folds": 3})
    # the fabricated table uses plain target labels and its own protein names (the real declaration assigns
    # cells from Frangieh's guide column; that selection has its own test in test_guide_assignment.py)
    decl["source"].update({"assign_from_guides": False, "control": ["NT"], "stratum": "IFNg",
                           "treated": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1"]})
    names = {"MHC1": ["HLA-A"], "PDL1": ["CD274"], "HLAE": ["HLA-E"]}
    for r in decl["readouts"]:
        r["features"] = names[r["id"]]
    decl["nodes"] = [n for n in decl["nodes"] if n["id"] != "MRNA_HLAE"]
    for n in decl["nodes"]:
        if n["id"] == "MRNA_HLA":
            n["op"] = "rna_of_gene:HLA-A,B2M"
    # on a 3-protein panel CLR would give HLA-E a compositional response (HLA-A and PD-L1 fall), so the
    # null readout uses log1p counts
    for r in decl["readouts"]:
        if r["id"] == "HLAE":
            r["statistic"] = "adt_log1p"
    adapter = PerturbCiteAdapter(paths, log=lambda *a: None)
    tok = adapter.export_tokens(decl, tmp_path / "tokens.npz")
    t = np.load(tok)
    # stand-in for the Mac TEDDY forward: random states with the right shapes
    m = len(t["L"])
    np.savez(tmp_path / "teddy_states.npz", z_last=rng.normal(size=(m, 512)).astype(np.float32),
             layer_means=rng.normal(size=(12, m, 512)).astype(np.float32), rows=t["rows"])
    adapter.paths["teddy_states"] = str(tmp_path / "teddy_states.npz")
    rec = run(decl, adapter, tmp_path / "out", log=lambda *a: None)
    assert set(rec["readouts"]) == {"MHC1", "PDL1", "HLAE"}
    for rid in ("MHC1", "PDL1"):
        assert rec["readouts"][rid]["states"][f"ADT_{'MHC1' if rid == 'MHC1' else 'PDL1'}"]["verdict"] == "support"
    assert rec["readouts"]["HLAE"]["root_cause"]["cause"] == "no_effect"
    assert rec["readouts"]["HLAE"]["states"] == {}
    # the random stand-in for TEDDY cannot carry the response that the measured RNA carries
    for rid in ("MHC1", "PDL1"):
        assert rec["readouts"][rid]["states"]["Z"]["verdict"] == "reject"
    assert (tmp_path / "out" / "anm_report.md").exists()
