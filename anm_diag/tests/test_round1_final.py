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


# ------------------------------------------------------------------ gate 0b: the choice of donors and cells (scripts/final_study.py)
import final_study  # noqa: E402


def _obs(seed=0):
    """A study with three centres: 30 healthy donors and 50 patients, some with two samples, one with too few cells,
    one whose status changes between samples, and a group that is not drawn."""
    g = np.random.default_rng(seed)
    donor, sample, status, centre = [], [], [], []

    def add(d, s, st, c, n):
        donor.extend([d] * n); sample.extend([s] * n); status.extend([st] * n); centre.extend([c] * n)

    for i in range(30):
        add(f"H{i:02d}", f"H{i:02d}_a", "Healthy", ["A", "B", "C"][i % 3], int(g.integers(1100, 1500)))
    for i in range(50):
        c = "A" if i < 35 else ("B" if i < 45 else "C")
        add(f"P{i:02d}", f"P{i:02d}_s1", "Covid", c, int(g.integers(1100, 1500)))
        if i % 7 == 0:
            add(f"P{i:02d}", f"P{i:02d}_s2", "Covid", c, 1200)
    add("SMALL", "SMALL_a", "Covid", "A", 400)
    add("MIXED", "MIXED_a", "Covid", "B", 1200); add("MIXED", "MIXED_b", "Healthy", "B", 1200)
    for i in range(5):
        add(f"L{i}", f"L{i}_a", "LPS", "A", 1200)
    order = g.permutation(len(donor))                              # rows of a sample are not contiguous
    return {"patient_id": np.array(donor)[order], "sample_id": np.array(sample)[order], "Status": np.array(status)[order], "Site": np.array(centre)[order]}


def test_donors_and_cells_are_chosen_by_the_rule():
    obs, cfg = _obs(), final_study.STUDIES["stephenson2021"]
    rows, chosen, cells = final_study.choose(obs, cfg)
    by = {r["donor"]: r for r in rows}
    assert len(chosen["healthy"]) == 20 and len(chosen["patient"]) == 20
    assert not {"SMALL", "MIXED"} & set(chosen["patient"]) and not any(d.startswith("L") for d in chosen["healthy"] + chosen["patient"])
    assert not by["SMALL"]["eligible"] and by["MIXED"]["group"] is None and by["L0"]["group"] is None
    # balanced over centres: taken in turn, so no centre is more than one ahead while all still have donors
    n = {c: sum(by[d]["centres"][0] == c for d in chosen["healthy"]) for c in "ABC"}
    assert max(n.values()) - min(n.values()) <= 1
    n = {c: sum(by[d]["centres"][0] == c for d in chosen["patient"]) for c in "ABC"}
    assert n == {"A": 8, "B": 7, "C": 5}                           # C has only 5 patients; A and B share the rest in turn
    for d, idx in cells.items():
        assert len(idx) == 1000 == len(set(idx.tolist())) and np.all(np.diff(idx) > 0)
        assert set(obs["patient_id"][idx]) == {d} and set(obs["sample_id"][idx]) == {by[d]["first_sample"]}
    rows2, chosen2, cells2 = final_study.choose(obs, cfg)          # the same file gives the same choice
    assert chosen2 == chosen and final_study.rows_digest(cells2) == final_study.rows_digest(cells)
    assert by["P00"]["samples"] == 2 and by["P00"]["first_sample"] == "P00_s1"


def test_the_choice_does_not_depend_on_other_columns():
    obs, cfg = _obs(), final_study.STUDIES["stephenson2021"]
    a = final_study.choose(obs, cfg)
    b = final_study.choose({**obs, "cell_type": np.random.default_rng(5).choice(["T", "B"], len(obs["Status"]))}, cfg)
    assert a[1] == b[1] and final_study.rows_digest(a[2]) == final_study.rows_digest(b[2])


def test_proteins_are_matched_by_name():
    ours = ["CD3", "CD4", "CD11a", "TCR"]
    m = final_study.protein_mapping(ours, ["CD3", "CD4", "CD8", "TCR_ab"], "stephenson2021")
    assert [bool(r["source"]) for r in m] == [True, True, False, False]
    final_study.ALIASES["stephenson2021"]["TCR"] = ["TCR_ab"]
    try:
        m = final_study.protein_mapping(ours, ["CD3", "CD4", "CD8", "TCR_ab"], "stephenson2021")
        assert m[3]["source"] == ["TCR_ab"] and m[3]["rule"].startswith("alias")
        with pytest.raises(SystemExit):
            final_study.protein_mapping(ours, ["CD3", "CD4"], "stephenson2021")
    finally:
        del final_study.ALIASES["stephenson2021"]["TCR"]


def test_the_pack_holds_the_counts_of_the_chosen_cells(tmp_path):
    h5py = pytest.importorskip("h5py")
    import gzip
    from scipy import sparse
    g = np.random.default_rng(3)
    study, cfg = "stephenson2021", final_study.STUDIES["stephenson2021"]
    n, n_gene, n_prot = 5000, 30, 6
    symbols = [f"G{i}" for i in range(n_gene)]
    symbols[7] = symbols[3] = "DUP"                                # a repeated symbol: anndata names the second DUP-1
    in_file = final_study.make_unique(symbols)
    assert in_file[3] == "DUP" and in_file[7] == "DUP-1"
    feats = in_file + ["AB_CD3", "AB_CD4", "AB_CD19", "AB_Isotype", "AB_TCR_ab", "AB_CD8"]
    C = sparse.random(n, n_gene + n_prot, density=0.3, random_state=4, format="csr", data_rvs=lambda k: g.integers(1, 9, k)).astype(np.float32)
    donors = np.repeat(["D1", "D2", "D3", "D4"], n // 4)
    order = g.permutation(n)
    obs = {"patient_id": donors[order], "sample_id": np.char.add(donors[order], "_s1"), "Status": np.where(np.isin(donors[order], ["D1", "D2"]), "Healthy", "Covid"),
           "Site": np.full(n, "A"), "initial_clustering": g.choice(["T", "B"], n)}
    path = tmp_path / "study.h5ad"
    with h5py.File(path, "w") as f:
        f.create_dataset("X", data=np.zeros((2, 2)))
        L = f.create_group("layers").create_group("raw")
        L.attrs["encoding-type"] = "csr_matrix"; L.attrs["shape"] = C.shape
        L.create_dataset("data", data=C.data); L.create_dataset("indices", data=C.indices); L.create_dataset("indptr", data=C.indptr)
        o = f.create_group("obs"); o.attrs["_index"] = "index"
        o.create_dataset("index", data=np.array([f"cell{i}" for i in range(n)], dtype="S"))
        for k, v in obs.items():
            cats, codes = np.unique(v, return_inverse=True)
            c = o.create_group(k); c.create_dataset("categories", data=cats.astype("S")); c.create_dataset("codes", data=codes.astype(np.int8))
        v = f.create_group("var"); v.attrs["_index"] = "index"
        v.create_dataset("index", data=np.array(feats, dtype="S"))
        v.create_dataset("feature_types", data=np.array(["Gene Expression"] * n_gene + ["Antibody Capture"] * n_prot, dtype="S"))
    ref = tmp_path / "features.tsv.gz"
    with gzip.open(ref, "wt") as f:
        for i, sname in enumerate(symbols):
            f.write(f"ENSG{i:04d}\t{sname}\tGene Expression\n")
    our_genes = ["ENSG0003", "ENSG0007", "ENSG9999", "ENSG0000", "ENSG0029"]          # one gene the study does not have
    ours = ["CD3", "CD4", "CD11a", "TCR", "CD8"]
    for d in ("train", "pbmc"):
        (tmp_path / d).mkdir()
    np.savez(tmp_path / "train" / "cite_arrays.npz", adt_names=np.array(ours), rna_names=np.array(our_genes), token_ids=np.arange(5))
    np.savez(tmp_path / "pbmc" / "cite_arrays.npz", adt_present=np.array([True, True, True, True, False]))
    final_study.ALIASES[study]["TCR"] = ["TCR_ab"]
    old = final_study.N_PER_GROUP, final_study.N_CELLS
    final_study.N_PER_GROUP, final_study.N_CELLS = 2, 200
    try:
        sel, cells, _, _ = final_study.select(path, study, tmp_path / "train", tmp_path / "pbmc")
        assert sel["proteins"]["among_the_112"] == ["CD3", "CD4", "TCR"] and sel["proteins"]["of_the_112_not_measured"] == ["CD11a"]
        assert [d["donor"] for d in sel["donors"]] == [f"{study}_D1", f"{study}_D2", f"{study}_D3", f"{study}_D4"]
        (tmp_path / "sel.json").write_text(__import__("json").dumps(sel))
        np.savez(tmp_path / "rows.npz", **cells)
        rc = final_study.main(["pack", str(path), "--selection", str(tmp_path / "sel.json"), "--rows", str(tmp_path / "rows.npz"), "--train-pack", str(tmp_path / "train"),
                               "--reference", str(ref), "--cache", str(tmp_path / "cache"), "--out", str(tmp_path / "pack")])
        assert rc == 0
    finally:
        del final_study.ALIASES[study]["TCR"]
        final_study.N_PER_GROUP, final_study.N_CELLS = old
    z = np.load(tmp_path / "pack" / "cite_arrays.npz")
    rows = z["source_row"]
    assert len(rows) == 800 and list(z["donors"][:200]) == [f"{study}_D1"] * 200 and set(z["status"][:400]) == {"healthy"} and set(z["status"][400:]) == {"patient"}
    assert np.array_equal(np.concatenate([cells[k] for k in ("D1", "D2", "D3", "D4")]), rows)
    X = sparse.csr_matrix((z["rna_data"], z["rna_indices"], z["rna_indptr"]), shape=tuple(z["rna_shape"])).toarray()
    D = C.toarray()
    assert np.array_equal(X[:, 0], D[rows, 3]) and np.array_equal(X[:, 1], D[rows, 7]) and not X[:, 2].any()
    assert np.array_equal(X[:, 3], D[rows, 0]) and np.array_equal(X[:, 4], D[rows, 29])
    P = D[rows, n_gene:]                                                               # CD3 CD4 CD19 Isotype TCR_ab CD8
    assert np.array_equal(z["adt_source_counts"], P) and list(z["adt_present"]) == [True, True, False, True, True]
    gm = np.exp(np.log1p(P).mean(axis=1))
    assert np.allclose(z["adt"][:, 0], np.log1p(P[:, 0] / gm), atol=1e-5) and np.allclose(z["adt"][:, 3], np.log1p(P[:, 4] / gm), atol=1e-5)
    assert np.isnan(z["adt"][:, 2]).all() and np.array_equal(z["n_count_rna"], D[rows, :n_gene].sum(axis=1))
    assert list(z["barcodes"][:3]) == [f"cell{r}" for r in rows[:3]] and set(z["split"]) == {"test"}
