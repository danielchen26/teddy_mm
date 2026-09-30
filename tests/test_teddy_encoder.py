"""Tests for the official TEDDY-G preprocessing in teddy_mm/teddy_encoder.py and the legacy path."""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT))

from teddy_mm.teddy_encoder import (  # noqa: E402
    load_teddy,
    median_factors,
    official_values,
    pool_hidden,
    rank_encode_matrix,
    rank_encode_official,
)

PAD = 43810


def _main_checkout() -> Path:
    """The main checkout (data/ and outputs/ live there when tests run from a git worktree)."""
    dotgit = ROOT / ".git"
    if dotgit.is_file():  # worktree: "gitdir: <main>/.git/worktrees/<name>"
        gitdir = dotgit.read_text().split("gitdir:", 1)[-1].strip()
        if "/.git/" in gitdir:
            return Path(gitdir.split("/.git/", 1)[0])
    return ROOT


CKPT = Path(os.environ.get("TEDDY_CKPT", str(_main_checkout().parent / "teddy_mwe/ckpt/teddy_g_70M")))
PROCESSED = Path(os.environ.get("TEDDY_PROCESSED", str(_main_checkout() / "data/processed/cite")))
HAVE_CKPT = (CKPT / "model.safetensors").exists() and (CKPT / "model.safetensors").stat().st_size > 1_000_000
HAVE_DATA = (PROCESSED / "cite_arrays.npz").exists() and (PROCESSED / "z_rna.npy").exists()


# ------------------------------------------------------------------ verbatim official snippets
def _official_tokenize_batch(X_batch: torch.Tensor, token_array: torch.Tensor, max_seq_len: int):
    """tokenization.py _build_batch_tensors (random_genes False, add_cls False) + the
    include_zero_genes False filter of tokenize(), copied verbatim."""
    seq_tokens = max_seq_len
    top_vals, top_indices = torch.topk(X_batch, k=min(seq_tokens, X_batch.shape[1]), largest=True, sorted=True)
    gene_ids = token_array[top_indices]
    gene_ids_batch, vals_batch = gene_ids, top_vals
    final_gene_list = []
    for row_idx in range(len(gene_ids_batch)):
        g_row = gene_ids_batch[row_idx]
        v_row = vals_batch[row_idx]
        nonzero_mask = v_row != 0
        g_row = g_row[nonzero_mask]
        v_row = v_row[nonzero_mask]
        final_gene_list.append(g_row)
    return [g.tolist() for g in final_gene_list]


def _official_collate(batch, pad_token_id: int, max_seq_len: int = 2048):
    """tutorial cell 13 collate_fn, copied verbatim (tokenizer.pad_token_id -> pad_token_id)."""
    batch_size = len(batch)
    max_len = max_seq_len
    input_ids = torch.full((batch_size, max_len), pad_token_id, dtype=torch.long)
    for i, sample in enumerate(batch):
        seq = sample["gene_ids"]
        input_ids[i, :len(seq)] = torch.tensor(seq, dtype=torch.long)
    attention_mask = (input_ids != pad_token_id).long()
    return {"gene_ids": input_ids, "attention_mask": attention_mask}


def _synthetic_rows(n=40, g=300, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.6, size=(n, g)).astype(np.float32)  # many zeros and many ties
    X[0] = 0.0                      # empty cell
    X[1] = 0.0; X[1, :5] = 3.0        # 5 genes, all tied
    X[2] = rng.integers(1, 3, g)    # every gene non-zero, heavy ties -> truncation inside a tie
    return X


# ------------------------------------------------------------------ (a) tokenisation
@pytest.mark.parametrize("max_len", [16, 64, 512])
def test_rank_encode_official_matches_verbatim_official_code(max_len):
    X = _synthetic_rows()
    fac = np.random.default_rng(1).choice([0.5, 1.0, 2.0], size=X.shape[1])  # keeps some ties
    vals = official_values(X, fac)
    token_ids = np.arange(1000, 1000 + X.shape[1], dtype=np.int64)
    ref_lists = _official_tokenize_batch(torch.tensor(vals, dtype=torch.float), torch.tensor(token_ids), max_len)
    ref = _official_collate([{"gene_ids": s} for s in ref_lists], PAD, max_seq_len=max_len)
    # full batch and awkward sub-batches must give the same tokens (tie order independent of the batch)
    for bs in (len(X), 7, 1):
        toks, attn = [], []
        for s in range(0, len(X), bs):
            t, a = rank_encode_official(vals[s:s + bs], token_ids, max_len=max_len, pad_id=PAD)
            toks.append(t); attn.append(a)
        toks, attn = np.concatenate(toks), np.concatenate(attn)
        assert np.array_equal(toks, ref["gene_ids"].numpy()), bs
        assert np.array_equal(attn, ref["attention_mask"].numpy()), bs
    assert attn[0].sum() == 0 and attn[1].sum() == 5 and attn[2].sum() == min(max_len, X.shape[1])


def test_rank_encode_official_pad_to_shorter_width_and_overflow():
    X = _synthetic_rows()
    vals = official_values(X, np.ones(X.shape[1]))
    tid = np.arange(X.shape[1])
    t_full, a_full = rank_encode_official(vals, tid, max_len=64, pad_id=PAD)
    t_wide, a_wide = rank_encode_official(vals, tid, max_len=64, pad_id=PAD, pad_to=80)
    assert np.array_equal(t_wide[:, :64], t_full) and np.all(t_wide[:, 64:] == PAD) and a_wide[:, 64:].sum() == 0
    with pytest.raises(ValueError):
        rank_encode_official(vals, tid, max_len=64, pad_id=PAD, pad_to=10)


def test_official_values_match_official_preprocess_ops():
    sp = pytest.importorskip("scipy.sparse")
    sf = pytest.importorskip("sklearn.utils.sparsefuncs")
    sff = pytest.importorskip("sklearn.utils.sparsefuncs_fast")
    X = _synthetic_rows()[1:]  # drop the empty cell (official L1 leaves it at 0, as do we)
    names = [f"G{i}" for i in range(X.shape[1])]
    med = {g: float(v) for g, v in zip(names[::2], np.random.default_rng(2).uniform(0.3, 3.0, len(names[::2])))}
    fac = median_factors(names, med)
    # preprocess.py normalize_data_inplace + scale_columns_by_median_dict, verbatim ops
    m = sp.csr_matrix(X.astype(np.float32))
    sff.inplace_csr_row_normalize_l1(m)
    sf.inplace_row_scale(m, np.array([10000.0] * m.shape[0]))
    factors = np.array([1.0 / med[g] if g in med else 1.0 for g in names])
    sf.inplace_csr_column_scale(m, factors)
    ref = m.toarray().astype(np.float32)
    ours = official_values(X, fac)
    # only the order reaches the model; values agree to float32 rounding and the order is identical
    assert np.allclose(ours, ref, rtol=1e-6, atol=0)
    tid = np.arange(X.shape[1])
    assert np.array_equal(rank_encode_official(ours, tid, 64, PAD)[0], rank_encode_official(ref, tid, 64, PAD)[0])


# ------------------------------------------------------------------ (b) medians
def test_median_factors_missing_gene_is_one():
    f = median_factors(["A", "B", "C"], {"A": 0.5, "C": 4.0})
    assert f.tolist() == [2.0, 1.0, 0.25]


def test_pool_hidden_definitions():
    h = torch.arange(2 * 4 * 3, dtype=torch.float32).reshape(2, 4, 3)
    m = torch.tensor([[1, 1, 0, 0], [1, 1, 1, 1]])
    assert torch.allclose(pool_hidden(h, m, "gene-mean")[0], h[0, :2].mean(0))
    assert torch.allclose(pool_hidden(h, m, "all-positions")[0], h[0].mean(0))
    assert torch.equal(pool_hidden(h, m, "first"), h[:, 0])
    with pytest.raises(ValueError):
        pool_hidden(h, m, "cls")


def test_legacy_rank_encode_unchanged():
    X = _synthetic_rows()
    t, a = rank_encode_matrix(X, np.arange(X.shape[1]), 16, PAD)
    for i in range(len(X)):
        o = np.argsort(-X[i], kind="stable"); o = o[X[i][o] > 0][:16]
        assert np.array_equal(t[i, :len(o)], o) and a[i].sum() == len(o)


def test_embed_cli_rejects_bad_combinations():
    spec = importlib.util.spec_from_file_location("embed03", ROOT / "scripts/03_embed_rna.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    a = mod.parse_args([])
    assert a.preprocessing == "official" and a.seq_len == 2048 and a.poolings == ["gene-mean"]
    assert mod.parse_args(["--preprocessing", "legacy"]).seq_len == 512
    with pytest.raises(SystemExit):
        mod.parse_args(["--length-buckets", "--save-poolings", "all-positions"])
    with pytest.raises(SystemExit):
        mod.parse_args(["--preprocessing", "legacy", "--save-poolings", "first"])


# ------------------------------------------------------------------ (c) layer loop vs encoder call
@pytest.mark.skipif(not HAVE_CKPT, reason="TEDDY-G weights not found")
def test_layer_loop_matches_encoder_call_gene_mean():
    torch.manual_seed(0)
    model = load_teddy(CKPT, torch.device("cpu"))
    L = 96
    ids = torch.full((3, L), PAD, dtype=torch.long)
    lens = [L, 60, 17]
    for i, n in enumerate(lens):
        ids[i, :n] = torch.randint(0, 43800, (n,))
    mask = (ids != PAD).long()
    with torch.no_grad():
        z_old = model(ids, mask)
        h = model.hidden_states(ids, mask)
        z_new = pool_hidden(h, mask, "gene-mean")
        # padding length does not change real-token states (bool mask): same cell, shorter width
        h2 = model.hidden_states(ids[1:2, :60], mask[1:2, :60])
        # official model.py mask (float {0,1} from torch.cat) on the CPU eval fast path == bool mask
        fmask = torch.cat([torch.empty(0), ~mask.bool()], dim=1)
        hf = model.embeddings(ids) + model.position_embeddings(torch.arange(L))
        for layer in model.encoder.layers:
            hf = layer(src=hf, src_key_padding_mask=fmask)
    cos = torch.nn.functional.cosine_similarity(z_old, z_new, dim=1)
    assert float(cos.min()) > 0.99999
    assert torch.allclose(h2[0], h[1, :60], atol=1e-4)
    assert torch.allclose(pool_hidden(hf, mask, "gene-mean"), z_new, atol=1e-5)


# ------------------------------------------------------------------ (d) legacy reproduces stored z
@pytest.mark.skipif(not (HAVE_CKPT and HAVE_DATA), reason="stored z_rna.npy / weights not found")
def test_legacy_reproduces_stored_rows():
    from teddy_mm.data import load_prepared
    from teddy_mm.teddy_encoder import load_pad_id, load_vocab

    pack = load_prepared(PROCESSED)
    z_ref = np.load(PROCESSED / "z_rna.npy", mmap_mode="r")
    idx = np.sort(np.random.default_rng(0).choice(pack["rna"].shape[0], 20, replace=False))
    model = load_teddy(CKPT, torch.device("cpu"))
    pad = load_pad_id(CKPT, load_vocab(CKPT))
    block = pack["rna"][idx].toarray()
    tot = block.sum(axis=1, keepdims=True); tot[tot == 0] = 1.0
    tok, attn = rank_encode_matrix(block / tot * 10000.0, pack["token_ids"], 512, pad)
    with torch.no_grad():
        z = model(torch.from_numpy(tok), torch.from_numpy(attn)).numpy()
    r = np.asarray(z_ref[idx])
    cos = (z * r).sum(1) / np.linalg.norm(z, axis=1) / np.linalg.norm(r, axis=1)
    assert cos.min() > 0.9999, cos.min()
