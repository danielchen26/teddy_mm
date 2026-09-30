from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from safetensors.torch import load_file


class TeddyGEncoder(nn.Module):
    def __init__(self, cfg: dict):
        super().__init__()
        self.d_model = int(cfg["d_model"])
        self.ntoken = int(cfg["ntoken"])
        layer = nn.TransformerEncoderLayer(
            d_model=self.d_model,
            nhead=int(cfg["nheads"]),
            dim_feedforward=int(cfg["d_hid"]),
            dropout=float(cfg.get("dropout", 0.0)),
            batch_first=True,
            norm_first=bool(cfg.get("pre_norm", False)),
            activation=cfg.get("layer_activation", "gelu"),
        )
        self.embeddings = nn.Embedding(self.ntoken, self.d_model)
        self.position_embeddings = nn.Embedding(int(cfg["max_position_embeddings"]), self.d_model)
        self.encoder = nn.TransformerEncoder(layer, num_layers=int(cfg["nlayers"]))
        self.decoder_head = nn.Linear(self.d_model, self.ntoken, bias=False)
        self.decoder_bias = nn.Parameter(torch.zeros(self.ntoken))

    def forward(self, gene_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        """Legacy teddy_mm embedding (kept unchanged): nn.TransformerEncoder call, then the mean
        over real tokens. On MPS this call needs PYTORCH_ENABLE_MPS_FALLBACK=1."""
        pos = torch.arange(gene_ids.size(1), device=gene_ids.device)
        hidden = self.embeddings(gene_ids) + self.position_embeddings(pos)
        hidden = self.encoder(hidden, src_key_padding_mask=~attention_mask.bool())
        mask = attention_mask.unsqueeze(-1).to(hidden.dtype)
        return (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1.0)

    def hidden_states(
        self,
        gene_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        return_layer_means: bool = False,
    ):
        """Official TEDDY-G forward without annotations (model.py:222-258): token + position
        embeddings, then each encoder layer in turn. Returns the last-layer states [B, L, d]
        (padding positions keep their computed values, as in model.py) and, optionally, the
        gene-mean of every layer's output [n_layers, B, d].

        Padding mask: model.py builds its key-padding mask as
        torch.cat([torch.empty(0), ~attention_mask.bool()], dim=1), which type-promotes to a
        float {0, 1} tensor. On the CPU/CUDA eval fast path that float mask is converted back to
        bool, so padding is masked (the official tutorial runs on CPU in eval/no_grad). On the
        slow path (MPS, or grad enabled) a float mask is *added* to the attention logits and
        padding is attended. We pass the bool mask so every device reproduces the tutorial's
        CPU result (checked: CPU fast path float mask == bool mask, max |diff| 0.0).
        """
        pos = torch.arange(gene_ids.size(1), device=gene_ids.device)
        hidden = self.embeddings(gene_ids) + self.position_embeddings(pos)
        pad_mask = ~attention_mask.bool()
        layer_means = []
        for layer in self.encoder.layers:
            hidden = layer(src=hidden, src_key_padding_mask=pad_mask)
            if return_layer_means:
                layer_means.append(pool_hidden(hidden, attention_mask, "gene-mean"))
        if return_layer_means:
            return hidden, torch.stack(layer_means, dim=0)
        return hidden


POOLINGS = ("gene-mean", "all-positions", "first")


def pool_hidden(hidden: torch.Tensor, attention_mask: torch.Tensor, how: str) -> torch.Tensor:
    """Cell embedding from last-layer states.

    gene-mean: mean over real gene tokens (paper §5.4, "averaging over output gene embeddings");
    all-positions: hidden.mean(1) over every position, padding included (official tutorial,
        cell 19; only meaningful when every cell is padded to the same length, 2048);
    first: hidden[:, 0] (model.py return_cell_embs_first_token; with no CLS or annotation
        token this is the cell's top-ranked gene).
    """
    if how == "gene-mean":
        mask = attention_mask.unsqueeze(-1).to(hidden.dtype)
        return (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1.0)
    if how == "all-positions":
        return hidden.mean(dim=1)
    if how == "first":
        return hidden[:, 0, :]
    raise ValueError(f"unknown pooling {how!r}; choose from {POOLINGS}")


def load_vocab(ckpt: Path) -> dict[str, int]:
    vocab = {}
    with (ckpt / "vocab.txt").open() as f:
        for i, line in enumerate(f):
            vocab[line.rstrip("\n")] = i
    return vocab


def load_pad_id(ckpt: Path, vocab: dict[str, int]) -> int:
    added = json.loads((ckpt / "added_tokens.json").read_text())
    return int(added.get("<pad>", vocab.get("<pad>", 43810)))


def load_teddy(ckpt: Path, device: torch.device) -> TeddyGEncoder:
    cfg = json.loads((ckpt / "config.json").read_text())
    weights = ckpt / "model.safetensors"
    if not weights.exists() or weights.stat().st_size < 1_000_000:
        raise FileNotFoundError(f"Need real TEDDY weights at {weights}. Run teddy_mwe/setup.sh.")
    model = TeddyGEncoder(cfg)
    state = load_file(str(weights), device="cpu")
    missing, _ = model.load_state_dict(state, strict=False)
    if missing:
        raise RuntimeError(f"missing TEDDY weights: {missing}")
    model.to(device)
    model.eval()
    return model


def load_gene_medians(path: Path) -> dict[str, float]:
    """TEDDY's gene-median file (teddy/data_processing/utils/medians/data/teddy_gene_medians.json):
    Ensembl id -> non-zero median of the gene after scaling each cell to 10,000 counts."""
    return {str(k): float(v) for k, v in json.loads(Path(path).read_text()).items()}


def median_factors(rna_names, medians: dict[str, float]) -> np.ndarray:
    """Per-gene factor 1/median, and 1.0 for a gene with no median
    (preprocess.py scale_columns_by_median_dict, :297-302)."""
    return np.array([1.0 / medians[g] if g in medians else 1.0 for g in map(str, rna_names)], dtype=np.float64)


def official_values(block: np.ndarray, factors: np.ndarray, normalized_total: float = 1e4) -> np.ndarray:
    """Official TEDDY-G values for a dense block of RNA rows [n_cells, n_mapped_genes]:
    L1 row normalisation times normalized_total (preprocess.py:272-281), then each gene divided
    by its median (:284-306); no log1p. Each step is computed in float64 and stored as float32,
    as the in-place sklearn/numpy operations on a float32 matrix do; torch.topk then reads float32.
    Only the order of the values reaches the model, so the per-cell total cancels (the stored
    matrix is restricted to vocab-mapped genes; the official total runs over all genes)."""
    x = np.asarray(block, dtype=np.float32)
    tot = x.sum(axis=1, keepdims=True, dtype=np.float64)
    tot[tot == 0] = 1.0
    x = (x.astype(np.float64) / tot).astype(np.float32)
    x = (x.astype(np.float64) * float(normalized_total)).astype(np.float32)
    return (x.astype(np.float64) * np.asarray(factors, dtype=np.float64)[None, :]).astype(np.float32)


def rank_encode_official(
    values: np.ndarray,
    token_ids: np.ndarray,
    max_len: int = 2048,
    pad_id: int = 43810,
    pad_to: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Official TEDDY-G tokenisation of a block [n_cells, n_mapped_genes] of float32 values.

    torch.topk(X, k=min(max_len, n_genes), largest=True, sorted=True) on CPU
    (tokenization.py:168; add_cls False so all max_len slots are genes), then zero-valued entries
    are dropped (include_zero_genes False, :261-264), then each cell is padded with pad_id to
    pad_to (default max_len, the tutorial collate_fn) with attention_mask = ids != pad.
    Tie order is whatever torch.topk returns, as in the official code.
    """
    X = torch.from_numpy(np.ascontiguousarray(values, dtype=np.float32))
    k = min(int(max_len), X.shape[1])
    top_vals, top_idx = torch.topk(X, k=k, largest=True, sorted=True)
    tok_arr = torch.as_tensor(np.asarray(token_ids, dtype=np.int64))
    gene_ids = tok_arr[top_idx]
    n = X.shape[0]
    width = int(pad_to if pad_to is not None else max_len)
    tokens = np.full((n, width), pad_id, dtype=np.int64)
    attn = np.zeros((n, width), dtype=np.int64)
    nz = (top_vals != 0).numpy()
    g = gene_ids.numpy()
    for i in range(n):
        row = g[i][nz[i]]
        if row.size > width:
            raise ValueError(f"cell {i} has {row.size} tokens > pad_to {width}")
        tokens[i, : row.size] = row
        attn[i, : row.size] = 1
    return tokens, attn


def rank_encode_matrix(
    expr: np.ndarray,
    token_ids: np.ndarray,
    seq_len: int,
    pad_id: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Legacy teddy_mm tokenisation (kept unchanged): stable argsort, stop at the first zero,
    truncate at seq_len; ties keep column (genomic) order.
    expr: [n_cells, n_mapped_genes], token_ids: [n_mapped_genes]."""
    n = expr.shape[0]
    tokens = np.full((n, seq_len), pad_id, dtype=np.int64)
    attn = np.zeros((n, seq_len), dtype=np.int64)
    for i in range(n):
        row = expr[i]
        order = np.argsort(-row, kind="stable")
        k = 0
        for j in order:
            if row[j] <= 0:
                break
            tokens[i, k] = int(token_ids[j])
            attn[i, k] = 1
            k += 1
            if k >= seq_len:
                break
    return tokens, attn
