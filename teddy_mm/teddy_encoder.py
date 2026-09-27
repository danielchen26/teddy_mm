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
        pos = torch.arange(gene_ids.size(1), device=gene_ids.device)
        hidden = self.embeddings(gene_ids) + self.position_embeddings(pos)
        hidden = self.encoder(hidden, src_key_padding_mask=~attention_mask.bool())
        mask = attention_mask.unsqueeze(-1).to(hidden.dtype)
        return (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1.0)


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


def rank_encode_matrix(
    expr: np.ndarray,
    token_ids: np.ndarray,
    seq_len: int,
    pad_id: int,
) -> tuple[np.ndarray, np.ndarray]:
    """expr: [n_cells, n_mapped_genes], token_ids: [n_mapped_genes]."""
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
