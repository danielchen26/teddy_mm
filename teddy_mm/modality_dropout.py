"""Modality dropout + latent fusion for bidirectional CITE (phase-2).

Modes (mutually exclusive per cell, RNA dropped ~15-20% by default):
  - rna_only:  RNA → ADT path (phase-1 style); z_obs = z_rna
  - adt_only:  ADT → latent / RNA-side; z_obs = z_adt
  - joint:     both present; z_obs = fuse(z_rna, z_adt)
Never drops both modalities.
"""
from __future__ import annotations

import torch
import torch.nn.functional as F


def clr_log_adt(adt: torch.Tensor, eps: float = 1e-6) -> torch.Tensor:
    """CLR on log1p counts: log(x) - mean(log(x)) with x = adt + eps."""
    x = adt.clamp(min=0) + eps
    log_x = torch.log(x)
    return log_x - log_x.mean(dim=-1, keepdim=True)


def sample_modality_masks(
    batch_size: int,
    *,
    p_drop_rna: float = 0.18,
    p_drop_adt: float = 0.15,
    device: torch.device | None = None,
    generator: torch.Generator | None = None,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Bernoulli masks; force at least one modality kept.

    Returns
    -------
    keep_rna, keep_adt : BoolTensor [B]
    """
    device = device or torch.device("cpu")
    u_rna = torch.rand(batch_size, device=device, generator=generator)
    u_adt = torch.rand(batch_size, device=device, generator=generator)
    keep_rna = u_rna >= p_drop_rna
    keep_adt = u_adt >= p_drop_adt
    both_drop = ~keep_rna & ~keep_adt
    # Prefer keeping RNA when both would drop (rarer ADT-only still sampled above).
    keep_rna = keep_rna | both_drop
    return keep_rna, keep_adt


def fuse_latents(
    z_rna: torch.Tensor,
    z_adt: torch.Tensor,
    keep_rna: torch.Tensor,
    keep_adt: torch.Tensor,
    mode: str = "mean",
) -> torch.Tensor:
    """Build z_obs from available modality latents.

    Parameters
    ----------
    z_rna, z_adt : [B, D]
    keep_rna, keep_adt : [B] bool
    mode : 'mean' | 'sum' | 'rna_priority'
    """
    kr = keep_rna.to(dtype=z_rna.dtype).unsqueeze(-1)
    ka = keep_adt.to(dtype=z_adt.dtype).unsqueeze(-1)
    if mode == "sum":
        return kr * z_rna + ka * z_adt
    if mode == "rna_priority":
        # Prefer RNA when present; else ADT.
        return torch.where(keep_rna.unsqueeze(-1), z_rna, z_adt)
    # mean over present modalities
    num = kr * z_rna + ka * z_adt
    den = (kr + ka).clamp(min=1e-6)
    return num / den


def masked_l2(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Mean L2 over rows where mask is True; 0 if none."""
    if mask.dtype != torch.bool:
        mask = mask.bool()
    if not mask.any():
        return pred.new_zeros(())
    diff = pred[mask] - target[mask]
    return (diff ** 2).mean()
