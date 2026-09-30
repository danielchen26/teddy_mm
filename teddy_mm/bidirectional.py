"""Phase-2 bidirectional CITE: ADT encoder + modality-dropout latent FM.

No ATAC / perturb-seq / 160M — CITE RNA↔ADT only. cond = z_obs from
present modalities after dropout.

ADT input transform (applied inside AdtEncoder, exactly once):
  The processed pack (teddy_mm.data.prepare_cite) stores the h5ad's ADT ``.X``
  as-is. For data/processed/cite that matrix is already normalized upstream
  (non-integer, non-negative values), so the encoder must not transform it
  again: ``adt_input_transform="none"`` (the default here). ``"clr"`` is the
  legacy behaviour (a second log-ratio on top of the stored values;
  ``clr_log_adt`` maps every stored zero to log(1e-6) minus the row mean).
  ``"log1p"`` is for packs that store raw integer counts.
  ``resolve_adt_input_transform("auto", adt)`` picks "none" for non-integer
  (already normalized) data and "log1p" for integer counts.
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from teddy_mm.modality_dropout import clr_log_adt, fuse_latents, sample_modality_masks
from teddy_mm.models import AdtDecoder, MLP, VelocityNet, integrate_fm, nb_nll

ADT_INPUT_TRANSFORMS = ("none", "clr", "log1p")
# Checkpoints written before the transform was recorded were trained with the
# encoder's old default (use_clr=True), i.e. clr on top of the stored ADT.
LEGACY_ADT_INPUT_TRANSFORM = "clr"


def transform_adt_input(adt: torch.Tensor, transform: str) -> torch.Tensor:
    """The single ADT transform applied on the way into the ADT encoder."""
    if transform == "none":
        return adt
    if transform == "clr":
        return clr_log_adt(adt)
    if transform == "log1p":
        return torch.log1p(adt.clamp(min=0))
    raise ValueError(f"unknown adt_input_transform {transform!r}; choose from {ADT_INPUT_TRANSFORMS}")


def adt_is_count_like(adt) -> bool:
    """True when every finite ADT value is a non-negative integer (raw counts)."""
    a = np.asarray(adt, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return False
    return bool(np.all(a >= 0) and np.all(a == np.round(a)))


def resolve_adt_input_transform(requested: str, adt) -> str:
    """Resolve 'auto' against the stored ADT so the transform happens exactly once.

    auto: 'log1p' when the stored ADT is raw integer counts, else 'none'
    (the pack already holds normalized values).
    """
    if requested == "auto":
        return "log1p" if adt_is_count_like(adt) else "none"
    if requested not in ADT_INPUT_TRANSFORMS:
        raise ValueError(f"unknown adt_input_transform {requested!r}")
    return requested


class AdtEncoder(nn.Module):
    """Shallow MLP on the (once-transformed) ADT → latent z_adt (same dim as TEDDY z_rna)."""

    def __init__(
        self,
        n_adt: int,
        z_dim: int,
        hidden: int = 256,
        n_layers: int = 2,
        dropout: float = 0.1,
        input_transform: str = "none",
    ):
        super().__init__()
        if input_transform not in ADT_INPUT_TRANSFORMS:
            raise ValueError(f"unknown input_transform {input_transform!r}")
        self.input_transform = input_transform
        self.mlp = MLP(n_adt, z_dim, hidden=hidden, n_layers=n_layers, dropout=dropout)

    def forward(self, adt: torch.Tensor) -> torch.Tensor:
        x = transform_adt_input(adt, self.input_transform)
        z = self.mlp(x)
        return F.normalize(z, dim=-1)


class BidirectionalCite(nn.Module):
    """Encode RNA and/or ADT → z_obs; latent FM with cond=z_obs; decode ADT.

    Optional RNA-side head reconstructs z_rna from z_obs when ADT-only / joint
    (ADT → z / RNA side).
    """

    def __init__(
        self,
        z_dim: int,
        n_adt: int,
        *,
        hidden: int = 512,
        adt_hidden: int = 256,
        adt_layers: int = 2,
        n_layers: int = 3,
        dropout: float = 0.1,
        p_drop_rna: float = 0.18,
        p_drop_adt: float = 0.15,
        fuse_mode: str = "mean",
        adt_input_transform: str = "none",
    ):
        super().__init__()
        self.z_dim = z_dim
        self.n_adt = n_adt
        self.p_drop_rna = p_drop_rna
        self.p_drop_adt = p_drop_adt
        self.fuse_mode = fuse_mode
        self.adt_input_transform = adt_input_transform
        self.config = {
            "z_dim": z_dim, "n_adt": n_adt, "hidden": hidden, "adt_hidden": adt_hidden,
            "adt_layers": adt_layers, "n_layers": n_layers, "dropout": dropout,
            "p_drop_rna": p_drop_rna, "p_drop_adt": p_drop_adt, "fuse_mode": fuse_mode,
            "adt_input_transform": adt_input_transform,
        }
        self.adt_enc = AdtEncoder(
            n_adt, z_dim, hidden=adt_hidden, n_layers=adt_layers, dropout=dropout,
            input_transform=adt_input_transform,
        )
        self.rna_proj = MLP(z_dim, z_dim, hidden=hidden, n_layers=2, dropout=dropout)
        self.v_net = VelocityNet(z_dim, cond_dim=z_dim, hidden=hidden, n_layers=n_layers, dropout=dropout)
        self.dec = AdtDecoder(z_dim, n_adt, hidden=hidden)
        # ADT → RNA-side: predict z_rna from z_obs
        self.rna_side = MLP(z_dim, z_dim, hidden=hidden, n_layers=2, dropout=dropout)

    def encode_obs(
        self,
        z_rna: torch.Tensor,
        adt: torch.Tensor,
        keep_rna: torch.Tensor | None = None,
        keep_adt: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Returns z_obs, z_rna_h, z_adt, (keep_rna, keep_adt)."""
        B = z_rna.size(0)
        if keep_rna is None or keep_adt is None:
            keep_rna, keep_adt = sample_modality_masks(
                B, p_drop_rna=self.p_drop_rna, p_drop_adt=self.p_drop_adt, device=z_rna.device
            )
        z_rna_h = F.normalize(self.rna_proj(z_rna), dim=-1)
        z_adt = self.adt_enc(adt)
        z_obs = fuse_latents(z_rna_h, z_adt, keep_rna, keep_adt, mode=self.fuse_mode)
        z_obs = F.normalize(z_obs, dim=-1)
        return z_obs, z_rna_h, z_adt, keep_rna, keep_adt

    def forward_train(
        self,
        z_rna: torch.Tensor,
        adt: torch.Tensor,
        size_factor: torch.Tensor,
        *,
        fm_steps: int = 10,
        lambda_nll: float = 1.0,
        lambda_fm: float = 1.0,
        lambda_rna_side: float = 0.5,
    ) -> dict[str, torch.Tensor]:
        z_obs, z_rna_h, z_adt, keep_rna, keep_adt = self.encode_obs(z_rna, adt)

        # Target latent for FM: joint rebuild toward fused full (both modalities).
        with torch.no_grad():
            z_full = fuse_latents(
                z_rna_h.detach(),
                z_adt.detach(),
                torch.ones_like(keep_rna),
                torch.ones_like(keep_adt),
                mode=self.fuse_mode,
            )
            z_full = F.normalize(z_full, dim=-1)

        # Direct ADT decode from z_obs (RNA→ADT / ADT auto / joint)
        mu_obs, th = self.dec(z_obs, size_factor)
        loss_nll = nb_nll(adt, mu_obs, th)

        # Latent FM: x0 ~ N(0,I) → x1 = z_full, cond = z_obs
        x0 = torch.randn_like(z_obs)
        t = torch.rand(z_obs.size(0), device=z_obs.device)
        x1 = z_full
        xt = (1.0 - t[:, None]) * x0 + t[:, None] * x1
        v_tgt = x1 - x0
        loss_fm = ((self.v_net(xt, t, z_obs) - v_tgt) ** 2).mean()

        z1 = integrate_fm(self.v_net, x0, z_obs, max(fm_steps, 4))
        mu_fm, th_fm = self.dec(z1, size_factor)
        loss_fm_nll = nb_nll(adt, mu_fm, th_fm)

        # ADT→RNA-side: when RNA dropped (or always lightly), predict z_rna_h
        z_rna_hat = F.normalize(self.rna_side(z_obs), dim=-1)
        # Stronger when RNA absent (must recover from ADT); light when joint/rna
        w = torch.where(keep_rna, 0.25, 1.0).to(z_rna_h.dtype)
        rna_err = ((z_rna_hat - z_rna_h.detach()) ** 2).mean(dim=-1)
        loss_rna_side = (rna_err * w).mean()

        loss = (
            lambda_nll * (loss_nll + 0.5 * loss_fm_nll)
            + lambda_fm * loss_fm
            + lambda_rna_side * loss_rna_side
        )
        return {
            "loss": loss,
            "loss_nll": loss_nll.detach(),
            "loss_fm": loss_fm.detach(),
            "loss_fm_nll": loss_fm_nll.detach(),
            "loss_rna_side": loss_rna_side.detach(),
            "frac_rna_kept": keep_rna.float().mean().detach(),
            "frac_adt_kept": keep_adt.float().mean().detach(),
        }

    @torch.no_grad()
    def predict_adt(
        self,
        z_rna: torch.Tensor | None,
        adt: torch.Tensor | None,
        size_factor: torch.Tensor,
        *,
        mode: str = "joint",
        fm_steps: int = 20,
        use_fm: bool = True,
    ) -> torch.Tensor:
        """Inference helper. mode: rna_only | adt_only | joint.

        rna_only never reads ``adt``: the encoder gets zeros and the fused latent
        masks it out, so a measured ADT passed by the caller cannot reach the
        prediction. adt_only / joint take the measured ADT as input by definition.
        ``size_factor`` multiplies the decoded mean; pass a constant (train median)
        at evaluation so the answer key's depth does not enter the prediction.
        """
        B = size_factor.size(0)
        device = size_factor.device
        if mode == "rna_only":
            assert z_rna is not None
            keep_rna = torch.ones(B, dtype=torch.bool, device=device)
            keep_adt = torch.zeros(B, dtype=torch.bool, device=device)
            adt_in = torch.zeros(B, self.n_adt, device=device, dtype=z_rna.dtype)
        elif mode == "adt_only":
            assert adt is not None
            keep_rna = torch.zeros(B, dtype=torch.bool, device=device)
            keep_adt = torch.ones(B, dtype=torch.bool, device=device)
            z_rna = z_rna if z_rna is not None else torch.zeros(B, self.z_dim, device=device)
            adt_in = adt
        else:
            assert z_rna is not None and adt is not None
            keep_rna = torch.ones(B, dtype=torch.bool, device=device)
            keep_adt = torch.ones(B, dtype=torch.bool, device=device)
            adt_in = adt
        z_obs, _, _, _, _ = self.encode_obs(z_rna, adt_in, keep_rna, keep_adt)
        if use_fm:
            z1 = integrate_fm(self.v_net, torch.randn_like(z_obs), z_obs, fm_steps)
            mu, _ = self.dec(z1, size_factor)
        else:
            mu, _ = self.dec(z_obs, size_factor)
        return mu


def load_bidirectional(
    blob: dict,
    device: torch.device | str = "cpu",
    *,
    adt_input_transform: str | None = None,
) -> tuple[BidirectionalCite, dict]:
    """Rebuild a BidirectionalCite from a phase-2 checkpoint dict.

    The ADT input transform comes from the checkpoint's config; checkpoints
    without one (written before the fix) were trained with the legacy 'clr'
    encoder, so they are rebuilt with it. ``adt_input_transform`` overrides.
    Returns (model in eval mode, info dict).
    """
    cfg = dict(blob.get("config") or {})
    recorded = cfg.get("adt_input_transform")
    cfg.setdefault("z_dim", int(blob["z_dim"]))
    cfg.setdefault("n_adt", int(blob["n_adt"]))
    cfg["adt_input_transform"] = adt_input_transform or recorded or LEGACY_ADT_INPUT_TRANSFORM
    z_dim, n_adt = int(cfg.pop("z_dim")), int(cfg.pop("n_adt"))
    model = BidirectionalCite(z_dim, n_adt, **cfg).to(device)
    model.load_state_dict(blob["model"])
    model.eval()
    info = {
        "adt_input_transform": model.adt_input_transform,
        "adt_input_transform_source": (
            "override" if adt_input_transform else ("checkpoint" if recorded else "legacy_default_clr")
        ),
        "train_median_size_factor": blob.get("train_median_size_factor"),
    }
    return model, info
