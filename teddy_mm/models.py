from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class MLP(nn.Module):
    def __init__(self, in_dim: int, out_dim: int, hidden: int = 512, n_layers: int = 3, dropout: float = 0.1):
        super().__init__()
        layers: list[nn.Module] = []
        d = in_dim
        for _ in range(max(n_layers - 1, 1)):
            layers += [nn.Linear(d, hidden), nn.SiLU(), nn.Dropout(dropout)]
            d = hidden
        layers.append(nn.Linear(d, out_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class VelocityNet(nn.Module):
    def __init__(self, dim: int, cond_dim: int, hidden: int = 512, n_layers: int = 3, dropout: float = 0.1):
        super().__init__()
        self.mlp = MLP(dim + 1 + cond_dim, dim, hidden=hidden, n_layers=n_layers, dropout=dropout)

    def forward(self, z_t: torch.Tensor, t: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        if t.ndim == 1:
            t = t[:, None]
        return self.mlp(torch.cat([z_t, t, cond], dim=-1))


class AdtDecoder(nn.Module):
    """Predict NB mean (and shared theta) from a latent."""

    def __init__(self, z_dim: int, n_adt: int, hidden: int = 512):
        super().__init__()
        self.mu = MLP(z_dim, n_adt, hidden=hidden, n_layers=3, dropout=0.1)
        self.log_theta = nn.Parameter(torch.zeros(n_adt))

    def forward(self, z: torch.Tensor, size_factor: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mu = F.softplus(self.mu(z)) * size_factor.unsqueeze(-1).clamp(min=1e-4)
        theta = F.softplus(self.log_theta).clamp(min=1e-4)
        return mu, theta


def nb_nll(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor) -> torch.Tensor:
    """Mean negative binomial NLL. x, mu: [B, P], theta: [P] or [B, P]."""
    eps = 1e-6
    mu = mu.clamp(min=eps)
    theta = theta.clamp(min=eps)
    t1 = torch.lgamma(theta + x) - torch.lgamma(theta) - torch.lgamma(x + 1.0)
    t2 = theta * torch.log(theta / (theta + mu))
    t3 = x * torch.log(mu / (theta + mu))
    return -(t1 + t2 + t3).mean()


@torch.no_grad()
def integrate_fm(v_net: VelocityNet, x0: torch.Tensor, cond: torch.Tensor, n_steps: int) -> torch.Tensor:
    x = x0
    dt = 1.0 / n_steps
    for k in range(n_steps):
        t = torch.full((x.size(0),), k * dt, device=x.device)
        x = x + dt * v_net(x, t, cond)
    return x
