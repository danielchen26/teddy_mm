"""Synthetic adapter with known ground truth, used to validate the runner.

Each cell has a latent program score p = effect * A * strength(target) + noise. The input x holds p in a
few genes plus nuisance; the readout O measures p with independent noise (its own measurement). The
chain is INPUT_VIEW -> KEEP (keeps p) -> LOSE (removes p) -> FINAL. With `loss_stage` the program is
dropped at that stage, so the runner must localise the loss there; with effect 0 it must report no_effect;
with loss_stage None it must report no_deficit.
"""
import numpy as np

from .base import Adapter


class SyntheticAdapter(Adapter):
    def __init__(self, n_targets=12, cells_per_target=60, n_controls=900, p=40, effect=1.0,
                 loss_stage="LOSE", seed=0):
        rng = np.random.default_rng(seed)
        self.loss_stage = loss_stage
        nT = n_targets * cells_per_target
        self.A = np.r_[np.ones(nT, int), np.zeros(n_controls, int)]
        self.cluster = np.r_[np.repeat([f"t{i}" for i in range(n_targets)], cells_per_target),
                             [f"c{i}" for i in range(n_controls)]]
        strength = np.r_[np.repeat(rng.uniform(0.5, 1.5, n_targets), cells_per_target), np.zeros(n_controls)]
        n = len(self.A)
        self.prog = effect * strength * self.A + rng.normal(0, 1.0, n)        # latent program
        nuis = rng.normal(0, 1, (n, p))
        nuis[:, :3] *= 0.05                                                    # program genes measured precisely
        self.x = nuis.copy()
        self.x[:, :3] += self.prog[:, None]                                   # program enters three genes
        self.x[:, 3:6] += 0.5 * self.A[:, None] * strength[:, None]           # an unrelated perturbation signal
        self.O = self.prog + rng.normal(0, 0.3, n)                             # independent measurement of p

    def cells(self, decl):
        return {"A": self.A, "cluster": self.cluster, "ids": np.arange(len(self.A))}

    def readout(self, rspec, cells):
        return self.O

    def _stage(self, node, cells, arg):
        nid = node["id"]
        order = ["INPUT_VIEW", "KEEP", "LOSE", "FINAL"]
        X = self.x.copy()
        if self.loss_stage is not None and order.index(nid) >= order.index(self.loss_stage):
            X[:, :3] = np.random.default_rng(1).normal(0, 1, (len(X), 3))     # program removed from here on
        return X

    def _pc(self, node, cells, arg):
        return self.x[:, :3]

    def _shuffle(self, node, cells, arg):
        return np.random.default_rng(2).permutation(self.x[:, :3].ravel()).reshape(-1, 3)

    ops = {"view": _stage, "stage": _stage, "program_genes": _pc, "shuffled": _shuffle}

    def provenance(self):
        return {"adapter": "synthetic", "loss_stage": self.loss_stage}


def synthetic_declaration(name="synthetic"):
    return {
        "anm_declaration": 1, "name": name,
        "system": {"id": "toy", "frozen": True, "inputs": ["x"]},
        "measurements": {"x": {"modality": "rna", "kind": "values"}, "o": {"modality": "protein", "kind": "values"}},
        "source": {"controlled": True, "kind": "simulated", "control": "control", "treated": "t*", "cluster": "target"},
        "readouts": [{"id": "O", "measurement": "o", "statistic": "value", "features": "all"}],
        "nodes": [
            {"id": "INPUT_VIEW", "from": "x", "op": "view", "operator": "analysis", "role": "chain"},
            {"id": "KEEP", "from": "INPUT_VIEW", "op": "stage", "operator": "system", "role": "chain"},
            {"id": "LOSE", "from": "KEEP", "op": "stage", "operator": "system", "role": "chain"},
            {"id": "FINAL", "from": "LOSE", "op": "stage", "operator": "system", "role": "chain"},
            {"id": "PROGRAM", "from": "x", "op": "program_genes", "operator": "analysis", "role": "positive_control"},
        ],
        "chain": ["INPUT_VIEW", "KEEP", "LOSE", "FINAL"],
        "inference": {"bootstrap": 300, "forward_bootstrap": 100},
    }
