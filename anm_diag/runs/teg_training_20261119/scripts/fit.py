"""Fit the head's MLP once per arm and held-out gene (DESIGN.md). Every arm has the same cells, split, optimiser, steps,
stopping rule, and the same two holding terms; the arms differ only in what the supervision of the PD-L1 output carries.

  python fit.py --check                 the head in torch against numpy; the reading of the head as trained against the record of
                                        runs/scope_20261116; 20 steps of every arm on training cells, for timing. Writes nothing.
  python fit.py <data set> [arm ...]    the fits of the named arms (default: all), one per held-out gene. Writes
                                        <out>/heads/<data set>__<gene>__<arm>.npz and <out>/heads/fits_<data set>__<arm>.json

No audit cell enters a fit or a stopping check. The audit cells are read by read.py, once.
"""
import json
import sys
import time

import numpy as np
import torch
import torch.nn.functional as Fn

from _common import (ARMS, BM_FIT, EVERY, F, FEATURES, GUARD_HEADS, KEEP, LAM, LAYERS, ACT, LR, MLP, OUT, PATIENCE, RUN, SEED, STEPS, TARGETS,
                     bm_score, forward_np, load, own_output_u, unit, weights)

torch.set_num_threads(max(torch.get_num_threads(), 1))
T32 = lambda a: torch.from_numpy(np.ascontiguousarray(a, dtype=F))


def head(P, x):
    """The head in torch: latent and output (log of the predicted mean)."""
    h, lat = x, None
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T + P[name + ".bias"]
        h = Fn.silu(h) if act else h
        if name == "mlp__net.6":
            lat = h
    return lat, torch.log(torch.clamp(Fn.softplus(h), min=1e-6))


def graph_similarity(genes):
    """What an embedding fitted to the graph holds about these genes. The graph a knowledge graph would carry here gives every
    gene the same two edges (it `regulates` PD-L1, it is `part_of` the IFN-gamma pathway), and an embedding fitted to a graph
    places entities with the same edges at the same point. So the similarity is the cosine of the genes' edge profiles: 1 for
    every pair. Nothing measured enters."""
    edges = [("regulates", "PD-L1"), ("part_of", "IFN-gamma pathway")]
    M = np.ones((len(genes), len(edges)))                       # every gene has both edges
    M = M / np.linalg.norm(M, axis=1, keepdims=True)
    return np.clip(M @ M.T, 0, 1)


def judge(x, y, A, fit, train_ctl):
    """A learned judgment model: a logistic regression from TEDDY's embedding to the judgment 'PD-L1 is low' (below the
    median of the training controls), fitted on the cells of the fit. Returns its probability for every cell."""
    from sklearn.linear_model import LogisticRegression
    cut = np.median(y[fit & train_ctl])
    m = LogisticRegression(C=1.0, max_iter=5000).fit(x[fit].astype(np.float64), (y[fit] < cut).astype(int))
    return m.predict_proba(x.astype(np.float64))[:, 1], float(cut), float(m.score(x[fit].astype(np.float64), (y[fit] < cut).astype(int)))


class Fit:
    """One data set with one gene held out: the cells of the fit, the set-aside cells, and the first term of every arm."""

    def __init__(self, ds, d, P0, j, held, bm):
        self.ds, self.d, self.j, self.held = ds, d, j, held
        self.genes = [g for g in TARGETS[ds] if g != held]
        A, tg = d["A"], d["tg"]
        train = d["train_ctl"] | ((A == 1) & (tg != held))
        self.sets = {"fit": train & ~d["val"], "val": train & d["val"]}
        self.P0 = {k: T32(v) for k, v in P0.items()}
        with torch.no_grad():
            lat0, out0 = head(self.P0, T32(d["x"]))
        o0 = out0.numpy(); c = d["train_ctl"]; y = d["y"]
        self.t = (o0[c, j].mean() + o0[c, j].std() * (y - y[c].mean()) / y[c].std()).astype(F)                  # the target of the repair already made
        fc = self.sets["fit"] & c
        self.s0 = float(o0[fc, j].std()); self.c0 = float(np.median(o0[fc, j])); self.tau = self.s0 * np.sqrt(3) / np.pi
        self.S = graph_similarity(TARGETS[ds])                                                                   # the graph knows all five genes
        self.p, self.cut, self.judge_acc = judge(d["x"], y, A, self.sets["fit"], c)
        self.xb, self.ob0, self.bval = bm
        self.batch = {}
        for name, rows in self.sets.items():
            r = np.where(rows)[0]
            self.batch[name] = dict(x=T32(d["x"][r]), out0=out0[r], t=T32(self.t[r]), y=torch.from_numpy(y[r]), A=torch.from_numpy(A[r].astype(np.float32)),
                                    p=T32(self.p[r]), ctl=torch.from_numpy(A[r] == 0), ev={g: torch.from_numpy((A[r] == 1) & (tg[r] == g)) for g in self.genes},
                                    xb=self.xb[~self.bval if name == "fit" else self.bval], ob0=self.ob0[~self.bval if name == "fit" else self.bval])

    def first(self, arm, b, out, lat, probe):
        """The supervision of the PD-L1 output. The control cells are fitted to their measured protein, cell by cell, in every
        arm (as in the repair already made). The arms differ in what they carry about the knockout cells."""
        oj, n = out[:, self.j], out.shape[0]
        ctl, trt = b["ctl"], ~b["ctl"]
        sq = (oj - b["t"]) ** 2
        control_side = (sq * ctl).sum() / n
        ev = [m for m in b["ev"].values() if int(m.sum()) > 1]

        def rows():                                          # for every event of the fit, u of the head's own output, squared
            pred = torch.exp(oj.double()); y = b["y"]
            pc, yc = pred[ctl], y[ctl]
            a = ((pc - pc.mean()) * (yc - yc.mean())).mean() / torch.clamp(((pc - pc.mean()) ** 2).mean(), min=1e-12)
            u = torch.stack([1 - a * (pred[e].mean() - pc.mean()) / (y[e].mean() - yc.mean()) for e in ev])
            return (u ** 2).mean().float()
        if arm == "data":                                    # the measured protein of every knockout cell (the repair already made)
            ko = (sq * trt).sum() / n
        elif arm == "class":                                 # the class of every cell, read from the latent by a linear probe
            ko = Fn.binary_cross_entropy_with_logits(lat @ probe["w"] + probe["b"], b["A"])
        elif arm.startswith("edge") or arm == "graph":       # a signed edge: the knockout lowers PD-L1, by at least a fixed margin
            m = {"edge05": 0.5, "edge1": 1.0, "edge2": 2.0, "graph": 1.0}[arm]
            sh = torch.stack([oj[e].mean() - oj[ctl].mean() for e in ev]) / self.s0
            ko = (torch.relu(sh + m) ** 2).mean()
            if arm == "graph":                               # and the geometry of the graph: events the graph holds alike move the latent alike
                dl = [lat[e].mean(0) - lat[ctl].mean(0) for e in ev]; gi = [TARGETS[self.ds].index(g) for g, e in b["ev"].items() if int(e.sum()) > 1]
                pen = [self.S[gi[a], gi[c]] * (1 - Fn.cosine_similarity(dl[a], dl[c], dim=0)) for a in range(len(dl)) for c in range(a + 1, len(dl))]
                ko = ko + torch.stack(pen).mean()
        elif arm == "judge":                                 # the probabilities of a learned judgment model, for the knockout cells
            ko = Fn.binary_cross_entropy_with_logits((self.c0 - oj[trt]) / self.tau, b["p"][trt])
        elif arm == "rows":                                  # the rows of the TEG
            ko = rows()
        elif arm == "data_rows":
            ko = (sq * trt).sum() / n + rows()
        else:
            raise ValueError(arm)
        return control_side + ko

    def loss(self, arm, P, probe, name):
        b = self.batch[name]
        lat, out = head(P, b["x"])
        other = torch.ones(out.shape[1], dtype=torch.bool); other[self.j] = False
        keep = ((out[:, other] - b["out0"][:, other]) ** 2).mean()
        hold = ((head(P, b["xb"])[1] - b["ob0"]) ** 2).mean()
        return self.first(arm, b, out, lat, probe) + KEEP * keep + LAM * hold

    def run(self, arm, steps=STEPS):
        """Adam on the three MLP layers only (and the probe of arm `class`); kept: the weights with the lowest loss on the set-aside cells."""
        P = {k: v.clone() for k, v in self.P0.items()}
        keys = [l + s for l in MLP for s in (".weight", ".bias")]
        for k in keys:
            P[k].requires_grad_(True)
        probe = {"w": torch.zeros(P["mlp__net.6.bias"].shape[0], requires_grad=True), "b": torch.zeros((), requires_grad=True)}
        opt = torch.optim.Adam([P[k] for k in keys] + (list(probe.values()) if arm == "class" else []), lr=LR, betas=(0.9, 0.999), eps=1e-8)
        best, best_P, curve, bad, t0 = np.inf, None, [], 0, time.time()
        for step in range(1, steps + 1):
            opt.zero_grad()
            L = self.loss(arm, P, probe, "fit")
            L.backward()
            opt.step()
            if step % EVERY == 0:
                with torch.no_grad():
                    Lv = float(self.loss(arm, P, probe, "val"))
                curve.append([step, float(L), Lv])
                if Lv < best:
                    best, best_P, bad = Lv, {k: v.detach().clone() for k, v in P.items()}, 0
                else:
                    bad += 1
                    if bad >= PATIENCE:
                        break
        out = {k: v.numpy().astype(F) for k, v in best_P.items()}
        assert all(np.array_equal(out[k], self.P0[k].numpy()) for k in out if k not in keys)             # the decoder is bit for bit the head as trained
        with torch.no_grad():
            moved = float(((head(best_P, self.xb)[1] - self.ob0) ** 2).mean().sqrt())
        log = {"arm": arm, "held_out": self.held, "events_of_the_fit": self.genes, "n_fit": int(self.sets["fit"].sum()), "n_set_aside": int(self.sets["val"].sum()),
               "curve": curve, "kept_step": int(min(curve, key=lambda c: c[2])[0]), "seconds": time.time() - t0, "bone_marrow_outputs_moved": moved}
        if arm == "judge":
            log["judge"] = {"cut": self.cut, "accuracy_on_the_cells_of_the_fit": self.judge_acc}
        if arm == "graph":
            log["graph_similarity_of_the_genes"] = {"genes": TARGETS[self.ds], "cosine": np.round(self.S, 4).tolist()}
        return out, log


def bone_marrow(P0):
    b = np.load(BM_FIT)
    xb = T32(unit(b["z_train"]))
    with torch.no_grad():
        ob0 = head({k: T32(v) for k, v in P0.items()}, xb)[1]
    bval = np.zeros(len(xb), bool); bval[np.random.default_rng(SEED).permutation(len(xb))[: len(xb) // 5]] = True
    return (xb, ob0, torch.from_numpy(bval)), b


def main():
    P0, names, j = weights()
    bm, b = bone_marrow(P0)
    if "--check" in sys.argv:
        for ds in TARGETS:
            d = load(ds)
            with torch.no_grad():
                lat, out = head({k: T32(v) for k, v in P0.items()}, T32(d["x"]))
            lat_n, out_n = forward_np(P0, d["x"])
            print(ds, "| torch against numpy, largest difference: latent", f"{np.abs(lat.numpy() - lat_n).max():.2e}", "output", f"{np.abs(out.numpy() - out_n).max():.2e}")
            use = (d["A"] == 1) | d["audit_ctl"]
            u, a = own_output_u(d["y"][use], np.exp(out_n[use, j]), d["A"][use], d["tg"][use])
            rec = json.load(open(RUN.parent / "scope_20261116" / "results" / "frozen_arms.json"))[{"frangieh2021": "melanoma", "papalexi2021": "thp1"}[ds]]["cells of the loop (half of the controls)"]["as_trained"]
            print(f"   the head as trained, own output: u {u[0]:.6f}, slope {a:.4f} | on record (runs/scope_20261116): u {rec['u']:.6f}, slope {rec['slope']:.4f}")
            assert abs(u[0] - rec["u"]) < 1e-4, "the reading does not reproduce the record"
        d = load("frangieh2021"); f = Fit("frangieh2021", d, P0, j, "IFNGR1", bm)
        print("graph similarity of the five genes, lowest:", round(float(f.S.min()), 3), "| judge accuracy on the cells of its fit:", round(f.judge_acc, 3))
        for arm in ARMS:
            _, log = f.run(arm, 20)
            print(f"   {arm:10s} 20 steps: {log['seconds']:.1f}s; loss {log['curve'][0][1]:.4f}, set-aside {log['curve'][0][2]:.4f}")
        return
    ds = sys.argv[1]; arms = sys.argv[2:] or ARMS
    d = load(ds); hdir = OUT / "heads"; hdir.mkdir(parents=True, exist_ok=True)
    xv, yv = unit(b["z_val"]), b["adt_val"].astype(np.float64)
    for arm in arms:
        logs = {}
        for held in TARGETS[ds]:
            P, log = Fit(ds, d, P0, j, held, bm).run(arm)
            log["validation_score"] = bm_score(P, xv, yv, names)["mean_pearson"]
            np.savez(hdir / f"{ds}__{held}__{arm}.npz", **P); logs[held] = log
            print(f"{ds}|{arm}|{held}: kept step {log['kept_step']}, set-aside loss {log['curve'][0][2]:.4f} to {min(c[2] for c in log['curve']):.4f}, "
                  f"bone-marrow outputs moved {log['bone_marrow_outputs_moved']:.4f}, {log['seconds']:.0f}s", flush=True)
        json.dump(logs, open(hdir / f"fits_{ds}__{arm}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
