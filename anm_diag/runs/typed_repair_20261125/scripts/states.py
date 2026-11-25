#!/usr/bin/env python3
"""Parts B and D: typed states of TEDDY-G 70M for one stratum, with the types chosen on the fitted half only.
Nothing is trained. On the owner's machine.

  python states.py --repo TEDDY_MM --ckpt TEDDY_G_70M --bundle BUNDLE --stratum thp1|mel_co|mel_ifn --tokens TOKENS.npz --labels LABELS.npz [--cache CACHE] --out OUT

BUNDLE: a copy of this repository's anm_diag folder (for the halves of the round). LABELS: for thp1 the arrays of
runs/encoder_stages_20261124 (ids, A, unit, L12); for a melanoma stratum the file of strata (keys <stratum>__ids, __A,
__unit, __Z). CACHE: for thp1, the output of the layer before the last (runs/last_layer_trial_20261124); without it
all twelve layers are run.

The halves are those of the round (anm_diag.cycle.halves on the whole stratum). For each half f, on the cells of that
half only, and with tokens and labels only (no protein value is opened):

  A        the 128 genes whose presence among a cell's tokens differs most between the knockout cells of the recorded
           events of the type and the control cells (genes among the tokens of at least 5% of the half's cells)
  A-g      the same with the cells of event g left out, for every event g of the type
  R        128 of the eligible genes at random: the control
  P1..P4   four programs of the control cells, with no record: the presence matrix of the half's control cells (genes
           present in 5% to 95% of them), each gene weighted by log(1 / its share), factorised into four non-negative
           programs; a gene belongs to the program where its loading is largest, and the 128 genes with the largest
           loading of a program form its block

For every cell of the stratum and every set of genes above: the sum of the last layer's output over the positions of
those genes, divided by the number of the cell's gene positions. Checks (the script stops if one fails): the mean over
gene positions equals the recorded state to 2e-3.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np

PAD, LMAX, LIMIT, SEED_T = 43810, 2048, 2e-3, 20261125
N_GENES, MIN_SHARE, MAX_SHARE_P, K = 128, 0.05, 0.95, 4
TYPES = {"thp1": ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], "mel_co": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"], "mel_ifn": ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]}
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--stratum", required=True, choices=list(TYPES))
    ap.add_argument("--tokens", type=Path, required=True)
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", type=int, default=32)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    import torch
    from sklearn.decomposition import NMF
    sys.path.insert(0, str(a.bundle))
    sys.path.insert(0, str(a.bundle / "runs" / "teg_training_20261119" / "scripts"))
    from anm_diag import cycle
    from _common import SEED
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    s, typ = a.stratum, TYPES[a.stratum]
    Lb = np.load(a.labels)
    if s == "thp1":
        ids, A, un, Zrec = Lb["ids"].astype(str), Lb["A"].astype(int) == 1, Lb["unit"].astype(str), Lb["L12"]
    else:
        ids, A, un, Zrec = Lb[f"{s}__ids"].astype(str), Lb[f"{s}__A"].astype(int) == 1, Lb[f"{s}__unit"].astype(str), Lb[f"{s}__Z"]
    T = np.load(a.tokens)
    tid = T["ids"].astype(str)
    pos = {c: i for i, c in enumerate(tid)}
    take = np.array([pos[c] for c in ids])
    tok, lens = T["tok"][take].astype(np.int64), T["L"][take].astype(np.int64)
    n = len(ids)
    live = np.arange(tok.shape[1])[None, :] < lens[:, None]
    assert (tok[live] != PAD).all()
    H = None
    if a.cache is not None:
        assert np.array_equal(np.load(a.cache / "index.npz")["ids"].astype(str), ids), "what is kept is of other cells"
        H = np.load(a.cache / "h.npy", mmap_mode="r")
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    d, nvoc = m.d_model, m.embeddings.weight.shape[0]
    names = {}
    if "map_symbols" in T.files:
        names = dict(zip(T["map_token_ids"].astype(np.int64).tolist(), T["map_symbols"].astype(str).tolist()))
    elif (a.ckpt / "vocab.txt").exists():
        names = {int(v): k for k, v in te.load_vocab(a.ckpt).items()}
    name = lambda g: names.get(int(g), f"token {int(g)}")

    # --- the halves of the round, and the sets of genes of each half; tokens and labels only
    first = cycle.halves(np.where(A, un, "control"), SEED)
    count = lambda rows: np.bincount(tok[rows][live[rows]], minlength=nvoc).astype(np.float64)
    sets, man_sets = [], {}
    for f, fit in ((0, first), (1, ~first)):
        ctl = fit & ~A
        share_all, share_ctl = count(fit) / fit.sum(), count(ctl) / ctl.sum()
        eligible = np.where(share_all >= MIN_SHARE)[0]
        rng = np.random.default_rng(SEED_T + f)

        def by_events(events):
            ko = fit & A & np.isin(un, events)
            diff = count(ko) / ko.sum() - share_ctl
            order = eligible[np.lexsort((eligible, -np.abs(diff[eligible])))]
            return np.sort(order[:min(N_GENES, len(eligible))])

        fold = {"A": by_events(typ), **{f"A-{g}": by_events([h for h in typ if h != g]) for g in typ}, "R": np.sort(rng.choice(eligible, min(N_GENES, len(eligible)), replace=False))}
        gp = np.where((share_ctl >= MIN_SHARE) & (share_ctl <= MAX_SHARE_P))[0]                 # the programs: control cells of the half, no record
        col = {g: j for j, g in enumerate(gp.tolist())}
        X = np.zeros((int(ctl.sum()), len(gp)), np.float32)
        for i, r in enumerate(np.where(ctl)[0]):
            js = [col[g] for g in tok[r, :lens[r]].tolist() if g in col]
            X[i, js] = 1.0
        X *= np.log(1.0 / share_ctl[gp])[None, :].astype(np.float32)
        nmf = NMF(n_components=min(K, len(gp)), init="nndsvda", random_state=SEED_T, max_iter=500)
        nmf.fit(X)
        Hk = nmf.components_
        owner = Hk.argmax(0)
        for k in range(Hk.shape[0]):
            mine = np.where(owner == k)[0]
            top = mine[np.argsort(-Hk[k, mine], kind="stable")[:N_GENES]]
            fold[f"P{k + 1}"] = np.sort(gp[top])
        for k_, g_ in fold.items():
            sets.append((f"{k_}|{f}", g_))
        man_sets[str(f)] = {"cells": int(fit.sum()), "control_cells": int(ctl.sum()), "eligible_genes": int(len(eligible)), "genes_of_the_programs_matrix": int(len(gp)), "nmf_iterations": int(nmf.n_iter_),
                            "nmf_error": float(nmf.reconstruction_err_), "sets": {k_: {"genes": [name(g) for g in g_][:N_GENES], "n": int(len(g_)), "positions_in_a_cell_mean": float((np.isin(tok, g_) & live).sum(1).mean())}
                                                                                 for k_, g_ in fold.items()},
                            "genes_of_A_shared_with": {k_: int(len(np.intersect1d(fold["A"], g_))) for k_, g_ in fold.items() if k_ != "A"}}
        say(f"half {f}: {int(fit.sum())} cells, {len(eligible)} eligible genes; programs matrix {X.shape}; NMF {nmf.n_iter_} iterations")
    both = {k_.split("|")[0] for k_, _g in sets}
    shared = {k_: int(len(np.intersect1d(dict(sets)[f"{k_}|0"], dict(sets)[f"{k_}|1"]))) for k_ in sorted(both)}

    # --- the encoder as released
    masks = torch.zeros((len(sets), nvoc), dtype=torch.bool)
    for j, (_k, g_) in enumerate(sets):
        masks[j, torch.from_numpy(g_)] = True
    masks = masks.to(dev)
    Zc, B = np.zeros((n, d), np.float32), np.zeros((n, len(sets), d), np.float32)
    worst = 0.0
    order_c = np.argsort(-lens, kind="stable")
    layers = m.encoder.layers
    with torch.no_grad():
        for s0 in range(0, n, a.batch):
            idx = np.sort(order_c[s0:s0 + a.batch])
            Lm = int(lens[idx].max())
            t = torch.from_numpy(tok[idx, :Lm]).to(dev)
            am = torch.from_numpy(live[idx, :Lm]).to(dev)
            t = torch.where(am, t, torch.full_like(t, PAD))
            if H is not None:
                h = layers[-1](src=torch.from_numpy(np.asarray(H[idx, :Lm], np.float32)).to(dev), src_key_padding_mask=~am)
            else:
                h = m.embeddings(t) + m.position_embeddings(torch.arange(Lm, device=dev))
                for layer in layers:
                    h = layer(src=h, src_key_padding_mask=~am)
            o = h.float() * am.unsqueeze(-1)
            N = am.sum(1, keepdim=True).float()
            Zc[idx] = (o.sum(1) / N).cpu().numpy()
            mt = masks[:, t]                                                                    # [sets, B, L]
            B[idx] = (torch.bmm(mt.float().permute(1, 0, 2), o) / N.unsqueeze(1)).cpu().numpy()
            worst = max(worst, float(np.abs(Zc[idx] - Zrec[idx]).max()))
            if (s0 // a.batch) % 20 == 0:
                say(f"{min(s0 + a.batch, n)}/{n} on {dev.type}; mean against the recorded state, largest difference so far {worst:.6f}")
    hold = bool(worst <= LIMIT)
    man = {"stratum": s, "cells": int(n), "d_model": int(d), "the_recorded_events_of_the_type": typ, "layers_run": "the last one, on what was kept of the layer before it" if H is not None else "all twelve",
           "mean_against_the_recorded_state": worst, "limit": LIMIT, "hold": hold, "halves": man_sets, "genes_shared_by_the_two_halves": shared, "sets": [k_ for k_, _g in sets],
           "sha256_of_the_halves": hashlib.sha256(first.tobytes()).hexdigest()[:16], "device": dev.type, "torch": torch.__version__, "seconds": round(time.time() - t0, 1)}
    (a.out / f"{s}_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    if not hold:
        raise SystemExit("the mean does not equal the recorded state: stopped, no state written")
    np.savez(a.out / f"{s}_typed.npz", ids=ids, first=first, Zc=Zc, sets=np.array([k_ for k_, _g in sets]), blocks=B)
    say("done:", json.dumps({"cells": n, "sets": len(sets), "mean_against_the_recorded_state": worst, "genes_shared_by_the_two_halves": shared}))


if __name__ == "__main__":
    main()
