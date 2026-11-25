#!/usr/bin/env python3
"""Typed states of TEDDY-G 70M for every cell of the THP-1 stratum. Nothing is trained. On the owner's machine.

  python typed.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --cache CACHE --recorded THP1_ENC.npz --out OUT [--batch 32]

Today a head reads one vector for a cell: the mean of the last layer's output over the cell's gene positions. Here
the same output is kept by type, so that the blocks of a typed state add up to that mean:

  TK8   8 types of genes: k-means on the model's own gene embedding (the genes of the screen that are in the
        vocabulary; unit length; 10 starts; seed 20261125). Block k is the sum of the last layer's output over the
        positions whose gene is of type k, divided by the number of the cell's gene positions.        [8 x d]
  TA    2 types: the 128 genes whose presence among a cell's tokens the events of the pathway change most (below),
        and the rest. Block 1 is the sum over the positions of those genes, divided as above; block 2 is the mean
        minus block 1.                                                                                    [2 x d]
  TR    the control of TA: the same with 128 genes drawn at random (seed 20261125).                       [2 x d]
  G     the last layer's output at the position of the gene of the measured protein (CD274 for PD-L1), zeros when
        that gene is not among the cell's tokens, and a last column that says whether it is.             [d + 1]
  EA    the typing of TA at the entry: the same two blocks of what enters the first layer (gene embedding plus
        position embedding).                                                                              [2 x d]

The genes of TA. Among the genes that are a token in at least 5% of the cells of the stratum: the difference between
the share of the knockout cells of the four pathway genes (IFNGR1, IFNGR2, JAK2, STAT1) that have the gene among
their tokens and the same share of the control cells; the 128 genes with the largest absolute difference (ties by
token id). It reads tokens and the labels of those events and of the control cells. **No protein value is read, and no
cell of IRF1 enters the choice.** The 128 of TR are drawn from the same eligible genes.

Checks, written to OUT/typed_manifest.json (the script stops if one fails):
  the mean over gene positions computed here equals the recorded state Z (THP1_ENC.npz, L12), and the mean of
  what enters the first layer the recorded one (E0, where the file has it)                               limit 2e-3
  the blocks of TK8, of TA, of TR and of EA add up to their mean                                          limit 1e-4
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

PAD, LMAX, LIMIT, SEED = 43810, 2048, 2e-3, 20261125
TYPE_A, K, N_GENES, MIN_SHARE = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], 8, 128, 0.05
PROTEIN_GENE = {"PDL1": ("CD274", "ENSG00000120217")}
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--tokens", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--recorded", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", type=int, default=32)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    import torch
    from sklearn.cluster import KMeans
    import sklearn
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    T = np.load(a.tokens)
    tok, lens, ids = T["tok"].astype(np.int64), T["L"].astype(np.int64), T["ids"].astype(str)
    sym, sym_id = T["map_symbols"].astype(str), T["map_token_ids"].astype(np.int64)
    R = np.load(a.recorded)                                       # only ids, A, unit, L12 and E0 are opened: no protein value
    assert np.array_equal(R["ids"].astype(str), ids), "the recorded states are of other cells, or in another order"
    A, un = R["A"].astype(int) == 1, R["unit"].astype(str)
    assert np.array_equal(np.load(a.cache / "index.npz")["ids"].astype(str), ids), "what is kept is of other cells"
    H = np.load(a.cache / "h.npy", mmap_mode="r")
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    n, d, lmax = len(ids), m.d_model, min(LMAX, tok.shape[1])
    nvoc = m.embeddings.weight.shape[0]
    live = np.arange(tok.shape[1])[None, :] < lens[:, None]
    assert (tok[live] != PAD).all()

    # --- the types of genes; tokens and labels only
    genes = np.unique(tok[live])
    assert np.isin(genes, sym_id).all(), "a token of a cell is not a gene of the screen's map"
    count = lambda rows: np.bincount(tok[rows][live[rows]], minlength=nvoc).astype(np.float64)
    srt = np.sort(np.where(live, tok, -1 - np.arange(tok.shape[1])[None, :]), axis=1)
    twice = int((srt[:, 1:] == srt[:, :-1]).sum())               # a gene is a token of a cell once: the rank encoding gives every gene one position
    share_all = count(np.ones(n, bool)) / n
    ctl, typ = ~A, A & np.isin(un, TYPE_A)
    diff = count(typ) / typ.sum() - count(ctl) / ctl.sum()
    eligible = np.where(share_all >= MIN_SHARE)[0]
    order = eligible[np.lexsort((eligible, -np.abs(diff[eligible])))]
    n_g = min(N_GENES, len(eligible))
    gA = np.sort(order[:n_g])
    gR = np.sort(np.random.default_rng(SEED).choice(eligible, n_g, replace=False))
    E = m.embeddings.weight.detach().float().cpu().numpy().astype(np.float64)[genes]
    E /= np.linalg.norm(E, axis=1, keepdims=True) + 1e-12
    km = KMeans(n_clusters=K, n_init=10, random_state=SEED).fit(E)
    type_k = np.full(nvoc, -1, np.int64)
    type_k[genes] = km.labels_
    in_a, in_r = np.zeros(nvoc, bool), np.zeros(nvoc, bool)
    in_a[gA], in_r[gR] = True, True
    name, ens = PROTEIN_GENE["PDL1"]
    own = sym_id[sym == name]
    assert len(own) == 1, f"{name} is not once in the screen's map"
    own = int(own[0])
    if (a.ckpt / "vocab.txt").exists():
        assert te.load_vocab(a.ckpt).get(ens) == own, f"{name} is not {ens} in the vocabulary"
    sym_of = dict(zip(sym_id.tolist(), sym.tolist()))

    # --- the last layer as released, on what is kept of the layer before it
    Zc = np.zeros((n, d), np.float32)
    TK, TA, TR, G = np.zeros((n, K * d), np.float32), np.zeros((n, 2 * d), np.float32), np.zeros((n, 2 * d), np.float32), np.zeros((n, d + 1), np.float32)
    EA, has_e0 = np.zeros((n, 2 * d), np.float32), "E0" in R.files
    tk, ta, tr, tg = (torch.from_numpy(x).to(dev) for x in (type_k, in_a, in_r, np.arange(nvoc) == own))
    last = m.encoder.layers[-1]
    worst = {"mean_against_the_recorded_state": 0.0, "entry_mean_against_the_recorded_one": 0.0, "blocks_add_up_to_the_mean": 0.0}
    order_c = np.argsort(-lens, kind="stable")
    with torch.no_grad():
        for s in range(0, n, a.batch):
            idx = np.sort(order_c[s:s + a.batch])
            Lb = int(lens[idx].max())
            t = torch.from_numpy(tok[idx, :Lb]).to(dev)
            am = torch.from_numpy(live[idx, :Lb]).to(dev)
            t = torch.where(am, t, torch.full_like(t, PAD))
            h = torch.from_numpy(np.asarray(H[idx, :Lb], np.float32)).to(dev)
            o = last(src=h, src_key_padding_mask=~am).float() * am.unsqueeze(-1)
            N = am.sum(1, keepdim=True).float()
            z = o.sum(1) / N
            kk = tk[t]
            bk = torch.stack([(o * (kk == k).unsqueeze(-1)).sum(1) / N for k in range(K)], 1)             # [B, K, d]
            sa, sr = (o * ta[t].unsqueeze(-1)).sum(1) / N, (o * tr[t].unsqueeze(-1)).sum(1) / N
            here = tg[t] & am
            g = (o * here.unsqueeze(-1)).sum(1)
            e = (m.embeddings(t) + m.position_embeddings(torch.arange(Lb, device=dev))).float() * am.unsqueeze(-1)   # what enters the first layer
            e0, ea = e.sum(1) / N, (e * ta[t].unsqueeze(-1)).sum(1) / N
            EA[idx] = torch.cat([ea, e0 - ea], 1).cpu().numpy()
            if has_e0:
                worst["entry_mean_against_the_recorded_one"] = max(worst["entry_mean_against_the_recorded_one"], float(np.abs(e0.cpu().numpy() - R["E0"][idx]).max()))
            worst["blocks_add_up_to_the_mean"] = max(worst["blocks_add_up_to_the_mean"], float((bk.sum(1) - z).abs().max()))
            Zc[idx], TK[idx] = z.cpu().numpy(), bk.reshape(len(idx), -1).cpu().numpy()
            TA[idx], TR[idx] = torch.cat([sa, z - sa], 1).cpu().numpy(), torch.cat([sr, z - sr], 1).cpu().numpy()
            G[idx] = torch.cat([g, here.any(1, keepdim=True).float()], 1).cpu().numpy()
            worst["mean_against_the_recorded_state"] = max(worst["mean_against_the_recorded_state"], float(np.abs(Zc[idx] - R["L12"][idx]).max()))
            if (s // a.batch) % 40 == 0:
                say(f"{min(s + a.batch, n)}/{n} on {dev.type}; largest differences so far {json.dumps({k: round(v, 6) for k, v in worst.items()})}")
    hold = bool(worst["mean_against_the_recorded_state"] <= LIMIT and worst["entry_mean_against_the_recorded_one"] <= LIMIT and worst["blocks_add_up_to_the_mean"] <= 1e-4)
    npos = lambda mask: float((mask[tok] & live).sum(1).mean())
    man = {"cells": int(n), "d_model": int(d), "genes_among_the_tokens": int(len(genes)), "times_a_gene_is_twice_among_a_cells_tokens": twice, "eligible_genes": int(len(eligible)), "min_share": MIN_SHARE,
           "largest_differences": worst, "limits": {"mean_against_the_recorded_state": LIMIT, "entry_mean_against_the_recorded_one": LIMIT, "blocks_add_up_to_the_mean": 1e-4},
           "the_entry_mean_is_recorded": bool(has_e0), "hold": hold,
           "TK8": {"genes_of_a_type": np.bincount(km.labels_, minlength=K).tolist(), "positions_of_a_type_in_a_cell_mean": [npos(type_k == k) for k in range(K)],
                   "sha256_of_the_assignment": hashlib.sha256(type_k.tobytes()).hexdigest()[:16], "inertia": float(km.inertia_)},
           "TA": {"genes": [{"gene": sym_of[int(g_)], "share_of_control_cells": float(count(ctl)[g_] / ctl.sum()), "share_of_pathway_knockout_cells": float(count(typ)[g_] / typ.sum())} for g_ in order[:n_g]],
                  "smallest_absolute_difference_chosen": float(np.abs(diff[order[n_g - 1]])), "positions_in_a_cell_mean": npos(in_a),
                  "holds_the_gene_of_the_protein": bool(in_a[own]), "pathway_knockout_cells": int(typ.sum()), "control_cells": int(ctl.sum())},
           "TR": {"genes": [sym_of[int(g_)] for g_ in gR], "positions_in_a_cell_mean": npos(in_r), "genes_shared_with_TA": int((in_a & in_r).sum())},
           "G": {"gene": name, "vocabulary_entry": ens, "share_of_all_cells_with_the_gene_among_their_tokens": float(G[:, -1].mean())},
           "device": dev.type, "torch": torch.__version__, "numpy": np.__version__, "sklearn": sklearn.__version__, "python": sys.version.split()[0], "seconds": round(time.time() - t0, 1)}
    (a.out / "typed_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    if not hold:
        raise SystemExit("a check does not hold: stopped, no state written")
    np.savez(a.out / "typed.npz", ids=ids, Zc=Zc, TK8=TK, TA=TA, TR=TR, G=G, EA=EA)
    say("done:", json.dumps({k: man[k] for k in ("cells", "largest_differences", "hold", "seconds")}), "| TK8 genes of a type", man["TK8"]["genes_of_a_type"],
        "| TA positions in a cell", round(man["TA"]["positions_in_a_cell_mean"], 1), "| the protein's gene among the tokens of", round(man["G"]["share_of_all_cells_with_the_gene_among_their_tokens"], 3))


if __name__ == "__main__":
    main()
