#!/usr/bin/env python3
"""A check added after the reading, not declared before it: the typing of TA without the knocked-out genes.

  python after_states.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --cache CACHE --recorded THP1_ENC.npz --typed OUT --out OUT_AFTER [--batch 32]

The 128 genes of TA hold IRF1 and JAK2: two genes that an event of this stratum knocks out. A reader fitted with the
cells of IRF1 could then read the knockout from the gene's own transcript, which is not the question. Here:

  TAx   the rule of TA (typed.py) with the five knocked-out genes of the stratum (IFNGR1, IFNGR2, JAK2, STAT1, IRF1)
        not eligible. Block 1 is the sum of the last layer's output over the positions of the 128 genes, divided by
        the number of the cell's gene positions; block 2 is the mean minus block 1.                       [2 x d]
  EAx   the same typing of what enters the first layer.                                                   [2 x d]

Nothing else changes. Nothing is trained. No protein value is read.

Checks, written to OUT_AFTER/after_manifest.json (the script stops if one fails):
  without the exclusion the rule gives the genes of TA that OUT/typed_manifest.json records;
  the TA and the EA computed here equal those of typed.py (OUT/typed.npz)                                  limit 1e-4
  the mean over gene positions equals the recorded state (THP1_ENC.npz, L12)                              limit 2e-3
"""
import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from typed import LIMIT, LMAX, MIN_SHARE, N_GENES, PAD, TYPE_A  # noqa: E402

WITHOUT = ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--tokens", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--recorded", type=Path, required=True)
    ap.add_argument("--typed", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--without", nargs="*", default=WITHOUT, help="other names only for a trial on made-up values")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    T = np.load(a.tokens)
    tok, lens, ids = T["tok"].astype(np.int64), T["L"].astype(np.int64), T["ids"].astype(str)
    sym, sym_id = T["map_symbols"].astype(str), T["map_token_ids"].astype(np.int64)
    R = np.load(a.recorded)                                       # only ids, A, unit and L12 are opened: no protein value
    assert np.array_equal(R["ids"].astype(str), ids), "the recorded states are of other cells, or in another order"
    A, un = R["A"].astype(int) == 1, R["unit"].astype(str)
    assert np.array_equal(np.load(a.cache / "index.npz")["ids"].astype(str), ids), "what is kept is of other cells"
    H = np.load(a.cache / "h.npy", mmap_mode="r")
    P = np.load(a.typed / "typed.npz")
    assert np.array_equal(P["ids"].astype(str), ids)
    was = json.loads((a.typed / "typed_manifest.json").read_text())
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    n, d = len(ids), m.d_model
    nvoc = m.embeddings.weight.shape[0]
    live = np.arange(tok.shape[1])[None, :] < lens[:, None]

    # --- the rule of TA, and the same with the knocked-out genes not eligible; tokens and labels only
    count = lambda rows: np.bincount(tok[rows][live[rows]], minlength=nvoc).astype(np.float64)
    share_all = count(np.ones(n, bool)) / n
    ctl, typ = ~A, A & np.isin(un, TYPE_A)
    diff = count(typ) / typ.sum() - count(ctl) / ctl.sum()
    eligible = np.where(share_all >= MIN_SHARE)[0]
    sym_of = dict(zip(sym_id.tolist(), sym.tolist()))
    out_ids = sym_id[np.isin(sym, a.without)]

    def choose(el):
        order = el[np.lexsort((el, -np.abs(diff[el])))]
        return order[:min(N_GENES, len(el))]

    oA, oX = choose(eligible), choose(eligible[~np.isin(eligible, out_ids)])
    same_rule = [sym_of[int(g)] for g in oA] == [g["gene"] for g in was["TA"]["genes"]]
    in_a, in_x = np.zeros(nvoc, bool), np.zeros(nvoc, bool)
    in_a[oA], in_x[oX] = True, True

    TA, TAx, EA, EAx = (np.zeros((n, 2 * d), np.float32) for _ in range(4))
    Zc = np.zeros((n, d), np.float32)
    ta, tx = torch.from_numpy(in_a).to(dev), torch.from_numpy(in_x).to(dev)
    last = m.encoder.layers[-1]
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
            e = (m.embeddings(t) + m.position_embeddings(torch.arange(Lb, device=dev))).float() * am.unsqueeze(-1)
            e0 = e.sum(1) / N
            for mask, To, Eo in ((ta, TA, EA), (tx, TAx, EAx)):
                sa, ea = (o * mask[t].unsqueeze(-1)).sum(1) / N, (e * mask[t].unsqueeze(-1)).sum(1) / N
                To[idx], Eo[idx] = torch.cat([sa, z - sa], 1).cpu().numpy(), torch.cat([ea, e0 - ea], 1).cpu().numpy()
            Zc[idx] = z.cpu().numpy()
            if (s // a.batch) % 40 == 0:
                say(f"{min(s + a.batch, n)}/{n} on {dev.type}")
    worst = {"TA_against_the_one_read": float(np.abs(TA - P["TA"]).max()), "EA_against_the_one_read": float(np.abs(EA - P["EA"]).max()),
             "mean_against_the_recorded_state": float(np.abs(Zc - R["L12"]).max())}
    hold = bool(same_rule and worst["TA_against_the_one_read"] <= 1e-4 and worst["EA_against_the_one_read"] <= 1e-4 and worst["mean_against_the_recorded_state"] <= LIMIT)
    npos = lambda mask: float((mask[tok] & live).sum(1).mean())
    left, came = [g for g in oA if g not in set(oX.tolist())], [g for g in oX if g not in set(oA.tolist())]
    row = lambda g: {"gene": sym_of[int(g)], "share_of_control_cells": float(count(ctl)[g] / ctl.sum()), "share_of_pathway_knockout_cells": float(count(typ)[g] / typ.sum())}
    man = {"note": "a check added after the reading; not declared before it", "cells": int(n), "without": list(a.without), "of_them_in_the_map_of_the_screen": [sym_of[int(g)] for g in out_ids],
           "the_rule_without_the_exclusion_gives_the_genes_of_TA": bool(same_rule), "largest_differences": worst, "hold": hold,
           "TAx": {"genes": [sym_of[int(g)] for g in oX], "genes_of_TA_left_out": [row(g) for g in left], "genes_that_took_their_place": [row(g) for g in came],
                   "positions_in_a_cell_mean": npos(in_x), "positions_in_a_cell_mean_of_TA": npos(in_a),
                   "share_of_cells_with_a_gene_left_out_among_their_tokens": {sym_of[int(g)]: {u: float(((tok == g) & live)[un == u].any(1).mean()) for u in sorted(set(un.tolist()))} for g in left}},
           "device": dev.type, "seconds": round(time.time() - t0, 1)}
    (a.out / "after_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    if not hold:
        raise SystemExit("a check does not hold: stopped, no state written")
    np.savez(a.out / "after.npz", ids=ids, TAx=TAx, EAx=EAx)
    say("done:", json.dumps(worst), "| left out", [sym_of[int(g)] for g in left], "| in their place", [sym_of[int(g)] for g in came])


if __name__ == "__main__":
    main()
