#!/usr/bin/env python3
"""The states of one checkpoint of TEDDY-G (70M, 160M or 400M) for the cells of a stratum, on the owner's machine.

  python sizes.py --repo TEDDY_MM --ckpt CHECKPOINT --tokens TOKENS.npz --out OUT_DIR [--batch 16]
                  [--official OFFICIAL_CODE --check 16] [--recorded LAYERS_B.npz]

Reads tokens only (TOKENS.npz of runs/encoder_stages_20261124/scripts/layers.py; the three checkpoints share one
vocabulary, so the tokens of a cell are the same for all of them). No protein value and no label of an event is opened.
The forward pass is that of runs/annotation_input_20261124/scripts/annot.py, which equals the official class, written
for any number of layers. Two ways of running, each cell once in each:

  N  no annotation positions: how the repository and the official tutorial run the model
  P  the four placeholders <disease> <tissue> <cell> <sex> in front of the genes

Kept for every cell (OUT_DIR/states.npz), each a mean over the cell's gene tokens unless said otherwise:

  E0   what enters the first layer (gene embedding plus position embedding), in N
  LP   the output of the layer before the last, in N
  Z    the output of the last layer, in N: the state a head would read
  ZP   the same with the placeholders in front of the genes
  C    the output of the last layer at the first annotation position (the disease placeholder), in P: what the
       official class returns as `cell_emb`
  top, lab   what the model's own output layer ranks first at the four placeholders, and its probabilities for the 43
             label tokens

--check K: only the checks, on K cells (the K/2 longest and the K/2 shortest); agreement numbers are printed and
nothing else. With --official the official class itself is loaded (it needs `transformers` below version 5) and its
output is compared position by position. With --recorded (the 70M checkpoint only) Z is compared with the recorded state.
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

PAD = 43810
PLACE = ["<disease>", "<tissue>", "<cell>", "<sex>"]              # the order of the official tokeniser: disease, tissue, cell type, sex
LABELS = (43761, 43804)                                          # label tokens: 11 diseases, 17 tissues, 12 kinds of cell, female, male, unknown
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def forward(m, te, t, am, ann, torch):
    """t [B, L] gene tokens (padding already PAD), am [B, L] 1 where a token is. ann: None or the four token ids.
    Returns the output of the last layer [B, k + L, d], the gene-means of what enters layer 1 and of the output of
    the layer before the last, and k."""
    dev = t.device
    h = m.embeddings(t) + m.position_embeddings(torch.arange(t.size(1), device=dev))
    e0 = te.pool_hidden(h.float(), am, "gene-mean")
    pad = ~am.bool()
    k = 0
    if ann is not None:
        k = int(ann.numel())
        h = torch.cat([m.embeddings(ann)[None].expand(t.size(0), -1, -1), h], dim=1)
        pad = torch.cat([torch.zeros(t.size(0), k, dtype=torch.bool, device=dev), pad], dim=1)
    n, lp = len(m.encoder.layers), None
    for i, layer in enumerate(m.encoder.layers):
        h = layer(src=h, src_key_padding_mask=pad)
        if i == n - 2:
            lp = te.pool_hidden(h[:, k:].float(), am, "gene-mean")
    return h, e0, lp, k


def batch(tok, lens, idx, torch, dev):
    Lb = int(lens[idx].max())
    t = torch.from_numpy(tok[idx, :Lb]).to(dev)
    am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[idx])[:, None]).long().to(dev)
    return torch.where(am.bool(), t, torch.full_like(t, PAD)), am


def checks(a, m, te, tok, lens, ids, ann, cfg, torch, dev):
    """Agreement on what is kept: the mean over a cell's gene tokens (limit 1e-4) and the annotation positions (limit
    1e-3). The largest difference at any single gene position is reported beside them and is not judged."""
    k2 = a.check // 2
    order = np.argsort(-lens, kind="stable")
    idx = np.concatenate([order[:k2], order[-k2:]])
    cpu = torch.device("cpu")
    res = {"cells": int(len(idx)), "tokens_of_the_cells": [int(lens[idx].min()), int(lens[idx].max())], "device": str(dev), "limit_mean_over_genes": 1e-4, "limit_annotation_positions": 1e-3,
           "layers": len(m.encoder.layers), "d_model": m.d_model}
    judged = []

    def agree(x, y, kk, amc, perm=None):
        real = amc.bool()
        zx, zy = te.pool_hidden(x[:, kk:], amc, "gene-mean"), te.pool_hidden(y[:, kk:], amc, "gene-mean")
        out = {"mean_over_genes": float((zx - zy).abs().max()), "any_single_gene_position": float((x[:, kk:] - y[:, kk:]).abs()[real].max())}
        judged.append(out["mean_over_genes"] / 1e-4)
        if kk:
            ya = y[:, :kk] if perm is None else y[:, :kk][:, np.argsort(perm)]
            out["annotation_positions"] = float((x[:, :kk] - ya).abs().max())
            judged.append(out["annotation_positions"] / 1e-3)
        return out

    with torch.inference_mode():
        t, am = batch(tok, lens, idx, torch, dev)
        amc, tc = am.cpu(), t.cpu()
        H = {"N": forward(m, te, t, am, None, torch)[0].float().cpu(), "P": forward(m, te, t, am, ann.to(dev), torch)[0].float().cpu()}
        if a.recorded:
            rec = np.load(a.recorded)
            pos = {i: j for j, i in enumerate(rec["ids"].astype(str))}
            R = rec["L12"][[pos[i] for i in ids[idx]]]
            res["N_against_the_recorded_state"] = float(np.abs(te.pool_hidden(H["N"], amc, "gene-mean").numpy() - R).max())
            judged.append(res["N_against_the_recorded_state"] / 1e-4)
        res["largest_value"] = {"mean_over_genes": float(te.pool_hidden(H["P"][:, 4:], amc, "gene-mean").abs().max()), "annotation_positions": float(H["P"][:, :4].abs().max())}
        perm = [3, 2, 1, 0]
        hp = forward(m, te, t, am, ann[perm].to(dev), torch)[0].float().cpu()
        res["order_of_the_annotations"] = agree(H["P"], hp, 4, amc, perm)
        mc = te.load_teddy(a.ckpt, cpu)
        C = {"N": forward(mc, te, tc, amc, None, torch)[0].float(), "P": forward(mc, te, tc, amc, ann, torch)[0].float()}
        res["this_device_against_cpu"] = {k: agree(C[k], H[k], 0 if k == "N" else 4, amc) for k in ("N", "P")}
        if a.official:
            sys.path.insert(0, str(a.official))
            from safetensors.torch import load_file
            from teddy.models.teddy_g.model import TeddyGConfig, TeddyGModel
            keep = ("ntoken", "max_position_embeddings", "nlayers", "nheads", "d_model", "d_hid", "layer_activation", "dropout", "initializer_range", "n_cls", "n_layers_cls", "pre_norm")
            om = TeddyGModel(TeddyGConfig(**{k: cfg[k] for k in keep if k in cfg}))
            missing, unexpected = om.load_state_dict(load_file(str(a.ckpt / "model.safetensors"), device="cpu"), strict=False)
            assert not missing and not unexpected, (missing, unexpected)
            om.eval()
            om.return_all_embs = True
            import transformers
            res["official_class_against_this_code_on_cpu"] = {"transformers": transformers.__version__}
            for k in ("N", "P"):
                kw = {} if k == "N" else dict(annotations=ann[None].expand(len(idx), -1), annotation_attention_mask=torch.ones(len(idx), 4, dtype=torch.long))
                out = om(gene_ids=tc, attention_mask=amc, **kw)["all_embs"].float()
                res["official_class_against_this_code_on_cpu"][k] = agree(out, C[k], 0 if k == "N" else 4, amc)
    res["hold"] = bool(max(judged) <= 1.0)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--tokens", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--official", type=Path, help="a folder holding the official code (teddy/models/teddy_g/model.py)")
    ap.add_argument("--recorded", type=Path, help="layers_b.npz of runs/encoder_stages_20261124 (ids, L12); the 70M checkpoint only")
    ap.add_argument("--check", type=int, default=0)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    T = np.load(a.tokens)
    tok, lens, ids = T["tok"].astype(np.int64), T["L"].astype(np.int64), T["ids"].astype(str)
    cfg = json.loads((a.ckpt / "config.json").read_text())
    words = [line.rstrip("\n") for line in open(a.ckpt / "vocab.txt")]
    voc = {w: i for i, w in enumerate(words)}
    assert [voc[w] for w in PLACE] == [43804, 43806, 43805, 43807] and words[LABELS[0]] == "brain_disease" and words[LABELS[1] - 1] == "unknown", "another vocabulary"
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    ann = torch.tensor([voc[w] for w in PLACE])
    weights = int(sum(p.numel() for p in m.parameters()))
    info = {"layers": len(m.encoder.layers), "d_model": m.d_model, "heads": int(cfg["nheads"]), "d_hid": int(cfg["d_hid"]), "weights": weights, "torch": torch.__version__, "numpy": np.__version__, "python": sys.version.split()[0]}
    if a.check:
        res = {**checks(a, m, te, tok, lens, ids, ann, cfg, torch, dev), **info, "seconds": round(time.time() - t0, 1)}
        (a.out / "check.json").write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res, indent=1))
        sys.exit(0 if res["hold"] else 2)
    f = a.out / "states.npz"
    if f.exists():
        say("the states are there")
        return
    n, d = len(lens), m.d_model
    lo, hi = LABELS
    order = np.argsort(-lens, kind="stable")
    E0, LP, Z, ZP, C = (np.zeros((n, d), np.float32) for _ in range(5))
    top, lab = np.zeros((n, 4), np.int32), np.zeros((n, 4, hi - lo), np.float32)
    annd = ann.to(dev)
    with torch.inference_mode():
        for s in range(0, n, a.batch):
            idx = order[s:s + a.batch]
            t, am = batch(tok, lens, idx, torch, dev)
            h, e0, lp, _k = forward(m, te, t, am, None, torch)
            E0[idx], LP[idx], Z[idx] = e0.cpu().numpy(), lp.cpu().numpy(), te.pool_hidden(h.float(), am, "gene-mean").cpu().numpy()
            h, _e, _l, kk = forward(m, te, t, am, annd, torch)
            ZP[idx] = te.pool_hidden(h[:, kk:].float(), am, "gene-mean").cpu().numpy()
            ha = h[:, :kk].float()
            p = torch.softmax(m.decoder_head(ha) + m.decoder_bias, dim=-1)
            C[idx], top[idx], lab[idx] = ha[:, 0].cpu().numpy(), p.argmax(-1).cpu().numpy(), p[..., lo:hi].cpu().numpy()
            if (s // a.batch) % 40 == 0:
                say(f"{s + len(idx)}/{n} on {dev}")
    np.savez(a.out / "tmp_states.npz", ids=ids, E0=E0, LP=LP, Z=Z, ZP=ZP, C=C, top=top, lab=lab)
    (a.out / "tmp_states.npz").rename(f)
    man = {"cells": int(n), **info, "device": str(dev), "placeholders": dict(zip(["disease", "tissue", "cell_type", "sex"], PLACE)), "label_tokens": [lo, hi - 1], "seconds": round(time.time() - t0, 1)}
    if a.recorded:
        rec = np.load(a.recorded)
        assert np.array_equal(rec["ids"].astype(str), ids), "the recorded states are for other cells"
        man["N_against_the_recorded_state"] = float(np.abs(Z - rec["L12"]).max())
        man["check_passed"] = bool(man["N_against_the_recorded_state"] <= 1e-4)
    (a.out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    say("done:", json.dumps({k: man[k] for k in ("cells", "layers", "d_model", "weights", "seconds")}))


if __name__ == "__main__":
    main()
