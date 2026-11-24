#!/usr/bin/env python3
"""The state of the encoder with and without TEDDY's annotation positions, on the owner's machine (TEDDY-G 70M).

  python annot.py --repo TEDDY_MM --ckpt TEDDY_G_70M --tokens TOKENS.npz --recorded LAYERS_B.npz --out OUT_DIR
                  [--official OFFICIAL_CODE] [--check 16]

Reads tokens only (TOKENS.npz of runs/encoder_stages_20261124/scripts/layers.py): no protein value and no label of an
event is opened here. The forward pass follows the official code (teddy/models/teddy_g/model.py, `forward`):

  genes        gene embedding + position embedding, positions 0, 1, ... in rank order
  annotations  four tokens placed BEFORE the genes, with no position embedding, never masked from attention
               (model.py: `torch.cat([annotation_embeddings, embeddings], dim=1)`), in the order of the official
               tokeniser: disease, tissue, cell type, sex
  layers       the twelve layers in turn, padding masked

Ways of running it (MODES). Every cell of the stratum gets the same four tokens:

  N  no annotation positions at all: how the repository and the official tutorial run the model
  P  the four placeholders <disease> <tissue> <cell> <sex>: what the paper says is put in front of the genes
  T  the labels of the THP-1 line by the official mapping files: cancer_disease, hematopoietic_tissue, precursor_cell, male
  D  what the official tokeniser writes when a study has no annotation columns: healthy, unknown, unknown, unknown
  W  labels that are wrong for these cells: brain_disease, central_nervous_tissue, neural_cell, female

Kept for every cell and mode (OUT_DIR/mode_<m>.npz):

  Z     mean over the cell's gene tokens of the output of layer 12 (the state the head reads; in N it is the record)
  L11   the same after layer 11
  A     the output of layer 12 at the four annotation positions (4 x 512)
  top   the token the model's own output layer ranks first at each annotation position (all 43,840 tokens)
  lab   its probability for each of the 43 label tokens at each annotation position, and their sum

--check K: only the checks, on K cells (the K/2 longest and the K/2 shortest). Agreement numbers are printed and
nothing else: no state, no label and no prediction is shown. With --official, the official class itself is loaded
(it needs the `transformers` package, a version below 5) and its output is compared position by position.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import importlib.util

import numpy as np

PAD = 43810
ORDER = ["disease", "tissue", "cell_type", "sex"]
MODES = {"N": None,
         "P": ["<disease>", "<tissue>", "<cell>", "<sex>"],
         "T": ["cancer_disease", "hematopoietic_tissue", "precursor_cell", "male"],
         "D": ["healthy", "unknown", "unknown", "unknown"],
         "W": ["brain_disease", "central_nervous_tissue", "neural_cell", "female"]}
LABELS = (43761, 43804)                                          # label tokens: 11 diseases, 17 tissues, 12 kinds of cell, female, male, unknown
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def forward(m, te, t, am, ann, torch):
    """t [B, L] gene tokens (padding already PAD), am [B, L] 1 where a token is. ann: None or the four token ids.
    Returns the output of layer 12 [B, k + L, d], the gene-mean after layer 11, and k."""
    dev = t.device
    h = m.embeddings(t) + m.position_embeddings(torch.arange(t.size(1), device=dev))
    pad = ~am.bool()
    k = 0
    if ann is not None:
        k = int(ann.numel())
        h = torch.cat([m.embeddings(ann)[None].expand(t.size(0), -1, -1), h], dim=1)
        pad = torch.cat([torch.zeros(t.size(0), k, dtype=torch.bool, device=dev), pad], dim=1)
    l11 = None
    for i, layer in enumerate(m.encoder.layers):
        h = layer(src=h, src_key_padding_mask=pad)
        if i == 10:
            l11 = te.pool_hidden(h[:, k:].float(), am, "gene-mean")
    return h, l11, k


def batch(tok, lens, idx, torch, dev):
    Lb = int(lens[idx].max())
    t = torch.from_numpy(tok[idx, :Lb]).to(dev)
    am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[idx])[:, None]).long().to(dev)
    return torch.where(am.bool(), t, torch.full_like(t, PAD)), am


def checks(a, m, te, tok, lens, ids, voc, torch, dev):
    """Agreement on what is kept: the mean over a cell's gene tokens (limit 1e-4, as for the recorded states) and the
    four annotation positions (single positions; limit 1e-3). The largest difference at any single gene position is
    reported beside them and is not judged."""
    k2 = a.check // 2
    order = np.argsort(-lens, kind="stable")
    idx = np.concatenate([order[:k2], order[-k2:]])
    cpu = torch.device("cpu")
    res = {"cells": int(len(idx)), "tokens_of_the_cells": [int(lens[idx].min()), int(lens[idx].max())], "device": str(dev), "limit_mean_over_genes": 1e-4, "limit_annotation_positions": 1e-3}
    rec = np.load(a.recorded)
    pos = {i: j for j, i in enumerate(rec["ids"].astype(str))}
    R = rec["L12"][[pos[i] for i in ids[idx]]]
    ann = {k: (None if v is None else torch.tensor([voc[x] for x in v])) for k, v in MODES.items()}
    judged = []

    def agree(x, y, kk, amc, perm=None):
        """x, y [B, kk + L, d] on the CPU. perm: y holds the annotation positions in this order."""
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
        amc = am.cpu()
        H = {}
        for k, v in ann.items():
            h, _l, kk = forward(m, te, t, am, None if v is None else v.to(dev), torch)
            H[k] = h.float().cpu()
        res["N_against_the_recorded_state"] = float(np.abs(te.pool_hidden(H["N"], amc, "gene-mean").numpy() - R).max())
        judged.append(res["N_against_the_recorded_state"] / 1e-4)
        res["largest_value"] = {"mean_over_genes": float(te.pool_hidden(H["P"][:, 4:], amc, "gene-mean").abs().max()), "annotation_positions": float(H["P"][:, :4].abs().max()),
                                "single_gene_position": float(H["P"][:, 4:].abs()[amc.bool()].max())}
        # the order of the four annotation tokens: no position is added to them, so it must not matter
        perm = [3, 2, 1, 0]
        hp, _l, _k = forward(m, te, t, am, ann["P"][perm].to(dev), torch)
        res["order_of_the_annotations"] = agree(H["P"], hp.float().cpu(), 4, amc, perm)
        # the same on the CPU (the device the official tutorial runs on)
        mc = te.load_teddy(a.ckpt, cpu)
        tc = t.cpu()
        C = {}
        for k in ("N", "P", "T"):
            C[k] = forward(mc, te, tc, amc, ann[k], torch)[0].float()
            res.setdefault("this_device_against_cpu", {})[k] = agree(C[k], H[k], 0 if ann[k] is None else 4, amc)
        if a.official:
            sys.path.insert(0, str(a.official))
            from safetensors.torch import load_file
            from teddy.models.teddy_g.model import TeddyGConfig, TeddyGModel
            cfg = json.loads((a.ckpt / "config.json").read_text())
            keep = ("ntoken", "max_position_embeddings", "nlayers", "nheads", "d_model", "d_hid", "layer_activation", "dropout", "initializer_range", "n_cls", "n_layers_cls")
            om = TeddyGModel(TeddyGConfig(**{k: cfg[k] for k in keep if k in cfg}))
            missing, unexpected = om.load_state_dict(load_file(str(a.ckpt / "model.safetensors"), device="cpu"), strict=False)
            assert not missing and not unexpected, (missing, unexpected)
            om.eval()
            om.return_all_embs = True
            import transformers
            res["official_class_against_this_code_on_cpu"] = {"transformers": transformers.__version__}
            for k in ("N", "P", "T"):
                kw = {} if ann[k] is None else dict(annotations=ann[k][None].expand(len(idx), -1), annotation_attention_mask=torch.ones(len(idx), 4, dtype=torch.long))
                out = om(gene_ids=tc, attention_mask=amc, **kw)["all_embs"].float()
                res["official_class_against_this_code_on_cpu"][k] = agree(out, C[k], 0 if ann[k] is None else 4, amc)
    res["hold"] = bool(max(judged) <= 1.0)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--tokens", type=Path, required=True)
    ap.add_argument("--recorded", type=Path, required=True, help="layers_b.npz of the run encoder_stages_20261124 (ids, L12)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--official", type=Path, help="a folder holding the official code (teddy/models/teddy_g/model.py)")
    ap.add_argument("--check", type=int, default=0)
    ap.add_argument("--modes", default="N,P,T,D,W")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    T = np.load(a.tokens)
    tok, lens, ids = T["tok"].astype(np.int64), T["L"].astype(np.int64), T["ids"].astype(str)
    words = [line.rstrip("\n") for line in open(a.ckpt / "vocab.txt")]
    voc = {w: i for i, w in enumerate(words)}
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    versions = {"torch": torch.__version__, "numpy": np.__version__, "python": sys.version.split()[0]}
    if a.check:
        res = {**checks(a, m, te, tok, lens, ids, voc, torch, dev), **versions, "seconds": round(time.time() - t0, 1)}
        (a.out / "check.json").write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res, indent=1))
        sys.exit(0 if res["hold"] else 2)
    modes = [x for x in a.modes.split(",") if x]
    n, d = len(lens), m.d_model
    lo, hi = LABELS
    order = np.argsort(-lens, kind="stable")
    man = {"cells": int(n), "modes": {k: (None if MODES[k] is None else dict(zip(ORDER, MODES[k]))) for k in modes}, "token_ids": {k: (None if MODES[k] is None else [voc[x] for x in MODES[k]]) for k in modes},
           "label_tokens": [lo, hi - 1], "device": str(dev), **versions}
    for k in modes:
        f = a.out / f"mode_{k}.npz"
        if f.exists():
            say("mode", k, "is there")
            continue
        ann = None if MODES[k] is None else torch.tensor([voc[x] for x in MODES[k]], device=dev)
        Z, L11 = np.zeros((n, d), np.float32), np.zeros((n, d), np.float32)
        A, top, lab = (np.zeros((n, 4, d), np.float32), np.zeros((n, 4), np.int32), np.zeros((n, 4, hi - lo), np.float32)) if ann is not None else (None, None, None)
        with torch.inference_mode():
            for s in range(0, n, 16):
                idx = order[s:s + 16]
                t, am = batch(tok, lens, idx, torch, dev)
                h, l11, kk = forward(m, te, t, am, ann, torch)
                Z[idx] = te.pool_hidden(h[:, kk:].float(), am, "gene-mean").cpu().numpy()
                L11[idx] = l11.cpu().numpy()
                if ann is not None:
                    ha = h[:, :kk].float()
                    p = torch.softmax(m.decoder_head(ha) + m.decoder_bias, dim=-1)
                    A[idx], top[idx], lab[idx] = ha.cpu().numpy(), p.argmax(-1).cpu().numpy(), p[..., lo:hi].cpu().numpy()
                if (s // 16) % 60 == 0:
                    say(f"mode {k}: {s + len(idx)}/{n} on {dev}")
        keep = dict(ids=ids, Z=Z, L11=L11) if ann is None else dict(ids=ids, Z=Z, L11=L11, A=A, top=top, lab=lab)
        np.savez(a.out / f"tmp_{k}.npz", **keep)
        (a.out / f"tmp_{k}.npz").rename(f)
        say("mode", k, "done")
    rec = np.load(a.recorded)
    assert np.array_equal(rec["ids"].astype(str), ids), "the recorded states are for other cells"
    if "N" in modes:
        man["N_against_the_recorded_state"] = float(np.abs(np.load(a.out / "mode_N.npz")["Z"] - rec["L12"]).max())
        man["check_passed"] = bool(man["N_against_the_recorded_state"] <= 1e-4)
    man["seconds"] = round(time.time() - t0, 1)
    (a.out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    say("done:", json.dumps({k: man[k] for k in ("cells", "N_against_the_recorded_state", "check_passed", "seconds") if k in man}))
    if man.get("check_passed") is False:
        sys.exit(2)


if __name__ == "__main__":
    main()
