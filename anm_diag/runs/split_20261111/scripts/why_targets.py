"""Why does the share lost in the head differ so much between targets? Exploratory, on the exported features and the head's weights.
For every target: how far the knockout moves Z, how strongly the head responds when control cells are moved by that same shift
(set against random directions of equal length), how well the shift lines up with the direction the head's PD-L1 output is most
sensitive to, and how broad the knockout's footprint is in the RNA."""
import json, os, numpy as np, scipy.sparse as sp
from scipy.special import expit
from _paths import FEATURES, RES
F = str(FEATURES) + "/"
W = dict(np.load(F + "head_weights.npz")); names = [str(x) for x in W["adt_names"]]; iPD = names.index("CD274")
silu = lambda a: a * expit(a)
def mlp(p, x):
    h = silu(x @ W[p + "net.0.weight"].T + W[p + "net.0.bias"]); h = silu(h @ W[p + "net.3.weight"].T + W[p + "net.3.bias"])
    return h @ W[p + "net.6.weight"].T + W[p + "net.6.bias"]
def head(z):
    z = np.asarray(z, np.float64); lat = mlp("mlp__", z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6))
    return lat, np.log(np.maximum(np.logaddexp(0.0, mlp("dec__mu.", lat)), 1e-6))
ISG = ["IRF1", "GBP1", "GBP2", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "STAT1", "TAP1", "PSMB9", "PSMB8", "HLA-DRA", "CIITA", "SOCS1", "ICAM1", "APOL6", "WARS", "PARP14", "IFI30", "SECTM1", "CD74"]
rng = np.random.default_rng(7); NR = 200
def load(sub):
    d = np.load(F + sub + "/dense.npz"); D = {k: d[k] for k in d.files}
    for k in ("Z", "HEAD_LAT", "HEAD_OUT"):
        if "node__" + k not in D: D["node__" + k] = np.load(F + sub + f"/node__{k}.npy")
    D["RAW"] = sp.load_npz(F + sub + "/node__RAW.npz").toarray().astype(np.float32); D["genes"] = json.load(open(F + sub + "/raw_genes.json"))
    return D
def analyse(sub, targets, readout, u_out, u_z, delta):
    D = load(sub); A = D["A"]; tg = D["target"].astype(str); Z = D["node__Z"].astype(np.float64); C = np.where(A == 0)[0]; y = D["readout__" + readout]
    lat0, out0 = head(Z); chk = (float(np.abs(lat0 - D["node__HEAD_LAT"]).max()), float(np.abs(out0 - D["node__HEAD_OUT"]).max()))
    Zc = Z[C]; latc, outc = lat0[C], out0[C]; G = D["genes"]; gi = {g: i for i, g in enumerate(G)}; raw = D["RAW"]; mc, sc = raw[C].mean(0), raw[C].std(0)
    sub_i = rng.choice(len(C), min(600, len(C)), replace=False); Zs = Zc[sub_i]; h = 0.5; grad = np.zeros(Z.shape[1])
    for j in range(Z.shape[1]):                                  # gradient of the mean PD-L1 output with respect to Z (central differences)
        e = np.zeros(Z.shape[1]); e[j] = h; grad[j] = (head(Zs + e)[1][:, iPD].mean() - head(Zs - e)[1][:, iPD].mean()) / (2 * h)
    res = {"dataset": sub, "check_head_reimplementation_maxdiff": chk, "n_control": int(len(C)), "grad_norm": float(np.linalg.norm(grad)), "targets": {}}
    dirs = {}
    for g in targets:
        T = np.where((A == 1) & (tg == g))[0]; d = Z[T].mean(0) - Zc.mean(0); n = float(np.linalg.norm(d)); dirs[g] = d / n
        lat1, out1 = head(Zc + d); r_pd = float((out1[:, iPD] - outc[:, iPD]).mean()); dl = float(np.linalg.norm(lat1 - latc, axis=1).mean()); do = float(np.linalg.norm(out1 - outc, axis=1).mean())
        rp, rl, ro = [], [], []
        for _ in range(NR):
            v = rng.standard_normal(Z.shape[1]); v *= n / np.linalg.norm(v); l2, o2 = head(Zs + v)
            rp.append(float((o2[:, iPD] - outc[sub_i][:, iPD]).mean())); rl.append(float(np.linalg.norm(l2 - latc[sub_i], axis=1).mean())); ro.append(float(np.linalg.norm(o2 - outc[sub_i], axis=1).mean()))
        own = float((out0[T][:, iPD].mean() - outc[:, iPD].mean()))      # the head's own predicted change, treated minus control
        cd = (raw[T].mean(0) - mc) / np.where(sc > 1e-6, sc, 1.0); isg = [cd[gi[x]] for x in ISG if x in gi]
        res["targets"][g] = dict(n_treated=int(len(T)), measured_delta=float(delta[g]), u_head_out=float(u_out[g]), u_z=float(u_z[g]), z_shift=n,
            pd_move=r_pd, pd_move_pct_random=float(100 * np.mean(np.abs(rp) < abs(r_pd))), pd_random_abs_median=float(np.median(np.abs(rp))),
            lat_move=dl, lat_vs_random=dl / float(np.median(rl)), out_move=do, out_vs_random=do / float(np.median(ro)),
            head_own_change=own, head_own_share=own / float(delta[g]), move_share=r_pd / float(delta[g]),
            cos_grad=float(d @ grad / (n * np.linalg.norm(grad))), linear_pred=float(d @ grad),
            rna_genes_d025=int((np.abs(cd) > 0.25).sum()), rna_genes_d05=int((np.abs(cd) > 0.5).sum()), own_gene_d=(float(cd[gi[g]]) if g in gi else None),
            cd274_rna_d=(float(cd[gi["CD274"]]) if "CD274" in gi else None), isg_mean_d=float(np.mean(isg)), isg_n=len(isg), isg_down=int(sum(x < -0.25 for x in isg)))
    res["cos_between_targets"] = {a: {b: float(dirs[a] @ dirs[b]) for b in targets} for a in targets}
    return res
sp_ = json.load(open(RES / "papalexi_split.json")); pt = json.load(open(RES / "frangieh_split.json"))
e = sp_["targets"]["PDL1"]; GT = sp_["targets_done"]
thp = analyse("papalexi2021", GT, "PDL1", {g: e["stages"]["HEAD_OUT"]["sites"][g]["u"] for g in GT}, {g: e["stages"]["Z"]["sites"][g]["u"] for g in GT}, {g: e["delta"][g]["v"] for g in GT})
e = pt["readouts"]["PDL1"]["evidence"]; GM = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
mel = analyse("frangieh2021", GM, "PDL1", {g: e["stages"]["HEAD_OUT"]["sites"][g]["u"] for g in GM}, {g: e["stages"]["Z"]["sites"][g]["u"] for g in GM}, {g: e["delta"][g]["v"] for g in GM})
json.dump(dict(thp=thp, mel=mel, isg_list=ISG, n_random=NR), open(RES / "why_targets.json", "w"), indent=1)
for R in (thp, mel):
    print("==", R["dataset"], "head check", [f"{x:.1e}" for x in R["check_head_reimplementation_maxdiff"]], "grad norm", round(R["grad_norm"], 4))
    for g, t in R["targets"].items():
        print(f"  {g:7s} Δ{t['measured_delta']:+.2f} u_out {t['u_head_out']:.2f} u_z {t['u_z']:+.2f} | |dZ| {t['z_shift']:.2f} | lat×rand {t['lat_vs_random']:.2f} out×rand {t['out_vs_random']:.2f} | PD-L1 move {t['pd_move']:+.3f} ({100*t['move_share']:.0f}% of Δ, pct {t['pd_move_pct_random']:.0f}) own {t['head_own_change']:+.3f} | cos∇ {t['cos_grad']:+.3f} | RNA genes>.25 {t['rna_genes_d025']} >.5 {t['rna_genes_d05']} own-gene d {t['own_gene_d'] if t['own_gene_d'] is None else round(t['own_gene_d'],2)} CD274 d {round(t['cd274_rna_d'],2) if t['cd274_rna_d'] is not None else None} ISG mean d {t['isg_mean_d']:+.2f} ({t['isg_down']}/{t['isg_n']} down)")
    print("  cos between targets:", {a: {b: round(v, 2) for b, v in r.items()} for a, r in R["cos_between_targets"].items()})
