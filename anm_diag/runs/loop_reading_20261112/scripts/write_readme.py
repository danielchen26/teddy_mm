"""Writes ../README.md from results/loop_reading.json; every number is read from that file."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = json.load(open(HERE.parent / "results" / "loop_reading.json"))
CEN = json.load(open(HERE.parent / "results" / "data_census.json"))
DS = [("frangieh2021", "melanoma"), ("papalexi2021", "THP-1")]
RO = {"PDL1": "PD-L1", "PDL2": "PD-L2", "HLAE": "HLA-E", "MHC1": "MHC-I"}
ARM = {"decoder_only": "control: decoder repaired", "mlp_only": "round 1: MLP repaired", "mlp_guarded": "round 2: MLP repaired, bone-marrow outputs held", "joint_s0": "round 3: one fit for both cell types"}
S, A, B = "Z->HEAD_OUT", "Z->HEAD_LAT", "HEAD_LAT->HEAD_OUT"
iv = lambda x: f"[{x[0]:+.2f}, {x[1]:+.2f}]"
iv3 = lambda x: f"[{x[0]:+.3f}, {x[1]:+.3f}]"
pc = lambda x: f"{100 * x:.0f}%"
pci = lambda x: f"[{100 * x[0]:.0f}%, {100 * x[1]:.0f}%]"
o = []
P = {ds: R["diagnosis"]["pooled"][ds] for ds, _ in DS}
m, t = P["frangieh2021"]["PDL1"][S], P["papalexi2021"]["PDL1"][S]
assert all(x["outcome"] == "located step confirmed" and x["paired"] for x in (m, t))
o.append(f"""# The loop on the teddy_mm head, re-read under the protocol — 12 November 2026

Post hoc and exploratory. No fit and no audit was made for this folder. It reads the records and bootstrap draws of the
three rounds ([`loop_20261112`](../loop_20261112/), [`loop_guard_20261112`](../loop_guard_20261112/),
[`repair_20261112`](../repair_20261112/)) with [`anm_diag/loop.py`](../../anm_diag/loop.py), under
[`LOOP_PROTOCOL.md`](../../LOOP_PROTOCOL.md). The outcomes those rounds recorded under their own rules stand as recorded.

Why it was made: the rounds decided by comparing a point value with 0.10. The audits of all heads of a data set use the
same cells, folds and resamples (the draws of the response are identical, checked), so the change between two stages and
the difference between two heads both have intervals that can be read from what was already computed.

The span is Z → head output: from the stage before the located step to what the head delivers.

## In short

- **The located step is confirmed, on both cell types.** Over Z → output the head as trained adds {m['baseline']['added']:.2f} (melanoma) and
  {t['baseline']['added']:.2f} (THP-1) to *u* for PD-L1. Repairing the MLP, the step the audit names, removes {pc(m['located']['share_removed']['point'])} {pci(m['located']['share_removed']['ci'])} and {pc(t['located']['share_removed']['point'])} {pci(t['located']['share_removed']['ci'])}
  of that. Repairing the decoder with the same cells, loss and budget removes {pc(m['control']['share_removed']['point'])} {pci(m['control']['share_removed']['ci'])} and {pc(t['control']['share_removed']['point'])} {pci(t['control']['share_removed']['ci'])}.
- **No head is shown to be repaired in the paper's sense.** No changed head has a supported end stage. On THP-1 the changed
  heads add no loss over the span, by their intervals ("relative" closure). On melanoma the loss is no longer resolved, and
  its absence is not shown either.
- **Round 3 is slightly worse than round 2 on melanoma**, in a direct paired comparison, and no different on THP-1.
- **HLA-E and MHC-I should not have been targets.** Their response is too small against its own interval to follow.

## 1. Admission (gate 0), on the audit cells of the loop

| cells | readout | response over the half-width of its interval | reference Z here | reference Z in the sealed run | head as trained, Z → output | half-width | diagnostic loop | closure attainable | clusters |
|---|---|---|---|---|---|---|---|---|---|""")
for ds, short in DS:
    for rid, a in R["admission"][ds].items():
        sz = R["reference_in_the_sealed_runs"][ds]["stages"].get(rid, {}).get("Z")
        base = a.get("baseline")
        o.append(f"| {short} | {RO[rid]} | {a['response']['snr']:.1f} | {a['reference']['verdict'] if 'reference' in a else '—'} | {sz['verdict'] if sz else '—'} | "
                 + (f"{base['added']:+.2f} {iv(base['ci'])} {base['reading']}" if base else "—") + f" | {a['precision']['span_half_width']:.2f} | {a['diagnostic']} | {a['certification']} | {a['clusters']['n']} |")
o.append(f"""
Reference verdicts are two-sided. The sealed runs used every control cell ({R['reference_in_the_sealed_runs']['frangieh2021']['n_cells']:,} and {R['reference_in_the_sealed_runs']['papalexi2021']['n_cells']:,} cells); the loop's audits keep half of
the controls for the fits. With 5 treated clusters every interval here is approximate, so the grade stays exploratory
whatever else holds. "Relative": a head that adds almost nothing over the span could be recognized. "Marginal": only just.

## 2. The diagnostic loop: is the MLP the step that matters?

Baseline: the head as trained. Located change: the MLP repaired (round 1). Control change: the decoder repaired with the
same cells, loss and budget. Each fit holds one target gene out, and each knockout cell is audited through the head that
never saw its gene.

| cells | readout | head as trained adds | MLP repaired: share removed | decoder repaired: share removed | MLP minus decoder, in *u* | outcome | selective |
|---|---|---|---|---|---|---|---|""")
for ds, short in DS:
    for rid, d in P[ds].items():
        x = d[S]
        o.append(f"| {short} | {RO[rid]} | {x['baseline']['added']:.2f} {iv(x['baseline']['ci'])} | {pc(x['located']['share_removed']['point'])} {pci(x['located']['share_removed']['ci'])} | "
                 f"{pc(x['control']['share_removed']['point'])} {pci(x['control']['share_removed']['ci'])} | {x['located_minus_control']['point']:+.2f} {iv(x['located_minus_control']['ci'])} | {x['outcome']} | {'yes' if x['selective'] else 'no'} |")
sm = R["steps"]["frangieh2021"]["PDL1"]
o.append(f"""
"Selective": whatever the decoder repair removes lies below the step margin. PD-L2 has no output in the head and was in no
loss; the MLP repair made for PD-L1 removes its loss as well.

**Where the loss shows and where it can be removed are not the same place** (melanoma, PD-L1). As trained the head adds a
resolved loss at both of its steps: {sm['as_trained'][A]['added']:+.2f} {iv(sm['as_trained'][A]['ci'])} from Z to the latent and {sm['as_trained'][B]['added']:+.2f} {iv(sm['as_trained'][B]['ci'])} from the latent to the output. With
the MLP repaired and the decoder untouched, the second step adds {sm['mlp_only'][B]['added']:+.2f} {iv(sm['mlp_only'][B]['ci'])}. With the decoder repaired and the MLP
untouched, it still adds {sm['decoder_only'][B]['added']:+.2f} {iv(sm['decoder_only'][B]['ci'])}. The loss at the decoder's step goes when the MLP is repaired and stays when
the decoder is.

One held-out gene at a time (PD-L1; the head fitted without that gene, audited on its knockouts):

| cells | gene | clusters | gate 0 | MLP repaired: share removed | decoder repaired: share removed | outcome |
|---|---|---|---|---|---|---|""")
n_ok = n_adm = 0
for ds, short in DS:
    for g, v in R["diagnosis"]["single_gene"][ds].items():
        d, a = v["diagnosis"], v["admission"]
        if a["diagnostic"] == "ready":
            n_adm += 1; n_ok += d["outcome"] == "located step confirmed"
            o.append(f"| {short} | {g} | {a['clusters']['n']} | ready | {pc(d['located']['share_removed']['point'])} {pci(d['located']['share_removed']['ci'])} | {pc(d['control']['share_removed']['point'])} {pci(d['control']['share_removed']['ci'])} | {d['outcome']} |")
        else:
            o.append(f"| {short} | {g} | {a['clusters']['n']} | not ready: response over half-width {a['response']['snr']:.1f} | — | — | not admitted |")
o.append(f"""
{n_ok} of the {n_adm} admitted genes read "located step confirmed". With about 200 knockout cells a share has a wide interval,
and shares above 100% mean the repaired head's *u* at the output came out below its *u* at Z.

## 3. The repair loop: is a changed head free of the loss?

Control: the head as trained, on the same cells. Grades as in the protocol. Bone marrow: the head's own task, mean
correlation over 134 proteins on 5,000 test cells, lowest to highest of the round's heads; {R['guard']['as_trained']:.3f} as trained, and the
guard allows a mean fall of {R['guard']['allowed_fall']:g}.

| cells | readout | head | added over Z → output | end stage, two-sided | closure | by step: Z → latent, latent → output | bone marrow | guard |
|---|---|---|---|---|---|---|---|---|""")
for ds, short in DS:
    for arm, lab in ARM.items():
        for rid, x in R["closure"][ds][arm].items():
            c = x[S]; g = R["guard"]["groups"].get(f"{ds}|{arm}") or R["guard"]["groups"].get(f"both|{arm}")
            o.append(f"| {short} | {RO[rid]} | {lab} | {c['arm_span']['added']:+.3f} {iv3(c['arm_span']['ci'])} | {c['arm_end']['verdict']} | **{c['grade']}** | {x[A]['grade']}, {x[B]['grade']} | "
                     f"{g['lowest']:.3f} to {g['highest']:.3f} | {'kept' if g['kept'] else 'not kept'} |")
cm = R["closure"]["frangieh2021"]
o.append(f"""
- **Absolute closure: nowhere.** The end stage of every changed head is unresolved. On these audit cells the reference
  stage Z is itself unresolved, so no later stage can be supported.
- **THP-1, PD-L1: relative closure in all three rounds.** The upper bound of the added loss is below 0.10 each time.
- **Melanoma, PD-L1: the loss is no longer resolved, in all three rounds.** Rounds 1 and 2 miss "relative" by their upper
  bounds ({cm['mlp_only']['PDL1'][S]['arm_span']['ci'][1]:.3f} and {cm['mlp_guarded']['PDL1'][S]['arm_span']['ci'][1]:.3f} against 0.10). Read one step at a time, round 2 adds no loss at either step.
- **The decoder repair closes nothing**, on either cell type.
- **Round 1 does not keep the guard.** Rounds 2 and 3 do.
""")
c3 = R["round3_against_round2"]
o.append(f"""## 4. Round 3 against round 2, on the same cells

Both are repairs of the MLP with the bone-marrow outputs held. Round 2 fits one cell type and PD-L1. Round 3 fits both
cell types and every measured protein with a head output. PD-L1, Z → output, paired:

| cells | round 2 adds | round 3 adds | round 3 minus round 2 | reading |
|---|---|---|---|---|""")
for ds, short in DS:
    c = c3[ds]["PDL1"]
    o.append(f"| {short} | {c['added_baseline']['added']:+.3f} {iv3(c['added_baseline']['ci'])} | {c['added_arm']['added']:+.3f} {iv3(c['added_arm']['ci'])} | {-c['removed']:+.3f} {iv3([-c['removed_ci'][1], -c['removed_ci'][0]])} | "
             + ("round 3 adds more" if c["reading"] == "adds to the loss" else ("no difference resolved" if c["reading"] in ("does not reduce the loss", "unresolved") else c["reading"])) + " |")
sp = R["round3_splits"]
o.append(f"""
On melanoma the joint head adds more than the per-cell-type head, by an interval that just excludes zero. The rule of
round 3 called this "regressed" from a point value 0.009 over the margin; the direct comparison supports the direction
of that reading and gives its size.

Round 3 on its three splits of the control cells (PD-L1, Z → output; no head as trained was audited on splits 1 and 2):

| cells | split 0 | split 1 | split 2 |
|---|---|---|---|""")
for ds, short in DS:
    o.append(f"| {short} | " + " | ".join(f"{sp[ds][s]['PDL1'][S]['added']:+.3f} {iv3(sp[ds][s]['PDL1'][S]['ci'])} {sp[ds][s]['PDL1'][S]['reading']}" for s in ("0", "1", "2")) + " |")
tr = R["transfer"]
o.append(f"""
## 5. A head repaired on the other cell type (round 1)

Audited on every cell of the data set, since the head never saw any of them. PD-L1, Z → output:

| audited on | head repaired on | added | reading | reference Z | end stage |
|---|---|---|---|---|---|
| melanoma ({tr['frangieh2021']['mlp_only']['n_cells']:,} cells) | THP-1, MLP | {tr['frangieh2021']['mlp_only']['PDL1'][S]['added']:+.2f} {iv(tr['frangieh2021']['mlp_only']['PDL1'][S]['ci'])} | {tr['frangieh2021']['mlp_only']['PDL1'][S]['reading']} | {tr['frangieh2021']['mlp_only']['PDL1']['reference']['verdict']} | {tr['frangieh2021']['mlp_only']['PDL1']['end']['verdict']} |
| THP-1 ({tr['papalexi2021']['mlp_only']['n_cells']:,} cells) | melanoma, MLP | {tr['papalexi2021']['mlp_only']['PDL1'][S]['added']:+.2f} {iv(tr['papalexi2021']['mlp_only']['PDL1'][S]['ci'])} | {tr['papalexi2021']['mlp_only']['PDL1'][S]['reading']} | {tr['papalexi2021']['mlp_only']['PDL1']['reference']['verdict']} | {tr['papalexi2021']['mlp_only']['PDL1']['end']['verdict']} |

A repair does not carry to the other cell type. With every control cell in the audit, the reference stage on melanoma is
supported: the loop's own audits lose that by keeping half of the controls for the fits.

## 6. What data are left: a census of the protein files

[`tools/data_census.py`](../../tools/data_census.py) reads only the protein files and lists, per condition, the knockouts that
change a measured protein clearly (at least {CEN['frangieh2021']['thresholds']['min_cells']} single-guide cells, a shift of at least {CEN['frangieh2021']['thresholds']['min_shift_sd']:g} control standard deviations, and
at least {CEN['frangieh2021']['thresholds']['min_shift_over_half_width']:g} times the half-width of its interval, with cells as units, so the ratio is optimistic). No model is involved
(`results/data_census.json`).

| data | cells | control-guide cells | target genes with enough cells | clear (gene, protein) responses | genes with one | genes that change PD-L1 |
|---|---|---|---|---|---|---|""")
for name, ds, conds in (("melanoma", "frangieh2021", ("IFNγ", "Co-culture", "Control")), ("THP-1", "papalexi2021", ("all",))):
    for c in conds:
        v = CEN[ds]["conditions"][c]; pd = [x for x in v["clear_responses"] if x["antibody"] in ("CD274", "PDL1")]
        o.append(f"| {name}" + ("" if c == "all" else f", {'IFN-γ' if c == 'IFNγ' else c.lower()}") + f" | {v['cells']:,} | {v['control_guide_cells']:,} | {v['targets_with_min_cells']} | {len(v['clear_responses'])} | {len(v['targets_with_a_clear_response'])} | "
                 + ", ".join(f"{x['target']} ({x['shift_sd']:+.2f})" for x in sorted(pd, key=lambda x: x["shift_sd"])) + " |")
co = CEN["frangieh2021"]["conditions"]["Co-culture"]; five = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
cop = {x["target"]: x for x in co["clear_responses"] if x["antibody"] == "CD274"}; coh = {x["target"]: x for x in co["clear_responses"] if x["antibody"] == "HLA_A"}
th = {x["target"]: x for x in CEN["papalexi2021"]["conditions"]["all"]["clear_responses"] if x["antibody"] == "PDL1"}; new = [g for g in th if g not in ("IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1")]
assert all(g in cop and g in coh for g in five)
o.append(f"""
The three rounds used the IFN-γ condition of the melanoma data and five genes of the THP-1 data. What they did not use:

- **Melanoma, co-culture.** The same five genes change PD-L1 there by {max(cop[g]['shift_sd'] for g in five):+.2f} to {min(cop[g]['shift_sd'] for g in five):+.2f} standard deviations, on {sum(cop[g]['cells'] for g in five):,} single-guide
  knockout cells, with {co['control_guide_cells']:,} control-guide cells. No fit and no audit of the loop saw these cells: they could serve as
  final data for the melanoma head (a new condition, the same cell type). Results of other heads on co-culture cells,
  from work not recorded in this repository, were known when the loop was designed (design of round 1). MHC-I, too weak
  to follow under IFN-γ, moves by {max(coh[g]['shift_sd'] for g in five):+.2f} to {min(coh[g]['shift_sd'] for g in five):+.2f} there.
- **THP-1, other genes.** {', '.join(f"{g} ({th[g]['shift_sd']:+.2f})" for g in new)} change PD-L1 clearly and were in no fit or audit: new source genes for a final audit.
  TEDDY's embedding would have to be computed for their cells.
- **More clusters.** For PD-L1 no data set has ten source genes ({len(CEN['frangieh2021']['conditions']['IFNγ']['targets_with_a_clear_response']) and sum(1 for x in CEN['frangieh2021']['conditions']['IFNγ']['clear_responses'] if x['antibody'] == 'CD274')} in melanoma, {len(th)} in THP-1). Gene clusters stay below the
  protocol's ten; guides, three or four per gene, give more units.

The census says where a response exists in the measurement. Whether TEDDY's embedding carries it (the reference stage)
is a question for an audit.

## What this changes in how the rounds are read

| the rounds' reports say | under the protocol |
|---|---|
| rounds 1 and 2: "closed" on both cell types | the located step is confirmed on both; the head is shown free of added loss on THP-1 only, and only relative to an unresolved reference |
| round 3: "regressed" on melanoma by 0.009 | on melanoma the joint head adds {-c3['frangieh2021']['PDL1']['removed']:.3f} {iv3([-c3['frangieh2021']['PDL1']['removed_ci'][1], -c3['frangieh2021']['PDL1']['removed_ci'][0]])} more than the round-2 head; on THP-1 no difference |
| round 3: HLA-E "still located in the head" | HLA-E does not pass gate 0 on these cells; it was not a target that could be followed |
| "no loss in the head" as a pass | three different things: no added loss shown (THP-1), loss no longer resolved (melanoma), and never a supported head |

## Limits

- Post hoc: the protocol and this reading were written after every result here was known.
- 5 treated clusters in the pooled audits and 4 in the THP-1 single-gene audits: intervals are approximate.
- The same audit cells served all three rounds. Nothing here is a final audit.
- Changes are fine-tunes of a trained head with one training seed, compared with the head as trained.
- Changes in *u* are taken in |*u*|. For PD-L1 and PD-L2 on these cells that is the same as the signed change.

## Files

`scripts/read_loop.py` reads the three rounds and writes `results/loop_reading.json` (every reading, with intervals) and
`results/ledger.jsonl` (the changes, audits and decisions as typed events). `scripts/write_readme.py` writes this file.
The bootstrap draws (`anm_bootstrap.npz` in each audit folder) are not kept in git; re-running an audit writes them again.
""")
(HERE.parent / "README.md").write_text("\n".join(o))
print("written", sum(len(x) for x in o))
