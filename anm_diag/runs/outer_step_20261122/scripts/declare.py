#!/usr/bin/env python3
"""The first outer step, declared before it is run: what package 2 changes against package 1 (the round of 20 November),
which open items each change answers, and a hypothesis with a forecast for every action (OUTER_LOOP.md, Sections 9, 10).

  python declare.py            (writes results/declared.json and results/declared.md; reads no cell)

Item ids are those of runs/outer_loop_design_20261120/results/open_items.json. Units are written gene@group.
A forecast is an interval for the estimate, and where an outcome has classes, a probability for each.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parents[1]))
from anm_diag import outer  # noqa: E402

G5 = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
T5 = ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]
U = {"mel_ifn": [f"{g}@mel_ifn" for g in G5], "mel_co": [f"{g}@mel_co" for g in G5], "thp1": [f"{g}@thp1" for g in T5]}
ALL = U["mel_ifn"] + U["mel_co"] + U["thp1"]
NEW = ["B2M@mel_ifn", "TAPBP@mel_ifn"]
SAME = "read again on the units of the proposal: exploratory, and the item can only become 'decided on development units'"

H = [
    outer.hypothesis("P1", ["O14"], "rule",
                     {"quantity": "|u of the repaired output| - |u of the reference stage| on the cells the repair was fitted on, and its interval",
                      "PDL1 mel_ifn": [-0.02, 0.04], "PDL1 mel_co": [-0.13, -0.07], "PDL1 thp1": [-0.11, -0.03], "MHC1 mel_co": [-0.07, -0.01],
                      "passes": {"each of the four": 0.95}},
                     "gate 3 passes when the interval lies below the step margin (0.10)",
                     "the repair does not reach the reference stage even on its own cells: a failed fit, and the target stays 'not graded'",
                     "the output as trained fails the same gate on the same cells (forecast: on all four)",
                     proposal_units=U["mel_co"], note=SAME),
    outer.hypothesis("P2", ["O15"], "rule",
                     {"quantity": "outcome of the test of the repair by rows with two controls (the repair of another group where there is one; the readout on control cells without rows)",
                      "PDL1 mel_ifn": {"unresolved": 0.85, "refuted": 0.10, "confirmed": 0.05}, "PDL1 mel_co": {"unresolved": 0.65, "confirmed": 0.35},
                      "PDL1 thp1": {"unresolved": 0.85, "confirmed": 0.10, "refuted": 0.05}, "MHC1 mel_co": {"unresolved": 0.70, "confirmed": 0.30},
                      "gain of rows minus gain of controls-only": {"PDL1 mel_ifn": [0.05, 0.10], "PDL1 mel_co": [0.13, 0.19], "PDL1 thp1": [0.30, 0.42], "MHC1 mel_co": [0.24, 0.31]},
                      "that difference lies above the step margin": {"PDL1 mel_ifn": 0.05, "PDL1 mel_co": 0.85, "PDL1 thp1": 0.97, "MHC1 mel_co": 0.99}},
                     "loop.test: confirmed if the rows arm's gain lies above the grid's ceiling minus the step margin and every control arm's below it",
                     "a readout fitted on the group's control cells already removes the loss: what was missing is a readout of cells of this type, and the rows add little",
                     "this amendment is the control; the difference rows minus controls-only is read on paired resamples",
                     proposal_units=U["mel_co"], note=SAME),
    outer.hypothesis("P3", ["O20"], "rule",
                     {"quantity": "status of the two groups that pass admission as a pool with no admitted event", "MHC1 mel_ifn (five pathway genes)": "not a target", "HLAE mel_co": "not a target"},
                     "a repair by rows needs an admitted row; without one the group is not a target and the reason is an open item (small response)",
                     "more cells per event would let single events pass: then it is a question of resolution (E), not of the rule",
                     "the pooled response passes A1 while no single event does (the round, replication run)",
                     proposal_units=U["mel_ifn"] + U["mel_co"], note="a definition: no reading changes"),
    outer.hypothesis("P4", ["O19"], "rule",
                     {"quantity": "the events admitted in the round, screened again on a random third of their cells: |change| / half-width",
                      "every event of PDL1 mel_ifn at least 3": 0.80, "every event of PDL1 mel_co at least 3": 0.97, "IRF1@thp1 at least 3": 0.40, "the other four of thp1 at least 3": 0.97,
                      "library knockouts other than the readout's own gene and the five that pass by chance (PDL1, of about 200)": [0, 1]},
                     "from this step on an event is admitted on its selection third only; the events of the round keep their admission and this screen is reported",
                     "the selection made on all cells favoured no event: every event passes on a third as well",
                     "the same screen on the library's other knockouts: how many pass on a third",
                     proposal_units=ALL, note=SAME),
    outer.hypothesis("S1", ["O18", "O02"], "rule",
                     {"quantity": "u of the repaired output, each cell read by the readout of its recorded type",
                      "PDL1 mel_ifn": [-0.01, 0.05], "PDL1 mel_ifn verdict": {"support": 0.60, "unresolved": 0.40},
                      "PDL1 mel_co": [0.10, 0.13], "PDL1 thp1 (five events pooled)": [0.14, 0.17],
                      "PDL1 mel_ifn status closed (absolute) on development cells": 0.55},
                     "the type is taken from the record whatever this reading shows: the reader of types gave 96% of the knockout cells of one group their own type, below the chosen 99%",
                     "the 0.04 that was left is not from mistyped cells: u does not move",
                     "the reader of types as a component: the share of held-out cells it gives their own type, per group",
                     proposal_units=U["mel_ifn"], note=SAME),
    outer.hypothesis("S2", ["O09", "O08", "O11"], "split",
                     {"quantity": "the THP-1 pool read as two types: A = IFNGR1, IFNGR2, JAK2, STAT1; B = IRF1",
                      "A reference (other genes of the type, at Z)": [0.02, 0.12], "A reference verdict": {"unresolved": 0.65, "support": 0.20, "reject": 0.15},
                      "A repaired output": [0.06, 0.15], "A repaired output verdict": {"unresolved": 0.70, "support": 0.10, "reject": 0.20}, "A adds no loss beyond its reference": 0.80,
                      "B repaired output (its own row)": [0.35, 0.55], "B verdict reject": 0.90, "B adds no loss beyond a reader of Z with its own cells": 0.60,
                      "B from the RNA view to Z, own cells": [0.10, 0.35], "B reading upstream (interval above the step margin)": 0.45},
                     "operated provisionally. Adopted only by its test on units not used for the proposal: within each type the repair with the type's rows adds no loss, "
                     "and the repair with the pooled rows adds a resolved loss on the same units. THP-1 has no such unit with states yet, so the test is not run in this step",
                     "the missing field is not the one proposed (the place of the gene in the pathway: signalling to STAT1, or the transcription factor it induces). Rivals from the record: "
                     "the replicate (IRF1 cells are 56% from replicate 1, the others 50%, controls 40%) and the size of the response, which is the endpoint and may not type an event",
                     "a split of the same sizes at random does not pass: no other single event of the pool lies outside the pool's interval (IFNGR1 0.10, IFNGR2 0.14, JAK2 0.17, STAT1 0.02; pool 0.16 [0.08, 0.29])",
                     proposal_units=U["thp1"], note="operated in the bridge reading beside the coarser schema; " + SAME),
    outer.hypothesis("V1", ["O05", "O06"], "view",
                     {"quantity": "PD-L1 in co-culture: u with the event's own cells at the input view, per reader; the gap from the best reader of the input to Z",
                      "RNA, 512 components, boosting": [0.15, 0.30], "RNA, 4,000 genes, ridge": [0.08, 0.20],
                      "the input view is no longer dominated (the gap's interval reaches zero)": 0.60, "some event or the pool is supported by some reader at some stage": 0.15},
                     "the input view can be read when its gap to Z no longer lies below zero and the positive control holds; until both, no reading of where the loss is",
                     "Z carries what no reader of RNA of these classes finds: the view stays dominated, and Z bounds it",
                     "the same view with its rows permuted across cells reads nothing (u at least 0.9)",
                     proposal_units=U["mel_co"], note=SAME),
    outer.hypothesis("V2", ["O13"], "view",
                     {"quantity": "MHC-I in co-culture: u with the event's own cells, per stage and reader (pool of five)",
                      "Z, boosting": [0.20, 0.40], "RNA, 512 components, boosting": [0.22, 0.42], "RNA, 4,000 genes, ridge": [0.22, 0.36],
                      "every reader at every stage rejects": 0.80, "some event or the pool is supported by some reader at some stage": 0.10,
                      "boosting reads at least 0.10 better than ridge at the same stage": 0.20},
                     "`unreadable` is drawn only if the positive control holds, two reader classes of which one is not linear reject at the input view with the event's own cells, "
                     "and the view is not dominated. Without a positive control the reading stays invalid and the next action is the owner's: a tolerance for this target",
                     "the relation between the readout and the state is not linear across the knockout (within the control cells a small change, across the knockout a large one): then boosting reads it",
                     "the same view with its rows permuted across cells reads nothing (u at least 0.9)",
                     proposal_units=U["mel_co"], note=SAME),
    outer.hypothesis("E1", ["O12"], "resolution",
                     {"quantity": "which knockouts of the IFN-gamma library are admitted for MHC-I on a selection third (A1: resolved and at least 3 half-widths; A2: not HLA-A, HLA-B, HLA-C)",
                      "B2M and TAPBP both admitted": 0.95, "exactly those two": 0.70, "SSR2, STAT3, CD44: each": 0.15, "any other": 0.05},
                     "an event is admitted on its selection third; those cells, and a third of the control cells, enter no reading of u",
                     "the screen on all cells picked its near-threshold events by chance: they fail on a third",
                     "the count of library knockouts that pass on a third against the count that passed on all cells (7, of which 5 are not the readout's own genes)",
                     proposal_units=["the screens of 13 and 22 November read the protein response of every library knockout on all its cells"], test_units=NEW,
                     note="the system was never read on these units for this readout; their protein response was"),
    outer.hypothesis("E2", ["O12"], "rows",
                     {"quantity": "MHC-I in melanoma under IFN-gamma, the admitted events (expected B2M and TAPBP), on the two thirds of their cells not used for admission",
                      "head as trained, B2M": [0.85, 1.0], "head as trained, TAPBP": [0.85, 1.0],
                      "Z, own cells, B2M": [0.30, 0.90], "Z, own cells, B2M verdict reject": 0.85, "Z, own cells, TAPBP": [0.40, 0.95], "Z, own cells, TAPBP verdict reject": 0.85,
                      "input view, own cells, B2M (best reader)": [0.20, 0.85], "reading upstream, B2M": 0.30, "reading upstream, TAPBP": 0.25,
                      "the group passes A3 (the head loses what other cells of this kind read)": 0.45,
                      "if it does: the repair by rows adds no loss beyond the reference": 0.70, "if it does: test outcome": {"unresolved": 0.5, "confirmed": 0.3, "refuted": 0.2},
                      "some event supported by some reader at some stage (positive control)": 0.15},
                     "loop.test with the three arms: rows of the admitted events; rows of the five pathway genes of the same type, which move MHC-I little (cells of a kind already "
                     "covered); the repair of co-culture applied here (another map). Gate 3 as amended. Read per event: two events give no statement about a family",
                     "these knockouts act on the measured complex itself (B2M is its light chain, TAPBP loads it), as CD274 does for PD-L1: no reader of RNA carries them, and the item "
                     "goes to the claim (a limit, or a system that is told the perturbation), not to rows. The field would be 'relation to the readout'",
                     "the three arms of the protocol for a coverage claim",
                     proposal_units=U["mel_ifn"] + U["mel_co"] + ["CD274@mel_ifn"], test_units=NEW,
                     note="forecasts read from the five pathway genes (MHC-I, head as trained 0.88) and from the CD274 probe (0.91 at Z with its own cells)"),
    outer.hypothesis("E3", ["O21"], "rule",
                     {"quantity": "units set aside as final: a third of the units of each stratum that the system has not been read on", "share per stratum": [0.30, 0.37]},
                     "drawn once by unit with a fixed seed, stratified by study and condition, before any state is computed for them; a run that touches them is refused unless it is the final audit",
                     "a later split along a field that was not stratified leaves one side without final units: then that side has no certified statement",
                     "the draw is by unit name only: no response and no state enters it",
                     proposal_units=[], note="no reading is made"),
]

DIFF = [
    {"lever": "P", "change": "3a gate 3 is stated against the reference stage", "items": ["O14"], "hypothesis": "P1"},
    {"lever": "P", "change": "3b the test of a repair by rows has two controls", "items": ["O15"], "hypothesis": "P2"},
    {"lever": "P", "change": "3c a group ready as a pool with no admitted event is not a target", "items": ["O20"], "hypothesis": "P3"},
    {"lever": "P", "change": "3d a selection split for every new admission", "items": ["O19"], "hypothesis": "P4"},
    {"lever": "S", "change": "3e the type of a cell comes from the record", "items": ["O18", "O02"], "hypothesis": "S1"},
    {"lever": "S", "change": "the THP-1 pool is read as two types, provisionally: four genes that signal to STAT1, and IRF1", "items": ["O09", "O08", "O11"], "hypothesis": "S2"},
    {"lever": "V", "change": "reader specification: the input view is read by ridge and boosting on 512 components and by ridge on 4,000 genes; Z by ridge and boosting; a permuted view as control",
     "items": ["O05", "O06"], "hypothesis": "V1"},
    {"lever": "V", "change": "the same reader specification for MHC-I in co-culture, with the positive control evaluated", "items": ["O13"], "hypothesis": "V2"},
    {"lever": "E", "change": "the events of MHC-I in melanoma under IFN-gamma are selected from the library on a selection third", "items": ["O12"], "hypothesis": "E1"},
    {"lever": "E", "change": "a new target: MHC-I in melanoma under IFN-gamma with the admitted events", "items": ["O12"], "hypothesis": "E2"},
    {"lever": "E", "change": "a census from the record, and a third of the unread units of each stratum set aside as final", "items": ["O21"], "hypothesis": "E3"},
]
NOT_ACTED = [
    {"items": ["O01"], "why": "the owner's decision: a limit for knockouts of the readout's own gene, or a system that is told the perturbation"},
    {"items": ["O03", "O04", "O10"], "why": "need units the system has not been read on; no such unit has states"},
    {"items": ["O07"], "why": "at the margin as it stands; read again after V1"},
    {"items": ["O16", "O17"], "why": "the owner's decision, or units of the pool without states"},
]


def md():
    f = lambda v: json.dumps(v, ensure_ascii=False).replace('"', "") if not isinstance(v, str) else v
    o = ["| change | lever | items | hypothesis |", "|---|---|---|---|"] + [f"| {d['change']} | {d['lever']} | {' '.join(d['items'])} | {d['hypothesis']} |" for d in DIFF]
    o += ["", "| items not acted on | why |", "|---|---|"] + [f"| {' '.join(d['items'])} | {d['why']} |" for d in NOT_ACTED]
    for h in H:
        o += ["", f"### {h['id']} ({h['claim']}; items {' '.join(h['item'])})", "", f"- **quantity**: {h['forecast']['quantity']}"]
        o += [f"- **forecast**: " + "; ".join(f"{k}: {f(v)}" for k, v in h["forecast"].items() if k != "quantity")]
        o += [f"- **rule**: {h['rule']}", f"- **rival**: {h['rival']}", f"- **control**: {h['control']}",
              f"- **units**: proposal {', '.join(h['proposal_units']) or 'none'}; test {', '.join(h['test_units']) or 'none'}" + (f". {h['note']}" if h["note"] else "")]
    return "\n".join(o) + "\n"


def main():
    items = json.loads((RUN.parent / "outer_loop_design_20261120" / "results" / "open_items.json").read_text())["items"]
    ids = {i["id"] for i in items}
    named = {k for d in DIFF + NOT_ACTED for k in d["items"]}
    assert named == ids, sorted(ids ^ named)                # every open item is acted on or has a reason why not
    findings = outer.outer_check([], items, diff=DIFF, hypotheses=H, roles=None)
    assert not [x for x in findings if x["level"] == "error"], findings
    (RUN / "results").mkdir(exist_ok=True)
    (RUN / "results" / "declared.json").write_text(json.dumps({"package": 2, "against": "package 1: runs/teg_loop_20261120", "diff": DIFF, "not_acted_on": NOT_ACTED, "hypotheses": H,
                                                                 "ledger_check": findings}, indent=1, ensure_ascii=False) + "\n")
    (RUN / "results" / "declared.md").write_text(md())
    print(md())


if __name__ == "__main__":
    main()
