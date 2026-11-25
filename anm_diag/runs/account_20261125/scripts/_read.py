"""What the readings of this run share: the cells of a stratum (as in runs/typed_repair_20261125), its typed states, and
the account of one event at several stages."""
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "typed_repair_20261125" / "scripts"))
from anm_diag import account as AC  # noqa: E402
import round_a  # noqa: E402
import round_b  # noqa: E402
from repair import SEED, unit  # noqa: E402,F401

cells, typed = round_a.cells, round_b.typed
KEYS = ("carried", "read", "others", "found_again")


def strip(r):
    return {k: v for k, v in r.items() if k != "_draws"}


def account(stages, c, event, others, inputs=("V",), state="Z", pairs=(), say=print, controls=True):
    """The readings of one event at every stage, their chain, the gains named in `pairs` (b over a), what a readout made
    to follow the mean adds to what the state carries, and the two controls at the state."""
    R = {}
    for name, xs in stages.items():
        R[name] = AC.reading(xs, c["y"], c["unit"], c["first"], event, others=others, o=c["o"] if name == state else None, seed=SEED)
        r = R[name]
        say(f"{event:7s} {name:14s} overlap {r['overlap']:.2f} ({r['carried_by']}) | " + " | ".join(f"{k} {v['share']:5.0%} [{v['ci'][0]:4.0%},{v['ci'][1]:4.0%}]" for k, v in r["shares"].items()))
    ch = AC.chain(R, list(inputs), state)
    say(f"        chain: " + " | ".join(f"{k} {v['share']:.0%} [{v['ci'][0]:.0%},{v['ci'][1]:.0%}] {v['reading']}" for k, v in ch["parts"].items()) + f"  -> {ch['says']}")
    out = {"stages": {n: strip(r) for n, r in R.items()}, "chain": ch,
           "gains": {f"{b} over {a}": {k: AC.gain(R[a], R[b], k) for k in KEYS if k in R[a]["_draws"] and k in R[b]["_draws"]} for a, b in pairs if a in R and b in R}}
    d = R[state]["_draws"]["found_again"] - R[state]["_draws"]["carried"]
    iv = AC.interval(d[1:])
    out["found_again_over_carried_at_the_state"] = {"gain": float(d[0]), "ci": iv, "reading": AC.step(iv)}
    if controls:
        out["must_fail"] = AC.must_fail(stages[state], c["y"], c["unit"], c["first"], event, seed=SEED)
    return out
