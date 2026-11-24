"""The outer step: what a round leaves, and which sealed thing the next iteration changes (OUTER_LOOP.md).

A round (LOOP_PROTOCOL.md) changes only the model. It holds five other things fixed: the events that are read (E), the
views (V), the schema (S), the claim (C) and the protocol (P). This module turns the state of the targets at the end of
a round into typed open items, checks whether each reading can be used, sends each item to a lever, and checks the
ledger of an outer step (rules O1 to O6).

Nothing here reads a cell or fits anything. The input is a list of targets in a plain format (`open_items`); a run
supplies an adapter from its own result files.
"""
import math

LEVERS = {"M": "the model (inner loop)", "E": "events: which units of the pool are read next, and in which role",
          "V": "views: the readers of a stage, the input, the readouts", "S": "schema: the fields of an event and what is pooled",
          "C": "claim: the targets, the grade and the tolerance asked for", "P": "protocol: the gates and the rules"}

# reading -> (lever, what the reading means, admissible actions in order of cost, what not to do).  OUTER_LOOP.md, Section 7.
ROUTES = {
    "invalid": ("V", "the reading cannot be used (Section 6)", ["repair the reading before any action, as Section 6 says for the check that failed"], "act on the reading"),
    "unresolved": ("E", "the interval does not decide between support and rejection",
                   ["if the units of the pool differ among themselves, the split comes first", "more independent units of the same type from the pool",
                    "if the pool cannot supply them, report the estimate as it is: at the margin"],
                   "change the model or the views"),
    "small response": ("E", "the response is too small against its interval to be followed",
                       ["more cells per event, or more events that move this readout", "otherwise the readout is not asked for this type"], "make it a target"),
    "unreadable": ("V", "a resolved loss at the input view, with the event's own cells, in a view that is not dominated",
                   ["a richer reader of the same input", "another input already measured on the same cells, never the readout itself", "a new measurement", "a stated limit of the claim"],
                   "train"),
    "upstream": ("M", "the input view reads it and the stage the head reads does not",
                 ["change what comes before the stage: an inner loop of another depth", "if that part cannot be changed, the item stays open under its owner"], "change the head"),
    "no rows": ("E", "cells of this kind are needed and the table has no admitted event of the type", ["events of the type from the pool"], "carry rows over from another type"),
    "pooled, and differs": ("S", "a distinction the schema pools is visible in the readings",
                            ["split along a field from the record, with rival fields, tested on new units"], "repair across the pool"),
    "pooled, may be the same": ("S", "a distinction the schema keeps looks invisible to the readout", ["merge after bounded support on new units, twice"], "merge on development readings"),
    "mistyped": ("S", "cells get a type from a reader of their state, and some get the wrong one", ["take the type from the record"], "tune the reader"),
    "self event": ("C", "the event acts on the readout itself: a system that is not told the assignment may have nothing to read (rule S3)",
                   ["a stated limit for this type", "or the assignment as an input, which is another system: the owner's decision"], "use it as a row"),
    "rule": ("P", "a gate or a test could not be evaluated as written", ["amend the protocol before the next iteration"], "re-grade the finished run"),
    "selection": ("P", "an event was selected by a measurement on the cells that were then audited", ["a selection split"], "keep selecting on audited cells"),
    "no final data": ("E", "no unit of the claim is unread", ["a census of the pool; units set aside as final before anything else is read"], "call a development reading a confirmation"),
    "confounded": ("E", "two fields vary together in the units that were read", ["units of the pool in which they vary apart", "otherwise the claim names their product"],
                   "attribute the difference to one of them"),
    "guard conflict": ("C", "every change that removes the loss breaks the guard on the system's own task",
                       ["a change that holds the guard: the inner loop again", "if none does, the owner decides between the target and the guard"], "widen the guard after the result"),
}
# the readings the inner loop keeps for itself (LOOP_PROTOCOL.md, version 2); they are statuses of a target, not open items of the outer step
INNER = ("function", "coverage", "unit")

# the two places where an item is worded for a system; a run passes its own words in `world["words"]`
WORDS = {"self event": "the event acts on the readout itself", "upstream": "a change of what comes before the stage the model reads"}
STATUSES = ("closed (absolute)", "the model's part closed (relative)", "open in the model", "undecided in the model", "open before the model", "undecided before the model", "not graded",
            "not a target")
STATES = ("open", "decided on development units", "closed", "limit", "at the margin", "not asked", "on the experiment list")
TERMINAL = STATES[2:]
NEEDS = {"invalid": "no new data", "rule": "no new data", "selection": "no new data", "mistyped": "no new data", "unreadable": "no new data first; then the owner's decision",
         "self event": "the owner's decision", "guard conflict": "the owner's decision", "small response": "units of the pool, or the owner's decision", "no final data": "units of the pool",
         "confounded": "units of the pool", "pooled, and differs": "units of the pool", "pooled, may be the same": "units of the pool", "no rows": "units of the pool",
         "upstream": WORDS["upstream"]}
DEFAULTS = {"margin": 0.10, "step_margin": 0.10, "min_snr": 3.0, "min_units": 10, "margin_factor": 4, "typed_right": 0.99, "reader_classes": 2}


def needed_units(u, ci, n, bound=DEFAULTS["margin"]):
    """Units needed for the interval to leave the bound on the side of the estimate, if the estimate holds and the
    half-width falls as 1/sqrt(n). None: the estimate sits on the bound. A planning number, not a promise."""
    q = abs(u)
    lo, hi = sorted((abs(ci[0]) if ci[0] * ci[1] > 0 else 0.0, max(abs(ci[0]), abs(ci[1]))))
    half, gap = (hi - q, bound - q) if q < bound else (q - lo, q - bound)
    if half <= 0 or gap >= half:
        return n
    if gap <= 0:
        return None
    return int(math.ceil(n * (half / gap) ** 2))


def reading_checks(view, d=DEFAULTS):
    """Section 6: the checks a reading of a view passes before anything is drawn from it. Returns the names of the checks that fail.

    view: {"gap_to_later": {"ci": [lo, hi]} (the later stage minus this one, same scope), "positive_control": bool (some
           reader reaches support for some unit of the family on this scale), "conflict": bool, "classes": [names],
           "nonlinear": bool, "readout_control": bool}; a missing key is not checked, so an adapter supplies every key it can."""
    bad = []
    g = view.get("gap_to_later")
    if g and g["ci"][1] < 0:
        bad.append("dominated view")
    if view.get("positive_control") is False:
        bad.append("positive control")
    if view.get("conflict"):
        bad.append("one reader specification")
    if view.get("readout_control") is False:
        bad.append("readout control")
    return bad


def target_status(t):
    """Exactly one status for a target at the end of an iteration (Section 5)."""
    if not t.get("admitted", True):
        return "not a target"
    if t.get("control_adds_a_loss") is False:                   # the response is followed, and neither the model nor its control adds a loss beyond the reference stage
        ev = (t.get("model") or {}).get("events") or {}
        at = [((t.get("reference") or {}).get("units") or {}).get("verdict")]
        at += [(((x.get("before") or {}).get("reference") or {}).get("units") or x.get("reference") or {}).get("verdict") for x in ev.values()]
        return "open before the model" if "reject" in at else "undecided before the model"       # a loss before the model is said only where a stage before it is rejected
    g = t.get("gates", {})
    if g.get("evaluable") is False:
        return "not graded"
    m = t["model"]
    if m["added"]["reading"] == "loss":
        return "open in the model"
    if m["added"]["reading"] == "unresolved":
        return "undecided in the model"
    return "closed (absolute)" if m["own"]["verdict"] == "support" else "the model's part closed (relative)"


def _item(items, t, reading, question, value=None, ci=None, n=None, needs=None, note=None, event="pool", evidence=None, standing=None, key=None):
    lever, meaning, actions, dont = ROUTES[reading]
    obs = {"readout": t.get("readout"), "group": t.get("group"), "event": event, "value": value, "ci": ci, "evidence": evidence, "note": note}
    if key:
        for it in items:
            if it.get("key") == key and it["reading"] == reading:
                if not any(o["group"] == obs["group"] and o["value"] == value for o in it["observations"]):
                    it["observations"].append(obs)
                return it
    it = {"id": None, "key": key, "readout": t.get("readout"), "group": t.get("group"), "type": t.get("type"), "event": event, "reading": reading, "lever": lever, "question": question,
          "value": value, "ci": ci, "resampling_units": n, "units_needed": needs, "more_units_in_pool": (t.get("pool") or {}).get("more_units_of_the_type"),
          "next_step_needs": NEEDS.get(reading), "after": [], "standing": standing, "actions": actions, "not": dont, "evidence": evidence, "note": note, "state": "open",
          "observations": [obs]}
    items.append(it)
    return it


def _failed(items, t, iv, bad, n, event="pool", key=None):
    """An item for every check of Section 6 that the reading of an input view fails."""
    for b in bad:
        g = iv.get("gap_to_later", {}) if b == "dominated view" else {}
        it = _item(items, t, "invalid", f"the reading of the input view fails the check: {b}", g.get("added"), g.get("ci"), n, event=event, evidence=iv.get("evidence"),
                   note=iv.get("note") if b == "dominated view" else iv.get(b.replace(" ", "_") + "_note"), key=key)
        if b == "positive control":                             # a richer reader may reach support; if none does, the tolerance for this target is the owner's to declare
            it["next_step_needs"] = "no new data first; then the owner's decision"


def _before(items, t, ref, iv, bad, n, d, event="pool"):
    """Where a loss before the model is: at the input view (unreadable), between the input view and the stage the model reads
    (upstream), not decided, or not placed because the view fails a check of Section 6 or was not read."""
    r = ref["units"]
    if r["verdict"] == "unresolved":
        _item(items, t, "unresolved", "is anything lost before the model", r["u"], r["ci"], n, event=event, evidence=ref.get("evidence"), needs=needed_units(r["u"], r["ci"], n, d["margin"]) if n else None)
    elif r["verdict"] == "reject":
        if iv is None:
            _item(items, t, "invalid", "a loss before the model, and the input view was not read", r["u"], r["ci"], n, event=event, evidence=ref.get("evidence"))
        elif not bad:                                           # a view that fails a check of Section 6 gives no reading of where the loss is
            o = iv["own"]
            gap = iv.get("gap_to_later")
            if o["verdict"] == "reject":
                weak = len(iv.get("classes", [])) < d["reader_classes"] or not iv.get("nonlinear")
                it = _item(items, t, "unreadable", "a loss at the input view, with the event's own cells", o["u"], o["ci"], n, event=event, evidence=iv.get("evidence"),
                           note="rests on one reader class or on linear readers only: a richer reader comes before any conclusion about the input" if weak else iv.get("note"))
                more = lambda x: (x.get("reads_more_than_the_input") or {}).get("ci", [0])[0] > d["step_margin"]      # on paired resamples, at least a step more than the model's input
                reads = sorted(k for k, x in (iv.get("other_inputs") or {}).items() if x["verdict"] == "support" or more(x))
                if reads and not weak:                          # the second step of the ladder is at hand: an input already measured that the model does not receive
                    it["next_step_needs"] = "no new data"
                    it["forecast"] = "another input already measured reads it: " + ", ".join(reads)
            if o["verdict"] == "support" or (gap and gap["ci"][0] > d["step_margin"]):       # a loss can sit in two places: at the input, and between the input and this stage
                where = "; and ".join(f"between {a} and {b}" for a, b in iv.get("lost_between") or [])      # where the run read stages inside the model
                _item(items, t, "upstream", "the input view reads it and the reference stage does not", gap["added"] if gap else None, gap["ci"] if gap else None, n, event=event, evidence=iv.get("evidence"),
                      note=("a step is lost " + where) if where else None)
            elif o["verdict"] == "unresolved":
                _item(items, t, "unresolved", "is the loss before the model already at the input view", o["u"], o["ci"], n, event=event, evidence=iv.get("evidence"),
                      needs=needed_units(o["u"], o["ci"], n, d["margin"]) if n else None)


def open_items(targets, world=None, d=DEFAULTS):
    """The open items of an iteration. Every reading is derived by rule from what the round recorded.

    targets: a list of dicts, one per (type, readout):
      readout, group, type {field: value}, admitted (bool), response {"snr": float}, events_n (resampling units),
      pool {"more_units_of_the_type": int or None, "admissible_now": [events]},
      reference {"units": {u, ci, verdict}, "own": {u, ci, verdict}}   the stage the model reads
      input_view {"own": {u, ci, verdict}, "gap_to_later": {added, ci}, "classes": [..], "nonlinear": bool, ...}  or None
      model {"own": {u, ci, verdict}, "added": {added, ci, reading}, "reading": one of INNER or None, "rows_of_type": int,
             "events": {site: {u, ci, verdict, "added": {reading}, "unseen": {u, ci, verdict, "own_row_adds": {gain_difference, ci, reading}}}}}
      attainable ("absolute" or "relative"),
      gates {"evaluable": bool, "why": [..], "guard_kept": bool, "guard_holds_under_some_change": bool},
      arms {"other": {added, ci, reading, "pair": str}, "rows_other": {...}},
      self_events [{site, u, ci}], confounded [(field, field)]
    world: {"typed_right": float, "typed_group": str, "selection_on_audited_cells": str or None, "final_units": int, "undeclared_cases": [str],
            "words": {"self event": str, "upstream": str}  (how two items are worded for this system; default WORDS)}
    """
    items = []
    words = {**WORDS, **((world or {}).get("words") or {})}
    for t in targets:
        n = t.get("events_n")
        iv = t.get("input_view")
        bad = reading_checks(iv, d) if iv else []
        r0 = t.get("reference")
        if not (r0 and r0["units"]["verdict"] == "support"):     # where the reference stage is supported nothing is lost before the model, and the input view is not consulted
            _failed(items, t, iv, bad, n)
        for e in t.get("self_events", []):
            _item(items, t, "self event", words["self event"], e.get("u"), e.get("ci"), None, event=e["site"], evidence=e.get("evidence"), note=e.get("note"))
        for pair in t.get("confounded", []):
            _item(items, t, "confounded", "two fields vary together: " + " and ".join(pair), key="confounded: " + " and ".join(sorted(pair)))
        if not t.get("admitted", True):
            snr = t["response"]["snr"]
            now = (t.get("pool") or {}).get("admissible_now")
            it = _item(items, t, "small response", "the change of the readout against the half-width of its interval", snr, None, n, evidence=t.get("evidence"), note=t.get("note"),
                       needs=None if snr >= d["min_snr"] else int(math.ceil(n * (d["min_snr"] / snr) ** 2)))
            if now:
                it["next_step_needs"], it["forecast"] = "units at hand", "within reach now: events with states at hand pass the screen: " + ", ".join(now)
            elif now is not None:
                it["next_step_needs"], it["forecast"] = "the owner's decision", "the pool at hand has no event that moves this readout: not asked, or on the experiment list"
            continue
        g = t.get("gates", {})
        for why in g.get("why", []) if g.get("evaluable") is False or g.get("why") else []:
            _item(items, t, "rule", why["question"], why.get("value"), None, None, evidence=why.get("evidence"), note=why.get("note"))
        if g.get("guard_kept") is False and g.get("guard_holds_under_some_change") is False:
            _item(items, t, "guard conflict", "every change that removes the loss breaks the guard", evidence=g.get("evidence"))
        m, ref = t["model"], t.get("reference")
        added = m["added"]
        if added["reading"] == "loss":
            if m.get("reading") == "coverage" and not m.get("rows_of_type", 1):
                _item(items, t, "no rows", "cells of this kind are needed and no event of the type is admitted", added["added"], added["ci"], n, evidence=m.get("evidence"))
            # any other loss the model still adds is the inner loop's: a status of the target, not an open item of the outer step
        elif added["reading"] == "unresolved":
            _item(items, t, "unresolved", "does the model add a loss beyond the reference stage", added["added"], added["ci"], n, evidence=m.get("evidence"),
                  needs=needed_units(added["added"], added["ci"], n, d["step_margin"]) if n else None)
        elif t.get("attainable") == "absolute" and m["own"]["verdict"] != "support":
            _item(items, t, "unresolved", "is the output of the model inside the margin", m["own"]["u"], m["own"]["ci"], n, evidence=m.get("evidence"),
                  note="the reference stage is supported: nothing is left before the model", needs=needed_units(m["own"]["u"], m["own"]["ci"], n, d["margin"]) if n else None)
        # what is left before the model: for the pool, and for an event that the run read on its own (its own reference and input view)
        # an output inside the margin has lost nothing on the way: no stage before it is consulted for that pool or event
        if ref and m["own"]["verdict"] != "support":            # a target read event by event has no reference of its pool
            _before(items, t, ref, iv, bad, n, d)
        for s, x in (m.get("events") or {}).items():
            if x.get("before") and x["verdict"] != "support":
                b = x["before"]
                bb = reading_checks(b["input_view"], d) if b.get("input_view") else []
                _failed(items, t, b.get("input_view"), [c for c in bb if c != "positive control"], None, event=s)
                if "positive control" in bb:                    # a check of the family: one item for the target, whichever event shows it
                    _failed(items, t, b["input_view"], ["positive control"], None, key=f"positive control: {t.get('readout')} {t.get('group')}")
                _before(items, t, b["reference"], b.get("input_view"), bb, None, d, event=s)
        # events that read differently from their pool: at the model's output, or already at the reference stage
        inpool = {s: x for s, x in (m.get("events") or {}).items() if x.get("in_pool", True)}      # an event read on its own is not pooled: the schema keeps it apart already
        held = [k for k, y in inpool.items() if (y.get("reference") or {}).get("verdict") == "support"]
        for s, x in inpool.items():
            at_output = x["verdict"] == "reject" and x["ci"][0] > m["own"]["ci"][1]
            at_reference = bool(held) and (x.get("reference") or {}).get("verdict") == "reject"
            if at_output or at_reference:
                note = [f"the pool reads {m['own']['u']:.2f} [{m['own']['ci'][0]:.2f}, {m['own']['ci'][1]:.2f}]; this event's interval lies above it"] if at_output else []
                if at_reference:
                    note.append(f"at the reference stage this event is rejected ({x['reference']['u']:.2f}) and {len(held)} of the pool are supported: the pool mixes a loss before the model with none")
                elif (x.get("added") or {}).get("reading") == "no added loss":
                    note.append("for this event the loss is already at the reference stage")
                _item(items, t, "pooled, and differs", f"one event inside the pool: {s}", x["u"], x["ci"], None, event=s, evidence=m.get("evidence"), note="; ".join(note))
        # a unit event whose repair needs the unit's own row: the rows of the other units of its type do not carry to it
        for s, x in (m.get("events") or {}).items():
            own = (x.get("unseen") or {}).get("own_row_adds")
            if own and own["reading"] == "a step":
                _item(items, t, "pooled, and differs", f"the rows of the other units do not carry to this unit: {s}", own["gain_difference"], own["ci"], None, event=s, evidence=m.get("evidence"),
                      note=f"with its own row the output leaves {x['u']:.2f}; with the rows of the other units alone, {x['unseen']['u']:.2f}")
        # what carries across types
        arms = t.get("arms") or {}
        o, ro = arms.get("other"), arms.get("rows_other")
        if o and o["reading"] == "loss":
            _item(items, t, "pooled, and differs", "the whole repair across: " + o["pair"], o["added"], o["ci"], n, evidence=o.get("evidence"), key="repair across: " + o["pair"],
                  standing="the schema already keeps them apart; seen on development cells, to be confirmed on new units")
        if ro and ro["reading"] == "no added loss":
            _item(items, t, "pooled, may be the same", "the rows across: " + ro["pair"], ro["added"], ro["ci"], n, evidence=ro.get("evidence"), key="rows across: " + ro["pair"])
        elif ro and ro["reading"] == "unresolved":
            _item(items, t, "unresolved", "can the rows be merged across: " + ro["pair"], ro["added"], ro["ci"], n, evidence=ro.get("evidence"), key="rows across: " + ro["pair"],
                  needs=needed_units(ro["added"], ro["ci"], n, d["step_margin"]) if n else None)
    w = world or {}
    for it in items:
        if it["reading"] == "upstream":
            it["next_step_needs"] = words["upstream"]
    wt = {"readout": "all", "group": w.get("typed_group")}
    if w.get("typed_right") is not None and w["typed_right"] < d["typed_right"]:
        _item(items, wt, "mistyped", "how a cell gets its type", w["typed_right"], evidence=w.get("typed_evidence"), note="the type is guessed from the state; the record has it")
    if w.get("selection_on_audited_cells"):
        _item(items, {"readout": "all", "group": None}, "selection", "events were screened on the cells that were then audited", note=w["selection_on_audited_cells"])
    for case in w.get("undeclared_cases", []):
        _item(items, {"readout": "all", "group": None}, "rule", case)
    if w.get("final_units") == 0:
        _item(items, {"readout": "all", "group": None}, "no final data", "units of the claim that no run has read", 0, note=w.get("final_note"))
    # ids, dependencies, and what the pool can do for an unresolved item
    for k, it in enumerate(items):
        it["id"] = f"O{k + 1:02d}"
    mist = [j["id"] for j in items if j["reading"] == "mistyped"]
    for it in items:
        it["approximate"] = bool(it["resampling_units"] is not None and it["resampling_units"] < d["min_units"])
        if it["reading"] != "unresolved":
            continue
        n, need, have = it["resampling_units"], it["units_needed"], it["more_units_in_pool"]
        it["after"] = [j["id"] for j in items if j["reading"] == "pooled, and differs" and j["group"] == it["group"] and j["readout"] == it["readout"] and j["event"] != "pool"]
        if need is None:
            it["forecast"], it["next_step_needs"] = "at the margin: the estimate sits on the bound", "at the margin"
        elif n is not None and (need - n) > (d["margin_factor"] - 1) * n:
            it["forecast"], it["next_step_needs"] = f"at the margin: more than {d['margin_factor']} times the units would be needed", "at the margin"
        elif have is not None and n is not None and (need - n) > have:
            it["forecast"], it["next_step_needs"] = "at the margin: the pool has no further unit that passes admission", "at the margin"
        else:
            it["forecast"], it["next_step_needs"] = "within reach if the pool holds units of the type: a census first", "units of the pool"
        # an estimate read on cells whose type was guessed is read again once the type comes from the record
        if mist and it["question"] == "is the output of the model inside the margin" and it["group"] == w.get("typed_group"):
            it["after"] += mist
            it["forecast"] = "first the bridge reading with the type from the record, which changes the cells this estimate is read on; then: " + it["forecast"]
            it["next_step_needs"] = "no new data"
    return items


def counts(items):
    return {"by_lever": {k: sum(i["lever"] == k for i in items) for k in LEVERS},
            "by_reading": {k: sum(i["reading"] == k for i in items) for k in ROUTES},
            "by_what_the_next_step_needs": {k: sum(i["next_step_needs"] == k for i in items) for k in sorted({i["next_step_needs"] for i in items})}}


def coverage(targets, items):
    """The table a campaign reports: every target with its status, every item with its state (Section 12)."""
    st = [target_status(t) for t in targets]
    return {"targets": {s: st.count(s) for s in STATUSES}, "items": {s: sum(i["state"] == s for i in items) for s in STATES},
            "closed": bool(targets) and all(s == "closed (absolute)" for s in st) and all(i["state"] == "closed" for i in items)}


CLAIMS = ("resolution", "view", "split", "merge", "rows", "limit", "rule")


def hypothesis(hid, item, claim, forecast, rule, rival, control, proposal_units=(), test_units=(), note=None):
    """An action of an outer step, written before it is taken (Section 10)."""
    if claim not in CLAIMS:
        raise ValueError(f"claim must be one of {CLAIMS}")
    return {"kind": "outer_hypothesis", "id": hid, "item": item, "claim": claim, "forecast": forecast, "rule": rule, "rival": rival, "control": control,
            "proposal_units": sorted(proposal_units), "test_units": sorted(test_units), "note": note}


def outer_check(targets, items, diff=(), hypotheses=(), roles=None, reads=None, claim_changes=(), cover=None):
    """Rules O1 to O6 of an outer step. Returns a list of {"rule", "level", "message"}; an empty list is a clean ledger.

    diff: [{"lever", "change", "items": [ids], "hypothesis": id or None}]   every difference between two packages
    roles: {"final": [units], ...};  reads: {unit: times read}
    claim_changes: [{"change", "limit": {"item", "reading_passed_checks": bool, "reopens": str}}]
    """
    out = []
    say = lambda rule, level, msg: out.append({"rule": rule, "level": level, "message": msg})
    keys = [(t.get("readout"), t.get("group")) for t in targets]
    for k in {k for k in keys if keys.count(k) > 1}:
        say("O1", "error", f"target {k} appears more than once")
    for t in targets:
        try:
            s = target_status(t)
        except (KeyError, TypeError) as e:
            say("O1", "error", f"target {(t.get('readout'), t.get('group'))} has no status: missing {e}")
            continue
        if s not in STATUSES:
            say("O1", "error", f"target {(t.get('readout'), t.get('group'))} has an unknown status")
    ids = {i["id"] for i in items}
    for i in items:
        if i["reading"] not in ROUTES:
            say("O2", "error", f"{i['id']}: reading {i['reading']!r} is not in the table")
        if i["state"] not in STATES:
            say("O2", "error", f"{i['id']}: unknown state {i['state']!r}")
        elif i["state"] in TERMINAL and not i.get("evidence") and not i.get("terminal_evidence"):
            say("O2", "error", f"{i['id']}: a terminal state without evidence")
        if i["state"] == "limit" and i["reading"] in ("unresolved", "invalid", "small response"):
            say("O2", "error", f"{i['id']}: an item read {i['reading']!r} cannot become a limit")
    hyp = {h["id"]: h for h in hypotheses}
    for c in diff:
        if not c.get("items"):
            say("O3", "error", f"a change to {c['lever']} names no open item: {c['change']}")
        for k in c.get("items", []):
            if k not in ids:
                say("O3", "error", f"a change to {c['lever']} names {k}, which is not an open item")
        h = hyp.get(c.get("hypothesis"))
        if h is None:
            say("O4", "error", f"a change to {c['lever']} has no hypothesis written before it: {c['change']}")
            continue
        for f in ("forecast", "rule", "rival", "control"):
            if not h.get(f):
                say("O4", "error", f"hypothesis {h['id']} has no {f}")
        both = set(h["proposal_units"]) & set(h["test_units"])
        if both:
            say("O4", "error", f"hypothesis {h['id']} is tested on units its proposal was read from: {sorted(both)}")
    final = set((roles or {}).get("final", []))
    for u in sorted(final & set(reads or {})):
        if (reads or {})[u]:
            say("O5", "error", f"final unit {u} has been read")
    for h in hypotheses:
        hit = final & (set(h["proposal_units"]) | set(h["test_units"]))
        if hit:
            say("O5", "error", f"hypothesis {h['id']} touches final units: {sorted(hit)}")
    if roles is not None and not final:
        say("O5", "warning", "no unit is set aside as final: every status is exploratory")
    for c in claim_changes:
        lim = c.get("limit")
        if not lim or not lim.get("reopens"):
            say("O6", "error", f"a change of the claim without a limit that says what would reopen it: {c['change']}")
        elif not lim.get("reading_passed_checks"):
            say("O6", "error", f"a limit resting on a reading that did not pass the checks: {c['change']}")
    if claim_changes and cover is None:
        say("O6", "error", "the claim changed and the coverage table was not computed again")
    return out
