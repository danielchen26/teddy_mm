# The design of the first training experiment the loop asks for

27 November 2026. **Not run.** [`DESIGN.md`](DESIGN.md) says what would be trained (four pooling vectors, 2,048
numbers; the encoder frozen), by what loss (a readout taught by control cells alone has to read the recorded
events), on which cells (development cells of the THP-1 stratum, read before), and how the loop judges it.
[`scripts/pool.py`](scripts/pool.py) is the training script; `python scripts/pool.py --selftest` runs it on made-up
cells with no model and no data.

It waits for the owner's word. Nothing was trained.
