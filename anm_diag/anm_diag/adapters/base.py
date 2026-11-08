"""Adapter protocol: turns a dataset and a system into the arrays a declaration names.

An adapter implements
  cells(decl)               -> {"A": int array (1 treated, 0 control), "cluster": str array, "ids": array}
  readout(rspec, cells)     -> float array (n,)        the observation, from a measurement only
  state(node, cells)        -> float array (n, p)      a candidate state or probe, full features
  provenance()              -> dict                    data/code hashes and versions for the record
  exact(rspec, cells, delta)-> dict or None            optional exact same-state test (system-specific)
Operations are looked up by the node's `op` in the adapter's registry; an unknown op is an error, so a
declaration can only name operations an adapter really implements.
"""
import hashlib
from pathlib import Path


class Adapter:
    ops = {}

    def cells(self, decl):
        raise NotImplementedError

    def readout(self, rspec, cells):
        raise NotImplementedError

    def state(self, node, cells):
        op = str(node["op"])
        name, _, arg = op.partition(":")
        fn = self.ops.get(name)
        if fn is None:
            raise KeyError(f"adapter {type(self).__name__} has no operation '{name}' (node '{node['id']}')")
        return fn(self, node, cells, arg)

    def provenance(self):
        return {}

    def exact(self, rspec, cells, delta):
        return None


def load_paths(paths):
    """A paths dict or YAML file; expands ~ in every string value (and in lists of strings)."""
    import os
    if isinstance(paths, (str, Path)):
        import yaml
        paths = yaml.safe_load(open(os.path.expanduser(str(paths))))
    out = {}
    for k, v in paths.items():
        if isinstance(v, str):
            out[k] = os.path.expanduser(v)
        elif isinstance(v, list):
            out[k] = [os.path.expanduser(x) if isinstance(x, str) else x for x in v]
        else:
            out[k] = v
    return out


def file_sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()
