"""Check downloaded scPerturb files against the MD5 checksums published on the Zenodo record.

  python tools/verify_downloads.py DIR [CHECKSUM_FILE]

Record: "scPerturb Single-Cell Perturbation Data: RNA and protein h5ad files", version 1.3,
https://zenodo.org/records/10044268 (checksums as listed on the record page).
"""
import hashlib
import sys
from pathlib import Path


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main(d, sums=None):
    sums = Path(sums) if sums else Path(__file__).with_name("scperturb_zenodo_10044268.md5")
    bad = 0
    for line in sums.read_text().split("\n"):
        if not line.strip():
            continue
        want, name = line.split()
        p = Path(d).expanduser() / name
        if not p.exists():
            print(f"MISSING   {name}")
            bad += 1
            continue
        got = md5(p)
        ok = got == want
        bad += not ok
        print(f"{'OK      ' if ok else 'MISMATCH'}  {name}  {p.stat().st_size / 1e6:.1f} MB")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
