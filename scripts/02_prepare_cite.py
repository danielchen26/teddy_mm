#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.data import prepare_cite
from teddy_mm.teddy_encoder import load_vocab


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw", type=Path, default=ROOT / "data/raw/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad")
    p.add_argument("--out", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--ckpt", type=Path, default=ROOT.parent / "teddy_mwe/ckpt/teddy_g_70M")
    p.add_argument("--test-sites", default="site4")
    p.add_argument("--val-fraction", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-cells", type=int, default=None)
    args = p.parse_args()
    if not args.raw.exists():
        raise SystemExit(f"missing {args.raw}\nRun: bash scripts/01_download_cite.sh")
    vocab = load_vocab(args.ckpt)
    meta = prepare_cite(
        raw_path=args.raw,
        out_dir=args.out,
        vocab=vocab,
        test_sites=[s.strip() for s in args.test_sites.split(",") if s.strip()],
        val_fraction=args.val_fraction,
        seed=args.seed,
        max_cells=args.max_cells,
    )
    print("prepared CITE:", meta)


if __name__ == "__main__":
    main()
