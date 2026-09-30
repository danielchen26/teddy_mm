#!/usr/bin/env bash
# Fetch TEDDY's gene-median file (used by the official preprocessing in scripts/03_embed_rna.py)
# from the Merck/TEDDY Hugging Face repo, pinned to commit 8c880aa6, and check its sha256.
# usage: bash scripts/01b_fetch_teddy_medians.sh [DEST]   (default data/reference/teddy_gene_medians.json)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${1:-$ROOT/data/reference/teddy_gene_medians.json}"
REV="8c880aa67674042ccdc254bccb5de081092a48bf"
URL="https://huggingface.co/Merck/TEDDY/resolve/${REV}/teddy/data_processing/utils/medians/data/teddy_gene_medians.json"
SHA="3d3637fa24ef81bf431aa42c8b85ffe656c22b999b2ef937428490704982c001"
check() { [ "$(shasum -a 256 "$1" | cut -d' ' -f1)" = "$SHA" ]; }
mkdir -p "$(dirname "$DEST")"
if [ -f "$DEST" ] && check "$DEST"; then
  echo "ok $DEST (sha256 $SHA)"; exit 0
fi
tmp="$DEST.part"
curl -fL --retry 3 -m 300 -o "$tmp" "$URL"
if ! check "$tmp"; then
  echo "sha256 mismatch for $URL" >&2; rm -f "$tmp"; exit 1
fi
mv "$tmp" "$DEST"
echo "ok $DEST (sha256 $SHA)"
