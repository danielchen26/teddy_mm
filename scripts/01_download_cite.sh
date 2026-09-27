#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST_DIR="$ROOT/data/raw"
DEST="$DEST_DIR/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz"
mkdir -p "$DEST_DIR"
if [[ -f "$DEST" ]]; then
  echo "already have $DEST"
  exit 0
fi
echo "downloading CITE BMMC (~587 MB gz) ..."
curl -L --fail --retry 3 -o "$DEST.gz" "$URL"
echo "decompressing ..."
gzip -d -f "$DEST.gz"
echo "ready: $DEST"
ls -lh "$DEST"
