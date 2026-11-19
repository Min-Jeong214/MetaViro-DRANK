#!/usr/bin/env bash
set -euo pipefail
FAA=${1:?Usage: $0 all_viral_proteins.faa out_prefix [threads]}
OUT=${2:?}
THREADS=${3:-32}
TMP="${OUT}.tmp"
mkdir -p "$TMP"
# Frozen pilot criterion: 30% sequence identity and 70% coverage.
mmseqs easy-cluster "$FAA" "$OUT" "$TMP" \
  --min-seq-id 0.30 -c 0.70 --cov-mode 0 --threads "$THREADS"
