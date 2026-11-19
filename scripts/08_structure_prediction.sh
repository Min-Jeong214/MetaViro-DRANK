#!/usr/bin/env bash
set -euo pipefail

QUERY=${1:?Usage: $0 query.faa outdir}
OUT=${2:?}

mkdir -p "$OUT"

###############################################################################
# Stage F representative single-sequence structure prediction
#
# Production configuration reconstructed from saved config.json files:
#
#   ColabFold              1.6.1
#   model_type             alphafold2_ptm
#   msa_mode               single_sequence
#   num_models             1
#   num_recycle            1
#   num_seeds              1
#   random_seed            0
#   rank                   plddt
#   templates              false
#   relax                  none
#   max_seq                1
#   max_extra_seq          1
#
# Production output model:
#   alphafold2_ptm_model_1_seed_000
#
# NOTE:
# Numerical reproduction should use ColabFold 1.6.1.
###############################################################################

echo "=================================================="
echo "STAGE F STRUCTURE PREDICTION"
echo "=================================================="

if command -v colabfold_batch >/dev/null 2>&1; then
    echo "COLABFOLD:"
    colabfold_batch --version 2>&1 || true
else
    echo "[ERROR] colabfold_batch not found"
    exit 1
fi

echo
echo "QUERY=$QUERY"
echo "OUT=$OUT"

colabfold_batch \
    --msa-mode single_sequence \
    --model-type alphafold2_ptm \
    --num-models 1 \
    --num-recycle 1 \
    --num-seeds 1 \
    --random-seed 0 \
    --rank plddt \
    --max-seq 1 \
    --max-extra-seq 1 \
    "$QUERY" \
    "$OUT"
