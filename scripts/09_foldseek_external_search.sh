#!/usr/bin/env bash
set -euo pipefail

# STAGE09_INTEGRATION_MODE_V1
#
# Usage:
#   bash scripts/09_foldseek_external_search.sh --integrate [OUTPUT_TSV]
#
# Generates the 191-family external evidence matrix
# from existing Foldseek raw results.
# Does not rerun Foldseek or modify canonical data.

if [[ "${1:-}" == "--integrate" ]]; then
    shift

    if (( $# > 1 )); then
        echo "Usage: $0 --integrate [OUTPUT_TSV]" >&2
        exit 2
    fi

    SCRIPT_DIR="$(
        cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd
    )"

    REPO_ROOT="$(
        cd -- "$SCRIPT_DIR/.." && pwd
    )"

    DEFAULT_OUT="$REPO_ROOT/Result/Annotation_saturation/Pilot_v2/Stage09_Reconstruction_DryRun/structure191_external_structure_matrix.tsv"

    OUT_TSV="${1:-$DEFAULT_OUT}"

    CANON_TSV="$REPO_ROOT/Result/Annotation_saturation/Pilot_v2/StageF_Structure/StageG3_ExternalDB/Integrated_ExternalStructure/structure191_external_structure_matrix.tsv"

    if [[ "$(realpath -m -- "$OUT_TSV")" == "$(realpath -m -- "$CANON_TSV")" ]]; then
        echo "ERROR: Refusing to overwrite canonical matrix" >&2
        exit 2
    fi

    echo "STAGE09_MODE=INTEGRATE"
    echo "REPOSITORY=$REPO_ROOT"
    echo "OUTPUT=$OUT_TSV"

    exec python3 \
        "$SCRIPT_DIR/09_external_integration.py" \
        --root "$REPO_ROOT" \
        --out "$OUT_TSV"
fi


###############################################################################
# Usage:
#
#   09_foldseek_external_search.sh \
#       queryDB targetDB resultDB tmp out.tsv \
#       evalue max_seqs [threads]
#
# Examples reproducing Stage G3:
#
# PDB initial:
#   evalue=1e-3  max_seqs=200
#
# PDB E1:
#   evalue=1     max_seqs=1000
#
# AF/Swiss-Prot:
#   evalue=1     max_seqs=1000
#
# BFVD:
#   evalue=1     max_seqs=500
#
###############################################################################

QDB=${1:?queryDB required}
TDB=${2:?targetDB required}
RDB=${3:?resultDB required}
TMP=${4:?tmp directory required}
OUT=${5:?output TSV required}

EVALUE=${6:?evalue required}
MAXSEQS=${7:?max_seqs required}
THREADS=${8:-16}

FOLDSEEK=${FOLDSEEK:-foldseek}

mkdir -p "$TMP"
mkdir -p "$(dirname "$OUT")"

rm -f "${RDB}"*
rm -rf "$TMP"
mkdir -p "$TMP"

echo "=================================================="
echo "FOLDSEEK EXTERNAL STRUCTURE SEARCH"
echo "=================================================="

echo "QUERY_DB=$QDB"
echo "TARGET_DB=$TDB"
echo "RESULT_DB=$RDB"
echo "TMP=$TMP"
echo "OUTPUT=$OUT"
echo "EVALUE=$EVALUE"
echo "MAX_SEQS=$MAXSEQS"
echo "THREADS=$THREADS"

"$FOLDSEEK" search \
    "$QDB" \
    "$TDB" \
    "$RDB" \
    "$TMP" \
    -s 9.5 \
    -e "$EVALUE" \
    --max-seqs "$MAXSEQS" \
    --threads "$THREADS" \
    -a

"$FOLDSEEK" convertalis \
    "$QDB" \
    "$TDB" \
    "$RDB" \
    "$OUT" \
    --format-output \
"query,target,evalue,bits,fident,alnlen,qlen,tlen,qcov,tcov,lddt,alntmscore,qtmscore,ttmscore,prob"

echo
echo "ROWS=$(wc -l < "$OUT")"

if [[ -s "$OUT" ]]; then
    echo "COLUMNS_FIRST_ROW=$(awk -F'\t' 'NR==1{print NF}' "$OUT")"
else
    echo "COLUMNS_FIRST_ROW=NA"
fi
