#!/usr/bin/env bash
set -euo pipefail

QUERY=${1:?Usage: $0 query.faa outdir swissprot_db phrog_seq_db pfam_hmm [threads]}
OUT=${2:?}
SWISSPROT=${3:?}
PHROG=${4:?}
PFAM=${5:?}
THREADS=${6:-32}

mkdir -p \
    "$OUT" \
    "$OUT/tmp_swiss" \
    "$OUT/phrog_work/tmp"


# ============================================================
# 1. Swiss-Prot
# ============================================================

mmseqs easy-search \
    "$QUERY" \
    "$SWISSPROT" \
    "$OUT/swissprot.tsv" \
    "$OUT/tmp_swiss" \
    -s 7.5 \
    -e 1e-3 \
    --max-seqs 50 \
    -a \
    --threads "$THREADS" \
    --format-output \
query,target,theader,pident,alnlen,qcov,tcov,evalue,bits


# ============================================================
# 2. PHROG — header-independent native MMseqs workflow
#
# Canonical DB:
#   phrogs_profile_db_seq
#
# easy-search is intentionally NOT used because this DB lacks
# the header information required by convertalis/easy-search.
# ============================================================

PHROG_QUERY="$OUT/phrog_work/queryDB"
PHROG_RESULT="$OUT/phrog_work/resultDB"
PHROG_TMP="$OUT/phrog_work/tmp"

rm -f \
    "${PHROG_QUERY}"* \
    "${PHROG_RESULT}"*

mmseqs createdb \
    "$QUERY" \
    "$PHROG_QUERY"

mmseqs search \
    "$PHROG_QUERY" \
    "$PHROG" \
    "$PHROG_RESULT" \
    "$PHROG_TMP" \
    -s 7.5 \
    -e 1e-5 \
    --max-seqs 50 \
    -a \
    --threads "$THREADS"

mmseqs createtsv \
    "$PHROG_QUERY" \
    "$PHROG" \
    "$PHROG_RESULT" \
    "$OUT/phrog.native.tsv"


# ============================================================
# 3. Pfam
#
# Canonical acceptance criterion:
# Pfam gathering thresholds (--cut_ga)
# ============================================================

hmmscan \
    --domtblout "$OUT/pfam.domtblout" \
    --cut_ga \
    --cpu "$THREADS" \
    "$PFAM" \
    "$QUERY" \
    > "$OUT/pfam.hmmscan.txt"


echo
echo "=============================================="
echo "ANNOTATION SEARCH COMPLETE"
echo "=============================================="

echo -n "SWISSPROT_ROWS="
wc -l < "$OUT/swissprot.tsv"

echo -n "PHROG_NATIVE_ROWS="
wc -l < "$OUT/phrog.native.tsv"

echo -n "PFAM_DOMAIN_ROWS="
grep -vc '^#' "$OUT/pfam.domtblout" || true
