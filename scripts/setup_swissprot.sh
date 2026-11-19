#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT/databases/swissprot"
FASTA="$DB/swissprot.faa"
GZ="$DB/swissprot.faa.gz"

URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

mkdir -p "$DB"

if [[ -s "$FASTA" ]]; then
    echo "SWISSPROT=ALREADY_INSTALLED"
    exit 0
fi

command -v curl >/dev/null || {
    echo "ERROR: curl is required"
    exit 1
}

echo "===== DOWNLOAD SWISS-PROT ====="

curl -fL --retry 3 -C - \
    -o "$GZ" "$URL"

echo "===== VERIFY ARCHIVE ====="
gzip -t "$GZ"

echo "===== EXTRACT FASTA ====="

TMP="$FASTA.tmp.$$"
trap 'rm -f "$TMP"' EXIT

gzip -dc "$GZ" > "$TMP"
test -s "$TMP"

mv "$TMP" "$FASTA"

echo "===== VERIFY FASTA ====="
grep -m 1 '^>' "$FASTA" || true

echo "SWISSPROT_INSTALL=PASS"
echo "DATABASE=$FASTA"
