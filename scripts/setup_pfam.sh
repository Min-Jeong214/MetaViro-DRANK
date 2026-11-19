#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DB="$ROOT/databases/pfam"
HMM="$DB/Pfam-A.hmm"
GZ="$DB/Pfam-A.hmm.gz"
URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

mkdir -p "$DB"

if [[ ! -s "$HMM" ]]; then
    echo "===== DOWNLOAD PFAM ====="
    curl -fL --retry 3 -C - -o "$GZ" "$URL"
    gzip -t "$GZ"
    gzip -dc "$GZ" > "$HMM.partial"
    mv "$HMM.partial" "$HMM"
fi

echo "===== BUILD HMM INDEX ====="
hmmpress -f "$HMM"

for ext in h3f h3i h3m h3p; do
    test -s "$HMM.$ext"
done

echo "PFAM_INSTALL=PASS"
echo "DATABASE=$HMM"
