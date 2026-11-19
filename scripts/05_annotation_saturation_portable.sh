#!/usr/bin/env bash
set -euo pipefail

QUERY=${1:?Usage: $0 query.faa outdir swissprot_db phrog_profile_db pfam_hmm [threads]}
OUT=${2:?}
SWISSPROT=${3:?}
PHROG=${4:?}
PFAM=${5:?}
THREADS=${6:-32}

# Maximum candidate sequences per PHROG profile.
# Increase for large datasets after checking resource usage.
PHROG_MAX_SEQS=${PHROG_MAX_SEQS:-1000}
[[ "$PHROG_MAX_SEQS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: invalid PHROG_MAX_SEQS" >&2
    exit 1
}
echo "PHROG_MAX_SEQS=$PHROG_MAX_SEQS"

mkdir -p \
    "$OUT" \
    "$OUT/tmp_swiss" \
    "$OUT/phrog_work/tmp"

# Prevent stale outputs from a previous run being accepted.
rm -f     "$OUT/search_run_complete.json"     "$OUT/swissprot.tsv"     "$OUT/phrog.native.tsv"     "$OUT/pfam.domtblout"     "$OUT/pfam.hmmscan.txt"


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
#   phrogs_profile_db
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
    "$PHROG" \
    "$PHROG_QUERY" \
    "$PHROG_RESULT" \
    "$PHROG_TMP" \
    -s 7.5 \
    -e 1e-5 \
    --max-seqs "$PHROG_MAX_SEQS" \
    -a \
    --threads "$THREADS"

mmseqs createtsv \
    "$PHROG" \
    "$PHROG_QUERY" \
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



# ============================================================
# 4. Verified completion manifest
# ============================================================

python - "$QUERY" "$OUT" "$SWISSPROT" "$PHROG" \
    "$PFAM" "$PHROG_MAX_SEQS" "$THREADS" <<'PY_MANIFEST'
import hashlib
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

query, out, swiss, phrog, pfam = [
    Path(x) for x in sys.argv[1:6]
]
max_seqs = int(sys.argv[6])
threads = int(sys.argv[7])

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def fingerprint(path):
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }

if not query.is_file() or query.stat().st_size == 0:
    raise SystemExit("ERROR: Empty or missing query FASTA")

sequence_count = sum(
    line.startswith(b">")
    for line in query.open("rb")
)

if sequence_count == 0:
    raise SystemExit("ERROR: No proteins in query FASTA")

outputs = {}

for name in (
    "swissprot.tsv",
    "phrog.native.tsv",
    "pfam.domtblout",
):
    path = out / name

    if not path.is_file():
        raise SystemExit(
            f"ERROR: Missing search output: {name}"
        )

    outputs[name] = {
        "size": path.stat().st_size,
        "sha256": sha256(path),
    }

manifest = {
    "schema_version": 1,
    "status": "COMPLETE",
    "finished_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "query": {
        "path": str(query.resolve()),
        "sha256": sha256(query),
        "sequence_count": sequence_count,
    },
    "databases": {
        "swissprot": fingerprint(swiss),
        "phrog": fingerprint(phrog),
        "pfam": fingerprint(pfam),
    },
    "parameters": {
        "threads": threads,
        "phrog_max_seqs": max_seqs,
    },
    "outputs": outputs,
}

out.mkdir(parents=True, exist_ok=True)

fd, tmp_name = tempfile.mkstemp(
    prefix=".search_run_complete.",
    suffix=".tmp",
    dir=out,
)

try:
    with os.fdopen(fd, "w") as f:
        json.dump(
            manifest,
            f,
            indent=2,
            sort_keys=True
        )
        f.write("\n")

    os.replace(
        tmp_name,
        out / "search_run_complete.json"
    )
finally:
    if os.path.exists(tmp_name):
        os.unlink(tmp_name)

print("SEARCH_MANIFEST=PASS")
print(f"VERIFIED_PROTEINS={sequence_count}")
print("SEARCH_STATUS=COMPLETE")
PY_MANIFEST

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
