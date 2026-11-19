#!/usr/bin/env python3

import csv
import hashlib
import io
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Result/Annotation_saturation/Pilot_v2"

CANON = (
    BASE / "StageI2_ConfidenceAwareRanking"
    / "recurrent191_confidence_aware_ranking.tsv"
)

OUT = (
    ROOT / "scripts/resources"
    / "10_i2_historical_tie_order.tsv"
)

README = (
    ROOT / "scripts/resources"
    / "10_i2_historical_tie_order.README.md"
)

with CANON.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 191
assert len({r["structure_id"] for r in rows}) == 191

rows.sort(key=lambda r: int(r["I2_integrated_rank"]))

assert [
    int(r["I2_integrated_rank"]) for r in rows
] == list(range(1, 192))

def pair(r):
    return (
        float(r["I2_integrated_score"]),
        int(r["noncandidate_context_sample_span"]),
    )

for a, b in zip(rows, rows[1:]):
    if pair(a) < pair(b):
        raise RuntimeError(
            "Canonical ranking does not obey score/span order"
        )

groups = defaultdict(list)

for row in rows:
    groups[pair(row)].append(row["structure_id"])

tied_groups = {
    key: ids
    for key, ids in groups.items()
    if len(ids) > 1
}

tied_ids = sum(len(ids) for ids in tied_groups.values())

print("TIE_GROUPS =", len(tied_groups))
print("FAMILIES_WITH_TIES =", tied_ids)

if len(tied_groups) != 34 or tied_ids != 168:
    raise RuntimeError("Unexpected historical tie structure")

position = {}

for group, ids in tied_groups.items():
    for rank, sid in enumerate(ids, 1):
        position[sid] = rank

columns = [
    "structure_id",
    "I2_integrated_score",
    "noncandidate_context_sample_span",
    "historical_tie_order",
]

buffer = io.StringIO(newline="")

writer = csv.DictWriter(
    buffer,
    fieldnames=columns,
    delimiter="\t",
    lineterminator="\n",
)

writer.writeheader()

for row in rows:
    sid = row["structure_id"]

    if sid not in position:
        continue

    writer.writerow({
        "structure_id": sid,
        "I2_integrated_score": row["I2_integrated_score"],
        "noncandidate_context_sample_span":
            row["noncandidate_context_sample_span"],
        "historical_tie_order": position[sid],
    })

content = buffer.getvalue().encode()

OUT.parent.mkdir(parents=True, exist_ok=True)

if OUT.exists():
    if OUT.read_bytes() != content:
        raise RuntimeError(
            "Existing tie manifest differs; refusing overwrite"
        )
    print("TIE_MANIFEST=EXISTING_EXACT")
else:
    OUT.write_bytes(content)
    print("TIE_MANIFEST=CREATED")

digest = hashlib.sha256(CANON.read_bytes()).hexdigest()

readme_text = f"""# Stage 10 historical tie-order provenance

This file accompanies the MetaViro Stage I2 ranking script.

## Provenance

Source canonical table:
`Result/Annotation_saturation/Pilot_v2/StageI2_ConfidenceAwareRanking/recurrent191_confidence_aware_ranking.tsv`

Source SHA256: `{digest}`

The manifest contains {tied_ids} structures belonging to
{len(tied_groups)} groups tied on both I2 integrated score and
noncandidate context sample span.

## Interpretation

- All I2 scores are calculated independently from Stage I and H3.
- Experimental ranks are calculated independently.
- Integrated rank uses descending integrated score, then
  descending noncandidate context sample span.
- The historical ordering within exact ties is retained
  from the frozen canonical output.
- The original historical tie-breaking implementation
  has NOT been independently recovered.
- This manifest is archival provenance, not biological evidence.

## Modes

Without `--historical-tie-order`, the generator uses
deterministic structure ID tie-breaking.

With `--historical-tie-order`, the generator reproduces
the frozen historical order using this manifest.

Do not use the canonical ranking table as a runtime input.
"""

if README.exists():
    if README.read_text() != readme_text:
        raise RuntimeError(
            "Existing provenance README differs"
        )
else:
    README.write_text(readme_text)

print("MANIFEST =", OUT)
print("MANIFEST_ROWS =", tied_ids)
print("CANONICAL_SHA256 =", digest)
print("HISTORICAL_TIE_FREEZE=PASS")
