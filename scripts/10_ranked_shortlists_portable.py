#!/usr/bin/env python3
"""Generate portable Stage I2 shortlists from verified rankings."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCORING = ROOT / "scripts/10_integrated_ranking.py"

REQUIRED = {
    "structure_id",
    "family_rep",
    "I2_integrated_rank",
    "I2_experimental_rank",
    "I2_integrated_score",
    "I2_experimental_score",
    "I2_assay_independent_sample_span",
    "noncandidate_context_sample_span",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path):
    if not path.is_file():
        raise ValueError(f"Missing ranking: {path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = reader.fieldnames or []
        rows = list(reader)

    if len(columns) != len(set(columns)):
        raise ValueError("Duplicate ranking columns")

    if not REQUIRED.issubset(columns):
        raise ValueError(
            f"Missing ranking fields: {sorted(REQUIRED-set(columns))}"
        )

    if not rows:
        raise ValueError("Empty ranking")

    return columns, rows


def write_table(path, columns, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def score(row, field):
    result = float(row[field])
    if not math.isfinite(result):
        raise ValueError(f"Invalid score: {field}")
    return result


def count(row, field):
    value = int(row[field])
    if value < 0:
        raise ValueError(f"Negative count: {field}")
    return value


def validate(rows, provenance):
    n = len(rows)

    ids = [r["structure_id"] for r in rows]
    reps = [r["family_rep"] for r in rows]

    if (
        not all(ids)
        or not all(reps)
        or len(set(ids)) != n
        or len(set(reps)) != n
    ):
        raise ValueError("Duplicate/empty structure or Family ID")

    for field in (
        "I2_integrated_rank",
        "I2_experimental_rank",
    ):
        ranks = [int(r[field]) for r in rows]

        if sorted(ranks) != list(range(1, n + 1)):
            raise ValueError(
                f"Invalid or incomplete ranking: {field}"
            )

    for row in rows:
        score(row, "I2_integrated_score")
        score(row, "I2_experimental_score")
        count(row, "noncandidate_context_sample_span")
        count(row, "I2_assay_independent_sample_span")

    integrated = sorted(
        rows,
        key=lambda r: int(r["I2_integrated_rank"]),
    )

    experimental = sorted(
        rows,
        key=lambda r: int(r["I2_experimental_rank"]),
    )

    # Preserve the historical Stage I2 ordering rules.
    if provenance.get("historical_tie_replay"):
        # Historical tie order may differ from structure ID order.
        integrated_keys = [
            (
                -score(r, "I2_integrated_score"),
                -count(r, "noncandidate_context_sample_span"),
            )
            for r in integrated
        ]
        if integrated_keys != sorted(integrated_keys):
            raise ValueError("Integrated ranking order mismatch")
    else:
        expected = sorted(
            rows,
            key=lambda r: (
                -score(r, "I2_integrated_score"),
                -count(r, "noncandidate_context_sample_span"),
                r["structure_id"],
            ),
        )
        if [r["structure_id"] for r in expected] != [
            r["structure_id"] for r in integrated
        ]:
            raise ValueError("Integrated ranking order mismatch")

    expected_exp = sorted(
        rows,
        key=lambda r: (
            -score(r, "I2_experimental_score"),
            -score(r, "I2_integrated_score"),
            -count(
                r, "I2_assay_independent_sample_span"
            ),
            r["structure_id"],
        ),
    )

    if [r["structure_id"] for r in expected_exp] != [
        r["structure_id"] for r in experimental
    ]:
        raise ValueError("Experimental ranking order mismatch")

    return integrated, experimental


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=Path, required=True)
    parser.add_argument(
        "--ranking-provenance", type=Path, required=True
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    ranking = args.ranking.resolve()
    record_file = args.ranking_provenance.resolve()
    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            f"Refusing existing output directory: {output}"
        )

    columns, rows = read_table(ranking)

    if not record_file.is_file():
        raise ValueError("Missing ranking provenance")

    provenance = json.loads(record_file.read_text())

    if provenance.get("algorithm") != "ORIGINAL_STAGE_I2":
        raise ValueError("Unexpected ranking algorithm")

    if provenance.get("ranking_sha256") != sha256(ranking):
        raise ValueError("Ranking SHA-256 mismatch")

    if provenance.get("scoring_script_sha256") != sha256(
        SCORING
    ):
        raise ValueError("Original I2 scoring version mismatch")

    if (
        provenance.get("families") != len(rows)
        or provenance.get("columns") != len(columns)
    ):
        raise ValueError("Ranking provenance shape mismatch")

    integrated, experimental = validate(rows, provenance)

    top20 = integrated[:20]
    top15 = experimental[:15]

    # Create outputs only after successful validation.
    output.mkdir(parents=True)

    top20_path = output / "top20_confidence_aware.tsv"
    top15_path = output / "top15_experimental_candidates.tsv"

    write_table(top20_path, columns, top20)
    write_table(top15_path, columns, top15)

    result = {
        "source": "VERIFIED_STAGE10_RANKING",
        "ranking_sha256": sha256(ranking),
        "ranking_provenance_sha256": sha256(record_file),
        "source_families": len(rows),
        "source_columns": len(columns),
        "top20_count": len(top20),
        "top15_count": len(top15),
        "top20_sha256": sha256(top20_path),
        "top15_sha256": sha256(top15_path),
        "historical_tie_replay": bool(
            provenance.get("historical_tie_replay")
        ),
        "final_top5_selected": False,
        "primary_top3_selected": False,
    }

    (output / "shortlist_provenance.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )

    print(f"INPUT_FAMILIES={len(rows)}")
    print(f"INPUT_COLUMNS={len(columns)}")
    print(f"TOP20_ROWS={len(top20)}")
    print(f"TOP15_ROWS={len(top15)}")
    print(
        "TOP20_IDS="
        + ",".join(r["structure_id"] for r in top20)
    )
    print(
        "TOP15_IDS="
        + ",".join(r["structure_id"] for r in top15)
    )
    print("FINAL_TOP5=NOT_SELECTED")
    print("PRIMARY_TOP3=NOT_SELECTED")
    print("STAGE10_PORTABLE_SHORTLIST=PASS")


if __name__ == "__main__":
    main()
