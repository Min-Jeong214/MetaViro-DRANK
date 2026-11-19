#!/usr/bin/env python3
"""Generate Stage I2 ranked shortlists from reconstructed ranking.

Does not read canonical I2 results.
"""

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CANON_DIR = (
    ROOT / "Result/Annotation_saturation/Pilot_v2"
    / "StageI2_ConfidenceAwareRanking"
)

OUTPUT_NAMES = (
    "top20_confidence_aware.tsv",
    "top15_experimental_candidates.tsv",
)


def read_table(path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def validate_ranks(rows, field):
    ranks = [int(row[field]) for row in rows]

    if sorted(ranks) != list(range(1, 192)):
        raise ValueError(
            f"{field}: ranks must be unique integers 1..191"
        )


def write_table(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ranking",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )

    args = parser.parse_args()

    ranking = args.ranking.resolve()
    out_dir = args.out_dir.resolve()
    canonical = CANON_DIR.resolve()

    if ranking == (
        canonical / "recurrent191_confidence_aware_ranking.tsv"
    ):
        raise ValueError(
            "Canonical I2 cannot be used as generator input"
        )

    if out_dir == canonical or canonical in out_dir.parents:
        raise ValueError(
            "Refusing to write into canonical I2 directory"
        )

    for filename in OUTPUT_NAMES:
        if out_dir / filename == ranking:
            raise ValueError(
                "Output would overwrite ranking input"
            )

    columns, rows = read_table(ranking)

    if len(rows) != 191 or len(columns) != 51:
        raise ValueError(
            f"Unexpected ranking shape: "
            f"{len(rows)} x {len(columns)}"
        )

    ids = [r["structure_id"] for r in rows]

    if len(set(ids)) != 191:
        raise ValueError("Duplicate structure IDs")

    validate_ranks(rows, "I2_integrated_rank")
    validate_ranks(rows, "I2_experimental_rank")

    integrated = sorted(
        rows,
        key=lambda r: int(r["I2_integrated_rank"])
    )

    experimental = sorted(
        rows,
        key=lambda r: int(r["I2_experimental_rank"])
    )

    top20 = integrated[:20]
    top15 = experimental[:15]

    top20_out = out_dir / OUTPUT_NAMES[0]
    top15_out = out_dir / OUTPUT_NAMES[1]

    write_table(top20_out, columns, top20)
    write_table(top15_out, columns, top15)

    print("===== GENERATED SHORTLISTS =====")

    print("TOP20_FILE =", top20_out)
    print("TOP20_ROWS =", len(top20))
    print("TOP20_COLUMNS =", len(columns))

    print(
        "TOP20_IDS =",
        ",".join(r["structure_id"] for r in top20)
    )

    print("TOP15_FILE =", top15_out)
    print("TOP15_ROWS =", len(top15))
    print("TOP15_COLUMNS =", len(columns))

    print(
        "TOP15_IDS =",
        ",".join(r["structure_id"] for r in top15)
    )

    print("STAGE10_SHORTLIST_GENERATION=COMPLETE")


if __name__ == "__main__":
    main()
