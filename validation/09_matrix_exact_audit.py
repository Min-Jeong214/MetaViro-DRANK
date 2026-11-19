#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BASE = (
    ROOT / "Result/Annotation_saturation/Pilot_v2"
)

CANON = (
    BASE / "StageF_Structure/StageG3_ExternalDB"
    / "Integrated_ExternalStructure"
    / "structure191_external_structure_matrix.tsv"
)

DEFAULT_CALC = (
    BASE / "Stage09_Reconstruction_DryRun"
    / "structure191_external_structure_matrix.tsv"
)

def read(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return reader.fieldnames, list(reader)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--calc",
        type=Path,
        default=DEFAULT_CALC
    )
    args = ap.parse_args()

    canon_cols, canon = read(CANON)
    calc_cols, calc = read(args.calc)

    print("===== STAGE09 FINAL AUDIT =====")

    print("CANON_ROWS=", len(canon))
    print("CALC_ROWS=", len(calc))
    print("CANON_COLS=", len(canon_cols))
    print("CALC_COLS=", len(calc_cols))

    schema = (
        len(canon) == len(calc) == 191
        and len(canon_cols) == len(calc_cols) == 39
        and canon_cols == calc_cols
    )

    print("SCHEMA_EXACT=", schema)

    if not schema:
        raise SystemExit("STAGE09_STATUS=FAIL")

    canon_ids = [r["structure_id"] for r in canon]
    calc_ids = [r["structure_id"] for r in calc]

    order_exact = canon_ids == calc_ids

    print("ORDER_EXACT=", order_exact)

    if (
        len(set(canon_ids)) != 191
        or len(set(calc_ids)) != 191
        or set(canon_ids) != set(calc_ids)
    ):
        raise SystemExit("STAGE09_ID_SET=FAIL")

    calc_map = {
        r["structure_id"]: r for r in calc
    }

    mismatch = Counter()

    for original in canon:
        generated = calc_map[original["structure_id"]]

        for column in canon_cols:
            if original[column] != generated[column]:
                mismatch[column] += 1

    print("MISMATCHED_FIELDS=", len(mismatch))
    print("MISMATCHED_CELLS=", sum(mismatch.values()))

    for field, count in mismatch.most_common():
        print("DIFF", field, count)

    byte_exact = (
        CANON.read_bytes() == args.calc.read_bytes()
    )

    print("BYTE_EXACT=", byte_exact)

    counts = Counter(
        r["external_status"] for r in calc
    )

    print("===== STATUS DISTRIBUTION =====")

    expected = {
        "CONVINCING_STRUCTURE_CONNECTION": 4,
        "WEAK_RAW_HIT_ONLY": 138,
        "NO_RAW_HIT_E1": 49,
    }

    for status, number in expected.items():
        print(
            status,
            "CALC=", counts[status],
            "CANON=", number
        )

    passed = (
        schema
        and order_exact
        and not mismatch
        and byte_exact
        and counts == expected
    )

    print("===== FINAL =====")
    print(
        "STAGE09_STATUS=",
        "PASS" if passed else "FAIL"
    )

    if not passed:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
