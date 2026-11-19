#!/usr/bin/env python3
"""Validate independently reconstructed Stage I2 outputs."""

import argparse
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BASE = (
    ROOT / "Result/Annotation_saturation/Pilot_v2"
)

CANON_DIR = BASE / "StageI2_ConfidenceAwareRanking"

FILES = {
    "MAIN_191": (
        "recurrent191_confidence_aware_ranking.tsv",
        191,
        "I2_integrated_rank",
    ),
    "TOP20": (
        "top20_confidence_aware.tsv",
        20,
        "I2_integrated_rank",
    ),
    "TOP15": (
        "top15_experimental_candidates.tsv",
        15,
        "I2_experimental_rank",
    ),
}


def read(path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_one(label, filename, expected_rows, rank_field, out_dir):
    canonical = CANON_DIR / filename
    generated = out_dir / filename

    print(f"\n===== {label} =====")

    print("CANON =", canonical)
    print("CALC  =", generated)

    if not canonical.is_file():
        raise FileNotFoundError(canonical)

    if not generated.is_file():
        raise FileNotFoundError(generated)

    canon_cols, canon_rows = read(canonical)
    calc_cols, calc_rows = read(generated)

    schema = (
        len(canon_cols) == len(calc_cols) == 51
        and len(canon_rows) == len(calc_rows) == expected_rows
        and canon_cols == calc_cols
    )

    print("CANON_ROWS =", len(canon_rows))
    print("CALC_ROWS =", len(calc_rows))
    print("CANON_COLUMNS =", len(canon_cols))
    print("CALC_COLUMNS =", len(calc_cols))
    print("SCHEMA_EXACT =", schema)

    if not schema:
        print("CANON_HEADER =", canon_cols)
        print("CALC_HEADER =", calc_cols)
        return False

    canon_ids = [r["structure_id"] for r in canon_rows]
    calc_ids = [r["structure_id"] for r in calc_rows]

    id_set = (
        len(set(canon_ids)) == expected_rows
        and len(set(calc_ids)) == expected_rows
        and set(canon_ids) == set(calc_ids)
    )

    order_exact = canon_ids == calc_ids

    print("ID_SET_EXACT =", id_set)
    print("ORDER_EXACT =", order_exact)

    differences = []

    if id_set:
        calc_map = {
            r["structure_id"]: r
            for r in calc_rows
        }

        for expected in canon_rows:
            sid = expected["structure_id"]
            observed = calc_map[sid]

            for column in canon_cols:
                if expected[column] != observed[column]:
                    differences.append(
                        (
                            sid,
                            column,
                            expected[column],
                            observed[column],
                        )
                    )

    print("MISMATCHED_CELLS =", len(differences))

    for sid, col, expected, actual in differences[:12]:
        print(
            "DIFF",
            sid,
            col,
            "CANON=", repr(expected),
            "CALC=", repr(actual),
        )

    if id_set:
        correct_ranks = (
            [int(r[rank_field]) for r in calc_rows]
            == list(range(1, expected_rows + 1))
        )
    else:
        correct_ranks = False

    print("RANK_ORDER_VALID =", correct_ranks)

    byte_exact = (
        canonical.read_bytes() == generated.read_bytes()
    )

    print("BYTE_EXACT =", byte_exact)

    print("CANON_SHA256 =", sha256(canonical))
    print("CALC_SHA256  =", sha256(generated))

    passed = (
        schema
        and id_set
        and order_exact
        and not differences
        and correct_ranks
        and byte_exact
    )

    print(
        f"{label}_STATUS =",
        "PASS" if passed else "FAIL"
    )

    return passed


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )

    args = parser.parse_args()

    out_dir = args.out_dir.resolve()

    if out_dir == CANON_DIR.resolve():
        raise ValueError(
            "Audit must use reconstructed output directory"
        )

    results = {}

    for label, (filename, nrows, rank_field) in FILES.items():
        results[label] = audit_one(
            label,
            filename,
            nrows,
            rank_field,
            out_dir,
        )

    print("\n===== STAGE 10 FINAL SUMMARY =====")

    for label, passed in results.items():
        print(
            label,
            "PASS" if passed else "FAIL"
        )

    final = all(results.values())

    print(
        "STAGE10_ALL_OUTPUTS_STATUS =",
        "PASS" if final else "FAIL"
    )

    print("CANONICAL_I2_UNMODIFIED=YES")

    if not final:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
