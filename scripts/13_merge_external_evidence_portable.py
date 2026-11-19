#!/usr/bin/env python3
"""Merge preliminary Stage 13 and verified post-ranking K2 evidence.

This merge does not calculate Stage 10 scores or select candidates.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


PRELIM_COLUMNS = [
    "structure_id",
    "priority_class",
    "experimental_score",
    "dominant_context",
    "assay_category",
    "independent_context_n",
    "independent_context_sample_span",
    "single_mean_plddt",
    "msa_mean_plddt",
    "delta_mean_plddt",
    "single_ptm",
    "msa_ptm",
    "delta_ptm",
    "external_status",
    "raw_sources",
    "convincing_sources",
    "BFVD_best_target",
    "BFVD_best_qtm",
    "BFVD_best_ttm",
    "very_conserved_positions",
    "priority_structural_positions",
    "priority_acidic",
    "priority_basic",
    "spatial_cluster_count",
    "spatial_clusters",
]

EXTERNAL_COLUMNS = [
    "structure_id",
    "external_status",
    "raw_sources",
    "convincing_sources",
    "BFVD_best_target",
    "BFVD_best_qtm",
    "BFVD_best_ttm",
    "BFVD_best_is_convincing",
    "msa_pdb_sha256",
]

EXTERNAL_FIELDS = [
    name for name in EXTERNAL_COLUMNS
    if name != "structure_id"
]

ADDED_FIELDS = [
    name for name in EXTERNAL_FIELDS
    if name not in PRELIM_COLUMNS
]

OUTPUT_COLUMNS = PRELIM_COLUMNS + ADDED_FIELDS


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("PROVENANCE_NOT_OBJECT")
    return data


def read_table(path, expected_fields):
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)

    if fields != expected_fields:
        raise ValueError(
            f"TABLE_SCHEMA_MISMATCH:{path.name}"
        )

    if not rows:
        raise ValueError(f"EMPTY_TABLE:{path.name}")

    if any(
        None in row or any(v is None for v in row.values())
        for row in rows
    ):
        raise ValueError(f"MALFORMED_TABLE:{path.name}")

    identifiers = [row["structure_id"] for row in rows]

    if (
        any(not sid.strip() for sid in identifiers)
        or len(set(identifiers)) != len(identifiers)
    ):
        raise ValueError(
            f"INVALID_OR_DUPLICATE_STRUCTURE_ID:{path.name}"
        )

    return rows


def validate_external(row):
    status = row["external_status"]
    raw = row["raw_sources"]
    convincing = row["convincing_sources"]

    if not status or status == "NOT_SEARCHED":
        raise ValueError("EXTERNAL_SEARCH_NOT_COMPLETED")

    if status == "NO_RAW_HIT":
        if raw or convincing:
            raise ValueError("INCONSISTENT_NO_RAW_HIT")

    if status == "WEAK_RAW_HIT_ONLY":
        if not raw or convincing:
            raise ValueError("INCONSISTENT_WEAK_RAW_HIT")

    flag = row["BFVD_best_is_convincing"]

    if flag not in ("0", "1"):
        raise ValueError("INVALID_BFVD_FLAG")

    if flag == "1" and "BFVD" not in convincing.upper():
        raise ValueError("INCONSISTENT_BFVD_EVIDENCE")

    digest = row["msa_pdb_sha256"]

    if (
        len(digest) != 64
        or any(x not in "0123456789abcdef" for x in digest)
    ):
        raise ValueError("INVALID_MSA_PDB_SHA256")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--preliminary", type=Path, required=True
    )
    ap.add_argument(
        "--preliminary-provenance",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--external", type=Path, required=True
    )
    ap.add_argument(
        "--external-provenance",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--out-dir", type=Path, required=True
    )

    args = ap.parse_args()

    input_paths = [
        args.preliminary.resolve(),
        args.preliminary_provenance.resolve(),
        args.external.resolve(),
        args.external_provenance.resolve(),
    ]

    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            f"Refusing existing output directory: {out}"
        )

    for path in input_paths:
        if not path.is_file():
            raise ValueError(f"MISSING_INPUT:{path}")

        if out == path or out in path.parents:
            raise ValueError("OUTPUT_OVERLAPS_INPUT")

    pre_rows = read_table(
        args.preliminary, PRELIM_COLUMNS
    )

    ext_rows = read_table(
        args.external, EXTERNAL_COLUMNS
    )

    pre_meta = read_json(args.preliminary_provenance)
    ext_meta = read_json(args.external_provenance)

    if ext_meta.get("control_only") is not False:
        raise ValueError("CONTROL_ONLY_EVIDENCE_REJECTED")

    if ext_meta.get("final_selection_performed") is not False:
        raise ValueError("UNEXPECTED_FINAL_SELECTION")

    if ext_meta.get("stage13_evidence_sha256") != sha(
        args.external
    ):
        raise ValueError(
            "EXTERNAL_TABLE_CHECKSUM_MISMATCH"
        )

    if pre_meta.get("matrix_sha256") != sha(args.preliminary):
        raise ValueError(
            "PRELIMINARY_CHECKSUM_MISMATCH:matrix_sha256"
        )

    for key in (
        "preliminary_evidence_sha256",
        "evidence_sha256",
    ):
        recorded = pre_meta.get(key)

        if recorded is not None:
            if recorded != sha(args.preliminary):
                raise ValueError(
                    f"PRELIMINARY_CHECKSUM_MISMATCH:{key}"
                )

    pre_ids = {
        row["structure_id"] for row in pre_rows
    }

    ext_ids = {
        row["structure_id"] for row in ext_rows
    }

    if pre_ids != ext_ids:
        raise ValueError(
            "FAMILY_ID_SET_MISMATCH:"
            f"missing_external={sorted(pre_ids-ext_ids)}:"
            f"extra_external={sorted(ext_ids-pre_ids)}"
        )

    ext_index = {
        row["structure_id"]: row
        for row in ext_rows
    }

    merged = []

    for row in pre_rows:
        sid = row["structure_id"]
        external = ext_index[sid]

        validate_external(external)

        if row["external_status"] != "NOT_SEARCHED":
            raise ValueError(
                f"PRELIMINARY_ALREADY_HAS_EXTERNAL_EVIDENCE:{sid}"
            )

        merged_row = dict(row)

        for field in EXTERNAL_FIELDS:
            merged_row[field] = external[field]

        if set(merged_row) != set(OUTPUT_COLUMNS):
            raise ValueError("MERGED_SCHEMA_MISMATCH")

        merged.append(merged_row)

    out.mkdir(parents=True, exist_ok=False)

    output = (
        out / "preliminary_candidate_evidence_with_external.tsv"
    )

    with output.open(
        "w", newline="", encoding="utf-8"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(merged)

    status_counts = Counter(
        row["external_status"] for row in merged
    )

    provenance = {
        "stage": "STAGE13_POST_RANKING_EXTERNAL_MERGE",
        "input_preliminary": str(
            args.preliminary.resolve()
        ),
        "input_preliminary_sha256": sha(
            args.preliminary
        ),
        "input_preliminary_provenance_sha256": sha(
            args.preliminary_provenance
        ),
        "input_external": str(
            args.external.resolve()
        ),
        "input_external_sha256": sha(
            args.external
        ),
        "input_external_provenance_sha256": sha(
            args.external_provenance
        ),
        "external_control_only": False,
        "stage10_ranking_recalculated": False,
        "final_candidate_selection_performed": False,
        "families": len(merged),
        "columns": len(OUTPUT_COLUMNS),
        "status_counts": dict(status_counts),
        "output_sha256": sha(output),
        "merge_policy": (
            "REPLACE_ONLY_NOT_SEARCHED_FIELDS;"
            "PRESERVE_NONEXTERNAL_PRELIMINARY_EVIDENCE"
        ),
    }

    (
        out / "preliminary_external_merge_provenance.json"
    ).write_text(
        json.dumps(
            provenance,
            indent=2,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )

    print(f"FAMILIES={len(merged)}")
    print(f"COLUMNS={len(OUTPUT_COLUMNS)}")
    print(f"EXTERNAL_STATUS_COUNTS={dict(status_counts)}")
    print(f"OUTPUT={output}")
    print("STAGE13_EXTERNAL_MERGE=PASS")


if __name__ == "__main__":
    main()
