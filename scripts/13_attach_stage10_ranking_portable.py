#!/usr/bin/env python3
"""Attach verified Stage 10 rankings to post-ranking Stage 13 evidence.

Keeps initial-structure Stage 10 external evidence separate
from MSA-structure Stage K2 external evidence.

Does not rerank candidates or perform final candidate selection.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ADD_COLUMNS = [
    "I2_integrated_rank",
    "I2_experimental_rank",
    "I2_integrated_score",
    "I2_initial_external_status",
    "I2_initial_raw_sources",
    "I2_initial_convincing_sources",
]

MAPPING = {
    "priority_class": "priority_class",
    "experimental_score": "I2_experimental_score",
    "assay_category": "I2_assay_category",
    "independent_context_n": "I2_assay_independent_members",
    "independent_context_sample_span":
        "I2_assay_independent_sample_span",
}

RANK_REQUIRED = set(MAPPING.values()) | {
    "structure_id",
    "family_rep",
    "I2_integrated_rank",
    "I2_experimental_rank",
    "I2_integrated_score",
    "external_status",
    "raw_sources",
    "convincing_sources",
}


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path):
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("INVALID_PROVENANCE_OBJECT")
    return obj


def load_tsv(path):
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if not fields or len(fields) != len(set(fields)):
        raise ValueError("INVALID_TABLE_HEADER")

    if not rows:
        raise ValueError("EMPTY_TABLE")

    if any(
        None in row or any(value is None for value in row.values())
        for row in rows
    ):
        raise ValueError("MALFORMED_TABLE")

    if "structure_id" not in fields:
        raise ValueError("MISSING_STRUCTURE_ID")

    ids = [row["structure_id"] for row in rows]

    if (
        any(not sid for sid in ids)
        or len(ids) != len(set(ids))
    ):
        raise ValueError("DUPLICATE_OR_EMPTY_STRUCTURE_ID")

    return fields, rows


def main():
    parser = argparse.ArgumentParser()

    for option in (
        "ranking",
        "ranking_provenance",
        "shortlist",
        "selection",
        "selection_provenance",
        "preliminary",
        "preliminary_provenance",
        "merged",
        "merged_provenance",
        "out_dir",
    ):
        parser.add_argument(
            "--" + option.replace("_", "-"),
            type=Path,
            required=True,
        )

    args = parser.parse_args()

    sources = [
        args.ranking,
        args.ranking_provenance,
        args.shortlist,
        args.selection,
        args.selection_provenance,
        args.preliminary,
        args.preliminary_provenance,
        args.merged,
        args.merged_provenance,
    ]

    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            f"Refusing existing output directory: {out}"
        )

    for path in sources:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"MISSING_INPUT:{path}")

        resolved = path.resolve()

        if resolved == out or out in resolved.parents:
            raise ValueError("OUTPUT_INPUT_OVERLAP")

    rank_meta = load_json(args.ranking_provenance)
    sel_meta = load_json(args.selection_provenance)
    pre_meta = load_json(args.preliminary_provenance)
    merge_meta = load_json(args.merged_provenance)

    rank_sha = sha(args.ranking)
    select_sha = sha(args.selection)
    prelim_sha = sha(args.preliminary)
    merged_sha = sha(args.merged)

    ranking_hashes = [
        value for key, value in rank_meta.items()
        if key.endswith("sha256")
        and isinstance(value, str)
    ]

    if rank_sha not in ranking_hashes:
        raise ValueError("RANKING_SOURCE_SHA256_MISMATCH")

    if (
        sel_meta.get("selection_sha256") != select_sha
        or sel_meta.get("source_ranking_sha256")
        != rank_sha
        or sel_meta.get("source_shortlist_sha256")
        != sha(args.shortlist)
        or sel_meta.get("selection_method")
        != "VERIFIED_STAGE10_TOP15"
    ):
        raise ValueError("SELECTION_PROVENANCE_MISMATCH")

    if (
        pre_meta.get("selection_sha256") != select_sha
        or pre_meta.get("matrix_sha256") != prelim_sha
        or pre_meta.get("stage10_ranking_available") is not False
    ):
        raise ValueError("STAGE13_PRELIMINARY_PROVENANCE_MISMATCH")

    if (
        merge_meta.get("output_sha256") != merged_sha
        or merge_meta.get("input_preliminary_sha256")
        != prelim_sha
        or merge_meta.get("external_control_only") is not False
        or merge_meta.get("stage10_ranking_recalculated") is not False
        or merge_meta.get("final_candidate_selection_performed") is not False
    ):
        raise ValueError("STAGE13_MERGE_PROVENANCE_MISMATCH")

    rank_fields, rankings = load_tsv(args.ranking)
    short_fields, shortlist = load_tsv(args.shortlist)
    selection_fields, selection = load_tsv(args.selection)
    pre_fields, prelim_rows = load_tsv(args.preliminary)
    merged_fields, merged_rows = load_tsv(args.merged)

    if not RANK_REQUIRED.issubset(rank_fields):
        raise ValueError("MISSING_RANKING_FIELDS")

    if rank_fields != short_fields:
        raise ValueError("SHORTLIST_SCHEMA_MISMATCH")

    if selection_fields != ["structure_id"]:
        raise ValueError("SELECTION_SCHEMA_MISMATCH")

    if len(merged_fields) != 27:
        raise ValueError("UNEXPECTED_STAGE13_INPUT_SCHEMA")

    if any(name in merged_fields for name in ADD_COLUMNS):
        raise ValueError("RANKING_FIELDS_ALREADY_PRESENT")

    rank_by_id = {
        row["structure_id"]: row for row in rankings
    }

    short_by_id = {
        row["structure_id"]: row for row in shortlist
    }

    selected_ids = [
        row["structure_id"] for row in selection
    ]

    if (
        set(selected_ids) != {
            row["structure_id"] for row in merged_rows
        }
        or set(selected_ids) != {
            row["structure_id"] for row in prelim_rows
        }
    ):
        raise ValueError("STAGE10_STAGE13_SELECTION_MISMATCH")

    if sel_meta.get("selected_ids") != selected_ids:
        raise ValueError("SELECTED_ID_PROVENANCE_MISMATCH")

    for sid, row in short_by_id.items():
        if sid not in rank_by_id:
            raise ValueError("SHORTLIST_ID_NOT_IN_RANKING")

        if row != rank_by_id[sid]:
            raise ValueError("SHORTLIST_RANKING_ROW_MISMATCH")

    for sid in selected_ids:
        if sid not in short_by_id:
            raise ValueError("SELECTED_ID_NOT_IN_SHORTLIST")

    output_rows = []

    for original in merged_rows:
        sid = original["structure_id"]
        ranking = rank_by_id[sid]

        for field in MAPPING:
            if field not in original:
                raise ValueError(
                    f"MISSING_STAGE13_FIELD:{field}"
                )

        for field, source in MAPPING.items():
            existing = original[field]
            expected = ranking[source]

            if existing not in ("", expected):
                raise ValueError(
                    f"CONFLICTING_STAGE13_RANKING_FIELD:{field}"
                )

        for name in (
            "I2_integrated_rank",
            "I2_experimental_rank",
        ):
            if int(ranking[name]) < 1:
                raise ValueError("INVALID_RANK")

        for name in (
            "I2_integrated_score",
            "I2_experimental_score",
        ):
            if not math.isfinite(float(ranking[name])):
                raise ValueError("INVALID_SCORE")

        updated = dict(original)

        for target, source in MAPPING.items():
            updated[target] = ranking[source]

        updated.update({
            "I2_integrated_rank":
                ranking["I2_integrated_rank"],
            "I2_experimental_rank":
                ranking["I2_experimental_rank"],
            "I2_integrated_score":
                ranking["I2_integrated_score"],
            "I2_initial_external_status":
                ranking["external_status"],
            "I2_initial_raw_sources":
                ranking["raw_sources"],
            "I2_initial_convincing_sources":
                ranking["convincing_sources"],
        })

        output_rows.append(updated)

    output_fields = merged_fields + ADD_COLUMNS

    out.mkdir(parents=True, exist_ok=False)

    output = out / "ranked_candidate_evidence.tsv"

    with output.open(
        "x", newline="", encoding="utf-8"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    provenance = {
        "stage": "STAGE13_VERIFIED_STAGE10_RANKING_HANDOFF",
        "ranking_sha256": rank_sha,
        "ranking_provenance_sha256":
            sha(args.ranking_provenance),
        "shortlist_sha256": sha(args.shortlist),
        "selection_sha256": select_sha,
        "selection_provenance_sha256":
            sha(args.selection_provenance),
        "preliminary_sha256": prelim_sha,
        "preliminary_provenance_sha256":
            sha(args.preliminary_provenance),
        "merged_input_sha256": merged_sha,
        "merged_provenance_sha256":
            sha(args.merged_provenance),
        "output_sha256": sha(output),
        "families": len(output_rows),
        "columns": len(output_fields),
        "initial_external_source": "STAGE09_SINGLE_STRUCTURE",
        "post_ranking_external_source": "STAGE_K2_MSA_STRUCTURE",
        "stage10_scores_recalculated": False,
        "final_candidate_selection_performed": False,
        "historical_top5_selected": False,
        "historical_top3_selected": False,
        "dominant_context_policy":
            "PRESERVED_UNCHANGED_WITHOUT_VERIFIED_MAPPING",
    }

    with (
        out / "ranked_candidate_evidence_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(provenance, fh, indent=2)
        fh.write("\n")

    print("RANKING_SOURCE_PROVENANCE=PASS")
    print("SHORTLIST_SELECTION_PROVENANCE=PASS")
    print("STAGE13_SOURCE_PROVENANCE=PASS")
    print(f"OUTPUT_FAMILIES={len(output_rows)}")
    print(f"OUTPUT_COLUMNS={len(output_fields)}")
    print(f"OUTPUT={output}")
    print("STAGE10_TO_STAGE13_RANKING_HANDOFF=PASS")


if __name__ == "__main__":
    main()
