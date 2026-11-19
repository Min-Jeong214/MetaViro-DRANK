#!/usr/bin/env python3
"""Prepare portable Stage 10 evidence without changing I2 scoring."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


EVIDENCE_FIELDS = [
    "structure_id",
    "family_rep",
    "priority_class",
    "family_size",
    "family_sample_span",
    "length_aa",
    "mean_plddt",
    "ptm",
    "frac_ge70",
    "external_status",
    "raw_sources",
    "convincing_sources",
    "noncandidate_context_count",
    "noncandidate_context_sample_span",
    "structural_dense_fraction",
]

VALID_STATUS = {
    "CONVINCING_STRUCTURE_CONNECTION",
    "WEAK_RAW_HIT_ONLY",
    "NO_RAW_HIT_E1",
}


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1024 * 1024), b""
        ):
            h.update(block)

    return h.hexdigest()


def read_tsv(path, required):
    if not path.is_file():
        raise ValueError(f"Missing input: {path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(
            fh, delimiter="\t"
        )
        fields = reader.fieldnames or []

        missing = set(required) - set(fields)

        if missing:
            raise ValueError(
                f"{path}: missing columns {sorted(missing)}"
            )

        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty input: {path}")

    return fields, rows


def unique_index(rows, key):
    result = {}

    for row in rows:
        value = row[key]

        if not value or value in result:
            raise ValueError(
                f"Empty/duplicate {key}: {value}"
            )

        result[value] = row

    return result


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_inputs(args):
    _, darkness = read_tsv(
        args.family_audit,
        (
            "representative",
            "family_size",
            "priority_class",
        ),
    )

    _, families = read_tsv(
        args.family_context,
        (
            "family_rep",
            "family_size",
            "sample_span",
            "noncandidate_count",
            "noncandidate_anchor_context_count",
            "noncandidate_context_sample_span",
            "structural_dense_ge2_count",
        ),
    )

    member_fields, members = read_tsv(
        args.member_context,
        (
            "family_rep",
            "member",
            "sample",
            "candidate_original",
            "anchor_within_window",
            "anchor_categories",
            "window",
        ),
    )

    _, structures = read_tsv(
        args.structures,
        (
            "structure_id",
            "family_rep",
            "family_size",
            "priority_class",
            "length_aa",
            "mean_plddt",
            "ptm",
            "pdb_file",
            "pdb_sha256",
            "prediction_status",
        ),
    )

    dark_by_rep = unique_index(
        darkness, "representative"
    )
    family_by_rep = unique_index(
        families, "family_rep"
    )
    struct_by_id = unique_index(
        structures, "structure_id"
    )

    struct_reps = [
        row["family_rep"] for row in structures
    ]

    if len(set(struct_reps)) != len(struct_reps):
        raise ValueError(
            "Multiple structures for one Family"
        )

    all_reps = set(family_by_rep)

    if not (
        all_reps
        == set(dark_by_rep)
        == set(struct_reps)
    ):
        raise ValueError(
            "Stage 06/07/08 Family sets differ"
        )

    member_ids = set()
    members_by_rep = Counter()
    samples_by_rep = {}

    for row in members:
        rep = row["family_rep"]
        member = row["member"]

        if rep not in all_reps:
            raise ValueError(
                f"Unknown member Family: {rep}"
            )

        if not member or member in member_ids:
            raise ValueError(
                f"Duplicate/empty member: {member}"
            )

        member_ids.add(member)
        members_by_rep[rep] += 1

        samples_by_rep.setdefault(
            rep, set()
        ).add(row["sample"])

        # Historical Stage I2 assay evidence
        # specifically requires a five-gene window.
        if row["window"] != "5":
            raise ValueError(
                "Historical I2 assay scoring requires "
                "Stage 07 window=5"
            )

        if row["anchor_within_window"] not in (
            "0", "1"
        ):
            raise ValueError(
                f"Invalid anchor state: {member}"
            )

        if row["candidate_original"] not in (
            "0", "1"
        ):
            raise ValueError(
                f"Invalid candidate state: {member}"
            )

    prepared = []

    for sid, structure in struct_by_id.items():
        rep = structure["family_rep"]

        family = family_by_rep[rep]
        darkness_row = dark_by_rep[rep]

        size = int(family["family_size"])
        span = int(family["sample_span"])

        if size < 1 or span < 1:
            raise ValueError(
                f"Invalid Family size/span: {rep}"
            )

        if not (
            size == members_by_rep[rep]
            == int(darkness_row["family_size"])
            == int(structure["family_size"])
        ):
            raise ValueError(
                f"Family size mismatch: {rep}"
            )

        if span != len(samples_by_rep[rep]):
            raise ValueError(
                f"Family sample span mismatch: {rep}"
            )

        priority = darkness_row["priority_class"]

        if structure["priority_class"] != priority:
            raise ValueError(
                f"Priority class mismatch: {rep}"
            )

        if structure["prediction_status"] != "PREDICTED":
            raise ValueError(
                f"Missing predicted structure: {sid}"
            )

        pdb = Path(
            structure["pdb_file"]
        ).resolve()

        if (
            not pdb.is_file()
            or sha256(pdb)
            != structure["pdb_sha256"]
        ):
            raise ValueError(
                f"PDB integrity mismatch: {sid}"
            )

        dense = int(
            family["structural_dense_ge2_count"]
        )

        noncandidate = int(
            family["noncandidate_count"]
        )

        anchored_noncandidate = int(
            family["noncandidate_anchor_context_count"]
        )

        nonspan = int(
            family["noncandidate_context_sample_span"]
        )

        if not (
            0 <= dense <= size
            and 0 <= anchored_noncandidate
            <= noncandidate <= size
            and 0 <= nonspan <= span
        ):
            raise ValueError(
                f"Invalid Family context counts: {rep}"
            )

        prepared.append({
            "structure_id": sid,
            "family_rep": rep,
            "priority_class": priority,
            "family_size": str(size),
            "family_sample_span": str(span),
            "length_aa": structure["length_aa"],
            "mean_plddt": structure["mean_plddt"],
            "ptm": structure["ptm"],
            "frac_ge70": "",
            "external_status": "",
            "raw_sources": "",
            "convincing_sources": "",
            "noncandidate_context_count": str(
                anchored_noncandidate
            ),
            "noncandidate_context_sample_span": str(
                nonspan
            ),
            "structural_dense_fraction": (
                f"{dense / size:.4f}"
            ),
        })

    return prepared, member_fields, members


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--family-audit", required=True, type=Path
    )
    parser.add_argument(
        "--family-context", required=True, type=Path
    )
    parser.add_argument(
        "--member-context", required=True, type=Path
    )
    parser.add_argument(
        "--structures", required=True, type=Path
    )
    parser.add_argument(
        "--external-matrix", type=Path
    )
    parser.add_argument(
        "--external-provenance", type=Path
    )
    parser.add_argument(
        "--out-dir", type=Path
    )
    parser.add_argument(
        "--preflight", action="store_true"
    )

    args = parser.parse_args()

    if args.preflight:
        if any((
            args.external_matrix,
            args.external_provenance,
            args.out_dir,
        )):
            parser.error(
                "Preflight cannot create ranking outputs"
            )
    elif not all((
        args.external_matrix,
        args.external_provenance,
        args.out_dir,
    )):
        parser.error(
            "Full mode requires external matrix, "
            "provenance, and output directory"
        )

    prepared, member_fields, members = (
        validate_inputs(args)
    )

    print(
        f"FAMILIES={len(prepared)}",
        flush=True,
    )
    print(
        f"MEMBERS={len(members)}",
        flush=True,
    )
    print(
        "STAGE06_07_08_FAMILY_CONSISTENCY=PASS",
        flush=True,
    )
    print(
        "HISTORICAL_WINDOW5_COMPATIBILITY=PASS",
        flush=True,
    )

    if args.preflight:
        print(
            "STAGE10_INPUT_PREFLIGHT=PASS",
            flush=True,
        )
        print(
            "STAGE10_RANKING=BLOCKED_EXTERNAL_SEARCH_NOT_PROVIDED",
            flush=True,
        )
        return

    matrix_fields, external_rows = read_tsv(
        args.external_matrix,
        (
            "structure_id",
            "priority_class",
            "length",
            "mean_plddt",
            "ptm",
            "frac_ge70",
            "external_status",
            "raw_sources",
            "convincing_sources",
        ),
    )

    if len(matrix_fields) != 39:
        raise ValueError(
            "Stage 09 matrix must have 39 columns"
        )

    record = json.loads(
        args.external_provenance.read_text()
    )

    if (
        record.get("evidence_origin")
        != "COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES"
        or record.get("matrix_sha256")
        != sha256(args.external_matrix)
        or record.get("family_count")
        != len(prepared)
    ):
        raise ValueError(
            "Incomplete/unverified external evidence"
        )

    searches = record.get(
        "external_searches", {}
    )

    if set(searches) != {
        "pdb", "afsp", "bfvd"
    }:
        raise ValueError(
            "All three historical external sources required"
        )

    external = unique_index(
        external_rows, "structure_id"
    )

    if set(external) != {
        row["structure_id"] for row in prepared
    }:
        raise ValueError(
            "Stage 09/10 structure IDs differ"
        )

    for row in prepared:
        sid = row["structure_id"]
        hit = external[sid]

        if row["priority_class"] != hit["priority_class"]:
            raise ValueError(
                f"Priority mismatch: {sid}"
            )

        if int(row["length_aa"]) != int(hit["length"]):
            raise ValueError(
                f"Structure length mismatch: {sid}"
            )

        if abs(
            float(row["mean_plddt"])
            - float(hit["mean_plddt"])
        ) > 0.05:
            raise ValueError(
                f"Structure confidence mismatch: {sid}"
            )

        if hit["external_status"] not in VALID_STATUS:
            raise ValueError(
                f"Unknown external status: {sid}"
            )

        row["frac_ge70"] = hit["frac_ge70"]
        row["external_status"] = hit["external_status"]
        row["raw_sources"] = hit["raw_sources"]
        row["convincing_sources"] = hit[
            "convincing_sources"
        ]

    # Existing output folders are never overwritten.
    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    if output == Path.cwd().resolve():
        raise ValueError(
            "Unsafe output directory"
        )

    output.mkdir(parents=True)

    evidence = output / "stage10_input_evidence.tsv"

    write_tsv(
        evidence,
        EVIDENCE_FIELDS,
        prepared,
    )

    # Translate the current window column into the
    # original I2 member-context interface.
    compatible_members = [
        {
            **row,
            "anchor_within5":
                row["anchor_within_window"],
        }
        for row in members
    ]

    member_out = output / "member_context_i2.tsv"

    write_tsv(
        member_out,
        member_fields + ["anchor_within5"],
        compatible_members,
    )

    (output / "evidence_provenance.json").write_text(
        json.dumps(
            {
                "status": "EXTERNAL_EVIDENCE_VERIFIED",
                "stage09_matrix_sha256": sha256(
                    args.external_matrix
                ),
                "stage09_provenance_sha256": sha256(
                    args.external_provenance
                ),
                "evidence_sha256": sha256(evidence),
                "member_context_sha256": sha256(
                    member_out
                ),
                "family_count": len(prepared),
                "member_count": len(members),
            },
            indent=2,
        ) + "\n"
    )

    print("STAGE10_INPUT_BRIDGE=PASS")
    print(f"STAGE10_EVIDENCE={evidence}")
    print(f"STAGE10_MEMBER_CONTEXT={member_out}")
    print("STAGE10_SCORING=NOT_RUN")


if __name__ == "__main__":
    main()
