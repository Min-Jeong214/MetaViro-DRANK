#!/usr/bin/env python3
"""Bridge verified Family-MSA Foldseek evidence to Stage 13."""

import argparse
import csv
import hashlib
import json
import math
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SOURCE = ROOT / "scripts/09_external_integration_portable.py"

SOURCE_HEADER = runpy.run_path(
    str(SOURCE),
    run_name="__stage_k2_schema__",
)["HEADER"]

FIELDS = [
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

PROFILES = {
    "pdb": "pdb_e1",
    "afsp": "afsp",
    "bfvd": "bfvd",
}


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path, required=()):
    if not path.is_file():
        raise ValueError(f"MISSING_INPUT:{path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = reader.fieldnames or []
        rows = list(reader)

    if (
        len(columns) != len(set(columns))
        or not set(required).issubset(columns)
        or not rows
    ):
        raise ValueError(f"INVALID_INPUT_SCHEMA:{path}")

    return columns, rows


def indexed(rows):
    result = {}

    for row in rows:
        sid = row["structure_id"]

        if not sid or sid in result:
            raise ValueError(f"DUPLICATE_STRUCTURE_ID:{sid}")

        result[sid] = row

    return result


def numeric(value, field):
    number = float(value)

    if not math.isfinite(number):
        raise ValueError(f"NONFINITE_METRIC:{field}")

    return number


def format_evidence(row, msa_structure):
    raw = int(row["bfvd_raw_count"])
    strong = int(row["bfvd_convincing_count"])

    if raw < 0 or strong < 0 or strong > raw:
        raise ValueError("INVALID_BFVD_HIT_COUNTS")

    if raw and not row["bfvd_best_target"]:
        raise ValueError("MISSING_BFVD_BEST_TARGET")

    if not raw and row["bfvd_best_target"]:
        raise ValueError("UNEXPECTED_BFVD_BEST_TARGET")

    status = row["external_status"]

    if status not in (
        "CONVINCING_STRUCTURE_CONNECTION",
        "WEAK_RAW_HIT_ONLY",
        "NO_RAW_HIT_E1",
    ):
        raise ValueError("INVALID_EXTERNAL_STATUS")

    return {
        "structure_id": row["structure_id"],
        "external_status": status,
        "raw_sources": row["raw_sources"],
        "convincing_sources": row["convincing_sources"],
        "BFVD_best_target": row["bfvd_best_target"],
        "BFVD_best_qtm": row["bfvd_best_qtm"],
        "BFVD_best_ttm": row["bfvd_best_ttm"],
        "BFVD_best_is_convincing": (
            "1" if strong > 0 else "0"
        ),
        "msa_pdb_sha256": msa_structure["pdb_sha256"],
    }


def self_test():
    assert len(SOURCE_HEADER) == 39

    for field in (
        "bfvd_best_target",
        "bfvd_best_qtm",
        "bfvd_best_ttm",
        "bfvd_convincing_count",
    ):
        assert field in SOURCE_HEADER

    mock = dict.fromkeys(SOURCE_HEADER, "")
    mock.update({
        "structure_id": "TEST001",
        "external_status": "WEAK_RAW_HIT_ONLY",
        "raw_sources": "BFVD",
        "convincing_sources": "",
        "bfvd_raw_count": "1",
        "bfvd_convincing_count": "0",
        "bfvd_best_target": "SIMULATED_TARGET",
        "bfvd_best_qtm": "0.40",
        "bfvd_best_ttm": "0.35",
    })

    source = {"pdb_sha256": "SYNTHETIC_HASH"}

    weak = format_evidence(mock, source)

    assert weak["BFVD_best_target"] == "SIMULATED_TARGET"
    assert weak["BFVD_best_is_convincing"] == "0"

    mock["external_status"] = (
        "CONVINCING_STRUCTURE_CONNECTION"
    )
    mock["convincing_sources"] = "BFVD"
    mock["bfvd_convincing_count"] = "1"

    strong = format_evidence(mock, source)

    assert strong["BFVD_best_is_convincing"] == "1"
    assert strong["convincing_sources"] == "BFVD"

    print("HISTORICAL_39_COLUMN_INPUT=PASS")
    print("STAGE13_BFVD_FIELD_MAPPING=PASS")
    print("WEAK_BFVD_HIT_LABELING=PASS")
    print("CONVINCING_BFVD_HIT_LABELING=PASS")
    print("STAGE_K2_BRIDGE_UNIT_TEST=PASS")
    print("SYNTHETIC_EXTERNAL_EVIDENCE_NOT_WRITTEN=YES")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--k2-input", type=Path)
    parser.add_argument("--k2-provenance", type=Path)
    parser.add_argument("--external-matrix", type=Path)
    parser.add_argument("--external-provenance", type=Path)
    parser.add_argument("--out-dir", type=Path)

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    required = (
        args.k2_input,
        args.k2_provenance,
        args.external_matrix,
        args.external_provenance,
        args.out_dir,
    )

    if any(path is None for path in required):
        parser.error("All five paths are required")

    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError("Refusing existing output directory")

    k2_meta = json.loads(args.k2_provenance.read_text())
    ext_meta = json.loads(args.external_provenance.read_text())

    if (
        k2_meta.get("query_origin") != "STAGE11_FAMILY_MSA"
        or k2_meta.get("purpose")
        != "STAGE_K2_EXTERNAL_STRUCTURE_SEARCH"
        or k2_meta.get("structure_table_sha256")
        != sha(args.k2_input)
    ):
        raise ValueError("UNVERIFIED_MSA_STRUCTURE_INPUT")

    if (
        ext_meta.get("evidence_origin")
        != "COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES"
        or ext_meta.get("structure_table_sha256")
        != sha(args.k2_input)
        or ext_meta.get("matrix_sha256")
        != sha(args.external_matrix)
    ):
        raise ValueError(
            "STAGE_K2_EXTERNAL_EVIDENCE_UNVERIFIED"
        )

    searches = ext_meta.get("external_searches", {})

    if set(searches) != set(PROFILES):
        raise ValueError("INCOMPLETE_EXTERNAL_DATABASE_SEARCH")

    for prefix, expected in PROFILES.items():
        search = searches[prefix]

        if search.get("source_profile") != expected:
            raise ValueError("INVALID_SEARCH_PROFILE")

        directory = Path(search["directory"]).resolve()
        info_path = directory / "search_provenance.json"
        hits_path = directory / "raw_hits.tsv"
        manifest_path = directory / "query_manifest.tsv"

        for path, field in (
            (info_path, "search_provenance_sha256"),
            (hits_path, "raw_hits_sha256"),
            (manifest_path, "query_manifest_sha256"),
        ):
            if not path.is_file() or sha(path) != search.get(field):
                raise ValueError(
                    f"SEARCH_PROVENANCE_MISMATCH:{prefix}"
                )

        info = json.loads(info_path.read_text())

        if (
            info.get("control_only") is not False
            or info.get("source_profile") != expected
        ):
            raise ValueError(
                f"CONTROL_OR_INVALID_SEARCH:{prefix}"
            )

    columns, external_rows = read_tsv(
        args.external_matrix,
        SOURCE_HEADER,
    )

    if columns != SOURCE_HEADER:
        raise ValueError("STAGE09_SCHEMA_MISMATCH")

    _, k2_rows = read_tsv(
        args.k2_input,
        (
            "structure_id",
            "pdb_file",
            "pdb_sha256",
            "priority_class",
            "length_aa",
            "mean_plddt",
            "ptm",
        ),
    )

    external = indexed(external_rows)
    structures = indexed(k2_rows)

    if (
        set(external) != set(structures)
        or ext_meta.get("family_count") != len(structures)
        or k2_meta.get("structure_count") != len(structures)
    ):
        raise ValueError("STAGE_K2_FAMILY_SET_MISMATCH")

    results = []

    for sid, structure in structures.items():
        row = external[sid]

        pdb = Path(structure["pdb_file"]).resolve()

        if (
            not pdb.is_file()
            or sha(pdb) != structure["pdb_sha256"]
        ):
            raise ValueError(f"MSA_PDB_INTEGRITY_FAILURE:{sid}")

        if (
            row["priority_class"]
            != structure["priority_class"]
            or int(row["length"])
            != int(structure["length_aa"])
            or abs(
                numeric(row["mean_plddt"], "mean_plddt")
                - numeric(
                    structure["mean_plddt"],
                    "structure_mean_plddt",
                )
            ) > 0.05
            or abs(
                numeric(row["ptm"], "ptm")
                - numeric(structure["ptm"], "structure_ptm")
            ) > 0.001
        ):
            raise ValueError(
                f"MSA_STRUCTURE_EVIDENCE_MISMATCH:{sid}"
            )

        results.append(format_evidence(row, structure))

    output.mkdir(parents=True, exist_ok=False)

    table = output / "msa_external_evidence_stage13.tsv"

    with table.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)

    provenance = {
        "source": "VERIFIED_STAGE_K2_EXTERNAL_SEARCH",
        "query_origin": "STAGE11_FAMILY_MSA",
        "k2_input_sha256": sha(args.k2_input),
        "k2_provenance_sha256": sha(args.k2_provenance),
        "external_matrix_sha256": sha(args.external_matrix),
        "external_provenance_sha256": sha(
            args.external_provenance
        ),
        "families": len(results),
        "stage13_evidence_sha256": sha(table),
        "control_only": False,
        "final_selection_performed": False,
    }

    (
        output / "msa_external_evidence_provenance.json"
    ).write_text(
        json.dumps(provenance, indent=2) + "\n"
    )

    print(f"EXTERNAL_EVIDENCE_FAMILIES={len(results)}")
    print(f"STAGE13_EXTERNAL_COLUMNS={len(FIELDS)}")
    print(f"STAGE13_EXTERNAL_TABLE={table}")
    print("STAGE_K2_EXTERNAL_EVIDENCE_BRIDGE=PASS")
    print("FINAL_CANDIDATE_SELECTION=NOT_RUN")


if __name__ == "__main__":
    main()
