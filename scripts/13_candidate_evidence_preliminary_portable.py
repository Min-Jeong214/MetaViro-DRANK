#!/usr/bin/env python3
"""Assemble preliminary candidate evidence without final selection."""

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


FIELDS = [
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


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(
            lambda: fh.read(1048576), b""
        ):
            h.update(chunk)
    return h.hexdigest()


def checked(path, expected):
    path = Path(path).resolve()

    if not path.is_file() or sha(path) != expected:
        raise ValueError(
            f"INPUT_INTEGRITY_MISMATCH:{path}"
        )

    return path


def read_table(path, required=(), allow_empty=False):
    if not path.is_file():
        raise ValueError(f"MISSING_INPUT:{path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        if (
            len(fields) != len(set(fields))
            or not set(required).issubset(fields)
        ):
            raise ValueError(
                f"INVALID_TABLE_SCHEMA:{path}"
            )

        rows = list(reader)

    if not rows and not allow_empty:
        raise ValueError(f"EMPTY_INPUT:{path}")

    return rows


def index(rows, key="structure_id"):
    result = {}

    for row in rows:
        sid = row[key]

        if not sid or sid in result:
            raise ValueError(
                f"DUPLICATE_OR_EMPTY_ID:{sid}"
            )

        result[sid] = row

    return result


def load_sources(args):
    selection = read_table(
        args.selection, ("structure_id",)
    )

    ids = [row["structure_id"] for row in selection]

    if (
        len(ids) != len(set(ids))
        or any(
            not re.fullmatch(
                r"[A-Za-z0-9][A-Za-z0-9_.-]*",
                sid,
            )
            for sid in ids
        )
    ):
        raise ValueError("INVALID_SELECTION_IDS")

    qc_rows = read_table(
        args.msa_qc,
        (
            "structure_id",
            "single_mean_plddt",
            "msa_mean_plddt",
            "delta_mean_plddt",
            "single_ptm",
            "msa_ptm",
            "delta_ptm",
            "msa_pdb_sha256",
        ),
    )

    qc_prov = json.loads(
        args.msa_qc_provenance.read_text()
    )

    if (
        qc_prov.get("algorithm")
        != "VERIFIED_SINGLE_VS_FAMILY_MSA"
        or qc_prov.get("msa_qc_sha256")
        != sha(args.msa_qc)
    ):
        raise ValueError("INVALID_MSA_QC_PROVENANCE")

    cons_dir = args.conservation_dir.resolve()
    cons_path = (
        cons_dir / "conservation_provenance.json"
    )

    cons_prov = json.loads(cons_path.read_text())

    if (
        cons_prov.get("algorithm")
        != "ORIGINAL_STAGE12_FUNCTIONS"
        or cons_prov.get("selection_sha256")
        != sha(args.selection)
    ):
        raise ValueError(
            "INVALID_CONSERVATION_PROVENANCE"
        )

    cons_summary = checked(
        cons_dir / "candidate_conservation_summary.tsv",
        cons_prov["summary_sha256"],
    )

    checked(
        cons_dir / "conserved_windows.tsv",
        cons_prov["windows_sha256"],
    )

    cons_rows = read_table(
        cons_summary,
        (
            "structure_id",
            "very_conserved_positions",
            "priority_structural_positions",
            "priority_acidic",
            "priority_basic",
        ),
    )

    spatial_dir = args.spatial_dir.resolve()
    spatial_path = (
        spatial_dir / "spatial_cluster_provenance.json"
    )

    spatial_prov = json.loads(
        spatial_path.read_text()
    )

    if (
        spatial_prov.get("analysis_mode")
        != "NEW_EXPLORATORY_L3"
        or spatial_prov.get(
            "historical_stage_l3_reproduction"
        ) is not False
        or spatial_prov.get(
            "biochemical_activity_inferred"
        ) is not False
        or spatial_prov.get(
            "stage12_provenance_sha256"
        ) != sha(cons_path)
    ):
        raise ValueError(
            "INVALID_EXPLORATORY_L3_PROVENANCE"
        )

    spatial_summary = checked(
        spatial_dir / "family_spatial_cluster_summary.tsv",
        spatial_prov["summary_table_sha256"],
    )

    spatial_clusters = checked(
        spatial_dir / "charged_spatial_clusters.tsv",
        spatial_prov["cluster_table_sha256"],
    )

    summary_rows = read_table(
        spatial_summary,
        ("structure_id", "cluster_count", "status"),
    )

    cluster_rows = read_table(
        spatial_clusters,
        ("structure_id", "cluster_rank", "residues"),
        allow_empty=True,
    )

    qc = index(qc_rows)
    conservation = index(cons_rows)
    spatial = index(summary_rows)

    selected = set(ids)

    if not (
        selected.issubset(qc)
        and selected == set(conservation)
        and selected == set(spatial)
    ):
        raise ValueError("CROSS_STAGE_FAMILY_SET_MISMATCH")

    cons_families = index(
        cons_prov.get("families", [])
    )

    qc_evidence = index(
        qc_prov.get("evidence", [])
    )

    spatial_inputs = index(
        spatial_prov.get("inputs", [])
    )

    if not (
        selected == set(cons_families)
        == set(spatial_inputs)
    ):
        raise ValueError("FAMILY_PROVENANCE_SET_MISMATCH")

    if not selected.issubset(qc_evidence):
        raise ValueError("MISSING_QC_SOURCE_EVIDENCE")

    if (
        cons_prov.get("family_count") != len(ids)
        or spatial_prov.get("family_count") != len(ids)
        or spatial_prov.get("cluster_count")
        != len(cluster_rows)
    ):
        raise ValueError("CROSS_STAGE_COUNT_MISMATCH")

    clusters = defaultdict(list)

    for row in cluster_rows:
        sid = row["structure_id"]

        if sid not in selected:
            raise ValueError(
                f"UNKNOWN_SPATIAL_CLUSTER:{sid}"
            )

        clusters[sid].append(row)

    for sid in ids:
        family = cons_families[sid]
        qc_source = qc_evidence[sid]
        spatial_source = spatial_inputs[sid]

        residue_table = checked(
            cons_dir / f"{sid}_residue_conservation.tsv",
            family["residue_table_sha256"],
        )

        if (
            family["prediction_provenance_sha256"]
            != qc_source["stage11_record_sha256"]
            or family["a3m_sha256"]
            != qc_source["a3m_sha256"]
            or family["pdb_sha256"]
            != qc[sid]["msa_pdb_sha256"]
            or spatial_source["residue_table_sha256"]
            != sha(residue_table)
            or spatial_source["msa_pdb_sha256"]
            != family["pdb_sha256"]
        ):
            raise ValueError(
                f"CROSS_STAGE_SOURCE_MISMATCH:{sid}"
            )

        grouped = sorted(
            clusters[sid],
            key=lambda row: int(row["cluster_rank"]),
        )

        if (
            len(grouped)
            != int(spatial[sid]["cluster_count"])
            or [
                int(row["cluster_rank"])
                for row in grouped
            ] != list(range(1, len(grouped) + 1))
        ):
            raise ValueError(
                f"INVALID_SPATIAL_CLUSTER_RANKS:{sid}"
            )

        clusters[sid] = grouped

    return (
        ids,
        qc,
        conservation,
        spatial,
        clusters,
        {
            "msa_qc_provenance": qc_prov,
            "conservation_provenance": cons_prov,
            "spatial_provenance": spatial_prov,
        },
    )


def assemble(ids, qc, conservation, spatial, clusters):
    results = []

    for sid in ids:
        q = qc[sid]
        c = conservation[sid]
        s = spatial[sid]

        row = dict.fromkeys(FIELDS, "")

        row["structure_id"] = sid

        for field in (
            "single_mean_plddt",
            "msa_mean_plddt",
            "delta_mean_plddt",
            "single_ptm",
            "msa_ptm",
            "delta_ptm",
        ):
            row[field] = q[field]

        for field in (
            "very_conserved_positions",
            "priority_structural_positions",
            "priority_acidic",
            "priority_basic",
        ):
            row[field] = c[field]

        row["external_status"] = "NOT_SEARCHED"
        row["spatial_cluster_count"] = s["cluster_count"]
        row["spatial_clusters"] = ";".join(
            r["residues"] for r in clusters[sid]
        )

        results.append(row)

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--msa-qc", type=Path, required=True)
    ap.add_argument(
        "--msa-qc-provenance",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--conservation-dir",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--spatial-dir",
        type=Path,
        required=True,
    )
    ap.add_argument("--out-dir", type=Path, required=True)

    args = ap.parse_args()
    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    ids, qc, cons, spatial, clusters, _ = load_sources(
        args
    )

    rows = assemble(
        ids, qc, cons, spatial, clusters
    )

    output.mkdir(parents=True, exist_ok=False)

    matrix = (
        output / "preliminary_candidate_evidence.tsv"
    )

    with matrix.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "status": "PRELIMINARY_NOT_FINALIZED",
        "selection_method": "EXPLICIT_INPUT_ONLY",
        "historical_stage13_reproduction": False,
        "stage10_ranking_available": False,
        "stage_k2_external_search_available": False,
        "spatial_analysis": "NEW_EXPLORATORY_L3",
        "final_top5_selected": False,
        "primary_top3_selected": False,
        "family_count": len(rows),
        "matrix_sha256": sha(matrix),
        "selection_sha256": sha(args.selection),
        "msa_qc_sha256": sha(args.msa_qc),
        "conservation_provenance_sha256": sha(
            args.conservation_dir /
            "conservation_provenance.json"
        ),
        "spatial_provenance_sha256": sha(
            args.spatial_dir /
            "spatial_cluster_provenance.json"
        ),
    }

    (
        output / "preliminary_evidence_provenance.json"
    ).write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    print(f"ASSEMBLED_FAMILIES={len(rows)}")
    print(f"MATRIX_COLUMNS={len(FIELDS)}")
    print("STAGE10_RANKING=NOT_PROVIDED")
    print("STAGE_K2_EXTERNAL_EVIDENCE=NOT_SEARCHED")
    print("SPATIAL_ANALYSIS=EXPLORATORY_ONLY")
    print("FINAL_TOP5=NOT_SELECTED")
    print("PRIMARY_TOP3=NOT_SELECTED")
    print("STAGE13_PRELIMINARY_ASSEMBLY=PASS")


if __name__ == "__main__":
    main()
