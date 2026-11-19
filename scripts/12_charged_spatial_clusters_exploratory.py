#!/usr/bin/env python3
"""Exploratory, non-historical C-alpha spatial charge clusters."""

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STAGE12 = ROOT / "scripts/12_msa_structure_conservation.py"

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N",
    "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

CHARGED = set("DEKRH")
ACIDIC = set("DE")
BASIC = set("KRH")

CLUSTER_FIELDS = [
    "structure_id",
    "cluster_rank",
    "residues",
    "residue_count",
    "acidic_count",
    "basic_count",
    "mean_plddt",
    "max_edge_ca_distance_angstrom",
]

SUMMARY_FIELDS = [
    "structure_id",
    "length",
    "eligible_charged_residues",
    "cluster_count",
    "status",
]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def checked(path, digest):
    path = Path(path).resolve()
    if not path.is_file() or sha(path) != digest:
        raise ValueError(f"INTEGRITY_FAILURE:{path}")
    return path


def read_tsv(path, required):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = reader.fieldnames or []

        if (
            len(columns) != len(set(columns))
            or not set(required).issubset(columns)
        ):
            raise ValueError(f"INVALID_TABLE_SCHEMA:{path}")

        rows = list(reader)

    if not rows:
        raise ValueError(f"EMPTY_RESIDUE_TABLE:{path}")

    return rows


def read_ca(path):
    atoms = {}

    with path.open(errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue

            if line[26] != " ":
                raise ValueError("PDB_INSERTION_CODES_UNSUPPORTED")

            position = int(line[22:26])

            if position in atoms:
                raise ValueError("DUPLICATE_PDB_CA")

            residue = AA3.get(line[17:20].strip())

            if residue is None:
                raise ValueError("NONSTANDARD_PDB_RESIDUE")

            xyz = tuple(
                float(line[start:start + 8])
                for start in (30, 38, 46)
            )

            confidence = float(line[60:66])

            if (
                not all(math.isfinite(x) for x in xyz)
                or not math.isfinite(confidence)
                or not 0 <= confidence <= 100
            ):
                raise ValueError("INVALID_PDB_COORDINATES")

            atoms[position] = {
                "position": position,
                "aa": residue,
                "xyz": xyz,
                "plddt": confidence,
            }

    if not atoms:
        raise ValueError("NO_PDB_CA_ATOMS")

    if sorted(atoms) != list(range(1, len(atoms) + 1)):
        raise ValueError("NONCONTIGUOUS_PDB_RESIDUES")

    return atoms


def distance(a, b):
    return math.sqrt(
        sum(
            (x - y) ** 2
            for x, y in zip(a["xyz"], b["xyz"])
        )
    )


def find_clusters(candidates, cutoff, minimum):
    """Connected components, not all-pairs spatial cliques."""

    graph = {r["position"]: set() for r in candidates}
    by_position = {
        r["position"]: r for r in candidates
    }
    edges = {}

    for i, first in enumerate(candidates):
        for second in candidates[i + 1:]:
            d = distance(first, second)

            if d <= cutoff:
                a = first["position"]
                b = second["position"]
                graph[a].add(b)
                graph[b].add(a)
                edges[(min(a, b), max(a, b))] = d

    visited = set()
    clusters = []

    for start in sorted(graph):
        if start in visited:
            continue

        stack = [start]
        component = set()

        while stack:
            position = stack.pop()

            if position in visited:
                continue

            visited.add(position)
            component.add(position)

            stack.extend(
                sorted(graph[position] - visited)
            )

        if len(component) < minimum:
            continue

        members = [
            by_position[p] for p in sorted(component)
        ]

        member_edges = [
            d for (a, b), d in edges.items()
            if a in component and b in component
        ]

        clusters.append({
            "positions": sorted(component),
            "members": members,
            "mean_plddt": sum(
                r["plddt"] for r in members
            ) / len(members),
            "max_edge_distance": max(member_edges),
        })

    clusters.sort(
        key=lambda c: (
            -len(c["members"]),
            -c["mean_plddt"],
            c["positions"][0],
        )
    )

    return clusters


def self_test():
    def point(position, aa, x, confidence=90):
        return {
            "position": position,
            "aa": aa,
            "xyz": (float(x), 0.0, 0.0),
            "plddt": confidence,
        }

    residues = [
        point(1, "D", 0),
        point(2, "K", 4),
        point(3, "E", 8),
        point(4, "R", 40),
    ]

    clusters = find_clusters(residues, 8.0, 3)

    assert len(clusters) == 1
    assert clusters[0]["positions"] == [1, 2, 3]
    assert clusters[0]["max_edge_distance"] == 8.0

    assert find_clusters(residues, 3.0, 3) == []
    assert find_clusters(residues, 8.0, 4) == []

    print("SYNTHETIC_3D_CLUSTER_DETECTION=PASS")
    print("DISTANCE_CUTOFF_BEHAVIOR=PASS")
    print("MINIMUM_RESIDUE_GUARD=PASS")
    print("EXPLORATORY_NOT_HISTORICAL=YES")
    print("STAGE_L3_GEOMETRY_UNIT_TEST=PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--conservation-dir", type=Path)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument(
        "--ca-cutoff",
        type=float,
        default=8.0,
    )
    ap.add_argument(
        "--min-residues",
        type=int,
        default=3,
    )
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    if args.conservation_dir is None or args.out_dir is None:
        ap.error("Conservation and output directories required")

    cutoff = args.ca_cutoff
    minimum = args.min_residues

    if (
        not math.isfinite(cutoff)
        or cutoff <= 0
        or minimum < 2
    ):
        raise ValueError("INVALID_CLUSTER_PARAMETERS")

    source = args.conservation_dir.resolve()
    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            f"Refusing existing output directory: {output}"
        )

    if (
        output == source
        or source in output.parents
        or output in source.parents
    ):
        raise ValueError("OUTPUT_OVERLAPS_INPUT")

    meta_path = source / "conservation_provenance.json"

    if not meta_path.is_file():
        raise ValueError("MISSING_STAGE12_PROVENANCE")

    meta = json.loads(meta_path.read_text())

    if (
        meta.get("algorithm")
        != "ORIGINAL_STAGE12_FUNCTIONS"
        or meta.get("original_script_sha256")
        != sha(STAGE12)
    ):
        raise ValueError("STAGE12_ALGORITHM_MISMATCH")

    checked(
        source / "candidate_conservation_summary.tsv",
        meta["summary_sha256"],
    )

    checked(
        source / "conserved_windows.tsv",
        meta["windows_sha256"],
    )

    families = meta.get("families", [])

    if (
        not families
        or len(families) != meta.get("family_count")
    ):
        raise ValueError("STAGE12_FAMILY_COUNT_MISMATCH")

    seen = set()
    cluster_rows = []
    summaries = []
    input_checks = []

    for family in families:
        sid = family["structure_id"]

        if (
            not re.fullmatch(
                r"[A-Za-z0-9][A-Za-z0-9_.-]*",
                sid,
            )
            or sid in seen
        ):
            raise ValueError("INVALID_OR_DUPLICATE_FAMILY_ID")

        seen.add(sid)

        residue_path = checked(
            source / f"{sid}_residue_conservation.tsv",
            family["residue_table_sha256"],
        )

        rows = read_tsv(
            residue_path,
            (
                "structure_id",
                "query_position",
                "query_aa",
                "msa_plddt",
                "high_conserved",
                "structure_supported",
                "priority_residue",
            ),
        )

        prediction_path = checked(
            family["prediction_provenance"],
            family["prediction_provenance_sha256"],
        )

        prediction = json.loads(
            prediction_path.read_text()
        )

        if (
            prediction.get("status") != "COMPLETED"
            or prediction.get("structure_id") != sid
        ):
            raise ValueError(f"INVALID_PREDICTION:{sid}")

        pdb = checked(
            prediction["pdb_file"],
            family["pdb_sha256"],
        )

        if prediction["pdb_sha256"] != sha(pdb):
            raise ValueError(f"PDB_RECORD_MISMATCH:{sid}")

        atoms = read_ca(pdb)

        if (
            len(rows) != len(atoms)
            or len(rows) != family["residue_count"]
        ):
            raise ValueError(f"RESIDUE_LENGTH_MISMATCH:{sid}")

        eligible = []
        seen_positions = set()

        for row in rows:
            pos = int(row["query_position"])

            if (
                row["structure_id"] != sid
                or pos not in atoms
                or pos in seen_positions
            ):
                raise ValueError(f"RESIDUE_MAPPING_ERROR:{sid}")

            seen_positions.add(pos)
            atom = atoms[pos]

            if row["query_aa"] != atom["aa"]:
                raise ValueError(
                    f"PDB_SEQUENCE_MISMATCH:{sid}:{pos}"
                )

            if abs(
                float(row["msa_plddt"])
                - atom["plddt"]
            ) > 0.05:
                raise ValueError(
                    f"PLDDT_MISMATCH:{sid}:{pos}"
                )

            flags = [
                row["high_conserved"],
                row["structure_supported"],
                row["priority_residue"],
            ]

            if any(flag not in ("0", "1") for flag in flags):
                raise ValueError("INVALID_CONSERVATION_FLAGS")

            if row["structure_supported"] != str(
                int(atom["plddt"] >= 70)
            ):
                raise ValueError("STRUCTURE_SUPPORT_MISMATCH")

            if row["priority_residue"] == "1":
                if (
                    row["high_conserved"] != "1"
                    or row["structure_supported"] != "1"
                ):
                    raise ValueError(
                        "INCONSISTENT_PRIORITY_RESIDUE"
                    )

                if atom["aa"] in CHARGED:
                    eligible.append(atom)

        if len(seen_positions) != len(atoms):
            raise ValueError("INCOMPLETE_RESIDUE_MAPPING")

        clusters = find_clusters(
            eligible, cutoff, minimum
        )

        for rank, cluster in enumerate(clusters, 1):
            members = cluster["members"]

            cluster_rows.append({
                "structure_id": sid,
                "cluster_rank": str(rank),
                "residues": ",".join(
                    f'{r["aa"]}{r["position"]}'
                    for r in members
                ),
                "residue_count": str(len(members)),
                "acidic_count": str(sum(
                    r["aa"] in ACIDIC for r in members
                )),
                "basic_count": str(sum(
                    r["aa"] in BASIC for r in members
                )),
                "mean_plddt": (
                    f'{cluster["mean_plddt"]:.2f}'
                ),
                "max_edge_ca_distance_angstrom": (
                    f'{cluster["max_edge_distance"]:.3f}'
                ),
            })

        summaries.append({
            "structure_id": sid,
            "length": str(len(atoms)),
            "eligible_charged_residues": str(len(eligible)),
            "cluster_count": str(len(clusters)),
            "status": (
                "CLUSTERS_DETECTED"
                if clusters
                else "NO_CLUSTER_DETECTED"
            ),
        })

        input_checks.append({
            "structure_id": sid,
            "residue_table_sha256": sha(residue_path),
            "msa_pdb_sha256": sha(pdb),
            "prediction_provenance_sha256": sha(
                prediction_path
            ),
        })

    # All input validation and clustering precede output creation.
    output.mkdir(parents=True, exist_ok=False)

    def write(path, fields, rows):
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    cluster_file = (
        output / "charged_spatial_clusters.tsv"
    )
    summary_file = (
        output / "family_spatial_cluster_summary.tsv"
    )

    write(cluster_file, CLUSTER_FIELDS, cluster_rows)
    write(summary_file, SUMMARY_FIELDS, summaries)

    provenance = {
        "analysis_mode": "NEW_EXPLORATORY_L3",
        "historical_stage_l3_reproduction": False,
        "geometry": "CA_DISTANCE_CONNECTED_COMPONENTS",
        "ca_distance_cutoff_angstrom": cutoff,
        "min_cluster_residues": minimum,
        "residue_filter": (
            "STAGE12_PRIORITY_RESIDUE_AND_DEKRH"
        ),
        "biochemical_activity_inferred": False,
        "stage12_provenance_sha256": sha(meta_path),
        "family_count": len(families),
        "cluster_count": len(cluster_rows),
        "cluster_table_sha256": sha(cluster_file),
        "summary_table_sha256": sha(summary_file),
        "inputs": input_checks,
    }

    (
        output / "spatial_cluster_provenance.json"
    ).write_text(
        json.dumps(provenance, indent=2) + "\n"
    )

    print(f"ANALYZED_FAMILIES={len(families)}")
    print(f"DETECTED_CLUSTERS={len(cluster_rows)}")
    print("ANALYSIS_MODE=NEW_EXPLORATORY_L3")
    print("HISTORICAL_L3_REPRODUCTION=NO")
    print("STAGE_L3_EXPLORATORY_ANALYSIS=PASS")


if __name__ == "__main__":
    main()
