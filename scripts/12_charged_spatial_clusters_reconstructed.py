#!/usr/bin/env python3
"""Reconstruct Stage L3 from historical PDB and TSV evidence.

This is an independently reconstructed algorithm, not recovered
historical source code.
"""

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path


TERMINAL = {
    "D": ("OD1", "OD2"),
    "E": ("OE1", "OE2"),
    "K": ("NZ",),
    "R": ("NH1", "NH2"),
    "H": ("ND1", "NE2"),
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N",
    "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

PAIR_COLUMNS = [
    "structure_id",
    "residue1",
    "residue2",
    "distance_A",
    "sequence_separation",
    "r1_conservation",
    "r2_conservation",
    "r1_plddt",
    "r2_plddt",
]

CLUSTER_COLUMNS = [
    "structure_id",
    "cluster_rank",
    "n_residues",
    "residues",
    "acidic",
    "basic",
    "sequence_span",
    "mean_conservation",
    "mean_plddt",
]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for part in iter(lambda: fh.read(1048576), b""):
            h.update(part)
    return h.hexdigest()


def read_residues(path, sid):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {
            "structure_id",
            "query_position",
            "query_aa",
            "query_residue_fraction",
            "msa_plddt",
            "high_conserved",
            "structure_supported",
            "priority_residue",
        }

        if not required.issubset(reader.fieldnames or []):
            raise ValueError("MISSING_STAGE12_COLUMNS")

        rows = list(reader)

    if not rows:
        raise ValueError("EMPTY_RESIDUE_TABLE")

    positions = [int(r["query_position"]) for r in rows]

    if positions != list(range(1, len(rows) + 1)):
        raise ValueError("NONCONTIGUOUS_RESIDUE_TABLE")

    for row in rows:
        if row["structure_id"] != sid:
            raise ValueError("STRUCTURE_ID_MISMATCH")

        for flag in (
            "high_conserved",
            "structure_supported",
            "priority_residue",
        ):
            if row[flag] not in ("0", "1"):
                raise ValueError(f"INVALID_FLAG:{flag}")

        conservation = float(row["query_residue_fraction"])
        confidence = float(row["msa_plddt"])

        if (
            not math.isfinite(conservation)
            or not 0 <= conservation <= 1
            or not math.isfinite(confidence)
            or not 0 <= confidence <= 100
        ):
            raise ValueError("INVALID_RESIDUE_METRIC")

        if row["priority_residue"] == "1":
            if (
                row["high_conserved"] != "1"
                or row["structure_supported"] != "1"
            ):
                raise ValueError("INCONSISTENT_PRIORITY_FLAG")

    return rows


def read_pdb(path):
    residues = defaultdict(dict)
    names = {}

    with path.open(errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue

            if line[21] != "A":
                continue

            if line[16] not in (" ", "A"):
                continue

            if line[26] != " ":
                raise ValueError("PDB_INSERTION_CODE_UNSUPPORTED")

            position = int(line[22:26])
            atom = line[12:16].strip()
            residue = AA3.get(line[17:20].strip())

            if residue is None:
                raise ValueError("NONSTANDARD_PDB_RESIDUE")

            if position in names and names[position] != residue:
                raise ValueError("INCONSISTENT_PDB_RESIDUE")

            names[position] = residue

            point = tuple(
                float(line[start:start + 8])
                for start in (30, 38, 46)
            )

            if not all(math.isfinite(x) for x in point):
                raise ValueError("INVALID_PDB_COORDINATE")

            residues[position][atom] = {
                "xyz": point,
                "bfactor": float(line[60:66]),
            }

    if not residues:
        raise ValueError("NO_PDB_CHAIN_A")

    if sorted(residues) != list(
        range(1, len(residues) + 1)
    ):
        raise ValueError("NONCONTIGUOUS_PDB")

    return dict(residues), names


def charge_center(atoms, position, aa):
    names = TERMINAL[aa]
    residue = atoms[position]

    if any(name not in residue for name in names):
        raise ValueError(
            f"MISSING_TERMINAL_CHARGE_ATOMS:{aa}{position}"
        )

    return tuple(
        sum(residue[name]["xyz"][axis] for name in names)
        / len(names)
        for axis in range(3)
    )


def connected_components(nodes, adjacency, minimum):
    visited = set()
    components = []

    for start in sorted(nodes):
        if start in visited:
            continue

        pending = [start]
        component = set()

        while pending:
            current = pending.pop()

            if current in visited:
                continue

            visited.add(current)
            component.add(current)

            pending.extend(
                sorted(adjacency[current] - visited)
            )

        if len(component) >= minimum:
            components.append(sorted(component))

    components.sort(
        key=lambda positions: (
            -len(positions),
            positions[0],
        )
    )

    return components


def calculate(sid, residue_rows, pdb_path, cutoff, minimum):
    atoms, names = read_pdb(pdb_path)

    if len(atoms) != len(residue_rows):
        raise ValueError("PDB_RESIDUE_COUNT_MISMATCH")

    eligible = []

    for row in residue_rows:
        position = int(row["query_position"])
        aa = row["query_aa"]

        if names.get(position) != aa:
            raise ValueError(
                f"PDB_SEQUENCE_MISMATCH:{position}"
            )

        if "CA" not in atoms[position]:
            raise ValueError(
                f"MISSING_CA_ATOM:{position}"
            )

        pdb_plddt = atoms[position]["CA"]["bfactor"]
        table_plddt = float(row["msa_plddt"])

        if abs(pdb_plddt - table_plddt) > 0.05:
            raise ValueError(
                f"PLDDT_MISMATCH:{position}"
            )

        if row["structure_supported"] != str(
            int(pdb_plddt >= 70)
        ):
            raise ValueError(
                f"STRUCTURE_SUPPORT_MISMATCH:{position}"
            )

        if (
            row["priority_residue"] == "1"
            and aa in TERMINAL
        ):
            eligible.append({
                "position": position,
                "aa": aa,
                "center": charge_center(
                    atoms, position, aa
                ),
                "conservation": float(
                    row["query_residue_fraction"]
                ),
                "plddt": pdb_plddt,
            })

    eligible.sort(key=lambda r: r["position"])

    graph = defaultdict(set)
    pairs = []

    for i, first in enumerate(eligible):
        for second in eligible[i + 1:]:
            distance = math.dist(
                first["center"],
                second["center"],
            )

            if distance > cutoff:
                continue

            a = first["position"]
            b = second["position"]

            graph[a].add(b)
            graph[b].add(a)

            pairs.append({
                "structure_id": sid,
                "residue1": f'{first["aa"]}{a}',
                "residue2": f'{second["aa"]}{b}',
                "distance_A": f"{distance:.2f}",
                "sequence_separation": str(b - a),
                "r1_conservation": (
                    f'{first["conservation"]:.3f}'
                ),
                "r2_conservation": (
                    f'{second["conservation"]:.3f}'
                ),
                "r1_plddt": f'{first["plddt"]:.2f}',
                "r2_plddt": f'{second["plddt"]:.2f}',
            })

    components = connected_components(
        graph.keys(),
        graph,
        minimum,
    )

    by_position = {
        r["position"]: r for r in eligible
    }

    clusters = []

    for rank, positions in enumerate(components, 1):
        members = [by_position[p] for p in positions]

        mean_cons = sum(
            m["conservation"] for m in members
        ) / len(members)

        mean_plddt = sum(
            m["plddt"] for m in members
        ) / len(members)

        clusters.append({
            "structure_id": sid,
            "cluster_rank": str(rank),
            "n_residues": str(len(members)),
            "residues": ",".join(
                f'{m["aa"]}{m["position"]}'
                for m in members
            ),
            "acidic": str(sum(
                m["aa"] in "DE" for m in members
            )),
            "basic": str(sum(
                m["aa"] in "KRH" for m in members
            )),
            "sequence_span": str(
                positions[-1] - positions[0]
            ),
            "mean_conservation": f"{mean_cons:.3f}",
            "mean_plddt": f"{mean_plddt:.2f}",
        })

    return eligible, pairs, clusters


def write_tsv(path, columns, rows):
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure-id", required=True)
    ap.add_argument(
        "--residue-table", type=Path, required=True
    )
    ap.add_argument("--pdb", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument(
        "--cutoff-angstrom",
        type=float,
        default=8.0,
    )
    ap.add_argument(
        "--min-residues",
        type=int,
        default=2,
    )

    args = ap.parse_args()

    sid = args.structure_id

    if not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid
    ):
        raise ValueError("INVALID_STRUCTURE_ID")

    if (
        not math.isfinite(args.cutoff_angstrom)
        or args.cutoff_angstrom <= 0
        or args.min_residues < 2
    ):
        raise ValueError("INVALID_CLUSTER_PARAMETERS")

    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            f"Refusing existing output directory: {out}"
        )

    residue_path = args.residue_table.resolve()
    pdb_path = args.pdb.resolve()

    for path in (residue_path, pdb_path):
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"MISSING_INPUT:{path}")

        if out in path.parents or path in out.parents:
            raise ValueError("OUTPUT_INPUT_PATH_OVERLAP")

    rows = read_residues(residue_path, sid)

    eligible, pairs, clusters = calculate(
        sid,
        rows,
        pdb_path,
        args.cutoff_angstrom,
        args.min_residues,
    )

    out.mkdir(parents=True, exist_ok=False)

    pairs_file = (
        out / "charged_spatial_pairs_le8A.tsv"
    )

    clusters_file = (
        out / "charged_spatial_clusters.tsv"
    )

    write_tsv(pairs_file, PAIR_COLUMNS, pairs)
    write_tsv(clusters_file, CLUSTER_COLUMNS, clusters)

    meta = {
        "algorithm": (
            "TERMINAL_CHARGE_CENTROID_CONNECTED_COMPONENTS"
        ),
        "method_origin": (
            "INDEPENDENT_RECONSTRUCTION_FROM_HISTORICAL_OUTPUTS"
        ),
        "original_source_code_recovered": False,
        "structure_id": sid,
        "candidate_filter": (
            "STAGE12_PRIORITY_RESIDUE_AND_DEKRH"
        ),
        "cutoff_angstrom": args.cutoff_angstrom,
        "min_residues": args.min_residues,
        "stage12_residue_table_sha256": digest(
            residue_path
        ),
        "pdb_sha256": digest(pdb_path),
        "eligible_residues": len(eligible),
        "pair_count": len(pairs),
        "cluster_count": len(clusters),
        "pairs_sha256": digest(pairs_file),
        "clusters_sha256": digest(clusters_file),
        "historical_output_equivalence": (
            "NOT_CHECKED_BY_THIS_SCRIPT"
        ),
    }

    (
        out / "spatial_reconstruction_provenance.json"
    ).write_text(
        json.dumps(meta, indent=2) + "\n"
    )

    print(f"STRUCTURE_ID={sid}")
    print(f"ELIGIBLE_CHARGED_RESIDUES={len(eligible)}")
    print(f"SPATIAL_PAIR_COUNT={len(pairs)}")
    print(f"SPATIAL_CLUSTER_COUNT={len(clusters)}")
    print("STAGE_L3_RECONSTRUCTED_CALCULATION=PASS")


if __name__ == "__main__":
    main()
