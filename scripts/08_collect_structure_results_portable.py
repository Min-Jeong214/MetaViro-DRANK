#!/usr/bin/env python3
"""Validate ColabFold structures and link them to viral families."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


FIELDS = [
    "structure_id",
    "family_rep",
    "source_member",
    "family_size",
    "priority_class",
    "length_aa",
    "sequence_sha256",
    "pdb_file",
    "pdb_sha256",
    "score_json",
    "score_json_sha256",
    "mean_plddt",
    "min_plddt",
    "max_plddt",
    "ptm",
    "confidence_class",
    "prediction_status",
]


def sha256_file(path):
    digest = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(4 * 1024 * 1024), b""
        ):
            digest.update(block)

    return digest.hexdigest()


def read_mapping(path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "structure_id",
            "family_rep",
            "source_member",
            "family_size",
            "length_aa",
            "sequence_sha256",
        }

        if not required <= set(reader.fieldnames or []):
            raise ValueError("Missing structure mapping columns")

        rows = list(reader)

    if not rows:
        raise ValueError("Empty structure mapping")

    ids = [r["structure_id"] for r in rows]

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate structure IDs")

    reps = [r["family_rep"] for r in rows]

    if len(reps) != len(set(reps)):
        raise ValueError("Duplicate family representatives")

    return rows


def select_pdb(directory, structure_id):
    matches = sorted(
        directory.glob(
            f"{structure_id}_*rank_001_*.pdb"
        )
    )

    relaxed = [
        p for p in matches
        if p.name.startswith(
            f"{structure_id}_relaxed_rank_001_"
        )
    ]

    if len(relaxed) == 1:
        return relaxed[0]

    if len(matches) == 1:
        return matches[0]

    raise ValueError(
        f"Expected one rank-001 PDB for {structure_id}; "
        f"found {len(matches)}"
    )


def read_pdb(path):
    structure = PDBParser(
        QUIET=True
    ).get_structure("prediction", path)

    models = list(structure.get_models())

    if len(models) != 1:
        raise ValueError("Expected one PDB model")

    chains = list(models[0].get_chains())

    if len(chains) != 1:
        raise ValueError("Expected one protein chain")

    residues = [
        residue
        for residue in chains[0].get_residues()
        if residue.id[0] == " "
    ]

    if not residues:
        raise ValueError("PDB contains no protein residues")

    if any("CA" not in residue for residue in residues):
        raise ValueError("PDB residue missing CA atom")

    sequence = "".join(
        seq1(residue.get_resname(), custom_map={"UNK": "X"})
        for residue in residues
    ).upper()

    plddt = [
        float(residue["CA"].get_bfactor())
        for residue in residues
    ]

    if any(
        not math.isfinite(value)
        or not 0 <= value <= 100
        for value in plddt
    ):
        raise ValueError("Invalid PDB pLDDT value")

    return sequence, plddt


def read_scores(directory, structure_id, length):
    matches = sorted(
        directory.glob(
            f"{structure_id}_scores_rank_001_*.json"
        )
    )

    if not matches:
        return "", "", None

    if len(matches) != 1:
        raise ValueError(
            f"Ambiguous score JSON: {structure_id}"
        )

    path = matches[0]
    scores = json.loads(path.read_text())

    ptm = scores.get("ptm")

    if ptm is not None:
        ptm = float(ptm)

        if not math.isfinite(ptm) or not 0 <= ptm <= 1:
            raise ValueError("Invalid pTM score")

    per_residue = scores.get("plddt")

    if per_residue is not None:
        if len(per_residue) != length:
            raise ValueError(
                "Score JSON pLDDT length mismatch"
            )

    return str(path.resolve()), sha256_file(path), ptm


def confidence_class(plddt):
    if plddt < 50:
        return "VERY_LOW"
    if plddt < 70:
        return "LOW"
    if plddt < 90:
        return "CONFIDENT"
    return "VERY_HIGH"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument(
        "--prediction-dir", type=Path, required=True
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.out_dir.exists():
        raise ValueError(
            "Refusing to overwrite existing output directory"
        )

    mapping = read_mapping(args.mapping)

    queries = {}

    for record in SeqIO.parse(args.queries, "fasta"):
        if record.id in queries:
            raise ValueError("Duplicate query FASTA ID")

        queries[record.id] = str(
            record.seq
        ).upper().rstrip("*")

    expected_ids = {
        row["structure_id"] for row in mapping
    }

    if set(queries) != expected_ids:
        raise ValueError("Query FASTA/mapping ID mismatch")

    output = []

    for row in mapping:
        sid = row["structure_id"]
        sequence = queries[sid]

        digest = hashlib.sha256(
            sequence.encode("ascii")
        ).hexdigest()

        if digest != row["sequence_sha256"]:
            raise ValueError(
                f"Sequence SHA-256 mismatch: {sid}"
            )

        if len(sequence) != int(row["length_aa"]):
            raise ValueError(
                f"Sequence length mismatch: {sid}"
            )

        if row["source_member"] != row["family_rep"]:
            raise ValueError(
                f"Unexpected representative mapping: {sid}"
            )

        pdb_path = select_pdb(
            args.prediction_dir, sid
        )

        pdb_sequence, plddt = read_pdb(pdb_path)

        if pdb_sequence != sequence:
            raise ValueError(
                f"PDB amino acid sequence mismatch: {sid}"
            )

        score_file, score_digest, ptm = read_scores(
            args.prediction_dir,
            sid,
            len(sequence),
        )

        mean = sum(plddt) / len(plddt)
        quality = confidence_class(mean)

        output.append({
            "structure_id": sid,
            "family_rep": row["family_rep"],
            "source_member": row["source_member"],
            "family_size": row["family_size"],
            "priority_class": row.get(
                "priority_class", ""
            ),
            "length_aa": len(sequence),
            "sequence_sha256": digest,
            "pdb_file": str(pdb_path.resolve()),
            "pdb_sha256": sha256_file(pdb_path),
            "score_json": score_file,
            "score_json_sha256": score_digest,
            "mean_plddt": f"{mean:.2f}",
            "min_plddt": f"{min(plddt):.2f}",
            "max_plddt": f"{max(plddt):.2f}",
            "ptm": (
                f"{ptm:.4f}"
                if ptm is not None else "NA"
            ),
            "confidence_class": quality,
            "prediction_status": "PREDICTED",
        })

    # Write only after all structures pass validation.
    args.out_dir.mkdir(parents=True)

    destination = args.out_dir / "structure_results.tsv"

    with destination.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output)

    low = sum(
        r["confidence_class"] in ("VERY_LOW", "LOW")
        for r in output
    )

    print("STAGE08_STRUCTURE_COLLECTION=PASS")
    print(f"VALIDATED_STRUCTURES={len(output)}")
    print(f"LOW_CONFIDENCE_STRUCTURES={low}")
    print(f"STRUCTURE_RESULTS={destination}")

    for row in output:
        print(
            f"STRUCTURE={row['structure_id']} "
            f"pLDDT={row['mean_plddt']} "
            f"pTM={row['ptm']} "
            f"QUALITY={row['confidence_class']}"
        )


if __name__ == "__main__":
    main()
