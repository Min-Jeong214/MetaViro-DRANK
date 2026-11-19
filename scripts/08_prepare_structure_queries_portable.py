#!/usr/bin/env python3
"""Prepare one structure-prediction query per viral family."""

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


def read_tsv(path, required):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path}: missing columns {sorted(missing)}"
            )
        return list(reader)


def read_fasta(path):
    sequences = {}
    current = None
    chunks = []

    def save():
        if current is None:
            return
        if current in sequences:
            raise ValueError(f"Duplicate FASTA ID: {current}")
        seq = "".join(chunks).upper().rstrip("*")
        if not seq or not seq.isalpha():
            raise ValueError(f"Invalid protein sequence: {current}")
        sequences[current] = seq

    with path.open() as fh:
        for raw in fh:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                save()
                current = line[1:].split()[0]
                if not current:
                    raise ValueError("Empty FASTA ID")
                chunks = []
            else:
                if current is None:
                    raise ValueError("Sequence before FASTA header")
                chunks.append(line)

    save()
    return sequences


def write_fasta(path, records):
    with path.open("w") as fh:
        for identifier, sequence in records:
            fh.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                fh.write(sequence[start:start + 80] + "\n")


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", type=Path, required=True)
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--family-audit", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    if args.out_dir.exists():
        raise ValueError(
            f"Refusing existing output directory: {args.out_dir}"
        )

    members = read_tsv(
        args.members,
        ("member", "representative"),
    )

    families = read_tsv(
        args.family_audit,
        ("representative", "family_size"),
    )

    sequences = read_fasta(args.fasta)

    if not members or not families:
        raise ValueError("Empty family or member table")

    groups = defaultdict(set)
    seen_members = set()

    for row in members:
        member = row["member"]
        rep = row["representative"]

        if not member or not rep:
            raise ValueError("Empty member or representative ID")

        if member in seen_members:
            raise ValueError(f"Duplicate member: {member}")

        seen_members.add(member)
        groups[rep].add(member)

    if set(sequences) != seen_members:
        missing = seen_members - set(sequences)
        extra = set(sequences) - seen_members
        raise ValueError(
            f"FASTA/member mismatch: missing={len(missing)}, "
            f"extra={len(extra)}"
        )

    audit_by_rep = {}

    for row in families:
        rep = row["representative"]
        if rep in audit_by_rep:
            raise ValueError(f"Duplicate family: {rep}")
        audit_by_rep[rep] = row

    if set(groups) != set(audit_by_rep):
        raise ValueError("Stage05/Stage06 family mismatch")

    records = []
    mapping = []

    for index, rep in enumerate(sorted(groups), start=1):
        group = groups[rep]
        audit = audit_by_rep[rep]

        expected_size = int(audit["family_size"])

        if expected_size != len(group):
            raise ValueError(
                f"Family size mismatch: {rep}"
            )

        if rep not in group:
            raise ValueError(
                f"Representative absent from family: {rep}"
            )

        sequence = sequences[rep]
        structure_id = f"MVF{index:06d}"

        sequence_sha256 = hashlib.sha256(
            sequence.encode("ascii")
        ).hexdigest()

        records.append((structure_id, sequence))

        mapping.append({
            "structure_id": structure_id,
            "family_rep": rep,
            "source_member": rep,
            "family_size": expected_size,
            "length_aa": len(sequence),
            "sequence_sha256": sequence_sha256,
            "priority_class": audit.get(
                "priority_class", ""
            ),
            "prediction_status": "NOT_RUN",
        })

    args.out_dir.mkdir(parents=True)

    write_fasta(
        args.out_dir / "representative_queries.faa",
        records,
    )

    write_tsv(
        args.out_dir / "structure_id_map.tsv",
        mapping,
        [
            "structure_id",
            "family_rep",
            "source_member",
            "family_size",
            "length_aa",
            "sequence_sha256",
            "priority_class",
            "prediction_status",
        ],
    )

    print("STAGE08_QUERY_PREPARATION=PASS")
    print(f"INPUT_MEMBERS={len(members)}")
    print(f"INPUT_FAMILIES={len(families)}")
    print(f"STRUCTURE_QUERIES={len(records)}")
    print("PREDICTION_STATUS=NOT_RUN")


if __name__ == "__main__":
    main()
