#!/usr/bin/env python3
"""Portable family-MSA preparation, based on the original Stage J."""

import argparse
import csv
import hashlib
import subprocess
from collections import defaultdict
from pathlib import Path


def table(path, required):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError(f"Missing columns: {path}")
        return list(reader)


def fasta(path):
    data = {}
    identifier = None
    chunks = []

    def save():
        if identifier is not None:
            if identifier in data:
                raise ValueError(f"Duplicate FASTA ID: {identifier}")
            data[identifier] = "".join(chunks)

    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                identifier = line[1:].split()[0]
                if not identifier:
                    raise ValueError("Empty FASTA identifier")
                chunks = []
            else:
                if identifier is None:
                    raise ValueError("FASTA sequence without header")
                chunks.append(line)
    save()
    return data


def write_fasta(path, rows):
    with path.open("w") as fh:
        for identifier, sequence in rows:
            fh.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                fh.write(sequence[start:start + 80] + "\n")


def write_table(path, fields, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def to_a3m(alignment):
    widths = {len(sequence) for _, sequence in alignment}
    if len(widths) != 1:
        raise ValueError("Inconsistent MAFFT alignment lengths")

    query = alignment[0][1]
    output = [[] for _ in alignment]

    for col, query_aa in enumerate(query):
        if query_aa != "-":
            for i, (_, seq) in enumerate(alignment):
                aa = seq[col]
                output[i].append(
                    "-" if aa == "-" else aa.upper()
                )
        else:
            for i in range(1, len(alignment)):
                aa = alignment[i][1][col]
                if aa != "-":
                    output[i].append(aa.lower())

    return [
        (identifier, "".join(output[i]))
        for i, (identifier, _) in enumerate(alignment)
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", type=Path, required=True)
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--structure-map", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    if args.threads < 1:
        raise ValueError("threads must be >= 1")

    if args.out_dir.exists():
        raise ValueError("Refusing existing output directory")

    members = table(
        args.members, ("member", "representative")
    )
    structure_map = table(
        args.structure_map,
        ("structure_id", "family_rep", "family_size",
         "sequence_sha256"),
    )

    raw = fasta(args.fasta)
    sequences = {}

    for identifier, sequence in raw.items():
        cleaned = sequence.upper().rstrip("*")
        if not cleaned or not cleaned.isalpha():
            raise ValueError(f"Invalid protein: {identifier}")
        sequences[identifier] = cleaned

    groups = defaultdict(list)
    seen = set()

    for row in members:
        member = row["member"]
        if member in seen:
            raise ValueError(f"Duplicate member: {member}")
        seen.add(member)
        groups[row["representative"]].append(member)

    if set(sequences) != seen:
        raise ValueError("Member FASTA does not match member table")

    if not structure_map:
        raise ValueError("No structure-mapped families")

    ids = [r["structure_id"] for r in structure_map]
    reps = [r["family_rep"] for r in structure_map]

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate structure ID")
    if len(reps) != len(set(reps)):
        raise ValueError("Duplicate representative")

    if not set(reps) <= set(groups):
        raise ValueError("Unknown structure-mapped family")

    version = subprocess.run(
        ["mafft", "--version"],
        capture_output=True,
        text=True,
        check=True,
    )
    mafft_version = (
        version.stdout.strip()
        or version.stderr.strip()
    ).splitlines()[0]

    root = args.out_dir
    folders = {
        name: root / name
        for name in (
            "family_fasta",
            "mafft",
            "a3m",
            "a3m_colabfold",
        )
    }

    for directory in folders.values():
        directory.mkdir(parents=True)

    summary = []
    header_map = []

    for row in sorted(
        structure_map, key=lambda x: x["structure_id"]
    ):
        sid = row["structure_id"]
        rep = row["family_rep"]
        family = groups[rep]

        if rep not in family:
            raise ValueError(f"Missing representative: {rep}")

        if len(family) != int(row["family_size"]):
            raise ValueError(f"Family size mismatch: {sid}")

        rep_sha = hashlib.sha256(
            sequences[rep].encode("ascii")
        ).hexdigest()

        if rep_sha != row["sequence_sha256"]:
            raise ValueError(
                f"Representative SHA-256 mismatch: {sid}"
            )

        # Historical Stage J:
        # 1. Collapse exact sequences in original member order.
        # 2. Move the representative to the first position.
        unique = []
        sequence_positions = {}

        for member in family:
            sequence = sequences[member]

            if sequence in sequence_positions:
                continue

            sequence_positions[sequence] = len(unique)
            unique.append((member, sequence))

        # Portable safety guard:
        # If an earlier member has the identical sequence,
        # retain the representative ID for the query.
        retained_ids = {member for member, _ in unique}

        if rep not in retained_ids:
            position = sequence_positions[sequences[rep]]
            unique[position] = (rep, sequences[rep])

        unique = (
            [record for record in unique if record[0] == rep]
            +
            [record for record in unique if record[0] != rep]
        )

        if not unique or unique[0][0] != rep:
            raise ValueError("Representative was not retained")

        family_faa = (
            folders["family_fasta"] / f"{sid}.family.faa"
        )
        aligned_faa = (
            folders["mafft"] / f"{sid}.aligned.faa"
        )
        mafft_log = (
            folders["mafft"] / f"{sid}.mafft.log"
        )

        write_fasta(family_faa, unique)

        if len(unique) == 1:
            write_fasta(aligned_faa, unique)
            mafft_log.write_text(
                "MAFFT skipped: one unique sequence\n"
            )
            status = "SINGLE_UNIQUE_SEQUENCE"
        else:
            with (
                aligned_faa.open("w") as out,
                mafft_log.open("w") as err,
            ):
                subprocess.run(
                    [
                        "mafft",
                        "--auto",
                        "--thread",
                        str(args.threads),
                        str(family_faa),
                    ],
                    stdout=out,
                    stderr=err,
                    check=True,
                )
            status = "FAMILY_MSA_READY"

        aligned = fasta(aligned_faa)

        if set(aligned) != {identifier for identifier, _ in unique}:
            raise ValueError(f"MAFFT member mismatch: {sid}")

        for identifier, sequence in unique:
            observed = aligned[identifier].replace(
                "-", ""
            ).upper()
            if observed != sequence:
                raise ValueError(
                    f"MAFFT sequence integrity failure: {identifier}"
                )

        # Restore representative-first ordering even if MAFFT
        # changes its output record order.
        aligned_rows = [
            (identifier, aligned[identifier])
            for identifier, _ in unique
        ]

        a3m = to_a3m(aligned_rows)

        query_a3m = a3m[0][1]

        if query_a3m != sequences[rep]:
            raise ValueError(
                f"A3M query reconstruction mismatch: {sid}"
            )

        for identifier, sequence in a3m:
            reconstructed = "".join(
                c for c in sequence if c != "-"
            ).upper()
            if reconstructed != sequences[identifier]:
                raise ValueError(
                    f"A3M homolog reconstruction mismatch: {identifier}"
                )

        write_fasta(
            folders["a3m"] / f"{sid}.a3m",
            a3m,
        )

        safe_rows = []

        for i, (original_id, sequence) in enumerate(a3m):
            safe_id = (
                sid if i == 0 else f"{sid}_h{i:04d}"
            )
            safe_rows.append((safe_id, sequence))
            header_map.append({
                "structure_id": sid,
                "new_header": safe_id,
                "original_header": original_id,
            })

        safe_path = (
            folders["a3m_colabfold"] / f"{sid}.a3m"
        )
        write_fasta(safe_path, safe_rows)

        summary.append({
            "structure_id": sid,
            "family_rep": rep,
            "family_size": len(family),
            "unique_sequences": len(unique),
            "collapsed_duplicates": len(family) - len(unique),
            "msa_status": status,
            "query_length_aa": len(sequences[rep]),
            "a3m_sha256": sha256(safe_path),
            "mafft_version": mafft_version,
        })

        print(
            f"FAMILY={sid} "
            f"members={len(family)} "
            f"unique={len(unique)} "
            f"collapsed={len(family)-len(unique)} "
            f"status={status}"
        )

    write_table(
        root / "family_msa_summary.tsv",
        [
            "structure_id",
            "family_rep",
            "family_size",
            "unique_sequences",
            "collapsed_duplicates",
            "msa_status",
            "query_length_aa",
            "a3m_sha256",
            "mafft_version",
        ],
        summary,
    )

    write_table(
        root / "a3m_colabfold" / "header_map.tsv",
        ["structure_id", "new_header", "original_header"],
        header_map,
    )

    print("STAGE11_FAMILY_MSA=PASS")
    print(f"PROCESSED_FAMILIES={len(summary)}")
    print(
        "FAMILIES_WITH_MULTIPLE_UNIQUE_SEQUENCES="
        + str(sum(
            r["unique_sequences"] > 1 for r in summary
        ))
    )
    print("COLABFOLD_PREDICTION=NOT_RUN")


if __name__ == "__main__":
    main()
