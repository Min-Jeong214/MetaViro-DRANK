#!/usr/bin/env python3
"""Extract all members of selected recurrent protein families."""
import argparse
import csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()

    for name in (
        "cluster-tsv", "family-tsv", "proteins-faa",
        "out-faa", "out-members"
    ):
        ap.add_argument("--" + name, required=True, type=Path)

    args = ap.parse_args()

    with args.family_tsv.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "representative" not in (reader.fieldnames or []):
            raise ValueError("Missing representative column")
        selected = [r["representative"].strip() for r in reader]

    if len(selected) != len(set(selected)) or not all(selected):
        raise ValueError("Invalid family representatives")

    wanted = set(selected)
    members = {}
    observed = set()

    with args.cluster_tsv.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("Invalid cluster TSV")
            rep, member = fields[:2]
            if rep not in wanted:
                continue
            observed.add(rep)
            if member in members and members[member] != rep:
                raise ValueError(f"Conflicting family: {member}")
            members[member] = rep

    if observed != wanted:
        raise ValueError("Selected families missing from clusters")

    args.out_faa.parent.mkdir(parents=True, exist_ok=True)
    args.out_members.parent.mkdir(parents=True, exist_ok=True)

    found = set()

    with args.proteins_faa.open() as src, \
         args.out_faa.open("w") as dst:

        include = False

        for line in src:
            if line.startswith(">"):
                parts = line[1:].split()
                if not parts:
                    raise ValueError("Empty FASTA identifier")
                pid = parts[0]
                include = pid in members
                if include:
                    if pid in found:
                        raise ValueError(f"Duplicate FASTA ID: {pid}")
                    found.add(pid)

            if include:
                dst.write(line if line.endswith("\n") else line + "\n")

    if found != set(members):
        raise ValueError("Selected proteins missing from FASTA")

    with args.out_members.open("w", newline="") as fh:
        writer = csv.writer(
            fh, delimiter="\t", lineterminator="\n"
        )
        writer.writerow(["member", "representative"])
        for member in sorted(members):
            writer.writerow([member, members[member]])

    print(f"SELECTED_FAMILIES={len(wanted)}")
    print(f"SELECTED_MEMBERS={len(members)}")
    print("STAGE05_QUERY_PREPARATION=PASS")


if __name__ == "__main__":
    main()
