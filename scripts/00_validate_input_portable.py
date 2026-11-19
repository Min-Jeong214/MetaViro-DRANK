#!/usr/bin/env python3
"""Validate MetaViro-DRANK sample input consistency."""

import argparse
import re
import sys
from pathlib import Path


def validate(sample):
    faa = sample / "proteins.faa"
    gff = sample / "genes.gff"
    eggnog = sample / "annotations.emapper.annotations"

    for path in (faa, gff, eggnog):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty: {path}")

    # --------------------------------------------------
    # 1. FASTA
    # --------------------------------------------------
    proteins = set()
    current = None
    length = 0

    with faa.open() as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current is not None and length == 0:
                    raise ValueError(
                        f"Empty protein sequence: {current}"
                    )

                pid = line[1:].split()

                if not pid:
                    raise ValueError("Empty FASTA identifier")

                current = pid[0]

                if current in proteins:
                    raise ValueError(
                        f"Duplicate FASTA ID: {current}"
                    )

                proteins.add(current)
                length = 0

            else:
                if current is None:
                    raise ValueError(
                        "FASTA sequence before first header"
                    )

                if not re.fullmatch(
                    r"[A-Za-z*]+", line
                ):
                    raise ValueError(
                        f"Invalid protein sequence: {current}"
                    )

                length += len(line)

    if current is not None and length == 0:
        raise ValueError(
            f"Empty protein sequence: {current}"
        )

    if not proteins:
        raise ValueError("No FASTA proteins")

    # --------------------------------------------------
    # 2. GFF / Prodigal ORF mapping
    # --------------------------------------------------
    gff_ids = set()

    with gff.open() as fh:
        for number, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"GFF line {number}: expected 9 columns"
                )

            if fields[2] != "CDS":
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                raise ValueError(
                    f"GFF line {number}: invalid coordinates"
                )

            if start < 1 or end < start:
                raise ValueError(
                    f"GFF line {number}: invalid interval"
                )

            match = re.search(
                r"(?:^|;)ID=([^;]+)", fields[8]
            )

            if not match:
                raise ValueError(
                    f"GFF line {number}: missing CDS ID"
                )

            raw_id = match.group(1)

            try:
                index = int(raw_id.rsplit("_", 1)[1])
            except (IndexError, ValueError):
                raise ValueError(
                    f"GFF line {number}: "
                    f"invalid Prodigal ID {raw_id}"
                )

            pid = f"{fields[0]}_{index}"

            if pid in gff_ids:
                raise ValueError(
                    f"Duplicate GFF protein ID: {pid}"
                )

            if pid not in proteins:
                raise ValueError(
                    f"GFF protein ID absent from FASTA: {pid}"
                )

            gff_ids.add(pid)

    missing = proteins - gff_ids

    if missing:
        raise ValueError(
            "FASTA proteins absent from GFF: "
            + ", ".join(sorted(missing)[:5])
        )

    # --------------------------------------------------
    # 3. eggNOG-mapper annotations
    # --------------------------------------------------
    required = {
        "query",
        "eggNOG_OGs",
        "Description",
        "Preferred_name",
    }

    header = None
    query_index = None
    annotated = set()

    with eggnog.open() as fh:
        for number, line in enumerate(fh, 1):
            if line.startswith("#query"):
                if header is not None:
                    raise ValueError(
                        "Duplicate eggNOG header"
                    )

                header = line[1:].rstrip("\n").split("\t")

                if not required.issubset(header):
                    raise ValueError(
                        "Missing eggNOG columns: "
                        + str(sorted(required - set(header)))
                    )

                query_index = header.index("query")
                continue

            if not line.strip() or line.startswith("#"):
                continue

            if header is None:
                raise ValueError(
                    "eggNOG data before #query header"
                )

            fields = line.rstrip("\n").split("\t")

            if len(fields) != len(header):
                raise ValueError(
                    f"eggNOG line {number}: "
                    "column count mismatch"
                )

            pid = fields[query_index]

            if pid not in proteins:
                raise ValueError(
                    f"eggNOG ID absent from FASTA: {pid}"
                )

            if pid in annotated:
                raise ValueError(
                    f"Duplicate eggNOG query: {pid}"
                )

            annotated.add(pid)

    if header is None:
        raise ValueError(
            "Missing eggNOG #query header"
        )

    print(f"SAMPLE={sample.name}")
    print(f"FASTA_PROTEINS={len(proteins)}")
    print(f"GFF_MATCHED_CDS={len(gff_ids)}")
    print(f"EGGNOG_ANNOTATED={len(annotated)}")
    print(
        f"EGGNOG_WITHOUT_HIT={len(proteins - annotated)}"
    )
    print("INPUT_VALIDATION=PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--sample-dir", required=True, type=Path
    )
    args = ap.parse_args()

    try:
        validate(args.sample_dir)
    except (ValueError, OSError) as exc:
        sys.exit(f"INPUT_VALIDATION=FAIL: {exc}")


if __name__ == "__main__":
    main()
