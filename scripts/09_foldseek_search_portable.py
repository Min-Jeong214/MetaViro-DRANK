#!/usr/bin/env python3
"""Portable and non-destructive Foldseek structure search."""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


FIELDS = (
    "query target evalue bits fident alnlen qlen tlen "
    "qcov tcov lddt alntmscore qtmscore ttmscore prob"
).split()

PRESETS = {
    "pdb_initial": (1e-3, 200),
    "pdb_e1": (1.0, 1000),
    "afsp": (1.0, 1000),
    "bfvd": (1.0, 500),
    "control_self": (1.0, 20),
}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {
            "structure_id", "pdb_file",
            "pdb_sha256", "prediction_status",
        }

        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Stage 08 structure table schema mismatch")

        rows = list(reader)

    if not rows:
        raise ValueError("Empty Stage 08 structure table")

    ids = [r["structure_id"] for r in rows]

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate structure IDs")

    for sid in ids:
        if not sid or not all(
            c.isalnum() or c in "_-" for c in sid
        ):
            raise ValueError(f"Unsafe structure ID: {sid}")

    return rows


def run(command, logfile):
    with logfile.open("w") as fh:
        result = subprocess.run(
            command,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}). "
            f"See {logfile}"
        )


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


def strong_hit(hit):
    # Preserve original Stage G3 evidence thresholds.
    return (
        float(hit["prob"]) >= 0.90
        and float(hit["lddt"]) >= 0.60
        and float(hit["qcov"]) >= 0.30
        and float(hit["tcov"]) >= 0.50
        and float(hit["ttmscore"]) >= 0.50
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--structures", type=Path, required=True
    )
    parser.add_argument(
        "--target-db", type=Path, required=True
    )
    parser.add_argument(
        "--source", choices=PRESETS, required=True
    )
    parser.add_argument(
        "--out-dir", type=Path, required=True
    )
    parser.add_argument(
        "--threads", type=int, default=2
    )

    args = parser.parse_args()

    if args.threads < 1:
        raise ValueError("Threads must be >= 1")

    output = args.out_dir.resolve()
    target = args.target_db.resolve()

    if output.exists():
        raise ValueError(
            f"Refusing existing output directory: {output}"
        )

    if not Path(str(target) + ".dbtype").is_file():
        raise ValueError(
            f"Foldseek target database missing: {target}"
        )

    structures = read_table(args.structures)

    # Validate every input before creating output files.
    checked = []

    for row in structures:
        sid = row["structure_id"]
        pdb = Path(row["pdb_file"]).resolve()

        if row["prediction_status"] != "PREDICTED":
            raise ValueError(
                f"Structure not predicted: {sid}"
            )

        if not pdb.is_file() or not pdb.stat().st_size:
            raise ValueError(
                f"PDB missing or empty: {pdb}"
            )

        actual = digest(pdb)

        if actual != row["pdb_sha256"]:
            raise ValueError(
                f"PDB SHA-256 mismatch: {sid}"
            )

        checked.append((sid, pdb, actual))

    version = subprocess.run(
        ["foldseek", "version"],
        check=True,
        capture_output=True,
        text=True,
    )

    foldseek_version = (
        version.stdout.strip()
        or version.stderr.strip()
    ).splitlines()[0]

    output.mkdir(parents=True)

    pdb_dir = output / "query_pdb"
    pdb_dir.mkdir()

    query_manifest = []

    for sid, pdb, sha in checked:
        copy = pdb_dir / f"{sid}.pdb"

        # Copy, never move or modify the original PDB.
        shutil.copyfile(pdb, copy)

        if digest(copy) != sha:
            raise RuntimeError(
                f"PDB copy integrity failure: {sid}"
            )

        query_manifest.append({
            "structure_id": sid,
            "source_pdb": str(pdb),
            "query_pdb": str(copy),
            "sha256": sha,
        })

    write_tsv(
        output / "query_manifest.tsv",
        [
            "structure_id",
            "source_pdb",
            "query_pdb",
            "sha256",
        ],
        query_manifest,
    )

    qdb = output / "query_db"
    rdb = output / "result_db"
    tmp = output / "tmp"
    raw = output / "raw_hits.tsv"

    evalue, max_seqs = PRESETS[args.source]

    print(
        f"FOLDSEEK_VERSION={foldseek_version}",
        flush=True,
    )
    print(
        f"SOURCE_PROFILE={args.source}",
        flush=True,
    )
    print(
        f"EVALUE={evalue} MAX_SEQS={max_seqs}",
        flush=True,
    )

    run(
        [
            "foldseek", "createdb",
            str(pdb_dir), str(qdb),
        ],
        output / "createdb.log",
    )

    print("FOLDSEEK_QUERY_DB=PASS", flush=True)

    run(
        [
            "foldseek", "search",
            str(qdb),
            str(target),
            str(rdb),
            str(tmp),
            "-s", "9.5",
            "-e", str(evalue),
            "--max-seqs", str(max_seqs),
            "--threads", str(args.threads),
            "-a",
        ],
        output / "search.log",
    )

    print("FOLDSEEK_SEARCH=PASS", flush=True)

    run(
        [
            "foldseek", "convertalis",
            str(qdb),
            str(target),
            str(rdb),
            str(raw),
            "--format-output",
            ",".join(FIELDS),
        ],
        output / "convertalis.log",
    )

    if not raw.is_file():
        raise RuntimeError(
            "Foldseek raw result file missing"
        )

    hits = defaultdict(list)
    normalized = []

    known_ids = [
        row["structure_id"] for row in structures
    ]

    with raw.open(newline="") as fh:
        for lineno, parts in enumerate(
            csv.reader(fh, delimiter="\t"), 1
        ):
            if not parts:
                continue

            if len(parts) != len(FIELDS):
                raise ValueError(
                    f"Foldseek output line {lineno}: "
                    f"expected {len(FIELDS)} columns, "
                    f"got {len(parts)}"
                )

            hit = dict(zip(FIELDS, parts))

            original_query = hit["query"]

            matches = [
                sid for sid in known_ids
                if original_query == sid
                or original_query.startswith(sid + "_")
                or original_query.startswith(sid + ".")
            ]

            if len(matches) != 1:
                raise ValueError(
                    "Unrecognized/ambiguous Foldseek query: "
                    + original_query
                )

            sid = matches[0]
            hit["query"] = sid

            # Check all numbers needed by Stage G3.
            for field in (
                "prob", "lddt", "qcov", "tcov",
                "ttmscore", "qtmscore", "evalue",
            ):
                float(hit[field])

            hits[sid].append(hit)
            normalized.append(hit)

    # Output matches the historical 15-column
    # Foldseek raw format, without a header line.
    with raw.open("w", newline="") as fh:
        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n",
        )

        for hit in normalized:
            writer.writerow([
                hit[field] for field in FIELDS
            ])

    summary = []

    for sid in known_ids:
        family_hits = hits[sid]

        summary.append({
            "structure_id": sid,
            "source": args.source,
            "raw_count": len(family_hits),
            "convincing_count": sum(
                strong_hit(h) for h in family_hits
            ),
            "status": (
                "CONTROL_ONLY"
                if args.source == "control_self"
                else (
                    "RAW_HITS_PRESENT"
                    if family_hits else "NO_RAW_HIT"
                )
            ),
        })

    write_tsv(
        output / "source_summary.tsv",
        [
            "structure_id", "source",
            "raw_count", "convincing_count",
            "status",
        ],
        summary,
    )

    provenance = {
        "source_profile": args.source,
        "foldseek_version": foldseek_version,
        "sensitivity": 9.5,
        "evalue": evalue,
        "max_seqs": max_seqs,
        "threads": args.threads,
        "target_db": str(target),
        "target_dbtype_sha256": digest(
            Path(str(target) + ".dbtype")
        ),
        "input_structure_table": str(
            args.structures.resolve()
        ),
        "query_manifest": query_manifest,
        "raw_hits_sha256": digest(raw),
        "control_only": args.source == "control_self",
    }

    (output / "search_provenance.json").write_text(
        json.dumps(
            provenance,
            indent=2,
            ensure_ascii=False,
        ) + "\n"
    )

    print(f"STRUCTURE_COUNT={len(known_ids)}")
    print(f"RAW_HIT_COUNT={len(normalized)}")
    print(f"RAW_HITS={raw}")
    print("STAGE09_PORTABLE_SEARCH=PASS")


if __name__ == "__main__":
    main()
