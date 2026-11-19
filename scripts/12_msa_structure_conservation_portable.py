#!/usr/bin/env python3
"""Portable Stage 12 using the original residue conservation functions."""

import argparse
import csv
import hashlib
import json
import re
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = ROOT / "scripts/12_msa_structure_conservation.py"

SUMMARY_FIELDS = [
    "structure_id",
    "length",
    "msa_nseq",
    "very_conserved_positions",
    "priority_structural_positions",
    "priority_acidic",
    "priority_basic",
    "priority_nucleophilic",
]

WINDOW_FIELDS = [
    "structure_id",
    "start",
    "end",
    "sequence",
    "conserved_structural",
    "charged",
    "acidic",
    "basic",
    "mean_query_fraction",
    "mean_plddt",
]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for part in iter(lambda: fh.read(1048576), b""):
            h.update(part)
    return h.hexdigest()


def read_table(path, required):
    if not path.is_file():
        raise ValueError(f"MISSING_INPUT:{path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        if len(fields) != len(set(fields)):
            raise ValueError(f"DUPLICATE_COLUMNS:{path}")

        if not set(required).issubset(fields):
            raise ValueError(f"MISSING_COLUMNS:{path}")

        rows = list(reader)

    if not rows:
        raise ValueError(f"EMPTY_INPUT:{path}")

    return rows


def fasta(path):
    if not path.is_file():
        raise ValueError(f"MISSING_ALIGNMENT:{path}")

    records = []
    name = None
    chunks = []

    for line in path.read_text().splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))

            name = line[1:].split()[0]
            chunks = []
        else:
            if name is None:
                raise ValueError(f"INVALID_FASTA:{path}")
            chunks.append(line)

    if name is not None:
        records.append((name, "".join(chunks)))

    if not records:
        raise ValueError(f"EMPTY_FASTA:{path}")

    return records


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


def validate_family(sid, aln_dir, a3m_dir, record_path):
    aln = aln_dir / f"{sid}.aligned.faa"
    a3m = a3m_dir / f"{sid}.a3m"

    if not record_path.is_file():
        raise ValueError(f"MISSING_PREDICTION_RECORD:{sid}")

    record = json.loads(record_path.read_text())

    if (
        record.get("status") != "COMPLETED"
        or record.get("structure_id") != sid
    ):
        raise ValueError(f"UNVERIFIED_PREDICTION:{sid}")

    if not aln.is_file() or not a3m.is_file():
        raise ValueError(f"MISSING_FAMILY_MSA:{sid}")

    if sha256(a3m) != record.get("a3m_sha256"):
        raise ValueError(f"A3M_SHA_MISMATCH:{sid}")

    aligned = fasta(aln)
    msa = fasta(a3m)

    if (
        len(aligned) != len(msa)
        or len(msa) != record["sequence_count"]
        or len(msa) < 2
    ):
        raise ValueError(f"MSA_SEQUENCE_COUNT_MISMATCH:{sid}")

    query = aligned[0][1].replace("-", "")

    if (
        msa[0][0] != sid
        or query != msa[0][1]
        or len(query) != record["query_length"]
    ):
        raise ValueError(f"MSA_QUERY_MISMATCH:{sid}")

    query_sha = hashlib.sha256(
        query.encode()
    ).hexdigest()

    if query_sha != record["query_sha256"]:
        raise ValueError(f"QUERY_SHA_MISMATCH:{sid}")

    pdb = Path(record["pdb_file"]).resolve()
    scores = Path(record["score_json"]).resolve()

    if (
        not pdb.is_file()
        or sha256(pdb) != record["pdb_sha256"]
    ):
        raise ValueError(f"PDB_SHA_MISMATCH:{sid}")

    if (
        not scores.is_file()
        or sha256(scores) != record["score_json_sha256"]
    ):
        raise ValueError(f"SCORE_SHA_MISMATCH:{sid}")

    return aln, a3m, pdb, record


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--aln-dir", type=Path, required=True)
    ap.add_argument("--a3m-dir", type=Path, required=True)
    ap.add_argument(
        "--prediction-index", type=Path, required=True
    )
    ap.add_argument("--out-dir", type=Path, required=True)

    args = ap.parse_args()

    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            f"Refusing existing output directory: {output}"
        )

    selection = read_table(
        args.selection, ("structure_id",)
    )

    selected = [r["structure_id"] for r in selection]

    if len(selected) != len(set(selected)):
        raise ValueError("DUPLICATE_SELECTED_FAMILIES")

    for sid in selected:
        if not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid
        ):
            raise ValueError(f"INVALID_STRUCTURE_ID:{sid}")

    index_rows = read_table(
        args.prediction_index,
        ("structure_id", "prediction_provenance"),
    )

    index = {}

    for row in index_rows:
        sid = row["structure_id"]

        if sid in index:
            raise ValueError(f"DUPLICATE_PREDICTION_ID:{sid}")

        index[sid] = Path(
            row["prediction_provenance"]
        ).resolve()

    original = runpy.run_path(
        str(ORIGINAL),
        run_name="__portable_stage12__",
    )

    build_residues = original["build_residue_table"]
    build_summary = original["build_summary"]
    make_windows = original["make_windows"]
    write_residues = original["write_residue_table"]

    analyzed = []
    summaries = []
    windows = []

    # Validate and calculate all Families before creating output.
    for sid in selected:
        if sid not in index:
            raise ValueError(
                f"MISSING_PREDICTION_INDEX:{sid}"
            )

        aln, a3m, pdb, record = validate_family(
            sid,
            args.aln_dir.resolve(),
            args.a3m_dir.resolve(),
            index[sid],
        )

        residues = build_residues(sid, aln, pdb)

        if len(residues) != record["query_length"]:
            raise ValueError(
                f"RESIDUE_LENGTH_MISMATCH:{sid}"
            )

        summary = build_summary(sid, residues)
        family_windows = make_windows(sid, residues)

        analyzed.append({
            "structure_id": sid,
            "alignment": str(aln),
            "alignment_sha256": sha256(aln),
            "a3m": str(a3m),
            "a3m_sha256": sha256(a3m),
            "prediction_provenance": str(index[sid]),
            "prediction_provenance_sha256": sha256(
                index[sid]
            ),
            "pdb_sha256": sha256(pdb),
            "sequence_count": len(fasta(a3m)),
            "residue_count": len(residues),
            "residues": residues,
        })

        summaries.append(summary)
        windows.extend(family_windows)

    output.mkdir(parents=True, exist_ok=False)

    for family in analyzed:
        sid = family["structure_id"]

        residue_file = (
            output / f"{sid}_residue_conservation.tsv"
        )

        write_residues(
            residue_file,
            family.pop("residues"),
        )

        family["residue_table_sha256"] = sha256(
            residue_file
        )

    summary_path = (
        output / "candidate_conservation_summary.tsv"
    )

    write_table(
        summary_path, SUMMARY_FIELDS, summaries
    )

    window_path = output / "conserved_windows.tsv"

    formatted_windows = []

    for row in windows:
        formatted_windows.append({
            **row,
            "mean_query_fraction":
                f'{row["mean_query_fraction"]:.3f}',
            "mean_plddt":
                f'{row["mean_plddt"]:.2f}',
        })

    write_table(
        window_path,
        WINDOW_FIELDS,
        formatted_windows,
    )

    provenance = {
        "algorithm": "ORIGINAL_STAGE12_FUNCTIONS",
        "original_script_sha256": sha256(ORIGINAL),
        "selection_sha256": sha256(args.selection),
        "prediction_index_sha256": sha256(
            args.prediction_index
        ),
        "selection_source": "USER_SUPPLIED_CANDIDATE_IDS",
        "final_top5_automatically_selected": False,
        "family_count": len(analyzed),
        "families": analyzed,
        "summary_sha256": sha256(summary_path),
        "windows_sha256": sha256(window_path),
    }

    (output / "conservation_provenance.json").write_text(
        json.dumps(
            provenance,
            ensure_ascii=False,
            indent=2,
        ) + "\n"
    )

    print(f"ANALYZED_FAMILIES={len(analyzed)}")
    print(f"RESIDUE_TABLES={len(analyzed)}")
    print(f"CONSERVATION_WINDOWS={len(windows)}")
    print(f"OUTPUT={output}")
    print("ORIGINAL_STAGE12_FUNCTIONS=REUSED")
    print("STAGE12_PORTABLE_CONSERVATION=PASS")


if __name__ == "__main__":
    main()
