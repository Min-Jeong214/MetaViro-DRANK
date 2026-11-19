#!/usr/bin/env python3
"""Compare verified Stage 08 and Stage 11 structures without GPU reruns."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


FIELDS = [
    "structure_id",
    "family_rep",
    "single_mean_plddt",
    "msa_mean_plddt",
    "delta_mean_plddt",
    "single_ptm",
    "msa_ptm",
    "delta_ptm",
    "msa_sequences",
    "query_length",
    "model_type",
    "single_pdb_sha256",
    "msa_pdb_sha256",
]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def table(path, required):
    if not path.is_file():
        raise ValueError(f"MISSING_FILE:{path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        if len(fields) != len(set(fields)):
            raise ValueError(f"DUPLICATE_COLUMNS:{path}")

        if not set(required).issubset(fields):
            raise ValueError(f"MISSING_COLUMNS:{path}")

        rows = list(reader)

    if not rows:
        raise ValueError(f"EMPTY_TABLE:{path}")

    return rows


def indexed(rows, key):
    found = {}

    for row in rows:
        value = row[key]

        if not value or value in found:
            raise ValueError(f"INVALID_OR_DUPLICATE_ID:{value}")

        found[value] = row

    return found


def checked_file(filename, checksum):
    path = Path(filename).resolve()

    if not path.is_file() or sha(path) != checksum:
        raise ValueError(f"FILE_INTEGRITY_FAILURE:{path}")

    return path


def numeric(value, label, low, high):
    number = float(value)

    if not math.isfinite(number) or not low <= number <= high:
        raise ValueError(f"INVALID_METRIC:{label}")

    return number


def pdb_mean(path, expected_length):
    scores = []
    residues = set()

    with path.open(errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue

            residue = (line[21], line[22:27])

            if residue in residues:
                raise ValueError(f"DUPLICATE_CA:{path}")

            residues.add(residue)
            scores.append(
                numeric(
                    line[60:66],
                    "pLDDT",
                    0,
                    100,
                )
            )

    if len(scores) != expected_length:
        raise ValueError(f"PDB_LENGTH_MISMATCH:{path}")

    return sum(scores) / len(scores)


def query_sha_from_a3m(path):
    header = None
    sequence = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()

            if not line or line.startswith("#"):
                # ColabFold A3M may start with a # metadata header.
                # Skip metadata without changing the actual A3M bytes.
                continue

            if line.startswith(">"):
                if header is not None:
                    break

                header = line[1:].split()[0]
            else:
                if header is None:
                    raise ValueError("INVALID_A3M_HEADER")
                sequence.append(line)

    query = "".join(sequence)

    if not query or not query.isalpha() or not query.isupper():
        raise ValueError("INVALID_A3M_QUERY")

    return (
        header,
        len(query),
        hashlib.sha256(query.encode()).hexdigest(),
    )


def collect(single, record_path):
    sid = single["structure_id"]

    if not record_path.is_file():
        raise ValueError(f"MISSING_PREDICTION:{sid}")

    record = json.loads(record_path.read_text())

    historical = (
        record.get("status") == "HISTORICAL_IMPORTED"
        and record.get("provenance_type") == "HISTORICAL_MSA_REUSE"
    )
    fresh = (
        record.get("status") == "COMPLETED"
        and record.get("provenance_type") not in
        ("HISTORICAL_MSA_REUSE", "HISTORICAL_IMPORTED")
    )
    if (
        not (historical or fresh)
        or record.get("structure_id") != sid
        or record.get("model_type") != "alphafold2_ptm"
    ):
        raise ValueError(f"INVALID_MSA_PREDICTION:{sid}")
    if historical:
        # Historical files are QC-checked, not certified as having been
        # generated from the freshly prepared A3M.
        if (
            record.get("gpu_prediction_performed") is not False
            or record.get("original_run_provenance_verified") is not False
            or record.get("historical_alignment_relation") not in
            ("EXACT_ALIGNED_RECORDS", "UNGAPPED_MULTISET_ONLY")
        ):
            raise ValueError(f"INVALID_HISTORICAL_ATTRIBUTION:{sid}")
        prepared_a3m = checked_file(
            record["prepared_a3m"], record["prepared_a3m_sha256"]
        )

    length = int(single["length_aa"])

    if (
        length < 1
        or int(record["query_length"]) != length
        or int(record["sequence_count"]) < 2
        or single["prediction_status"] != "PREDICTED"
    ):
        raise ValueError(f"INVALID_SEQUENCE_METADATA:{sid}")

    if (
        single["sequence_sha256"]
        != record["query_sha256"]
    ):
        raise ValueError(f"QUERY_SEQUENCE_MISMATCH:{sid}")

    a3m = checked_file(
        record["source_a3m"],
        record["a3m_sha256"],
    )

    query_id, query_length, query_sha = (
        query_sha_from_a3m(a3m)
    )

    if (
        query_id != sid
        or query_length != length
        or query_sha != single["sequence_sha256"]
    ):
        raise ValueError(f"A3M_QUERY_MISMATCH:{sid}")

    single_pdb = checked_file(
        single["pdb_file"],
        single["pdb_sha256"],
    )

    msa_pdb = checked_file(
        record["pdb_file"],
        record["pdb_sha256"],
    )

    single_score_file = checked_file(
        single["score_json"],
        single["score_json_sha256"],
    )

    msa_score_file = checked_file(
        record["score_json"],
        record["score_json_sha256"],
    )

    single_scores = json.loads(
        single_score_file.read_text()
    )

    msa_scores = json.loads(
        msa_score_file.read_text()
    )

    single_plddt = pdb_mean(single_pdb, length)
    msa_plddt = pdb_mean(msa_pdb, length)

    if (
        abs(
            single_plddt
            - numeric(
                single["mean_plddt"],
                "stage08_mean_plddt",
                0,
                100,
            )
        ) > 0.05
        or abs(
            msa_plddt
            - numeric(
                record["mean_plddt"],
                "stage11_mean_plddt",
                0,
                100,
            )
        ) > 0.05
    ):
        raise ValueError(f"PLDDT_RECORD_MISMATCH:{sid}")

    single_ptm = numeric(
        single_scores["ptm"], "single_ptm", 0, 1
    )

    msa_ptm = numeric(
        msa_scores["ptm"], "msa_ptm", 0, 1
    )

    if (
        abs(single_ptm - float(single["ptm"])) > 0.011
        or abs(msa_ptm - float(record["ptm"])) > 0.001
    ):
        raise ValueError(f"PTM_RECORD_MISMATCH:{sid}")

    row = {
        "structure_id": sid,
        "family_rep": single["family_rep"],
        "single_mean_plddt": f"{single_plddt:.2f}",
        "msa_mean_plddt": f"{msa_plddt:.2f}",
        "delta_mean_plddt": (
            f"{msa_plddt-single_plddt:+.2f}"
        ),
        "single_ptm": f"{single_ptm:.4f}",
        "msa_ptm": f"{msa_ptm:.4f}",
        "delta_ptm": f"{msa_ptm-single_ptm:+.4f}",
        "msa_sequences": str(record["sequence_count"]),
        "query_length": str(length),
        "model_type": record["model_type"],
        "single_pdb_sha256": sha(single_pdb),
        "msa_pdb_sha256": sha(msa_pdb),
    }

    provenance = {
        "structure_id": sid,
        "prediction_provenance_type": (
            "HISTORICAL_MSA_REUSE" if historical else "FRESH_COLABFOLD"
        ),
        "historical_alignment_relation": (
            record.get("historical_alignment_relation") if historical else None
        ),
        "prepared_a3m_sha256": (
            sha(prepared_a3m) if historical else None
        ),
        "stage11_record": str(record_path.resolve()),
        "stage11_record_sha256": sha(record_path),
        "a3m_sha256": sha(a3m),
        "single_score_sha256": sha(single_score_file),
        "msa_score_sha256": sha(msa_score_file),
    }

    return row, provenance


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--structures",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--prediction-index",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )
    args = ap.parse_args()

    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            f"Refusing existing output directory: {out}"
        )

    structures = indexed(
        table(
            args.structures,
            (
                "structure_id",
                "family_rep",
                "sequence_sha256",
                "length_aa",
                "prediction_status",
                "mean_plddt",
                "ptm",
                "pdb_file",
                "pdb_sha256",
                "score_json",
                "score_json_sha256",
            ),
        ),
        "structure_id",
    )

    prediction_rows = table(
        args.prediction_index,
        ("structure_id", "prediction_provenance"),
    )

    prediction_index = indexed(
        prediction_rows, "structure_id"
    )

    if not set(prediction_index).issubset(structures):
        raise ValueError("UNKNOWN_MSA_STRUCTURE_ID")

    results = []
    evidence = []

    for sid, item in prediction_index.items():
        record_path = Path(
            item["prediction_provenance"]
        ).resolve()

        row, source = collect(
            structures[sid],
            record_path,
        )

        results.append(row)
        evidence.append(source)

    # All inputs are verified before output creation.
    out.mkdir(parents=True, exist_ok=False)

    result_path = out / "msa_structure_qc.tsv"

    with result_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)

    metadata = {
        "algorithm": "VERIFIED_SINGLE_VS_FAMILY_MSA",
        "source_stage08_sha256": sha(args.structures),
        "prediction_index_sha256": sha(
            args.prediction_index
        ),
        "families": len(results),
        "msa_qc_sha256": sha(result_path),
        "evidence": evidence,
    }

    (
        out / "msa_structure_qc_provenance.json"
    ).write_text(
        json.dumps(
            metadata,
            indent=2,
            ensure_ascii=False,
        ) + "\n"
    )

    print(f"MSA_QC_FAMILIES={len(results)}")

    for row in results:
        print(
            "STRUCTURE_COMPARISON="
            f"{row['structure_id']} "
            f"pLDDT:{row['single_mean_plddt']}"
            f"->{row['msa_mean_plddt']} "
            f"delta:{row['delta_mean_plddt']} "
            f"pTM:{row['single_ptm']}"
            f"->{row['msa_ptm']}"
        )

    print(f"MSA_QC_TABLE={result_path}")
    print("STAGE11_MSA_STRUCTURE_QC=PASS")
    print("GPU_PREDICTION=NOT_RUN")


if __name__ == "__main__":
    main()
