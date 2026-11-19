#!/usr/bin/env python3
"""Convert verified Family-MSA PDBs to portable Foldseek input."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path, required):
    if not path.is_file():
        raise ValueError(f"MISSING_INPUT:{path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = reader.fieldnames or []
        rows = list(reader)

    if (
        len(columns) != len(set(columns))
        or not set(required).issubset(columns)
        or not rows
    ):
        raise ValueError(f"INVALID_INPUT_SCHEMA:{path}")

    return columns, rows


def index(rows):
    result = {}
    for row in rows:
        sid = row["structure_id"]
        if not sid or sid in result:
            raise ValueError(f"DUPLICATE_STRUCTURE_ID:{sid}")
        result[sid] = row
    return result


def verified(path, expected):
    path = Path(path).resolve()
    if not path.is_file() or sha(path) != expected:
        raise ValueError(f"FILE_INTEGRITY_FAILURE:{path}")
    return path


def query_a3m(path):
    header = None
    chunks = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].split()[0]
        elif header is not None:
            chunks.append(line)
        else:
            raise ValueError("INVALID_A3M")

    sequence = "".join(chunks)
    if not sequence or not sequence.isalpha() or not sequence.isupper():
        raise ValueError("INVALID_A3M_QUERY")

    return header, sequence


def pdb_confidence(path, length):
    values = []
    residues = set()

    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[12:16].strip() != "CA":
            continue
        if line[16] not in (" ", "A"):
            continue

        residue = (line[21], line[22:27])

        if residue in residues:
            raise ValueError("DUPLICATE_PDB_RESIDUE")

        residues.add(residue)
        score = float(line[60:66])

        if not math.isfinite(score) or not 0 <= score <= 100:
            raise ValueError("INVALID_PDB_PLDDT")

        values.append(score)

    if len(values) != length:
        raise ValueError("MSA_PDB_LENGTH_MISMATCH")

    return min(values), max(values), sum(values) / length


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", type=Path, required=True)
    ap.add_argument("--prediction-index", type=Path, required=True)
    ap.add_argument("--msa-qc", type=Path, required=True)
    ap.add_argument("--msa-qc-provenance", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError("Refusing existing output directory")

    stage08_fields, stage08_rows = read_tsv(
        args.structures,
        (
            "structure_id", "sequence_sha256",
            "length_aa", "prediction_status",
            "pdb_file", "pdb_sha256",
        ),
    )

    _, prediction_rows = read_tsv(
        args.prediction_index,
        ("structure_id", "prediction_provenance"),
    )

    _, qc_rows = read_tsv(
        args.msa_qc,
        (
            "structure_id", "single_pdb_sha256",
            "msa_pdb_sha256", "msa_mean_plddt",
            "msa_ptm", "msa_sequences",
        ),
    )

    qc_meta = json.loads(
        args.msa_qc_provenance.read_text()
    )

    if (
        qc_meta.get("algorithm")
        != "VERIFIED_SINGLE_VS_FAMILY_MSA"
        or qc_meta.get("source_stage08_sha256")
        != sha(args.structures)
        or qc_meta.get("prediction_index_sha256")
        != sha(args.prediction_index)
        or qc_meta.get("msa_qc_sha256")
        != sha(args.msa_qc)
    ):
        raise ValueError("MSA_QC_PROVENANCE_MISMATCH")

    originals = index(stage08_rows)
    predictions = index(prediction_rows)
    quality = index(qc_rows)

    if (
        set(predictions) != set(quality)
        or not set(predictions).issubset(originals)
        or qc_meta.get("families") != len(predictions)
    ):
        raise ValueError("FAMILY_SET_MISMATCH")

    evidence = {
        item["structure_id"]: item
        for item in qc_meta.get("evidence", [])
    }

    if set(evidence) != set(predictions):
        raise ValueError("QC_EVIDENCE_SET_MISMATCH")

    converted = []
    provenance = []

    for sid, item in predictions.items():
        original = originals[sid]
        qc = quality[sid]

        record_path = Path(
            item["prediction_provenance"]
        ).resolve()

        verified(
            record_path,
            evidence[sid]["stage11_record_sha256"],
        )

        record = json.loads(record_path.read_text())

        if (
            record.get("status") != "COMPLETED"
            or record.get("structure_id") != sid
            or record.get("model_type") != "alphafold2_ptm"
            or original["prediction_status"] != "PREDICTED"
        ):
            raise ValueError(f"INVALID_PREDICTION:{sid}")

        if (
            qc["single_pdb_sha256"]
            != original["pdb_sha256"]
        ):
            raise ValueError(f"SINGLE_PDB_MISMATCH:{sid}")

        verified(
            original["pdb_file"],
            original["pdb_sha256"],
        )

        length = int(original["length_aa"])

        if (
            length < 1
            or length != int(record["query_length"])
            or int(record["sequence_count"]) < 2
            or int(qc["msa_sequences"])
            != int(record["sequence_count"])
            or original["sequence_sha256"]
            != record["query_sha256"]
        ):
            raise ValueError(f"MSA_SEQUENCE_MISMATCH:{sid}")

        a3m = verified(
            record["source_a3m"],
            record["a3m_sha256"],
        )

        query_id, query = query_a3m(a3m)

        if (
            query_id != sid
            or len(query) != length
            or hashlib.sha256(
                query.encode()
            ).hexdigest() != original["sequence_sha256"]
        ):
            raise ValueError(f"MSA_QUERY_MISMATCH:{sid}")

        pdb = verified(
            record["pdb_file"],
            record["pdb_sha256"],
        )

        score_file = verified(
            record["score_json"],
            record["score_json_sha256"],
        )

        if qc["msa_pdb_sha256"] != sha(pdb):
            raise ValueError(f"MSA_QC_PDB_MISMATCH:{sid}")

        scores = json.loads(score_file.read_text())
        ptm = float(scores["ptm"])

        if not math.isfinite(ptm) or not 0 <= ptm <= 1:
            raise ValueError(f"INVALID_PTM:{sid}")

        minimum, maximum, mean = pdb_confidence(
            pdb, length
        )

        if (
            abs(mean - float(qc["msa_mean_plddt"])) > 0.05
            or abs(mean - float(record["mean_plddt"])) > 0.05
            or abs(ptm - float(qc["msa_ptm"])) > 0.001
            or abs(ptm - float(record["ptm"])) > 0.001
        ):
            raise ValueError(f"MSA_QC_METRIC_MISMATCH:{sid}")

        row = {
            field: original[field]
            for field in stage08_fields
        }

        row.update({
            "pdb_file": str(pdb),
            "pdb_sha256": sha(pdb),
            "score_json": str(score_file),
            "score_json_sha256": sha(score_file),
            "mean_plddt": f"{mean:.2f}",
            "min_plddt": f"{minimum:.2f}",
            "max_plddt": f"{maximum:.2f}",
            "ptm": f"{ptm:.4f}",
            "confidence_class": "",
            "prediction_status": "PREDICTED",
        })

        converted.append(row)

        provenance.append({
            "structure_id": sid,
            "source": "FAMILY_MSA_PREDICTION",
            "prediction_record_sha256": sha(record_path),
            "a3m_sha256": sha(a3m),
            "msa_pdb_sha256": sha(pdb),
            "msa_score_sha256": sha(score_file),
        })

    # No output is created until all Families pass validation.
    output.mkdir(parents=True, exist_ok=False)

    result = output / "msa_structure_results.tsv"

    with result.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=stage08_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(converted)

    metadata = {
        "query_origin": "STAGE11_FAMILY_MSA",
        "purpose": "STAGE_K2_EXTERNAL_STRUCTURE_SEARCH",
        "external_evidence": "NOT_SEARCHED",
        "stage08_sha256": sha(args.structures),
        "msa_qc_sha256": sha(args.msa_qc),
        "msa_qc_provenance_sha256": sha(
            args.msa_qc_provenance
        ),
        "structure_count": len(converted),
        "structure_table_sha256": sha(result),
        "structures": provenance,
    }

    (
        output / "msa_foldseek_input_provenance.json"
    ).write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    print(f"MSA_STRUCTURE_COUNT={len(converted)}")
    print(f"FOLDSEEK_INPUT_COLUMNS={len(stage08_fields)}")
    print(f"FOLDSEEK_INPUT={result}")
    print("STAGE_K2_MSA_QUERY_CONVERSION=PASS")
    print("EXTERNAL_EVIDENCE=NOT_SEARCHED")


if __name__ == "__main__":
    main()
