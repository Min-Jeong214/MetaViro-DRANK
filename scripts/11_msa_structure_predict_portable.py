#!/usr/bin/env python3
"""Non-destructive ColabFold prediction from a prepared family A3M."""

import argparse
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def validate_a3m(path):
    if not path.is_file() or path.suffix.lower() != ".a3m":
        raise ValueError("A3M_INPUT_MISSING_OR_INVALID")

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
                raise ValueError("A3M_MISSING_HEADER")
            chunks.append(line)

    if name is not None:
        records.append((name, "".join(chunks)))

    if len(records) < 2:
        raise ValueError("FAMILY_MSA_REQUIRES_MULTIPLE_SEQUENCES")

    ids = [name for name, _ in records]

    if len(set(ids)) != len(ids):
        raise ValueError("DUPLICATE_A3M_IDS")

    query_id, query = records[0]

    if query_id != path.stem:
        raise ValueError("QUERY_ID_FILENAME_MISMATCH")

    if not re.fullmatch(r"[A-Z]+", query):
        raise ValueError("INVALID_UNGAPPED_QUERY")

    for member_id, sequence in records:
        if not re.fullmatch(r"[A-Za-z-]+", sequence):
            raise ValueError(f"INVALID_A3M_SEQUENCE:{member_id}")

        columns = sum(
            char == "-" or char.isupper()
            for char in sequence
        )

        if columns != len(query):
            raise ValueError(
                f"A3M_COLUMN_MISMATCH:{member_id}"
            )

        if not any(char.isalpha() for char in sequence):
            raise ValueError(f"EMPTY_MEMBER_SEQUENCE:{member_id}")

    return {
        "structure_id": query_id,
        "query_length": len(query),
        "sequence_count": len(records),
        "query_sha256": hashlib.sha256(
            query.encode()
        ).hexdigest(),
        "a3m_sha256": digest(path),
    }


def get_version():
    try:
        return version("colabfold")
    except PackageNotFoundError:
        return "UNKNOWN"


def get_command(binary, input_path, results_path, args):
    return [
        binary,
        "--model-type", "alphafold2_ptm",
        "--num-models", "1",
        "--num-recycle", "1",
        "--num-seeds", "1",
        "--random-seed", "0",
        "--rank", "plddt",
        "--max-seq", str(args.max_seq),
        "--max-extra-seq", str(args.max_extra_seq),
        "--num-relax", "0",
        str(input_path),
        str(results_path),
    ]


def validate_prediction(folder, query_length):
    pdb_files = sorted(
        folder.glob("*unrelaxed_rank_001_*.pdb")
    )
    score_files = sorted(
        folder.glob("*scores_rank_001_*.json")
    )

    if len(pdb_files) != 1 or len(score_files) != 1:
        raise ValueError(
            "EXPECTED_ONE_RANK1_PDB_AND_SCORE_JSON"
        )

    pdb = pdb_files[0]
    scores_file = score_files[0]

    confidences = []
    residues = set()

    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[12:16].strip() != "CA":
            continue
        if line[16] not in (" ", "A"):
            continue

        residue = (line[21], line[22:27])

        if residue in residues:
            raise ValueError("DUPLICATE_PDB_CA")

        residues.add(residue)
        confidences.append(float(line[60:66]))

    if len(confidences) != query_length:
        raise ValueError("PDB_QUERY_LENGTH_MISMATCH")

    scores = json.loads(scores_file.read_text())

    if "ptm" not in scores:
        raise ValueError("MISSING_PTM_SCORE")

    return {
        "pdb_file": str(pdb.resolve()),
        "pdb_sha256": digest(pdb),
        "score_json": str(scores_file.resolve()),
        "score_json_sha256": digest(scores_file),
        "mean_plddt": round(
            sum(confidences) / len(confidences), 2
        ),
        "ptm": round(float(scores["ptm"]), 4),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--a3m", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--max-seq", type=int, default=512)
    parser.add_argument("--max-extra-seq", type=int, default=1024)

    args = parser.parse_args()

    source = args.a3m.resolve()
    output = args.out_dir.resolve()

    if args.max_seq < 2 or args.max_extra_seq < 1:
        raise ValueError("INVALID_MSA_CAPACITY")

    if output.exists():
        raise ValueError(
            f"Refusing existing output directory: {output}"
        )

    if output == source or output in source.parents:
        raise ValueError("OUTPUT_OVERLAPS_MSA_INPUT")

    metadata = validate_a3m(source)

    binary = shutil.which("colabfold_batch")

    if binary is None:
        raise RuntimeError(
            "COLABFOLD_NOT_FOUND_IN_ACTIVE_ENVIRONMENT"
        )

    input_copy = output / "input" / source.name
    predictions = output / "predictions"

    command = get_command(
        binary, input_copy, predictions, args
    )

    metadata.update({
        "source_a3m": str(source),
        "colabfold_version": get_version(),
        "model_type": "alphafold2_ptm",
        "num_models": 1,
        "num_recycle": 1,
        "num_seeds": 1,
        "random_seed": 0,
        "ranking_metric": "plddt",
        "max_seq": args.max_seq,
        "max_extra_seq": args.max_extra_seq,
        "num_relax": 0,
        "external_msa_search": False,
        "command": command,
    })

    print(f"STRUCTURE_ID={metadata['structure_id']}")
    print(f"QUERY_LENGTH={metadata['query_length']}")
    print(f"MSA_SEQUENCE_COUNT={metadata['sequence_count']}")
    print(f"COLABFOLD_VERSION={metadata['colabfold_version']}")
    print(f"A3M_SHA256={metadata['a3m_sha256']}")
    print("INPUT_A3M_VALIDATION=PASS")
    print("COLABFOLD_COMMAND=" + shlex.join(command))

    if not args.execute:
        print("GPU_PREDICTION=NOT_RUN")
        print("OUTPUT_DIRECTORY=NOT_CREATED")
        print("STAGE11_MSA_PREDICTION_PREFLIGHT=PASS")
        return

    output.mkdir(parents=True, exist_ok=False)
    input_copy.parent.mkdir()
    predictions.mkdir()

    shutil.copy2(source, input_copy)

    if digest(input_copy) != metadata["a3m_sha256"]:
        raise ValueError("A3M_COPY_INTEGRITY_FAILURE")

    record_path = output / "prediction_provenance.json"
    metadata["status"] = "RUNNING"

    record_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    log_path = output / "colabfold.log"

    print("FAMILY_MSA_GPU_PREDICTION=STARTED", flush=True)

    with log_path.open("w") as log:
        with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        ) as process:
            for line in process.stdout:
                print(line, end="", flush=True)
                log.write(line)

            returncode = process.wait()

    metadata["colabfold_returncode"] = returncode
    metadata["colabfold_log_sha256"] = digest(log_path)

    if returncode != 0:
        metadata["status"] = "FAILED"
        record_path.write_text(
            json.dumps(metadata, indent=2) + "\n"
        )
        raise RuntimeError(
            f"COLABFOLD_FAILED_RETURN_CODE_{returncode}"
        )

    if digest(source) != metadata["a3m_sha256"]:
        raise ValueError("SOURCE_A3M_CHANGED_DURING_RUN")

    prediction = validate_prediction(
        predictions, metadata["query_length"]
    )

    metadata.update(prediction)
    metadata["status"] = "COMPLETED"

    record_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    print(f"MSA_MEAN_PLDDT={metadata['mean_plddt']}")
    print(f"MSA_PTM={metadata['ptm']}")
    print(f"MSA_PDB={metadata['pdb_file']}")
    print("MSA_PDB_INTEGRITY=PASS")
    print("STAGE11_MSA_GPU_PREDICTION=PASS")


if __name__ == "__main__":
    main()
