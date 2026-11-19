#!/usr/bin/env python3
"""Prepare selected family MSAs and ColabFold A3M inputs.

Uses all constituent proteins of each selected family.
Never performs GPU prediction or final candidate selection.
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MSA_SCRIPT = ROOT / "scripts/11_family_msa_portable.py"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_OR_EMPTY_INPUT:{path}")
    return path


def table(path, required):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not rows
        or len(fields) != len(set(fields))
        or not set(required).issubset(fields)
        or any(
            None in row
            or any(value is None for value in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    return fields, rows


def fasta_ids(path):
    ids = []

    with require(path).open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                ids.append(identifier)

    if not ids or len(ids) != len(set(ids)):
        raise ValueError("INVALID_MEMBER_FASTA_IDS")

    return set(ids)


def validate(args):
    up = args.stage08_output.resolve()
    down = args.stage09_10_output.resolve()
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    for source in (up, down):
        if not source.is_dir():
            raise ValueError(f"MISSING_INPUT_DIRECTORY:{source}")

        if (
            out == source
            or out in source.parents
            or source in out.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    if args.threads < 1:
        raise ValueError("INVALID_THREADS")

    members_file = require(
        up / "05_annotation/recurrent_members.tsv"
    )
    fasta_file = require(
        up / "05_annotation/recurrent_members.faa"
    )

    shortlist = require(
        down /
        "10_shortlists/top15_experimental_candidates.tsv"
    )

    ranking = require(
        down /
        "10_ranking/confidence_aware_ranking.tsv"
    )

    ranking_prov_file = require(
        down / "10_ranking/ranking_provenance.json"
    )

    manifest = json.loads(
        require(
            down / "orchestration_provenance.json"
        ).read_text()
    )

    rank_meta = json.loads(
        ranking_prov_file.read_text()
    )

    artifacts = manifest.get("artifacts", {})

    for key, path in (
        ("top15", shortlist),
        ("ranking", ranking),
    ):
        item = artifacts.get(key, {})

        if (
            item.get("sha256") != sha(path)
            or Path(item.get("path", "/")).resolve()
            != path
        ):
            raise ValueError(
                f"STAGE10_ARTIFACT_MISMATCH:{key}"
            )

    if sha(ranking) not in [
        value for key, value in rank_meta.items()
        if key.endswith("sha256")
        and isinstance(value, str)
    ]:
        raise ValueError("RANKING_PROVENANCE_MISMATCH")

    structure_candidates = sorted(
        (
            up / "08_structure"
        ).glob(
            "run_*/collected/structure_results.tsv"
        )
    )

    structures = [
        path for path in structure_candidates
        if sha(path) == manifest.get(
            "stage08_structure_table_sha256"
        )
    ]

    if len(structures) != 1:
        raise ValueError(
            "STAGE08_STRUCTURE_SOURCE_NOT_UNIQUE"
        )

    structures_file = structures[0]
    mapping_file = require(
        structures_file.parent.parent /
        "prepared/structure_id_map.tsv"
    )

    map_fields, map_rows = table(
        mapping_file,
        (
            "structure_id", "family_rep",
            "family_size", "length_aa",
            "sequence_sha256",
        ),
    )

    _, structure_rows = table(
        structures_file,
        (
            "structure_id", "family_rep",
            "sequence_sha256",
        ),
    )

    _, members = table(
        members_file,
        ("member", "representative"),
    )

    _, top15 = table(
        shortlist,
        ("structure_id", "family_rep"),
    )

    _, ranks = table(
        ranking,
        ("structure_id", "family_rep"),
    )

    rank_index = {
        row["structure_id"]: row
        for row in ranks
    }

    map_index = {
        row["structure_id"]: row
        for row in map_rows
    }

    structure_index = {
        row["structure_id"]: row
        for row in structure_rows
    }

    if (
        len(rank_index) != len(ranks)
        or len(map_index) != len(map_rows)
        or len(structure_index) != len(structure_rows)
    ):
        raise ValueError("DUPLICATE_STRUCTURE_ID")

    members_in_fasta = fasta_ids(fasta_file)
    member_ids = [row["member"] for row in members]

    if (
        len(member_ids) != len(set(member_ids))
        or set(member_ids) != members_in_fasta
    ):
        raise ValueError("MEMBER_FASTA_TABLE_MISMATCH")

    counts = Counter(
        row["representative"] for row in members
    )

    selected = []
    seen = set()

    for row in top15:
        sid = row["structure_id"]

        if (
            sid in seen
            or sid not in map_index
            or sid not in rank_index
            or sid not in structure_index
        ):
            raise ValueError(
                f"INVALID_SELECTED_STRUCTURE:{sid}"
            )

        if rank_index[sid] != row:
            raise ValueError(
                f"SHORTLIST_RANKING_MISMATCH:{sid}"
            )

        mapped = map_index[sid]
        structure = structure_index[sid]
        rep = mapped["family_rep"]

        if (
            rep != row["family_rep"]
            or rep != structure["family_rep"]
            or mapped["sequence_sha256"]
            != structure["sequence_sha256"]
            or int(mapped["family_size"])
            != counts[rep]
        ):
            raise ValueError(
                f"FAMILY_OR_SEQUENCE_MISMATCH:{sid}"
            )

        if not re.fullmatch(
            r"[A-Za-z0-9_.-]+", sid
        ):
            raise ValueError("UNSAFE_STRUCTURE_ID")

        seen.add(sid)
        selected.append(mapped)

    print(
        "SELECTED_FAMILIES=",
        len(selected),
    )
    print(
        "SELECTED_FAMILY_MEMBER_COUNTS=",
        {
            row["structure_id"]:
                counts[row["family_rep"]]
            for row in selected
        },
    )
    print("STAGE10_SELECTION_PROVENANCE=PASS")
    print("STAGE08_STRUCTURE_MAP=PASS")
    print("FAMILY_MEMBER_FASTA_CONSISTENCY=PASS")

    return {
        "out": out,
        "map_fields": map_fields,
        "selected": selected,
        "members_file": members_file,
        "fasta_file": fasta_file,
        "mapping_file": mapping_file,
        "structures_file": structures_file,
        "shortlist": shortlist,
        "ranking": ranking,
    }


def execute(args, state):
    out = state["out"]
    out.mkdir(parents=True, exist_ok=False)

    selected_map = (
        out / "selected_structure_map.tsv"
    )

    with selected_map.open(
        "x", newline="", encoding="utf-8"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=state["map_fields"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(state["selected"])

    msa_dir = out / "msa"
    log = out / "stage11_msa.log"

    command = [
        sys.executable, "-B", str(MSA_SCRIPT),
        "--members", str(state["members_file"]),
        "--fasta", str(state["fasta_file"]),
        "--structure-map", str(selected_map),
        "--out-dir", str(msa_dir),
        "--threads", str(args.threads),
    ]

    print("[RUN] 11 selected family MSA", flush=True)

    with log.open("x", encoding="utf-8") as fh:
        result = subprocess.run(
            command,
            cwd=ROOT,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=False,
        )

    if result.returncode:
        raise RuntimeError(
            f"STAGE11_MSA_FAILED:SEE_LOG:{log}"
        )

    print("[PASS] 11 selected family MSA", flush=True)

    summary_file = require(
        msa_dir / "family_msa_summary.tsv"
    )

    _, summary = table(
        summary_file,
        (
            "structure_id", "family_rep",
            "family_size", "unique_sequences",
            "query_length_aa", "a3m_sha256",
        ),
    )

    summary_index = {
        row["structure_id"]: row
        for row in summary
    }

    expected = {
        row["structure_id"]: row
        for row in state["selected"]
    }

    if (
        len(summary_index) != len(summary)
        or set(summary_index) != set(expected)
    ):
        raise ValueError("MSA_OUTPUT_FAMILY_MISMATCH")

    a3m_inputs = {}
    eligibility_rows = []

    for sid, record in summary_index.items():
        mapping = expected[sid]

        a3m = require(
            msa_dir / "a3m_colabfold" /
            f"{sid}.a3m"
        )

        require(
            msa_dir / "a3m" / f"{sid}.a3m"
        )

        if (
            record["family_rep"]
            != mapping["family_rep"]
            or int(record["family_size"])
            != int(mapping["family_size"])
            or int(record["query_length_aa"])
            != int(mapping["length_aa"])
            or not (
                1 <= int(
                    record["unique_sequences"]
                ) <= int(record["family_size"])
            )
            or record["a3m_sha256"] != sha(a3m)
        ):
            raise ValueError(
                f"MSA_PROVENANCE_MISMATCH:{sid}"
            )

        a3m_inputs[sid] = {
            "path": str(a3m.resolve()),
            "sha256": sha(a3m),
            "family_members":
                int(mapping["family_size"]),
            "unique_sequences":
                int(record["unique_sequences"]),
        }

        unique_n = int(record["unique_sequences"])

        status = (
            "MSA_PREDICTION_ELIGIBLE"
            if unique_n >= 2
            else "SINGLE_UNIQUE_SEQUENCE_NO_MSA_REPREDICTION"
        )

        a3m_inputs[sid]["msa_prediction_eligible"] = (
            unique_n >= 2
        )
        a3m_inputs[sid]["msa_prediction_status"] = status

        eligibility_rows.append({
            "structure_id": sid,
            "family_rep": mapping["family_rep"],
            "family_members": mapping["family_size"],
            "unique_sequences": str(unique_n),
            "msa_prediction_status": status,
        })

    eligible_count = sum(
        row["msa_prediction_status"]
        == "MSA_PREDICTION_ELIGIBLE"
        for row in eligibility_rows
    )

    eligibility_file = out / "prediction_eligibility.tsv"

    with eligibility_file.open(
        "x", newline="", encoding="utf-8"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "structure_id",
                "family_rep",
                "family_members",
                "unique_sequences",
                "msa_prediction_status",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(eligibility_rows)

    provenance = {
        "stage": "STAGE11_PREPARE_SELECTED_FAMILIES",
        "mode": "MSA_AND_A3M_PREPARE_ONLY",
        "stage10_shortlist_sha256":
            sha(state["shortlist"]),
        "stage10_ranking_sha256":
            sha(state["ranking"]),
        "stage08_structures_sha256":
            sha(state["structures_file"]),
        "family_members_sha256":
            sha(state["members_file"]),
        "member_fasta_sha256":
            sha(state["fasta_file"]),
        "original_structure_map_sha256":
            sha(state["mapping_file"]),
        "selected_structure_map_sha256":
            sha(selected_map),
        "msa_summary_sha256":
            sha(summary_file),
        "families": len(summary),
        "msa_prediction_eligible_families": eligible_count,
        "single_unique_sequence_families":
            len(summary) - eligible_count,
        "prediction_eligibility_sha256":
            sha(eligibility_file),
        "colabfold_inputs": a3m_inputs,
        "gpu_prediction_performed": False,
        "old_prediction_reuse_verified": False,
        "final_candidate_selection_performed": False,
    }

    with (
        out / "preparation_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(
            provenance,
            fh,
            indent=2,
        )
        fh.write("\n")

    print(f"MSA_GENERATED_FAMILIES={len(summary)}")
    print(
        f"MSA_PREDICTION_ELIGIBLE_FAMILIES="
        f"{eligible_count}"
    )
    print(
        f"SINGLE_UNIQUE_SEQUENCE_FAMILIES="
        f"{len(summary) - eligible_count}"
    )
    print(f"COLABFOLD_A3M_DIRECTORY={msa_dir / 'a3m_colabfold'}")
    print("GPU_PREDICTION=NOT_RUN")
    print("OLD_PREDICTION_REUSE=NOT_ASSUMED")
    print("STAGE11_PREPARE_ORCHESTRATION=PASS")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage08-output",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--stage09-10-output",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--mode",
        choices=("check", "run"),
        default="check",
    )

    args = parser.parse_args()
    state = validate(args)

    if args.mode == "check":
        print("CHECK_MODE=PASS")
        print("OUTPUT_CREATED=NO")
        return

    execute(args, state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
