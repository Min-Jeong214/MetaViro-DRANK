#!/usr/bin/env python3
"""Reject MSA reuse if the actual Stage 08 family differs from its source.

Historical predictions require their source member table and FASTA.
Newer prepared predictions may use hashed preparation provenance.
No prediction, search, ranking, or output directory is created.
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(path):
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_INPUT:{path}")
    return path


def tsv(path, fields):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)

    if (
        not rows
        or len(header) != len(set(header))
        or not set(fields).issubset(header)
        or any(
            None in row
            or any(value is None for value in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    return rows


def json_record(path):
    obj = json.loads(require(path).read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return obj


def proteins(path):
    records = {}
    name = None
    chunks = []

    def finish():
        if name is not None:
            sequence = "".join(chunks).upper()
            if not sequence or name in records:
                raise ValueError(f"INVALID_FASTA:{path}")
            records[name] = sequence

    with require(path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                finish()
                parts = line[1:].split()
                if not parts:
                    raise ValueError("EMPTY_FASTA_HEADER")
                name = parts[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError("SEQUENCE_WITHOUT_HEADER")
                chunks.append(line)

    finish()

    if not records:
        raise ValueError("EMPTY_FASTA")

    return records


def membership(members_file, fasta_file, representative):
    rows = tsv(
        members_file, ("member", "representative")
    )
    sequences = proteins(fasta_file)

    member_ids = [r["member"] for r in rows]

    if (
        not all(member_ids)
        or len(member_ids) != len(set(member_ids))
        or set(member_ids) != set(sequences)
    ):
        raise ValueError("MEMBER_TABLE_FASTA_MISMATCH")

    selected = {
        r["member"]: sequences[r["member"]]
        for r in rows
        if r["representative"] == representative
    }

    if not selected:
        raise ValueError(
            f"FAMILY_NOT_PRESENT:{representative}"
        )

    return selected


def a3m_sequences(path):
    records = proteins(path)
    sequences = []

    for seq in records.values():
        # A3M lowercase letters are insertions; '-' are gaps.
        normalized = "".join(
            character.upper()
            for character in seq
            if character not in "-."
        )
        if not normalized:
            raise ValueError("EMPTY_A3M_SEQUENCE")

        sequences.append(normalized)

    return records, sequences


def check(args):
    upstream = args.stage08_output.resolve()
    qc_dir = args.msa_qc_dir.resolve()

    member_file = require(
        upstream /
        "05_annotation/recurrent_members.tsv"
    )
    fasta_file = require(
        upstream /
        "05_annotation/recurrent_members.faa"
    )

    candidates = sorted(
        (
            upstream / "08_structure"
        ).glob(
            "run_*/collected/structure_results.tsv"
        )
    )

    if len(candidates) != 1:
        raise ValueError("STRUCTURE_TABLE_NOT_UNIQUE")

    structure_file = require(candidates[0])

    structures = tsv(
        structure_file,
        (
            "structure_id", "family_rep",
            "sequence_sha256", "pdb_sha256",
        ),
    )
    qc_rows = tsv(
        qc_dir / "msa_structure_qc.tsv",
        (
            "structure_id", "family_rep",
            "msa_pdb_sha256",
        ),
    )
    qc_prov = json_record(
        qc_dir / "msa_structure_qc_provenance.json"
    )

    if (
        qc_prov.get("source_stage08_sha256")
        != sha(structure_file)
        or qc_prov.get("msa_qc_sha256")
        != sha(qc_dir / "msa_structure_qc.tsv")
    ):
        raise ValueError("STAGE08_QC_PROVENANCE_MISMATCH")

    evidence = qc_prov.get("evidence")

    if not isinstance(evidence, list) or not evidence:
        raise ValueError("MSA_QC_EVIDENCE_MISSING")

    stage08_index = {
        r["structure_id"]: r
        for r in structures
    }
    qc_index = {
        r["structure_id"]: r
        for r in qc_rows
    }
    evidence_index = {
        r["structure_id"]: r
        for r in evidence
    }

    if (
        len(stage08_index) != len(structures)
        or len(qc_index) != len(qc_rows)
        or len(evidence_index) != len(evidence)
        or set(qc_index) != set(evidence_index)
    ):
        raise ValueError("QC_EVIDENCE_ID_MISMATCH")

    if args.stage09_10_output:
        shortlist = tsv(
            args.stage09_10_output.resolve() /
            "10_shortlists/top15_experimental_candidates.tsv",
            ("structure_id", "family_rep"),
        )
        selected = {
            row["structure_id"]: row["family_rep"]
            for row in shortlist
        }

        if (
            len(selected) != len(shortlist)
            or set(selected) != set(qc_index)
        ):
            raise ValueError(
                "SELECTED_MSA_QC_COVERAGE_MISMATCH"
            )

        for sid, row in qc_index.items():
            if selected[sid] != row["family_rep"]:
                raise ValueError(
                    f"SELECTED_FAMILY_MISMATCH:{sid}"
                )

    for sid, qc_row in qc_index.items():
        if sid not in stage08_index:
            raise ValueError(
                f"QC_NOT_IN_STAGE08:{sid}"
            )

        source = stage08_index[sid]
        detail = evidence_index[sid]

        if source["family_rep"] != qc_row["family_rep"]:
            raise ValueError(
                f"QC_FAMILY_REP_MISMATCH:{sid}"
            )

        record_file = require(
            detail["stage11_record"]
        )

        if (
            sha(record_file)
            != detail["stage11_record_sha256"]
        ):
            raise ValueError(
                f"MSA_RECORD_TAMPERED:{sid}"
            )

        record = json_record(record_file)

        if (
            record.get("structure_id") != sid
            or record.get("status") != "COMPLETED"
            or record.get("colabfold_returncode") != 0
            or record.get("query_sha256")
            != source["sequence_sha256"]
            or record.get("pdb_sha256")
            != qc_row["msa_pdb_sha256"]
            or sha(require(record["pdb_file"]))
            != record["pdb_sha256"]
        ):
            raise ValueError(
                f"MSA_PREDICTION_SOURCE_MISMATCH:{sid}"
            )

        source_a3m = require(
            record["source_a3m"]
        )

        if (
            sha(source_a3m)
            != record.get("a3m_sha256")
            or sha(source_a3m)
            != detail.get("a3m_sha256")
        ):
            raise ValueError(
                f"MSA_A3M_SHA256_MISMATCH:{sid}"
            )

        a3m_records, a3m_proteins = a3m_sequences(
            source_a3m
        )

        if (
            len(a3m_records) < 2
            or len(a3m_records)
            != int(record["sequence_count"])
        ):
            raise ValueError(
                f"MSA_SEQUENCE_COUNT_MISMATCH:{sid}"
            )

        msa_dir = source_a3m.parent.parent
        summary_path = require(
            msa_dir / "family_msa_summary.tsv"
        )

        matches = [
            r for r in tsv(
                summary_path,
                (
                    "structure_id", "family_rep",
                    "family_size",
                    "unique_sequences", "a3m_sha256",
                ),
            )
            if r["structure_id"] == sid
        ]

        if len(matches) != 1:
            raise ValueError(
                f"MSA_SUMMARY_NOT_UNIQUE:{sid}"
            )

        summary = matches[0]
        rep = source["family_rep"]
        current = membership(
            member_file, fasta_file, rep
        )

        if (
            summary["family_rep"] != rep
            or summary["a3m_sha256"] != sha(source_a3m)
            or len(current) != int(summary["family_size"])
            or len(set(current.values()))
            != int(summary["unique_sequences"])
            or set(current.values()) != set(a3m_proteins)
            or source["sequence_sha256"]
            != hashlib.sha256(
                current[rep].encode()
            ).hexdigest()
            if rep in current
            else True
        ):
            raise ValueError(
                f"MSA_REUSE_FAMILY_MEMBERSHIP_MISMATCH:{sid}"
            )

        # Verify the precise source membership whenever historical
        # source files are available; this includes duplicate members.
        historical_members = (
            msa_dir.parent / "members.tsv"
        )
        historical_fasta = (
            msa_dir.parent / "members.faa"
        )

        prepared_prov = (
            msa_dir.parent /
            "preparation_provenance.json"
        )

        if (
            historical_members.is_file()
            and historical_fasta.is_file()
        ):
            prior = membership(
                historical_members,
                historical_fasta,
                rep,
            )

            if prior != current:
                raise ValueError(
                    f"MSA_SOURCE_MEMBER_IDS_OR_SEQUENCES_MISMATCH:{sid}"
                )

            source_type = "HISTORICAL_MEMBER_TABLE_AND_FASTA"

        elif prepared_prov.is_file():
            prepared = json_record(prepared_prov)

            if (
                prepared.get("family_members_sha256")
                != sha(member_file)
                or prepared.get("member_fasta_sha256")
                != sha(fasta_file)
                or prepared.get("colabfold_inputs", {})
                    .get(sid, {})
                    .get("sha256")
                != sha(source_a3m)
            ):
                raise ValueError(
                    f"MSA_PREPARED_MEMBERSHIP_SHA256_MISMATCH:{sid}"
                )

            source_type = "PREPARATION_PROVENANCE"

        else:
            raise ValueError(
                f"MSA_SOURCE_MEMBERSHIP_PROVENANCE_MISSING:{sid}"
            )

        print(
            f"MSA_MEMBERSHIP_VERIFIED={sid} "
            f"members={len(current)} "
            f"unique={len(set(current.values()))} "
            f"source={source_type}"
        )

    print(
        f"MSA_REUSE_VERIFIED_FAMILIES={len(qc_index)}"
    )
    print("MSA_FAMILY_MEMBERSHIP_REUSE_GUARD=PASS")
    print("OUTPUT_CREATED=NO")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--stage08-output", type=Path, required=True
    )
    ap.add_argument(
        "--msa-qc-dir", type=Path, required=True
    )
    ap.add_argument(
        "--stage09-10-output", type=Path
    )

    check(ap.parse_args())


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, KeyError,
        TypeError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
