#!/usr/bin/env python3
"""Collect eligible MSA predictions or retain Stage 08-only evidence.

Never substitutes a single-sequence structure for an MSA prediction.
Only verified predictions proceed to Stage 11 MSA structure QC.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
QC_SCRIPT = ROOT / "scripts/11_msa_structure_qc_portable.py"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(path):
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"MISSING_OR_EMPTY_INPUT:{path}")
    return path


def metadata(path):
    obj = json.loads(require(path).read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"INVALID_METADATA:{path}")
    return obj


def table(path, required):
    with require(path).open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not rows
        or len(fields) != len(set(fields))
        or not set(required).issubset(fields)
        or any(
            None in row
            or any(v is None for v in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    if "structure_id" in fields:
        ids = [r["structure_id"] for r in rows]
        if (
            any(not sid for sid in ids)
            or len(ids) != len(set(ids))
        ):
            raise ValueError(f"INVALID_IDS:{path}")

    return rows


def verify(args):
    up = args.stage08_output.resolve()
    down = args.stage09_10_output.resolve()
    prep = args.prepared_dir.resolve()
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError("Refusing existing output directory")

    for source in (up, down, prep):
        if not source.is_dir():
            raise ValueError(f"INPUT_DIRECTORY_MISSING:{source}")
        if (
            out == source
            or source in out.parents
            or out in source.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    if args.prediction_index:
        idx_path = args.prediction_index.resolve()
        if (
            out == idx_path
            or idx_path in out.parents
            or out in idx_path.parents
        ):
            raise ValueError("PREDICTION_INDEX_OUTPUT_OVERLAP")

    prep_meta = metadata(
        prep / "preparation_provenance.json"
    )

    selected_path = require(
        prep / "selected_structure_map.tsv"
    )
    summary_path = require(
        prep / "msa/family_msa_summary.tsv"
    )
    policy_path = require(
        prep / "prediction_eligibility.tsv"
    )

    expected_hashes = (
        ("selected_structure_map_sha256", selected_path),
        ("msa_summary_sha256", summary_path),
        ("prediction_eligibility_sha256", policy_path),
        (
            "stage10_shortlist_sha256",
            down /
            "10_shortlists/top15_experimental_candidates.tsv",
        ),
        (
            "stage10_ranking_sha256",
            down /
            "10_ranking/confidence_aware_ranking.tsv",
        ),
        (
            "family_members_sha256",
            up / "05_annotation/recurrent_members.tsv",
        ),
        (
            "member_fasta_sha256",
            up / "05_annotation/recurrent_members.faa",
        ),
    )

    for key, path in expected_hashes:
        if prep_meta.get(key) != sha(require(path)):
            raise ValueError(f"PREPARATION_SHA256_MISMATCH:{key}")

    structures = [
        p for p in (
            up / "08_structure"
        ).glob("run_*/collected/structure_results.tsv")
        if sha(p) == prep_meta.get(
            "stage08_structures_sha256"
        )
    ]

    if len(structures) != 1:
        raise ValueError("STAGE08_STRUCTURE_SOURCE_NOT_UNIQUE")

    structures_path = structures[0]

    mapping = table(
        selected_path,
        (
            "structure_id", "family_rep",
            "family_size", "sequence_sha256",
        ),
    )
    summary = table(
        summary_path,
        (
            "structure_id", "family_rep",
            "family_size", "unique_sequences",
            "a3m_sha256",
        ),
    )
    policy = table(
        policy_path,
        (
            "structure_id", "family_rep",
            "family_members", "unique_sequences",
            "msa_prediction_status",
        ),
    )
    single = table(
        structures_path,
        (
            "structure_id", "family_rep",
            "sequence_sha256", "pdb_file",
            "pdb_sha256",
        ),
    )
    shortlist = table(
        down /
        "10_shortlists/top15_experimental_candidates.tsv",
        ("structure_id", "family_rep"),
    )

    smap = {r["structure_id"]: r for r in summary}
    pmap = {r["structure_id"]: r for r in policy}
    singlemap = {r["structure_id"]: r for r in single}
    shortmap = {r["structure_id"]: r for r in shortlist}

    ids = [r["structure_id"] for r in mapping]

    if not (
        set(ids) == set(smap) == set(pmap) == set(shortmap)
    ):
        raise ValueError("SELECTED_FAMILY_ID_MISMATCH")

    inputs = prep_meta.get("colabfold_inputs", {})
    if set(inputs) != set(ids):
        raise ValueError("PREPARED_A3M_ID_MISMATCH")

    eligible = []
    status_rows = []

    for row in mapping:
        sid = row["structure_id"]

        if sid not in singlemap:
            raise ValueError(f"STAGE08_ID_MISSING:{sid}")

        s = smap[sid]
        p = pmap[sid]
        original = singlemap[sid]
        info = inputs[sid]

        n = int(s["unique_sequences"])
        members = int(s["family_size"])

        expected_status = (
            "MSA_PREDICTION_ELIGIBLE"
            if n >= 2
            else "SINGLE_UNIQUE_SEQUENCE_NO_MSA_REPREDICTION"
        )

        a3m = require(
            prep / "msa/a3m_colabfold" / f"{sid}.a3m"
        )
        pdb = require(Path(original["pdb_file"]))

        if not (
            1 <= n <= members
            and members == int(row["family_size"])
            and p["family_members"] == str(members)
            and p["unique_sequences"] == str(n)
            and p["msa_prediction_status"] == expected_status
            and info.get("msa_prediction_status") == expected_status
            and info.get("msa_prediction_eligible") is (n >= 2)
            and info.get("sha256") == sha(a3m)
            and s["a3m_sha256"] == sha(a3m)
            and row["sequence_sha256"]
                == original["sequence_sha256"]
            and row["family_rep"]
                == original["family_rep"]
                == s["family_rep"]
                == p["family_rep"]
                == shortmap[sid]["family_rep"]
            and original["pdb_sha256"] == sha(pdb)
        ):
            raise ValueError(
                f"STAGE11_PREPARATION_OR_STRUCTURE_MISMATCH:{sid}"
            )

        if n >= 2:
            eligible.append(sid)

        status_rows.append({
            "structure_id": sid,
            "family_rep": row["family_rep"],
            "family_members": str(members),
            "unique_sequences": str(n),
            "collection_status": (
                "MSA_PREDICTION_QC_REQUIRED"
                if n >= 2
                else "STAGE08_ONLY_SINGLE_UNIQUE_SEQUENCE"
            ),
            "stage08_pdb_file": str(pdb.resolve()),
            "stage08_pdb_sha256": sha(pdb),
            "msa_prediction_provenance": "",
        })

    if (
        len(eligible)
        != prep_meta.get("msa_prediction_eligible_families")
        or len(ids) - len(eligible)
        != prep_meta.get("single_unique_sequence_families")
    ):
        raise ValueError("ELIGIBILITY_COUNT_MISMATCH")

    prediction_rows = []

    if args.prediction_index:
        if not eligible:
            raise ValueError(
                "PREDICTION_INDEX_NOT_ALLOWED_WITH_ZERO_ELIGIBLE"
            )

        prediction_rows = table(
            args.prediction_index,
            ("structure_id", "prediction_provenance"),
        )

        by_id = {
            r["structure_id"]: r
            for r in prediction_rows
        }

        if (
            not set(by_id).issubset(eligible)
            or (
                not args.allow_partial
                and set(by_id) != set(eligible)
            )
        ):
            raise ValueError(
                "PREDICTION_INDEX_ELIGIBLE_IDS_MISMATCH"
            )

        for sid in by_id:
            record_path = require(
                Path(by_id[sid]["prediction_provenance"])
            )
            record = metadata(record_path)

            a3m = (
                prep / "msa/a3m_colabfold" / f"{sid}.a3m"
            )

            historical = (
                record.get("provenance_type") == "HISTORICAL_MSA_REUSE"
                and record.get("status") == "HISTORICAL_IMPORTED"
            )
            if historical:
                if (
                    not args.allow_partial
                    or record.get("prepared_a3m_sha256") != sha(a3m)
                    or Path(record.get("prepared_a3m", "/")).resolve()
                    != a3m.resolve()
                    or record.get("gpu_prediction_performed") is not False
                    or record.get("original_run_provenance_verified") is not False
                    or record.get("a3m_sha256") == sha(a3m)
                    and record.get("historical_alignment_relation")
                    == "UNGAPPED_MULTISET_ONLY"
                ):
                    raise ValueError(f"INVALID_HISTORICAL_REUSE:{sid}")
                from collections import Counter

                def aligned_sequences(path):
                    result = []
                    with path.open(encoding="utf-8") as fh:
                        for line in fh:
                            line = line.strip()
                            if not line or line.startswith("#"):
                                continue
                            if line.startswith(">"):
                                result.append("")
                            elif result:
                                result[-1] += line
                            else:
                                raise ValueError(f"INVALID_A3M:{path}")
                    return result

                historical_a3m = require(
                    Path(record["source_a3m"])
                )
                if sha(historical_a3m) != record["a3m_sha256"]:
                    raise ValueError(f"HISTORICAL_A3M_SHA_MISMATCH:{sid}")
                seq_old = aligned_sequences(historical_a3m)
                seq_new = aligned_sequences(a3m)
                normalize = lambda x: x.replace("-", "").replace(".", "").upper()
                if (
                    not seq_old or not seq_new
                    or Counter(map(normalize, seq_old))
                    != Counter(map(normalize, seq_new))
                ):
                    raise ValueError(f"HISTORICAL_MEMBERSHIP_MISMATCH:{sid}")
            elif record.get("a3m_sha256") != sha(a3m):
                raise ValueError(
                    f"GPU_INPUT_A3M_SHA256_MISMATCH:{sid}"
                )

            returncode = record.get("colabfold_returncode")

            if returncode is not None and returncode != 0:
                raise ValueError(
                    f"GPU_PREDICTION_NOT_SUCCESSFUL:{sid}"
                )

            for row in status_rows:
                if row["structure_id"] == sid:
                    row["msa_prediction_provenance"] = (
                        str(record_path.resolve())
                    )
                    row["collection_status"] = (
                        "MSA_HISTORICAL_PREDICTION_QC_REQUIRED"
                        if historical else "MSA_PREDICTION_QC_REQUIRED"
                    )

    registered = {r["structure_id"] for r in prediction_rows}
    if args.allow_partial:
        for row in status_rows:
            if (
                row["structure_id"] in eligible
                and row["structure_id"] not in registered
            ):
                row["collection_status"] = "MSA_PREDICTION_PENDING"

    print(f"SELECTED_FAMILIES={len(ids)}")
    print(f"REGISTERED_PREDICTIONS={len(registered)}")
    print(f"PENDING_PREDICTIONS={len(eligible) - len(registered)}")
    print(f"MSA_ELIGIBLE_FAMILIES={len(eligible)}")
    print(
        "SINGLE_SEQUENCE_FALLBACK_FAMILIES="
        f"{len(ids) - len(eligible)}"
    )
    print("STAGE08_PDB_SHA256=PASS")
    print("STAGE11_PREPARATION_SHA256=PASS")

    return {
        "out": out,
        "prep": prep,
        "structures": structures_path,
        "eligible": eligible,
        "statuses": status_rows,
        "prediction_rows": prediction_rows,
        "prepared_provenance":
            prep / "preparation_provenance.json",
        "prediction_index": args.prediction_index,
        "allow_partial": args.allow_partial,
    }


def collect(state):
    out = state["out"]

    if (
        state["eligible"] and not state["prediction_rows"]
        and not state["allow_partial"]
    ):
        raise ValueError(
            "PREDICTION_INDEX_REQUIRED_FOR_ELIGIBLE_FAMILIES"
        )

    out.mkdir(parents=True, exist_ok=False)

    qc_file = None
    qc_prov = None

    if state["prediction_rows"]:
        index_path = out / "eligible_prediction_index.tsv"

        with index_path.open(
            "x", newline="", encoding="utf-8"
        ) as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "structure_id",
                    "prediction_provenance",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in state["prediction_rows"]:
                writer.writerow({
                    "structure_id": row["structure_id"],
                    "prediction_provenance":
                        row["prediction_provenance"],
                })

        qc_dir = out / "msa_qc"
        log = out / "msa_qc.log"

        command = [
            sys.executable, "-B", str(QC_SCRIPT),
            "--structures", str(state["structures"]),
            "--prediction-index", str(index_path),
            "--out-dir", str(qc_dir),
        ]

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
                f"MSA_STRUCTURE_QC_FAILED:SEE_LOG:{log}"
            )

        qc_file = require(
            qc_dir / "msa_structure_qc.tsv"
        )
        qc_prov = require(
            qc_dir / "msa_structure_qc_provenance.json"
        )

        qc_rows = table(
            qc_file,
            ("structure_id", "msa_pdb_sha256"),
        )

        if {
            r["structure_id"] for r in qc_rows
        } != {r["structure_id"] for r in state["prediction_rows"]}:
            raise ValueError("MSA_QC_FAMILY_COVERAGE_MISMATCH")

        for row in state["statuses"]:
            if row["msa_prediction_provenance"]:
                if row["collection_status"] == "MSA_HISTORICAL_PREDICTION_QC_REQUIRED":
                    row["collection_status"] = (
                        "MSA_HISTORICAL_PREDICTION_COLLECTED_QC_VERIFIED"
                    )
                else:
                    row["collection_status"] = (
                        "MSA_PREDICTION_COLLECTED_QC_VERIFIED"
                    )

        print("MSA_STRUCTURE_QC=PASS")
    else:
        print("MSA_STRUCTURE_QC=NOT_APPLICABLE")
        print("SINGLE_SEQUENCE_STAGE08_FALLBACK=PASS")

    status_path = out / "collection_status.tsv"

    with status_path.open(
        "x", newline="", encoding="utf-8"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "structure_id",
                "family_rep",
                "family_members",
                "unique_sequences",
                "collection_status",
                "stage08_pdb_file",
                "stage08_pdb_sha256",
                "msa_prediction_provenance",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(state["statuses"])

    record = {
        "stage": "STAGE11_COLLECT_OR_STAGE08_FALLBACK",
        "preparation_provenance_sha256":
            sha(state["prepared_provenance"]),
        "source_stage08_sha256":
            sha(state["structures"]),
        "collection_status_sha256":
            sha(status_path),
        "selected_families": len(state["statuses"]),
        "msa_collected_families":
            len(state["prediction_rows"]),
        "historical_imported_families": sum(
            r["collection_status"] ==
            "MSA_HISTORICAL_PREDICTION_COLLECTED_QC_VERIFIED"
            for r in state["statuses"]
        ),
        "msa_prediction_pending_families": sum(
            r["collection_status"] == "MSA_PREDICTION_PENDING"
            for r in state["statuses"]
        ),
        "stage08_fallback_families": sum(
            r["collection_status"] == "STAGE08_ONLY_SINGLE_UNIQUE_SEQUENCE"
            for r in state["statuses"]
        ),
        "msa_qc_sha256":
            sha(qc_file) if qc_file else None,
        "msa_qc_provenance_sha256":
            sha(qc_prov) if qc_prov else None,
        "gpu_prediction_performed": False,
        "msa_results_fabricated": False,
        "stage12_analysis_performed": False,
        "final_candidate_selection_performed": False,
    }

    with (
        out / "collection_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")

    print(
        "STAGE11_COLLECTION_STATUS="
        "VERIFIED_WITH_EXPLICIT_FALLBACK"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage08-output", type=Path, required=True
    )
    parser.add_argument(
        "--stage09-10-output", type=Path, required=True
    )
    parser.add_argument(
        "--prepared-dir", type=Path, required=True
    )
    parser.add_argument(
        "--prediction-index", type=Path
    )
    parser.add_argument(
        "--out-dir", type=Path, required=True
    )
    parser.add_argument(
        "--allow-partial", action="store_true",
        help="Permit a verified subset of eligible predictions; others remain pending",
    )
    parser.add_argument(
        "--mode",
        choices=("check", "collect"),
        default="check",
    )

    args = parser.parse_args()
    state = verify(args)

    if args.mode == "check":
        print("COLLECT_CHECK=PASS")
        print("OUTPUT_CREATED=NO")
        return

    collect(state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, csv.Error, TypeError
    ) as exc:
        sys.exit(f"ERROR: {exc}")
