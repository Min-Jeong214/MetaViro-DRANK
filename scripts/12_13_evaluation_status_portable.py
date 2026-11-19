#!/usr/bin/env python3
"""Preserve Stage 10 ranking and record Stage 11-13 evaluation states.

This is a status layer, not completed Stage 13 biological evidence.
No MSA prediction, Foldseek search or final selection is performed.
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


STATUS_COLUMNS = [
    "stage11_collection_status",
    "stage11_msa_qc_status",
    "stage12_conservation_status",
    "stage12_spatial_status",
    "stage13_msa_external_status",
    "stage13_integrated_evidence_status",
]

FALLBACK = "STAGE08_ONLY_SINGLE_UNIQUE_SEQUENCE"
COLLECTED = "MSA_PREDICTION_COLLECTED_QC_VERIFIED"
HISTORICAL = "MSA_HISTORICAL_PREDICTION_COLLECTED_QC_VERIFIED"
PENDING = "MSA_PREDICTION_PENDING"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_OR_EMPTY_INPUT:{path}")
    return path


def metadata(path):
    data = json.loads(require(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f"INVALID_METADATA:{path}")
    return data


def table(path):
    with require(path).open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not fields
        or len(fields) != len(set(fields))
        or "structure_id" not in fields
        or not rows
        or any(
            None in row
            or any(v is None for v in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    ids = [r["structure_id"] for r in rows]
    if (
        any(not x for x in ids)
        or len(ids) != len(set(ids))
    ):
        raise ValueError(f"INVALID_STRUCTURE_IDS:{path}")

    return fields, rows


def verify(args):
    up = args.stage08_output.resolve()
    down = args.stage09_10_output.resolve()
    prep = args.prepared_dir.resolve()
    collected = args.collected_dir.resolve()
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError("Refusing existing output directory")

    for directory in (up, down, prep, collected):
        if not directory.is_dir():
            raise ValueError(f"MISSING_INPUT_DIRECTORY:{directory}")
        if (
            out == directory
            or out in directory.parents
            or directory in out.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    rank = require(
        down / "10_ranking/confidence_aware_ranking.tsv"
    )
    shortlist = require(
        down / "10_shortlists/top15_experimental_candidates.tsv"
    )
    rank_prov = require(
        down / "10_ranking/ranking_provenance.json"
    )

    down_meta = metadata(
        down / "orchestration_provenance.json"
    )
    rank_meta = metadata(rank_prov)

    for name, path in (
        ("ranking", rank),
        ("top15", shortlist),
    ):
        entry = down_meta.get("artifacts", {}).get(name, {})

        if (
            entry.get("sha256") != sha(path)
            or Path(entry.get("path", "/")).resolve() != path
        ):
            raise ValueError(
                f"STAGE10_MANIFEST_MISMATCH:{name}"
            )

    rank_hashes = [
        value
        for key, value in rank_meta.items()
        if key.endswith("sha256") and isinstance(value, str)
    ]

    if sha(rank) not in rank_hashes:
        raise ValueError("STAGE10_RANKING_PROVENANCE_MISMATCH")

    rank_fields, ranks = table(rank)
    short_fields, short_rows = table(shortlist)

    if rank_fields != short_fields:
        raise ValueError("STAGE10_SHORTLIST_HEADER_MISMATCH")

    if any(field in rank_fields for field in STATUS_COLUMNS):
        raise ValueError("STATUS_COLUMNS_ALREADY_PRESENT")

    rank_by_id = {r["structure_id"]: r for r in ranks}

    for row in short_rows:
        sid = row["structure_id"]
        if sid not in rank_by_id or row != rank_by_id[sid]:
            raise ValueError("STAGE10_SHORTLIST_ROW_MISMATCH")

    prep_prov_path = require(
        prep / "preparation_provenance.json"
    )
    collect_prov_path = require(
        collected / "collection_provenance.json"
    )
    status_path = require(
        collected / "collection_status.tsv"
    )

    prep_meta = metadata(prep_prov_path)
    collect_meta = metadata(collect_prov_path)

    selected_map = require(
        prep / "selected_structure_map.tsv"
    )
    eligibility = require(
        prep / "prediction_eligibility.tsv"
    )

    for key, path in (
        ("stage10_ranking_sha256", rank),
        ("stage10_shortlist_sha256", shortlist),
        ("selected_structure_map_sha256", selected_map),
        ("prediction_eligibility_sha256", eligibility),
        (
            "family_members_sha256",
            up / "05_annotation/recurrent_members.tsv",
        ),
        (
            "member_fasta_sha256",
            up / "05_annotation/recurrent_members.faa",
        ),
    ):
        if prep_meta.get(key) != sha(require(path)):
            raise ValueError(f"PREPARATION_PROVENANCE_MISMATCH:{key}")

    structures = [
        p for p in (
            up / "08_structure"
        ).glob("run_*/collected/structure_results.tsv")
        if p.is_file()
        and sha(p) == prep_meta.get("stage08_structures_sha256")
    ]

    if len(structures) != 1:
        raise ValueError("STAGE08_STRUCTURE_SOURCE_MISMATCH")

    structures_file = structures[0]

    if (
        down_meta.get("stage08_structure_table_sha256")
        != sha(structures_file)
        or collect_meta.get("source_stage08_sha256")
        != sha(structures_file)
        or collect_meta.get("preparation_provenance_sha256")
        != sha(prep_prov_path)
    ):
        raise ValueError("CROSS_STAGE_SOURCE_PROVENANCE_MISMATCH")

    if collect_meta.get("collection_status_sha256") != sha(status_path):
        raise ValueError("COLLECTION_STATUS_SHA256_MISMATCH")

    if (
        collect_meta.get("gpu_prediction_performed") is not False
        or collect_meta.get("msa_results_fabricated") is not False
        or collect_meta.get("stage12_analysis_performed") is not False
        or collect_meta.get(
            "final_candidate_selection_performed"
        ) is not False
    ):
        raise ValueError("UNEXPECTED_COMPLETED_ANALYSIS_CLAIM")

    _, statuses = table(status_path)
    _, selected = table(selected_map)
    _, policy = table(eligibility)
    _, single_structures = table(structures_file)

    statuses_by_id = {
        r["structure_id"]: r for r in statuses
    }
    selected_by_id = {
        r["structure_id"]: r for r in selected
    }
    policy_by_id = {
        r["structure_id"]: r for r in policy
    }
    structures_by_id = {
        r["structure_id"]: r for r in single_structures
    }

    selected_ids = [
        r["structure_id"] for r in short_rows
    ]

    if not (
        set(selected_ids)
        == set(statuses_by_id)
        == set(selected_by_id)
        == set(policy_by_id)
    ):
        raise ValueError("STAGE10_STAGE11_SELECTED_ID_MISMATCH")

    output = []
    msa_ids = []
    historical_ids = []
    pending_ids = []
    fallback_count = 0

    for rank_row in short_rows:
        sid = rank_row["structure_id"]
        status = statuses_by_id[sid]
        selection = selected_by_id[sid]
        policy_row = policy_by_id[sid]

        if sid not in structures_by_id:
            raise ValueError(f"STAGE08_STRUCTURE_MISSING:{sid}")

        original = structures_by_id[sid]
        pdb = require(Path(original["pdb_file"]))

        if (
            sha(pdb) != original["pdb_sha256"]
            or sha(pdb) != status["stage08_pdb_sha256"]
            or Path(status["stage08_pdb_file"]).resolve()
            != pdb.resolve()
        ):
            raise ValueError(f"STAGE08_PDB_MISMATCH:{sid}")

        family = rank_row["family_rep"]
        n = int(status["unique_sequences"])
        members = int(status["family_members"])

        if not (
            family
            == status["family_rep"]
            == selection["family_rep"]
            == policy_row["family_rep"]
            == original["family_rep"]
            and n >= 1
            and members >= n
            and int(selection["family_size"]) == members
            and int(policy_row["family_members"]) == members
            and int(policy_row["unique_sequences"]) == n
        ):
            raise ValueError(f"FAMILY_MEMBERSHIP_MISMATCH:{sid}")

        if status["collection_status"] == FALLBACK:
            if (
                n != 1
                or policy_row["msa_prediction_status"]
                != "SINGLE_UNIQUE_SEQUENCE_NO_MSA_REPREDICTION"
                or status["msa_prediction_provenance"]
            ):
                raise ValueError(f"INVALID_FALLBACK:{sid}")

            stage11_qc = "NOT_EVALUATED"
            stage12 = "NOT_EVALUATED"
            stage13 = "NOT_EVALUATED"
            integrated = "STAGE10_ONLY_SINGLE_UNIQUE_SEQUENCE"
            fallback_count += 1

        elif status["collection_status"] == PENDING:
            if (
                n < 2
                or policy_row["msa_prediction_status"] != "MSA_PREDICTION_ELIGIBLE"
                or status["msa_prediction_provenance"]
            ):
                raise ValueError(f"INVALID_PENDING_PREDICTION:{sid}")
            stage11_qc = "PENDING_PREDICTION"
            stage12 = "PENDING_PREDICTION"
            stage13 = "PENDING_PREDICTION"
            integrated = "STAGE10_ONLY_PREDICTION_PENDING"
            pending_ids.append(sid)

        elif status["collection_status"] in (COLLECTED, HISTORICAL):
            if (
                n < 2
                or policy_row["msa_prediction_status"]
                != "MSA_PREDICTION_ELIGIBLE"
                or not status["msa_prediction_provenance"]
            ):
                raise ValueError(f"INVALID_MSA_COLLECTED_STATUS:{sid}")

            require(Path(status["msa_prediction_provenance"]))
            record_path = Path(status["msa_prediction_provenance"])
            record = metadata(record_path)
            is_historical = status["collection_status"] == HISTORICAL
            if is_historical != (
                record.get("provenance_type") == "HISTORICAL_MSA_REUSE"
                and record.get("status") == "HISTORICAL_IMPORTED"
            ):
                raise ValueError(f"PREDICTION_PROVENANCE_TYPE_MISMATCH:{sid}")
            stage11_qc = (
                "QC_VERIFIED_HISTORICAL" if is_historical else "QC_VERIFIED"
            )
            if is_historical:
                historical_ids.append(sid)
            stage12 = "PENDING_ANALYSIS"
            stage13 = "PENDING_SEARCH"
            integrated = "AWAITING_STAGE12_13_EVIDENCE"
            msa_ids.append(sid)

        else:
            raise ValueError(
                f"UNSUPPORTED_COLLECTION_STATUS:{sid}"
            )

        result = dict(rank_row)

        result.update({
            "stage11_collection_status":
                status["collection_status"],
            "stage11_msa_qc_status": stage11_qc,
            "stage12_conservation_status": stage12,
            "stage12_spatial_status": stage12,
            "stage13_msa_external_status": stage13,
            "stage13_integrated_evidence_status": integrated,
        })

        output.append(result)

    if (
        collect_meta.get("selected_families") != len(output)
        or collect_meta.get("msa_collected_families")
        != len(msa_ids)
        or collect_meta.get("stage08_fallback_families")
        != fallback_count
        or collect_meta.get("msa_prediction_pending_families", 0)
        != len(pending_ids)
        or collect_meta.get("historical_imported_families", 0)
        != len(historical_ids)
    ):
        raise ValueError("COLLECTION_COUNTS_MISMATCH")

    qc_file = collected / "msa_qc/msa_structure_qc.tsv"

    if msa_ids:
        qc_file = require(qc_file)
        qc_prov = require(
            collected / "msa_qc/msa_structure_qc_provenance.json"
        )
        qc_meta = metadata(qc_prov)
        index = require(
            collected / "eligible_prediction_index.tsv"
        )

        _, qc_rows = table(qc_file)

        if (
            {r["structure_id"] for r in qc_rows}
            != set(msa_ids)
            or collect_meta.get("msa_qc_sha256") != sha(qc_file)
            or collect_meta.get("msa_qc_provenance_sha256")
            != sha(qc_prov)
            or qc_meta.get("msa_qc_sha256") != sha(qc_file)
            or qc_meta.get("source_stage08_sha256")
            != sha(structures_file)
            or qc_meta.get("prediction_index_sha256") != sha(index)
        ):
            raise ValueError("MSA_QC_PROVENANCE_MISMATCH")
    else:
        if (
            qc_file.exists()
            or collect_meta.get("msa_qc_sha256") is not None
        ):
            raise ValueError("UNEXPECTED_MSA_QC_FOR_FALLBACK")

    print(f"SELECTED_FAMILIES={len(output)}")
    print(f"MSA_QC_VERIFIED_FAMILIES={len(msa_ids)}")
    print(f"STAGE08_ONLY_FAMILIES={fallback_count}")
    print(f"HISTORICAL_MSA_QC_FAMILIES={len(historical_ids)}")
    print(f"MSA_PREDICTION_PENDING_FAMILIES={len(pending_ids)}")
    print("STAGE10_RANKING_PRESERVED=PASS")
    print("STAGE11_COLLECTION_PROVENANCE=PASS")
    print("STAGE12_13_STATUS_POLICY=PASS")

    return {
        "out": out,
        "fields": rank_fields + STATUS_COLUMNS,
        "rows": output,
        "ranking": rank,
        "shortlist": shortlist,
        "preparation_provenance": prep_prov_path,
        "collection_provenance": collect_prov_path,
        "collection_status": status_path,
        "source_structures": structures_file,
        "msa_count": len(msa_ids),
        "fallback_count": fallback_count,
        "historical_count": len(historical_ids),
        "pending_count": len(pending_ids),
    }


def execute(state):
    out = state["out"]
    out.mkdir(parents=True, exist_ok=False)

    matrix = out / "stage10_ranked_evaluation_status.tsv"

    with matrix.open(
        "x", newline="", encoding="utf-8"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=state["fields"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(state["rows"])

    record = {
        "stage": "STAGE12_13_EVALUATION_STATUS",
        "ranking_sha256": sha(state["ranking"]),
        "shortlist_sha256": sha(state["shortlist"]),
        "stage08_structures_sha256":
            sha(state["source_structures"]),
        "preparation_provenance_sha256":
            sha(state["preparation_provenance"]),
        "collection_provenance_sha256":
            sha(state["collection_provenance"]),
        "collection_status_sha256":
            sha(state["collection_status"]),
        "status_matrix_sha256": sha(matrix),
        "families": len(state["rows"]),
        "columns": len(state["fields"]),
        "msa_qc_verified": state["msa_count"],
        "stage08_only": state["fallback_count"],
        "historical_msa_qc_verified": state["historical_count"],
        "msa_prediction_pending": state["pending_count"],
        "stage12_analysis_performed": False,
        "stage13_msa_external_search_performed": False,
        "stage10_scores_recalculated": False,
        "final_top5_top3_selected": False,
    }

    with (
        out / "evaluation_status_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")

    print(f"STATUS_COLUMNS={len(state['fields'])}")
    print(f"STATUS_FAMILIES={len(state['rows'])}")
    print("STAGE12_13_EVALUATION_STATUS=PASS")


def main():
    ap = argparse.ArgumentParser()

    for name in (
        "stage08_output",
        "stage09_10_output",
        "prepared_dir",
        "collected_dir",
        "out_dir",
    ):
        ap.add_argument(
            "--" + name.replace("_", "-"),
            type=Path,
            required=True,
        )

    ap.add_argument(
        "--mode",
        choices=("check", "run"),
        default="check",
    )

    args = ap.parse_args()
    state = verify(args)

    if args.mode == "check":
        print("CHECK_MODE=PASS")
        print("OUTPUT_CREATED=NO")
        return

    execute(state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, TypeError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
