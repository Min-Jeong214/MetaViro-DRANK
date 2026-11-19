#!/usr/bin/env python3
"""Run Stage 09-10 and adaptive Stage 11-13 from collected Stage 08.

This implementation uses verified Foldseek search reuse.
It never runs GPU prediction, fresh Foldseek searches, or
automatic final Top5/Top3 selection.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1048576), b""
        ):
            h.update(block)
    return h.hexdigest()


def require(path):
    path = Path(path).resolve()
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"MISSING_OR_EMPTY_FILE:{path}")
    return path


def read_json(path):
    record = json.loads(require(path).read_text())
    if not isinstance(record, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return record


def read_tsv(path):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as fh:
        reader = csv.DictReader(
            fh, delimiter="\t"
        )
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not fields
        or len(fields) != len(set(fields))
        or "structure_id" not in fields
        or not rows
        or any(
            None in row
            or any(value is None for value in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    ids = [
        row["structure_id"] for row in rows
    ]

    if (
        any(not sid for sid in ids)
        or len(ids) != len(set(ids))
    ):
        raise ValueError(
            f"INVALID_STRUCTURE_IDS:{path}"
        )

    return fields, rows


def invoke(label, script, arguments, log=None):
    command = [
        sys.executable,
        "-B",
        str(SCRIPTS / script),
        *map(str, arguments),
    ]

    print(f"[RUN] {label}", flush=True)

    if log is None:
        result = subprocess.run(
            command,
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode:
            raise RuntimeError(
                f"{label}: "
                + (result.stderr or result.stdout)
            )
    else:
        log.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with log.open(
            "x", encoding="utf-8"
        ) as fh:
            result = subprocess.run(
                command,
                cwd=ROOT,
                stdout=fh,
                stderr=subprocess.STDOUT,
                check=False,
            )

        if result.returncode:
            tail = log.read_text(
                encoding="utf-8",
                errors="replace",
            )[-5000:]

            raise RuntimeError(
                f"{label} FAILED\n"
                f"LOG={log}\n{tail}"
            )

    print(f"[PASS] {label}", flush=True)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--stage08-output",
        type=Path,
        required=True,
    )
    search_group = ap.add_mutually_exclusive_group(
        required=True
    )
    search_group.add_argument(
        "--reuse-search-root", type=Path
    )
    search_group.add_argument(
        "--db-root", type=Path
    )
    ap.add_argument(
        "--prediction-index",
        type=Path,
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=2,
    )
    ap.add_argument(
        "--mode",
        choices=("check", "run"),
        default="check",
    )

    args = ap.parse_args()

    up = args.stage08_output.resolve()
    search = (
        args.reuse_search_root or args.db_root
    ).resolve()
    search_option = (
        "--reuse-search-root"
        if args.reuse_search_root
        else "--db-root"
    )
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    if args.threads < 1:
        raise ValueError("INVALID_THREADS")

    for source in (up, search):
        if not source.is_dir():
            raise ValueError(
                f"MISSING_INPUT_DIRECTORY:{source}"
            )

        if (
            out == source
            or source in out.parents
            or out in source.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    require(
        up / "05_annotation/recurrent_members.tsv"
    )
    require(
        up / "05_annotation/recurrent_members.faa"
    )

    prediction_index = None

    if args.prediction_index:
        prediction_index = require(
            args.prediction_index
        )

        if (
            prediction_index == out
            or out in prediction_index.parents
            or prediction_index in out.parents
        ):
            raise ValueError(
                "PREDICTION_INDEX_OUTPUT_OVERLAP"
            )

    stage09_dir = out / "09_10"
    adaptive_dir = out / "11_13_adaptive"

    stage09_args = [
        "--stage08-output", up,
        search_option, search,
        "--out-dir", stage09_dir,
        "--threads", args.threads,
    ]

    invoke(
        "09-10 verified search preflight",
        "09_10_orchestrator_portable.py",
        stage09_args + ["--mode", "check"],
    )

    if args.mode == "check":
        print("STAGE09_10_PREFLIGHT=PASS")
        print(
            "STAGE11_13_CHECK="
            "DEFERRED_UNTIL_STAGE10_RANKING"
        )
        print("OUTPUT_CREATED=NO")
        print("STAGE08_13_ADAPTIVE_CHECK=PASS")
        return

    out.mkdir(
        parents=True, exist_ok=False
    )

    invoke(
        "09-10 external evidence and ranking",
        "09_10_orchestrator_portable.py",
        stage09_args + ["--mode", "run"],
        out / "logs/09_10.log",
    )

    rank = require(
        stage09_dir /
        "10_ranking/confidence_aware_ranking.tsv"
    )
    shortlist = require(
        stage09_dir /
        "10_shortlists/top15_experimental_candidates.tsv"
    )
    stage09_manifest = require(
        stage09_dir / "orchestration_provenance.json"
    )

    m09 = read_json(stage09_manifest)

    for name, path in (
        ("ranking", rank),
        ("top15", shortlist),
    ):
        entry = m09.get(
            "artifacts", {}
        ).get(name, {})

        if (
            entry.get("sha256") != sha(path)
            or Path(
                entry.get("path", "/")
            ).resolve() != path
        ):
            raise ValueError(
                f"STAGE09_10_ARTIFACT_MISMATCH:{name}"
            )

    rank_fields, rank_rows = read_tsv(rank)
    top_fields, selected_rows = read_tsv(shortlist)

    if (
        len(rank_fields) != 30
        or top_fields != rank_fields
        or m09.get("family_count")
        != len(rank_rows)
    ):
        raise ValueError(
            "STAGE10_RANKING_FORMAT_MISMATCH"
        )

    rank_index = {
        row["structure_id"]: row
        for row in rank_rows
    }

    for row in selected_rows:
        if rank_index.get(
            row["structure_id"]
        ) != row:
            raise ValueError(
                "STAGE10_SHORTLIST_RANKING_MISMATCH"
            )

    print(
        f"STAGE10_RANKED_FAMILIES={len(rank_rows)}"
    )
    print(
        f"STAGE10_SELECTED_FAMILIES={len(selected_rows)}"
    )

    adaptive_args = [
        "--stage08-output", up,
        "--stage09-10-output", stage09_dir,
        "--out-dir", adaptive_dir,
        "--threads", args.threads,
        "--mode", "run",
    ]

    if prediction_index:
        adaptive_args.extend([
            "--prediction-index", prediction_index
        ])

    invoke(
        "11-13 adaptive preparation and status",
        "11_13_adaptive_orchestrator_portable.py",
        adaptive_args,
        out / "logs/11_13_adaptive.log",
    )

    adaptive_manifest = require(
        adaptive_dir / "orchestration_provenance.json"
    )

    m11 = read_json(adaptive_manifest)

    if (
        m11.get("stage10_ranking_sha256")
        != sha(rank)
        or m11.get("selected_families")
        != len(selected_rows)
        or m11.get("gpu_prediction_performed")
        is not False
        or m11.get("foldseek_search_performed")
        is not False
        or m11.get("stage10_scores_recalculated")
        is not False
        or m11.get("final_top5_top3_selected")
        is not False
    ):
        raise ValueError(
            "ADAPTIVE_STAGE10_HANDOFF_MISMATCH"
        )

    for name, entry in m11.get(
        "artifacts", {}
    ).items():
        path = require(entry["path"])

        if entry.get("sha256") != sha(path):
            raise ValueError(
                f"ADAPTIVE_ARTIFACT_SHA256_MISMATCH:{name}"
            )

        if (
            path != adaptive_dir
            and adaptive_dir not in path.parents
        ):
            raise ValueError(
                f"ADAPTIVE_ARTIFACT_OUTSIDE_OUTPUT:{name}"
            )

    state = m11.get("workflow_state")

    valid_states = {
        "AWAITING_MSA_PREDICTION",
        "STAGE10_ONLY_SINGLE_SEQUENCE_FALLBACK",
        "STAGE11_QC_VERIFIED_STAGE12_13_PENDING",
    }

    if state not in valid_states:
        raise ValueError(
            f"UNKNOWN_ADAPTIVE_STATE:{state}"
        )

    result_file = (
        adaptive_dir /
        "12_13_status/"
        "stage10_ranked_evaluation_status.tsv"
    )

    if state == "AWAITING_MSA_PREDICTION":
        if (
            result_file.exists()
            or not m11.get(
                "msa_prediction_eligible_families"
            )
        ):
            raise ValueError(
                "INVALID_GPU_PENDING_STATE"
            )

        status_sha = None
        status_columns = None

    else:
        result_file = require(result_file)

        status_fields, status_rows = read_tsv(
            result_file
        )

        if (
            len(status_fields) != 36
            or status_fields[:30] != rank_fields
            or len(status_rows)
            != len(selected_rows)
        ):
            raise ValueError(
                "STAGE12_13_STATUS_FORMAT_MISMATCH"
            )

        selected_index = {
            row["structure_id"]: row
            for row in selected_rows
        }

        if (
            set(selected_index)
            != {
                row["structure_id"]
                for row in status_rows
            }
        ):
            raise ValueError(
                "STAGE12_13_SELECTED_FAMILIES_MISMATCH"
            )

        for row in status_rows:
            original = selected_index[
                row["structure_id"]
            ]

            if any(
                row[field] != original[field]
                for field in rank_fields
            ):
                raise ValueError(
                    "STAGE10_RANKING_CHANGED_IN_STATUS"
                )

        status_sha = sha(result_file)
        status_columns = len(status_fields)

    # STAGE12_CONSERVATION_INTEGRATION_STEP126
    # Only QC-verified MSA families undergo conservation analysis.
    # Families awaiting GPU prediction remain pending.
    conservation_output = out / "12_conservation"
    conservation_manifest_sha = None
    conservation_status_sha = None
    conservation_executed = False
    conservation_families = 0

    if state != "AWAITING_MSA_PREDICTION":
        invoke(
            "12 adaptive conservation analysis",
            "12_conservation_adaptive_orchestrator_portable.py",
            [
                "--stage09-10-output", stage09_dir,
                "--adaptive-output", adaptive_dir,
                "--out-dir", conservation_output,
                "--mode", "run",
            ],
            out / "logs/12_conservation.log",
        )

        conservation_manifest = require(
            conservation_output /
            "orchestration_provenance.json"
        )
        conservation_status = require(
            conservation_output /
            "stage12_evaluation_status.tsv"
        )

        cm = read_json(conservation_manifest)
        cfields, crows = read_tsv(conservation_status)
        pfields, prows = read_tsv(result_file)

        if (
            cfields != pfields
            or len(cfields) != 36
            or len(crows) != len(prows)
        ):
            raise ValueError(
                "INTEGRATED_CONSERVATION_STATUS_SCHEMA_MISMATCH"
            )

        prior = {
            row["structure_id"]: row
            for row in prows
        }

        if (
            len(prior) != len(prows)
            or set(prior) != {
                row["structure_id"]
                for row in crows
            }
        ):
            raise ValueError(
                "INTEGRATED_CONSERVATION_FAMILY_MISMATCH"
            )

        completed = 0

        for row in crows:
            old = prior[row["structure_id"]]

            for field in pfields:
                if field == "stage12_conservation_status":
                    continue

                if row[field] != old[field]:
                    raise ValueError(
                        "STAGE10_OR_STAGE11_STATUS_MODIFIED:"
                        + row["structure_id"]
                        + ":"
                        + field
                    )

            before = old["stage12_conservation_status"]
            after = row["stage12_conservation_status"]

            if before == "NOT_EVALUATED":
                if after != "NOT_EVALUATED":
                    raise ValueError(
                        "FALSE_SINGLE_SEQUENCE_CONSERVATION"
                    )
            elif before == "PENDING_ANALYSIS":
                if after != "COMPLETED":
                    raise ValueError(
                        "ELIGIBLE_CONSERVATION_NOT_COMPLETED"
                    )
                completed += 1
            else:
                raise ValueError(
                    "UNEXPECTED_PREVIOUS_CONSERVATION_STATUS"
                )

        if (
            cm.get("stage11_manifest_sha256")
            != sha(adaptive_manifest)
            or cm.get("stage10_ranking_sha256")
            != sha(rank)
            or cm.get("stage11_status_sha256")
            != sha(result_file)
            or cm.get("stage12_status_sha256")
            != sha(conservation_status)
            or cm.get("selected_families")
            != len(crows)
            or cm.get("conservation_completed_families")
            != completed
            or cm.get(
                "stage12_conservation_analysis_performed"
            ) is not bool(completed)
            or cm.get("stage10_scores_recalculated")
            is not False
            or cm.get("gpu_prediction_performed")
            is not False
        ):
            raise ValueError(
                "INTEGRATED_CONSERVATION_PROVENANCE_MISMATCH"
            )

        conservation_manifest_sha = sha(
            conservation_manifest
        )
        conservation_status_sha = sha(
            conservation_status
        )
        conservation_executed = bool(completed)
        conservation_families = completed

    # STAGE12_SPATIAL_INTEGRATION_STEP128
    # Geometry is exploratory evidence, never biochemical proof.
    spatial_output = out / "12_spatial"
    spatial_manifest_sha = None
    spatial_status_sha = None
    spatial_evaluated = 0
    exploratory_cluster_count = 0

    if state != "AWAITING_MSA_PREDICTION":
        invoke(
            "12 adaptive exploratory spatial analysis",
            "12_spatial_adaptive_orchestrator_portable.py",
            [
                "--conservation-output",
                conservation_output,
                "--out-dir",
                spatial_output,
                "--mode",
                "run",
            ],
            out / "logs/12_spatial.log",
        )

        spatial_manifest = require(
            spatial_output /
            "orchestration_provenance.json"
        )
        spatial_status = require(
            spatial_output /
            "stage12_spatial_evaluation_status.tsv"
        )
        conservation_status_file = require(
            conservation_output /
            "stage12_evaluation_status.tsv"
        )
        conservation_manifest_file = require(
            conservation_output /
            "orchestration_provenance.json"
        )

        sm = read_json(spatial_manifest)

        sfields, srows = read_tsv(spatial_status)
        cfields, crows = read_tsv(
            conservation_status_file
        )

        if (
            sfields != cfields
            or len(sfields) != 36
            or len(srows) != len(crows)
        ):
            raise ValueError(
                "SPATIAL_EVALUATION_SCHEMA_MISMATCH"
            )

        old_index = {
            row["structure_id"]: row
            for row in crows
        }

        if (
            len(old_index) != len(crows)
            or set(old_index) != {
                row["structure_id"]
                for row in srows
            }
        ):
            raise ValueError(
                "SPATIAL_FAMILY_COVERAGE_MISMATCH"
            )

        evaluated = 0

        for row in srows:
            sid = row["structure_id"]
            old = old_index[sid]

            for field in cfields:
                if field == "stage12_spatial_status":
                    continue

                if row[field] != old[field]:
                    raise ValueError(
                        "SPATIAL_UNEXPECTED_FIELD_MODIFICATION:"
                        + sid + ":" + field
                    )

            previous = old["stage12_spatial_status"]
            current = row["stage12_spatial_status"]

            if previous == "NOT_EVALUATED":
                if current != "NOT_EVALUATED":
                    raise ValueError(
                        "SINGLE_SEQUENCE_FALSE_SPATIAL_EVIDENCE:"
                        + sid
                    )

            elif previous == "PENDING_ANALYSIS":
                if current != "COMPLETED_EXPLORATORY":
                    raise ValueError(
                        "SPATIAL_ANALYSIS_NOT_COMPLETED:"
                        + sid
                    )

                if (
                    old["stage12_conservation_status"]
                    != "COMPLETED"
                ):
                    raise ValueError(
                        "SPATIAL_WITHOUT_CONSERVATION:"
                        + sid
                    )

                evaluated += 1

            else:
                raise ValueError(
                    "UNEXPECTED_SPATIAL_INPUT_STATUS:"
                    + sid
                )

        if (
            sm.get(
                "stage12_conservation_manifest_sha256"
            ) != sha(conservation_manifest_file)
            or sm.get(
                "stage12_conservation_status_sha256"
            ) != sha(conservation_status_file)
            or sm.get(
                "stage12_spatial_status_sha256"
            ) != sha(spatial_status)
            or sm.get("selected_families")
            != len(srows)
            or sm.get("spatial_evaluated_families")
            != evaluated
            or sm.get(
                "spatial_not_evaluated_families"
            ) != len(srows) - evaluated
            or sm.get("biochemical_activity_validated")
            is not False
            or sm.get("stage13_analysis_performed")
            is not False
            or sm.get("gpu_prediction_performed")
            is not False
            or sm.get("foldseek_search_performed")
            is not False
            or sm.get("stage10_scores_recalculated")
            is not False
            or sm.get("final_top5_top3_selected")
            is not False
        ):
            raise ValueError(
                "SPATIAL_PROVENANCE_MISMATCH"
            )

        clusters = sm.get(
            "exploratory_cluster_count"
        )

        if (
            type(clusters) is not int
            or clusters < 0
        ):
            raise ValueError(
                "INVALID_EXPLORATORY_CLUSTER_COUNT"
            )

        if evaluated:
            actual_spatial_prov = require(
                spatial_output /
                "spatial/spatial_cluster_provenance.json"
            )

            if (
                sm.get("spatial_provenance_sha256")
                != sha(actual_spatial_prov)
            ):
                raise ValueError(
                    "SPATIAL_ANALYSIS_ARTIFACT_MISMATCH"
                )
        else:
            if (
                (spatial_output / "spatial").exists()
                or clusters != 0
            ):
                raise ValueError(
                    "SPATIAL_FALLBACK_CREATED_FALSE_OUTPUT"
                )

        spatial_manifest_sha = sha(
            spatial_manifest
        )
        spatial_status_sha = sha(
            spatial_status
        )
        spatial_evaluated = evaluated
        exploratory_cluster_count = clusters

        print(
            "STAGE12_SPATIAL_EVALUATED_FAMILIES="
            f"{spatial_evaluated}"
        )
        print(
            "STAGE12_EXPLORATORY_CLUSTERS="
            f"{exploratory_cluster_count}"
        )

    # STAGE13_INTEGRATION_STEP133
    # No historical external evidence is automatically attached.
    # Missing MSA Foldseek evidence stays PENDING_SEARCH.
    stage13_output = out / "13_adaptive"
    stage13_manifest_sha = None
    stage13_status_sha = None
    stage13_state = (
        "AWAITING_MSA_PREDICTION"
        if state == "AWAITING_MSA_PREDICTION"
        else None
    )

    if state != "AWAITING_MSA_PREDICTION":
        invoke(
            "13 adaptive external evidence handoff",
            "13_adaptive_orchestrator_portable.py",
            [
                "--stage09-10-output", stage09_dir,
                "--adaptive-output", adaptive_dir,
                "--conservation-output",
                conservation_output,
                "--spatial-output", spatial_output,
                "--out-dir", stage13_output,
                "--mode", "run",
            ],
            out / "logs/13_adaptive.log",
        )

        stage13_manifest = require(
            stage13_output /
            "orchestration_provenance.json"
        )
        stage13_status = require(
            stage13_output /
            "stage13_evaluation_status.tsv"
        )
        spatial_status_file = require(
            spatial_output /
            "stage12_spatial_evaluation_status.tsv"
        )

        m13 = read_json(stage13_manifest)
        f13, r13 = read_tsv(stage13_status)
        fs, rs = read_tsv(spatial_status_file)

        if (
            f13 != fs
            or len(f13) != 36
            or r13 != rs
            or sha(stage13_status)
            != sha(spatial_status_file)
        ):
            raise ValueError(
                "STAGE13_PENDING_STATUS_HANDOFF_MISMATCH"
            )

        expected_state = (
            "AWAITING_MSA_EXTERNAL_SEARCH"
            if spatial_evaluated > 0
            else "STAGE13_NOT_EVALUATED_SINGLE_UNIQUE"
        )

        if (
            m13.get("workflow_state")
            != expected_state
            or m13.get("stage10_ranking_sha256")
            != sha(rank)
            or m13.get(
                "stage12_spatial_manifest_sha256"
            ) != sha(spatial_manifest)
            or m13.get(
                "stage12_spatial_status_sha256"
            ) != sha(spatial_status_file)
            or m13.get("stage13_status_sha256")
            != sha(stage13_status)
            or m13.get("selected_families")
            != len(selected_rows)
            or m13.get(
                "msa_external_eligible_families"
            ) != spatial_evaluated
            or m13.get("verified_external_families")
            != 0
            or m13.get("external_mode")
            != "NOT_EXECUTED"
            or m13.get("final_33col_evidence_sha256")
            is not None
            or m13.get("new_foldseek_search_performed")
            is not False
            or m13.get("gpu_prediction_performed")
            is not False
            or m13.get("stage10_scores_recalculated")
            is not False
            or m13.get("final_top5_top3_selected")
            is not False
        ):
            raise ValueError(
                "STAGE13_HANDOFF_PROVENANCE_MISMATCH"
            )

        if (
            stage13_output / "13_ranked"
        ).exists():
            raise ValueError(
                "STAGE13_UNVERIFIED_FINAL_EVIDENCE_CREATED"
            )

        stage13_manifest_sha = sha(
            stage13_manifest
        )
        stage13_status_sha = sha(
            stage13_status
        )
        stage13_state = expected_state

        print(
            "STAGE13_ADAPTIVE_STATE="
            f"{stage13_state}"
        )

    else:
        if stage13_output.exists():
            raise ValueError(
                "GPU_PENDING_STAGE13_OUTPUT_UNEXPECTED"
            )

    record = {
        "stage": "STAGE08_13_ADAPTIVE_HANDOFF",
        "search_mode": m09.get("search_mode"),
        "stage09_10_manifest_sha256":
            sha(stage09_manifest),
        "stage11_13_manifest_sha256":
            sha(adaptive_manifest),
        "stage10_ranking_sha256": sha(rank),
        "stage10_shortlist_sha256": sha(shortlist),
        "evaluation_status_sha256": status_sha,
        "stage12_conservation_manifest_sha256":
            conservation_manifest_sha,
        "stage12_conservation_status_sha256":
            conservation_status_sha,
        "stage12_conservation_analysis_performed":
            conservation_executed,
        "stage12_conservation_completed_families":
            conservation_families,
        "stage12_spatial_manifest_sha256":
            spatial_manifest_sha,
        "stage12_spatial_status_sha256":
            spatial_status_sha,
        "stage12_spatial_evaluated_families":
            spatial_evaluated,
        "stage12_exploratory_cluster_count":
            exploratory_cluster_count,
        "biochemical_activity_validated": False,
        "stage13_manifest_sha256":
            stage13_manifest_sha,
        "stage13_status_sha256":
            stage13_status_sha,
        "stage13_workflow_state":
            stage13_state,
        "stage13_verified_external_families": 0,
        "stage13_final_33col_sha256": None,
        "ranked_families": len(rank_rows),
        "selected_families": len(selected_rows),
        "evaluation_status_columns": status_columns,
        "workflow_state": state,
        "stage12_13_analysis_completed": False,
        "gpu_prediction_performed": False,
        "foldseek_search_performed":
            args.db_root is not None,
        "final_top5_top3_selected": False,
    }

    with (
        out / "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")

    print(f"ADAPTIVE_WORKFLOW_STATE={state}")
    print(
        f"STAGE10_SELECTED_FAMILIES={len(selected_rows)}"
    )
    print(f"OUTPUT_DIRECTORY={out}")
    print("STAGE08_TO_13_ADAPTIVE_HANDOFF=PASS")


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, TypeError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
