#!/usr/bin/env python3
"""MetaViro-DRANK unified entrypoint.

Existing mode:
    Reuse completed Stage 08, verified searches, and adaptive Stage 11-13.

Fresh mode:
    Execute the existing Stage 00-08 runner, then invoke downstream.
    Requires compatible precomputed Foldseek search provenance.

Does not automatically execute new MSA GPU prediction,
new Foldseek search, Stage 12-13 biological analysis, or Top5/Top3.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UPSTREAM = ROOT / "run_pipeline.py"
DOWNSTREAM = (
    ROOT / "scripts/08_13_adaptive_orchestrator_portable.py"
)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(
            lambda: stream.read(1048576), b""
        ):
            h.update(chunk)
    return h.hexdigest()


def require(path):
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_OR_EMPTY_FILE:{path}")
    return path


def record(path):
    obj = json.loads(require(path).read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"INVALID_MANIFEST:{path}")
    return obj


def table(path):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as stream:
        reader = csv.DictReader(
            stream, delimiter="\t"
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

    ids = [row["structure_id"] for row in rows]

    if (
        any(not sid for sid in ids)
        or len(ids) != len(set(ids))
    ):
        raise ValueError(f"DUPLICATE_STRUCTURE_ID:{path}")

    return fields, rows


def run(label, command, log=None):
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
            parents=True, exist_ok=True
        )

        with log.open("x", encoding="utf-8") as stream:
            result = subprocess.run(
                command,
                cwd=ROOT,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )

        if result.returncode:
            tail = log.read_text(
                encoding="utf-8",
                errors="replace",
            )[-5000:]

            raise RuntimeError(
                f"{label}:FAILED\nLOG={log}\n{tail}"
            )

    print(f"[PASS] {label}", flush=True)


def stage08_inputs(upstream):
    if not upstream.is_dir():
        raise ValueError(
            f"STAGE08_OUTPUT_DIRECTORY_MISSING:{upstream}"
        )

    members = require(
        upstream /
        "05_annotation/recurrent_members.tsv"
    )
    proteins = require(
        upstream /
        "05_annotation/recurrent_members.faa"
    )

    structure_tables = sorted(
        (upstream / "08_structure").glob(
            "run_*/collected/structure_results.tsv"
        )
    )

    if len(structure_tables) != 1:
        raise ValueError(
            "STAGE08_COLLECTED_STRUCTURE_TABLE_NOT_UNIQUE"
        )

    structures = require(structure_tables[0])
    _, rows = table(structures)

    return {
        "members": members,
        "proteins": proteins,
        "structures": structures,
        "families": len(rows),
    }


def fresh_output(config):
    try:
        import yaml
    except ImportError as exc:
        raise ValueError(
            "PYYAML_REQUIRED_FOR_CONFIG_MODE"
        ) from exc

    settings = yaml.safe_load(
        require(config).read_text(encoding="utf-8")
    )

    if not isinstance(settings, dict):
        raise ValueError("INVALID_CONFIG")

    value = settings.get("output_dir")

    if not isinstance(value, str) or not value.strip():
        raise ValueError("CONFIG_OUTPUT_DIR_MISSING")

    output = Path(value).expanduser()

    if not output.is_absolute():
        output = ROOT / output

    return output.resolve()


def validate_disjoint(output, sources):
    if output.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    for source in sources:
        source = source.resolve()

        if (
            output == source
            or output in source.parents
            or source in output.parents
        ):
            raise ValueError(
                f"INPUT_OUTPUT_OVERLAP:{source}"
            )


def main():
    ap = argparse.ArgumentParser()

    group = ap.add_mutually_exclusive_group(
        required=True
    )

    group.add_argument(
        "--existing-stage08-output",
        type=Path,
    )
    group.add_argument(
        "--config",
        type=Path,
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
        "--structure-mode",
        choices=("collect", "predict"),
        default="collect",
    )
    ap.add_argument(
        "--structure-prediction-dir",
        type=Path,
    )
    ap.add_argument(
        "--structure-env",
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

    output = args.out_dir.resolve()
    search = (
        args.reuse_search_root or args.db_root
    ).resolve()
    search_option = (
        "--reuse-search-root"
        if args.reuse_search_root
        else "--db-root"
    )

    if args.threads < 1:
        raise ValueError("INVALID_THREADS")

    if not search.is_dir():
        raise ValueError(
            "SEARCH_REUSE_DIRECTORY_MISSING"
        )

    config = None

    if args.existing_stage08_output:
        upstream = (
            args.existing_stage08_output.resolve()
        )
        source_mode = "EXISTING_STAGE08"
    else:
        config = require(args.config)
        upstream = fresh_output(config)
        source_mode = "RUN_STAGE00_08"

        if args.reuse_search_root:
            raise ValueError(
                "FRESH_MODE_REQUIRES_DB_ROOT"
            )

        if (
            args.structure_mode == "collect"
            and not args.structure_prediction_dir
        ):
            raise ValueError(
                "STRUCTURE_PREDICTION_DIR_REQUIRED"
            )

    sources = [upstream, search]

    if config:
        sources.append(config)

    if args.prediction_index:
        sources.append(
            require(args.prediction_index)
        )

    if args.structure_prediction_dir:
        prediction_dir = (
            args.structure_prediction_dir.resolve()
        )
        if not prediction_dir.is_dir():
            raise ValueError(
                "STRUCTURE_PREDICTION_DIR_MISSING"
            )
        sources.append(prediction_dir)

    validate_disjoint(output, sources)

    if config:
        if upstream.exists():
            raise ValueError(
                "Refusing existing upstream output directory"
            )

        upstream_command = [
            sys.executable,
            "-B",
            str(UPSTREAM),
            "--config",
            str(config),
        ]

        # Build identical Stage 08 options for check and run.
        # Previously --mode check silently used Stage 02 defaults
        # and skipped the Stage 05-08 database preflight.
        upstream_command.extend([
            "--through", "08",
            "--structure-mode", args.structure_mode,
        ])

        if args.structure_prediction_dir:
            upstream_command.extend([
                "--structure-prediction-dir",
                str(prediction_dir),
            ])

        if args.structure_env:
            upstream_command.extend([
                "--structure-env", args.structure_env,
            ])

        if args.mode == "check":
            run(
                "Stage 00-08 configuration preflight",
                upstream_command + ["--check"],
            )

            print("UPSTREAM_PREFLIGHT=PASS")
            print(
                "DOWNSTREAM_PREFLIGHT="
                "DEFERRED_UNTIL_STAGE08_COLLECTION"
            )
            print("OUTPUT_CREATED=NO")
            print("FULL_PIPELINE_CHECK=PASS")
            return

    else:
        stage08_inputs(upstream)

    downstream_output = output / "downstream"

    downstream_command = [
        sys.executable,
        "-B",
        str(DOWNSTREAM),
        "--stage08-output",
        str(upstream),
        search_option,
        str(search),
        "--out-dir",
        str(downstream_output),
        "--threads",
        str(args.threads),
    ]

    if args.prediction_index:
        downstream_command.extend([
            "--prediction-index",
            str(args.prediction_index.resolve()),
        ])

    if args.mode == "check":
        run(
            "Stage 08-13 downstream preflight",
            downstream_command + [
                "--mode", "check"
            ],
        )

        print("UPSTREAM_MODE=EXISTING_STAGE08")
        print("OUTPUT_CREATED=NO")
        print("FULL_PIPELINE_CHECK=PASS")
        return

    output.mkdir(
        parents=True,
        exist_ok=False,
    )

    if config:
        run(
            "Stage 00-08 upstream pipeline",
            upstream_command,
            output / "logs/stage00_08.log",
        )

    inputs = stage08_inputs(upstream)

    run(
        "Stage 09-13 adaptive pipeline",
        downstream_command + [
            "--mode", "run"
        ],
        output / "logs/stage09_13.log",
    )

    downstream_manifest = require(
        downstream_output /
        "orchestration_provenance.json"
    )
    down_meta = record(
        downstream_manifest
    )

    if (
        down_meta.get("ranked_families")
        != inputs["families"]
        or down_meta.get(
            "gpu_prediction_performed"
        ) is not False
        or down_meta.get(
            "foldseek_search_performed"
        ) != (args.db_root is not None)
        or down_meta.get(
            "final_top5_top3_selected"
        ) is not False
    ):
        raise ValueError(
            "FULL_PIPELINE_MANIFEST_MISMATCH"
        )

    ranking = require(
        downstream_output /
        "09_10/10_ranking/"
        "confidence_aware_ranking.tsv"
    )

    rank_fields, rank_rows = table(ranking)

    if (
        len(rank_fields) != 30
        or len(rank_rows) != inputs["families"]
    ):
        raise ValueError(
            "FULL_PIPELINE_RANKING_MISMATCH"
        )

    evaluation = (
        downstream_output /
        "11_13_adaptive/12_13_status/"
        "stage10_ranked_evaluation_status.tsv"
    )

    workflow_state = down_meta[
        "workflow_state"
    ]

    if workflow_state == "AWAITING_MSA_PREDICTION":
        if evaluation.exists():
            raise ValueError(
                "UNEXPECTED_EVALUATION_STATUS_FILE"
            )
        evaluation_sha = None
    else:
        evaluation = require(evaluation)
        fields, rows = table(evaluation)

        if (
            len(fields) != 36
            or len(rows) != down_meta.get(
                "selected_families"
            )
        ):
            raise ValueError(
                "FULL_PIPELINE_STATUS_MISMATCH"
            )

        evaluation_sha = sha(evaluation)

    # STAGE13_INTEGRATION_STEP133
    stage13_status_sha = None
    stage13_state = None

    stage13_manifest_path = (
        downstream_output /
        "13_adaptive/orchestration_provenance.json"
    )
    stage13_status_path = (
        downstream_output /
        "13_adaptive/stage13_evaluation_status.tsv"
    )

    if workflow_state == "AWAITING_MSA_PREDICTION":
        if (
            stage13_manifest_path.exists()
            or stage13_status_path.exists()
            or down_meta.get("stage13_status_sha256")
            is not None
        ):
            raise ValueError(
                "STAGE13_EVIDENCE_WHILE_GPU_PENDING"
            )

        stage13_state = (
            down_meta.get("stage13_workflow_state")
        )

        if stage13_state != "AWAITING_MSA_PREDICTION":
            raise ValueError(
                "STAGE13_GPU_PENDING_STATE_MISMATCH"
            )

    else:
        stage13_manifest_path = require(
            stage13_manifest_path
        )
        stage13_status_path = require(
            stage13_status_path
        )

        stage13_meta = record(
            stage13_manifest_path
        )

        final_fields, final_rows = table(
            stage13_status_path
        )

        if (
            len(final_fields) != 36
            or final_fields[:30] != rank_fields
            or len(final_rows)
            != down_meta["selected_families"]
            or stage13_meta.get(
                "stage13_status_sha256"
            ) != sha(stage13_status_path)
            or stage13_meta.get(
                "stage10_ranking_sha256"
            ) != sha(ranking)
            or down_meta.get(
                "stage13_manifest_sha256"
            ) != sha(stage13_manifest_path)
            or down_meta.get(
                "stage13_status_sha256"
            ) != sha(stage13_status_path)
            or stage13_meta.get(
                "verified_external_families"
            ) != 0
            or stage13_meta.get(
                "final_33col_evidence_sha256"
            ) is not None
        ):
            raise ValueError(
                "FULL_STAGE13_PROVENANCE_MISMATCH"
            )

        rank_by_id = {
            row["structure_id"]: row
            for row in rank_rows
        }

        for row in final_rows:
            original = rank_by_id.get(
                row["structure_id"]
            )

            if original is None:
                raise ValueError(
                    "STAGE13_STRUCTURE_NOT_IN_RANKING"
                )

            for field in rank_fields:
                if row[field] != original[field]:
                    raise ValueError(
                        "STAGE13_MODIFIED_STAGE10_RANKING"
                    )

        stage13_state = stage13_meta.get(
            "workflow_state"
        )

        expected_stage13_state = (
            "STAGE13_NOT_EVALUATED_SINGLE_UNIQUE"
            if workflow_state
            == "STAGE10_ONLY_SINGLE_SEQUENCE_FALLBACK"
            else "AWAITING_MSA_EXTERNAL_SEARCH"
        )

        if (
            stage13_state != expected_stage13_state
            or down_meta.get(
                "stage13_workflow_state"
            ) != stage13_state
        ):
            raise ValueError(
                "FULL_STAGE13_WORKFLOW_STATE_MISMATCH"
            )

        stage13_status_sha = sha(
            stage13_status_path
        )

    provenance = {
        "stage": "METAVIRO_DRANK_UNIFIED_ENTRYPOINT",
        "source_mode": source_mode,
        "upstream_config_sha256":
            sha(config) if config else None,
        "stage08_structure_results_sha256":
            sha(inputs["structures"]),
        "stage08_family_members_sha256":
            sha(inputs["members"]),
        "stage08_member_fasta_sha256":
            sha(inputs["proteins"]),
        "downstream_manifest_sha256":
            sha(downstream_manifest),
        "stage10_ranking_sha256":
            sha(ranking),
        "evaluation_status_sha256":
            evaluation_sha,
        "stage13_workflow_state":
            stage13_state,
        "stage13_status_sha256":
            stage13_status_sha,
        "stage13_final_33col_sha256": None,
        "families": inputs["families"],
        "selected_families":
            down_meta["selected_families"],
        "workflow_state": workflow_state,
        "search_mode": down_meta.get("search_mode"),
        "new_foldseek_search_performed":
            args.db_root is not None,
        "msa_gpu_prediction_performed":
            False,
        "stage12_13_biological_analysis_completed":
            False,
        "final_top5_top3_selected": False,
    }

    with (
        output /
        "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as stream:
        json.dump(
            provenance, stream, indent=2
        )
        stream.write("\n")

    print(
        f"FULL_PIPELINE_SOURCE_MODE={source_mode}"
    )
    print(
        f"FULL_PIPELINE_FAMILIES={inputs['families']}"
    )
    print(
        f"FULL_PIPELINE_STATE={workflow_state}"
    )
    print(
        f"FULL_PIPELINE_OUTPUT={output}"
    )
    print("FULL_PIPELINE_ENTRYPOINT=PASS")


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, TypeError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
