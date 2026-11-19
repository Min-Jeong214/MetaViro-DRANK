#!/usr/bin/env python3
"""Adaptive Stage 11-13 preparation and evaluation-status workflow.

Reuses verified Stage 10 rankings. Does not perform ColabFold GPU
prediction, Foldseek searching, or final candidate selection.
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

PREP = "11_prepare_orchestrator_portable.py"
COLLECT = "11_collect_orchestrator_portable.py"
STATUS = "12_13_evaluation_status_portable.py"
IMPORT = "11_import_historical_reuse_portable.py"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def require(path):
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"MISSING_OR_EMPTY_FILE:{path}")
    return path


def read_tsv(path):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not fields
        or "structure_id" not in fields
        or len(fields) != len(set(fields))
        or not rows
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    ids = [r["structure_id"] for r in rows]

    if (
        any(not sid for sid in ids)
        or len(ids) != len(set(ids))
        or any(
            None in row
            or any(v is None for v in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"MALFORMED_TABLE:{path}")

    return fields, rows


def meta(path):
    record = json.loads(require(path).read_text())
    if not isinstance(record, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return record


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
            parents=True, exist_ok=True
        )

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
                f"{label} failed; inspect {log}"
            )

    print(f"[PASS] {label}", flush=True)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--stage08-output",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--stage09-10-output",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--prediction-index",
        type=Path,
    )
    ap.add_argument(
        "--reuse-register",
        type=Path,
        help="Historical prediction reuse register with verified artifact hashes",
    )
    ap.add_argument(
        "--prepared-dir",
        type=Path,
        help="Reuse a previously generated Stage 11 preparation (no MAFFT rerun)",
    )
    ap.add_argument(
        "--allow-partial-predictions",
        action="store_true",
        help="With a prediction index, permit missing eligible predictions to remain pending",
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
    down = args.stage09_10_output.resolve()
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    for source in (up, down):
        if not source.is_dir():
            raise ValueError(
                f"MISSING_SOURCE_DIRECTORY:{source}"
            )
        if (
            out == source
            or source in out.parents
            or out in source.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    if args.threads < 1:
        raise ValueError("INVALID_THREADS")
    if args.prediction_index and args.reuse_register:
        raise ValueError("PREDICTION_INDEX_AND_REUSE_REGISTER_CONFLICT")

    prediction_index = None

    if args.prediction_index:
        prediction_index = require(
            args.prediction_index.resolve()
        )

        if (
            prediction_index == out
            or out in prediction_index.parents
            or prediction_index in out.parents
        ):
            raise ValueError(
                "PREDICTION_INDEX_OUTPUT_OVERLAP"
            )

    reuse_register = None
    if args.reuse_register:
        reuse_register = require(args.reuse_register.resolve())
        if (
            reuse_register == out or out in reuse_register.parents
            or reuse_register in out.parents
        ):
            raise ValueError("REUSE_REGISTER_OUTPUT_OVERLAP")

    prep_dir = (
        args.prepared_dir.resolve()
        if args.prepared_dir else out / "11_prepare"
    )
    if args.prepared_dir and (
        not prep_dir.is_dir() or out == prep_dir
        or out in prep_dir.parents or prep_dir in out.parents
    ):
        raise ValueError("PREPARED_DIR_MISSING_OR_OVERLAP")
    collect_dir = out / "11_collect"
    status_dir = out / "12_13_status"
    import_dir = out / "11_historical_import"

    prep_args = [
        "--stage08-output", up,
        "--stage09-10-output", down,
        "--out-dir", prep_dir,
        "--threads", args.threads,
    ]

    if args.prepared_dir:
        require(prep_dir / "preparation_provenance.json")
        require(prep_dir / "prediction_eligibility.tsv")
        print("PREPARATION_SOURCE=EXISTING_VERIFIED_INPUT")
    else:
        invoke(
            "11 preparation preflight",
            PREP,
            prep_args + ["--mode", "check"],
        )

    if args.mode == "check":
        if reuse_register and args.prepared_dir:
            invoke("historical import preflight", IMPORT, [
                "--reuse-register", reuse_register,
                "--prepared-dir", prep_dir,
                "--stage08-output", up,
                "--out-dir", import_dir,
                "--mode", "check",
            ])
        print("PREPARATION_PREFLIGHT=PASS")
        print("COLLECT_AND_STATUS_PREFLIGHT=DEFERRED_UNTIL_COLLECTION")
        print("CHECK_MODE_NO_WRITE=PASS")
        return

    out.mkdir(parents=True, exist_ok=False)

    if not args.prepared_dir:
        invoke(
            "11 family MSA preparation",
            PREP,
            prep_args + ["--mode", "run"],
            out / "logs/11_prepare.log",
        )

    eligibility_path = require(
        prep_dir / "prediction_eligibility.tsv"
    )
    prepare_prov_path = require(
        prep_dir / "preparation_provenance.json"
    )
    prepare_meta = meta(prepare_prov_path)

    _, eligibility = read_tsv(
        eligibility_path
    )

    valid_statuses = {
        "MSA_PREDICTION_ELIGIBLE",
        "SINGLE_UNIQUE_SEQUENCE_NO_MSA_REPREDICTION",
    }

    if any(
        row["msa_prediction_status"]
        not in valid_statuses
        for row in eligibility
    ):
        raise ValueError(
            "UNEXPECTED_ELIGIBILITY_STATUS"
        )

    eligible = sum(
        row["msa_prediction_status"]
        == "MSA_PREDICTION_ELIGIBLE"
        for row in eligibility
    )

    if (
        prepare_meta.get("families")
        != len(eligibility)
        or prepare_meta.get(
            "msa_prediction_eligible_families"
        ) != eligible
        or prepare_meta.get(
            "single_unique_sequence_families"
        ) != len(eligibility) - eligible
        or prepare_meta.get(
            "prediction_eligibility_sha256"
        ) != sha(eligibility_path)
    ):
        raise ValueError(
            "PREPARATION_ELIGIBILITY_PROVENANCE_MISMATCH"
        )

    artifacts = {
        "preparation_provenance": prepare_prov_path,
        "prediction_eligibility": eligibility_path,
    }

    if reuse_register:
        invoke("11 historical artifact import", IMPORT, [
            "--reuse-register", reuse_register,
            "--prepared-dir", prep_dir,
            "--stage08-output", up,
            "--out-dir", import_dir,
            "--mode", "run",
        ], out / "logs/11_historical_import.log")
        prediction_index = require(import_dir / "prediction_index.tsv")
        artifacts["historical_import_provenance"] = require(
            import_dir / "import_provenance.json"
        )
        artifacts["historical_import_index"] = prediction_index

    state = "AWAITING_MSA_PREDICTION"
    partial = bool(reuse_register or args.allow_partial_predictions)

    if eligible == 0 and prediction_index:
        raise ValueError(
            "PREDICTION_INDEX_NOT_ALLOWED_WITH_ZERO_ELIGIBLE"
        )

    if eligible and prediction_index is None:
        print(
            "ELIGIBLE_FAMILIES_REQUIRE_GPU_RESULTS="
            f"{eligible}"
        )
        print(
            "NEXT_ACTION=PROVIDE_VERIFIED_PREDICTION_INDEX"
        )

    else:
        collect_args = [
            "--stage08-output", up,
            "--stage09-10-output", down,
            "--prepared-dir", prep_dir,
            "--out-dir", collect_dir,
            "--mode", "collect",
        ]

        if prediction_index:
            collect_args.extend(["--prediction-index", prediction_index])
        if partial:
            collect_args.append("--allow-partial")

        invoke(
            "11 verified structure collection",
            COLLECT,
            collect_args,
            out / "logs/11_collect.log",
        )

        status_args = [
            "--stage08-output", up,
            "--stage09-10-output", down,
            "--prepared-dir", prep_dir,
            "--collected-dir", collect_dir,
            "--out-dir", status_dir,
            "--mode", "run",
        ]

        invoke(
            "12-13 evaluation status",
            STATUS,
            status_args,
            out / "logs/12_13_status.log",
        )

        collect_prov = require(
            collect_dir /
            "collection_provenance.json"
        )
        status_prov = require(
            status_dir /
            "evaluation_status_provenance.json"
        )
        status_table = require(
            status_dir /
            "stage10_ranked_evaluation_status.tsv"
        )

        collect_meta = meta(collect_prov)
        status_meta = meta(status_prov)
        columns, rows = read_tsv(status_table)

        if (
            len(columns) != 36
            or len(rows) != len(eligibility)
            or status_meta.get(
                "status_matrix_sha256"
            ) != sha(status_table)
            or status_meta.get(
                "collection_provenance_sha256"
            ) != sha(collect_prov)
            or collect_meta.get(
                "preparation_provenance_sha256"
            ) != sha(prepare_prov_path)
            or status_meta.get(
                "stage10_scores_recalculated"
            ) is not False
        ):
            raise ValueError(
                "ADAPTIVE_CHAIN_PROVENANCE_MISMATCH"
            )

        artifacts.update({
            "collection_provenance": collect_prov,
            "evaluation_status_provenance": status_prov,
            "evaluation_status_table": status_table,
        })

        verified = collect_meta.get("msa_collected_families", 0)
        pending = collect_meta.get("msa_prediction_pending_families", 0)
        state = (
            "PARTIAL_HISTORICAL_QC_VERIFIED_STAGE12_13_PENDING"
            if partial and verified and pending
            else "PARTIAL_PREDICTION_PENDING" if partial and pending
            else "STAGE11_QC_VERIFIED_STAGE12_13_PENDING"
            if verified else "STAGE10_ONLY_SINGLE_SEQUENCE_FALLBACK"
        )

    manifest = {
        "stage": "STAGE11_13_ADAPTIVE_WORKFLOW",
        "workflow_state": state,
        "selected_families": len(eligibility),
        "msa_prediction_eligible_families":
            eligible,
        "stage08_fallback_families":
            len(eligibility) - eligible,
        "historical_imported_families": (
            collect_meta.get("historical_imported_families", 0)
            if "collect_meta" in locals() else 0
        ),
        "msa_prediction_pending_families": (
            collect_meta.get("msa_prediction_pending_families", eligible)
            if "collect_meta" in locals() else eligible
        ),
        "historical_original_execution_verified": False,
        "stage10_ranking_sha256":
            prepare_meta["stage10_ranking_sha256"],
        "artifacts": {
            key: {
                "path": str(path.resolve()),
                "sha256": sha(path),
            }
            for key, path in artifacts.items()
        },
        "gpu_prediction_performed": False,
        "foldseek_search_performed": False,
        "stage12_13_analysis_completed": False,
        "stage10_scores_recalculated": False,
        "final_top5_top3_selected": False,
    }

    manifest_path = (
        out / "orchestration_provenance.json"
    )

    with manifest_path.open(
        "x", encoding="utf-8"
    ) as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    print(
        f"SELECTED_FAMILIES={len(eligibility)}"
    )
    print(f"MSA_ELIGIBLE_FAMILIES={eligible}")
    print(f"ADAPTIVE_WORKFLOW_STATE={state}")
    print(f"ADAPTIVE_OUTPUT={out}")
    print("STAGE11_13_ADAPTIVE_ORCHESTRATION=PASS")


if __name__ == "__main__":
    try:
        main()
    except (
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        TypeError,
        csv.Error,
    ) as exc:
        sys.exit(f"ERROR: {exc}")
