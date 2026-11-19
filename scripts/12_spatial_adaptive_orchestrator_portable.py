#!/usr/bin/env python3
"""Adaptive exploratory spatial clustering after Stage 12 conservation.

Only QC-supported, conservation-completed families are analyzed.
No geometry result is treated as biochemical activity validation.
"""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ENGINE = (
    Path(__file__).resolve().parent /
    "12_charged_spatial_clusters_exploratory.py"
)


def digest(path):
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
        raise ValueError(f"MISSING_INPUT:{path}")
    return path


def metadata(path):
    result = json.loads(require(path).read_text())
    if not isinstance(result, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return result


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
        or not rows
        or any(
            None in row
            or any(value is None for value in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    if "structure_id" not in fields:
        raise ValueError(
            f"MISSING_STRUCTURE_ID_COLUMN:{path}"
        )

    ids = [r["structure_id"] for r in rows]

    if (
        any(not sid for sid in ids)
        or len(ids) != len(set(ids))
    ):
        raise ValueError(
            f"DUPLICATE_STRUCTURE_IDS:{path}"
        )

    return fields, rows


def validate(args):
    source = args.conservation_output.resolve()
    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    if not source.is_dir():
        raise ValueError(
            f"MISSING_CONSERVATION_OUTPUT:{source}"
        )

    if (
        output == source
        or source in output.parents
        or output in source.parents
    ):
        raise ValueError("INPUT_OUTPUT_OVERLAP")

    source_manifest = require(
        source / "orchestration_provenance.json"
    )
    source_status = require(
        source / "stage12_evaluation_status.tsv"
    )

    manifest = metadata(source_manifest)
    fields, rows = table(source_status)

    if (
        len(fields) != 36
        or manifest.get("stage12_status_sha256")
        != digest(source_status)
        or manifest.get("selected_families")
        != len(rows)
        or manifest.get(
            "stage10_scores_recalculated"
        ) is not False
        or manifest.get("gpu_prediction_performed")
        is not False
    ):
        raise ValueError(
            "STAGE12_STATUS_PROVENANCE_MISMATCH"
        )

    eligible = []

    for row in rows:
        conservation = row[
            "stage12_conservation_status"
        ]
        spatial = row["stage12_spatial_status"]

        if conservation == "NOT_EVALUATED":
            if spatial != "NOT_EVALUATED":
                raise ValueError(
                    "INVALID_SINGLE_SEQUENCE_FALLBACK"
                )

        elif conservation == "COMPLETED":
            if (
                spatial != "PENDING_ANALYSIS"
                or row["stage11_msa_qc_status"]
                != "QC_VERIFIED"
            ):
                raise ValueError(
                    "INVALID_MULTISEQUENCE_ELIGIBILITY"
                )
            eligible.append(row["structure_id"])

        else:
            raise ValueError(
                "UNSUPPORTED_CONSERVATION_STATE:"
                + conservation
            )

    if (
        manifest.get(
            "conservation_completed_families"
        ) != len(eligible)
    ):
        raise ValueError(
            "CONSERVATION_FAMILY_COUNT_MISMATCH"
        )

    conservation_dir = source / "conservation"
    conservation_prov = None
    summary = None

    if eligible:
        conservation_prov = require(
            conservation_dir /
            "conservation_provenance.json"
        )
        summary = require(
            conservation_dir /
            "candidate_conservation_summary.tsv"
        )

        if (
            manifest.get(
                "conservation_provenance_sha256"
            ) != digest(conservation_prov)
            or manifest.get(
                "conservation_summary_sha256"
            ) != digest(summary)
        ):
            raise ValueError(
                "CONSERVATION_ARTIFACT_SHA256_MISMATCH"
            )

        _, summary_rows = table(summary)

        if {
            r["structure_id"] for r in summary_rows
        } != set(eligible):
            raise ValueError(
                "CONSERVATION_FAMILY_COVERAGE_MISMATCH"
            )

    return {
        "source": source,
        "output": output,
        "manifest": source_manifest,
        "status": source_status,
        "fields": fields,
        "rows": rows,
        "eligible": eligible,
        "conservation_dir": conservation_dir,
        "conservation_prov": conservation_prov,
    }


def run(state):
    output = state["output"]
    output.mkdir(
        parents=True, exist_ok=False
    )

    eligible = state["eligible"]
    spatial_prov = None
    spatial_summary = None
    spatial_clusters = None
    cluster_count = 0

    if eligible:
        spatial_output = output / "spatial"
        log = output / "spatial_analysis.log"

        command = [
            sys.executable,
            "-B",
            str(ENGINE),
            "--conservation-dir",
            str(state["conservation_dir"]),
            "--out-dir",
            str(spatial_output),
        ]

        print(
            "[RUN] Stage 12 exploratory spatial clustering",
            flush=True
        )

        with log.open(
            "x", encoding="utf-8"
        ) as stream:
            process = subprocess.run(
                command,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=False,
            )

        if process.returncode:
            tail = log.read_text(
                encoding="utf-8",
                errors="replace",
            )[-6000:]

            raise RuntimeError(
                "SPATIAL_ANALYSIS_FAILED\n"
                f"LOG={log}\n{tail}"
            )

        spatial_prov = require(
            spatial_output /
            "spatial_cluster_provenance.json"
        )
        spatial_summary = require(
            spatial_output /
            "family_spatial_cluster_summary.tsv"
        )
        spatial_clusters = require(
            spatial_output /
            "charged_spatial_clusters.tsv"
        )

        record = metadata(spatial_prov)
        _, family_rows = table(spatial_summary)

        if (
            record.get("stage12_provenance_sha256")
            != digest(
                state["conservation_prov"]
            )
            or record.get("summary_table_sha256")
            != digest(spatial_summary)
            or record.get("cluster_table_sha256")
            != digest(spatial_clusters)
            or record.get("family_count")
            != len(eligible)
            or {
                row["structure_id"]
                for row in family_rows
            } != set(eligible)
            or record.get(
                "biochemical_activity_inferred"
            ) is not False
        ):
            raise ValueError(
                "SPATIAL_RESULT_PROVENANCE_MISMATCH"
            )

        cluster_count = record["cluster_count"]

        if (
            not isinstance(cluster_count, int)
            or cluster_count < 0
        ):
            raise ValueError(
                "INVALID_SPATIAL_CLUSTER_COUNT"
            )

        print(
            f"SPATIAL_CLUSTER_COUNT={cluster_count}",
            flush=True
        )
        print(
            "SPATIAL_ANALYSIS_PROVENANCE=VERIFIED",
            flush=True
        )

    else:
        print(
            "SPATIAL_ANALYSIS=NOT_EVALUATED",
            flush=True
        )

    status_out = (
        output /
        "stage12_spatial_evaluation_status.tsv"
    )

    if not eligible:
        shutil.copyfile(
            state["status"], status_out
        )
    else:
        modified_rows = []

        for row in state["rows"]:
            updated = dict(row)

            if row["structure_id"] in eligible:
                updated[
                    "stage12_spatial_status"
                ] = "COMPLETED_EXPLORATORY"

            modified_rows.append(updated)

        with status_out.open(
            "x", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=state["fields"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(modified_rows)

    provenance = {
        "stage": "STAGE12_ADAPTIVE_SPATIAL_EXPLORATORY",
        "stage12_conservation_manifest_sha256":
            digest(state["manifest"]),
        "stage12_conservation_status_sha256":
            digest(state["status"]),
        "stage12_spatial_status_sha256":
            digest(status_out),
        "spatial_provenance_sha256":
            digest(spatial_prov)
            if spatial_prov else None,
        "spatial_summary_sha256":
            digest(spatial_summary)
            if spatial_summary else None,
        "spatial_clusters_sha256":
            digest(spatial_clusters)
            if spatial_clusters else None,
        "selected_families":
            len(state["rows"]),
        "spatial_evaluated_families":
            len(eligible),
        "spatial_not_evaluated_families":
            len(state["rows"]) - len(eligible),
        "exploratory_cluster_count":
            cluster_count,
        "biochemical_activity_validated": False,
        "stage13_analysis_performed": False,
        "gpu_prediction_performed": False,
        "foldseek_search_performed": False,
        "stage10_scores_recalculated": False,
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
        "SPATIAL_EVALUATED_FAMILIES="
        f"{len(eligible)}"
    )
    print("STAGE12_ADAPTIVE_SPATIAL=PASS")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--conservation-output",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--mode",
        choices=("check", "run"),
        default="check",
    )

    args = ap.parse_args()
    state = validate(args)

    if args.mode == "check":
        print(
            "SPATIAL_ELIGIBLE_FAMILIES="
            f"{len(state['eligible'])}"
        )
        print("OUTPUT_CREATED=NO")
        print("STAGE12_SPATIAL_CHECK=PASS")
        return

    run(state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, TypeError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
