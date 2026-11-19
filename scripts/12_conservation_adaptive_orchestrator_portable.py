#!/usr/bin/env python3
"""Run Stage 12 conservation for QC-verified MSA families only.

Preserves original Stage 10 ranking and status data.
Never executes GPU prediction, Foldseek or spatial analysis.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "12_msa_structure_conservation_portable.py"

QC_OK = "MSA_PREDICTION_COLLECTED_QC_VERIFIED"
FALLBACK = "STAGE08_ONLY_SINGLE_UNIQUE_SEQUENCE"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(path):
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_INPUT:{path}")
    return path


def metadata(path):
    obj = json.loads(require(path).read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return obj


def table(path):
    with require(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not fields
        or len(fields) != len(set(fields))
        or "structure_id" not in fields
        or not rows
        or any(None in row or any(
            value is None for value in row.values()
        ) for row in rows)
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    ids = [row["structure_id"] for row in rows]

    if any(not x for x in ids) or len(ids) != len(set(ids)):
        raise ValueError(f"DUPLICATE_STRUCTURE_ID:{path}")

    return fields, rows


def validate(args):
    adaptive = args.adaptive_output.resolve()
    down = args.stage09_10_output.resolve()
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError("Refusing existing output directory")

    for folder in (adaptive, down):
        if not folder.is_dir():
            raise ValueError(f"MISSING_DIRECTORY:{folder}")
        if (
            out == folder
            or folder in out.parents
            or out in folder.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    manifest_path = require(
        adaptive / "orchestration_provenance.json"
    )
    manifest = metadata(manifest_path)

    ranking = require(
        down / "10_ranking/confidence_aware_ranking.tsv"
    )
    shortlist = require(
        down / "10_shortlists/top15_experimental_candidates.tsv"
    )

    if manifest.get("stage10_ranking_sha256") != sha(ranking):
        raise ValueError("STAGE10_RANKING_PROVENANCE_MISMATCH")

    for name, item in manifest.get("artifacts", {}).items():
        artifact = require(item["path"])
        if (
            item.get("sha256") != sha(artifact)
            or (
                artifact != adaptive
                and adaptive not in artifact.parents
            )
        ):
            raise ValueError(
                f"ADAPTIVE_ARTIFACT_INTEGRITY_MISMATCH:{name}"
            )

    fields, statuses = table(
        adaptive / "12_13_status/"
        "stage10_ranked_evaluation_status.tsv"
    )
    rank_fields, rank_rows = table(ranking)
    _, selected = table(shortlist)

    if (
        len(fields) != 36
        or len(rank_fields) != 30
        or fields[:30] != rank_fields
        or len(selected) != len(statuses)
    ):
        raise ValueError("STATUS_OR_RANKING_FORMAT_MISMATCH")

    selected_by_id = {
        r["structure_id"]: r for r in selected
    }
    ranking_by_id = {
        r["structure_id"]: r for r in rank_rows
    }

    if (
        len(selected_by_id) != len(selected)
        or set(selected_by_id) != {
            r["structure_id"] for r in statuses
        }
    ):
        raise ValueError("STAGE10_STAGE11_SELECTION_MISMATCH")

    prep_dir = adaptive / "11_prepare"
    collect_dir = adaptive / "11_collect"

    prep_path = require(
        prep_dir / "preparation_provenance.json"
    )
    collect_path = require(
        collect_dir / "collection_provenance.json"
    )

    prep = metadata(prep_path)
    collected = metadata(collect_path)

    if (
        prep.get("stage10_ranking_sha256") != sha(ranking)
        or prep.get("stage10_shortlist_sha256") != sha(shortlist)
        or collected.get("preparation_provenance_sha256")
        != sha(prep_path)
    ):
        raise ValueError("STAGE11_HANDOFF_SHA256_MISMATCH")

    eligible = []
    output = []

    for row in statuses:
        sid = row["structure_id"]

        if selected_by_id[sid] != {
            key: row[key] for key in rank_fields
        }:
            raise ValueError(
                f"STAGE10_RANKING_CHANGED:{sid}"
            )

        if ranking_by_id.get(sid) != selected_by_id[sid]:
            raise ValueError(
                f"STAGE10_SELECTION_CHANGED:{sid}"
            )

        record = dict(row)
        status = row["stage11_collection_status"]

        if status == FALLBACK:
            if (
                row["stage12_conservation_status"]
                != "NOT_EVALUATED"
                or row["stage11_msa_qc_status"]
                != "NOT_EVALUATED"
            ):
                raise ValueError(
                    f"INVALID_FALLBACK_STATUS:{sid}"
                )

        elif status == QC_OK:
            if (
                row["stage11_msa_qc_status"]
                != "QC_VERIFIED"
                or row["stage12_conservation_status"]
                != "PENDING_ANALYSIS"
            ):
                raise ValueError(
                    f"INVALID_MSA_STATUS:{sid}"
                )
            eligible.append(sid)

        else:
            raise ValueError(
                f"UNSUPPORTED_COLLECTION_STATUS:{sid}"
            )

        output.append(record)

    if (
        manifest.get("selected_families") != len(statuses)
        or manifest.get("msa_prediction_eligible_families")
        != len(eligible)
        or collected.get("msa_collected_families")
        != len(eligible)
    ):
        raise ValueError("ADAPTIVE_FAMILY_COUNT_MISMATCH")

    index = None
    aln_dir = None
    a3m_dir = None

    if eligible:
        index = require(
            collect_dir / "eligible_prediction_index.tsv"
        )
        _, index_rows = table(index)
        _, qc_rows = table(
            collect_dir / "msa_qc/msa_structure_qc.tsv"
        )

        if (
            set(r["structure_id"] for r in index_rows)
            != set(eligible)
            or set(r["structure_id"] for r in qc_rows)
            != set(eligible)
        ):
            raise ValueError("MSA_QC_COVERAGE_MISMATCH")

        aln_dir = prep_dir / "msa/mafft"
        # Stage 12 must use the exact A3M used for GPU prediction.
        a3m_dir = prep_dir / "msa/a3m_colabfold"

        prepared_inputs = prep.get("colabfold_inputs", {})

        for row in index_rows:
            sid = row["structure_id"]
            a3m = require(
                prep_dir / "msa/a3m_colabfold" /
                f"{sid}.a3m"
            )
            require(aln_dir / f"{sid}.aligned.faa")
            require(a3m_dir / f"{sid}.a3m")

            prediction = metadata(
                row["prediction_provenance"]
            )

            if (
                prediction.get("status") != "COMPLETED"
                or prediction.get("colabfold_returncode") != 0
                or prediction.get("a3m_sha256") != sha(a3m)
                or prepared_inputs.get(sid, {}).get("sha256")
                != sha(a3m)
            ):
                raise ValueError(
                    f"MSA_PREDICTION_INPUT_MISMATCH:{sid}"
                )

    return {
        "out": out,
        "eligible": eligible,
        "rows": output,
        "fields": fields,
        "index": index,
        "aln_dir": aln_dir,
        "a3m_dir": a3m_dir,
        "stage11_manifest": manifest_path,
        "stage10_ranking": ranking,
        "source_status": (
            adaptive / "12_13_status/"
            "stage10_ranked_evaluation_status.tsv"
        ),
    }


def execute(state):
    out = state["out"]
    out.mkdir(parents=True, exist_ok=False)

    eligible = state["eligible"]

    selection = None
    provenance = None
    summary = None

    if eligible:
        selection = out / "msa_selection.tsv"

        with selection.open(
            "x", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(
                stream,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writerow(["structure_id"])
            for sid in eligible:
                writer.writerow([sid])

        output = out / "conservation"

        command = [
            sys.executable,
            "-B",
            str(ANALYSIS),
            "--selection",
            str(selection),
            "--aln-dir",
            str(state["aln_dir"]),
            "--a3m-dir",
            str(state["a3m_dir"]),
            "--prediction-index",
            str(state["index"]),
            "--out-dir",
            str(output),
        ]

        log = out / "stage12_conservation.log"

        print(
            "[RUN] Stage 12 conservation analysis",
            flush=True,
        )

        with log.open("x", encoding="utf-8") as stream:
            result = subprocess.run(
                command,
                cwd=HERE.parent,
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
                "STAGE12_CONSERVATION_FAILED\n"
                f"LOG={log}\n{tail}"
            )

        provenance = require(
            output / "conservation_provenance.json"
        )
        summary = require(
            output / "candidate_conservation_summary.tsv"
        )

        meta = metadata(provenance)
        _, summary_rows = table(summary)

        if (
            meta.get("selection_sha256")
            != sha(selection)
            or meta.get("prediction_index_sha256")
            != sha(state["index"])
            or set(r["structure_id"] for r in summary_rows)
            != set(eligible)
        ):
            raise ValueError(
                "STAGE12_CONSERVATION_RESULT_MISMATCH"
            )

        for row in state["rows"]:
            if row["structure_id"] in eligible:
                row["stage12_conservation_status"] = (
                    "COMPLETED"
                )

        print(
            "STAGE12_CONSERVATION_RESULT=VERIFIED",
            flush=True,
        )
    else:
        print(
            "STAGE12_CONSERVATION=NOT_EVALUATED",
            flush=True,
        )

    status = out / "stage12_evaluation_status.tsv"

    with status.open(
        "x", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=state["fields"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(state["rows"])

    result = {
        "stage": "STAGE12_ADAPTIVE_CONSERVATION",
        "stage11_manifest_sha256":
            sha(state["stage11_manifest"]),
        "stage10_ranking_sha256":
            sha(state["stage10_ranking"]),
        "stage11_status_sha256":
            sha(state["source_status"]),
        "stage12_status_sha256": sha(status),
        "selection_sha256":
            sha(selection) if selection else None,
        "conservation_provenance_sha256":
            sha(provenance) if provenance else None,
        "conservation_summary_sha256":
            sha(summary) if summary else None,
        "selected_families": len(state["rows"]),
        "conservation_completed_families":
            len(eligible),
        "conservation_not_evaluated_families":
            len(state["rows"]) - len(eligible),
        "stage12_conservation_analysis_performed":
            bool(eligible),
        "stage12_spatial_analysis_performed": False,
        "stage13_external_search_performed": False,
        "gpu_prediction_performed": False,
        "stage10_scores_recalculated": False,
        "final_top5_top3_selected": False,
    }

    with (
        out / "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")

    print(
        "STAGE12_CONSERVATION_COMPLETED_FAMILIES="
        f"{len(eligible)}"
    )
    print(
        "STAGE12_CONSERVATION_NOT_EVALUATED_FAMILIES="
        f"{len(state['rows']) - len(eligible)}"
    )
    print("STAGE12_ADAPTIVE_CONSERVATION=PASS")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--stage09-10-output",
        type=Path,
        required=True,
    )
    ap.add_argument(
        "--adaptive-output",
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
            "MSA_ANALYSIS_ELIGIBLE_FAMILIES="
            f"{len(state['eligible'])}"
        )
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
