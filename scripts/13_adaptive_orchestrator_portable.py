#!/usr/bin/env python3
"""Stage 13 adaptive external evidence and ranking handoff.

Without a verified MSA external search, preserve PENDING_SEARCH.
With verified legacy inputs, regenerate the 25/27/33-column
evidence tables and explicitly label the result as verified reuse.

No GPU prediction, Foldseek search or final Top5/Top3 selection.
"""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(path):
    path = Path(path).resolve()
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"MISSING_FILE:{path}")
    return path


def meta(path):
    obj = json.loads(require(path).read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return obj


def table(path):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not fields
        or len(fields) != len(set(fields))
        or not rows
        or "structure_id" not in fields
        or any(
            None in row
            or any(v is None for v in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    ids = [r["structure_id"] for r in rows]

    if (
        any(not sid for sid in ids)
        or len(ids) != len(set(ids))
    ):
        raise ValueError(f"INVALID_IDS:{path}")

    return fields, rows


def execute(label, name, options, log=None):
    command = [
        sys.executable, "-B",
        str(SCRIPTS / name),
        *map(str, options),
    ]

    result = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("x", encoding="utf-8") as stream:
            stream.write(result.stdout)
            stream.write(result.stderr)

    if result.returncode:
        raise RuntimeError(
            f"{label}:FAILED\n"
            + (result.stderr or result.stdout)[-5000:]
        )

    print(f"[PASS] {label}", flush=True)


def validate(args):
    out = args.out_dir.resolve()
    down = args.stage09_10_output.resolve()
    adaptive = args.adaptive_output.resolve()
    conservation = args.conservation_output.resolve()
    spatial = args.spatial_output.resolve()

    if out.exists():
        raise ValueError("Refusing existing output directory")

    for directory in (down, adaptive, conservation, spatial):
        if not directory.is_dir():
            raise ValueError(f"MISSING_DIRECTORY:{directory}")

        if (
            out == directory
            or out in directory.parents
            or directory in out.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    rank = require(
        down / "10_ranking/confidence_aware_ranking.tsv"
    )
    rank_prov = require(
        down / "10_ranking/ranking_provenance.json"
    )
    shortlist = require(
        down / "10_shortlists/top15_experimental_candidates.tsv"
    )

    adaptive_manifest = require(
        adaptive / "orchestration_provenance.json"
    )
    conservation_manifest = require(
        conservation / "orchestration_provenance.json"
    )
    spatial_manifest = require(
        spatial / "orchestration_provenance.json"
    )

    conservation_status = require(
        conservation / "stage12_evaluation_status.tsv"
    )
    spatial_status = require(
        spatial / "stage12_spatial_evaluation_status.tsv"
    )

    am = meta(adaptive_manifest)
    cm = meta(conservation_manifest)
    sm = meta(spatial_manifest)

    if (
        am.get("stage10_ranking_sha256") != sha(rank)
        or cm.get("stage10_ranking_sha256") != sha(rank)
        or cm.get("stage11_manifest_sha256")
        != sha(adaptive_manifest)
        or sm.get("stage12_conservation_manifest_sha256")
        != sha(conservation_manifest)
        or sm.get("stage12_conservation_status_sha256")
        != sha(conservation_status)
        or sm.get("stage12_spatial_status_sha256")
        != sha(spatial_status)
    ):
        raise ValueError("ADAPTIVE_STAGE11_12_PROVENANCE_MISMATCH")

    fields, rows = table(spatial_status)
    rank_fields, rank_rows = table(rank)
    _, selected = table(shortlist)

    if (
        len(fields) != 36
        or len(rank_fields) != 30
        or fields[:30] != rank_fields
        or len(rows) != len(selected)
    ):
        raise ValueError("RANKING_STATUS_SCHEMA_MISMATCH")

    selected_index = {
        row["structure_id"]: row for row in selected
    }

    if (
        len(selected_index) != len(selected)
        or set(selected_index) != {
            row["structure_id"] for row in rows
        }
    ):
        raise ValueError("SELECTED_FAMILY_COVERAGE_MISMATCH")

    eligible = []

    for row in rows:
        sid = row["structure_id"]

        if any(
            row[field] != selected_index[sid][field]
            for field in rank_fields
        ):
            raise ValueError(
                f"STAGE10_VALUES_CHANGED:{sid}"
            )

        msa_status = row["stage11_msa_qc_status"]
        conservation_status_value = row[
            "stage12_conservation_status"
        ]
        spatial_status_value = row[
            "stage12_spatial_status"
        ]
        external_status = row[
            "stage13_msa_external_status"
        ]

        if msa_status == "NOT_EVALUATED":
            if (
                conservation_status_value != "NOT_EVALUATED"
                or spatial_status_value != "NOT_EVALUATED"
                or external_status != "NOT_EVALUATED"
            ):
                raise ValueError(
                    f"INVALID_SINGLE_SEQUENCE_FALLBACK:{sid}"
                )

        elif msa_status == "QC_VERIFIED":
            if (
                conservation_status_value != "COMPLETED"
                or spatial_status_value != "COMPLETED_EXPLORATORY"
                or external_status != "PENDING_SEARCH"
            ):
                raise ValueError(
                    f"INVALID_STAGE13_READY_STATUS:{sid}"
                )
            eligible.append(sid)

        else:
            raise ValueError(
                f"UNSUPPORTED_MSA_STATUS:{sid}"
            )

    has_external = args.external_dir is not None

    external_args = [
        args.stage08_output,
        args.external_dir,
        args.k2_structures,
        args.legacy_qc_dir,
        args.legacy_conservation_dir,
        args.legacy_spatial_dir,
        args.selection_dir,
    ]

    if has_external:
        if not all(value is not None for value in external_args):
            raise ValueError(
                "VERIFIED_EXTERNAL_INPUTS_INCOMPLETE"
            )

        if (
            not eligible
            or len(eligible) != len(rows)
        ):
            raise ValueError(
                "VERIFIED_EXTERNAL_REUSE_REQUIRES_ALL_SELECTED_MSA"
            )

        current_qc = require(
            adaptive /
            "11_collect/msa_qc/msa_structure_qc.tsv"
        )
        legacy_qc = require(
            args.legacy_qc_dir /
            "msa_structure_qc.tsv"
        )

        if sha(current_qc) != sha(legacy_qc):
            raise ValueError(
                "CURRENT_LEGACY_MSA_QC_MISMATCH"
            )

        selection = require(
            args.selection_dir /
            "stage10_selected_families.tsv"
        )
        _, selected_rows = table(selection)

        if {
            r["structure_id"] for r in selected_rows
        } != set(eligible):
            raise ValueError(
                "EXTERNAL_SELECTION_FAMILY_MISMATCH"
            )

        # This established workflow validates full Stage08 family
        # membership, GPU PDB, conservation, spatial and K2 sources.
        check_arguments = [
            "--stage08-output", args.stage08_output,
            "--stage09-10-output", down,
            "--msa-qc-dir", args.legacy_qc_dir,
            "--conservation-dir",
            args.legacy_conservation_dir,
            "--spatial-dir", args.legacy_spatial_dir,
            "--k2-structures", args.k2_structures,
            "--external-dir", args.external_dir,
            "--out-dir", out / "legacy_check_only",
            "--mode", "check",
        ]

        execute(
            "Verified MSA external provenance preflight",
            "11_13_orchestrator_portable.py",
            check_arguments,
        )

        if (out / "legacy_check_only").exists():
            raise ValueError("CHECK_MODE_WROTE_FILES")

    elif any(value is not None for value in external_args):
        raise ValueError(
            "EXTERNAL_OPTIONS_REQUIRE_EXTERNAL_DIR"
        )

    return {
        "out": out,
        "down": down,
        "adaptive": adaptive,
        "conservation": conservation,
        "spatial": spatial,
        "rank": rank,
        "rank_prov": rank_prov,
        "shortlist": shortlist,
        "spatial_status": spatial_status,
        "spatial_manifest": spatial_manifest,
        "fields": fields,
        "rows": rows,
        "eligible": eligible,
        "external": has_external,
        "args": args,
    }


def run(state):
    out = state["out"]
    out.mkdir(parents=True, exist_ok=False)

    rows = [dict(row) for row in state["rows"]]
    evidence = None
    preliminary = None
    merged = None

    if state["external"]:
        args = state["args"]
        selected = args.selection_dir
        qc = (
            state["adaptive"] /
            "11_collect/msa_qc"
        )
        cons = (
            state["conservation"] /
            "conservation"
        )
        spatial = (
            state["spatial"] /
            "spatial"
        )

        preliminary_dir = out / "13_preliminary"

        execute(
            "13 new preliminary evidence",
            "13_candidate_evidence_preliminary_portable.py",
            [
                "--selection",
                selected / "stage10_selected_families.tsv",
                "--msa-qc", qc / "msa_structure_qc.tsv",
                "--msa-qc-provenance",
                qc / "msa_structure_qc_provenance.json",
                "--conservation-dir", cons,
                "--spatial-dir", spatial,
                "--out-dir", preliminary_dir,
            ],
            out / "logs/13_preliminary.log",
        )

        preliminary = require(
            preliminary_dir /
            "preliminary_candidate_evidence.tsv"
        )
        preliminary_prov = require(
            preliminary_dir /
            "preliminary_evidence_provenance.json"
        )

        merged_dir = out / "13_external"

        execute(
            "13 external evidence merge",
            "13_merge_external_evidence_portable.py",
            [
                "--preliminary", preliminary,
                "--preliminary-provenance",
                preliminary_prov,
                "--external",
                args.external_dir /
                "msa_external_evidence_stage13.tsv",
                "--external-provenance",
                args.external_dir /
                "msa_external_evidence_provenance.json",
                "--out-dir", merged_dir,
            ],
            out / "logs/13_merge.log",
        )

        merged = require(
            merged_dir /
            "preliminary_candidate_evidence_with_external.tsv"
        )
        merged_prov = require(
            merged_dir /
            "preliminary_external_merge_provenance.json"
        )

        ranked_dir = out / "13_ranked"

        execute(
            "13 attach preserved Stage10 ranking",
            "13_attach_stage10_ranking_portable.py",
            [
                "--ranking", state["rank"],
                "--ranking-provenance",
                state["rank_prov"],
                "--shortlist", state["shortlist"],
                "--selection",
                selected / "stage10_selected_families.tsv",
                "--selection-provenance",
                selected / "selection_provenance.json",
                "--preliminary", preliminary,
                "--preliminary-provenance",
                preliminary_prov,
                "--merged", merged,
                "--merged-provenance", merged_prov,
                "--out-dir", ranked_dir,
            ],
            out / "logs/13_ranked.log",
        )

        evidence = require(
            ranked_dir / "ranked_candidate_evidence.tsv"
        )
        ef, er = table(evidence)

        if (
            len(ef) != 33
            or len(er) != len(state["eligible"])
            or {
                row["structure_id"] for row in er
            } != set(state["eligible"])
        ):
            raise ValueError(
                "STAGE13_FINAL_EVIDENCE_FORMAT_MISMATCH"
            )

        final_prov = meta(
            ranked_dir /
            "ranked_candidate_evidence_provenance.json"
        )

        if final_prov.get("output_sha256") != sha(evidence):
            raise ValueError(
                "STAGE13_FINAL_EVIDENCE_SHA256_MISMATCH"
            )

        for row in rows:
            if row["structure_id"] in state["eligible"]:
                row["stage13_msa_external_status"] = (
                    "COMPLETED_VERIFIED_REUSE"
                )

        workflow = "STAGE13_VERIFIED_EXTERNAL_REUSE"

    else:
        workflow = (
            "AWAITING_MSA_EXTERNAL_SEARCH"
            if state["eligible"]
            else "STAGE13_NOT_EVALUATED_SINGLE_UNIQUE"
        )

    status_path = out / "stage13_evaluation_status.tsv"

    if not state["external"]:
        shutil.copyfile(
            state["spatial_status"], status_path
        )
    else:
        with status_path.open(
            "x", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=state["fields"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    provenance = {
        "stage": "STAGE13_ADAPTIVE",
        "workflow_state": workflow,
        "stage10_ranking_sha256":
            sha(state["rank"]),
        "stage12_spatial_manifest_sha256":
            sha(state["spatial_manifest"]),
        "stage12_spatial_status_sha256":
            sha(state["spatial_status"]),
        "stage13_status_sha256":
            sha(status_path),
        "preliminary_evidence_sha256":
            sha(preliminary) if preliminary else None,
        "merged_external_evidence_sha256":
            sha(merged) if merged else None,
        "final_33col_evidence_sha256":
            sha(evidence) if evidence else None,
        "selected_families": len(rows),
        "msa_external_eligible_families":
            len(state["eligible"]),
        "verified_external_families":
            len(state["eligible"])
            if evidence else 0,
        "external_mode":
            "VERIFIED_REUSE" if evidence else "NOT_EXECUTED",
        "new_foldseek_search_performed": False,
        "gpu_prediction_performed": False,
        "stage10_scores_recalculated": False,
        "final_top5_top3_selected": False,
        "biochemical_activity_validated": False,
    }

    with (
        out / "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as stream:
        json.dump(
            provenance, stream, indent=2
        )
        stream.write("\n")

    print(f"STAGE13_WORKFLOW_STATE={workflow}")
    print(
        "STAGE13_VERIFIED_EXTERNAL_FAMILIES="
        f"{provenance['verified_external_families']}"
    )
    print("STAGE13_ADAPTIVE_ORCHESTRATION=PASS")


def main():
    parser = argparse.ArgumentParser()

    for name in (
        "stage09-10-output",
        "adaptive-output",
        "conservation-output",
        "spatial-output",
        "out-dir",
    ):
        parser.add_argument(
            "--" + name,
            type=Path,
            required=True,
        )

    for name in (
        "stage08-output",
        "external-dir",
        "k2-structures",
        "legacy-qc-dir",
        "legacy-conservation-dir",
        "legacy-spatial-dir",
        "selection-dir",
    ):
        parser.add_argument(
            "--" + name,
            type=Path,
        )

    parser.add_argument(
        "--mode",
        choices=("check", "run"),
        default="check",
    )

    args = parser.parse_args()
    state = validate(args)

    if args.mode == "check":
        print(
            "STAGE13_ELIGIBLE_FAMILIES="
            f"{len(state['eligible'])}"
        )
        print(
            "VERIFIED_EXTERNAL_INPUT="
            f"{state['external']}"
        )
        print("OUTPUT_CREATED=NO")
        print("STAGE13_ADAPTIVE_CHECK=PASS")
    else:
        run(state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, TypeError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
