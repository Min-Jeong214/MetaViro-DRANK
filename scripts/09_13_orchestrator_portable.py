#!/usr/bin/env python3
"""Unified Stage 09-13 portable orchestrator.

Stage 09-10 supports Foldseek search or verified reuse.
Stage 11-13 currently requires independently verified reuse.
Does not perform GPU prediction or final Top5/Top3 selection.
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
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_OR_EMPTY_FILE:{path}")
    return path


def read_json(path):
    value = json.loads(require(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"INVALID_JSON:{path}")
    return value


def read_tsv(path):
    with require(path).open(
        newline="", encoding="utf-8"
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)

    if (
        not fields
        or len(fields) != len(set(fields))
        or "structure_id" not in fields
        or not rows
    ):
        raise ValueError(f"INVALID_TABLE:{path}")

    ids = [row["structure_id"] for row in rows]

    if (
        len(ids) != len(set(ids))
        or any(not sid for sid in ids)
        or any(
            None in row
            or any(v is None for v in row.values())
            for row in rows
        )
    ):
        raise ValueError(f"MALFORMED_TABLE:{path}")

    return fields, rows


def invoke(label, script, arguments, log=None):
    command = [
        sys.executable, "-B",
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


def prepare(args):
    out = args.out_dir.resolve()
    up = args.stage08_output.resolve()

    if out.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    inputs = [
        up,
        args.msa_qc_dir.resolve(),
        args.conservation_dir.resolve(),
        args.spatial_dir.resolve(),
        args.k2_structures.resolve(),
        args.external_dir.resolve(),
    ]

    if args.reuse_search_root:
        inputs.append(
            args.reuse_search_root.resolve()
        )
    else:
        inputs.append(args.db_root.resolve())

    for source in inputs:
        if (
            out == source
            or source in out.parents
            or out in source.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    if not up.is_dir():
        raise ValueError("STAGE08_DIRECTORY_MISSING")

    candidates = sorted(
        (up / "08_structure").glob(
            "run_*/collected/structure_results.tsv"
        )
    )

    if args.structure_results is not None:
        structures = (
            args.structure_results.resolve()
        )

        if structures not in [
            path.resolve()
            for path in candidates
        ]:
            raise ValueError(
                "INVALID_SELECTED_STRUCTURE_TABLE"
            )
    else:
        if len(candidates) != 1:
            raise ValueError(
                "STAGE08_STRUCTURE_RUN_NOT_UNIQUE"
            )

        structures = candidates[0].resolve()

    _, stage08_rows = read_tsv(structures)

    qcdir = args.msa_qc_dir.resolve()
    consdir = args.conservation_dir.resolve()
    spatdir = args.spatial_dir.resolve()
    extdir = args.external_dir.resolve()

    qc_file = require(
        qcdir / "msa_structure_qc.tsv"
    )
    qc_meta = read_json(
        qcdir / "msa_structure_qc_provenance.json"
    )

    if (
        qc_meta.get("source_stage08_sha256")
        != sha(structures)
        or qc_meta.get("msa_qc_sha256")
        != sha(qc_file)
    ):
        raise ValueError(
            "STAGE08_11_PROVENANCE_MISMATCH"
        )

    cons_file = require(
        consdir / "conservation_provenance.json"
    )
    cons_meta = read_json(cons_file)

    if (
        cons_meta.get("prediction_index_sha256")
        != qc_meta.get("prediction_index_sha256")
    ):
        raise ValueError(
            "STAGE11_12_PREDICTION_INDEX_MISMATCH"
        )

    spat_meta = read_json(
        spatdir / "spatial_cluster_provenance.json"
    )

    if (
        spat_meta.get("stage12_provenance_sha256")
        != sha(cons_file)
    ):
        raise ValueError(
            "STAGE12_SPATIAL_PROVENANCE_MISMATCH"
        )

    k2_file = require(
        args.k2_structures.resolve()
    )

    _, k2_rows = read_tsv(k2_file)
    _, qc_rows = read_tsv(qc_file)

    ext_file = require(
        extdir / "msa_external_evidence_stage13.tsv"
    )
    ext_meta = read_json(
        extdir / "msa_external_evidence_provenance.json"
    )
    _, ext_rows = read_tsv(ext_file)

    if (
        ext_meta.get("stage13_evidence_sha256")
        != sha(ext_file)
        or ext_meta.get("control_only") is not False
        or ext_meta.get(
            "final_selection_performed"
        ) is not False
    ):
        raise ValueError(
            "K2_EXTERNAL_EVIDENCE_SHA256_MISMATCH"
        )

    stage08 = {
        row["structure_id"]: row
        for row in stage08_rows
    }
    qc = {
        row["structure_id"]: row
        for row in qc_rows
    }
    k2 = {
        row["structure_id"]: row
        for row in k2_rows
    }
    ext = {
        row["structure_id"]: row
        for row in ext_rows
    }

    if set(k2) != set(ext):
        raise ValueError(
            "K2_EXTERNAL_STRUCTURE_ID_MISMATCH"
        )

    for sid, msa in k2.items():
        if sid not in stage08 or sid not in qc:
            raise ValueError(
                f"K2_STRUCTURE_NOT_IN_UPSTREAM:{sid}"
            )

        single = stage08[sid]
        checked = qc[sid]
        external = ext[sid]

        if (
            single["family_rep"]
            != msa["family_rep"]
            or single["sequence_sha256"]
            != msa["sequence_sha256"]
            or checked["single_pdb_sha256"]
            != single["pdb_sha256"]
            or checked["msa_pdb_sha256"]
            != msa["pdb_sha256"]
            or external["msa_pdb_sha256"]
            != msa["pdb_sha256"]
        ):
            raise ValueError(
                f"CROSS_STAGE_STRUCTURE_SHA256_MISMATCH:{sid}"
            )

    print("STAGE08_TO_STAGE11_PROVENANCE=PASS")
    print("STAGE11_TO_STAGE12_PROVENANCE=PASS")
    print("STAGE12_SPATIAL_PROVENANCE=PASS")
    print("STAGE_K2_EXTERNAL_PROVENANCE=PASS")

    return out


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage08-output",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--structure-results",
        type=Path,
    )

    group = parser.add_mutually_exclusive_group(
        required=True
    )

    group.add_argument(
        "--reuse-search-root",
        type=Path,
    )

    group.add_argument(
        "--db-root",
        type=Path,
    )

    for name in (
        "msa_qc_dir",
        "conservation_dir",
        "spatial_dir",
        "k2_structures",
        "external_dir",
    ):
        parser.add_argument(
            "--" + name.replace("_", "-"),
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
    # MSA_MEMBERSHIP_GUARD_STEP119
    import subprocess as _guard_subprocess
    import sys as _guard_sys
    from pathlib import Path as _guard_Path
    _guard_command = [
        _guard_sys.executable, "-B",
        str(_guard_Path(__file__).with_name(
            "11_msa_reuse_membership_guard_portable.py"
        )),
        "--stage08-output", str(args.stage08_output),
        "--msa-qc-dir", str(args.msa_qc_dir),
    ]
    _guard_result = _guard_subprocess.run(
        _guard_command,
        capture_output=True, text=True, check=False,
    )
    if _guard_result.returncode != 0:
        raise ValueError(
            "MSA_MEMBERSHIP_REUSE_REJECTED: "
            + (_guard_result.stderr or _guard_result.stdout)
        )
    print("MSA_MEMBERSHIP_REUSE_GUARD=PASS")

    out = prepare(args)

    first = [
        "--stage08-output",
        args.stage08_output,
        "--out-dir",
        out / "09_10",
        "--threads",
        args.threads,
    ]

    if args.structure_results:
        first.extend([
            "--structure-results",
            args.structure_results,
        ])

    if args.reuse_search_root:
        first.extend([
            "--reuse-search-root",
            args.reuse_search_root,
        ])
    else:
        first.extend([
            "--db-root",
            args.db_root,
        ])

    invoke(
        "09-10 input check",
        "09_10_orchestrator_portable.py",
        first + ["--mode", "check"],
    )

    if args.mode == "check":
        print(
            "STAGE11_13_SELECTION_DEPENDENT_CHECK="
            "DEFERRED_UNTIL_RANKING"
        )
        print("UNIFIED_CHECK_MODE=PASS")
        print("OUTPUT_CREATED=NO")
        return

    out.mkdir(
        parents=True,
        exist_ok=False,
    )

    invoke(
        "09-10 automated analysis",
        "09_10_orchestrator_portable.py",
        first + ["--mode", "run"],
        out / "logs/09_10.log",
    )

    second = [
        "--stage08-output",
        args.stage08_output,
        "--stage09-10-output",
        out / "09_10",
        "--msa-qc-dir",
        args.msa_qc_dir,
        "--conservation-dir",
        args.conservation_dir,
        "--spatial-dir",
        args.spatial_dir,
        "--k2-structures",
        args.k2_structures,
        "--external-dir",
        args.external_dir,
        "--out-dir",
        out / "11_13",
        "--mode",
        "run",
    ]

    invoke(
        "11-13 verified reuse analysis",
        "11_13_orchestrator_portable.py",
        second,
        out / "logs/11_13.log",
    )

    rank_file = require(
        out /
        "09_10/10_ranking/"
        "confidence_aware_ranking.tsv"
    )
    top15_file = require(
        out /
        "09_10/10_shortlists/"
        "top15_experimental_candidates.tsv"
    )
    final_file = require(
        out /
        "11_13/13_ranked/"
        "ranked_candidate_evidence.tsv"
    )

    m09 = read_json(
        out /
        "09_10/orchestration_provenance.json"
    )
    m13 = read_json(
        out /
        "11_13/orchestration_provenance.json"
    )

    if (
        m13.get("source_stage10_ranking_sha256")
        != sha(rank_file)
        or m13.get(
            "source_stage10_shortlist_sha256"
        ) != sha(top15_file)
    ):
        raise ValueError(
            "ORCHESTRATOR_CHAIN_SHA256_MISMATCH"
        )

    rank_fields, rank_rows = read_tsv(
        rank_file
    )
    _, top15_rows = read_tsv(top15_file)
    final_fields, final_rows = read_tsv(
        final_file
    )

    selected = {
        row["structure_id"]
        for row in top15_rows
    }

    if (
        selected
        != {
            row["structure_id"]
            for row in final_rows
        }
        or len(rank_fields) != 30
        or len(final_fields) != 33
    ):
        raise ValueError(
            "ORCHESTRATOR_FINAL_TABLE_MISMATCH"
        )

    if (
        m09.get("family_count")
        != len(rank_rows)
        or m13.get("result_families")
        != len(final_rows)
        or m13.get(
            "final_candidate_selection_performed"
        ) is not False
    ):
        raise ValueError(
            "ORCHESTRATOR_MANIFEST_MISMATCH"
        )

    manifest = {
        "stage": "STAGE09_13_UNIFIED_ORCHESTRATION",
        "mode": "STAGE11_13_VERIFIED_REUSE",
        "stage09_10_manifest_sha256": sha(
            out /
            "09_10/orchestration_provenance.json"
        ),
        "stage11_13_manifest_sha256": sha(
            out /
            "11_13/orchestration_provenance.json"
        ),
        "ranking_sha256": sha(rank_file),
        "top15_sha256": sha(top15_file),
        "stage13_evidence_sha256": sha(final_file),
        "ranking_families": len(rank_rows),
        "selected_families": len(final_rows),
        "stage13_columns": len(final_fields),
        "gpu_prediction_performed": False,
        "final_top5_top3_selected": False,
    }

    with (
        out / "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    print(
        f"RANKED_FAMILIES={len(rank_rows)}"
    )
    print(
        f"STAGE13_SELECTED_FAMILIES={len(final_rows)}"
    )
    print(
        f"STAGE13_COLUMNS={len(final_fields)}"
    )
    print(
        f"UNIFIED_OUTPUT={out}"
    )
    print("STAGE09_TO_STAGE13_UNIFIED_ORCHESTRATION=PASS")


if __name__ == "__main__":
    try:
        main()
    except (
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        csv.Error,
    ) as exc:
        sys.exit(f"ERROR: {exc}")
