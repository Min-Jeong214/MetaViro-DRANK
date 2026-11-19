#!/usr/bin/env python3
"""Record reproducible pre-run provenance."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

SCRIPT_NAMES = [
    "run_pipeline.py",
    "scripts/00_validate_input_portable.py",
    "scripts/00_capture_run_provenance.py",
    "scripts/01_define_unresolved_candidates.py",
    "scripts/02_extract_replication_context.py",
    "scripts/03_global_family_clustering.sh",
    "scripts/04_select_recurrent_families_portable.py",
    "scripts/05_extract_recurrent_members.py",
    "scripts/05_annotation_saturation_portable.sh",
    "scripts/05_build_member_evidence.py",
    "scripts/06_family_darkness_audit_portable.py",
    "scripts/07_familywide_context_portable.py",
    "scripts/08_prepare_structure_queries_portable.py",
    "scripts/08_collect_structure_results_portable.py",
    "scripts/08_structure_prediction.sh",
]

INPUT_NAMES = [
    "proteins.faa",
    "genes.gff",
    "annotations.emapper.annotations",
]


def resolve(value):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else ROOT / path).resolve()


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(4 * 1024 * 1024), b""
        ):
            digest.update(block)

    return digest.hexdigest()


def fingerprint(path):
    path = Path(path)

    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty file: {path}")

    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def command_output(command):
    result = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--through",
        choices=["02", "03", "04", "05", "06", "07", "08"],
        default="06",
    )
    parser.add_argument(
        "--structure-mode",
        choices=["prepare", "collect", "predict"],
        default="prepare",
    )
    parser.add_argument(
        "--structure-env",
        default="metaviro-structure",
    )
    parser.add_argument(
        "--structure-prediction-dir",
        default="",
    )
    args = parser.parse_args()

    config_path = resolve(args.config)
    output_path = resolve(args.out)

    if output_path.exists():
        raise ValueError(
            "Refusing to overwrite existing provenance: "
            + str(output_path)
        )

    cfg = yaml.safe_load(config_path.read_text())

    if not isinstance(cfg, dict):
        raise ValueError("Invalid YAML configuration")

    input_dir = resolve(cfg["input_dir"])
    samples = sorted(
        p for p in input_dir.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )

    if not samples:
        raise ValueError("No input samples")

    inputs = {}

    for sample in samples:
        inputs[sample.name] = {
            name: fingerprint(sample / name)
            for name in INPUT_NAMES
        }

    scripts = {
        name: fingerprint(ROOT / name)
        for name in SCRIPT_NAMES
    }

    databases = {}
    db_cfg = cfg.get("databases", {})

    if args.through in ("05", "06", "07", "08"):
        for name in (
            "swissprot",
            "phrog",
            "pfam",
            "phrog_annotations",
        ):
            databases[name] = fingerprint(
                resolve(db_cfg[name])
            )

        phrog_base = resolve(db_cfg["phrog"])
        pfam_base = resolve(db_cfg["pfam"])

        for suffix in (
            ".dbtype",
            ".index",
            "_h",
            "_h.dbtype",
            "_h.index",
        ):
            databases["phrog" + suffix] = fingerprint(
                Path(str(phrog_base) + suffix)
            )

        for suffix in (
            ".h3f",
            ".h3i",
            ".h3m",
            ".h3p",
        ):
            databases["pfam" + suffix] = fingerprint(
                Path(str(pfam_base) + suffix)
            )

    # Executable prerequisites depend on the chosen stage.  Stage 00-02
    # must remain runnable on a clean Python + PyYAML installation; MMseqs2
    # begins at Stage 03, while HMMER is first used at Stage 05.
    if args.through in ("03", "04", "05", "06", "07", "08"):
        mmseqs_version = command_output(["mmseqs", "version"]).splitlines()[0]
    else:
        mmseqs_version = "NOT_REQUIRED_FOR_STAGE_02"

    if args.through in ("05", "06", "07", "08"):
        hmmer_help = command_output(["hmmscan", "-h"])
        hmmer_version = next(
            (line.strip() for line in hmmer_help.splitlines() if "HMMER" in line),
            "VERSION_NOT_IDENTIFIED",
        )
    else:
        hmmer_version = "NOT_REQUIRED_FOR_THIS_STAGE"

    # Git metadata are recorded for cloned checkouts, but a release source
    # archive is also a supported installation path. Absence is explicit.
    try:
        git_head = command_output(["git", "rev-parse", "HEAD"])
        git_dirty = bool(command_output([
            "git", "status", "--porcelain", "--untracked-files=no",
        ]))
    except (FileNotFoundError, subprocess.CalledProcessError):
        git_head = None
        git_dirty = None

    max_seqs = os.environ.get(
        "PHROG_MAX_SEQS", "1000"
    )

    if not max_seqs.isdigit() or int(max_seqs) < 1:
        raise ValueError("Invalid PHROG_MAX_SEQS")

    record = {
        "schema_version": 1,
        "record_type": "PRE_RUN_PROVENANCE",
        "analysis_status": "NOT_STARTED",
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "intended_through": args.through,
        "config": fingerprint(config_path),
        "parameters": {
            "threads": int(cfg.get("threads", 8)),
            "min_aa": int(cfg["min_aa"]),
            "window": int(cfg["window"]),
            "phrog_max_seqs": int(max_seqs),
            "structure_mode": (
                args.structure_mode
                if args.through == "08" else None
            ),
            "structure_env": (
                args.structure_env
                if args.through == "08" else None
            ),
            "structure_prediction_dir": (
                args.structure_prediction_dir or None
            ),
        },
        "software": {
            "python": sys.version.split()[0],
            "mmseqs2": mmseqs_version,
            "hmmer": hmmer_version,
        },
        "git": {
            "head": git_head,
            "tracked_files_modified": git_dirty,
        },
        "samples": inputs,
        "databases": databases,
        "pipeline_scripts": scripts,
    }

    output_path.parent.mkdir(
        parents=True, exist_ok=True
    )

    fd, temporary = tempfile.mkstemp(
        prefix=".run_provenance.",
        suffix=".tmp",
        dir=output_path.parent,
    )

    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(
                record, fh,
                indent=2,
                sort_keys=True,
            )
            fh.write("\n")

        os.replace(temporary, output_path)

    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

    print("RUN_PROVENANCE=PASS")
    print(f"RECORDED_SAMPLES={len(samples)}")
    print(f"RECORDED_DATABASE_FILES={len(databases)}")
    print(f"RECORDED_SCRIPTS={len(scripts)}")
    print("ANALYSIS_STATUS=NOT_STARTED")
    print(f"PROVENANCE_FILE={output_path}")


if __name__ == "__main__":
    try:
        main()
    except (
        ValueError,
        OSError,
        KeyError,
        subprocess.CalledProcessError,
    ) as exc:
        sys.exit(f"RUN_PROVENANCE=FAIL: {exc}")
