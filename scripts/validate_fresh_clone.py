#!/usr/bin/env python3
"""Repeatable fresh-checkout verification; no research samples are bundled.

Default smoke: independently generated 1- and 3-sample cohorts, real Stage
00-02 execution, result/provenance checks, and a deliberate invalid GFF case.
Database mode: preflight (or optionally execute) Stage 00-07 plus CPU context
screen on caller-owned inputs and database paths.  This does not assert that
optional structure searches or biological functions have been validated.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def call(*argv: str | Path, expect_ok: bool = True) -> str:
    args = [str(arg) for arg in argv]
    p = subprocess.run(args, cwd=ROOT, text=True, capture_output=True, check=False)
    if (p.returncode == 0) != expect_ok:
        raise RuntimeError(
            f"WRONG_EXIT={p.returncode}: {' '.join(args)}\n"
            f"STDOUT:\n{p.stdout[-5000:]}\nSTDERR:\n{p.stderr[-5000:]}"
        )
    return p.stdout + p.stderr


def tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"MISSING_OUTPUT:{path}")
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def verify_cohort(directory: Path, count: int, prefix: str) -> dict:
    fixture = directory / f"fixture_{count}"
    created = call(sys.executable, "-B", ROOT / "examples/generate_demo_input.py",
                   "--out-dir", fixture, "--sample-count", str(count),
                   "--prefix", prefix)
    if "DEMO_GENERATED=PASS" not in created:
        raise RuntimeError("FIXTURE_GENERATION_FAILED")
    config = fixture / "config.yaml"
    check = call(sys.executable, "-B", ROOT / "run_pipeline.py",
                 "--config", config, "--through", "02", "--check")
    if "INPUT_PREFLIGHT=PASS" not in check or f"SAMPLES={count}" not in check:
        raise RuntimeError("PREFLIGHT_FAILED")
    executed = call(sys.executable, "-B", ROOT / "run_pipeline.py",
                    "--config", config, "--through", "02")
    if "STAGES_01_02=COMPLETE" not in executed:
        raise RuntimeError("STAGE02_COMPLETION_MISSING")
    output = fixture / "analysis_output"
    events = list((output / "provenance").glob("run_*.json"))
    if len(events) != 1:
        raise RuntimeError("INVALID_PROVENANCE_COUNT")
    meta = json.loads(events[0].read_text())
    if (len(meta["samples"]) != count
            or meta["software"]["mmseqs2"] != "NOT_REQUIRED_FOR_STAGE_02"
            or meta["software"]["hmmer"] != "NOT_REQUIRED_FOR_THIS_STAGE"):
        raise RuntimeError("INVALID_PROVENANCE_CONTENT")
    names = sorted(p.name for p in (fixture / "input").iterdir())
    ids = set()
    for name in names:
        unresolved = tsv_rows(output / "01_unresolved" / f"{name}.tsv")
        context = tsv_rows(output / "02_context" / f"{name}.tsv")
        anchors = tsv_rows(output / "02_context" / f"{name}_anchors.tsv")
        if (len(unresolved) != 3 or len(context) != 2 or len(anchors) != 1
                or sum(row["strict_unresolved"] == "1" for row in unresolved) != 2):
            raise RuntimeError(f"WRONG_STAGE02_COUNTS:{name}")
        local = {row["candidate_id"] for row in context}
        if ids & local:
            raise RuntimeError("CROSS_SAMPLE_ID_COLLISION")
        ids.update(local)
    print(f"FRESH_CLONE_STAGE00_02_{count}_COHORTS=PASS", flush=True)
    return {"samples": count, "status": "PASS", "protein_count": 3 * count,
            "context_candidate_count": 2 * count}


def test_invalid_input(directory: Path) -> None:
    fixture = directory / "invalid"
    call(sys.executable, "-B", ROOT / "examples/generate_demo_input.py",
         "--out-dir", fixture, "--sample-count", "2", "--prefix", "BAD")
    gff = fixture / "input" / "BAD_A" / "genes.gff"
    gff.write_text(gff.read_text().replace("ID=1_1", "ID=1_999"))
    output = call(sys.executable, "-B", ROOT / "run_pipeline.py",
                  "--config", fixture / "config.yaml", "--through", "02",
                  "--check", expect_ok=False)
    if "Stage 00 validation failed" not in output:
        raise RuntimeError("GFF_ID_MISMATCH_NOT_REJECTED")
    print("FRESH_CLONE_INVALID_GFF_REJECTED=PASS", flush=True)


def smoke() -> dict:
    for runner in ("run_metaviro.py", "run_pipeline.py", "run_pipeline_full.py"):
        call(sys.executable, "-B", ROOT / runner, "--help")
    with tempfile.TemporaryDirectory(prefix="metaviro_fresh_clone_") as directory:
        root = Path(directory)
        results = [verify_cohort(root, 1, "ONE"), verify_cohort(root, 3, "THREE")]
        test_invalid_input(root)
    print("FRESH_CLONE_SMOKE=PASS", flush=True)
    return {"mode": "smoke", "status": "PASS", "cohorts": results,
            "invalid_gff_rejection": "PASS"}


def database(args: argparse.Namespace) -> dict:
    if not args.config or not args.screen_out_dir:
        raise ValueError("DATABASE_MODE_REQUIRES_CONFIG_AND_SCREEN_OUT_DIR")
    config = args.config.resolve()
    output = args.screen_out_dir.resolve()
    argv = [sys.executable, "-B", ROOT / "run_metaviro.py", "run",
            "--config", config, "--screen-out-dir", output]
    if args.normalize_ids:
        argv.append("--normalize-ids")
    if args.execute and output.exists():
        raise ValueError(f"OUTPUT_ALREADY_EXISTS:{output}")
    checked = call(*argv, "--mode", "check")
    if "STEP181_NEW_COHORT_PREFLIGHT=PASS" not in checked:
        raise RuntimeError("STAGE07_PREFLIGHT_NOT_CONFIRMED")
    print("STAGE00_07_DATABASE_PREFLIGHT=PASS", flush=True)
    if not args.execute:
        return {"mode": "database", "status": "PREFLIGHT_PASS", "executed": False}
    actual = call(*argv, "--mode", "run")
    if "STEP181_MAIN_CLI=PASS" not in actual:
        raise RuntimeError("STAGE07_EXECUTION_NOT_CONFIRMED")
    manifest = json.loads((output / "run_metaviro_provenance.json").read_text())
    if not (manifest["stage07_completed"] and manifest["context_screen_completed"]):
        raise RuntimeError("INCOMPLETE_STAGE07_MANIFEST")
    screen = output / "review/10_context_screen/all_family_context_screen.tsv"
    if not screen.is_file():
        raise RuntimeError("MISSING_CONTEXT_REVIEW")
    print("STAGE00_07_DATABASE_EXECUTION=PASS", flush=True)
    return {"mode": "database", "status": "PASS", "executed": True,
            "family_count": manifest["family_count"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("smoke", "database"), default="smoke")
    ap.add_argument("--config", type=Path, help="User-owned config with real DB paths for database mode")
    ap.add_argument("--screen-out-dir", type=Path, help="Nonexistent run workspace (database mode)")
    ap.add_argument("--normalize-ids", action="store_true")
    ap.add_argument("--execute", action="store_true", help="Run the database-backed Stage 00-07 analysis")
    ap.add_argument("--report", type=Path, help="Optional JSON report; no outputs written by default")
    args = ap.parse_args()
    if args.mode == "smoke" and (args.config or args.screen_out_dir or args.execute or args.normalize_ids):
        ap.error("Database-specific arguments require --mode database")
    result = smoke() if args.mode == "smoke" else database(args)
    if args.report:
        target = args.report.resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise ValueError(f"REPORT_ALREADY_EXISTS:{target}")
        target.write_text(json.dumps(result, indent=2) + "\n")
        print(f"VALIDATION_REPORT={target}")
    print("RELEASE_VALIDATION=" + result["status"])


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        sys.exit(f"RELEASE_VALIDATION=FAIL: {exc}")
