#!/usr/bin/env python3
"""Stage 09-10 orchestrator with verified zero-hit support.

Check or run using existing Foldseek searches or fresh searches.
Stage 08, GPU prediction, and final selection are not modified.
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

SOURCES = {
    "pdb_e1": "PDB",
    "afsp": "AF_SwissProt",
    "bfvd": "BFVD",
}


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def require(path, allow_empty=False):
    if not path.is_file():
        raise ValueError(f"MISSING_FILE:{path}")
    if not allow_empty and path.stat().st_size == 0:
        raise ValueError(f"EMPTY_FILE:{path}")
    return path


def read_tsv(path):
    with require(path).open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    if not fields or not rows:
        raise ValueError(f"EMPTY_OR_INVALID_TABLE:{path}")
    return fields, rows


def read_json(path):
    obj = json.loads(require(path).read_text())
    if not isinstance(obj, dict):
        raise ValueError(f"INVALID_JSON_PROVENANCE:{path}")
    return obj


def ensure_separate(out, source):
    source = source.resolve()
    if (
        out == source
        or source in out.parents
        or out in source.parents
    ):
        raise ValueError("INPUT_OUTPUT_OVERLAP")


def run(label, script, arguments, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable, "-B",
        str(SCRIPTS / script),
        *map(str, arguments),
    ]
    print(f"[RUN] {label}", flush=True)
    with log.open("x", encoding="utf-8") as fh:
        result = subprocess.run(
            command, cwd=ROOT,
            stdout=fh, stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(
            f"MODULE_FAILED:{label}:SEE_LOG:{log}"
        )
    print(f"[PASS] {label}", flush=True)


def verify_search(directory, source, structures):
    """Allow empty raw TSV only with matching completed-search provenance."""
    directory = directory.resolve()

    raw = require(
        directory / "raw_hits.tsv",
        allow_empty=True,
    )
    meta_path = require(
        directory / "search_provenance.json"
    )
    require(directory / "source_summary.tsv")
    require(directory / "query_manifest.tsv")
    require(directory / "search.log")
    require(directory / "convertalis.log")

    meta = read_json(meta_path)

    if meta.get("source_profile") != source:
        raise ValueError(
            f"SEARCH_SOURCE_MISMATCH:{source}"
        )

    if meta.get("control_only") is not False:
        raise ValueError(
            f"CONTROL_SEARCH_REJECTED:{source}"
        )

    recorded_input = meta.get("input_structure_table")

    if (
        not isinstance(recorded_input, str)
        or Path(recorded_input).resolve()
        != structures.resolve()
    ):
        raise ValueError(
            f"SEARCH_STRUCTURE_INPUT_MISMATCH:{source}"
        )

    if meta.get("raw_hits_sha256") != sha(raw):
        raise ValueError(
            f"RAW_HITS_SHA256_MISMATCH:{source}"
        )

    # Search settings required for the Stage 10 evidence profile.
    if (
        float(meta.get("evalue", -1)) != 1.0
        or int(meta.get("max_seqs", -1))
        != (500 if source == "bfvd" else 1000)
    ):
        raise ValueError(
            f"SEARCH_PROFILE_PARAMETERS_MISMATCH:{source}"
        )

    target = meta.get("target_db")
    target_digest = meta.get("target_dbtype_sha256")

    if not isinstance(target, str) or not isinstance(
        target_digest, str
    ):
        raise ValueError(
            f"SEARCH_DATABASE_PROVENANCE_MISSING:{source}"
        )

    dbtype = Path(target + ".dbtype")

    if dbtype.is_file():
        if sha(dbtype) != target_digest:
            raise ValueError(
                f"TARGET_DBTYPE_SHA256_MISMATCH:{source}"
            )
        db_status = "LOCAL_DBTYPE_VERIFIED"
    else:
        db_status = "LOCAL_DB_NOT_AVAILABLE"

    hit_count = sum(
        1 for line in raw.open(encoding="utf-8")
        if line.strip()
    )

    if raw.stat().st_size == 0:
        assert hit_count == 0
        print(
            f"ZERO_HIT_VERIFIED={source}:"
            "EMPTY_SHA256_MATCHED"
        )
    else:
        print(
            f"RAW_HITS_VERIFIED={source}:{hit_count}"
        )

    print(f"DATABASE_STATUS={source}:{db_status}")

    return {
        "source": source,
        "raw_hits_sha256": sha(raw),
        "provenance_sha256": sha(meta_path),
        "raw_hit_count": hit_count,
        "database_status": db_status,
    }


def validate(args):
    upstream = args.stage08_output.resolve()
    out = args.out_dir.resolve()

    if not upstream.is_dir():
        raise ValueError("STAGE08_OUTPUT_MISSING")

    if out.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    ensure_separate(out, upstream)

    audit = require(
        upstream / "06_darkness/family_darkness_audit.tsv"
    )
    family = require(
        upstream / "07_family_context/family_context.tsv"
    )
    member = require(
        upstream / "07_family_context/member_context.tsv"
    )

    candidates = sorted(
        (
            upstream / "08_structure"
        ).glob("run_*/collected/structure_results.tsv")
    )

    if args.structure_results:
        structures = args.structure_results.resolve()
        if structures not in [x.resolve() for x in candidates]:
            raise ValueError("INVALID_STAGE08_STRUCTURE_TABLE")
    else:
        if len(candidates) != 1:
            raise ValueError(
                "SPECIFY_STRUCTURE_RESULTS_FOR_MULTIPLE_RUNS"
            )
        structures = candidates[0].resolve()

    fields, rows = read_tsv(structures)

    for name in (
        "structure_id", "family_rep",
        "pdb_file", "pdb_sha256",
    ):
        if name not in fields:
            raise ValueError(
                f"MISSING_STRUCTURE_FIELD:{name}"
            )

    ids = [r["structure_id"] for r in rows]
    if len(ids) != len(set(ids)) or any(not x for x in ids):
        raise ValueError("INVALID_STAGE08_STRUCTURE_IDS")

    for row in rows:
        pdb = require(Path(row["pdb_file"]))
        if sha(pdb) != row["pdb_sha256"]:
            raise ValueError("STAGE08_PDB_SHA256_MISMATCH")

    if args.threads < 1:
        raise ValueError("INVALID_THREADS")

    searches = None
    db = None

    if args.reuse_search_root:
        searches = args.reuse_search_root.resolve()
        if not searches.is_dir():
            raise ValueError("REUSE_SEARCH_DIRECTORY_MISSING")
        ensure_separate(out, searches)

        for source in SOURCES:
            verify_search(
                searches / source,
                source,
                structures,
            )

    else:
        db = args.db_root.resolve()
        if not db.is_dir():
            raise ValueError("FOLDSEEK_DB_ROOT_MISSING")
        ensure_separate(out, db)

        for name in SOURCES.values():
            prefix = db / name
            for suffix in ("", ".dbtype", ".index"):
                require(Path(str(prefix) + suffix))

    command = [
        sys.executable, "-B",
        str(
            SCRIPTS /
            "10_prepare_evidence_portable.py"
        ),
        "--family-audit", str(audit),
        "--family-context", str(family),
        "--member-context", str(member),
        "--structures", str(structures),
        "--preflight",
    ]

    result = subprocess.run(
        command, cwd=ROOT,
        capture_output=True, text=True,
        check=False,
    )

    if (
        result.returncode
        or "STAGE06_07_08_FAMILY_CONSISTENCY=PASS"
        not in result.stdout
    ):
        raise ValueError(
            "STAGE10_UPSTREAM_PREFLIGHT_FAILED:"
            + (result.stderr or result.stdout)
        )

    print("STAGE06_07_08_CONSISTENCY=PASS")

    return {
        "upstream": upstream,
        "out": out,
        "audit": audit,
        "family": family,
        "member": member,
        "structures": structures,
        "reuse": searches,
        "db": db,
    }


def execute(args, state):
    out = state["out"]
    out.mkdir(parents=True, exist_ok=False)
    logs = out / "logs"

    search_dirs = {}

    if state["reuse"]:
        print("SEARCH_MODE=VERIFIED_REUSE")
        for source in SOURCES:
            search_dirs[source] = (
                state["reuse"] / source
            )
    else:
        print("SEARCH_MODE=FRESH_FOLDSEEK")
        for source, db_name in SOURCES.items():
            destination = out / "09_search" / source

            run(
                f"09 Foldseek {source}",
                "09_foldseek_search_portable.py",
                [
                    "--structures", state["structures"],
                    "--target-db", state["db"] / db_name,
                    "--source", source,
                    "--out-dir", destination,
                    "--threads", args.threads,
                ],
                logs / f"09_{source}.log",
            )

            search_dirs[source] = destination

    search_manifest = {}

    for source, directory in search_dirs.items():
        search_manifest[source] = verify_search(
            directory, source, state["structures"]
        )

    stage09 = out / "09_external"

    run(
        "09 integration",
        "09_external_integration_portable.py",
        [
            "--structures", state["structures"],
            "--pdb-e1-dir", search_dirs["pdb_e1"],
            "--afsp-dir", search_dirs["afsp"],
            "--bfvd-dir", search_dirs["bfvd"],
            "--out-dir", stage09,
        ],
        logs / "09_integration.log",
    )

    matrix = require(
        stage09 / "external_structure_matrix.tsv"
    )
    external_prov = require(
        stage09 / "integration_provenance.json"
    )

    evidence_dir = out / "10_evidence"

    run(
        "10 evidence preparation",
        "10_prepare_evidence_portable.py",
        [
            "--family-audit", state["audit"],
            "--family-context", state["family"],
            "--member-context", state["member"],
            "--structures", state["structures"],
            "--external-matrix", matrix,
            "--external-provenance", external_prov,
            "--out-dir", evidence_dir,
        ],
        logs / "10_evidence.log",
    )

    evidence = require(
        evidence_dir / "stage10_input_evidence.tsv"
    )
    member_i2 = require(
        evidence_dir / "member_context_i2.tsv"
    )
    bridge_prov = require(
        evidence_dir / "evidence_provenance.json"
    )

    rank_dir = out / "10_ranking"

    run(
        "10 integrated ranking",
        "10_integrated_ranking_portable.py",
        [
            "--evidence", evidence,
            "--member-context", member_i2,
            "--bridge-provenance", bridge_prov,
            "--stage09-matrix", matrix,
            "--stage09-provenance", external_prov,
            "--out-dir", rank_dir,
        ],
        logs / "10_ranking.log",
    )

    ranking = require(
        rank_dir / "confidence_aware_ranking.tsv"
    )
    ranking_prov = require(
        rank_dir / "ranking_provenance.json"
    )

    shortlist_dir = out / "10_shortlists"

    run(
        "10 shortlists",
        "10_ranked_shortlists_portable.py",
        [
            "--ranking", ranking,
            "--ranking-provenance", ranking_prov,
            "--out-dir", shortlist_dir,
        ],
        logs / "10_shortlists.log",
    )

    top20 = require(
        shortlist_dir / "top20_confidence_aware.tsv"
    )
    top15 = require(
        shortlist_dir / "top15_experimental_candidates.tsv"
    )

    matrix_fields, matrix_rows = read_tsv(matrix)
    ranking_fields, ranking_rows = read_tsv(ranking)
    _, top20_rows = read_tsv(top20)
    _, top15_rows = read_tsv(top15)

    _, structure_rows = read_tsv(state["structures"])

    expected_ids = {
        row["structure_id"] for row in structure_rows
    }

    if (
        {r["structure_id"] for r in matrix_rows}
        != expected_ids
        or {r["structure_id"] for r in ranking_rows}
        != expected_ids
    ):
        raise ValueError("STAGE08_09_10_ID_MISMATCH")

    if not (
        {r["structure_id"] for r in top20_rows}
        .issubset(expected_ids)
        and {r["structure_id"] for r in top15_rows}
        .issubset(expected_ids)
    ):
        raise ValueError("INVALID_SHORTLIST_ID")

    meta09 = read_json(external_prov)
    meta10 = read_json(bridge_prov)

    if (
        meta09.get("matrix_sha256") != sha(matrix)
        or meta10.get("evidence_sha256")
        != sha(evidence)
        or meta10.get("member_context_sha256")
        != sha(member_i2)
    ):
        raise ValueError("INTEGRATED_PROVENANCE_MISMATCH")

    artifacts = {
        "stage09_matrix": matrix,
        "stage09_provenance": external_prov,
        "stage10_evidence": evidence,
        "stage10_member_context": member_i2,
        "stage10_bridge_provenance": bridge_prov,
        "ranking": ranking,
        "ranking_provenance": ranking_prov,
        "top20": top20,
        "top15": top15,
    }

    manifest = {
        "stage": "STAGE09_10_PORTABLE_ORCHESTRATION",
        "search_mode": (
            "VERIFIED_REUSE"
            if state["reuse"] else "FRESH_FOLDSEEK"
        ),
        "stage08_structure_table_sha256":
            sha(state["structures"]),
        "searches": search_manifest,
        "artifacts": {
            name: {
                "path": str(path.resolve()),
                "sha256": sha(path),
            }
            for name, path in artifacts.items()
        },
        "family_count": len(ranking_rows),
        "ranking_columns": len(ranking_fields),
        "stage09_columns": len(matrix_fields),
        "top20_count": len(top20_rows),
        "top15_count": len(top15_rows),
        "gpu_prediction_performed": False,
        "final_candidate_selection_performed": False,
    }

    with (
        out / "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    print(f"RANKED_FAMILIES={len(ranking_rows)}")
    print(f"RANKING_COLUMNS={len(ranking_fields)}")
    print(f"TOP20_ROWS={len(top20_rows)}")
    print(f"TOP15_ROWS={len(top15_rows)}")
    print("STAGE09_10_ORCHESTRATION=PASS")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage08-output", type=Path, required=True
    )
    parser.add_argument(
        "--structure-results", type=Path
    )

    group = parser.add_mutually_exclusive_group(
        required=True
    )
    group.add_argument(
        "--reuse-search-root", type=Path
    )
    group.add_argument(
        "--db-root", type=Path
    )

    parser.add_argument(
        "--out-dir", type=Path, required=True
    )
    parser.add_argument(
        "--threads", type=int, default=2
    )
    parser.add_argument(
        "--mode", choices=("check", "run"),
        default="check",
    )

    args = parser.parse_args()
    state = validate(args)

    if args.mode == "check":
        print("CHECK_MODE=PASS")
        print("OUTPUT_CREATED=NO")
        return

    execute(args, state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
