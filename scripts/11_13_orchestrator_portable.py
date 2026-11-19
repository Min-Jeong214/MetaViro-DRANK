#!/usr/bin/env python3
"""Stage 10-to-13 verified-reuse orchestrator.

Reuses existing MSA, prediction, conservation, spatial, and K2
evidence. Does not run GPU predictions, Foldseek searches,
Stage 10 scoring, or final Top5/Top3 candidate selection.
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
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def require(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"MISSING_OR_EMPTY_FILE:{path}")
    return path


def read_json(path):
    value = json.loads(require(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"INVALID_PROVENANCE:{path}")
    return value


def read_table(path):
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

    if any(
        None in row
        or any(value is None for value in row.values())
        for row in rows
    ):
        raise ValueError(f"MALFORMED_TABLE:{path}")

    ids = [row["structure_id"] for row in rows]

    if (
        len(ids) != len(set(ids))
        or any(not sid for sid in ids)
    ):
        raise ValueError(f"INVALID_STRUCTURE_IDS:{path}")

    return fields, rows


def indexed(rows):
    return {
        row["structure_id"]: row
        for row in rows
    }


def verify_sources(args):
    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    sources = [
        args.stage08_output,
        args.stage09_10_output,
        args.msa_qc_dir,
        args.conservation_dir,
        args.spatial_dir,
        args.k2_structures,
        args.external_dir,
    ]

    for source in sources:
        source = source.resolve()

        if (
            out == source
            or source in out.parents
            or out in source.parents
        ):
            raise ValueError("INPUT_OUTPUT_OVERLAP")

    up = args.stage08_output.resolve()
    down = args.stage09_10_output.resolve()
    qcdir = args.msa_qc_dir.resolve()
    consdir = args.conservation_dir.resolve()
    spatdir = args.spatial_dir.resolve()
    extdir = args.external_dir.resolve()

    manifest_file = require(
        down / "orchestration_provenance.json"
    )
    manifest = read_json(manifest_file)

    rank = require(
        down / "10_ranking/confidence_aware_ranking.tsv"
    )
    rank_prov = require(
        down / "10_ranking/ranking_provenance.json"
    )
    shortlist = require(
        down / "10_shortlists/top15_experimental_candidates.tsv"
    )

    artifacts = manifest.get("artifacts", {})

    for name, path in (
        ("ranking", rank),
        ("top15", shortlist),
    ):
        item = artifacts.get(name, {})

        if (
            item.get("sha256") != sha(path)
            or Path(item.get("path", "/")).resolve()
            != path.resolve()
        ):
            raise ValueError(
                f"STAGE10_MANIFEST_MISMATCH:{name}"
            )

    rank_meta = read_json(rank_prov)
    ranking_hashes = [
        value
        for key, value in rank_meta.items()
        if key.endswith("sha256")
        and isinstance(value, str)
    ]

    if sha(rank) not in ranking_hashes:
        raise ValueError(
            "RANKING_PROVENANCE_MISMATCH"
        )

    rank_fields, rank_rows = read_table(rank)
    short_fields, short_rows = read_table(shortlist)

    if rank_fields != short_fields:
        raise ValueError("SHORTLIST_SCHEMA_MISMATCH")

    rank_by_id = indexed(rank_rows)

    for row in short_rows:
        if (
            row["structure_id"] not in rank_by_id
            or rank_by_id[row["structure_id"]] != row
        ):
            raise ValueError(
                "SHORTLIST_RANKING_ROW_MISMATCH"
            )

    ids = [row["structure_id"] for row in short_rows]

    if not ids:
        raise ValueError("NO_SELECTED_CANDIDATES")

    selection_bytes = (
        "structure_id\n"
        + "".join(sid + "\n" for sid in ids)
    ).encode("utf-8")

    selection_sha = hashlib.sha256(
        selection_bytes
    ).hexdigest()

    structures_found = sorted(
        (
            up / "08_structure"
        ).glob("run_*/collected/structure_results.tsv")
    )

    matching_structures = [
        path for path in structures_found
        if sha(path)
        == manifest.get(
            "stage08_structure_table_sha256"
        )
    ]

    if len(matching_structures) != 1:
        raise ValueError(
            "STAGE08_STRUCTURE_SOURCE_NOT_UNIQUE"
        )

    structures = matching_structures[0]
    _, stage08_rows = read_table(structures)
    stage08 = indexed(stage08_rows)

    qc_file = require(
        qcdir / "msa_structure_qc.tsv"
    )
    qc_meta_file = require(
        qcdir / "msa_structure_qc_provenance.json"
    )
    qc_meta = read_json(qc_meta_file)
    _, qc_rows = read_table(qc_file)
    qc = indexed(qc_rows)

    if (
        qc_meta.get("source_stage08_sha256")
        != sha(structures)
        or qc_meta.get("msa_qc_sha256")
        != sha(qc_file)
    ):
        raise ValueError(
            "STAGE08_TO_STAGE11_PROVENANCE_MISMATCH"
        )

    cons_meta_file = require(
        consdir / "conservation_provenance.json"
    )
    cons_meta = read_json(cons_meta_file)

    spat_meta_file = require(
        spatdir / "spatial_cluster_provenance.json"
    )
    spat_meta = read_json(spat_meta_file)

    if cons_meta.get("selection_sha256") != selection_sha:
        raise ValueError(
            "STAGE12_SELECTION_SHA256_MISMATCH:"
            "EXISTING_CONSERVATION_CANNOT_BE_REUSED"
        )

    if (
        cons_meta.get("prediction_index_sha256")
        != qc_meta.get("prediction_index_sha256")
        or spat_meta.get("stage12_provenance_sha256")
        != sha(cons_meta_file)
    ):
        raise ValueError(
            "STAGE11_12_SPATIAL_PROVENANCE_MISMATCH"
        )

    k2_file = require(args.k2_structures.resolve())
    _, k2_rows = read_table(k2_file)
    k2 = indexed(k2_rows)

    external_file = require(
        extdir / "msa_external_evidence_stage13.tsv"
    )
    external_meta_file = require(
        extdir / "msa_external_evidence_provenance.json"
    )

    external_meta = read_json(external_meta_file)
    _, external_rows = read_table(external_file)
    external = indexed(external_rows)

    if (
        external_meta.get("stage13_evidence_sha256")
        != sha(external_file)
        or external_meta.get("control_only") is not False
        or external_meta.get(
            "final_selection_performed"
        ) is not False
    ):
        raise ValueError(
            "K2_EXTERNAL_PROVENANCE_MISMATCH"
        )

    if set(external) != set(ids):
        raise ValueError(
            "K2_EXTERNAL_SELECTION_MISMATCH"
        )

    for sid in ids:
        if not all(
            sid in group
            for group in (stage08, qc, k2, external)
        ):
            raise ValueError(
                f"CROSS_STAGE_MISSING_STRUCTURE_ID:{sid}"
            )

        single = stage08[sid]
        checked = qc[sid]
        msa = k2[sid]
        ext = external[sid]

        if (
            single["family_rep"]
            != msa["family_rep"]
            or checked["family_rep"]
            != single["family_rep"]
            or single["sequence_sha256"]
            != msa["sequence_sha256"]
            or single["pdb_sha256"]
            != checked["single_pdb_sha256"]
            or msa["pdb_sha256"]
            != checked["msa_pdb_sha256"]
            or msa["pdb_sha256"]
            != ext["msa_pdb_sha256"]
        ):
            raise ValueError(
                f"STRUCTURE_PROVENANCE_MISMATCH:{sid}"
            )

        for row in (single, msa):
            pdb = require(Path(row["pdb_file"]))

            if sha(pdb) != row["pdb_sha256"]:
                raise ValueError(
                    f"PDB_SHA256_MISMATCH:{sid}"
                )

    print("STAGE10_SHORTLIST_PROVENANCE=PASS")
    print("STAGE08_TO_STAGE11_PDB_PROVENANCE=PASS")
    print("STAGE11_TO_STAGE12_PROVENANCE=PASS")
    print("STAGE12_TO_SPATIAL_PROVENANCE=PASS")
    print("STAGE_K2_EXTERNAL_PROVENANCE=PASS")
    print("STAGE12_REUSE_SELECTION_SHA256=PASS")
    print(f"SELECTED_FAMILIES={len(ids)}")

    return {
        "out": out,
        "rank": rank,
        "rank_prov": rank_prov,
        "shortlist": shortlist,
        "selection_bytes": selection_bytes,
        "selection_sha": selection_sha,
        "ids": ids,
        "structures": structures,
        "qc_file": qc_file,
        "qc_meta_file": qc_meta_file,
        "consdir": consdir,
        "spatdir": spatdir,
        "ext_file": external_file,
        "ext_meta_file": external_meta_file,
        "k2_file": k2_file,
    }


def run_module(label, script, arguments, log):
    log.parent.mkdir(
        parents=True, exist_ok=True
    )

    command = [
        sys.executable, "-B",
        str(SCRIPTS / script),
        *map(str, arguments),
    ]

    print(f"[RUN] {label}", flush=True)

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


def execute(state):
    out = state["out"]
    out.mkdir(parents=True, exist_ok=False)

    selection_dir = out / "selection"
    selection_dir.mkdir()

    selection = (
        selection_dir /
        "stage10_selected_families.tsv"
    )

    with selection.open("xb") as fh:
        fh.write(state["selection_bytes"])

    if sha(selection) != state["selection_sha"]:
        raise ValueError(
            "GENERATED_SELECTION_HASH_MISMATCH"
        )

    selection_prov = (
        selection_dir / "selection_provenance.json"
    )

    record = {
        "stage": "STAGE10_TO_STAGE11_13_SELECTION",
        "selection_method": "VERIFIED_STAGE10_TOP15",
        "source_shortlist": str(state["shortlist"]),
        "source_shortlist_sha256":
            sha(state["shortlist"]),
        "source_ranking_sha256":
            sha(state["rank"]),
        "source_ranking_provenance_sha256":
            sha(state["rank_prov"]),
        "source_stage08_sha256":
            sha(state["structures"]),
        "source_stage11_qc_sha256":
            sha(state["qc_file"]),
        "source_k2_structure_table_sha256":
            sha(state["k2_file"]),
        "selected_ids": state["ids"],
        "selection_sha256": sha(selection),
        "biological_priority_validated": False,
        "historical_top5_selected": False,
        "historical_top3_selected": False,
    }

    with selection_prov.open(
        "x", encoding="utf-8"
    ) as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")

    logs = out / "logs"

    prelim_dir = out / "13_preliminary"

    run_module(
        "13 preliminary evidence",
        "13_candidate_evidence_preliminary_portable.py",
        [
            "--selection", selection,
            "--msa-qc", state["qc_file"],
            "--msa-qc-provenance",
            state["qc_meta_file"],
            "--conservation-dir",
            state["consdir"],
            "--spatial-dir",
            state["spatdir"],
            "--out-dir", prelim_dir,
        ],
        logs / "13_preliminary.log",
    )

    preliminary = (
        prelim_dir / "preliminary_candidate_evidence.tsv"
    )
    preliminary_prov = (
        prelim_dir / "preliminary_evidence_provenance.json"
    )

    merged_dir = out / "13_external"

    run_module(
        "13 K2 external evidence merge",
        "13_merge_external_evidence_portable.py",
        [
            "--preliminary", preliminary,
            "--preliminary-provenance",
            preliminary_prov,
            "--external", state["ext_file"],
            "--external-provenance",
            state["ext_meta_file"],
            "--out-dir", merged_dir,
        ],
        logs / "13_external.log",
    )

    merged = (
        merged_dir /
        "preliminary_candidate_evidence_with_external.tsv"
    )
    merged_prov = (
        merged_dir /
        "preliminary_external_merge_provenance.json"
    )

    ranked_dir = out / "13_ranked"

    run_module(
        "13 attach verified Stage 10 ranking",
        "13_attach_stage10_ranking_portable.py",
        [
            "--ranking", state["rank"],
            "--ranking-provenance",
            state["rank_prov"],
            "--shortlist", state["shortlist"],
            "--selection", selection,
            "--selection-provenance",
            selection_prov,
            "--preliminary", preliminary,
            "--preliminary-provenance",
            preliminary_prov,
            "--merged", merged,
            "--merged-provenance",
            merged_prov,
            "--out-dir", ranked_dir,
        ],
        logs / "13_ranking_handoff.log",
    )

    ranked = (
        ranked_dir / "ranked_candidate_evidence.tsv"
    )
    ranked_prov = (
        ranked_dir /
        "ranked_candidate_evidence_provenance.json"
    )

    fields, rows = read_table(ranked)

    if (
        len(fields) != 33
        or {row["structure_id"] for row in rows}
        != set(state["ids"])
    ):
        raise ValueError(
            "FINAL_STAGE13_EVIDENCE_MISMATCH"
        )

    artifacts = {
        "selection": selection,
        "selection_provenance": selection_prov,
        "preliminary": preliminary,
        "preliminary_provenance": preliminary_prov,
        "merged": merged,
        "merged_provenance": merged_prov,
        "ranked": ranked,
        "ranked_provenance": ranked_prov,
    }

    manifest = {
        "stage": "STAGE11_13_VERIFIED_REUSE",
        "mode": "REUSE_EXISTING_VALIDATED_ANALYSES",
        "source_stage10_ranking_sha256":
            sha(state["rank"]),
        "source_stage10_shortlist_sha256":
            sha(state["shortlist"]),
        "source_stage11_qc_sha256":
            sha(state["qc_file"]),
        "source_stage12_provenance_sha256":
            sha(
                state["consdir"] /
                "conservation_provenance.json"
            ),
        "selected_ids": state["ids"],
        "result_columns": len(fields),
        "result_families": len(rows),
        "artifacts": {
            name: {
                "path": str(path.resolve()),
                "sha256": sha(path),
            }
            for name, path in artifacts.items()
        },
        "gpu_prediction_performed": False,
        "foldseek_search_performed": False,
        "stage10_scores_recalculated": False,
        "final_candidate_selection_performed": False,
    }

    with (
        out / "orchestration_provenance.json"
    ).open("x", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    print(f"STAGE13_COLUMNS={len(fields)}")
    print(f"STAGE13_FAMILIES={len(rows)}")
    print(f"STAGE13_RESULT={ranked}")
    print("STAGE11_13_VERIFIED_REUSE_ORCHESTRATION=PASS")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--stage08-output", type=Path, required=True
    )
    parser.add_argument(
        "--stage09-10-output", type=Path, required=True
    )
    parser.add_argument(
        "--msa-qc-dir", type=Path, required=True
    )
    parser.add_argument(
        "--conservation-dir", type=Path, required=True
    )
    parser.add_argument(
        "--spatial-dir", type=Path, required=True
    )
    parser.add_argument(
        "--k2-structures", type=Path, required=True
    )
    parser.add_argument(
        "--external-dir", type=Path, required=True
    )
    parser.add_argument(
        "--out-dir", type=Path, required=True
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
        "--stage09-10-output",
        str(args.stage09_10_output),
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
    state = verify_sources(args)

    if args.mode == "check":
        print("CHECK_MODE=PASS")
        print("OUTPUT_CREATED=NO")
        return

    execute(state)


if __name__ == "__main__":
    try:
        main()
    except (
        OSError, ValueError, RuntimeError,
        KeyError, csv.Error
    ) as exc:
        sys.exit(f"ERROR: {exc}")
