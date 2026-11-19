#!/usr/bin/env python3
"""MetaViro-DRANK portable pipeline (Stages 01-02)."""

import argparse
import csv
import re
import uuid
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
REQUIRED = (
    "proteins.faa",
    "genes.gff",
    "annotations.emapper.annotations",
)


def resolve(value):
    p = Path(value).expanduser()
    return (p if p.is_absolute() else ROOT / p).resolve()


def run(label, command, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"[RUN] {label}", flush=True)
    with log.open("w") as fh:
        result = subprocess.run(
            command, stdout=fh, stderr=subprocess.STDOUT,
            cwd=ROOT, check=False
        )
    if result.returncode:
        raise RuntimeError(f"{label} failed. See {log}")
    print(f"[PASS] {label}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--through",
        choices=["02", "03", "04", "05", "06", "07", "08"],
        default="02",
    )
    parser.add_argument(
        "--structure-mode",
        choices=["prepare", "collect", "predict"],
        default="prepare",
    )
    parser.add_argument(
        "--structure-prediction-dir",
        default=None,
        help="Existing ColabFold PDB directory for collect mode",
    )
    parser.add_argument(
        "--structure-env",
        default="metaviro-structure",
        help="Conda environment for ColabFold and BioPython",
    )
    args = parser.parse_args()

    if args.through != "08":
        if (
            args.structure_mode != "prepare"
            or args.structure_prediction_dir is not None
        ):
            raise ValueError(
                "Structure options require --through 08"
            )

    if args.through == "08":
        if (
            args.structure_mode == "collect"
            and not args.structure_prediction_dir
        ):
            raise ValueError(
                "Collect mode requires --structure-prediction-dir"
            )

        if (
            args.structure_mode != "collect"
            and args.structure_prediction_dir is not None
        ):
            raise ValueError(
                "Prediction directory is only valid in collect mode"
            )

        if args.structure_mode == "collect":
            prediction_source = resolve(
                args.structure_prediction_dir
            )
            if (
                not prediction_source.is_dir()
                or not any(prediction_source.glob("*.pdb"))
            ):
                raise ValueError(
                    "Missing predicted PDB files: "
                    + str(prediction_source)
                )

    with resolve(args.config).open() as fh:
        cfg = yaml.safe_load(fh)

    inp = resolve(cfg["input_dir"])
    out = resolve(cfg["output_dir"])
    min_aa = int(cfg["min_aa"])
    window = int(cfg["window"])

    if inp == out or inp in out.parents or out in inp.parents:
        raise ValueError("Input and output directories overlap")
    if min_aa < 1 or window < 0:
        raise ValueError("Invalid min_aa or window")
    if not inp.is_dir():
        raise ValueError(f"Input directory missing: {inp}")

    samples = sorted(
        p for p in inp.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )

    seen = set()
    for sample in samples:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample.name):
            raise ValueError(f"Invalid sample name: {sample.name}")
        for filename in REQUIRED:
            path = sample / filename
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing or empty: {path}")
        with (sample / "proteins.faa").open() as fh:
            for line in fh:
                if line.startswith(">"):
                    parts = line[1:].split()
                    if not parts or parts[0] in seen:
                        raise ValueError(
                            f"Empty or duplicate protein ID: {line.strip()}"
                        )
                    seen.add(parts[0])

    print(f"SAMPLES={len(samples)}")
    print(f"PROTEINS={len(seen)}")

    # Stage 00: validate all inputs before analysis.
    for sample in samples:
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts/00_validate_input_portable.py"),
                "--sample-dir",
                str(sample),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            details = (
                result.stderr or result.stdout
            ).strip()

            raise ValueError(
                f"Stage 00 validation failed "
                f"for {sample.name}: {details}"
            )

        print(
            f"[PASS] 00 input validation {sample.name}",
            flush=True
        )

    if samples:
        print("INPUT_CONSISTENCY=PASS")



    db_paths = {}

    if args.through in ("05", "06", "07", "08"):
        databases = cfg.get("databases")

        if not isinstance(databases, dict):
            raise ValueError("Missing database configuration")

        for name in (
            "swissprot",
            "phrog",
            "pfam",
            "phrog_annotations",
        ):
            if name not in databases:
                raise ValueError(
                    f"Missing database configuration: {name}"
                )

            path = resolve(databases[name])

            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(
                    f"Missing or empty database: {path}"
                )

            db_paths[name] = path

        for key in ("swissprot", "phrog"):
            db = db_paths[key]
            if not Path(str(db) + ".dbtype").is_file():
                raise ValueError(f"{key} MMseqs database incomplete: missing .dbtype")
        for suffix in (".index", "_h", "_h.dbtype", "_h.index"):
            if not Path(str(db_paths["phrog"]) + suffix).is_file():
                raise ValueError(f"PHROG profile incomplete: missing {suffix}")
        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
            if not Path(str(db_paths["pfam"]) + suffix).is_file():
                raise ValueError(f"Pfam HMM index missing: {suffix}")

        print("DATABASE_PREFLIGHT=PASS")

    if args.check:
        print("INPUT_PREFLIGHT=" + ("PASS" if samples else "NO_SAMPLES"))
        return

    if not samples:
        raise ValueError("No sample folders in input/")

    # Capture immutable pre-run information before Stage 01.
    # Each execution receives a separate provenance file.
    run_id = uuid.uuid4().hex
    provenance_file = (
        out / "provenance" / f"run_{run_id}.json"
    )

    run("00 capture run provenance", [
        sys.executable,
        str(ROOT / "scripts/00_capture_run_provenance.py"),
        "--config", str(resolve(args.config)),
        "--out", str(provenance_file),
        "--through", args.through,
        "--structure-mode", args.structure_mode,
        "--structure-env", args.structure_env,
        "--structure-prediction-dir",
        (
            str(resolve(args.structure_prediction_dir))
            if args.structure_prediction_dir else ""
        ),
    ], out / "logs" / f"provenance_{run_id}.log")

    if (
        not provenance_file.is_file()
        or provenance_file.stat().st_size == 0
    ):
        raise RuntimeError(
            "Run provenance was not generated"
        )

    print(
        f"RUN_PROVENANCE_FILE={provenance_file}",
        flush=True
    )

    for sample in samples:
        name = sample.name
        s1 = out / "01_unresolved" / f"{name}.tsv"
        s2 = out / "02_context" / f"{name}.tsv"
        anchors = out / "02_context" / f"{name}_anchors.tsv"
        s1.parent.mkdir(parents=True, exist_ok=True)
        s2.parent.mkdir(parents=True, exist_ok=True)

        run(f"01 {name}", [
            sys.executable, str(ROOT / "scripts/01_define_unresolved_candidates.py"),
            "--annotations", str(sample / REQUIRED[2]),
            "--faa", str(sample / REQUIRED[0]),
            "--out", str(s1), "--min-aa", str(min_aa)
        ], out / "logs" / f"{name}_01.log")

        run(f"02 {name}", [
            sys.executable, str(ROOT / "scripts/02_extract_replication_context.py"),
            "--sample", name,
            "--gff", str(sample / REQUIRED[1]),
            "--faa", str(sample / REQUIRED[0]),
            "--annotations", str(sample / REQUIRED[2]),
            "--unresolved-tsv", str(s1),
            "--out", str(s2),
            "--anchors-out", str(anchors),
            "--window", str(window),
            "--min-aa", str(min_aa)
        ], out / "logs" / f"{name}_02.log")

    if args.through == "02":
        print("STAGES_01_02=COMPLETE")
        return

    threads = int(cfg.get("threads", 8))
    if threads < 1:
        raise ValueError("threads must be >= 1")

    stage3 = out / "03_clustering"
    stage3.mkdir(parents=True, exist_ok=True)

    combined = stage3 / "all_viral_proteins.faa"
    prefix = stage3 / "viral_families"

    with combined.open("w") as dst:
        for sample in samples:
            faa = sample / "proteins.faa"
            with faa.open() as src:
                for line in src:
                    dst.write(line.rstrip("\r\n") + "\n")

    run("03 MMseqs2 clustering", [
        "bash",
        str(ROOT / "scripts/03_global_family_clustering.sh"),
        str(combined),
        str(prefix),
        str(threads)
    ], out / "logs" / "stage03.log")

    cluster_tsv = Path(str(prefix) + "_cluster.tsv")

    if not cluster_tsv.is_file():
        raise RuntimeError(
            f"Missing clustering result: {cluster_tsv}"
        )

    print(f"CLUSTER_TABLE={cluster_tsv}")
    if args.through == "03":
        print("STAGES_01_03=COMPLETE")
        return

    stage4 = out / "04_recurrent"
    stage4.mkdir(parents=True, exist_ok=True)
    candidates_file = stage4 / "candidate_provenance.tsv"

    count = 0
    with candidates_file.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["candidate_id", "sample"],
            delimiter="\t"
        )
        writer.writeheader()

        for sample in samples:
            source = (
                out / "02_context" / f"{sample.name}.tsv"
            )
            with source.open(newline="") as src:
                reader = csv.DictReader(
                    src, delimiter="\t"
                )
                required = {"candidate_id", "sample"}
                if not required.issubset(
                    reader.fieldnames or []
                ):
                    raise ValueError(
                        f"Invalid Stage 02 columns: {source}"
                    )
                for row in reader:
                    pid = (row["candidate_id"] or "").strip()
                    if not pid:
                        continue
                    if row["sample"] != sample.name:
                        raise ValueError(
                            f"Sample mismatch: {pid}"
                        )
                    writer.writerow({
                        "candidate_id": pid,
                        "sample": sample.name
                    })
                    count += 1

    unresolved_args = []
    for sample in samples:
        table = (
            out / "01_unresolved" / f"{sample.name}.tsv"
        )
        unresolved_args.extend([
            "--unresolved",
            f"{sample.name}={table}"
        ])

    prefix = stage4 / "recurrent_families"

    run("04 recurrent family selection", [
        sys.executable,
        str(ROOT / "scripts/04_select_recurrent_families_portable.py"),
        "--cluster-tsv", str(cluster_tsv),
        "--candidates", str(candidates_file),
        *unresolved_args,
        "--out-prefix", str(prefix)
    ], out / "logs" / "stage04.log")

    for suffix in (
        "_all.tsv", "_ge2.tsv", "_ge3.tsv",
        "_ge4.tsv", "_ge5.tsv"
    ):
        result = Path(str(prefix) + suffix)
        if not result.is_file():
            raise RuntimeError(
                f"Missing Stage 04 output: {result}"
            )

    print(f"STAGE04_CANDIDATES={count}")

    print("STAGES_01_04=COMPLETE")

    if args.through == "04":
        return

    # --------------------------------------------------------
    # Stage 05: recurrent family member annotation
    # --------------------------------------------------------

    family_tsv = stage4 / "recurrent_families_ge2.tsv"

    with family_tsv.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if "representative" not in (reader.fieldnames or []):
            raise ValueError(
                f"Invalid recurrent family table: {family_tsv}"
            )

        selected = sum(
            bool((row["representative"] or "").strip())
            for row in reader
        )

    if selected == 0:
        raise ValueError(
            "No recurrent families in ge2 table. "
            "Stage 05/06 cannot be evaluated."
        )

    print(f"STAGE05_SELECTED_FAMILIES={selected}")

    stage5 = out / "05_annotation"
    stage5.mkdir(parents=True, exist_ok=True)

    query_faa = stage5 / "recurrent_members.faa"
    members_tsv = stage5 / "recurrent_members.tsv"
    search_dir = stage5 / "search"
    search_log = out / "logs" / "stage05_search.log"
    evidence_tsv = stage5 / "member_level_evidence.tsv"

    run("05 extract family members", [
        sys.executable,
        str(ROOT / "scripts/05_extract_recurrent_members.py"),
        "--cluster-tsv", str(cluster_tsv),
        "--family-tsv", str(family_tsv),
        "--proteins-faa", str(combined),
        "--out-faa", str(query_faa),
        "--out-members", str(members_tsv),
    ], out / "logs" / "stage05_members.log")

    if not query_faa.is_file() or query_faa.stat().st_size == 0:
        raise RuntimeError("Empty Stage 05 protein FASTA")

    run("05 functional annotation", [
        "bash",
        str(ROOT / "scripts/05_annotation_saturation_portable.sh"),
        str(query_faa),
        str(search_dir),
        str(db_paths["swissprot"]),
        str(db_paths["phrog"]),
        str(db_paths["pfam"]),
        str(threads),
    ], search_log)

    run("05 build member evidence", [
        sys.executable,
        str(ROOT / "scripts/05_build_member_evidence.py"),
        "--members", str(members_tsv),
        "--proteins-faa", str(query_faa),
        "--search-out", str(search_dir),
        "--search-log", str(search_log),
        "--phrog-annotations",
        str(db_paths["phrog_annotations"]),
        "--out", str(evidence_tsv),
    ], out / "logs" / "stage05_evidence.log")

    if not evidence_tsv.is_file() or evidence_tsv.stat().st_size == 0:
        raise RuntimeError("Stage 05 member evidence missing")

    print(f"STAGE05_EVIDENCE={evidence_tsv}")
    print("STAGES_01_05=COMPLETE")

    if args.through == "05":
        return

    # --------------------------------------------------------
    # Stage 06: family darkness audit
    # --------------------------------------------------------

    stage6 = out / "06_darkness"
    stage6.mkdir(parents=True, exist_ok=True)

    family_audit = stage6 / "family_darkness_audit.tsv"

    run("06 family darkness audit", [
        sys.executable,
        str(ROOT / "scripts/06_family_darkness_audit_portable.py"),
        "--member-evidence", str(evidence_tsv),
        "--out", str(family_audit),
    ], out / "logs" / "stage06.log")

    if not family_audit.is_file() or family_audit.stat().st_size == 0:
        raise RuntimeError("Stage 06 family audit missing")

    print(f"STAGE06_FAMILY_AUDIT={family_audit}")

    print("STAGES_01_06=COMPLETE")

    if args.through == "06":
        return

    # --------------------------------------------------------
    # Stage 07: family-wide genomic neighborhood analysis
    # --------------------------------------------------------

    stage7 = out / "07_family_context"
    stage7.mkdir(parents=True, exist_ok=True)

    run("07 family-wide context", [
        sys.executable,
        str(
            ROOT
            / "scripts/07_familywide_context_portable.py"
        ),
        "--input-dir", str(inp),
        "--stage01-dir", str(out / "01_unresolved"),
        "--stage02-dir", str(out / "02_context"),
        "--members", str(members_tsv),
        "--family-audit", str(family_audit),
        "--out-dir", str(stage7),
        "--window", str(window),
        "--min-aa", str(min_aa),
    ], out / "logs" / "stage07.log")

    member_context = stage7 / "member_context.tsv"
    family_context = stage7 / "family_context.tsv"

    for result in (member_context, family_context):
        if (
            not result.is_file()
            or result.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing Stage 07 output: {result}"
            )

    print(f"STAGE07_MEMBER_CONTEXT={member_context}")
    print(f"STAGE07_FAMILY_CONTEXT={family_context}")
    print("STAGES_01_07=COMPLETE")

    if args.through == "07":
        return

    # Stage 08: prepare, optionally predict, then collect.
    stage8 = out / "08_structure" / f"run_{run_id}"
    prepared = stage8 / "prepared"

    run("08 prepare structure queries", [
        sys.executable,
        str(
            ROOT
            / "scripts/08_prepare_structure_queries_portable.py"
        ),
        "--members", str(members_tsv),
        "--fasta", str(query_faa),
        "--family-audit", str(family_audit),
        "--out-dir", str(prepared),
    ], out / "logs" / f"stage08_prepare_{run_id}.log")

    query_path = prepared / "representative_queries.faa"
    mapping_path = prepared / "structure_id_map.tsv"

    for item in (query_path, mapping_path):
        if not item.is_file() or item.stat().st_size == 0:
            raise RuntimeError(
                f"Missing Stage 08 prepared output: {item}"
            )

    print(f"STAGE08_MODE={args.structure_mode}")
    print(f"STAGE08_QUERY_FASTA={query_path}")
    print(f"STAGE08_STRUCTURE_MAP={mapping_path}")

    if args.structure_mode == "prepare":
        print("STAGE08_STATUS=PREPARED_ONLY")
        print("STAGE08_GPU_PREDICTION=NOT_RUN")
        return

    environment_command = [
        "conda", "run", "--no-capture-output",
        "-n", args.structure_env,
    ]

    if args.structure_mode == "predict":
        prediction_dir = stage8 / "prediction"

        run("08 ColabFold prediction",
            environment_command + [
                "bash",
                str(ROOT / "scripts/08_structure_prediction.sh"),
                str(query_path),
                str(prediction_dir),
            ],
            out / "logs" / f"stage08_predict_{run_id}.log",
        )
    else:
        prediction_dir = resolve(
            args.structure_prediction_dir
        )

    collected = stage8 / "collected"

    run("08 collect structure results",
        environment_command + [
            "python",
            str(
                ROOT
                / "scripts/08_collect_structure_results_portable.py"
            ),
            "--mapping", str(mapping_path),
            "--queries", str(query_path),
            "--prediction-dir", str(prediction_dir),
            "--out-dir", str(collected),
        ],
        out / "logs" / f"stage08_collect_{run_id}.log",
    )

    results_tsv = collected / "structure_results.tsv"

    if (
        not results_tsv.is_file()
        or results_tsv.stat().st_size == 0
    ):
        raise RuntimeError(
            "Stage 08 structure results missing"
        )

    print(f"STAGE08_STRUCTURE_RESULTS={results_tsv}")
    print("STAGES_01_08=COMPLETE")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, RuntimeError) as exc:
        sys.exit(f"ERROR: {exc}")
