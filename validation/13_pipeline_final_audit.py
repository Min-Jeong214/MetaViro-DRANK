#!/usr/bin/env python3
"""Final MetaViro source and canonical-output integrity audit."""

import csv
import hashlib
import re
import subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Result/Annotation_saturation/Pilot_v2"

OUTPUTS = [
    (
        "07_MEMBER",
        "StageH3_CorrectedFamilyWideContext/"
        "recurrent191_member_context_corrected.tsv",
        "Stage07_Reconstruction_DryRun/"
        "recurrent191_member_context_corrected.tsv",
        10327, 19,
    ),
    (
        "07_FAMILY",
        "StageH3_CorrectedFamilyWideContext/"
        "recurrent191_family_context_corrected.tsv",
        "Stage07_Reconstruction_DryRun/"
        "recurrent191_family_context_corrected.tsv",
        191, 24,
    ),
    (
        "09_EXTERNAL",
        "StageF_Structure/StageG3_ExternalDB/"
        "Integrated_ExternalStructure/"
        "structure191_external_structure_matrix.tsv",
        "Stage09_Reconstruction_DryRun/"
        "structure191_external_structure_matrix.tsv",
        191, 39,
    ),
    (
        "10_MAIN",
        "StageI2_ConfidenceAwareRanking/"
        "recurrent191_confidence_aware_ranking.tsv",
        "Stage10_Reconstruction_DryRun/"
        "recurrent191_confidence_aware_ranking.tsv",
        191, 51,
    ),
    (
        "10_TOP20",
        "StageI2_ConfidenceAwareRanking/"
        "top20_confidence_aware.tsv",
        "Stage10_Reconstruction_DryRun/"
        "top20_confidence_aware.tsv",
        20, 51,
    ),
    (
        "10_TOP15",
        "StageI2_ConfidenceAwareRanking/"
        "top15_experimental_candidates.tsv",
        "Stage10_Reconstruction_DryRun/"
        "top15_experimental_candidates.tsv",
        15, 51,
    ),
]

def read_tsv(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

print("===== 1. STAGE 01-13 SOURCE COVERAGE =====")

scripts_dir = ROOT / "scripts"
stage_scripts = defaultdict(list)

pattern = re.compile(r"^(\d{2})_.+\.(py|sh)$")

for path in sorted(scripts_dir.iterdir()):
    if not path.is_file():
        continue

    match = pattern.fullmatch(path.name)
    if not match:
        continue

    stage = int(match.group(1))
    if 1 <= stage <= 13:
        stage_scripts[stage].append(path)

coverage_ok = True
syntax_ok = True

for stage in range(1, 14):
    files = stage_scripts[stage]

    if not files:
        coverage_ok = False
        print(f"STAGE{stage:02d}=MISSING")
        continue

    print(
        f"STAGE{stage:02d}=FOUND",
        ",".join(p.name for p in files)
    )

    for path in files:
        try:
            if path.suffix == ".py":
                compile(
                    path.read_bytes(),
                    str(path),
                    "exec"
                )
            elif path.suffix == ".sh":
                subprocess.run(
                    ["bash", "-n", str(path)],
                    check=True,
                    capture_output=True,
                    text=True
                )
        except Exception as e:
            syntax_ok = False
            print("SYNTAX_FAIL=", path, str(e))

print("SCRIPT_STAGE_COVERAGE=", "PASS" if coverage_ok else "FAIL")
print("SCRIPT_SYNTAX=", "PASS" if syntax_ok else "FAIL")

print("\n===== 2. CANONICAL OUTPUT AUDIT =====")

output_ok = True

for label, gold_rel, calc_rel, nrows, ncols in OUTPUTS:
    gold_path = BASE / gold_rel
    calc_path = BASE / calc_rel

    print(f"\nOUTPUT={label}")

    if not gold_path.is_file() or not calc_path.is_file():
        output_ok = False
        print("STATUS=FAIL_MISSING_FILE")
        print("CANON_EXISTS=", gold_path.is_file())
        print("CALC_EXISTS=", calc_path.is_file())
        continue

    gold_cols, gold_rows = read_tsv(gold_path)
    calc_cols, calc_rows = read_tsv(calc_path)

    schema = (
        len(gold_rows) == len(calc_rows) == nrows
        and len(gold_cols) == len(calc_cols) == ncols
        and gold_cols == calc_cols
    )

    byte_exact = (
        gold_path.read_bytes() == calc_path.read_bytes()
    )

    print("EXPECTED_SHAPE=", (nrows, ncols))
    print("CANON_SHAPE=", (len(gold_rows), len(gold_cols)))
    print("CALC_SHAPE=", (len(calc_rows), len(calc_cols)))
    print("SCHEMA_EXACT=", schema)
    print("BYTE_EXACT=", byte_exact)
    print("SHA256=", digest(calc_path))

    passed = schema and byte_exact
    output_ok &= passed

    print("STATUS=", "PASS" if passed else "FAIL")

print("\n===== 3. HISTORICAL TIE MANIFEST =====")

manifest = (
    ROOT / "scripts/resources/"
    "10_i2_historical_tie_order.tsv"
)

main_canon = (
    BASE / "StageI2_ConfidenceAwareRanking/"
    "recurrent191_confidence_aware_ranking.tsv"
)

tie_ok = False

if manifest.is_file() and main_canon.is_file():
    columns, tie_rows = read_tsv(manifest)
    _, canon_rows = read_tsv(main_canon)

    canon_rows.sort(
        key=lambda r: int(r["I2_integrated_rank"])
    )

    groups = defaultdict(list)

    for r in canon_rows:
        key = (
            r["I2_integrated_score"],
            r["noncandidate_context_sample_span"],
        )
        groups[key].append(r["structure_id"])

    expected = {}

    for key, ids in groups.items():
        if len(ids) <= 1:
            continue
        for index, sid in enumerate(ids, 1):
            expected[sid] = (key, index)

    actual_ids = set()
    mismatches = 0

    for r in tie_rows:
        sid = r["structure_id"]

        if sid in actual_ids or sid not in expected:
            mismatches += 1
            continue

        actual_ids.add(sid)

        key = (
            r["I2_integrated_score"],
            r["noncandidate_context_sample_span"],
        )

        if (
            key != expected[sid][0]
            or int(r["historical_tie_order"])
            != expected[sid][1]
        ):
            mismatches += 1

    tie_ok = (
        len(tie_rows) == 168
        and len(expected) == 168
        and len([
            ids for ids in groups.values()
            if len(ids) > 1
        ]) == 34
        and actual_ids == set(expected)
        and mismatches == 0
    )

    print("MANIFEST_ROWS=", len(tie_rows))
    print("EXPECTED_TIED_MEMBERS=", len(expected))
    print("MISMATCHES=", mismatches)
else:
    print("MANIFEST_OR_CANONICAL_MISSING")

print("HISTORICAL_MANIFEST=", "PASS" if tie_ok else "FAIL")

print("\n===== 4. GITHUB CANDIDATE SCRIPTS =====")

for stage in range(1, 14):
    for path in stage_scripts[stage]:
        print("SOURCE=", path.relative_to(ROOT))

for name in [
    "scripts/resources/10_i2_historical_tie_order.tsv",
    "scripts/resources/10_i2_historical_tie_order.README.md",
    "validation/09_matrix_exact_audit.py",
    "validation/10_full_i2_outputs_audit.py",
    "validation/13_pipeline_final_audit.py",
]:
    path = ROOT / name
    print(
        "SUPPORT=",
        name,
        "EXISTS=", path.is_file()
    )

print("\n===== FINAL =====")

passed = (
    coverage_ok
    and syntax_ok
    and output_ok
    and tie_ok
)

print(
    "CANONICAL_OUTPUTS_STATUS=",
    "PASS" if output_ok else "FAIL"
)
print(
    "FINAL_PIPELINE_AUDIT=",
    "PASS" if passed else "FAIL"
)

print("CANONICAL_DATA_UNMODIFIED=YES")
print("PRODUCTION_SCRIPTS_UNMODIFIED=YES")
print("GIT_COMMIT_OR_PUSH=NOT_EXECUTED")

print(
    "NOTE=Stages 01-06, 08, 11-13 were "
    "syntax-checked but not re-executed in this audit"
)

if not passed:
    raise SystemExit(1)
