#!/usr/bin/env python3
"""Portable Stage I2: reuse original scoring with verified external evidence."""

import argparse
import csv
import hashlib
import json
import math
import runpy
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

ORIGINAL = ROOT / "scripts/10_integrated_ranking.py"
BRIDGE = ROOT / "scripts/10_prepare_evidence_portable.py"

legacy = runpy.run_path(
    str(ORIGINAL), run_name="_legacy_scoring_"
)
bridge = runpy.run_path(
    str(BRIDGE), run_name="_portable_evidence_"
)

calculate = legacy["calculate"]
build_assay_evidence = legacy["build_assay_evidence"]
load_historical_tie_order = legacy["load_historical_tie_order"]
I2_FIELDS = legacy["I2_FIELDS"]
EVIDENCE_FIELDS = bridge["EVIDENCE_FIELDS"]

MEMBER_REQUIRED = {
    "family_rep", "member", "sample",
    "candidate_original", "anchor_within5",
    "anchor_categories", "window",
}


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(
            lambda: fh.read(1048576), b""
        ):
            h.update(chunk)
    return h.hexdigest()


def read_table(path):
    if not path.is_file():
        raise ValueError(f"Missing file: {path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate TSV columns: {path}")

        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty table: {path}")

    return fields, rows


def read_json(path):
    if not path.is_file():
        raise ValueError(f"Missing provenance: {path}")
    return json.loads(path.read_text())


def validate_provenance(args, family_count, member_count):
    b = read_json(args.bridge_provenance)
    s = read_json(args.stage09_provenance)

    if b.get("status") != "EXTERNAL_EVIDENCE_VERIFIED":
        raise ValueError("UNVERIFIED_STAGE10_BRIDGE")

    checks = (
        (
            "evidence_sha256",
            args.evidence,
        ),
        (
            "member_context_sha256",
            args.member_context,
        ),
        (
            "stage09_matrix_sha256",
            args.stage09_matrix,
        ),
        (
            "stage09_provenance_sha256",
            args.stage09_provenance,
        ),
    )

    for key, path in checks:
        if not path.is_file() or b.get(key) != sha(path):
            raise ValueError(
                f"BRIDGE_PROVENANCE_MISMATCH: {key}"
            )

    if (
        b.get("family_count") != family_count
        or b.get("member_count") != member_count
    ):
        raise ValueError("BRIDGE_COUNT_MISMATCH")

    if s.get("evidence_origin") != (
        "COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES"
    ):
        raise ValueError("EXTERNAL_SEARCH_NOT_VERIFIED")

    if s.get("matrix_sha256") != sha(
        args.stage09_matrix
    ):
        raise ValueError("STAGE09_MATRIX_SHA_MISMATCH")

    if s.get("family_count") != family_count:
        raise ValueError("STAGE09_FAMILY_COUNT_MISMATCH")

    if set(s.get("external_searches", {})) != {
        "pdb", "afsp", "bfvd"
    }:
        raise ValueError("INCOMPLETE_EXTERNAL_DATABASES")

    matrix_fields, matrix = read_table(
        args.stage09_matrix
    )

    if len(matrix_fields) != 39:
        raise ValueError("INVALID_STAGE09_SCHEMA")

    matrix_by_id = {
        r["structure_id"]: r
        for r in matrix
    }

    if (
        len(matrix_by_id) != family_count
        or len(matrix) != family_count
    ):
        raise ValueError("STAGE09_STRUCTURE_ID_MISMATCH")

    return matrix_by_id


def validate_inputs(evidence, members, matrix):
    ids = [r["structure_id"] for r in evidence]
    reps = [r["family_rep"] for r in evidence]

    if (
        len(set(ids)) != len(ids)
        or len(set(reps)) != len(reps)
    ):
        raise ValueError("DUPLICATE_FAMILY_OR_STRUCTURE")

    if set(ids) != set(matrix):
        raise ValueError("STAGE09_STAGE10_ID_MISMATCH")

    member_ids = set()
    counts = Counter()
    samples = {}

    for m in members:
        rep = m["family_rep"]
        member = m["member"]

        if not member or member in member_ids:
            raise ValueError("DUPLICATE_MEMBER")

        member_ids.add(member)
        counts[rep] += 1
        samples.setdefault(rep, set()).add(m["sample"])

        if m["window"] != "5":
            raise ValueError("HISTORICAL_WINDOW5_REQUIRED")

        if m["candidate_original"] not in ("0", "1"):
            raise ValueError("INVALID_CANDIDATE_STATUS")

        if m["anchor_within5"] not in ("0", "1"):
            raise ValueError("INVALID_ANCHOR_STATUS")

    if set(counts) != set(reps):
        raise ValueError("FAMILY_MEMBER_SET_MISMATCH")

    for r in evidence:
        sid = r["structure_id"]
        rep = r["family_rep"]
        original = matrix[sid]

        if counts[rep] != int(r["family_size"]):
            raise ValueError(f"FAMILY_SIZE_MISMATCH: {rep}")

        if len(samples[rep]) != int(
            r["family_sample_span"]
        ):
            raise ValueError(f"SAMPLE_SPAN_MISMATCH: {rep}")

        for key in (
            "priority_class",
            "external_status",
            "raw_sources",
            "convincing_sources",
        ):
            if r[key] != original[key]:
                raise ValueError(
                    f"EXTERNAL_EVIDENCE_MISMATCH: {sid}:{key}"
                )

        if int(r["length_aa"]) != int(
            original["length"]
        ):
            raise ValueError("STRUCTURE_LENGTH_MISMATCH")

        for key in (
            "mean_plddt", "frac_ge70", "ptm"
        ):
            value = float(r[key])
            reference = float(original[key])

            if (
                not math.isfinite(value)
                or abs(value - reference) > 0.0001
            ):
                raise ValueError(
                    f"STRUCTURE_METRIC_MISMATCH: {sid}:{key}"
                )

        for key in (
            "frac_ge70", "structural_dense_fraction"
        ):
            value = float(r[key])
            if not 0 <= value <= 1:
                raise ValueError(
                    f"INVALID_FRACTION: {sid}:{key}"
                )


def rank(evidence, members, tie_order=None):
    counts, samples = build_assay_evidence(members)

    rows = [
        calculate(r, counts, samples)
        for r in evidence
    ]

    historical = (
        load_historical_tie_order(
            tie_order, rows
        )
        if tie_order is not None
        else None
    )

    integrated = sorted(
        rows,
        key=lambda r: (
            -float(r["I2_integrated_score"]),
            -int(
                r["noncandidate_context_sample_span"]
            ),
            (
                historical.get(
                    r["structure_id"], 0
                )
                if historical is not None
                else r["structure_id"]
            ),
        ),
    )

    for i, row in enumerate(integrated, 1):
        row["I2_integrated_rank"] = str(i)

    experimental = sorted(
        rows,
        key=lambda r: (
            -float(r["I2_experimental_score"]),
            -float(r["I2_integrated_score"]),
            -int(
                r["I2_assay_independent_sample_span"]
            ),
            r["structure_id"],
        ),
    )

    for i, row in enumerate(experimental, 1):
        row["I2_experimental_rank"] = str(i)

    return integrated


def self_test():
    def evidence(sid, rep, priority, external, plddt,
                 fraction, span, context_span,
                 context_count, length):
        return {
            "structure_id": sid,
            "family_rep": rep,
            "priority_class": priority,
            "family_size": "1",
            "family_sample_span": str(span),
            "length_aa": str(length),
            "mean_plddt": str(plddt),
            "ptm": "0.5",
            "frac_ge70": str(fraction),
            "external_status": external,
            "raw_sources": "",
            "convincing_sources": "",
            "noncandidate_context_count": str(
                context_count
            ),
            "noncandidate_context_sample_span": str(
                context_span
            ),
            "structural_dense_fraction": "0.05",
        }

    a = evidence(
        "T001", "R1", "A_DEEP_DARK",
        "NO_RAW_HIT_E1", 80, 0.8,
        3, 2, 2, 200,
    )

    b = evidence(
        "T002", "R2", "D_KNOWN_EVIDENCE",
        "CONVINCING_STRUCTURE_CONNECTION",
        34, 0, 2, 0, 0, 400,
    )

    members = [
        {
            "family_rep": "R1",
            "member": "M1",
            "sample": "S1",
            "candidate_original": "0",
            "anchor_within5": "1",
            "anchor_categories": "nuclease",
        },
        {
            "family_rep": "R2",
            "member": "M2",
            "sample": "S2",
            "candidate_original": "1",
            "anchor_within5": "0",
            "anchor_categories": "",
        },
    ]

    rows = rank([b, a], members)

    by_id = {
        r["structure_id"]: r for r in rows
    }

    assert [
        r["structure_id"] for r in rows
    ] == ["T001", "T002"]

    assert by_id["T001"]["I2_dark"] == "4.00"
    assert by_id["T001"]["I2_external"] == "4.00"
    assert by_id["T001"]["I2_structure"] == "3.00"
    assert by_id["T001"]["I2_context"] == "2.00"
    assert by_id["T001"]["I2_recurrence"] == "0.50"

    assert (
        by_id["T001"]["I2_integrated_score"]
        == "13.50"
    )
    assert by_id["T001"]["I2_assay_bonus"] == "0.60"
    assert (
        by_id["T001"]["I2_experimental_score"]
        == "15.10"
    )

    assert by_id["T002"]["I2_integrated_score"] == "0.00"
    assert by_id["T002"]["I2_experimental_score"] == "0.50"

    assert by_id["T001"]["I2_integrated_rank"] == "1"
    assert by_id["T002"]["I2_integrated_rank"] == "2"

    assert by_id["T001"]["I2_experimental_rank"] == "1"
    assert by_id["T002"]["I2_experimental_rank"] == "2"

    assert len(I2_FIELDS) == 13
    assert len(EVIDENCE_FIELDS) == 15

    print("ORIGINAL_I2_SCORE_FUNCTIONS=REUSED")
    print("DARKNESS_AND_EXTERNAL_SCORING=PASS")
    print("STRUCTURE_CONTEXT_RECURRENCE_SCORING=PASS")
    print("ASSAY_AND_EXPRESSION_SCORING=PASS")
    print("DYNAMIC_FAMILY_RANKING=PASS")
    print("STAGE10_PORTABLE_UNIT_TEST=PASS")
    print("ACTUAL_RESEARCH_RANKING=NOT_RUN")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--evidence", type=Path)
    ap.add_argument("--member-context", type=Path)
    ap.add_argument("--bridge-provenance", type=Path)
    ap.add_argument("--stage09-matrix", type=Path)
    ap.add_argument("--stage09-provenance", type=Path)
    ap.add_argument("--historical-tie-order", type=Path)
    ap.add_argument("--out-dir", type=Path)

    args = ap.parse_args()

    if args.self_test:
        self_test()
        return

    required = (
        args.evidence,
        args.member_context,
        args.bridge_provenance,
        args.stage09_matrix,
        args.stage09_provenance,
        args.out_dir,
    )

    if any(x is None for x in required):
        ap.error(
            "Ranking requires evidence, member context, "
            "both provenance records, Stage09 matrix "
            "and output directory"
        )

    output = args.out_dir.resolve()

    if output.exists():
        raise ValueError(
            "Refusing existing output directory"
        )

    fields, evidence_rows = read_table(args.evidence)
    member_fields, member_rows = read_table(
        args.member_context
    )

    if fields != EVIDENCE_FIELDS:
        raise ValueError("STAGE10_EVIDENCE_SCHEMA_MISMATCH")

    if not MEMBER_REQUIRED.issubset(member_fields):
        raise ValueError("STAGE10_MEMBER_SCHEMA_MISMATCH")

    matrix = validate_provenance(
        args,
        len(evidence_rows),
        len(member_rows),
    )

    validate_inputs(
        evidence_rows,
        member_rows,
        matrix,
    )

    rows = rank(
        evidence_rows,
        member_rows,
        args.historical_tie_order,
    )

    output.mkdir(parents=True)

    result = output / "confidence_aware_ranking.tsv"

    header = [
        "I2_integrated_rank",
        "I2_experimental_rank",
    ] + fields + I2_FIELDS

    with result.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    record = {
        "algorithm": "ORIGINAL_STAGE_I2",
        "scoring_script_sha256": sha(ORIGINAL),
        "evidence_sha256": sha(args.evidence),
        "member_context_sha256": sha(
            args.member_context
        ),
        "bridge_provenance_sha256": sha(
            args.bridge_provenance
        ),
        "stage09_provenance_sha256": sha(
            args.stage09_provenance
        ),
        "historical_tie_replay": (
            args.historical_tie_order is not None
        ),
        "families": len(rows),
        "columns": len(header),
        "ranking_sha256": sha(result),
    }

    (output / "ranking_provenance.json").write_text(
        json.dumps(
            record, ensure_ascii=False, indent=2
        ) + "\n"
    )

    print(f"RANKED_FAMILIES={len(rows)}")
    print(f"RANKING_COLUMNS={len(header)}")
    print(f"RANKING_OUTPUT={result}")
    print("STAGE10_PORTABLE_RANKING=PASS")


if __name__ == "__main__":
    main()
