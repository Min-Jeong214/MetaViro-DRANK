#!/usr/bin/env python3
"""Portable Stage G3 integration with mandatory complete-search evidence."""

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

HIT_FIELDS = (
    "query target evalue bits fident alnlen qlen tlen "
    "qcov tcov lddt alntmscore qtmscore ttmscore prob"
).split()

SOURCES = (
    ("pdb", "PDB", "pdb_e1", 1000),
    ("afsp", "AFSP", "afsp", 1000),
    ("bfvd", "BFVD", "bfvd", 500),
)

CORE = (
    "structure_id priority_class length mean_plddt ptm "
    "frac_ge70 external_status raw_sources convincing_sources"
).split()

DETAIL = (
    "raw_count convincing_count best_target best_prob "
    "best_qtm best_ttm best_lddt best_qcov best_tcov "
    "best_evalue"
).split()

HEADER = CORE + [
    f"{prefix}_{field}"
    for prefix, _, _, _ in SOURCES
    for field in DETAIL
]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def tsv(path, required=()):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"Missing columns in {path}")
        return list(reader)


def numeric(value, label):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite numeric field: {label}")
    return result


def convincing(hit):
    return (
        numeric(hit["prob"], "prob") >= 0.90
        and numeric(hit["lddt"], "lddt") >= 0.60
        and numeric(hit["qcov"], "qcov") >= 0.30
        and numeric(hit["tcov"], "tcov") >= 0.50
        and numeric(hit["ttmscore"], "ttmscore") >= 0.50
    )


def best_hit(hits):
    strong = [hit for hit in hits if convincing(hit)]
    pool = strong if strong else hits
    return max(
        pool,
        key=lambda h: (
            float(h["prob"]),
            max(float(h["qtmscore"]), float(h["ttmscore"])),
        ),
    )


def pdb_confidence(path, expected_length, expected_mean, score_json=None, score_json_sha256=None):
    scores = []
    residues = set()

    with path.open(errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            if line[12:16].strip() != "CA":
                continue
            if line[16] not in (" ", "A"):
                continue

            residue = (line[21], line[22:27])

            if residue in residues:
                raise ValueError(f"Duplicate CA residue: {path}")

            residues.add(residue)
            scores.append(
                numeric(line[60:66], "PDB CA pLDDT")
            )

    if len(scores) != expected_length:
        if not score_json or not score_json_sha256:
            raise ValueError(
                f"Missing JSON provenance for incomplete PDB: {path}"
            )

        json_path = Path(score_json)

        if (
            not json_path.is_file()
            or sha(json_path) != score_json_sha256
        ):
            raise ValueError(
                f"Prediction JSON integrity mismatch: {path}"
            )

        values = json.loads(
            json_path.read_text(encoding="utf-8")
        ).get("plddt")

        if not isinstance(values, list) or len(values) != expected_length:
            raise ValueError(
                f"Prediction JSON residue count mismatch: {path}"
            )

        full_scores = [
            numeric(v, "Prediction JSON pLDDT")
            for v in values
        ]

        if any(not 0 <= v <= 100 for v in full_scores):
            raise ValueError(f"Invalid JSON pLDDT: {path}")

        positions = set()
        chains = set()

        with path.open(errors="replace") as fh:
            for line in fh:
                if (
                    not line.startswith("ATOM  ")
                    or line[12:16].strip() != "CA"
                    or line[16] not in (" ", "A")
                ):
                    continue

                chains.add(line[21])
                position = int(line[22:26])

                if (
                    position < 1
                    or position > expected_length
                    or position in positions
                    or line[26] != " "
                ):
                    raise ValueError(
                        f"Invalid PDB residue mapping: {path}"
                    )

                positions.add(position)

                observed = numeric(
                    line[60:66], "PDB CA pLDDT"
                )

                if abs(
                    observed - full_scores[position - 1]
                ) > 0.051:
                    raise ValueError(
                        f"PDB/JSON residue mismatch: {path}"
                    )

        if len(chains) != 1 or len(positions) != len(scores):
            raise ValueError(
                f"Inconsistent PDB residue mapping: {path}"
            )

        scores = full_scores

    mean = sum(scores) / len(scores)

    if abs(mean - expected_mean) > 0.05:
        raise ValueError(
            f"PDB/Stage08 confidence mismatch: {path}"
        )

    return sum(v >= 70 for v in scores) / len(scores)


def load_structures(path):
    required = (
        "structure_id", "priority_class", "length_aa",
        "mean_plddt", "ptm", "pdb_file",
        "pdb_sha256", "prediction_status",
    )
    rows = tsv(path, required)

    if not rows:
        raise ValueError("No Stage 08 structures")

    observed = set()
    metadata = {}

    for row in rows:
        sid = row["structure_id"]

        if not sid or sid in observed:
            raise ValueError(f"Duplicate/empty structure ID: {sid}")

        observed.add(sid)

        if row["prediction_status"] != "PREDICTED":
            raise ValueError(f"Unpredicted structure: {sid}")

        pdb = Path(row["pdb_file"]).resolve()

        if not pdb.is_file() or sha(pdb) != row["pdb_sha256"]:
            raise ValueError(f"PDB integrity mismatch: {sid}")

        length = int(row["length_aa"])
        mean = numeric(row["mean_plddt"], "mean_plddt")
        ptm = numeric(row["ptm"], "ptm")

        if length < 1 or not (0 <= ptm <= 1):
            raise ValueError(f"Invalid Stage 08 metrics: {sid}")

        fraction = pdb_confidence(
            pdb, length, mean,
            row.get("score_json"),
            row.get("score_json_sha256"),
        )

        metadata[sid] = {
            "structure_id": sid,
            "priority_class": row["priority_class"],
            "length": str(length),
            "mean_plddt": f"{mean:.2f}",
            "ptm": f"{ptm:.4f}",
            "frac_ge70": f"{fraction:.4f}",
            "pdb_sha256": row["pdb_sha256"],
        }

    return metadata


def load_search(directory, expected_profile, max_seqs, metadata):
    directory = directory.resolve()
    info_path = directory / "search_provenance.json"
    manifest_path = directory / "query_manifest.tsv"
    hits_path = directory / "raw_hits.tsv"

    if not all(
        p.is_file()
        for p in (info_path, manifest_path, hits_path)
    ):
        raise ValueError(
            f"Incomplete Foldseek search directory: {directory}"
        )

    info = json.loads(info_path.read_text())

    if info.get("control_only") is not False:
        raise ValueError(
            f"Control/unverified search cannot be external evidence: "
            f"{directory}"
        )

    if info.get("source_profile") != expected_profile:
        raise ValueError(
            f"Incorrect Foldseek search profile: {directory}"
        )

    if (
        numeric(info.get("evalue"), "evalue") != 1.0
        or numeric(info.get("sensitivity"), "sensitivity") != 9.5
        or int(info.get("max_seqs")) != max_seqs
    ):
        raise ValueError(
            f"Historical search parameter mismatch: {directory}"
        )

    if sha(hits_path) != info.get("raw_hits_sha256"):
        raise ValueError(f"Raw hit SHA mismatch: {directory}")

    manifest = tsv(
        manifest_path,
        ("structure_id", "sha256"),
    )

    found = {}

    for item in manifest:
        sid = item["structure_id"]
        if sid in found:
            raise ValueError(f"Duplicate query in manifest: {sid}")
        found[sid] = item["sha256"]

    if set(found) != set(metadata):
        raise ValueError(
            f"Incomplete query coverage: {directory}"
        )

    for sid in metadata:
        if found[sid] != metadata[sid]["pdb_sha256"]:
            raise ValueError(
                f"Query PDB differs between stages: {sid}"
            )

    hits = defaultdict(list)

    with hits_path.open(newline="") as fh:
        for line_number, parts in enumerate(
            csv.reader(fh, delimiter="\t"), 1
        ):
            if not parts:
                continue

            if len(parts) != len(HIT_FIELDS):
                raise ValueError(
                    f"Malformed Foldseek row {line_number}: "
                    f"{directory}"
                )

            hit = dict(zip(HIT_FIELDS, parts))
            sid = hit["query"]

            if sid not in metadata:
                raise ValueError(
                    f"Unknown Foldseek query: {sid}"
                )

            for field in (
                "evalue", "prob", "lddt",
                "qcov", "tcov", "qtmscore", "ttmscore",
            ):
                value = numeric(hit[field], field)

                if field != "evalue" and not 0 <= value <= 1:
                    raise ValueError(
                        f"Out-of-range Foldseek {field}: {value}"
                    )

            hits[sid].append(hit)

    return hits, {
        "directory": str(directory),
        "source_profile": expected_profile,
        "raw_hits_sha256": sha(hits_path),
        "search_provenance_sha256": sha(info_path),
        "query_manifest_sha256": sha(manifest_path),
        "foldseek_version": info.get("foldseek_version"),
    }


def integrate(metadata, sources):
    output = []

    for sid, meta in metadata.items():
        row = {
            key: meta[key]
            for key in CORE[:6]
        }

        raw_labels = []
        convincing_labels = []

        for prefix, label, _, _ in SOURCES:
            hits = sources[prefix].get(sid, [])
            strong = [hit for hit in hits if convincing(hit)]

            row[f"{prefix}_raw_count"] = str(len(hits))
            row[f"{prefix}_convincing_count"] = str(len(strong))

            if hits:
                raw_labels.append(label)
                best = best_hit(hits)

                for field, source in (
                    ("target", "target"),
                    ("prob", "prob"),
                    ("qtm", "qtmscore"),
                    ("ttm", "ttmscore"),
                    ("lddt", "lddt"),
                    ("qcov", "qcov"),
                    ("tcov", "tcov"),
                    ("evalue", "evalue"),
                ):
                    row[f"{prefix}_best_{field}"] = best[source]
            else:
                for field in (
                    "target", "prob", "qtm", "ttm",
                    "lddt", "qcov", "tcov", "evalue"
                ):
                    row[f"{prefix}_best_{field}"] = ""

            if strong:
                convincing_labels.append(label)

        row["raw_sources"] = ",".join(raw_labels)
        row["convincing_sources"] = ",".join(
            convincing_labels
        )

        if convincing_labels:
            row["external_status"] = (
                "CONVINCING_STRUCTURE_CONNECTION"
            )
        elif raw_labels:
            row["external_status"] = "WEAK_RAW_HIT_ONLY"
        else:
            row["external_status"] = "NO_RAW_HIT_E1"

        if set(row) != set(HEADER):
            raise ValueError("Stage G3 output schema mismatch")

        output.append(row)

    return output


def self_test():
    assert len(HEADER) == 39

    sid = "TEST001"

    meta = {
        sid: {
            "structure_id": sid,
            "priority_class": "A_DEEP_DARK",
            "length": "160",
            "mean_plddt": "80.00",
            "ptm": "0.7000",
            "frac_ge70": "0.8000",
            "pdb_sha256": "SIMULATED",
        }
    }

    empty = {
        prefix: defaultdict(list)
        for prefix, _, _, _ in SOURCES
    }

    result = integrate(meta, empty)

    assert result[0]["external_status"] == "NO_RAW_HIT_E1"
    assert result[0]["raw_sources"] == ""
    assert len(result[0]) == 39

    print("SIMULATED_COMPLETE_NO_HIT_RULE=PASS")

    weak = dict.fromkeys(HIT_FIELDS, "1")
    weak.update({
        "query": sid,
        "target": "SIMULATED_WEAK",
        "prob": "0.80",
        "lddt": "0.80",
        "qcov": "0.80",
        "tcov": "0.80",
        "ttmscore": "0.80",
        "qtmscore": "0.80",
    })

    empty["pdb"][sid] = [weak]

    result = integrate(meta, empty)

    assert result[0]["external_status"] == "WEAK_RAW_HIT_ONLY"
    assert result[0]["pdb_raw_count"] == "1"
    assert result[0]["pdb_convincing_count"] == "0"

    print("SIMULATED_WEAK_HIT_RULE=PASS")

    strong = {
        **weak,
        "target": "SIMULATED_STRONG",
        "prob": "0.95",
    }

    empty["bfvd"][sid] = [strong]

    result = integrate(meta, empty)

    assert (
        result[0]["external_status"]
        == "CONVINCING_STRUCTURE_CONNECTION"
    )
    assert result[0]["convincing_sources"] == "BFVD"
    assert result[0]["bfvd_best_target"] == (
        "SIMULATED_STRONG"
    )

    print("SIMULATED_CONVINCING_HIT_RULE=PASS")
    print("HISTORICAL_39_COLUMN_SCHEMA=PASS")
    print("STAGE09_INTEGRATION_UNIT_TEST=PASS")
    print("SIMULATED_RESULTS_NOT_WRITTEN=YES")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--structures", type=Path)
    parser.add_argument("--pdb-e1-dir", type=Path)
    parser.add_argument("--afsp-dir", type=Path)
    parser.add_argument("--bfvd-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    required = (
        args.structures, args.pdb_e1_dir,
        args.afsp_dir, args.bfvd_dir, args.out_dir,
    )

    if any(x is None for x in required):
        parser.error("All five input/output paths are required")

    dirs = {
        "pdb": args.pdb_e1_dir.resolve(),
        "afsp": args.afsp_dir.resolve(),
        "bfvd": args.bfvd_dir.resolve(),
    }

    if len(set(dirs.values())) != 3:
        raise ValueError("Source search directories must differ")

    out = args.out_dir.resolve()

    if out.exists():
        raise ValueError(
            f"Refusing existing output directory: {out}"
        )

    if any(
        out == d or d in out.parents or out in d.parents
        for d in dirs.values()
    ):
        raise ValueError("Output overlaps source search directory")

    metadata = load_structures(args.structures)

    sources = {}
    provenance = {}

    for prefix, _, profile, max_seqs in SOURCES:
        sources[prefix], provenance[prefix] = load_search(
            dirs[prefix], profile, max_seqs, metadata
        )

    rows = integrate(metadata, sources)

    # Do not create any output until all evidence is validated.
    out.mkdir(parents=True)

    matrix = out / "external_structure_matrix.tsv"

    with matrix.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=HEADER,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    record = {
        "evidence_origin": "COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES",
        "historical_evidence_rules": True,
        "structure_table": str(args.structures.resolve()),
        "structure_table_sha256": sha(args.structures),
        "external_searches": provenance,
        "family_count": len(rows),
        "status_counts": dict(
            Counter(row["external_status"] for row in rows)
        ),
        "matrix_sha256": sha(matrix),
    }

    (out / "integration_provenance.json").write_text(
        json.dumps(
            record, indent=2, ensure_ascii=False
        ) + "\n"
    )

    print(f"FAMILIES={len(rows)}")
    print(f"COLUMNS={len(HEADER)}")
    print(f"EXTERNAL_STATUS_COUNTS={record['status_counts']}")
    print(f"OUTPUT={matrix}")
    print("STAGE09_PORTABLE_INTEGRATION=PASS")


if __name__ == "__main__":
    main()
