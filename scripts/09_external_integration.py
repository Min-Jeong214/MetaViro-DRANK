#!/usr/bin/env python3
"""Reconstruct the Stage G3 191-family external structure matrix.

Inputs: Stage F structure metadata and three raw Foldseek outputs.
The canonical integrated matrix is NOT an input.
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

HIT_FIELDS = (
    "query target evalue bits fident alnlen qlen tlen "
    "qcov tcov lddt alntmscore qtmscore ttmscore prob"
).split()

SOURCES = (
    ("pdb", "PDB", "PDB_Search_E1/pdb_hits_e1.tsv"),
    ("afsp", "AFSP", "AF_SwissProt_Search/af_swissprot_hits.tsv"),
    ("bfvd", "BFVD", "BFVD_Search/bfvd_hits.tsv"),
)

CORE_COLUMNS = (
    "structure_id priority_class length mean_plddt ptm frac_ge70 "
    "external_status raw_sources convincing_sources"
).split()

SOURCE_COLUMNS = (
    "raw_count convincing_count best_target best_prob best_qtm "
    "best_ttm best_lddt best_qcov best_tcov best_evalue"
).split()

HEADER = CORE_COLUMNS + [
    f"{prefix}_{field}"
    for prefix, _, _ in SOURCES
    for field in SOURCE_COLUMNS
]


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def passing(hit):
    return (
        float(hit["prob"]) >= 0.90
        and float(hit["lddt"]) >= 0.60
        and float(hit["qcov"]) >= 0.30
        and float(hit["tcov"]) >= 0.50
        and float(hit["ttmscore"]) >= 0.50
    )


def select_best(hits):
    strong = [h for h in hits if passing(h)]
    pool = strong if strong else hits

    return max(
        pool,
        key=lambda h: (
            float(h["prob"]),
            max(
                float(h["qtmscore"]),
                float(h["ttmscore"])
            )
        )
    )


def format_num(value, decimals):
    return f"{float(value):.{decimals}f}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(".")
    )
    parser.add_argument(
        "--out", type=Path, required=True
    )
    args = parser.parse_args()

    base = (
        args.root / "Result/Annotation_saturation/Pilot_v2"
    )
    g3 = base / "StageF_Structure/StageG3_ExternalDB"

    mapping = read_tsv(
        base / "StageF_Structure/structure_id_map.tsv"
    )
    qc = read_tsv(
        base / "StageF_Structure/structure191_qc.tsv"
    )

    qcm = {x["structure_id"]: x for x in qc}
    ids = [r["structure_id"] for r in mapping]

    if not (
        len(ids) == len(set(ids)) == len(qcm) == 191
        and set(ids) == set(qcm)
    ):
        raise ValueError(
            "Structure QC/map must contain 191 matching unique IDs"
        )

    by_source = {}
    totals = Counter()

    for prefix, label, suffix in SOURCES:
        hits = defaultdict(list)

        with (g3 / suffix).open(newline="") as fh:
            for lineno, parts in enumerate(
                csv.reader(fh, delimiter="\t"), 1
            ):
                if not parts:
                    continue

                if len(parts) != len(HIT_FIELDS):
                    raise ValueError(
                        f"{suffix}:{lineno}: {len(parts)} fields"
                    )

                hit = dict(zip(HIT_FIELDS, parts))

                if hit["query"] not in qcm:
                    raise ValueError(
                        f"Unknown Foldseek query {hit['query']}"
                    )

                hits[hit["query"]].append(hit)
                totals[(prefix, "raw")] += 1
                totals[(prefix, "convincing")] += int(
                    passing(hit)
                )

        by_source[prefix] = hits

    output = []

    for meta in mapping:
        sid = meta["structure_id"]
        quality = qcm[sid]

        row = {
            "structure_id": sid,
            "priority_class": quality["priority_class"],
            "length": meta["length"],
            "mean_plddt": format_num(
                quality["mean_plddt"], 2
            ),
            "ptm": format_num(quality["ptm"], 4),
            "frac_ge70": format_num(
                quality["frac_ge70"], 4
            ),
        }

        raw_labels = []
        convincing_labels = []

        for prefix, label, _ in SOURCES:
            hits = by_source[prefix][sid]
            n_strong = sum(passing(h) for h in hits)

            row[f"{prefix}_raw_count"] = str(len(hits))
            row[f"{prefix}_convincing_count"] = str(
                n_strong
            )

            if hits:
                raw_labels.append(label)
                best = select_best(hits)

                for key, source in (
                    ("target", "target"),
                    ("prob", "prob"),
                    ("qtm", "qtmscore"),
                    ("ttm", "ttmscore"),
                    ("lddt", "lddt"),
                    ("qcov", "qcov"),
                    ("tcov", "tcov"),
                    ("evalue", "evalue"),
                ):
                    row[f"{prefix}_best_{key}"] = best[
                        source
                    ]
            else:
                for key in (
                    "target", "prob", "qtm", "ttm",
                    "lddt", "qcov", "tcov", "evalue"
                ):
                    row[f"{prefix}_best_{key}"] = ""

            if n_strong:
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

        output.append(row)

    args.out.parent.mkdir(
        parents=True, exist_ok=True
    )

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=HEADER,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output)

    print(
        "WROTE", args.out,
        "ROWS=", len(output),
        "COLS=", len(HEADER)
    )

    for prefix, _, _ in SOURCES:
        print(
            prefix.upper(),
            "RAW=", totals[(prefix, "raw")],
            "CONVINCING=", totals[(prefix, "convincing")]
        )

    print(
        "EXTERNAL_STATUS_COUNTS=",
        dict(Counter(
            r["external_status"] for r in output
        ))
    )


if __name__ == "__main__":
    main()
