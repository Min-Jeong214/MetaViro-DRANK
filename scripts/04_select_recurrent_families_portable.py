#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def sample_from_prefixed_id(pid):
    m = re.match(r"(V\d+)_Virus_Ch\|", pid)
    if m:
        return m.group(1)
    return None


def load_candidates(path):
    """
    Accepts candidate tables containing one of:
      candidate
      candidate_id
      protein_id
      id

    Returns:
      candidate_sample[protein_id] = sample
    """
    candidate_sample = {}

    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")

        fields = r.fieldnames or []

        idcol = None

        for x in [
            "candidate",
            "candidate_id",
            "protein_id",
            "id",
        ]:
            if x in fields:
                idcol = x
                break

        if idcol is None:
            raise SystemExit(
                "[ERROR] Candidate ID column not found."
            )

        for row in r:
            pid = (row.get(idcol) or "").strip()

            if not pid:
                continue

            sample = (
                (row.get("sample") or "").strip()
                or sample_from_prefixed_id(pid)
            )

            candidate_sample[pid] = sample

    return candidate_sample


def load_unresolved(paths):
    """
    validation/01 tables contain raw protein IDs.

    Protein IDs are used unchanged across all input files.

    Family-level filter uses strict_unresolved only.
    ge_min_aa is NOT applied to all family members.
    """
    strict = {}

    for item in paths:

        if "=" not in item:
            raise SystemExit(
                "[ERROR] --unresolved must be SAMPLE=FILE"
            )

        sample, path = item.split("=", 1)

        # Preserve original FASTA protein identifiers.

        with open(path) as fh:
            r = csv.DictReader(fh, delimiter="\t")

            for row in r:
                raw = row["protein_id"]
                pid = raw

                strict[pid] = (
                    row["strict_unresolved"] == "1"
                )

    return strict


def load_clusters(path):
    fam_members = defaultdict(list)
    member_rep = {}

    with open(path) as fh:
        for line in fh:

            if not line.strip():
                continue

            rep, member = (
                line.rstrip("\n")
                .split("\t")[:2]
            )

            fam_members[rep].append(member)
            member_rep[member] = rep

    return fam_members, member_rep


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--cluster-tsv",
        required=True,
        help="MMseqs cluster TSV: representative<TAB>member",
    )

    ap.add_argument(
        "--candidates",
        required=True,
        help="Candidate/provenance TSV",
    )

    ap.add_argument(
        "--unresolved",
        action="append",
        required=True,
        help=(
            "Per-sample strict-unresolved table as "
            "SAMPLE=FILE. Repeat for each sample."
        ),
    )

    ap.add_argument(
        "--out-prefix",
        required=True,
    )

    args = ap.parse_args()

    candidate_sample = load_candidates(
        args.candidates
    )

    strict_unresolved = load_unresolved(
        args.unresolved
    )

    fam_members, member_rep = load_clusters(
        args.cluster_tsv
    )

    missing = [
        pid
        for pid in candidate_sample
        if pid not in member_rep
    ]

    if missing:
        raise SystemExit(
            f"[ERROR] {len(missing)} candidates "
            "missing from cluster TSV."
        )

    family_candidates = defaultdict(list)

    for pid in candidate_sample:
        rep = member_rep[pid]
        family_candidates[rep].append(pid)

    candidate_families = set(
        family_candidates
    )

    all_unresolved = {}

    for rep in candidate_families:
        all_unresolved[rep] = all(
            strict_unresolved.get(member, False)
            for member in fam_members[rep]
        )

    eligible = {
        rep
        for rep in candidate_families
        if all_unresolved[rep]
    }

    rows = []

    for rep in sorted(candidate_families):

        candidates = family_candidates[rep]

        samples = sorted({
            candidate_sample[pid]
            for pid in candidates
            if candidate_sample[pid]
        })

        rows.append({
            "representative": rep,
            "family_size": len(
                fam_members[rep]
            ),
            "candidate_count": len(
                candidates
            ),
            "all_members_unresolved":
                int(all_unresolved[rep]),
            "sample_span": len(samples),
            "samples": ";".join(samples),
        })

    out_all = Path(
        args.out_prefix + "_all.tsv"
    )

    with out_all.open(
        "w",
        newline="",
    ) as fh:

        fields = [
            "representative",
            "family_size",
            "candidate_count",
            "all_members_unresolved",
            "sample_span",
            "samples",
        ]

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        w.writeheader()
        w.writerows(rows)

    for threshold in [2, 3, 4, 5]:

        out = Path(
            f"{args.out_prefix}_ge{threshold}.tsv"
        )

        with out.open(
            "w",
            newline="",
        ) as fh:

            fields = [
                "representative",
                "family_size",
                "candidate_count",
                "sample_span",
                "samples",
            ]

            w = csv.DictWriter(
                fh,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )

            w.writeheader()

            for row in rows:

                if not row[
                    "all_members_unresolved"
                ]:
                    continue

                if row["sample_span"] < threshold:
                    continue

                w.writerow({
                    "representative":
                        row["representative"],
                    "family_size":
                        row["family_size"],
                    "candidate_count":
                        row["candidate_count"],
                    "sample_span":
                        row["sample_span"],
                    "samples":
                        row["samples"],
                })

    spans = {}

    for rep in eligible:

        samples = {
            candidate_sample[x]
            for x in family_candidates[rep]
            if candidate_sample[x]
        }

        spans[rep] = len(samples)

    print(
        f"CANDIDATE_ORFS="
        f"{len(candidate_sample)}"
    )

    print(
        f"CANDIDATE_FAMILIES="
        f"{len(candidate_families)}"
    )

    print(
        "ALL_MEMBERS_UNRESOLVED_FAMILIES="
        f"{len(eligible)}"
    )

    for threshold in [2, 3, 4, 5]:

        n = sum(
            span >= threshold
            for span in spans.values()
        )

        print(
            f"RECURRENT_GE{threshold}={n}"
        )


if __name__ == "__main__":
    main()
