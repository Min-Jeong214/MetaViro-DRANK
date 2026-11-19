#!/usr/bin/env python3
"""Conservative family-level darkness audit."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

STATUSES = {
    "KNOWN_SEQUENCE_EVIDENCE",
    "STRICT_DARK",
    "DUF_ONLY",
    "PHROG_UNKNOWN",
    "UNRESOLVED_HOMOLOGY",
}

REQUIRED = {
    "member",
    "representative",
    "search_status",
    "sequence_status",
    "phrog_structural",
    "phrog_nucleic",
}

FIELDS = [
    "representative",
    "family_size",
    "search_complete_members",
    "known_members",
    "strict_dark_members",
    "duf_only_members",
    "phrog_unknown_members",
    "unresolved_homology_members",
    "function_unresolved_members",
    "all_members_function_unresolved",
    "all_members_strict_dark",
    "family_sequence_status",
    "priority_class",
    "members_with_phrog_structural",
    "members_with_phrog_nucleic",
]


def audit(input_path):
    families = defaultdict(list)
    seen_members = set()

    with input_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Missing columns: " + ", ".join(sorted(missing))
            )

        for line_number, row in enumerate(reader, 2):
            member = (row["member"] or "").strip()
            rep = (row["representative"] or "").strip()
            status = (row["sequence_status"] or "").strip()
            search = (row["search_status"] or "").strip()

            if not member or not rep:
                raise ValueError(
                    f"Line {line_number}: empty member/representative"
                )

            if member in seen_members:
                raise ValueError(
                    f"Line {line_number}: duplicate member {member}"
                )
            seen_members.add(member)

            if search != "COMPLETE":
                raise ValueError(
                    f"Line {line_number}: search not complete "
                    f"for {member}: {search}"
                )

            if status not in STATUSES:
                raise ValueError(
                    f"Line {line_number}: invalid sequence_status "
                    f"for {member}: {status}"
                )

            for key in ("phrog_structural", "phrog_nucleic"):
                if row[key] not in ("0", "1"):
                    raise ValueError(
                        f"Line {line_number}: invalid {key} "
                        f"for {member}: {row[key]}"
                    )

            families[rep].append(row)

    if not families:
        raise ValueError("No member evidence records")

    results = []

    for rep, members in sorted(families.items()):
        counts = Counter(
            row["sequence_status"] for row in members
        )

        n = len(members)
        known = counts["KNOWN_SEQUENCE_EVIDENCE"]
        strict = counts["STRICT_DARK"]
        duf = counts["DUF_ONLY"]
        unknown = counts["PHROG_UNKNOWN"]
        homology = counts["UNRESOLVED_HOMOLOGY"]

        unresolved = strict + duf + unknown + homology

        if strict == n:
            family_status = "ALL_MEMBER_STRICT_DARK"
            priority = "A_DEEP_DARK"
        elif known == 0:
            family_status = "ALL_MEMBER_FUNCTION_UNRESOLVED"
            priority = "B_FUNCTION_UNRESOLVED"
        elif known == n:
            family_status = "KNOWN_EVIDENCE_PRESENT"
            priority = "D_KNOWN_EVIDENCE"
        else:
            family_status = "MIXED_KNOWN_UNRESOLVED"
            priority = "C_KNOWN_DARK_BRIDGE"

        results.append({
            "representative": rep,
            "family_size": n,
            "search_complete_members": n,
            "known_members": known,
            "strict_dark_members": strict,
            "duf_only_members": duf,
            "phrog_unknown_members": unknown,
            "unresolved_homology_members": homology,
            "function_unresolved_members": unresolved,
            "all_members_function_unresolved": int(known == 0),
            "all_members_strict_dark": int(strict == n),
            "family_sequence_status": family_status,
            "priority_class": priority,
            "members_with_phrog_structural": sum(
                int(row["phrog_structural"])
                for row in members
            ),
            "members_with_phrog_nucleic": sum(
                int(row["phrog_nucleic"])
                for row in members
            ),
        })

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--member-evidence", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    results = audit(args.member_evidence)

    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)

    counts = Counter(row["priority_class"] for row in results)

    print("FAMILY_DARKNESS_AUDIT=PASS")
    print(f"TOTAL_FAMILIES={len(results)}")
    print(
        "TOTAL_MEMBERS="
        + str(sum(row["family_size"] for row in results))
    )

    for key in (
        "A_DEEP_DARK",
        "B_FUNCTION_UNRESOLVED",
        "C_KNOWN_DARK_BRIDGE",
        "D_KNOWN_EVIDENCE",
    ):
        print(f"{key}={counts[key]}")


if __name__ == "__main__":
    main()
