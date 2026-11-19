#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict, Counter


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--member-evidence",
        required=True,
        help=(
            "recurrent191_all_members_sequence_evidence.tsv "
            "or equivalent member-level evidence table"
        ),
    )

    ap.add_argument(
        "--out",
        required=True,
    )

    args = ap.parse_args()

    families = defaultdict(list)

    with open(args.member_evidence) as fh:
        r = csv.DictReader(fh, delimiter="\t")

        required = {
            "member",
            "representative",
            "sequence_status",
            "phrog_structural",
            "phrog_nucleic",
        }

        missing = required - set(r.fieldnames or [])

        if missing:
            raise SystemExit(
                "[ERROR] Missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in r:
            families[row["representative"]].append(row)

    rows = []

    for rep, members in families.items():

        c = Counter(
            x["sequence_status"]
            for x in members
        )

        n = len(members)

        known = c["KNOWN_SEQUENCE_EVIDENCE"]
        strict = c["STRICT_DARK"]
        duf = c["DUF_ONLY"]
        phrog_unknown = c["PHROG_UNKNOWN"]

        unresolved = (
            strict
            + duf
            + phrog_unknown
        )

        if strict == n:
            family_status = (
                "ALL_MEMBER_STRICT_DARK"
            )
            priority_class = "A_DEEP_DARK"

        elif known == 0:
            family_status = (
                "ALL_MEMBER_FUNCTION_UNRESOLVED"
            )
            priority_class = (
                "B_FUNCTION_UNRESOLVED"
            )

        elif known == n:
            family_status = (
                "KNOWN_EVIDENCE_PRESENT"
            )
            priority_class = (
                "D_KNOWN_EVIDENCE"
            )

        else:
            family_status = (
                "MIXED_KNOWN_UNRESOLVED"
            )
            priority_class = (
                "C_KNOWN_DARK_BRIDGE"
            )

        rows.append({
            "representative": rep,
            "family_size": n,
            "known_members": known,
            "strict_dark_members": strict,
            "duf_only_members": duf,
            "phrog_unknown_members":
                phrog_unknown,
            "function_unresolved_members":
                unresolved,
            "all_members_function_unresolved":
                int(known == 0),
            "all_members_strict_dark":
                int(strict == n),
            "family_sequence_status":
                family_status,
            "priority_class":
                priority_class,
            "members_with_phrog_structural":
                sum(
                    int(x["phrog_structural"])
                    for x in members
                ),
            "members_with_phrog_nucleic":
                sum(
                    int(x["phrog_nucleic"])
                    for x in members
                ),
        })

    rows.sort(
        key=lambda x: x["representative"]
    )

    fields = [
        "representative",
        "family_size",
        "known_members",
        "strict_dark_members",
        "duf_only_members",
        "phrog_unknown_members",
        "function_unresolved_members",
        "all_members_function_unresolved",
        "all_members_strict_dark",
        "family_sequence_status",
        "priority_class",
        "members_with_phrog_structural",
        "members_with_phrog_nucleic",
    ]

    with open(
        args.out,
        "w",
        newline="",
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        w.writeheader()
        w.writerows(rows)

    counts = Counter(
        x["priority_class"]
        for x in rows
    )

    print(
        f"TOTAL_FAMILIES={len(rows)}"
    )

    for cls in [
        "A_DEEP_DARK",
        "B_FUNCTION_UNRESOLVED",
        "C_KNOWN_DARK_BRIDGE",
        "D_KNOWN_EVIDENCE",
    ]:
        print(
            f"{cls}={counts[cls]}"
        )

    print(
        "ALL_MEMBERS_FUNCTION_UNRESOLVED="
        + str(
            sum(
                x[
                    "all_members_function_unresolved"
                ]
                for x in rows
            )
        )
    )

    print(
        "ALL_MEMBERS_STRICT_DARK="
        + str(
            sum(
                x["all_members_strict_dark"]
                for x in rows
            )
        )
    )


if __name__ == "__main__":
    main()
