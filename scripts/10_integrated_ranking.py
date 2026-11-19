#!/usr/bin/env python3
"""MetaViro pilot Stage I2 confidence-aware ranking.

Inputs:
    Stage I integrated evidence (36 columns)
    Stage H3 all-member context (19 columns)

The canonical Stage I2 table is not used for calculations.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


DARK = {
    "A_DEEP_DARK": 4.0,
    "B_FUNCTION_UNRESOLVED": 2.5,
    "C_KNOWN_DARK_BRIDGE": 0.5,
    "D_KNOWN_EVIDENCE": 0.0,
}

ASSAY_BASE = {
    "nuclease": 3.0,
    "restriction": 2.5,
    "ssdna_binding": 2.0,
    "recombination": 2.0,
    "ligase": 1.75,
    "helicase": 1.5,
    "primase": 1.25,
    "polymerase": 1.0,
    "dna_methylation": 0.75,
}

I2_FIELDS = [
    "I2_dark",
    "I2_external",
    "I2_structure",
    "I2_context",
    "I2_recurrence",
    "I2_structural_penalty",
    "I2_integrated_score",
    "I2_assay_category",
    "I2_assay_independent_members",
    "I2_assay_independent_sample_span",
    "I2_assay_bonus",
    "I2_expression_bonus",
    "I2_experimental_score",
]


def read_table(path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def fmt(value):
    return f"{value:.2f}"


def structure_score(plddt, frac):
    if plddt >= 70:
        return 3.0
    if plddt >= 60:
        return 2.0
    if plddt >= 50:
        return 1.0
    if frac >= 0.10:
        return 0.5
    return 0.0


def external_score(status, plddt, frac):
    if status == "CONVINCING_STRUCTURE_CONNECTION":
        return 0.0

    if status == "WEAK_RAW_HIT_ONLY":
        return 1.0 if plddt >= 50 else 0.0

    if status == "NO_RAW_HIT_E1":
        if plddt >= 70:
            return 4.0
        if plddt >= 60:
            return 2.5
        if plddt >= 50:
            return 1.5
        if frac >= 0.10:
            return 0.5
        return 0.0

    raise ValueError(f"Unknown external status: {status}")


def context_score(span, count):
    if span <= 0:
        return 0.0
    if span == 1:
        return 0.75
    if span == 2:
        return 2.0
    if span == 3:
        return 3.5 if count >= 5 else 3.0
    if span == 4:
        return 4.5
    if span >= 5:
        return 5.0 if count >= 20 else 4.5
    raise ValueError(span)


def recurrence_score(span):
    return {
        3: 0.5,
        4: 1.0,
        5: 1.5,
    }.get(span, 0.0)


def structural_penalty(dense):
    if dense >= 0.30:
        return -2.0
    if dense >= 0.20:
        return -1.5
    if dense >= 0.10:
        return -0.75
    return 0.0


def expression_bonus(length):
    if length <= 350:
        return 1.0
    if length <= 600:
        return 0.5
    return 0.0


def assay_factor(span, members):
    if span >= 3:
        return 1.0
    if span == 2:
        return 0.75
    if span == 1:
        return min(0.20 * members, 0.40)
    return 0.0


def build_assay_evidence(member_rows):
    """Count noncandidate members and independent samples.

    Commas separate anchor combinations.
    Semicolons within one combination are preserved.
    """

    counts = defaultdict(lambda: defaultdict(set))
    samples = defaultdict(lambda: defaultdict(set))

    for row in member_rows:
        if row["candidate_original"] != "0":
            continue

        if row["anchor_within5"] != "1":
            continue

        rep = row["family_rep"]
        member = row["member"]
        sample = row["sample"]

        categories = {
            category
            for category in row["anchor_categories"].split(",")
            if category
        }

        for category in categories:
            if category not in ASSAY_BASE:
                continue

            counts[rep][category].add(member)
            samples[rep][category].add(sample)

    return counts, samples


def best_assay(rep, counts, samples):
    """Select the category with the highest assay bonus."""

    best = ("", 0, 0, 0.0)
    best_score = 0.0

    # ASSAY_BASE insertion order resolves exact bonus ties.
    for category, base_score in ASSAY_BASE.items():
        members = len(counts[rep][category])
        span = len(samples[rep][category])

        if members == 0:
            continue

        bonus = round(
            base_score * assay_factor(span, members),
            2
        )

        if bonus > best_score:
            best_score = bonus
            best = (
                category,
                members,
                span,
                bonus,
            )

    return best


def calculate(source, counts, samples):
    rep = source["family_rep"]

    plddt = float(source["mean_plddt"])
    frac = float(source["frac_ge70"])

    family_span = int(source["family_sample_span"])

    context_span = int(
        source["noncandidate_context_sample_span"]
    )
    context_count = int(
        source["noncandidate_context_count"]
    )

    dense = float(source["structural_dense_fraction"])
    length = int(source["length_aa"])

    priority = source["priority_class"]

    if priority not in DARK:
        raise ValueError(
            f"Unknown priority class: {priority}"
        )

    dark = DARK[priority]

    ext = external_score(
        source["external_status"],
        plddt,
        frac
    )

    structure = structure_score(plddt, frac)

    context = context_score(
        context_span,
        context_count
    )

    recurrence = recurrence_score(family_span)
    penalty = structural_penalty(dense)

    integrated = round(
        dark
        + ext
        + structure
        + context
        + recurrence
        + penalty,
        2
    )

    category, members, span, assay = best_assay(
        rep, counts, samples
    )

    expression = expression_bonus(length)

    experimental = round(
        integrated + assay + expression,
        2
    )

    output = dict(source)

    scores = {
        "I2_dark": dark,
        "I2_external": ext,
        "I2_structure": structure,
        "I2_context": context,
        "I2_recurrence": recurrence,
        "I2_structural_penalty": penalty,
        "I2_integrated_score": integrated,
        "I2_assay_bonus": assay,
        "I2_expression_bonus": expression,
        "I2_experimental_score": experimental,
    }

    for field, value in scores.items():
        output[field] = fmt(value)

    output["I2_assay_category"] = category

    output["I2_assay_independent_members"] = str(
        members
    )
    output["I2_assay_independent_sample_span"] = str(
        span
    )

    return output



# STAGE10_HISTORICAL_TIE_REPLAY_V1
def load_historical_tie_order(path, rows):
    """Replay archived tie order without reading canonical I2."""

    from collections import defaultdict

    columns, archived = read_table(path)

    required = [
        "structure_id",
        "I2_integrated_score",
        "noncandidate_context_sample_span",
        "historical_tie_order",
    ]

    if columns != required:
        raise ValueError(
            "Historical tie manifest schema mismatch"
        )

    groups = defaultdict(list)
    current = {}

    for row in rows:
        sid = row["structure_id"]

        key = (
            row["I2_integrated_score"],
            int(row["noncandidate_context_sample_span"]),
        )

        if sid in current:
            raise ValueError("Duplicate structure ID")

        current[sid] = key
        groups[key].append(sid)

    expected_tied = {
        sid
        for ids in groups.values()
        if len(ids) > 1
        for sid in ids
    }

    observed = {}
    positions = defaultdict(list)

    for item in archived:
        sid = item["structure_id"]

        if sid in observed:
            raise ValueError(
                "Duplicate historical tie record: " + sid
            )

        if sid not in expected_tied:
            raise ValueError(
                "Unexpected historical tie structure: " + sid
            )

        expected_score, expected_span = current[sid]

        if item["I2_integrated_score"] != expected_score:
            raise ValueError(
                "Historical tie score mismatch: " + sid
            )

        if int(
            item["noncandidate_context_sample_span"]
        ) != expected_span:
            raise ValueError(
                "Historical tie context-span mismatch: " + sid
            )

        rank = int(item["historical_tie_order"])

        if rank < 1:
            raise ValueError("Invalid historical tie rank")

        observed[sid] = rank
        positions[current[sid]].append(rank)

    if set(observed) != expected_tied:
        missing = expected_tied - set(observed)
        extra = set(observed) - expected_tied

        raise ValueError(
            "Incomplete historical tie manifest: "
            f"missing={len(missing)}, extra={len(extra)}"
        )

    for key, ranks in positions.items():
        if sorted(ranks) != list(range(1, len(ranks) + 1)):
            raise ValueError(
                f"Nonconsecutive historical ranks: {key}"
            )

    print(
        "HISTORICAL_TIE_GROUPS =",
        len(positions)
    )
    print(
        "HISTORICAL_TIE_MEMBERS =",
        len(observed)
    )

    return observed


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--evidence",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--member-context",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--out",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--historical-tie-order",
        type=Path,
        default=None,
        help=(
            "Optional archived order for exact integrated-"
            "rank ties; not used for score calculation"
        ),
    )

    args = parser.parse_args()

    if not args.evidence.is_file():
        raise FileNotFoundError(args.evidence)

    if not args.member_context.is_file():
        raise FileNotFoundError(args.member_context)

    input_columns, evidence = read_table(
        args.evidence
    )

    _, member_rows = read_table(
        args.member_context
    )

    if len(evidence) != 191:
        raise ValueError(
            f"Expected 191 families: {len(evidence)}"
        )

    if len(member_rows) != 10327:
        raise ValueError(
            f"Expected 10327 members: {len(member_rows)}"
        )

    if len(input_columns) != 36:
        raise ValueError(
            f"Expected 36 input columns: {len(input_columns)}"
        )

    ids = [r["structure_id"] for r in evidence]

    if len(set(ids)) != 191:
        raise ValueError("Duplicate structure IDs")

    reps = [r["family_rep"] for r in evidence]

    if len(set(reps)) != 191:
        raise ValueError("Duplicate family representatives")

    context_reps = {
        r["family_rep"] for r in member_rows
    }

    if set(reps) != context_reps:
        raise ValueError(
            "Family representatives do not match H3"
        )

    counts, samples = build_assay_evidence(
        member_rows
    )

    rows = [
        calculate(r, counts, samples)
        for r in evidence
    ]

    # Reconstructed numeric ranking rules.
    # The optional manifest replays only historical exact ties.
    historical_order = (
        load_historical_tie_order(
            args.historical_tie_order,
            rows
        )
        if args.historical_tie_order is not None
        else None
    )

    integrated_order = sorted(
        rows,
        key=lambda r: (
            -float(r["I2_integrated_score"]),
            -int(
                r["noncandidate_context_sample_span"]
            ),
            (
                historical_order.get(
                    r["structure_id"], 0
                )
                if historical_order is not None
                else r["structure_id"]
            ),
        )
    )

    for rank, row in enumerate(
        integrated_order, 1
    ):
        row["I2_integrated_rank"] = str(rank)

    # Experimental ranking:
    # score DESC, integrated score DESC,
    # assay sample span DESC, structure ID ASC.
    experimental_order = sorted(
        rows,
        key=lambda r: (
            -float(r["I2_experimental_score"]),
            -float(r["I2_integrated_score"]),
            -int(
                r["I2_assay_independent_sample_span"]
            ),
            r["structure_id"],
        )
    )

    for rank, row in enumerate(
        experimental_order, 1
    ):
        row["I2_experimental_rank"] = str(rank)

    header = (
        [
            "I2_integrated_rank",
            "I2_experimental_rank",
        ]
        + input_columns
        + I2_FIELDS
    )

    if len(header) != 51:
        raise ValueError(
            f"Expected 51 columns: {len(header)}"
        )

    out_resolved = args.out.resolve()

    for protected in (
        args.evidence.resolve(),
        args.member_context.resolve(),
    ):
        if out_resolved == protected:
            raise ValueError(
                "Refusing to overwrite an input"
            )

    if (
        args.out.name
        == "recurrent191_confidence_aware_ranking.tsv"
        and args.out.parent.name
        == "StageI2_ConfidenceAwareRanking"
    ):
        raise ValueError(
            "Refusing to overwrite canonical Stage I2"
        )

    args.out.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(integrated_order)

    print("WROTE =", args.out)
    print("ROWS =", len(rows))
    print("COLUMNS =", len(header))
    print("ASSAY_FAMILIES =", sum(
        bool(r["I2_assay_category"])
        for r in rows
    ))

    print(
        "INTEGRATED_TOP10 =",
        ",".join(
            r["structure_id"]
            for r in integrated_order[:10]
        )
    )

    print(
        "EXPERIMENTAL_TOP10 =",
        ",".join(
            r["structure_id"]
            for r in experimental_order[:10]
        )
    )

    print("STAGE10_RECONSTRUCTION_GENERATED=YES")


if __name__ == "__main__":
    main()
