#!/usr/bin/env python3
"""Portable family-wide genomic neighborhood analysis."""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

STRUCTURAL = re.compile(
    r"\b(?:virion|tail|capsid|portal|terminase|"
    r"baseplate|sheath|fiber|tape measure)\b",
    re.I,
)


def read_tsv(path, required):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path}: missing columns {sorted(missing)}"
            )
        return list(reader)


def read_faa(path):
    result = {}
    pid = None
    seq = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()

            if line.startswith(">"):
                if pid is not None:
                    result[pid] = len("".join(seq).rstrip("*"))

                pid = line[1:].split()[0]
                if pid in result:
                    raise ValueError(f"Duplicate FASTA ID: {pid}")
                seq = []

            elif line:
                if pid is None:
                    raise ValueError("Sequence before FASTA header")
                seq.append(line)

        if pid is not None:
            if pid in result:
                raise ValueError(f"Duplicate FASTA ID: {pid}")
            result[pid] = len("".join(seq).rstrip("*"))

    return result


def structural_ids(path):
    ids = set()
    header = None

    with path.open() as fh:
        for line in fh:
            if line.startswith("#query"):
                header = line[1:].rstrip("\n").split("\t")

                for col in (
                    "query", "Description", "Preferred_name"
                ):
                    if col not in header:
                        raise ValueError(
                            f"Missing eggNOG column: {col}"
                        )
                continue

            if not line.strip() or line.startswith("#"):
                continue

            if header is None:
                raise ValueError("Missing eggNOG header")

            values = line.rstrip("\n").split("\t")
            if len(values) != len(header):
                raise ValueError("Invalid eggNOG column count")

            row = dict(zip(header, values))

            description = " ".join((
                row["Description"],
                row["Preferred_name"],
            ))

            if STRUCTURAL.search(description):
                ids.add(row["query"])

    return ids


def read_gff(path):
    contigs = defaultdict(list)
    ids = set()

    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "CDS":
                continue

            match = re.search(
                r"(?:^|;)ID=([^;]+)", cols[8]
            )

            if not match:
                raise ValueError("GFF CDS missing ID")

            raw = match.group(1)
            index = raw.rsplit("_", 1)[-1]

            if not index.isdigit():
                raise ValueError(f"Invalid Prodigal ID: {raw}")

            # Same native ID convention as Stage 02.
            pid = f"{cols[0]}_{index}"

            if pid in ids:
                raise ValueError(f"Duplicate GFF ID: {pid}")
            ids.add(pid)

            start, end = int(cols[3]), int(cols[4])

            if start < 1 or end < start:
                raise ValueError(f"Invalid GFF interval: {pid}")

            contigs[cols[0]].append(
                (start, end, int(index), pid)
            )

    lookup = {}

    for contig, genes in contigs.items():
        genes.sort(key=lambda x: (x[0], x[1], x[2]))

        for i, gene in enumerate(genes):
            lookup[gene[3]] = (contig, i, genes)

    return lookup


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def fraction(value, total):
    return f"{value / total:.4f}" if total else "0.0000"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--stage01-dir", type=Path, required=True)
    ap.add_argument("--stage02-dir", type=Path, required=True)
    ap.add_argument("--members", type=Path, required=True)
    ap.add_argument("--family-audit", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--min-aa", type=int, default=150)
    args = ap.parse_args()

    if args.window < 0 or args.min_aa < 1:
        raise ValueError("Invalid window or min-aa")

    members = read_tsv(
        args.members, ("member", "representative")
    )
    audits = read_tsv(
        args.family_audit, ("representative", "family_size")
    )

    expected = {
        r["representative"]: int(r["family_size"])
        for r in audits
    }

    if not members or len(expected) != len(audits):
        raise ValueError("Empty or duplicate family audit")

    sample_data = {}
    ownership = {}

    samples = sorted(
        p for p in args.input_dir.iterdir() if p.is_dir()
    )

    for sample_dir in samples:
        sample = sample_dir.name

        lengths = read_faa(
            sample_dir / "proteins.faa"
        )

        for pid in lengths:
            if pid in ownership:
                raise ValueError(
                    f"Protein ID used in multiple samples: {pid}"
                )
            ownership[pid] = sample

        gff = read_gff(sample_dir / "genes.gff")

        unresolved_rows = read_tsv(
            args.stage01_dir / f"{sample}.tsv",
            ("protein_id", "strict_unresolved", "eggnog_status"),
        )

        unresolved = {
            r["protein_id"]: r for r in unresolved_rows
        }

        candidates = {
            r["candidate_id"]
            for r in read_tsv(
                args.stage02_dir / f"{sample}.tsv",
                ("candidate_id",),
            )
        }

        anchor_rows = read_tsv(
            args.stage02_dir / f"{sample}_anchors.tsv",
            ("anchor_id", "anchor_categories"),
        )

        anchors = {
            r["anchor_id"]: set(
                x for x in r["anchor_categories"].split(";") if x
            )
            for r in anchor_rows
        }

        if len(anchors) != len(anchor_rows):
            raise ValueError(f"Duplicate anchors: {sample}")

        structural = structural_ids(
            sample_dir / "annotations.emapper.annotations"
        )

        if set(lengths) != set(gff):
            raise ValueError(f"FASTA/GFF mismatch: {sample}")

        if set(lengths) != set(unresolved):
            raise ValueError(
                f"Stage01/FASTA mismatch: {sample}"
            )

        if not candidates <= set(lengths):
            raise ValueError(
                f"Stage02 candidate mismatch: {sample}"
            )

        if not set(anchors) <= set(gff):
            raise ValueError(
                f"Stage02 anchor/GFF mismatch: {sample}"
            )

        sample_data[sample] = (
            lengths, gff, unresolved,
            candidates, anchors, structural
        )

    context_rows = []
    groups = defaultdict(list)
    seen = set()

    for item in members:
        pid = item["member"]
        rep = item["representative"]

        if pid in seen:
            raise ValueError(f"Duplicate member: {pid}")
        seen.add(pid)

        if rep not in expected:
            raise ValueError(f"Unknown family: {rep}")

        if pid not in ownership:
            raise ValueError(f"Member absent from FASTA: {pid}")

        sample = ownership[pid]

        (
            lengths, gff, unresolved,
            candidates, anchors, structural
        ) = sample_data[sample]

        contig, position, genes = gff[pid]

        left = max(0, position - args.window)
        right = min(len(genes), position + args.window + 1)

        neighbors = [
            gene[3]
            for k, gene in enumerate(genes[left:right], left)
            if k != position
        ]

        nearby_anchors = [
            p for p in neighbors if p in anchors
        ]

        categories = sorted({
            category
            for p in nearby_anchors
            for category in anchors[p]
        })

        anchor_positions = [
            k for k, gene in enumerate(genes)
            if gene[3] in nearby_anchors
        ]

        distance = min(
            (abs(position - k) for k in anchor_positions),
            default=None,
        )

        structural_count = sum(
            p in structural for p in neighbors
        )

        status = unresolved[pid]
        strict = status["strict_unresolved"] == "1"
        length = lengths[pid]

        reconstructed = (
            strict and length >= args.min_aa
            and bool(nearby_anchors)
        )

        original_candidate = pid in candidates

        reasons = []
        if length < args.min_aa:
            reasons.append("BELOW_MIN_AA")
        if not strict:
            reasons.append("NOT_STRICT_UNRESOLVED")
        if not nearby_anchors:
            reasons.append("NO_ANCHOR_IN_WINDOW")

        row = {
            "family_rep": rep,
            "member": pid,
            "sample": sample,
            "length_aa": length,
            "candidate_original": int(original_candidate),
            "eggnog_status": status["eggnog_status"],
            "strict_unresolved_eggnog": int(strict),
            "gff_seqid": contig,
            "gff_local_index": genes[position][2],
            "gff_order": position + 1,
            "window": args.window,
            "neighbors_available": len(neighbors),
            "anchor_within_window": int(bool(nearby_anchors)),
            "anchor_min_gene_distance":
                distance if distance is not None else "",
            "anchor_categories": ";".join(categories),
            "anchor_ids": ";".join(sorted(nearby_anchors)),
            "structural_neighbors_within_window":
                structural_count,
            "structural_dense_ge2":
                int(structural_count >= 2),
            "reconstructed_candidate": int(reconstructed),
            "noncandidate_exclusion_reason":
                ";".join(reasons) if not reconstructed else "",
        }

        context_rows.append(row)
        groups[rep].append(row)

    if set(groups) != set(expected):
        raise ValueError("Family sets differ from Stage06")

    family_rows = []

    for rep, group in sorted(groups.items()):
        n = len(group)

        if n != expected[rep]:
            raise ValueError(
                f"Family size mismatch: {rep}"
            )

        cand = [r for r in group if r["candidate_original"]]
        noncand = [
            r for r in group if not r["candidate_original"]
        ]
        anchor = [
            r for r in group if r["anchor_within_window"]
        ]
        c_anchor = [
            r for r in cand if r["anchor_within_window"]
        ]
        nc_anchor = [
            r for r in noncand if r["anchor_within_window"]
        ]

        counts = Counter()

        for r in anchor:
            for cat in set(
                filter(None, r["anchor_categories"].split(";"))
            ):
                counts[cat] += 1

        dominant, support = (
            sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0]
            if counts else ("", 0)
        )

        nonspan = len({
            r["sample"] for r in nc_anchor
        })

        if nonspan == 0:
            tier = "NO_NONCANDIDATE_CONTEXT_SUPPORT"
        elif nonspan == 1:
            tier = "NONCANDIDATE_CONTEXT_1_SAMPLE"
        elif nonspan == 2:
            tier = "NONCANDIDATE_CONTEXT_2_SAMPLES"
        else:
            tier = "NONCANDIDATE_CONTEXT_GE3_SAMPLES"

        structural_any = sum(
            int(r["structural_neighbors_within_window"] > 0)
            for r in group
        )

        dense = sum(
            r["structural_dense_ge2"] for r in group
        )

        family_rows.append({
            "family_rep": rep,
            "family_size": n,
            "sample_span": len({
                r["sample"] for r in group
            }),
            "samples": ";".join(sorted({
                r["sample"] for r in group
            })),
            "candidate_count": len(cand),
            "noncandidate_count": len(noncand),
            "all_anchor_context_count": len(anchor),
            "all_anchor_context_fraction":
                fraction(len(anchor), n),
            "candidate_anchor_context_count":
                len(c_anchor),
            "candidate_anchor_context_fraction":
                fraction(len(c_anchor), len(cand)),
            "noncandidate_anchor_context_count":
                len(nc_anchor),
            "noncandidate_anchor_context_fraction":
                fraction(len(nc_anchor), len(noncand)),
            "all_context_sample_span": len({
                r["sample"] for r in anchor
            }),
            "noncandidate_context_sample_span": nonspan,
            "dominant_anchor_category": dominant,
            "dominant_anchor_supported_members": support,
            "dominant_anchor_consistency":
                fraction(support, len(anchor)),
            "all_anchor_categories":
                ";".join(sorted(counts)),
            "structural_context_any_count": structural_any,
            "structural_dense_ge2_count": dense,
            "family_context_tier": tier,
        })

    write_tsv(
        args.out_dir / "member_context.tsv",
        context_rows,
        list(context_rows[0]),
    )
    write_tsv(
        args.out_dir / "family_context.tsv",
        family_rows,
        list(family_rows[0]),
    )

    print("STAGE07_FAMILY_CONTEXT=PASS")
    print(f"TOTAL_MEMBERS={len(context_rows)}")
    print(f"TOTAL_FAMILIES={len(family_rows)}")
    print(
        "MEMBERS_WITH_ANCHOR="
        + str(sum(r["anchor_within_window"] for r in context_rows))
    )
    print(f"WINDOW={args.window}")


if __name__ == "__main__":
    main()
