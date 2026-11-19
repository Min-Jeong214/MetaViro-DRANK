#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


MISSING = {
    "",
    "-",
    "--",
    "NA",
    "N/A",
    "na",
    "n/a",
    "None",
    "none",
    "nan",
    "NaN",
}


def blank(x):
    return str(x or "").strip() in MISSING


def read_fasta_lengths(path):
    lengths = {}
    current = None
    seq = []

    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith(">"):
                if current is not None:
                    aa = "".join(seq).rstrip("*")
                    lengths[current] = len(aa)

                current = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

        if current is not None:
            aa = "".join(seq).rstrip("*")
            lengths[current] = len(aa)

    return lengths


def parse_gff(path):
    """
    Reproduce the original Prodigal GFF -> FAA mapping.

    Example:
      contig = NODE_xxx||full
      GFF ID = 3_2
      protein ID = NODE_xxx||full_2

    The last numeric component of the Prodigal ID is the
    local ORF index used in the FAA identifier.
    """
    location = {}
    orf_lookup = {}
    missing_id = 0

    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            p = line.rstrip("\n").split("\t")

            if len(p) < 9 or p[2] != "CDS":
                continue

            contig = p[0]
            attrs = p[8]

            m = re.search(
                r"(?:^|;)ID=([^;]+)",
                attrs
            )

            if not m:
                missing_id += 1
                continue

            gid = m.group(1)

            try:
                local_idx = int(
                    gid.rsplit("_", 1)[1]
                )
            except Exception:
                missing_id += 1
                continue

            qid = f"{contig}_{local_idx}"

            location[qid] = (
                contig,
                local_idx,
            )

            orf_lookup[
                (contig, local_idx)
            ] = qid

    return location, orf_lookup


# ============================================================
# Final Pilot-v2 anchor definitions
# ============================================================

PATTERNS = {
    "polymerase":
        re.compile(r"\bpolymerase\b", re.I),

    "primase":
        re.compile(r"\bprimase\b", re.I),

    "helicase":
        re.compile(r"\bhelicase\b", re.I),

    "nuclease":
        re.compile(
            r"nuclease|"
            r"\brnase\b",
            re.I,
        ),

    "ligase":
        re.compile(r"\bligase\b", re.I),

    "recombination":
        re.compile(
            r"\brecombinase\b|"
            r"\brecombination\b|"
            r"\bresolvase\b|"
            r"\btopoisomerase\b",
            re.I,
        ),

    "dna_methylation":
        re.compile(
            r"\bdna\b.*\bmethylase\b|"
            r"\bdna\b.*\bmethyltransferase\b",
            re.I,
        ),

    "ssdna_binding":
        re.compile(
            r"single[- ]stranded(?:[- ]dna)?[- ]binding|"
            r"\bssb\b|"
            r"\bgp32\b",
            re.I,
        ),

    "restriction":
        re.compile(
            r"restriction enzyme|"
            r"restriction.*subunit|"
            r"\bres subunit\b",
            re.I,
        ),

    "dna_repair":
        re.compile(
            r"dna[- ]repair|"
            r"\bphotolyase\b|"
            r"\bgp59\b|"
            r"\buvsy\b",
            re.I,
        ),

    "reverse_transcriptase":
        re.compile(
            r"reverse transcriptase",
            re.I,
        ),
}


ATPASE_PATTERN = re.compile(r"\batpase\b", re.I)

DNA_CONTEXT_PATTERN = re.compile(
    r"\bdna\b|"
    r"replication|"
    r"repair|"
    r"recombination|"
    r"helicase|"
    r"primase|"
    r"polymerase",
    re.I,
)


def classify_anchor(text):
    text = text or ""

    hits = []

    for category, rx in PATTERNS.items():
        if rx.search(text):
            hits.append(category)

    # DNA-associated ATPase only.
    if (
        ATPASE_PATTERN.search(text)
        and DNA_CONTEXT_PATTERN.search(text)
    ):
        hits.append("dna_associated_atpase")

    return sorted(set(hits))


def load_unresolved(path):
    unresolved = set()

    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")

        for row in r:
            if row.get("ge_min_aa") == "1":
                unresolved.add(row["protein_id"])

    return unresolved


def load_eggnog(path):
    records = {}

    required = [
        "query",
        "eggNOG_OGs",
        "Description",
        "Preferred_name",
    ]

    with open(path, errors="replace") as fh:
        header = None
        idx = None

        for line in fh:
            if line.startswith("#query"):
                header = (
                    line[1:]
                    .rstrip("\n")
                    .split("\t")
                )

                missing = [
                    x for x in required
                    if x not in header
                ]

                if missing:
                    raise SystemExit(
                        f"[ERROR] Missing eggNOG columns: {missing}"
                    )

                idx = {
                    x: header.index(x)
                    for x in required
                }

                continue

            if (
                line.startswith("#")
                or not line.strip()
                or idx is None
            ):
                continue

            p = line.rstrip("\n").split("\t")

            if len(p) <= max(idx.values()):
                continue

            q = p[idx["query"]]

            records[q] = {
                "ogs": p[idx["eggNOG_OGs"]],
                "desc": p[idx["Description"]],
                "pref": p[idx["Preferred_name"]],
            }

    return records


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--sample", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--annotations", required=True)
    ap.add_argument("--unresolved-tsv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--anchors-out", required=True)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--min-aa", type=int, default=150)

    args = ap.parse_args()

    prefix = f"{args.sample}_Virus_Ch"

    lengths_raw = read_fasta_lengths(args.faa)
    unresolved_raw = load_unresolved(args.unresolved_tsv)
    eggnog = load_eggnog(args.annotations)

    location_raw, orf_lookup_raw = parse_gff(args.gff)

    # GFF and eggNOG use raw protein IDs.
    # Validation unresolved table also currently uses raw IDs.
    anchors = []

    for q, rec in eggnog.items():

        # Final Pilot-v2 viral support
        if "10239|Viruses" not in rec["ogs"]:
            continue

        text = " ".join([
            rec["desc"] if not blank(rec["desc"]) else "",
            rec["pref"] if not blank(rec["pref"]) else "",
        ]).strip()

        categories = classify_anchor(text)

        if not categories:
            continue

        if q not in location_raw:
            continue

        contig, idx = location_raw[q]

        anchors.append({
            "anchor_id": q,
            "sample": args.sample,
            "contig": contig,
            "orf_index": idx,
            "anchor_categories": ";".join(categories),
            "anchor_text": text,
        })

    candidates = {}

    for anchor in anchors:
        contig = anchor["contig"]
        idx = anchor["orf_index"]

        for j in range(
            idx - args.window,
            idx + args.window + 1,
        ):
            pid = orf_lookup_raw.get((contig, j))

            if pid is None:
                continue

            if pid not in unresolved_raw:
                continue

            # unresolved file already applies >=150 aa,
            # but retain explicit canonical check.
            if lengths_raw.get(pid, 0) < args.min_aa:
                continue

            key = pid

            if key not in candidates:
                candidates[key] = {
                    "candidate_id": pid,
                    "sample": args.sample,
                    "contig": contig,
                    "anchor_ids": set(),
                    "anchor_categories": set(),
                    "offsets": set(),
                }

            candidates[key]["anchor_ids"].add(
                anchor["anchor_id"]
            )

            candidates[key]["anchor_categories"].update(
                anchor["anchor_categories"].split(";")
            )

            candidates[key]["offsets"].add(
                j - idx
            )

    with open(args.anchors_out, "w", newline="") as out:
        fields = [
            "anchor_id",
            "sample",
            "contig",
            "orf_index",
            "anchor_categories",
            "anchor_text",
        ]

        w = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        w.writeheader()
        w.writerows(anchors)

    with open(args.out, "w", newline="") as out:
        fields = [
            "candidate_id",
            "sample",
            "contig",
            "anchor_count",
            "anchor_ids",
            "anchor_categories",
            "orf_offsets",
        ]

        w = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )

        w.writeheader()

        for pid in sorted(candidates):
            r = candidates[pid]

            w.writerow({
                "candidate_id": pid,
                "sample": r["sample"],
                "contig": r["contig"],
                "anchor_count": len(r["anchor_ids"]),
                "anchor_ids":
                    ";".join(sorted(r["anchor_ids"])),
                "anchor_categories":
                    ";".join(sorted(r["anchor_categories"])),
                "orf_offsets":
                    ";".join(
                        str(x)
                        for x in sorted(r["offsets"])
                    ),
            })

    print(
        f"SAMPLE={args.sample}"
    )
    print(
        f"VIRAL_ANCHORS={len(anchors)}"
    )
    print(
        f"CANDIDATE_ORFS={len(candidates)}"
    )


if __name__ == "__main__":
    main()
