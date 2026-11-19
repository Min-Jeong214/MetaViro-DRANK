#!/usr/bin/env python3

###############################################################################
# Stage L
# MSA conservation + MSA-structure-supported residue prioritization
#
# Canonically reconstructed from frozen Stage L outputs.
#
# Residue-level rules
# -------------------
# Alignment statistics are calculated only from the 20 standard amino acids.
# Gaps and non-standard residues such as X are excluded.
#
# msa_coverage
#   = number of standard-AA residues in alignment column / total MSA sequences
#
# conservation_entropy_score
#   = 1 - H / ln(20)
#   where H is Shannon entropy across standard amino acids only.
#
# query_residue_fraction
#   = frequency of the query residue among standard-AA residues in the column
#
# high_conserved
#   = msa_coverage >= 0.70
#     AND query_residue_fraction >= 0.70
#
# very_conserved
#   = msa_coverage >= 0.70
#     AND query_residue_fraction >= 0.85
#
# structure_supported
#   = residue pLDDT >= 70
#
# priority_residue
#   = high_conserved
#     AND structure_supported
#     AND query residue in DEKRHSTCYNQ
#
#
# Conserved structural windows
# ----------------------------
# Window size = 9 aa
#
# conserved_structural
#   = number of residues satisfying:
#       high_conserved == 1
#       AND structure_supported == 1
#
# charged
#   = number of priority residues in D/E/K/R/H
#
# Candidate window:
#   conserved_structural >= 4
#   AND charged >= 2
#
# Candidates are ranked by:
#   conserved_structural DESC
#   charged DESC
#   mean pLDDT DESC
#   mean query fraction DESC
#   start ASC
#
# Greedy redundancy removal:
#   reject a window if it overlaps an already selected window
#   from the same family by >= 5 residues.
###############################################################################

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

ACIDIC = set("DE")
BASIC = set("KRH")

# Residues allowed for priority_residue.
FUNCTIONAL = set("DEKRHSTCYNQ")

# Canonical "nucleophilic" summary category:
# residues whose residue_class contains nucleophilic.
NUCLEOPHILIC = set("STCY")

TOP5 = [
    "F135",
    "F005",
    "F008",
    "F025",
    "F023",
]


###############################################################################
# FASTA
###############################################################################

def read_fasta(path):
    records = []
    header = None
    seq = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        (
                            header,
                            "".join(seq),
                        )
                    )

                header = line[1:].split()[0]
                seq = []

            else:
                seq.append(line.strip())

        if header is not None:
            records.append(
                (
                    header,
                    "".join(seq),
                )
            )

    return records


###############################################################################
# Residue class
###############################################################################

def residue_class(aa):

    if aa in ACIDIC:
        return "acidic"

    if aa in BASIC:
        return "basic"

    if aa in set("ST"):
        return "polar,nucleophilic"

    if aa in set("NQ"):
        return "polar"

    if aa in set("CY"):
        return "nucleophilic"

    if aa in set("GP"):
        return "gly_pro"

    return "other"


###############################################################################
# Shannon conservation score
###############################################################################

def conservation_score(chars):

    counts = Counter(chars)
    total = sum(counts.values())

    if total == 0:
        return 0.0

    entropy = 0.0

    for n in counts.values():
        p = n / total
        entropy -= p * math.log(p)

    return 1.0 - (
        entropy / math.log(20.0)
    )


###############################################################################
# pLDDT from AlphaFold/ColabFold PDB
#
# Canonical residue confidence is taken from CA-atom B factors.
###############################################################################

def read_plddt_from_pdb(path):

    plddt = {}

    with open(path) as fh:

        for line in fh:

            if not line.startswith("ATOM"):
                continue

            atom = line[12:16].strip()

            if atom != "CA":
                continue

            try:
                resseq = int(
                    line[22:26].strip()
                )

                bfactor = float(
                    line[60:66].strip()
                )

            except ValueError:
                continue

            # Single-chain query structure.
            # First CA for each residue number.
            if resseq not in plddt:
                plddt[resseq] = bfactor

    return plddt


###############################################################################
# Locate PDB
###############################################################################

def find_pdb(pdb_dir, sid):

    patterns = [
        f"{sid}_unrelaxed_rank_001_*.pdb",
        f"{sid}_relaxed_rank_001_*.pdb",
        f"{sid}*.pdb",
    ]

    for pattern in patterns:

        hits = sorted(
            pdb_dir.glob(pattern)
        )

        if hits:
            return hits[0]

    raise FileNotFoundError(
        f"No PDB found for {sid} in {pdb_dir}"
    )


###############################################################################
# Build residue table
###############################################################################

def build_residue_table(
    sid,
    aln_path,
    pdb_path,
):

    records = read_fasta(
        aln_path
    )

    if not records:
        raise RuntimeError(
            f"{sid}: empty alignment"
        )

    query = records[0][1]
    nseq = len(records)

    width = len(query)

    if any(
        len(seq) != width
        for _, seq in records
    ):
        raise RuntimeError(
            f"{sid}: inconsistent alignment widths"
        )

    plddt = read_plddt_from_pdb(
        pdb_path
    )

    rows = []
    query_pos = 0

    for col in range(width):

        query_aa = query[col]

        # Query-centered coordinate system.
        if query_aa == "-":
            continue

        query_pos += 1

        # Canonical conservation statistics:
        # standard 20 aa only.
        valid = [
            seq[col]
            for _, seq in records
            if seq[col] in STANDARD_AA
        ]

        if not valid:
            raise RuntimeError(
                f"{sid}: no standard AA at query position {query_pos}"
            )

        counts = Counter(
            valid
        )

        dominant_aa, dominant_n = (
            counts.most_common(1)[0]
        )

        coverage = (
            len(valid) / nseq
        )

        query_fraction = (
            counts.get(query_aa, 0)
            / len(valid)
        )

        dominant_fraction = (
            dominant_n
            / len(valid)
        )

        entropy_score = (
            conservation_score(
                valid
            )
        )

        residue_plddt = plddt.get(
            query_pos
        )

        if residue_plddt is None:
            raise RuntimeError(
                f"{sid}: PDB missing CA/pLDDT for residue {query_pos}"
            )

        high_conserved = int(
            coverage >= 0.70
            and
            query_fraction >= 0.70
        )

        very_conserved = int(
            coverage >= 0.70
            and
            query_fraction >= 0.85
        )

        structure_supported = int(
            residue_plddt >= 70.0
        )

        priority_residue = int(
            high_conserved == 1
            and
            structure_supported == 1
            and
            query_aa in FUNCTIONAL
        )

        rows.append(
            {
                "structure_id": sid,
                "query_position": query_pos,
                "query_aa": query_aa,
                "residue_class":
                    residue_class(query_aa),
                "msa_sequences":
                    nseq,
                "msa_coverage":
                    coverage,
                "conservation_entropy_score":
                    entropy_score,
                "query_residue_fraction":
                    query_fraction,
                "dominant_aa":
                    dominant_aa,
                "dominant_fraction":
                    dominant_fraction,
                "msa_plddt":
                    residue_plddt,
                "high_conserved":
                    high_conserved,
                "very_conserved":
                    very_conserved,
                "structure_supported":
                    structure_supported,
                "priority_residue":
                    priority_residue,
            }
        )

    return rows


###############################################################################
# TSV output
###############################################################################

RESIDUE_COLUMNS = [
    "structure_id",
    "query_position",
    "query_aa",
    "residue_class",
    "msa_sequences",
    "msa_coverage",
    "conservation_entropy_score",
    "query_residue_fraction",
    "dominant_aa",
    "dominant_fraction",
    "msa_plddt",
    "high_conserved",
    "very_conserved",
    "structure_supported",
    "priority_residue",
]


def write_residue_table(
    path,
    rows,
):

    with open(
        path,
        "w",
        newline="",
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            RESIDUE_COLUMNS
        )

        for r in rows:

            writer.writerow(
                [
                    r["structure_id"],
                    r["query_position"],
                    r["query_aa"],
                    r["residue_class"],
                    r["msa_sequences"],
                    f'{r["msa_coverage"]:.4f}',
                    f'{r["conservation_entropy_score"]:.4f}',
                    f'{r["query_residue_fraction"]:.4f}',
                    r["dominant_aa"],
                    f'{r["dominant_fraction"]:.4f}',
                    f'{r["msa_plddt"]:.2f}',
                    r["high_conserved"],
                    r["very_conserved"],
                    r["structure_supported"],
                    r["priority_residue"],
                ]
            )


###############################################################################
# Summary
###############################################################################

def build_summary(
    sid,
    rows,
):

    priority = [
        r
        for r in rows
        if r["priority_residue"] == 1
    ]

    return {
        "structure_id":
            sid,

        "length":
            len(rows),

        "msa_nseq":
            rows[0]["msa_sequences"],

        "very_conserved_positions":
            sum(
                r["very_conserved"]
                for r in rows
            ),

        "priority_structural_positions":
            len(priority),

        "priority_acidic":
            sum(
                r["query_aa"] in ACIDIC
                for r in priority
            ),

        "priority_basic":
            sum(
                r["query_aa"] in BASIC
                for r in priority
            ),

        "priority_nucleophilic":
            sum(
                r["query_aa"] in NUCLEOPHILIC
                for r in priority
            ),
    }


###############################################################################
# Conserved windows
###############################################################################

def window_overlap(a, b):

    return max(
        0,
        min(a["end"], b["end"])
        -
        max(a["start"], b["start"])
        + 1
    )


def make_windows(
    sid,
    rows,
    size=9,
):

    candidates = []

    for i in range(
        len(rows) - size + 1
    ):

        w = rows[
            i:i + size
        ]

        start = w[0][
            "query_position"
        ]

        end = w[-1][
            "query_position"
        ]

        sequence = "".join(
            r["query_aa"]
            for r in w
        )

        conserved_structural = sum(
            (
                r["high_conserved"] == 1
                and
                r["structure_supported"] == 1
            )
            for r in w
        )

        charged = sum(
            (
                r["priority_residue"] == 1
                and
                r["query_aa"] in (
                    ACIDIC | BASIC
                )
            )
            for r in w
        )

        acidic = sum(
            (
                r["priority_residue"] == 1
                and
                r["query_aa"] in ACIDIC
            )
            for r in w
        )

        basic = sum(
            (
                r["priority_residue"] == 1
                and
                r["query_aa"] in BASIC
            )
            for r in w
        )

        mean_qfrac = (
            sum(
                r["query_residue_fraction"]
                for r in w
            )
            / size
        )

        mean_plddt = (
            sum(
                r["msa_plddt"]
                for r in w
            )
            / size
        )

        if (
            conserved_structural >= 4
            and
            charged >= 2
        ):

            candidates.append(
                {
                    "structure_id":
                        sid,

                    "start":
                        start,

                    "end":
                        end,

                    "sequence":
                        sequence,

                    "conserved_structural":
                        conserved_structural,

                    "charged":
                        charged,

                    "acidic":
                        acidic,

                    "basic":
                        basic,

                    "mean_query_fraction":
                        mean_qfrac,

                    "mean_plddt":
                        mean_plddt,
                }
            )

    # Canonical ranking sufficient to reproduce
    # the frozen selected-window set.
    candidates.sort(
        key=lambda r: (
            -r["conserved_structural"],
            -r["charged"],
            -r["mean_plddt"],
            -r["mean_query_fraction"],
            r["start"],
        )
    )

    selected = []

    for cand in candidates:

        redundant = False

        for previous in selected:

            if (
                window_overlap(
                    cand,
                    previous,
                )
                >= 5
            ):
                redundant = True
                break

        if not redundant:
            selected.append(
                cand
            )

    # Canonical output order:
    # strongest window first within each family.
    selected.sort(
        key=lambda r: (
            -r["conserved_structural"],
            -r["charged"],
            -r["mean_plddt"],
            -r["mean_query_fraction"],
            r["start"],
        )
    )

    return selected


###############################################################################
# Main
###############################################################################

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--aln-dir",
        required=True,
        help="StageJ_FamilyMSA/mafft",
    )

    ap.add_argument(
        "--pdb-dir",
        required=True,
        help="StageK_MSARePrediction/predictions",
    )

    ap.add_argument(
        "--out-dir",
        required=True,
    )

    args = ap.parse_args()

    aln_dir = Path(
        args.aln_dir
    )

    pdb_dir = Path(
        args.pdb_dir
    )

    out_dir = Path(
        args.out_dir
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    all_rows = {}
    summaries = []
    windows = []

    for sid in TOP5:

        aln_path = (
            aln_dir /
            f"{sid}.aligned.faa"
        )

        if not aln_path.exists():
            raise FileNotFoundError(
                aln_path
            )

        pdb_path = find_pdb(
            pdb_dir,
            sid,
        )

        rows = build_residue_table(
            sid,
            aln_path,
            pdb_path,
        )

        all_rows[sid] = rows

        write_residue_table(
            out_dir /
            f"{sid}_residue_conservation.tsv",
            rows,
        )

        summaries.append(
            build_summary(
                sid,
                rows,
            )
        )

        windows.extend(
            make_windows(
                sid,
                rows,
            )
        )

        print(
            sid,
            f"length={len(rows)}",
            f"nseq={rows[0]['msa_sequences']}",
            f"very={sum(r['very_conserved'] for r in rows)}",
            f"priority={sum(r['priority_residue'] for r in rows)}",
            sep="\t",
        )


    ###########################################################################
    # summary
    ###########################################################################

    summary_cols = [
        "structure_id",
        "length",
        "msa_nseq",
        "very_conserved_positions",
        "priority_structural_positions",
        "priority_acidic",
        "priority_basic",
        "priority_nucleophilic",
    ]

    with open(
        out_dir /
        "five_candidate_conservation_summary.tsv",
        "w",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=summary_cols,
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(
            summaries
        )


    ###########################################################################
    # conserved windows
    ###########################################################################

    window_cols = [
        "structure_id",
        "start",
        "end",
        "sequence",
        "conserved_structural",
        "charged",
        "acidic",
        "basic",
        "mean_query_fraction",
        "mean_plddt",
    ]

    with open(
        out_dir /
        "conserved_windows.tsv",
        "w",
        newline="",
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            window_cols
        )

        for r in windows:

            writer.writerow(
                [
                    r["structure_id"],
                    r["start"],
                    r["end"],
                    r["sequence"],
                    r["conserved_structural"],
                    r["charged"],
                    r["acidic"],
                    r["basic"],
                    f'{r["mean_query_fraction"]:.3f}',
                    f'{r["mean_plddt"]:.2f}',
                ]
            )


if __name__ == "__main__":
    main()
