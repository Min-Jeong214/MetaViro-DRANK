#!/usr/bin/env python3

###############################################################################
# Stage M
# Final integrated candidate evidence assembly
#
# Canonical reconstruction from frozen Stage M outputs.
#
# Important:
#   This stage does NOT re-rank all 191 families.
#
#   The final selected sets are frozen integrated selections:
#
#       Top5:
#       F135, F005, F008, F025, F023
#
#       Primary Top3:
#       F135, F008, F005
#
#   Stage M integrates evidence from:
#
#     - Stage I2 confidence-aware ranking
#     - Stage K family-MSA structure prediction
#     - Stage K2 external structure search
#     - Stage L conservation analysis
#     - Stage L3 charged spatial clusters
#
#   It writes:
#
#     final5_evidence_matrix.tsv
#     final_top5_representatives.faa
#     primary_top3_representatives.faa
###############################################################################

import argparse
import csv
from collections import defaultdict
from pathlib import Path


TOP5 = [
    "F135",
    "F005",
    "F008",
    "F025",
    "F023",
]

TOP3 = [
    "F135",
    "F008",
    "F005",
]


###############################################################################
# Generic TSV
###############################################################################

def read_tsv(path):

    with open(path) as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t",
            )
        )


def index_by(rows, column="structure_id"):

    return {
        r[column]: r
        for r in rows
        if r.get(column)
    }


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
                            "".join(seq).rstrip("*"),
                        )
                    )

                header = line[1:].split()[0]
                seq = []

            else:

                seq.append(
                    line.strip()
                )

        if header is not None:

            records.append(
                (
                    header,
                    "".join(seq).rstrip("*"),
                )
            )

    return records


def write_fasta(path, ids, sequences):

    with open(path, "w") as fh:

        for sid in ids:

            seq = sequences[sid]

            fh.write(
                f">{sid}\n"
            )

            # Canonical output uses ordinary wrapped FASTA.
            # 60 aa per line.
            for i in range(
                0,
                len(seq),
                60,
            ):
                fh.write(
                    seq[i:i+60]
                    + "\n"
                )


###############################################################################
# Spatial cluster aggregation
###############################################################################

def aggregate_spatial_clusters(rows):

    grouped = defaultdict(list)

    for r in rows:

        sid = r.get(
            "structure_id",
            "",
        )

        if sid:
            grouped[sid].append(r)

    result = {}

    for sid, rr in grouped.items():

        rr = sorted(
            rr,
            key=lambda x:
            int(
                x.get(
                    "cluster_rank",
                    "999999",
                )
            ),
        )

        result[sid] = {
            "count":
                str(len(rr)),

            "clusters":
                ";".join(
                    r.get(
                        "residues",
                        "",
                    )
                    for r in rr
                ),
        }

    return result


###############################################################################
# Main
###############################################################################

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--i2-ranking",
        required=True,
    )

    ap.add_argument(
        "--msa-qc",
        required=True,
    )

    ap.add_argument(
        "--external-matrix",
        required=True,
    )

    ap.add_argument(
        "--conservation",
        required=True,
    )

    ap.add_argument(
        "--spatial-clusters",
        required=True,
    )

    ap.add_argument(
        "--structure-fasta",
        required=True,
        help=(
            "StageF structure_all191.faa "
            "with Fxxx headers"
        ),
    )

    ap.add_argument(
        "--out-dir",
        required=True,
    )

    args = ap.parse_args()

    out = Path(
        args.out_dir
    )

    out.mkdir(
        parents=True,
        exist_ok=True,
    )


    ###########################################################################
    # Load inputs
    ###########################################################################

    i2 = index_by(
        read_tsv(
            args.i2_ranking
        )
    )

    msa = index_by(
        read_tsv(
            args.msa_qc
        )
    )

    ext = index_by(
        read_tsv(
            args.external_matrix
        )
    )

    cons = index_by(
        read_tsv(
            args.conservation
        )
    )

    spatial = aggregate_spatial_clusters(
        read_tsv(
            args.spatial_clusters
        )
    )

    structure_sequences = dict(
        read_fasta(
            args.structure_fasta
        )
    )


    ###########################################################################
    # Validate selected IDs
    ###########################################################################

    for sid in TOP5:

        for name, idx in [
            ("I2", i2),
            ("MSA", msa),
            ("EXTERNAL", ext),
            ("CONSERVATION", cons),
        ]:

            if sid not in idx:

                raise RuntimeError(
                    f"{sid}: missing from {name}"
                )

        if sid not in structure_sequences:

            raise RuntimeError(
                f"{sid}: missing from structure FASTA"
            )


    ###########################################################################
    # Build evidence matrix
    ###########################################################################

    columns = [
        "structure_id",
        "priority_class",
        "experimental_score",
        "dominant_context",
        "assay_category",
        "independent_context_n",
        "independent_context_sample_span",
        "single_mean_plddt",
        "msa_mean_plddt",
        "delta_mean_plddt",
        "single_ptm",
        "msa_ptm",
        "delta_ptm",
        "external_status",
        "raw_sources",
        "convincing_sources",
        "BFVD_best_target",
        "BFVD_best_qtm",
        "BFVD_best_ttm",
        "very_conserved_positions",
        "priority_structural_positions",
        "priority_acidic",
        "priority_basic",
        "spatial_cluster_count",
        "spatial_clusters",
    ]

    matrix = []

    for sid in TOP5:

        ir = i2[sid]
        mr = msa[sid]
        er = ext[sid]
        cr = cons[sid]

        sr = spatial.get(
            sid,
            {
                "count": "0",
                "clusters": "",
            },
        )

        row = {
            "structure_id":
                sid,

            "priority_class":
                ir.get(
                    "priority_class",
                    "",
                ),

            # These fields are intentionally blank in
            # the frozen canonical Stage M matrix.
            "experimental_score":
                "",

            "dominant_context":
                "",

            "assay_category":
                "",

            "independent_context_n":
                "",

            # Canonical Stage M uses the family-wide
            # noncandidate context sample span.
            "independent_context_sample_span":
                ir.get(
                    "noncandidate_context_sample_span",
                    "",
                ),

            "single_mean_plddt":
                mr.get(
                    "single_mean_plddt",
                    "",
                ),

            "msa_mean_plddt":
                mr.get(
                    "msa_mean_plddt",
                    "",
                ),

            "delta_mean_plddt":
                mr.get(
                    "delta_mean_plddt",
                    "",
                ),

            "single_ptm":
                mr.get(
                    "single_ptm",
                    "",
                ),

            "msa_ptm":
                mr.get(
                    "msa_ptm",
                    "",
                ),

            "delta_ptm":
                mr.get(
                    "delta_ptm",
                    "",
                ),

            "external_status":
                er.get(
                    "external_status",
                    "",
                ),

            "raw_sources":
                er.get(
                    "raw_sources",
                    "",
                ),

            "convincing_sources":
                er.get(
                    "convincing_sources",
                    "",
                ),

            "BFVD_best_target":
                er.get(
                    "BFVD_best_target",
                    "",
                ),

            "BFVD_best_qtm":
                er.get(
                    "BFVD_best_qtm",
                    "",
                ),

            "BFVD_best_ttm":
                er.get(
                    "BFVD_best_ttm",
                    "",
                ),

            "very_conserved_positions":
                cr.get(
                    "very_conserved_positions",
                    "",
                ),

            "priority_structural_positions":
                cr.get(
                    "priority_structural_positions",
                    "",
                ),

            "priority_acidic":
                cr.get(
                    "priority_acidic",
                    "",
                ),

            "priority_basic":
                cr.get(
                    "priority_basic",
                    "",
                ),

            "spatial_cluster_count":
                sr["count"],

            "spatial_clusters":
                sr["clusters"],
        }

        matrix.append(
            row
        )


    ###########################################################################
    # Write final5 evidence matrix
    ###########################################################################

    with open(
        out /
        "final5_evidence_matrix.tsv",
        "w",
        newline="",
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        writer.writerows(
            matrix
        )


    ###########################################################################
    # FASTA
    ###########################################################################

    write_fasta(
        out /
        "final_top5_representatives.faa",
        TOP5,
        structure_sequences,
    )

    write_fasta(
        out /
        "primary_top3_representatives.faa",
        TOP3,
        structure_sequences,
    )


    ###########################################################################
    # Summary
    ###########################################################################

    print(
        "TOP5 =",
        ",".join(TOP5),
    )

    print(
        "TOP3 =",
        ",".join(TOP3),
    )

    print(
        "FINAL5_ROWS =",
        len(matrix),
    )

    print(
        "SPATIAL_CLUSTER_FAMILIES =",
        ",".join(
            sorted(
                spatial
            )
        )
        if spatial
        else "NONE",
    )


if __name__ == "__main__":
    main()
