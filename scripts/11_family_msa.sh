#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Stage J
# Family-level MSA preparation for the frozen Top10 recurrent families
#
# Reconstructed canonical workflow:
#
#   H3 recurrent191 member table
#        ↓
#   terminal "*" removal
#        ↓
#   exact amino-acid sequence collapse within each family
#   (first occurrence in H3 member order retained)
#        ↓
#   family representative moved to first position
#        ↓
#   MAFFT v7.526 --auto --thread 16
#        ↓
#   query-centered A3M conversion
#        ↓
#   ColabFold-safe header conversion
#
# Frozen Top10:
# F023 F008 F025 F048 F031 F005 F065 F135 F009 F001
#
# Output family order is structure-ID order for deterministic file generation:
# F001 F005 F008 F009 F023 F025 F031 F048 F065 F135
###############################################################################

BASE="${1:-Result/Annotation_saturation/Pilot_v2}"
THREADS="${2:-16}"

H3_MEMBER="$BASE/StageH3_CorrectedFamilyWideContext/recurrent191_member_context_corrected.tsv"
H3_FAMILY="$BASE/StageH3_CorrectedFamilyWideContext/recurrent191_family_context_corrected.tsv"

ALL_MEMBERS="$BASE/StageD_Recurrent191_FamilyAudit/recurrent191_all_members.faa"

OUT="$BASE/StageJ_FamilyMSA"
FAMILY_FASTA="$OUT/family_fasta"
MAFFT_OUT="$OUT/mafft"
A3M_OUT="$OUT/a3m"
CF_OUT="$OUT/a3m_colabfold"

mkdir -p \
    "$FAMILY_FASTA" \
    "$MAFFT_OUT" \
    "$A3M_OUT" \
    "$CF_OUT"

for f in \
    "$H3_MEMBER" \
    "$H3_FAMILY" \
    "$ALL_MEMBERS"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing/empty input: $f" >&2
        exit 1
    fi
done

TOP10=(
    F001
    F005
    F008
    F009
    F023
    F025
    F031
    F048
    F065
    F135
)

###############################################################################
# 1. Build canonical family FASTA files
###############################################################################

python - \
    "$H3_MEMBER" \
    "$H3_FAMILY" \
    "$ALL_MEMBERS" \
    "$FAMILY_FASTA" <<'PY'
from pathlib import Path
import sys
import pandas as pd

member_path = Path(sys.argv[1])
family_path = Path(sys.argv[2])
faa_path = Path(sys.argv[3])
outdir = Path(sys.argv[4])

TOP10 = [
    "F001",
    "F005",
    "F008",
    "F009",
    "F023",
    "F025",
    "F031",
    "F048",
    "F065",
    "F135",
]


def read_fasta(path):
    records = []
    header = None
    seq = []

    with path.open() as fh:
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
                seq.append(line.strip())

        if header is not None:
            records.append(
                (
                    header,
                    "".join(seq).rstrip("*"),
                )
            )

    return records


def write_fasta(records, path):
    with path.open("w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


members = pd.read_csv(
    member_path,
    sep="\t",
    dtype=str,
)

families = pd.read_csv(
    family_path,
    sep="\t",
    dtype=str,
)

all_records = read_fasta(
    faa_path
)

seq_by_id = {
    header: seq
    for header, seq in all_records
}

sid_to_rep = dict(
    zip(
        families["structure_id"],
        families["family_rep"],
    )
)

for sid in TOP10:

    if sid not in sid_to_rep:
        raise RuntimeError(
            f"{sid}: structure_id missing from H3 family table"
        )

    rep = sid_to_rep[sid]

    family_members = (
        members.loc[
            members["family_rep"] == rep,
            "member",
        ]
        .tolist()
    )

    if rep not in family_members:
        raise RuntimeError(
            f"{sid}: representative missing from member table: {rep}"
        )

    missing = [
        x
        for x in family_members
        if x not in seq_by_id
    ]

    if missing:
        raise RuntimeError(
            f"{sid}: {len(missing)} members missing from FASTA; "
            f"first={missing[:5]}"
        )

    # Exact AA sequence collapse.
    #
    # Important:
    # preserve H3 physical member order and retain the first
    # member encountered for each exact sequence.
    seen_sequences = set()
    unique_records = []

    for member in family_members:

        seq = seq_by_id[member].rstrip("*")

        if seq in seen_sequences:
            continue

        seen_sequences.add(seq)

        unique_records.append(
            (
                member,
                seq,
            )
        )

    retained_ids = [
        h
        for h, _ in unique_records
    ]

    if rep not in retained_ids:
        raise RuntimeError(
            f"{sid}: representative sequence was not retained "
            f"during exact-sequence collapse"
        )

    # Representative must be query / first sequence.
    unique_records = (
        [
            x
            for x in unique_records
            if x[0] == rep
        ]
        +
        [
            x
            for x in unique_records
            if x[0] != rep
        ]
    )

    outfile = (
        outdir /
        f"{sid}.family.faa"
    )

    write_fasta(
        unique_records,
        outfile,
    )

    print(
        sid,
        f"family={len(family_members)}",
        f"unique={len(unique_records)}",
        f"collapsed={len(family_members)-len(unique_records)}",
        f"rep={rep}",
        sep="\t",
    )
PY


###############################################################################
# 2. MAFFT
#
# Canonical production:
#   MAFFT v7.526
#   --auto
#   --thread 16
#
# Saved logs show L-INS-i for all ten families.
###############################################################################

if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: mafft not found in PATH" >&2
    echo "Activate the environment containing MAFFT v7.526." >&2
    exit 1
fi

echo "MAFFT:"
mafft --version 2>&1 || true

for SID in "${TOP10[@]}"
do
    FAA="$FAMILY_FASTA/${SID}.family.faa"
    ALN="$MAFFT_OUT/${SID}.aligned.faa"
    LOG="$MAFFT_OUT/${SID}.mafft.log"

    echo "[MAFFT] $SID"

    mafft \
        --auto \
        --thread "$THREADS" \
        "$FAA" \
        > "$ALN" \
        2> "$LOG"
done


###############################################################################
# 3. Convert aligned FASTA to query-centered A3M
###############################################################################

python - \
    "$MAFFT_OUT" \
    "$A3M_OUT" <<'PY'
from pathlib import Path
import sys

indir = Path(sys.argv[1])
outdir = Path(sys.argv[2])

TOP10 = [
    "F001",
    "F005",
    "F008",
    "F009",
    "F023",
    "F025",
    "F031",
    "F048",
    "F065",
    "F135",
]


def read_fasta(path):
    records = []
    header = None
    seq = []

    with path.open() as fh:
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


def write_fasta(records, path):
    with path.open("w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            out.write(seq + "\n")


def msa_to_a3m(records):

    if not records:
        return []

    seqs = [
        seq
        for _, seq in records
    ]

    width = len(seqs[0])

    if any(
        len(seq) != width
        for seq in seqs
    ):
        raise RuntimeError(
            "MSA sequences have inconsistent alignment widths"
        )

    converted = [
        [
            header,
            [],
        ]
        for header, _ in records
    ]

    query = seqs[0]

    for col in range(width):

        q = query[col]

        if q != "-":

            for i, seq in enumerate(seqs):

                aa = seq[col]

                if aa == "-":
                    converted[i][1].append("-")
                else:
                    converted[i][1].append(
                        aa.upper()
                    )

        else:

            # Query-gap alignment columns become lowercase
            # insertions in homolog sequences.
            for i in range(
                1,
                len(seqs),
            ):

                aa = seqs[i][col]

                if aa == "-":
                    continue

                converted[i][1].append(
                    aa.lower()
                )

    return [
        (
            header,
            "".join(seq),
        )
        for header, seq in converted
    ]


for sid in TOP10:

    src = (
        indir /
        f"{sid}.aligned.faa"
    )

    dst = (
        outdir /
        f"{sid}.a3m"
    )

    records = read_fasta(
        src
    )

    converted = msa_to_a3m(
        records
    )

    write_fasta(
        converted,
        dst,
    )

    print(
        sid,
        f"n={len(converted)}",
        f"query_length={len(converted[0][1])}",
        sep="\t",
    )
PY


###############################################################################
# 4. Make ColabFold-safe A3M headers
#
# query:
#   F001
#
# homologs:
#   F001_h0001
#   F001_h0002
#   ...
###############################################################################

python - \
    "$A3M_OUT" \
    "$CF_OUT" <<'PY'
from pathlib import Path
import sys
import csv

indir = Path(sys.argv[1])
outdir = Path(sys.argv[2])

TOP10 = [
    "F001",
    "F005",
    "F008",
    "F009",
    "F023",
    "F025",
    "F031",
    "F048",
    "F065",
    "F135",
]


def read_fasta(path):
    records = []
    header = None
    seq = []

    with path.open() as fh:
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


header_map = []

for sid in TOP10:

    src = (
        indir /
        f"{sid}.a3m"
    )

    dst = (
        outdir /
        f"{sid}.a3m"
    )

    records = read_fasta(
        src
    )

    with dst.open("w") as out:

        for i, (
            original_header,
            seq,
        ) in enumerate(records):

            if i == 0:
                new_header = sid
            else:
                new_header = (
                    f"{sid}_h{i:04d}"
                )

            out.write(
                f">{new_header}\n"
            )

            out.write(
                seq + "\n"
            )

            header_map.append(
                (
                    sid,
                    new_header,
                    original_header,
                )
            )

map_path = (
    outdir /
    "header_map.tsv"
)

with map_path.open(
    "w",
    newline="",
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        [
            "structure_id",
            "new_header",
            "original_header",
        ]
    )

    writer.writerows(
        header_map
    )

print(
    "HEADER_MAP_ROWS =",
    len(header_map)
)
PY


echo
echo "=================================================="
echo "STAGE J COMPLETE"
echo "=================================================="
echo "OUTPUT=$OUT"
