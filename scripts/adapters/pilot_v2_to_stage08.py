#!/usr/bin/env python3
"""Read-only Pilot_v2 to Stage08 compatibility adapter."""
import argparse
import ast
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

def rows(path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as f:
        data = list(csv.DictReader(f, delimiter="\t"))
    if not data:
        raise ValueError(f"Empty TSV: {path}")
    return data

def index(data, key):
    result = {}
    for row in data:
        k = row[key]
        if not k or k in result:
            raise ValueError(f"Missing/duplicate {key}: {k}")
        result[k] = row
    return result

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def fasta(path):
    result = {}
    key, seq = None, []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if key is not None:
                    result[key] = "".join(seq).upper()
                key, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if key is not None:
        result[key] = "".join(seq).upper()
    return result

def ca_range(path):
    values = []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
                values.append(float(line[60:66]))
    if not values:
        raise ValueError(f"No CA B-factors: {path}")
    return min(values), max(values)

def original_confidence(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    fn = next(
        (n for n in tree.body
         if isinstance(n, ast.FunctionDef)
         and n.name == "confidence_class"),
        None
    )
    if fn is None:
        raise ValueError(f"Missing confidence_class: {path}")
    namespace = {}
    module = ast.fix_missing_locations(
        ast.Module(body=[fn], type_ignores=[])
    )
    exec(compile(module, str(path), "exec"), namespace)
    return namespace["confidence_class"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot-root", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    pilot = args.pilot_root.resolve()
    out = args.out_dir.resolve()
    root = Path(__file__).resolve().parents[2]

    if out == pilot or pilot in out.parents or out in pilot.parents:
        raise ValueError("Output overlaps source research data")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be empty")

    paths = {
        "family": pilot / "StageE_Recurrent191_AllMemberSaturation/recurrent191_family_sequence_saturation.tsv",
        "members": pilot / "StageE_Recurrent191_AllMemberSaturation/recurrent191_all_members_sequence_evidence.tsv",
        "fmap": pilot / "StageF_Structure/structure_id_map.tsv",
        "qc": pilot / "StageF_Structure/structure191_qc.tsv",
        "fctx": pilot / "StageH3_CorrectedFamilyWideContext/recurrent191_family_context_corrected.tsv",
        "mctx": pilot / "StageH3_CorrectedFamilyWideContext/recurrent191_member_context_corrected.tsv",
        "seq": pilot / "StageF_Structure/structure_all191.faa",
        "allseq": pilot / "StageD_Recurrent191_FamilyAudit/recurrent191_all_members.faa",
    }

    families = index(rows(paths["family"]), "representative")
    members = rows(paths["members"])
    fmap = index(rows(paths["fmap"]), "structure_id")
    qc = index(rows(paths["qc"]), "structure_id")
    contexts = index(rows(paths["fctx"]), "family_rep")
    member_context = rows(paths["mctx"])
    seq = fasta(paths["seq"])
    allseq = fasta(paths["allseq"])

    if set(families) != set(contexts):
        raise ValueError("Family/context keys differ")
    if set(fmap) != set(qc) or set(fmap) != set(seq):
        raise ValueError("Structure IDs differ")

    grouped = defaultdict(list)
    for m in members:
        grouped[m["representative"]].append(m)

    if set(grouped) != set(families):
        raise ValueError("Family membership mismatch")
    if {m["member"] for m in members} != {
        m["member"] for m in member_context
    }:
        raise ValueError("Member context mismatch")

    qc_rep = {r["representative"]: r for r in qc.values()}
    if set(qc_rep) != set(families):
        raise ValueError("Structure family coverage mismatch")

    darkness = []
    allowed = {
        "KNOWN_SEQUENCE", "KNOWN_SEQUENCE_EVIDENCE", "STRICT_DARK", "DUF_ONLY",
        "PHROG_UNKNOWN", "UNRESOLVED_HOMOLOGY",
        "HOMOLOGY_UNRESOLVED"
    }

    for rep, row in families.items():
        group = grouped[rep]
        n = int(row["family_size"])
        statuses = Counter(m["sequence_status"] for m in group)

        if set(statuses) - allowed:
            raise ValueError(
                f"Unknown member status: {rep}: {statuses}"
            )

        known = int(row["known_members"])
        strict = int(row["strict_dark_members"])
        unresolved = sum(
            int(m["function_unresolved"]) for m in group
        )

        if (
            len(group) != n
            or int(contexts[rep]["family_size"]) != n
            or unresolved != n - known
        ):
            raise ValueError(f"Family count mismatch: {rep}")

        if (
            (statuses["KNOWN_SEQUENCE"] + statuses["KNOWN_SEQUENCE_EVIDENCE"]) != known
            or statuses["STRICT_DARK"] != strict
        ):
            raise ValueError(
                f"Annotation classification mismatch: {rep}"
            )

        if int(row["all_members_strict_dark"]) != int(strict == n):
            raise ValueError(f"Strict-dark flag mismatch: {rep}")

        converted = dict(row)
        converted.update(
            search_complete_members=str(n),
            unresolved_homology_members=str(
                statuses["UNRESOLVED_HOMOLOGY"]
                + statuses["HOMOLOGY_UNRESOLVED"]
            ),
            function_unresolved_members=str(unresolved),
            priority_class=qc_rep[rep]["priority_class"],
        )
        darkness.append(converted)

    converted_context = []
    for row in member_context:
        converted = dict(row)
        converted["window"] = "5"
        converted["anchor_within_window"] = row["anchor_within5"]
        converted["structural_neighbors_within_window"] = (
            row["structural_neighbors_within5"]
        )

        if converted["anchor_within_window"] not in ("0", "1"):
            raise ValueError(f"Invalid anchor: {row['member']}")

        converted_context.append(converted)

    confidence = original_confidence(
        root / "scripts/08_collect_structure_results_portable.py"
    )

    structures = []
    for sid, row in fmap.items():
        q = qc[sid]
        rep = row["representative"]
        sequence = seq[sid]

        if q["representative"] != rep or allseq.get(rep, "").rstrip("*") != sequence.rstrip("*"):
            raise ValueError(f"Structure sequence mismatch: {sid}")

        if (
            int(row["family_size"]) != int(families[rep]["family_size"])
            or int(q["family_size"]) != int(row["family_size"])
        ):
            raise ValueError(f"Structure family size mismatch: {sid}")

        if (
            q["priority_class"] != row["priority_class"]
            or q["length_match"] != "PASS"
        ):
            raise ValueError(f"Historical structure QC mismatch: {sid}")

        if (
            len(sequence) != int(row["length"])
            or len(sequence) != int(q["length_structure"])
        ):
            raise ValueError(f"Sequence length mismatch: {sid}")

        pdb = Path(q["pdb"])
        score = Path(q["score_json"])

        if not pdb.is_file() or not score.is_file():
            raise ValueError(f"Missing structure files: {sid}")

        lo, hi = ca_range(pdb)
        mean = float(q["mean_plddt"])
        if not lo - 0.1 <= mean <= hi + 0.1:
            raise ValueError(f"Structure confidence mismatch: {sid}")

        structures.append(dict(
            structure_id=sid,
            family_rep=rep,
            source_member=rep,
            family_size=row["family_size"],
            priority_class=q["priority_class"],
            length_aa=len(sequence),
            sequence_sha256=hashlib.sha256(
                sequence.encode()
            ).hexdigest(),
            pdb_file=str(pdb),
            pdb_sha256=sha(pdb),
            score_json=str(score),
            score_json_sha256=sha(score),
            mean_plddt=q["mean_plddt"],
            min_plddt=f"{lo:.2f}",
            max_plddt=f"{hi:.2f}",
            ptm=q["ptm"],
            confidence_class=confidence(mean),
            prediction_status="PREDICTED",
        ))

    schemas = {
        "06_darkness/family_darkness_audit.tsv":
            "representative family_size search_complete_members known_members strict_dark_members duf_only_members phrog_unknown_members unresolved_homology_members function_unresolved_members all_members_function_unresolved all_members_strict_dark family_sequence_status priority_class members_with_phrog_structural members_with_phrog_nucleic".split(),
        "07_family_context/family_context.tsv":
            "family_rep family_size sample_span samples candidate_count noncandidate_count all_anchor_context_count all_anchor_context_fraction candidate_anchor_context_count candidate_anchor_context_fraction noncandidate_anchor_context_count noncandidate_anchor_context_fraction all_context_sample_span noncandidate_context_sample_span dominant_anchor_category dominant_anchor_supported_members dominant_anchor_consistency all_anchor_categories structural_context_any_count structural_dense_ge2_count family_context_tier".split(),
        "07_family_context/member_context.tsv":
            "family_rep member sample length_aa candidate_original eggnog_status strict_unresolved_eggnog gff_seqid gff_local_index gff_order window neighbors_available anchor_within_window anchor_min_gene_distance anchor_categories anchor_ids structural_neighbors_within_window structural_dense_ge2 reconstructed_candidate noncandidate_exclusion_reason".split(),
        "08_structure/run_pilot_import/collected/structure_results.tsv":
            "structure_id family_rep source_member family_size priority_class length_aa sequence_sha256 pdb_file pdb_sha256 score_json score_json_sha256 mean_plddt min_plddt max_plddt ptm confidence_class prediction_status".split(),
    }

    datasets = [
        darkness,
        list(contexts.values()),
        converted_context,
        structures,
    ]

    for (name, columns), data in zip(schemas.items(), datasets):
        missing = set(columns) - set(data[0])
        if missing:
            raise ValueError(f"Missing fields in {name}: {missing}")

        dest = out / name
        dest.parent.mkdir(parents=True, exist_ok=True)

        with dest.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=columns,
                delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(data)

        print(f"CREATED={dest} ROWS={len(data)}")

    provenance = {
        "source_type": "historical Pilot_v2 reuse",
        "families": len(families),
        "members": len(members),
        "structures": len(structures),
        "window_mapping": "anchor_within5 -> window=5",
        "priority_source": "historical StageF structure QC",
        "new_gpu_predictions": False,
        "new_foldseek_searches": False,
        "source_sha256": {
            str(path.relative_to(pilot)): sha(path)
            for path in paths.values()
        },
    }

    (out / "historical_source_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n",
        encoding="utf-8"
    )

    print("STEP143_REAL_STAGE08_ADAPTER=PASS")

if __name__ == "__main__":
    main()
