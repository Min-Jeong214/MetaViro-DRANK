# Output directory and interpretation

Current runner creates output under YAML `output_dir`, following these main paths:

| Stages | Important paths | What to inspect |
|---|---|---|
| 00 | `provenance/run_<id>.json`, `logs/` | CLI input hashes, parameters, tools, logs |
| 01 | `01_unresolved/<sample>.tsv` | `protein_id`, length, `eggnog_status`, `strict_unresolved`, `ge_min_aa` |
| 02 | `02_context/<sample>.tsv`, `<sample>_anchors.tsv` | candidate protein, viral replication anchors, ORF offsets |
| 03 | `03_clustering/viral_families_cluster.tsv` | representative-member MMseqs relationship |
| 04 | `04_recurrent/recurrent_families_{all,ge2,ge3,ge4,ge5}.tsv` | repeated families; Stage 05 uses ge2 |
| 05 | `05_annotation/recurrent_members.{faa,tsv}`, `member_level_evidence.tsv`, `search/` | member sequences and Swiss-Prot/PHROG/Pfam results |
| 06 | `06_darkness/family_darkness_audit.tsv` | family-level evidence status |
| 07 | `07_family_context/member_context.tsv`, `07_family_context/family_context.tsv` | family/member genome context tables |
| 08 | `08_structure/run_<id>/prepared/`, `/collected/structure_results.tsv` | structure queries and collected prediction integrity |
| 09–10 | *full-run output* `/downstream/09_10/09_external/`, `/10_ranking/confidence_aware_ranking.tsv`, `/10_shortlists/` | Foldseek evidence status and ranked shortlist |
| 11–13 | *full-run output* `/downstream/11_13_adaptive/`, `/13_adaptive/` | prediction/evidence statuses, may be pending |

The unified full runner wraps the upstream directory with a separate `--out-dir`; its `downstream` and `orchestration_provenance.json` are written under that output folder. Existing output paths are intentionally refused; use a new output directory for each run.

Do **not** interpret:
- `strict_unresolved` as a guaranteed biologically novel protein;
- `nuclease` assay category as confirmed nuclease activity;
- Foldseek structural similarity as validated enzyme function;
- empty or pending Stage 11–13 fields as zero biological signal;
- a previously published, cohort-specific shortlist as a new independent selection.

A reportable candidate should include protein/family ID, stage and method provenance, all evidence statuses, sequence/structure identifiers and proposed validation assay. Results must be accompanied by the relevant software/database versions and input fingerprints.
