# Methods — implementation-linked overview

This is a description of algorithms actually present in the `scripts/` directory, not an independent proof of biological performance.

| Stage | Implementation and decision rule | Interpretation |
|---|---|---|
| 00 | `00_validate_input_portable.py`; FASTA/GFF/eggNOG schema and Prodigal ID cross-check, `00_capture_run_provenance.py` | data and runtime integrity |
| 01 | `01_define_unresolved_candidates.py`; missing functional annotation fields + `min_aa` threshold | *annotation unresolved*, not necessarily novel |
| 02 | `02_extract_replication_context.py`; textual viral replication/NA metabolism anchor and +/- `window` ORF indices | hypothesis-generating neighborhood |
| 03 | `03_global_family_clustering.sh`: MMseqs easy-cluster at min sequence identity **0.30**, coverage **0.70**, cov-mode **0** | sequence similarity family groups |
| 04 | `04_select_recurrent_families_portable.py`; unresolved Family selection at >=2 members (ge2 TSV) | within-input recurrence; not proof of cross-sample recurrence |
| 05 | `05_annotation_saturation_portable.sh`, `05_build_member_evidence.py`; Swiss-Prot / PHROG profile / Pfam | sequence/domain functional evidence |
| 06 | `06_family_darkness_audit_portable.py`; integrates member annotations | unresolved family evidence |
| 07 | `07_familywide_context_portable.py`; family-wide genome context | contextual association |
| 08 | `08_prepare_structure_queries_portable.py`; ColabFold optional external run, verified collection | predicted structure, not experimentally resolved |
| 09 | `09_foldseek_search_portable.py` and `09_external_integration_portable.py`; Foldseek PDB, AF_SwissProt and BFVD | external structural similarity |
| 10 | `10_integrated_ranking_portable.py`, `10_ranked_shortlists_portable.py`; confidence-aware ranking + shortlists | prioritization, not diagnosis |
| 11 | `11_family_msa_portable.py`, prediction/QC and adaptive wrappers | optional Family MSA/structure support |
| 12 | `12_msa_structure_conservation*.py`, `12_charged_spatial_clusters_exploratory.py` | conserved residues and exploratory charged surface clusters; methodology depends on workflow version |
| 13 | `13_*portable.py`; preliminary/final-evidence bridge and explicit pending state | incomplete external evidence must stay pending |

Stages 00–08 are driven by `run_pipeline.py`; `run_pipeline_full.py` hands Stage 08 outputs to the Stage 09–13 adaptive orchestrator, but currently **does not automatically run all MSA/Stage 12–13 biological analyses**. Historical Pilot-only scripts remain for compatibility and should not be substituted for portable modules in a new run.

## Scientific cautions

- A conserved Lys/Arg/Asp/Glu/His cluster can reflect nucleic-acid binding, global electrostatic surface, or other roles; it is **not** a catalytic-site call.
- MSA structures can differ from single-sequence structures; keep search provenance separate.
- The historical Pilot spatial-cluster definition is **not identical** to the current exploratory Stage 12 cluster definition. Cluster counts are not directly numerically reproducible across these algorithms.
- Biochemical nuclease activity, enzyme identity and substrate specificity require experimental validation.
