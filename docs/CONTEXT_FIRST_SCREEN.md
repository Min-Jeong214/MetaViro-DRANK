# Context-first screening without full structural predictions

MetaViro-DRANK is intended for **new marine viral protein cohorts**, not for reproducing a fixed Pilot family list. This optional CPU-only screening path lets a user continue after validated Stages 00–07 even if none of the families has a predicted structure. It does **not** execute Stage 09 Foldseek, the published Pilot I2 integrated scoring, or Stage 13 final candidate selection.

## Prerequisites

- Finished workflow output from `run_pipeline.py --config <config> --through 07` containing nonempty:
  `05_annotation/recurrent_members.tsv`, `05_annotation/recurrent_members.faa`, `06_darkness/family_darkness_audit.tsv`, `07_family_context/{family_context,member_context}.tsv`.
- Existing `scripts/08_prepare_structure_queries_portable.py`, `scripts/08_partial_structure_plan_portable.py` and `scripts/10_context_first_screen_portable.py`.
- Python 3.11 and no extra database or GPU dependency for this step.

## One-command analysis

```bash
python -B run_context_screen.py \
  --workflow-output /analysis/my_cohort/workflow \
  --out-dir /analysis/my_cohort/context_screen_01 \
  --limit 30 \
  --mode check

python -B run_context_screen.py \
  --workflow-output /analysis/my_cohort/workflow \
  --out-dir /analysis/my_cohort/context_screen_01 \
  --limit 30 \
  --mode run
```

The output folder must not exist. `--mode check` writes no files; `run` creates a new directory. The optional `--historical-structures /path/to/collected/structure_results.tsv` allows SHA/length-verified reuse **assessment**, not automatic reuse/annotation transfer. Omit it entirely for a never-seen dataset. The family count is determined from the input, never pinned to a Pilot number.

## Output files

```
context_screen_01/
  08_prepared/structure_id_map.tsv
  08_prepared/representative_queries.faa
  08_partial_plan/partial_structure_inventory.tsv
  08_partial_plan/historical_reuse_candidates.tsv
  08_partial_plan/new_prediction_review_shortlist.tsv
  08_partial_plan/new_prediction_review_queries.faa
  08_partial_plan/partial_structure_plan_provenance.json
  10_context_screen/all_family_context_screen.tsv
  10_context_screen/dark_replication_context_review.tsv
  10_context_screen/context_first_provenance.json
  context_screen_run_provenance.json
```

`all_family_context_screen.tsv`: **every family**, including those without structural predictions and those without replication-anchor support.

`dark_replication_context_review.tsv`: subset with darkness class A/B/C and an identified nearby replication/nucleotide-metabolism category supported by noncandidate family members in Stage 07. The `best_noncandidate_anchor_category`, `category_supported_members`, and `category_supported_samples` columns show specific *neighboring gene context*, not function assignment to the dark protein.

`context_review_order` is only a deterministic **review order** (independent sample support, supported members, overall context span and sample recurrence); it is **not** a calibrated activity score, the original I2 integrated/experimental ranking, or a final enzyme prediction.

`structure_status=VERIFIED_HISTORICAL_REUSE_CANDIDATE` only means the previous PDB and score JSON match a new representative sequence exactly and their file hashes pass. No structure is considered incorporated without subsequent verified Stage 08 collection.

## Scientific and release limitations

- **No hard-coded family identifiers, fixed cohort-specific family counts, or external historical-project input paths in the general screening workflow.**
- An anchor in the +/- N gene neighborhood supports association with a replication/nucleotide-metabolism region; a nearby nuclease does **not** establish nuclease catalytic activity for an uncharacterized protein.
- Missing structure, no Foldseek hit, and no noncandidate anchor are not identical to negative biochemical evidence. Missing downstream evidence stays missing/pending.
- The Stage 10 structure-aware I2 rank requires full verified Stage 08/09 structures; this CPU screening does not create fake `mean_plddt`, TM-score, or `external_status`.
- Full end-to-end Stage 08–13 partial-prediction/partial-external integration is still a **separate incomplete development task**. GitHub release must describe the limitation.
- For authorizations/license, follow the institute's disclosure and rights policy; this code does not include private raw FASTQ/FASTA, databases, or Pilot paper results.
