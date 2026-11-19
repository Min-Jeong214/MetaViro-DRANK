# MetaViro-DRANK official CPU-first entrypoint

This guide is for **new marine viral protein cohorts**, not the historical Pilot results. `run_metaviro.py` is the recommended entrypoint for a reproducible **structure-optional preliminary screen**. The legacy `run_pipeline_full.py` remains a distinct, structure-dependent advanced analysis and is *not* silently invoked.

## Required inputs

For each sample, provide `proteins.faa`, `genes.gff` (Prodigal CDS coordinates), and `annotations.emapper.annotations` (eggNOG-mapper). The CLI requires the `metaviro-drank` Conda environment and separately installed Swiss-Prot MMseqs2, PHROG MMseqs2 profile and annotation TSV, and pressed Pfam HMM databases. These paths are supplied in `config.local.yaml`, following `config/example.yaml`.

The `--normalize-ids` option **copies and namespaces** sample protein IDs, GFF coordinates and eggNOG query IDs together, preserving unmodified source files. It rejects invalid within-sample relationships; it does not fix biological mapping mismatches.

## From raw protein/GFF/eggNOG cohort to contextual candidates

```bash
# Fill config.local.yaml and set an UNUSED output_dir and an UNUSED screen-out-dir.
conda activate metaviro-drank
python -B run_metaviro.py run \
  --config config.local.yaml \
  --screen-out-dir results/run01_context \
  --normalize-ids --limit 30 --mode check

python -B run_metaviro.py run \
  --config config.local.yaml \
  --screen-out-dir results/run01_context \
  --normalize-ids --limit 30 --mode run
```

This runs Stage 00–07 (including Swiss-Prot, PHROG and Pfam annotation) exactly once, then Stage 08 representative preparation, structure-reuse inventory, and CPU-only context-first candidate review. `--normalize-ids` is optional when all input protein IDs are already unique across samples. For large cohorts Stage 05 HMMER can take substantial time. Stop on failure and inspect the current output; do not rerun the same target paths without first deciding how to resume safely.

## From existing Stage 07 results (no repeat searches)

```bash
python -B run_metaviro.py resume \
  --workflow-output /path/to/completed/workflow \
  --screen-out-dir results/resume01_context \
  --limit 30 --mode check

python -B run_metaviro.py resume \
  --workflow-output /path/to/completed/workflow \
  --screen-out-dir results/resume01_context \
  --limit 30 --mode run
```

Both modes support optional `--historical-structures /path/to/verified/structure_results.tsv`. This is **never required** for a new dataset; where given, exact-sequence SHA and source PDB integrity determine potential reuse. Family identifiers are generated per cohort, not mapped by their numeric suffixes to the Pilot.

## Outputs and evidential boundaries

- `run_metaviro_provenance.json`: completion state, counts, and provenance hash.
- `review/08_prepared/structure_id_map.tsv`: family IDs and representative sequences.
- `review/08_partial_plan/partial_structure_inventory.tsv`: all families with optional structure reuse status.
- `review/08_partial_plan/new_prediction_review_shortlist.tsv`: review shortlist, **not** predicted structures.
- `review/10_context_screen/all_family_context_screen.tsv`: all families and structure status.
- `review/10_context_screen/dark_replication_context_review.tsv`: context-supported unresolved families.
- `review/context_screen_run_provenance.json`: evidence source hashes and partial-analysis state.

**No GPU, ColabFold or Foldseek is run here.** These outputs are **preliminary contextual priorities**, *not* confirmed nuclease/primase/helicase function, Stage 10's fully structure-informed confidence rank, Stage 13 completed evaluation, or biochemical validation. Families without structure remain in results and must not be treated as negative structural evidence.

For Stage 09–13 advanced structural analysis see [METHODS.md](METHODS.md), [USAGE.md](USAGE.md), and [VALIDATION_STATUS.md](VALIDATION_STATUS.md). Study-specific comparisons and validation results belong in the associated scientific manuscript and Zenodo record, not as required run inputs.
