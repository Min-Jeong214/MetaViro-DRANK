# Execute MetaViro-DRANK on your own marine viral dataset

## A. Check sample format without external search DBs

```bash
python run_pipeline.py --config config.local.yaml --through 02 --check
```

See [INPUT_FORMAT.md](INPUT_FORMAT.md) and [config/example.yaml](../config/example.yaml). Input sample directory contains proteins+GFF+eggNOG annotations; raw reads and bare contig FASTA are not accepted.

## B. Run upstream stages 00–08

```bash
# Check all stage-05 databases and sample consistency first
python run_pipeline.py --config config.local.yaml --through 08 --check

# Run only through Stage 07, if you want Family/context results before GPU work
python run_pipeline.py --config config.local.yaml --through 07

# For the initial Stage 08 query preparation, use a NEW output_dir in YAML
python run_pipeline.py --config config.local.yaml --through 08 --structure-mode prepare
```

`output_dir` already written by an earlier run will contain provenance and logs; do not manually mix stages or reuse stale files to imply a new run. Stage 08 `prepare` **does not** produce collected PDBs.

Stage 08 collected structure results require either `--structure-mode predict` (configured external ColabFold environment) or `--structure-mode collect --structure-prediction-dir /path/to/matching_predictions`. Confirm input family→structure mapping before collection; do not reuse unrelated predictions from another dataset.

## C. Unified Stage 00–13 entrypoint (current supported automation boundary)

```bash
# New samples; this command can do fresh Foldseek searches after Stage 08 collection.
python run_pipeline_full.py \
  --config config.local.yaml \
  --db-root /absolute/path/to/foldseek_dbs \
  --structure-mode collect \
  --structure-prediction-dir /absolute/path/to/matching_colabfold_results \
  --out-dir /absolute/path/to/new_full_output \
  --mode check

# After reviewing check results and ensuring all dependencies exist,
# use a NEW --out-dir and --mode run.
```

`--mode check` validates new-data Stage 00–08 inputs/databases and structural prediction source, but downstream search validation is deferred until structure collection. A check PASS is **not** proof that the complete run will succeed.

Current full entrypoint automatically produces Stage 09–10 ranking and adaptive Stage 11–13 status/hand-off, but **not** all of the optional MSA GPU predictions, new MSA-structure external searches or Stage 12–13 residue biology. For verified Stage 08 reuse use `--existing-stage08-output /path` with compatible `--reuse-search-root /path`; the match is SHA/provenance-sensitive.

## D. Outputs and documentation

See [OUTPUTS.md](OUTPUTS.md) for directory structure and table semantics, [DATABASES.md](DATABASES.md) for DB layouts and [METHODS.md](METHODS.md) for algorithms. Do not claim completed final candidate ranking when predictions remain pending.

## E. Limitations of the first deployable release

- `environment.yml` has not been fresh-install accepted or fully locked across Ubuntu platforms.
- Bare contigs/raw reads require upstream gene prediction and annotation; this tool currently starts from compatible protein/GFF/eggNOG files.
- Structure predictions and optional Stage 11–13 functions have separate requirements and may remain pending.
- No biochemical nuclease activity is predicted with validated assay accuracy.
