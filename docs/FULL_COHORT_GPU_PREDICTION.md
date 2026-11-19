# Optional full-cohort MSA and GPU structure prediction

This is a **generic, checkpointed workflow** for predicting protein-family representative structures. No sample IDs, local research paths, or fixed number of families are built into the scripts. The number of families is taken from the user's Stage 08 structure map. Biological databases, proprietary data and generated PDBs are not included in this repository.

## Requirements

- Completed Stage 05–07 outputs: `recurrent_members.tsv` and `recurrent_members.faa`.
- Prepared Stage 08 outputs: `structure_id_map.tsv` and `representative_queries.faa`, derived from the **same cohort**.
- Main Conda environment with Python and MAFFT for family MSA preparation.
- Separate working GPU ColabFold environment. A tested example was Python 3.11.17, ColabFold 1.6.1 and JAX 0.5.3 on an NVIDIA RTX 3050 (8 GiB). These versions/hardware are **not guaranteed to support every query**; longer models may exceed available VRAM.

## 1. Set paths for your own analysis

```bash
REPO=/path/to/MetaViro-DRANK
MEMBERS=/path/to/your/results/recurrent_members.tsv
PROTEINS=/path/to/your/results/recurrent_members.faa
MAP=/path/to/your/results/structure_id_map.tsv
QUERIES=/path/to/your/results/representative_queries.faa
MSA=/path/to/new/full_cohort_msa_directory
GPU=/path/to/new/gpu_predictions_directory
LOGS=/path/to/your/logs
mkdir -p "$LOGS"
cd "$REPO"
```

Each output path should be distinct from all inputs and any existing experiment's results. Do not copy data into the source repository.

## 2. Check and prepare MSAs for all mapped families

```bash
python -B scripts/11_full_cohort_msa_batch_portable.py \
  --members "$MEMBERS" --fasta "$PROTEINS" \
  --structure-map "$MAP" --out-dir "$MSA" \
  --chunk-size 10 --threads 8 --mode check

python -B scripts/11_full_cohort_msa_batch_portable.py \
  --members "$MEMBERS" --fasta "$PROTEINS" \
  --structure-map "$MAP" --out-dir "$MSA" \
  --chunk-size 10 --threads 8 --mode run \
  2>&1 | tee "$LOGS/full_cohort_msa.log"
```

Completed chunks are reused on rerun, with SHA-256 input matching. Failed/partial chunks are **not overwritten automatically**: diagnose `chunks/chunk_NNNN/run.log`, archive the *incomplete chunk* outside the active output tree, and rerun. The final `a3m_colabfold/`, `mafft/`, `family_msa_summary.tsv` and `header_map.tsv` cover every mapped family, even those with a single unique sequence after deduplication.

## 3. Check GPU inputs, then test a small batch

```bash
# Optional strict reuse only when A3M and query SHA-256 match exactly:
REUSE_ARGS=()
# REUSE_ARGS=(--reuse-prediction /path/to/existing/completed/prediction_directory)

conda run --no-capture-output -n metaviro-structure \
  python -B scripts/11_gpu_prediction_batch_portable.py \
  --structure-map "$MAP" --queries "$QUERIES" \
  --msa-dir "$MSA" --out-dir "$GPU" \
  "${REUSE_ARGS[@]}" --mode check

conda run --no-capture-output -n metaviro-structure \
  python -u -B scripts/11_gpu_prediction_batch_portable.py \
  --structure-map "$MAP" --queries "$QUERIES" \
  --msa-dir "$MSA" --out-dir "$GPU" \
  "${REUSE_ARGS[@]}" --max-jobs 3 --mode run \
  2>&1 | tee "$LOGS/gpu_pilot.log"
```

If strict reuse fails because the A3M differs (`REUSE_HASH_OR_PDB_INVALID`), **do not bypass the integrity check**. The simplest safe remedy is to omit `--reuse-prediction` and generate a fresh model with the current A3M. Alternatively add `--skip-incompatible-reuse` to the existing command: it prints `REUSE_INCOMPATIBLE_SKIPPED=<structure_id>:<reason>`, ignores only that incompatible reuse request, and schedules a fresh prediction from the current cohort MSA while preserving the old prediction. The default remains fail-closed. It never falsely marks an incompatible model as `REUSED`. Check or run mode uses the same policy. Only one prediction of a family should be launched in a given batch.

## 4. Continue remaining GPU predictions

After the pilot, verify GPU memory, model length and output correctness before starting all pending jobs:

```bash
nohup conda run --no-capture-output -n metaviro-structure \
  python -u -B scripts/11_gpu_prediction_batch_portable.py \
  --structure-map "$MAP" --queries "$QUERIES" \
  --msa-dir "$MSA" --out-dir "$GPU" \
  "${REUSE_ARGS[@]}" --mode run \
  > "$LOGS/gpu_batch_driver.log" 2>&1 < /dev/null &
echo "GPU_BATCH_PID=$!"
tail -f "$LOGS/gpu_batch_driver.log"
```

To continue after interruptions, run the same command again. The runner skips verified `COMPLETED` outputs and valid explicit `REUSED` outputs; output folders, per-family logs and failed attempts are preserved. An explicit `--retry-failed` uses a new `attempt_NNN` directory instead of overwriting previous files. The runner acquires a lock to prevent simultaneous work in the same GPU output directory.

**Single-unique-sequence handling:** families without at least two distinct sequences are recorded as `SINGLE_UNIQUE_SEQUENCE` and **not** given artificial homologs. They require a separately verified single-sequence prediction method. Extremely long queries can exhaust GPU VRAM; a failed query is retained as `FAILED` and does not invalidate completed predictions.

## 5. Check the manifest and downstream readiness

```bash
python - "$GPU/prediction_manifest.tsv" <<'PY'
import csv,sys
from collections import Counter
with open(sys.argv[1]) as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
print('FAMILIES=',len(rows))
print('STATUS=',dict(Counter(r['status'] for r in rows)))
PY
```

This extension does **not** claim that Foldseek searching, spatial analysis, biochemical candidate selection or final family-catalog integration has run for all families. Those stages must consume the verified PDBs through their own input/integrity contracts. A predicted structure or a charged cluster does not prove protein function.

## Repository release boundaries

GitHub holds general-purpose executable code, documentation and small synthetic tests. Real research cohort files, statistical validation results, structure models, raw hits and local installation paths remain outside the software repository. Scientific analyses can be reported separately in the manuscript and dedicated validation archive.
