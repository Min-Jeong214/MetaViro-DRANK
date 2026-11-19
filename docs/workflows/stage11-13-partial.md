<!-- MetaViro-DRANK release snapshot (2026-10-27) -->
# Stage 11–13: partial MSA-evidence workflow

The current development pipeline supports selective re-use of **provenance-checked historical** MSA predictions. Stage 11 must retain both the actual source A3M and source predicted-PDB identity; a newly computed alignment is not retroactively the input to a historical prediction.

## Components and source-code entry points

| Task | Entry point | Scope |
| --- | --- | --- |
| Historical MSA prediction import | `scripts/11_import_historical_reuse_portable.py` | Validate prediction index, A3M, PDB, associated hashes |
| Stage11 adaptive collect | `scripts/11_13_adaptive_orchestrator_portable.py` | 9 verified historical predictions + 6 pending |
| Stage12/13 preliminary partial evidence | `scripts/12_13_partial_historical_portable.py` | Conservation + exploratory 3D evidence for verified subset |
| Stage K2 Foldseek | `scripts/09_foldseek_search_portable.py` | Run on verified MSA-based PDB input, separately from Stage09 single-sequence PDB input |
| K2 evidence bridge | `scripts/13_k2_external_evidence_bridge_portable.py` | Validate external structure matrix and source files before partial merger |
| Stage13 K2 evidence merge | `scripts/13_partial_k2_merge_portable.py` | Verify source hashes and transfer `k2_` evidence without overwriting Stage10 |

The STEP164/165 Pilot K2 query preparation and search orchestration were executed as verified, dataset-specific adapters; a generic one-command K2 orchestration has **not** passed fresh-input acceptance. Do not treat the one-off WSL Pilot shell script as a portable module.

Inspect installed CLI arguments with `python -B <script> --help`; database paths, candidate sequences, prediction provenance, and output directories are user-supplied, not bundled. Do not use historical paths as defaults for new datasets.

## Partial integration CLI

```bash
python -B scripts/13_partial_k2_merge_portable.py \
  --stage162-output /path/to/valid/stage162/downstream \
  --stage165-output /path/to/valid/msa-k2-search \
  --stage09-10-output /path/to/valid/ranked-stage09-10 \
  --out-dir /path/to/NEW_OUTPUT_DIRECTORY \
  --mode check

# Following a successful read-only check, run with a fresh output directory:
python -B scripts/13_partial_k2_merge_portable.py \
  --stage162-output /path/to/valid/stage162/downstream \
  --stage165-output /path/to/valid/msa-k2-search \
  --stage09-10-output /path/to/valid/ranked-stage09-10 \
  --out-dir /path/to/ANOTHER_NEW_OUTPUT_DIRECTORY \
  --mode run
```

Key outputs:

- `stage13_partial_evidence_15.tsv`: all 15 candidate IDs, 9 verified with `k2_` columns plus 6 explicitly pending.
- `stage13_local_k2_verified_subset.tsv`: nine structures with both local Stage12 and K2 structural evidence.
- `stage13_partial_evidence_provenance.json`: source checksums and completion flags.

The old Stage10 `external_status` is preserved for **single-sequence structures**. The `k2_external_status` is the **MSA-based structure** assessment, not an overwrite of the original column.

## Validation gates

1. All listed input files exist and match the recorded SHA-256 provenance. Any mismatch is a hard failure, not a prompt to edit the source hash.
2. Selected-family IDs remain consistent (15 = 9 analyzed + 6 pending); Stage10 rank, candidate identifiers and assay categories are not silently recomputed.
3. The K2 evidence is linked to the exact MSA predicted PDB, not a similarly named single-sequence structure.
4. Run Stage159, Stage162 and Stage166 synthetic tests; inspect `PASS` output. The synthetic tests validate pipeline logic, **not biological function**.
5. Historical Pilot Top5 and Top3 are audit reference points, not newly selected final outputs.

## Out of scope

- No automatic 6-family GPU prediction, reassembly, or additional Foldseek is implied by partial-mode operation.
- No default interpretation of active-site chemistry from conserved charged clusters.
- No function confirmation from BFVD best-hit ID alone. Always retrieve the exact DB release annotation when needed.
- Full Stage13 final comparison and release acceptance remain outstanding.
