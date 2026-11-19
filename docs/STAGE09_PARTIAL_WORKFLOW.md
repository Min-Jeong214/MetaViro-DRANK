# Stage 09: structure-optional Foldseek workflow

This is a **general-purpose** extension for new marine viral datasets, not a Pilot-only procedure.

## Preconditions

- Stage 05–07 complete; Stage 08 structure planning completed by `run_context_screen.py` or `run_metaviro.py resume`.
- `08_partial_plan/partial_structure_inventory.tsv`, its provenance JSON and `08_prepared/representative_queries.faa` must exist.
- Optional existing PDBs from independently verified Stage 08 results (`--predicted-structures`) or exact-sequence-matched historical reuse candidates. Neither is mandatory.

## 1. Prepare the structure-verified subset

```bash
python scripts/09_partial_structure_evidence_portable.py prepare \
  --plan-dir results/my_screen/review/08_partial_plan \
  --queries results/my_screen/review/08_prepared/representative_queries.faa \
  --out-dir results/my_stage09/prepare --mode check

python scripts/09_partial_structure_evidence_portable.py prepare \
  --plan-dir results/my_screen/review/08_partial_plan \
  --queries results/my_screen/review/08_prepared/representative_queries.faa \
  --out-dir results/my_stage09/prepare --mode run
```

To use new structures collected via Stage 08, add `--predicted-structures /absolute/path/to/structure_results.tsv` to both commands. If verified historical models exist in the Stage179 inventory, they are optional; the output marks their origin and does not call them new predictions.

`structure_results_for_search.tsv` is **only** the verified PDB subset. `all_family_structure_status.tsv` keeps every Family, including all pending predictions.

## 2. Perform external searches only when PDBs are available

The DB paths are **user-configured**, never built into the script:

```bash
INPUT=results/my_stage09/prepare/structure_results_for_search.tsv
python scripts/09_foldseek_search_portable.py --structures "$INPUT" \
  --target-db /database/path/PDB --source pdb_e1 \
  --out-dir results/my_stage09/search/pdb_e1 --threads 8
python scripts/09_foldseek_search_portable.py --structures "$INPUT" \
  --target-db /database/path/AF_SwissProt --source afsp \
  --out-dir results/my_stage09/search/afsp --threads 8
python scripts/09_foldseek_search_portable.py --structures "$INPUT" \
  --target-db /database/path/BFVD --source bfvd \
  --out-dir results/my_stage09/search/bfvd --threads 8
python scripts/09_external_integration_portable.py \
  --structures "$INPUT" --pdb-e1-dir results/my_stage09/search/pdb_e1 \
  --afsp-dir results/my_stage09/search/afsp \
  --bfvd-dir results/my_stage09/search/bfvd \
  --out-dir results/my_stage09/external
```

Do not re-execute completed Foldseek searches on error; inspect their provenance and resume from the missing stage. A database search is not a biochemical function assay.

## 3. Verify Foldseek integration and merge with *all* Family statuses

```bash
python scripts/09_partial_structure_evidence_portable.py merge \
  --prepared-dir results/my_stage09/prepare \
  --external-dir results/my_stage09/external \
  --out-dir results/my_stage09/merged --mode check

python scripts/09_partial_structure_evidence_portable.py merge \
  --prepared-dir results/my_stage09/prepare \
  --external-dir results/my_stage09/external \
  --out-dir results/my_stage09/merged --mode run
```

For an entirely novel dataset with **zero available structures**, omit `--external-dir` and the search commands: the valid merge result retains all families as `NOT_EVALUATED_NO_STRUCTURE`.

## Outputs and interpretation

- `prepare/all_family_structure_status.tsv`: all families with `PREDICTED_VERIFIED` or `PREDICTION_PENDING`.
- `prepare/structure_results_for_search.tsv`: only PDB-verified search inputs. This table may have a header but no data rows.
- `merged/all_family_stage09_evidence.tsv`: all families and all 38 Stage09 evidence columns (prefixed `stage09_`), with search result status.
- `merged/stage182_merge_provenance.json`: counts and SHAs.

**Missing PDB never means NO_RAW_HIT_E1.** The external status is empty and the `stage09_evidence_state` is `NOT_EVALUATED_NO_STRUCTURE`. `NO_RAW_HIT_E1` is permissible only after all three database searches have actually completed for that structure.

Step182 does not generate Stage10 final ranks, Stage11 MSA, Stage12 conserved residues, Stage13 functional confirmation, or wet-lab evidence. Separate downstream logic will consume the full family evidence table.
