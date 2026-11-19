# Stage 186: structure-supported residue evidence in the Stage 13 catalog

This utility joins a **fresh Stage 185 catalog generated from the latest Stage 10 review** with the verified structure-supported residue analysis from `12_msa_structure_conservation_portable.py`.

It does **not** label catalytic residues or confirm enzymatic activities. The Stage 12 `conserved_windows.tsv` records *sequence windows of conserved structure-supported residues*, not 3D spatial clusters; hence `stage12_spatial_state` remains `NOT_EVALUATED`.

## Inputs

- `--catalog-dir`: new `Stage 185` output generated from the **same updated Stage 10 review** as the new Stage 184 prepare/collect outputs. `stage185_provisional_provenance.json` required.
- `--stage12-dir`: validated output of `12_msa_structure_conservation_portable.py`, including `conservation_provenance.json`, `candidate_conservation_summary.tsv`, `conserved_windows.tsv`, and residue tables.
- `--out-dir`: NEW directory; tool refuses to overwrite.

## Execution

```
python -B scripts/13_integrate_structure_residue_evidence_portable.py \
    --catalog-dir "$WORK/results/stage11_13_structure_refresh/catalog" \
    --stage12-dir "$WORK/results/stage12_MVF000472_structure" \
    --out-dir "$WORK/results/stage13_structure_integrated" --mode check

python -B scripts/13_integrate_structure_residue_evidence_portable.py \
    --catalog-dir "$WORK/results/stage11_13_structure_refresh/catalog" \
    --stage12-dir "$WORK/results/stage12_MVF000472_structure" \
    --out-dir "$WORK/results/stage13_structure_integrated" --mode run
```

## Outputs

- `all_family_stage13_structure_integrated.tsv`: Stage 185 full Family catalog plus Stage 12 structure-supported residue count columns. All Family rows are preserved.
- `structure_evaluated_stage13_review.tsv`: evaluated subset only.
- `stage186_structure_integration_provenance.json`: source and output SHA-256 checksums; explicit `spatial_3d_cluster_analysis_completed = false` and `final_biochemical_selection_completed = false`.

## Protection and validation

The program verifies Stage 185 source/output SHA checksums, Stage 12 summary and per-Family provenance, residue-table hashes, cross-stage A3M sequence content hashes, and matching Family IDs. It refuses to overwrite files. For a new experiment, repeat the earlier stages; do not reuse unmatched A3Ms or provenance across cohorts.

Synthetic test:

```
python -B tests/test_stage186_structure_integration.py
```
