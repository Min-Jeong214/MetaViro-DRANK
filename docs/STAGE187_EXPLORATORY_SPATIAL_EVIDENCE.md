# Stage 187 — Exploratory 3D spatial evidence integration

This portable step joins exploratory spatial clustering from `12_charged_spatial_clusters_exploratory.py` to the structure-enriched Stage 186 catalog. It is **not** a historical reconstruction, catalytic-site prediction, confirmed nuclease activity, or final biochemical ranking.

## Inputs

- Stage 186: `all_family_stage13_structure_integrated.tsv` and `stage186_structure_integration_provenance.json`
- Exploratory Stage 12: `family_spatial_cluster_summary.tsv`, `charged_spatial_clusters.tsv`, and `spatial_cluster_provenance.json`

Stage 187 verifies file hashes, lineage to the same original Stage 12 conservation provenance, Family identity, PDB/residue-input fingerprints, cluster counts and classification consistency **before writing any output**.

## Command

```bash
python -B scripts/13_integrate_exploratory_spatial_evidence_portable.py \
  --stage186-dir /path/to/stage13_structure_integrated \
  --spatial-dir /path/to/stage12_spatial \
  --out-dir /new/path/to/stage13_spatial_integrated \
  --mode check

# Only if check passed:
python -B scripts/13_integrate_exploratory_spatial_evidence_portable.py \
  --stage186-dir /path/to/stage13_structure_integrated \
  --spatial-dir /path/to/stage12_spatial \
  --out-dir /new/path/to/stage13_spatial_integrated \
  --mode run
```

Never provide an existing output directory. `check` writes no files.

## Results and interpretation

- `all_family_stage13_spatial_integrated.tsv`: Every original Family, plus spatial status and parameter metadata.
- `spatial_evaluated_stage13_review.tsv`: Only Families actually evaluated in the exploratory 3D analysis.
- `spatial_cluster_details.tsv`: Preserves full cluster details; header-only if zero clusters detected.
- `stage187_spatial_integration_provenance.json`: Input/output hashes and counts.

`EXPLORATORY_NO_CLUSTER_DETECTED` means the Family **was assessed** but no qualifying cluster was found under the stated exploratory criteria. `NOT_EVALUATED` means no spatial analysis was performed. Do not set unevaluated Family counts to zero. The procedure uses C-alpha distance-connected components, not strict all-pairs cliques. Zero clusters under one distance/minimum-size choice does not rule out enzyme function. It does not experimentally confirm function or produce a final I2 ranking.
