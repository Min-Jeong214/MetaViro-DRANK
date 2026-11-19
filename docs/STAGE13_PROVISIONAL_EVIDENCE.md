# Stage 13: provisional, structure-optional evidence catalog

A standalone, cohort-agnostic bridge for MetaViro-DRANK **Stage183 (all-family context + verified Stage09 structural status)** and **Stage184 (subset MAFFT and sequence-only conservation)**. It covers *all* families, independent of the number of input samples, regardless of whether a protein structure exists. No Pilot family lists, laboratory-validated nuclease assumptions, or historical Top5/Top3 IDs are used.

**Important:** This is **NOT** the complete Stage13 `final_candidate_selection` and **NOT** full structure-based Stage12. It is deliberately named a provisional catalog. It does not run GPU structure prediction, spatial residue clustering, Foldseek, catalytic-site identification, an assay, historical I2 ranking, or a confirmed biochemical target nomination. Existing Stage13 full evaluation scripts remain separate and are not overwritten.

## Usage

```bash
python scripts/13_partial_evidence_catalog_portable.py \
  --review-dir "$STAGE183_REVIEW" \
  --conservation-dir "$STAGE184_CONSERVATION" \
  --out-dir "$NEW_DIR" --limit 30 --mode check

python scripts/13_partial_evidence_catalog_portable.py \
  --review-dir "$STAGE183_REVIEW" \
  --conservation-dir "$STAGE184_CONSERVATION" \
  --out-dir "$NEW_DIR" --limit 30 --mode run
```

`--mode check` does not create output files. `--mode run` refuses an existing directory; never delete prior results merely to repeat analysis. Checks Stage183/184 SHA-256, the Stage184 prepare source and output SHA-256 chain, exact family IDs and sizes, MSA coverage, per-residue coordinate completeness and MSA summary counts. Selected family MSA status is distinct from all other families. Files are re-checked before writing.

## Outputs

* `all_family_stage13_evidence.tsv`: every Family in the verified Stage183 review order; Stage09 state, Stage11/12 MSA status, and explicit `NOT_EVALUATED` spatial status.
* `context_supported_stage13_review.tsv`: T1/T2 Family rows with supported contextual evidence (not an activity finding).
* `msa_sequence_stage13_evidence.tsv`: all selected sequence-only MSA Family rows.
* `provisional_stage13_review_queue.tsv`: up to `--limit` A/B/C dark Family rows from the existing Stage183 review order. This is *not* a functional score or final rank.
* `stage185_provisional_provenance.json`: source/output SHA-256 plus explicit completion boundaries.

Sequence-conserved positions indicate exploratory homologous conservation, not catalytic residues. Strong sequence conservation, structural similarity and genomic neighborhood are all hypothesis-generating, not experimental validation.

## Not yet integrated / future GitHub release gate

The new program is a **separate CLI module**, not yet in `run_metaviro.py advanced`. To complete full Stage11-13, integrate optional structure prediction from family MSA, Stage12 PDB/pLDDT and spatial analysis for verified matching sequences, and confidence-aware Stage13 evidence integration, then add proper resume modes and full new-cohort tests. The GitHub release should report exactly which functional modules have and have not been exercised.
