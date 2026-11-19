# Stage 10 preliminary evidence triage (structure-optional)

`10_partial_evidence_triage_portable.py` attaches **verified** Stage 09 evidence to the **full** Stage 180 independent-noncandidate-genomic-context review, preserving all protein families.

It produces `all_family_stage10_review.tsv`, `dark_family_stage10_review.tsv`, `context_supported_stage10_shortlist.tsv`, and SHA-256 source/output provenance. Source results are not modified; `check` is read-only; `run` will not overwrite an existing output directory.

```bash
python scripts/10_partial_evidence_triage_portable.py \
  --context-dir RESULTS/context_screen/10_context_screen \
  --stage09-dir RESULTS/STAGE09/merged \
  --prepared-dir RESULTS/STAGE09/prepare \
  --out-dir RESULTS/stage10_review \
  --mode check
# Use --mode run only after check succeeds.
```

## Output semantics

- **T1**: Dark family with independent noncandidate replication/nucleotide-metabolism context **and** a convincing structural connection by configured structural criteria. Still not a confirmed enzyme activity.
- **T2**: Dark family with independently supported context, without convincing structure. Includes both weak structural matches and unsearched families; the latter are NOT treated as negative.
- **T3**: Dark family with convincing structural connection but no independently supported context. May warrant review; structural homology alone is not biochemical activity.
- **T4**: Other dark/unresolved families, retained for future evidence.
- **T5**: Known-evidence reference families, preserved outside the dark-priority selection.

Review order is deterministic: tier, independently supported samples, number of members supporting best context category, noncandidate context sample span, overall sample span, then structure ID. It is **not** the historical I2 score, a probability, an enzymatic activity predictor, nor final Stage 13 selection. In particular, `WEAK_RAW_HIT_ONLY` is not promoted to a convincing functional connection, and `NO_RAW_HIT_E1` (searched/no hit) is distinct from `NOT_EVALUATED_NO_STRUCTURE` (not searched).

Conclusions are contingent on viral gene calling, dark-family selection criteria, profile searches, context-window configuration, and quality of structural predictions. The script does not itself run Foldseek, GPU predictions, MSA, or Stage 11–13 analysis. It is a release-facing *intermediate* Stage 10 module for partially observed cohorts, not proof that the full pipeline is completed.
