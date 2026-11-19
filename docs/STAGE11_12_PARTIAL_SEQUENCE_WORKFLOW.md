# MetaViro-DRANK Stage 11~12: optional MSA and sequence-only conservation

This is a **preliminary sequence-only analysis** on a user-configurable subset of Stage 183 families. It does not claim that structure prediction, spatial clustering, Foldseek reanalysis, the original Stage 12 structural conservation, biochemical activity validation, or full Stage 13 final candidate selection have been completed.

Inputs must come from the **same cohort**: Stage 183 `review/all_family_stage10_review.tsv` and manifest; Stage 08 `structure_id_map.tsv`; Stage 05 `recurrent_members.tsv` and `recurrent_members.faa`. All member/family IDs, FASTA query sequences, family sizes, and the SHA-256 of representative sequences are checked; Stage 183 provenance is checked. The generic CLI has no Pilot path or constant number of families.

## Typical command sequence

```bash
# first verify inputs read-only
python scripts/11_12_partial_msa_conservation_portable.py prepare \
    --review-dir "$STAGE183_REVIEW" \
    --structure-map "$STRUCTURE_MAP" \
    --members "$STAGE05/recurrent_members.tsv" \
    --fasta "$STAGE05/recurrent_members.faa" \
    --limit 10 --out-dir "$RUN/prepare" --mode check

# materialize top review subset (not a biochemical function ranking)
python scripts/11_12_partial_msa_conservation_portable.py prepare \
    --review-dir "$STAGE183_REVIEW" \
    --structure-map "$STRUCTURE_MAP" \
    --members "$STAGE05/recurrent_members.tsv" \
    --fasta "$STAGE05/recurrent_members.faa" \
    --limit 10 --out-dir "$RUN/prepare" --mode run

# existing validated Stage11 family MSA script, using only the selected subset
python scripts/11_family_msa_portable.py \
    --members "$RUN/prepare/selected_members.tsv" \
    --fasta "$RUN/prepare/selected_proteins.faa" \
    --structure-map "$RUN/prepare/selected_structure_map.tsv" \
    --out-dir "$RUN/msa" --threads 4

# Stage12 sequence-only signal, no PDB or structure prediction required
python scripts/11_12_partial_msa_conservation_portable.py collect \
    --prepare-dir "$RUN/prepare" --msa-dir "$RUN/msa" \
    --min-homologs 2 --identity 0.8 \
    --out-dir "$RUN/conservation" --mode check
python scripts/11_12_partial_msa_conservation_portable.py collect \
    --prepare-dir "$RUN/prepare" --msa-dir "$RUN/msa" \
    --min-homologs 2 --identity 0.8 \
    --out-dir "$RUN/conservation" --mode run
```

The alignment script processes exact-duplicate collapsed homologs and requires MAFFT. A conserved position counts only when at least `--min-homologs` **non-gap homologous sequences other than the query** cover that position, and the fraction matching the query residue is at least `--identity`. It is an exploratory conservation signal, not evidence of a catalytic residue. Small or redundant families may yield `INSUFFICIENT_HOMOLOG_SUPPORT` or `SINGLE_UNIQUE_SEQUENCE` rather than being marked nonfunctional. All unselected families are retained in `all_family_stage11_12_status.tsv` as `NOT_SELECTED_FOR_MSA`; no structure confidence/3D catalytic spatial assertions are generated.

Output: `all_family_stage11_12_status.tsv`, `selected_family_msa_conservation.tsv`, `selected_query_residue_conservation.tsv` and SHA-256 manifest. All output directories must be new. Completed MSA files are never deleted or reprocessed on failure; run `collect` independently after resolving a downstream issue.
