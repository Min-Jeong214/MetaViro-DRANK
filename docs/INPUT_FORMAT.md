# New dataset input format

## Required layout

```text
input_dir/
  sample_A/
    proteins.faa
    genes.gff
    annotations.emapper.annotations
  sample_B/
    proteins.faa
    genes.gff
    annotations.emapper.annotations
```

`input_dir` points to a directory containing sample subdirectories, **not** directly to FASTA files. Sample folder names use ASCII letters, digits, `_`, `.` and `-` and are used in output tables.

- `proteins.faa`: amino-acid FASTA. Sequence IDs (first header token) must be **globally unique across all samples**. All sequences must be nonempty. Terminal `*` is accepted.
- `genes.gff`: 9-column GFF, `CDS` feature, positive 1-based coordinates, `ID=...` attributes. **Current code assumes Prodigal-style ID with a trailing `_N`; mapped protein ID = `<contig_column_1>_<N>`.** Every FASTA protein must match one CDS under this convention. It is *not* a generic GFF3 protein-ID parser.
- `annotations.emapper.annotations`: eggNOG-mapper tab-delimited output with a literal `#query` header. Mandatory columns: `query`, `eggNOG_OGs`, `Description`, `Preferred_name`; additional `KEGG_ko`, `PFAMs` are used by unresolved classification when present. Each reported ID must be in FASTA. Proteins without hits may be absent from data rows.

## Biological filtering used today

An unresolved protein has missing `Description`, `Preferred_name`, `KEGG_ko`, `PFAMs` and `PFAM` values (or no row); candidate `ge_min_aa` also requires length >= `min_aa` (default in provided example: 150 aa). Genomic replication/nucleotide-metabolism anchors are recognized from functional text **only for eggNOG records whose `eggNOG_OGs` contains `10239|Viruses`**; candidates are unresolved ORFs within `window` ORF indices of a matching anchor on the same contig. Family recurrence threshold in this runner: >=2 family members, based on MMseqs 30% identity / 70% coverage clustering.

**Important:** this is a contextual enrichment heuristic. It is not direct enzyme-function prediction, nor is it a viral-origin validator; avoid claims that all selected proteins are bona-fide viral or biochemical nucleases without independent evidence.

## Input validation

```bash
python run_pipeline.py --config config.local.yaml --through 02 --check
```

For full Stage 08 DB preflight:

```bash
python run_pipeline.py --config config.local.yaml --through 08 --check
```

No Stage 00 input check generates annotation, GFF or protein IDs; prepare those beforehand using an upstream ORF caller/annotation workflow. See generated illustrative files from `examples/generate_demo_input.py`.
