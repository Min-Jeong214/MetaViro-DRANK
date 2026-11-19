# Database setup and provenance

The repository does **not** distribute biological reference databases, pretrained model weights or the underlying research datasets. Download or prepare resources from the original maintainers under their terms, and record their exact releases/checksums and dates.

## Main sequence/domain annotation databases (Stage 05–07)

| YAML field | Required object | How to prepare / known limitation |
|---|---|---|
| `swissprot` | Existing MMseqs2 protein search DB prefix | `bash scripts/setup_swissprot.sh` downloads UniProtKB/Swiss-Prot FASTA, then use `mmseqs createdb` as below |
| `phrog` | Existing **PHROG profile** MMseqs DB prefix with `.dbtype` | **Manual PHROG profile preparation required**; do not substitute simple FASTA or HMM |
| `pfam` | Pfam-A HMM with `hmmpress` indexes (`.h3f/.h3i/.h3m/.h3p`) | `bash scripts/setup_pfam.sh` |
| `phrog_annotations` | PHROG TSV with `phrog`, `annot`, `category` headers | Prepare from the same PHROG release as the profile DB; must be ID-compatible |

### Swiss-Prot

```bash
bash scripts/setup_swissprot.sh
mmseqs createdb databases/swissprot/swissprot.faa \
  databases/swissprot/swissprot_mmseqs
# YAML: swissprot: /absolute/repository/path/databases/swissprot/swissprot_mmseqs
```

`setup_swissprot.sh` downloads the **current** UniProtKB/Swiss-Prot release, which is not a fixed historical version. It **does not** create the MMseqs2 DB by itself. Record download date, URL and accession-release details for reproducibility.

### Pfam

```bash
bash scripts/setup_pfam.sh
# YAML: pfam: /absolute/repository/path/databases/pfam/Pfam-A.hmm
```

The script downloads `Pfam-A.hmm.gz` from the current EBI release location and presses it with HMMER. Record the exact Pfam release used.

### PHROG

The Stage 05 script runs `mmseqs search` **with the PHROG profile DB as its first input**, followed by `mmseqs createtsv`. As a result, an arbitrary FASTA-derived MMseqs DB is **not a drop-in replacement** for the validated PHROG profile resource. The software currently requires users to provide the profile DB and its annotation mapping; a release-independent automatic preparation recipe has **not** been verified. Check the profile format before starting a long run.

### Command-line validation

```bash
python -B run_pipeline.py --config config.local.yaml --through 08 --check
```

This checks path/format preconditions, but it **cannot** establish biological completeness, correct PHROG release mapping, or a valid end-to-end scientific result.

## Optional Foldseek structure-reference DBs (Stage 09–10)

```text
/absolute/path/to/foldseek_db_root/
  PDB             (+ Foldseek DB index and related files)
  AF_SwissProt    (+ Foldseek DB index and related files)
  BFVD            (+ Foldseek DB index and related files)
```

Run `python -B run_pipeline_full.py --help` for `--db-root` or compatible `--reuse-search-root` usage. `--db-root` executes new searches (if databases and structures are provided); reuse is appropriate **only** when structural inputs, IDs, search settings and provenance agree. Stages 09/10 employ PDB and AF_SwissProt e-value `1.0`, `--max-seqs 1000`; BFVD `--max-seqs 500` (documented current profile). See [Methods](METHODS.md).

**Interpretation:** a Foldseek hit is **not** direct evidence of enzymatic function, and no hit is **not** proof of functional inactivity. Preserve source release and structure mapping information for any scientific report.

## Provenance checklist

- Reference resource name, version/release and provider URL
- Download date, local path and checksum when available
- Format conversion and indexing commands plus software versions
- License/reuse terms of source sequences, models and structures

Do **not** commit DB files, personal data or large generated outputs into GitHub.
