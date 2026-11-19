# MetaViro-DRANK

**Marine Viral Dark-Protein Discovery and Functional Prioritization**

MetaViro-DRANK is a research pipeline for family-level analysis and functional prioritization of **unresolved proteins from user-provided viral metagenomic datasets**. It integrates clustering, functional annotation, genomic context, and optional structural evidence to identify candidates for investigating **viral replication and nucleic-acid metabolism**. Outputs include family-level evidence tables, genome-neighborhood summaries, and hypothesis-generating candidate lists.

## 1. Overview and workflow

MetaViro-DRANK is a **downstream** workflow. It requires precomputed protein FASTA sequences, compatible gene coordinates, and eggNOG-mapper annotations, or a compatible completed family/context analysis for resumption. It does **not** assemble FASTQ/contigs, identify viral sequences from raw assemblies, or provide experimentally validated nuclease or clinical/diagnostic classification.

The analysis proceeds through these stages:

```
Input: protein FASTA + Prodigal-style GFF + eggNOG-mapper annotations
  │
  ├─ 00–02  Input validation → unresolved proteins → viral replication-context screen
  ├─ 03–04  MMseqs2 clustering → recurrent unresolved protein families
  ├─ 05–07  Swiss-Prot / PHROG / Pfam evidence → darkness audit → family-wide context
  ├─ 08     Representative structure queries / optional external structure prediction
  ├─ 09–10 Optional Foldseek comparisons → confidence-aware prioritization
  └─ 11–13 Optional/advanced family MSA, conservation and integrated evidence
                         ↓
           Family catalogs, evidence status and candidate shortlists

```

**Interpretation:** An *unresolved* protein is not necessarily novel. MMseqs2 clustering uses **≥30% sequence identity and ≥70% coverage**; a recurrent family here has at least two members **within the input cohort**, not necessarily in different samples. Functional/domain matches and genomic proximity to replication or nucleotide-metabolism genes provide complementary evidence. Predicted structure matches, conserved residues, and charged surface clusters support **hypotheses**, not confirmed catalytic activity.

For decision rules and input requirements, see [Methods](docs/METHODS.md) and [Input format](docs/INPUT_FORMAT.md).

## 2. System requirements and software versions

**Platform:** 64-bit Linux (Ubuntu recommended) or Ubuntu under **WSL2**; Conda or Mamba; Bash and Git. Other operating systems have not been established as supported deployment targets.

**Hardware planning guidance (not validated minimum specifications):** for CPU-first analyses, plan for approximately **8 CPU threads and 32 GB RAM**, with disk space for inputs, temporary files and reference databases. Requirements vary with dataset size. **GPU acceleration and a separately installed ColabFold environment are optional** for new structure predictions; CPU-only context analysis does not require a GPU.

| **Software**                    | **Role**                                       | **Version information in this repository**                                                                   |
| ------------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| Python                          | Pipeline execution                             | **3.11** specified in `environment.yml`                                                                      |
| PyYAML, Biopython               | Configuration, protein/structure handling      | PyYAML unpinned; Biopython **1.88** pinned                                                                   |
| MMseqs2                         | Protein clustering; Swiss-Prot/PHROG searches  | Conda dependency; **version not pinned**                                                                     |
| HMMER (`hmmscan`, `hmmpress`)   | Pfam domain analysis                           | Conda dependency; **version not pinned**                                                                     |
| MAFFT                           | Optional family multiple-sequence alignment    | Conda dependency; **version not pinned**                                                                     |
| Foldseek                        | Optional protein structure similarity searches | Conda dependency; **version not pinned**                                                                     |
| eggNOG-mapper                   | **Upstream** protein annotation                | External input producer; **no version pinned here**                                                          |
| Prodigal-compatible gene caller | **Upstream** proteins/GFF                      | External input producer; **no version pinned here**                                                          |
| ColabFold                       | Optional GPU structure prediction              | **1.6.1 recorded in a historical prediction script**; separate installation, not pinned by `environment.yml` |

Unpinned does **not** mean a specific version was validated. Record actual executable and database versions before publishing an analysis; see [Software versions and reproducibility](docs/SOFTWARE_VERSIONS.md).

## 3. Installation and analysis

### 3.1. Install the main environment

```
git clone https://github.com/Min-Jeong214/MetaViro-DRANK.git
cd MetaViro-DRANK
conda env create -f environment.yml
conda activate metaviro-drank
python -B run_metaviro.py --help
python -B run_pipeline.py --help
```

Public availability does not change the [permitted use conditions](#5-copyright-and-permitted-use). ColabFold and its GPU dependencies are **not** installed by the main environment file.

### 3.2. Download or prepare reference databases

The following resources are **not bundled** with this code. Stage 05–07 searches require: (1) a Swiss-Prot **MMseqs2 DB** including its `.dbtype`, (2) a complete PHROG **MMseqs2 profile DB** with index/header files, (3) a pressed Pfam HMM including `.h3f/.h3i/.h3m/.h3p`, and (4) a PHROG annotation table with `phrog`, `annot`, `category` columns. An optional advanced Foldseek search additionally requires PDB, AF_SwissProt and BFVD Foldseek DBs.

```
# Downloads Swiss-Prot FASTA; it is NOT yet an MMseqs2 search DB.
bash scripts/setup_swissprot.sh
mmseqs createdb databases/swissprot/swissprot.faa \
  databases/swissprot/swissprot_mmseqs

# Downloads Pfam HMM and creates HMMER pressed indexes.
bash scripts/setup_pfam.sh
```

**Important:** PHROG profiles/annotations and the three Foldseek databases require separate preparation. The repository does **not** provide a one-command installation for these resources. Check formats, sources, and usage terms before analysis. See [Database setup](docs/DATABASES.md) and [Installation details](docs/INSTALLATION.md).

### 3.3. Prepare input data and configuration

```
input/my_dataset/
  sample_A/proteins.faa
  sample_A/genes.gff
  sample_A/annotations.emapper.annotations
  sample_B/proteins.faa
  sample_B/genes.gff
  sample_B/annotations.emapper.annotations

```

Input files must have mutually consistent protein IDs and **Prodigal-style CDS identifiers**; arbitrary GFF3 is not supported. The `annotations.emapper.annotations` file must retain the `#query` column header. Copy [config/example.yaml](config/example.yaml), then set your own input, output and DB paths:

```
cp config/example.yaml config.local.yaml
# Edit config.local.yaml: input_dir, output_dir and database prefixes
```

The configuration is local and ignored by Git. Use **new, non-overlapping output directories** for each run.

### 3.4. Recommended CPU-first analysis

```
# Validate paths, input consistency and required references first.
python -B run_metaviro.py run \
  --config config.local.yaml \
  --screen-out-dir results/my_context_review \
  --normalize-ids --mode check

# Execute Stage 00–07, prepare Stage 08 representatives,
# and perform the preliminary context-first review.
python -B run_metaviro.py run \
  --config config.local.yaml \
  --screen-out-dir results/my_context_review \
  --normalize-ids --mode run
```

`--normalize-ids` makes a consistent *copy* of sample IDs and related annotations; omit it when identifiers are already unique across all samples. **This command does not run ColabFold, Foldseek, or final Stage 13 biochemical target selection.**

To resume a compatible completed upstream workflow:

```
python -B run_metaviro.py resume \
  --workflow-output /path/to/completed/workflow \
  --screen-out-dir results/resumed_review \
  --mode check
python -B run_metaviro.py resume \
  --workflow-output /path/to/completed/workflow \
  --screen-out-dir results/resumed_review \
  --mode run
```

**Advanced structural analysis (separate execution paths):** `run_pipeline.py` prepares or collects Stage 08 predictions; `run_pipeline_full.py` and dedicated Stage 09–13 modules support structure search, alignment, residue analysis and evidence integration, subject to the required inputs and reference resources. These steps are **not** automatically completed by the CPU-first command. See the [CLI manual](docs/CLI_USER_WORKFLOW.md), [Advanced usage](docs/USAGE.md), and [Structure/sequence evidence](docs/STAGE11_12_PARTIAL_SEQUENCE_WORKFLOW.md).

### 3.5. Fresh-clone verification (recommended before real-data analysis)

From a newly cloned GitHub checkout, first run the built-in synthetic verification:

```
conda activate metaviro-drank
python -B scripts/validate_fresh_clone.py --mode smoke
```

This creates **independent synthetic datasets with one and three samples**, executes Stage 00–02 (not just preflight), checks output counts and provenance, and confirms that a deliberately incorrect GFF identifier is rejected. Temporary fixtures are removed afterward; research datasets are not bundled. Stage 00–02 requires Python and PyYAML but does **not** require MMseqs2/HMMER.

Then verify the actual database paths and CLI for your own dataset:

```
python -B scripts/validate_fresh_clone.py --mode database \
  --config config.local.yaml \
  --screen-out-dir results/validation_context \
  --normalize-ids
```

The command above is a **non-writing database/Stage 00–07 preflight**. To execute Stage 00–07 and the CPU context review with the configured input and references, append `--execute` and use unused output directories. Database-backed execution may take substantial time and storage. For automated checks, see [Validation and tests](docs/VALIDATION_STATUS.md).

## 4. Outputs and interpretation

| **Stage/group**  | **Main outputs**                                                                                              | **Interpretation**                                                               |
| ---------------- | ------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| 00–02            | Input/provenance logs, `01_unresolved/*.tsv`, `02_context/*.tsv`                                              | Annotation-unresolved ORFs and nearby replication-related gene context           |
| 03–04            | `03_clustering/viral_families_cluster.tsv`, `04_recurrent/recurrent_families_*.tsv`                           | Protein families and within-cohort recurrence                                    |
| 05–07            | `05_annotation/member_level_evidence.tsv`, `06_darkness/family_darkness_audit.tsv`, `07_family_context/*.tsv` | Functional/domain annotations, family darkness and genomic-neighborhood evidence |
| 08               | `08_structure/.../prepared/` and, if available, `collected/structure_results.tsv`                             | Structure queries and verified structure collection status                       |
| 09–10 (optional) | `09_10/09_external/`, `09_10/10_ranking/confidence_aware_ranking.tsv`, `09_10/10_shortlists/`                 | External structure matches, ranked candidate evidence                            |
| 11–13 (advanced) | MSA, conservation and integrated evidence tables                                                              | Optional evidence; evaluation status depends on compatible inputs and execution  |

The **official CPU-first entrypoint** also writes `run_metaviro_provenance.json`, `review/10_context_screen/all_family_context_screen.tsv`, and `review/10_context_screen/dark_replication_context_review.tsv` under `--screen-out-dir`. Stage-specific paths depend on the entrypoint and run mode. See [Complete output guide](docs/OUTPUTS.md).

**Scientific boundary:** a candidate shortlist expresses a testable functional hypothesis, **not confirmation of nuclease activity, enzyme specificity, or novelty**. Missing predicted structures or database matches must remain *not evaluated / pending* rather than negative functional evidence.

## 5. Copyright and permitted use

Copyright © 2026 Korea Institute of Ocean Science and Technology (KIOST).

MetaViro-DRANK is registered with the Korea Copyright Commission as a computer software work.

- **Registration No.:** C-2026-060878
- **Creation and first publication:** November 9, 2026
- **Registration date:** November 18, 2026

Non-commercial research and educational use is permitted. Commercial use, modification, and redistribution require prior written permission from KIOST.

**All rights reserved.**

Third-party software and database rights remain with their respective owners. See [Data and release policy](docs/DATA_AND_RELEASE_POLICY.md).

## 6. Manuscript, data availability and validation

### Associated scientific manuscript

**Title:** *An Evidence-Integrated Framework for Prioritizing Recurrent Viral Protein Families of Unknown Function in Marine Metaviromes*

**Authors:** Min-Jeong Kim, Yu Jin Kim, Hyun Ju Ha, Joon Sang Park, and Taek-Kyun Lee.

**Journal:** Marine Life Science & Technology

**Status:** Submitted; under editorial consideration.

The manuscript has not yet been published. Bibliographic details and the article DOI will be added when available.

### Scientific data availability

Study-specific sequence inputs and computational validation evidence are available in the [Zenodo dataset (DOI: 10.5281/zenodo.24184202)](https://doi.org/10.5281/zenodo.24184202), separately from the source-code repository.

### Software verification

Synthetic fresh-clone tests, regression tests, and GitHub Actions checks assess software behavior, **not experimentally validated biological functions**. Test coverage, database-dependent procedures, and pending evaluations are described in [Validation and tests](docs/VALIDATION_STATUS.md).

For reproducibility, record the Git commit, software and database versions, configuration parameters, and input checksums.

## 7. Documentation

- **Getting started:** [Installation](docs/INSTALLATION.md), [Database preparation](docs/DATABASES.md), [Input schema](docs/INPUT_FORMAT.md), [Command-line workflow](docs/CLI_USER_WORKFLOW.md)
- **Methods and outputs:** [Methods](docs/METHODS.md), [Output guide](docs/OUTPUTS.md), [Software versions](docs/SOFTWARE_VERSIONS.md)
- **Advanced analyses:** [Advanced usage](docs/USAGE.md), [Stage 11–12 sequence/structure evidence](docs/STAGE11_12_PARTIAL_SEQUENCE_WORKFLOW.md), [Full-cohort GPU prediction](docs/FULL_COHORT_GPU_PREDICTION.md)
- **Verification and release:** [Validation status](docs/VALIDATION_STATUS.md), [Data and release policy](docs/DATA_AND_RELEASE_POLICY.md)

Optional full-cohort MSA preparation, restartable single-GPU ColabFold prediction, integrity checks, and retries require a **separate ColabFold environment**; research dataset outputs are not bundled.
