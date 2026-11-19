# Installation and preparation — Linux / Ubuntu / WSL2

> MetaViro-DRANK requires upstream protein calling and annotation. The main Conda environment does **not** install ColabFold, all biological reference databases, or upstream eggNOG-mapper databases automatically.

## 1. System preparation

Recommended environment: 64-bit Linux / Ubuntu (or WSL2 Ubuntu), Conda/Mamba, Bash, Git, curl, gzip and a working network connection for acquiring reference databases. Hardware requirements depend on input size and database versions; see the main [README](../README.md).

```bash
# Clone requires GitHub permission while the repository remains private.
git clone https://github.com/Min-Jeong214/MetaViro-DRANK.git
cd MetaViro-DRANK
conda env create -f environment.yml
conda activate metaviro-drank

python --version
python -c 'import yaml, Bio; print("Python dependencies ready")'
mmseqs version
hmmscan -h | head -n 4
mafft --version
foldseek version
python -B run_metaviro.py --help
```

`environment.yml` pins Python 3.11 but **does not pin all tool versions**. Record installed versions using [Software versions](SOFTWARE_VERSIONS.md); this environment has not been accepted as a clean install on every supported Linux configuration.

## 2. Input prerequisite (performed outside this pipeline)

For **every** input sample create `proteins.faa`, `genes.gff`, `annotations.emapper.annotations` in a named sample folder. The GFF must follow the Prodigal-style CDS ID mapping required by the current parser; protein and eggNOG `query` IDs must match. Check [Input schema](INPUT_FORMAT.md) before attempting database installation.

## 3. Download and prepare databases

### 3.1 Swiss-Prot (sequence search)

```bash
bash scripts/setup_swissprot.sh
mmseqs createdb databases/swissprot/swissprot.faa \
  databases/swissprot/swissprot_mmseqs
```

**Do not** point the YAML `swissprot` key at the downloaded raw FASTA. The correct prefix after these commands is `<repo>/databases/swissprot/swissprot_mmseqs`.

### 3.2 Pfam (domain search)

```bash
bash scripts/setup_pfam.sh
# Expected: databases/pfam/Pfam-A.hmm and .h3f/.h3i/.h3m/.h3p
```

### 3.3 PHROG (protein-profile search)

**Manual preparation is currently required.** Provide a validated **MMseqs2 profile** DB and a tab-delimited PHROG annotation mapping with headers `phrog`, `annot`, `category` from a properly sourced PHROG release. A raw PHROG FASTA file or an arbitrary HMM file is **not** interchangeable with the profile DB expected by `scripts/05_annotation_saturation_portable.sh`. The repository does not yet provide an end-to-end downloader/converter for this DB.

### 3.4 Advanced Foldseek databases (optional)

For Stage 09–10 structure comparisons, configure a Foldseek DB root containing the exact prefixes `PDB`, `AF_SwissProt`, `BFVD`. These large databases are **not bundled** and require release-specific preparation; there is currently no all-in-one verified download command for them. For sources/formats/checks read [Database documentation](DATABASES.md).

## 4. Point the local configuration at the installed data

```bash
cp config/example.yaml config.local.yaml
# Edit local paths in config.local.yaml.
# 'swissprot': .../databases/swissprot/swissprot_mmseqs
# 'pfam': .../databases/pfam/Pfam-A.hmm
# 'phrog' and 'phrog_annotations': your prepared external paths
python -B run_pipeline.py --config config.local.yaml --through 08 --check
```

`config.local.yaml` is ignored by Git. If only the sample parser needs checking, use `--through 02 --check` without search databases.

## 5. Optional structure prediction environment

ColabFold is **separately installed and versioned**, not part of `environment.yml`. Prepare a supported CUDA/GPU and compatible ColabFold runtime according to the upstream installation instructions. The historical `scripts/08_structure_prediction.sh` notes version 1.6.1 for reproduction of its settings; this is **not** a newly verified clean-install specification. Alternatively generate structures elsewhere and use Stage 08 collection with their verified mapping:

```bash
python -B run_pipeline.py --config config.local.yaml --through 08 \
  --structure-mode collect \
  --structure-prediction-dir /absolute/path/to/matching_pdbs
```

Stage 08 `prepare` needs no GPU but creates **queries only**, not PDB predictions. See [Advanced usage](USAGE.md).

## 6. Validate a fresh GitHub checkout

```bash
# Synthetic datasets of one and three samples; genuine Stage 00–02 execution.
python -B scripts/validate_fresh_clone.py --mode smoke

# Optional independent regression test (runs the same smoke in unittest).
python -m unittest discover -s tests -p 'test_release_smoke.py' -v
```

This synthetic verification requires PyYAML and Python, and works without MMseqs2, HMMER, real sample files or database downloads. Stage 00–02 provenance marks unused programs as `NOT_REQUIRED_FOR_THIS_STAGE` (or `NOT_REQUIRED_FOR_STAGE_02`). A Git checkout records the commit; a source archive explicitly records unavailable Git metadata.

## 7. Check installed reference databases and run the analysis

```bash
cp config/example.yaml config.local.yaml
# Configure absolute input, output and reference database paths.
python -B scripts/validate_fresh_clone.py --mode database \
  --config config.local.yaml \
  --screen-out-dir results/validation_context \
  --normalize-ids

# After successful preflight, run with new, unused directories:
python -B scripts/validate_fresh_clone.py --mode database \
  --config config.local.yaml \
  --screen-out-dir results/validation_context \
  --normalize-ids --execute
```

Stage 03 onward requires MMseqs2, while Stage 05 onward requires HMMER and complete Swiss-Prot, PHROG and Pfam files. Stage 08–13 and GPU structure generation have separate runtime and reference dependencies described in [Advanced usage](USAGE.md). Database-backed verification is not part of lightweight GitHub Actions because it requires separately licensed/large reference data.
