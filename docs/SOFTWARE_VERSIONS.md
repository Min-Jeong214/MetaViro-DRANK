# Software versions and reproducibility

MetaViro-DRANK does **not** currently ship a fully locked or cross-platform validated environment. The difference between a declared dependency, a historical analysis version, and a version installed on your machine is important.

| Dependency | Installation source | Declared/pinned version | Record at runtime |
|---|---|---|---|
| Python | `environment.yml` | `3.11` | `python --version` |
| PyYAML | Conda | Unpinned | `python -c 'import yaml; print(yaml.__version__)'` |
| Biopython | Conda | Unpinned | `python -c 'import Bio; print(Bio.__version__)'` |
| MMseqs2 | Bioconda/Conda; required from Stage 03 | Unpinned | `mmseqs version` |
| HMMER | Bioconda/Conda; required from Stage 05 | Unpinned | `hmmscan -h` (HMMER header) |
| MAFFT | Bioconda/Conda | Unpinned | `mafft --version` |
| Foldseek | Bioconda/Conda | Unpinned | `foldseek version` |
| eggNOG-mapper | External upstream | Unspecified | `emapper.py --version` if installed |
| Prodigal or compatible tool | External upstream | Unspecified | `prodigal -v` if installed |
| ColabFold | Separate optional structure environment | Historical prediction script documents `1.6.1`; **not a general environment pin** | `colabfold_batch --version` inside its dedicated environment |

Recommended capture on WSL/Linux after `conda activate metaviro-drank`:

```bash
python --version
python -c 'import yaml, Bio; print("PyYAML", yaml.__version__, "Biopython", Bio.__version__)'
mmseqs version
hmmscan -h | head -n 4
mafft --version
foldseek version
conda list --explicit > installed-conda-spec.txt
# Keep installed-conda-spec.txt with your run provenance; avoid committing local paths.
```

For reproducibility also record `git rev-parse HEAD`, the exact Swiss-Prot, PHROG, Pfam, PDB, AF_SwissProt and BFVD releases, resource hashes where possible, and the source of protein/GFF/eggNOG input annotations. A current remote `latest` database download is **not** a permanently pinned reference release.

**Hardware sizing** is input- and database-dependent; recommendations in README are planning guidance, not measured acceptance specifications. Validate on the intended target machine before production use.
