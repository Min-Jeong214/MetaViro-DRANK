# Validation and reproducibility tests

## Repeatable software checks

The repository provides the Stage 00–07 workflow, a CPU context-review entrypoint, and separate structure/sequence evidence modules. Its tests distinguish parser-level execution, integration fixtures and database-dependent analyses. Tests are independent of the research datasets used in the associated publication.

| Verification path | Command | Requirement |
|---|---|---|
| Fresh GitHub checkout | `python -B scripts/validate_fresh_clone.py --mode smoke` | Python 3.11 + PyYAML; no biological databases |
| Synthetic regression | `python -m unittest discover -s tests -p 'test_release_smoke.py' -v` | Python 3.11 + PyYAML |
| Database/Stage 00–07 preflight | `python -B scripts/validate_fresh_clone.py --mode database --config config.local.yaml --screen-out-dir results/validation_context --normalize-ids` | Local inputs and full reference DBs |
| Database-backed execution | The previous command with `--execute` and unused output paths | MMseqs2, HMMER, search DBs, CPU/RAM/disk |
| Stage 08–13 module fixtures | `python -B tests/test_stage185_partial_evidence_catalog.py` and related test scripts | Module-specific synthetic/mock evidence fixtures |

GitHub Actions performs the synthetic fresh-clone smoke on each push and pull request, separately from Python source syntax checks. No biological reference databases or user research datasets are downloaded by CI.

## What tests can establish

- Python source parsing and interface compatibility at the test-fixture level.
- Input schema/ID checks, deterministic parameter handling and safeguards in covered paths.
- Actual synthetic Stage 00–02 execution for **both one and three sample inputs**, with output count/provenance assertions and a negative ID-mismatch test.
- Functional evidence integration *only within the explicitly tested scope* of the supplied code/test fixtures.

## Scientific and external-dependency verification

- The synthetic smoke proves fresh-source usability with Python/PyYAML; complete Conda/DB installation must also be tested on the intended Ubuntu/WSL2 system.
- Complete uninterrupted end-to-end execution of every optional Stage 08–13 process on an arbitrary new research dataset.
- A validated nuclease activity classifier, substrate-specific assay, catalytic-site assignment or experimentally verified enzyme function.
- Direct input support for FASTQ or bare contig FASTA.

## Location of study-specific evidence

Scientific cohort descriptions, figures, benchmarking, evidence counts and supplementary validation are documented in the associated manuscript (**submitted; under editorial consideration**) and the **Zenodo** validation dataset (**currently a draft; DOI/URL to be added after publication**). They are intentionally excluded from the general README and repository-facing validation summary. No particular research cohort is a fixed prerequisite for using the pipeline.

**Reporting rule:** distinguish `PASS` of software tests from scientific validation, and preserve `NOT_EVALUATED` / `PENDING` for unavailable structural or biochemical evidence rather than interpreting it as a negative result.
