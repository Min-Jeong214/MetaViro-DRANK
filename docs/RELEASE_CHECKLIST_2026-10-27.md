# MetaViro-DRANK research-software release checklist

**Release maturity:** `PARTIAL_VALIDATED`. The recommended CPU-first workflow is distinct from optional and partly automated structure/MSA/evidence integration. No experimentally validated nuclease function is claimed.

## Software quality gate

- [ ] Fresh checkout and Conda environment installation on a supported Ubuntu/WSL2 platform
- [ ] CLI `--help` and source-level automated tests pass
- [ ] Input fixture, schema, and Stage 00–02 smoke tests pass
- [ ] Independent arbitrary-cohort Stage 00–07 and context screening execution pass
- [ ] Correct and versioned Swiss-Prot, PHROG profile, Pfam resources verified
- [ ] Optional ColabFold/Foldseek/Stage 08–13 workflows explicitly tested where claimed
- [ ] Status/provenance clearly distinguishes results, pending analyses and missing evidence

## Documentation, research and rights gate

- [x] General-purpose README covering principles, requirements, install, usage and outputs
- [x] Scientific cohort-specific numbers omitted from general-purpose README
- [ ] Corresponding manuscript metadata added after peer review
- [ ] Zenodo record DOI/URL and licensing of validation datasets confirmed
- [ ] KIOST copyright registration and institutional rights/permissions established
- [ ] External software and database licenses checked
- [ ] Official version tag or release issued only after the preceding checks

## Pre-push review

```bash
git status --short
git diff --check
git diff --stat
```

Do not stage sequencing data, local inputs, reference DBs, generated outputs, unpublished structures, credentials or local configurations. Historical study-specific validation details are retained with the research materials and manuscript/Zenodo record rather than release-facing repository documentation.
