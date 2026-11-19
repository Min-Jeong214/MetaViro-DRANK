# Example: reproducible synthetic input-format smoke test

```bash
python examples/generate_demo_input.py --out-dir results/demo_step173
python run_pipeline.py --config results/demo_step173/config.yaml --through 02 --check
```

This generator makes two **artificial sample directories** with matching Prodigal-style protein/GFF IDs and eggNOG tabular headers. It does not contain natural viral sequences, research data or verified enzymes, and it is not a Stage 03–13 biological test.

`--through 02` (no `--check`) can run when system provenance tools MMseqs2 and HMMER are installed. Avoid treating demo results as evidence for candidate discovery or external databases.
