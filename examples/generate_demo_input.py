#!/usr/bin/env python3
"""Generate a non-biological, ID-consistent smoke fixture for Stage 00–02."""
import argparse
from pathlib import Path
import json
import re


def generate(root, sample_count=2, prefix="DEMO"):
    if not 1 <= sample_count <= 26:
        raise ValueError("sample_count must be from 1 to 26")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", prefix):
        raise ValueError("prefix must contain letters, numbers and underscores")
    root = root.expanduser().resolve()
    if root.exists():
        raise ValueError(f"Refusing existing directory: {root}")
    if root == Path.cwd().resolve() or root in Path.cwd().resolve().parents:
        raise ValueError("Unsafe output directory")
    # The runner resolves config relative paths against its own repository;
    # absolute paths avoid unintended location changes.
    out = root / "analysis_output"
    source = root / "input"
    for sample in (f"{prefix}_{chr(65+i)}" for i in range(sample_count)):
        path = source / sample
        path.mkdir(parents=True)
        contig = f"{sample}_contig"
        names = [f"{contig}_{i}" for i in (1, 2, 3)]
        (path / "proteins.faa").write_text(
            "".join(f">{n}\n{'M' * (155+i)}\n" for i, n in enumerate(names)), encoding="utf-8"
        )
        (path / "genes.gff").write_text(
            "##gff-version 3\n" + "".join(
                f"{contig}\tDemo\tCDS\t{1+600*i}\t{480+600*i}\t.\t+\t0\tID=1_{i+1}\n"
                for i in range(3)
            ), encoding="utf-8"
        )
        # Anchor protein at index 2; 1 and 3 are intentionally unannotated.
        (path / "annotations.emapper.annotations").write_text(
            "#query\teggNOG_OGs\tDescription\tPreferred_name\tKEGG_ko\tPFAMs\n"
            + f"{names[1]}\t10239|Viruses\tDNA helicase\t-\t-\t-\n",
            encoding="utf-8"
        )
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required; conda env create -f environment.yml") from exc
    config = {
        "input_dir": str(source), "output_dir": str(out),
        "min_aa": 150, "window": 2, "threads": 2,
    }
    (root / "config.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf8"
    )
    (root / "FIXTURE_NOTE.json").write_text(
        json.dumps({"kind": "synthetic_nonbiological_smoke", "expected_samples":sample_count,
                    "expected_unresolved_per_sample":2, "expected_anchors_per_sample":1}, indent=2) + "\n"
    )
    print(f"DEMO_INPUT={source}")
    print(f"DEMO_CONFIG={root / 'config.yaml'}")
    print("DEMO_GENERATED=PASS")
    return root


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--sample-count", type=int, default=2)
    ap.add_argument("--prefix", default="DEMO")
    args = ap.parse_args()
    generate(args.out_dir, sample_count=args.sample_count, prefix=args.prefix)
