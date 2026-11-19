#!/usr/bin/env python3
"""Regression tests for real multi-sample ID collisions."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts' / '00_namespace_sample_ids.py'
CURRENT_VALIDATOR = ROOT / 'scripts' / '00_validate_input_portable.py'
PIPELINE = ROOT / 'run_pipeline.py'


def run(*args, good=True):
    p = subprocess.run(list(map(str, args)), capture_output=True, text=True, check=False)
    if good and p.returncode:
        raise AssertionError(f'FAILED COMMAND: {args}\n{p.stdout}\n{p.stderr}')
    if not good and not p.returncode:
        raise AssertionError(f'SHOULD REJECT: {args}\n{p.stdout}\n{p.stderr}')
    return p


def make(sample):
    sample.mkdir(parents=True)
    # Same protein ID and GFF contig in two independent samples.
    contig = 'k127_999||lt2gene'
    pid = contig + '_1'
    (sample / 'proteins.faa').write_text('>' + pid + ' # 1 # 1200 # 1 # ID=294_1;partial=11\n' + 'MKTL' * 30 + '\n', encoding='utf-8')
    (sample / 'genes.gff').write_text('##gff-version 3\n' + contig + '\tProdigal\tCDS\t1\t360\t.\t+\t0\tID=294_1;partial=11\n', encoding='utf-8')
    (sample / 'annotations.emapper.annotations').write_text('#query\teggNOG_OGs\tDescription\tPreferred_name\n' + pid + '\t10239|Viruses\tDNA polymerase\t-\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    with tempfile.TemporaryDirectory(prefix='metaviro_stage174_namespace_') as tmp:
        base = Path(tmp)
        source = base / 'input'
        for name in ('V695','V696'):
            make(source / name)
        original = {str(p): sha(p) for p in source.rglob('*') if p.is_file()}
        target = base / 'normalized'
        c = run(sys.executable, SCRIPT, '--input-dir', source, '--output-dir', target, '--mode','check')
        assert 'CROSS_SAMPLE_DUPLICATE_IDS=1' in c.stdout
        assert not target.exists()
        run(sys.executable, SCRIPT, '--input-dir', source, '--output-dir', target, '--mode','run')
        j = json.loads((target / 'namespace_manifest.json').read_text())
        assert j['cross_sample_duplicate_count'] == 1
        for sample in ('V695','V696'):
            sd = target / sample
            out = run(sys.executable, CURRENT_VALIDATOR, '--sample-dir', sd)
            assert 'INPUT_VALIDATION=PASS' in out.stdout
            assert (sd / 'proteins.faa').read_text().startswith('>' + sample + '|k127_999||lt2gene_1')
            assert (sample + '|k127_999||lt2gene\t') in (sd / 'genes.gff').read_text()
            assert (sample + '|k127_999||lt2gene_1\t') in (sd / 'annotations.emapper.annotations').read_text()
        assert all(sha(Path(path)) == expected for path, expected in original.items())
        print('STEP174_NAMESPACE_CROSS_SAMPLE_DUPLICATE_RESOLVED=PASS')
        print('STEP174_SOURCE_SHA256_UNCHANGED=PASS')

        # Pipeline Stage 00 on newly normalized input, not precomputed Pilot output.
        config = base / 'config.yaml'
        config.write_text('input_dir: ' + str(target) + '\noutput_dir: ' + str(base/'workflow') + '\nmin_aa: 150\nwindow: 5\nthreads: 2\n')
        result = run(sys.executable, PIPELINE, '--config', config, '--through','04','--check')
        assert 'INPUT_PREFLIGHT=PASS' in result.stdout
        print('STEP174_NORMALIZED_PIPELINE_STAGE00_PREFLIGHT=PASS')

        # No silent overwrite; sources with invalid cross-reference must fail.
        run(sys.executable, SCRIPT, '--input-dir', source, '--output-dir', target, '--mode','run', good=False)
        bad = base / 'bad'
        make(bad / 'V999')
        (bad/'V999'/'genes.gff').write_text('##gff-version 3\nWRONG\tProdigal\tCDS\t1\t360\t.\t+\t0\tID=294_1\n')
        bad_target = base / 'bad_out'
        err = run(sys.executable, SCRIPT, '--input-dir', bad, '--output-dir', bad_target, '--mode','run', good=False)
        assert 'FASTA/GFF mismatch' in err.stderr
        assert not bad_target.exists()
        print('STEP174_INVALID_GFF_MAPPING_REJECTED=PASS')
        print('STEP174_OUTPUT_OVERWRITE_REJECTED=PASS')
    print('ALL_TESTS_PASS')

if __name__ == '__main__':
    main()
