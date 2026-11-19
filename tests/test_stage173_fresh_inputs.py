"""New dataset interface smoke: real validation and Stage 01–02 parsers."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(*args, success=True):
    p = subprocess.run([sys.executable, '-B', *map(str, args)], cwd=ROOT,
                       capture_output=True, text=True)
    if success and p.returncode:
        raise AssertionError(f"FAILED: {' '.join(map(str,args))}\n{p.stdout}\n{p.stderr}")
    if not success and p.returncode == 0:
        raise AssertionError("Unexpectedly passed negative test")
    return p


def read(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def test():
    with tempfile.TemporaryDirectory(prefix='metaviro_step173_') as temporary:
        home = Path(temporary)
        fixture = home / 'demo'
        run(ROOT/'examples/generate_demo_input.py', '--out-dir', fixture)
        check = run(ROOT/'run_pipeline.py', '--config', fixture/'config.yaml',
                    '--through','02','--check')
        assert 'INPUT_PREFLIGHT=PASS' in check.stdout
        for sample in ('DEMO_A','DEMO_B'):
            source = fixture/'input'/sample
            one = home/f'{sample}.01.tsv'
            two = home/f'{sample}.02.tsv'
            anchors = home/f'{sample}.anchors.tsv'
            run(ROOT/'scripts/01_define_unresolved_candidates.py', '--annotations', source/'annotations.emapper.annotations',
                '--faa',source/'proteins.faa','--out',one,'--min-aa','150')
            run(ROOT/'scripts/02_extract_replication_context.py','--sample',sample,'--gff',source/'genes.gff',
                '--faa',source/'proteins.faa','--annotations',source/'annotations.emapper.annotations',
                '--unresolved-tsv',one,'--out',two,'--anchors-out',anchors,'--window','2','--min-aa','150')
            rows=read(one); found=read(two); detected=read(anchors)
            assert len(rows)==3 and sum(r['strict_unresolved']=='1' for r in rows)==2
            assert len(found)==2 and len(detected)==1
            assert set(r['anchor_categories'] for r in detected)=={'helicase'}
        # Tamper protein accession => GFF/FASTA contract must fail.
        bad = fixture/'input'/'DEMO_A'/'genes.gff'
        bad.write_text(bad.read_text().replace('ID=1_1','ID=1_99'))
        invalid = run(ROOT/'run_pipeline.py','--config',fixture/'config.yaml','--through','02','--check',success=False)
        assert 'Stage 00 validation failed' in invalid.stderr
    print('STEP173_NEW_INPUT_VALIDATION=PASS')
    print('STEP173_STAGE01_02_REAL_PARSERS=PASS')
    print('STEP173_INVALID_MAPPING_REJECTION=PASS')
    print('STEP173_STAGE03_13=NOT_TESTED')

if __name__=='__main__':
    test()
