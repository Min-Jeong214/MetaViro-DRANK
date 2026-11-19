#!/usr/bin/env python3
"""Official structure-optional CPU workflow for new marine viral protein cohorts.

New run: raw/canonical FASTA+GFF+eggNOG -> Stage00-07 -> context-first review.
Resume: a verified Stage05-07 workflow -> context-first review.
No ColabFold, GPU, Foldseek, or biochemical function assignment is done here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
UPSTREAM = ROOT / "run_pipeline.py"
NAMESPACE = ROOT / "scripts/00_namespace_sample_ids.py"
SCREEN = ROOT / "run_context_screen.py"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def resolve(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else ROOT / path).resolve()


def check_file(path: Path, purpose: str) -> None:
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f'MISSING_{purpose}:{path}')


def disjoint(one: Path, another: Path) -> bool:
    return not (one == another or one in another.parents or another in one.parents)


def command(label: str, argv: list[str | Path]) -> None:
    print(f'[RUN] {label}', flush=True)
    result = subprocess.run([str(item) for item in argv], cwd=ROOT, check=False)
    if result.returncode:
        raise RuntimeError(f'COMMAND_FAILED:{label}:exit={result.returncode}')
    print(f'[PASS] {label}', flush=True)


def config_info(file: Path):
    check_file(file, 'CONFIG')
    data = yaml.safe_load(file.read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise ValueError('INVALID_CONFIG_YAML')
    for key in ('input_dir', 'output_dir', 'min_aa', 'window'):
        if key not in data:
            raise ValueError(f'CONFIG_KEY_MISSING:{key}')
    source = resolve(data['input_dir'])
    upstream = resolve(data['output_dir'])
    if not source.is_dir() or not disjoint(source, upstream):
        raise ValueError('INVALID_INPUT_OR_OVERLAP')
    if int(data['min_aa']) < 1 or int(data['window']) < 0:
        raise ValueError('INVALID_CONFIG_THRESHOLDS')
    dbs = data.get('databases')
    if not isinstance(dbs, dict):
        raise ValueError('MISSING_DATABASE_CONFIGURATION')
    for key in ('swissprot','phrog','pfam','phrog_annotations'):
        if not isinstance(dbs.get(key), str):
            raise ValueError(f'DATABASE_KEY_MISSING:{key}')
        db = resolve(dbs[key])
        check_file(db, f'DATABASE_{key}')
        if key in ('swissprot', 'phrog'):
            check_file(Path(str(db)+'.dbtype'), key.upper()+'_DBTYPE')
        if key == 'phrog':
            for ext in ('.index', '_h', '_h.dbtype', '_h.index'):
                check_file(Path(str(db)+ext), 'PHROG_INDEX')
        if key == 'pfam':
            for ext in ('.h3f','.h3i','.h3m','.h3p'):
                check_file(Path(str(db)+ext), 'PFAM_INDEX')
    return data, source, upstream


def ensure_destination(path: Path, sources: list[Path]):
    if path.exists():
        raise ValueError(f'OUTPUT_ALREADY_EXISTS:{path}')
    for src in sources:
        if not disjoint(src, path):
            raise ValueError(f'OUTPUT_OVERLAPS_SOURCE:{src}:{path}')


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='operation', required=True)
    new = sub.add_parser('run', help='New FASTA/GFF/eggNOG cohort through CPU screening')
    new.add_argument('--config', type=Path, required=True)
    new.add_argument('--screen-out-dir', type=Path, required=True)
    new.add_argument('--normalize-ids', action='store_true', help='Namespace FAA/GFF/eggNOG IDs by sample in safe copies')
    resume = sub.add_parser('resume', help='Screen already completed Stage05-07 workflow without rerunning searches')
    resume.add_argument('--workflow-output', type=Path, required=True)
    resume.add_argument('--screen-out-dir', type=Path, required=True)
    for ap in (new,resume):
        ap.add_argument('--historical-structures', type=Path, help='Optional verified structure table (same-sequence SHA only)')
        ap.add_argument('--limit', type=int, default=30)
        ap.add_argument('--mode', choices=('check','run'), default='check')
    return parser.parse_args()


def main():
    args = cli()
    screen_out = resolve(args.screen_out_dir)
    if args.limit < 1:
        raise ValueError('INVALID_REVIEW_LIMIT')
    hist = resolve(args.historical_structures) if args.historical_structures else None
    if hist:
        check_file(hist, 'HISTORICAL_STRUCTURE_TABLE')
    for script in (SCREEN,):
        check_file(script, 'INSTALLED_SCRIPT')
    if args.operation == 'resume':
        workflow = resolve(args.workflow_output)
        if not workflow.is_dir():
            raise ValueError(f'WORKFLOW_MISSING:{workflow}')
        ensure_destination(screen_out, [workflow] + ([hist] if hist else []))
        screen_args = [sys.executable, '-B', SCREEN, '--workflow-output', workflow,
                       '--out-dir', screen_out / 'review', '--limit', str(args.limit)]
        if hist:
            screen_args += ['--historical-structures', hist]
        command('Existing Stage05-07 CPU screening preflight', screen_args + ['--mode', 'check'])
        if args.mode == 'check':
            print('STEP181_RESUME_PREFLIGHT=PASS\nFILES_WRITTEN=NO')
            return
        screen_out.mkdir(parents=True, exist_ok=False)
        command('CPU context-first screening', screen_args + ['--mode', 'run'])
        mode = 'RESUME_VERIFIED_STAGE05_07'
        input_config_hash = None
    else:
        check_file(UPSTREAM, 'UPSTREAM_RUNNER')
        conf = resolve(args.config)
        settings, raw_input, workflow = config_info(conf)
        ensure_destination(workflow, [raw_input, conf, screen_out] + ([hist] if hist else []))
        ensure_destination(screen_out, [raw_input, workflow, conf] + ([hist] if hist else []))
        input_config_hash = digest(conf)
        effective_config = conf
        if args.normalize_ids:
            check_file(NAMESPACE, 'NAMESPACE_MODULE')
            normal_args = [sys.executable, '-B', NAMESPACE,
                           '--input-dir', raw_input,
                           '--output-dir', screen_out / 'normalized_input']
            command('Safe input ID namespace preflight', normal_args + ['--mode', 'check'])
        else:
            command('New cohort Stage00-07 preflight',
                    [sys.executable, '-B', UPSTREAM, '--config', conf,
                     '--through', '07', '--check'])
        if args.mode == 'check':
            print('STEP181_NEW_COHORT_PREFLIGHT=PASS')
            if args.normalize_ids:
                print('UPSTREAM_PREFLIGHT=DEFERRED_UNTIL_NAMESPACE_CREATED')
            print('FILES_WRITTEN=NO')
            return
        screen_out.mkdir(parents=True, exist_ok=False)
        if args.normalize_ids:
            command('Namespace input copies', normal_args + ['--mode', 'run'])
            settings['input_dir'] = str(screen_out / 'normalized_input')
            settings['output_dir'] = str(workflow)
            settings['databases'] = {key:str(resolve(path)) for key,path in settings['databases'].items()}
            effective_config = screen_out / 'effective_config.yaml'
            with effective_config.open('x', encoding='utf-8') as stream:
                yaml.safe_dump(settings, stream, sort_keys=False, allow_unicode=True)
            command('New cohort normalized Stage00-07 preflight',
                    [sys.executable, '-B', UPSTREAM, '--config', effective_config,
                     '--through', '07', '--check'])
        command('New cohort Stage00-07 analysis',
                [sys.executable, '-B', UPSTREAM, '--config', effective_config,
                 '--through', '07'])
        screen_args = [sys.executable, '-B', SCREEN, '--workflow-output', workflow,
                       '--out-dir', screen_out / 'review', '--limit', str(args.limit)]
        if hist:
            screen_args += ['--historical-structures', hist]
        command('New cohort structure-optional CPU screen preflight',
                screen_args + ['--mode', 'check'])
        command('New cohort structure-optional CPU screen',
                screen_args + ['--mode', 'run'])
        mode = 'NEW_INPUT_STAGE00_07_AND_CPU_SCREEN'
    child = screen_out / 'review/context_screen_run_provenance.json'
    check_file(child, 'CPU_SCREEN_PROVENANCE')
    result = json.loads(child.read_text(encoding='utf-8'))
    if result.get('family_coverage_verified') is not True or result.get('stage10_full_ranking_completed') is not False:
        raise ValueError('INVALID_CHILD_RUN_PROVENANCE')
    meta = {
        'workflow_mode':mode,
        'stage07_workflow':str(workflow),
        'stage07_completed':True,
        'context_screen_completed':True,
        'family_count':result['source_families'],
        'member_count':result['source_members'],
        'context_supported_dark_families':result['context_supported_dark_families'],
        'historical_structure_reuse_optional':True,
        'structure_prediction_executed':False,
        'external_structure_search_executed':False,
        'stage10_final_ranking_completed':False,
        'stage13_full_integration_completed':False,
        'original_config_sha256':input_config_hash,
        'cpu_screen_manifest_sha256':digest(child),
    }
    (screen_out/'run_metaviro_provenance.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print('STEP181_MAIN_CLI=PASS')
    print('MODE='+mode)
    print('FAMILY_COUNT='+str(meta['family_count']))
    print('CONTEXT_SUPPORTED_DARK='+str(meta['context_supported_dark_families']))
    print('GPU_FOLDSEEK=NOT_RUN')
    print('STAGE10_FINAL_RANKING=NOT_RUN')
    print('GITHUB_PUSH=NOT_PERFORMED')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as exc:
        sys.exit(f'ERROR: {exc}')
