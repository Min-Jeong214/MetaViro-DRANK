#!/usr/bin/env python3
"""Single entrypoint: verified Stage05-07 -> Stage08 planning -> CPU context review.

Optional historical structure index is for reuse assessment, never required.
This runner NEVER invokes GPU prediction, Foldseek or final Stage10/13 ranking.
"""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
PREP=ROOT/'scripts/08_prepare_structure_queries_portable.py'
PLAN=ROOT/'scripts/08_partial_structure_plan_portable.py'
SCREEN=ROOT/'scripts/10_context_first_screen_portable.py'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def tsv(p,keys):
    if not p.is_file() or not p.stat().st_size:raise ValueError(f'MISSING_OR_EMPTY:{p}')
    with p.open(encoding='utf-8-sig',newline='') as f:
        rd=csv.DictReader(f,delimiter='\t'); fields=rd.fieldnames or []
        if not set(keys)<=set(fields) or len(fields)!=len(set(fields)):
            raise ValueError(f'INVALID_TSV_SCHEMA:{p}')
        rows=list(rd)
    if not rows:raise ValueError(f'EMPTY_TSV:{p}')
    return rows

def check_stage57(workflow):
    files={
        'members': workflow/'05_annotation/recurrent_members.tsv',
        'fasta': workflow/'05_annotation/recurrent_members.faa',
        'audit': workflow/'06_darkness/family_darkness_audit.tsv',
        'family': workflow/'07_family_context/family_context.tsv',
        'member': workflow/'07_family_context/member_context.tsv',
    }
    for p in (PREP,PLAN,SCREEN):
        if not p.is_file():raise ValueError(f'SCRIPT_NOT_INSTALLED:{p}')
    m=tsv(files['members'],('member','representative'))
    a=tsv(files['audit'],('representative','family_size','priority_class'))
    f=tsv(files['family'],('family_rep','family_size'))
    c=tsv(files['member'],('member','family_rep','sample','anchor_within_window','candidate_original','anchor_categories','window'))
    if len({r['member'] for r in m})!=len(m) or len({r['member'] for r in c})!=len(c):raise ValueError('DUPLICATE_MEMBER')
    if {r['member'] for r in m}!={r['member'] for r in c}:raise ValueError('MEMBER_ID_SET_MISMATCH')
    if len({r['representative'] for r in a})!=len(a) or len({r['family_rep'] for r in f})!=len(f):raise ValueError('DUPLICATE_FAMILY')
    if {r['representative'] for r in a}!={r['family_rep'] for r in f}:raise ValueError('FAMILY_ID_SET_MISMATCH')
    if not files['fasta'].is_file() or not files['fasta'].stat().st_size:raise ValueError('MISSING_PROTEIN_FASTA')
    print(f'STAGE05_07_FAMILIES={len(a)}')
    print(f'STAGE05_07_MEMBERS={len(m)}')
    print('STAGE05_07_PRECONDITIONS=PASS')
    return files,len(a),len(m)

def go(label,command):
    print(f'[RUN] {label}',flush=True)
    result=subprocess.run([sys.executable,'-B',*map(str,command)],cwd=ROOT)
    if result.returncode:raise RuntimeError(f'SUBPROCESS_FAILED:{label}:exit={result.returncode}')
    print(f'[PASS] {label}',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workflow-output',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    p.add_argument('--historical-structures',type=Path)
    p.add_argument('--limit',type=int,default=30)
    p.add_argument('--mode',choices=('check','run'),default='check')
    a=p.parse_args()
    workflow=a.workflow_output.resolve(); out=a.out_dir.resolve()
    if a.limit<1:raise ValueError('INVALID_REVIEW_LIMIT')
    if out.exists():raise ValueError(f'REFUSE_EXISTING_OUTPUT:{out}')
    if out==workflow or workflow in out.parents or out in workflow.parents:
        # Allow sibling stage workspaces, not nested inside source workflow.
        raise ValueError('INPUT_OUTPUT_OVERLAP')
    files,nf,nm=check_stage57(workflow)
    if a.historical_structures and not a.historical_structures.is_file():
        raise ValueError('HISTORICAL_STRUCTURE_TABLE_MISSING')
    if a.mode=='check':
        print('STAGE180_UNIFIED_PREFLIGHT=PASS');print('FILES_WRITTEN=NO');return
    out.mkdir(parents=True,exist_ok=False)
    prepared=out/'08_prepared'; plan=out/'08_partial_plan';screen=out/'10_context_screen'
    go('Stage08 structure queries', [PREP,'--members',files['members'],'--fasta',files['fasta'],'--family-audit',files['audit'],'--out-dir',prepared])
    options=[PLAN,'--mapping',prepared/'structure_id_map.tsv','--queries',prepared/'representative_queries.faa','--family-audit',files['audit'],'--family-context',files['family'],'--limit',a.limit,'--out-dir',plan,'--mode','run']
    if a.historical_structures:options.extend(['--historical-structures',a.historical_structures])
    go('Partial structure plan',options)
    go('Context-first review', [SCREEN,'--partial-plan-dir',plan,'--member-context',files['member'],'--family-context',files['family'],'--out-dir',screen,'--mode','run'])
    inf=tsv(plan/'partial_structure_inventory.tsv',('structure_id',))
    screens=tsv(screen/'all_family_context_screen.tsv',('structure_id',))
    if len(inf)!=len(screens) or len(inf)!=nf:raise ValueError('FAMILY_COUNT_MISMATCH_AFTER_SCREEN')
    status=json.loads((screen/'context_first_provenance.json').read_text())
    mainfest={
        'workflow':'STRUCTURE_OPTIONAL_CPU_FIRST_PRELIMINARY_SCREEN',
        'source_workflow':str(workflow),
        'source_families':nf,'source_members':nm,
        'family_coverage_verified':True,
        'structure_predictions_executed':False,'external_structure_search_executed':False,
        'stage10_full_ranking_completed':False,'stage13_completed':False,
        'historical_reuse_only_optional':True,
        'review_limit':a.limit,
        'context_supported_dark_families':status['context_supported_dark_families'],
        'input_sha256':{k:sha(v) for k,v in files.items()},
        'child_manifests_sha256':{
            'stage179':sha(plan/'partial_structure_plan_provenance.json'),
            'stage180':sha(screen/'context_first_provenance.json')},
    }
    (out/'context_screen_run_provenance.json').write_text(json.dumps(mainfest,indent=2)+'\n')
    print('STAGE180_UNIFIED_CONTEXT_SCREEN=PASS')
    print(f'ALL_FAMILIES={nf}')
    print(f'CONTEXT_SUPPORTED_DARK_FAMILIES={status["context_supported_dark_families"]}')
    print(f'OUTPUT_DIR={out}')
    print('FOLDSEEK_GPU_REEXECUTED=NO')
    print('STAGE10_FINAL_RANKING=NOT_RUN')

if __name__=='__main__':main()
