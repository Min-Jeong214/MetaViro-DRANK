#!/usr/bin/env python3
"""Portable single-command software regression with mock new viral cohort."""
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PKG=Path(__file__).resolve().parents[1]
STAGE08=PKG/'scripts/08_prepare_structure_queries_portable.py'
STAGE179=PKG/'scripts/08_partial_structure_plan_portable.py'
# Isolated package QA fallback; an installed GitHub checkout uses PKG/scripts.
if not STAGE08.is_file():
    STAGE08=PKG.parent/'metaviro_step172_review/scripts/08_prepare_structure_queries_portable.py'
if not STAGE179.is_file():
    STAGE179=PKG.parent/'stage179_release/scripts/08_partial_structure_plan_portable.py'

def tab(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)

def run(path,*args):
    return subprocess.run([sys.executable,'-B',str(path),*map(str,args)],text=True,capture_output=True,timeout=35)

with tempfile.TemporaryDirectory(prefix='metaviro_stage180_unified_') as tmp:
    root=Path(tmp)
    repo=root/'repo';repo.mkdir();(repo/'scripts').mkdir()
    originals=(
        ('run_context_screen.py',PKG/'run_context_screen.py'),
        ('scripts/08_prepare_structure_queries_portable.py',STAGE08),
        ('scripts/08_partial_structure_plan_portable.py',STAGE179),
        ('scripts/10_context_first_screen_portable.py',PKG/'scripts/10_context_first_screen_portable.py'),
    )
    for dst,src in originals:
        assert src.is_file(),src
        shutil.copyfile(src,repo/dst)
    work=root/'new_cohort_workflow'
    members=[];audit=[];fm=[];context=[];seq=[]
    for idx,klass in ((1,'A_DEEP_DARK'),(2,'D_KNOWN_EVIDENCE')):
        rep=f'S{idx}_protein1';ids=[rep,f'S{idx}_protein2',f'S{idx}_protein3']
        for j,m in enumerate(ids):
            members.append({'member':m,'representative':rep})
            context.append({'family_rep':rep,'member':m,'sample':('V1','V2','V1')[j],'candidate_original':('1' if j==2 else '0'),'anchor_within_window':('1' if j<2 else '0'),'anchor_categories':('nuclease' if j<2 else ''),'window':'5'})
            seq.append((m,'M'+ 'ACDEFGHIKLMNPQRSTVWY'[idx:idx+7]*3+str(idx).replace('1','A').replace('2','C')))
        audit.append({'representative':rep,'family_size':'3','priority_class':klass})
        fm.append({'family_rep':rep,'family_size':'3','sample_span':'2','noncandidate_context_sample_span':'2','all_context_sample_span':'2','noncandidate_anchor_context_count':'2','all_anchor_context_count':'2','dominant_anchor_category':'nuclease'})
    tab(work/'05_annotation/recurrent_members.tsv',members)
    fa=work/'05_annotation/recurrent_members.faa'
    fa.write_text(''.join('>'+sid+'\n'+aa+'\n' for sid,aa in seq))
    tab(work/'06_darkness/family_darkness_audit.tsv',audit)
    tab(work/'07_family_context/family_context.tsv',fm)
    tab(work/'07_family_context/member_context.tsv',context)
    out=root/'result'
    r=run(repo/'run_context_screen.py','--workflow-output',work,'--out-dir',out,'--mode','check')
    assert r.returncode==0,(r.stdout,r.stderr)
    assert not out.exists()
    print('STAGE180_UNIFIED_PREFLIGHT_NO_WRITE=PASS')
    r=run(repo/'run_context_screen.py','--workflow-output',work,'--out-dir',out,'--mode','run')
    assert r.returncode==0,(r.stdout,r.stderr)
    assert 'STAGE180_UNIFIED_CONTEXT_SCREEN=PASS' in r.stdout
    with (out/'10_context_screen/all_family_context_screen.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==2 and {x['family_rep'] for x in rows}=={'S1_protein1','S2_protein1'}
    assert (out/'08_partial_plan/new_prediction_review_queries.faa').exists()
    meta=json.loads((out/'context_screen_run_provenance.json').read_text())
    assert meta['source_families']==2 and meta['context_supported_dark_families']==1
    assert meta['structure_predictions_executed'] is False
    print('STAGE180_UNIFIED_NEW_COHORT_E2E=PASS')
    assert run(repo/'run_context_screen.py','--workflow-output',work,'--out-dir',out,'--mode','run').returncode!=0
    print('STAGE180_UNIFIED_OUTPUT_OVERWRITE_REJECTED=PASS')
print('ALL_TESTS_PASS')
