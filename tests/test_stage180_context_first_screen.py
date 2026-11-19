#!/usr/bin/env python3
"""Synthetic contract tests for portable context-first screening (no DB/GPU)."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/10_context_first_screen_portable.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
def read(p):
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def execute(*args):
    return subprocess.run([sys.executable,'-B',str(SCRIPT),*map(str,args)],capture_output=True,text=True)
with tempfile.TemporaryDirectory(prefix='metaviro_stage180_') as td:
    root=Path(td); plan=root/'plan';plan.mkdir()
    inventory=[]; family=[]; members=[]
    classes=['A_DEEP_DARK','B_FUNCTION_UNRESOLVED','C_KNOWN_DARK_BRIDGE','D_KNOWN_EVIDENCE']
    for i,klass in enumerate(classes,1):
        rep=f'S{i}|rep'; sid=f'MVF{i:06d}'
        tag='VERIFIED_HISTORICAL_REUSE_CANDIDATE' if i==1 else 'NO_STRUCTURE_AVAILABLE'
        categories={1:['nuclease;helicase','nuclease'],2:['helicase',''],3:['',''],4:['nuclease','nuclease']}[i]
        for j in range(3):
            candidate='1' if j==2 else '0'
            cat=categories[j] if j<2 else ''
            members.append(dict(family_rep=rep,member=f'{rep}_gene{j}',sample=('V1','V2','V1')[j],candidate_original=candidate,anchor_within_window=('1' if cat else '0'),anchor_categories=cat,window='5'))
        relevant=[r for r in members if r['family_rep']==rep]
        anchors=[r for r in relevant if r['anchor_within_window']=='1']
        samp=len({r['sample'] for r in anchors})
        family.append(dict(family_rep=rep,family_size='3',sample_span='2',all_context_sample_span=str(samp),noncandidate_context_sample_span=str(samp),noncandidate_anchor_context_count=str(len(anchors)),all_anchor_context_count=str(len(anchors)),dominant_anchor_category=('nuclease' if i in (1,4) else 'helicase' if i==2 else '')))
        inventory.append(dict(structure_id=sid,family_rep=rep,priority_class=klass,family_size='3',sample_span='2',noncandidate_context_sample_span=str(samp),all_context_sample_span=str(samp),noncandidate_anchor_context_count=str(len(anchors)),all_anchor_context_count=str(len(anchors)),dominant_anchor_category=family[-1]['dominant_anchor_category'],structure_status=tag,prediction_selection='REUSE_CANDIDATE_NOT_NEW_PREDICTION' if i==1 else 'DEFERRED_NO_STRUCTURE'))
    inv=plan/'partial_structure_inventory.tsv'; fc=root/'family_context.tsv'; mc=root/'member_context.tsv'
    write(inv,inventory);write(fc,family);write(mc,members)
    metadata=dict(input_sha256={'family_context':sha(fc)},output_sha256={'partial_structure_inventory.tsv':sha(inv)},gpu_prediction_performed=False,external_search_performed=False)
    manifest=plan/'partial_structure_plan_provenance.json';manifest.write_text(json.dumps(metadata))
    common=['--partial-plan-dir',plan,'--family-context',fc,'--member-context',mc]
    output=root/'output'
    result=execute(*common,'--out-dir',output,'--mode','check')
    assert result.returncode==0,result.stderr
    assert 'DARK_CONTEXT_SUPPORTED_FAMILIES=2' in result.stdout
    assert not output.exists()
    print('STAGE180_CHECK_MODE_NO_WRITE=PASS')
    result=execute(*common,'--out-dir',output,'--mode','run')
    assert result.returncode==0,result.stderr
    allrows=read(output/'all_family_context_screen.tsv'); screen=read(output/'dark_replication_context_review.tsv')
    assert len(allrows)==4 and len(screen)==2
    assert [r['structure_id'] for r in screen]==['MVF000001','MVF000002']
    assert screen[0]['best_noncandidate_anchor_category']=='nuclease'
    assert screen[0]['category_supported_samples']=='2'
    assert allrows[0]['structure_status']=='VERIFIED_HISTORICAL_REUSE_CANDIDATE'
    assert allrows[0]['screening_interpretation']=='CONTEXT_ASSOCIATION_NOT_CONFIRMED_ENZYME_ACTIVITY'
    status_by_id={r['structure_id']:r['screening_status'] for r in allrows}
    assert status_by_id['MVF000003']=='DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY'
    assert status_by_id['MVF000004']=='KNOWN_EVIDENCE_FAMILY'
    manifest_out=json.loads((output/'context_first_provenance.json').read_text())
    assert manifest_out['not_stage10_i2_ranking'] and not manifest_out['gpu_prediction_performed']
    assert sha(output/'all_family_context_screen.tsv')==manifest_out['output_sha256']['all_family_context_screen.tsv']
    print('STAGE180_VALID_DARK_SCREEN_WITHOUT_STRUCTURE=PASS')
    assert execute(*common,'--out-dir',output,'--mode','run').returncode!=0
    print('STAGE180_OUTPUT_OVERWRITE_REJECTED=PASS')
    # Historical plan evidence tampered after compilation: detect by recorded SHA.
    inv.write_text(inv.read_text()+'X\n')
    result=execute(*common,'--out-dir',root/'bad1','--mode','check')
    assert result.returncode!=0 and 'STAGE179_INVENTORY_SHA_MISMATCH' in result.stderr,result.stderr
    print('STAGE180_STAGE179_TAMPER_REJECTED=PASS')
    write(inv,inventory)
    # Source Stage07 family table changed: detect before downstream join.
    fc.write_text(fc.read_text()+'X\n')
    result=execute(*common,'--out-dir',root/'bad2','--mode','check')
    assert result.returncode!=0 and 'STAGE179_STAGE07_FAMILY_SHA_MISMATCH' in result.stderr,result.stderr
    print('STAGE180_SOURCE_CONTEXT_TAMPER_REJECTED=PASS')
    write(fc,family)
    broken=[dict(r) for r in members]
    broken[0]['anchor_within_window']='0'
    write(mc,broken)
    result=execute(*common,'--out-dir',root/'bad3','--mode','check')
    assert result.returncode!=0 and 'CONTEXT_COUNT_MISMATCH' in result.stderr,result.stderr
    print('STAGE180_MEMBER_CONTEXT_INCONSISTENCY_REJECTED=PASS')
print('ALL_TESTS_PASS')
