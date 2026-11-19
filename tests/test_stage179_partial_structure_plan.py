#!/usr/bin/env python3
"""No GPU or database needed; validate strict Stage179 local input contracts."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/08_partial_structure_plan_portable.py'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def table(p, fields, rows):
    with p.open('w', encoding='utf-8',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def run(*args):
    return subprocess.run([sys.executable,'-B',str(SCRIPT),*map(str,args)],capture_output=True,text=True)

with tempfile.TemporaryDirectory(prefix='metaviro_stage179_') as td:
    root=Path(td)
    mapping,query,audit,context,old=(root/n for n in ('map.tsv','q.faa','audit.tsv','context.tsv','old.tsv'))
    labels=['A_DEEP_DARK','A_DEEP_DARK','B_FUNCTION_UNRESOLVED','C_KNOWN_DARK_BRIDGE','D_KNOWN_EVIDENCE','D_KNOWN_EVIDENCE']
    m=[]; a=[]; c=[]; sequences={}
    for i,priority in enumerate(labels,1):
        sid=f'MVF{i:06d}'; rep=f'S{i}|rep'
        seq='M'+('ACDEFGHIKLMNPQRSTVWY'[i:i+7])*2
        sequences[sid]=seq
        m.append(dict(structure_id=sid,family_rep=rep,family_size='3',length_aa=str(len(seq)),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),priority_class=priority))
        a.append(dict(representative=rep,family_size='3',priority_class=priority))
        c.append(dict(family_rep=rep,family_size='3',sample_span='5',noncandidate_context_sample_span=str(5-i if i<4 else 1),all_context_sample_span='5',all_anchor_context_count=str(10+i),noncandidate_anchor_context_count=str(5+i),dominant_anchor_category='nuclease'))
    table(mapping,list(m[0]),m)
    table(audit,list(a[0]),a)
    table(context,list(c[0]),c)
    query.write_text(''.join(f'>{k}\n{v}\n' for k,v in sequences.items()))
    pdb=root/'old.pdb'; score=root/'old.json'
    pdb.write_text('MODEL   1\nENDMDL\n'); score.write_text('{"ptm": 0.8}\n')
    table(old,['structure_id','sequence_sha256','length_aa','prediction_status','pdb_file','pdb_sha256','score_json','score_json_sha256','mean_plddt','ptm'], [
        dict(structure_id='FOLD',sequence_sha256=m[0]['sequence_sha256'],length_aa=m[0]['length_aa'],prediction_status='PREDICTED',pdb_file=str(pdb),pdb_sha256=sha(pdb),score_json=str(score),score_json_sha256=sha(score),mean_plddt='95.0',ptm='0.8')
    ])
    base=['--mapping',mapping,'--queries',query,'--family-audit',audit,'--family-context',context,'--historical-structures',old,'--limit','2','--expect-families','6']
    out=root/'run'
    t=run(*base,'--out-dir',out,'--mode','check')
    assert t.returncode==0, t.stderr
    assert 'STRUCTURE_REUSE_CANDIDATES=1' in t.stdout
    assert not out.exists()
    print('STAGE179_CHECK_MODE_NO_WRITE=PASS')
    t=run(*base,'--out-dir',out,'--mode','run')
    assert t.returncode==0, t.stderr
    with (out/'partial_structure_inventory.tsv').open() as fh:
        rows=list(csv.DictReader(fh,delimiter='\t'))
    assert len(rows)==6
    assert sum(r['structure_status']=='VERIFIED_HISTORICAL_REUSE_CANDIDATE' for r in rows)==1
    with (out/'new_prediction_review_shortlist.tsv').open() as fh:
        top=list(csv.DictReader(fh,delimiter='\t'))
    assert [r['structure_id'] for r in top]==['MVF000002','MVF000003'],[r['structure_id'] for r in top]
    assert len((out/'new_prediction_review_queries.faa').read_text().split('>'))-1==2
    assert json.loads((out/'partial_structure_plan_provenance.json').read_text())['gpu_prediction_performed'] is False
    assert run(*base,'--out-dir',out,'--mode','run').returncode!=0
    print('STAGE179_CONSISTENT_PARTIAL_INVENTORY_AND_SHORTLIST=PASS')
    print('STAGE179_OVERWRITE_REJECTION=PASS')
    pdb.write_text('TAMPER\n')
    t=run(*base,'--out-dir',root/'bad','--mode','check')
    assert t.returncode!=0 and 'OLD_SOURCE_SHA_MISMATCH' in t.stderr,t.stderr
    print('STAGE179_HISTORICAL_PDB_TAMPER_REJECTED=PASS')
    pdb.write_text('MODEL   1\nENDMDL\n')
    # Source-less new dataset: allow zero historical reuse.
    args_no_hist=['--mapping',mapping,'--queries',query,'--family-audit',audit,'--family-context',context,'--limit','2','--expect-families','6']
    t=run(*args_no_hist,'--out-dir',root/'nohist','--mode','check')
    assert t.returncode==0 and 'STRUCTURE_REUSE_CANDIDATES=0' in t.stdout,t.stderr
    print('STAGE179_GENERAL_NEW_DATASET_NO_HISTORICAL_REQUIRED=PASS')
    query.write_text(query.read_text().replace(sequences['MVF000001'], 'AAAAAAAAAAAAAAAA'))
    t=run(*base,'--out-dir',root/'bad2','--mode','check')
    assert t.returncode!=0 and 'QUERY_' in t.stderr,t.stderr
    print('STAGE179_QUERY_SEQUENCE_MISMATCH_REJECTED=PASS')
print('ALL_TESTS_PASS')
