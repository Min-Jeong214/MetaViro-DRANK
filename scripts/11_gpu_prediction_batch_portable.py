#!/usr/bin/env python3
"""Sequential, resumable ColabFold batch for all cohort family representatives.

Run from an independently installed ColabFold environment (tested: Python 3.11,
ColabFold 1.6.1, JAX 0.5.3). No unverified structure is marked COMPLETED.
"""
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def tsv(path, cols):
    with path.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        names=reader.fieldnames or []
        if len(names)!=len(set(names)) or not set(cols)<=set(names):
            raise ValueError('BAD_TSV_COLUMNS:'+str(path))
        rows=list(reader)
    if not rows or any(None in r or None in r.values() for r in rows):
        raise ValueError('BAD_TSV:'+str(path))
    return rows


def fasta(path):
    data={}; sid=None
    for raw in path.read_text().splitlines():
        x=raw.strip()
        if not x:continue
        if x.startswith('>'):
            sid=x[1:].split()[0]
            if not sid or sid in data:raise ValueError('DUPLICATE_FASTA_ID:'+str(sid))
            data[sid]=''
        else:
            if sid is None:raise ValueError('FASTA_MISSING_HEADER')
            data[sid]+=x
    return {k:v.upper().rstrip('*') for k,v in data.items()}


def first_a3m(path):
    entries=0; sid=None; seq=''; started=False
    for raw in path.read_text().splitlines():
        x=raw.strip()
        if not x: continue
        if x.startswith('>'):
            entries+=1
            if not started:
                sid=x[1:].split()[0];started=True
        else:
            if not started: raise ValueError('A3M_MISSING_HEADER:'+str(path))
            if entries==1: seq+=x
    return sid,seq,entries


def completed(path, sid, a3m_hash, query_hash):
    record=path/'prediction_provenance.json'
    if not record.is_file():return False
    try:
        m=json.loads(record.read_text())
        if m.get('status')!='COMPLETED' or m.get('structure_id')!=sid:
            return False
        if m.get('a3m_sha256')!=a3m_hash or m.get('query_sha256')!=query_hash:
            return False
        for k,h in (('pdb_file','pdb_sha256'),('score_json','score_json_sha256')):
            p=Path(m[k]);
            if not p.is_file() or not p.stat().st_size or sha(p)!=m[h]:
                return False
        return True
    except (KeyError,ValueError,OSError,TypeError):
        return False


def report(path, rows):
    fields=('structure_id','status','attempt','path')
    tmp=path.with_name(path.name+'.partial.'+str(os.getpid()))
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
    os.replace(tmp,path)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structure-map',type=Path,required=True)
    p.add_argument('--queries',type=Path,required=True)
    p.add_argument('--msa-dir',type=Path,required=True,help='Directory containing a3m_colabfold/ and family_msa_summary.tsv')
    p.add_argument('--out-dir',type=Path,required=True)
    p.add_argument('--reuse-prediction',type=Path,action='append',default=[])
    p.add_argument('--max-jobs',type=int,default=0,help='0: all pending; positive: at most this many per invocation')
    p.add_argument('--max-query-length',type=int,default=0,help='0: no length cutoff')
    p.add_argument('--max-seq',type=int,default=64)
    p.add_argument('--max-extra-seq',type=int,default=128)
    p.add_argument('--retry-failed',action='store_true')
    p.add_argument('--skip-incompatible-reuse',action='store_true',
                   help='Do not reuse an existing prediction if its input/hash integrity differs; schedule a new prediction instead. Default: fail closed.')
    p.add_argument('--mode',choices=('check','run'),default='check')
    a=p.parse_args(argv)
    if a.max_jobs<0 or a.max_query_length<0 or a.max_seq<2 or a.max_extra_seq<1:
        raise ValueError('INVALID_BATCH_OPTIONS')
    inputs=(a.structure_map.resolve(),a.queries.resolve(),(a.msa_dir/'family_msa_summary.tsv').resolve())
    for file in inputs:
        if not file.is_file():raise ValueError('MISSING_INPUT:'+str(file))
    for file in inputs:
        if a.out_dir.resolve() in file.parents or file in a.out_dir.resolve().parents:
            raise ValueError('INPUT_OUTPUT_OVERLAP')
    m=tsv(a.structure_map,('structure_id','family_rep','sequence_sha256','length_aa'))
    summaries=tsv(a.msa_dir/'family_msa_summary.tsv',('structure_id','a3m_sha256','unique_sequences'))
    index={r['structure_id']:r for r in summaries}
    ids=[r['structure_id'] for r in m]
    if len(set(ids))!=len(ids) or len(index)!=len(summaries) or set(ids)!=set(index):
        raise ValueError('MSA_MAPPING_COHORT_MISMATCH')
    queries=fasta(a.queries)
    if set(queries)!=set(ids):raise ValueError('QUERY_MAPPING_COHORT_MISMATCH')
    validated={}
    for r in m:
        sid=r['structure_id'];s=queries[sid]; sm=index[sid]
        if not sid.replace('_','').isalnum() or '/' in sid or '..' in sid:
            raise ValueError('UNSAFE_STRUCTURE_ID:'+sid)
        if len(s)!=int(r['length_aa']) or hashlib.sha256(s.encode('ascii')).hexdigest()!=r['sequence_sha256']:
            raise ValueError('QUERY_MAPPING_SHA_MISMATCH:'+sid)
        path=a.msa_dir/'a3m_colabfold'/(sid+'.a3m')
        if not path.is_file() or sha(path)!=sm['a3m_sha256']:
            raise ValueError('MSA_SHA_MISMATCH:'+sid)
        qid,qseq,count=first_a3m(path)
        if qid!=sid or qseq!=s or count!=int(sm['unique_sequences']):
            raise ValueError('MSA_QUERY_OR_COUNT_MISMATCH:'+sid)
        validated[sid]={'path':path,'sha':sm['a3m_sha256'], 'query_sha':r['sequence_sha256'],
                        'count':count,'length':len(s)}
    # Strict reuse: only matching query, identical MSA and verified PDB/score hashes.
    # Explicit skip means fresh prediction; it does NOT reinterpret old models as
    # predictions from the newly prepared MSA.
    reused={}; skipped_reuse=[];seen_reuse_ids=set()
    for d in a.reuse_prediction:
        pth=d.resolve();file=pth/'prediction_provenance.json'
        if not file.is_file():raise ValueError('REUSE_PROVENANCE_MISSING:'+str(pth))
        meta=json.loads(file.read_text());sid=meta.get('structure_id')
        if sid not in validated:raise ValueError('REUSE_ID_OUT_OF_COHORT:'+str(sid))
        if sid in seen_reuse_ids:raise ValueError('DUPLICATE_REUSE:'+sid)
        seen_reuse_ids.add(sid)
        data=validated[sid]
        if not completed(pth,sid,data['sha'],data['query_sha']):
            why=[]
            if meta.get('status')!='COMPLETED':why.append('status')
            if meta.get('query_sha256')!=data['query_sha']:why.append('query_sha256')
            if meta.get('a3m_sha256')!=data['sha']:why.append('a3m_sha256')
            if not why:why.append('PDB_or_score_integrity')
            reason=','.join(why)
            if a.skip_incompatible_reuse:
                skipped_reuse.append(sid)
                print('REUSE_INCOMPATIBLE_SKIPPED='+sid+':'+reason,flush=True)
                continue
            raise ValueError('REUSE_HASH_OR_PDB_INVALID:'+sid+':'+reason)
        reused[sid]=pth
    rows=[];pending=[];failures=[]
    for sid in sorted(ids):
        d=validated[sid]
        fam=a.out_dir/'families'/sid
        previous=sorted(p for p in fam.glob('attempt_*') if p.is_dir()) if fam.is_dir() else []
        done=next((q for q in previous if completed(q,sid,d['sha'],d['query_sha'])),None)
        if done:
            status='COMPLETED';attempt=done
        elif sid in reused:
            status='REUSED';attempt=reused[sid]
        elif d['count']<2:
            status='SINGLE_UNIQUE_SEQUENCE';attempt=None
        elif a.max_query_length and d['length']>a.max_query_length:
            status='EXCEEDS_LENGTH_GUARD';attempt=None
        elif previous and not a.retry_failed:
            status='NEEDS_EXPLICIT_RETRY';attempt=previous[-1];failures.append(sid)
        else:
            status='PENDING';attempt=None;pending.append(sid)
        rows.append({'structure_id':sid,'status':status,'attempt':str(attempt.name if attempt else ''),
                     'path':str(attempt or '')})
    counts={k:sum(r['status']==k for r in rows) for k in set(r['status'] for r in rows)}
    print('TOTAL_FAMILIES='+str(len(rows)))
    print('BATCH_STATUS='+json.dumps(counts,sort_keys=True))
    print('GPU_PENDING='+str(len(pending)))
    print('GPU_EXISTING_REUSED='+str(len(reused)))
    print('REUSE_INCOMPATIBLE_SKIPPED_COUNT='+str(len(skipped_reuse)))
    if a.mode=='check':
        print('STAGE188_GPU_BATCH_PREFLIGHT=PASS');print('FILES_WRITTEN=NO');return
    if shutil.which('colabfold_batch') is None:
        raise ValueError('COLABFOLD_BATCH_NOT_INSTALLED_IN_ACTIVE_ENV')
    try:
        cf=version('colabfold')
    except PackageNotFoundError:
        cf='UNKNOWN'
    print('RUNNING_COLABFOLD_VERSION='+cf,flush=True)
    print('NOTE_VALIDATED_ENVIRONMENT=ColabFold_1.6.1_JAX_0.5.3_Python_3.11',flush=True)
    wrapper=Path(__file__).resolve().parent/'11_msa_structure_predict_portable.py'
    if not wrapper.is_file():raise ValueError('STRUCTURE_WRAPPER_NOT_FOUND')
    a.out_dir.mkdir(parents=True,exist_ok=True)
    lock=a.out_dir/'.batch.lock'
    with lock.open('a+') as f:
        try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('GPU_BATCH_ALREADY_RUNNING')
        manifest=a.out_dir/'gpu_batch_plan.json'
        plan={'schema':1,'mapping_sha256':sha(a.structure_map),'queries_sha256':sha(a.queries),
              'msa_summary_sha256':sha(a.msa_dir/'family_msa_summary.tsv'),
              'ids':sorted(ids),'max_seq':a.max_seq,'max_extra_seq':a.max_extra_seq}
        if manifest.is_file():
            if json.loads(manifest.read_text())!=plan:
                raise ValueError('GPU_BATCH_INPUT_OR_PARAMETERS_CHANGED')
        else:manifest.write_text(json.dumps(plan,indent=2)+'\n')
        report(a.out_dir/'prediction_manifest.tsv',rows)
        launched=0
        for sid in pending:
            if a.max_jobs and launched>=a.max_jobs:break
            fam=a.out_dir/'families'/sid
            fam.mkdir(parents=True,exist_ok=True)
            attempts=list(fam.glob('attempt_*'))
            attempt=fam/f'attempt_{len(attempts)+1:03d}'
            if attempt.exists():raise ValueError('ATTEMPT_ALREADY_EXISTS:'+str(attempt))
            d=validated[sid]
            cmd=[sys.executable,'-u','-B',str(wrapper),'--a3m',str(d['path']),
                 '--out-dir',str(attempt),'--max-seq',str(a.max_seq),
                 '--max-extra-seq',str(a.max_extra_seq),'--execute']
            log=fam/(attempt.name+'.driver.log')
            env=os.environ.copy();env['XLA_PYTHON_CLIENT_PREALLOCATE']='false';env['PYTHONUNBUFFERED']='1'
            print('GPU_START='+sid,flush=True)
            with log.open('x') as output:
                proc=subprocess.run(cmd,env=env,stdout=output,stderr=subprocess.STDOUT,check=False)
            good=proc.returncode==0 and completed(attempt,sid,d['sha'],d['query_sha'])
            row=next(r for r in rows if r['structure_id']==sid)
            row.update(status='COMPLETED' if good else 'FAILED',attempt=attempt.name,path=str(attempt))
            report(a.out_dir/'prediction_manifest.tsv',rows)
            launched+=1
            print('GPU_RESULT='+sid+':' + ('PASS' if good else 'FAIL')+':exit='+str(proc.returncode),flush=True)
        print('GPU_LAUNCHED_THIS_INVOCATION='+str(launched))
        print('BATCH_RESUME_SAFE=YES')
        print('STAGE188_GPU_BATCH_RUN_FINISHED=YES')


if __name__=='__main__':
    main()
