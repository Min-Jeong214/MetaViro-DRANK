#!/usr/bin/env python3
"""Prepare MSAs for *every* structure-mapped family in restartable chunks.

Uses the validated 11_family_msa_portable.py on deterministic chunks. Does
not select only dark/context-supported families or silently drop singletons.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_tsv(path, cols):
    with path.open(newline='', encoding='utf-8-sig') as f:
        r = csv.DictReader(f, delimiter='\t')
        names = r.fieldnames or []
        if len(names) != len(set(names)) or not set(cols) <= set(names):
            raise ValueError('INVALID_TABLE_COLUMNS:' + str(path))
        rows = list(r)
    if not rows or any(None in r or None in r.values() for r in rows):
        raise ValueError('INVALID_TABLE:' + str(path))
    return names, rows


def read_fasta(path):
    result = {}
    sid = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            sid = line[1:].split()[0]
            if sid in result or not sid:
                raise ValueError('DUPLICATE_FASTA:' + str(sid))
            result[sid] = ''
        else:
            if sid is None:
                raise ValueError('FASTA_HEADER_MISSING')
            result[sid] += line
    return {k:v.upper().rstrip('*') for k,v in result.items()}


def put_table(path, cols, rows):
    with path.open('x', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def put_fasta(path, records):
    with path.open('x') as f:
        for identifier, seq in records:
            f.write('>' + identifier + '\n' + seq + '\n')


def copy_verified(src, dst):
    if not src.is_file():
        raise ValueError('MISSING_CHUNK_FILE:' + str(src))
    if dst.exists():
        if sha(src) != sha(dst):
            raise ValueError('CONFLICTING_AGGREGATE_FILE:' + str(dst))
        return
    tmp = dst.with_name(dst.name + '.partial.' + str(os.getpid()))
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--members', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--structure-map', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    p.add_argument('--chunk-size', type=int, default=10)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--mode', choices=('check','run'), default='check')
    a = p.parse_args(argv)
    if a.chunk_size < 1 or a.threads < 1:
        raise ValueError('INVALID_CHUNK_SIZE_OR_THREADS')
    paths = [a.members.resolve(), a.fasta.resolve(), a.structure_map.resolve()]
    for source in paths:
        if not source.is_file() or not source.stat().st_size:
            raise ValueError('MISSING_INPUT:' + str(source))
        if source == a.out_dir.resolve() or a.out_dir.resolve() in source.parents or source in a.out_dir.resolve().parents:
            raise ValueError('OUTPUT_INPUT_OVERLAP')
    mapcols, mapping = read_tsv(a.structure_map, ('structure_id','family_rep','family_size','sequence_sha256'))
    membercols, members = read_tsv(a.members, ('member','representative'))
    proteins = read_fasta(a.fasta)
    ids = [r['structure_id'] for r in mapping]
    reps = [r['family_rep'] for r in mapping]
    if len(set(ids)) != len(ids) or len(set(reps)) != len(reps):
        raise ValueError('DUPLICATE_STRUCTURE_OR_FAMILY_ID')
    byrep = {}
    seen = set()
    for r in members:
        rep, mid = r['representative'], r['member']
        if mid in seen:
            raise ValueError('DUPLICATE_MEMBER:' + mid)
        seen.add(mid)
        byrep.setdefault(rep, []).append(mid)
    if set(proteins) != seen or set(byrep) != set(reps):
        raise ValueError('COHORT_MEMBERSHIP_MISMATCH')
    for r in mapping:
        rep = r['family_rep']
        if rep not in byrep[rep] or len(byrep[rep]) != int(r['family_size']):
            raise ValueError('FAMILY_SIZE_MISMATCH:' + rep)
        if hashlib.sha256(proteins[rep].encode('ascii')).hexdigest() != r['sequence_sha256']:
            raise ValueError('QUERY_SHA_MISMATCH:' + rep)
    mapping.sort(key=lambda r:r['structure_id'])
    inputs = {str(p):sha(p) for p in paths}
    payload = dict(schema=1, inputs=inputs, chunk_size=a.chunk_size,
                   total_families=len(mapping), ordered_ids=[r['structure_id'] for r in mapping])
    manifest = a.out_dir / 'full_cohort_msa_plan.json'
    print('FULL_COHORT_FAMILIES=' + str(len(mapping)))
    print('MSA_CHUNKS=' + str((len(mapping)+a.chunk_size-1)//a.chunk_size))
    if a.mode == 'check':
        if manifest.exists() and json.loads(manifest.read_text()) != payload:
            raise ValueError('EXISTING_MSA_PLAN_MISMATCH')
        print('FULL_COHORT_MSA_PREFLIGHT=PASS')
        print('FILES_WRITTEN=NO')
        return
    if shutil.which('mafft') is None:
        raise ValueError('MAFFT_NOT_FOUND')
    script = Path(__file__).resolve().parent / '11_family_msa_portable.py'
    if not script.is_file():
        raise ValueError('FAMILY_MSA_SCRIPT_MISSING')
    a.out_dir.mkdir(parents=True, exist_ok=True)
    if manifest.exists():
        if json.loads(manifest.read_text()) != payload:
            raise ValueError('EXISTING_MSA_PLAN_MISMATCH')
    else:
        manifest.write_text(json.dumps(payload, indent=2)+'\n')
    for section in range(0,len(mapping),a.chunk_size):
        selected = mapping[section:section+a.chunk_size]
        chunk = a.out_dir / 'chunks' / f'chunk_{section//a.chunk_size+1:04d}'
        chunk.mkdir(parents=True, exist_ok=True)
        selection = set(r['family_rep'] for r in selected)
        filtered = [r for r in members if r['representative'] in selection]
        for name, produce in [
            ('structure_map.tsv', lambda x:put_table(x,mapcols,selected)),
            ('members.tsv', lambda x:put_table(x,membercols,filtered)),
            ('proteins.faa', lambda x:put_fasta(x,[(r['member'],proteins[r['member']]) for r in filtered])),
        ]:
            dest = chunk / name
            if not dest.exists():
                produce(dest)
        chunk_out = chunk / 'msa'
        summary = chunk_out / 'family_msa_summary.tsv'
        if not summary.exists():
            if chunk_out.exists():
                raise ValueError('PARTIAL_CHUNK_REQUIRES_INSPECTION:' + str(chunk_out))
            cmd = [sys.executable,'-B',str(script),'--members',str(chunk/'members.tsv'),
                   '--fasta',str(chunk/'proteins.faa'),'--structure-map',str(chunk/'structure_map.tsv'),
                   '--out-dir',str(chunk_out),'--threads',str(a.threads)]
            with (chunk/'run.log').open('a') as log:
                proc = subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=False)
            if proc.returncode:
                raise RuntimeError('CHUNK_MSA_FAILED:'+str(chunk)+':'+str(proc.returncode))
        _, rows = read_tsv(summary, ('structure_id','msa_status','a3m_sha256'))
        if set(r['structure_id'] for r in rows) != set(r['structure_id'] for r in selected):
            raise ValueError('CHUNK_ID_MISMATCH:'+str(chunk))
        for r in rows:
            a3m = chunk_out / 'a3m_colabfold' / (r['structure_id']+'.a3m')
            if not a3m.is_file() or sha(a3m)!=r['a3m_sha256']:
                raise ValueError('CHUNK_A3M_INTEGRITY_ERROR:'+str(a3m))
        print(f'CHUNK_MSA_PASS={section//a.chunk_size+1}',flush=True)
    # Aggregate only once all chunks have validated summaries.
    agg = {name:a.out_dir/name for name in ('a3m_colabfold','mafft','a3m','family_fasta')}
    for folder in agg.values():
        folder.mkdir(exist_ok=True)
    summed=[]; head=[]
    for i in range((len(mapping)+a.chunk_size-1)//a.chunk_size):
        src = a.out_dir/'chunks'/f'chunk_{i+1:04d}'/'msa'
        _,rows=read_tsv(src/'family_msa_summary.tsv',('structure_id','a3m_sha256'))
        _,headers=read_tsv(src/'a3m_colabfold'/'header_map.tsv',('structure_id','new_header','original_header'))
        summed.extend(rows); head.extend(headers)
        for name in agg:
            for file in (src/name).iterdir():
                if file.name == 'header_map.tsv':
                    continue
                copy_verified(file,agg[name]/file.name)
    if len(summed)!=len(mapping) or len(set(r['structure_id'] for r in summed))!=len(mapping):
        raise ValueError('AGGREGATED_FAMILY_COUNT_MISMATCH')
    sumfields=list(summed[0])
    for p,cols,rows in ((a.out_dir/'family_msa_summary.tsv',sumfields,summed),
                        (a.out_dir/'a3m_colabfold'/'header_map.tsv',list(head[0]),head)):
        if p.exists():
            _,current=read_tsv(p,cols)
            if current!=rows:
                raise ValueError('AGGREGATE_TABLE_MISMATCH:'+str(p))
        else:
            put_table(p,cols,rows)
    print('ALL_MSA_FAMILIES='+str(len(summed)))
    print('MSA_MULTIPLE_UNIQUE='+str(sum(int(r['unique_sequences'])>=2 for r in summed)))
    print('MSA_SINGLE_UNIQUE='+str(sum(int(r['unique_sequences'])==1 for r in summed)))
    print('STAGE188_FULL_COHORT_MSA=PASS')


if __name__=='__main__':
    main()
