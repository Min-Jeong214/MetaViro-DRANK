#!/usr/bin/env python3
"""Stage 11/12 structure-optional MSA selection and sequence-only conservation.

No structure prediction is performed. Full Stage12 (PDB + pLDDT and spatial
residue analysis) is *not* represented by the sequence-only result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

DARK_TIERS = {
    'T1_CONTEXT_AND_CONVINCING_STRUCTURE',
    'T2_CONTEXT_SUPPORTED_DARK',
    'T3_CONVINCING_STRUCTURE_NO_SUPPORTED_CONTEXT',
}
REVIEW_COLUMNS = ('structure_id','family_rep','family_size','review_order',
                  'review_tier','priority_class','stage09_evidence_state')
MAP_COLUMNS = ('structure_id','family_rep','family_size','sequence_sha256')
MEMBER_COLUMNS = ('member','representative')
MSA_COLUMNS = ('structure_id','family_rep','family_size','unique_sequences',
               'collapsed_duplicates','msa_status','query_length_aa',
               'a3m_sha256','mafft_version')


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()


def table(path, required):
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f'MISSING_TABLE:{path}')
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        reader=csv.DictReader(f, delimiter='\t')
        cols=reader.fieldnames or []
        if len(cols)!=len(set(cols)) or not set(required).issubset(cols):
            raise ValueError(f'INVALID_COLUMNS:{path}:{sorted(set(required)-set(cols))}')
        rows=list(reader)
    if not rows or any(None in r or None in r.values() for r in rows):
        raise ValueError(f'INVALID_TABLE_ROWS:{path}')
    return rows, cols


def uniq(rows,key):
    obj={}
    for r in rows:
        name=r[key]
        if not name or name in obj:
            raise ValueError(f'DUPLICATE_OR_EMPTY:{key}:{name}')
        obj[name]=r
    return obj


def fasta(path):
    if not path.is_file():
        raise ValueError(f'MISSING_FASTA:{path}')
    data={}; key=None; chunks=[]
    def flush():
        if key is not None:
            seq=''.join(chunks).upper().rstrip('*')
            if not seq or not seq.isalpha() or key in data:
                raise ValueError(f'INVALID_FASTA_RECORD:{key}')
            data[key]=seq
    with path.open(encoding='utf-8') as f:
        for raw in f:
            line=raw.strip()
            if not line:continue
            if line.startswith('>'):
                flush(); key=line[1:].split()[0]; chunks=[]
                if not key:raise ValueError('EMPTY_FASTA_ID')
            else:
                if key is None:raise ValueError('FASTA_MISSING_HEADER')
                chunks.append(line)
    flush()
    if not data:raise ValueError('EMPTY_FASTA')
    return data


def aligned_fasta(path, *, require_equal_width=True):
    records=[]; sid=None; chunks=[]
    if not path.is_file():raise ValueError(f'MISSING_MSA:{path}')
    def flush():
        if sid is not None:
            seq=''.join(chunks).upper()
            if not seq or not set(seq)<=set('ABCDEFGHIJKLMNOPQRSTUVWXYZ-'):
                raise ValueError(f'INVALID_MSA_SEQUENCE:{sid}')
            records.append((sid,seq))
    with path.open() as f:
        for raw in f:
            line=raw.strip()
            if not line:continue
            if line.startswith('>'):
                flush();sid=line[1:].split()[0];chunks=[]
            else:
                if sid is None:raise ValueError('MSA_MISSING_HEADER')
                chunks.append(line)
    flush()
    if not records or len(set(k for k,_ in records))!=len(records):
        raise ValueError('MISSING_OR_DUPLICATED_ALIGNMENT_RECORDS')
    if require_equal_width and len(set(len(v) for _,v in records))!=1:
        raise ValueError('UNEQUAL_MSA_WIDTH')
    return records


def dump_table(path,fields,rows):
    with path.open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def dump_fasta(path,records):
    with path.open('x',encoding='utf-8') as f:
        for sid,seq in records:
            f.write(f'>{sid}\n')
            for start in range(0,len(seq),80):f.write(seq[start:start+80]+'\n')


def jload(path):
    if not path.is_file():raise ValueError(f'MISSING_PROVENANCE:{path}')
    obj=json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(obj,dict):raise ValueError('INVALID_PROVENANCE')
    return obj


def refuse_overlap(out, sources):
    if out.exists():raise ValueError(f'REFUSE_EXISTING_OUTPUT:{out}')
    target=out.resolve()
    for path in sources:
        source=path.resolve()
        if target==source or source in target.parents or target in source.parents:
            raise ValueError(f'OUTPUT_OVERLAPS_SOURCE:{path}')


def stage183_valid(review_dir):
    reviewfile=review_dir/'all_family_stage10_review.tsv'
    manifestfile=review_dir/'stage183_review_provenance.json'
    meta=jload(manifestfile)
    if meta.get('stage')!='STAGE10_PARTIAL_EVIDENCE_REVIEW' or meta.get('stage10_historical_I2_completed') is not False:
        raise ValueError('NOT_PRELIMINARY_STAGE183')
    if meta.get('output_sha256',{}).get(reviewfile.name)!=sha(reviewfile):
        raise ValueError('STAGE183_REVIEW_SHA_MISMATCH')
    rows,_=table(reviewfile,REVIEW_COLUMNS)
    ids=uniq(rows,'structure_id')
    if len(ids)!=meta.get('all_families'):
        raise ValueError('STAGE183_FAMILY_COUNT_MISMATCH')
    for i,row in enumerate(rows,1):
        if row['review_order']!=str(i) or row['review_tier'] not in {
            'T1_CONTEXT_AND_CONVINCING_STRUCTURE','T2_CONTEXT_SUPPORTED_DARK',
            'T3_CONVINCING_STRUCTURE_NO_SUPPORTED_CONTEXT',
            'T4_DARK_CONTEXT_UNRESOLVED','T5_KNOWN_EVIDENCE_REFERENCE'}:
            raise ValueError('INVALID_REVIEW_ORDER_OR_TIER')
    return rows,reviewfile,manifestfile


def prepare(args):
    if args.limit<1:raise ValueError('LIMIT_MUST_BE_POSITIVE')
    sources=(args.review_dir,args.structure_map,args.members,args.fasta)
    refuse_overlap(args.out_dir,sources)
    reviews,reviewfile,reviewmeta=stage183_valid(args.review_dir)
    mappings,mapfields=table(args.structure_map,MAP_COLUMNS)
    memberrows,memfields=table(args.members,MEMBER_COLUMNS)
    seq=fasta(args.fasta)
    m=uniq(mappings,'structure_id')
    if len(m)!=len(reviews) or set(m)!={r['structure_id'] for r in reviews}:
        raise ValueError('STAGE08_STAGE183_FAMILY_SET_MISMATCH')
    group=defaultdict(list);member_ids=set()
    for r in memberrows:
        member,rep=r['member'],r['representative']
        if not member or not rep or member in member_ids:
            raise ValueError('DUPLICATE_OR_EMPTY_MEMBER')
        member_ids.add(member);group[rep].append(member)
    if member_ids!=set(seq):raise ValueError('MEMBER_FASTA_DIFFERENCE')
    for r in reviews:
        mr=m[r['structure_id']]
        if r['family_rep']!=mr['family_rep'] or r['family_size']!=mr['family_size']:
            raise ValueError('STAGE08_STAGE183_FAMILY_ID_SIZE_MISMATCH')
        rep=mr['family_rep']
        if rep not in group or rep not in group[rep] or len(group[rep])!=int(r['family_size']):
            raise ValueError('FAMILY_MEMBERS_MISMATCH')
        if hashlib.sha256(seq[rep].encode('ascii')).hexdigest()!=mr['sequence_sha256']:
            raise ValueError('REPRESENTATIVE_SEQUENCE_SHA_MISMATCH')
    if set(group)!={r['family_rep'] for r in reviews} or len({r['family_rep'] for r in reviews})!=len(reviews):
        raise ValueError('FAMILY_REP_SET_MISMATCH')
    chosen=[r for r in reviews if r['review_tier'] in DARK_TIERS][:args.limit]
    if not chosen:raise ValueError('NO_SUPPORTED_DARK_FAMILIES')
    selection={r['structure_id'] for r in chosen};reps={r['family_rep'] for r in chosen}
    selected_members=[r for r in memberrows if r['representative'] in reps]
    selected_faa=[(r['member'],seq[r['member']]) for r in selected_members]
    inventory=[{'structure_id':r['structure_id'], 'family_rep':r['family_rep'],
                'review_order':r['review_order'],'review_tier':r['review_tier'],
                'priority_class':r['priority_class'],
                'structure_state':r.get('structure_state',''),
                'msa_plan_state':'SELECTED_FOR_SEQUENCE_MSA' if r['structure_id'] in selection else 'NOT_SELECTED_FOR_SEQUENCE_MSA',
                'stage12_structure_state':'NOT_EVALUATED_BY_STAGE184'} for r in reviews]
    print(f'ALL_FAMILIES={len(reviews)}')
    print(f'SELECTED_MSA_FAMILIES={len(chosen)}')
    print(f'SELECTED_PROTEINS={len(selected_members)}')
    print('FULL_STAGE12_SPATIAL=NOT_RUN')
    if args.mode=='check':
        print('STEP184_PREPARE_CHECK=PASS');print('FILES_WRITTEN=NO');return
    args.out_dir.mkdir(parents=True,exist_ok=False)
    dump_table(args.out_dir/'selected_structure_map.tsv',mapfields,[m[r['structure_id']] for r in chosen])
    dump_table(args.out_dir/'selected_members.tsv',memfields,selected_members)
    dump_fasta(args.out_dir/'selected_proteins.faa',selected_faa)
    dump_table(args.out_dir/'all_family_msa_plan.tsv',list(inventory[0]),inventory)
    manifest={'stage':'STAGE11_SELECTED_MSA_PREPARATION','schema_version':1,
              'status':'SELECTED_ONLY_NOT_ALL_FAMILIES','all_families':len(reviews),
              'selected_structure_ids':[r['structure_id'] for r in chosen],
              'selected_protein_count':len(selected_members),
              'source_sha256':{str(p.resolve()):sha(p) for p in (reviewfile,reviewmeta,args.structure_map,args.members,args.fasta)},
              'output_sha256':{p.name:sha(p) for p in sorted(args.out_dir.glob('*')) if p.is_file()},
              'not_full_stage12':True,'gpu_prediction_run':False}
    (args.out_dir/'stage184_prepare_provenance.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    print('STEP184_PREPARE=PASS');print('OUTPUT_DIR='+str(args.out_dir))


def collect(args):
    if not (0<args.identity<=1) or args.min_homologs<1:
        raise ValueError('BAD_CONSERVATION_THRESHOLDS')
    refuse_overlap(args.out_dir,(args.prepare_dir,args.msa_dir))
    prep=args.prepare_dir;msa=args.msa_dir
    meta_file=prep/'stage184_prepare_provenance.json';meta=jload(meta_file)
    if meta.get('stage')!='STAGE11_SELECTED_MSA_PREPARATION' or meta.get('not_full_stage12') is not True:
        raise ValueError('NOT_STAGE184_PREPARATION')
    for name,digest in meta['output_sha256'].items():
        if sha(prep/name)!=digest:raise ValueError('STAGE184_PREPARATION_TAMPER')
    for fname,digest in meta['source_sha256'].items():
        if sha(Path(fname))!=digest:raise ValueError('UPSTREAM_SOURCE_TAMPER')
    mapping,_=table(prep/'selected_structure_map.tsv',MAP_COLUMNS)
    plans,_=table(prep/'all_family_msa_plan.tsv',('structure_id','msa_plan_state','family_rep'))
    members,_=table(prep/'selected_members.tsv',MEMBER_COLUMNS)
    seq=fasta(prep/'selected_proteins.faa')
    summary_file=msa/'family_msa_summary.tsv'
    summary,_=table(summary_file,MSA_COLUMNS)
    mi=uniq(mapping,'structure_id');si=uniq(summary,'structure_id')
    selected=meta['selected_structure_ids']
    if set(mi)!=set(si) or set(mi)!=set(selected) or len(selected)!=len(mi):
        raise ValueError('MSA_SELECTED_STRUCTURE_SET_MISMATCH')
    if len(plans)!=meta['all_families'] or len(uniq(plans,'structure_id'))!=len(plans):
        raise ValueError('MSA_PLAN_TOTAL_FAMILY_MISMATCH')
    if {r['structure_id'] for r in plans if r['msa_plan_state']=='SELECTED_FOR_SEQUENCE_MSA'}!=set(selected):
        raise ValueError('MSA_PLAN_SELECTION_MISMATCH')
    groups=defaultdict(list)
    for r in members:groups[r['representative']].append(r['member'])
    if len({r['member'] for r in members})!=len(members) or set(seq)!={r['member'] for r in members}:
        raise ValueError('SELECTED_MEMBER_FASTA_MISMATCH')
    residues=[];summaries=[];hashes={str(summary_file.resolve()):sha(summary_file)};stats={}
    for sid in selected:
        mr=mi[sid];sr=si[sid];rep=mr['family_rep'];fam=groups[rep]
        if rep not in fam or len(fam)!=int(mr['family_size']) or sr['family_rep']!=rep or sr['family_size']!=mr['family_size']:
            raise ValueError('MSA_FAMILY_MEMBERSHIP_MISMATCH')
        aligned_file=msa/'mafft'/f'{sid}.aligned.faa'
        a3m_file=msa/'a3m_colabfold'/f'{sid}.a3m'
        if sha(a3m_file)!=sr['a3m_sha256']:
            raise ValueError('MAFFT_A3M_SHA_MISMATCH')
        hashes[str(aligned_file.resolve())]=sha(aligned_file)
        hashes[str(a3m_file.resolve())]=sha(a3m_file)
        aln=aligned_fasta(aligned_file)
        # The first ColabFold A3M record must be the exact representative
        # sequence, not simply a self-consistent checksum supplied by MSA TSV.
        a3m_records=aligned_fasta(a3m_file, require_equal_width=False)
        if (len(a3m_records)!=len(aln) or a3m_records[0][0]!=sid or
            a3m_records[0][1]!=seq[rep] or
            len({name for name,_ in a3m_records})!=len(a3m_records)):
            raise ValueError('A3M_REPRESENTATIVE_ID_OR_SEQUENCE_MISMATCH')
        for idx,(name,sequence) in enumerate(a3m_records[1:],1):
            if name!=f'{sid}_h{idx:04d}' or not sequence.replace('-',''):
                raise ValueError('A3M_HOMOLOG_HEADER_OR_SEQUENCE_INVALID')
        if aln[0][0]!=rep or len({k for k,_ in aln})!=len(aln) or len(aln)!=int(sr['unique_sequences']):
            raise ValueError('MSA_ALIGNED_MEMBERS_MISMATCH')
        if aln[0][1].replace('-','')!=seq[rep] or len(seq[rep])!=int(sr['query_length_aa']):
            raise ValueError('MSA_REP_SEQUENCE_MISMATCH')
        if len({v for _,v in aln})!=len(aln):
            raise ValueError('UNEXPECTED_DUPLICATE_ALIGNED_SEQUENCE')
        for member,align in aln:
            if member not in fam or align.replace('-','')!=seq[member]:
                raise ValueError('MSA_ALIGNMENT_OR_SEQUENCE_TAMPER')
        if {seq[k] for k in fam}!={v.replace('-','') for _,v in aln}:
            raise ValueError('MSA_UNIQUE_SEQUENCE_CONTENT_MISMATCH')
        if {v.replace('-','') for _,v in a3m_records}!={seq[k] for k in fam}:
            raise ValueError('A3M_UNIQUE_SEQUENCE_CONTENT_MISMATCH')
        if int(sr['collapsed_duplicates'])!=len(fam)-len(aln):
            raise ValueError('MSA_COLLAPSED_DUPLICATE_COUNT_MISMATCH')
        ref=aln[0][1];homologs=[s for _,s in aln[1:]]
        n_supported=0;n_conserved=0;querypos=0
        for col,aa in enumerate(ref):
            if aa=='-':continue
            querypos+=1
            column=[s[col] for s in homologs if s[col]!='-']
            count=len(column);match=sum(v==aa for v in column)
            reliable=count>=args.min_homologs
            fraction=match/count if count else None
            conserved=reliable and fraction>=args.identity
            if reliable:n_supported+=1
            if conserved:n_conserved+=1
            residues.append({'structure_id':sid,'query_position':querypos,'query_residue':aa,
                             'non_gap_homologs':count,'matching_homologs':match,
                             'identity_fraction':f'{fraction:.6f}' if fraction is not None else '',
                             'minimum_homolog_support_met':'YES' if reliable else 'NO',
                             'sequence_conserved':'YES' if conserved else 'NO'})
        if querypos!=len(seq[rep]):raise ValueError('QUERY_RESIDUE_COORDINATE_MISMATCH')
        state=('SINGLE_UNIQUE_SEQUENCE' if len(aln)==1 else
               'INSUFFICIENT_HOMOLOG_SUPPORT' if n_supported==0 else
               'SEQUENCE_ONLY_CONSERVATION_AVAILABLE')
        summaries.append({'structure_id':sid,'family_rep':rep,'family_size':len(fam),
                          'msa_unique_sequences':len(aln),'query_length':querypos,
                          'min_homologs_per_position':args.min_homologs,'minimum_identity_fraction':f'{args.identity:.4f}',
                          'analyzable_positions':n_supported,'conserved_positions':n_conserved,
                          'msa_conservation_state':state,'structure_residue_analysis':'NOT_EVALUATED_NO_MATCHED_PDB_PREDICTION_RECORD'})
        stats[sid]=summaries[-1]
    out=[]
    for r in plans:
        sid=r['structure_id'];s=stats.get(sid)
        out.append({**r,'msa_conservation_state':s['msa_conservation_state'] if s else 'NOT_SELECTED_FOR_MSA',
                    'analyzable_positions':s['analyzable_positions'] if s else '',
                    'conserved_positions':s['conserved_positions'] if s else '',
                    'full_stage12_structure_analysis':'NOT_EVALUATED'})
    print('ALL_FAMILIES='+str(len(out)))
    print('MAFFT_FAMILIES='+str(len(summaries)))
    print('SEQUENCE_ONLY_CONSERVATION_AVAILABLE='+str(sum(x['msa_conservation_state']=='SEQUENCE_ONLY_CONSERVATION_AVAILABLE' for x in summaries)))
    print('STAGE12_FULL_STRUCTURE_ANALYSIS=NOT_RUN')
    if args.mode=='check':
        print('STEP184_COLLECT_CHECK=PASS');print('FILES_WRITTEN=NO');return
    args.out_dir.mkdir(parents=True,exist_ok=False)
    dump_table(args.out_dir/'selected_family_msa_conservation.tsv',list(summaries[0]),summaries)
    dump_table(args.out_dir/'selected_query_residue_conservation.tsv',list(residues[0]),residues)
    dump_table(args.out_dir/'all_family_stage11_12_status.tsv',list(out[0]),out)
    result={'stage':'STAGE11_STAGE12_SEQUENCE_ONLY_PARTIAL','schema_version':1,
            'all_families':len(out),'msa_families':len(summaries),
            'structure_based_stage12_completed':False,
            'biochemical_function_confirmed':False,
            'thresholds':{'min_homologs':args.min_homologs,'identity_fraction':args.identity},
            'source_sha256':{**hashes,str(meta_file.resolve()):sha(meta_file)},
            'output_sha256':{p.name:sha(p) for p in sorted(args.out_dir.glob('*.tsv'))}}
    (args.out_dir/'stage184_conservation_provenance.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('STEP184_SEQUENCE_ONLY_CONSERVATION=PASS')
    print('OUTPUT_DIR='+str(args.out_dir))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    first=sub.add_parser('prepare',help='Select Family subset and build portable Stage11 inputs')
    first.add_argument('--review-dir',type=Path,required=True)
    first.add_argument('--structure-map',type=Path,required=True)
    first.add_argument('--members',type=Path,required=True)
    first.add_argument('--fasta',type=Path,required=True)
    first.add_argument('--limit',type=int,default=10)
    first.add_argument('--out-dir',type=Path,required=True)
    first.add_argument('--mode',choices=('check','run'),default='check')
    second=sub.add_parser('collect',help='Validate Stage11 MSA and generate sequence-only conservation')
    second.add_argument('--prepare-dir',type=Path,required=True)
    second.add_argument('--msa-dir',type=Path,required=True)
    second.add_argument('--min-homologs',type=int,default=2)
    second.add_argument('--identity',type=float,default=0.8)
    second.add_argument('--out-dir',type=Path,required=True)
    second.add_argument('--mode',choices=('check','run'),default='check')
    args=p.parse_args()
    if args.command=='prepare':prepare(args)
    else:collect(args)


if __name__=='__main__':main()
