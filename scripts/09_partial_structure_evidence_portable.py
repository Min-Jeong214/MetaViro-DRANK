#!/usr/bin/env python3
"""Stage09 structure-optional adapter for ANY marine viral cohort.

prepare: build a verified PDB subset from Stage179 inventory and/or new Stage08
results; preserve every family as searched-ready or missing-prediction. It does
not search or assign biochemical activity.
merge: audit the three genuine Foldseek searches with the existing Stage09
integration implementation, then combine searched and pending families without
conflating absence of a structure with absence of a structural match.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

STAGE08_COLUMNS = [
    'structure_id','family_rep','source_member','family_size','priority_class',
    'length_aa','sequence_sha256','pdb_file','pdb_sha256','score_json',
    'score_json_sha256','mean_plddt','min_plddt','max_plddt','ptm',
    'confidence_class','prediction_status',
]
INVENTORY_KEYS = ('structure_id','family_rep','priority_class','family_size',
                  'length_aa','sequence_sha256','structure_status')
STATE_KEYS = ('structure_id','family_rep','priority_class','family_size',
              'length_aa','sequence_sha256','structure_state','structure_origin',
              'pdb_sha256','external_search_state')


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for chunk in iter(lambda:fh.read(1048576),b''):
            h.update(chunk)
    return h.hexdigest()


def table(path: Path, required=(), allow_empty=False):
    if not path.is_file():
        raise ValueError(f'MISSING_TABLE:{path}')
    with path.open(encoding='utf-8-sig',newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)) or not set(required)<=set(fields):
            raise ValueError(f'INVALID_TABLE_SCHEMA:{path}')
        rows=list(reader)
    if not rows and not allow_empty:
        raise ValueError(f'EMPTY_TABLE:{path}')
    return fields,rows


def unique(rows,key):
    result={}
    for row in rows:
        value=row.get(key,'')
        if not value or value in result:
            raise ValueError(f'EMPTY_OR_DUPLICATE_{key}:{value}')
        result[value]=row
    return result


def fasta(path):
    if not path.is_file():
        raise ValueError(f'MISSING_FASTA:{path}')
    values={};ident=None;chunks=[]
    def commit():
        if ident is None:return
        seq=''.join(chunks).upper().rstrip('*')
        if not seq or not seq.isalpha() or not seq.isascii() or ident in values:
            raise ValueError(f'INVALID_FASTA_RECORD:{ident}')
        values[ident]=seq
    with path.open(encoding='utf-8') as fh:
        for raw in fh:
            line=raw.strip()
            if not line:continue
            if line.startswith('>'):
                commit();ident=line[1:].split()[0];chunks=[]
                if not ident:raise ValueError('EMPTY_FASTA_ID')
            else:
                if ident is None:raise ValueError('FASTA_SEQUENCE_BEFORE_HEADER')
                chunks.append(line)
    commit()
    return values


def write_tsv(path, fields, rows):
    with path.open('x',encoding='utf-8',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def checked_model(sid,model,inventory):
    if model['prediction_status']!='PREDICTED':
        raise ValueError(f'UNPREDICTED_MODEL:{sid}')
    if sid != model['structure_id']:
        raise ValueError(f'MODEL_ID_MISMATCH:{sid}')
    if model['sequence_sha256']!=inventory['sequence_sha256']:
        raise ValueError(f'MODEL_SEQUENCE_SHA_MISMATCH:{sid}')
    if int(model['length_aa'])!=int(inventory['length_aa']):
        raise ValueError(f'MODEL_LENGTH_MISMATCH:{sid}')
    if model['priority_class']!=inventory['priority_class']:
        raise ValueError(f'MODEL_CLASS_MISMATCH:{sid}')
    if model['family_rep']!=inventory['family_rep']:
        raise ValueError(f'MODEL_FAMILY_MISMATCH:{sid}')
    if int(model['family_size'])!=int(inventory['family_size']):
        raise ValueError(f'MODEL_SIZE_MISMATCH:{sid}')
    if model['source_member']!=inventory['family_rep']:
        raise ValueError(f'MODEL_REP_MISMATCH:{sid}')
    for file_key,hash_key in (('pdb_file','pdb_sha256'),('score_json','score_json_sha256')):
        filename=model.get(file_key,'')
        if not filename:
            if file_key=='score_json':
                if model.get(hash_key,''):raise ValueError(f'ORPHAN_SCORE_SHA:{sid}')
                continue
            raise ValueError(f'MISSING_PDB_PATH:{sid}')
        file=Path(filename)
        if not file.is_file() or not file.stat().st_size or sha(file)!=model.get(hash_key):
            raise ValueError(f'SOURCE_SHA256_MISMATCH:{sid}:{file_key}')
    # The existing Stage09 integration validates PDB confidence/CA coordinates.
    return dict(model)


def load_inventory(path):
    _,rows=table(path,INVENTORY_KEYS)
    return rows,unique(rows,'structure_id')


def check_plan_provenance(plan_dir,inventory):
    p=plan_dir/'partial_structure_plan_provenance.json'
    if not p.is_file():raise ValueError('STAGE179_PROVENANCE_MISSING')
    meta=json.loads(p.read_text(encoding='utf-8'))
    if meta.get('output_sha256',{}).get('partial_structure_inventory.tsv')!=sha(inventory):
        raise ValueError('STAGE179_INVENTORY_PROVENANCE_MISMATCH')
    return p


def build_preparation(args):
    inventory_file=args.plan_dir/'partial_structure_inventory.tsv'
    plan_prov=check_plan_provenance(args.plan_dir,inventory_file)
    inventory,by_id=load_inventory(inventory_file)
    sequences=fasta(args.queries)
    if set(sequences)!=set(by_id):raise ValueError('FASTA_INVENTORY_ID_MISMATCH')
    for sid,r in by_id.items():
        seq=sequences[sid]
        if len(seq)!=int(r['length_aa']) or hashlib.sha256(seq.encode('ascii')).hexdigest()!=r['sequence_sha256']:
            raise ValueError(f'QUERY_SEQUENCE_SHA_MISMATCH:{sid}')
    new={}
    if args.predicted_structures:
        _,rows=table(args.predicted_structures,STAGE08_COLUMNS)
        new=unique(rows,'structure_id')
        if not set(new)<=set(by_id):raise ValueError('UNEXPECTED_PREDICTED_IDS')
    models=[];states=[]
    for r in inventory:
        sid=r['structure_id']; reuse=r['structure_status']=='VERIFIED_HISTORICAL_REUSE_CANDIDATE'
        if r['structure_status'] not in ('VERIFIED_HISTORICAL_REUSE_CANDIDATE','NO_STRUCTURE_AVAILABLE'):
            raise ValueError(f'UNRECOGNIZED_STRUCTURE_STATE:{sid}')
        if sid in new and reuse:raise ValueError(f'AMBIGUOUS_OLD_AND_NEW_MODEL:{sid}')
        model=None;origin=''
        if sid in new:
            model=checked_model(sid,new[sid],r)
            origin='NEW_STAGE08_PREDICTED'
        elif reuse:
            old_id=r['historical_structure_id']
            if not old_id:raise ValueError(f'MISSING_OLD_ID:{sid}')
            model={
                'structure_id':sid,'family_rep':r['family_rep'],
                'source_member':r['family_rep'],'family_size':r['family_size'],
                'priority_class':r['priority_class'],'length_aa':r['length_aa'],
                'sequence_sha256':r['sequence_sha256'],
                'pdb_file':r['historical_pdb_file'],'pdb_sha256':r['historical_pdb_sha256'],
                'score_json':r['historical_score_json'],
                'score_json_sha256':r['historical_score_json_sha256'],
                'mean_plddt':r['historical_mean_plddt'],
                'min_plddt':'','max_plddt':'',
                'ptm':r['historical_ptm'],'confidence_class':'',
                'prediction_status':'PREDICTED',
            }
            model=checked_model(sid,model,r)
            origin=f'HISTORICAL_EXACT_SEQUENCE_REUSE:{old_id}'
        if model:models.append(model)
        states.append({
            'structure_id':sid,'family_rep':r['family_rep'],
            'priority_class':r['priority_class'],'family_size':r['family_size'],
            'length_aa':r['length_aa'],'sequence_sha256':r['sequence_sha256'],
            'structure_state':'PREDICTED_VERIFIED' if model else 'PREDICTION_PENDING',
            'structure_origin':origin,'pdb_sha256':model['pdb_sha256'] if model else '',
            'external_search_state':'READY_FOR_SEARCH' if model else 'NOT_EVALUATED_NO_STRUCTURE',
        })
    return models,states,plan_prov,inventory_file


def refuse_overlap(destination, *sources):
    dest=destination.resolve()
    for source in sources:
        if source is None:continue
        src=source.resolve()
        if dest==src or dest in src.parents or src in dest.parents:
            raise ValueError(f'INPUT_OUTPUT_OVERLAP:{destination}:{source}')


def prepare(args):
    if args.out_dir.exists():raise ValueError(f'REFUSE_EXISTING_OUTPUT:{args.out_dir}')
    refuse_overlap(args.out_dir,args.plan_dir,args.queries,args.predicted_structures)
    models,states,plan_prov,inventory=build_preparation(args)
    print(f'ALL_FAMILIES={len(states)}')
    print(f'VERIFIED_SEARCH_READY={len(models)}')
    print(f'PREDICTION_PENDING={len(states)-len(models)}')
    if args.mode=='check':
        print('STAGE182_PREPARE_CHECK=PASS\nFILES_WRITTEN=NO');return
    args.out_dir.mkdir(parents=True)
    write_tsv(args.out_dir/'structure_results_for_search.tsv',STAGE08_COLUMNS,models)
    write_tsv(args.out_dir/'all_family_structure_status.tsv',STATE_KEYS,states)
    record={
        'stage':'STAGE09_PARTIAL_PREPARE', 'schema_version':1,
        'source_inventory_sha256':sha(inventory),
        'source_plan_provenance_sha256':sha(plan_prov),
        'source_fasta_sha256':sha(args.queries),
        'source_new_structure_table_sha256':sha(args.predicted_structures) if args.predicted_structures else None,
        'all_family_count':len(states),'search_ready_count':len(models),
        'prediction_pending_count':len(states)-len(models),
        'origin_counts':dict(Counter('OLD_REUSE' if r['structure_origin'].startswith('HISTORICAL_') else 'NEW_PREDICTION' for r in states if r['structure_origin'])),
        'searched':False,'stage10_full_ranking_completed':False,
        'input_structure_sha256':sha(args.out_dir/'structure_results_for_search.tsv'),
        'all_family_status_sha256':sha(args.out_dir/'all_family_structure_status.tsv'),
    }
    (args.out_dir/'stage182_preparation_provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print('STEP182_PREPARE=PASS')
    print(f'OUTPUT_DIR={args.out_dir}')
    print('EXTERNAL_SEARCH=NOT_EXECUTED')


def stage09_module():
    module_path=Path(__file__).with_name('09_external_integration_portable.py')
    if not module_path.is_file():raise ValueError(f'MISSING_EXISTING_STAGE09_MODULE:{module_path}')
    spec=importlib.util.spec_from_file_location('stage09_portable_real',module_path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_merge(args):
    prep=args.prepared_dir
    manifest=prep/'stage182_preparation_provenance.json'
    if not manifest.is_file():raise ValueError('STAGE182_PREPARATION_MISSING')
    info=json.loads(manifest.read_text())
    states_file=prep/'all_family_structure_status.tsv'
    models_file=prep/'structure_results_for_search.tsv'
    if sha(states_file)!=info['all_family_status_sha256'] or sha(models_file)!=info['input_structure_sha256']:
        raise ValueError('STAGE182_PREPARATION_TAMPERED')
    _,states=table(states_file,STATE_KEYS)
    _,model_rows=table(models_file,STAGE08_COLUMNS,allow_empty=True)
    full=unique(states,'structure_id');ids=unique(model_rows,'structure_id')
    expected={sid for sid,r in full.items() if r['external_search_state']=='READY_FOR_SEARCH'}
    if set(ids)!=expected:raise ValueError('SEARCH_READY_ID_MISMATCH')
    if len(full)!=info['all_family_count'] or len(ids)!=info['search_ready_count']:
        raise ValueError('PROVENANCE_COUNT_MISMATCH')
    module=stage09_module()
    if not ids:
        if args.external_dir:raise ValueError('NO_MODELS_CANNOT_HAVE_EXTERNAL')
        ext_fields=list(module.HEADER);matrix=[]
        external_prov=None
    else:
        if not args.external_dir:raise ValueError('EXTERNAL_RESULTS_REQUIRED_FOR_MODELS')
        matrix_path=args.external_dir/'external_structure_matrix.tsv'
        prov_path=args.external_dir/'integration_provenance.json'
        if not prov_path.is_file():raise ValueError('EXTERNAL_PROVENANCE_MISSING')
        p=json.loads(prov_path.read_text())
        if (p.get('structure_table_sha256')!=sha(models_file) or
            p.get('matrix_sha256')!=sha(matrix_path) or
            p.get('family_count')!=len(ids) or
            p.get('evidence_origin')!='COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES'):
            raise ValueError('EXTERNAL_INTEGRATION_PROVENANCE_MISMATCH')
        # Re-audit the raw Foldseek results and full query/PDB hashes by the original
        # Stage09 integration code. Does not execute Foldseek or write any files.
        metadata=module.load_structures(models_file)
        sources={}
        for prefix,_,profile,maximum in module.SOURCES:
            sourceinfo=p['external_searches'][prefix]
            folder=Path(sourceinfo['directory'])
            hits,validated=module.load_search(folder,profile,maximum,metadata)
            if validated['raw_hits_sha256']!=sourceinfo['raw_hits_sha256'] or validated['query_manifest_sha256']!=sourceinfo['query_manifest_sha256'] or validated['search_provenance_sha256']!=sourceinfo['search_provenance_sha256']:
                raise ValueError(f'EXTERNAL_SEARCH_PROVENANCE_MISMATCH:{prefix}')
            sources[prefix]=hits
        recomputed=module.integrate(metadata,sources)
        ext_fields,matrix=table(matrix_path,module.HEADER)
        if len(matrix)!=len(ids) or unique(matrix,'structure_id')!=unique(recomputed,'structure_id'):
            raise ValueError('EXTERNAL_MATRIX_CONTENT_MISMATCH')
        if set(ext_fields)!=set(module.HEADER):raise ValueError('EXTERNAL_MATRIX_SCHEMA_MISMATCH')
        external_prov=sha(prov_path)
    by_ext=unique(matrix,'structure_id')
    extra=['stage09_'+c for c in module.HEADER if c!='structure_id']
    output=[]
    for row in states:
        sid=row['structure_id'];hit=by_ext.get(sid)
        merged=dict(row)
        merged['stage09_evidence_state']='SEARCH_VERIFIED' if hit else 'NOT_EVALUATED_NO_STRUCTURE'
        for c in module.HEADER:
            if c!='structure_id':merged['stage09_'+c]=hit[c] if hit else ''
        output.append(merged)
    return output,list(STATE_KEYS)+['stage09_evidence_state']+extra,info,external_prov


def merge(args):
    if args.out_dir.exists():raise ValueError(f'REFUSE_EXISTING_OUTPUT:{args.out_dir}')
    refuse_overlap(args.out_dir,args.prepared_dir,args.external_dir)
    output,fields,info,external_sha=build_merge(args)
    verified=sum(row['stage09_evidence_state']=='SEARCH_VERIFIED' for row in output)
    print(f'ALL_FAMILIES={len(output)}')
    print(f'EXTERNAL_EVIDENCE_VERIFIED={verified}')
    print(f'PREDICTION_PENDING={len(output)-verified}')
    if args.mode=='check':
        print('STAGE182_MERGE_CHECK=PASS\nFILES_WRITTEN=NO');return
    args.out_dir.mkdir(parents=True)
    dest=args.out_dir/'all_family_stage09_evidence.tsv'
    write_tsv(dest,fields,output)
    provenance={
        'stage':'STAGE09_PARTIAL_EVIDENCE_MERGED','schema_version':1,
        'source_preparation_manifest_sha256':sha(args.prepared_dir/'stage182_preparation_provenance.json'),
        'source_external_integration_sha256':external_sha,
        'all_family_count':len(output),'externally_evaluated_family_count':verified,
        'prediction_pending_count':len(output)-verified,
        'not_evaluated_is_not_no_hit':True,
        'gpu_prediction_executed':False,'stage10_full_ranking_completed':False,
        'output_sha256':sha(dest),
        'stage09_status_counts':dict(Counter(row['stage09_external_status'] for row in output if row['stage09_evidence_state']=='SEARCH_VERIFIED')),
    }
    (args.out_dir/'stage182_merge_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('STEP182_PARTIAL_STAGE09_MERGE=PASS')
    print(f'OUTPUT_DIR={args.out_dir}')
    print('STAGE10_FINAL_RANKING=NOT_RUN')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare',help='Validate and collect predicted PDB subset')
    a.add_argument('--plan-dir',type=Path,required=True)
    a.add_argument('--queries',type=Path,required=True)
    a.add_argument('--predicted-structures',type=Path)
    a.add_argument('--out-dir',type=Path,required=True)
    a.add_argument('--mode',choices=('check','run'),default='check')
    b=sub.add_parser('merge',help='Verify real Foldseek search and preserve pending families')
    b.add_argument('--prepared-dir',type=Path,required=True)
    b.add_argument('--external-dir',type=Path)
    b.add_argument('--out-dir',type=Path,required=True)
    b.add_argument('--mode',choices=('check','run'),default='check')
    opts=p.parse_args()
    if opts.action=='prepare':prepare(opts)
    else:merge(opts)

if __name__=='__main__':
    main()
