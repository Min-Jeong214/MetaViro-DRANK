#!/usr/bin/env python3
"""Integrity-checked Stage13 partial join of local Stage162 and K2 Stage165.

This joins only observed evidence, NEVER recomputes Stage10 ranking and NEVER
claims completed final integration or selects top candidates. Strictly keep
Stage10 single-sequence external evidence distinct from MSA K2 evidence.
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

HISTORICAL = 'MSA_HISTORICAL_PREDICTION_COLLECTED_QC_VERIFIED'
FRESH = 'MSA_PREDICTION_COLLECTED_QC_VERIFIED'
PENDING = 'MSA_PREDICTION_PENDING'
FALLBACK = 'STAGE08_ONLY_SINGLE_UNIQUE_SEQUENCE'
EXTERNAL_COLUMNS = [
    'structure_id','external_status','raw_sources','convincing_sources',
    'BFVD_best_target','BFVD_best_qtm','BFVD_best_ttm',
    'BFVD_best_is_convincing','msa_pdb_sha256',
]
NEW_COLUMNS = [
    'k2_search_status', 'k2_external_status', 'k2_raw_sources',
    'k2_convincing_sources', 'k2_BFVD_best_target', 'k2_BFVD_best_qtm',
    'k2_BFVD_best_ttm','k2_BFVD_best_is_convincing',
    'k2_msa_pdb_sha256','stage13_partial_evidence_status',
]
K2_STATUS = {'NO_RAW_HIT_E1', 'WEAK_RAW_HIT_ONLY', 'CONVINCING_STRUCTURE_CONNECTION'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for piece in iter(lambda: stream.read(1 << 20), b''):
            h.update(piece)
    return h.hexdigest()


def check_file(path):
    path = Path(path).resolve()
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f'MISSING_OR_EMPTY_INPUT:{path}')
    return path


def check_hash(path, digest):
    path = check_file(path)
    if sha(path) != digest:
        raise ValueError(f'SOURCE_SHA256_MISMATCH:{path}')
    return path


def obj(path):
    item = json.loads(check_file(path).read_text(encoding='utf-8'))
    if not isinstance(item, dict):
        raise ValueError(f'PROVENANCE_NOT_OBJECT:{path}')
    return item


def table(path, required=(), columns=None):
    with check_file(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields, rows = reader.fieldnames or [], list(reader)
    if (not fields or len(fields) != len(set(fields)) or
        not set(required).issubset(fields) or not rows or
        (columns is not None and fields != columns) or
        any(None in r or any(v is None for v in r.values()) for r in rows)):
        raise ValueError(f'INVALID_TABLE:{path}')
    ids = [r['structure_id'] for r in rows]
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', x) for x in ids) or len(set(ids)) != len(ids):
        raise ValueError(f'INVALID_OR_DUPLICATE_IDS:{path}')
    return fields, rows


def index(rows):
    return {r['structure_id']: r for r in rows}


def preflight(args):
    local = Path(args.stage162_output).resolve()
    k2 = Path(args.stage165_output).resolve()
    ranked = Path(args.stage09_10_output).resolve()
    out = Path(args.out_dir).resolve()
    if out.exists():
        raise ValueError('REFUSE_EXISTING_OUTPUT')
    for source in (local,k2,ranked):
        if not source.is_dir() or out == source or out in source.parents or source in out.parents:
            raise ValueError(f'INVALID_INPUT_OUTPUT_DIRECTORY:{source}')

    local_prov_path = local/'orchestration_provenance.json'
    local_prov = obj(local_prov_path)
    if (local_prov.get('stage') != 'STAGE12_13_PARTIAL_HISTORICAL_LOCAL_ANALYSIS'
        or local_prov.get('stage13_full_integration_completed') is not False
        or local_prov.get('stage10_scores_recalculated') is not False
        or local_prov.get('final_top5_top3_selected') is not False):
        raise ValueError('STAGE162_UNEXPECTED_COMPLETION_STATE')
    artifact_hashes = local_prov.get('artifacts', {})
    if not isinstance(artifact_hashes, dict):
        raise ValueError('INVALID_STAGE162_ARTIFACT_MANIFEST')
    local_status_path = local/'stage12_13_partial_status.tsv'
    joined_path = local/'13_preliminary/stage10_plus_partial_evidence.tsv'
    local_pre_path = local/'13_preliminary/preliminary_candidate_evidence.tsv'
    local_pre_prov_path = local/'13_preliminary/preliminary_evidence_provenance.json'
    cons_prov_path = local/'12_conservation/conservation_provenance.json'
    for path in (local_status_path,joined_path,local_pre_path,local_pre_prov_path,cons_prov_path):
        key = path.relative_to(local).as_posix()
        if key not in artifact_hashes:
            raise ValueError(f'STAGE162_ARTIFACT_MISSING:{key}')
        check_hash(path, artifact_hashes[key]['sha256'])
    local_pre_prov = obj(local_pre_prov_path)
    if (local_pre_prov.get('preliminary_matrix_sha256') != sha(local_pre_path)
        or local_pre_prov.get('joined_matrix_sha256') != sha(joined_path)
        or local_pre_prov.get('status_matrix_sha256') != sha(local_status_path)
        or local_pre_prov.get('stage13_external_msa_search_performed') is not False
        or local_pre_prov.get('new_final_top5_selected') is not False):
        raise ValueError('STAGE162_PRELIMINARY_PROVENANCE_MISMATCH')
    cons_prov = obj(cons_prov_path)
    if (cons_prov.get('algorithm') != 'ORIGINAL_STAGE12_FUNCTIONS'
        or cons_prov.get('family_count') != local_prov.get('stage12_conservation_completed_families')):
        raise ValueError('STAGE162_CONSERVATION_PROVENANCE_MISMATCH')
    cons_by = index(cons_prov['families'])
    if len(cons_by) != len(cons_prov['families']):
        raise ValueError('DUPLICATE_CONSERVATION_FAMILY')

    ranking_path = ranked/'10_ranking/confidence_aware_ranking.tsv'
    shortlist_path = ranked/'10_shortlists/top15_experimental_candidates.tsv'
    if local_prov.get('input_stage10_ranking_sha256') != sha(check_file(ranking_path)):
        raise ValueError('STAGE10_RANKING_SHA_MISMATCH')
    ranking_fields, ranking_rows = table(ranking_path, ('structure_id',))
    shortlist_fields, shortlist_rows = table(shortlist_path, ('structure_id',))
    if shortlist_fields != ranking_fields or len(ranking_fields) != 30:
        raise ValueError('STAGE10_RANKING_SCHEMA_MISMATCH')
    rank_by = index(ranking_rows)
    for r in shortlist_rows:
        if rank_by.get(r['structure_id']) != r:
            raise ValueError(f'STAGE10_SHORTLIST_CHANGED:{r["structure_id"]}')

    status_fields, status_rows = table(local_status_path, ('structure_id','stage11_collection_status',
                      'stage12_conservation_status','stage12_spatial_status','stage13_integrated_evidence_status'))
    joined_fields, joined_rows = table(joined_path, ('structure_id','preliminary_very_conserved_positions'))
    prelim_fields, pre_rows = table(local_pre_path, ('structure_id',))
    if (status_fields[:30] != ranking_fields or joined_fields[:len(status_fields)] != status_fields
        or len(status_rows) != len(joined_rows) != len(shortlist_rows)
        or [r['structure_id'] for r in status_rows] != [r['structure_id'] for r in joined_rows]
        or set(index(status_rows)) != set(index(shortlist_rows))):
        raise ValueError('STAGE162_SELECTION_SCHEMA_MISMATCH')
    selected = index(shortlist_rows)
    prelim = index(pre_rows)
    for s, j in zip(status_rows, joined_rows):
        sid=s['structure_id']
        if {k:s[k] for k in ranking_fields} != selected[sid] or any(j[k] != s[k] for k in status_fields):
            raise ValueError(f'STAGE10_RANKING_OR_STAGE162_STATUS_CHANGED:{sid}')
        for field in prelim_fields:
            if field != 'structure_id' and f'preliminary_{field}' in j:
                if j[f'preliminary_{field}'] != prelim.get(sid,{}).get(field,''):
                    raise ValueError(f'STAGE162_PRELIMINARY_VALUE_MISMATCH:{sid}:{field}')
    computed = set()
    pending = set()
    for row in status_rows:
        sid=row['structure_id']
        if row['stage11_collection_status'] in (HISTORICAL,FRESH):
            if (row['stage12_conservation_status'] != 'COMPLETED'
                or row['stage12_spatial_status'] != 'COMPLETED_EXPLORATORY'
                or row['stage13_integrated_evidence_status'] != 'PRELIMINARY_LOCAL_EVIDENCE_ONLY'):
                raise ValueError(f'INVALID_COMPLETED_LOCAL_STATE:{sid}')
            computed.add(sid)
        elif row['stage11_collection_status'] in (PENDING,FALLBACK):
            if (sid in prelim or row['stage12_conservation_status'] == 'COMPLETED'):
                raise ValueError(f'INVALID_PENDING_LOCAL_STATE:{sid}')
            pending.add(sid)
        else:
            raise ValueError(f'UNEXPECTED_STAGE11_STATE:{sid}')
    if (set(prelim) != computed or set(cons_by) != computed
        or len(computed) != local_prov.get('stage13_preliminary_evidence_families')
        or len(pending) != local_prov.get('pending_or_fallback_families')
        or len(status_rows) != local_prov.get('selected_families')
        or len(computed)+len(pending) != len(status_rows)):
        raise ValueError('LOCAL_EVIDENCE_COVERAGE_MISMATCH')

    bridge_path=k2/'13_k2_bridge/msa_external_evidence_stage13.tsv'
    bridge_prov_path=k2/'13_k2_bridge/msa_external_evidence_provenance.json'
    bridge_meta=obj(bridge_prov_path)
    if (bridge_meta.get('source') != 'VERIFIED_STAGE_K2_EXTERNAL_SEARCH'
        or bridge_meta.get('query_origin') != 'STAGE11_FAMILY_MSA'
        or bridge_meta.get('control_only') is not False
        or bridge_meta.get('final_selection_performed') is not False
        or bridge_meta.get('stage13_evidence_sha256') != sha(check_file(bridge_path))):
        raise ValueError('UNVERIFIED_STAGE165_BRIDGE')
    _, external_rows=table(bridge_path, columns=EXTERNAL_COLUMNS)
    ext=index(external_rows)
    if (set(ext) != computed or len(ext) != bridge_meta.get('families')):
        raise ValueError(f'K2_LOCAL_COVERAGE_MISMATCH:missing={sorted(computed-set(ext))}:extra={sorted(set(ext)-computed)}')
    k2_input_path=k2/'msa_structure_results.tsv'
    k2_source_meta_path=k2/'msa_foldseek_input_provenance.json'
    k2_source=obj(k2_source_meta_path)
    if (bridge_meta.get('k2_input_sha256') != sha(check_file(k2_input_path))
        or bridge_meta.get('k2_provenance_sha256') != sha(k2_source_meta_path)
        or k2_source.get('structure_table_sha256') != sha(k2_input_path)
        or k2_source.get('query_origin') != 'STAGE11_FAMILY_MSA'
        or k2_source.get('structure_count') != len(computed)
        or k2_source.get('original_gpu_execution_provenance_verified') is not False):
        raise ValueError('K2_SEARCH_INPUT_PROVENANCE_MISMATCH')
    _, input_rows=table(k2_input_path, ('structure_id','pdb_file','pdb_sha256'))
    k2_input=index(input_rows)
    if set(k2_input) != computed:
        raise ValueError('K2_STRUCTURE_SET_MISMATCH')
    for sid in computed:
        s=cons_by[sid]; item=k2_input[sid]; e=ext[sid]
        if not (e['msa_pdb_sha256'] == item['pdb_sha256'] == s['pdb_sha256']):
            raise ValueError(f'K2_LOCAL_PDB_SHA_MISMATCH:{sid}')
        check_hash(item['pdb_file'], item['pdb_sha256'])
        if e['external_status'] not in K2_STATUS or e['BFVD_best_is_convincing'] not in ('0','1'):
            raise ValueError(f'INVALID_K2_EVIDENCE_STATUS:{sid}')
        if e['external_status']=='NO_RAW_HIT_E1' and (e['raw_sources'] or e['convincing_sources']):
            raise ValueError(f'INVALID_K2_NO_HIT:{sid}')
        if e['external_status']=='WEAK_RAW_HIT_ONLY' and (not e['raw_sources'] or e['convincing_sources']):
            raise ValueError(f'INVALID_K2_WEAK_HIT:{sid}')
        if e['external_status']=='CONVINCING_STRUCTURE_CONNECTION' and not e['convincing_sources']:
            raise ValueError(f'INVALID_K2_CONVINCING_STATUS:{sid}')
    # Bridge must be bound to the actual verified Foldseek matrix and its three search logs.
    ext_matrix_path=k2/'09_external/external_structure_matrix.tsv'
    if bridge_meta.get('external_matrix_sha256') != sha(check_file(ext_matrix_path)):
        raise ValueError('K2_MATRIX_SHA_MISMATCH')
    matches=[]
    for path in (k2/'09_external').glob('*.json'):
        if sha(path)==bridge_meta.get('external_provenance_sha256'):
            matches.append(path)
    if len(matches)!=1:
        raise ValueError('K2_EXTERNAL_PROVENANCE_UNAVAILABLE')
    ext_meta=obj(matches[0])
    if (ext_meta.get('evidence_origin') != 'COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES'
        or ext_meta.get('matrix_sha256') != sha(ext_matrix_path)
        or ext_meta.get('structure_table_sha256') != sha(k2_input_path)
        or set(ext_meta.get('external_searches',{})) != {'pdb','afsp','bfvd'}):
        raise ValueError('K2_EXTERNAL_SEARCH_PROVENANCE_MISMATCH')
    for prefix, profile in (('pdb','pdb_e1'),('afsp','afsp'),('bfvd','bfvd')):
        src=ext_meta['external_searches'][prefix]
        if src.get('source_profile') != profile:
            raise ValueError(f'K2_SEARCH_PROFILE_MISMATCH:{prefix}')
        folder=Path(src['directory'])
        for name,key in (('search_provenance.json','search_provenance_sha256'),
                         ('raw_hits.tsv','raw_hits_sha256'),('query_manifest.tsv','query_manifest_sha256')):
            check_hash(folder/name,src[key])
    return dict(local=local, k2=k2, ranked=ranked, out=out, joined_fields=joined_fields,
                joined_rows=joined_rows, ext=ext, computed=computed, pending=pending,
                local_prov_path=local_prov_path, local_pre_prov_path=local_pre_prov_path,
                bridge_path=bridge_path, bridge_prov_path=bridge_prov_path,
                k2_input_path=k2_input_path, ext_matrix_path=ext_matrix_path,
                shortlist_path=shortlist_path, ranking_path=ranking_path)


def save_tsv(path, fields, rows):
    with path.open('x',encoding='utf-8',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def execute(state):
    out=state['out']
    result=[]
    for original in state['joined_rows']:
        sid=original['structure_id']
        row=dict(original)
        external=state['ext'].get(sid)
        if external:
            row.update(k2_search_status='SEARCH_COMPLETED',
                k2_external_status=external['external_status'],
                k2_raw_sources=external['raw_sources'],
                k2_convincing_sources=external['convincing_sources'],
                k2_BFVD_best_target=external['BFVD_best_target'],
                k2_BFVD_best_qtm=external['BFVD_best_qtm'],
                k2_BFVD_best_ttm=external['BFVD_best_ttm'],
                k2_BFVD_best_is_convincing=external['BFVD_best_is_convincing'],
                k2_msa_pdb_sha256=external['msa_pdb_sha256'],
                stage13_partial_evidence_status='LOCAL_AND_K2_EXTERNAL_EVIDENCE_AVAILABLE')
        else:
            row.update({c:'' for c in NEW_COLUMNS})
            row['k2_search_status']='NOT_SEARCHED_PREDICTION_PENDING'
            row['stage13_partial_evidence_status']='PREDICTION_PENDING'
        result.append(row)
    fields=state['joined_fields'] + NEW_COLUMNS
    out.mkdir(parents=True,exist_ok=False)
    full=out/'stage13_partial_evidence_15.tsv'
    verified=out/'stage13_local_k2_verified_subset.tsv'
    save_tsv(full,fields,result)
    save_tsv(verified,fields,[r for r in result if r['structure_id'] in state['computed']])
    source_sha={
        'stage162_orchestration_provenance_sha256':sha(state['local_prov_path']),
        'stage162_preliminary_evidence_provenance_sha256':sha(state['local_pre_prov_path']),
        'stage165_bridge_provenance_sha256':sha(state['bridge_prov_path']),
        'stage165_external_evidence_sha256':sha(state['bridge_path']),
        'stage165_input_table_sha256':sha(state['k2_input_path']),
        'stage165_external_matrix_sha256':sha(state['ext_matrix_path']),
        'stage10_rank_sha256':sha(state['ranking_path']),
        'stage10_shortlist_sha256':sha(state['shortlist_path']),
    }
    metadata={
        'stage':'STAGE13_PARTIAL_LOCAL_PLUS_K2_EXTERNAL_EVIDENCE',
        'selected_families':len(result),
        'local_and_k2_evidence_families':len(state['computed']),
        'prediction_pending_families':len(state['pending']),
        'external_status_counts_for_verified_families':dict(Counter(
            state['ext'][sid]['external_status'] for sid in state['computed'])),
        'historical_gpu_execution_provenance_verified':False,
        'stage10_ranking_recalculated':False,
        'stage13_full_integration_completed':False,
        'final_top5_top3_selected':False,
        'biochemical_function_experimentally_verified':False,
        'single_sequence_stage10_external_evidence_preserved':True,
        'msa_k2_external_evidence_namespaced':True,
        'msa_k2_external_search_reexecuted':False,
        'source_sha256':source_sha,
        'partial_table_sha256':sha(full),
        'verified_subset_sha256':sha(verified),
        'output_columns':len(fields),
        'notes':'Partial integration of existing results; not a final candidate selection.'
    }
    (out/'stage13_partial_evidence_provenance.json').write_text(
        json.dumps(metadata,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('STAGE10_RANKING_PRESERVED=PASS')
    print('STAGE12_LOCAL_EVIDENCE_VERIFIED='+str(len(state['computed'])))
    print('STAGE_K2_EXTERNAL_EVIDENCE_VERIFIED='+str(len(state['ext'])))
    print('PREDICTION_PENDING='+str(len(state['pending'])))
    print('K2_EXTERNAL_STATUS_COUNTS='+str(metadata['external_status_counts_for_verified_families']))
    print('STAGE13_FULL_INTEGRATION=NOT_COMPLETED')
    print('FINAL_TOP5_TOP3_SELECTED=NO')
    print('STEP166_PARTIAL_K2_MERGE=PASS')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage162-output',type=Path,required=True)
    ap.add_argument('--stage165-output',type=Path,required=True)
    ap.add_argument('--stage09-10-output',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True)
    ap.add_argument('--mode',choices=('check','run'),default='check')
    args=ap.parse_args()
    state=preflight(args)
    print(f'SELECTED_FAMILIES={len(state["joined_rows"])}')
    print(f'STAGE12_COMPLETED_AND_K2_VERIFIED={len(state["computed"])}')
    print(f'PENDING_PREDICTIONS={len(state["pending"])}')
    if args.mode=='check':
        print('CHECK_MODE_NO_WRITE=PASS')
    else:
        execute(state)

if __name__=='__main__':
    try:
        main()
    except (OSError,ValueError,KeyError,TypeError,csv.Error) as exc:
        sys.exit('ERROR: '+str(exc))
