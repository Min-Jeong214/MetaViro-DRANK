#!/usr/bin/env python3
"""Synthetic Stage166 integration and tamper tests; no Foldseek, no GPU."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/13_partial_k2_merge_portable.py'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(value, encoding='utf-8')
    return p

def js(p, value):
    return write(p, json.dumps(value, indent=2) + '\n')

def tab(p, cols, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    return p

def run(loc, ext, ranked, out, mode='check', success=True):
    proc = subprocess.run([sys.executable, '-B', str(SCRIPT),
        '--stage162-output', str(loc), '--stage165-output', str(ext),
        '--stage09-10-output', str(ranked), '--out-dir', str(out),
        '--mode', mode], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if success and proc.returncode:
        raise AssertionError(proc.stdout)
    if not success and not proc.returncode:
        raise AssertionError('Tampered source was accepted')
    return proc.stdout

with tempfile.TemporaryDirectory(prefix='metaviro_stage166_') as dirname:
    root=Path(dirname); loc=root/'local'; ext=root/'k2'; rank=root/'rank'
    topfields=['structure_id','family_rep'] + [f'field{i:02d}' for i in range(28)]
    base=[dict(zip(topfields,[sid,f'rep{sid}',*['v']*28])) for sid in ('F001','F002')]
    ranking=tab(rank/'10_ranking/confidence_aware_ranking.tsv',topfields,base)
    shortlist=tab(rank/'10_shortlists/top15_experimental_candidates.tsv',topfields,base)
    statusfields=topfields+['stage11_collection_status','stage11_msa_qc_status','stage12_conservation_status',
            'stage12_spatial_status','stage13_msa_external_status','stage13_integrated_evidence_status']
    statuses=[dict(base[0],stage11_collection_status='MSA_HISTORICAL_PREDICTION_COLLECTED_QC_VERIFIED',
        stage11_msa_qc_status='QC_VERIFIED_HISTORICAL',stage12_conservation_status='COMPLETED',
        stage12_spatial_status='COMPLETED_EXPLORATORY',stage13_msa_external_status='PENDING',
        stage13_integrated_evidence_status='PRELIMINARY_LOCAL_EVIDENCE_ONLY'),
        dict(base[1],stage11_collection_status='MSA_PREDICTION_PENDING',
        stage11_msa_qc_status='PENDING_PREDICTION',stage12_conservation_status='PENDING_PREDICTION',
        stage12_spatial_status='PENDING_PREDICTION',stage13_msa_external_status='PENDING',
        stage13_integrated_evidence_status='STAGE10_ONLY_PREDICTION_PENDING')]
    statustable=tab(loc/'stage12_13_partial_status.tsv',statusfields,statuses)
    prelim=tab(loc/'13_preliminary/preliminary_candidate_evidence.tsv',
        ['structure_id','very_conserved_positions'],[{'structure_id':'F001','very_conserved_positions':'12'}])
    joined=tab(loc/'13_preliminary/stage10_plus_partial_evidence.tsv',
        statusfields+['preliminary_very_conserved_positions'],
        [dict(row, preliminary_very_conserved_positions=('12' if row['structure_id']=='F001' else '')) for row in statuses])
    pp=js(loc/'13_preliminary/preliminary_evidence_provenance.json',{
        'preliminary_matrix_sha256':sha(prelim),'joined_matrix_sha256':sha(joined),
        'status_matrix_sha256':sha(statustable),'stage13_external_msa_search_performed':False,
        'new_final_top5_selected':False})
    pdb=write(root/'F001.pdb','PDB SYNTHETIC\n')
    cons=js(loc/'12_conservation/conservation_provenance.json',{'algorithm':'ORIGINAL_STAGE12_FUNCTIONS',
        'family_count':1,'families':[{'structure_id':'F001','pdb_sha256':sha(pdb)}]})
    watched=[statustable,joined,prelim,pp,cons]
    localmeta=js(loc/'orchestration_provenance.json',{'stage':'STAGE12_13_PARTIAL_HISTORICAL_LOCAL_ANALYSIS',
        'stage13_full_integration_completed':False,'stage10_scores_recalculated':False,'final_top5_top3_selected':False,
        'input_stage10_ranking_sha256':sha(ranking),'stage12_conservation_completed_families':1,
        'stage13_preliminary_evidence_families':1,'pending_or_fallback_families':1,'selected_families':2,
        'artifacts':{str(x.relative_to(loc)):{'sha256':sha(x)} for x in watched}})
    k2in=tab(ext/'msa_structure_results.tsv',['structure_id','pdb_file','pdb_sha256'],[
        {'structure_id':'F001','pdb_file':str(pdb),'pdb_sha256':sha(pdb)}])
    k2meta=js(ext/'msa_foldseek_input_provenance.json',{
        'structure_table_sha256':sha(k2in),'query_origin':'STAGE11_FAMILY_MSA','structure_count':1,
        'original_gpu_execution_provenance_verified':False})
    matrix=write(ext/'09_external/external_structure_matrix.tsv','structure_id\tstatus\nF001\tOK\n')
    searches={}
    for prefix,profile in (('pdb','pdb_e1'),('afsp','afsp'),('bfvd','bfvd')):
        folder=ext/'search'/profile
        searchmeta=js(folder/'search_provenance.json', {'source_profile':profile,'control_only':False})
        hits=write(folder/'raw_hits.tsv','query\tresult\nF001\tX\n')
        qmanifest=write(folder/'query_manifest.tsv','structure_id\nF001\n')
        searches[prefix]={'source_profile':profile,'directory':str(folder),
            'search_provenance_sha256':sha(searchmeta),'raw_hits_sha256':sha(hits),
            'query_manifest_sha256':sha(qmanifest)}
    extmeta=js(ext/'09_external/external_provenance.json',{
        'evidence_origin':'COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES',
        'structure_table_sha256':sha(k2in),'matrix_sha256':sha(matrix),
        'external_searches':searches})
    extcols=['structure_id','external_status','raw_sources','convincing_sources','BFVD_best_target',
        'BFVD_best_qtm','BFVD_best_ttm','BFVD_best_is_convincing','msa_pdb_sha256']
    evidence=tab(ext/'13_k2_bridge/msa_external_evidence_stage13.tsv', extcols,[{
        'structure_id':'F001','external_status':'CONVINCING_STRUCTURE_CONNECTION',
        'raw_sources':'BFVD','convincing_sources':'BFVD','BFVD_best_target':'synthetic_target',
        'BFVD_best_qtm':'0.6','BFVD_best_ttm':'0.6','BFVD_best_is_convincing':'1','msa_pdb_sha256':sha(pdb)}])
    bridge=js(ext/'13_k2_bridge/msa_external_evidence_provenance.json',{
        'source':'VERIFIED_STAGE_K2_EXTERNAL_SEARCH','query_origin':'STAGE11_FAMILY_MSA',
        'control_only':False,'final_selection_performed':False,'stage13_evidence_sha256':sha(evidence),
        'families':1,'k2_input_sha256':sha(k2in),'k2_provenance_sha256':sha(k2meta),
        'external_matrix_sha256':sha(matrix),'external_provenance_sha256':sha(extmeta)})
    check_output=root/'check_only'
    print(run(loc,ext,rank,check_output))
    assert not check_output.exists()
    out=root/'merged'
    print(run(loc,ext,rank,out,mode='run'))
    with (out/'stage13_partial_evidence_15.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==2
    assert rows[0]['k2_external_status']=='CONVINCING_STRUCTURE_CONNECTION'
    assert rows[1]['k2_search_status']=='NOT_SEARCHED_PREDICTION_PENDING'
    assert rows[1]['k2_external_status']==''
    assert all(a==b for a,b in zip([r['field03'] for r in rows],[r['field03'] for r in base]))
    assert json.loads((out/'stage13_partial_evidence_provenance.json').read_text())['final_top5_top3_selected'] is False
    print('STEP166_PARTIAL_INTEGRATION_SYNTHETIC=PASS')
    pdb.write_text('PDB SYNTHETIC TAMPERED\n')
    print(run(loc,ext,rank,root/'tampered',success=False))
    assert not (root/'tampered').exists()
    print('STEP166_TAMPER_REJECTION=PASS')
    print('ALL_TESTS_PASS')
