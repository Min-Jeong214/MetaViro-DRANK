#!/usr/bin/env python3
"""Stage 12 conservation/spatial and preliminary Stage 13 for partial Stage 11 QC.

Read-only source stages. Historical predictions use THEIR ACTUAL source A3M,
not newly-prepared MSA. Existing original biological algorithms are invoked
unaltered. The Stage13 result is PRELIMINARY; NO MSA Foldseek is performed.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import runpy
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ORIGINAL = SCRIPTS / '12_msa_structure_conservation.py'
SPATIAL = SCRIPTS / '12_charged_spatial_clusters_exploratory.py'
PRELIM = SCRIPTS / '13_candidate_evidence_preliminary_portable.py'
COLLECTED = 'MSA_PREDICTION_COLLECTED_QC_VERIFIED'
HISTORICAL = 'MSA_HISTORICAL_PREDICTION_COLLECTED_QC_VERIFIED'
PENDING = 'MSA_PREDICTION_PENDING'
FALLBACK = 'STAGE08_ONLY_SINGLE_UNIQUE_SEQUENCE'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def file(path):
    path = Path(path).resolve()
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f'MISSING_OR_EMPTY_SOURCE:{path}')
    return path


def checked(path, digest):
    path = file(path)
    if sha(path) != digest:
        raise ValueError(f'SOURCE_SHA256_MISMATCH:{path}')
    return path


def meta(path):
    value = json.loads(file(path).read_text(encoding='utf-8'))
    if not isinstance(value, dict):
        raise ValueError(f'INVALID_PROVENANCE:{path}')
    return value


def tsv(path):
    with file(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        rows = list(reader)
    if not fields or len(set(fields)) != len(fields) or not rows or any(
        None in r or any(v is None for v in r.values()) for r in rows
    ):
        raise ValueError(f'INVALID_TABLE:{path}')
    return fields, rows


def indexed(rows):
    result = {}
    for row in rows:
        sid = row['structure_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or sid in result:
            raise ValueError(f'INVALID_OR_DUPLICATE_ID:{sid}')
        result[sid] = row
    return result


def save_tsv(path, columns, rows):
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def save_meta(path, data):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def aligned_a3m(path, structure_id, prediction):
    """Remove A3M insertion states; use query-anchored match columns only."""
    records = []
    with path.open(encoding='utf-8') as stream:
        for raw in stream:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
                if not name:
                    raise ValueError('EMPTY_A3M_HEADER')
                records.append([name, ''])
            elif not records:
                raise ValueError('A3M_SEQUENCE_BEFORE_HEADER')
            else:
                records[-1][1] += line
    if not records or records[0][0] != structure_id:
        raise ValueError(f'A3M_QUERY_ID_MISMATCH:{structure_id}')
    aligned = []
    for name, raw in records:
        if not raw or any(not(c.isascii() and (c.isalpha() or c in '-.')) for c in raw):
            raise ValueError(f'INVALID_A3M_CHARACTERS:{structure_id}')
        # Lowercase letters and dots are insert states in A3M, not match columns.
        seq = ''.join(c for c in raw if c != '.' and not c.islower())
        if any(c not in 'ACDEFGHIKLMNPQRSTVWYXBUZOJ-' for c in seq):
            raise ValueError(f'INVALID_A3M_MATCH_COLUMN:{structure_id}')
        aligned.append((name, seq))
    query = aligned[0][1]
    if '-' in query or not query or len(query) != prediction['query_length']:
        raise ValueError(f'INVALID_QUERY_ALIGNMENT:{structure_id}')
    if hashlib.sha256(query.encode()).hexdigest() != prediction['query_sha256']:
        raise ValueError(f'QUERY_SHA_MISMATCH:{structure_id}')
    if any(len(seq) != len(query) for _, seq in aligned):
        raise ValueError(f'A3M_MATCH_WIDTH_MISMATCH:{structure_id}')
    if len(aligned) != prediction['sequence_count']:
        raise ValueError(f'A3M_SEQUENCE_COUNT_MISMATCH:{structure_id}')
    return aligned


def validate_inputs(args):
    source = args.adaptive_output.resolve()
    ranking_dir = args.stage09_10_output.resolve()
    out = args.out_dir.resolve()
    if out.exists():
        raise ValueError('REFUSE_EXISTING_OUTPUT_DIRECTORY')
    for folder in (source, ranking_dir):
        if not folder.is_dir() or out == folder or out in folder.parents or folder in out.parents:
            raise ValueError(f'SOURCE_DIRECTORY_INVALID_OR_OVERLAP:{folder}')
    if not all(p.is_file() for p in (ORIGINAL, SPATIAL, PRELIM)):
        raise ValueError('ORIGINAL_STAGE12_OR_STAGE13_SOURCE_MISSING')
    adaptive_manifest_path = file(source / 'orchestration_provenance.json')
    manifest = meta(adaptive_manifest_path)
    if (manifest.get('stage12_13_analysis_completed') is not False or
        manifest.get('stage10_scores_recalculated') is not False or
        manifest.get('gpu_prediction_performed') is not False):
        raise ValueError('STAGE11_SOURCE_STATE_NOT_VERIFIED')
    ranking = file(ranking_dir / '10_ranking/confidence_aware_ranking.tsv')
    shortlist = file(ranking_dir / '10_shortlists/top15_experimental_candidates.tsv')
    if sha(ranking) != manifest['stage10_ranking_sha256']:
        raise ValueError('STAGE10_RANKING_HASH_MISMATCH')
    for label, entry in manifest.get('artifacts', {}).items():
        checked(entry['path'], entry['sha256'])
    ranking_fields, ranking_rows = tsv(ranking)
    short_fields, selected = tsv(shortlist)
    status_path = file(source / '12_13_status/stage10_ranked_evaluation_status.tsv')
    status_fields, rows = tsv(status_path)
    if not (len(ranking_fields) == 30 and ranking_fields == short_fields and
            len(status_fields) == 36 and status_fields[:30] == ranking_fields and
            len(rows) == len(selected) == manifest['selected_families']):
        raise ValueError('SELECTION_OR_STATUS_SCHEMA_MISMATCH')
    ranked = indexed(ranking_rows)
    selection = indexed(selected)
    records = indexed(rows)
    if set(selection) != set(records):
        raise ValueError('STAGE10_SELECTION_SET_MISMATCH')
    for sid, row in records.items():
        if {k:row[k] for k in ranking_fields} != selection[sid] or ranked[sid] != selection[sid]:
            raise ValueError(f'STAGE10_RANKING_CHANGED:{sid}')
    qc_dir = source / '11_collect/msa_qc'
    qc_path = file(qc_dir / 'msa_structure_qc.tsv')
    qc_prov = meta(qc_dir / 'msa_structure_qc_provenance.json')
    if qc_prov.get('algorithm') != 'VERIFIED_SINGLE_VS_FAMILY_MSA' or qc_prov.get('msa_qc_sha256') != sha(qc_path):
        raise ValueError('STAGE11_QC_HASH_MISMATCH')
    _, qc_rows = tsv(qc_path)
    qc_by_id = indexed(qc_rows)
    qc_evidence = indexed(qc_prov['evidence'])
    index_path = file(source / '11_collect/eligible_prediction_index.tsv')
    _, indices = tsv(index_path)
    prediction_by_id = indexed(indices)
    if qc_prov.get('prediction_index_sha256') != sha(index_path):
        raise ValueError('QC_PREDICTION_INDEX_HASH_MISMATCH')
    evaluated, pending = [], []
    for sid in [r['structure_id'] for r in selected]:
        row = records[sid]
        status = row['stage11_collection_status']
        if status in (HISTORICAL, COLLECTED):
            target_qc = 'QC_VERIFIED_HISTORICAL' if status == HISTORICAL else 'QC_VERIFIED'
            if (row['stage11_msa_qc_status'] != target_qc or
                row['stage12_conservation_status'] != 'PENDING_ANALYSIS' or
                row['stage12_spatial_status'] != 'PENDING_ANALYSIS' or
                sid not in qc_by_id or sid not in qc_evidence or sid not in prediction_by_id):
                raise ValueError(f'INVALID_VERIFIED_STAGE11_STATUS:{sid}')
            path = checked(prediction_by_id[sid]['prediction_provenance'],
                           qc_evidence[sid]['stage11_record_sha256'])
            prediction = meta(path)
            historical = status == HISTORICAL
            if prediction.get('structure_id') != sid or prediction.get('model_type') != 'alphafold2_ptm':
                raise ValueError(f'INVALID_PREDICTION_MODEL:{sid}')
            if historical:
                if (prediction.get('status') != 'HISTORICAL_IMPORTED' or
                    prediction.get('provenance_type') != 'HISTORICAL_MSA_REUSE' or
                    prediction.get('original_run_provenance_verified') is not False or
                    prediction.get('gpu_prediction_performed') is not False):
                    raise ValueError(f'INVALID_HISTORICAL_PREDICTION_ATTRIBUTION:{sid}')
                checked(prediction['prepared_a3m'], prediction['prepared_a3m_sha256'])
            elif prediction.get('status') != 'COMPLETED':
                raise ValueError(f'INVALID_NEW_PREDICTION_ATTRIBUTION:{sid}')
            a3m = checked(prediction['source_a3m'], prediction['a3m_sha256'])
            pdb = checked(prediction['pdb_file'], prediction['pdb_sha256'])
            score = checked(prediction['score_json'], prediction['score_json_sha256'])
            if (qc_evidence[sid]['a3m_sha256'] != sha(a3m) or
                qc_by_id[sid]['msa_pdb_sha256'] != sha(pdb)):
                raise ValueError(f'STAGE11_QC_SOURCE_MISMATCH:{sid}')
            scores = meta(score)
            if (len(scores.get('plddt', [])) != int(prediction['query_length']) or
                not all(isinstance(v, (int,float)) and math.isfinite(v) and 0 <= v <= 100
                        for v in scores['plddt'])):
                raise ValueError(f'INVALID_RESIDUE_CONFIDENCE:{sid}')
            aligned = aligned_a3m(a3m, sid, prediction)
            evaluated.append(dict(sid=sid, status=status, prediction=prediction,
                                  record_path=path, a3m=a3m, pdb=pdb, score=score,
                                  aligned=aligned, qc=qc_by_id[sid]))
        elif status == PENDING:
            if (row['stage11_msa_qc_status'] != 'PENDING_PREDICTION' or
                row['stage12_conservation_status'] != 'PENDING_PREDICTION' or
                row['stage12_spatial_status'] != 'PENDING_PREDICTION' or
                sid in qc_by_id or sid in prediction_by_id):
                raise ValueError(f'INVALID_PENDING_STATE:{sid}')
            pending.append(sid)
        elif status == FALLBACK:
            if (row['stage11_msa_qc_status'] != 'NOT_EVALUATED' or
                row['stage12_conservation_status'] != 'NOT_EVALUATED' or
                sid in qc_by_id):
                raise ValueError(f'INVALID_FALLBACK_STATE:{sid}')
            pending.append(sid)
        else:
            raise ValueError(f'UNSUPPORTED_STAGE11_STATE:{sid}:{status}')
    if set(qc_by_id) != {r['sid'] for r in evaluated} or set(prediction_by_id) != set(qc_by_id) or set(qc_evidence) != set(qc_by_id):
        raise ValueError('STAGE11_QC_COVERAGE_MISMATCH')
    if (manifest.get('historical_imported_families') != sum(r['status'] == HISTORICAL for r in evaluated)
        or manifest.get('msa_prediction_pending_families') != sum(records[s]['stage11_collection_status']==PENDING for s in pending)):
        raise ValueError('STAGE11_MANIFEST_COUNT_MISMATCH')
    return dict(source=source, down=ranking_dir, output=out, manifest=adaptive_manifest_path,
                ranking=ranking, shortlist=shortlist, status=status_path, status_fields=status_fields,
                rows=rows, evaluated=evaluated, pending=pending, qc=qc_path,
                qc_prov=qc_dir / 'msa_structure_qc_provenance.json')


def compute(state):
    original = runpy.run_path(str(ORIGINAL), run_name='__stage162_original__')
    spatial = runpy.run_path(str(SPATIAL), run_name='__stage162_spatial__')
    preliminary = runpy.run_path(str(PRELIM), run_name='__stage162_stage13__')
    with tempfile.TemporaryDirectory(prefix='metaviro_stage162_alignment_') as d:
        tmp = Path(d)
        computed = []
        for item in state['evaluated']:
            sid = item['sid']
            path = tmp / f'{sid}.aligned.faa'
            with path.open('x', encoding='utf-8') as stream:
                for name, seq in item['aligned']:
                    stream.write(f'>{name}\n{seq}\n')
            residues = original['build_residue_table'](sid, path, item['pdb'])
            if len(residues) != item['prediction']['query_length']:
                raise ValueError(f'ORIGINAL_STAGE12_RESIDUE_COUNT_MISMATCH:{sid}')
            atom = spatial['read_ca'](item['pdb'])
            if len(atom) != len(residues):
                raise ValueError(f'SPATIAL_RESIDUE_COUNT_MISMATCH:{sid}')
            eligible = []
            for row in residues:
                pos = row['query_position']
                a = atom.get(pos)
                if (a is None or a['aa'] != row['query_aa'] or
                    abs(a['plddt'] - row['msa_plddt']) > 0.051):
                    raise ValueError(f'SPATIAL_CA_MAPPING_MISMATCH:{sid}:{pos}')
                if row['priority_residue'] and a['aa'] in spatial['CHARGED']:
                    eligible.append(a)
            clusters = spatial['find_clusters'](eligible, 8.0, 3)
            computed.append(dict(**item, residues=residues,
                 summary=original['build_summary'](sid, residues),
                 windows=original['make_windows'](sid, residues), clusters=clusters,
                 eligible_count=len(eligible)))
    return original, spatial, preliminary, computed


def execute(state, original, spatial, preliminary, computed):
    out = state['output']
    out.mkdir(parents=True, exist_ok=False)
    cons = out / '12_conservation'; cons.mkdir()
    align_dir = cons / 'prediction_input_aligned'; align_dir.mkdir()
    families = []; summaries=[]; windows=[]; cluster_rows=[]; spatial_summaries=[]; cluster_map=defaultdict(list)
    for item in computed:
        sid=item['sid']
        alignment=align_dir / f'{sid}.aligned.faa'
        with alignment.open('x', encoding='utf-8') as stream:
            for name, sequence in item['aligned']:
                stream.write(f'>{name}\n{sequence}\n')
        residue=cons / f'{sid}_residue_conservation.tsv'
        original['write_residue_table'](residue, item['residues'])
        families.append(dict(structure_id=sid, alignment=str(alignment),
             alignment_sha256=sha(alignment), a3m=str(item['a3m']), a3m_sha256=sha(item['a3m']),
             prediction_provenance=str(item['record_path']),
             prediction_provenance_sha256=sha(item['record_path']),
             pdb_sha256=sha(item['pdb']), sequence_count=len(item['aligned']),
             residue_count=len(item['residues']), residue_table_sha256=sha(residue),
             prediction_type=('HISTORICAL' if item['status']==HISTORICAL else 'FRESH')))
        summaries.append(item['summary']); windows.extend(item['windows'])
        for rank, cluster in enumerate(item['clusters'],1):
            members=cluster['members']; record=dict(
                structure_id=sid, cluster_rank=str(rank),
                residues=','.join(f"{r['aa']}{r['position']}" for r in members),
                residue_count=str(len(members)),
                acidic_count=str(sum(r['aa'] in spatial['ACIDIC'] for r in members)),
                basic_count=str(sum(r['aa'] in spatial['BASIC'] for r in members)),
                mean_plddt=f"{cluster['mean_plddt']:.2f}",
                max_edge_ca_distance_angstrom=f"{cluster['max_edge_distance']:.3f}")
            cluster_rows.append(record); cluster_map[sid].append(record)
        spatial_summaries.append(dict(structure_id=sid,
            length=str(len(item['residues'])), eligible_charged_residues=str(item['eligible_count']),
            cluster_count=str(len(item['clusters'])), status=('CLUSTERS_DETECTED' if item['clusters'] else 'NO_CLUSTER_DETECTED')))
    summary_path=cons/'candidate_conservation_summary.tsv'
    save_tsv(summary_path, original['SUMMARY_FIELDS'] if 'SUMMARY_FIELDS' in original else [
        'structure_id','length','msa_nseq','very_conserved_positions','priority_structural_positions',
        'priority_acidic','priority_basic','priority_nucleophilic'], summaries)
    windows_path=cons/'conserved_windows.tsv'
    formatted=[{**r,'mean_query_fraction':f"{r['mean_query_fraction']:.3f}",
                    'mean_plddt':f"{r['mean_plddt']:.2f}"} for r in windows]
    save_tsv(windows_path,['structure_id','start','end','sequence','conserved_structural','charged',
        'acidic','basic','mean_query_fraction','mean_plddt'],formatted)
    cons_prov=cons/'conservation_provenance.json'
    save_meta(cons_prov, dict(algorithm='ORIGINAL_STAGE12_FUNCTIONS',
        original_script_sha256=sha(ORIGINAL), alignment_source='ACTUAL_PREDICTION_INPUT_A3M',
        alignment_conversion='DROP_A3M_INSERT_STATE_LOWERCASE_AND_DOTS',
        new_msa_substituted_for_historical=False, selection_source='STAGE11_QC_VERIFIED_ONLY',
        final_top5_automatically_selected=False,
        stage11_manifest_sha256=sha(state['manifest']),stage10_ranking_sha256=sha(state['ranking']),
        family_count=len(families), families=families,
        summary_sha256=sha(summary_path),windows_sha256=sha(windows_path)))
    spatial_dir=out/'12_spatial'; spatial_dir.mkdir()
    spatial_clusters=spatial_dir/'charged_spatial_clusters.tsv'
    spatial_summary=spatial_dir/'family_spatial_cluster_summary.tsv'
    save_tsv(spatial_clusters,spatial['CLUSTER_FIELDS'],cluster_rows)
    save_tsv(spatial_summary,spatial['SUMMARY_FIELDS'],spatial_summaries)
    spatial_prov=spatial_dir/'spatial_cluster_provenance.json'
    save_meta(spatial_prov, dict(analysis_mode='NEW_EXPLORATORY_L3',
        historical_stage_l3_reproduction=False, geometry='CA_DISTANCE_CONNECTED_COMPONENTS',
        ca_distance_cutoff_angstrom=8.0,min_cluster_residues=3,
        residue_filter='STAGE12_PRIORITY_RESIDUE_AND_DEKRH',biochemical_activity_inferred=False,
        original_spatial_script_sha256=sha(SPATIAL),stage12_provenance_sha256=sha(cons_prov),
        family_count=len(families),cluster_count=len(cluster_rows),
        cluster_table_sha256=sha(spatial_clusters),summary_table_sha256=sha(spatial_summary)))
    stage13=out/'13_preliminary';stage13.mkdir()
    ids=[item['sid'] for item in computed]
    qc_map=indexed(tsv(state['qc'])[1]); cons_map=indexed(summaries); spatial_map=indexed(spatial_summaries)
    prelim_rows=preliminary['assemble'](ids,qc_map,cons_map,spatial_map,cluster_map)
    preliminary_file=stage13/'preliminary_candidate_evidence.tsv'
    save_tsv(preliminary_file,preliminary['FIELDS'],prelim_rows)
    preliminary_by_id=indexed(prelim_rows)
    # The original 36-column source is never changed; append preliminary fields.
    extra=[f for f in preliminary['FIELDS'] if f!='structure_id' and f not in state['status_fields']]
    status_updates=[];joined=[]
    for row in state['rows']:
        sid=row['structure_id']; r=dict(row)
        if sid in preliminary_by_id:
            r['stage12_conservation_status']='COMPLETED'
            r['stage12_spatial_status']='COMPLETED_EXPLORATORY'
            r['stage13_integrated_evidence_status']='PRELIMINARY_LOCAL_EVIDENCE_ONLY'
        status_updates.append(r)
        extended=dict(r)
        evidence=preliminary_by_id.get(sid,{})
        for f in extra:
            extended['preliminary_'+f]=evidence.get(f,'')
        joined.append(extended)
    status_result=out/'stage12_13_partial_status.tsv'
    save_tsv(status_result,state['status_fields'],status_updates)
    joined_path=stage13/'stage10_plus_partial_evidence.tsv'
    save_tsv(joined_path,state['status_fields']+['preliminary_'+f for f in extra],joined)
    save_meta(stage13/'preliminary_evidence_provenance.json',dict(
        status='PRELIMINARY_LOCAL_EVIDENCE_ONLY',stage13_external_msa_search_performed=False,
        original_stage13_assembler_sha256=sha(PRELIM),stage10_scores_recalculated=False,
        new_final_top5_selected=False,biochemical_activity_verified=False,
        selected_families=len(state['rows']),evidence_families=len(ids),
        pending_or_fallback_families=len(state['rows'])-len(ids),
        stage11_status_sha256=sha(state['status']),
        conservation_provenance_sha256=sha(cons_prov),
        spatial_provenance_sha256=sha(spatial_prov),
        preliminary_matrix_sha256=sha(preliminary_file),
        joined_matrix_sha256=sha(joined_path),status_matrix_sha256=sha(status_result)))
    save_meta(out/'orchestration_provenance.json',dict(
        stage='STAGE12_13_PARTIAL_HISTORICAL_LOCAL_ANALYSIS',
        input_stage11_manifest_sha256=sha(state['manifest']),
        input_stage10_ranking_sha256=sha(state['ranking']),
        input_stage11_status_sha256=sha(state['status']),
        conservation_original_source_sha256=sha(ORIGINAL),
        spatial_original_source_sha256=sha(SPATIAL),
        selected_families=len(state['rows']),stage12_conservation_completed_families=len(ids),
        stage12_spatial_completed_families=len(ids),stage13_preliminary_evidence_families=len(ids),
        pending_or_fallback_families=len(state['rows'])-len(ids),
        historical_reuse_families=sum(x['status']==HISTORICAL for x in computed),
        gpu_prediction_performed=False,foldseek_search_performed=False,
        stage13_msa_external_search_performed=False,
        stage13_full_integration_completed=False,final_top5_top3_selected=False,
        stage10_scores_recalculated=False,
        artifacts={str(path.relative_to(out)):{'sha256':sha(path)} for path in out.rglob('*') if path.is_file()}))
    print(f'SELECTED_FAMILIES={len(state["rows"])}')
    print(f'STAGE12_CONSERVATION_COMPLETED={len(ids)}')
    print(f'STAGE12_SPATIAL_COMPLETED={len(ids)}')
    print(f'STAGE13_PRELIMINARY_EVIDENCE={len(ids)}')
    print(f'PENDING_OR_FALLBACK={len(state["rows"])-len(ids)}')
    print('STAGE13_EXTERNAL_MSA_SEARCH=NOT_PERFORMED')
    print('FINAL_TOP5_TOP3_SELECTED=NO')
    print('STAGE162_PARTIAL_LOCAL_EVIDENCE=PASS')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--adaptive-output',type=Path,required=True)
    ap.add_argument('--stage09-10-output',type=Path,required=True)
    ap.add_argument('--out-dir',type=Path,required=True)
    ap.add_argument('--mode',choices=['check','run'],default='check')
    args=ap.parse_args()
    state=validate_inputs(args)
    if args.mode=='check':
        print(f'SELECTED_FAMILIES={len(state["rows"])}')
        print(f'STAGE12_QC_ELIGIBLE={len(state["evaluated"])}')
        print(f'PENDING_OR_FALLBACK={len(state["pending"])}')
        print('CHECK_MODE_NO_WRITE=PASS')
        return
    original,spatial,preliminary,computed=compute(state)
    execute(state,original,spatial,preliminary,computed)


if __name__=='__main__':
    try:
        main()
    except (OSError,ValueError,RuntimeError,KeyError,TypeError,csv.Error) as exc:
        sys.exit('ERROR: '+str(exc))
