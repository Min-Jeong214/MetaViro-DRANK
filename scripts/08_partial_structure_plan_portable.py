#!/usr/bin/env python3
"""Plan structure-aware follow-up for new cohorts without requiring all PDBs.

Read-only source validation; writes only a new output directory in --mode run.
No biological activity is inferred, and this is not Stage 10 final ranking.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = [
    "structure_id", "family_rep", "priority_class", "family_size",
    "sample_span", "noncandidate_context_sample_span", "all_context_sample_span",
    "all_anchor_context_count", "noncandidate_anchor_context_count",
    "dominant_anchor_category", "length_aa", "sequence_sha256",
    "structure_status", "historical_structure_id", "historical_pdb_file",
    "historical_pdb_sha256", "historical_score_json",
    "historical_score_json_sha256", "historical_mean_plddt", "historical_ptm",
    "historical_match_count", "prediction_selection", "prediction_priority_order",
]


def sha_file(p):
    h = hashlib.sha256()
    with p.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path, required):
    if not path.is_file():
        raise ValueError(f"MISSING_FILE:{path}")
    with path.open(encoding='utf-8-sig', newline='') as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        fields = rd.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"DUPLICATE_COLUMNS:{path}")
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"MISSING_COLUMNS:{path}:{sorted(missing)}")
        rows = list(rd)
    if not rows:
        raise ValueError(f"EMPTY_TABLE:{path}")
    return rows


def unique(rows, field):
    by = {}
    for r in rows:
        key = r[field]
        if not key or key in by:
            raise ValueError(f"EMPTY_OR_DUPLICATE_{field}:{key}")
        by[key] = r
    return by


def positive_int(value, key):
    try:
        x = int(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f"INVALID_INTEGER:{key}:{value}") from e
    if x < 0:
        raise ValueError(f"NEGATIVE_INTEGER:{key}")
    return x


def fasta(path):
    if not path.is_file():
        raise ValueError(f"MISSING_FASTA:{path}")
    seqs = {}
    pid = None
    chunks = []
    def emit():
        if pid is not None:
            if pid in seqs:
                raise ValueError(f"DUPLICATE_FASTA_ID:{pid}")
            sequence = ''.join(chunks).upper().rstrip('*')
            if not sequence or any(not ('A' <= c <= 'Z') for c in sequence):
                raise ValueError(f"INVALID_SEQUENCE:{pid}")
            seqs[pid] = sequence
    with path.open(encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                emit()
                pid = line[1:].split()[0]
                if not pid:
                    raise ValueError('EMPTY_FASTA_ID')
                chunks = []
            else:
                if pid is None:
                    raise ValueError('SEQUENCE_BEFORE_FASTA_HEADER')
                chunks.append(line)
    emit()
    return seqs


def matching_structures(rows):
    index = defaultdict(list)
    if rows is None:
        return index
    ids = unique(rows, 'structure_id')
    for sid, r in ids.items():
        if r['prediction_status'] != 'PREDICTED':
            continue
        seqsha = r['sequence_sha256']
        if len(seqsha) != 64:
            raise ValueError(f"INVALID_OLD_SEQUENCE_SHA:{sid}")
        index[seqsha].append(r)
    return index


def verify_old(row, new):
    sid = row['structure_id']
    if row['sequence_sha256'] != new['sequence_sha256']:
        raise ValueError(f"OLD_SEQUENCE_MISMATCH:{sid}")
    if positive_int(row['length_aa'], 'old_length') != positive_int(new['length_aa'], 'new_length'):
        raise ValueError(f"OLD_LENGTH_MISMATCH:{sid}")
    for file_key, hash_key in (('pdb_file','pdb_sha256'),('score_json','score_json_sha256')):
        p = Path(row[file_key])
        if not p.is_file() or not p.stat().st_size:
            raise ValueError(f"OLD_SOURCE_MISSING:{sid}:{file_key}")
        if sha_file(p) != row[hash_key]:
            raise ValueError(f"OLD_SOURCE_SHA_MISMATCH:{sid}:{file_key}")


def write_tsv(path, rows, fields):
    with path.open('x', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=fields, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def make_plan(args):
    required_map = ('structure_id', 'family_rep', 'family_size', 'length_aa', 'sequence_sha256', 'priority_class')
    required_audit = ('representative','family_size','priority_class')
    required_context = (
        'family_rep', 'family_size', 'sample_span', 'noncandidate_context_sample_span',
        'all_context_sample_span','all_anchor_context_count', 'noncandidate_anchor_context_count',
        'dominant_anchor_category',
    )
    mapping = read_tsv(args.mapping, required_map)
    audit = unique(read_tsv(args.family_audit, required_audit), 'representative')
    context = unique(read_tsv(args.family_context, required_context), 'family_rep')
    seqs = fasta(args.queries)
    mids = unique(mapping, 'structure_id')
    reps = unique(mapping, 'family_rep')
    if set(seqs) != set(mids) or set(reps) != set(audit) or set(reps) != set(context):
        raise ValueError('STAGE06_07_08_COVERAGE_MISMATCH')
    if args.expect_families is not None and len(mapping) != args.expect_families:
        raise ValueError('UNEXPECTED_FAMILY_COUNT')
    old_rows = None
    if args.historical_structures:
        old_rows = read_tsv(args.historical_structures, (
            'structure_id', 'sequence_sha256', 'length_aa', 'prediction_status',
            'pdb_file','pdb_sha256','score_json','score_json_sha256',
            'mean_plddt','ptm',
        ))
    old_by_sha = matching_structures(old_rows)
    output = []
    for sid in sorted(mids):
        m = mids[sid]
        rep = m['family_rep']
        a, c = audit[rep], context[rep]
        for r in (a, c):
            if positive_int(r['family_size'], 'family_size') != positive_int(m['family_size'], 'mapping_size'):
                raise ValueError(f"FAMILY_SIZE_MISMATCH:{sid}")
        if a['priority_class'] != m['priority_class']:
            raise ValueError(f"PRIORITY_CLASS_MISMATCH:{sid}")
        if a['priority_class'] not in {'A_DEEP_DARK','B_FUNCTION_UNRESOLVED','C_KNOWN_DARK_BRIDGE','D_KNOWN_EVIDENCE'}:
            raise ValueError(f"UNKNOWN_PRIORITY_CLASS:{sid}")
        seq = seqs[sid]
        if len(seq) != positive_int(m['length_aa'], 'length_aa'):
            raise ValueError(f"QUERY_LENGTH_MISMATCH:{sid}")
        if hashlib.sha256(seq.encode('ascii')).hexdigest() != m['sequence_sha256']:
            raise ValueError(f"QUERY_SEQUENCE_SHA_MISMATCH:{sid}")
        for k in ('sample_span','noncandidate_context_sample_span','all_context_sample_span',
                  'all_anchor_context_count', 'noncandidate_anchor_context_count'):
            positive_int(c[k], k)
        if not (0 <= int(c['noncandidate_context_sample_span']) <= int(c['all_context_sample_span']) <= int(c['sample_span'])):
            raise ValueError(f"INVALID_SAMPLE_SPAN:{sid}")
        if int(c['noncandidate_anchor_context_count']) > int(c['all_anchor_context_count']):
            raise ValueError(f"INVALID_ANCHOR_COUNT:{sid}")
        matches = sorted(old_by_sha.get(m['sequence_sha256'], []), key=lambda r:r['structure_id'])
        for old in matches:
            verify_old(old,m)
        reuse = matches[0] if matches else None
        output.append({
            'structure_id':sid, 'family_rep':rep, 'priority_class':a['priority_class'],
            'family_size':m['family_size'], 'sample_span':c['sample_span'],
            'noncandidate_context_sample_span':c['noncandidate_context_sample_span'],
            'all_context_sample_span':c['all_context_sample_span'],
            'all_anchor_context_count':c['all_anchor_context_count'],
            'noncandidate_anchor_context_count':c['noncandidate_anchor_context_count'],
            'dominant_anchor_category':c['dominant_anchor_category'],
            'length_aa':m['length_aa'], 'sequence_sha256':m['sequence_sha256'],
            'structure_status':'VERIFIED_HISTORICAL_REUSE_CANDIDATE' if reuse else 'NO_STRUCTURE_AVAILABLE',
            'historical_structure_id':reuse['structure_id'] if reuse else '',
            'historical_pdb_file':str(Path(reuse['pdb_file']).resolve()) if reuse else '',
            'historical_pdb_sha256':reuse['pdb_sha256'] if reuse else '',
            'historical_score_json':str(Path(reuse['score_json']).resolve()) if reuse else '',
            'historical_score_json_sha256':reuse['score_json_sha256'] if reuse else '',
            'historical_mean_plddt':reuse['mean_plddt'] if reuse else '',
            'historical_ptm':reuse['ptm'] if reuse else '',
            'historical_match_count':len(matches),
            'prediction_selection':'', 'prediction_priority_order':'',
        })

    dark_classes = {'A_DEEP_DARK','B_FUNCTION_UNRESOLVED','C_KNOWN_DARK_BRIDGE'}
    missing_dark = [r for r in output if r['priority_class'] in dark_classes and r['structure_status']=='NO_STRUCTURE_AVAILABLE']
    # Transparent heuristic for resource allocation; NOT a calibrated functional score.
    missing_dark.sort(key=lambda r: (
        -int(r['noncandidate_context_sample_span']),
        -int(r['all_context_sample_span']),
        -int(r['noncandidate_anchor_context_count']),
        -int(r['all_anchor_context_count']),
        -int(r['sample_span']),
        r['family_rep'],
    ))
    shortlist = missing_dark[:args.limit]
    priority = {r['structure_id']:i+1 for i,r in enumerate(shortlist)}
    for r in output:
        sid = r['structure_id']
        if r['structure_status']=='VERIFIED_HISTORICAL_REUSE_CANDIDATE':
            r['prediction_selection']='REUSE_CANDIDATE_NOT_NEW_PREDICTION'
        elif sid in priority:
            r['prediction_selection']='SELECTED_FOR_REVIEW_NOT_AUTOMATIC_GPU_RUN'
            r['prediction_priority_order']=priority[sid]
        elif r['priority_class']=='D_KNOWN_EVIDENCE':
            r['prediction_selection']='DEFERRED_KNOWN_FUNCTION_EVIDENCE'
        else:
            r['prediction_selection']='DEFERRED_NO_STRUCTURE'
    by_id = {r['structure_id']:r for r in output}
    return output, [by_id[r['structure_id']] for r in shortlist], seqs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mapping',type=Path,required=True)
    p.add_argument('--queries',type=Path,required=True)
    p.add_argument('--family-audit',type=Path,required=True)
    p.add_argument('--family-context',type=Path,required=True)
    p.add_argument('--historical-structures',type=Path)
    p.add_argument('--limit',type=int,default=30)
    p.add_argument('--expect-families',type=int)
    p.add_argument('--out-dir',type=Path,required=True)
    p.add_argument('--mode',choices=('check','run'),default='check')
    args=p.parse_args()
    if args.limit < 1:
        p.error('--limit must be positive')
    if args.out_dir.exists():
        raise ValueError(f"REFUSE_OVERWRITE:{args.out_dir}")
    output,shortlist,seqs=make_plan(args)
    counts=Counter(r['structure_status'] for r in output)
    dark = sum(r['priority_class']!='D_KNOWN_EVIDENCE' for r in output)
    remaining_dark=sum(r['priority_class']!='D_KNOWN_EVIDENCE' and r['structure_status']=='NO_STRUCTURE_AVAILABLE' for r in output)
    print(f"ALL_FAMILIES={len(output)}")
    print(f"STRUCTURE_REUSE_CANDIDATES={counts['VERIFIED_HISTORICAL_REUSE_CANDIDATE']}")
    print(f"NO_STRUCTURE_AVAILABLE={counts['NO_STRUCTURE_AVAILABLE']}")
    print(f"A_TO_C_FAMILIES={dark}")
    print(f"A_TO_C_NO_STRUCTURE={remaining_dark}")
    print(f"NEW_PREDICTION_REVIEW_SHORTLIST={len(shortlist)}")
    if args.mode=='check':
        print('STEP179_PARTIAL_STRUCTURE_PLAN_CHECK=PASS')
        print('FILES_WRITTEN=NO')
        return
    args.out_dir.mkdir(parents=True)
    write_tsv(args.out_dir/'partial_structure_inventory.tsv',output,FIELDS)
    write_tsv(args.out_dir/'new_prediction_review_shortlist.tsv',shortlist,FIELDS)
    write_tsv(args.out_dir/'historical_reuse_candidates.tsv',[
        r for r in output if r['structure_status']=='VERIFIED_HISTORICAL_REUSE_CANDIDATE'
    ],FIELDS)
    with (args.out_dir/'new_prediction_review_queries.faa').open('x', encoding='utf-8') as fh:
        for row in shortlist:
            sid=row['structure_id']; seq=seqs[sid]
            fh.write(f'>{sid}\n')
            for i in range(0,len(seq),80):
                fh.write(seq[i:i+80]+'\n')
    provenance={
        'mode':'READ_ONLY_SOURCE_TRIAGE',
        'not_final_function_ranking':True,
        'external_search_performed':False,
        'gpu_prediction_performed':False,
        'structure_reuse_status':'CANDIDATE_REQUIRES_DOWNSTREAM_STAGE08_IMPORT',
        'ranking_rule':'noncandidate_context_sample_span desc; all_context_sample_span desc; noncandidate_anchor_context_count desc; all_anchor_context_count desc; sample_span desc; family_rep asc',
        'count_all_families':len(output),
        'count_verified_historical_reuse_candidates':counts['VERIFIED_HISTORICAL_REUSE_CANDIDATE'],
        'count_no_structure':counts['NO_STRUCTURE_AVAILABLE'],
        'count_dark_without_structure':remaining_dark,
        'count_prediction_review_shortlist':len(shortlist),
        'input_sha256':{key:sha_file(val) for key,val in (
            ('mapping',args.mapping),('queries',args.queries),('family_audit',args.family_audit),('family_context',args.family_context)
        )},
        'historical_index_sha256':sha_file(args.historical_structures) if args.historical_structures else None,
        'output_sha256':{fn:sha_file(args.out_dir/fn) for fn in (
            'partial_structure_inventory.tsv','new_prediction_review_shortlist.tsv',
            'historical_reuse_candidates.tsv','new_prediction_review_queries.faa'
        )},
    }
    (args.out_dir/'partial_structure_plan_provenance.json').write_text(
        json.dumps(provenance,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('STEP179_PARTIAL_STRUCTURE_PLAN=PASS')
    print(f'OUTPUT_DIR={args.out_dir}')
    print('GPU_PREDICTION=NOT_RUN')
    print('FOLDSEEK_SEARCH=NOT_RUN')
    print('FINAL_FUNCTION_RANKING=NOT_RUN')

if __name__=='__main__':
    main()
