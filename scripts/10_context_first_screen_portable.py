#!/usr/bin/env python3
"""Stage 10A: cohort-independent, structure-optional context-first review table.

This is *not* Stage 10 I2 scoring, a biochemical function prediction, or a
completed Stage 09/10 structural ranking. Candidate inclusion is based on
existing Stage 06 darkness and independent Stage 07 noncandidate anchors.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

CATS=("nuclease","restriction","ssdna_binding","recombination","ligase","helicase","primase","polymerase","dna_methylation")
DARK=("A_DEEP_DARK","B_FUNCTION_UNRESOLVED","C_KNOWN_DARK_BRIDGE")
ALLCLASSES=set(DARK)|{"D_KNOWN_EVIDENCE"}
OUT_FIELDS=("context_review_order","structure_id","family_rep","priority_class","family_size","sample_span","best_noncandidate_anchor_category","category_supported_members","category_supported_samples","all_noncandidate_anchor_categories","noncandidate_context_sample_span","noncandidate_anchor_context_count","all_context_sample_span","all_anchor_context_count","dominant_anchor_category","structure_status","prediction_selection","screening_status","screening_interpretation")


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20), b''):
            h.update(block)
    return h.hexdigest()


def read(path, needed):
    if not path.is_file(): raise ValueError(f'MISSING_INPUT:{path}')
    with path.open('r',encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)) or not set(needed)<=set(fields):
            raise ValueError(f'INVALID_SCHEMA:{path}:missing={sorted(set(needed)-set(fields))}')
        rows=list(reader)
    if not rows or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(f'INVALID_ROWS:{path}')
    return rows


def unique(rows,key):
    d={}
    for r in rows:
        if not r[key] or r[key] in d:raise ValueError(f'DUPLICATE_OR_EMPTY:{key}:{r[key]}')
        d[r[key]]=r
    return d


def num(v,field):
    try:n=int(v)
    except (ValueError,TypeError) as e:raise ValueError(f'INVALID_NUMBER:{field}:{v}') from e
    if n<0:raise ValueError(f'NEGATIVE_NUMBER:{field}')
    return n


def build(args):
    plan=args.partial_plan_dir
    inventory=plan/'partial_structure_inventory.tsv'
    manifest=plan/'partial_structure_plan_provenance.json'
    if not manifest.is_file():raise ValueError('PARTIAL_PLAN_MANIFEST_MISSING')
    meta=json.loads(manifest.read_text(encoding='utf-8'))
    if meta.get('gpu_prediction_performed') is not False or meta.get('external_search_performed') is not False:
        raise ValueError('PARTIAL_PLAN_STATUS_UNEXPECTED')
    expected=meta.get('output_sha256',{}).get('partial_structure_inventory.tsv')
    if expected!=sha(inventory):raise ValueError('STAGE179_INVENTORY_SHA_MISMATCH')
    source_family_sha=meta.get('input_sha256',{}).get('family_context')
    if source_family_sha!=sha(args.family_context):
        raise ValueError('STAGE179_STAGE07_FAMILY_SHA_MISMATCH')
    inventory_rows=read(inventory,('structure_id','family_rep','priority_class','family_size','sample_span','noncandidate_context_sample_span','all_context_sample_span','noncandidate_anchor_context_count','all_anchor_context_count','dominant_anchor_category','structure_status','prediction_selection'))
    member_rows=read(args.member_context,('family_rep','member','sample','candidate_original','anchor_within_window','anchor_categories','window'))
    family_rows=read(args.family_context,('family_rep','family_size','sample_span','all_context_sample_span','noncandidate_context_sample_span','noncandidate_anchor_context_count','all_anchor_context_count','dominant_anchor_category'))
    inv=unique(inventory_rows,'family_rep'); fam=unique(family_rows,'family_rep')
    if set(inv)!=set(fam) or len({r['structure_id'] for r in inventory_rows})!=len(inv):
        raise ValueError('COVERAGE_OR_ID_MISMATCH')
    by_family=defaultdict(list); member_seen=set(); windows=set()
    for r in member_rows:
        if r['member'] in member_seen or not r['member']:raise ValueError('DUPLICATE_MEMBER')
        member_seen.add(r['member'])
        if r['family_rep'] not in inv:raise ValueError('UNKNOWN_MEMBER_FAMILY')
        if r['candidate_original'] not in ('0','1') or r['anchor_within_window'] not in ('0','1'):
            raise ValueError('INVALID_MEMBER_STATUS')
        windows.add(num(r['window'],'window'))
        by_family[r['family_rep']].append(r)
    if len(windows)!=1:raise ValueError('INCONSISTENT_CONTEXT_WINDOWS')
    output=[]
    for rep,r in inv.items():
        members=by_family[rep]; f=fam[rep]
        if not members:raise ValueError(f'MISSING_MEMBERS:{rep}')
        if r['priority_class'] not in ALLCLASSES:raise ValueError('UNKNOWN_DARKNESS_CLASS')
        n=num(r['family_size'],'family_size')
        s={x['sample'] for x in members}
        if n!=len(members) or n!=num(f['family_size'],'family_size'):
            raise ValueError(f'FAMILY_SIZE_MISMATCH:{rep}')
        if len(s)!=num(r['sample_span'],'sample_span') or len(s)!=num(f['sample_span'],'sample_span'):
            raise ValueError(f'SAMPLE_SPAN_MISMATCH:{rep}')
        anchor=[m for m in members if m['anchor_within_window']=='1']
        non=[m for m in anchor if m['candidate_original']=='0']
        non_s={m['sample'] for m in non}; anchor_s={m['sample'] for m in anchor}
        for key,val in (('noncandidate_context_sample_span',len(non_s)),('all_context_sample_span',len(anchor_s)),('noncandidate_anchor_context_count',len(non)),('all_anchor_context_count',len(anchor))):
            if num(r[key],key)!=val or num(f[key],key)!=val:
                raise ValueError(f'CONTEXT_COUNT_MISMATCH:{rep}:{key}')
        if r['dominant_anchor_category']!=f['dominant_anchor_category']:
            raise ValueError('DOMINANT_ANCHOR_MISMATCH')
        catmembers=defaultdict(set); catsamples=defaultdict(set)
        for m in non:
            # The portable Stage07 format is semicolon-separated; unlike
            # historical H3's comma-separated group notation.
            categories=set(filter(None,m['anchor_categories'].split(';')))
            for category in categories:
                if category not in CATS:continue
                catmembers[category].add(m['member']);catsamples[category].add(m['sample'])
        # No biological activity probability/score assigned.
        best=sorted((c for c in CATS if c in catmembers),key=lambda c:(-len(catsamples[c]),-len(catmembers[c]),CATS.index(c)))
        category=best[0] if best else ''
        is_dark=r['priority_class'] in DARK
        if is_dark and category:
            status='DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT'
            interpretation='CONTEXT_ASSOCIATION_NOT_CONFIRMED_ENZYME_ACTIVITY'
        elif is_dark:
            status='DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY'
            interpretation='KEEP_UNRESOLVED_NOT_NEGATIVE'
        else:
            status='KNOWN_EVIDENCE_FAMILY'
            interpretation='NOT_DARK_PRIORITY_SCREEN'
        record={
            'context_review_order':'','structure_id':r['structure_id'],'family_rep':rep,'priority_class':r['priority_class'],
            'family_size':str(n),'sample_span':str(len(s)),
            'best_noncandidate_anchor_category':category,
            'category_supported_members':str(len(catmembers[category])) if category else '',
            'category_supported_samples':str(len(catsamples[category])) if category else '',
            'all_noncandidate_anchor_categories':';'.join(best),
            'noncandidate_context_sample_span':r['noncandidate_context_sample_span'],
            'noncandidate_anchor_context_count':r['noncandidate_anchor_context_count'],
            'all_context_sample_span':r['all_context_sample_span'],
            'all_anchor_context_count':r['all_anchor_context_count'],
            'dominant_anchor_category':r['dominant_anchor_category'],
            'structure_status':r['structure_status'],'prediction_selection':r['prediction_selection'],
            'screening_status':status,'screening_interpretation':interpretation,
        }
        output.append(record)
    # Sorting for expert review only, not a calibrated function or activity score.
    def sortkey(r):
        active=r['screening_status']=='DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT'
        return (not active,-num(r['category_supported_samples'] or '0','catspan'),-num(r['category_supported_members'] or '0','catmembers'),-num(r['noncandidate_context_sample_span'],'nonspan'),-num(r['sample_span'],'samplespan'),r['structure_id'])
    output.sort(key=sortkey)
    for i,r in enumerate(output,1):r['context_review_order']=str(i)
    candidate=[r for r in output if r['screening_status']=='DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT']
    return output,candidate,windows.pop(),meta


def write(path,rows):
    with path.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=OUT_FIELDS,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--partial-plan-dir',type=Path,required=True)
    p.add_argument('--member-context',type=Path,required=True)
    p.add_argument('--family-context',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    p.add_argument('--mode',choices=('check','run'),default='check')
    a=p.parse_args()
    if a.out_dir.exists():raise ValueError(f'REFUSE_EXISTING_OUTPUT:{a.out_dir}')
    if not a.partial_plan_dir.is_dir():raise ValueError('PLAN_DIRECTORY_MISSING')
    output,candidate,window,meta=build(a)
    by_status=Counter(r['screening_status'] for r in output)
    by_class=Counter(r['priority_class'] for r in output)
    if len(output)!=len({r['structure_id'] for r in output}):raise ValueError('DUPLICATE_ID_OUTPUT')
    print(f'INPUT_FAMILIES={len(output)}')
    print(f'DARK_CONTEXT_SUPPORTED_FAMILIES={len(candidate)}')
    print(f'DARK_WITHOUT_ENZYME_CONTEXT={by_status["DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY"]}')
    print(f'KNOWN_EVIDENCE_FAMILIES={by_status["KNOWN_EVIDENCE_FAMILY"]}')
    print('DARKNESS_CLASSES='+json.dumps(dict(sorted(by_class.items())),ensure_ascii=False))
    print(f'GENOMIC_CONTEXT_WINDOW={window}')
    if a.mode=='check':
        print('STAGE180_CONTEXT_FIRST_CHECK=PASS');print('FILES_WRITTEN=NO');return
    a.out_dir.mkdir(parents=True)
    write(a.out_dir/'all_family_context_screen.tsv',output)
    write(a.out_dir/'dark_replication_context_review.tsv',candidate)
    outputs={k:sha(a.out_dir/k) for k in ('all_family_context_screen.tsv','dark_replication_context_review.tsv')}
    provenance={
        'analysis_type':'STRUCTURE_OPTIONAL_CONTEXT_FIRST_SCREEN',
        'not_stage10_i2_ranking':True,'not_confirmed_enzyme_function':True,
        'ranking_is_review_order_not_activity_score':True,
        'all_families':len(output),'context_supported_dark_families':len(candidate),
        'darkness_class_counts':dict(sorted(by_class.items())),
        'context_window':window,
        'input_sha256':{'partial_inventory':sha(a.partial_plan_dir/'partial_structure_inventory.tsv'),'partial_plan_provenance':sha(a.partial_plan_dir/'partial_structure_plan_provenance.json'),'member_context':sha(a.member_context),'family_context':sha(a.family_context)},
        'output_sha256':outputs,
        'gpu_prediction_performed':False,'foldseek_search_performed':False,
        'stage10_i2_completed':False,'stage13_completed':False,
        'sorting_rule':'dark-supported first; best-category independent-sample span descending; supported-member count descending; noncandidate context sample span desc; family sample span desc; structure_id asc',
        'annotation_note':'Using independent noncandidate anchors within Stage07 window; this does not prove the protein encodes the nearby enzymatic function',
    }
    (a.out_dir/'context_first_provenance.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False)+'\n')
    print('STAGE180_CONTEXT_FIRST_SCREEN=PASS')
    print(f'OUTPUT_DIR={a.out_dir}')
    print('GPU_OR_FOLDSEEK=NOT_EXECUTED')
    print('STAGE10_FINAL_I2_RANKING=NOT_EXECUTED')

if __name__=='__main__':main()
