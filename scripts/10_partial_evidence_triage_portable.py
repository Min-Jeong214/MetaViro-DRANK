#!/usr/bin/env python3
"""Structure-optional Stage10 review: all families, no invented I2/biochemical scores.

Grounded in Stage180 independent noncandidate context and verified Stage182
Foldseek evidence. The output is a deterministic *review order*, not a calibrated
probability, nuclease activity claim, or complete historical I2 ranking.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

DARK = {'A_DEEP_DARK', 'B_FUNCTION_UNRESOLVED', 'C_KNOWN_DARK_BRIDGE'}
KNOWN = 'D_KNOWN_EVIDENCE'
CTX_YES = 'DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT'
CTX_NO = 'DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY'
CTX_KNOWN = 'KNOWN_EVIDENCE_FAMILY'
STAGE09_STATES = {'SEARCH_VERIFIED', 'NOT_EVALUATED_NO_STRUCTURE'}
SEARCH_RESULTS = {'CONVINCING_STRUCTURE_CONNECTION', 'WEAK_RAW_HIT_ONLY', 'NO_RAW_HIT_E1'}

CTX_FIELDS = (
    'context_review_order', 'structure_id', 'family_rep', 'priority_class',
    'family_size', 'sample_span', 'best_noncandidate_anchor_category',
    'category_supported_members', 'category_supported_samples',
    'noncandidate_context_sample_span', 'noncandidate_anchor_context_count',
    'screening_status', 'screening_interpretation',
)
STAGE09_FIELDS = (
    'structure_id', 'family_rep', 'priority_class', 'family_size',
    'structure_state', 'external_search_state', 'stage09_evidence_state',
    'stage09_external_status', 'stage09_raw_sources',
    'stage09_convincing_sources',
)
OUT_FIELDS = (
    'review_order', 'review_tier', 'structure_id', 'family_rep',
    'priority_class', 'family_size', 'sample_span',
    'genomic_context_state', 'best_noncandidate_anchor_category',
    'category_supported_samples', 'category_supported_members',
    'noncandidate_context_sample_span', 'noncandidate_anchor_context_count',
    'structure_state', 'stage09_evidence_state', 'stage09_external_status',
    'stage09_raw_sources', 'stage09_convincing_sources',
    'stage09_interpretation', 'biochemical_activity_confirmed',
    'historical_I2_score', 'review_interpretation',
)
TIER_ORDER = (
    'T1_CONTEXT_AND_CONVINCING_STRUCTURE',
    'T2_CONTEXT_SUPPORTED_DARK',
    'T3_CONVINCING_STRUCTURE_NO_SUPPORTED_CONTEXT',
    'T4_DARK_CONTEXT_UNRESOLVED',
    'T5_KNOWN_EVIDENCE_REFERENCE',
)


def sha(file):
    h = hashlib.sha256()
    with Path(file).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_table(file, required):
    if not file.is_file() or not file.stat().st_size:
        raise ValueError(f'MISSING_INPUT:{file}')
    with file.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(required).issubset(fields):
            raise ValueError(f'INVALID_SCHEMA:{file}:missing={sorted(set(required)-set(fields))}')
        rows = list(reader)
    if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f'INVALID_ROWS:{file}')
    return rows


def unique(rows, name):
    out = {}
    for row in rows:
        value = row[name]
        if not value or value in out:
            raise ValueError(f'EMPTY_OR_DUPLICATE:{name}:{value}')
        out[value] = row
    return out


def positive(row, key, allow_zero=False):
    try:
        value = int(row[key])
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError(f'INVALID_INTEGER:{key}:{row.get(key)}') from exc
    if value < (0 if allow_zero else 1):
        raise ValueError(f'INVALID_INTEGER_RANGE:{key}:{value}')
    return value


def load_json(file):
    if not file.is_file():
        raise ValueError(f'MISSING_PROVENANCE:{file}')
    data = json.loads(file.read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise ValueError(f'INVALID_PROVENANCE:{file}')
    return data


def check_origin(args, context_file, evidence_file):
    context_manifest_file = args.context_dir / 'context_first_provenance.json'
    stage09_manifest_file = args.stage09_dir / 'stage182_merge_provenance.json'
    prep_manifest_file = args.prepared_dir / 'stage182_preparation_provenance.json'
    context = load_json(context_manifest_file)
    stage09 = load_json(stage09_manifest_file)
    prep = load_json(prep_manifest_file)
    if context.get('analysis_type') != 'STRUCTURE_OPTIONAL_CONTEXT_FIRST_SCREEN' or not context.get('not_stage10_i2_ranking'):
        raise ValueError('STAGE180_NOT_PRELIMINARY_CONTEXT')
    if context.get('output_sha256', {}).get(context_file.name) != sha(context_file):
        raise ValueError('STAGE180_CONTEXT_SHA_MISMATCH')
    if stage09.get('stage') != 'STAGE09_PARTIAL_EVIDENCE_MERGED' or stage09.get('output_sha256') != sha(evidence_file):
        raise ValueError('STAGE182_EVIDENCE_SHA_MISMATCH')
    if stage09.get('source_preparation_manifest_sha256') != sha(prep_manifest_file):
        raise ValueError('STAGE182_PREPARATION_PROVENANCE_MISMATCH')
    inventory_sha = context.get('input_sha256', {}).get('partial_inventory')
    if not inventory_sha or inventory_sha != prep.get('source_inventory_sha256'):
        raise ValueError('CONTEXT_STAGE09_DIFFERENT_INVENTORY')
    if prep.get('source_plan_provenance_sha256') != context.get('input_sha256', {}).get('partial_plan_provenance'):
        raise ValueError('CONTEXT_STAGE09_DIFFERENT_PLAN')
    status_file = args.prepared_dir / 'all_family_structure_status.tsv'
    models_file = args.prepared_dir / 'structure_results_for_search.tsv'
    if prep.get('all_family_status_sha256') != sha(status_file) or prep.get('input_structure_sha256') != sha(models_file):
        raise ValueError('PREPARED_STRUCTURE_INPUT_SHA_MISMATCH')
    return context, stage09, prep, (context_manifest_file, stage09_manifest_file, prep_manifest_file, status_file, models_file)


def build(args):
    ctxfile = args.context_dir / 'all_family_context_screen.tsv'
    evfile = args.stage09_dir / 'all_family_stage09_evidence.tsv'
    context_meta, evidence_meta, prep_meta, input_files = check_origin(args, ctxfile, evfile)
    ctxrows = read_table(ctxfile, CTX_FIELDS)
    evrows = read_table(evfile, STAGE09_FIELDS)
    c = unique(ctxrows, 'structure_id')
    e = unique(evrows, 'structure_id')
    if set(c) != set(e) or len(c) != context_meta.get('all_families') or len(e) != evidence_meta.get('all_family_count'):
        raise ValueError('CONTEXT_STAGE09_FAMILY_SET_MISMATCH')
    if len({r['family_rep'] for r in ctxrows}) != len(ctxrows):
        raise ValueError('DUPLICATE_CONTEXT_FAMILY_REP')
    if len({r['family_rep'] for r in evrows}) != len(evrows):
        raise ValueError('DUPLICATE_STAGE09_FAMILY_REP')
    statusrows = read_table(args.prepared_dir / 'all_family_structure_status.tsv',
                            ('structure_id', 'family_rep', 'structure_state', 'external_search_state'))
    statusby = unique(statusrows, 'structure_id')
    if set(statusby) != set(c):
        raise ValueError('PREPARATION_STATUS_COVERAGE_MISMATCH')
    if set(r['stage09_evidence_state'] for r in evrows) - STAGE09_STATES:
        raise ValueError('UNKNOWN_STAGE09_STATE')
    result = []
    for sid, x in c.items():
        y = e[sid]
        z = statusby[sid]
        for key in ('family_rep', 'priority_class', 'family_size'):
            if x[key] != y[key]:
                raise ValueError(f'CROSS_STAGE_FIELD_MISMATCH:{sid}:{key}')
        for key in ('family_rep', 'structure_state', 'external_search_state'):
            if y[key] != z[key]:
                raise ValueError(f'PREPARATION_FIELD_MISMATCH:{sid}:{key}')
        n = positive(x, 'family_size')
        sp = positive(x, 'sample_span')
        if sp > n:
            raise ValueError(f'SPAN_EXCEEDS_FAMILY_SIZE:{sid}')
        is_dark = x['priority_class'] in DARK
        if not is_dark and x['priority_class'] != KNOWN:
            raise ValueError(f'UNKNOWN_DARKNESS_CLASS:{sid}')
        ctx_supported = x['screening_status'] == CTX_YES
        expected = (CTX_YES if ctx_supported else CTX_NO) if is_dark else CTX_KNOWN
        if x['screening_status'] != expected:
            raise ValueError(f'DARKNESS_CONTEXT_INCONSISTENT:{sid}')
        cat = x['best_noncandidate_anchor_category']
        cat_n = positive(x, 'category_supported_samples', True) if x['category_supported_samples'] else 0
        mem_n = positive(x, 'category_supported_members', True) if x['category_supported_members'] else 0
        if ctx_supported and (not cat or not cat_n or not mem_n or cat_n > sp):
            raise ValueError(f'CONTEXT_SUPPORT_COUNT_INVALID:{sid}')
        if not cat and (cat_n or mem_n):
            raise ValueError(f'CONTEXT_CATEGORY_COUNT_WITHOUT_CATEGORY:{sid}')
        searched = y['stage09_evidence_state'] == 'SEARCH_VERIFIED'
        ext = y['stage09_external_status']
        if searched:
            if y['structure_state'] != 'PREDICTED_VERIFIED' or y['external_search_state'] != 'READY_FOR_SEARCH' or ext not in SEARCH_RESULTS:
                raise ValueError(f'INVALID_SEARCH_VERIFIED_STATUS:{sid}')
            structure_explanation = {
                'CONVINCING_STRUCTURE_CONNECTION': 'VERIFIED_STRUCTURAL_CONNECTION_NOT_BIOCHEMICAL_FUNCTION',
                'WEAK_RAW_HIT_ONLY': 'WEAK_STRUCTURAL_SIMILARITY_NOT_FUNCTION_PROOF',
                'NO_RAW_HIT_E1': 'SEARCH_PERFORMED_NO_RAW_HIT_AT_E1',
            }[ext]
        else:
            if (y['structure_state'] != 'PREDICTION_PENDING' or
                y['external_search_state'] != 'NOT_EVALUATED_NO_STRUCTURE' or
                ext or y['stage09_raw_sources'] or y['stage09_convincing_sources']):
                raise ValueError(f'PENDING_STRUCTURAL_EVIDENCE_NOT_BLANK:{sid}')
            structure_explanation = 'NOT_EVALUATED_NO_STRUCTURE_NOT_NEGATIVE'
        if searched and ext == 'CONVINCING_STRUCTURE_CONNECTION' and not y['stage09_convincing_sources']:
            raise ValueError(f'CONVINCING_MISSING_SOURCE:{sid}')
        strong = searched and ext == 'CONVINCING_STRUCTURE_CONNECTION'
        if not is_dark:
            tier = TIER_ORDER[4]
        elif ctx_supported and strong:
            tier = TIER_ORDER[0]
        elif ctx_supported:
            tier = TIER_ORDER[1]
        elif strong:
            tier = TIER_ORDER[2]
        else:
            tier = TIER_ORDER[3]
        result.append({
            'review_order': '', 'review_tier': tier, 'structure_id': sid,
            'family_rep': x['family_rep'], 'priority_class': x['priority_class'],
            'family_size': str(n), 'sample_span': str(sp),
            'genomic_context_state': x['screening_status'],
            'best_noncandidate_anchor_category': cat,
            'category_supported_samples': x['category_supported_samples'],
            'category_supported_members': x['category_supported_members'],
            'noncandidate_context_sample_span': x['noncandidate_context_sample_span'],
            'noncandidate_anchor_context_count': x['noncandidate_anchor_context_count'],
            'structure_state': y['structure_state'],
            'stage09_evidence_state': y['stage09_evidence_state'],
            'stage09_external_status': ext,
            'stage09_raw_sources': y['stage09_raw_sources'],
            'stage09_convincing_sources': y['stage09_convincing_sources'],
            'stage09_interpretation': structure_explanation,
            'biochemical_activity_confirmed': 'NOT_TESTED',
            'historical_I2_score': '',
            'review_interpretation': (
                'CONTEXT_PLUS_STRUCTURAL_HOMOLOGY_NOT_FUNCTION_PROOF' if tier == TIER_ORDER[0] else
                'INDEPENDENT_GENOMIC_CONTEXT_REVIEW' if tier == TIER_ORDER[1] else
                'STRUCTURAL_HOMOLOGY_ONLY_NOT_FUNCTION_PROOF' if tier == TIER_ORDER[2] else
                'DARK_UNRESOLVED_KEEP_FOR_FURTHER_EVIDENCE' if tier == TIER_ORDER[3] else
                'REFERENCE_KNOWN_FAMILY_NOT_DARK_PRIORITY'
            ),
        })
    # No numerical weight is assigned to weak/raw hits, missing structures,
    # or lack of a match. Deterministic *review* order, not scientific I2 score.
    def key(row):
        return (
            TIER_ORDER.index(row['review_tier']),
            -positive(row, 'category_supported_samples', True) if row['category_supported_samples'] else 0,
            -positive(row, 'category_supported_members', True) if row['category_supported_members'] else 0,
            -positive(row, 'noncandidate_context_sample_span', True),
            -positive(row, 'sample_span'),
            row['structure_id'],
        )
    result.sort(key=key)
    for i, row in enumerate(result, 1):
        row['review_order'] = str(i)
    nchecked = sum(r['stage09_evidence_state'] == 'SEARCH_VERIFIED' for r in result)
    nctx = sum(r['genomic_context_state'] == CTX_YES for r in result)
    if nchecked != evidence_meta.get('externally_evaluated_family_count') or nctx != context_meta.get('context_supported_dark_families'):
        raise ValueError('PROVENANCE_COUNTS_INCONSISTENT')
    if nchecked != prep_meta.get('search_ready_count'):
        raise ValueError('SEARCH_READY_COUNT_MISMATCH')
    return result, ctxfile, evfile, input_files


def write_tsv(file, rows):
    with file.open('x', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--context-dir', type=Path, required=True)
    p.add_argument('--stage09-dir', type=Path, required=True)
    p.add_argument('--prepared-dir', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    p.add_argument('--mode', choices=('check', 'run'), default='check')
    args = p.parse_args()
    if args.out_dir.exists():
        raise ValueError(f'REFUSE_EXISTING_OUTPUT:{args.out_dir}')
    target = args.out_dir.resolve()
    for source in (args.context_dir, args.stage09_dir, args.prepared_dir):
        src = source.resolve()
        if target == src or target in src.parents or src in target.parents:
            raise ValueError('SOURCE_OUTPUT_OVERLAP')
    rows, context, evidence, other = build(args)
    counts = dict(Counter(r['review_tier'] for r in rows))
    print(f'ALL_FAMILIES={len(rows)}')
    print('TIER_COUNTS=' + json.dumps(counts, sort_keys=True))
    print(f'FOLDSEEK_VERIFIED={sum(r["stage09_evidence_state"] == "SEARCH_VERIFIED" for r in rows)}')
    print(f'PREDICTION_PENDING={sum(r["stage09_evidence_state"] == "NOT_EVALUATED_NO_STRUCTURE" for r in rows)}')
    print(f'CONTEXT_SUPPORTED_DARK={sum(r["genomic_context_state"] == CTX_YES for r in rows)}')
    if args.mode == 'check':
        print('STEP183_STAGE10_REVIEW_CHECK=PASS')
        print('FILES_WRITTEN=NO')
        return
    args.out_dir.mkdir(parents=True, exist_ok=False)
    allout = args.out_dir / 'all_family_stage10_review.tsv'
    darkout = args.out_dir / 'dark_family_stage10_review.tsv'
    focusout = args.out_dir / 'context_supported_stage10_shortlist.tsv'
    write_tsv(allout, rows)
    write_tsv(darkout, [r for r in rows if r['priority_class'] in DARK])
    write_tsv(focusout, [r for r in rows if r['genomic_context_state'] == CTX_YES])
    manifest = {
        'stage': 'STAGE10_PARTIAL_EVIDENCE_REVIEW', 'schema_version': 1,
        'status': 'PRELIMINARY_EVIDENCE_REVIEW_NOT_I2_FINAL',
        'all_families': len(rows),
        'tier_counts': counts,
        'structural_evidence_not_evaluated_is_not_negative': True,
        'weak_raw_hit_is_not_positive_function_evidence': True,
        'biochemical_function_confirmed': False,
        'stage10_historical_I2_completed': False,
        'stage13_completed': False,
        'sort_order': 'tier; category_supported_samples desc; category_supported_members desc; noncandidate_context_sample_span desc; sample_span desc; structure_id asc',
        'source_sha256': {str(file): sha(file) for file in (context, evidence, *other)},
        'output_sha256': {file.name: sha(file) for file in (allout, darkout, focusout)},
    }
    (args.out_dir / 'stage183_review_provenance.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print('STEP183_STAGE10_PARTIAL_REVIEW=PASS')
    print(f'OUTPUT_DIR={args.out_dir}')
    print('I2_FINAL_RANKING=NOT_RUN')
    print('GPU_FOLDSEEK_REEXECUTED=NO')


if __name__ == '__main__':
    main()
