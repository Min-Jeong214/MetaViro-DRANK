#!/usr/bin/env python3
"""Stage 13: structure-optional, *provisional* evidence catalog for any cohort.

Merge verified Stage183 all-family review and Stage184 selected sequence-only
conservation. Does not infer biochemical function, catalytic sites, I2 scores,
structure-based Stage12, or historical Top5/Top3. No GPU or Foldseek jobs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

STAGE183_COLUMNS = (
    'review_order', 'review_tier', 'structure_id', 'family_rep', 'family_size',
    'priority_class', 'sample_span', 'genomic_context_state',
    'best_noncandidate_anchor_category', 'category_supported_samples',
    'category_supported_members', 'stage09_evidence_state',
    'stage09_external_status', 'biochemical_activity_confirmed',
    'historical_I2_score',
)
STAGE184_COLUMNS = (
    'structure_id', 'family_rep', 'priority_class', 'msa_plan_state',
    'msa_conservation_state', 'analyzable_positions', 'conserved_positions',
    'full_stage12_structure_analysis',
)
SUMMARY_COLUMNS = (
    'structure_id', 'family_rep', 'family_size', 'msa_unique_sequences',
    'query_length', 'analyzable_positions', 'conserved_positions',
    'msa_conservation_state', 'structure_residue_analysis',
)
RESIDUE_COLUMNS = (
    'structure_id', 'query_position', 'query_residue', 'non_gap_homologs',
    'matching_homologs', 'minimum_homolog_support_met', 'sequence_conserved',
)
DARK_TIERS = {
    'T1_CONTEXT_AND_CONVINCING_STRUCTURE',
    'T2_CONTEXT_SUPPORTED_DARK',
    'T3_CONVINCING_STRUCTURE_NO_SUPPORTED_CONTEXT',
}
TIERS = DARK_TIERS | {'T4_DARK_CONTEXT_UNRESOLVED', 'T5_KNOWN_EVIDENCE_REFERENCE'}
MSA_STATES = {'SEQUENCE_ONLY_CONSERVATION_AVAILABLE',
              'INSUFFICIENT_HOMOLOG_SUPPORT', 'SINGLE_UNIQUE_SEQUENCE'}
OUT_FIELDS = (
    'review_order', 'review_tier', 'structure_id', 'family_rep',
    'priority_class', 'family_size', 'sample_span',
    'genomic_context_state', 'best_noncandidate_anchor_category',
    'category_supported_samples', 'category_supported_members',
    'stage09_evidence_state', 'stage09_external_status',
    'msa_plan_state', 'msa_conservation_state', 'msa_unique_sequences',
    'msa_query_length', 'msa_analyzable_positions', 'msa_conserved_positions',
    'stage12_spatial_state', 'biochemical_activity_confirmed',
    'historical_I2_score', 'stage13_evidence_interpretation',
    'stage13_selection_state',
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def require(path: Path) -> Path:
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f'MISSING_OR_EMPTY:{path}')
    return path


def provenance(path: Path) -> dict:
    obj = json.loads(require(path).read_text(encoding='utf-8'))
    if not isinstance(obj, dict):
        raise ValueError(f'INVALID_PROVENANCE:{path}')
    return obj


def table(path: Path, columns: tuple | list):
    with require(path).open('r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        names = reader.fieldnames or []
        if len(set(names)) != len(names) or not set(columns).issubset(names):
            raise ValueError(f'INVALID_COLUMNS:{path}:{sorted(set(columns)-set(names))}')
        rows = list(reader)
    if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f'MALFORMED_TABLE:{path}')
    return rows


def index(rows, field):
    out = {}
    for row in rows:
        key = row[field]
        if not key or key in out:
            raise ValueError(f'DUPLICATE_OR_BLANK:{field}:{key}')
        out[key] = row
    return out


def nonneg(value, field):
    try:
        if not str(value).isdigit():
            raise ValueError()
        x = int(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f'INVALID_INTEGER:{field}:{value}') from e
    return x


def assert_sha(path, expected, label):
    if not isinstance(expected, str) or len(expected) != 64 or digest(require(path)) != expected:
        raise ValueError(f'{label}_SHA_MISMATCH:{path}')


def source_preflight(review_dir: Path, conservation_dir: Path):
    rfile = review_dir / 'all_family_stage10_review.tsv'
    rmeta_file = review_dir / 'stage183_review_provenance.json'
    status_file = conservation_dir / 'all_family_stage11_12_status.tsv'
    summary_file = conservation_dir / 'selected_family_msa_conservation.tsv'
    residues_file = conservation_dir / 'selected_query_residue_conservation.tsv'
    cmeta_file = conservation_dir / 'stage184_conservation_provenance.json'
    rmeta = provenance(rmeta_file)
    cmeta = provenance(cmeta_file)
    if rmeta.get('stage') != 'STAGE10_PARTIAL_EVIDENCE_REVIEW' or rmeta.get('stage10_historical_I2_completed') is not False:
        raise ValueError('INVALID_STAGE183_STATUS')
    if cmeta.get('stage') != 'STAGE11_STAGE12_SEQUENCE_ONLY_PARTIAL' or cmeta.get('structure_based_stage12_completed') is not False or cmeta.get('biochemical_function_confirmed') is not False:
        raise ValueError('INVALID_STAGE184_STATUS')
    assert_sha(rfile, rmeta.get('output_sha256', {}).get(rfile.name), 'STAGE183_OUTPUT')
    for f in (status_file, summary_file, residues_file):
        assert_sha(f, cmeta.get('output_sha256', {}).get(f.name), 'STAGE184_OUTPUT')
    src = cmeta.get('source_sha256')
    if not isinstance(src, dict) or not src:
        raise ValueError('MISSING_STAGE184_SOURCE_PROVENANCE')
    for name, checksum in src.items():
        assert_sha(Path(name), checksum, 'STAGE184_SOURCE')
    # The Stage184 prepare manifest anchors its own source hashes to Stage183;
    # verify its reference chain rather than trusting a merely self-consistent TSV.
    prep_files = [Path(name) for name in src if Path(name).name == 'stage184_prepare_provenance.json']
    if len(prep_files) != 1:
        raise ValueError('STAGE184_PREPARE_MANIFEST_MISSING')
    prep = provenance(prep_files[0])
    if prep.get('stage') != 'STAGE11_SELECTED_MSA_PREPARATION' or prep.get('not_full_stage12') is not True:
        raise ValueError('INVALID_STAGE184_PREPARE_MANIFEST')
    originals = prep.get('source_sha256')
    if not isinstance(originals, dict) or originals.get(str(rfile.resolve())) != digest(rfile) or originals.get(str(rmeta_file.resolve())) != digest(rmeta_file):
        raise ValueError('STAGE183_STAGE184_COHORT_PROVENANCE_MISMATCH')
    for name, checksum in originals.items():
        assert_sha(Path(name), checksum, 'STAGE184_PREPARE_SOURCE')
    for name, checksum in prep.get('output_sha256', {}).items():
        assert_sha(prep_files[0].parent / name, checksum, 'STAGE184_PREPARE_OUTPUT')
    return (rfile, rmeta_file, status_file, summary_file, residues_file, cmeta_file), rmeta, cmeta, prep


def assemble(args):
    if args.limit <= 0:
        raise ValueError('LIMIT_MUST_BE_POSITIVE')
    out = args.out_dir.resolve()
    if out.exists():
        raise ValueError(f'REFUSE_EXISTING_OUTPUT:{out}')
    for source in (args.review_dir, args.conservation_dir):
        src = source.resolve()
        if src == out or src in out.parents or out in src.parents:
            raise ValueError('OUTPUT_OVERLAPS_SOURCE')
    sources, rm, cm, pm = source_preflight(args.review_dir, args.conservation_dir)
    reviews = table(sources[0], STAGE183_COLUMNS)
    all_status = table(sources[2], STAGE184_COLUMNS)
    msa = table(sources[3], SUMMARY_COLUMNS)
    residues = table(sources[4], RESIDUE_COLUMNS)
    rindex = index(reviews, 'structure_id')
    sindex = index(all_status, 'structure_id')
    mindex = index(msa, 'structure_id')
    if len(reviews) != rm.get('all_families') or len(reviews) != cm.get('all_families') or len(reviews) != pm.get('all_families'):
        raise ValueError('PROVENANCE_FAMILY_COUNT_MISMATCH')
    if len(msa) != cm.get('msa_families') or len(msa) != len(pm.get('selected_structure_ids', [])):
        raise ValueError('PROVENANCE_MSA_COUNT_MISMATCH')
    if set(rindex) != set(sindex) or set(mindex) != set(pm['selected_structure_ids']):
        raise ValueError('CROSS_STAGE_FAMILY_SET_MISMATCH')
    if len({r['family_rep'] for r in reviews}) != len(reviews):
        raise ValueError('DUPLICATE_FAMILY_REP')
    count_res = Counter()
    positions = defaultdict(set)
    for row in residues:
        sid = row['structure_id']
        if sid not in mindex:
            raise ValueError(f'UNSELECTED_RESIDUE:{sid}')
        p = nonneg(row['query_position'], 'query_position')
        if p < 1 or p > nonneg(mindex[sid]['query_length'], 'query_length') or p in positions[sid]:
            raise ValueError('RESIDUE_POSITION_INVALID_OR_DUPLICATE')
        positions[sid].add(p)
        n = nonneg(row['non_gap_homologs'], 'non_gap_homologs')
        matches = nonneg(row['matching_homologs'], 'matching_homologs')
        if matches > n or row['minimum_homolog_support_met'] not in {'YES', 'NO'} or row['sequence_conserved'] not in {'YES', 'NO'}:
            raise ValueError('INVALID_RESIDUE_COUNT_OR_STATUS')
        if row['sequence_conserved'] == 'YES' and row['minimum_homolog_support_met'] != 'YES':
            raise ValueError('CONSERVATION_WITHOUT_SUPPORT')
        count_res[(sid, 'analyzable')] += row['minimum_homolog_support_met'] == 'YES'
        count_res[(sid, 'conserved')] += row['sequence_conserved'] == 'YES'
    rows = []
    for i, review in enumerate(reviews, 1):
        sid = review['structure_id']
        stage = sindex[sid]
        if review['review_order'] != str(i) or review['review_tier'] not in TIERS:
            raise ValueError('STAGE183_REVIEW_ORDER_OR_TIER_INVALID')
        if any(review[key] != stage[key] for key in ('structure_id', 'family_rep', 'priority_class')):
            raise ValueError(f'CROSS_STAGE_FAMILY_ID_MISMATCH:{sid}')
        if stage['full_stage12_structure_analysis'] != 'NOT_EVALUATED':
            raise ValueError('STAGE184_FALSE_SPATIAL_STATUS')
        selected = sid in mindex
        if selected:
            sm = mindex[sid]
            if stage['msa_plan_state'] != 'SELECTED_FOR_SEQUENCE_MSA' or stage['msa_conservation_state'] not in MSA_STATES:
                raise ValueError('SELECTED_MSA_STATUS_INCONSISTENT')
            if any(review[k] != sm[k] for k in ('structure_id', 'family_rep', 'family_size')) or sm['msa_conservation_state'] != stage['msa_conservation_state']:
                raise ValueError('MSA_FAMILY_METADATA_MISMATCH')
            if sm['structure_residue_analysis'] != 'NOT_EVALUATED_NO_MATCHED_PDB_PREDICTION_RECORD':
                raise ValueError('FALSE_STRUCTURE_BASED_RESIDUE_EVIDENCE')
            qlen = nonneg(sm['query_length'], 'query_length')
            if len(positions[sid]) != qlen or not positions[sid] == set(range(1, qlen + 1)):
                raise ValueError('MISSING_RESIDUE_POSITIONS')
            usable = nonneg(sm['analyzable_positions'], 'analyzable_positions')
            conserved = nonneg(sm['conserved_positions'], 'conserved_positions')
            if (conserved > usable or usable > qlen or
                usable != count_res[(sid, 'analyzable')] or conserved != count_res[(sid, 'conserved')] or
                str(usable) != stage['analyzable_positions'] or str(conserved) != stage['conserved_positions']):
                raise ValueError('MSA_RESIDUE_SUMMARY_INCONSISTENT')
            msa_unique = str(nonneg(sm['msa_unique_sequences'], 'msa_unique_sequences'))
        else:
            if (stage['msa_plan_state'] != 'NOT_SELECTED_FOR_SEQUENCE_MSA' or stage['msa_conservation_state'] != 'NOT_SELECTED_FOR_MSA' or stage['analyzable_positions'] or stage['conserved_positions']):
                raise ValueError('UNSELECTED_MSA_STATUS_INCONSISTENT')
            msa_unique = ''
        ext = review['stage09_external_status']
        search_state = review['stage09_evidence_state']
        if search_state == 'NOT_EVALUATED_NO_STRUCTURE' and ext:
            raise ValueError('PENDING_STRUCTURE_WITH_HIT_STATUS')
        if search_state not in {'SEARCH_VERIFIED', 'NOT_EVALUATED_NO_STRUCTURE'}:
            raise ValueError('UNSUPPORTED_STAGE09_STATUS')
        if review['biochemical_activity_confirmed'] != 'NOT_TESTED' or review['historical_I2_score']:
            raise ValueError('STAGE183_FALSE_EXPERIMENTAL_OR_I2_EVIDENCE')
        if review['review_tier'] in DARK_TIERS:
            interpretation = ('PRELIMINARY_CONTEXT_AND_VERIFIED_STRUCTURE_HOMOLOGY' if review['review_tier']=='T1_CONTEXT_AND_CONVINCING_STRUCTURE'
                              else 'PRELIMINARY_CONTEXT_SUPPORTED_DARK' if review['review_tier']=='T2_CONTEXT_SUPPORTED_DARK'
                              else 'PRELIMINARY_STRUCTURE_HOMOLOGY_ONLY')
        elif review['review_tier']=='T4_DARK_CONTEXT_UNRESOLVED':
            interpretation = 'UNRESOLVED_DARK_INVENTORY_NOT_NEGATIVE'
        else:
            interpretation = 'KNOWN_REFERENCE_NOT_DARK_TARGET'
        rows.append({
            'review_order': review['review_order'], 'review_tier': review['review_tier'],
            'structure_id': sid, 'family_rep': review['family_rep'],
            'priority_class': review['priority_class'], 'family_size': review['family_size'],
            'sample_span': review['sample_span'], 'genomic_context_state': review['genomic_context_state'],
            'best_noncandidate_anchor_category': review['best_noncandidate_anchor_category'],
            'category_supported_samples': review['category_supported_samples'],
            'category_supported_members': review['category_supported_members'],
            'stage09_evidence_state': search_state, 'stage09_external_status': ext,
            'msa_plan_state': stage['msa_plan_state'],
            'msa_conservation_state': stage['msa_conservation_state'],
            'msa_unique_sequences': msa_unique,
            'msa_query_length': sm['query_length'] if selected else '',
            'msa_analyzable_positions': stage['analyzable_positions'],
            'msa_conserved_positions': stage['conserved_positions'],
            'stage12_spatial_state': 'NOT_EVALUATED',
            'biochemical_activity_confirmed': 'NOT_TESTED',
            'historical_I2_score': '',
            'stage13_evidence_interpretation': interpretation,
            'stage13_selection_state': 'PROVISIONAL_REVIEW_ONLY_NOT_VALIDATED',
        })
    context_dark = [r for r in rows if r['review_tier'] in {'T1_CONTEXT_AND_CONVINCING_STRUCTURE', 'T2_CONTEXT_SUPPORTED_DARK'}]
    selected_msa = [r for r in rows if r['structure_id'] in mindex]
    # Shortlist is an ordered *review queue*, never an efficacy or activity rank.
    review_queue = [r for r in rows if r['review_tier'] in DARK_TIERS][:args.limit]
    print('ALL_FAMILIES='+str(len(rows)))
    print('CONTEXT_SUPPORTED_DARK='+str(len(context_dark)))
    print('MSA_SEQUENCE_EVALUATED='+str(len(selected_msa)))
    print('STRUCTURE_SEARCHED='+str(sum(r['stage09_evidence_state']=='SEARCH_VERIFIED' for r in rows)))
    print('REVIEW_QUEUE='+str(len(review_queue)))
    print('STRUCTURE_BASED_STAGE12=NOT_RUN')
    print('FINAL_BIOCHEMICAL_SELECTION=NOT_RUN')
    if args.mode == 'check':
        print('STEP185_STAGE13_CHECK=PASS');print('FILES_WRITTEN=NO');return
    out.mkdir(parents=True, exist_ok=False)
    outputs = {
        'all_family_stage13_evidence.tsv': rows,
        'context_supported_stage13_review.tsv': context_dark,
        'msa_sequence_stage13_evidence.tsv': selected_msa,
        'provisional_stage13_review_queue.tsv': review_queue,
    }
    for name, data in outputs.items():
        with (out / name).open('x', encoding='utf-8', newline='') as f:
            wr=csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter='\t', lineterminator='\n')
            wr.writeheader();wr.writerows(data)
    manifest = {
        'stage':'STAGE13_PROVISIONAL_EVIDENCE_CATALOG', 'schema_version':1,
        'all_families':len(rows), 'context_supported_dark':len(context_dark),
        'msa_sequence_evaluated':len(selected_msa), 'review_queue_count':len(review_queue),
        'stage12_structure_based_completed':False, 'stage10_historical_i2_completed':False,
        'biochemical_activity_confirmed':False, 'final_candidate_selection_completed':False,
        'source_sha256':{str(f.resolve()):digest(f) for f in sources},
        'output_sha256':{p.name:digest(p) for p in sorted(out.glob('*.tsv'))},
    }
    (out/'stage185_provisional_provenance.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print('STEP185_STAGE13_PROVISIONAL_CATALOG=PASS')
    print('OUTPUT_DIR='+str(out))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--review-dir', type=Path, required=True)
    p.add_argument('--conservation-dir', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    p.add_argument('--limit', type=int, default=30,
                   help='Number of dark family entries in the exploratory review queue')
    p.add_argument('--mode', choices=('check','run'), default='check')
    a=p.parse_args()
    assemble(a)


if __name__=='__main__':
    main()
