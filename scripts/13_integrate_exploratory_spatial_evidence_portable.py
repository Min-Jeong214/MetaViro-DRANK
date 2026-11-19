#!/usr/bin/env python3
"""Safely merge exploratory Stage12 C-alpha cluster results into a Stage186 catalog.

Zero detected clusters means EVALUATED and NO_CLUSTER_DETECTED, not NOT_EVALUATED.
This is exploratory spatial evidence, not a validated catalytic-site or enzyme claim.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

CATALOG_REQUIRED = (
    'structure_id', 'stage12_structure_residue_state', 'stage12_spatial_state',
    'stage12_pdb_sha256', 'stage12_residue_table_sha256',
    'biochemical_activity_confirmed', 'historical_I2_score',
)
SUMMARY_REQUIRED = ('structure_id', 'length', 'eligible_charged_residues', 'cluster_count', 'status')
CLUSTER_REQUIRED = (
    'structure_id', 'cluster_rank', 'residues', 'residue_count',
    'acidic_count', 'basic_count', 'mean_plddt', 'max_edge_ca_distance_angstrom',
)
ADDED = (
    'stage12_spatial_eligible_charged_residues', 'stage12_spatial_cluster_count',
    'stage12_spatial_ca_cutoff_angstrom', 'stage12_spatial_min_cluster_residues',
    'stage12_spatial_analysis_mode',
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def good_file(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f'MISSING_OR_EMPTY:{path}')
    return path


def verify(path: Path, digest: str) -> None:
    if not isinstance(digest, str) or len(digest) != 64 or sha(good_file(path)) != digest:
        raise ValueError(f'SHA256_MISMATCH:{path}')


def jsonfile(path: Path) -> dict:
    obj = json.loads(good_file(path).read_text(encoding='utf-8'))
    if not isinstance(obj, dict):
        raise ValueError(f'INVALID_JSON:{path}')
    return obj


def tsv(path: Path, fields, *, allow_empty=False):
    with good_file(path).open('r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not set(fields).issubset(names):
            raise ValueError(f'INVALID_TSV_SCHEMA:{path}')
        records = list(reader)
    if (not records and not allow_empty) or any(
        None in record or any(value is None for value in record.values())
        for record in records
    ):
        raise ValueError(f'INVALID_TSV_ROWS:{path}')
    return names, records


def unique(rows, label='structure_id'):
    result = {}
    for row in rows:
        key = row[label]
        if not key or key in result:
            raise ValueError(f'DUPLICATE_OR_EMPTY_{label}:{key}')
        result[key] = row
    return result


def integer(s, label, *, minimum=0):
    if not isinstance(s, str) or not s.isdigit() or int(s) < minimum:
        raise ValueError(f'INVALID_INTEGER:{label}:{s}')
    return int(s)


def positive_number(s, label):
    try:
        n = float(s)
    except (TypeError, ValueError):
        raise ValueError(f'INVALID_FLOAT:{label}') from None
    if not math.isfinite(n) or n <= 0:
        raise ValueError(f'INVALID_FLOAT:{label}')
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage186-dir', required=True, type=Path)
    ap.add_argument('--spatial-dir', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--mode', choices=('check', 'run'), default='check')
    args = ap.parse_args()
    base, spatial, out = (p.expanduser().resolve() for p in
                          (args.stage186_dir, args.spatial_dir, args.out_dir))
    if out.exists():
        raise ValueError(f'OUTPUT_ALREADY_EXISTS:{out}')
    if any(out == p or out in p.parents or p in out.parents for p in (base, spatial)):
        raise ValueError('OUTPUT_OVERLAPS_INPUT')

    base_meta_file = base / 'stage186_structure_integration_provenance.json'
    base_catalog = base / 'all_family_stage13_structure_integrated.tsv'
    spatial_meta_file = spatial / 'spatial_cluster_provenance.json'
    spatial_summary = spatial / 'family_spatial_cluster_summary.tsv'
    spatial_details = spatial / 'charged_spatial_clusters.tsv'

    bm = jsonfile(base_meta_file)
    sm = jsonfile(spatial_meta_file)
    if bm.get('stage') != 'STAGE13_STRUCTURE_RESIDUE_EVIDENCE_INTEGRATION' or \
       bm.get('spatial_3d_cluster_analysis_completed') is not False or \
       bm.get('final_biochemical_selection_completed') is not False:
        raise ValueError('INVALID_STAGE186_PROVENANCE')
    if sm.get('analysis_mode') != 'NEW_EXPLORATORY_L3' or \
       sm.get('historical_stage_l3_reproduction') is not False or \
       sm.get('biochemical_activity_inferred') is not False or \
       sm.get('geometry') != 'CA_DISTANCE_CONNECTED_COMPONENTS' or \
       sm.get('residue_filter') != 'STAGE12_PRIORITY_RESIDUE_AND_DEKRH':
        raise ValueError('INVALID_SPATIAL_PROVENANCE')

    for filename, digest in bm.get('output_sha256', {}).items():
        verify(base / filename, digest)
    verify(base_catalog, bm.get('output_sha256', {}).get(base_catalog.name))
    for filename, digest in (
        (spatial_summary, sm.get('summary_table_sha256')),
        (spatial_details, sm.get('cluster_table_sha256')),
    ):
        verify(filename, digest)

    stage12_src = [
        (Path(filename), digest)
        for filename, digest in bm.get('source_sha256', {}).items()
        if Path(filename).name == 'conservation_provenance.json'
    ]
    if len(stage12_src) != 1:
        raise ValueError('STAGE12_LINEAGE_MISSING_OR_AMBIGUOUS')
    stage12_file, stage12_digest = stage12_src[0]
    verify(stage12_file, stage12_digest)
    if sm.get('stage12_provenance_sha256') != stage12_digest:
        raise ValueError('STAGE12_LINEAGE_SHA256_MISMATCH')

    cutoff = positive_number(sm.get('ca_distance_cutoff_angstrom'), 'ca-cutoff')
    minimum = sm.get('min_cluster_residues')
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 2:
        raise ValueError('INVALID_MIN_CLUSTER_RESIDUES')
    fields, all_rows = tsv(base_catalog, CATALOG_REQUIRED)
    _, summaries = tsv(spatial_summary, SUMMARY_REQUIRED)
    _, clusters = tsv(spatial_details, CLUSTER_REQUIRED, allow_empty=True)
    by_id = unique(all_rows)
    spatial_by_id = unique(summaries)
    if len(all_rows) != bm.get('all_families') or \
       len(summaries) != sm.get('family_count') or \
       len(clusters) != sm.get('cluster_count'):
        raise ValueError('PROVENANCE_ROW_COUNTS_MISMATCH')
    if not set(spatial_by_id).issubset(by_id):
        raise ValueError('SPATIAL_FAMILY_NOT_IN_CATALOG')

    inputs = sm.get('inputs')
    if not isinstance(inputs, list) or len(inputs) != len(summaries):
        raise ValueError('INVALID_SPATIAL_INPUTS')
    input_by_id = unique(inputs)
    if set(input_by_id) != set(spatial_by_id):
        raise ValueError('SPATIAL_INPUT_FAMILY_MISMATCH')

    per_family = defaultdict(list)
    ranks = defaultdict(set)
    for c in clusters:
        sid = c['structure_id']
        if sid not in spatial_by_id:
            raise ValueError('CLUSTER_FAMILY_NOT_IN_SUMMARY')
        rank = integer(c['cluster_rank'], 'cluster_rank', minimum=1)
        if rank in ranks[sid]:
            raise ValueError('DUPLICATE_CLUSTER_RANK')
        ranks[sid].add(rank)
        n = integer(c['residue_count'], 'residue_count', minimum=minimum)
        acid = integer(c['acidic_count'], 'acidic_count')
        basic = integer(c['basic_count'], 'basic_count')
        if acid + basic != n or len(c['residues'].split(',')) != n:
            raise ValueError('CLUSTER_RESIDUE_COUNT_MISMATCH')
        plddt = float(c['mean_plddt'])
        edge = positive_number(c['max_edge_ca_distance_angstrom'], 'maxedge')
        if not math.isfinite(plddt) or not 0 <= plddt <= 100 or edge > cutoff + 0.001:
            raise ValueError('INVALID_CLUSTER_METRICS')
        per_family[sid].append(c)

    statuses = Counter()
    for sid, s in spatial_by_id.items():
        cat = by_id[sid]
        checks = input_by_id[sid]
        if cat['stage12_structure_residue_state'] != 'STRUCTURE_SUPPORTED_RESIDUES_AVAILABLE' or \
           cat['stage12_spatial_state'] != 'NOT_EVALUATED':
            raise ValueError(f'SPATIAL_REQUIRES_EVALUATED_STAGE12:{sid}')
        if cat['stage12_residue_table_sha256'] != checks.get('residue_table_sha256') or \
           cat['stage12_pdb_sha256'] != checks.get('msa_pdb_sha256'):
            raise ValueError(f'SPATIAL_INPUT_DIGEST_MISMATCH:{sid}')
        n = integer(s['cluster_count'], 'cluster_count')
        eligible = integer(s['eligible_charged_residues'], 'eligible_charged_residues')
        integer(s['length'], 'length', minimum=1)
        if n != len(per_family[sid]) or \
           sum(integer(c['residue_count'], 'residue_count') for c in per_family[sid]) > eligible or \
           (n > 0 and set(ranks[sid]) != set(range(1, n + 1))):
            raise ValueError(f'CLUSTER_SUMMARY_MISMATCH:{sid}')
        expected = 'CLUSTERS_DETECTED' if n else 'NO_CLUSTER_DETECTED'
        if s['status'] != expected:
            raise ValueError(f'SPATIAL_STATUS_INCONSISTENT:{sid}')
        statuses[expected] += 1

    enriched = []
    for row in all_rows:
        value = dict(row)
        value.update({key: '' for key in ADDED})
        sid = row['structure_id']
        if sid in spatial_by_id:
            s = spatial_by_id[sid]
            n = integer(s['cluster_count'], 'cluster_count')
            value.update({
                'stage12_spatial_state': ('EXPLORATORY_CLUSTERS_DETECTED' if n
                                          else 'EXPLORATORY_NO_CLUSTER_DETECTED'),
                'stage12_spatial_eligible_charged_residues': s['eligible_charged_residues'],
                'stage12_spatial_cluster_count': s['cluster_count'],
                'stage12_spatial_ca_cutoff_angstrom': str(cutoff),
                'stage12_spatial_min_cluster_residues': str(minimum),
                'stage12_spatial_analysis_mode': 'NEW_EXPLORATORY_L3',
            })
        if value['biochemical_activity_confirmed'] != 'NOT_TESTED' or value['historical_I2_score']:
            raise ValueError('UNEXPECTED_BIOCHEMICAL_OR_FINAL_SCORE')
        enriched.append(value)

    if len(enriched) != len(all_rows):
        raise ValueError('FAMILY_COVERAGE_MISMATCH')
    print(f'ALL_FAMILIES={len(all_rows)}')
    print(f'SPATIAL_EVALUATED={len(spatial_by_id)}')
    print(f'NO_CLUSTER_DETECTED_FAMILIES={statuses["NO_CLUSTER_DETECTED"]}')
    print(f'CLUSTERS_DETECTED_FAMILIES={statuses["CLUSTERS_DETECTED"]}')
    print(f'TOTAL_SPATIAL_CLUSTERS={len(clusters)}')
    print(f'SPATIAL_NOT_EVALUATED={len(all_rows)-len(spatial_by_id)}')
    print('FINAL_BIOCHEMICAL_SELECTION=NOT_RUN')
    if args.mode == 'check':
        print('STAGE187_SPATIAL_INTEGRATION_CHECK=PASS')
        print('FILES_WRITTEN=NO')
        return

    out.mkdir(parents=True, exist_ok=False)
    final_fields = fields + list(ADDED)
    produced = {
        'all_family_stage13_spatial_integrated.tsv': enriched,
        'spatial_evaluated_stage13_review.tsv': [
            r for r in enriched if r['structure_id'] in spatial_by_id
        ],
    }
    for filename, rows in produced.items():
        with (out / filename).open('x', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=final_fields, delimiter='\t',
                                    lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    with (out / 'spatial_cluster_details.tsv').open('x', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CLUSTER_REQUIRED, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        writer.writerows(clusters)
    manifest = {
        'stage': 'STAGE13_EXPLORATORY_SPATIAL_EVIDENCE_INTEGRATION',
        'schema_version': 1,
        'analysis_mode': 'NEW_EXPLORATORY_L3',
        'all_families': len(all_rows),
        'spatial_evaluated_families': len(spatial_by_id),
        'spatial_no_cluster_families': statuses['NO_CLUSTER_DETECTED'],
        'spatial_cluster_positive_families': statuses['CLUSTERS_DETECTED'],
        'spatial_not_evaluated_families': len(all_rows) - len(spatial_by_id),
        'total_spatial_clusters': len(clusters),
        'final_biochemical_selection_completed': False,
        'source_sha256': {str(p): sha(p) for p in
                          (base_meta_file, base_catalog, spatial_meta_file,
                           spatial_summary, spatial_details, stage12_file)},
        'output_sha256': {p.name: sha(p) for p in sorted(out.glob('*.tsv'))},
    }
    (out / 'stage187_spatial_integration_provenance.json').write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print('STAGE187_SPATIAL_INTEGRATION=PASS')
    print(f'OUTPUT_DIR={out}')


if __name__ == '__main__':
    main()
