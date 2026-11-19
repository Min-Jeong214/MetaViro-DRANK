#!/usr/bin/env python3
"""Add verified structure-supported residue evidence to an existing Stage185 catalog.

Does NOT calculate spatial residue clusters, catalytic activity or final I2 ranking.
Stage185 must be rebuilt from the same Stage183 review that includes real
external structure evidence. This script preserves all Family rows.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

STAGE12_SUMMARY = (
    'structure_id', 'length', 'msa_nseq', 'very_conserved_positions',
    'priority_structural_positions', 'priority_acidic', 'priority_basic',
    'priority_nucleophilic',
)
WINDOW_COLS = ('structure_id', 'start', 'end', 'sequence', 'conserved_structural',
               'charged', 'acidic', 'basic', 'mean_query_fraction', 'mean_plddt')
CATALOG_COLS = (
    'structure_id', 'review_tier', 'stage09_evidence_state',
    'stage09_external_status', 'msa_conservation_state', 'stage12_spatial_state',
    'biochemical_activity_confirmed', 'historical_I2_score',
)
NEW_FIELDS = (
    'stage12_structure_residue_state', 'stage12_very_conserved_positions',
    'stage12_priority_structural_positions', 'stage12_priority_acidic',
    'stage12_priority_basic', 'stage12_priority_nucleophilic',
    'stage12_conserved_window_count', 'stage12_pdb_sha256',
    'stage12_a3m_sha256', 'stage12_residue_table_sha256',
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def required(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'MISSING_OR_EMPTY:{path}')
    return path


def check_sha(path: Path, digest: str) -> None:
    if not isinstance(digest, str) or len(digest) != 64 or sha(required(path)) != digest:
        raise ValueError(f'SHA256_MISMATCH:{path}')


def load_json(path: Path) -> dict:
    data = json.loads(required(path).read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise ValueError(f'INVALID_JSON:{path}')
    return data


def load_tsv(path: Path, columns: tuple, *, allow_empty: bool = False):
    with required(path).open('r', encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(columns).issubset(fields):
            raise ValueError(f'INVALID_TSV_COLUMNS:{path}')
        rows = list(reader)
    if (not rows and not allow_empty) or any(None in r or None in r.values() for r in rows):
        raise ValueError(f'INVALID_TSV_ROWS:{path}')
    return fields, rows


def unique(rows, key):
    result = {}
    for r in rows:
        val = r[key]
        if not val or val in result:
            raise ValueError(f'DUPLICATE_OR_EMPTY:{key}:{val}')
        result[val] = r
    return result


def nonneg(value, label):
    if not str(value).isdigit():
        raise ValueError(f'INVALID_INTEGER:{label}:{value}')
    return int(value)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--catalog-dir', type=Path, required=True,
                    help='Stage185 rebuilt from structure-updated Stage10')
    ap.add_argument('--stage12-dir', type=Path, required=True,
                    help='Verified Stage12 structure-residue output directory')
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--mode', choices=('check', 'run'), default='check')
    args = ap.parse_args()
    cdir, sdir, out = (p.expanduser().resolve() for p in
                       (args.catalog_dir, args.stage12_dir, args.out_dir))
    if out.exists():
        raise ValueError(f'REFUSE_EXISTING_OUTPUT:{out}')
    for p in (cdir, sdir):
        if out == p or out in p.parents or p in out.parents:
            raise ValueError('OUTPUT_OVERLAPS_INPUT')

    catalog_file = cdir / 'all_family_stage13_evidence.tsv'
    catalog_meta_file = cdir / 'stage185_provisional_provenance.json'
    stage12_meta_file = sdir / 'conservation_provenance.json'
    summary_file = sdir / 'candidate_conservation_summary.tsv'
    windows_file = sdir / 'conserved_windows.tsv'

    cmeta = load_json(catalog_meta_file)
    smeta = load_json(stage12_meta_file)
    if cmeta.get('stage') != 'STAGE13_PROVISIONAL_EVIDENCE_CATALOG' or \
       cmeta.get('stage12_structure_based_completed') is not False or \
       cmeta.get('biochemical_activity_confirmed') is not False:
        raise ValueError('UNEXPECTED_STAGE185_MANIFEST')
    if smeta.get('algorithm') != 'ORIGINAL_STAGE12_FUNCTIONS' or \
       smeta.get('final_top5_automatically_selected') is not False:
        raise ValueError('UNEXPECTED_STAGE12_MANIFEST')

    for filename, digest in cmeta.get('output_sha256', {}).items():
        check_sha(cdir / filename, digest)
    check_sha(catalog_file, cmeta.get('output_sha256', {}).get(catalog_file.name))
    for filename, digest in cmeta.get('source_sha256', {}).items():
        check_sha(Path(filename), digest)
    check_sha(summary_file, smeta.get('summary_sha256'))
    check_sha(windows_file, smeta.get('windows_sha256'))

    # Compare Family-level MSA evidence by content SHA, not by path. Both stages
    # can be run in separate output directories using the same aligned sequences.
    stage184_sources = {}
    for path in cmeta.get('source_sha256', {}):
        if Path(path).name == 'stage184_conservation_provenance.json':
            meta = load_json(Path(path))
            if meta.get('stage') == 'STAGE11_STAGE12_SEQUENCE_ONLY_PARTIAL':
                stage184_sources = meta.get('source_sha256', {})
    if not stage184_sources:
        raise ValueError('STAGE184_PROVENANCE_NOT_FOUND')

    base_fields, base_rows = load_tsv(catalog_file, CATALOG_COLS)
    _, summary_rows = load_tsv(summary_file, STAGE12_SUMMARY)
    _, windows = load_tsv(windows_file, WINDOW_COLS, allow_empty=True)
    cat = unique(base_rows, 'structure_id')
    summ = unique(summary_rows, 'structure_id')
    if len(base_rows) != cmeta.get('all_families'):
        raise ValueError('STAGE185_FAMILY_COUNT_MISMATCH')
    if len(summary_rows) != smeta.get('family_count'):
        raise ValueError('STAGE12_FAMILY_COUNT_MISMATCH')
    families = smeta.get('families')
    if not isinstance(families, list) or len(families) != len(summ):
        raise ValueError('INVALID_STAGE12_FAMILY_PROVENANCE')

    window_count = Counter()
    for w in windows:
        sid = w['structure_id']
        if sid not in summ:
            raise ValueError(f'WINDOW_FOR_UNEVALUATED_FAMILY:{sid}')
        st, en = nonneg(w['start'], 'start'), nonneg(w['end'], 'end')
        if st < 1 or en < st or en > nonneg(summ[sid]['length'], 'length'):
            raise ValueError(f'INVALID_CONSERVATION_WINDOW:{sid}')
        window_count[sid] += 1

    seen = set()
    for family in families:
        sid = family.get('structure_id')
        if sid not in summ or sid in seen or sid not in cat:
            raise ValueError(f'FAMILY_ID_MISMATCH:{sid}')
        seen.add(sid)
        row = cat[sid]
        if row['msa_conservation_state'] != 'SEQUENCE_ONLY_CONSERVATION_AVAILABLE':
            raise ValueError(f'STAGE185_MSA_NOT_EVALUATED:{sid}')
        if row['stage12_spatial_state'] != 'NOT_EVALUATED':
            raise ValueError(f'FALSE_SPATIAL_STATUS:{sid}')
        if row['biochemical_activity_confirmed'] != 'NOT_TESTED' or row['historical_I2_score']:
            raise ValueError('FALSE_BIOCHEMICAL_OR_FINAL_SCORE')
        a3m = Path(family['a3m']).resolve()
        aligned = Path(family['alignment']).resolve()
        pfile = Path(family['prediction_provenance']).resolve()
        pmeta = load_json(pfile)
        pdb = Path(pmeta['pdb_file']).resolve()
        seq = pmeta.get('a3m_sha256')
        if pmeta.get('status') != 'COMPLETED' or pmeta.get('structure_id') != sid:
            raise ValueError(f'INVALID_PREDICTION_PROVENANCE:{sid}')
        for file, digest in ((a3m, family['a3m_sha256']),
                             (aligned, family['alignment_sha256']),
                             (pfile, family['prediction_provenance_sha256']),
                             (pdb, family['pdb_sha256']),
                             (sdir / f'{sid}_residue_conservation.tsv', family['residue_table_sha256'])):
            check_sha(file, digest)
        if seq != family['a3m_sha256'] or pmeta.get('pdb_sha256') != family['pdb_sha256']:
            raise ValueError(f'PREDICTION_INPUT_MISMATCH:{sid}')
        if not any(Path(path).name == f'{sid}.a3m' and digest == seq
                   for path, digest in stage184_sources.items()):
            raise ValueError(f'CROSS_STAGE_A3M_SHA_MISMATCH:{sid}')
        residue_fields, resrows = load_tsv(sdir / f'{sid}_residue_conservation.tsv', ('structure_id',))
        if len(resrows) != nonneg(summ[sid]['length'], 'length') or \
           len(resrows) != family.get('residue_count') or \
           nonneg(summ[sid]['msa_nseq'], 'msa_nseq') != family.get('sequence_count'):
            raise ValueError(f'STAGE12_RESIDUE_SUMMARY_MISMATCH:{sid}')
        if any(r['structure_id'] != sid for r in resrows):
            raise ValueError(f'RESIDUE_ID_MISMATCH:{sid}')
        if nonneg(summ[sid]['priority_structural_positions'], 'priority_structural_positions') > len(resrows):
            raise ValueError(f'IMPOSSIBLE_PRIORITY_COUNT:{sid}')

    outputs = []
    for row in base_rows:
        sid = row['structure_id']
        enriched = dict(row)
        enriched.update({k: '' for k in NEW_FIELDS})
        enriched['stage12_structure_residue_state'] = 'NOT_EVALUATED'
        if sid in summ:
            sr = summ[sid]
            pr = next(f for f in families if f['structure_id'] == sid)
            enriched.update({
                'stage12_structure_residue_state': 'STRUCTURE_SUPPORTED_RESIDUES_AVAILABLE',
                'stage12_very_conserved_positions': sr['very_conserved_positions'],
                'stage12_priority_structural_positions': sr['priority_structural_positions'],
                'stage12_priority_acidic': sr['priority_acidic'],
                'stage12_priority_basic': sr['priority_basic'],
                'stage12_priority_nucleophilic': sr['priority_nucleophilic'],
                'stage12_conserved_window_count': str(window_count[sid]),
                'stage12_pdb_sha256': pr['pdb_sha256'],
                'stage12_a3m_sha256': pr['a3m_sha256'],
                'stage12_residue_table_sha256': pr['residue_table_sha256'],
            })
        outputs.append(enriched)

    if len(outputs) != len(base_rows) or sum(r['stage12_structure_residue_state'] ==
                                            'STRUCTURE_SUPPORTED_RESIDUES_AVAILABLE' for r in outputs) != len(summ):
        raise ValueError('INTEGRATION_FAMILY_COVERAGE_MISMATCH')

    print(f'ALL_FAMILIES={len(outputs)}')
    print(f'STRUCTURE_RESIDUE_EVALUATED={len(summ)}')
    print(f'CONSERVATION_WINDOWS={len(windows)}')
    print('SPATIAL_3D_CLUSTER_ANALYSIS=NOT_RUN')
    print('FINAL_BIOCHEMICAL_SELECTION=NOT_RUN')
    if args.mode == 'check':
        print('STAGE186_STRUCTURE_INTEGRATION_CHECK=PASS')
        print('FILES_WRITTEN=NO')
        return

    out.mkdir(parents=True, exist_ok=False)
    fields = list(base_fields) + list(NEW_FIELDS)
    paths = {'all_family_stage13_structure_integrated.tsv': outputs,
             'structure_evaluated_stage13_review.tsv': [r for r in outputs
                                                       if r['stage12_structure_residue_state'] ==
                                                       'STRUCTURE_SUPPORTED_RESIDUES_AVAILABLE']}
    for filename, rows in paths.items():
        with (out / filename).open('x', encoding='utf-8', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    manifest = {
        'stage': 'STAGE13_STRUCTURE_RESIDUE_EVIDENCE_INTEGRATION',
        'schema_version': 1,
        'all_families': len(outputs),
        'structure_residue_evaluated': len(summ),
        'conservation_windows': len(windows),
        'spatial_3d_cluster_analysis_completed': False,
        'final_biochemical_selection_completed': False,
        'source_sha256': {str(p.resolve()): sha(p) for p in
                          (catalog_file, catalog_meta_file, stage12_meta_file,
                           summary_file, windows_file)},
        'output_sha256': {p.name: sha(p) for p in sorted(out.glob('*.tsv'))},
    }
    (out / 'stage186_structure_integration_provenance.json').write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print('STAGE186_STRUCTURE_INTEGRATION=PASS')
    print(f'OUTPUT_DIR={out}')


if __name__ == '__main__':
    main()
