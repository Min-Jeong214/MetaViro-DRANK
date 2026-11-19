#!/usr/bin/env python3
"""Validate an external historical ColabFold artifact register.

This does NOT authenticate an original GPU execution. It verifies source
bytes, sequence composition, model metadata and local PDB/score consistency.
Differences from the freshly prepared A3M remain explicitly attributed.
"""

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def nonempty(path):
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'MISSING_OR_EMPTY_FILE:{path}')
    return path


def read_table(path, required):
    with nonempty(path).open(encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(required).issubset(fields):
            raise ValueError(f'INVALID_TABLE_COLUMNS:{path}')
        records = list(reader)
    if not records or any(None in r or any(v is None for v in r.values()) for r in records):
        raise ValueError(f'INVALID_TABLE_ROWS:{path}')
    ids = [r['structure_id'] for r in records]
    if any(not i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError(f'DUPLICATE_OR_EMPTY_IDS:{path}')
    return records


def sequences(path):
    entries = []
    with nonempty(path).open(encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                entry_id = line[1:].split()[0]
                if not entry_id:
                    raise ValueError(f'INVALID_A3M_HEADER:{path}')
                entries.append([entry_id, ''])
            elif entries:
                entries[-1][1] += line
            else:
                raise ValueError(f'A3M_SEQUENCE_BEFORE_HEADER:{path}')
    if not entries or any(not seq for _, seq in entries):
        raise ValueError(f'INVALID_A3M_SEQUENCES:{path}')
    return entries


def normalized(seq):
    return seq.replace('-', '').replace('.', '').upper().rstrip('*')


def finite(value, low, high, name):
    x = float(value)
    if not math.isfinite(x) or x < low or x > high:
        raise ValueError(f'INVALID_{name}:{value}')
    return x


AA3 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E',
    'GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K',
    'MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
    'TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O',
}


def validate_pdb_score(pdb, score, expected_length, sid, query):
    data = json.loads(score.read_text(encoding='utf-8'))
    values = data.get('plddt')
    if not isinstance(values, list) or len(values) != expected_length:
        raise ValueError(f'PREDICTION_PLDDT_LENGTH_MISMATCH:{sid}')
    values = [finite(v, 0, 100, 'PLDDT') for v in values]
    ptm = finite(data.get('ptm'), 0, 1, 'PTM')
    seen = set()
    chains = set()
    for line in pdb.read_text(errors='replace').splitlines():
        if not line.startswith('ATOM  ') or line[12:16].strip() != 'CA' or line[16] not in (' ', 'A'):
            continue
        pos = int(line[22:26])
        if pos in seen or pos < 1 or pos > expected_length or line[26] != ' ':
            raise ValueError(f'PDB_RESIDUE_MAPPING_INVALID:{sid}')
        seen.add(pos)
        chains.add(line[21])
        amino = AA3.get(line[17:20].strip())
        if amino is None or amino != query[pos - 1]:
            raise ValueError(f'PDB_QUERY_RESIDUE_MISMATCH:{sid}:{pos}')
        if abs(finite(line[60:66], 0, 100, 'PDB_PLDDT') - values[pos - 1]) > 0.051:
            raise ValueError(f'PDB_SCORE_RESIDUE_MISMATCH:{sid}:{pos}')
    if len(chains) != 1 or len(seen) != expected_length:
        # The existing Stage 11 QC requires a complete CA trace.
        raise ValueError(f'HISTORICAL_PDB_INCOMPLETE_FOR_QC:{sid}:{len(seen)}')
    return round(sum(values)/len(values), 2), round(ptm, 4)


def prepare(args):
    prep = args.prepared_dir.resolve()
    up = args.stage08_output.resolve()
    out = args.out_dir.resolve()
    if out.exists():
        raise ValueError(f'EXISTING_OUTPUT_DIRECTORY:{out}')
    if any(out == d or out in d.parents or d in out.parents for d in (prep, up)):
        raise ValueError('INPUT_OUTPUT_OVERLAP')
    reg = nonempty(args.reuse_register)
    if out == reg or out in reg.parents or reg in out.parents:
        raise ValueError('REGISTER_OUTPUT_OVERLAP')
    register = read_table(reg, (
        'structure_id','reuse_class','new_a3m_sha256','historical_a3m',
        'historical_a3m_sha256','historical_pdb','historical_pdb_sha256',
        'historical_score_json','historical_score_sha256'
    ))
    mapping = read_table(nonempty(prep/'selected_structure_map.tsv'), (
        'structure_id', 'family_rep','length_aa','sequence_sha256'
    ))
    policy = read_table(nonempty(prep/'prediction_eligibility.tsv'), (
        'structure_id','msa_prediction_status'
    ))
    eligible = {r['structure_id'] for r in policy if r['msa_prediction_status']=='MSA_PREDICTION_ELIGIBLE'}
    mmap = {r['structure_id']:r for r in mapping}
    if set(mmap) != {r['structure_id'] for r in policy} or {r['structure_id'] for r in register} != set(mmap):
        raise ValueError('REGISTER_SELECTED_FAMILY_SET_MISMATCH')
    statuses = {}
    records = []
    for r in register:
        sid = r['structure_id']
        new_a3m = nonempty(prep/'msa/a3m_colabfold'/f'{sid}.a3m')
        if sha(new_a3m) != r['new_a3m_sha256']:
            raise ValueError(f'PREPARED_A3M_SHA_MISMATCH:{sid}')
        status = r['reuse_class']
        if status == 'PREDICTION_PENDING':
            if sid not in eligible or any(r.get(k) for k in (
                'historical_a3m','historical_pdb','historical_score_json'
            )):
                raise ValueError(f'INVALID_PENDING_ROW:{sid}')
            statuses[sid] = 'PREDICTION_PENDING'
            continue
        if status not in ('HISTORICAL_MSA_REUSE_CANDIDATE', 'EXACT_MSA_REUSE_CANDIDATE') or sid not in eligible:
            raise ValueError(f'REUSE_CLASS_NOT_ALLOWED:{sid}:{status}')
        a3m = nonempty(r['historical_a3m'])
        pdb = nonempty(r['historical_pdb'])
        score = nonempty(r['historical_score_json'])
        for source, checksum, label in (
            (a3m, r['historical_a3m_sha256'], 'A3M'),
            (pdb, r['historical_pdb_sha256'], 'PDB'),
            (score, r['historical_score_sha256'], 'SCORE'),
        ):
            if sha(source) != checksum:
                raise ValueError(f'HISTORICAL_{label}_SHA_MISMATCH:{sid}')
        fresh = sequences(new_a3m)
        old = sequences(a3m)
        seqsha = mmap[sid]['sequence_sha256']
        length = int(mmap[sid]['length_aa'])
        for label, seqs in (('NEW', fresh),('OLD',old)):
            first_id, first_seq = seqs[0]
            if first_id != sid or len(normalized(first_seq)) != length or hashlib.sha256(normalized(first_seq).encode()).hexdigest() != seqsha:
                raise ValueError(f'{label}_QUERY_MISMATCH:{sid}')
        if len(old) < 2 or Counter(map(normalized, (s for _,s in old))) != Counter(map(normalized,(s for _,s in fresh))):
            raise ValueError(f'HISTORICAL_SEQUENCE_MULTISET_MISMATCH:{sid}')
        exact = [s for _,s in old] == [s for _,s in fresh]
        alignment_relation = 'EXACT_ALIGNED_RECORDS' if exact else 'UNGAPPED_MULTISET_ONLY'
        if status == 'EXACT_MSA_REUSE_CANDIDATE' and not exact:
            raise ValueError(f'FALSE_EXACT_MSA_REUSE_CLASS:{sid}')
        mean, ptm = validate_pdb_score(pdb, score, length, sid, normalized(old[0][1]))
        record = {
            'structure_id': sid,
            'status': 'HISTORICAL_IMPORTED',
            'provenance_type': 'HISTORICAL_MSA_REUSE',
            'historical_alignment_relation': alignment_relation,
            'source_a3m': str(a3m), 'a3m_sha256': sha(a3m),
            'prepared_a3m': str(new_a3m), 'prepared_a3m_sha256': sha(new_a3m),
            'pdb_file': str(pdb), 'pdb_sha256': sha(pdb),
            'score_json': str(score), 'score_json_sha256': sha(score),
            'model_type': 'alphafold2_ptm',
            'query_length': length, 'query_sha256': seqsha,
            'sequence_count': len(old),
            'mean_plddt': mean, 'ptm': ptm,
            'gpu_prediction_performed': False,
            'original_run_provenance_verified': False,
            'historical_attribution': 'SOURCE_FILES_SUPPLIED_BY_USER_NOT_ORIGINAL_EXECUTION_LOG',
            'reuse_register_sha256': sha(reg),
        }
        records.append(record)
        statuses[sid] = alignment_relation
    return out, register, records, statuses


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--reuse-register', type=Path, required=True)
    parser.add_argument('--prepared-dir', type=Path, required=True)
    parser.add_argument('--stage08-output', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--mode', choices=['check','run'], default='check')
    args=parser.parse_args()
    out, register, records, statuses = prepare(args)
    if args.mode == 'run':
        out.mkdir(parents=True, exist_ok=False)
        index_path = out/'prediction_index.tsv'
        with index_path.open('x', encoding='utf-8', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=['structure_id','prediction_provenance'],delimiter='\t',lineterminator='\n')
            writer.writeheader()
            for record in records:
                path = out/f"{record['structure_id']}.prediction_provenance.json"
                path.write_text(json.dumps(record, indent=2, ensure_ascii=False)+'\n',encoding='utf-8')
                writer.writerow({'structure_id':record['structure_id'], 'prediction_provenance':str(path)})
        (out/'import_provenance.json').write_text(json.dumps({
            'status': 'HISTORICAL_ARTIFACT_IMPORT_NOT_GPU_RUN',
            'imported_families': len(records),
            'pending_families': sum(v=='PREDICTION_PENDING' for v in statuses.values()),
            'reuse_register_sha256':sha(args.reuse_register),
            'prediction_index_sha256': sha(index_path),
            'original_execution_provenance_verified':False,
        },indent=2)+'\n',encoding='utf-8')
        print('PREDICTION_INDEX='+str(index_path))
    print('HISTORICAL_IMPORT_CANDIDATES_QC_PASS='+str(len(records)))
    print('PENDING_PREDICTIONS='+str(sum(v=='PREDICTION_PENDING' for v in statuses.values())))
    print('HISTORICAL_ALIGNMENT_MODIFIED='+str(sum(v=='UNGAPPED_MULTISET_ONLY' for v in statuses.values())))
    print('ORIGINAL_GPU_EXECUTION_PROVENANCE_VERIFIED=NO')
    print('IMPORT_'+args.mode.upper()+'=PASS')


if __name__=='__main__':
    try:
        main()
    except (OSError,ValueError,KeyError,TypeError,AssertionError,csv.Error) as exc:
        sys.exit('ERROR: '+str(exc))
