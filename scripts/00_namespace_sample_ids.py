#!/usr/bin/env python3
"""Copy sample FASTA/GFF/eggNOG inputs into a unique sample-ID namespace.

Read-only with respect to input datasets. This is an explicit preprocessing
operation, NOT a biological inference or a correction of invalid GFF mappings.
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

NEEDED = ('proteins.faa', 'genes.gff', 'annotations.emapper.annotations')
NAME_RX = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
ID_RX = re.compile(r'(?:^|;)ID=([^;]+)')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for buf in iter(lambda: fh.read(1 << 20), b''):
            h.update(buf)
    return h.hexdigest()


def probe(sample_dir):
    paths = [sample_dir / name for name in NEEDED]
    for p in paths:
        if not p.is_file() or p.stat().st_size == 0:
            raise ValueError(f'Missing/empty source input: {p}')

    faa, gff, annotations = paths
    ids = set()
    current = None
    seq_length = 0
    with faa.open(encoding='utf-8-sig') as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith('>'):
                if current is not None and seq_length == 0:
                    raise ValueError(f'Empty protein {sample_dir.name}:{current}')
                match = re.match(r'^>(\S+)', line)
                if not match:
                    raise ValueError(f'Empty FASTA header at {faa}:{lineno}')
                current = match.group(1)
                if current in ids:
                    raise ValueError(f'Duplicate protein WITHIN sample {sample_dir.name}: {current}')
                ids.add(current)
                seq_length = 0
            elif line.strip():
                if current is None:
                    raise ValueError(f'FASTA sequence before header at {faa}:{lineno}')
                if not re.fullmatch(r'[A-Za-z*]+', line.strip()):
                    raise ValueError(f'Invalid protein sequence at {faa}:{lineno}')
                seq_length += len(line.strip())
    if not ids or seq_length == 0:
        raise ValueError(f'Empty last protein in {faa}')

    gff_ids = set()
    contigs = set()
    with gff.open(encoding='utf-8-sig') as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith('##FASTA'):
                raise ValueError(f'Embedded GFF FASTA unsupported: {gff}')
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\r\n').split('\t')
            if len(cols) != 9:
                raise ValueError(f'GFF column count !=9: {gff}:{lineno}')
            contigs.add(cols[0])
            if cols[2] != 'CDS':
                continue
            match = ID_RX.search(cols[8])
            if not match:
                raise ValueError(f'GFF CDS lacks ID: {gff}:{lineno}')
            try:
                i = int(match.group(1).rsplit('_', 1)[1])
                a, b = int(cols[3]), int(cols[4])
                if a < 1 or b < a or i < 0:
                    raise ValueError('Invalid CDS index/range')
            except (IndexError, ValueError) as exc:
                raise ValueError(f'GFF invalid Prodigal ID/coordinates: {gff}:{lineno}') from exc
            pid = f'{cols[0]}_{i}'
            if pid in gff_ids:
                raise ValueError(f'Duplicate GFF-derived protein ID: {pid}')
            gff_ids.add(pid)
    if gff_ids != ids:
        miss = sorted(ids - gff_ids)[:3]
        extra = sorted(gff_ids - ids)[:3]
        raise ValueError(f'FASTA/GFF mismatch for {sample_dir.name}; FAA-only={miss}; GFF-only={extra}')

    header = None
    annotated = set()
    with annotations.open(encoding='utf-8-sig') as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith('#query'):
                if header is not None:
                    raise ValueError('Duplicate eggNOG header')
                header = line[1:].rstrip('\r\n').split('\t')
                for required in ('query', 'eggNOG_OGs', 'Description', 'Preferred_name'):
                    if required not in header:
                        raise ValueError(f'Missing eggNOG {required}')
                continue
            if not line.strip() or line.startswith('#'):
                continue
            if header is None:
                raise ValueError(f'eggNOG data before header: {annotations}:{lineno}')
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != len(header):
                raise ValueError(f'eggNOG column count mismatch: {annotations}:{lineno}')
            pid = fields[header.index('query')]
            if pid not in ids:
                raise ValueError(f'eggNOG protein ID not found: {sample_dir.name}:{pid}')
            if pid in annotated:
                raise ValueError(f'Duplicate eggNOG query: {pid}')
            annotated.add(pid)
    if header is None:
        raise ValueError(f'Missing eggNOG #query header in {annotations}')
    return {'ids': ids, 'protein_count': len(ids), 'gff_cds_count': len(gff_ids),
            'annotation_count': len(annotated), 'contig_count': len(contigs),
            'files': {name: {'sha256': digest(path), 'bytes': path.stat().st_size}
                      for name, path in zip(NEEDED, paths)}}


def transform(faa, gff, ann, sample, dest):
    """Transform all three cross-reference surfaces; source remains untouched."""
    prefix = sample + '|'
    dest.mkdir(parents=True, exist_ok=False)
    with faa.open(encoding='utf-8-sig') as inp, (dest / NEEDED[0]).open('w', encoding='utf-8', newline='') as out:
        for line in inp:
            if line.startswith('>'):
                m = re.match(r'^>(\S+)(.*)$', line.rstrip('\r\n'))
                assert m is not None
                out.write('>' + prefix + m.group(1) + m.group(2) + '\n')
            else:
                out.write(line)
    with gff.open(encoding='utf-8-sig') as inp, (dest / NEEDED[1]).open('w', encoding='utf-8', newline='') as out:
        for line in inp:
            if line.startswith('##sequence-region '):
                parts = line.rstrip('\r\n').split(' ', 2)
                assert len(parts) == 3
                out.write(parts[0] + ' ' + prefix + parts[1] + ' ' + parts[2] + '\n')
            elif not line.strip() or line.startswith('#'):
                out.write(line)
            else:
                cols = line.rstrip('\r\n').split('\t')
                cols[0] = prefix + cols[0]
                out.write('\t'.join(cols) + '\n')
    with ann.open(encoding='utf-8-sig') as inp, (dest / NEEDED[2]).open('w', encoding='utf-8', newline='') as out:
        fields = None
        for line in inp:
            if line.startswith('#query'):
                fields = line[1:].rstrip('\r\n').split('\t')
                out.write(line)
            elif not line.strip() or line.startswith('#'):
                out.write(line)
            else:
                cols = line.rstrip('\r\n').split('\t')
                cols[fields.index('query')] = prefix + cols[fields.index('query')]
                out.write('\t'.join(cols) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir', required=True, type=Path)
    ap.add_argument('--output-dir', required=True, type=Path)
    ap.add_argument('--mode', choices=('check', 'run'), default='check')
    args = ap.parse_args()
    source, target = args.input_dir.resolve(), args.output_dir.resolve()
    if not source.is_dir():
        raise ValueError(f'Input directory missing: {source}')
    if target == source or target in source.parents or source in target.parents:
        raise ValueError('Input and output directory cannot overlap')
    if target.exists():
        raise ValueError(f'Destination already exists (will not overwrite): {target}')
    samples = sorted(p for p in source.iterdir() if p.is_dir() and not p.name.startswith('.'))
    if not samples:
        raise ValueError('No sample directories')
    report = {}
    id_owner = {}
    collisions = set()
    for sample in samples:
        if not NAME_RX.fullmatch(sample.name):
            raise ValueError(f'Invalid sample name: {sample.name}')
        info = probe(sample)
        report[sample.name] = info
        for pid in info['ids']:
            if pid in id_owner:
                collisions.add(pid)
            else:
                id_owner[pid] = sample.name
        print(f'INPUT_SAMPLE={sample.name} PROTEINS={info["protein_count"]} GFF_CDS={info["gff_cds_count"]} EGGNOG_ANNOTATED={info["annotation_count"]}', flush=True)
    print(f'CROSS_SAMPLE_DUPLICATE_IDS={len(collisions)}', flush=True)
    print('EXAMPLE_DUPLICATE_IDS=' + ','.join(sorted(collisions)[:5]), flush=True)
    print('SOURCE_FASTA_GFF_EGGNOG_MAPPING=PASS', flush=True)
    if args.mode == 'check':
        print('NAMESPACE_CHECK=PASS NO_FILES_WRITTEN', flush=True)
        return
    # Only after all source inputs pass strict integrity checks, write copies.
    target.mkdir(parents=True, exist_ok=False)
    for sample in samples:
        transform(*(sample / name for name in NEEDED), sample.name, target / sample.name)
        result = probe(target / sample.name)
        assert result['protein_count'] == report[sample.name]['protein_count']
        assert result['gff_cds_count'] == report[sample.name]['gff_cds_count']
        assert result['annotation_count'] == report[sample.name]['annotation_count']
        assert result['ids'] == {sample.name + '|' + pid for pid in report[sample.name]['ids']}
    manifest = {
        'purpose': 'Sample-namespace normalization for portable MetaViro-DRANK',
        'source_directory': str(source), 'output_directory': str(target),
        'method': 'Prefix each protein ID, each GFF seqid and each eggNOG query with sample|; retain sequence bodies, coordinates, annotations and GFF CDS attributes',
        'cross_sample_duplicate_count': len(collisions),
        'original_inputs': {
            name: {k: v for k, v in info.items() if k != 'ids'}
            for name, info in report.items()
        },
        'normalized_inputs': {
            sample.name: probe(target / sample.name)['files'] for sample in samples
        }
    }
    (target / 'namespace_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(f'NORMALIZED_INPUT_DIRECTORY={target}', flush=True)
    print(f'NORMALIZED_PROTEINS={sum(x["protein_count"] for x in report.values())}', flush=True)
    print('NAMESPACE_RUN=PASS', flush=True)
    print('ORIGINAL_FILES_MODIFIED=NO', flush=True)

if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, AssertionError) as exc:
        sys.exit('NAMESPACE_VALIDATION=FAIL: ' + str(exc))
