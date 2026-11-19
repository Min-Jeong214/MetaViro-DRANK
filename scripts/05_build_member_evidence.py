#!/usr/bin/env python3
"""Join Stage05 evidence for all selected family members."""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

UNKNOWN = re.compile(
    r'\b(uncharacteri[sz]ed|hypothetical|unknown function|DUF\d+)\b',
    re.I
)
DUF = re.compile(r'^(DUF\d+|UPF\d+)$', re.I)
STRUCTURAL = {'head and packaging', 'tail', 'connector'}
NUCLEIC = 'dna, rna and nucleotide metabolism'


def fasta_ids(path):
    ids, aliases = set(), {}

    with path.open() as fh:
        for line in fh:
            if not line.startswith('>'):
                continue

            token = line[1:].split()[0]

            if token in ids:
                raise ValueError(f'Duplicate FASTA ID: {token}')

            ids.add(token)
            candidates = {token}

            match = re.fullmatch(
                r'(?:sp|tr)\|([^|]+)\|[^|]+', token
            )
            if match:
                candidates.add(match.group(1))

            for alias in candidates:
                old = aliases.get(alias)
                if old is not None and old != token:
                    raise ValueError(
                        f'Ambiguous alias: {alias}: {old}, {token}'
                    )
                aliases[alias] = token

    if not ids:
        raise ValueError(f'Empty FASTA: {path}')

    return ids, aliases


def tsv_rows(path):
    with path.open() as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip('\r\n')
            if line and not line.startswith('#'):
                yield n, line.split('\t')



def sha256_file(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1024 * 1024), b""
        ):
            h.update(block)

    return h.hexdigest()


def verify_search_manifest(args):
    path = args.search_out / "search_run_complete.json"

    if not path.is_file():
        raise ValueError(
            "Missing verified Stage05 completion manifest: "
            + str(path)
        )

    manifest = json.loads(path.read_text())

    if (
        manifest.get("schema_version") != 1
        or manifest.get("status") != "COMPLETE"
    ):
        raise ValueError(
            "Invalid Stage05 completion status or schema"
        )

    query = manifest.get("query", {})

    if query.get("path") != str(
        args.proteins_faa.resolve()
    ):
        raise ValueError(
            "Input FASTA path differs from search manifest"
        )

    if query.get("sha256") != sha256_file(
        args.proteins_faa
    ):
        raise ValueError(
            "Input FASTA SHA256 mismatch"
        )

    sequence_count = sum(
        line.startswith(b">")
        for line in args.proteins_faa.open("rb")
    )

    if (
        sequence_count == 0
        or query.get("sequence_count") != sequence_count
    ):
        raise ValueError(
            "Input FASTA protein count mismatch"
        )

    for name in (
        "swissprot.tsv",
        "phrog.native.tsv",
        "pfam.domtblout",
    ):
        output = args.search_out / name
        expected = manifest.get(
            "outputs", {}
        ).get(name, {})

        if not output.is_file():
            raise ValueError(
                f"Missing search output: {name}"
            )

        if (
            expected.get("size") != output.stat().st_size
            or expected.get("sha256")
            != sha256_file(output)
        ):
            raise ValueError(
                f"Search result fingerprint mismatch: {name}"
            )

    for name in ("swissprot", "phrog", "pfam"):
        info = manifest.get(
            "databases", {}
        ).get(name, {})

        if not info.get("path"):
            raise ValueError(
                f"Missing database provenance: {name}"
            )

        db = Path(info["path"])

        if not db.is_file():
            raise ValueError(
                f"Database missing: {name}"
            )

        stat = db.stat()

        if (
            str(db.resolve()) != info["path"]
            or stat.st_size != info.get("size")
            or stat.st_mtime_ns != info.get("mtime_ns")
        ):
            raise ValueError(
                f"Database fingerprint mismatch: {name}"
            )

    print("SEARCH_MANIFEST_VALIDATION=PASS")
    print("QUERY_SHA256=PASS")
    print("SEARCH_OUTPUT_SHA256=PASS")
    print("DATABASE_PROVENANCE=PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--members', type=Path, required=True)
    ap.add_argument('--proteins-faa', type=Path, required=True)
    ap.add_argument('--search-out', type=Path, required=True)
    ap.add_argument('--search-log', type=Path, required=True)
    ap.add_argument(
        '--phrog-annotations', type=Path, required=True
    )
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()

    for path in (
        args.members,
        args.proteins_faa,
        args.search_log,
        args.phrog_annotations
    ):
        if not path.is_file():
            raise ValueError(f'Missing input: {path}')

    log = args.search_log.read_text(errors='replace')

    if 'ANNOTATION SEARCH COMPLETE' not in log:
        raise ValueError(
            'Stage05 completion marker missing from search log'
        )

    for name in (
        'swissprot.tsv',
        'phrog.native.tsv',
        'pfam.domtblout'
    ):
        if not (args.search_out / name).is_file():
            raise ValueError(f'Missing search output: {name}')

    verify_search_manifest(args)

    seqs, aliases = fasta_ids(args.proteins_faa)

    members = {}

    with args.members.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')

        if not {'member', 'representative'} <= set(
            reader.fieldnames or []
        ):
            raise ValueError('Invalid members TSV columns')

        for row in reader:
            member = row['member'].strip()
            rep = row['representative'].strip()

            if not member or not rep or member in members:
                raise ValueError(
                    f'Duplicate or empty member: {row}'
                )

            members[member] = rep

    if set(members) != seqs:
        raise ValueError(
            'FASTA/member mismatch: '
            f'missing_in_fasta={sorted(set(members)-seqs)[:5]}, '
            f'extra_in_fasta={sorted(seqs-set(members))[:5]}'
        )

    def resolve(raw, source):
        if raw not in aliases:
            raise ValueError(
                f'{source}: unmatched protein ID: {raw}'
            )
        return aliases[raw]

    annotations = {}

    with args.phrog_annotations.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')

        required = {'phrog', 'annot', 'category'}

        if not required <= set(reader.fieldnames or []):
            raise ValueError('Invalid PHROG annotation columns')

        for row in reader:
            phrog_id = (
                'phrog_' + row['phrog'].removeprefix('phrog_')
            )

            if phrog_id in annotations:
                raise ValueError(
                    f'Duplicate PHROG annotation: {phrog_id}'
                )

            annotations[phrog_id] = row

    swiss = defaultdict(list)
    phrog = defaultdict(list)
    pfam = defaultdict(list)

    # Swiss-Prot: 9 columns
    for n, cols in tsv_rows(
        args.search_out / 'swissprot.tsv'
    ):
        if len(cols) != 9:
            raise ValueError(
                f'Swiss-Prot line {n}: expected 9 columns'
            )

        member = resolve(cols[0], 'Swiss-Prot')
        desc = (
            cols[2].split(' OS=', 1)[0]
            .partition(' ')[2].strip()
        )

        swiss[member].append((
            float(cols[7]),
            -float(cols[8]),
            cols[1],
            desc,
            not bool(UNKNOWN.search(desc))
        ))

    # PHROG: first column is PHROG profile ID
    #        second column is target protein ID
    for n, cols in tsv_rows(
        args.search_out / 'phrog.native.tsv'
    ):
        if len(cols) < 5:
            raise ValueError(
                f'PHROG line {n}: expected >=5 columns'
            )

        phrog_id, raw_member = cols[:2]
        member = resolve(raw_member, 'PHROG')

        if phrog_id not in annotations:
            raise ValueError(
                f'Missing PHROG annotation: {phrog_id}'
            )

        ann = annotations[phrog_id]
        cat = (ann['category'] or '').strip()
        desc = (ann['annot'] or '').strip()

        known = (
            cat.lower() not in {
                'unknown function', 'unknown', ''
            }
            and (
                cat.lower() != 'other'
                or (
                    desc
                    and desc.lower() not in {
                        'na', 'null', 'unknown function'
                    }
                )
            )
        )

        phrog[member].append((
            float(cols[4]),
            phrog_id,
            cat,
            desc,
            known
        ))

    # Pfam domtblout: target HMM, query protein
    with (args.search_out / 'pfam.domtblout').open() as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith('#'):
                continue

            cols = line.split()

            if len(cols) < 13:
                raise ValueError(
                    f'Pfam line {n}: expected >=13 columns'
                )

            member = resolve(cols[3], 'Pfam')

            pfam[member].append((
                float(cols[12]),
                cols[0],
                cols[1],
                not bool(DUF.fullmatch(cols[0]))
            ))

    fields = [
        'member',
        'representative',
        'search_status',
        'sequence_status',
        'swissprot_hits',
        'swissprot_known_hits',
        'swissprot_best_target',
        'swissprot_best_description',
        'phrog_hits',
        'phrog_known_hits',
        'phrog_unknown_hits',
        'phrog_ids',
        'phrog_categories',
        'phrog_structural',
        'phrog_nucleic',
        'pfam_domains',
        'pfam_known_hits',
        'pfam_duf_hits'
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    counts = defaultdict(int)

    with args.out.open('w', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter='\t',
            lineterminator='\n'
        )
        writer.writeheader()

        for member in sorted(members):
            s = sorted(swiss[member])
            p = sorted(phrog[member])
            f = sorted(pfam[member])

            s_known = sum(x[4] for x in s)
            p_known = sum(x[4] for x in p)
            f_known = sum(x[3] for x in f)

            f_duf = len(f) - f_known
            p_unknown = len(p) - p_known

            if s_known or p_known or f_known:
                status = 'KNOWN_SEQUENCE_EVIDENCE'
            elif f_duf:
                status = 'DUF_ONLY'
            elif p_unknown:
                status = 'PHROG_UNKNOWN'
            elif s or p or f:
                status = 'UNRESOLVED_HOMOLOGY'
            else:
                status = 'STRICT_DARK'

            counts[status] += 1

            writer.writerow({
                'member': member,
                'representative': members[member],
                'search_status': 'COMPLETE',
                'sequence_status': status,
                'swissprot_hits': len(s),
                'swissprot_known_hits': s_known,
                'swissprot_best_target':
                    s[0][2] if s else '',
                'swissprot_best_description':
                    s[0][3] if s else '',
                'phrog_hits': len(p),
                'phrog_known_hits': p_known,
                'phrog_unknown_hits': p_unknown,
                'phrog_ids': ';'.join(
                    dict.fromkeys(x[1] for x in p)
                ),
                'phrog_categories': ';'.join(
                    sorted(set(x[2] for x in p))
                ),
                'phrog_structural': int(
                    any(x[2].lower() in STRUCTURAL for x in p)
                ),
                'phrog_nucleic': int(
                    any(x[2].lower() == NUCLEIC for x in p)
                ),
                'pfam_domains': ';'.join(
                    dict.fromkeys(x[1] for x in f)
                ),
                'pfam_known_hits': f_known,
                'pfam_duf_hits': f_duf
            })

    print('MEMBER_EVIDENCE=PASS')
    print(f'TOTAL_MEMBERS={len(members)}')

    for status in (
        'KNOWN_SEQUENCE_EVIDENCE',
        'DUF_ONLY',
        'PHROG_UNKNOWN',
        'UNRESOLVED_HOMOLOGY',
        'STRICT_DARK'
    ):
        print(f'{status}={counts[status]}')

    if counts['UNRESOLVED_HOMOLOGY']:
        print(
            'NOTE=Stage06 requires '
            'UNRESOLVED_HOMOLOGY support'
        )


if __name__ == '__main__':
    main()
