#!/usr/bin/env python3
"""Rebuild H3 family-wide ORF context from Stage D/E/F and raw GFF/eggNOG.

Prototype: writes ONLY to a separate output directory. It never uses H3 as input.
"""
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

SAMPLES = ("V695", "V696", "V699", "V710", "V733")
STRUCTURAL = re.compile(
    r"\b(?:virion|tail|capsid|portal|terminase|baseplate|sheath|fiber|tape measure|virion structural)\b", re.I
)
ANCHOR_REGEX = (
    ("polymerase", r"\bpolymerase\b"),
    ("primase", r"\bprimase\b"),
    ("helicase", r"\bhelicase\b"),
    ("nuclease", r"nuclease|\brnase\b"),
    ("ligase", r"\bligase\b"),
    ("recombination", r"recombinase|recombination|resolvase|topoisomerase"),
    ("dna_methylation", r"\bdna[ -]+methyl(?:ase|transferase)\b"),
    ("ssdna_binding", r"\bsingle[- ]stranded(?:[- ]dna)?[- ]binding\b|\bsingle[- ]stranded\s+DNA\s+binding\b|\bssb\b|\bgp32\b"),
    ("restriction", r"\brestriction\b.{0,35}\b(?:enzyme|subunit)\b|\bres\s+subunit\b"),
    ("dna_repair", r"\bdna[- ]repair\b|\bphotolyase\b|\bgp59\b|\buvsy\b"),
    ("reverse_transcriptase", r"\breverse[- ]transcriptase\b"),
)
ANCHORS = [(name, re.compile(pat, re.I)) for name, pat in ANCHOR_REGEX]
MEMBER_FIELDS = [
    'family_rep','member','sample','length_aa','candidate_original','eggnog_status',
    'strict_unresolved_eggnog','gff_seqid','gff_local_index','gff_order',
    'neighbors_available','anchor_within5','anchor_min_gene_distance',
    'anchor_categories','anchor_ids','structural_neighbors_within5',
    'structural_dense_ge2','reconstructed_candidate','noncandidate_exclusion_reason',
]
FAMILY_FIELDS = [
    'structure_id','family_rep','family_size','sample_span','samples','candidate_count',
    'noncandidate_count','all_anchor_context_count','all_anchor_context_fraction',
    'candidate_anchor_context_count','candidate_anchor_context_fraction',
    'noncandidate_anchor_context_count','noncandidate_anchor_context_fraction',
    'all_context_sample_span','noncandidate_context_sample_span',
    'dominant_anchor_category','dominant_anchor_supported_members',
    'dominant_anchor_consistency','all_anchor_categories',
    'structural_context_any_count','structural_context_any_fraction',
    'structural_dense_ge2_count','structural_dense_ge2_fraction','family_context_tier',
]


def tsv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def fasta(path):
    found = {}
    name, seq = None, []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    if name in found: raise ValueError(f'duplicate FASTA ID: {name}')
                    found[name] = ''.join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        if name in found: raise ValueError(f'duplicate FASTA ID: {name}')
        found[name] = ''.join(seq)
    return found


def sample_of(member):
    m = re.match(r'^(V\d+)_Virus_Ch\|', member)
    if not m: raise ValueError(f'Invalid member ID: {member}')
    return m.group(1)


def ontology(desc, pref):
    text = ' '.join(x for x in (desc, pref) if x and x != '-')
    cats = {name for name, rx in ANCHORS if rx.search(text)}
    # The H3 ontology uses UvsY as a multifunctional anchor.
    if re.search(r'\buvsy\b', text, re.I):
        cats.update(('dna_repair', 'recombination', 'ssdna_binding'))
    if re.search(r'\batpase\b', text, re.I) and re.search(
        r'\bdna\b|\breplication\b|\brepair\b|\brecombination\b|\bhelicase\b|\bprimase\b|\bpolymerase\b', text, re.I
    ):
        cats.add('dna_associated_atpase')  # Audit will verify the exact H3 label.
    return ';'.join(sorted(cats)), bool(STRUCTURAL.search(text))


def read_eggnog(path):
    hits = {}
    names = None
    with path.open() as f:
        for line in f:
            if line.startswith('#query'):
                names = line[1:].rstrip('\n').split('\t')
                if 'Description' not in names or 'Preferred_name' not in names:
                    raise ValueError(f'eggNOG columns unavailable in {path}: {names}')
                continue
            if not line.strip() or line.startswith('#'):
                continue
            if names is None: raise ValueError(f'No eggNOG header in {path}')
            vals = line.rstrip('\n').split('\t')
            if len(vals) != len(names):
                raise ValueError(f'Bad eggNOG row {path}: {len(vals)} != {len(names)}')
            row = dict(zip(names, vals))
            query = row['query']
            hits[query] = ontology(row['Description'], row['Preferred_name'])
    return hits


def read_gff(path, sample, eggnog, official_anchors):
    contigs = defaultdict(list)
    seen = set()
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'CDS': continue
            attrs = dict(z.split('=',1) for z in p[8].split(';') if '=' in z)
            gene_idx = attrs.get('ID','').rsplit('_',1)[-1]
            if not gene_idx.isdigit():
                raise ValueError(f'Cannot derive gene index: {path} {p[:3]} {attrs.get("ID")}')
            seqid = p[0]
            local_id = f'{seqid}_{gene_idx}'
            full_id = f'{sample}_Virus_Ch|{local_id}'
            if full_id in seen: raise ValueError(f'Duplicate GFF ID: {full_id}')
            seen.add(full_id)
            _, structural = eggnog.get(local_id, ('', False))
            label = official_anchors.get(full_id, '')
            contigs[seqid].append((int(p[3]), int(p[4]), len(contigs[seqid]), full_id, gene_idx, label, structural))
    lookup = {}
    for seqid, genes in contigs.items():
        genes.sort(key=lambda r:(r[0],r[1],r[2]))
        for order, gene in enumerate(genes, 1):
            lookup[gene[3]] = (seqid, order, genes)
    return lookup


def member_context(entry, sequence, status, gl, window):
    mid = entry['member']
    rep = entry['representative']
    sample = sample_of(mid)
    length = len(sequence.rstrip('*'))
    seqid, order, genes = gl
    i = order - 1
    left, right = max(0,i-window), min(len(genes),i+window+1)
    neighbors = [genes[j] for j in range(left,right) if j != i]
    anchor_hits = [g for g in neighbors if g[5]]
    cats = ','.join(sorted({g[5] for g in anchor_hits}))
    aids = ','.join(sorted(g[3] for g in anchor_hits))
    min_dist = min((abs(order-(genes.index(g)+1)) for g in anchor_hits),default=None)
    struct_count = sum(1 for g in neighbors if g[6])
    cand = int(length >= 150 and bool(anchor_hits) and status in ('EGGNOG_NO_HIT','EGGNOG_UNRESOLVED_HIT'))
    reasons = []
    if not cand:
        if length < 150: reasons.append('LENGTH_LT150')
        if not anchor_hits: reasons.append('NO_ANCHOR_WITHIN5')
    return dict(zip(MEMBER_FIELDS, [
        rep,mid,sample,str(length),str(int(entry['is_candidate'])),
        {'EGGNOG_NO_HIT':'NO_HIT','EGGNOG_UNRESOLVED_HIT':'HIT_STRICT_UNRESOLVED'}[status],
        '1',seqid,genes[i][4],str(order),str(len(neighbors)),str(int(bool(anchor_hits))),
        '' if min_dist is None else str(min_dist),cats,aids,str(struct_count),
        str(int(struct_count>=2)),str(cand),','.join(reasons),
    ]))


def frac(n, total):
    return f'{n/total:.4f}' if total else '0.0000'


def family_row(rep, group, structure):
    size = len(group)
    cand = [r for r in group if r['candidate_original']=='1']
    noncand = [r for r in group if r['candidate_original']=='0']
    anchors = [r for r in group if r['anchor_within5']=='1']
    c_anchor = [r for r in cand if r['anchor_within5']=='1']
    nc_anchor = [r for r in noncand if r['anchor_within5']=='1']
    counts = Counter()
    for r in anchors:
        # Count each distinct category-combination at most once per member.
        for cat in dict.fromkeys(x for x in r['anchor_categories'].split(',') if x):
            counts[cat] += 1
    dom, support = counts.most_common(1)[0] if counts else ('',0)
    nonspan = len({r['sample'] for r in nc_anchor})
    tier = ('NO_NONCANDIDATE_CONTEXT_SUPPORT' if nonspan==0 else
            'NONCANDIDATE_CONTEXT_1_SAMPLE' if nonspan==1 else
            'NONCANDIDATE_CONTEXT_2_SAMPLES' if nonspan==2 else
            'NONCANDIDATE_CONTEXT_GE3_SAMPLES')
    structural = [r for r in group if int(r['structural_neighbors_within5']) > 0]
    dense = [r for r in group if r['structural_dense_ge2']=='1']
    samples = sorted({r['sample'] for r in group})
    return dict(zip(FAMILY_FIELDS, [
        structure,rep,str(size),str(len(samples)),','.join(samples),str(len(cand)),str(len(noncand)),
        str(len(anchors)),frac(len(anchors),size),str(len(c_anchor)),frac(len(c_anchor),len(cand)),
        str(len(nc_anchor)),frac(len(nc_anchor),len(noncand)),
        str(len({r['sample'] for r in anchors})),str(nonspan),dom,str(support),frac(support,len(anchors)),
        ','.join(sorted(counts)),str(len(structural)),frac(len(structural),size),
        str(len(dense)),frac(len(dense),size),tier,
    ]))


def write_table(path, columns, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(data)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--out-dir',type=Path,required=True)
    ap.add_argument(
        '--stage-h-order',
        action='store_true',
        help='Reproduce historical Stage H row ordering'
    )
    ap.add_argument('--window',type=int,default=5)
    a=ap.parse_args()
    if a.window!=5: raise ValueError('Stage H3 canonical requires window=5')
    base=a.root/'Result/Annotation_saturation/Pilot_v2'
    d=base/'StageD_Recurrent191_FamilyAudit'
    e=base/'StageE_Recurrent191_AllMemberSaturation'
    sf=base/'StageF_Structure/structure_id_map.tsv'
    entries=tsv(d/'recurrent191_member_map.tsv')
    statuses={r['member']:r['eggnog_status'] for r in tsv(e/'recurrent191_all_members_original_status.tsv')}
    seqs=fasta(d/'recurrent191_all_members.faa')
    struct_map={r['representative']:r['structure_id'] for r in tsv(sf)}

    anchor_path = base/'pilot_v2_viral_anchors.tsv'
    anchor_rows = tsv(anchor_path)

    if len(anchor_rows) != 11338:
        raise ValueError(
            f'Unexpected official anchor count: {len(anchor_rows)}'
        )

    official_anchors = {}

    for row in anchor_rows:
        anchor_id = row['anchor']
        categories = row['categories']

        if (
            not anchor_id
            or not categories
            or anchor_id in official_anchors
        ):
            raise ValueError(
                f'Invalid or duplicate official anchor: {anchor_id!r}'
            )

        official_anchors[anchor_id] = categories

    print(
        'OFFICIAL_ANCHORS_LOADED=',
        len(official_anchors),
        flush=True
    )

    if len(entries)!=10327 or len(seqs)!=10327 or len(statuses)!=10327 or len(struct_map)!=191:
        raise ValueError('Frozen pilot counts unexpected')
    print('READING RAW ANNOTATION + GFF',flush=True)
    gff={}
    for sample in SAMPLES:
        egg=read_eggnog(a.root/f'Analysis/{sample}_Virus_Ch.emapper.annotations')
        lookup=read_gff(a.root/f'Analysis/{sample}_Virus_Ch.gff',sample,egg,official_anchors)
        gff.update(lookup)
        print(sample,'eggnog',len(egg),'gff CDS',len(lookup),flush=True)
    groups=defaultdict(list)
    members=[]
    for row in entries:
        mid=row['member']
        if mid not in gff or mid not in statuses or mid not in seqs:
            raise ValueError(f'Unmatched member {mid}: GFF={mid in gff}, status={mid in statuses}, FAA={mid in seqs}')
        rec=member_context(row,seqs[mid],statuses[mid],gff[mid],a.window)
        members.append(rec)
        groups[row['representative']].append(rec)

    # Stage E order determines first-seen category tie resolution.
    e_order = [
        r['member']
        for r in tsv(
            e/'recurrent191_all_members_original_status.tsv'
        )
    ]

    member_set = {r['member'] for r in members}

    if (
        len(e_order) != len(members)
        or len(set(e_order)) != len(e_order)
        or set(e_order) != member_set
    ):
        raise ValueError("Stage E member ordering mismatch")

    e_rank = {
        member: i
        for i, member in enumerate(e_order)
    }

    for group in groups.values():
        group.sort(
            key=lambda r: e_rank[r['member']]
        )

    families = [
        family_row(rep, group, struct_map[rep])
        for rep, group in groups.items()
    ]

    # Stage H supplies only historical serialization order.
    # No H3 canonical values are used in reconstruction.
    if a.stage_h_order:

        h = base/'StageH_FamilyWideContext'

        m_order = [
            r['member']
            for r in tsv(
                h/'recurrent191_familywide_member_context.tsv'
            )
        ]

        f_order = [
            r['family_rep']
            for r in tsv(
                h/'recurrent191_familywide_context_summary.tsv'
            )
        ]

        member_by_id = {
            r['member']: r for r in members
        }

        family_by_rep = {
            r['family_rep']: r for r in families
        }

        if (
            len(m_order) != len(members)
            or len(set(m_order)) != len(m_order)
            or set(m_order) != set(member_by_id)
        ):
            raise ValueError(
                "Historical member order mismatch"
            )

        if (
            len(f_order) != len(families)
            or len(set(f_order)) != len(f_order)
            or set(f_order) != set(family_by_rep)
        ):
            raise ValueError(
                "Historical family order mismatch"
            )

        members = [
            member_by_id[mid] for mid in m_order
        ]

        families = [
            family_by_rep[rep] for rep in f_order
        ]

        print(
            'STAGE_H_HISTORICAL_ORDER_APPLIED=YES',
            flush=True
        )

    else:

        members.sort(
            key=lambda r: (
                struct_map[r['family_rep']],
                r['member']
            )
        )

        families.sort(
            key=lambda r: r['structure_id']
        )

        print(
            'STAGE_H_HISTORICAL_ORDER_APPLIED=NO',
            flush=True
        )

    if len(members)!=10327 or len(families)!=191:
        raise ValueError('Output counts unexpected')
    member_out=a.out_dir/'recurrent191_member_context_corrected.tsv'
    family_out=a.out_dir/'recurrent191_family_context_corrected.tsv'
    write_table(member_out,MEMBER_FIELDS,members)
    write_table(family_out,FAMILY_FIELDS,families)
    print('WROTE',member_out, len(members))
    print('WROTE',family_out, len(families))
    print('RECONSTRUCTION_GENERATED=YES',flush=True)

if __name__ == '__main__': main()
