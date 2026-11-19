#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

MISSING = {"", "-", "NA", "N/A", "nan", "None"}

def blank(v):
    return (v or "").strip() in MISSING

def fasta_lengths(path):
    d, name, seq = {}, None, []
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    d[name] = len(''.join(seq).rstrip('*'))
                name=line[1:].split()[0]; seq=[]
            else:
                seq.append(line.strip())
        if name is not None:
            d[name] = len(''.join(seq).rstrip('*'))
    return d

ap=argparse.ArgumentParser()
ap.add_argument('--annotations', required=True, help='eggNOG-mapper .emapper.annotations')
ap.add_argument('--faa', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--min-aa', type=int, default=150)
args=ap.parse_args()

lens=fasta_lengths(args.faa)
rows={}
with open(args.annotations) as fh:
    header=None
    for line in fh:
        if line.startswith('##'): continue
        if line.startswith('#query'):
            header=line.lstrip('#').rstrip('\n').split('\t'); continue
        if line.startswith('#') or not line.strip(): continue
        if header is None: continue
        vals=line.rstrip('\n').split('\t')
        rec=dict(zip(header, vals))
        rows[rec['query']]=rec

fields=['Description','Preferred_name','KEGG_ko','PFAMs','PFAM']
with open(args.out,'w',newline='') as out:
    w=csv.writer(out, delimiter='\t')
    w.writerow(['protein_id','length_aa','eggnog_status','strict_unresolved','ge_min_aa'])
    for pid,L in lens.items():
        rec=rows.get(pid)
        if rec is None:
            status='NO_HIT'; unresolved=True
        else:
            status='HIT'
            unresolved=all(blank(rec.get(f,'')) for f in fields)
        w.writerow([pid,L,status,int(unresolved),int(unresolved and L>=args.min_aa)])
