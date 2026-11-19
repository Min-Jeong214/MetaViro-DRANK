#!/usr/bin/env python3
"""Structure-informed Stage13 catalog join smoke tests (synthetic cohort)."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/13_integrate_structure_residue_evidence_portable.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path, rows):
    with path.open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


class Stage186Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mv186_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.catalog = self.root / 'stage185'
        self.stage12 = self.root / 'stage12'
        self.out = self.root / 'stage186'
        self.review = self.root / 'review'
        self.review.mkdir()
        self.prepare = self.root / 'prepare'
        self.msa = self.root / 'msa'
        self.msa.mkdir()
        (self.msa / 'a3m_colabfold').mkdir()
        fasta = self.root / 'all.faa'
        fasta.write_text('>s1|a\nACDE\n>s1|b\nACDF\n>s1|c\nACDG\n')
        member = self.root / 'members.tsv'
        write_tsv(member, [dict(member=a, representative='s1|a') for a in ('s1|a','s1|b','s1|c')])
        mapfile = self.root / 'map.tsv'
        write_tsv(mapfile,[dict(structure_id='X001',family_rep='s1|a',family_size='3',sequence_sha256=hashlib.sha256(b'ACDE').hexdigest())])
        rfile = self.review / 'all_family_stage10_review.tsv'
        write_tsv(rfile,[dict(review_order='1',review_tier='T1_CONTEXT_AND_CONVINCING_STRUCTURE',structure_id='X001',family_rep='s1|a',family_size='3',priority_class='A_DEEP_DARK',sample_span='2',genomic_context_state='DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT',best_noncandidate_anchor_category='nuclease',category_supported_samples='2',category_supported_members='3',stage09_evidence_state='SEARCH_VERIFIED',stage09_external_status='CONVINCING_STRUCTURE_CONNECTION',biochemical_activity_confirmed='NOT_TESTED',historical_I2_score='')])
        (self.review/'stage183_review_provenance.json').write_text(json.dumps({'stage':'STAGE10_PARTIAL_EVIDENCE_REVIEW','all_families':1,'stage10_historical_I2_completed':False,'output_sha256':{rfile.name:sha(rfile)}}))
        stage184 = ROOT/'scripts/11_12_partial_msa_conservation_portable.py'
        prep = subprocess.run([sys.executable,str(stage184),'prepare','--review-dir',str(self.review),'--structure-map',str(mapfile),'--members',str(member),'--fasta',str(fasta),'--out-dir',str(self.prepare),'--limit','1','--mode','run'],capture_output=True,text=True)
        self.assertEqual(prep.returncode, 0, prep.stdout + prep.stderr)
        (self.msa/'mafft').mkdir()
        aligned = self.msa/'mafft/X001.aligned.faa'
        aligned.write_text('>s1|a\nACDE\n>s1|b\nACDF\n>s1|c\nACDG\n')
        a3m = self.msa/'a3m_colabfold/X001.a3m'
        a3m.write_text('>X001\nACDE\n>X001_h0001\nACDF\n>X001_h0002\nACDG\n')
        write_tsv(self.msa/'family_msa_summary.tsv',[dict(structure_id='X001',family_rep='s1|a',family_size='3',unique_sequences='3',collapsed_duplicates='0',msa_status='FAMILY_MSA_READY',query_length_aa='4',a3m_sha256=sha(a3m),mafft_version='fixture')])
        cons = self.root / 'sequence'
        r = subprocess.run([sys.executable,str(stage184),'collect','--prepare-dir',str(self.prepare),'--msa-dir',str(self.msa),'--out-dir',str(cons),'--mode','run'],capture_output=True,text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        stage185 = ROOT/'scripts/13_partial_evidence_catalog_portable.py'
        r = subprocess.run([sys.executable,str(stage185),'--review-dir',str(self.review),'--conservation-dir',str(cons),'--out-dir',str(self.catalog),'--mode','run'],capture_output=True,text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.stage12.mkdir()
        pdb = self.root / 'X001.pdb';pdb.write_text('PDB fixture\n')
        pred = self.root / 'prediction_provenance.json'
        pred.write_text(json.dumps({'status':'COMPLETED','structure_id':'X001','pdb_file':str(pdb),'a3m_sha256':sha(a3m),'pdb_sha256':sha(pdb)}))
        residues = self.stage12 / 'X001_residue_conservation.tsv'
        write_tsv(residues,[dict(structure_id='X001',position=str(i)) for i in range(1,5)])
        summary = self.stage12/'candidate_conservation_summary.tsv'
        write_tsv(summary,[dict(structure_id='X001',length='4',msa_nseq='3',very_conserved_positions='3',priority_structural_positions='1',priority_acidic='1',priority_basic='0',priority_nucleophilic='0')])
        windows = self.stage12/'conserved_windows.tsv'
        write_tsv(windows,[dict(structure_id='X001',start='1',end='2',sequence='AC',conserved_structural='2',charged='0',acidic='0',basic='0',mean_query_fraction='0.9',mean_plddt='80')])
        family = {'structure_id':'X001','alignment':str(aligned),'alignment_sha256':sha(aligned),'a3m':str(a3m),'a3m_sha256':sha(a3m),'prediction_provenance':str(pred),'prediction_provenance_sha256':sha(pred),'pdb_sha256':sha(pdb),'residue_table_sha256':sha(residues),'residue_count':4,'sequence_count':3}
        meta = {'algorithm':'ORIGINAL_STAGE12_FUNCTIONS','final_top5_automatically_selected':False,'family_count':1,'families':[family],'summary_sha256':sha(summary),'windows_sha256':sha(windows)}
        (self.stage12/'conservation_provenance.json').write_text(json.dumps(meta))

    def call(self,mode='check'):
        return subprocess.run([sys.executable,str(SCRIPT),'--catalog-dir',str(self.catalog),'--stage12-dir',str(self.stage12),'--out-dir',str(self.out),'--mode',mode],text=True,capture_output=True)

    def test_success_and_preserves_spatial_status(self):
        x = self.call()
        self.assertEqual(x.returncode, 0, x.stdout + x.stderr)
        self.assertFalse(self.out.exists())
        x = self.call('run')
        self.assertEqual(x.returncode, 0, x.stdout + x.stderr)
        with (self.out/'all_family_stage13_structure_integrated.tsv').open() as f:
            rows=list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['stage12_spatial_state'],'NOT_EVALUATED')
        self.assertEqual(rows[0]['stage12_structure_residue_state'],'STRUCTURE_SUPPORTED_RESIDUES_AVAILABLE')
        self.assertEqual(rows[0]['stage12_conserved_window_count'],'1')
        print('STAGE186_STRUCTURE_INTEGRATION_SMOKE=PASS')

    def test_reject_changed_stage12_summary(self):
        path=self.stage12/'candidate_conservation_summary.tsv'
        path.write_text(path.read_text()+'\n')
        self.assertNotEqual(self.call().returncode,0)

    def test_reject_mismatched_stage12_msa(self):
        meta=self.stage12/'conservation_provenance.json'
        d=json.loads(meta.read_text())
        d['families'][0]['a3m_sha256']='0'*64
        meta.write_text(json.dumps(d))
        self.assertNotEqual(self.call().returncode,0)

    def test_reject_existing_output(self):
        self.out.mkdir()
        self.assertNotEqual(self.call('run').returncode,0)

if __name__ == '__main__':
    unittest.main(verbosity=2)
