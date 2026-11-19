#!/usr/bin/env python3
"""Unit regression: new-cohort MSA selection, checksums, all-family status."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/11_12_partial_msa_conservation_portable.py'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def tsv(file,columns,rows):
    file.parent.mkdir(parents=True,exist_ok=True)
    with file.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


class Checks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='mv184_')
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.review=self.root/'stage10';self.review.mkdir()
        self.rows=[
            dict(structure_id='MVF000001',family_rep='S1|a',family_size='3',review_order='1',review_tier='T2_CONTEXT_SUPPORTED_DARK',priority_class='A_DEEP_DARK',stage09_evidence_state='NOT_EVALUATED_NO_STRUCTURE',structure_state='PREDICTION_PENDING'),
            dict(structure_id='MVF000002',family_rep='S2|a',family_size='1',review_order='2',review_tier='T5_KNOWN_EVIDENCE_REFERENCE',priority_class='D_KNOWN_EVIDENCE',stage09_evidence_state='NOT_EVALUATED_NO_STRUCTURE',structure_state='PREDICTION_PENDING'),
        ]
        reviewfile=self.review/'all_family_stage10_review.tsv'
        tsv(reviewfile,list(self.rows[0]),self.rows)
        (self.review/'stage183_review_provenance.json').write_text(json.dumps({
            'stage':'STAGE10_PARTIAL_EVIDENCE_REVIEW',
            'stage10_historical_I2_completed':False,
            'all_families':2,
            'output_sha256':{reviewfile.name:digest(reviewfile)}}))
        self.seqs={'S1|a':'ACDE','S1|b':'ACDF','S1|c':'ACDG','S2|a':'TTR'}
        with (self.root/'members.faa').open('w') as f:
            for k,v in self.seqs.items():f.write(f'>{k}\n{v}\n')
        self.members=self.root/'members.tsv'
        tsv(self.members,['member','representative'],[
            {'member':member,'representative':('S2|a' if member=='S2|a' else 'S1|a')}
            for member in self.seqs])
        self.mapping=self.root/'structure_id_map.tsv'
        tsv(self.mapping,['structure_id','family_rep','family_size','sequence_sha256'],[
            dict(structure_id='MVF000001',family_rep='S1|a',family_size='3',sequence_sha256=hashlib.sha256(b'ACDE').hexdigest()),
            dict(structure_id='MVF000002',family_rep='S2|a',family_size='1',sequence_sha256=hashlib.sha256(b'TTR').hexdigest()),
        ])
        self.prepare=self.root/'prepared';self.mafft=self.root/'msa';self.conservation=self.root/'conservation'

    def execute(self,command,**kw):
        args=([command,'--review-dir',str(self.review),'--structure-map',str(self.mapping),
               '--members',str(self.members),'--fasta',str(self.root/'members.faa'),
               '--out-dir',str(self.prepare),'--limit','1'] if command=='prepare' else
              [command,'--prepare-dir',str(self.prepare),'--msa-dir',str(self.mafft),
               '--out-dir',str(self.conservation)])
        args+=sum(([k.replace('_','-').join(['--','']),str(v)] for k,v in kw.items()),[])
        return subprocess.run([sys.executable,str(SCRIPT),*args],capture_output=True,text=True)

    def assert_pass(self,result):
        self.assertEqual(result.returncode,0,msg=result.stdout+'\n'+result.stderr)

    def prepared(self):
        self.assert_pass(self.execute('prepare',mode='run'))

    def mock_mafft(self,mutate=False):
        (self.mafft/'mafft').mkdir(parents=True)
        (self.mafft/'a3m_colabfold').mkdir(parents=True)
        aligned=self.mafft/'mafft/MVF000001.aligned.faa'
        aligned.write_text('>S1|a\nACDE\n>S1|b\nACDF\n>S1|c\nACDG\n' if not mutate else '>S1|a\nACDE\n>S1|b\nACDF\n>S1|c\nACDA\n')
        a3m=self.mafft/'a3m_colabfold/MVF000001.a3m'
        a3m.write_text('>MVF000001\nACDE\n>MVF000001_h0001\nACDF\n>MVF000001_h0002\nACDG\n')
        tsv(self.mafft/'family_msa_summary.tsv',
            ['structure_id','family_rep','family_size','unique_sequences','collapsed_duplicates','msa_status','query_length_aa','a3m_sha256','mafft_version'],
            [dict(structure_id='MVF000001',family_rep='S1|a',family_size='3',unique_sequences='3',collapsed_duplicates='0',msa_status='FAMILY_MSA_READY',query_length_aa='4',a3m_sha256=digest(a3m),mafft_version='MAFFT synthetic fixture')])

    def test_new_cohort_without_pilot(self):
        self.assert_pass(self.execute('prepare',mode='check'))
        self.assertFalse(self.prepare.exists())
        self.prepared();self.mock_mafft()
        self.assert_pass(self.execute('collect',mode='check'))
        self.assertFalse(self.conservation.exists())
        self.assert_pass(self.execute('collect',mode='run'))
        with (self.conservation/'all_family_stage11_12_status.tsv').open() as f:
            data=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(data),2)
        self.assertEqual(data[0]['msa_conservation_state'],'SEQUENCE_ONLY_CONSERVATION_AVAILABLE')
        self.assertEqual(data[1]['msa_conservation_state'],'NOT_SELECTED_FOR_MSA')
        with (self.conservation/'selected_family_msa_conservation.tsv').open() as f:
            family=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(family[0]['analyzable_positions'],'4')
        self.assertEqual(family[0]['conserved_positions'],'3')
        self.assertEqual(family[0]['structure_residue_analysis'],'NOT_EVALUATED_NO_MATCHED_PDB_PREDICTION_RECORD')
        print('STAGE184_NEW_DATASET_NO_PILOT=PASS')

    def test_provenance_corruption_rejected(self):
        p=self.review/'all_family_stage10_review.tsv'
        p.write_text(p.read_text().replace('T2_CONTEXT_SUPPORTED_DARK','T4_DARK_CONTEXT_UNRESOLVED'))
        self.assertNotEqual(self.execute('prepare',mode='check').returncode,0)
        print('STAGE184_STAGE183_PROVENANCE_TAMPER_REJECTED=PASS')

    def test_mapping_sequence_tamper_rejected(self):
        p=self.root/'members.faa';p.write_text(p.read_text().replace('ACDE','FCDE'))
        self.assertNotEqual(self.execute('prepare',mode='check').returncode,0)
        print('STAGE184_SEQUENCE_SHA_MISMATCH_REJECTED=PASS')

    def test_msa_alignment_tamper_rejected(self):
        self.prepared();self.mock_mafft(mutate=True)
        self.assertNotEqual(self.execute('collect',mode='check').returncode,0)
        self.assertFalse(self.conservation.exists())
        print('STAGE184_MSA_ALIGNMENT_TAMPER_REJECTED=PASS')

    def test_no_overwrite(self):
        self.prepared();self.assertNotEqual(self.execute('prepare',mode='run').returncode,0)
        self.mock_mafft();self.assert_pass(self.execute('collect',mode='run'))
        self.assertNotEqual(self.execute('collect',mode='run').returncode,0)
        print('STAGE184_OUTPUT_OVERWRITE_REJECTED=PASS')

    def test_source_changed_after_preparation_rejected(self):
        self.prepared();self.mock_mafft()
        with self.mapping.open('a') as f:f.write('\n')
        self.assertNotEqual(self.execute('collect',mode='check').returncode,0)
        print('STAGE184_UPSTREAM_SOURCE_TAMPER_REJECTED=PASS')

    def test_a3m_insert_lengths_are_not_forced_equal(self):
        import runpy
        mod=runpy.run_path(str(SCRIPT),run_name='__test__')
        a3m=self.root/'insertions.a3m'
        a3m.write_text('>MVF000001\nACDE\n>MVF000001_h0001\nAgCDE\n')
        records=mod['aligned_fasta'](a3m,require_equal_width=False)
        self.assertEqual(len(records),2)
        with self.assertRaises(ValueError):
            mod['aligned_fasta'](a3m)
        print('STAGE184_VARIABLE_A3M_INSERTIONS_SUPPORTED=PASS')

    def test_zero_limit_rejected(self):
        self.assertNotEqual(self.execute('prepare',mode='check',limit=0).returncode,0)
        print('STAGE184_INVALID_SELECTION_LIMIT_REJECTED=PASS')


if __name__=='__main__':unittest.main(verbosity=1)
