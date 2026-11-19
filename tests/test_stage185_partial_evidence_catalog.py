#!/usr/bin/env python3
"""Integration regressions for generic Stage183 -> Stage184 -> Stage185 catalog."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STAGE184=ROOT/'scripts/11_12_partial_msa_conservation_portable.py'
STAGE185=ROOT/'scripts/13_partial_evidence_catalog_portable.py'


def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()


def write_tsv(file,rows):
    file.parent.mkdir(parents=True,exist_ok=True)
    with file.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


class Stage185Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='mv185_')
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.review=self.root/'review';self.review.mkdir()
        self.prep=self.root/'prep';self.msa=self.root/'msa'
        self.cons=self.root/'conservation';self.out=self.root/'out'
        self.reviews=[
            dict(review_order='1',review_tier='T2_CONTEXT_SUPPORTED_DARK',structure_id='X001',family_rep='S1|a',family_size='3',priority_class='A_DEEP_DARK',sample_span='2',genomic_context_state='DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT',best_noncandidate_anchor_category='nuclease',category_supported_samples='2',category_supported_members='3',stage09_evidence_state='NOT_EVALUATED_NO_STRUCTURE',stage09_external_status='',biochemical_activity_confirmed='NOT_TESTED',historical_I2_score=''),
            dict(review_order='2',review_tier='T4_DARK_CONTEXT_UNRESOLVED',structure_id='X002',family_rep='S2|a',family_size='1',priority_class='B_FUNCTION_UNRESOLVED',sample_span='1',genomic_context_state='DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY',best_noncandidate_anchor_category='',category_supported_samples='',category_supported_members='',stage09_evidence_state='NOT_EVALUATED_NO_STRUCTURE',stage09_external_status='',biochemical_activity_confirmed='NOT_TESTED',historical_I2_score=''),
            dict(review_order='3',review_tier='T5_KNOWN_EVIDENCE_REFERENCE',structure_id='X003',family_rep='S3|a',family_size='1',priority_class='D_KNOWN_EVIDENCE',sample_span='1',genomic_context_state='KNOWN_EVIDENCE_FAMILY',best_noncandidate_anchor_category='',category_supported_samples='',category_supported_members='',stage09_evidence_state='NOT_EVALUATED_NO_STRUCTURE',stage09_external_status='',biochemical_activity_confirmed='NOT_TESTED',historical_I2_score=''),
        ]
        rf=self.review/'all_family_stage10_review.tsv';write_tsv(rf,self.reviews)
        (self.review/'stage183_review_provenance.json').write_text(json.dumps({'stage':'STAGE10_PARTIAL_EVIDENCE_REVIEW','all_families':3,'stage10_historical_I2_completed':False,'output_sha256':{rf.name:sha(rf)}}))
        self.faa=self.root/'proteins.faa';self.faa.write_text('>S1|a\nACDE\n>S1|b\nACDF\n>S1|c\nACDG\n>S2|a\nTTR\n>S3|a\nGGA\n')
        self.members=self.root/'members.tsv'
        write_tsv(self.members,[dict(member=k,representative=rep) for k,rep in [('S1|a','S1|a'),('S1|b','S1|a'),('S1|c','S1|a'),('S2|a','S2|a'),('S3|a','S3|a')]])
        self.mapping=self.root/'mapping.tsv'
        write_tsv(self.mapping,[dict(structure_id=sid,family_rep=rep,family_size=n,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()) for sid,rep,n,seq in [('X001','S1|a','3','ACDE'),('X002','S2|a','1','TTR'),('X003','S3|a','1','GGA')]])
        self.good(self.exec184('prepare','--mode','run'))
        (self.msa/'mafft').mkdir(parents=True)
        (self.msa/'a3m_colabfold').mkdir(parents=True)
        (self.msa/'mafft/X001.aligned.faa').write_text('>S1|a\nACDE\n>S1|b\nACDF\n>S1|c\nACDG\n')
        a3m=self.msa/'a3m_colabfold/X001.a3m'
        a3m.write_text('>X001\nACDE\n>X001_h0001\nACDF\n>X001_h0002\nACDG\n')
        write_tsv(self.msa/'family_msa_summary.tsv',[dict(structure_id='X001',family_rep='S1|a',family_size='3',unique_sequences='3',collapsed_duplicates='0',msa_status='FAMILY_MSA_READY',query_length_aa='4',a3m_sha256=sha(a3m),mafft_version='test fixture')])
        self.good(self.exec184('collect','--mode','run'))

    def run_cmd(self,script,*args):
        return subprocess.run([sys.executable,'-B',str(script),*map(str,args)],text=True,capture_output=True)

    def exec184(self,stage,*extra):
        if stage=='prepare':
            return self.run_cmd(STAGE184,stage,'--review-dir',self.review,'--structure-map',self.mapping,'--members',self.members,'--fasta',self.faa,'--out-dir',self.prep,'--limit',1,*extra)
        return self.run_cmd(STAGE184,stage,'--prepare-dir',self.prep,'--msa-dir',self.msa,'--out-dir',self.cons,*extra)

    def execute(self,mode='check',limit=2):
        return self.run_cmd(STAGE185,'--review-dir',self.review,'--conservation-dir',self.cons,'--out-dir',self.out,'--limit',limit,'--mode',mode)

    def good(self,result):
        self.assertEqual(result.returncode,0,msg=result.stdout+'\n'+result.stderr)

    def refresh_conservation_meta(self,path):
        meta=self.cons/'stage184_conservation_provenance.json';obj=json.loads(meta.read_text())
        obj['output_sha256'][path.name]=sha(path)
        meta.write_text(json.dumps(obj))

    def test_new_cohort_full_family_catalog(self):
        self.good(self.execute('check'))
        self.assertFalse(self.out.exists())
        self.good(self.execute('run'))
        with (self.out/'all_family_stage13_evidence.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(rows),3)
        self.assertEqual(rows[0]['msa_conserved_positions'],'3')
        self.assertEqual(rows[0]['msa_conservation_state'],'SEQUENCE_ONLY_CONSERVATION_AVAILABLE')
        self.assertEqual(rows[1]['msa_conservation_state'],'NOT_SELECTED_FOR_MSA')
        self.assertEqual(rows[1]['stage09_evidence_state'],'NOT_EVALUATED_NO_STRUCTURE')
        self.assertTrue(all(r['stage12_spatial_state']=='NOT_EVALUATED' for r in rows))
        self.assertTrue(all(r['biochemical_activity_confirmed']=='NOT_TESTED' for r in rows))
        with (self.out/'context_supported_stage13_review.tsv').open() as f:
            self.assertEqual(len(list(csv.DictReader(f,delimiter='\t'))),1)
        meta=json.loads((self.out/'stage185_provisional_provenance.json').read_text())
        self.assertFalse(meta['final_candidate_selection_completed'])
        print('STAGE185_NEW_COHORT_NO_PILOT=PASS')

    def test_output_overwrite_rejected(self):
        self.good(self.execute('run'))
        self.assertNotEqual(self.execute('run').returncode,0)
        print('STAGE185_OUTPUT_OVERWRITE_REJECTED=PASS')

    def test_stage183_sha_tampering_rejected(self):
        rf=self.review/'all_family_stage10_review.tsv'
        rf.write_text(rf.read_text().replace('nuclease','helicase'))
        self.assertNotEqual(self.execute().returncode,0)
        print('STAGE185_UPSTREAM_STAGE183_TAMPER_REJECTED=PASS')

    def test_stage184_sha_tampering_rejected(self):
        path=self.cons/'all_family_stage11_12_status.tsv'
        path.write_text(path.read_text()+'\n')
        self.assertNotEqual(self.execute().returncode,0)
        print('STAGE185_STAGE184_SHA_TAMPER_REJECTED=PASS')

    def test_cross_cohort_id_mismatch_rejected_even_if_rehashed(self):
        path=self.cons/'all_family_stage11_12_status.tsv'
        path.write_text(path.read_text().replace('S2|a','different|protein'))
        self.refresh_conservation_meta(path)
        self.assertNotEqual(self.execute().returncode,0)
        print('STAGE185_CROSS_COHORT_ID_REJECTED=PASS')

    def test_fake_conserved_position_rejected_even_if_rehashed(self):
        path=self.cons/'selected_family_msa_conservation.tsv'
        path.write_text(path.read_text().replace('\t3\tSEQUENCE_ONLY_CONSERVATION_AVAILABLE','\t4\tSEQUENCE_ONLY_CONSERVATION_AVAILABLE'))
        self.refresh_conservation_meta(path)
        self.assertNotEqual(self.execute().returncode,0)
        print('STAGE185_RESIDUE_SUMMARY_MISMATCH_REJECTED=PASS')

    def test_final_claim_and_invalid_limit_rejected(self):
        self.assertNotEqual(self.execute('check',limit=0).returncode,0)
        p=self.cons/'stage184_conservation_provenance.json';meta=json.loads(p.read_text())
        meta['structure_based_stage12_completed']=True
        p.write_text(json.dumps(meta))
        self.assertNotEqual(self.execute().returncode,0)
        print('STAGE185_FALSE_FULL_STRUCTURE_STATUS_REJECTED=PASS')

    def test_stage184_prepare_source_tamper_rejected(self):
        self.mapping.write_text(self.mapping.read_text()+'\n')
        self.assertNotEqual(self.execute().returncode,0)
        print('STAGE185_STAGE184_PREPARE_SOURCE_TAMPER_REJECTED=PASS')


if __name__=='__main__':unittest.main(verbosity=1)
