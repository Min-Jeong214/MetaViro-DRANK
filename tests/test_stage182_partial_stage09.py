#!/usr/bin/env python3
"""Offline Stage182 contract tests; no Foldseek, GPU, or Pilot data."""
import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]

def module(filename,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/filename)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

s182=module('09_partial_structure_evidence_portable.py','s182')
s09=module('09_external_integration_portable.py','s09')


def write_tsv(path,rows,columns):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)

class Stage182Test(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='mv_182_')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.plan=self.root/'plan';self.plan.mkdir()
        self.queries=self.root/'queries.faa'
        self.seqs={'MVF000001':'ACD','MVF000002':'MEK'}
        self.queries.write_text(''.join(f'>{k}\n{v}\n' for k,v in self.seqs.items()))
        self.model_pdb=self.root/'original_model.pdb'
        self.model_score=self.root/'score.json'
        self.model_score.write_text(json.dumps({'plddt':[80,80,80],'ptm':0.5}))
        self.model_pdb.write_text(''.join(
            f'ATOM  {i:5d}  CA  {aa} A{i:4d}    {i:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{80.0:6.2f}           C\n'
            for i,aa in enumerate(['ALA','CYS','ASP'],1))+'END\n')
        self.rows=[]
        for i,(sid,seq) in enumerate(self.seqs.items(),1):
            self.rows.append({
                'structure_id':sid,'family_rep':f'R{i}','priority_class':'A_DEEP_DARK',
                'family_size':'2','length_aa':str(len(seq)),
                'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                'structure_status':'NO_STRUCTURE_AVAILABLE',
                'historical_structure_id':'','historical_pdb_file':'',
                'historical_pdb_sha256':'','historical_score_json':'',
                'historical_score_json_sha256':'','historical_mean_plddt':'',
                'historical_ptm':'',
            })
        self.inventory=self.plan/'partial_structure_inventory.tsv'
        write_tsv(self.inventory,self.rows,list(self.rows[0]))
        self.planprov=self.plan/'partial_structure_plan_provenance.json'
        self.update_plan_prov()
        self.pred=self.root/'new_stage08.tsv'
        self.new_model={
            'structure_id':'MVF000001','family_rep':'R1','source_member':'R1',
            'family_size':'2','priority_class':'A_DEEP_DARK','length_aa':'3',
            'sequence_sha256':self.rows[0]['sequence_sha256'],
            'pdb_file':str(self.model_pdb),'pdb_sha256':s182.sha(self.model_pdb),
            'score_json':str(self.model_score),'score_json_sha256':s182.sha(self.model_score),
            'mean_plddt':'80.0','min_plddt':'80.0','max_plddt':'80.0',
            'ptm':'0.5','confidence_class':'CONFIDENT','prediction_status':'PREDICTED',
        }
        write_tsv(self.pred,[self.new_model],s182.STAGE08_COLUMNS)

    def update_plan_prov(self):
        self.planprov.write_text(json.dumps({'output_sha256':{
            'partial_structure_inventory.tsv':s182.sha(self.inventory)
        }}))

    def call_prepare(self,out,mode='run',new=True):
        argv=['cmd','prepare','--plan-dir',str(self.plan),'--queries',str(self.queries),
              '--out-dir',str(out),'--mode',mode]
        if new:argv+=['--predicted-structures',str(self.pred)]
        with patch('sys.argv',argv):s182.main()

    def call_merge(self,prep,out,ext=None,mode='run'):
        argv=['cmd','merge','--prepared-dir',str(prep),'--out-dir',str(out),'--mode',mode]
        if ext:argv+=['--external-dir',str(ext)]
        with patch('sys.argv',argv):s182.main()

    def fake_external(self,prep):
        external=self.root/'external';external.mkdir()
        model=s09.load_structures(prep/'structure_results_for_search.tsv')
        search={}
        for pre,_,profile,max_seqs in s09.SOURCES:
            directory=self.root/'search'/profile;directory.mkdir(parents=True)
            hits=[]
            if pre=='bfvd':
                hit={field:'1' for field in s09.HIT_FIELDS}
                hit.update(query='MVF000001',target='ARTIFICIAL_BFVD_TARGET',
                           evalue='1e-20',prob='0.99',lddt='0.9',qcov='1.0',tcov='0.9',
                           qtmscore='0.8',ttmscore='0.8',bits='250',fident='0.45',
                           alnlen='3',qlen='3',tlen='3',alntmscore='0.8')
                hits=[hit]
            raw=directory/'raw_hits.tsv'
            with raw.open('w',newline='') as f:
                w=csv.writer(f,delimiter='\t',lineterminator='\n')
                for h in hits:w.writerow([h[k] for k in s09.HIT_FIELDS])
            manifest=directory/'query_manifest.tsv'
            write_tsv(manifest,[{'structure_id':'MVF000001',
               'sha256':self.new_model['pdb_sha256']}],['structure_id','sha256'])
            prov=directory/'search_provenance.json'
            prov.write_text(json.dumps({'source_profile':profile,'control_only':False,
                'evalue':1.0,'sensitivity':9.5,'max_seqs':max_seqs,
                'raw_hits_sha256':s182.sha(raw),'foldseek_version':'fixture'}))
            search[pre]=defaultdict(list,{'MVF000001':hits})
        matrix=s09.integrate(model,search)
        mpath=external/'external_structure_matrix.tsv'
        write_tsv(mpath,matrix,s09.HEADER)
        pv={'evidence_origin':'COMPLETED_EXTERNAL_FOLDSEEK_SEARCHES',
            'structure_table_sha256':s182.sha(prep/'structure_results_for_search.tsv'),
            'matrix_sha256':s182.sha(mpath),'family_count':1,
            'external_searches':{}}
        for pre,_,profile,max_seqs in s09.SOURCES:
            folder=self.root/'search'/profile
            pv['external_searches'][pre]={
                'directory':str(folder),'source_profile':profile,
                'raw_hits_sha256':s182.sha(folder/'raw_hits.tsv'),
                'query_manifest_sha256':s182.sha(folder/'query_manifest.tsv'),
                'search_provenance_sha256':s182.sha(folder/'search_provenance.json'),
            }
        (external/'integration_provenance.json').write_text(json.dumps(pv))
        return external

    def test_complete_mixed_and_search_verification(self):
        prep=self.root/'prepare'
        self.call_prepare(prep,'check')
        self.assertFalse(prep.exists())
        self.call_prepare(prep)
        _,initial=s182.table(prep/'all_family_structure_status.tsv')
        self.assertEqual([r['structure_state'] for r in initial],['PREDICTED_VERIFIED','PREDICTION_PENDING'])
        external=self.fake_external(prep)
        dest=self.root/'merged'
        self.call_merge(prep,dest,external,'check')
        self.assertFalse(dest.exists())
        self.call_merge(prep,dest,external)
        _,result=s182.table(dest/'all_family_stage09_evidence.tsv')
        self.assertEqual(len(result),2)
        self.assertEqual(result[0]['stage09_external_status'],'CONVINCING_STRUCTURE_CONNECTION')
        self.assertEqual(result[1]['stage09_evidence_state'],'NOT_EVALUATED_NO_STRUCTURE')
        self.assertEqual(result[1]['stage09_external_status'],'')
        self.assertEqual(json.loads((dest/'stage182_merge_provenance.json').read_text())['prediction_pending_count'],1)
        print('STAGE182_MIXED_NEW_PDB_AND_PENDING=PASS')

    def test_missing_structure_all_pending(self):
        prep=self.root/'none'
        self.call_prepare(prep,new=False)
        dest=self.root/'all_pending'
        self.call_merge(prep,dest)
        _,rows=s182.table(dest/'all_family_stage09_evidence.tsv')
        self.assertTrue(all(not r['stage09_external_status'] for r in rows))
        print('STAGE182_NO_PILOT_NO_PDB_WORKS=PASS')

    def test_historical_exact_reuse(self):
        row=self.rows[0]
        row.update(structure_status='VERIFIED_HISTORICAL_REUSE_CANDIDATE',
                   historical_structure_id='OLDX',historical_pdb_file=str(self.model_pdb),
                   historical_pdb_sha256=s182.sha(self.model_pdb),
                   historical_score_json=str(self.model_score),
                   historical_score_json_sha256=s182.sha(self.model_score),
                   historical_mean_plddt='80.0',historical_ptm='0.5')
        write_tsv(self.inventory,self.rows,list(self.rows[0]))
        self.update_plan_prov()
        prep=self.root/'old'
        self.call_prepare(prep,new=False)
        _,rows=s182.table(prep/'structure_results_for_search.tsv')
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['structure_id'],'MVF000001')
        self.assertEqual(rows[0]['pdb_sha256'],s182.sha(self.model_pdb))
        print('STAGE182_HISTORICAL_REUSE_IS_EXACT_SEQUENCE_ONLY=PASS')

    def test_malformed_model_and_tampered_sources_rejected(self):
        wrong=dict(self.new_model,sequence_sha256='0'*64)
        write_tsv(self.pred,[wrong],s182.STAGE08_COLUMNS)
        with self.assertRaisesRegex(ValueError,'MODEL_SEQUENCE_SHA_MISMATCH'):
            self.call_prepare(self.root/'wrong')
        self.assertFalse((self.root/'wrong').exists())
        write_tsv(self.pred,[self.new_model],s182.STAGE08_COLUMNS)
        self.model_pdb.write_text('TAMPERED')
        with self.assertRaisesRegex(ValueError,'SOURCE_SHA256_MISMATCH'):
            self.call_prepare(self.root/'tampered')
        self.assertFalse((self.root/'tampered').exists())
        print('STAGE182_SEQUENCE_AND_PDB_TAMPER_REJECTED=PASS')

    def test_external_raw_tamper_and_overwrite_rejected(self):
        prep=self.root/'prepare';self.call_prepare(prep)
        external=self.fake_external(prep)
        (self.root/'search/bfvd/raw_hits.tsv').write_text('TAMPER')
        with self.assertRaisesRegex(ValueError,'Raw hit SHA mismatch'):
            self.call_merge(prep,self.root/'merged',external)
        self.assertFalse((self.root/'merged').exists())
        self.call_prepare(self.root/'overwrite')
        with self.assertRaisesRegex(ValueError,'REFUSE_EXISTING_OUTPUT'):
            self.call_prepare(self.root/'overwrite')
        print('STAGE182_EXTERNAL_RAW_TAMPER_AND_OVERWRITE_REJECTED=PASS')

if __name__=='__main__':
    unittest.main(verbosity=0)
    print('ALL_TESTS_PASS')
