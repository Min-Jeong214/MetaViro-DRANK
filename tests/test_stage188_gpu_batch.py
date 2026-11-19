import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

BUNDLE=Path(__file__).resolve().parents[1]
SOURCE=BUNDLE/'scripts'


def tsv(path, cols, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)


class TestFullCohortGPU(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.scripts=self.root/'scripts';self.scripts.mkdir()
        for x in ('11_full_cohort_msa_batch_portable.py','11_gpu_prediction_batch_portable.py'):
            shutil.copy2(BUNDLE/'scripts'/x,self.scripts/x)
        for x in ('11_family_msa_portable.py','11_msa_structure_predict_portable.py'):
            shutil.copy2(SOURCE/x,self.scripts/x)
        self.bin=self.root/'bin';self.bin.mkdir();self.env=os.environ.copy();self.env['PATH']=str(self.bin)+os.pathsep+self.env.get('PATH','')
        mafft=self.bin/'mafft';mafft.write_text('#!/bin/sh\nif [ "$1" = "--version" ]; then echo "mock-mafft 1" >&2; exit 0; fi\nfor arg; do last="$arg"; done\ncat "$last"\n')
        mafft.chmod(0o755)
        colab=self.bin/'colabfold_batch'
        colab.write_text('''#!/usr/bin/env python3
import sys, pathlib, json
inp=pathlib.Path(sys.argv[-2]);out=pathlib.Path(sys.argv[-1]);out.mkdir(exist_ok=True)
sid=inp.stem;seq=''
for line in inp.read_text().splitlines()[1:]:
    if line.startswith('>'):break
    seq+=line.strip()
pdb=out/f'{sid}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb'
with pdb.open('w') as f:
    for i,aa in enumerate(seq,1):
        f.write(f'ATOM  {i:5d}  CA  ALA A{i:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{80.:6.2f}           C\\n')
(out/f'{sid}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json').write_text(json.dumps({'ptm':0.71}))
''')
        colab.chmod(0o755)
        sequences={'R1':'AAAA','H1':'AAAT','R2':'GGGG','H2':'GGGG','R3':'CCCC','H3':'CCCT'}
        self.mapping=[];self.members=[]
        for i,rep in enumerate(('R1','R2','R3'),1):
            self.mapping.append(dict(structure_id=f'MVF{i:06d}',family_rep=rep,family_size='2',length_aa='4',sequence_sha256=hashlib.sha256(sequences[rep].encode()).hexdigest()))
            self.members.extend([{'member':rep,'representative':rep},{'member':'H'+str(i),'representative':rep}])
        self.map=self.root/'map.tsv';self.mem=self.root/'members.tsv';self.faa=self.root/'members.faa';self.queries=self.root/'queries.faa'
        tsv(self.map,list(self.mapping[0]),self.mapping);tsv(self.mem,list(self.members[0]),self.members)
        self.faa.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sequences.items()))
        self.queries.write_text(''.join('>'+r['structure_id']+'\n'+sequences[r['family_rep']]+'\n' for r in self.mapping))
        self.msa=self.root/'msa';self.gpu=self.root/'gpu'

    def cmd(self,name,*args,expect=0):
        p=subprocess.run([sys.executable,str(self.scripts/name),*map(str,args)],env=self.env,capture_output=True,text=True)
        self.assertEqual(p.returncode,expect,p.stdout+'\n'+p.stderr)
        return p.stdout+p.stderr

    def prep(self):
        args=['--members',self.mem,'--fasta',self.faa,'--structure-map',self.map,'--out-dir',self.msa,'--chunk-size','2','--threads','2']
        self.assertIn('FILES_WRITTEN=NO',self.cmd('11_full_cohort_msa_batch_portable.py',*args,'--mode','check'))
        self.assertFalse(self.msa.exists())
        self.assertIn('STAGE188_FULL_COHORT_MSA=PASS',self.cmd('11_full_cohort_msa_batch_portable.py',*args,'--mode','run'))
        self.assertIn('STAGE188_FULL_COHORT_MSA=PASS',self.cmd('11_full_cohort_msa_batch_portable.py',*args,'--mode','run'))

    def gpu_args(self):
        return ['--structure-map',self.map,'--queries',self.queries,'--msa-dir',self.msa,'--out-dir',self.gpu]

    def test_full_msa_and_gpu_resume(self):
        self.prep()
        x=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),'--mode','check')
        self.assertIn('TOTAL_FAMILIES=3',x);self.assertIn('SINGLE_UNIQUE_SEQUENCE',x)
        self.assertFalse(self.gpu.exists())
        x=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),'--mode','run','--max-jobs','1')
        self.assertIn('GPU_LAUNCHED_THIS_INVOCATION=1',x)
        x=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),'--mode','run')
        self.assertIn('GPU_LAUNCHED_THIS_INVOCATION=1',x)
        x=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),'--mode','run')
        self.assertIn('GPU_LAUNCHED_THIS_INVOCATION=0',x)
        with (self.gpu/'prediction_manifest.tsv').open() as manifest:
            m=list(csv.DictReader(manifest,delimiter='\t'))
        self.assertEqual([r['status'] for r in m],['COMPLETED','SINGLE_UNIQUE_SEQUENCE','COMPLETED'])

    def test_strict_reuse_then_explicit_safe_skip(self):
        self.prep()
        self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),'--mode','run','--max-jobs','1')
        old=self.gpu/'families'/'MVF000001'/'attempt_001'
        self.assertTrue((old/'prediction_provenance.json').is_file())
        self.gpu=self.root/'second_gpu'
        normal=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),
                        '--reuse-prediction',old,'--mode','check')
        self.assertIn('GPU_EXISTING_REUSED=1',normal)
        self.assertFalse(self.gpu.exists())
        # A distinct MSA with the exact same query but a different homolog ID.
        a3m=self.msa/'a3m_colabfold'/'MVF000001.a3m'
        orig=a3m.read_text()
        self.assertIn('>MVF000001_h0001',orig)
        a3m.write_text(orig.replace('>MVF000001_h0001','>MVF000001_h0001_renamed'))
        summary=self.msa/'family_msa_summary.tsv'
        with summary.open(newline='') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        for row in rows:
            if row['structure_id']=='MVF000001':
                row['a3m_sha256']=hashlib.sha256(a3m.read_bytes()).hexdigest()
        tsv(summary,list(rows[0]),rows)
        p=subprocess.run([sys.executable,str(self.scripts/'11_gpu_prediction_batch_portable.py'),
                *map(str,self.gpu_args()),'--reuse-prediction',str(old),'--mode','check'],
                env=self.env,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('REUSE_HASH_OR_PDB_INVALID:MVF000001:a3m_sha256',p.stderr)
        x=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),
                '--reuse-prediction',old,'--skip-incompatible-reuse','--mode','check')
        self.assertIn('REUSE_INCOMPATIBLE_SKIPPED=MVF000001:a3m_sha256',x)
        self.assertIn('REUSE_INCOMPATIBLE_SKIPPED_COUNT=1',x)
        self.assertIn('GPU_EXISTING_REUSED=0',x)
        self.assertIn('GPU_PENDING=2',x)
        self.assertFalse(self.gpu.exists())
        x=self.cmd('11_gpu_prediction_batch_portable.py',*self.gpu_args(),
                '--reuse-prediction',old,'--skip-incompatible-reuse',
                '--max-jobs','1','--mode','run')
        self.assertIn('GPU_LAUNCHED_THIS_INVOCATION=1',x)
        self.assertTrue((old/'prediction_provenance.json').is_file())

    def test_hash_tamper_blocks_run(self):
        self.prep()
        q=self.msa/'a3m_colabfold'/'MVF000001.a3m';q.write_text(q.read_text()+'A')
        p=subprocess.run([sys.executable,str(self.scripts/'11_gpu_prediction_batch_portable.py'),*map(str,self.gpu_args()),'--mode','check'],env=self.env,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('MSA_SHA_MISMATCH',p.stderr)
        self.assertFalse(self.gpu.exists())

    def test_chunk_input_tamper_blocks_resume(self):
        self.prep();self.faa.write_text(self.faa.read_text()+'\n')
        p=subprocess.run([sys.executable,str(self.scripts/'11_full_cohort_msa_batch_portable.py'),
              '--members',str(self.mem),'--fasta',str(self.faa),'--structure-map',str(self.map),
              '--out-dir',str(self.msa),'--mode','check'],env=self.env,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('EXISTING_MSA_PLAN_MISMATCH',p.stderr)


if __name__=='__main__': unittest.main()
