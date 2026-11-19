import csv, hashlib, json, subprocess, tempfile, shutil
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
SCR=BASE/'scripts'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,t):
 p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(t)
 return p
def tsv(p, fields, records):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f, fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
 return p
def js(p,obj):write(p,json.dumps(obj,indent=2)+'\n');return p
def call(name,args,expect=0):
 cmd=['python','-B',str(SCR/name),*map(str,args)]
 p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,cwd=BASE)
 print('CALL',name,'EXIT',p.returncode, '\n'+p.stdout[-1200:])
 if expect == 0 and p.returncode:raise RuntimeError(p.stdout)
 if expect != 0 and p.returncode == 0:raise AssertionError('Expected failure')
 return p.stdout

def pdbtext(n):
 lines=[]
 for i in range(1,n+1):
  line=f"ATOM  {i:5d}  CA  ALA A{i:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{72.:6.2f}           C\n"
  assert line[12:16].strip()=='CA', repr(line)
  assert line[21]=='A' and int(line[22:26])==i
  assert float(line[60:66])==72.
  lines.append(line)
 return ''.join(lines)+'END\n'

with tempfile.TemporaryDirectory(prefix='mvstage158_') as t:
 root=Path(t)
 up=root/'stage08'; down=root/'ranking'; prep=root/'prepared'; hist=root/'historical'
 members=tsv(up/'05_annotation/recurrent_members.tsv',['member','representative'],[
  {'member':'rep1','representative':'rep1'}, {'member':'m1','representative':'rep1'},
  {'member':'rep2','representative':'rep2'}, {'member':'m2','representative':'rep2'},
 ])
 faa=write(up/'05_annotation/recurrent_members.faa','>rep1\nAAAA\n>m1\nAAA\n>rep2\nACDE\n>m2\nACE\n')
 structures=[]; maprows=[]; prep_summary=[]; policy=[]; inp={}; selected=[]; reg=[]
 for sid,rep,seq,other in [('F001','rep1','AAAA','AAA'),('F002','rep2','ACDE','ACE')]:
  seqsha=hashlib.sha256(seq.encode()).hexdigest()
  sf=write(up/f'08_structure/run_test/collected/{sid}.pdb',pdbtext(4))
  score=js(up/f'08_structure/run_test/collected/{sid}.json',{'plddt':[72.]*4,'ptm':0.3})
  structures.append(dict(structure_id=sid,family_rep=rep,sequence_sha256=seqsha,length_aa='4',prediction_status='PREDICTED',mean_plddt='72.00',ptm='0.3',pdb_file=str(sf),pdb_sha256=sha(sf),score_json=str(score),score_json_sha256=sha(score)))
  maprows.append(dict(structure_id=sid,family_rep=rep,family_size='2',length_aa='4',sequence_sha256=seqsha))
  # query and member with valid A3M columns count 4
  other_aligned='A-AA' if sid=='F001' else 'AC-E'
  prepared_a3m=write(prep/f'msa/a3m_colabfold/{sid}.a3m',f'>{sid}\n{seq}\n>m_{sid}\n{other_aligned}\n')
  prep_summary.append(dict(structure_id=sid,family_rep=rep,family_size='2',unique_sequences='2',a3m_sha256=sha(prepared_a3m)))
  policy.append(dict(structure_id=sid,family_rep=rep,family_members='2',unique_sequences='2',msa_prediction_status='MSA_PREDICTION_ELIGIBLE'))
  inp[sid]=dict(sha256=sha(prepared_a3m), msa_prediction_eligible=True,msa_prediction_status='MSA_PREDICTION_ELIGIBLE')
  selected.append({'structure_id':sid,'family_rep':rep,**{f'field{i}':i for i in range(28)}})
  if sid=='F001':
   old_a3m=write(hist/f'{sid}.a3m',f'#4,2\n>{sid}\n{seq}\n>m_{sid}\nAA-A\n')
   old_pdb=write(hist/f'{sid}.pdb',pdbtext(4))
   old_score=js(hist/f'{sid}.json',{'plddt':[72.]*4,'ptm':0.3})
   reg.append(dict(structure_id=sid,reuse_class='HISTORICAL_MSA_REUSE_CANDIDATE',new_a3m_sha256=sha(prepared_a3m),historical_a3m=str(old_a3m),historical_a3m_sha256=sha(old_a3m),historical_pdb=str(old_pdb),historical_pdb_sha256=sha(old_pdb),historical_score_json=str(old_score),historical_score_sha256=sha(old_score)))
  else:
   reg.append(dict(structure_id=sid,reuse_class='PREDICTION_PENDING',new_a3m_sha256=sha(prepared_a3m),historical_a3m='',historical_a3m_sha256='',historical_pdb='',historical_pdb_sha256='',historical_score_json='',historical_score_sha256=''))
 structpath=tsv(up/'08_structure/run_test/collected/structure_results.tsv',list(structures[0]),structures)
 selpath=tsv(prep/'selected_structure_map.tsv',list(maprows[0]),maprows)
 summpath=tsv(prep/'msa/family_msa_summary.tsv',list(prep_summary[0]),prep_summary)
 eligpath=tsv(prep/'prediction_eligibility.tsv',list(policy[0]),policy)
 fields=list(selected[0]); rank=tsv(down/'10_ranking/confidence_aware_ranking.tsv',fields,selected)
 short=tsv(down/'10_shortlists/top15_experimental_candidates.tsv',fields,selected)
 js(down/'10_ranking/ranking_provenance.json',{'ranking_sha256':sha(rank)})
 js(down/'orchestration_provenance.json',{'stage08_structure_table_sha256':sha(structpath),'artifacts':{'ranking':{'sha256':sha(rank),'path':str(rank)},'top15':{'sha256':sha(short),'path':str(short)}}})
 js(prep/'preparation_provenance.json',dict(stage10_ranking_sha256=sha(rank),stage10_shortlist_sha256=sha(short),selected_structure_map_sha256=sha(selpath),msa_summary_sha256=sha(summpath),prediction_eligibility_sha256=sha(eligpath),family_members_sha256=sha(members),member_fasta_sha256=sha(faa),stage08_structures_sha256=sha(structpath),families=2,msa_prediction_eligible_families=2,single_unique_sequence_families=0,colabfold_inputs=inp))
 register=tsv(root/'register.tsv',list(reg[0]),reg)
 call('11_import_historical_reuse_portable.py',['--reuse-register',register,'--prepared-dir',prep,'--stage08-output',up,'--out-dir',root/'check','--mode','check'])
 run=root/'workflow'
 call('11_13_adaptive_orchestrator_portable.py',['--stage08-output',up,'--stage09-10-output',down,'--prepared-dir',prep,'--reuse-register',register,'--out-dir',run,'--mode','run'])
 with (run/'12_13_status/stage10_ranked_evaluation_status.tsv').open() as fh:
  output={r['structure_id']:r for r in csv.DictReader(fh,delimiter='\t')}
 assert output['F001']['stage11_msa_qc_status']=='QC_VERIFIED_HISTORICAL',output
 assert output['F002']['stage11_collection_status']=='MSA_PREDICTION_PENDING',output
 assert output['F002']['stage13_integrated_evidence_status']=='STAGE10_ONLY_PREDICTION_PENDING',output
 assert len(output)==2
 assert json.loads((run/'orchestration_provenance.json').read_text())['stage12_13_analysis_completed'] is False
 print('PARTIAL_E2E_PASS')

 # STEP162: verified historical A3M drives original conservation and exploratory 3D clusters.
 target=root/'stage162_output'
 call('12_13_partial_historical_portable.py',[
  '--adaptive-output',run,'--stage09-10-output',down,'--out-dir',target,'--mode','check'])
 assert not target.exists()
 call('12_13_partial_historical_portable.py',[
  '--adaptive-output',run,'--stage09-10-output',down,'--out-dir',target,'--mode','run'])
 with (target/'stage12_13_partial_status.tsv').open() as fh:
  stage162={r['structure_id']:r for r in csv.DictReader(fh,delimiter='\t')}
 assert stage162['F001']['stage12_conservation_status']=='COMPLETED'
 assert stage162['F001']['stage12_spatial_status']=='COMPLETED_EXPLORATORY'
 assert stage162['F002']['stage12_conservation_status']=='PENDING_PREDICTION'
 assert stage162['F002']['stage13_integrated_evidence_status']=='STAGE10_ONLY_PREDICTION_PENDING'
 assert len(list(csv.DictReader((target/'12_conservation/F001_residue_conservation.tsv').open(),delimiter='\t'))) == 4
 with (target/'13_preliminary/stage10_plus_partial_evidence.tsv').open() as fh:
  evidence={r['structure_id']:r for r in csv.DictReader(fh,delimiter='\t')}
 assert evidence['F002']['preliminary_very_conserved_positions']==''
 assert evidence['F001']['preliminary_very_conserved_positions']!=''
 assert json.loads((target/'orchestration_provenance.json').read_text())['stage13_full_integration_completed'] is False
 print('STAGE162_SYNTHETIC_PARTIAL_TEST=PASS')
 # Confirm old strict full-index flow remains operational without --allow-partial.
 freshrows=[]
 for r in structures:
  sid=r['structure_id']
  a3m=prep/f'msa/a3m_colabfold/{sid}.a3m'
  provenance=js(root/f'{sid}.fresh.json',dict(
   status='COMPLETED',structure_id=sid,model_type='alphafold2_ptm',
   query_length=4,sequence_count=2,query_sha256=r['sequence_sha256'],
   source_a3m=str(a3m),a3m_sha256=sha(a3m),
   pdb_file=r['pdb_file'],pdb_sha256=r['pdb_sha256'],
   score_json=r['score_json'],score_json_sha256=r['score_json_sha256'],
   mean_plddt=72.00,ptm=0.3,colabfold_returncode=0,
  ))
  freshrows.append({'structure_id':sid,'prediction_provenance':str(provenance)})
 freshindex=tsv(root/'fresh_index.tsv',['structure_id','prediction_provenance'],freshrows)
 full=root/'full_workflow'
 call('11_13_adaptive_orchestrator_portable.py',[
  '--stage08-output',up,'--stage09-10-output',down,
  '--prepared-dir',prep,'--prediction-index',freshindex,
  '--out-dir',full,'--mode','run'])
 with (full/'12_13_status/stage10_ranked_evaluation_status.tsv').open() as fh:
  statuses=list(csv.DictReader(fh,delimiter='\t'))
 assert len(statuses)==2 and all(row['stage11_msa_qc_status']=='QC_VERIFIED' for row in statuses)
 print('STRICT_FULL_E2E_PASS')
 # tampered historical artifact must be detected before creating an output.
 old_a3m.write_text(old_a3m.read_text()+'\n')
 call('11_import_historical_reuse_portable.py',['--reuse-register',register,'--prepared-dir',prep,'--stage08-output',up,'--out-dir',root/'bad','--mode','check'],expect=1)
 assert not (root/'bad').exists()
 print('TAMPER_NEGATIVE_TEST_PASS')
 # Stage162 must refuse stale/tampered historical predictions as well.
 call('12_13_partial_historical_portable.py',[
  '--adaptive-output',run,'--stage09-10-output',down,
  '--out-dir',root/'tampered_stage162','--mode','check'],expect=1)
 assert not (root/'tampered_stage162').exists()
 print('STAGE162_TAMPER_NEGATIVE_TEST=PASS')

 # subset not accepted by default (compatibility with prior strict behavior).
 partialidx=run/'11_historical_import/prediction_index.tsv'
 call('11_collect_orchestrator_portable.py',['--stage08-output',up,'--stage09-10-output',down,'--prepared-dir',prep,'--prediction-index',partialidx,'--out-dir',root/'strict_bad','--mode','check'],expect=1)
 print('STRICT_MODE_NEGATIVE_TEST_PASS')
print('ALL_TESTS_PASS')
