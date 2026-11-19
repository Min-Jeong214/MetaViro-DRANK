#!/usr/bin/env python3
"""Stage183 portable review regression tests (no GPU, Foldseek or user data)."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/10_partial_evidence_triage_portable.py'


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def write(file, rows):
    file.parent.mkdir(parents=True, exist_ok=True)
    with file.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read(file):
    with file.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


class Stage183(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='metaviro183_')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ctx = self.root / 'context'
        self.ev = self.root / 'merged'
        self.prep = self.root / 'prepare'
        for folder in (self.ctx, self.ev, self.prep):
            folder.mkdir()
        self.ctx_rows = []
        self.ev_rows = []
        self.status_rows = []
        specs = [
            ('A_DEEP_DARK', 'DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT', 'SEARCH_VERIFIED', 'CONVINCING_STRUCTURE_CONNECTION'),
            ('B_FUNCTION_UNRESOLVED', 'DARK_WITH_NONCANDIDATE_REPLICATION_CONTEXT', 'NOT_EVALUATED_NO_STRUCTURE', ''),
            ('A_DEEP_DARK', 'DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY', 'SEARCH_VERIFIED', 'WEAK_RAW_HIT_ONLY'),
            ('C_KNOWN_DARK_BRIDGE', 'DARK_WITHOUT_SUPPORTED_ENZYME_CATEGORY', 'SEARCH_VERIFIED', 'NO_RAW_HIT_E1'),
            ('D_KNOWN_EVIDENCE', 'KNOWN_EVIDENCE_FAMILY', 'NOT_EVALUATED_NO_STRUCTURE', ''),
        ]
        for idx, (klass, context_state, search_state, external_status) in enumerate(specs, 1):
            sid = f'MVF{idx:06d}'
            supported = context_state.endswith('REPLICATION_CONTEXT')
            searched = search_state == 'SEARCH_VERIFIED'
            self.ctx_rows.append(dict(
                context_review_order=str(idx), structure_id=sid,
                family_rep=f'V{idx}|rep', priority_class=klass,
                family_size='3', sample_span='2',
                best_noncandidate_anchor_category='nuclease' if supported else '',
                category_supported_members='2' if supported else '',
                category_supported_samples='2' if supported else '',
                noncandidate_context_sample_span='2' if supported else '0',
                noncandidate_anchor_context_count='2' if supported else '0',
                screening_status=context_state,
                screening_interpretation='CONTEXT_ASSOCIATION_NOT_CONFIRMED_ENZYME_ACTIVITY',
            ))
            self.ev_rows.append(dict(
                structure_id=sid, family_rep=f'V{idx}|rep', priority_class=klass,
                family_size='3',
                structure_state='PREDICTED_VERIFIED' if searched else 'PREDICTION_PENDING',
                external_search_state='READY_FOR_SEARCH' if searched else 'NOT_EVALUATED_NO_STRUCTURE',
                stage09_evidence_state=search_state,
                stage09_external_status=external_status,
                stage09_raw_sources='BFVD' if external_status in ('CONVINCING_STRUCTURE_CONNECTION','WEAK_RAW_HIT_ONLY') else '',
                stage09_convincing_sources='BFVD' if external_status == 'CONVINCING_STRUCTURE_CONNECTION' else '',
            ))
            self.status_rows.append(dict(
                structure_id=sid, family_rep=f'V{idx}|rep',
                structure_state='PREDICTED_VERIFIED' if searched else 'PREDICTION_PENDING',
                external_search_state='READY_FOR_SEARCH' if searched else 'NOT_EVALUATED_NO_STRUCTURE',
            ))
        self.ctxfile = self.ctx / 'all_family_context_screen.tsv'
        self.evfile = self.ev / 'all_family_stage09_evidence.tsv'
        self.statusfile = self.prep / 'all_family_structure_status.tsv'
        self.modelfile = self.prep / 'structure_results_for_search.tsv'
        write(self.ctxfile, self.ctx_rows)
        write(self.evfile, self.ev_rows)
        write(self.statusfile, self.status_rows)
        write(self.modelfile, [dict(structure_id='placeholder')])
        self.update_provenance()

    def update_provenance(self):
        (self.ctx / 'context_first_provenance.json').write_text(json.dumps(dict(
            analysis_type='STRUCTURE_OPTIONAL_CONTEXT_FIRST_SCREEN',
            not_stage10_i2_ranking=True,
            output_sha256={'all_family_context_screen.tsv':sha(self.ctxfile)},
            input_sha256={'partial_inventory':'a'*64,'partial_plan_provenance':'b'*64},
            all_families=5, context_supported_dark_families=2,
        )))
        self.prepprov = self.prep / 'stage182_preparation_provenance.json'
        self.prepprov.write_text(json.dumps(dict(
            source_inventory_sha256='a'*64,
            source_plan_provenance_sha256='b'*64,
            all_family_status_sha256=sha(self.statusfile),
            input_structure_sha256=sha(self.modelfile), search_ready_count=3,
        )))
        (self.ev / 'stage182_merge_provenance.json').write_text(json.dumps(dict(
            stage='STAGE09_PARTIAL_EVIDENCE_MERGED',
            output_sha256=sha(self.evfile),
            source_preparation_manifest_sha256=sha(self.prepprov),
            all_family_count=5, externally_evaluated_family_count=3,
        )))

    def execute(self, name='out', mode='check'):
        return subprocess.run([
            sys.executable, '-B', str(SCRIPT),
            '--context-dir', str(self.ctx), '--stage09-dir', str(self.ev),
            '--prepared-dir', str(self.prep), '--out-dir', str(self.root / name),
            '--mode', mode,
        ], capture_output=True, text=True)

    def test_all_families_and_statuses(self):
        result = self.execute(mode='check')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.root / 'out').exists())
        self.assertIn('FILES_WRITTEN=NO', result.stdout)
        result = self.execute(mode='run')
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = read(self.root/'out/all_family_stage10_review.tsv')
        self.assertEqual(len(rows), 5)
        self.assertEqual([r['review_tier'] for r in rows], [
            'T1_CONTEXT_AND_CONVINCING_STRUCTURE','T2_CONTEXT_SUPPORTED_DARK',
            'T4_DARK_CONTEXT_UNRESOLVED','T4_DARK_CONTEXT_UNRESOLVED',
            'T5_KNOWN_EVIDENCE_REFERENCE'
        ])
        self.assertEqual(len(read(self.root/'out/dark_family_stage10_review.tsv')), 4)
        self.assertEqual(len(read(self.root/'out/context_supported_stage10_shortlist.tsv')), 2)
        self.assertEqual(rows[1]['stage09_interpretation'],'NOT_EVALUATED_NO_STRUCTURE_NOT_NEGATIVE')
        self.assertEqual(rows[2]['stage09_interpretation'],'WEAK_STRUCTURAL_SIMILARITY_NOT_FUNCTION_PROOF')
        self.assertEqual(rows[3]['stage09_interpretation'],'SEARCH_PERFORMED_NO_RAW_HIT_AT_E1')
        self.assertTrue(all(not r['historical_I2_score'] for r in rows))
        self.assertTrue(all(r['biochemical_activity_confirmed']=='NOT_TESTED' for r in rows))
        print('STAGE183_ALL_FAMILIES_MIXED_EVIDENCE=PASS')

    def test_reject_context_modified(self):
        self.ctxfile.write_text(self.ctxfile.read_text() + 'tamper\n')
        result = self.execute()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('STAGE180_CONTEXT_SHA_MISMATCH',result.stderr)
        print('STAGE183_CONTEXT_PROVENANCE_TAMPER_REJECTED=PASS')

    def test_reject_stage09_modified(self):
        self.ev_rows[0]['stage09_external_status'] = 'NO_RAW_HIT_E1'
        write(self.evfile, self.ev_rows)
        result = self.execute()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('STAGE182_EVIDENCE_SHA_MISMATCH',result.stderr)
        print('STAGE183_STAGE09_PROVENANCE_TAMPER_REJECTED=PASS')

    def test_no_false_negative_pending(self):
        self.ev_rows[1]['stage09_external_status'] = 'NO_RAW_HIT_E1'
        write(self.evfile, self.ev_rows)
        self.update_provenance()
        result = self.execute()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('PENDING_STRUCTURAL_EVIDENCE_NOT_BLANK',result.stderr)
        print('STAGE183_PENDING_IS_NOT_NO_HIT=PASS')

    def test_reject_cross_cohort_and_overwrite(self):
        self.ev_rows[0]['family_rep'] = 'OTHER_REP'
        write(self.evfile, self.ev_rows)
        self.update_provenance()
        result = self.execute()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('CROSS_STAGE_FIELD_MISMATCH',result.stderr)
        write(self.evfile, self.ev_rows[:0]+[dict(r, family_rep=f'V{i}|rep') for i,r in enumerate(self.ev_rows,1)])
        self.update_provenance()
        self.assertEqual(self.execute(mode='run').returncode, 0)
        result = self.execute(mode='run')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('REFUSE_EXISTING_OUTPUT',result.stderr)
        print('STAGE183_CROSS_COHORT_AND_OVERWRITE_REJECTED=PASS')

    def test_new_cohort_no_historical_structure(self):
        # Use only two dark families, with no known PDBs and no Foldseek searches.
        self.ctx_rows = self.ctx_rows[1:3]
        self.ev_rows = self.ev_rows[1:3]
        self.status_rows = self.status_rows[1:3]
        for row in self.ev_rows:
            row.update(structure_state='PREDICTION_PENDING',external_search_state='NOT_EVALUATED_NO_STRUCTURE',stage09_evidence_state='NOT_EVALUATED_NO_STRUCTURE',stage09_external_status='',stage09_raw_sources='',stage09_convincing_sources='')
        for row in self.status_rows:
            row.update(structure_state='PREDICTION_PENDING', external_search_state='NOT_EVALUATED_NO_STRUCTURE')
        write(self.ctxfile, self.ctx_rows)
        write(self.evfile, self.ev_rows)
        write(self.statusfile, self.status_rows)
        self.update_provenance()
        c=json.loads((self.ctx / 'context_first_provenance.json').read_text())
        c.update(all_families=2,context_supported_dark_families=1)
        (self.ctx / 'context_first_provenance.json').write_text(json.dumps(c))
        p=json.loads(self.prepprov.read_text());p['search_ready_count']=0
        self.prepprov.write_text(json.dumps(p))
        e=json.loads((self.ev/'stage182_merge_provenance.json').read_text())
        e.update(all_family_count=2,externally_evaluated_family_count=0,source_preparation_manifest_sha256=sha(self.prepprov))
        (self.ev/'stage182_merge_provenance.json').write_text(json.dumps(e))
        result = self.execute(mode='run')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(len(read(self.root/'out/all_family_stage10_review.tsv')),2)
        print('STAGE183_NEW_COHORT_NO_PILOT_NO_STRUCTURE=PASS')


if __name__ == '__main__':
    unittest.main()
