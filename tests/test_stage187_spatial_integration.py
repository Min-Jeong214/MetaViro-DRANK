#!/usr/bin/env python3
"""Generic, synthetic Stage187 checks: zero is evaluated, missing is untested."""
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/13_integrate_exploratory_spatial_evidence_portable.py'


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def save(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


class TestStage187(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='stage187_')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.base, self.spatial, self.out = (
            self.root / 'stage186', self.root / 'spatial', self.root / 'stage187')
        self.base.mkdir()
        self.spatial.mkdir()
        self.stage12 = self.root / 'conservation_provenance.json'
        self.stage12.write_text(json.dumps({'algorithm': 'ORIGINAL_STAGE12_FUNCTIONS'}))
        catalog = self.base / 'all_family_stage13_structure_integrated.tsv'
        self.catalog_rows = [
            dict(structure_id=sid, stage12_structure_residue_state=state,
                 stage12_spatial_state='NOT_EVALUATED',
                 stage12_pdb_sha256='1' * 64 if sid == 'X1' else '',
                 stage12_residue_table_sha256='2' * 64 if sid == 'X1' else '',
                 biochemical_activity_confirmed='NOT_TESTED',
                 historical_I2_score='')
            for sid, state in (
                ('X1', 'STRUCTURE_SUPPORTED_RESIDUES_AVAILABLE'),
                ('X2', 'NOT_EVALUATED'),
            )
        ]
        save(catalog, self.catalog_rows)
        (self.base / 'stage186_structure_integration_provenance.json').write_text(json.dumps({
            'stage': 'STAGE13_STRUCTURE_RESIDUE_EVIDENCE_INTEGRATION',
            'all_families': 2,
            'spatial_3d_cluster_analysis_completed': False,
            'final_biochemical_selection_completed': False,
            'output_sha256': {catalog.name: sha(catalog)},
            'source_sha256': {str(self.stage12): sha(self.stage12)},
        }))
        self.summary_file = self.spatial / 'family_spatial_cluster_summary.tsv'
        self.cluster_file = self.spatial / 'charged_spatial_clusters.tsv'
        self.summaries = [dict(structure_id='X1', length='146',
                               eligible_charged_residues='6', cluster_count='0',
                               status='NO_CLUSTER_DETECTED')]
        self.clusters = []
        self.refresh_spatial()

    def refresh_spatial(self):
        save(self.summary_file, self.summaries)
        save(self.cluster_file, self.clusters, fields=(
            'structure_id', 'cluster_rank', 'residues', 'residue_count',
            'acidic_count', 'basic_count', 'mean_plddt',
            'max_edge_ca_distance_angstrom'))
        (self.spatial / 'spatial_cluster_provenance.json').write_text(json.dumps({
            'analysis_mode': 'NEW_EXPLORATORY_L3',
            'historical_stage_l3_reproduction': False,
            'geometry': 'CA_DISTANCE_CONNECTED_COMPONENTS',
            'residue_filter': 'STAGE12_PRIORITY_RESIDUE_AND_DEKRH',
            'biochemical_activity_inferred': False,
            'ca_distance_cutoff_angstrom': 8.0,
            'min_cluster_residues': 3,
            'stage12_provenance_sha256': sha(self.stage12),
            'family_count': 1,
            'cluster_count': len(self.clusters),
            'cluster_table_sha256': sha(self.cluster_file),
            'summary_table_sha256': sha(self.summary_file),
            'inputs': [dict(structure_id='X1',
                            residue_table_sha256='2' * 64,
                            msa_pdb_sha256='1' * 64,
                            prediction_provenance_sha256='3' * 64)],
        }))

    def run_cli(self, mode='check'):
        return subprocess.run(
            [sys.executable, '-B', str(SCRIPT),
             '--stage186-dir', str(self.base), '--spatial-dir', str(self.spatial),
             '--out-dir', str(self.out), '--mode', mode],
            text=True, capture_output=True,
        )

    def test_zero_clusters_is_evaluated_and_other_family_not_evaluated(self):
        check = self.run_cli()
        self.assertEqual(check.returncode, 0, check.stderr)
        self.assertFalse(self.out.exists())
        self.assertIn('NO_CLUSTER_DETECTED_FAMILIES=1', check.stdout)
        go = self.run_cli('run')
        self.assertEqual(go.returncode, 0, go.stderr)
        with (self.out / 'all_family_stage13_spatial_integrated.tsv').open() as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['stage12_spatial_state'], 'EXPLORATORY_NO_CLUSTER_DETECTED')
        self.assertEqual(rows[0]['stage12_spatial_cluster_count'], '0')
        self.assertEqual(rows[0]['stage12_spatial_eligible_charged_residues'], '6')
        self.assertEqual(rows[1]['stage12_spatial_state'], 'NOT_EVALUATED')
        self.assertEqual(rows[1]['stage12_spatial_cluster_count'], '')
        self.assertEqual(rows[0]['biochemical_activity_confirmed'], 'NOT_TESTED')
        print('STAGE187_ZERO_IS_EVALUATED=PASS')

    def test_detected_cluster_preserved(self):
        self.summaries[0]['cluster_count'] = '1'
        self.summaries[0]['status'] = 'CLUSTERS_DETECTED'
        self.clusters = [dict(structure_id='X1', cluster_rank='1', residues='D1,E2,K3',
                              residue_count='3', acidic_count='2', basic_count='1',
                              mean_plddt='85.00', max_edge_ca_distance_angstrom='7.500')]
        self.refresh_spatial()
        go = self.run_cli('run')
        self.assertEqual(go.returncode, 0, go.stderr)
        with (self.out / 'all_family_stage13_spatial_integrated.tsv').open() as fh:
            row = list(csv.DictReader(fh, delimiter='\t'))[0]
        self.assertEqual(row['stage12_spatial_state'], 'EXPLORATORY_CLUSTERS_DETECTED')
        self.assertEqual(row['stage12_spatial_cluster_count'], '1')

    def test_reject_altered_spatial_summary(self):
        self.summary_file.write_text(self.summary_file.read_text() + 'BAD\n')
        self.assertNotEqual(self.run_cli('run').returncode, 0)
        self.assertFalse(self.out.exists())

    def test_reject_wrong_stage12_provenance(self):
        meta = self.spatial / 'spatial_cluster_provenance.json'
        d = json.loads(meta.read_text())
        d['stage12_provenance_sha256'] = '0' * 64
        meta.write_text(json.dumps(d))
        self.assertNotEqual(self.run_cli().returncode, 0)

    def test_reject_inconsistent_zero_status(self):
        self.summaries[0]['status'] = 'CLUSTERS_DETECTED'
        self.refresh_spatial()
        self.assertNotEqual(self.run_cli().returncode, 0)

    def test_reject_existing_output(self):
        self.out.mkdir()
        self.assertNotEqual(self.run_cli('run').returncode, 0)

    def test_reject_nonexistent_family(self):
        self.summaries[0]['structure_id'] = 'Z999'
        self.refresh_spatial()
        self.assertNotEqual(self.run_cli().returncode, 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
