#!/usr/bin/env python3
"""Self-contained tests for official CPU-first wrapper routing and refusal guards."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / 'run_metaviro.py'
spec = importlib.util.spec_from_file_location('run_metaviro', MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class MainCliTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mv_step181_')
        self.addCleanup(self.temp.cleanup)
        self.p = Path(self.temp.name)
        self.input = self.p / 'input'
        self.input.mkdir()
        self.upstream = self.p / 'workflow'
        self.screen = self.p / 'screen'
        self.dbs = self.p / 'db'
        self.dbs.mkdir()
        for name in ('swissprot','phrog','pfam','phrog_annotations'):
            (self.dbs / name).write_text('fixture\n')
        (self.dbs / 'swissprot.dbtype').write_text('1')
        for ext in ('.dbtype', '.index', '_h', '_h.dbtype', '_h.index'):
            (self.dbs / ('phrog'+ext)).write_text('1')
        for ext in ('.h3f','.h3i','.h3m','.h3p'):
            (self.dbs / ('pfam'+ext)).write_text('1')
        self.conf = self.p / 'config.yaml'
        self.conf.write_text(
            f'input_dir: {self.input}\noutput_dir: {self.upstream}\n'
            'min_aa: 150\nwindow: 5\nthreads: 2\n'
            'databases:\n'
            + ''.join(f'  {key}: {self.dbs / key}\n' for key in
                      ('swissprot','phrog','pfam','phrog_annotations')))
        self.call_log = []
        self.patch_root = patch.multiple(mod, ROOT=self.p,
                                         UPSTREAM=self.p/'run_pipeline.py',
                                         NAMESPACE=self.p/'namespace.py',
                                         SCREEN=self.p/'run_context_screen.py')
        self.patch_root.start()
        self.addCleanup(self.patch_root.stop)
        for name in ('run_pipeline.py','namespace.py','run_context_screen.py'):
            (self.p / name).write_text('# mock runner\n')

    def invoke(self, *args):
        with patch('sys.argv', ['run_metaviro.py', *args]):
            mod.main()

    def test_new_config_no_normalizer_preflight(self):
        with patch.object(mod, 'command', side_effect=lambda label,argv: self.call_log.append((label, list(map(str,argv))))):
            self.invoke('run','--config',str(self.conf),'--screen-out-dir',str(self.screen),'--mode','check')
        self.assertEqual(len(self.call_log), 1)
        self.assertIn('--through',self.call_log[0][1])
        self.assertIn('07',self.call_log[0][1])
        self.assertFalse(self.screen.exists())
        self.assertFalse(self.upstream.exists())
        print('STEP181_NEW_CONFIG_PREFLIGHT_NO_WRITE=PASS')

    def test_new_normalized_input_is_optional_and_safe(self):
        with patch.object(mod, 'command', side_effect=lambda label,argv: self.call_log.append((label,list(map(str,argv))))):
            self.invoke('run','--config',str(self.conf),'--screen-out-dir',str(self.screen),
                        '--normalize-ids','--mode','check')
        self.assertEqual(len(self.call_log),1)
        self.assertIn('namespace.py',self.call_log[0][1][2])
        self.assertFalse(self.screen.exists())
        print('STEP181_ID_NAMESPACE_PREFLIGHT_NO_WRITE=PASS')

    def test_missing_db_and_overlap_are_rejected(self):
        (self.dbs/'pfam.h3p').unlink()
        with self.assertRaisesRegex(ValueError, 'PFAM_INDEX'):
            self.invoke('run','--config',str(self.conf),'--screen-out-dir',str(self.screen))
        (self.dbs/'pfam.h3p').write_text('1')
        with self.assertRaisesRegex(ValueError,'OUTPUT_OVERLAPS_SOURCE'):
            self.invoke('run','--config',str(self.conf),'--screen-out-dir',str(self.input/'nested'))
        print('STEP181_DB_AND_PATH_REJECTION=PASS')

    def test_resume_check_and_run_record_state(self):
        self.upstream.mkdir()
        def fake_cmd(label,argv):
            self.call_log.append((label,list(map(str,argv))))
            if 'CPU context-first screening' == label:
                child=self.screen/'review/context_screen_run_provenance.json'
                child.parent.mkdir(parents=True)
                child.write_text(json.dumps({'family_coverage_verified':True,
                      'stage10_full_ranking_completed':False,
                      'source_families':2,'source_members':12,
                      'context_supported_dark_families':1}))
        with patch.object(mod,'command',side_effect=fake_cmd):
            self.invoke('resume','--workflow-output',str(self.upstream),
                        '--screen-out-dir',str(self.screen), '--mode','check')
            self.assertFalse(self.screen.exists())
            self.invoke('resume','--workflow-output',str(self.upstream),
                        '--screen-out-dir',str(self.screen), '--mode','run')
        report=json.loads((self.screen/'run_metaviro_provenance.json').read_text())
        self.assertEqual(report['family_count'],2)
        self.assertFalse(report['stage10_final_ranking_completed'])
        self.assertFalse(report['external_structure_search_executed'])
        self.assertEqual(len(self.call_log),3)
        print('STEP181_RESUME_PARTIAL_STATUS_PROVENANCE=PASS')

    def test_resume_existing_output_rejected(self):
        self.upstream.mkdir()
        self.screen.mkdir()
        with self.assertRaisesRegex(ValueError, 'OUTPUT_ALREADY_EXISTS'):
            self.invoke('resume','--workflow-output',str(self.upstream),'--screen-out-dir',str(self.screen))
        print('STEP181_RESUME_OVERWRITE_REJECTION=PASS')


if __name__=='__main__':
    unittest.main(verbosity=0)
    print('ALL_TESTS_PASS')
