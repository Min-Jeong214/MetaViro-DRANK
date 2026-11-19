"""Fresh clone checks using varying synthetic cohort sizes and invalid inputs."""
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class FreshCloneTests(unittest.TestCase):
    def test_synthetic_multi_cohort_stage02_and_negative_preflight(self):
        result = subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts/validate_fresh_clone.py")],
            cwd=ROOT, text=True, capture_output=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("FRESH_CLONE_STAGE00_02_1_COHORTS=PASS", result.stdout)
        self.assertIn("FRESH_CLONE_STAGE00_02_3_COHORTS=PASS", result.stdout)
        self.assertIn("FRESH_CLONE_INVALID_GFF_REJECTED=PASS", result.stdout)
        self.assertIn("RELEASE_VALIDATION=PASS", result.stdout)


if __name__ == "__main__":
    unittest.main()
