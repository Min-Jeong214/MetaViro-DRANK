"""Regression: fresh unified --check must request through 08, not through 02."""
import ast
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1]/'run_pipeline_full.py'
code = SOURCE.read_text(encoding='utf-8')
ast.parse(code)
assert '"--through", "08"' in code
# The upstream command must acquire --through 08 before the check branch.
ix=code.index('upstream_command.extend([\n            "--through", "08"')
check=code.index('if args.mode == "check":',ix)
assert ix < check
print('STEP173_FULL_CONFIG_STAGE08_PREFLIGHT_REGRESSION=PASS')

# Dynamic check against a new synthetic sample: missing Stage 05 DB config
# must be rejected by the unified check (historically this incorrectly passed).
import subprocess,sys,tempfile
with tempfile.TemporaryDirectory(prefix='metaviro_full_stage08_preflight_') as temp:
    root=Path(temp)
    fixture=root/'demo'
    gen=subprocess.run([sys.executable, str(SOURCE.parent/'examples/generate_demo_input.py'), '--out-dir',str(fixture)],cwd=SOURCE.parent,capture_output=True,text=True)
    assert gen.returncode==0,(gen.stdout,gen.stderr)
    db=root/'db';db.mkdir()
    proc=subprocess.run([sys.executable,str(SOURCE), '--config',str(fixture/'config.yaml'), '--db-root',str(db), '--structure-mode','predict', '--out-dir',str(root/'full'), '--mode','check'],cwd=SOURCE.parent,capture_output=True,text=True)
    assert proc.returncode!=0,'Missing Stage 05 DBs must fail preflight'
    assert 'Missing database configuration' in (proc.stderr + proc.stdout), (proc.stdout,proc.stderr)
print('STEP173_FULL_CONFIG_MISSING_DB_REJECTED=PASS')
