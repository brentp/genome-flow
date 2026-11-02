"""Check the actual WDL resource expressions and low-memory guard."""
from pathlib import Path
import subprocess
import tempfile
import unittest

try:
    import WDL
except ImportError:
    WDL = None

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipIf(WDL is None, 'requires miniwdl')
class GLnexusMemoryTests(unittest.TestCase):
    def test_memory_budget_and_guard(self):
        task = WDL.load(str(ROOT/'wdl/tasks/glnexus.wdl')).tasks[0]
        class LocalStdLib(WDL.StdLib.Base):
            def _virtualize_filename(self, name): return name
        with tempfile.TemporaryDirectory() as tmp:
            stdlib = LocalStdLib('1.0',write_dir=tmp)
            for memory, budget in ((64,48),(80,60),(3,2),(1,0)):
                with self.subTest(memory=memory):
                    env = WDL.Env.Bindings()
                    inputs = dict(memory_gb=memory, cpu=1, model_type='WGS', output_prefix='cohort', fasta='unused.fa', fasta_index='unused.fa.fai', gvcfs=['unused.g.vcf.gz'], gvcf_indexes=['unused.g.vcf.gz.tbi'], docker='unused', disk_gb=1)
                    for decl in task.inputs:
                        value = WDL.Value.from_json(decl.type,inputs[decl.name]) if decl.name in inputs else decl.expr.eval(env,stdlib)
                        env = env.bind(decl.name,value)
                    for decl in task.postinputs:
                        env = env.bind(decl.name,decl.expr.eval(env,stdlib))
                    self.assertEqual(env.resolve('memory_budget_gb').value,budget)
                    if memory == 1:
                        result = subprocess.run(['bash','-c',task.command.eval(env,stdlib).value],cwd=tmp,capture_output=True,timeout=5)
                        self.assertNotEqual(result.returncode,0)
                        self.assertIn(b'at least 2',result.stderr)
                        self.assertFalse((Path(tmp)/'cohort.bcf').exists())
