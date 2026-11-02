"""Check karyotype and haploid-contig handling in the input generator."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class HaploidDefaultTests(unittest.TestCase):
    def generate(self, samples, *extra, check=True):
        with tempfile.TemporaryDirectory() as tmp:
            sheet = Path(tmp)/'samples.json'
            sheet.write_text(json.dumps(samples))
            # Inference must not require a local copy of a remote FAI.
            cmd = [sys.executable,str(ROOT/'scripts/sheet_to_inputs.py'),str(sheet),
                   '--fasta','s3://reference/genome.fa',
                   '--fasta-index','s3://reference/custom.fai',*extra]
            return subprocess.run(cmd,check=check,capture_output=True,text=True,timeout=5)

    def test_input_generator_preserves_explicit_overrides(self):
        sample = dict(sample_id='sample1',fastq_r1=['r1.fq.gz'],fastq_r2=['r2.fq.gz'])
        for override in (None,'','chrM','X,Y'):
            with self.subTest(override=override):
                extra = [] if override is None else ['--haploid-contigs',override]
                inputs = json.loads(self.generate([sample],*extra).stdout)
                self.assertEqual(inputs['GenomeFlow.fasta_index'],'s3://reference/custom.fai')
                self.assertFalse(inputs['GenomeFlow.detailed_alignment_stats'])
                detailed = json.loads(self.generate([sample],*extra,'--detailed-alignment-stats').stdout)
                self.assertTrue(detailed['GenomeFlow.detailed_alignment_stats'])
                if override is None:
                    self.assertNotIn('GenomeFlow.haploid_contigs',inputs)
                else:
                    self.assertEqual(inputs['GenomeFlow.haploid_contigs'],override)

    def test_input_generator_karyotypes(self):
        base = dict(fastq_r1=['r1.fq.gz'],fastq_r2=['r2.fq.gz'])
        samples = [dict(base,sample_id='a',karyotype='xy'),dict(base,sample_id='b',karyotype='XX'),dict(base,sample_id='c')]
        inputs = json.loads(self.generate(samples,'--par-bed','par.bed').stdout)
        self.assertEqual([s.get('karyotype') for s in inputs['GenomeFlow.samples']],['XY','XX',None])
        # XY samples default to the bundled GRCh38 PAR BED; explicit choices win.
        bundled = str(ROOT/'references/GRCh38.par.bed')
        self.assertEqual(json.loads(self.generate(samples).stdout)['GenomeFlow.par_bed'],bundled)
        self.assertNotIn('GenomeFlow.par_bed',json.loads(self.generate(samples,'--haploid-contigs','').stdout))
        self.assertNotIn('GenomeFlow.par_bed',json.loads(self.generate(samples[1:]).stdout))
        self.assertEqual(inputs['GenomeFlow.par_bed'],'par.bed')
        bad = self.generate([dict(base,sample_id='a',karyotype='XXY')],check=False)
        self.assertNotEqual(bad.returncode,0)
        self.assertIn('karyotype',bad.stderr)


if __name__ == '__main__':
    unittest.main()
