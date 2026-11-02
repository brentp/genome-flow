"""Run the real WDL commands against repeat indels and multiallelic records."""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

try:
    import WDL
except ImportError:
    WDL = None

ROOT = Path(__file__).resolve().parents[1]


def rows(text):
    return [line for line in text.splitlines() if not line.startswith('#')]


@unittest.skipUnless(WDL and shutil.which('bcftools') and shutil.which('samtools'), 'requires miniwdl, bcftools, samtools')
class NormalizationTests(unittest.TestCase):
    def setUp(self):
        if '--write-index' not in subprocess.run(['bcftools','norm','--help'],capture_output=True,text=True).stderr:
            self.skipTest('requires bcftools norm with --write-index')
        self.tmp = tempfile.TemporaryDirectory(prefix='norm test ')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.ref = self.base/'reference.fa'
        self.ref.write_text('>chr1\nCAAAAAAAAGCCCCCCCCCC\n')
        self.cmd(['samtools','faidx',str(self.ref)])
        self.header = '''##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##INFO=<ID=TAG,Number=1,Type=String,Description="Annotation">
##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttiny
'''
        self.records = [
            'chr1\t4\tsnp\tA\tT\t50\tPASS\tTAG=one\tGT:DP:AD:PL\t0/1:12:6,6:80,0,90',
            'chr1\t7\tdel\tAA\tA\t60\tPASS\tTAG=two\tGT:DP:AD:PL\t0/1:13:7,6:90,0,100',
            'chr1\t15\tmulti\tC\tA,T\t70\tPASS\tTAG=three\tGT:DP:AD:PL\t1/2:14:2,5,7:0,10,20,30,40,50',
        ]

    def cmd(self, args, **kwargs):
        return subprocess.run(args, check=True, capture_output=True, text=True, timeout=20, **kwargs)

    def prepare(self, mismatch=False):
        records = list(self.records)
        if mismatch:
            records[0] = records[0].replace('\tA\tT\t', '\tG\tT\t')
        source = self.base/'source.vcf'
        source.write_text(self.header+'\n'.join(records)+'\n')
        self.cmd(['bcftools','view','-Oz','-o',str(self.base/'input with space.vcf.gz'),str(source)])
        self.cmd(['bcftools','index','-t',str(self.base/'input with space.vcf.gz')])
        self.cmd(['bcftools','view','-Ob','-o',str(self.base/'input.bcf'),str(source)])
        block = self.base/'block.vcf'
        block.write_text(self.header+'chr1\t1\t.\tC\t<*>\t.\t.\tEND=20\tGT:DP\t0/0:20\n')
        self.cmd(['bcftools','view','-Oz','-o',str(self.base/'input.g.vcf.gz'),str(block)])
        self.cmd(['bcftools','index','-t',str(self.base/'input.g.vcf.gz')])

    def run_task(self, task_name, mismatch=False):
        self.prepare(mismatch)
        work = self.base/task_name
        work.mkdir()
        path = ROOT/'wdl/tasks'/('deepvariant.wdl' if task_name=='MergeSampleCalls' else 'glnexus.wdl')
        task = next(t for t in WDL.load(str(path)).tasks if t.name == task_name)
        values = dict(fasta=str(self.ref), fasta_index=str(self.ref)+'.fai', docker='unused')
        if task_name == 'MergeSampleCalls':
            values.update(sample_id='tiny',vcfs=[str(self.base/'input with space.vcf.gz')],vcf_indexes=[str(self.base/'input with space.vcf.gz.tbi')],gvcfs=[str(self.base/'input.g.vcf.gz')],gvcf_indexes=[str(self.base/'input.g.vcf.gz.tbi')])
        else:
            values.update(gvcfs=[str(self.base/'input.g.vcf.gz')],gvcf_indexes=[str(self.base/'input.g.vcf.gz.tbi')],output_prefix='tiny',model_type='WGS',cpu=1,memory_gb=64,disk_gb=1)
            stub = work/'glnexus_cli'
            stub.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > glnexus.args\ncat "$TEST_BCF"\n')
            stub.chmod(0o755)
        class LocalStdLib(WDL.StdLib.Base):
            def _virtualize_filename(self, name): return name
            def _devirtualize_filename(self, name): return name
        stdlib = LocalStdLib('1.0',write_dir=str(work))
        env = WDL.Env.Bindings()
        for decl in task.inputs:
            value = WDL.Value.from_json(decl.type,values[decl.name]) if decl.name in values else decl.expr.eval(env,stdlib)
            env = env.bind(decl.name,value)
        for decl in task.postinputs:
            env = env.bind(decl.name,decl.expr.eval(env,stdlib))
        command = task.command.eval(env,stdlib).value
        result = subprocess.run(['bash','-c',command],cwd=work,capture_output=True,text=True,timeout=30,env=dict(os.environ, PATH=str(work)+':'+os.environ['PATH'], TEST_BCF=str(self.base/'input.bcf')))
        if mismatch:
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Reference allele mismatch',result.stderr)
            self.assertFalse((work/'tiny.bcf.csi').exists())
            return
        self.assertEqual(result.returncode,0,result.stderr)
        output = work/'tiny.bcf'
        actual = rows(self.cmd(['bcftools','view','-Ov',str(output)]).stdout)
        self.assertEqual([int(r.split('\t')[1]) for r in actual],[1,4,15])
        expected = {r.split('\t')[2]:r.split('\t') for r in self.records}
        for row in actual:
            fields = row.split('\t')
            original = expected[fields[2]]
            self.assertEqual(fields[5:],original[5:])  # QUAL, FILTER, INFO, FORMAT, all sample data
            if fields[2]=='multi': self.assertEqual(fields[3:5],['C','A,T'])
            if fields[2]=='del': self.assertEqual(fields[3:5],['CA','C'])
        normalized_again = rows(self.cmd(['bcftools','norm','-f',str(self.ref),'-c','e','-Ov',str(output)]).stdout)
        self.assertEqual(actual,normalized_again)
        self.assertEqual(self.cmd(['bcftools','query','-r','chr1:1-2','-f','%ID\n',str(output)]).stdout,'del\n')
        self.assertTrue((work/'tiny.bcftools.stats.txt').stat().st_size)
        if task_name=='GLnexus':
            args = (work/'glnexus.args').read_text().splitlines()
            self.assertEqual(args[args.index('--mem-gbytes')+1],'48')
        self.assertTrue((work/'tiny.bcf.csi').is_file())
        self.assertFalse(list(work.glob('*.vcf.gz')) if task_name=='GLnexus' else list(work.glob('tiny.vcf.gz')))
        if task_name=='MergeSampleCalls':
            before = rows(self.cmd(['bcftools','view','-Ov',str(self.base/'input.g.vcf.gz')]).stdout)
            after = rows(self.cmd(['bcftools','view','-Ov',str(work/'tiny.g.bcf')]).stdout)
            self.assertEqual(before,after)

    def test_deepvariant_reuses_or_repairs_shard_indexes(self):
        self.prepare()
        for source,name in [('input with space.vcf.gz','tiny.part.vcf.gz'),('input.g.vcf.gz','tiny.part.g.vcf.gz')]:
            shutil.copyfile(self.base/source,self.base/name)
            shutil.copyfile(self.base/(source+'.tbi'),self.base/(name+'.tbi'))
        source = (ROOT/'wdl/tasks/deepvariant.wdl').read_text()
        command = source.split('    # DeepVariant normally writes',1)[1].split('>>>',1)[0]
        command = '# DeepVariant normally writes'+command
        command = command.replace('~{sample_id}','tiny').replace('~{shard_id}','part')
        indexes = [self.base/'tiny.part.vcf.gz.tbi',self.base/'tiny.part.g.vcf.gz.tbi']
        mtimes = [p.stat().st_mtime_ns for p in indexes]
        self.cmd(['bash','-euc',command],cwd=self.base)
        self.assertEqual([p.stat().st_mtime_ns for p in indexes],mtimes)
        indexes[0].write_bytes(b'')
        indexes[1].unlink()
        self.cmd(['bash','-euc',command],cwd=self.base)
        for index in indexes:
            self.assertGreater(index.stat().st_size,0)
            self.cmd(['bcftools','view','-H','-r','chr1',str(index)[:-4]])

    def test_sample_merge(self): self.run_task('MergeSampleCalls')
    def test_cohort_stream(self): self.run_task('GLnexus')
    def test_sample_ref_mismatch(self): self.run_task('MergeSampleCalls',mismatch=True)
    def test_cohort_ref_mismatch(self): self.run_task('GLnexus',mismatch=True)
