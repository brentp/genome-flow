import gzip
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stream_align
import sheet_to_inputs


class ValidationTests(unittest.TestCase):
    def test_json_rejects_mismatched_lanes_strings_empty_and_duplicates(self):
        good = dict(sample_id='tiny', library='lib1', platform='ILLUMINA', fastq_r1=['a'], fastq_r2=['b'])
        for data in [[], [dict(good, fastq_r2=[])], [dict(good, fastq_r1='a')], [good, good]]:
            with self.subTest(data=data), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / 'sheet.json'
                path.write_text(json.dumps(data))
                with self.assertRaises(SystemExit):
                    sheet_to_inputs.validate_samples(sheet_to_inputs.load_json_sheet(path))

    def test_pipeline_failure_stops_blocked_process(self):
        p = stream_align.Pipeline()
        blocked = p.start([sys.executable, '-c', 'import time; time.sleep(60)'])
        p.start([sys.executable, '-c', 'raise SystemExit(7)'])
        try:
            with self.assertRaisesRegex(RuntimeError, 'status 7'):
                p.wait()
        finally:
            p.close()
        self.assertIsNotNone(blocked.poll())

    def test_successful_early_consumer_does_not_hang(self):
        p = stream_align.Pipeline()
        p.start([sys.executable, '-c', 'import time; time.sleep(60)'])
        p.start([sys.executable, '-c', 'pass'])
        try:
            with self.assertRaisesRegex(RuntimeError, 'downstream exit'):
                p.wait()
        finally:
            p.close()

    def test_cancellation_terminates_children(self):
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / 'supervise.py'
            pid_file = Path(tmp) / 'child.pid'
            script.write_text(
                'import sys, signal\n'
                f'sys.path.insert(0, {str(ROOT / "scripts")!r})\n'
                'from stream_align import Pipeline\n'
                'p = Pipeline()\n'
                'def stop(*args): raise RuntimeError("cancelled")\n'
                'signal.signal(signal.SIGTERM, stop)\n'
                'try:\n'
                '    child = p.start([sys.executable, "-c", "import time; time.sleep(60)"])\n'
                f'    open({str(pid_file)!r}, "w").write(str(child.pid))\n'
                '    p.wait()\n'
                'finally: p.close()\n')
            proc = subprocess.Popen([sys.executable, str(script)], stderr=subprocess.DEVNULL)
            try:
                deadline = time.monotonic() + 5
                while not pid_file.exists() and time.monotonic() < deadline:
                    time.sleep(.02)
                self.assertTrue(pid_file.exists())
                child_pid = int(pid_file.read_text())
                proc.terminate()
                proc.wait(timeout=5)
                with self.assertRaises(ProcessLookupError):
                    os.kill(child_pid, 0)
            finally:
                if proc.poll() is None:
                    proc.kill()
                proc.wait()


@unittest.skipUnless(os.environ.get('GENOME_FLOW_INTEGRATION') == '1', 'set GENOME_FLOW_INTEGRATION=1 inside alignment image')
class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='genome flow ')
        cls.base = Path(cls.tmp.name)
        spec = importlib.util.spec_from_file_location('tiny', ROOT / 'testdata/tiny/make_tiny.py')
        tiny = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(tiny)
        tiny.__file__ = str(cls.base / 'make_tiny.py')
        tiny.main()
        cls.ref = cls.base / 'ref.fa'
        cls.cmd(['samtools', 'faidx', str(cls.ref)])
        cls.cmd(['bwa-mem3', 'index', str(cls.ref)])

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def cmd(args, **kw):
        return subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90, **kw)

    def config(self):
        return dict(sample_id='tiny', library='lib1', platform='ILLUMINA',
                    fastq_r1=[str(self.base / 'illumina.R1.fq.gz')],
                    fastq_r2=[str(self.base / 'illumina.R2.fq.gz')], bam=None,
                    fasta=str(self.ref), fasta_index=str(self.ref)+'.fai',
                    cpu=2, sequencing_duplicate_detection=True, detailed_stats=False,
                    bwa_mem3_indexes=[str(self.ref)+s for s in ('.amb','.ann','.pac','.bwt.2bit.64')])

    def align(self, config, name, env=None, success=True):
        work = self.base / name
        work.mkdir()
        manifest = work / 'inputs.json'
        manifest.write_text(json.dumps(config))
        r = subprocess.run([sys.executable, str(ROOT / 'scripts/stream_align.py'), str(manifest)], cwd=work, env=env, capture_output=True, timeout=60)
        (work/'align.stderr.txt').write_bytes(r.stderr)
        if success:
            self.assertEqual(r.returncode, 0, r.stderr.decode())
            self.cmd(['samtools', 'quickcheck', str(work/'tiny.cram')])
            self.assertTrue((work/'tiny.cram.crai').is_file())
            self.assertTrue(self.cmd(['samtools','view','-T',str(self.ref),str(work/'tiny.cram'),'chrTiny1:1-1000']).stdout)
            # In-pipe QC must match a rescan of the written CRAM.
            rescan = self.cmd(['samtools','flagstat','-O','json','--input-fmt-option','reference='+str(self.ref),str(work/'tiny.cram')]).stdout
            self.assertEqual((work/'tiny.flagstat.json').read_bytes(), rescan)
            stats = work/'tiny.samtools.stats.txt'
            self.assertEqual(stats.exists(), bool(config.get('detailed_stats')))
            if stats.exists():
                rescan = self.cmd(['samtools','stats','-r',str(self.ref),'--reference',str(self.ref),str(work/'tiny.cram')]).stdout
                sn = lambda text: [l for l in text.splitlines() if l.startswith(b'SN')]
                self.assertEqual(sn(stats.read_bytes()), sn(rescan))
        else:
            self.assertNotEqual(r.returncode, 0)
        self.assertFalse(list(work.glob('*.fastq*')))
        self.assertFalse(list(work.glob('stream-*')))
        return work

    def records(self, work):
        return sorted(self.cmd(['samtools','view','-T',str(self.ref),str(work/'tiny.cram')]).stdout.splitlines())

    def test_single_multilane_and_bam_cram_parity(self):
        c = self.config()
        # Split on whole FASTQ records, preserving input order.
        for mate in (1, 2):
            lines = gzip.decompress(Path(c[f'fastq_r{mate}'][0]).read_bytes()).splitlines(keepends=True)
            split = (len(lines)//8)*4
            lanes = []
            for i, part in enumerate((lines[:split], lines[split:])):
                path = self.base / f'lane.{mate}.{i}.fq.gz'
                path.write_bytes(gzip.compress(b''.join(part)))
                lanes.append(str(path))
            c[f'fastq_r{mate}'] = lanes
        multi = self.align(c, 'multi')
        plain = self.base/'plain.fq'
        plain.write_bytes(gzip.decompress(Path(c['fastq_r1'][1]).read_bytes()))
        failed = self.align(dict(c, fastq_r1=[c['fastq_r1'][0], str(plain)]), 'mixed', success=False)
        self.assertIn(b'all gzip or all uncompressed', (failed/'align.stderr.txt').read_bytes())
        single = self.align(self.config(), 'single')
        self.assertGreaterEqual(len(self.records(single)), 8560)
        self.assertEqual(self.records(multi), self.records(single))
        self.assertEqual((multi/'tiny.dupblaster.duplicate-metrics.tsv').read_bytes(), (single/'tiny.dupblaster.duplicate-metrics.tsv').read_bytes())
        # dupblaster relies on bwa-mem3's own mate tags.
        self.assertTrue(all(b'\tMC:Z:' in r and b'\tMQ:i:' in r for r in self.records(single)))
        for ext in ('bam', 'cram'):
            source = single/'tiny.cram'
            if ext == 'bam':
                source = self.base / 'input.bam'
                self.cmd(['samtools','view','-T',str(self.ref),'-b','-o',str(source),str(single/'tiny.cram')])
            # Compare with a materialized baseline using the identical collation order.
            bam = self.cmd(['samtools','collate','-Ou','-f','-r','1000000','--reference',str(self.ref),str(source)])
            r1, r2 = [self.base/f'{ext}.R{m}.fq.gz' for m in (1,2)]
            self.cmd(['samtools','fastq','-n','-1',str(r1),'-2',str(r2),'-0','/dev/null','-s','/dev/null','-'], input=bam.stdout)
            baseline = self.align(dict(c,fastq_r1=[str(r1)],fastq_r2=[str(r2)]),ext+'-baseline')
            streamed = self.align(dict(c,fastq_r1=[],fastq_r2=[],bam=str(source),detailed_stats=True),ext+'-stream')
            self.assertEqual(self.records(streamed), self.records(baseline))
            self.assertEqual((streamed/'tiny.dupblaster.duplicate-metrics.tsv').read_bytes(), (baseline/'tiny.dupblaster.duplicate-metrics.tsv').read_bytes())

    def test_query_grouped_input_skips_collate_and_false_claims_fail(self):
        c = self.config()
        ubam = self.base/'unaligned.bam'
        self.cmd(['samtools','import','-1',c['fastq_r1'][0],'-2',c['fastq_r2'][0],'-o',str(ubam)])
        self.assertIn(b'GO:query', self.cmd(['samtools','view','-H',str(ubam)]).stdout)
        baseline = self.align(c, 'ubam-baseline')
        streamed = self.align(dict(c,fastq_r1=[],fastq_r2=[],bam=str(ubam)), 'ubam')
        self.assertEqual(self.records(streamed), self.records(baseline))
        # A coordinate-sorted file falsely claiming query grouping yields
        # unpaired reads; that must fail rather than produce a near-empty CRAM.
        sorted_bam = self.base/'sorted.bam'
        self.cmd(['samtools','view','-b','-T',str(self.ref),'-o',str(sorted_bam),str(baseline/'tiny.cram')])
        header = self.cmd(['samtools','view','-H',str(sorted_bam)]).stdout.replace(b'SO:coordinate',b'SO:queryname')
        false_claim = self.base/'false-claim.bam'
        with open(false_claim,'wb') as out:
            subprocess.run(['samtools','reheader','-',str(sorted_bam)],input=header,stdout=out,check=True,timeout=30)
        failed = self.align(dict(c,fastq_r1=[],fastq_r2=[],bam=str(false_claim)), 'false-claim', success=False)
        self.assertIn(b'unpaired', (failed/'align.stderr.txt').read_bytes())

    @unittest.skipUnless(os.environ.get('GENOME_FLOW_SHM_TEST') == '1', 'requires disposable private IPC namespace')
    def test_shared_memory_reused_across_alignment_processes(self):
        # run_pipe.sh opts in inside a fresh Docker IPC namespace. Refuse to
        # touch a namespace containing staged indexes owned by someone else.
        self.assertFalse(self.cmd(['bwa-mem3','shm','-l']).stdout.strip())
        c = self.config()
        disk = self.align(c, 'shm-disk')
        self.cmd(['bwa-mem3','shm',str(self.ref)])
        try:
            for i in range(2):
                shared = self.align(c, f'shm-shared-{i}')
                self.assertIn('attached from shm', (shared/'align.stderr.txt').read_text())
                self.assertEqual(self.records(disk), self.records(shared))
                self.assertEqual((disk/'tiny.dupblaster.duplicate-metrics.tsv').read_bytes(), (shared/'tiny.dupblaster.duplicate-metrics.tsv').read_bytes())
                self.assertTrue(self.cmd(['bwa-mem3','shm','-l']).stdout.strip())
        finally:
            self.cmd(['bwa-mem3','shm','-d'])

    def test_task_scratch_does_not_depend_on_system_tmp(self):
        blocked_tmp = self.base/'not-a-directory'
        blocked_tmp.write_text('System temporary storage is unavailable in this test.')
        c = self.config()
        baseline = self.align(c, 'scratch-baseline')
        bam = self.base/'scratch.bam'
        self.cmd(['samtools','view','-T',str(self.ref),'-b','-o',str(bam),str(baseline/'tiny.cram')])
        self.align(dict(c,fastq_r1=[],fastq_r2=[],bam=str(bam)), 'scratch-stream',
                   env=dict(os.environ,TMPDIR=str(blocked_tmp)))

    def test_corrupt_lanes_fail(self):
        c = self.config()
        data = Path(c['fastq_r1'][0]).read_bytes()
        bad_crc = bytearray(data)
        bad_crc[-8] ^= 0xFF
        for name, content in (('truncated', data[:-8]), ('crc', bytes(bad_crc))):
            with self.subTest(corruption=name):
                bad = self.base/f'{name}.fq.gz'
                bad.write_bytes(content)
                self.align(dict(c, fastq_r1=[c['fastq_r1'][0], str(bad)], fastq_r2=c['fastq_r2']*2), name, success=False)

    def test_simulated_names_detection_policy(self):
        c = self.config()
        c['fastq_r1'] = [str(self.base/'simulated.R1.fq.gz')]
        c['fastq_r2'] = [str(self.base/'simulated.R2.fq.gz')]
        self.align(c, 'sim-on', success=False)
        self.align(dict(c, sequencing_duplicate_detection=False), 'sim-off')

    def test_early_aligner_exit_stops_lane_producers(self):
        c = self.config()
        c['fastq_r1'] *= 2
        c['fastq_r2'] *= 2
        fake = self.base/'bin'
        fake.mkdir()
        tool = fake/'bwa-mem3'
        tool.write_text('#!/bin/sh\nexit 3\n')
        tool.chmod(0o755)
        self.align(c, 'early-exit', env=dict(os.environ, PATH=str(fake)+':'+os.environ['PATH']),success=False)


if __name__ == '__main__':
    unittest.main()
