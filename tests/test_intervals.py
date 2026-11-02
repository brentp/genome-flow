"""Run scripts/make_intervals.py on small references and capture BEDs."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/make_intervals.py'


def read_bed(path):
    return [(chrom, int(start), int(end)) for chrom, start, end in (line.split() for line in path.read_text().splitlines())]


class IntervalTests(unittest.TestCase):
    def scatter(self, bed, scatter_count=3, contigs=None, haploid_contigs=None, karyotypes=(), par_bed=None):
        """Return (process result, shards, calling regions, per-sample haploid contigs)."""
        with tempfile.TemporaryDirectory(prefix='interval test ') as tmp:
            base = Path(tmp)
            fai = base / 'ref.fa.fai'
            lengths = contigs if isinstance(contigs,dict) else {chrom:1000 for chrom in (contigs if contigs is not None else ['chr1','chr2','chr10'])}
            fai.write_text(''.join(f'{chrom}\t{length}\n' for chrom,length in lengths.items()))
            capture = base / 'capture.bed'
            if bed is not None:
                capture.write_text(bed)
            par = base / 'par.bed'
            if par_bed is not None:
                par.write_text(par_bed)
            manifest = base / 'inputs.json'
            manifest.write_text(json.dumps(dict(fasta_index=str(fai), capture_bed=str(capture) if bed is not None else None, scatter_count=scatter_count,
                                                haploid_override=haploid_contigs, karyotypes=list(karyotypes), par_bed=str(par) if par_bed is not None else None)))
            result = subprocess.run([sys.executable, str(SCRIPT), str(manifest)], text=True, capture_output=True, cwd=base, timeout=5)
            listed = (base/'shard_list.txt').read_text().split() if result.returncode == 0 else []
            shards = [read_bed(base/path) for path in listed]
            calling = read_bed(base/'calling.bed') if (base/'calling.bed').exists() else None
            haploid_file = base/'sample_haploid_contigs.json'
            haploid = json.loads(haploid_file.read_text()) if haploid_file.exists() else None
            return result, shards, calling, haploid

    def test_reference_order_and_exact_coverage(self):
        result, shards, calling, _ = self.scatter('chr10\t0\t150\nchr2\t0\t150\nchr1\t0\t150\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(shards, [[('chr1',0,150)], [('chr2',0,150)], [('chr10',0,150)]])
        self.assertEqual(calling, [('chr1',0,150), ('chr2',0,150), ('chr10',0,150)])

    def test_overlaps_adjacent_duplicates_and_containment_merge(self):
        result, shards, calling, _ = self.scatter('chr1\t150\t250\nchr1\t100\t200\nchr1\t100\t200\nchr1\t120\t140\nchr1\t250\t300\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(calling, [('chr1', 100, 300)])
        self.assertEqual(shards, [[('chr1',100,167)], [('chr1',167,234)], [('chr1',234,300)]])

    def test_more_shards_than_bases(self):
        result, shards, _, _ = self.scatter('chr1\t0\t2\n', scatter_count=5)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(shards, [[('chr1',0,1)], [('chr1',1,2)]])

    def test_whole_reference(self):
        result, shards, _, _ = self.scatter(None)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(shards, [[('chr1',0,1000)], [('chr2',0,1000)], [('chr10',0,1000)]])

    def test_exact_balanced_shards_preserve_coverage_and_order(self):
        cases = [
            {'chr1':1000,'chr2':1001,'chr10':1002},
            {f'chr{i+1}':length for i,length in enumerate([
                248956422,242193529,198295559,190214555,181538259,170805979,
                159345973,145138636,138394717,133797422,135086622,133275309,
                114364328,107043718,101991189,90338345,83257441,80373285,
                58617616,64444167,46709983,50818468,156040895,57227415])},
        ]
        for lengths in cases:
            with self.subTest(lengths=lengths):
                result,shards,_,_ = self.scatter(None,scatter_count=20,contigs=lengths)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertEqual(len(shards),20)
                sizes = [sum(e-s for _,s,e in shard) for shard in shards]
                self.assertLessEqual(max(sizes)-min(sizes),1)
                rank = {chrom:i for i,chrom in enumerate(lengths)}
                prior_rank = -1
                ends = {chrom:0 for chrom in lengths}
                for shard in shards:
                    for chrom,start,end in shard:
                        self.assertGreaterEqual(rank[chrom],prior_rank)
                        self.assertEqual(start,ends[chrom])
                        self.assertLess(start,end)
                        ends[chrom] = end
                        prior_rank = rank[chrom]
                self.assertEqual(ends,lengths)

    def test_haploid_contigs_from_reference_not_capture_regions(self):
        cases = [(['chr1','chrY','chrX'], 'chrX,chrY'),
                 (['chr1','chrX'], 'chrX'), (['chr1','chrY'], 'chrY'),
                 (['chr1','X','Y','chrX_random'], '')]
        for contigs, expected in cases:
            with self.subTest(contigs=contigs):
                result, _, _, haploid = self.scatter('chr1\t0\t100\n', contigs=contigs, karyotypes=['XY'], par_bed='chr1\t0\t10\n')
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertEqual(haploid,[expected])

    def test_haploid_resolution_requires_matching_par_for_xy(self):
        par = 'chrX\t10\t20\nchrY\t10\t20\n'
        ok = [(None, ['XX', ''], None, ['', '']),
              (None, ['XY'], par, ['chrX,chrY']),
              ('', ['XY'], None, ['']),
              ('chrX', ['XY', 'XX'], par, ['chrX', ''])]
        for override, karyotypes, par_bed, expected in ok:
            with self.subTest(override=override, karyotypes=karyotypes, par=par_bed):
                result, _, _, resolved = self.scatter('chr1\t0\t100\n', contigs=['chr1','chrX','chrY'], haploid_contigs=override, karyotypes=karyotypes, par_bed=par_bed)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(resolved, expected)
        # A PAR BED for another naming scheme or build must fail, not be ignored.
        for override, karyotypes, par_bed, message in [(None, ['XX', 'XY'], None, 'par_bed'),
                                                       ('chrZ', ['XY'], par, 'chrZ'),
                                                       (None, ['xy'], par, 'karyotype'),
                                                       (None, ['XY'], 'X\t10\t20\n', 'unknown reference contig X'),
                                                       (None, ['XY'], 'chrY\t10\t2000\n', 'outside contig bounds')]:
            with self.subTest(override=override, karyotypes=karyotypes, par=par_bed):
                result, shards, _, _ = self.scatter('chr1\t0\t100\n', contigs=['chr1','chrX','chrY'], haploid_contigs=override, karyotypes=karyotypes, par_bed=par_bed)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertEqual(shards, [])

    def test_invalid_inputs_fail_promptly(self):
        for bed in ('chrOther\t0\t10\n', 'chr1\t-1\t10\n', 'chr1\t0\t1001\n', 'chr1\t10\t10\n', 'chr1\tbad\t20\n', 'chr1\t10\n', ''):
            with self.subTest(bed=bed):
                result, shards, _, _ = self.scatter(bed)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(shards)
        for scatter_count in (0, -1):
            with self.subTest(scatter_count=scatter_count):
                result, _, _, _ = self.scatter(None, scatter_count)
                self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
