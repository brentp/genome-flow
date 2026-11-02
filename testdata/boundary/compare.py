#!/usr/bin/env python3
"""Compare calls and per-base gVCF genotype coverage on the 8 kb fixture."""
import json
from pathlib import Path
import sys
import subprocess


def records(path):
    text = subprocess.check_output(['bcftools','view','-H',str(path)],text=True)
    return [line.split('\t') for line in text.splitlines()]


def calls(path):
    return [(r[0], int(r[1]), r[3], r[4], r[6], tuple(dict(zip(r[8].split(':'),sample.split(':')))['GT'] for sample in r[9:])) for r in records(path)]


def coverage(path):
    bases = [set() for _ in range(8000)]
    seen = set()
    for r in records(path):
        start = int(r[1])-1
        info = dict(item.split('=',1) for item in r[7].split(';') if '=' in item)
        end = int(info.get('END', start+len(r[3])))
        gt = dict(zip(r[8].split(':'),r[9].split(':')))['GT']
        kind = 'reference' if 'END' in info else 'variant'
        key = (start, end, r[3], r[4], gt)
        assert key not in seen, f'duplicate gVCF record in {path}: {key}'
        seen.add(key)
        for pos in range(start, end):
            # Normalization can shift a deletion into a reference block even
            # without scattering. Compare both tracks, allowing that overlap.
            label = (kind, gt)
            assert kind != 'reference' or label not in bases[pos], f'overlapping reference blocks at {pos+1}'
            bases[pos].add(label)
    assert all(bases), f'gVCF coverage gaps in {path}'
    return bases


def compare(left, right):
    a = json.loads(Path(left).read_text())
    b = json.loads(Path(right).read_text())
    assert a['GenomeFlow.alignment_stats'] == [], 'lean QC should omit samtools stats'
    assert len(b['GenomeFlow.alignment_stats']) == 2, 'detailed QC should report both samples'
    for key in ('sample_bcfs', 'sample_gbcfs'):
        for p, q in zip(a['GenomeFlow.'+key], b['GenomeFlow.'+key]):
            if key == 'sample_bcfs':
                x, y = calls(p), calls(q)
                assert len(x) == 3, f'expected SNP and two indels in {p}: {x}'
                assert x == y, f'variant/genotype mismatch:\n{x}\n{y}'
                assert records(p) == records(q), f'VCF quality/annotation mismatch: {p}, {q}'
            else:
                assert coverage(p) == coverage(q), f'gVCF genotype coverage differs: {p}, {q}'
    assert records(a['GenomeFlow.cohort_bcf']) == records(b['GenomeFlow.cohort_bcf']), 'joint-call mismatch'
    for output in (a,b):
        for path, index in zip(output['GenomeFlow.sample_bcfs']+output['GenomeFlow.sample_gbcfs']+[output['GenomeFlow.cohort_bcf']], output['GenomeFlow.sample_bcf_indexes']+output['GenomeFlow.sample_gbcf_indexes']+[output['GenomeFlow.cohort_bcf_index']]):
            assert path.endswith('.bcf') and index.endswith('.bcf.csi')
            assert Path(index).stat().st_size
            indexed = subprocess.check_output(['bcftools','view','-H','-r','chrTiny',path+'##idx##'+index],text=True)
            assert [line.split('\t') for line in indexed.splitlines()] == records(path)
        for key in ('flagstats','sample_variant_stats'):
            assert len(output['GenomeFlow.'+key]) == 2
            assert all(Path(p).stat().st_size for p in output['GenomeFlow.'+key])
        for path in output['GenomeFlow.alignment_stats']:
            assert Path(path).stat().st_size
        assert len(output['GenomeFlow.alignment_stats']) in (0,2)
        assert len(output['GenomeFlow.coverage_reports']) == 2
        for files in output['GenomeFlow.coverage_reports']:
            assert any(p.endswith('.mosdepth.summary.txt') for p in files)
            assert all(Path(p).stat().st_size for p in files)
        for key in ('somalier_samples','somalier_pairs','somalier_report','cohort_variant_stats'):
            assert Path(output['GenomeFlow.'+key]).stat().st_size
    print('Variant/genotype calls, gVCF coverage, joint calls, and QC outputs agree.')


if __name__ == '__main__':
    compare(*sys.argv[1:])
