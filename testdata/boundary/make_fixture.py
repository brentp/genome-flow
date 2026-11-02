#!/usr/bin/env python3
"""Deterministic ~70x synthetic paired reads with SNP/indels near a 4 kb boundary."""
import gzip
import json
from pathlib import Path
import random


def main():
    out = Path(__file__).resolve().parent / 'runs/fixture'
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(73)
    ref = ''.join(rng.choice('ACGT') for _ in range(8000))
    # SNP at 3500, insertion at the 4000 split, deletion at 4500 (1-based anchors).
    alt = ref[:3499] + next(b for b in 'ACGT' if b != ref[3499]) + ref[3500:4000] + 'TTGACC' + ref[4000:4500] + ref[4505:]
    # The second sample has a deletion whose REF span crosses the split.
    alt2 = ref[:3499] + next(b for b in 'ACGT' if b != ref[3499]) + ref[3500:3997] + ref[4003:4500] + 'TTGACC' + ref[4500:]
    (out/'ref.fa').write_text('>chrTiny\n' + ref + '\n')
    samples = []
    for sid, haps in [('boundary1', (ref, alt)), ('boundary2', (ref, alt2))]:
        records = [[], []]
        for h, seq in enumerate(haps):
            for start in range(h*4, len(seq)-350, 8):
                fragment = seq[start:start+350]
                pair = (fragment[:150], fragment[-150:].translate(str.maketrans('ACGT','TGCA'))[::-1])
                for mate, bases in enumerate(pair):
                    records[mate].append(f'@{sid}.{h}.{start}\n{bases}\n+\n'+ 'I'*150 +'\n')
        paths = []
        for mate, rows in enumerate(records, 1):
            p = out/f'{sid}.R{mate}.fq.gz'
            with gzip.open(p, 'wt') as fh:
                fh.write(''.join(rows))
            paths.append(str(p))
        samples.append(dict(sample_id=sid, library='lib1', platform='ILLUMINA', fastq_r1=paths[:1], fastq_r2=paths[1:]))
    header = '##fileformat=VCFv4.2\n##contig=<ID=chrTiny,length=8000>\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    (out/'sites.vcf').write_text(header+''.join(f'chrTiny\t{p}\t.\t{ref[p-1]}\t{next(b for b in "ACGT" if b != ref[p-1])}\t.\tPASS\tAF=0.5\n' for p in range(50,7950,50)))
    inputs = {'GenomeFlow.samples':samples,'GenomeFlow.fasta':str(out/'ref.fa'), 'GenomeFlow.fasta_index':str(out/'ref.fa.fai'), 'GenomeFlow.bwa_mem3_indexes':[str(out/'ref.fa')+s for s in ('.amb','.ann','.pac','.bwt.2bit.64')], 'GenomeFlow.sequencing_duplicate_detection':False,'GenomeFlow.somalier_sites':str(out/'sites.vcf'), 'GenomeFlow.align_cpu':2,'GenomeFlow.align_memory_gb':4,'GenomeFlow.dv_cpu':2,'GenomeFlow.dv_memory_gb':4,'GenomeFlow.glnexus_cpu':2,'GenomeFlow.glnexus_memory_gb':4}
    for n in (1,2):
        (out/f'inputs.{n}.json').write_text(json.dumps(dict(inputs, **{'GenomeFlow.scatter_count':n, 'GenomeFlow.detailed_alignment_stats':n == 2}), indent=2)+'\n')
    print(out)


if __name__ == '__main__':
    main()
