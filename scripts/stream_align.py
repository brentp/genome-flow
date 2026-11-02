#!/usr/bin/env python3
"""Align one sample to a duplicate-marked, coordinate-sorted CRAM in one streaming pass.

    reads -> bwa-mem3 -> dupblaster -> mako (sort) -> tee -> samtools view   (CRAM + CRAI)
                                                         -> samtools flagstat (+ optional stats)

Reads come from one FASTQ lane pair (passed to bwa-mem3 directly), several lane
pairs (concatenated into two FIFOs), or a BAM/CRAM (converted by samtools fastq).
Nothing between the tools is written to disk except mako's sort spill and, for
BAM/CRAM input that is not grouped by read name, a small samtools collate buffer.
Both live in a scratch directory in the task's working directory.

Usage: stream_align.py manifest.json  (written by wdl/tasks/align.wdl). Paths
from the manifest are passed to tools as argv, never through a shell.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import time

NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*')
# Required bwa-mem3 index files, then optional ones (.alt enables ALT-aware mapping).
BWA_REQUIRED = ('.bwt.2bit.64', '.amb', '.ann', '.pac')
BWA_OPTIONAL = ('.alt', '.hdr', '.dict')


class Pipeline:
    """Start connected processes; fail as soon as any exits nonzero; stop all on exit."""

    def __init__(self, env=None):
        self.jobs = []
        self.sinks = []
        self.env = env

    def start(self, args, stdin=None, stdout=None, stderr=None, sink=False):
        # Own process group, so close() also stops any children a tool spawns.
        job = subprocess.Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, start_new_session=True, env=self.env)
        self.jobs.append((args[0], job))
        if sink:
            self.sinks.append(job)
        return job

    def wait(self):
        # Sinks are the final readers (CRAM writer and tee'd QC). Without explicit
        # sinks the last started process is the only one.
        sinks = self.sinks or [self.jobs[-1][1]]
        finished_at = None
        while True:
            pending = False
            for name, job in self.jobs:
                code = job.poll()
                if code is None:
                    pending = True
                elif code != 0:
                    raise RuntimeError(f"{name} exited with status {code}")
            if not pending:
                return
            # Successful sink exits must not leave producers blocked in FIFO
            # open/write forever (including consumers that exit zero early).
            if all(job.poll() == 0 for job in sinks):
                finished_at = finished_at or time.monotonic()
                if time.monotonic() - finished_at > 2:
                    raise RuntimeError('upstream process remained running after downstream exit')
            time.sleep(0.05)

    def close(self):
        for _, job in self.jobs:
            if job.poll() is None:
                try:
                    os.killpg(job.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic() + 2
        while any(job.poll() is None for _, job in self.jobs) and time.monotonic() < deadline:
            time.sleep(0.05)
        for _, job in self.jobs:
            if job.poll() is None:
                try:
                    os.killpg(job.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            job.wait()


def validate(c):
    # scripts/sheet_to_inputs.py checks sample sheets fully; this guards direct WDL inputs.
    for key in ('sample_id', 'library', 'platform'):
        if not isinstance(c.get(key), str) or not NAME.fullmatch(c[key]):
            raise ValueError(f"{key}: use letters, digits, underscores, dots or hyphens")
    if len(c['fastq_r1']) != len(c['fastq_r2']) or bool(c['fastq_r1']) == bool(c.get('bam')):
        raise ValueError('provide equal nonempty FASTQ lane lists or BAM/CRAM, exclusively')


def link(source, target):
    Path(target).symlink_to(Path(source).resolve(strict=True))


def stage_reference(c):
    """Link the reference and its bwa-mem3 index files as ref.fa*, as bwa-mem3 expects."""
    link(c['fasta'], 'ref.fa')
    link(c['fasta_index'], 'ref.fa.fai')
    found = set()
    for path in c['bwa_mem3_indexes']:
        suffix = next((s for s in BWA_REQUIRED + BWA_OPTIONAL if path.endswith(s)), None)
        if suffix is None or suffix in found:
            raise ValueError(f'unrecognized or duplicate bwa-mem3 index suffix: {path}')
        found.add(suffix)
        link(path, 'ref.dict' if suffix == '.dict' else 'ref.fa' + suffix)
    missing = set(BWA_REQUIRED) - found
    if missing:
        raise ValueError(f'missing bwa-mem3 indexes: {sorted(missing)}')


def is_gzip(path):
    with open(path, 'rb') as fh:
        return fh.read(2) == b'\x1f\x8b'


def query_grouped(path):
    # An @HD grouping claim lets samtools fastq pair mates in input order with no
    # collate spill; check_pairing() catches false claims after alignment.
    r = subprocess.run(['samtools', 'view', '-H', '--no-PG', path], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f'cannot read header of {path}: {r.stderr.strip()}')
    for line in r.stdout.splitlines():
        if line.startswith('@HD'):
            tags = line.split('\t')[1:]
            return 'SO:queryname' in tags or 'GO:query' in tags
    return False


def start_reads(pipeline, c, scratch):
    """Start any read producers; return (bwa-mem3 read arguments, bwa-mem3 stdin)."""
    if c.get('bam'):
        # Interleaved pairs on bwa-mem3's stdin. samtools fastq keeps paired
        # primary reads only; its summary in fastq.log is used by check_pairing().
        # collate -f buffers mates in memory and spills only pairs whose mate is
        # not found within -r reads.
        bam = c['bam']
        fastq = ['samtools', 'fastq', '--reference', 'ref.fa', '-n', '-0', '/dev/null', '-s', '/dev/null']
        r, w = os.pipe()
        with os.fdopen(w, 'wb') as out, open(Path(scratch, 'fastq.log'), 'wb') as log:
            if query_grouped(bam):
                pipeline.start(fastq + ['-@', '2', bam], stdout=out, stderr=log)
            else:
                collate = pipeline.start(['samtools', 'collate', '-Ou', '-f', '-r', '1000000', '--reference', 'ref.fa',
                                          '--threads', '2', bam, scratch + '/collate'], stdout=subprocess.PIPE)
                pipeline.start(fastq + ['-'], stdin=collate.stdout, stdout=out, stderr=log)
                collate.stdout.close()
        return ['-p', 'ref.fa', '-'], os.fdopen(r, 'rb')

    if len(c['fastq_r1']) == 1:
        return ['ref.fa', c['fastq_r1'][0], c['fastq_r2'][0]], None

    # Several lanes: concatenate each mate's lanes into a FIFO. Concatenated gzip
    # files are one valid gzip stream, and bwa-mem3 rejects truncated or corrupt
    # gzip input itself, so the bytes are passed through undecoded.
    fifos = []
    for mate, paths in (('r1', c['fastq_r1']), ('r2', c['fastq_r2'])):
        if len({is_gzip(p) for p in paths}) > 1:
            raise ValueError('FASTQ lanes of one sample must be all gzip or all uncompressed')
        fifo = f'{scratch}/{mate}.fifo'
        os.mkfifo(fifo)
        pipeline.start(['sh', '-c', 'exec cat -- "$@" > "$0"', fifo, *paths])
        fifos.append(fifo)
    return ['ref.fa', *fifos], None


def start_alignment(pipeline, c, scratch):
    """Start bwa-mem3 | dupblaster | mako; return the process emitting sorted uncompressed BAM."""
    sample, cpu = c['sample_id'], str(c['cpu'])
    reads, stdin = start_reads(pipeline, c, scratch)
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tLB:{c['library']}\\tPL:{c['platform']}"
    # -K fixes the batch size so output does not depend on the thread count.
    # bwa-mem3 writes MC/MQ mate tags itself, which dupblaster needs.
    align = pipeline.start(['bwa-mem3', 'mem', '-K', '100000000', '-t', cpu, '-R', rg, *reads],
                           stdin=stdin, stdout=subprocess.PIPE)
    if stdin is not None:
        stdin.close()
    detection = 'on' if c['sequencing_duplicate_detection'] else 'off'
    mark = pipeline.start(['dupblaster', '--metrics-prefix', sample + '.dupblaster', '--sample', sample,
                           '--sequencing-duplicate-detection', detection, '-o', '-'],
                          stdin=align.stdout, stdout=subprocess.PIPE)
    align.stdout.close()
    # A coordinate sort emits nothing until its input ends; it spills to scratch.
    sort = pipeline.start(['mako', '-i', '-', '-o', '/dev/stdout', '--order', 'coordinate', '--threads', cpu,
                           '--compression-level', '0', '--tmp-dir', scratch],
                          stdin=mark.stdout, stdout=subprocess.PIPE)
    mark.stdout.close()
    return sort


def start_outputs(pipeline, c, scratch, sorted_bam):
    """tee the sorted stream into the CRAM writer and QC, so QC never re-reads the CRAM."""
    sample = c['sample_id']
    qc = {sample + '.flagstat.json': ['samtools', 'flagstat', '-O', 'json']}
    if c['detailed_stats']:
        qc[sample + '.samtools.stats.txt'] = ['samtools', 'stats', '-r', 'ref.fa']
    fifos = []
    for i, (path, args) in enumerate(qc.items()):
        fifo = f'{scratch}/qc{i}.fifo'
        os.mkfifo(fifo)
        with open(path, 'wb') as out:
            pipeline.start(args + [fifo], stdout=out, sink=True)
        fifos.append(fifo)
    tee = pipeline.start(['tee', *fifos], stdin=sorted_bam.stdout, stdout=subprocess.PIPE)
    sorted_bam.stdout.close()
    pipeline.start(['samtools', 'view', '--write-index', '-@', '4', '-C', '-T', 'ref.fa', '-o', sample + '.cram'],
                   stdin=tee.stdout, sink=True)
    tee.stdout.close()


def check_pairing(scratch, sample):
    # Each read given to the aligner yields exactly one primary record, so any
    # shortfall is singleton or unpaired (single-end, READ_OTHER) input that
    # samtools fastq discarded, e.g. a false query-grouped @HD claim.
    summary = re.search(r'processed (\d+) reads', Path(scratch, 'fastq.log').read_text(errors='replace'))
    if not summary:
        raise RuntimeError('samtools fastq summary not found')
    processed = int(summary.group(1))
    stats = json.loads(Path(sample + '.flagstat.json').read_text())
    aligned = sum(stats[k]['primary'] for k in ('QC-passed reads', 'QC-failed reads'))
    dropped = processed - aligned
    if dropped:
        print(f'genome-flow alignment: {dropped} of {processed} input reads were unpaired and not aligned', file=sys.stderr)
    if dropped > processed * 0.01:
        raise RuntimeError(f'{dropped} of {processed} input reads were unpaired; '
                           'input must be paired-end and its @HD grouping claim, if any, true')


def cancelled(signum, frame):
    raise RuntimeError(f'cancelled by signal {signum}')


def run(c):
    validate(c)
    stage_reference(c)
    signal.signal(signal.SIGTERM, cancelled)
    signal.signal(signal.SIGINT, cancelled)
    with tempfile.TemporaryDirectory(prefix='stream-', dir='.') as scratch:
        # Tools that use TMPDIR write to task scratch, not the container's /tmp.
        pipeline = Pipeline(env=dict(os.environ, TMPDIR=str(Path(scratch).resolve())))
        try:
            start_outputs(pipeline, c, scratch, start_alignment(pipeline, c, scratch))
            pipeline.wait()
            if c.get('bam'):
                check_pairing(scratch, c['sample_id'])
        finally:
            # Stop every tool before the scratch directory is removed; a second
            # signal must not interrupt that.
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
            signal.signal(signal.SIGINT, signal.SIG_IGN)
            pipeline.close()
            log = Path(scratch, 'fastq.log')
            if log.exists():
                sys.stderr.write(log.read_text(errors='replace'))


if __name__ == '__main__':
    try:
        run(json.loads(Path(sys.argv[1]).read_text()))
    except (ValueError, OSError, RuntimeError) as exc:
        print(f'genome-flow alignment: {exc}', file=sys.stderr)
        sys.exit(1)
