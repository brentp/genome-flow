#!/usr/bin/env python3
"""Split the calling regions into DeepVariant shards and resolve haploid contigs per sample.

Usage: make_intervals.py manifest.json  (written by wdl/tasks/intervals.wdl), with keys
  fasta_index      reference .fai; its contig order is the output order
  capture_bed      BED of regions to call, or null for the whole reference
  scatter_count    number of shards wanted
  haploid_override comma-separated haploid contigs for XY samples, or null for
                   chrX/chrY when present in the reference ("" disables)
  karyotypes       one of "", "XX", "XY" per sample
  par_bed          pseudoautosomal regions BED, required when XY samples are
                   called haploid, or null

Writes, in the working directory:
  calling.bed                  merged calling regions (also used for coverage QC)
  shards/shard_NNNN.bed        contiguous shards with equal base counts (within one base)
  shard_list.txt               the shard paths, in reference order
  sample_haploid_contigs.json  per sample: the haploid contigs for XY samples, else ""
"""
import json
import os
import sys


def load_fai(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            name, length, *_ = line.split()
            if name in lengths or int(length) < 1:
                sys.exit(f"invalid or duplicate FASTA index contig: {name}")
            lengths[name] = int(length)
    return lengths


def load_bed(path, lengths):
    """Validated intervals sorted in reference order, with overlapping/adjacent ones merged."""
    intervals = []
    rank = {chrom: i for i, chrom in enumerate(lengths)}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            try:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            except (IndexError, ValueError):
                sys.exit(f"{path}:{lineno}: expected BED chrom, start, end")
            if chrom not in lengths:
                sys.exit(f"{path}:{lineno}: unknown reference contig {chrom}")
            if not 0 <= start < end <= lengths[chrom]:
                sys.exit(f"{path}:{lineno}: interval outside contig bounds or empty")
            intervals.append((chrom, start, end))
    intervals.sort(key=lambda x: (rank[x[0]], x[1], x[2]))
    merged = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def haploid_contigs(override, lengths):
    if override is None:
        return [chrom for chrom in ("chrX", "chrY") if chrom in lengths]
    contigs = [chrom.strip() for chrom in override.split(",") if chrom.strip()]
    unknown = [chrom for chrom in contigs if chrom not in lengths]
    if unknown:
        sys.exit(f"haploid contigs not in the FASTA index: {', '.join(unknown)}")
    return contigs


def split_evenly(regions, n):
    """Cut regions into at most n contiguous shards whose base counts differ by at most one."""
    total = sum(end - start for _, start, end in regions)
    n = min(n, total)
    size, extra = divmod(total, n)
    shards = [[]]
    remaining = size + (extra > 0)
    for chrom, start, end in regions:
        while start < end:
            if remaining == 0:
                shards.append([])
                remaining = size + (len(shards) <= extra)
            stop = min(end, start + remaining)
            shards[-1].append((chrom, start, stop))
            remaining -= stop - start
            start = stop
    return shards


def write_bed(path, intervals):
    with open(path, "w") as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")


def main(args):
    if args["scatter_count"] < 1:
        sys.exit("scatter_count must be positive")
    lengths = load_fai(args["fasta_index"])

    # DeepVariant: haploid contigs apply only to XY samples ("for samples with XX
    # karyotype --haploid_contigs flag should not be used"), and XY samples need
    # PAR regions so the PARs stay diploid.
    karyotypes = args["karyotypes"]
    invalid = sorted({k for k in karyotypes if k not in ("", "XX", "XY")})
    if invalid:
        sys.exit(f"karyotype must be XX, XY or empty, got: {', '.join(invalid)}")
    haploid = haploid_contigs(args["haploid_override"], lengths)
    if haploid and "XY" in karyotypes:
        if not args["par_bed"]:
            sys.exit("XY samples need par_bed for haploid calling of "
                     f"{','.join(haploid)}; supply par_bed or set haploid_contigs to \"\"")
        # Contig names and bounds reject a PAR BED from another build or naming scheme.
        load_bed(args["par_bed"], lengths)

    if args["capture_bed"]:
        regions = load_bed(args["capture_bed"], lengths)
    else:
        regions = [(chrom, 0, length) for chrom, length in lengths.items()]
    if not regions:
        sys.exit("no intervals to scatter")

    os.makedirs("shards", exist_ok=True)
    paths = []
    for i, shard in enumerate(split_evenly(regions, args["scatter_count"])):
        paths.append(f"shards/shard_{i:04d}.bed")
        write_bed(paths[-1], shard)
    write_bed("calling.bed", regions)
    with open("shard_list.txt", "w") as fh:
        fh.write("\n".join(paths) + "\n")
    with open("sample_haploid_contigs.json", "w") as fh:
        json.dump([",".join(haploid) if k == "XY" else "" for k in karyotypes], fh)


if __name__ == "__main__":
    with open(sys.argv[1]) as fh:
        main(json.load(fh))
