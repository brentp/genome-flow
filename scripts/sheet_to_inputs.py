#!/usr/bin/env python3
"""Convert a sample sheet (TSV or JSON) into Cromwell/miniwdl inputs JSON.

Task sizes and container images keep the workflow defaults; edit the
generated JSON to change them (e.g. "GenomeFlow.align_cpu": 32).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

GRCH38_PAR_BED = Path(__file__).resolve().parents[1] / "references" / "GRCh38.par.bed"
NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def parse_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def load_tsv(path: Path) -> list[dict[str, Any]]:
    """One row per sample; lanes are comma-separated. Empty and #-prefixed rows are skipped."""
    rows = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"{path}: empty TSV")
        reader.fieldnames = [name.lstrip("#").strip() for name in reader.fieldnames]
        if "sample_id" not in reader.fieldnames:
            raise SystemExit(f"{path}: missing sample_id column")
        for row in reader:
            row = {key: (value or "").strip() for key, value in row.items() if key}
            if not row.get("sample_id") or row["sample_id"].startswith("#"):
                continue
            row["fastq_r1"] = parse_list(row.get("fastq_r1", ""))
            row["fastq_r2"] = parse_list(row.get("fastq_r2", ""))
            rows.append(row)
    return rows


def load_json_sheet(path: Path) -> list[dict[str, Any]]:
    """A list of sample objects, or {"samples": [...]}."""
    data = json.loads(path.read_text())
    if isinstance(data, dict) and "samples" in data:
        data = data["samples"]
    if not isinstance(data, list) or not all(isinstance(item, dict) for item in data):
        raise SystemExit(f"{path}: expected a list of sample objects or {{'samples': [...]}}")
    return data


def validate_samples(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Fill defaults and check every sample; the result is GenomeFlow.samples."""
    if not rows:
        raise SystemExit("sample sheet has no samples")
    samples = []
    seen = set()
    for row in rows:
        sid = row.get("sample_id")
        sample = {
            "sample_id": sid,
            "library": row.get("library") or sid,
            "platform": row.get("platform") or "ILLUMINA",
            "fastq_r1": row.get("fastq_r1") or [],
            "fastq_r2": row.get("fastq_r2") or [],
        }
        for field in ("sample_id", "library", "platform"):
            if not isinstance(sample[field], str) or not NAME.fullmatch(sample[field]):
                raise SystemExit(f"invalid {field} {sample[field]!r}: use letters, digits, _, . or -")
        if sid in seen:
            raise SystemExit(f"duplicate sample_id: {sid}")
        seen.add(sid)
        r1, r2, bam = sample["fastq_r1"], sample["fastq_r2"], row.get("bam") or None
        if not isinstance(r1, list) or not isinstance(r2, list) or not all(isinstance(p, str) and p for p in r1 + r2):
            raise SystemExit(f"{sid}: fastq_r1/fastq_r2 must be lists of nonempty paths")
        if len(r1) != len(r2):
            raise SystemExit(f"{sid}: fastq_r1 and fastq_r2 need the same number of lanes")
        if bam is not None and not isinstance(bam, str):
            raise SystemExit(f"{sid}: bam must be a path")
        if bool(r1) == bool(bam):
            raise SystemExit(f"{sid}: provide FASTQ lanes or bam, not both or neither")
        if bam:
            sample["bam"] = bam
        karyotype = row.get("karyotype") or ""
        if not isinstance(karyotype, str) or karyotype.upper() not in ("", "XX", "XY"):
            raise SystemExit(f"{sid}: karyotype must be XX, XY or empty, got {karyotype!r}")
        if karyotype:
            sample["karyotype"] = karyotype.upper()
        samples.append(sample)
    return samples


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("sheet", type=Path, help="TSV or JSON sample sheet")
    p.add_argument("--fasta", required=True, help="Reference FASTA")
    p.add_argument("--fasta-index", help="Default: <fasta>.fai")
    p.add_argument("--bwa-mem3-indexes", nargs="+",
                   help="Index files, including optional .alt/.hdr/.dict. Default: <fasta>.{amb,ann,pac,bwt.2bit.64}")
    p.add_argument("--model-type", default="WGS", choices=["WGS", "WES"])
    p.add_argument("--capture-bed", help="Regions to call and report coverage for. Default: whole reference")
    p.add_argument("--par-bed", help="PAR regions for XY samples. Default with XY samples: bundled references/GRCh38.par.bed")
    p.add_argument("--somalier-sites", help="Matching sites VCF with INFO/AF; enables identity QC")
    p.add_argument(
        "--haploid-contigs",
        help='Haploid contigs for XY samples. Default: chrX/chrY present in the FAI '
        '(inferred by the workflow); use "" to call XY samples diploid',
    )
    p.add_argument("--use-gpu", action="store_true", help="GPU DeepVariant (requires a GPU backend)")
    p.add_argument("--detailed-alignment-stats", action="store_true", help="Also write samtools stats reports")
    p.add_argument("--scatter-count", type=int, default=20, help="DeepVariant tasks per sample")
    p.add_argument("--cohort-id", default="cohort", help="Name for the joint-called output")
    p.add_argument("-o", "--output", type=Path, help="Default: stdout")
    args = p.parse_args()
    if args.scatter_count < 1:
        p.error("--scatter-count must be positive")

    if args.sheet.suffix.lower() == ".json":
        samples = validate_samples(load_json_sheet(args.sheet))
    else:
        samples = validate_samples(load_tsv(args.sheet))
    if (any(s.get("karyotype") == "XY" for s in samples) and not args.par_bed
            and args.haploid_contigs != ""):
        # The bundled BED is GRCh38/hg38 with chr names; the workflow rejects it
        # if its contigs or bounds do not fit the reference.
        args.par_bed = str(GRCH38_PAR_BED)

    inputs: dict[str, Any] = {
        "GenomeFlow.samples": samples,
        "GenomeFlow.fasta": args.fasta,
        "GenomeFlow.fasta_index": args.fasta_index or args.fasta + ".fai",
        "GenomeFlow.bwa_mem3_indexes": args.bwa_mem3_indexes
        or [args.fasta + suffix for suffix in (".amb", ".ann", ".pac", ".bwt.2bit.64")],
        "GenomeFlow.model_type": args.model_type,
        "GenomeFlow.use_gpu": args.use_gpu,
        "GenomeFlow.detailed_alignment_stats": args.detailed_alignment_stats,
        "GenomeFlow.scatter_count": args.scatter_count,
        "GenomeFlow.cohort_id": args.cohort_id,
    }
    optional = {
        "GenomeFlow.capture_bed": args.capture_bed,
        "GenomeFlow.somalier_sites": args.somalier_sites,
        "GenomeFlow.haploid_contigs": args.haploid_contigs,
        "GenomeFlow.par_bed": args.par_bed,
    }
    inputs.update({key: value for key, value in optional.items() if value is not None})

    text = json.dumps(inputs, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
