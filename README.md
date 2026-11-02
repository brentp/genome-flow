WIP that may be abandoned. Trying to write an efficient, but simple workflow to go from FASTQ -> jointly called VCF using new tools like bwa-mem3.

# genome-flow [![test](https://github.com/brentp/genome-flow/actions/workflows/test.yml/badge.svg)](https://github.com/brentp/genome-flow/actions/workflows/test.yml)

Germline short-read variant calling for whole genomes and exomes:
**bwa-mem3 → duplicate marking → CRAM → DeepVariant → GLnexus**, with coverage,
alignment and sample-identity QC.

You need **Python 3, [miniwdl](https://miniwdl.readthedocs.io), and Docker on
Linux x86-64** (AVX2 CPU). Container images are pulled automatically.

## Quick start

1. Index the reference once (tools are in the container):

   ```sh
   docker run --rm -v /ref:/ref brentp/genome-flow-align:0.8.0 \
     sh -c 'samtools faidx /ref/GRCh38.fa && bwa-mem3 index /ref/GRCh38.fa'
   ```

2. Write a [sample sheet](#sample-sheet) and generate workflow inputs:

   ```sh
   python3 scripts/sheet_to_inputs.py samples.tsv \
     --fasta /ref/GRCh38.fa -o inputs.json
   ```

3. Run:

   ```sh
   miniwdl run wdl/genome_flow.wdl -i inputs.json --dir runs
   ```

For exomes, add `--model-type WES --capture-bed targets.bed`. For identity and
relatedness QC, add `--somalier-sites sites.hg38.vcf.gz` (a sites VCF with
INFO/AF matching your reference, e.g. from the
[Somalier releases](https://github.com/brentp/somalier/releases)).

## Sample sheet

Tab-separated, one row per sample ([example](examples/samples.tsv)). JSON is
also accepted ([example](examples/samples.json)).

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Unique name: letters, digits, `_`, `.`, `-`. |
| `library` | No | Library name; defaults to `sample_id`. |
| `platform` | No | Defaults to `ILLUMINA`. |
| `fastq_r1`, `fastq_r2` | FASTQ input | Comma-separated lanes, R2 in the same order as R1. |
| `bam` | BAM/CRAM input | Aligned or unaligned BAM/CRAM; reads are realigned. |
| `karyotype` | No | `XY` for haploid chrX/chrY calling; `XX` or empty calls everything diploid. |

Give each sample either FASTQ lanes or one BAM/CRAM. Lanes of one sample must
be all gzip-compressed or all uncompressed. Reads must be paired-end.

## Options

`scripts/sheet_to_inputs.py --help` lists all options. The most common:

| Option | Default | Purpose |
|---|---|---|
| `--fasta FILE` | required | Reference FASTA; index files are found next to it. |
| `--model-type WGS\|WES` | `WGS` | DeepVariant model. |
| `--bwa-mem3-indexes FILE ...` | next to FASTA | Index files elsewhere, plus optional `.alt` for ALT-aware alignment. |
| `--capture-bed FILE` | whole genome | Regions to call and report coverage for. |
| `--somalier-sites FILE` | off | Enables Somalier identity/relatedness QC. |
| `--par-bed FILE` | bundled GRCh38 | Pseudoautosomal regions, called diploid in `XY` samples. |
| `--haploid-contigs LIST` | `chrX,chrY` | Contigs called haploid in `XY` samples; `""` disables. |
| `--detailed-alignment-stats` | off | Also write `samtools stats` reports. |
| `--use-gpu` | off | GPU DeepVariant (requires a GPU backend). |
| `--cohort-id NAME` | `cohort` | Name for the joint-called output. |

**Sex chromosomes.** `XY` samples are called haploid on chrX/chrY outside the
pseudoautosomal regions. The bundled [`references/GRCh38.par.bed`](references/GRCh38.par.bed)
(DeepVariant's GRCh38 PAR file) is used automatically for `chr`-named GRCh38;
supply `--par-bed` for other references. The run stops immediately if the PAR
file does not match the reference.

**Resources.** Edit the generated JSON to change task sizes or container images,
for example `"GenomeFlow.align_cpu": 32`; [`wdl/genome_flow.wdl`](wdl/genome_flow.wdl)
lists every setting with its default. Defaults: alignment 16 CPUs / 48 GiB /
400 GB SSD; DeepVariant 8 CPUs / 32 GiB / 50 GB per shard, `scatter_count` 20
shards per sample; coverage QC 2 CPUs / 8 GiB / 50 GB; GLnexus 16 CPUs / 64 GiB.
Disk defaults leave room for a localized copy of a 30× CRAM. Set
`"GenomeFlow.sequencing_duplicate_detection": false` if read names are not
Illumina/Element format (e.g. SRA downloads).

## Outputs

| Output | Workflow keys |
|---|---|
| Alignments (CRAM + CRAI) | `crams`, `cram_indexes` |
| Variant calls per sample (BCF + CSI) | `sample_bcfs`, `sample_bcf_indexes` |
| gVCF per sample (BCF + CSI) | `sample_gbcfs`, `sample_gbcf_indexes` |
| Joint calls, two or more samples (BCF + CSI) | `cohort_bcf`, `cohort_bcf_index` |
| Duplicate metrics, flagstat, mosdepth coverage | `dup_metrics`, `flagstats`, `coverage_reports` |
| samtools stats (if requested) | `alignment_stats` |
| Variant stats | `sample_variant_stats`, `cohort_variant_stats` |
| Somalier (if sites given) | `somalier_sketches`, `somalier_samples`, `somalier_pairs`, `somalier_report` |
| Calling regions | `calling_bed` |

Variant calls are left-aligned and normalized against the reference. Somalier
reports metrics; review sex and relatedness against your sample sheet.

## Running on a cluster

For Slurm with Singularity, adapt [conf/cromwell.slurm.conf](conf/cromwell.slurm.conf)
and run:

```sh
java -Dconfig.file=conf/cromwell.slurm.conf -jar cromwell.jar run \
  wdl/genome_flow.wdl -i inputs.json
```

Other Cromwell backends (cloud) work with your existing configuration. Give
tasks fast local disk; alignment and DeepVariant write temporary files there.

## Development

```sh
make -C docker align                   # rebuild the tools image
python3 -m unittest discover -s tests  # unit tests
bash testdata/tiny/run_pipe.sh         # alignment integration tests in Docker
bash testdata/boundary/run.sh          # full workflow on a synthetic fixture
```

Set `IMAGE=...` to run both integration scripts against a locally built image.
[CI](.github/workflows/test.yml) runs all three test commands on every push and
pull request, using an image built from the commit.
