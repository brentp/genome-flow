version 1.0

import "structs.wdl"
import "tasks/align.wdl" as align
import "tasks/intervals.wdl" as intervals
import "tasks/deepvariant.wdl" as dv
import "tasks/glnexus.wdl" as jnt
import "tasks/qc.wdl" as qc

# Germline short-read WGS/WES:
#   per sample: bwa-mem3 | dupblaster | mako -> CRAM -> DeepVariant (scattered) -> BCF + gVCF BCF
#   per cohort: GLnexus joint calling (two or more samples), Somalier relatedness (optional)
# Task defaults live here only; tasks take every setting as an input.
workflow GenomeFlow {
  input {
    Array[Sample] samples

    File fasta
    File fasta_index
    # .amb .ann .pac .bwt.2bit.64, plus optional .alt for ALT-aware alignment.
    Array[File] bwa_mem3_indexes

    String model_type = "WGS"
    # Regions to call and report coverage for; default is the whole reference.
    File? capture_bed
    # Pseudoautosomal regions, called diploid in XY samples. Required when
    # XY samples are called haploid.
    File? par_bed
    # Haploid contigs for XY samples. Default: chrX/chrY present in the
    # reference; "" calls XY samples diploid.
    String? haploid_contigs
    File? somalier_sites
    Boolean use_gpu = false
    # Turn off when read names are not Illumina/Element format (e.g. SRA).
    Boolean sequencing_duplicate_detection = true
    Boolean detailed_alignment_stats = false
    String cohort_id = "cohort"

    # DeepVariant tasks per sample; each covers an equal share of the calling
    # regions. Each task costs ~20 s of fixed startup (HG003 35x chr20, 8 CPUs),
    # small next to a 1/20 WGS share; on backends that copy inputs, each task
    # also copies the whole CRAM.
    Int scatter_count = 20

    String align_docker = "brentp/genome-flow-align:0.8.0"
    String deepvariant_docker = "google/deepvariant:1.10.0@sha256:962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6"
    String deepvariant_gpu_docker = "google/deepvariant:1.10.0-gpu@sha256:78fb1bb960f3b74a4308b62ba8b6232723d4c75e16a4a50f0d810d7086467dd1"

    Int align_cpu = 16
    Int align_memory_gb = 48
    Int align_disk_gb = 400

    Int qc_cpu = 2
    Int qc_memory_gb = 8
    # Holds a copy of the CRAM (~15-20 GB at 30x) on backends that copy inputs.
    Int qc_disk_gb = 50

    Int dv_cpu = 8
    Int dv_memory_gb = 32
    # Holds a copy of the CRAM on backends that copy inputs. Intermediate files
    # are ~3 MB per Mb called at 35x, i.e. under 1 GB for a 1/20 WGS shard.
    Int dv_disk_gb = 50
    Int dv_gpu_count = 1

    Int glnexus_cpu = 16
    Int glnexus_memory_gb = 64
    Int glnexus_disk_gb = 200
  }

  scatter (s in samples) {
    String karyotype = select_first([s.karyotype, ""])
  }

  call intervals.MakeIntervals {
    input:
      fasta_index = fasta_index,
      capture_bed = capture_bed,
      scatter_count = scatter_count,
      haploid_override = haploid_contigs,
      karyotypes = karyotype,
      par_bed = par_bed,
      docker = align_docker
  }

  scatter (i in range(length(samples))) {
    Sample sample = samples[i]

    call align.AlignMarkSort {
      input:
        sample_id = sample.sample_id,
        library = sample.library,
        platform = sample.platform,
        fastq_r1 = sample.fastq_r1,
        fastq_r2 = sample.fastq_r2,
        bam = sample.bam,
        fasta = fasta,
        fasta_index = fasta_index,
        bwa_mem3_indexes = bwa_mem3_indexes,
        sequencing_duplicate_detection = sequencing_duplicate_detection,
        detailed_stats = detailed_alignment_stats,
        docker = align_docker,
        cpu = align_cpu,
        memory_gb = align_memory_gb,
        disk_gb = align_disk_gb
    }

    call qc.AlignmentQC {
      input:
        sample_id = sample.sample_id,
        cram = AlignMarkSort.cram,
        cram_index = AlignMarkSort.cram_index,
        fasta = fasta,
        fasta_index = fasta_index,
        calling_bed = MakeIntervals.calling_bed,
        somalier_sites = somalier_sites,
        docker = align_docker,
        cpu = qc_cpu,
        memory_gb = qc_memory_gb,
        disk_gb = qc_disk_gb
    }

    scatter (bed in MakeIntervals.interval_beds) {
      call dv.DeepVariant {
        input:
          sample_id = sample.sample_id,
          shard_id = basename(bed, ".bed"),
          cram = AlignMarkSort.cram,
          cram_index = AlignMarkSort.cram_index,
          fasta = fasta,
          fasta_index = fasta_index,
          intervals = bed,
          model_type = model_type,
          par_bed = par_bed,
          haploid_contigs = MakeIntervals.sample_haploid_contigs[i],
          docker = if use_gpu then deepvariant_gpu_docker else deepvariant_docker,
          cpu = dv_cpu,
          memory_gb = dv_memory_gb,
          disk_gb = dv_disk_gb,
          gpu_count = if use_gpu then dv_gpu_count else 0
      }
    }

    call dv.MergeSampleCalls {
      input:
        sample_id = sample.sample_id,
        fasta = fasta,
        fasta_index = fasta_index,
        vcfs = DeepVariant.vcf,
        vcf_indexes = DeepVariant.vcf_index,
        gvcfs = DeepVariant.gvcf,
        gvcf_indexes = DeepVariant.gvcf_index,
        docker = align_docker
    }
  }

  if (defined(somalier_sites)) {
    call qc.SomalierRelate {
      input:
        sketches = flatten(AlignmentQC.somalier_sketch),
        sites = select_first([somalier_sites]),
        docker = align_docker
    }
  }

  if (length(samples) >= 2) {
    call jnt.GLnexus {
      input:
        fasta = fasta,
        fasta_index = fasta_index,
        gvcfs = MergeSampleCalls.gbcf,
        gvcf_indexes = MergeSampleCalls.gbcf_index,
        model_type = model_type,
        output_prefix = cohort_id,
        docker = align_docker,
        cpu = glnexus_cpu,
        memory_gb = glnexus_memory_gb,
        disk_gb = glnexus_disk_gb
    }
  }

  output {
    Array[File] crams = AlignMarkSort.cram
    Array[File] cram_indexes = AlignMarkSort.cram_index
    Array[File] dup_metrics = AlignMarkSort.dup_metrics
    Array[File] flagstats = AlignMarkSort.flagstat
    Array[File] alignment_stats = flatten(AlignMarkSort.alignment_stats)
    Array[Array[File]] coverage_reports = AlignmentQC.coverage
    Array[File] sample_variant_stats = MergeSampleCalls.stats
    Array[File] somalier_sketches = flatten(AlignmentQC.somalier_sketch)
    File? somalier_samples = SomalierRelate.samples
    File? somalier_pairs = SomalierRelate.pairs
    File? somalier_report = SomalierRelate.report
    File? cohort_variant_stats = GLnexus.stats
    File calling_bed = MakeIntervals.calling_bed
    Array[File] sample_bcfs = MergeSampleCalls.bcf
    Array[File] sample_bcf_indexes = MergeSampleCalls.bcf_index
    Array[File] sample_gbcfs = MergeSampleCalls.gbcf
    Array[File] sample_gbcf_indexes = MergeSampleCalls.gbcf_index
    File? cohort_bcf = GLnexus.bcf
    File? cohort_bcf_index = GLnexus.bcf_index
  }
}
