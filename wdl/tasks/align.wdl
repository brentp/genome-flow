version 1.0

# Stream reads through bwa-mem3, duplicate marking and sorting into a CRAM.
# flagstat (and optional samtools stats) read the sorted stream, not the CRAM.
# See scripts/stream_align.py.
task AlignMarkSort {
  input {
    String sample_id
    String library
    String platform
    Array[File] fastq_r1
    Array[File] fastq_r2
    File? bam
    File fasta
    File fasta_index
    Array[File] bwa_mem3_indexes
    Boolean sequencing_duplicate_detection
    Boolean detailed_stats
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  File manifest = write_json(object {
    sample_id: sample_id, library: library, platform: platform,
    fastq_r1: fastq_r1, fastq_r2: fastq_r2, bam: bam,
    fasta: fasta, fasta_index: fasta_index, bwa_mem3_indexes: bwa_mem3_indexes,
    sequencing_duplicate_detection: sequencing_duplicate_detection,
    detailed_stats: detailed_stats, cpu: cpu
  })

  command <<<
    set -euo pipefail
    python3 /usr/local/bin/genome-flow-stream-align "~{manifest}"
  >>>

  output {
    File cram = "~{sample_id}.cram"
    File cram_index = "~{sample_id}.cram.crai"
    File dup_metrics = "~{sample_id}.dupblaster.duplicate-metrics.tsv"
    File flagstat = "~{sample_id}.flagstat.json"
    Array[File] alignment_stats = glob("~{sample_id}.samtools.stats.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    memory_mb: memory_gb * 1024
    # The sort spills here and its final merge reads every spill file at once.
    disks: "local-disk ~{disk_gb} SSD"
    runtime_minutes: 1440
  }
}
