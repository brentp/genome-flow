version 1.0

task GLnexus {
  input {
    File fasta
    File fasta_index
    Array[File] gvcfs
    Array[File] gvcf_indexes
    String model_type
    String output_prefix
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb
  }

  # Reserve 25% for pipeline buffers and overhead; keep one memory input.
  Int memory_budget_gb = floor(memory_gb * 0.75)
  String config = if model_type == "WES" then "DeepVariantWES" else "DeepVariantWGS"

  command <<<
    set -euo pipefail
    if (( ~{memory_gb} < 2 )); then
      echo "GLnexus memory_gb must be at least 2" >&2
      exit 1
    fi
    ln -s "~{fasta}" ref.fa
    ln -s "~{fasta_index}" ref.fa.fai
    # norm re-sorts records it moves within --site-win; index-on-write fails on
    # any residual disorder instead of emitting an unsorted file.
    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2 glnexus_cli \
      --config ~{config} \
      --threads ~{cpu} \
      --mem-gbytes ~{memory_budget_gb} \
      --list "~{write_lines(gvcfs)}" \
      | bcftools norm --fasta-ref ref.fa --check-ref e --site-win 100000 \
          --write-index -Ob -o "~{output_prefix}.bcf"
    bcftools stats "~{output_prefix}.bcf" > "~{output_prefix}.bcftools.stats.txt"
  >>>

  output {
    File bcf = "~{output_prefix}.bcf"
    File bcf_index = "~{output_prefix}.bcf.csi"
    File stats = "~{output_prefix}.bcftools.stats.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    memory_mb: memory_gb * 1024
    disks: "local-disk ~{disk_gb} HDD"
    runtime_minutes: 1440
  }
}
