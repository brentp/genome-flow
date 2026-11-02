version 1.0

task DeepVariant {
  input {
    String sample_id
    String shard_id
    File cram
    File cram_index
    File fasta
    File fasta_index
    File intervals
    String model_type
    File? par_bed
    String haploid_contigs
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb
    Int gpu_count
    String gpu_type = "nvidia-tesla-t4"
  }

  String gpu = if gpu_count > 0 then "--gres=gpu:" + gpu_count else ""

  command <<<
    set -euo pipefail

    ln -sf "~{fasta}" ref.fa
    ln -sf "~{fasta_index}" ref.fa.fai
    ln -sf "~{cram}" "~{sample_id}.cram"
    ln -sf "~{cram_index}" "~{sample_id}.cram.crai"

    # haploid_contigs is empty except for XY samples (resolved by MakeIntervals,
    # which also requires par_bed for them).
    extra=()
    if [[ -n "~{haploid_contigs}" ]]; then
      extra+=(--haploid_contigs="~{haploid_contigs}")
      par_bed_path="~{default="" par_bed}"
      if [[ -n "${par_bed_path}" ]]; then
        extra+=(--par_regions_bed="${par_bed_path}")
      fi
    fi

    # Keep make_examples/call_variants intermediates on the task disk, not the
    # container's /tmp; they are not outputs.
    intermediate="$PWD/dv-intermediate"
    mkdir "${intermediate}"
    trap 'rm -rf "${intermediate}"' EXIT

    /opt/deepvariant/bin/run_deepvariant \
      --model_type="~{model_type}" \
      --ref=ref.fa \
      --reads="~{sample_id}.cram" \
      --regions="~{intervals}" \
      --output_vcf="~{sample_id}.~{shard_id}.vcf.gz" \
      --output_gvcf="~{sample_id}.~{shard_id}.g.vcf.gz" \
      --intermediate_results_dir="${intermediate}" \
      --num_shards=~{cpu} \
      --vcf_stats_report=false \
      "${extra[@]}"

    # DeepVariant normally writes these indexes; only build missing/empty ones.
    if [[ ! -s "~{sample_id}.~{shard_id}.vcf.gz.tbi" ]]; then
      bcftools index -f -t "~{sample_id}.~{shard_id}.vcf.gz"
    fi
    if [[ ! -s "~{sample_id}.~{shard_id}.g.vcf.gz.tbi" ]]; then
      bcftools index -f -t "~{sample_id}.~{shard_id}.g.vcf.gz"
    fi
  >>>

  output {
    File vcf = "~{sample_id}.~{shard_id}.vcf.gz"
    File vcf_index = "~{sample_id}.~{shard_id}.vcf.gz.tbi"
    File gvcf = "~{sample_id}.~{shard_id}.g.vcf.gz"
    File gvcf_index = "~{sample_id}.~{shard_id}.g.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    memory_mb: memory_gb * 1024
    disks: "local-disk ~{disk_gb} HDD"
    gpuCount: gpu_count
    gpuType: gpu_type
    gpu: gpu_count > 0
    gres: gpu
    runtime_minutes: 1440
  }
}

task MergeSampleCalls {
  input {
    String sample_id
    File fasta
    File fasta_index
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[File] gvcfs
    Array[File] gvcf_indexes
    String docker
    Int cpu = 2
    Int memory_gb = 8
    Int disk_gb = 50
  }

  command <<<
    set -euo pipefail

    ln -s "~{fasta}" ref.fa
    ln -s "~{fasta_index}" ref.fa.fai
    # Uncompressed BCF between stages; only the final BCF touches disk. norm
    # re-sorts records it moves within --site-win; index-on-write fails on any
    # residual disorder instead of emitting an unsorted file.
    bcftools concat -a -Ou --file-list "~{write_lines(vcfs)}" \
      | bcftools norm --fasta-ref ref.fa --check-ref e --site-win 100000 \
          --write-index -Ob -o "~{sample_id}.bcf"
    bcftools concat -a --write-index --threads ~{cpu} -Ob -o "~{sample_id}.g.bcf" --file-list "~{write_lines(gvcfs)}"
    bcftools stats "~{sample_id}.bcf" > "~{sample_id}.bcftools.stats.txt"
  >>>

  output {
    File stats = "~{sample_id}.bcftools.stats.txt"
    File bcf = "~{sample_id}.bcf"
    File bcf_index = "~{sample_id}.bcf.csi"
    File gbcf = "~{sample_id}.g.bcf"
    File gbcf_index = "~{sample_id}.g.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    memory_mb: memory_gb * 1024
    disks: "local-disk ~{disk_gb} HDD"
    runtime_minutes: 180
  }
}
