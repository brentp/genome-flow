version 1.0

# Split the calling regions into DeepVariant shards and resolve each sample's
# haploid contigs. See scripts/make_intervals.py for the rules.
task MakeIntervals {
  input {
    File fasta_index
    File? capture_bed
    Int scatter_count
    String? haploid_override
    Array[String] karyotypes
    File? par_bed
    String docker
  }

  File manifest = write_json(object {
    fasta_index: fasta_index, capture_bed: capture_bed, scatter_count: scatter_count,
    haploid_override: haploid_override, karyotypes: karyotypes, par_bed: par_bed
  })

  command <<<
    set -euo pipefail
    python3 /usr/local/bin/genome-flow-make-intervals "~{manifest}"
  >>>

  output {
    Array[File] interval_beds = read_lines("shard_list.txt")
    File calling_bed = "calling.bed"
    Array[String] sample_haploid_contigs = read_json("sample_haploid_contigs.json")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GiB"
    memory_mb: 2048
    disks: "local-disk 10 HDD"
    runtime_minutes: 60
  }
}
