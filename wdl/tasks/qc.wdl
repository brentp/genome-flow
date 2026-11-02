version 1.0

# Separate QC from alignment so a reporting failure never requires realignment.
# flagstat/samtools stats are computed inside the alignment pipe; this task is
# the only extra CRAM read: mosdepth coverage plus optional indexed Somalier sites.
task AlignmentQC {
  input {
    String sample_id
    File cram
    File cram_index
    File fasta
    File fasta_index
    File calling_bed
    File? somalier_sites
    String docker
    Int cpu
    Int memory_gb
    Int disk_gb
  }
  command <<<
    set -euo pipefail
    ln -s "~{fasta}" ref.fa
    ln -s "~{fasta_index}" ref.fa.fai
    ln -s "~{cram}" sample.cram
    ln -s "~{cram_index}" sample.cram.crai
    sites="~{default="" somalier_sites}"
    if [[ -n "${sites}" ]]; then
      # Sketch is named from the CRAM's @RG SM, which alignment sets to sample_id.
      somalier extract --sites "${sites}" --fasta ref.fa sample.cram &
      somalier_pid=$!
    fi
    mosdepth --threads ~{cpu} --fasta ref.fa --no-per-base \
      --by "~{calling_bed}" --thresholds 1,10,20,30 "~{sample_id}.coverage" sample.cram
    if [[ -n "${sites}" ]]; then
      wait "${somalier_pid}"
      test -s "~{sample_id}.somalier"
    fi
  >>>
  output {
    Array[File] coverage = glob("~{sample_id}.coverage*")
    Array[File] somalier_sketch = glob("~{sample_id}.somalier")
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    memory_mb: memory_gb * 1024
    disks: "local-disk ~{disk_gb} HDD"
    runtime_minutes: 360
  }
}

task SomalierRelate {
  input {
    Array[File] sketches
    File sites
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s "~{write_lines(sketches)}" samples.list
    somalier relate --sites "~{sites}" --output-prefix cohort samples.list
  >>>
  output {
    File samples = "cohort.samples.tsv"
    File pairs = "cohort.pairs.tsv"
    File report = "cohort.html"
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    memory_mb: 4096
    disks: "local-disk 10 HDD"
    runtime_minutes: 180
  }
}
