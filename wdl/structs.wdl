version 1.0

# One sample. Provide either FASTQ R1/R2 (one or more lane pairs) or a BAM/CRAM.
# BAM/CRAM is always converted to FASTQ and realigned.
# karyotype: "XY" enables haploid sex-chromosome calling; "XX" or omitted is diploid.
struct Sample {
  String sample_id
  String library
  String platform
  Array[File] fastq_r1
  Array[File] fastq_r2
  File? bam
  String? karyotype
}
