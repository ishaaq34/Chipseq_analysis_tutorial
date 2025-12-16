#!/bin/bash
set -euo pipefail

mkdir -p bowalign

while read -r sample; do
  echo "Aligning $sample..."

  bowtie2 -x index \
  -U ${sample}.clean.fq.gz \
  -p 6 --no-unal \
  2> bowalign/${sample}.log \
| samtools sort -@ 6 -o bowalign/${sample}.sorted.bam

  samtools index bowalign/${sample}.sorted.bam
  
done < sample_id.txt

2nd)

#!/bin/bash

mkdir -p bowalign

while read -r sample; do
  echo "Aligning $sample"

  bowtie2 -x hg38_index \
    -U clean_${sample}.fastq.gz \
    -p 6 --no-unal \
  | samtools sort -@ 6 -o bowalign/${sample}.sorted.bam

  samtools index bowalign/${sample}.sorted.bam

done < sample_id.txt

