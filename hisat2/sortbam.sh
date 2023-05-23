#!/bin/bash
#hisat2 -p 12 -x ./hg19_idx/hg19 -1 ./sra/SRR3589956_1.fastq -2 ./sra/SRR3589956_2.fastq -S ./sam/SRR3589956.sam
for i in $(ls ../data | grep SRR1)
do
	nohup samtools sort -@8 ./bam/${i}/${i}.bam -o ./bam/${i}/${i}_sort.bam &
done
