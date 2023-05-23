#!/bin/bash
#hisat2 -p 12 -x ./hg19_idx/hg19 -1 ./sra/SRR3589956_1.fastq -2 ./sra/SRR3589956_2.fastq -S ./sam/SRR3589956.sam
mkdir bam
for i in $(ls ../data | grep SRR1)
do
	mkdir bam/${i}
	samtools view -@8 -b ./sam/${i}/${i}.sam > ./bam/${i}/${i}.bam
done
