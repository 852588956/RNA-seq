#!/bin/bash
#hisat2 -p 12 -x ./hg19_idx/hg19 -1 ./sra/SRR3589956_1.fastq -2 ./sra/SRR3589956_2.fastq -S ./sam/SRR3589956.sam
mkdir gtf
for i in $(ls ../data | grep SRR1)
do
	mkdir gtf/${i}
	nohup stringtie -p 4 -e \
	       -G ../index/hg38_tran/hg38_ucsc.annotated.gtf \
	       -o gtf/${i}/${i}.gtf \
	       ../hisat2/bam/${i}/${i}_sort.bam &

done
