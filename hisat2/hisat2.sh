#!/bin/bash
#hisat2 -p 12 -x ./hg19_idx/hg19 -1 ./sra/SRR3589956_1.fastq -2 ./sra/SRR3589956_2.fastq -S ./sam/SRR3589956.sam
mkdir sam
for i in $(ls ../data | grep SRR1)
do
	mkdir sam/${i}
	nohup hisat2 -p 12 -x ../index/hg38_tran/genome_tran -1 ../trimmomatic/${i}/${i}_paired_R1.fq.gz -2 ../trimmomatic/${i}/${i}_paired_R2.fq.gz -S ./sam/${i}/${i}.sam > ${i}_hisat2.out &

done
