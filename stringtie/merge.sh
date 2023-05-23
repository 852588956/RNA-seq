#!/bin/bash
#hisat2 -p 12 -x ./hg19_idx/hg19 -1 ./sra/SRR3589956_1.fastq -2 ./sra/SRR3589956_2.fastq -S ./sam/SRR3589956.sam
stringtie --merge -p 8 -G ../index/hg38_tran/hg38_ucsc.annotated.gtf -o stringtie_merged.gtf merglist.txt
python prepDE.py -i sample_lst.txt -g gene_count_matrix.csv -t transcript_count_matrix.csv
