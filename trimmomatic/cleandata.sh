for i in $(ls ../data/ | grep SRR1)
do
	mkdir ${i}
        nohup trimmomatic PE \
	    ../data/${i}/${i}_1.fastq ../data/${i}/${i}_2.fastq \
            ${i}/${i}_paired_R1.fq.gz ${i}/${i}_unpairedd_R1.fq.gz \
	    ${i}/${i}_paired_R2.fq.gz ${i}/${i}_unpairedd_R2.fq.gz \
	    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
	    LEADING:3 \
	    TRAILING:3 \
	    SLIDINGWINDOW:4:15 \
	    MINLEN:36 > ${i}_trimmomatic.out &

done
