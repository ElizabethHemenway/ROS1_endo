#!/bin/bash
#	sort and index bam files
#	04/03/2024 EAH
#	/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/BAMindex.sh


FILES=$(find . -name "*.bam")
for file in $FILES ; do
	JID_JOB1=`sbatch -p 20 --mem=10gb --job-name=samsort --wrap " samtools sort "$file" -o "${file%.*}.sort.bam" " | cut -d " " -f 4`
 	sbatch -p 20 --mem=5gb --dependency=afterok:$JID_JOB1 --job-name=samindex --wrap " samtools index "${file%.*}.sort.bam" "
done
