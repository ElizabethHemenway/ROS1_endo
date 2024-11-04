#!/bin/bash
# EAH 04/02/2024
# use to make files for IGV browsing with negative -0.1 instead of 0 so you can see where coverage is

#chromosome size file for generating .bw files to browse data using IGV - replace with your own if needed.
chrom_sizes=/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/TAIR10.chrom.sizes 

FILES=$(find . -name "*.min5.bedGraph")
for file in $FILES ; do
	echo "input file = ""$file"
	forIGV="${file%.*}.forIGV.bedGraph"
	rm -f $forIGV "${forIGV%.*}.bw"
	echo "generate file = ""$forIGV"
	awk -F'\t' '{ $4 = ($4 == 0 ? -10 : $4) } 1' OFS="\t" "$file" > $forIGV
	bedGraphToBigWig $forIGV "$chrom_sizes" "${forIGV%.*}.bw"
done
