#!/bin/bash
##07/18/2024 EAH
##sum methylation over DMRs 
##v2 allows iteration over set files, have a data file named "data" and a regions file named "regions" in the directory you run the script from


path='/lab/solexa_gehring/elizabeth/allelic_emseq/sumby_2/'
srnas=$path'srna_data/'
outpath=$path'slurmout/'
data=$path'data/'
DMRs=$path'DMRs/'

cd $path

FILES=$(find $srnas -name "*.bed")
for file in $FILES ; do
	base=$(basename -- "$file" .bed)
	
	DMR=$(find $DMRs -type f ! -name "*.merge.bed")
	for region in $DMR ; do
		regbase=$(basename -- "$region" .bed )
		intout="$base""_int_""$regbase.bed"
		echo $intout
		sumout="$base""_sum_""$regbase.bed"
		echo $sumout
		
		#bedtools intersect -wao -a $region -b $file > $path$intout
		#bedtools groupby -i $intout -g 1,2,3 -c 7 -o sum > $path$sumout
	done
	DME=$(find $DMRs -type f -name "*.merge.bed")
	for region in $DME ; do
		regbase=$(basename -- "$region" .merge.bed )
		echo $intout
		sumout="$base""_sum_""$regbase.bed"
		echo $sumout
		
		#bedtools intersect -wao -a $region -b $file > $path$intout
		#bedtools groupby -i $intout -g 1,2,3 -c 7 -o sum > $path$sumout
	done
done
