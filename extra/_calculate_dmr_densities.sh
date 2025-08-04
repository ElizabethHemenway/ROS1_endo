#!/bin/bash
# 06/11/24 - calculate "coverage" of DMRs across chromosome windows to plot as distribution

path=./
DMRs=$path'DMRs/'
outpath=$path'slurmout/'

# make genome windows
window_len=100000
gwin="TAIR10_"$window_len"bp0step.bed"

bedtools makewindows -g TAIR10.chrom.sizes -w $window_len > $gwin


# Loop over all files to calculate coverage
for fileB in "$DMRs"/*; do
    baseB=$(basename "$fileB")
    output_file="$path${baseB%.bed}_cov_"$window_len"bpwindows.bed"
	echo $output_file
	# compare to genome windows
	bedtools coverage -a $gwin -b "$fileB" > $output_file 
        
done
    
    
#also sumby mCG over same windows

data=$path'mC_data/'

FILES=$(find $data -name "*.bed")
for file in $FILES ; do
	base=$(basename -- "$file" .bed)
	output="$base""_sumby_"$window_len"bpwindows.bed"
	echo $output
	sumby=$path'sumByFeature.sh'
	sbatch -p 20 --job-name=sumby_mC --output $outpath"$output"'.out' --mem=8gb --wrap " $sumby -i "$file" -r "$gwin" -o "$output" -m 1 -x 5 "
done