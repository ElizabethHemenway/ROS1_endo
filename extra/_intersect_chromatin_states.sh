#!/bin/bash
# 06/03/2025 - compare DMRs to chromatin features of interest

path=./
DMRs=$path"DMRs/"
genes=$path"genes/"

states9=$path'2compare/Chromatin_states_2014.bed'

# Loop over all files in folder A
for fileA in "$DMRs"/*; do
    baseA=$(basename "$fileA")
  
    baseB='chromatin_states_2014'
    # Construct output file name
    output_file="$path${baseA%.bed}_chromatin_states_2014_intersect.bed"
    echo $output_file
    #rm -f $output_file
    # Run bedtools intersect
    bedtools intersect -wa -wb -a "$fileA" -b $states9 > "$output_file"
    
done
