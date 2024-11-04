#!/bin/bash
# 04/08/2024 EAH
# make pseudogenome using bedtools maskfasta
# usage: make_pseudogenome.sh ref_genome ref_name alt_name output_directory scripts_directory 
# example: make_pseudogenome.sh ./TAIR10.fa T10 C24 /lab/solexa_gehring/elizabeth/pseudoG/ /lab/solexa_gehring/elizabeth/scripts/

syri_out="$1" # syri output file (if you hav a snp bed file, change columns=TSV to columns=bed in syri_snps2bed.py line 22)
ref_name="$2" # name of reference
alt_name="$3" # name of alt
reference_genome="$4" # fasta file of reference genome
outdir="$4" # output directory
scripts="$5" # must include helper script syri_snps2bed.py

name=$ref_name'_'$alt_name

cd $outdir

# run python script to make output bed files for each SNP base
$scripts'syri_snps2bed.py' $syri_out $ref_name $alt_name $outdir

# use bedtools maskfasta to change each SNP from reference to alt base 
bedtools maskfasta -fi $reference_genome -bed $name'_2A.bed' -fo $name'_2A.fa' -mc A
bedtools maskfasta -fi $name'_2A.fa' -bed $name'_2T.bed' -fo $name'_2AT.fa' -mc T
bedtools maskfasta -fi $name'_2AT.fa' -bed $name'_2C.bed' -fo $name'_2ATC.fa' -mc C
bedtools maskfasta -fi $name'_2ATC.fa' -bed $name'_2G.bed' -fo $name'_2ATCG.fa' -mc G

# make an index for pseudogenome
samtools faidx $name'_2ATCG.fa'

# remove intermediate files
rm -f $name'_2AT.fa' $name'_2ATC.fa' $name'_2A.fa'

# rename pseudogenome and index files
mv $name'_2ATCG.fa' $name'_pseudogenome.fa'
mv $name'_2ATCG.fa.fai' $name'_pseudogenome.fa.fai'


