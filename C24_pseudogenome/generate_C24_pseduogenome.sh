#!/bin/bash
# Generate C24 pseudogenome and make it Bismark-ready
# EAH 2024


# Get C24 genome files
wget -r -np -nH http://1001genomes.org/data/MPIPZ/MPIPZJiao2020/releases/current/strains/C24/ ./


# 04/08/2024 EAH
# Used python jupyter lab notebook ./make_snp_beds.ipynb to generate SNP beds for each base
# (I wrote a python script to do the same thing from command line in the future)


# make C24 pseudogenome using bedtools maskfasta


bedtools maskfasta -fi TAIR10.fa -bed T10_C24_2A.bed -fo T10_C24_2A.fa -mc A
bedtools maskfasta -fi T10_C24_2A.fa -bed T10_C24_2T.bed -fo T10_C24_2AT.fa -mc T
bedtools maskfasta -fi T10_C24_2AT.fa -bed T10_C24_2C.bed -fo T10_C24_2ATC.fa -mc C
bedtools maskfasta -fi T10_C24_2ATC.fa -bed T10_C24_2G.bed -fo T10_C24_2ATCG.fa -mc G

samtools faidx T10_C24_2ATCG.fa

rm -f T10_C24_2AT.fa T10_C24_2ATC.fa T10_C24_2A.fa

mv T10_C24_2ATCG.fa T10_C24_pseduoG.fa
mv T10_C24_2ATCG.fa.fai T10_C24_pseduoG.fa.fai


# Generate bismark genome

sbatch -p 20 --mem 20G --job-name bismarkprep --output C24pseudo_bismarkprep.out --wrap " bismark_genome_preparation ./C24_pseudo/ "


