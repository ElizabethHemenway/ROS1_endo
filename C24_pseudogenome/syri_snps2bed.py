#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# To generate Col-C24 pseudogenome file using syri TSV output file, start by making a bed file for each SNP base
# Usage: syri_snps2bed.py syri_out_file 

for arg in sys.argv:
    print(arg)
syri_out=sys.argv[1] # syri file (if you hav a snp bed file, change columns=TSV to columns=bed in line 22)
ref_name=sys.argv[2] # name of reference genome (ex=T10)
alt_name=sys.argv[3] # name of alternate genome (ex=C24)
output_dir=sys.argv[4] # location of output directory for saving snp bed files

tab="\t"

TSV=['Chr','Start','End','Ref','Alt','ChrAlt','StartAlt','EndAlt','ID','Type','Variant','Copy']
bed=['Chr','Start','End','Ref','Alt']
tab="\t"

columns=TSV ##Change to columns=bed if you want to use a snp bed file as input.

syri=pd.read_csv(syri_out, header=None, sep=tab, names=columns)

snps=syri[syri['Variant']=="SNP"]

# Generate bed files with individual SNP positions to replace with each nucleotide

snps_bed=pd.DataFrame(columns=bed)

snps_bed['Chr']=snps['Chr']
snps_bed['Start']=snps['Start']
snps_bed['End']=snps['End']
snps_bed['Ref']=snps['Ref']
snps_bed['Alt']=snps['Alt']

snps2A=snps_bed[snps_bed['Alt']=='A']
snps2T=snps_bed[snps_bed['Alt']=='T']
snps2C=snps_bed[snps_bed['Alt']=='C']
snps2G=snps_bed[snps_bed['Alt']=='G']

snps2A.to_csv(output_dir+ref_name+alt_name+'2a.bed', header=False, index=False, sep=tab)
snps2T.to_csv(output_dir+ref_name+alt_name+'2T.bed', header=False, index=False, sep=tab)
snps2C.to_csv(output_dir+ref_name+alt_name+'2C.bed', header=False, index=False, sep=tab)
snps2G.to_csv(output_dir+ref_name+alt_name+'2G.bed', header=False, index=False, sep=tab)
