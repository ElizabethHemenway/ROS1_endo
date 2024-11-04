#!/usr/bin/env python
# coding: utf-8
##EAH 04/12/2024
##Filter final output bedfiles before doing analyses. 

import pandas as pd
import numpy as np
import sys

for arg in sys.argv:
    print(arg)
methyldir=sys.argv[1] # where inputs and outputs should live (for EMmap sampledir/merge_a2a)
output=sys.argv[2] # sample prefix (ex. Col_1)
gnome1=sys.argv[3] # genome name
gnome2=sys.argv[4] # genome name
minimum=int(sys.argv[5]) # minimum reads for inclusion in browsing file
outdir=sys.argv[6] # where output files will go (string)
SNPs=sys.argv[7]

bedcol=['chr','start','end','count_unmethylated','count_methylated','methylation_percent']
sumcolumns=['chr','start','end','name','avg_methy','nC']
tab="\t"
minstr=str(minimum)

gnome1_CG=output+"_"+gnome1+"_CG.bed"
gnome2_CG=output+"_"+gnome2+"_CG.bed"
gnome1_CHG=output+"_"+gnome1+"_CHG.bed"
gnome2_CHG=output+"_"+gnome2+"_CHG.bed"
gnome1_CHH=output+"_"+gnome1+"_CHH.bed"
gnome2_CHH=output+"_"+gnome2+"_CHH.bed"


snps=pd.read_csv(SNPs, header=None, sep=tab, names=['chr','start','end','snp'])


gnome1_CG_data=pd.read_csv(methyldir+gnome1_CG, header=None, sep=tab, names=bedcol)
gnome2_CG_data=pd.read_csv(methyldir+gnome2_CG, header=None, sep=tab, names=bedcol)

gnome1_CHG_data=pd.read_csv(methyldir+gnome1_CHG, header=None, sep=tab, names=bedcol)
gnome2_CHG_data=pd.read_csv(methyldir+gnome2_CHG, header=None, sep=tab, names=bedcol)

gnome1_CHH_data=pd.read_csv(methyldir+gnome1_CHH, header=None, sep=tab, names=bedcol)
gnome2_CHH_data=pd.read_csv(methyldir+gnome2_CHH, header=None, sep=tab, names=bedcol)


# Filter coverage file by minimum reads
def coverage_filter (df):
    print("Filter data from "+df+" by minimum: "+minstr+" reads minimum")
    df_min=df[df['count_methylated'].astype('int64') + df['count_unmethylated'].astype('int64') >= minimum]
    return (df_min)


gnome1_CG_min5=coverage_filter(gnome1_CG_data)
gnome2_CG_min5=coverage_filter(gnome2_CG_data)

gnome1_CHG_min5=coverage_filter(gnome1_CHG_data)
gnome2_CHG_min5=coverage_filter(gnome2_CHG_data)

gnome1_CHH_min5=coverage_filter(gnome1_CHH_data)
gnome2_CHH_min5=coverage_filter(gnome2_CHH_data)


#Are there min reads on both genomes of a sample?
def coverage_on_both_genomes(g1_df, g2_df):
    g1_context_merge=g1_df.merge(g2_df, on=['chr','start','end'], how="left", suffixes=('_G1','_G2'), indicator=True)
    g1_context_merge_pass=g1_context_merge[g1_context_merge['_merge']=='both']
    g1_context_merge_pass=g1_context_merge_pass.drop(['count_unmethylated_G2','count_methylated_G2','methylation_percent_G2','_merge'], axis=1)
    g1_context_merge_pass=g1_context_merge_pass.set_axis(['chr','start','end','count_unmethylated','count_methylated','methylation_percent'], axis=1)
    
    g2_context_merge=g2_df.merge(g1_df, on=['chr','start','end'], how="left", suffixes=('_G2','_G1'), indicator=True)
    g2_context_merge_pass=g2_context_merge[g2_context_merge['_merge']=='both']
    g2_context_merge_pass=g2_context_merge_pass.drop(['count_unmethylated_G1','count_methylated_G1','methylation_percent_G1','_merge'], axis=1)
    g2_context_merge_pass=g2_context_merge_pass.set_axis(['chr','start','end','count_unmethylated','count_methylated','methylation_percent'], axis=1)

    return(g1_context_merge_pass, g2_context_merge_pass)
    
gnome1_CG_min5_covpass, gnome2_CG_min5_covpass=coverage_on_both_genomes(gnome1_CG_min5, gnome2_CG_min5)
gnome1_CHG_min5_covpass, gnome2_CHG_min5_covpass=coverage_on_both_genomes(gnome1_CHG_min5, gnome2_CHG_min5)
gnome1_CHH_min5_covpass, gnome2_CHH_min5_covpass=coverage_on_both_genomes(gnome1_CHH_min5, gnome2_CHH_min5)

# Remove SNPs from data files
def remove_snps(df):
    df_snps_merge=df.merge(snps, on=['chr','start','end'], how="left", indicator=True)
    not_snp=df_snps_merge[df_snps_merge['_merge']=='left_only']
    not_snp=not_snp.drop(['snp','_merge'], axis=1)
    return(not_snp)

gnome1_CG_min5_covpass_notsnp=remove_snps(gnome1_CG_min5_covpass)
gnome2_CG_min5_covpass_notsnp=remove_snps(gnome2_CG_min5_covpass)

gnome1_CHG_min5_covpass_notsnp=remove_snps(gnome1_CHG_min5_covpass)
gnome2_CHG_min5_covpass_notsnp=remove_snps(gnome2_CHG_min5_covpass)

gnome1_CHH_min5_covpass_notsnp=remove_snps(gnome1_CHH_min5_covpass)
gnome2_CHH_min5_covpass_notsnp=remove_snps(gnome2_CHH_min5_covpass)

gnome1_CG_min5_covpass_notsnp.to_csv(methyldir+output+"_"+gnome1+"_CG_min5inbothgnome_noSNPs.bed", header=False, index=False, sep=tab)
gnome2_CG_min5_covpass_notsnp.to_csv(methyldir+output+"_"+gnome2+"_CG_min5inbothgnome_noSNPs.bed", header=False, index=False, sep=tab)

gnome1_CHG_min5_covpass_notsnp.to_csv(methyldir+output+"_"+gnome1+"_CHG_min5inbothgnome_noSNPs.bed", header=False, index=False, sep=tab)
gnome2_CHG_min5_covpass_notsnp.to_csv(methyldir+output+"_"+gnome2+"_CHG_min5inbothgnome_noSNPs.bed", header=False, index=False, sep=tab)

gnome1_CHH_min5_covpass_notsnp.to_csv(methyldir+output+"_"+gnome1+"_CHH_min5inbothgnome_noSNPs.bed", header=False, index=False, sep=tab)
gnome2_CHH_min5_covpass_notsnp.to_csv(methyldir+output+"_"+gnome2+"_CHH_min5inbothgnome_noSNPs.bed", header=False, index=False, sep=tab)

print("filtered methyl data files saved to "+methyldir)


