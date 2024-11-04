#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# combine methylation data for a genome from both mapping instances and generate filtered bedgraph files for browsing

for arg in sys.argv:
    print(arg)
methyldir_1=sys.argv[1] # where methylextract outputs live
methyldir_2=sys.argv[2] # where methylextract outputs live
output=sys.argv[3] # sample prefix (ex. Col_1)
methyl1_1=sys.argv[4] # name of unzipped CX summary file
methyl1_2=sys.argv[5] # name of unzipped CX summary file
methyl2_1=sys.argv[6] # name of unzipped CX summary file
methyl2_2=sys.argv[7] # name of unzipped CX summary file
gnome1=sys.argv[8] # name of genome1
gnome2=sys.argv[9] # name of genome2
outdir=sys.argv[10] # where output files will go (string)
minimum=int(sys.argv[11]) # minimum reads for inclusion in browsing file
#snps=sys.argv[12]

tab="\t"
minstr=str(minimum)

# Bismark CX report header: <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
methylCol=['chr','pos','strand','count_methylated','count_unmethylated','C_context','trinucleotide_context']

methyl1_1_data=pd.read_csv(methyl1_1, sep=tab, header=None, names=methylCol)
methyl1_2_data=pd.read_csv(methyl1_2, sep=tab, header=None, names=methylCol)
methyl2_1_data=pd.read_csv(methyl2_1, sep=tab, header=None, names=methylCol)
methyl2_2_data=pd.read_csv(methyl2_2, sep=tab, header=None, names=methylCol)

# SNPs dataframe: [chr, start, end, ref>alt]
snpscol=['chr','start','end','ref>alt']
#SNPs=pd.read_csv(snps, sep=tab, header=None, names=snpscol)

## remove SNPs here, spit out two files one with snps (browsing) and one without (analysis). 
def make_bed (df):
    bed=pd.DataFrame()
    bed['chr']=df['chr']
    bed['start']=df['pos']-1
    bed['end']=df['pos']
    bed['count_unmethylated']=df['count_unmethylated_x']+df['count_unmethylated_y']
    bed['count_methylated']=df['count_methylated_x']+df['count_methylated_y']
    bed['methylation_percent']=bed['count_methylated']/(bed['count_methylated']+bed['count_unmethylated'])
    bed=bed.fillna(value=0, axis=1)
    return (bed)

def make_merged_data_bed (df1, df2, context):
    
	df=df1.merge(df2, on=['chr','pos','strand','C_context','trinucleotide_context'], how="outer")
	df=df[df['chr']!="methylated_CG_pUC19_NEB"]
	df_context=df[df['C_context']==context]
	#df_context_nosnps=df_context[~df_context['chr','pos'].isin(SNPs['chr','end'])]

	bed_context= make_bed(df_context)
	#bed_context_nosnps= make_bed(df_context_nosnps)

	return (bed_context)
    
methyl1CG=make_merged_data_bed(methyl1_1_data, methyl1_2_data, 'CG')
methyl2CG =make_merged_data_bed(methyl2_1_data, methyl2_2_data, 'CG')

methyl1CHG=make_merged_data_bed(methyl1_1_data, methyl1_2_data, 'CHG')
methyl2CHG=make_merged_data_bed(methyl2_1_data, methyl2_2_data, 'CHG')
 
methyl1CHH=make_merged_data_bed(methyl1_1_data, methyl1_2_data, 'CHH')
methyl2CHH=make_merged_data_bed(methyl2_1_data, methyl2_2_data, 'CHH')


methyl1CG.to_csv(outdir+output+"_"+gnome1+"_CG.bed", header=False, index=False, sep=tab)
methyl2CG.to_csv(outdir+output+"_"+gnome2+"_CG.bed", header=False, index=False, sep=tab)

methyl1CHG.to_csv(outdir+output+"_"+gnome1+"_CHG.bed", header=False, index=False, sep=tab)
methyl2CHG.to_csv(outdir+output+"_"+gnome2+"_CHG.bed", header=False, index=False, sep=tab)

methyl1CHH.to_csv(outdir+output+"_"+gnome1+"_CHH.bed", header=False, index=False, sep=tab)
methyl2CHH.to_csv(outdir+output+"_"+gnome2+"_CHH.bed", header=False, index=False, sep=tab)

# Filter coverage file by minimum reads
print("Filter data from file by minimum: "+minstr+" reads minimum")
def coverage_filter (df):
    df_min=df[df['count_methylated'].astype('int64') + df['count_unmethylated'].astype('int64') >= minimum]
    return (df_min)
    
methyl1CG_pass=coverage_filter(methyl1CG)
methyl2CG_pass=coverage_filter(methyl2CG)

methyl1CHG_pass=coverage_filter(methyl1CHG)
methyl2CHG_pass=coverage_filter(methyl2CHG)

methyl1CHH_pass=coverage_filter(methyl1CHH)
methyl2CHH_pass=coverage_filter(methyl2CHH)


# Generate filtered bedgraph and .bed files by context from filtered methyl data

def bedgraph_by_context(df):
    bedgraph_context=pd.DataFrame()
    bedgraph_context['chr']=df['chr']
    bedgraph_context['start']=df['start']
    bedgraph_context['end']=df['end']
    bedgraph_context['methylation_percent']=df['methylation_percent']
    return (bedgraph_context)
    
methyl1CG_bg=bedgraph_by_context(methyl1CG_pass)
methyl1CHG_bg=bedgraph_by_context(methyl1CHG_pass)
methyl1CHH_bg=bedgraph_by_context(methyl1CHH_pass)

methyl2CG_bg=bedgraph_by_context(methyl2CG_pass)
methyl2CHG_bg=bedgraph_by_context(methyl2CHG_pass)
methyl2CHH_bg=bedgraph_by_context(methyl2CHH_pass)

# Save filtered output files to methylextract directory
methyl1CG_bg.to_csv(outdir+output+'_'+gnome1+'_CG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
methyl1CHG_bg.to_csv(outdir+output+'_'+gnome1+'_CHG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
methyl1CHH_bg.to_csv(outdir+output+'_'+gnome1+'_CHH.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)

methyl2CG_bg.to_csv(outdir+output+'_'+gnome2+'_CG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
methyl2CHG_bg.to_csv(outdir+output+'_'+gnome2+'_CHG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
methyl2CHH_bg.to_csv(outdir+output+'_'+gnome2+'_CHH.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)

# Get a summary file

def methylfilter (df1, df2):
	df=df1.merge(df2, on=['chr','pos','strand','C_context','trinucleotide_context'], how="outer")
	df_filtered = df[~df['chr'].isin(['ChrC', 'ChrM', 'methylated_CG_pUC19_NEB'])]
	bed=pd.DataFrame()
	bed['chr']=df_filtered['chr']
	bed['start']=df_filtered['pos']-1
	bed['end']=df_filtered['pos']
	bed['count_unmethylated']=df_filtered['count_unmethylated_x']+df_filtered['count_unmethylated_y']
	bed['count_methylated']=df_filtered['count_methylated_x']+df_filtered['count_methylated_y']
	bed['C_context']=df_filtered['C_context']
	return (bed)
    
def make_methyl_summary (df1, df2, out):
	df=methylfilter(df1, df2)
	dfsize=df.groupby(['C_context']).count()
	dfsum=df.groupby(['C_context']).sum()
	##dfsum=dfsum.join(dfsize)
	dfsum['fraction_methylated']=dfsum['count_methylated']/(dfsum['count_methylated']+dfsum['count_unmethylated'])
	##dfsum['n_count/total']=(dfsum['count_methylated']+dfsum['count_unmethylated'])/dfsum['counts']
	dfsum['n_count']=(dfsum['count_methylated']+dfsum['count_unmethylated'])
	dfout=dfsum.drop(labels=['start','end','count_methylated','count_unmethylated'], axis=1)
	dfout.to_csv(out, header=True, index=True, sep=tab)
	print("length "+out+" = ")
	print(dfsize)
	return (dfout)   	

genome1summary=make_methyl_summary(methyl1_1_data, methyl1_2_data, outdir+output+"_"+gnome1+"_mCsummary.txt") 
genome2summary=make_methyl_summary(methyl2_1_data, methyl2_2_data, outdir+output+"_"+gnome2+"_mCsummary.txt") 


bigsum=pd.concat([genome1summary, genome2summary], keys=[gnome1, gnome2])
bigsum.to_csv(outdir+output+"_genome_mC_summary.txt", header=True, index=True, sep=tab)

