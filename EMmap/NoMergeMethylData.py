#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# make methylation data files ready for downstream analysis (no merging)
for arg in sys.argv:
    print(arg)
methyldir_1=sys.argv[1] # where methylextract outputs live
output=sys.argv[2] # sample prefix (ex. Col_1)
methyl1=sys.argv[3] # name of unzipped CX summary file
gnome1=sys.argv[4] # name of genome
outdir=sys.argv[5] # where output files will go (string)
minimum=int(sys.argv[6]) # minimum reads for inclusion in browsing file
#snps=sys.argv[12]

tab="\t"
minstr=str(minimum)

# Bismark CX report header: <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
methylCol=['chr','pos','strand','count_methylated','count_unmethylated','C_context','trinucleotide_context']

methyl1_1_data=pd.read_csv(methyl1, sep=tab, header=None, names=methylCol)

## remove SNPs here, spit out two files one with snps (browsing) and one without (analysis). 
def make_bed (df):
    bed=pd.DataFrame()
    bed['chr']=df['chr']
    bed['start']=df['pos']-1
    bed['end']=df['pos']
    bed['count_unmethylated']=df['count_unmethylated']
    bed['count_methylated']=df['count_methylated']
    bed['methylation_percent']=bed['count_methylated']/(bed['count_methylated']+bed['count_unmethylated'])
    bed=bed.fillna(value=0, axis=1)
    return (bed)

def make_data_bed (df, context):
	df=df[df['chr']!="methylated_CG_pUC19_NEB"]
	df_context=df[df['C_context']==context]

	bed_context= make_bed(df_context)

	return (bed_context)
    
methyl1CG=make_data_bed(methyl1_1_data, 'CG')

methyl1CHG=make_data_bed(methyl1_1_data, 'CHG')
 
methyl1CHH=make_data_bed(methyl1_1_data, 'CHH')


methyl1CG.to_csv(outdir+output+"_"+gnome1+"_CG.bed", header=False, index=False, sep=tab)

methyl1CHG.to_csv(outdir+output+"_"+gnome1+"_CHG.bed", header=False, index=False, sep=tab)

methyl1CHH.to_csv(outdir+output+"_"+gnome1+"_CHH.bed", header=False, index=False, sep=tab)

# Filter coverage file by minimum reads
print("Filter data from file by minimum: "+minstr+" reads minimum")
def coverage_filter (df):
    df_min=df[df['count_methylated'].astype('int64') + df['count_unmethylated'].astype('int64') >= minimum]
    #return (df_min)
    
#methyl1CG_pass=coverage_filter(methyl1CG)

#methyl1CHG_pass=coverage_filter(methyl1CHG)

#methyl1CHH_pass=coverage_filter(methyl1CHH)

# Generate filtered bedgraph and .bed files by context from filtered methyl data

def bedgraph_by_context(df):
    bedgraph_context=pd.DataFrame()
    bedgraph_context['chr']=df['chr']
    bedgraph_context['start']=df['start']
    bedgraph_context['end']=df['end']
    bedgraph_context['methylation_percent']=df['methylation_percent']
    #return (bedgraph_context)
    
#methyl1CG_bg=bedgraph_by_context(methyl1CG_pass)
#methyl1CHG_bg=bedgraph_by_context(methyl1CHG_pass)
#methyl1CHH_bg=bedgraph_by_context(methyl1CHH_pass)

# Save filtered output files to methylextract directory
#methyl1CG_bg.to_csv(outdir+output+'_'+gnome1+'_CG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
#methyl1CHG_bg.to_csv(outdir+output+'_'+gnome1+'_CHG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
#methyl1CHH_bg.to_csv(outdir+output+'_'+gnome1+'_CHH.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
# Get a summary file

def methylfilter (df):
	df_filtered = df[~df['chr'].isin(['ChrC', 'ChrM', 'methylated_CG_pUC19_NEB'])]
	bed=pd.DataFrame()
	bed['chr']=df_filtered['chr']
	bed['start']=df_filtered['pos']-1
	bed['end']=df_filtered['pos']
	bed['count_unmethylated']=df_filtered['count_unmethylated']
	bed['count_methylated']=df_filtered['count_methylated']
	bed['C_context']=df_filtered['C_context']
	#return (bed)
    
def make_methyl_summary (df1, out):
	df=methylfilter(df1)
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
	#return (dfout)   	

#genome1summary=make_methyl_summary(methyl1_1_data, outdir+output+"_"+gnome1+"_mCsummary.txt") 
