#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# summarize methylation from Bismark methylextract CX report by chromosome and context (ex: Chr1 CG%, CHG%, CHH%) 
# and make browsing files for IGV

for arg in sys.argv:
    print(arg)
methyldir=sys.argv[1] # where methylextract outputs live
output=sys.argv[2] # sample prefix (ex. Col_1)
methylout=sys.argv[3] # name of unzipped CX summary file
coverageout=sys.argv[4] # name of unzipped bismark coverage file
minimum=int(sys.argv[5]) # minimum reads for inclusion in browsing file

tab="\t"
minstr=str(minimum)

# Bismark CX report header: <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
methylCol=['chr','pos','strand','count_methylated','count_unmethylated','C_context','trinucleotide_context']
methylinfo=pd.read_csv(methyldir+methylout, sep=tab, header=None, names=methylCol)

# Bismark coverage header: <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
coverageCol=['chr','start','end','methylation_percent','count_methylated','count_unmethylated']
coverageinfo=pd.read_csv(methyldir+coverageout, sep=tab, header=None, names=coverageCol)

# Fraction mC/total C per chromosome by sequence context 
def methylChr (df):
	dfsum=df.groupby(['chr','C_context']).sum()
	dfsum['fraction_methylated']=dfsum['count_methylated']/(dfsum['count_methylated']+dfsum['count_unmethylated'])
	dfout=dfsum.drop(labels='pos', axis=1)
	return (dfout)
        
methyl_summary=methylChr(methylinfo)

# Save summary output dataframe to methylextract directory
methyl_summary.to_csv(methyldir+output+'_mC_summary_by_chr.txt', header=True, index=True, sep=tab)

print("Find summary table here: "+methyldir+output+"_mC_summary_by_chr.txt")


def methylContext (df):
	df_filtered = df[~df['chr'].isin(['ChrC', 'ChrM', 'methylated_CG_pUC19_NEB'])]
	dfsum=df_filtered.groupby(['C_context']).sum()
	dfsum['fraction_methylated']=dfsum['count_methylated']/(dfsum['count_methylated']+dfsum['count_unmethylated'])
	dfout=dfsum.drop(labels='pos', axis=1)
	return (dfout)
    
total_summary=methylContext(methylinfo)
    
# Save summary output dataframe to methylextract directory
methyl_summary.to_csv(methyldir+output+'_mC_summary_by_context.txt', header=True, index=True, sep=tab)

print("Total methylation (excluding ChrC, ChrM, or pUC19) info summarized for "+methylout+".")
print("Find summary table here: "+methyldir+output+"_mC_summary_by_context.txt")
    
# Merge unfiltered coverage and methyl summary files to incorporate coverage filter with methyl call info
unfilteredinfo=coverageinfo.merge(methylinfo, how='left', left_on=['chr','end','count_methylated','count_unmethylated'], 
                   right_on=['chr','pos','count_methylated','count_unmethylated'])

# Generate unfiltered .bed file for DSS dmrs and other analysis
def bed_by_context(df, context):
	df_context=df[df['C_context']==context]
	bed_context=pd.DataFrame()
	bed_context['chr']=df_context['chr']
	bed_context['start']=df_context['start']
	bed_context['end']=df_context['end']
	bed_context['count_unmethylated']=df_context['count_unmethylated']
	bed_context['count_methylated']=df_context['count_methylated']
	bed_context['methylation_percent']=df_context['methylation_percent']
	return (bed_context)

CG_allbed=bed_by_context(unfilteredinfo, 'CG')
CHG_allbed=bed_by_context(unfilteredinfo, 'CHG')
CHH_allbed=bed_by_context(unfilteredinfo, 'CHH')

# Save unfiltered bed files
CG_allbed.to_csv(methyldir+output+'_CG.bed', header=False, index=False, sep=tab)
CHG_allbed.to_csv(methyldir+output+'_CHG.bed', header=False, index=False, sep=tab)
CHH_allbed.to_csv(methyldir+output+'_CHH.bed', header=False, index=False, sep=tab)

#print("Unfiltered methylation info by context for "+methylout+" found here: "+methyldir+output+"_Cx.bed")

# Filter coverage file by minimum reads
print("Filter data from coverage file by minimum: "+minstr+" reads minimum")
def coverage_filter (df):
	df_min=df[df['count_methylated'].astype('int64') + df['count_unmethylated'].astype('int64') >= minimum]
	return (df_min)
    
coveragepass=coverage_filter(unfilteredinfo)

# Generate filtered bedgraph and .bed files by context from filtered methyl data

def bedgraph_by_context(df, context):
	df_context=df[df['C_context']==context]
	bedgraph_context=pd.DataFrame()
	bedgraph_context['chr']=df_context['chr']
	bedgraph_context['start']=df_context['start']
	bedgraph_context['end']=df_context['end']
	bedgraph_context['methylation_percent']=df_context['methylation_percent']
	return (bedgraph_context)
    
CG_bg=bedgraph_by_context(coveragepass, 'CG')
CHG_bg=bedgraph_by_context(coveragepass, 'CHG')
CHH_bg=bedgraph_by_context(coveragepass, 'CHH')

CG_bed=bed_by_context(coveragepass, 'CG')
CHG_bed=bed_by_context(coveragepass, 'CHG')
CHH_bed=bed_by_context(coveragepass, 'CHH')

# Save filtered output files to methylextract directory
CG_bg.to_csv(methyldir+output+'_CG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
CHG_bg.to_csv(methyldir+output+'_CHG.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)
CHH_bg.to_csv(methyldir+output+'_CHH.min'+minstr+'.bedGraph', header=False, index=False, sep=tab)

CG_bed.to_csv(methyldir+output+'_CG.min'+minstr+'.bed', header=False, index=False, sep=tab)
CHG_bed.to_csv(methyldir+output+'_CHG.min'+minstr+'.bed', header=False, index=False, sep=tab)
CHH_bed.to_csv(methyldir+output+'_CHH.min'+minstr+'.bed', header=False, index=False, sep=tab)

print("Filtered methylation info by context for "+methylout+" found here: "+methyldir+output+"_Cx.min"+minstr+".bed")

