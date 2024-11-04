#!/usr/bin/env python
##EAH 07/31/24
## Remove Ws regions from ros1-3 and rdd methylation data files
##Usage: remove_Ws.py inputfile outputfile (r3 or rdd)

import pandas as pd
import numpy as np
import sys

for arg in sys.argv:
    print(arg)
infile=sys.argv[1] # input file
output=sys.argv[2] # prefix for output file (ex. path...r3_1_leaf)
genotype=sys.argv[3] # remove just Chr2 region (r3) or Chr2 and Chr3 region (rdd)

bedcol=['chr','start','end','count_unmethylated','count_methylated','methylation_percent']
tab="\t"
		
data=pd.read_csv(infile, header=None, sep=tab, names=bedcol)

def exclude_ws_r3 (bed):
    
    df_chr1 = pd.DataFrame(bed[bed['chr']=='Chr1'])
    df_chr2 = pd.DataFrame(bed[bed['chr']=='Chr2'])
    df_chr3 = pd.DataFrame(bed[bed['chr']=='Chr3'])
    df_chr4 = pd.DataFrame(bed[bed['chr']=='Chr4'])
    df_chr5 = pd.DataFrame(bed[bed['chr']=='Chr5'])
    df_chrC = pd.DataFrame(bed[bed['chr']=='ChrC'])
    df_chrM = pd.DataFrame(bed[bed['chr']=='ChrM'])

    chr2Col = pd.DataFrame(df_chr2[(df_chr2['start']<=8802496) | (df_chr2['start']>=15397296)])
    
    df_Col0 = pd.DataFrame(columns=bedcol)
    df_Col0 =df_Col0.append(df_chr1)
    df_Col0 =df_Col0.append(chr2Col)
    df_Col0 =df_Col0.append(df_chr3)
    df_Col0 =df_Col0.append(df_chr4)
    df_Col0 =df_Col0.append(df_chr5)
    df_Col0 =df_Col0.append(df_chrC)
    df_Col0 =df_Col0.append(df_chrM)
    
    return(df_Col0)
        
def exclude_ws_rdd (bed):
    
    df_chr1 = pd.DataFrame(bed[bed['chr']=='Chr1'])
    df_chr2 = pd.DataFrame(bed[bed['chr']=='Chr2'])
    df_chr3 = pd.DataFrame(bed[bed['chr']=='Chr3'])
    df_chr4 = pd.DataFrame(bed[bed['chr']=='Chr4'])
    df_chr5 = pd.DataFrame(bed[bed['chr']=='Chr5'])
    df_chrC = pd.DataFrame(bed[bed['chr']=='ChrC'])
    df_chrM = pd.DataFrame(bed[bed['chr']=='ChrM'])

    chr2Col = pd.DataFrame(df_chr2[(df_chr2['start']<=8802496) | (df_chr2['start']>=15397296)])
    chr3Col = pd.DataFrame(df_chr3[(df_chr3['start']<=677340) | (df_chr3['start']>=5117803)])
    
    df_Col0 = pd.DataFrame(columns=bedcol)
    df_Col0 =df_Col0.append(df_chr1)
    df_Col0 =df_Col0.append(chr2Col)
    df_Col0 =df_Col0.append(chr3Col)
    df_Col0 =df_Col0.append(df_chr4)
    df_Col0 =df_Col0.append(df_chr5)
    df_Col0 =df_Col0.append(df_chrC)
    df_Col0 =df_Col0.append(df_chrM)
    
    return(df_Col0)
    
if genotype=="r3":
	print("Excluding coordinates chr2: 8,802,496-15,397,296")
	outfile=output+"_exWs_Chr2.bed"
	data_filt=exclude_ws_r3(data)
	
if genotype=="rdd":
	print("Excluding coordinates chr2: 8,802,496-15,397,296 and chr3: 677,340-5,117,803")
	outfile=output+"_exWs_Chr2Chr3.bed"
	data_filt=exclude_ws_rdd(data)

print("Saving filtered data file here: "+outfile)	
data_filt.to_csv(outfile, sep=tab, header=False, index=False)


  
