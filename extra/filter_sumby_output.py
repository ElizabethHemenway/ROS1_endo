#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

# to filter sumby dataframes for threshold mC
# ./filter_sumby_output.py /lab/solexa_gehring/elizabeth/TE_spreading/leaf_spreading/sumby/ wt_1_leaf_sumby.bed 0.1 ">"or"<"
for arg in sys.argv:
    print(arg)
file=sys.argv[1] # file to process
name=sys.argv[2]
thresh=float(sys.argv[3]) # filter by this as threshold fraction mC
threshstr=str(thresh)
direction=sys.argv[4] #return greater than or less than threshold

col=['chr','start','end','feature','frac_mC','numC'] #'strand', 'DMR_chr','DMR_start','DMR_end']
tab="\t"

df=pd.read_csv(file, sep=tab, header=None, names=col)

if direction == ">":
	pass_df=df[df['frac_mC']>=thresh]
	
	pass_df.to_csv(name+"_greaterthanequalto_"+threshstr+".bed", sep=tab, header=None, index=None)
	print("saved output file as: "+name+"_greaterthanequalto_"+threshstr+".bed")
elif direction == "<":
	pass_df=df[df['frac_mC']<=thresh]
	
	pass_df.to_csv(name+"_lessthanequalto_"+threshstr+".bed", sep=tab, header=None, index=None)
	print("saved output file as: "+name+"_lessthanequalto_"+threshstr+".bed")

print("done!")