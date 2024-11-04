#!/bin/bash

###EAH 08/06/2024


cd /lab/solexa_gehring/elizabeth/demethylase_em_remap2024/DMRs

rm -f *.sorted.bed
rm -f *.hyper*
rm -f *.hypo*
rm -f *allC.merge.bed

/lab/solexa_gehring/elizabeth/scripts_and_code/EAH_scripts/hyper_hypo.sh
wait

makeallC () {
	CG=$1
	CHG=$2
	CHH=$3
	name=$4
	allC=$name"_allC"
	cat $CG $CHG $CHH > $allC".cat.bed"
	sort -k1,1 -k2,2n $allC".cat.bed" > $allC".sorted.bed"
	bedtools merge -i $allC".sorted.bed" > $allC".merge.bed"
		
}

sortandmerge () {
	filename=$1".bed"
	sort -k1,1 -k2,2n $filename > $1".sorted.bed"
	bedtools merge -i $1".sorted.bed" > $1".merge.bed"
		
}


makeallC r7_v_wt_leaf_CG.DMRs.hyper.bed r7_v_wt_leaf_CHG.DMRs.hyper.bed r7_v_wt_leaf_CHH.DMRs.hyper.bed r7_v_wt_leaf_allC_hyper
makeallC r3_v_wt_leaf_CG.DMRs.hyper.bed r3_v_wt_leaf_CHG.DMRs.hyper.bed r3_v_wt_leaf_CHH.DMRs.hyper.bed r3_v_wt_leaf_allC_hyper
makeallC r7_v_wt_CG.DMRs.hyper.bed r7_v_wt_CHG.DMRs.hyper.bed r7_v_wt_CHH.DMRs.hyper.bed r7_v_wt_allC_hyper
makeallC r3_v_wt_CG.DMRs.hyper.bed r3_v_wt_CHG.DMRs.hyper.bed r3_v_wt_CHH.DMRs.hyper.bed r3_v_wt_allC_hyper
makeallC rdd_v_wt_CG.DMRs.hyper.bed rdd_v_wt_CHG.DMRs.hyper.bed rdd_v_wt_CHH.DMRs.hyper.bed rdd_v_wt_allC_hyper


rm -f *.sorted.bed
rm -f *cat.bed
