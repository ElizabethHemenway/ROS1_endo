#!/bin/bash

###EAH 4/21/24 bedtools intersect -u -f 0.1 -a


cd /lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/

rm -f recip*
rm -f *.sorted.bed
rm -r -f processed 
rm -f *.hyper*
rm -f *.hypo*

/lab/solexa_gehring/elizabeth/scripts_and_code/EAH_scripts/hyper_hypo.sh
wait

# sort bed files
FILES=$(find . -type f -name "*.bed") 
for file in $FILES ; do
	name="${file%.*}"
	sort -k1,1 -k2,2n $file > $name".sorted.bed"
done

# intersect -u -f 0.1 -a mat v pat DMRs from both cross directions. (just keeping a,  I really don't think it matters.)
bedtools intersect -u -f 0.1 -a C24xCol_mat_v_pat_CG.DMRs.hyper.sorted.bed -b ColxC24_mat_v_pat_CG.DMRs.hyper.sorted.bed > recipA_WT_matVpat_CGr.bed
bedtools intersect -u -f 0.1 -a ColxC24_mat_v_pat_CG.DMRs.hyper.sorted.bed -b C24xCol_mat_v_pat_CG.DMRs.hyper.sorted.bed > recipB_WT_matVpat_CGr.bed
cat recipA_WT_matVpat_CGr.bed recipB_WT_matVpat_CGr.bed > recipCat_WT_matVpat_CGr.bed

bedtools intersect -u -f 0.1 -a C24xCol_mat_v_pat_CG.DMRs.hypo.sorted.bed -b ColxC24_mat_v_pat_CG.DMRs.hypo.sorted.bed > recipA_WT_matVpat_CGo.bed
bedtools intersect -u -f 0.1 -a ColxC24_mat_v_pat_CG.DMRs.hypo.sorted.bed -b C24xCol_mat_v_pat_CG.DMRs.hypo.sorted.bed > recipB_WT_matVpat_CGo.bed
cat recipA_WT_matVpat_CGo.bed recipB_WT_matVpat_CGo.bed > recipCat_WT_matVpat_CGo.bed

bedtools intersect -u -f 0.1 -a r1xr3_mat_v_pat_CG.DMRs.hyper.sorted.bed -b r3xr1_mat_v_pat_CG.DMRs.hyper.sorted.bed > recipA_ros1_matVpat_CGr.bed
bedtools intersect -u -f 0.1 -a r3xr1_mat_v_pat_CG.DMRs.hyper.sorted.bed -b r1xr3_mat_v_pat_CG.DMRs.hyper.sorted.bed > recipB_ros1_matVpat_CGr.bed
cat recipA_ros1_matVpat_CGr.bed recipB_ros1_matVpat_CGr.bed > recipCat_ros1_matVpat_CGr.bed

bedtools intersect -u -f 0.1 -a r1xr3_mat_v_pat_CG.DMRs.hypo.sorted.bed -b r3xr1_mat_v_pat_CG.DMRs.hypo.sorted.bed > recipA_ros1_matVpat_CGo.bed
bedtools intersect -u -f 0.1 -a r3xr1_mat_v_pat_CG.DMRs.hypo.sorted.bed -b r1xr3_mat_v_pat_CG.DMRs.hypo.sorted.bed > recipB_ros1_matVpat_CGo.bed
cat recipA_ros1_matVpat_CGo.bed recipB_ros1_matVpat_CGo.bed > recipCat_ros1_matVpat_CGo.bed

bedtools intersect -u -f 0.1 -a C24xCol_mat_v_pat_CHG.DMRs.hyper.sorted.bed -b ColxC24_mat_v_pat_CHG.DMRs.hyper.sorted.bed > recipA_WT_matVpat_CHGr.bed
bedtools intersect -u -f 0.1 -a ColxC24_mat_v_pat_CHG.DMRs.hyper.sorted.bed -b C24xCol_mat_v_pat_CHG.DMRs.hyper.sorted.bed > recipB_WT_matVpat_CHGr.bed
cat recipA_WT_matVpat_CHGr.bed recipB_WT_matVpat_CHGr.bed > recipCat_WT_matVpat_CHGr.bed

bedtools intersect -u -f 0.1 -a C24xCol_mat_v_pat_CHG.DMRs.hypo.sorted.bed -b ColxC24_mat_v_pat_CHG.DMRs.hypo.sorted.bed > recipA_WT_matVpat_CHGo.bed
bedtools intersect -u -f 0.1 -a ColxC24_mat_v_pat_CHG.DMRs.hypo.sorted.bed -b C24xCol_mat_v_pat_CHG.DMRs.hypo.sorted.bed > recipB_WT_matVpat_CHGo.bed
cat recipA_WT_matVpat_CHGo.bed recipB_WT_matVpat_CHGo.bed > recipCat_WT_matVpat_CHGo.bed

bedtools intersect -u -f 0.1 -a r1xr3_mat_v_pat_CHG.DMRs.hyper.sorted.bed -b r3xr1_mat_v_pat_CHG.DMRs.hyper.sorted.bed > recipA_ros1_matVpat_CHGr.bed
bedtools intersect -u -f 0.1 -a r3xr1_mat_v_pat_CHG.DMRs.hyper.sorted.bed -b r1xr3_mat_v_pat_CHG.DMRs.hyper.sorted.bed > recipB_ros1_matVpat_CHGr.bed
cat recipA_ros1_matVpat_CHGr.bed recipB_ros1_matVpat_CHGr.bed > recipCat_ros1_matVpat_CHGr.bed

bedtools intersect -u -f 0.1 -a r1xr3_mat_v_pat_CHG.DMRs.hypo.sorted.bed -b r3xr1_mat_v_pat_CHG.DMRs.hypo.sorted.bed > recipA_ros1_matVpat_CHGo.bed
bedtools intersect -u -f 0.1 -a r3xr1_mat_v_pat_CHG.DMRs.hypo.sorted.bed -b r1xr3_mat_v_pat_CHG.DMRs.hypo.sorted.bed > recipB_ros1_matVpat_CHGo.bed
cat recipA_ros1_matVpat_CHGo.bed recipB_ros1_matVpat_CHGo.bed > recipCat_ros1_matVpat_CHGo.bed


bedtools intersect -u -f 0.1 -a C24xCol_mat_v_pat_CHH.DMRs.hyper.sorted.bed -b ColxC24_mat_v_pat_CHH.DMRs.hyper.sorted.bed > recipA_WT_matVpat_CHHr.bed
bedtools intersect -u -f 0.1 -a ColxC24_mat_v_pat_CHH.DMRs.hyper.sorted.bed -b C24xCol_mat_v_pat_CHH.DMRs.hyper.sorted.bed > recipB_WT_matVpat_CHHr.bed
cat recipA_WT_matVpat_CHHr.bed recipB_WT_matVpat_CHHr.bed > recipCat_WT_matVpat_CHHr.bed

bedtools intersect -u -f 0.1 -a C24xCol_mat_v_pat_CHH.DMRs.hypo.sorted.bed -b ColxC24_mat_v_pat_CHH.DMRs.hypo.sorted.bed > recipA_WT_matVpat_CHHo.bed
bedtools intersect -u -f 0.1 -a ColxC24_mat_v_pat_CHH.DMRs.hypo.sorted.bed -b C24xCol_mat_v_pat_CHH.DMRs.hypo.sorted.bed > recipB_WT_matVpat_CHHo.bed
cat recipA_WT_matVpat_CHHo.bed recipB_WT_matVpat_CHHo.bed > recipCat_WT_matVpat_CHHo.bed


bedtools intersect -u -f 0.1 -a r1xr3_mat_v_pat_CHH.DMRs.hyper.sorted.bed -b r3xr1_mat_v_pat_CHH.DMRs.hyper.sorted.bed > recipA_ros1_matVpat_CHHr.bed
bedtools intersect -u -f 0.1 -a r3xr1_mat_v_pat_CHH.DMRs.hyper.sorted.bed -b r1xr3_mat_v_pat_CHH.DMRs.hyper.sorted.bed > recipB_ros1_matVpat_CHHr.bed
cat recipA_ros1_matVpat_CHHr.bed recipB_ros1_matVpat_CHHr.bed > recipCat_ros1_matVpat_CHHr.bed


bedtools intersect -u -f 0.1 -a r1xr3_mat_v_pat_CHH.DMRs.hypo.sorted.bed -b r3xr1_mat_v_pat_CHH.DMRs.hypo.sorted.bed > recipA_ros1_matVpat_CHHo.bed
bedtools intersect -u -f 0.1 -a r3xr1_mat_v_pat_CHH.DMRs.hypo.sorted.bed -b r1xr3_mat_v_pat_CHH.DMRs.hypo.sorted.bed > recipB_ros1_matVpat_CHHo.bed
cat recipA_ros1_matVpat_CHHo.bed recipB_ros1_matVpat_CHHo.bed > recipCat_ros1_matVpat_CHHo.bed


# intersect -u -f 0.1 -a ros1 v wt by genome from both cross directions. 
bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_pat_Col_spiked_CG.DMRs.hyper.sorted.bed -b r3xr1_v_ColxC24_pat_C24_pseudo_CG.DMRs.hyper.sorted.bed > recipA_r1Vwt_pat_CGr.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_pat_C24_pseudo_CG.DMRs.hyper.sorted.bed -b r1xr3_v_C24xCol_pat_Col_spiked_CG.DMRs.hyper.sorted.bed > recipB_r1Vwt_pat_CGr.bed
cat recipA_r1Vwt_pat_CGr.bed recipB_r1Vwt_pat_CGr.bed > recipCat_ros1_v_wt_pat_CGr.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_pat_Col_spiked_CG.DMRs.hypo.sorted.bed -b r3xr1_v_ColxC24_pat_C24_pseudo_CG.DMRs.hypo.sorted.bed > recipA_r1Vwt_pat_CGo.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_pat_C24_pseudo_CG.DMRs.hypo.sorted.bed -b r1xr3_v_C24xCol_pat_Col_spiked_CG.DMRs.hypo.sorted.bed > recipB_r1Vwt_pat_CGo.bed
cat recipA_r1Vwt_pat_CGo.bed recipB_r1Vwt_pat_CGo.bed > recipCat_ros1_v_wt_pat_CGo.bed


bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_mat_C24_pseudo_CG.DMRs.hyper.sorted.bed -b r3xr1_v_ColxC24_mat_Col_spiked_CG.DMRs.hyper.sorted.bed > recipA_r1Vwt_mat_CGr.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_mat_Col_spiked_CG.DMRs.hyper.sorted.bed -b r1xr3_v_C24xCol_mat_C24_pseudo_CG.DMRs.hyper.sorted.bed > recipB_r1Vwt_mat_CGr.bed
cat recipA_r1Vwt_mat_CGr.bed recipB_r1Vwt_mat_CGr.bed > recipCat_ros1_v_wt_mat_CGr.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_mat_C24_pseudo_CG.DMRs.hypo.sorted.bed -b r3xr1_v_ColxC24_mat_Col_spiked_CG.DMRs.hypo.sorted.bed > recipA_r1Vwt_mat_CGo.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_mat_Col_spiked_CG.DMRs.hypo.sorted.bed -b r1xr3_v_C24xCol_mat_C24_pseudo_CG.DMRs.hypo.sorted.bed > recipB_r1Vwt_mat_CGo.bed
cat recipA_r1Vwt_mat_CGo.bed recipB_r1Vwt_mat_CGo.bed > recipCat_ros1_v_wt_mat_CGo.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_pat_Col_spiked_CHG.DMRs.hyper.sorted.bed -b r3xr1_v_ColxC24_pat_C24_pseudo_CHG.DMRs.hyper.sorted.bed > recipA_r1Vwt_pat_CHGr.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_pat_C24_pseudo_CHG.DMRs.hyper.sorted.bed -b r1xr3_v_C24xCol_pat_Col_spiked_CHG.DMRs.hyper.sorted.bed > recipB_r1Vwt_pat_CHGr.bed
cat recipA_r1Vwt_pat_CHGr.bed recipB_r1Vwt_pat_CHGr.bed > recipCat_ros1_v_wt_pat_CHGr.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_pat_Col_spiked_CHG.DMRs.hypo.sorted.bed -b r3xr1_v_ColxC24_pat_C24_pseudo_CHG.DMRs.hypo.sorted.bed > recipA_r1Vwt_pat_CHGo.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_pat_C24_pseudo_CHG.DMRs.hypo.sorted.bed -b r1xr3_v_C24xCol_pat_Col_spiked_CHG.DMRs.hypo.sorted.bed > recipB_r1Vwt_pat_CHGo.bed
cat recipA_r1Vwt_pat_CHGo.bed recipB_r1Vwt_pat_CHGo.bed > recipCat_ros1_v_wt_pat_CHGo.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_mat_C24_pseudo_CHG.DMRs.hyper.sorted.bed -b r3xr1_v_ColxC24_mat_Col_spiked_CHG.DMRs.hyper.sorted.bed > recipA_r1Vwt_mat_CHGr.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_mat_Col_spiked_CHG.DMRs.hyper.sorted.bed -b r1xr3_v_C24xCol_mat_C24_pseudo_CHG.DMRs.hyper.sorted.bed > recipB_r1Vwt_mat_CHGr.bed
cat recipA_r1Vwt_mat_CHGr.bed recipB_r1Vwt_mat_CHGr.bed > recipCat_ros1_v_wt_mat_CHGr.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_mat_C24_pseudo_CHG.DMRs.hypo.sorted.bed -b r3xr1_v_ColxC24_mat_Col_spiked_CHG.DMRs.hypo.sorted.bed > recipA_r1Vwt_mat_CHGo.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_mat_Col_spiked_CHG.DMRs.hypo.sorted.bed -b r1xr3_v_C24xCol_mat_C24_pseudo_CHG.DMRs.hypo.sorted.bed > recipB_r1Vwt_mat_CHGo.bed
cat recipA_r1Vwt_mat_CHGo.bed recipB_r1Vwt_mat_CHGo.bed > recipCat_ros1_v_wt_mat_CHGo.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_pat_Col_spiked_CHH.DMRs.hyper.sorted.bed -b r3xr1_v_ColxC24_pat_C24_pseudo_CHH.DMRs.hyper.sorted.bed > recipA_r1Vwt_pat_CHHr.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_pat_C24_pseudo_CHH.DMRs.hyper.sorted.bed -b r1xr3_v_C24xCol_pat_Col_spiked_CHH.DMRs.hyper.sorted.bed > recipB_r1Vwt_pat_CHHr.bed
cat recipA_r1Vwt_pat_CHHr.bed recipB_r1Vwt_pat_CHHr.bed > recipCat_ros1_v_wt_pat_CHHr.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_pat_Col_spiked_CHH.DMRs.hypo.sorted.bed -b r3xr1_v_ColxC24_pat_C24_pseudo_CHH.DMRs.hypo.sorted.bed > recipA_r1Vwt_pat_CHHo.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_pat_C24_pseudo_CHH.DMRs.hypo.sorted.bed -b r1xr3_v_C24xCol_pat_Col_spiked_CHH.DMRs.hypo.sorted.bed > recipB_r1Vwt_pat_CHHo.bed
cat recipA_r1Vwt_pat_CHHo.bed recipB_r1Vwt_pat_CHHo.bed > recipCat_ros1_v_wt_pat_CHHo.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_mat_C24_pseudo_CHH.DMRs.hyper.sorted.bed -b r3xr1_v_ColxC24_mat_Col_spiked_CHH.DMRs.hyper.sorted.bed > recipA_r1Vwt_mat_CHHr.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_mat_Col_spiked_CHH.DMRs.hyper.sorted.bed -b r1xr3_v_C24xCol_mat_C24_pseudo_CHH.DMRs.hyper.sorted.bed > recipB_r1Vwt_mat_CHHr.bed
cat recipA_r1Vwt_mat_CHHr.bed recipB_r1Vwt_mat_CHHr.bed > recipCat_ros1_v_wt_mat_CHHr.bed

bedtools intersect -u -f 0.1 -a r1xr3_v_C24xCol_mat_C24_pseudo_CHH.DMRs.hypo.sorted.bed -b r3xr1_v_ColxC24_mat_Col_spiked_CHH.DMRs.hypo.sorted.bed > recipA_r1Vwt_mat_CHHo.bed
bedtools intersect -u -f 0.1 -a r3xr1_v_ColxC24_mat_Col_spiked_CHH.DMRs.hypo.sorted.bed -b r1xr3_v_C24xCol_mat_C24_pseudo_CHH.DMRs.hypo.sorted.bed > recipB_r1Vwt_mat_CHHo.bed
cat recipA_r1Vwt_mat_CHHo.bed recipB_r1Vwt_mat_CHHo.bed > recipCat_ros1_v_wt_mat_CHHo.bed

mkdir -p processed

CAT=$(find . -type f -name "recipCat*") 
for file in $CAT ; do
	name="${file%.*}"
	sort -k1,1 -k2,2n $file > $name".sorted.bed"
	mv $name".sorted.bed" ./processed/$name".bed"
done

for file in ./processed/*.bed; do
	name="${file%.*}"
	bedtools merge -i $file > $name".merge.bed"
done

cd ./processed/

# Presumed DME only/primary targets (mat hypomethylated in wt and in ros1)
bedtools intersect -u -f 0.1 -a recipCat_WT_matVpat_CGo.merge.bed -b recipCat_ros1_matVpat_CGo.merge.bed > WT_and_ros1_matVpat_CGo.bed
bedtools intersect -u -f 0.1 -a recipCat_WT_matVpat_CHGo.merge.bed -b recipCat_ros1_matVpat_CHGo.merge.bed > WT_and_ros1_matVpat_CHGo.bed
bedtools intersect -u -f 0.1 -a recipCat_WT_matVpat_CHHo.merge.bed -b recipCat_ros1_matVpat_CHHo.merge.bed > WT_and_ros1_matVpat_CHHo.bed

# Presumed "novel" DME and ROS1 (low meth both alleles in wildtype,  mat only hypomethylated in ros1)
bedtools intersect -v -f 0.1 -a recipCat_ros1_matVpat_CGo.merge.bed -b recipCat_WT_matVpat_CGo.merge.bed > ros1_only_matVpat_CGo.bed
bedtools intersect -v -f 0.1 -a recipCat_ros1_matVpat_CHGo.merge.bed -b recipCat_WT_matVpat_CHGo.merge.bed > ros1_only_matVpat_CHGo.bed
bedtools intersect -v -f 0.1 -a recipCat_ros1_matVpat_CHHo.merge.bed -b recipCat_WT_matVpat_CHHo.merge.bed > ros1_only_matVpat_CHHo.bed

# maternal hyper (?, mat hypermethylated in ros1 – confusing but possible)
bedtools intersect -u -f 0.1 -a recipCat_WT_matVpat_CGr.merge.bed -b recipCat_ros1_matVpat_CGr.merge.bed > WT_and_ros1_matVpat_CGr.bed
bedtools intersect -u -f 0.1 -a recipCat_WT_matVpat_CHGr.merge.bed -b recipCat_ros1_matVpat_CHGr.merge.bed > WT_and_ros1_matVpat_CHGr.bed
bedtools intersect -u -f 0.1 -a recipCat_WT_matVpat_CHHr.merge.bed -b recipCat_ros1_matVpat_CHHr.merge.bed > WT_and_ros1_matVpat_CHHr.bed

bedtools intersect -v -f 0.1 -a recipCat_ros1_matVpat_CGr.merge.bed -b recipCat_WT_matVpat_CGr.merge.bed > ros1_only_matVpat_CGr.bed
bedtools intersect -v -f 0.1 -a recipCat_ros1_matVpat_CHGr.merge.bed -b recipCat_WT_matVpat_CHGr.merge.bed > ros1_only_matVpat_CHGr.bed
bedtools intersect -v -f 0.1 -a recipCat_ros1_matVpat_CHHr.merge.bed -b recipCat_WT_matVpat_CHHr.merge.bed > ros1_only_matVpat_CHHr.bed

# gain of methylation on paternal allele in both genotypes of ros1, in any sequence context (ros1 paternal targets)

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


# gain of methylation on allele in both genotypes of ros1, in any sequence context (ros1 paternal or maternal targets)
makeallC recipCat_ros1_v_wt_pat_CGr.bed recipCat_ros1_v_wt_pat_CHGr.bed recipCat_ros1_v_wt_pat_CHHr.bed "ros1_v_wt_pat_hyper"
makeallC recipCat_ros1_v_wt_pat_CGo.bed recipCat_ros1_v_wt_pat_CHGo.bed recipCat_ros1_v_wt_pat_CHHo.bed "ros1_v_wt_pat_hypo"

makeallC recipCat_ros1_v_wt_mat_CGr.bed recipCat_ros1_v_wt_mat_CHGr.bed recipCat_ros1_v_wt_mat_CHHr.bed "ros1_v_wt_mat_hyper"
makeallC recipCat_ros1_v_wt_mat_CGo.bed recipCat_ros1_v_wt_mat_CHGo.bed recipCat_ros1_v_wt_mat_CHHo.bed "ros1_v_wt_mat_hypo"

#biallelic gain of mC? 
bedtools intersect -u -f 0.1 -a ros1_v_wt_pat_hyper_allC.merge.bed -b ros1_v_wt_mat_hyper_allC.merge.bed > ros1_v_wt_PatMat_hyper_allC_A.bed
bedtools intersect -u -f 0.1 -a ros1_v_wt_mat_hyper_allC.merge.bed -b ros1_v_wt_pat_hyper_allC.merge.bed > ros1_v_wt_PatMat_hyper_allC_B.bed
cat ros1_v_wt_PatMat_hyper_allC_A.bed ros1_v_wt_PatMat_hyper_allC_A.bed > ros1_v_wt_PatMat_hyper_allC_C.bed
sortandmerge ros1_v_wt_PatMat_hyper_allC_C



# allC DME only or ros1 only targets
makeallC WT_and_ros1_matVpat_CGo.bed WT_and_ros1_matVpat_CHGo.bed WT_and_ros1_matVpat_CHHo.bed "WT_and_ros1_matVpat_hypo"
makeallC ros1_only_matVpat_CGo.bed ros1_only_matVpat_CHGo.bed ros1_only_matVpat_CHHo.bed "ros1_only_matVpat_hypo"


###EAH 4/23/24 extra intersects and processing 


cd /lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/


makeallC C24xCol_mat_v_pat_CG.DMRs.hypo.bed C24xCol_mat_v_pat_CHG.DMRs.hypo.bed C24xCol_mat_v_pat_CHH.DMRs.hypo.bed "C24xCol_mat_v_pat_hypo"
makeallC ColxC24_mat_v_pat_CG.DMRs.hypo.bed ColxC24_mat_v_pat_CHG.DMRs.hypo.bed ColxC24_mat_v_pat_CHH.DMRs.hypo.bed "ColxC24_mat_v_pat_hypo"


# only ros1 NEITHER wt 
bedtools intersect -v -f 0.1 -a ./processed/ros1_only_matVpat_hypo_allC.merge.bed -b C24xCol_mat_v_pat_hypo_allC.merge.bed ColxC24_mat_v_pat_hypo_allC.merge.bed > ./processed/ros1_only_NEITHERcheck_matVpat_hypo_allC.bed

cd /lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/processed/

#biallelic gain of mC? 
bedtools intersect -v -f 0.1 -a ros1_v_wt_mat_hyper_allC.merge.bed -b ros1_v_wt_PatMat_hyper_allC_C.merge.bed > ros1_v_wt_matonly_hyper_allC.bed
bedtools intersect -v -f 0.1 -a ros1_v_wt_pat_hyper_allC.merge.bed -b ros1_v_wt_PatMat_hyper_allC_C.merge.bed > ros1_v_wt_patonly_hyper_allC.bed



