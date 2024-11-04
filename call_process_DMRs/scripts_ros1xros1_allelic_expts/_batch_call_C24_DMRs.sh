#!/bin/bash
##04/19/2024 EAH
##batch call DMRs in allelic data 
##/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs/_batch_call_C24_DMRs.sh
#USAGE: ./DSS_3rep_v2.r output.file context samplelist (Ax3 Bx3)

 
#prep

script_path='/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/'
DSS=$script_path'DSS_3rep_v2.r'

path='/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/r1_vs_C24/'
outpath=$path'slurmout/'

bismarkpath='/lab/solexa_gehring/elizabeth/allelic_emseq/bismark/C24_pseudo/'


mkdir -p $outpath



name="r1_v_C24"

sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=120gb --wrap " $DSS $path$name"_CG" "CG" $bismarkpath"r1_1/r1_1_C24_pseudo_CG.bed" $bismarkpath"r1_2/r1_2_C24_pseudo_CG.bed" $bismarkpath"r1_3/r1_3_C24_pseudo_CG.bed" \
$bismarkpath"C24_1/C24_1_C24_pseudo_CG.bed" $bismarkpath"C24_2/C24_2_C24_pseudo_CG.bed" $bismarkpath"C24_3/C24_3_C24_pseudo_CG.bed" "
	
sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=120gb --wrap " $DSS $path$name"_CHG" "CHG" $bismarkpath"r1_1/r1_1_C24_pseudo_CHG.bed" $bismarkpath"r1_2/r1_2_C24_pseudo_CHG.bed" $bismarkpath"r1_3/r1_3_C24_pseudo_CHG.bed" \
$bismarkpath"C24_1/C24_1_C24_pseudo_CHG.bed" $bismarkpath"C24_2/C24_2_C24_pseudo_CHG.bed" $bismarkpath"C24_3/C24_3_C24_pseudo_CHG.bed" "

sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=120gb --wrap " $DSS $path$name"_CHH" "CHH" $bismarkpath"r1_1/r1_1_C24_pseudo_CHH.bed" $bismarkpath"r1_2/r1_2_C24_pseudo_CHH.bed" $bismarkpath"r1_3/r1_3_C24_pseudo_CHH.bed" \
$bismarkpath"C24_1/C24_1_C24_pseudo_CHH.bed" $bismarkpath"C24_2/C24_2_C24_pseudo_CHH.bed" $bismarkpath"C24_3/C24_3_C24_pseudo_CHH.bed" "
