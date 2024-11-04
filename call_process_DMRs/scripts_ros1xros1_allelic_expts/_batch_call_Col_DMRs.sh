#!/bin/bash
##04/19/2024 EAH
##batch call DMRs in allelic data 
##/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs/_batch_call_Col_DMRs.sh
#USAGE: ./DSS_3rep_v2.r output.file context samplelist (Ax2 Bx2)

 
#prep

script_path='/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/'
DSS=$script_path'DSS_2rep_v2.r'

path='/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/r3_vs_Col/'
outpath=$path'slurmout/'

bismarkpath='/lab/solexa_gehring/elizabeth/allelic_emseq/bismark/Col_spiked/'


mkdir -p $outpath

cd $path

name="r3_vs_Col"

sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=120gb --wrap " $DSS $path$name"_CG" "CG" $bismarkpath"r3_1/r3_1_Col_spiked_CG.bed" $bismarkpath"r3_2/r3_2_Col_spiked_CG.bed" \
$bismarkpath"Col_1/Col_1_Col_spiked_CG.bed" $bismarkpath"Col_2/Col_2_Col_spiked_CG.bed" "
	
sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=120gb --wrap " $DSS $path$name"_CHG" "CHG" $bismarkpath"r3_1/r3_1_Col_spiked_CHG.bed" $bismarkpath"r3_2/r3_2_Col_spiked_CHG.bed" \
$bismarkpath"Col_1/Col_1_Col_spiked_CHG.bed" $bismarkpath"Col_2/Col_2_Col_spiked_CHG.bed" "

sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=120gb --wrap " $DSS $path$name"_CHH" "CHH" $bismarkpath"r3_1/r3_1_Col_spiked_CHH.bed" $bismarkpath"r3_2/r3_2_Col_spiked_CHH.bed" \
$bismarkpath"Col_1/Col_1_Col_spiked_CHH.bed" $bismarkpath"Col_2/Col_2_Col_spiked_CHH.bed" "
