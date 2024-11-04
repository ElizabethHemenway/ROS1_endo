#!/bin/bash
##07/18/2024 EAH
##batch call DMRs in allelic data - REDO because I found a typo in DSS script, working backwards through analysis. 
##/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs/_batch_callDMRs.sh
#USAGE: ./DSS_3rep_v2.r output.file context samplelist (Ax3 Bx3)

 
#prep

script_path='/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/'
DSS=$script_path'DSS_3rep_v2.r'

path='/lab/solexa_gehring/elizabeth/demethylase_em_remap2024/DMRs/'

outpath=$path'slurmout/'

bismarkpath='/lab/solexa_gehring/elizabeth/demethylase_em_remap2024/bismark/Col_spiked/'
mergea2a="merge_a2a/"

mkdir -p $outpath

cd $path

declare -a mutants=("r3" "r7")

#leaf ros1 v wt DMRs
for i in "${mutants[@]}" ; do
	name=$i
	data1=$bismarkpath$name"_1_leaf/"
	data2=$bismarkpath$name"_2_leaf/"
	data3=$bismarkpath$name"_3_leaf/"
	
	wt1=$bismarkpath"wt_1_leaf/"
	wt2=$bismarkpath"wt_2_leaf/"
	wt3=$bismarkpath"wt_3_leaf/"

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=80gb --wrap " $DSS $name"_v_wt_leaf_CG" "CG" $data1$name"_1_leaf_Col_spiked_CG_exWs_Chr2.bed" $data2$name"_2_leaf_Col_spiked_CG_exWs_Chr2.bed" $data3$name"_3_leaf_Col_spiked_CG_exWs_Chr2.bed" \
	#$wt1"wt_1_leaf_Col_spiked_CG_exWs_Chr2.bed" $wt2"wt_2_leaf_Col_spiked_CG_exWs_Chr2.bed" $wt3"wt_3_leaf_Col_spiked_CG_exWs_Chr2.bed" "
	
	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=100gb --wrap " $DSS $name"_v_wt_leaf_CHG" "CHG" $data1$name"_1_leaf_Col_spiked_CHG_exWs_Chr2.bed" $data2$name"_2_leaf_Col_spiked_CHG_exWs_Chr2.bed" $data3$name"_3_leaf_Col_spiked_CHG_exWs_Chr2.bed" \
	#$wt1"wt_1_leaf_Col_spiked_CHG_exWs_Chr2.bed" $wt2"wt_2_leaf_Col_spiked_CHG_exWs_Chr2.bed" $wt3"wt_3_leaf_Col_spiked_CHG_exWs_Chr2.bed" "

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=238gb --cpus-per-task=42 --wrap " $DSS $name"_v_wt_leaf_CHH" "CHH" $data1$name"_1_leaf_Col_spiked_CHH_exWs_Chr2.bed" $data2$name"_2_leaf_Col_spiked_CHH_exWs_Chr2.bed" $data3$name"_3_leaf_Col_spiked_CHH_exWs_Chr2.bed" \
	#$wt1"wt_1_leaf_Col_spiked_CHH_exWs_Chr2.bed" $wt2"wt_2_leaf_Col_spiked_CHH_exWs_Chr2.bed" $wt3"wt_3_leaf_Col_spiked_CHH_exWs_Chr2.bed" "

done

#leaf ros1 v wt DMRs - no exWS difference test
for i in "${mutants[@]}" ; do
	name=$i
	data1=$bismarkpath$name"_1_leaf/"
	data2=$bismarkpath$name"_2_leaf/"
	data3=$bismarkpath$name"_3_leaf/"
	
	wt1=$bismarkpath"wt_1_leaf/"
	wt2=$bismarkpath"wt_2_leaf/"
	wt3=$bismarkpath"wt_3_leaf/"

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=80gb --cpus-per-task=20 --wrap " $DSS "./DMRs_not_exWS_test/"$name"_v_wt_leaf_CG_notex" "CG" $data1$name"_1_leaf_Col_spiked_CG.bed" $data2$name"_2_leaf_Col_spiked_CG.bed" $data3$name"_3_leaf_Col_spiked_CG.bed" \
	#$wt1"wt_1_leaf_Col_spiked_CG.bed" $wt2"wt_2_leaf_Col_spiked_CG.bed" $wt3"wt_3_leaf_Col_spiked_CG.bed" "
	
	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=100gb --wrap " $DSS $name"_v_wt_leaf_CHG" "CHG" $data1$name"_1_leaf_Col_spiked_CHG_exWs_Chr2.bed" $data2$name"_2_leaf_Col_spiked_CHG_exWs_Chr2.bed" $data3$name"_3_leaf_Col_spiked_CHG_exWs_Chr2.bed" \
	#$wt1"wt_1_leaf_Col_spiked_CHG_exWs_Chr2.bed" $wt2"wt_2_leaf_Col_spiked_CHG_exWs_Chr2.bed" $wt3"wt_3_leaf_Col_spiked_CHG_exWs_Chr2.bed" "

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=238gb --cpus-per-task=42 --wrap " $DSS $name"_v_wt_leaf_CHH" "CHH" $data1$name"_1_leaf_Col_spiked_CHH_exWs_Chr2.bed" $data2$name"_2_leaf_Col_spiked_CHH_exWs_Chr2.bed" $data3$name"_3_leaf_Col_spiked_CHH_exWs_Chr2.bed" \
	#$wt1"wt_1_leaf_Col_spiked_CHH_exWs_Chr2.bed" $wt2"wt_2_leaf_Col_spiked_CHH_exWs_Chr2.bed" $wt3"wt_3_leaf_Col_spiked_CHH_exWs_Chr2.bed" "

done

#endo ros1 v wt DMRs
for i in "${mutants[@]}" ; do
	name=$i
	data1=$bismarkpath$name"_1/"
	data2=$bismarkpath$name"_2/"
	data3=$bismarkpath$name"_3/"
	
	wt1=$bismarkpath"wt_1/"
	wt2=$bismarkpath"wt_2/"
	wt3=$bismarkpath"wt_3/"

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=80gb --wrap " $DSS $name"_v_wt_CG" "CG" $data1$name"_1_Col_spiked_CG_exWs_Chr2.bed" $data2$name"_2_Col_spiked_CG_exWs_Chr2.bed" $data3$name"_3_Col_spiked_CG_exWs_Chr2.bed" \
	#$wt1"wt_1_Col_spiked_CG_exWs_Chr2.bed" $wt2"wt_2_Col_spiked_CG_exWs_Chr2.bed" $wt3"wt_3_Col_spiked_CG_exWs_Chr2.bed" "
	
	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=100gb --wrap " $DSS $name"_v_wt_CHG" "CHG" $data1$name"_1_Col_spiked_CHG_exWs_Chr2.bed" $data2$name"_2_Col_spiked_CHG_exWs_Chr2.bed" $data3$name"_3_Col_spiked_CHG_exWs_Chr2.bed" \
	#$wt1"wt_1_Col_spiked_CHG_exWs_Chr2.bed" $wt2"wt_2_Col_spiked_CHG_exWs_Chr2.bed" $wt3"wt_3_Col_spiked_CHG_exWs_Chr2.bed" "

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=120gb --wrap " $DSS $name"_v_wt_CHH" "CHH" $data1$name"_1_Col_spiked_CHH_exWs_Chr2.bed" $data2$name"_2_Col_spiked_CHH_exWs_Chr2.bed" $data3$name"_3_Col_spiked_CHH_exWs_Chr2.bed" \
	#$wt1"wt_1_Col_spiked_CHH_exWs_Chr2.bed" $wt2"wt_2_Col_spiked_CHH_exWs_Chr2.bed" $wt3"wt_3_Col_spiked_CHH_exWs_Chr2.bed" "

done

declare -a rdd=("rdd")

#endo ros1 v wt DMRs
for i in "${rdd[@]}" ; do
	name=$i
	data1=$bismarkpath$name"_1/"
	data2=$bismarkpath$name"_2/"
	data3=$bismarkpath$name"_3/"
	
	wt1=$bismarkpath"wt_1/"
	wt2=$bismarkpath"wt_2/"
	wt3=$bismarkpath"wt_3/"

	sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=80gb --wrap " $DSS $name"_v_wt_CG" "CG" $data1$name"_1_CG_exWs_Chr2Chr3.bed" $data2$name"_2_CG_exWs_Chr2Chr3.bed" $data3$name"_3_CG_exWs_Chr2Chr3.bed" \
	$wt1"wt_1_CG_exWs_Chr2Chr3.bed" $wt2"wt_2_CG_exWs_Chr2Chr3.bed" $wt3"wt_3_CG_exWs_Chr2Chr3.bed" "
	
	sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=100gb --wrap " $DSS $name"_v_wt_CHG" "CHG" $data1$name"_1_CHG_exWs_Chr2Chr3.bed" $data2$name"_2_CHG_exWs_Chr2Chr3.bed" $data3$name"_3_CHG_exWs_Chr2Chr3.bed" \
	$wt1"wt_1_CHG_exWs_Chr2Chr3.bed" $wt2"wt_2_CHG_exWs_Chr2Chr3.bed" $wt3"wt_3_CHG_exWs_Chr2Chr3.bed" "

	sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=120gb --wrap " $DSS $name"_v_wt_CHH" "CHH" $data1$name"_1_CHH_exWs_Chr2Chr3.bed" $data2$name"_2_CHH_exWs_Chr2Chr3.bed" $data3$name"_3_CHH_exWs_Chr2Chr3.bed" \
	$wt1"wt_1_CHH_exWs_Chr2Chr3.bed" $wt2"wt_2_CHH_exWs_Chr2Chr3.bed" $wt3"wt_3_CHH_exWs_Chr2Chr3.bed" "

done






