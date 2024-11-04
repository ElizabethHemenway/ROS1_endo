#!/bin/bash
##07/18/2024 EAH
##batch call DMRs in allelic data - REDO because I found a typo in DSS script, working backwards through analysis. 
##/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs/_batch_callDMRs.sh
#USAGE: ./DSS_3rep_v2.r output.file context samplelist (Ax3 Bx3)

 
#prep

script_path='/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/'
DSS=$script_path'DSS_3rep_v2.r'

path='/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/'

outpath=$path'slurmout/'

bismarkpath='/lab/solexa_gehring/elizabeth/allelic_emseq/bismark/Col_spiked/'
mergea2a="merge_a2a/"

#mkdir -p $outpath

cd $path

declare -a c24mats=("r1xr3" "C24xCol")

#self mat vs pat
for i in "${c24mats[@]}" ; do
	name=$i
	data1="$bismarkpath$name"_1/"$mergea2a"
	data2="$bismarkpath$name"_2/"$mergea2a"
	data3="$bismarkpath$name"_3/"$mergea2a"

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=120gb --wrap " $DSS $name"_mat_v_pat_CG" "CG" $data1$name"_1_C24_pseudo_CG.bed" $data2$name"_2_C24_pseudo_CG.bed" $data3$name"_3_C24_pseudo_CG.bed" \
	#$data1$name"_1_Col_spiked_CG.bed" $data2$name"_2_Col_spiked_CG.bed" $data3$name"_3_Col_spiked_CG.bed" "
	
	sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=120gb --wrap " $DSS $name"_mat_v_pat_CHG" "CHG" $data1$name"_1_C24_pseudo_CHG.bed" $data2$name"_2_C24_pseudo_CHG.bed" $data3$name"_3_C24_pseudo_CHG.bed" \
	$data1$name"_1_Col_spiked_CHG.bed" $data2$name"_2_Col_spiked_CHG.bed" $data3$name"_3_Col_spiked_CHG.bed" "

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=120gb --wrap " $DSS $name"_mat_v_pat_CHH" "CHH" $data1$name"_1_C24_pseudo_CHH.bed" $data2$name"_2_C24_pseudo_CHH.bed" $data3$name"_3_C24_pseudo_CHH.bed" \
	#$data1$name"_1_Col_spiked_CHH.bed" $data2$name"_2_Col_spiked_CHH.bed" $data3$name"_3_Col_spiked_CHH.bed" "


done

declare -a colmats=("r3xr1" "ColxC24")

for i in "${colmats[@]}" ; do
	name=$i
	data1="$bismarkpath$name"_1/"$mergea2a"
	data2="$bismarkpath$name"_2/"$mergea2a"
	data3="$bismarkpath$name"_3/"$mergea2a"

	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CG.txt' --mem=120gb --wrap " $DSS $name"_mat_v_pat_CG" "CG" $data1$name"_1_Col_spiked_CG.bed" $data2$name"_2_Col_spiked_CG.bed" $data3$name"_3_Col_spiked_CG.bed" \
	#$data1$name"_1_C24_pseudo_CG.bed" $data2$name"_2_C24_pseudo_CG.bed" $data3$name"_3_C24_pseudo_CG.bed" "
	
	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHG.txt' --mem=120gb --wrap " $DSS $name"_mat_v_pat_CHG" "CHG" $data1$name"_1_Col_spiked_CHG.bed" $data2$name"_2_Col_spiked_CHG.bed" $data3$name"_3_Col_spiked_CHG.bed" \
	#$data1$name"_1_C24_pseudo_CHG.bed" $data2$name"_2_C24_pseudo_CHG.bed" $data3$name"_3_C24_pseudo_CHG.bed" "
	
	#sbatch -p 20 --job-name=DSS --output $outpath$name'_CHH.txt' --mem=120gb --wrap " $DSS $name"_mat_v_pat_CHH" "CHH" $data1$name"_1_Col_spiked_CHH.bed" $data2$name"_2_Col_spiked_CHH.bed" $data3$name"_3_Col_spiked_CHH.bed" \
	#$data1$name"_1_C24_pseudo_CHH.bed" $data2$name"_2_C24_pseudo_CHH.bed" $data3$name"_3_C24_pseudo_CHH.bed" "

done


# ros1 vs wt DMRs by parent of origin


callDMRs () {
	name1=$1
	name2=$2
	echo 'sample 1 = $name1'
	echo 'sample 2 = $name2'
	
	if [[ $3 == "C24mat" ]] ; then
  		echo "$name1 in c24mats"
  		matG="C24_pseudo"
  		patG="Col_spiked"
	else
 		echo "$name1 in colmats"
 		matG="Col_spiked"
 		patG="C24_pseudo"
	fi
	
	data1="$bismarkpath$name1"_1/"$mergea2a"
	data2="$bismarkpath$name1"_2/"$mergea2a"
	data3="$bismarkpath$name1"_3/"$mergea2a"

	data4="$bismarkpath$name2"_1/"$mergea2a"
	data5="$bismarkpath$name2"_2/"$mergea2a"
	data6="$bismarkpath$name2"_3/"$mergea2a"
	
	sbatch -p 20 --job-name=DSS --output $outpath$name1'_v_'$name2'_CG.txt' --mem=120gb --wrap " $DSS $name1"_v_"$name2"_mat_"$matG"_CG" "CG" $data1$name1"_1_"$matG"_CG.bed" $data2$name1"_2_"$matG"_CG.bed" $data3$name1"_3_"$matG"_CG.bed" \
	$data4$name2"_1_"$matG"_CG.bed" $data5$name2"_2_"$matG"_CG.bed" $data6$name2"_3_"$matG"_CG.bed" "
	
	sbatch -p 20 --job-name=DSS --output $outpath$name1'_v_'$name2'_CHG.txt' --mem=120gb --wrap " $DSS $name1"_v_"$name2"_mat_"$matG"_CHG" "CHG" $data1$name1"_1_"$matG"_CHG.bed" $data2$name1"_2_"$matG"_CHG.bed" $data3$name1"_3_"$matG"_CHG.bed" \
	$data4$name2"_1_"$matG"_CHG.bed" $data5$name2"_2_"$matG"_CHG.bed" $data6$name2"_3_"$matG"_CHG.bed" "

	sbatch -p 20 --job-name=DSS --output $outpath$name1'_v_'$name2'_CHH.txt' --mem=120gb --wrap " $DSS $name1"_v_"$name2"_mat_"$matG"_CHH" "CHH" $data1$name1"_1_"$matG"_CHH.bed" $data2$name1"_2_"$matG"_CHH.bed" $data3$name1"_3_"$matG"_CHH.bed" \
	$data4$name2"_1_"$matG"_CHH.bed" $data5$name2"_2_"$matG"_CHH.bed" $data6$name2"_3_"$matG"_CHH.bed" "

	

	sbatch -p 20 --job-name=DSS --output $outpath$name1'_v_'$name2'_CG.txt' --mem=120gb --wrap " $DSS $name1"_v_"$name2"_pat_"$patG"_CG" "CG" $data1$name1"_1_"$patG"_CG.bed" $data2$name1"_2_"$patG"_CG.bed" $data3$name1"_3_"$patG"_CG.bed" \
	$data4$name2"_1_"$patG"_CG.bed" $data5$name2"_2_"$patG"_CG.bed" $data6$name2"_3_"$patG"_CG.bed" "
	
	sbatch -p 20 --job-name=DSS --output $outpath$name1'_v_'$name2'_CHG.txt' --mem=120gb --wrap " $DSS $name1"_v_"$name2"_pat_"$patG"_CHG" "CHG" $data1$name1"_1_"$patG"_CHG.bed" $data2$name1"_2_"$patG"_CHG.bed" $data3$name1"_3_"$patG"_CHG.bed" \
	$data4$name2"_1_"$patG"_CHG.bed" $data5$name2"_2_"$patG"_CHG.bed" $data6$name2"_3_"$patG"_CHG.bed" "

	sbatch -p 20 --job-name=DSS --output $outpath$name1'_v_'$name2'_CHH.txt' --mem=120gb --wrap " $DSS $name1"_v_"$name2"_pat_"$patG"_CHH" "CHH" $data1$name1"_1_"$patG"_CHH.bed" $data2$name1"_2_"$patG"_CHH.bed" $data3$name1"_3_"$patG"_CHH.bed" \
	$data4$name2"_1_"$patG"_CHH.bed" $data5$name2"_2_"$patG"_CHH.bed" $data6$name2"_3_"$patG"_CHH.bed" "

}

#r1_3 pat vs Col pat

#callDMRs "r1xr3" "C24xCol" "C24mat"

#callDMRs "r3xr1" "ColxC24" "Colmat"


