#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=window     # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1             # number of cpus/threads requested.
#SBATCH --mem=42gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

# 01/29/2025 - make random windows to compare to DMRs

genome=/lab/solexa_gehring/elizabeth/genomez/TAIR10.chrom.sizes
regions=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/random_regions/regions/
data=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/ends_hmaps/data/
scripts=/lab/solexa_gehring/elizabeth/ros1_endo_review_analysis/random_regions/scripts/


# make genome windows to test
bedtools makewindows -g $genome -w 200 -s 50 > $regions"TAIR10_200bp_step50_windows.bed"

#intersect with ros1 dmrs (not ros1 dmr)

r3_hyper=$regions"r3_v_wt_allC_hyper_allC.merge.bed"
r7_hyper=$regions"r7_v_wt_allC_hyper_allC.merge.bed"



cd $regions

bedtools intersect -v -a $regions"TAIR10_200bp_step50_windows.bed" -b $r3_hyper > $regions"not_r3hyper.bed"
bedtools intersect -v -a $regions"TAIR10_200bp_step50_windows.bed" -b $r7_hyper > $regions"not_r7hyper.bed"

bedtools intersect -f 1 -a $regions"not_r3hyper.bed" -b $regions"not_r7hyper.bed" > $regions"not_ros1_windows.bed"

#sumby WT methylation over regions
FILES=$(find $data -type f -name '*min5.bed' -name '*wt*')
for file in $FILES ; do
	echo $file
	base=$(basename -- "$file" .bed)
	output="./sumby_out/"$base"_sumby_not_ros1_windows.bed"
	#echo $output
	$scripts'sumByFeature.sh' -i "$file" -r $regions"not_ros1_windows.bed" -o $path$output -m 1 -x 5 
	#needs to be unmethylated in wild-type so you could see gain of mC
	#filter by less than 50% methylated all context
	$scripts"filter_sumby_output.py" "$output" "$base" 0.5 "<"

done

wt1=wt_1_Col_spiked_
wt2=wt_2_Col_spiked_
wt3=wt_3_Col_spiked_

CG=CG_min5_lessthanequalto_0.5.bed
CHG=CHG_min5_lessthanequalto_0.5.bed
CHH=CHH_min5_lessthanequalto_0.5.bed

sumout=./sumby_out/

cd $regions

bedtools intersect -u -f 1 -a $wt1$CG -b $wt1$CHG > wt_1_not_r3_thresh_CG_HG.bed
bedtools intersect -u -f 1 -a wt_1_not_r3_thresh_CG_HG.bed -b $wt1$CHH > wt_1_not_r3_thresh.bed

bedtools intersect -u -f 1 -a $wt2$CG -b $wt2$CHG > wt_2_not_r3_thresh_CG_HG.bed
bedtools intersect -u -f 1 -a wt_2_not_r3_thresh_CG_HG.bed -b $wt2$CHH > wt_2_not_r3_thresh.bed

bedtools intersect -u -f 1 -a $wt3$CG -b $wt3$CHG > wt_3_not_r3_thresh_CG_HG.bed
bedtools intersect -u -f 1 -a wt_3_not_r3_thresh_CG_HG.bed -b $wt3$CHH > wt_3_not_r3_thresh.bed

bedtools intersect -u -f 1 -a wt_1_not_r3_thresh.bed -b wt_2_not_r3_thresh.bed > wt_12_not_r3_thresh.bed
bedtools intersect -u -f 1 -a wt_12_not_r3_thresh.bed -b wt_3_not_r3_thresh.bed > wt_not_r3_thresh.bed

#select random subset
N=1000
shuf -n $N wt_not_r3_thresh.bed > wt_not_r3_thresh_random1.bed
shuf -n $N wt_not_r3_thresh.bed > wt_not_r3_thresh_random2.bed
shuf -n $N wt_not_r3_thresh.bed > wt_not_r3_thresh_random3.bed


#bedtools window with genes and TEs. 

#bedtools window -a $regions"genes.bed" -b $regions"wt_not_r3_thresh_random1.bed" -w 1000 > genes_1kb_wt_not_r3_thresh_random1.bed
#bedtools window -a $regions"TE_fragments.bed" -b $regions"wt_not_r3_thresh_random1.bed" -w 1000 > TEs_1kb_wt_not_r3_thresh_random1.bed


#bedtools window -a $regions"genes.bed" -b $regions"wt_not_r3_thresh_random2.bed" -w 1000 > genes_1kb_wt_not_r3_thresh_random2.bed
#bedtools window -a $regions"TE_fragments.bed" -b $regions"wt_not_r3_thresh_random2.bed" -w 1000 > TEs_1kb_wt_not_r3_thresh_random2.bed

#bedtools window -a $regions"genes.bed" -b $regions"wt_not_r3_thresh_random3.bed" -w 1000 > genes_1kb_wt_not_r3_thresh_random3.bed
#bedtools window -a $regions"TE_fragments.bed" -b $regions"wt_not_r3_thresh_random3.bed" -w 1000 > TEs_1kb_wt_not_r3_thresh_random3.bed


