#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=merge_genomes   # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1             # number of cpus/threads requested.
#SBATCH --mem=30gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

sampledir="$1" ##output directory for whole sample (where mapped reads should be also)
output="$2" ##output prefix ex: Col_1 or C24xCol_1_C24 if allelic 
gnome1="$3"
gnome2="$4"
scripts="$5"
chrom_sizes="$6"

minimum=5 ##minimum reads for inclusion in browsing file - adjust here if you want to. 


methyldir_1="$sampledir""assign2allele/methylextract/"  ##methylextract output directory
methyldir_2="$sampledir""remap2"$gnome2"/assign2allele/methylextract/"  ##methylextract output directory
outdir=$sampledir"merge_a2a/"
rm -r -f $outdir
mkdir -p $outdir

methyl1_1=$methyldir_1$output"_"$gnome1".CX_report.txt"
methyl1_2=$methyldir_2$output"_"$gnome1".CX_report.txt"

methyl2_1=$methyldir_1$output"_"$gnome2".CX_report.txt"
methyl2_2=$methyldir_2$output"_"$gnome2".CX_report.txt"

start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")
echo -e "Start time for allelic data merge and cleanup:\nCurrent date: $start_date\nCurrent time: $start_time" 

# Unzip methylation reports
echo -e "Unzip CX reports..." 
gunzip $methyl1_1'.gz' $methyl1_2'.gz' $methyl2_1'.gz' $methyl2_2'.gz'

cd $outdir

# Run python script for dataframe manipulation - outputs a methylation summary table by chromosome and context. 
echo -e "Generating merged methylation data for each genome and bedgraph files..." 

$scripts'MergeMethylData.py' $methyldir_1 $methyldir_2 $output $methyl1_1 $methyl1_2 $methyl2_1 $methyl2_2 $gnome1 $gnome2 $outdir $minimum

wait

# Make IGV browsing files
echo "Making .bw files for browsing in IGV..."
FILES=$(find . -name "*.min5.bedGraph")
for file in $FILES ; do
	sort -k1,1 -k2,2n "$file" > "${file%.*}"".sort.bedGraph"
	sorted="${file%.*}"".sort.bedGraph"
	forIGV="${sorted%.*}.forIGV.bedGraph"
	rm -f $forIGV "${forIGV%.*}.bw" #remove IGV files if they exist to avoid error
	echo "generate file = ""$forIGV"
	awk -F'\t' '{ $4 = ($4 == 0 ? -10 : $4) } 1' OFS="\t" "$sorted" > $forIGV
	bedGraphToBigWig $forIGV "$chrom_sizes" "${forIGV%.*}.bw"
	rm -f "$file" 
done

echo -e "Rezip files..." 

gzip $methyl1_1 $methyl1_2 $methyl2_1 $methyl2_2
gzip *.bedGraph

##Add more details to the output log
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time: Current date: $current_date\nCurrent time: $current_time"
echo "Allelic merge and clean up done!"
