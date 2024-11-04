#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=tidy_methyldata   # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1             # number of cpus/threads requested.
#SBATCH --mem=30gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use
 
sampledir="$1" ##output directory for whole sample (where mapped reads should be also)
output="$2" ##output prefix ex: Col_1 or C24xCol_1_C24 if allelic 
gnome1="$3"
scripts="/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/"
chrom_sizes="/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/TAIR10.chrom.sizes" ##needed for generating .bw files to browse data using IGV

minimum=5 ##minimum reads for inclusion in browsing file - adjust here if you want to. 


methyldir_1="$sampledir""methylextract/"  ##methylextract output directory
outdir=$sampledir

methyl1_1=$methyldir_1$output".deduplicated.CX_report.txt"

start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")
echo -e "Start time for data tidy:\nCurrent date: $start_date\nCurrent time: $start_time" 

# Unzip methylation reports
echo -e "Unzip CX report..." 
gunzip $methyl1_1'.gz'

cd $outdir

# Run python script for dataframe manipulation - outputs a methylation summary table by chromosome and context. 
echo -e "Generating useful methylation dataframe and bedgraph files..." 

$scripts'NoMergeMethylData.py' $methyldir_1 $output $methyl1_1 $gnome1 $outdir $minimum

wait

# Make IGV browsing files
#echo "Making .bw files for browsing in IGV..."
#FILES=$(find . -name "*.min5.bedGraph")
#for file in $FILES ; do
	#sort -k1,1 -k2,2n "$file" > "${file%.*}"".sort.bedGraph"
	#sorted="${file%.*}"".sort.bedGraph"
	#forIGV="${sorted%.*}.forIGV.bedGraph"
	#rm -f $forIGV "${forIGV%.*}.bw" #remove IGV files if they exist to avoid error
	#echo "generate file = ""$forIGV"
	#awk -F'\t' '{ $4 = ($4 == 0 ? -10 : $4) } 1' OFS="\t" "$sorted" > $forIGV
	#bedGraphToBigWig $forIGV "$chrom_sizes" "${forIGV%.*}.bw"
	#rm -f "$file" 
#done

#echo -e "Rezip files..." 

#gzip $methyl1_1
#gzip *.bedGraph

##Add more details to the output log
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time: Current date: $current_date\nCurrent time: $current_time"
echo "Data clean up done!"
