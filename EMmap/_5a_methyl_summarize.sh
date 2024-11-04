#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=cleanup_methyl    # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=2             # number of cpus/threads requested.
#SBATCH --mem=30gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

sampledir="$1" ##output directory for sample (where mapped reads should be also)
methyldir="$sampledir""methylextract/"  ##methylextract output directory
output="$2" ##output prefix ex: Col_1 or C24xCol_1_C24 if allelic 
scripts="$3"
chrom_sizes="$4"
 
minimum=5 ##minimum reads for inclusion in browsing file - adjust here if you want to. 

start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")

echo "methyldir = $methyldir"
cd $methyldir

echo -e "Unzip CX report and coverage file..." 

# Unzip chromosome methylation reports and coverage files
gunzip *.CX_report.txt.gz *.bismark.cov.gz

methylout=$( find . \( -name "*.CX_report.txt" -and -name "*$output*" \))
coverageout=$( find . \( -name "*.bismark.cov" -and -name "*$output*" \))

echo "methylout = $methylout"
echo "coverageout = $coverageout"

# Run python script for dataframe manipulation - outputs a methylation summary table by chromosome and context. 
echo -e "Generating chromosome summary output and bedgraph files..." 

$scripts'make_browsing_and_chrSummary_2.py' $methyldir $output $methylout $coverageout $minimum

wait

# Make IGV browsing files
echo "Making .bw files for browsing in IGV..."
#FILES=$(find . -name "*.min5.bedGraph")
#for file in $FILES ; do
	sort -k1,1 -k2,2n "$file" > "${file%.*}"".sort.bedGraph"
	sorted="${file%.*}"".sort.bedGraph"
	forIGV="${sorted%.*}.forIGV.bedGraph"
	#rm -f $forIGV "${forIGV%.*}.bw" #remove IGV files if they exist to avoid error
	echo "generate file = ""$forIGV"
	awk -F'\t' '{ $4 = ($4 == 0 ? -10 : $4) } 1' OFS="\t" "$sorted" > $forIGV
	bedGraphToBigWig $forIGV "$chrom_sizes" "${forIGV%.*}.bw"
	#rm -f "$file" 
	##bedGraphToBigWig "${file%.*}.sort.bedGraph" "$chrom_sizes" "${file%.*}.bw"
#done

#echo -e "Rezip files..." 

#gzip $methylout $coverageout
#gzip *.bedGraph

##Add more details to the output log
echo -e "Start time for cleanup and summary file generation:\nCurrent date: $start_date\nCurrent time: $start_time" 
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time: Current date: $current_date\nCurrent time: $current_time"
echo "Clean up and summarize done!"
