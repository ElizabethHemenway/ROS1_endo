#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=deduplicate    # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=10             # number of cpus/threads requested.
#SBATCH --mem=100gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

mapped="$1"  ##mapped reads 
outdir="$2"  ##output directory (sample directory in _0_script)
output="$3" ##output prefix ex: Col_1

cd $outdir
dedupedir=$outdir"deduped/"
mkdir -p $dedupedir

# sort forward+reverse mapped read files
start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")

echo "Sorting mapped reads by name..."
echo -e "Start time for sort: Current date: $start_date\nCurrent time: $start_time" 

sortout="${mapped%%.*}_sorted.bam"

samtools sort -o $sortout -n $mapped
wait

end_date=$(date +"%Y-%m-%d")
end_time=$(date +"%H:%M:%S")

##log
echo -e "mapped read bam file sorted by read name.\nOutput found here: $dedupedir$output"_sorted.bam""
echo -e "End time for sort: Current date: $end_date\nCurrent time: $end_time" 

echo "Deduplicating mapped reads using bismark deduplicate..."
echo -e "Start time for deduplication: Current date: $dedupestart_date\nCurrent time: $dedupestart_time" 

##Deduplicate
dedupestart_date=$(date +"%Y-%m-%d")
dedupestart_time=$(date +"%H:%M:%S")

echo "Deduplicating mapped and sorted reads..."
deduplicate_bismark -p --outfile $output --output_dir $dedupedir $sortout

wait 

# Delete some redundant files to save space
echo "checking that bam files are not truncated before deleting..."
samtools quickcheck -v *.bam > bad_bams.txt
if [ -s bad_bams.txt ]; then
	echo "some files seem bad, check"
    # The file is not empty?
	# keep the extra files for more checking
elif [ ! -s bad_bams.txt ]; then
    # The file is empty.
	rm -f $mapped bad_bams.txt
	echo "all files seem ok. deleted files $mapped, instead keep sorted $sortout for space saving"
fi


##Add more details to the output log
echo "find deduped outputs here: $dedupedir"
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time for deduplication: Current date: $current_date\nCurrent time: $current_time"
echo "deduplication done!"




