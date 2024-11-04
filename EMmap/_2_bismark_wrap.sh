#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=bismark    # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=10             # number of cpus/threads requested.
#SBATCH --mem=230gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

forward="$1" ##trimmed reads forward
reverse="$2" ##trimmed reads reverse
outdir="$3" ##output directory (sample directory in _0_script)
output="$4" ##output prefix
genome="$5" ##path to BS genome
namepair1="$6"
namepair2="$7"

cd "$outdir"

start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")

echo -e "Mapping now..." 

bismark -N 1 -X 600 --un -o "$outdir" --genome $genome -1 $forward -2 $reverse 


wait 

# Rename unmapped read files for later
mv $namepair1"_val_1.fq.gz_unmapped_reads_1.fq.gz" $namepair1"_val_1_unmapped_reads_1.fq.gz"
mv $namepair2"_val_2.fq.gz_unmapped_reads_2.fq.gz" $namepair2"_val_2_unmapped_reads_2.fq.gz"

# Add more details to the output log
echo "find mapped reads here: $outdir" 
echo -e "Start time for bismark mapping: Current date: $start_date\nCurrent time: $start_time" 
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time for bismark mapping: Current date: $current_date\nCurrent time: $current_time"
echo "mapping done!"



