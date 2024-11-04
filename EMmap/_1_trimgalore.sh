#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=trimgalore     # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1             # number of cpus/threads requested.
#SBATCH --mem=8gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

forward="$1" ##raw data forward reads
reverse="$2" ##raw data reverse reads
outdir="$3" ##output directory

start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")

# Trim reads (default Phred score=20)
echo -e "Running bismark trimgalore..." 

trim_galore --phred33 --paired --illumina -o $outdir $forward $reverse

wait 

# Add more details to the output log
echo "trimmed reads will go here: $outdir" 
echo -e "Start time for trim_galore: Current date: $start_date\nCurrent time: $start_time" 
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time for trim_galore: Current date: $current_date\nCurrent time: $current_time"
echo "trimming done!"
