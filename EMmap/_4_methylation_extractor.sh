#!/bin/bash
# Configuration values for SLURM job submission.
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=methylextract    # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=30             # number of cpus/threads requested.
#SBATCH --mem=80gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

deduped="$1" ##deduplicated paired-end reads
sampledir="$2" ##output directory for sample (where mapped reads should be also)
outdir="$sampledir""methylextract/"  ##methylextract output directory
output="$3" ##output prefix ex: Col_1
genome="$4" ##path to genome used for mapping

cd $sampledir
#rm -r -f $outdir
mkdir -p $outdir

start_date=$(date +"%Y-%m-%d")
start_time=$(date +"%H:%M:%S")

echo -e "Running bismark methyl extractor..." 

bismark_methylation_extractor -p --no_header --multicore 10 --gzip --bedGraph --CX_context --comprehensive --cytosine_report --output $outdir --genome_folder $genome $deduped

##Add more details to the output log
echo "find methyl extract outputs here: $outdir"
echo -e "Start time for methyl extract: Current date: $start_date\nCurrent time: $start_time" 
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time for methyl extract: Current date: $current_date\nCurrent time: $current_time"
echo "methyl extract done!"

