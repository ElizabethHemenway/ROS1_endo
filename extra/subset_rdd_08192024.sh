#!/bin/bash
# assign to allele slurm wrapper
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=subset    # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=4             # number of cpus/threads requested.
#SBATCH --mem=230gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use


##08/19/2024
##subsample rdd 2 before mapping

endo_data="/archive/gehring/2022.07.21-31497/lab/solexa_gehring/elizabeth/demethylase_endosperm_emseq/demethylase_endosperm_emseq_FASTQ/" ##location of raw reads
path="/lab/solexa_gehring/elizabeth/demethylase_em_remap2024/"
fastout=$path"FASTQ/"
slurmout=$path"slurmout/"

#rdd2_1=$endo_data"rdd_2_f.fastq.gz"
#rdd2_2=$endo_data"rdd_2_r.fastq.gz"

#cp $rdd2_1 $rdd2_2 $fastout 
#wait 

rdd2_1=$fastout"rdd_2_f.fastq.gz"
rdd2_2=$fastout"rdd_2_r.fastq.gz"

#seqtk sample -s100 $rdd2_1 150000000 > $fastout"rdd_2_150m_f.fastq"
#seqtk sample -s100 $rdd2_2 150000000 > $fastout"rdd_2_150m_r.fastq"


gzip $fastout"rdd_2_150m_f.fastq" $fastout"rdd_2_150m_r.fastq"