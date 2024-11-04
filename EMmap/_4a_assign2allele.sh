#!/bin/bash
# assign to allele slurm wrapper
# One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=allele    # friendly name for job.
#SBATCH --nodes=1                      # ensure cpus are on one node
#SBATCH --ntasks=1                     # run a single task
#SBATCH --cpus-per-task=1             # number of cpus/threads requested.
#SBATCH --mem=20gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use

current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "Start time for assign2allele: Current date: $current_date\nCurrent time: $current_time" 

deduped="$1" ##deduplicated paired-end reads
sampledir="$2" ##output directory for sample (where mapped reads should be also)
outdir="$sampledir""assign2allele/"  ##methylextract output directory
output="$3" ##output prefix ex: Col_1
scripts="$4" ##path to script directory (including assign_to_allele)
SNPs="$5" ##needs format "the SNP file must be in .bed format with the fourth field = ref>alt (e.g. A>T, where A is the reference base and T is the alt base)"
Gnome1="$6" ##reference genome/genome 1
Gnome2="$7" ##alt genome/genome 2


mkdir -p $outdir
cd $outdir

echo "sorting deduplicated bam file by coordinate..."
# samtools sort
sortout=$sampledir"deduped/"$output".deduplicated.sort.bam"
samtools sort -o $sortout $deduped

echo "converting bam to sam..."
# bam to sam
sortsam=$sampledir"deduped/"$output".deduplicated.sort.sam"
samtools view -h -o $sortsam $sortout

echo "running assign_to_allele..."
# run assign2allele script
$scripts"assign_to_allele_3.py" $SNPs $sortsam $output --refname $Gnome1 --altname $Gnome2 --bisulfite --allow_sub

wait 

echo "convert sams to bams..."
# convert sams to bams
samtools view -o $output"_"$Gnome1.bam $output"_"$Gnome1.sam
samtools view -o $output"_"$Gnome2.bam $output"_"$Gnome2.sam
samtools view -o $outpu"_confl".bam $output"_confl".sam
samtools view -o $output"_none".bam $output"_none".sam


# Delete some redundant files to save space
echo "checking that bam files are not truncated before deleting redundancy..."
samtools quickcheck -v *.bam > bad_bams.txt   
if [ -s bad_bams.txt ]; then
	echo "some files seem bad, check"
    # The file is not-empty.
	# keep the extra files for more checking
elif [ ! -s bad_bams.txt ]; then
     # The file is empty.
	echo "all files seem ok. deleted files $sortsam $deduped and assign2allele output sam files for space saving"
    rm -f $sortsam $deduped $output"_"$Gnome1.sam $output"_"$Gnome2.sam $output"_confl".sam $output"_none".sam bad_bams.txt
    #rename sorted reads file to avoid confusion later.
	mv $sortout $sampledir"deduped/"$output".deduplicated.bam"
fi



# Add more details to the output log
echo "find allele assignments here : $outdir"
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "End time for assign2allele: Current date: $current_date\nCurrent time: $current_time"
echo "Assign2allele done!"



