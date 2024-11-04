#!/bin/bash
#	Wrapper script for EMmap pipeline

# Version info and usage/
read -d '' usage <<"EOF"
------------------------------------------------------------------------------------
v3.1 by Elizabeth Hemenway, 03/29/2024
This script is the wrapper for running a Bismark-based processing and mapping of EM-seq
data. You can use it to trim, map, and extract methylation data from paired-end BS/EM-seq 
GZIPPED FASTQ data files. You'll need to modifiy a set of paths at the top of the script
to match your own destinations, I have not built these in as options. 

Additional options available for allelic data:
- remap unmapped reads from TAIR10 to a second genome, deduplicate/methylextract
- run AssignToAllele script written by Colette Picard
- extract mC info from outputs of AssignToAllele

Usage: _0_EMmap_wrap.sh 1.fastq 2.fastq prefix map2TAIR10spikedpUC19? map2genome2? trim? 
	map? deduplicate? methylextract? allelic? remap2unmapped? allelicmethylextract? g2methylextract?
Example usage: sbatch -p 20  --mem=5gb --output=$slurmout"Col_2_T10_EMmap.out" 
	--job-name=Col_2_EMmap --wrap=" "$scripts"_0_EMmap_wrap.sh Col_2_F.fastq.gz 
	Col_2_R.fastq.gz Col_2"
If you leave all questions blank instead of yes/no, script will default to mapping to 
pUC19-spiked TAIR10 genome (genome1), and run all non-allelic steps. 
Helper scripts need to be in the directory $scripts, defined below. 
------------------------------------------------------------------------------------
EOF
[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# required inputs
forward="$1" ## pair 1 out of 2 raw reads fastq, can be gzipped
reverse="$2" ## pair 2 out of 2 raw reads fastq, can be gzipped
output="$3" ## output prefix; ex "C24_1"
# optional genome selection, defaults to TAIR10 with pUC19 spike-in - plan to re-write eventually as a suite of genome options. 
Col=${4:-"yes"}
G2=${5:-"no"} ## map to second genome of choice?

# these are settings for testing and troubleshooting purposes, default is yes for non-allelic steps so not required arguments if you want to run everything.
trim=${6:-"yes"} 
map=${7:-"yes"}
deduplicate=${8:-"yes"}
methylextract=${9:-"yes"}
allele=${10:-"no"}
unmap=${11:-"no"}
allelicmethy=${12:-"no"}
g2methy=${13:-"no"}

# make and set paths to useful places if needed - for now replace these with your own paths, I'll make them arguments eventually. 
genome2_name="C24_pseudo" ##name of genome2 if you want to reference it
scripts=/lab/solexa_gehring/elizabeth/ros1_parental_emseq/scripts/ ##location of all referenced scripts (_0 - _4)
data=/lab/solexa_gehring/elizabeth/ros1_parental_emseq/FASTQ/ ##location of raw reads
outdir=/lab/solexa_gehring/elizabeth/ros1_parental_emseq/ ##primary output directory - outermost directory and where logfiles will go
slurmout=/lab/solexa_gehring/elizabeth/ros1_parental_emseq/slurmout/ ##where the slurm-generated log files will go
genome1=/lab/solexa_gehring/elizabeth/genomez/TAIR10_spiked/ ##Col-0 Bismark-ready genome (has pUC19 spiked in)
genome2=/lab/solexa_gehring/elizabeth/allelic_emseq/C24_genome_prep/C24_pseudo/ ##C24 Bismark-ready pseudogenome, can replace with another Bismark-prepped genome
SNPs=/lab/solexa_gehring/elizabeth/allelic_emseq/C24_genome_prep/C24/T10_C24_snps_checkDup.bed ##Col-0 to C24 snps, replace with SNP file you need
##SNPs needs format: "the SNP file must be in .bed format with the fourth field = ref>alt (e.g. A>T, where A is the reference base and T is the alt base)"
chrom_sizes=/lab/solexa_gehring/elizabeth/allelic_emseq/scripts/TAIR10.chrom.sizes ##needed for generating .bw files to browse data using IGV

# Make a log file and update with input information.
logfile="$outdir$output"_EMmap.txt

echo -e "\n------------------------------------------------------------------------------------" >> $logfile
echo -e "Bismark-based processing and mapping wrapper for allelic BS/EM-seq:\nEAH 03/29/2024 - v3.1" >> $logfile
current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
echo -e "Current date: $current_date\nCurrent time: $current_time" >> $logfile
echo -e "Forward : $forward\nReverse : $reverse\nOutput : $output" >> $logfile

# Which steps are going to be run?
echo -e "\nWhich steps to do?" >> $logfile
echo "Trim reads? - $trim" >> $logfile
echo "Map? - $map" >> $logfile
echo "Deduplicate? - $deduplicate" >> $logfile
echo "Methyl extract (all reads mapped to genome1)? - $methylextract" >> $logfile
echo "Assign2allele? - $allele" >> $logfile
echo "Remap unmapped reads to other genome? - $unmap" >> $logfile
echo "allelic methyl extract? - $allelicmethy" >> $logfile
echo "Methyl extract (all reads mapped to alt/genome2)? - $g2methy" >> $logfile


# If trim=yes, trim reads. 
if [ $trim = "yes" ]; then
	trimdir=$outdir"trimgalore/"
	mkdir -p $trimdir
    echo -e "\nTrimming reads using trim-galore (default Phred-20)\nTrimmed reads will go here: "$outdir"trimgalore/" >> $logfile
	# Trim reads using trim-galore (default Phred33 score=20)
	JID_JOB1=`sbatch --output "$slurmout""$output""trimgalore.out" --wait "$scripts"_1_trimgalore.sh $forward $reverse $trimdir | cut -d " " -f 4`
	trimdepend=afterok:$JID_JOB1
	echo "Trim-galore slurm job ID = $JID_JOB1" >> $logfile
	echo "Trim-galore run details found here: $slurmout$output"trimgalore.out"" >> $logfile
elif [ $trim = "no" ]; then 
	echo -e "\nNot trimming reads" >> $logfile
	trimdir=$outdir"trimgalore/"
	trimdepend=""
fi

namepair1="$(basename $forward .fastq.gz)"
namepair2="$(basename $reverse .fastq.gz)"
echo "basename of most files moving forward: $namepair1 or $namepair2" >> $logfile

# Trimmed reads PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz) 
trim_for=$trimdir"${namepair1%%.*}_val_1.fq.gz"
trim_rev=$trimdir"${namepair2%%.*}_val_2.fq.gz"

# Set directories for specific genome and sample - may need to personalize. need the linux for "if this directory doesn't exist, make it"
if [ $map = "yes" ]; then
	if [ "$Col" = "yes" ]; then
		genome=$genome1
		Gnome="Col_spiked"
		Gnome2=$genome2_name
		sampledir=$outdir"bismark/$Gnome/"$output"/"
		mkdir -p $sampledir
	elif [ "$Col" = "no" ]; then
		echo -e "\nnot mapping to Col genome" >> $logfile
	fi
	if [ "$G2" = "yes" ]; then
		genome=$genome2
		Gnome=$genome2_name
		sampledir=$outdir"bismark/$Gnome/"$output"/"
		mkdir -p $sampledir
	elif [ "$G2" = "no" ]; then
		echo -e "\nnot mapping to $genome2_name genome" >> $logfile
	fi
		
	# Mapping using bismark
	echo -e "\nMapping trimmed reads to $Gnome \nPair 1 : $trim_for \nPair 2 : $trim_rev" >> $logfile
	echo "Outputs for mapping to $Gnome go here: $sampledir" >>$logfile
	#	bismark_wrap needs in order: forward reverse outputdirectory outputprefix genome
	JID_JOB2=`sbatch --output $slurmout$output$Gnome"_bismark.out" --dependency=$trimdepend --wait "$scripts"_2_bismark_wrap.sh $trim_for $trim_rev $sampledir $output $genome $namepair1 $namepair2 | cut -d " " -f 4`
	mapdepend=afterok:$JID_JOB2
	echo "bismark mapping slurm job ID = $JID_JOB2" >> $logfile
	echo "$Gnome mapping details found here: $slurmout$output$Gnome"_bismark.out"" >> $logfile
elif [ $map = "no" ]; then
	echo -e "\nNot mapping reads" >> $logfile
	mapdepend=""
	if [ "$Col" = "yes" ]; then
		genome=$genome1
		Gnome="Col_spiked"
		Gnome2=$genome2_name
		sampledir=$outdir"bismark/$Gnome/"$output"/"
		mkdir -p $sampledir
		echo -e "\nAlready mapped to $genome" >> $logfile
	elif [ "$Col" = "no" ]; then
		echo -e "\nnot mapping to $genome" >> $logfile
	fi
	if [ "$G2" = "yes" ]; then
		genome=$genome2
		Gnome=$genome2_name
		sampledir=$outdir"bismark/$Gnome/"$output"/"
		mkdir -p $sampledir
		echo -e "\nAlready mapped to $genome" >> $logfile
	elif [ "$G2" = "no" ]; then
		echo -e "\nnot mapping to $genome2_name genome" >> $logfile
	fi
fi

# Mapped reads output to use  (simulated_1_bismark_bt2_pe.bam)
mapped=$sampledir"${namepair1%%.*}_val_1_bismark_bt2_pe.bam"
echo "$mapped" >> $logfile

# 3 Sort mapped bam file and deduplicate
if [ $deduplicate = "yes" ]; then
	echo -e "\nSort mapped bam file: $mapped then deduplicate" >> $logfile
	JID_JOB3=`sbatch --output $slurmout$output$Gnome"_sort_and_deduplicate.out" --dependency=$mapdepend --wait "$scripts"_3_sort_deduplicate.sh $mapped $sampledir $output | cut -d " " -f 4`
	dedupdepend=afterok:$JID_JOB3
	echo "sort and deduplicate slurm job ID = $JID_JOB3" >> $logfile
	echo "sort and deduplicate details found here: $slurmout$output$Gnome"_sort_and_deduplicate.out"" >> $logfile
elif [ $deduplicate = "no" ]; then
	echo -e "\n\nNot sorting and deduplicating" >> $logfile
	dedupdepend=""
fi

# Deduplicated reads output to use (bams/{sample}.deduplicated.bam)
deduped=$sampledir"deduped/"$output".deduplicated.bam"

###Allele-specific steps###
### -------------- meanwhile...--------------------- ###
# 2.1 Map unmapped reads to other genome in F1 hybrid data
if [ $unmap = "yes" ]; then
	unmapdir=$sampledir"remap2"$Gnome2"/"
	mkdir -p $unmapdir
	
	unmap_for=$sampledir"${namepair1%%.*}_val_1_unmapped_reads_1.fq.gz"
	unmap_rev=$sampledir"${namepair2%%.*}_val_2_unmapped_reads_2.fq.gz"

	# Mapping using bismark
	echo -e "\n\nMapping unmapped trimmed reads to other genome --> $Gnome2 \nPair 1 : $unmap_for \nPair 2 : $unmap_rev" >> $logfile
	echo "Re-mapped reads to $Gnome2 go here: $unmapdir" >>$logfile
	#	bismark_wrap needs in order: forward reverse outputdirectory outputprefix genome
	JID_JOB2a=`sbatch --output $slurmout$output$Gnome2"_remapunmapped.out" --dependency=$mapdepend "$scripts"_2_bismark_wrap.sh $unmap_for $unmap_rev $unmapdir $output $genome2 $namepair1 $namepair2 | cut -d " " -f 4`
	unmapdepend="" #afterok:$JID_JOB2a
	echo "$Gnome2 bismark remmapping of unmapped reads slurm job ID = $JID_JOB2a" >> $logfile
	echo "$Gnome2 remapping of unmapped reads details found here: $slurmout$output$Gnome2"_unmapped.out"" >> $logfile

	remapped=$unmapdir"${namepair1%%.*}_val_1_unmapped_reads_1_bismark_bt2_pe.bam"

	# 3.1 deduplicate re-mapped (unmapped) reads 
	echo -e "\nSort re-mapped bam file: $remapped then deduplicate" >> $logfile
	JID_JOB3a=`sbatch --output $slurmout$output$Gnome2"_unmapped_dedup.out" --dependency=$unmapdepend "$scripts"_3_sort_deduplicate.sh $remapped $unmapdir $output | cut -d " " -f 4`
	undedupdepend=$JID_JOB3a
	echo "sort and deduplicate slurm job ID = $JID_JOB3a" >> $logfile
	echo "sort and deduplicate details found here: $slurmout$output$Gnome2"_unmapped_dedup.out"" >> $logfile
	undeduped=$unmapdir"deduped/"$output".deduplicated.bam"
	

	# 3.15 run a2a on deduplicated remapped reads
	echo -e "\nRun assign to allele script on file: $undeduped" >> $logfile
	JID_JOB3b=`sbatch --output "$slurmout""$output""$Gnome2""_unmapped_a2a.out" --dependency=$undedupdepend "$scripts"_4a_assign2allele.sh $undeduped $unmapdir $output $scripts $SNPs $Gnome $Gnome2 | cut -d " " -f 4`
	echo "assign2allele slurm job ID = $JID_JOB4" >> $logfile
	unalleledepend=$JID_JOB3b
	echo "assign2allele details found here: $slurmout$output$Gnome2"_unmapped_a2a.out"" >> $logfile
	
elif [ $unmap = "no" ]; then
	echo -e "\nNot mapping unmapped trimmed reads to $Gnome2" >> $logfile
fi

### ------------------------------------------------ ###	
###-----3.2 Option to run assign to allele script on deduplicated reads mapped to genome 1-----###
if [ $allele = "yes" ]; then
	echo -e "\nRun assign to allele script on file: $deduped" >> $logfile
	JID_JOB4=`sbatch --output "$slurmout""$output""$Gnome""_assign2allele.out" --dependency=$dedupdepend "$scripts"_4a_assign2allele.sh $deduped $sampledir $output $scripts $SNPs $Gnome $Gnome2 | cut -d " " -f 4`
	echo "assign2allele slurm job ID = $JID_JOB4" >> $logfile
	alleledepend=$JID_JOB4
	echo "assign2allele details found here: $slurmout$output$Gnome"_assign2allele.out"" >> $logfile
	
elif [ $allele = "no" ]; then
	echo -e "\nNot running assign2allele" >> $logfile
	alleledepend=""
fi
### ------------------------------------------------ ###	
###End allele-specific steps###

# 4. Methylation Extractor - have pipeline wait here, should not finish completely until this job is done.
if [ $methylextract = "yes" ]; then
	#methyl extract whole experiment - all reads mapped to genome 1.  
	echo -e "\nRun methylation extractor on file: $deduped" >> $logfile
	JID_JOB5=`sbatch --output "$slurmout""$output""$Gnome""_methylextract.out" --dependency=$dedupdepend "$scripts"_4_methylation_extractor.sh $deduped $sampledir $output $genome | cut -d " " -f 4`
	echo "methylation extractor slurm job ID = $JID_JOB5" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome"_methylextract.out"" >> $logfile
	##methyldepend=afterok:$JID_JOB5
	methyldepend=""
	
	# 5. Summarize methyl extract by chromosome and context - note that this CX output doesn't filter by coverage. 
	echo -e "\nGet chromosome methylation summaries from $deduped" >> $logfile
	JID_JOB6=`sbatch --output "$slurmout""$output""$Gnome""_methylsummarize.out" --dependency=$methyldepend "$scripts"_5a_methyl_summarize.sh $sampledir $output $scripts $chrom_sizes | cut -d " " -f 4`
	echo "methylation summarize slurm job ID = $JID_JOB6" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome"_methylsummarize.out"" >> $logfile
	
elif [ $methylextract = "no" ]; then
	echo  -e "\nNot running methylation extractor." >> $logfile
	methyldepend=""
fi
###-----4.1 option to run methylation extractor on reads specific to each genome independently.-----###
if [ $allelicmethy = "yes" ]; then
	
	a2a=$sampledir"assign2allele/"
	remap=$sampledir"remap2"$Gnome2"/assign2allele/"
		
	genome1reads_strict=$a2a$output"_"$Gnome".bam" ## Genome 1 reads = ColxC24_2_Col.bam $output"_"$Gnome.bam 
	genome1_strict_output=$output"_"$Gnome
	genome2reads_strict=$a2a$output"_"$Gnome2".bam" ## Genome 2 reads = ColxC24_2_C24.bam $output"_"$Gnome2.bam 
	genome2_strict_output=$output"_"$Gnome2
	
	genome1remapreads_strict=$remap$output"_"$Gnome".bam" ## Genome 1 reads from remapping= ColxC24_2_Col.bam $output"_"$Gnome.bam 
	genome1remap_strict_output=$output"_"$Gnome
	genome2remapreads_strict=$remap$output"_"$Gnome2".bam" ## Genome 2 reads from remapping = ColxC24_2_C24.bam $output"_"$Gnome2.bam 
	genome2remap_strict_output=$output"_"$Gnome2
	#genome2remapreads_mapped2gnome2=$remap"deduped/"$output".deduplicated.bam" ## and unmapped-->remapped/deduped to C24 reads (mapped to what I think is pseudogenome)
	
	echo -e "\nRun methylation extractor on file: $genome1reads_strict" >> $logfile
	JID_JOB7=`sbatch --output "$slurmout""$output""$Gnome""strict_methylextract.out" --dependency=$alleledepend "$scripts"_4_methylation_extractor.sh $genome1reads_strict $a2a $output $genome | cut -d " " -f 4`
	echo "methylation extractor $Gnome strict reads slurm job ID = $JID_JOB7" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome"_a2a_methylextract.out"" >> $logfile
	methyldepend_a1=afterok:$JID_JOB7
	echo -e "\nGet chromosome methylation summaries from $genome1_strict_output" >> $logfile
	JID_JOB7a=`sbatch --output "$slurmout""$output""$Gnome""_methylsummarize.out" --dependency=$methyldepend_a1 "$scripts"_5a_methyl_summarize.sh $a2a $genome1_strict_output $scripts $chrom_sizes | cut -d " " -f 4`
	echo "methylation summarize slurm job ID = $JID_JOB7a" >> $logfile
	echo "methylation extractor details found here: $slurmout$genome1_strict_output"_methylsummarize.out"" >> $logfile
	
	#######
	echo -e "\nRun methylation extractor on file: $genome2reads_strict" >> $logfile
	JID_JOB8=`sbatch --output "$slurmout""$output""$Gnome2""strict_methylextract.out" --dependency=$alleledepend "$scripts"_4_methylation_extractor.sh $genome2reads_strict $a2a $output $genome | cut -d " " -f 4`
	echo "methylation extractor $Gnome2 strict reads slurm job ID = $JID_JOB8" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome2"_a2a_methylextract.out"" >> $logfile
	methyldepend_a2=afterok:$JID_JOB8
	echo -e "\nGet chromosome methylation summaries from $genome2_strict_output" >> $logfile
	JID_JOB8a=`sbatch --output "$slurmout""$output""$Gnome2""_methylsummarize.out" --dependency=$methyldepend_a2 "$scripts"_5a_methyl_summarize.sh $a2a $genome2_strict_output $scripts $chrom_sizes | cut -d " " -f 4`
	echo "methylation summarize slurm job ID = $JID_JOB8a" >> $logfile
	echo "methylation extractor details found here: $slurmout$genome2_strict_output"_methylsummarize.out"" >> $logfile
	
	
	######
	echo -e "\nRun methylation extractor on file: $genome1remapreads_strict" >> $logfile
	JID_JOB9=`sbatch --output "$slurmout""$output""$Gnome""strict_methylextract.out" --dependency=$unalleledepend "$scripts"_4_methylation_extractor.sh $genome1remapreads_strict $remap $output $genome | cut -d " " -f 4`
	echo "methylation extractor $Gnome strict reads slurm job ID = $JID_JOB9" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome"_remapa2a_methylextract.out"" >> $logfile
	remethyldepend_a1=afterok:$JID_JOB9
	
	echo -e "\nGet chromosome methylation summaries from $genome1remap_strict_output" >> $logfile
	JID_JOB9a=`sbatch --output "$slurmout""$output""$Gnome""_methylsummarize.out" --dependency=$remethyldepend_a1 "$scripts"_5a_methyl_summarize.sh $remap $genome1remap_strict_output $scripts $chrom_sizes | cut -d " " -f 4`
	echo "methylation summarize slurm job ID = $JID_JOB9a" >> $logfile
	echo "methylation extractor details found here: $slurmout$genome1remap_strict_output"_remap_methylsummarize.out"" >> $logfile
	
	#######
	echo -e "\nRun methylation extractor on file: $genome2remapreads_strict" >> $logfile
	JID_JOB10=`sbatch --output "$slurmout""$output""$Gnome2""strict_methylextract.out" --dependency=$unalleledepend "$scripts"_4_methylation_extractor.sh $genome2remapreads_strict $remap $output $genome | cut -d " " -f 4`
	echo "methylation extractor $Gnome strict reads slurm job ID = $JID_JOB10" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome2"_remapa2a_methylextract.out"" >> $logfile
	remethyldepend_a2=afterok:$JID_JOB10

	echo -e "\nGet chromosome methylation summaries from $genome2remap_strict_output" >> $logfile
	JID_JOB10a=`sbatch --output "$slurmout""$output""$Gnome""_methylsummarize.out" --dependency=$remethyldepend_a2 "$scripts"_5a_methyl_summarize.sh $remap $genome2remap_strict_output $scripts $chrom_sizes | cut -d " " -f 4`
	echo "methylation summarize slurm job ID = $JID_JOB10a" >> $logfile
	echo "methylation extractor details found here: $slurmout$genome2remap_strict_output"_remap_methylsummarize.out"" >> $logfile
	
	#######
	#gnomedepend=""
	gnomedepend=afterok:$JID_JOB7,$JID_JOB8,$JID_JOB9,$JID_JOB10
	# 5b step to merge allelic data and re-extract methyl info
	echo -e "\nMerge methylextract output and clean up data from from $Gnome" >> $logfile
	JID_JOB11=`sbatch --output $slurmout$output"_merge_genome_data.out" --dependency=$gnomedepend "$scripts"_5b_allelic_data_merge.sh $sampledir $output $Gnome $Gnome2 $scripts $chrom_sizes | cut -d " " -f 4`
	echo "merge genome data slurm job ID = $JID_JOB11" >> $logfile
	echo "Genome methyl data merge details found here: $slurmout$output"_merge_genome_data.out"" >> $logfile
	
fi

###-----4.2 run methyl extract in bulk on second-genome mapped reads.-----###
if [ $g2methy = "yes" ]; then
	
	remap=$sampledir"remap2"$Gnome2"/"
	deduped2=$remap"deduped/"$output".deduplicated.bam"
	## test ##
	echo -e "test sorting reads by name before g2 methy extract"
	samtools sort -n $deduped2
	wait
	##
	echo -e "\nRun methylation extractor on file: $deduped2" >> $logfile
	JID_JOB12=`sbatch --output "$slurmout""$output""$Gnome""_methylextract.out" --dependency=$undedupdepend --wait "$scripts"_4_methylation_extractor.sh $deduped2 $remap $output $genome2 | cut -d " " -f 4`
	echo "methylation extractor slurm job ID = $JID_JOB12" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome2"_methylextract.out"" >> $logfile
	methyldepend2=afterok:$JID_JOB12
	
	# 5. Summarize methyl extract by chromosome and context - note that this CX output doesn't filter by coverage. 
	echo -e "\nGet chromosome methylation summaries from $deduped2" >> $logfile
	JID_JOB12=`sbatch --output "$slurmout""$output""$Gnome""_methylsummarize.out" --dependency=$methyldepend2 "$scripts"_5a_methyl_summarize.sh $remap $output $scripts $chrom_sizes | cut -d " " -f 4`
	echo "methylation summarize slurm job ID = $JID_JOB6" >> $logfile
	echo "methylation extractor details found here: $slurmout$output$Gnome2"_methylsummarize.out"" >> $logfile
	
elif [ $g2methy = "no" ]; then
	echo  -e "\nNot running methylation extractor $Gnome2 mapped reads." >> $logfile
	methyldepend=""
fi	

wait



# Finish log file and close out
end_date=$(date +"%Y-%m-%d")
end_time=$(date +"%H:%M:%S")
echo -e "Current date: $end_date\nCurrent time: $end_time" >> $logfile
echo -e "\nEMmap is done!" >> $logfile













