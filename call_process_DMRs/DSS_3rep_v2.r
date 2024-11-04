#!/usr/bin/env Rscript


###
###  For a methyl-seq experiment, compare methylation levels
###    of 2 samples, 2 replicates each
###    using the function DMLtest() of the R package DSS.
###
###  Version 1.0: 19 July 2021
###  Version 1.1: 22 July 2021 -- Add option of filtering sites by difference before calculating "filtered FDRs"
###                               Add option of minimum coverage (to filter out low-coverage sites before trying to analyze)
###  Version 1.2: 23 July 2021 -- Add first-version of DMR calling (with default parameters)
###                               Round methylation differences in bedgraph files (to avoid scientific notation, which genome browsers don't like) 
###  Version 2.0: 27 April 2022 -- Modify code to expect 2 samples, 3 biological replicates of each
###
###  George Bell, Whitehead Bioinformatics and Research Computing
###  
###  Reference for method and R implementation
###  Park Y, Wu H (2016). ìDifferential methylation analysis for BS-seq data 
###  under general experimental design.î Bioinformatics. doi: 10.1093/bioinformatics/btw026.
###
###  Details of the DSS package and the DMLtest() function
###    https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html#3_Using_DSS_for_BS-seq_differential_methylation_analysis
###    https://www.rdocumentation.org/packages/DSS/versions/2.12.0/topics/DMLtest
###
###  Version 2.1: 5 April 2024 -- Modify code to take DMR-calling paramaters as arguments. 


#####  Analysis parameters (can be modified by the user)  #####

# Should methylation levels be smoothed?
# smoothing = TRUE
smoothing = FALSE
# How wide a window should be smooth (if smoothing=TRUE)?
smoothing.span = 500
# How many cores should we use (to speed up analysis)?
num.cores = 4
# Generic sample names (2 samples, 3 biological replicates of each)
sample.names = c("SampleA.rep1", "SampleA.rep2", "SampleA.rep3", "SampleB.rep1", "SampleB.rep2", "SampleB.rep3")
# How many reads are required across a site?  If this minimum isn't met, counts for that site are removed
minimum.coverage = 5
# If the methylation fraction difference is less than this, ignore that site and don't use it to calculate FDR
min.meth.difference = 0

#####  End of analysis parameters  #####


# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./compare_methylation_counts.R" }

if (length(commandArgs()) < 10)
{
	message("\nFor a methyl-seq experiment, compare methylation levels")
	message("  of 2 samples, 3 replicates each")
	message("  using the function DMLtest() from the R package DSS.\n")
	message("USAGE: ", this.script, " output.file sampleA.rep1.file sampleA.rep2.file sampleA.rep3.file sampleB.rep1.file sampleB.rep2.file sampleB.rep3.file\n")
	
	message("Details: Expected input file format: chr, start, end, unmeth readcount, meth readcount, percent meth")
	message("         smoothing => ", smoothing)
	if (smoothing == TRUE)
	{
		message("         smoothing span => ", smoothing.span)
	}
	message("         number of cores => ", num.cores)
	message("         (These parameters can be modified at the top of the script.)\n")
	q()
}

# Arguments are the input BED files
output.filename = commandArgs()[6]
context = commandArgs()[7]
input.filenames = commandArgs()[8 : length(commandArgs())]

### DMR-calling parameters
if (context == "CG")
{
	DMR.delta = 0.3
}
if (context == "CHG")
{
	DMR.delta = 0.2
}
if (context == "CHH")
{
	DMR.delta = 0.1
}
DMR.p.threshold = 1e-02
DMR.minlen = 50
DMR.minCG = 5
DMR.distance.merge = 50
DMR.pct.sig = 0.20

# Expected input format:
# 	(1) Chr
# 	(2) Start
# 	(3) End
# 	(4) unmethylated reads
# 	(5) methylated reads
# 	(6) %methylation

# DSS expects this data structure: The data frame MUST contain following columns in correct order: 
# 	(1) Chromosome number
# 	(2) Genomic coordinates
# 	(3) Read coverage of the position from BS-seq data
# 	(4) Number of reads showing methylation of the position.

message("\nThe main output file will be ", output.filename)
if (minimum.coverage > 0)
{
	message("In any sample, removing site if it has less than ", minimum.coverage, " reads.\n")
}

number.of.input.files = length(input.filenames)
if (number.of.input.files != 6)
{
	message("\nERROR: This code currently requires exactly 2 different samples, with three replicates each.")
	message("       Change the number of input files to 6 or ask BaRC to modify the code.\n")
	q()
}

message("Loading DSS library ... ")
library(DSS)

message("Reading ", number.of.input.files, " input methylation BED files:")

# Create an empty list to hold our methylation data
this.methyl.dataset = list()

# Read one input file at a time and add the data to this.methyl.dataset
for (i in 1:length(input.filenames))
{
	message("  ", input.filenames[i], " ==> ", sample.names[i])
	methyl.data.all = read.delim(input.filenames[i], header = FALSE)
	total.num.reads = methyl.data.all[,4] + methyl.data.all[,5]
	
	# Temporarily make matrix with chr in place of total reads 
	# Use first genome coordinate as coordinate
	methyl.dataset = methyl.data.all[,c(1,2,1,5)]
	colnames(methyl.dataset) = c("chr", "pos", "N", "X")

	# Replace total reads in the third column
	methyl.dataset[,3] = total.num.reads
	
	methyl.dataset = methyl.dataset[methyl.dataset$N >= minimum.coverage,]
	
	# Add this dataset to our big list of methylation data
	this.methyl.dataset[[i]] = methyl.dataset
}
# Combine all of the datasets into one big object
message("Making BSseq object to hold all data ....")
BSobj = makeBSseqData(this.methyl.dataset, sample.names)

message("Analyzing dataset of ", nrow(BSobj), " positions.")
##  DML test

message("\nPerforming statistical tests for differential methylation (SampleA vs SampleB) at each site ...")

if (smoothing == FALSE)
{
	message("Smoothing is turned off.")
} else {
	message("Smoothing is turned on, with a smoothing span of ", smoothing.span)
}

dmlTest <- DMLtest(BSobj, group1=c("SampleA.rep1", "SampleA.rep2", "SampleA.rep3"), 
                          group2=c("SampleB.rep1", "SampleB.rep2", "SampleB.rep3"),
                          smoothing=smoothing, smoothing.span=smoothing.span, ncores=num.cores)


# Filter sites and only run FDR on sites that have at least a preset methylation difference
if (min.meth.difference > 0)
{
	message("To get \"filtered FDRs\", we're ignoring sites that have a methylation fraction difference less than ", min.meth.difference)
	pvals.filtered = rep(NA, nrow(dmlTest))
	for (i in 1:nrow(dmlTest))
	{
		if (abs(dmlTest$diff[i]) >= min.meth.difference)
		{
			pvals.filtered[i] = dmlTest$pval[i]
		}
	}
	fdr.filtered = p.adjust(pvals.filtered, "fdr")
} else
{
	fdr.filtered = dmlTest$pval
}

# Create output BED files
chrs = as.vector(dmlTest[,1])
start.pos = dmlTest[,2]
end.pos = dmlTest[,2] + 1
meth.difference = signif(dmlTest[,5], 4)

# Can't have scientific notation as bedgraph score
output.bedgraph = cbind(chrs, sprintf("%.0f",start.pos),sprintf("%.0f",end.pos), sprintf("%.5f",meth.difference))
# head(output.bedgraph)

message("\nPrinting results ...")

write.table(cbind(dmlTest[,c(1:5,10,11)], fdr.filtered), file=output.filename, sep="\t", quote=F, row.names=F)

bedgraph.output.file = paste(output.filename, ".bedgraph", sep="")
bedgraph.output.file = gsub(".txt.", ".", bedgraph.output.file)
write.table(output.bedgraph, file=bedgraph.output.file, sep="\t", quote=F, row.names=F, col.names=F)

message("\nCalling DMRs ...")

dmrs = callDMR(dmlTest, delta=DMR.delta, p.threshold=DMR.p.threshold, minlen=DMR.minlen, minCG=DMR.minCG, dis.merge=DMR.distance.merge, pct.sig=DMR.pct.sig)

if (exists("dmrs"))
{
	message("Found ", nrow(dmrs), " DMRs.")

	dmr.output.file = paste(output.filename, ".DMRs.txt", sep="")
	dmr.output.file = gsub(".txt.", ".", dmr.output.file)
	write.table(dmrs, file=dmr.output.file, sep="\t", quote=F, row.names=F)

	dmr.bed.output.file = paste(output.filename, ".DMRs.bed", sep="")
	dmr.bed.output.file = gsub(".txt.", ".", dmr.bed.output.file)
	dmr.output.bed = dmrs[,1:3]
	write.table(dmr.output.bed, file=dmr.bed.output.file, sep="\t", quote=F, row.names=F, col.names=F)
} else {
	message("No DMRs found with the given parameters.")
}

message("All done -- for output, see ", output.filename)
message("                            ", bedgraph.output.file)
if (exists("dmrs"))
{
	message("                            ", dmr.output.file)
	message("                            ", dmr.bed.output.file)
}
message()

# Uncomment these commands if you want to save the data structures created during the session
#   (to further explore them interactively, for example)
# RData.file = paste(output.filename, ".RData", sep="")
# RData.file = gsub(".txt.", ".", RData.file)
# message(paste("Saving data structure as ", RData.file, " ..."))
# save(file = RData.file)
# message("All done!")
