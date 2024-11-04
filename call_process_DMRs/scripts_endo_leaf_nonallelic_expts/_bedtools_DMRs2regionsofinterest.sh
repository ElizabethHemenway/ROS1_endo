#!/bin/bash

###EAH 07/19/24 compare DMR lists to other regions of interest


cd /lab/solexa_gehring/elizabeth/demethylase_em_remap2024/DMRs/make_comparisons/

processed="/lab/solexa_gehring/elizabeth/demethylase_em_remap2024/DMRs/"
comp="/lab/solexa_gehring/elizabeth/demethylase_em_remap2024/DMRs/make_comparisons/"
tocomp=$comp"regions2compare/"


# get genes and TEs of interest

##bedtools window -a $tocomp"genes.bed" -b $processed"r3_v_wt_leaf_allC_hyper_allC.merge.bed" -w 1000 > genes_1kb_r3leaf_hyper.bed
##bedtools window -a $tocomp"TE_fragments.bed" -b $processed"r3_v_wt_leaf_allC_hyper_allC.merge.bed" -w 1000 > TE_fragments_1kb_r3leaf_hyper.bed
##bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"r3_v_wt_leaf_allC_hyper_allC.merge.bed" -w 1000 > TEmerge_1kb_r3leaf_hyper.bed


##bedtools window -a $tocomp"genes.bed" -b $processed"r7_v_wt_leaf_allC_hyper_allC.merge.bed" -w 1000 > genes_1kb_r7leaf_hyper.bed
##bedtools window -a $tocomp"TE_fragments.bed" -b $processed"r7_v_wt_leaf_allC_hyper_allC.merge.bed" -w 1000 > TE_fragments_1kb_r7leaf_hyper.bed
##bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"r7_v_wt_leaf_allC_hyper_allC.merge.bed" -w 1000 > TEmerge_1kb_r7leaf_hyper.bed

#bedtools window -a $tocomp"genes.bed" -b $processed"r3_v_wt_allC_hyper_allC.merge.bed" -w 1000 > genes_1kb_r3_hyper.bed
#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"r3_v_wt_allC_hyper_allC.merge.bed" -w 1000 > TE_fragments_1kb_r3_hyper.bed
#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"r3_v_wt_allC_hyper_allC.merge.bed" -w 1000 > TEmerge_1kb_r3_hyper.bed


#bedtools window -a $tocomp"genes.bed" -b $processed"r7_v_wt_allC_hyper_allC.merge.bed" -w 1000 > genes_1kb_r7_hyper.bed
#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"r7_v_wt_allC_hyper_allC.merge.bed" -w 1000 > TE_fragments_1kb_r7_hyper.bed
#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"r7_v_wt_allC_hyper_allC.merge.bed" -w 1000 > TEmerge_1kb_r7_hyper.bed


# use bedtools nearest to calculate distance from a TE

#bedtools closest -t first -d -a $processed"r3_v_wt_allC_hyper_allC.merge.bed" -b $tocomp"TE_fragments.bed" > r3_v_wt_allC_hyper_distance_to_nearestTEfrag.bed
#bedtools closest -t first -d -a $processed"r3_v_wt_allC_hyper_allC.merge.bed" -b $tocomp"genes.bed" > r3_v_wt_allC_hyper_distance_to_nearestGene.bed

# not ros1 target TEs
#bedtools intersect -v -a $tocomp"TE_fragments.bed" -b $comp"TE_fragments_1kb_r3leaf_hyper.bed" > TE_fragments_NOT_1kb_r3leaf_hyper.bed
#bedtools intersect -v -a $tocomp"TE_fragments.bed" -b $comp"TE_fragments_1kb_r7leaf_hyper.bed" > TE_fragments_NOT_1kb_r7leaf_hyper.bed

#bedtools intersect -v -a $tocomp"TEmerge_named.bed" -b $comp"TEmerge_1kb_r3leaf_hyper.bed" > TEmerge_NOT_1kb_r3leaf_hyper.bed
#bedtools intersect -v -a $tocomp"TEmerge_named.bed" -b $comp"TEmerge_1kb_r7leaf_hyper.bed" > TEmerge_NOT_1kb_r7leaf_hyper.bed

#bedtools intersect -v -a $tocomp"TE_fragments.bed" -b $comp"TE_fragments_1kb_r3_hyper.bed" > TE_fragments_NOT_1kb_r3_hyper.bed
#bedtools intersect -v -a $tocomp"TE_fragments.bed" -b $comp"TE_fragments_1kb_r7_hyper.bed" > TE_fragments_NOT_1kb_r7_hyper.bed

#bedtools intersect -v -a $tocomp"TEmerge_named.bed" -b $comp"TEmerge_1kb_r3_hyper.bed" > TEmerge_NOT_1kb_r3_hyper.bed
#bedtools intersect -v -a $tocomp"TEmerge_named.bed" -b $comp"TEmerge_1kb_r7_hyper.bed" > TEmerge_NOT_1kb_r7_hyper.bed


# not ros1 target genes
#bedtools intersect -v -a $tocomp"genes.bed" -b $comp"genes_1kb_r3leaf_hyper.bed" > genes_NOT_1kb_r3leaf_hyper.bed
#bedtools intersect -v -a $tocomp"genes.bed" -b $comp"genes_1kb_r7leaf_hyper.bed" > genes_NOT_1kb_r7leaf_hyper.bed

#bedtools intersect -v -a $tocomp"genes.bed" -b $comp"genes_1kb_r3_hyper.bed" > genes_NOT_1kb_r3_hyper.bed
#bedtools intersect -v -a $tocomp"genes.bed" -b $comp"genes_1kb_r7_hyper.bed" > genes_NOT_1kb_r7_hyper.bed


# get genes and TEs of interest

processed="/lab/solexa_gehring/elizabeth/demethylase_em_remap2024/DMRs/hypermethylation_limited/"
comp=$processed"make_comparisons/"
tocomp=$processed"regions2compare/"
mkdir -p $comp
cd $comp

bedtools intersect -u -a $processed"r3_CGr_limit50.bed" -b $tocomp"genes.bed" > r3CGr_int_genes_limit50.bed
bedtools intersect -u -a $processed"r3_CGr_limit50.bed" -b $tocomp"TEmerge_named.bed" > r3CGr_int_TEmerge_limit50.bed
bedtools intersect -u -a $processed"r3_CGr_notlimit50.bed" -b $tocomp"genes.bed" > r3CGr_int_genes_notlimit50.bed
bedtools intersect -u -a $processed"r3_CGr_notlimit50.bed" -b $tocomp"TEmerge_named.bed" > r3CGr_int_TEmerge_notlimit50.bed

bedtools intersect -u -a $processed"r7_CGr_limit50.bed" -b $tocomp"genes.bed" > r7CGr_int_genes_limit50.bed
bedtools intersect -u -a $processed"r7_CGr_limit50.bed" -b $tocomp"TEmerge_named.bed" > r7CGr_int_TEmerge_limit50.bed
bedtools intersect -u -a $processed"r7_CGr_notlimit50.bed" -b $tocomp"genes.bed" > r7CGr_int_genes_notlimit50.bed
bedtools intersect -u -a $processed"r7_CGr_notlimit50.bed" -b $tocomp"TEmerge_named.bed" > r7CGr_int_TEmerge_notlimit50.bed
