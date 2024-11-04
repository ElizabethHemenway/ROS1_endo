#!/bin/bash

###EAH 07/19/24 compare DMR lists to other regions of interest


cd /lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/make_comparisons/

processed="/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/processed/"
comp="/lab/solexa_gehring/elizabeth/allelic_emseq/DMRs_2/make_comparisons/"
tocomp=$comp"regions2compare/"

#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"dme_genes_gehring.bed" -w 1000 > ros1_needed_matVpathypo_1kb_dmegenes.bed
#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"dme_genes_gehring.bed" -w 1000 > WT_and_ros1_matVpathypo_1kb_dmegenes.bed

#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"Col_endo_vs_embryo_CpG_neg_DMRs.sorted.uniq.bed" -w 1000 > WT_and_ros1_matVpathypo_1kb_endohypo.bed
#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"Col_endo_vs_embryo_CpG_neg_DMRs.sorted.uniq.bed" -w 1000 > ros1_needed_matVpathypo_1kb_endohypo.bed

#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"Col_endo_vs_embryo_CpG_pos_DMRs.sorted.uniq.bed" -w 1000 > ros1_needed_matVpathypo_1kb_endohyper.bed
#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"Col_endo_vs_embryo_CpG_pos_DMRs.sorted.uniq.bed" -w 1000 > WT_and_ros1_matVpathypo_1kb_endohyper.bed


# get genes and TEs of interest
#bedtools window -a $tocomp"genes.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > genes_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"genes.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > genes_500bp_ros1_needed_matVpathypo.bed
#bedtools window -a $tocomp"genes.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 500 > genes_500bp_ros1_patonly_hyper.bed

#bedtools window -a $tocomp"genes.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > genes_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"genes.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > genes_1kb_ros1_needed_matVpathypo.bed
#bedtools window -a $tocomp"genes.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 1000 > genes_1kb_ros1_patonly_hyper.bed


#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > TE_fragments_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > TE_fragments_500bp_ros1_needed_matVpathypo.bed
#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 500 > TE_fragments_500bp_ros1_patonly_hyper.bed

#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > TE_fragments_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > TE_fragments_1kb_ros1_needed_matVpathypo.bed
#bedtools window -a $tocomp"TE_fragments.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 1000 > TE_fragments_1kb_ros1_patonly_hyper.bed


#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > TEmerge_named_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > TEmerge_named_500bp_ros1_needed_matVpathypo.bed
#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 500 > TEmerge_named_500bp_ros1_patonly_hyper.bed

#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > TEmerge_named_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > TEmerge_named_1kb_ros1_needed_matVpathypo.bed
#bedtools window -a $tocomp"TEmerge_named.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 1000 > TEmerge_named_1kb_ros1_patonly_hyper.bed



#imprinted genes near DMRs

#bedtools window -a $tocomp"pignatta_picard_megs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > pignatta_picard_megs_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta_picard_megs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > pignatta_picard_megs_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"pignatta_picard_megs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > pignatta_picard_megs_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta_picard_megs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > pignatta_picard_megs_1kb_ros1_needed_matVpathypo.bed


#bedtools window -a $tocomp"pignatta_picard_pegs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > pignatta_picard_pegs_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta_picard_pegs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > pignatta_picard_pegs_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"pignatta_picard_pegs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > pignatta_picard_pegs_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta_picard_pegs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > pignatta_picard_pegs_1kb_ros1_needed_matVpathypo.bed

#just pignatta 2014 genes
#bedtools window -a $tocomp"pignatta2014_megs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > pignatta2014_megs_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta2014_megs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > pignatta2014_megs_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"pignatta2014_megs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > pignatta2014_megs_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta2014_megs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > pignatta2014_megs_1kb_ros1_needed_matVpathypo.bed


#bedtools window -a $tocomp"pignatta2014_pegs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > pignatta2014_pegs_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta2014_pegs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > pignatta2014_pegs_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"pignatta2014_pegs.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > pignatta2014_pegs_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"pignatta2014_pegs.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > pignatta2014_pegs_1kb_ros1_needed_matVpathypo.bed


#Col/ler ISRs
#bedtools window -a $tocomp"Col_Ler_matISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > Col_Ler_matISR_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Ler_matISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > Col_Ler_matISR_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"Col_Ler_matISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > Col_Ler_matISR_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Ler_matISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > Col_Ler_matISR_1kb_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"Col_Ler_patISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > Col_Ler_patISR_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Ler_patISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > Col_Ler_patISR_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"Col_Ler_patISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > Col_Ler_patISR_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Ler_patISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > Col_Ler_patISR_1kb_ros1_needed_matVpathypo.bed


#Col/Cvi ISRs
#bedtools window -a $tocomp"Col_Cvi_matISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > Col_Cvi_matISR_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Cvi_matISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > Col_Cvi_matISR_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"Col_Cvi_matISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > Col_Cvi_matISR_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Cvi_matISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > Col_Cvi_matISR_1kb_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"Col_Cvi_patISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 500 > Col_Cvi_patISR_500bp_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Cvi_patISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 500 > Col_Cvi_patISR_500bp_ros1_needed_matVpathypo.bed

#bedtools window -a $tocomp"Col_Cvi_patISR.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > Col_Cvi_patISR_1kb_WT_and_ros1_matVpathypo.bed
#bedtools window -a $tocomp"Col_Cvi_patISR.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > Col_Cvi_patISR_1kb_ros1_needed_matVpathypo.bed


#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"Col_Cvi_matISR.bed" -w 1000 > ros1_needed_1kb_ColCvimatISR.bed
#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"Col_Cvi_patISR.bed" -w 1000 > ros1_needed_1kb_ColCvipatISR.bed
#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"Col_Cvi_matISR.bed" -w 1000 > dme_dom_1kb_ColCvimatISR.bed
#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"Col_Cvi_patISR.bed" -w 1000 > dme_dom_1kb_ColCvipatISR.bed

#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"Col_Ler_matISR.bed" -w 1000 > ros1_needed_1kb_ColLermatISR.bed
#bedtools window -a $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -b $tocomp"Col_Ler_patISR.bed" -w 1000 > ros1_needed_1kb_ColLerpatISR.bed
#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"Col_Ler_matISR.bed" -w 1000 > dme_dom_1kb_ColLermatISR.bed
#bedtools window -a $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed" -b $tocomp"Col_Ler_patISR.bed" -w 1000 > dme_dom_1kb_ColLerpatISR.bed




bedtools window -a $tocomp"drdd_CG_hyperdmrs_CGonly.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  -w 1000 > drdd_CG_hyperdmrs_CGonly_1kb_WT_and_ros1_matVpathypo.bed
bedtools window -a $tocomp"drdd_CG_hyperdmrs_CGonly.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" -w 1000 > drdd_CG_hyperdmrs_CGonly_1kb_ros1_needed_matVpathypo.bed
bedtools window -a $tocomp"drdd_CG_hyperdmrs_CGonly.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" -w 1000 > drdd_CG_hyperdmrs_CGonly_1kb_ros1_patonly_hyper.bed

bedtools intersect -wb -a $tocomp"drdd_CG_hyperdmrs_CGonly.bed" -b $processed"WT_and_ros1_matVpat_hypo_allC.merge.bed"  > drdd_CG_hyperdmrs_CGonly_int_WT_and_ros1_matVpathypo.bed
bedtools intersect -wb -a $tocomp"drdd_CG_hyperdmrs_CGonly.bed" -b $processed"ros1_only_NEITHERcheck_matVpat_hypo_allC.bed" > drdd_CG_hyperdmrs_CGonly_int_ros1_needed_matVpathypo.bed
bedtools intersect -wb -a $tocomp"drdd_CG_hyperdmrs_CGonly.bed" -b $processed"ros1_v_wt_patonly_hyper_allC.bed" > drdd_CG_hyperdmrs_CGonly_int_ros1_patonly_hyper.bed







