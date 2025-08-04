#!/usr/bin/env Rscript

### written by: Elizabeth Hemenway
### 06/17/2025 - filter DME TEs by methylation difference in ros1 endosperm

setwd("./")
Colpat1=$path'dmedom_tes_CG_in2kbout2kb_C24xCol_1_Col_spiked_CG_min5.bed_mat.txt'
Colpat2=$path'dmedom_tes_CG_in2kbout2kb_C24xCol_2_Col_spiked_CG_min5.bed_mat.txt' 
Colpat3=$path'dmedom_tes_CG_in2kbout2kb_C24xCol_3_Col_spiked_CG_min5.bed_mat.txt' 
r3pat1=$path'dmedom_tes_CG_in2kbout2kb_r1xr3_1_Col_spiked_CG_min5.bed_mat.txt' 
r3pat2=$path'dmedom_tes_CG_in2kbout2kb_r1xr3_2_Col_spiked_CG_min5.bed_mat.txt' 
r3pat3=$path'dmedom_tes_CG_in2kbout2kb_r1xr3_3_Col_spiked_CG_min5.bed_mat.txt'

Col1=$path'dmedom_tes_CG_in2kbout2kb_wt_1_Col_spiked_CG_min5.bed_mat.txt'
Col2=$path'dmedom_tes_CG_in2kbout2kb_wt_2_Col_spiked_CG_min5.bed_mat.txt' 
Col3=$path'dmedom_tes_CG_in2kbout2kb_wt_3_Col_spiked_CG_min5.bed_mat.txt' 
r31=$path'dmedom_tes_CG_in2kbout2kb_r3_1_Col_spiked_CG_min5.bed_mat.txt'
r32=$path'dmedom_tes_CG_in2kbout2kb_r3_2_Col_spiked_CG_min5.bed_mat.txt' 
r33=$path'dmedom_tes_CG_in2kbout2kb_r3_3_Col_spiked_CG_min5.bed_mat.txt' 


cd $regions'ends/'
$regions'make_diff_hmaps.r' $regions'ends/' dmedom_TEs_CG $r31 $r32 $r33 $Col1 $Col2 $Col3 30



r3_TEs_endo <- read.csv("./data/TE_fragments_1kb_r3_hyper_TE_24ntsRNA_endo24nt_mat.txt", header=TRUE, sep="\t")
r3_TEs_em <- read.csv("./data/TE_fragments_1kb_r3_hyper_TE_24ntsRNA_embryo24nt_mat.txt", header=TRUE, sep="\t")

# Now, ensure that only numeric columns are included in the matrix - but I want to keep NAs as place holder? 
lastcol=ncol(r3_TEs_endo)
r3_TEs_endo[,2:lastcol] <- sapply(r3_TEs_endo[,2:lastcol], as.numeric)

lastcol=ncol(r3_TEs_em)
r3_TEs_em[,2:lastcol] <- sapply(r3_TEs_em[,2:lastcol], as.numeric)

#print(r3_TEs_em)

compare_flank_averages <- function(df1, df2, id_col = 1, min_diff = 0.01) {
  col_names <- names(df1)
  id_col_name <- col_names[id_col]

  if (!"feature_start" %in% col_names) stop("Column 'feature_start' not found.")
  if (!"feature_end" %in% col_names) stop("Column 'feature_end' not found.")

  start_idx <- which(col_names == "feature_start")[1]
  end_idx <- which(col_names == "feature_end")[1]

  cols_before_start <- col_names[(id_col + 1):(start_idx - 1)]
  cols_after_end <- col_names[(end_idx + 1):length(col_names)]

  # Compute averages
  df1$avg_before <- rowMeans(df1[ , cols_before_start], na.rm = TRUE)
  df1$avg_after  <- rowMeans(df1[ , cols_after_end], na.rm = TRUE)
  df2$avg_before <- rowMeans(df2[ , cols_before_start], na.rm = TRUE)
  df2$avg_after  <- rowMeans(df2[ , cols_after_end], na.rm = TRUE)

  # Differences
  diff_before <- df1$avg_before - df2$avg_before
  diff_after  <- df1$avg_after - df2$avg_after

  # Prepare subsets
  df1_sub <- df1[ , c(id_col_name, "avg_before", "avg_after")]
  df2_sub <- df2[ , c(id_col_name, "avg_before", "avg_after")]

  merged <- merge(df1_sub, df2_sub, by = id_col_name, suffixes = c("_df1", "_df2"))

  # Add difference columns
  merged$diff_before <- diff_before[match(merged[[id_col_name]], df1[[id_col_name]])]
  merged$diff_after  <- diff_after[match(merged[[id_col_name]], df1[[id_col_name]])]

  # Logical vector for rows meeting the condition
  condition_met <- (merged$diff_before > min_diff) | (merged$diff_after > min_diff)

  matched <- merged[condition_met, ]
  unmatched <- merged[!condition_met, ]

  return(list(matched = matched, unmatched = unmatched))
}
results <- compare_flank_averages(r3_TEs_endo, r3_TEs_em)

write.csv(results$matched, "mindiff_rows.csv", row.names = FALSE)
write.csv(results$unmatched, "notdiff_rows.csv", row.names = FALSE)



compare_flank_sums <- function(df1, df2, id_col = 1, min_diff = 0.5) {
  col_names <- names(df1)
  id_col_name <- col_names[id_col]

  if (!"feature_start" %in% col_names) stop("Column 'feature_start' not found.")
  if (!"feature_end" %in% col_names) stop("Column 'feature_end' not found.")

  start_idx <- which(col_names == "feature_start")[1]
  end_idx <- which(col_names == "feature_end")[1]

  cols_before_start <- col_names[(id_col + 1):(start_idx - 1)]
  cols_after_end    <- col_names[(end_idx + 1):length(col_names)]

  # Compute sums
  df1$sum_before <- rowSums(df1[ , cols_before_start], na.rm = TRUE)
  df1$sum_after  <- rowSums(df1[ , cols_after_end], na.rm = TRUE)
  df2$sum_before <- rowSums(df2[ , cols_before_start], na.rm = TRUE)
  df2$sum_after  <- rowSums(df2[ , cols_after_end], na.rm = TRUE)

  # Differences
  diff_before <- df1$sum_before - df2$sum_before
  diff_after  <- df1$sum_after  - df2$sum_after

  # Subset to ID and sums
  df1_sub <- df1[ , c(id_col_name, "sum_before", "sum_after")]
  df2_sub <- df2[ , c(id_col_name, "sum_before", "sum_after")]

  # Merge on ID
  merged <- merge(df1_sub, df2_sub, by = id_col_name, suffixes = c("_df1", "_df2"))

  # Add difference columns
  merged$diff_before <- diff_before[match(merged[[id_col_name]], df1[[id_col_name]])]
  merged$diff_after  <- diff_after[match(merged[[id_col_name]], df1[[id_col_name]])]

  # Logical vector for rows meeting the condition
  condition_met <- (merged$diff_before > min_diff) | (merged$diff_after > min_diff)

  matched   <- merged[condition_met, ]
  unmatched <- merged[!condition_met, ]

  return(list(matched = matched, unmatched = unmatched))
}
results <- compare_flank_sums(r3_TEs_endo, r3_TEs_em)

write.csv(results$matched, "mindiff_sumrows.csv", row.names = FALSE)
write.csv(results$unmatched, "notdiff_sumrows.csv", row.names = FALSE)


compare_flank_diff_vs_sd <- function(df1, df2, id_col = 1, sd_factor = 1) {
  col_names <- names(df1)
  id_col_name <- col_names[id_col]

  if (!"feature_start" %in% col_names) stop("Column 'feature_start' not found.")
  if (!"feature_end" %in% col_names) stop("Column 'feature_end' not found.")

  start_idx <- which(col_names == "feature_start")[1]
  end_idx   <- which(col_names == "feature_end")[1]

  cols_before <- col_names[(id_col + 1):(start_idx - 1)]
  cols_after  <- col_names[(end_idx + 1):length(col_names)]

  # Compute row sums
  df1$sum_before <- rowSums(df1[ , cols_before], na.rm = TRUE)
  df1$sum_after  <- rowSums(df1[ , cols_after], na.rm = TRUE)
  df2$sum_before <- rowSums(df2[ , cols_before], na.rm = TRUE)
  df2$sum_after  <- rowSums(df2[ , cols_after], na.rm = TRUE)

  # Subset for merging
  df1_sub <- df1[ , c(id_col_name, "sum_before", "sum_after")]
  df2_sub <- df2[ , c(id_col_name, "sum_before", "sum_after")]

  # Merge on ID
  merged <- merge(df1_sub, df2_sub, by = id_col_name, suffixes = c("_df1", "_df2"))

  # Compute differences
  merged$diff_before <- merged$sum_before_df1 - merged$sum_before_df2
  merged$diff_after  <- merged$sum_after_df1  - merged$sum_after_df2

  # Compute mean and standard deviation of differences
  mean_diff_before <- mean(merged$diff_before, na.rm = TRUE)
  mean_diff_after  <- mean(merged$diff_after, na.rm = TRUE)

  sd_diff_before <- sd(merged$diff_before, na.rm = TRUE)
  sd_diff_after  <- sd(merged$diff_after, na.rm = TRUE)

  # Thresholds
  threshold_before <- mean_diff_before + sd_factor * sd_diff_before
  threshold_after  <- mean_diff_after  + sd_factor * sd_diff_after

  # Filter based on thresholds
  condition_met <- (merged$diff_before > threshold_before) |
                   (merged$diff_after  > threshold_after)

  matched   <- merged[condition_met, ]
  unmatched <- merged[!condition_met, ]

  return(list(matched = matched, unmatched = unmatched))
}

compare_flank_diff_vs_sd <- function(df1, df2, id_col = 1, sd_factor = 1) {
  col_names <- names(df1)
  id_col_name <- col_names[id_col]

  if (!"feature_start" %in% col_names) stop("Column 'feature_start' not found.")
  if (!"feature_end" %in% col_names) stop("Column 'feature_end' not found.")

  start_idx <- which(col_names == "feature_start")[1]
  end_idx   <- which(col_names == "feature_end")[1]

  cols_before <- col_names[(id_col + 1):(start_idx - 1)]
  cols_after  <- col_names[(end_idx + 1):length(col_names)]

  # Compute row sums
  df1$sum_before <- rowSums(df1[ , cols_before], na.rm = TRUE)
  df1$sum_after  <- rowSums(df1[ , cols_after], na.rm = TRUE)
  df2$sum_before <- rowSums(df2[ , cols_before], na.rm = TRUE)
  df2$sum_after  <- rowSums(df2[ , cols_after], na.rm = TRUE)

  # Subset for merging
  df1_sub <- df1[ , c(id_col_name, "sum_before", "sum_after")]
  df2_sub <- df2[ , c(id_col_name, "sum_before", "sum_after")]

  # Merge on ID
  merged <- merge(df1_sub, df2_sub, by = id_col_name, suffixes = c("_df1", "_df2"))

  # Compute differences
  merged$diff_before <- merged$sum_before_df1 - merged$sum_before_df2
  merged$diff_after  <- merged$sum_after_df1  - merged$sum_after_df2

  # Compute mean and standard deviation of differences
  mean_diff_before <- mean(merged$diff_before, na.rm = TRUE)
  mean_diff_after  <- mean(merged$diff_after, na.rm = TRUE)

  sd_diff_before <- sd(merged$diff_before, na.rm = TRUE)
  sd_diff_after  <- sd(merged$diff_after, na.rm = TRUE)

  # Thresholds
  threshold_before <- mean_diff_before + sd_factor * sd_diff_before
  threshold_after  <- mean_diff_after  + sd_factor * sd_diff_after

  # Filter based on thresholds
  condition_met <- (merged$diff_before > threshold_before) |
                   (merged$diff_after  > threshold_after)

  matched   <- merged[condition_met, ]
  unmatched <- merged[!condition_met, ]

  return(list(matched = matched, unmatched = unmatched))
}
results <- compare_flank_diff_vs_sd(r3_TEs_endo, r3_TEs_em)

write.csv(results$matched, "mindiff_1SD_sumrows.csv", row.names = FALSE)
write.csv(results$unmatched, "notdiff_1SD_sumrows.csv", row.names = FALSE)

plot_flank_diff_histogram <- function(df1, df2, id_col = 1, bins = 30) {
  library(ggplot2)

  col_names <- names(df1)
  id_col_name <- col_names[id_col]

  if (!"feature_start" %in% col_names) stop("Column 'feature_start' not found.")
  if (!"feature_end" %in% col_names) stop("Column 'feature_end' not found.")

  start_idx <- which(col_names == "feature_start")[1]
  end_idx   <- which(col_names == "feature_end")[1]

  cols_before <- col_names[(id_col + 1):(start_idx - 1)]
  cols_after  <- col_names[(end_idx + 1):length(col_names)]

  # Compute row sums
  df1$sum_before <- rowSums(df1[ , cols_before], na.rm = TRUE)
  df1$sum_after  <- rowSums(df1[ , cols_after], na.rm = TRUE)
  df2$sum_before <- rowSums(df2[ , cols_before], na.rm = TRUE)
  df2$sum_after  <- rowSums(df2[ , cols_after], na.rm = TRUE)

  # Subset for merging
  df1_sub <- df1[ , c(id_col_name, "sum_before", "sum_after")]
  df2_sub <- df2[ , c(id_col_name, "sum_before", "sum_after")]

  # Merge on ID
  merged <- merge(df1_sub, df2_sub, by = id_col_name, suffixes = c("_df1", "_df2"))

  # Compute differences
  merged$diff_before <- merged$sum_before_df1 - merged$sum_before_df2
  merged$diff_after  <- merged$sum_after_df1  - merged$sum_after_df2

  # Plot helper
  plot_diff_histogram <- function(diff_vector, label) {
    mu <- mean(diff_vector, na.rm = TRUE)
    sd_val <- sd(diff_vector, na.rm = TRUE)

    ggplot(data.frame(diff = diff_vector), aes(x = diff)) +
      geom_histogram(bins = bins, fill = "grey70", color = "black") +
      geom_vline(xintercept = mu, color = "blue", linetype = "solid", size = 1.2) +
      geom_vline(xintercept = c(mu + 0.5 * sd_val, mu - 0.5 * sd_val),
                 color = "orange", linetype = "dashed", size = 1) +
      geom_vline(xintercept = c(mu + 1 * sd_val, mu - 1 * sd_val),
                 color = "red", linetype = "dashed", size = 1) +
      labs(title = paste("Histogram of", label, "Differences"),
           x = "Difference",
           y = "Count") +
      theme_minimal()
  }

  # Plot both
  p1 <- plot_diff_histogram(merged$diff_before, "Before")
  p2 <- plot_diff_histogram(merged$diff_after, "After")

  return(list(before_plot = p1, after_plot = p2))
}
plots <- plot_flank_diff_histogram(r3_TEs_endo, r3_TEs_em)
print(plots$before_plot)
print(plots$after_plot)
