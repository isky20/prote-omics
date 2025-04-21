# Load required libraries
library(rbioapi)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(openxlsx)

# Main function to perform enrichment analysis and generate comparative results
process_enrichment <- function(obj, spid, int, output_file) {
  
  # STEP 1: Load the CSV file into a dataframe
  df <- read.csv(obj)
  
  # STEP 2: Aggregate rows by Gene name and take the max intensity per group
  df <- aggregate(df[], list(df$Gene), FUN = max)
  
  # STEP 3: Clean up column names after aggregation
  df["Gene"] <- NULL
  colnames(df)[1] <- "Gene"
  
  # STEP 4: Initialize a list to store enrichment results for each condition
  df_list <- list()
  
  # STEP 5: Loop through each sample/condition column (excluding 'Gene')
  for (i in colnames(df)[2:length(colnames(df))]) {
    # Select genes with intensity > 0
    protlist <- df$Gene[df[[i]] > 0]
    
    # Perform functional enrichment using STRING API
    enrichment <- rbioapi::rba_string_enrichment(protlist, spid, split_df = FALSE)
    
    # Filter for selected pathway/ontology categories
    categories_to_match <- c("Process", "Function", "Component", "KEGG", "RCTM", "WikiPathways")
    subset_df <- enrichment[enrichment$category %in% categories_to_match,
                            c("term", "description", "category", "preferredNames", "number_of_genes")]
    
    # Rename columns to reflect the current condition
    colnames(subset_df)[colnames(subset_df) == "number_of_genes"] <- paste0("count_", i)
    colnames(subset_df)[colnames(subset_df) == "preferredNames"] <- paste0("genes_", i)
    
    # Store the processed enrichment result in the list
    df_list[[i]] <- subset_df
  }
  
  # STEP 6: Merge all enrichment dataframes into one
  final_df <- Reduce(function(x, y) merge(x, y, by = c("term", "description", "category"), all = TRUE), df_list)
  final_df[is.na(final_df)] <- 0  # Replace NAs with 0
  
  # STEP 7: Process gene lists for each term
  gene_columns <- grep("^(term|genes_)", names(final_df), value = TRUE)
  gene_df <- final_df[, gene_columns]
  
  # Helper function to combine unique gene names across columns
  find_unique <- function(lists) {
    combined_list <- unlist(strsplit(unlist(lists), ","))
    combined_list <- combined_list[combined_list != "0"]
    return(paste(unique(combined_list), collapse = ","))
  }
  
  # Apply function to rows to generate a unique list of genes
  gene_df$geneid <- apply(gene_df[, grepl("genes_", names(gene_df))], 1, find_unique)
  
  # Keep only relevant gene info
  gene_df <- gene_df[, c("term", "geneid")]
  
  # STEP 8: Prepare count data (excluding gene columns)
  counts_df <- final_df[, !names(final_df) %in% grep("genes_", names(final_df), value = TRUE)]
  
  # Calculate total counts per condition (optional)
  numeric_cols <- sapply(counts_df, is.numeric)
  wsum.col.of.subjs <- colSums(counts_df[, numeric_cols])
  
  # Extract numeric part (counts per condition)
  psm <- counts_df[, 4:ncol(counts_df)]
  Data.psm <- as.data.frame(psm)
  
  # STEP 9: Compute mean expression per group (e.g., Control, ALlambda, ATTR)
  mean_list <- list()
  for (i in int) {
    mean_list[[i]] <- rowMeans(Data.psm[, grepl(i, names(Data.psm))])
  }
  mean_df <- as.data.frame(mean_list)
  Data.psm <- cbind(Data.psm, mean_df)
  
  # STEP 10: Compute difference ratios between each group pair (DAvs_X_Y)
  combinations <- t(combn(colnames(mean_df), 2))
  sub_dataframes <- list()
  
  DAvs_X_Y <- function(row) {
    if ((row[1] + row[2]) != 0) {
      return(((row[1] - row[2]) / (row[1] + row[2])) / 0.5)
    } else {
      return(0)
    }
  }
  
  for (i in 1:nrow(combinations)) {
    col1 <- combinations[i, 1]
    col2 <- combinations[i, 2]
    h <- paste0("DAvs_", col1, "_", col2)
    sub_df <- mean_df[, c(col1, col2)]
    sub_dataframes[[h]] <- apply(sub_df, 1, DAvs_X_Y)
  }
  
  sub_dataframes_df <- as.data.frame(sub_dataframes)
  Data.psm <- cbind(Data.psm, sub_dataframes_df)
  
  # STEP 11: Perform MANOVA test
  psm.norm.t <- t(psm)
  
  counts <- sapply(int, function(subword) sum(grepl(subword, colnames(df))))
  labels <- rep(int, counts)
  final_psm_norm_t <- data.frame(psm.norm.t, labels)
  
  dependent_vars <- as.matrix(final_psm_norm_t[, -ncol(final_psm_norm_t)])
  independent_var <- as.factor(final_psm_norm_t$labels)
  
  manova_model <- manova(dependent_vars ~ independent_var, data = final_psm_norm_t)
  summary <- summary.aov(manova_model)
  
  # Extract p-values from MANOVA result
  output_LDA <- data.frame()
  for (i in 1:length(summary)) {
    output_LDA <- rbind(output_LDA, summary[[i]][["Pr(>F)"]][1])
  }
  
  output_LDA$description <- final_df$description
  output_LDA$term <- final_df$term
  colnames(output_LDA)[1] <- "Pvalue"
  
  # STEP 12: Keep only significant terms (p < 0.05)
  selected_sig <- output_LDA[output_LDA$Pvalue < 0.05, ]
  
  # Merge significant terms with the full data
  Data.psm$description <- final_df$description
  Data.psm$term <- final_df$term
  merged_df <- merge(Data.psm, selected_sig, by = c("description", "term"), all.x = TRUE)
  
  # Add back category and merge with gene info
  the_data <- merge(merged_df, final_df[, c("description", "category", "term")], by = c("description", "term"))
  no.na <- na.omit(the_data)
  df <- subset(no.na, !duplicated(no.na))
  
  # Reorder columns and merge final gene list
  new_column_order <- c("category", names(df)[!names(df) %in% "category"])
  thedf <- merge(df[, new_column_order], gene_df, by = "term")
  
  # STEP 13: Save results to Excel
  write.xlsx(thedf, output_file, overwrite = TRUE)
  
  return(thedf)
}

# Example input to run the function
obj <- "test3.csv"                  # Input file
spid <- 9606                       # Species ID for Homo sapiens
int <- c("Control", "ALlambda", "ATTR")  # Experimental conditions
output_file <- "enrichment_analysis_results_test3.xlsx"

# Run the function
result <- process_enrichment(obj, spid, int, output_file)
