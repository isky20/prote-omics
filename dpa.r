# ---------------------------------------------------------------
# Title: Differential Protein Abundance Analysis via Protein-Level Aggregation and MANOVA
# ---------------------------------------------------------------

# Load required libraries
library(readxl)
library(dplyr)
library(stringr)
library(tibble)

# Step 1: Load the input Excel file containing protein-level data
rdata <- read_excel("your_proteomics_file.xlsx")

# Step 2: Extract protein names (gene symbols) from the 'Description' column using regex
rdata$Protein <- str_extract(rdata$Description, "(?<=GN=)[^ ]+")

# Step 3: Filter out rows where protein names could not be extracted (i.e., NA)
rdata <- rdata[!is.na(rdata$Protein), ]

# Step 4: Aggregate data by protein using the maximum intensity across duplicates
# This ensures one row per protein with the highest observed values
protein_agg <- aggregate(rdata[, sapply(rdata, is.numeric)],
                         list(rdata$Protein), FUN = max)

# Step 5: Prepare the expression matrix by selecting only intensity columns
# The first 4 columns are assumed to be metadata and are skipped
expr_matrix <- protein_agg[, 5:ncol(protein_agg)]

# Step 6: Set protein names as row names for the expression matrix
rownames(expr_matrix) <- protein_agg$Group.1

# Step 7: Normalize the data by dividing each column by its column sum
# This makes sample intensities comparable
expr_norm <- apply(expr_matrix, 2, function(x) x / sum(x))

# Step 8: Transpose the matrix so that rows = samples, columns = proteins
expr_transposed <- t(expr_norm)

# Step 9: Create group labels for the samples (e.g., control vs treated)
# Here: First 3 samples are control, next 3 are treated
group <- data.frame(group = rep(c("control", "treat"), each = 3))

# Step 10: Run MANOVA to find proteins differentially abundant between groups
fit.manova <- manova(expr_transposed ~ group$group)

# Step 11: Extract p-values from the MANOVA result
p_values <- summary.aov(fit.manova)

# Step 12: Identify significant proteins with p-value < 0.01
p.value <- unlist(lapply(p_values, function(x) x[[1]][["Pr(>F)"]][1]))
sig_proteins <- names(p.value[p.value < 0.01])

# Step 13: Extract normalized expression values for significant proteins only
res <- expr_norm[sig_proteins, ]

# Step 14: Save the result to a CSV file
write.csv(res, "MANOVA_significant_proteins.csv")
