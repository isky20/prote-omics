# Load required packages
library(MSstats)
library(tidyverse)

# Step 1: Load your protein-level data
data_raw <- read.csv("protein_counts.csv")

# Remove duplicated Pro_id names (like repeated C4A)
data_clean <- data_raw[!duplicated(data_raw$Pro_id), ]
# Remove rows with all zeros or NAs across replicates
data_clean <- data_clean %>%
  filter(if_any(-Pro_id, ~ !is.na(.) & . != 0))  # keep rows with at least one non-zero


data_clean <- data_clean %>%
  group_by(Pro_id) %>%
  summarise(across(where(is.numeric), max, na.rm = TRUE))

# Step 2: Reshape data to long format
data_long <- data_clean %>%
  pivot_longer(cols = -Pro_id, names_to = "Replicate", values_to = "Intensity") %>%
  mutate(
    Condition = ifelse(grepl("^E_", Replicate), "E", "F"),
    BioReplicate = Replicate,
    ProteinName = Pro_id
  )

# Step 3: Add dummy but unique values required by MSstats
msstats_input_fixed <- data_long %>%
  mutate(
    PeptideSequence = paste0("Peptide_", ProteinName),                 # unique dummy peptide per protein
    PeptideModifiedSequence = paste0("Peptide_", ProteinName),
    PrecursorCharge = 2,
    FragmentIon = "NA",
    ProductCharge = 1,
    IsotopeLabelType = "L",
    Run = BioReplicate
  ) %>%
  dplyr::select(ProteinName, PeptideSequence, PeptideModifiedSequence,
                PrecursorCharge, FragmentIon, ProductCharge,
                IsotopeLabelType, Condition, BioReplicate, Run, Intensity)

# Step 4: Run MSstats data processing
processed <- dataProcess(msstats_input_fixed,
                         normalization = "equalizeMedians",
                         summaryMethod = "TMP",
                         censoredInt = "NA")

# Step 5: Define comparison (F vs E)
comparison <- matrix(c(-1, 1), nrow = 1)
colnames(comparison) <- c("E", "F")
rownames(comparison) <- "F_vs_E"
# Step 6: Differential expression
results <- groupComparison(contrast.matrix = comparison,
                           data = processed)

# Step 7: View and save results
head(results$ComparisonResult)
write.csv(results$ComparisonResult, "Differential_Proteins_MSstats.csv", row.names = FALSE)

# Optional: Volcano Plot
groupComparisonPlots(data = results$ComparisonResult, type = "VolcanoPlot")
