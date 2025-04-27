# 🧬 Proteomics Analysis Pipelines: Differential Abundance and Functional Enrichment

This project provides two connected workflows for analyzing proteomics data:
- Identifying differentially abundant proteins across groups
- Performing functional enrichment analysis based on STRING database, followed by MANOVA analysis

---

## 🧬 Differential Protein Abundance Analysis via Protein-Level Aggregation and MANOVA

This pipeline identifies differentially abundant proteins between experimental groups using protein-level intensity aggregation and multivariate analysis of variance (MANOVA).

---

## 📂 Input Requirements
- Proteomics data in Excel format (with a `Description` column containing gene symbols)
- Experimental design: sample grouping information (e.g., control vs treated)

---

## 🔹 Pipeline Steps

### 📁 1. Load Protein-Level Data
- Load the input Excel file containing protein intensity measurements.

### 🔎 2. Extract Protein Names
- Parse the `Description` column to extract gene symbols using regular expressions (`GN=...` pattern).

### 🧹 3. Filter Missing Protein Names
- Remove any rows where a protein name could not be extracted.

### 📊 4. Aggregate Intensities by Protein
- Aggregate data by taking the maximum intensity value for each protein across duplicates.

### 📈 5. Prepare the Expression Matrix
- Keep only intensity measurement columns, skipping metadata columns.
- Set protein names as row names for the expression matrix.

### ⚖️ 6. Normalize Intensities
- Normalize each sample by dividing by the column sum to allow comparability across samples.

### 🔁 7. Transpose the Expression Matrix
- Transpose the matrix so that rows represent samples and columns represent proteins.

### 🧪 8. Create Sample Group Labels
- Define experimental groups (e.g., control vs treated).

### 🧬 9. Perform MANOVA
- Apply multivariate analysis of variance (MANOVA) to detect proteins with differential abundance between groups.

### 📉 10. Extract p-Values
- Extract p-values for each protein from the MANOVA results.

### 🎯 11. Select Significant Proteins
- Identify proteins with p-value < 0.01 as significantly differentially abundant.

### 💾 12. Save Results
- Save the normalized expression values of significant proteins into a CSV file (`MANOVA_significant_proteins.csv`).

---

## 🧬 Functional Enrichment and MANOVA-Based Analysis of Proteomics Data

This pipeline performs functional enrichment analysis for proteins across multiple conditions using the STRING database, followed by multivariate analysis (MANOVA) to identify significant biological processes and pathways.

---

## 📂 Input Requirements
- Proteomics data file (CSV format) with gene names and intensity values
- STRING taxonomy ID for the organism (e.g., `9606` for Homo sapiens)
- A list of condition/group names (e.g., `Control`, `ALlambda`, `ATTR`)

---

## 🔹 Pipeline Steps

### 📁 1. Load Proteomics Data
- Load the input CSV file containing protein intensities across multiple conditions.

### 🧹 2. Aggregate Protein Intensities
- Aggregate data by gene name, keeping the maximum observed intensity per protein.

### ✍️ 3. Clean and Rename Columns
- Adjust column names to ensure correct structure after aggregation.

### 🔁 4. Loop Through Conditions for Enrichment
- For each condition:
  - Select proteins with intensity > 0.
  - Perform functional enrichment using the **STRING database**.

### 🔎 5. Filter Enrichment Categories
- Keep relevant biological categories:
  - Biological Process
  - Molecular Function
  - Cellular Component
  - KEGG, Reactome, WikiPathways

### 🧬 6. Merge Enrichment Results
- Merge all enrichment outputs into one table, replacing NAs with zeros.

### 🧬 7. Merge Unique Gene Lists per Term
- Merge unique gene lists across all conditions for each biological term.

### 📊 8. Prepare Data for Comparative Analysis
- Calculate total protein counts per group.
- Compute mean expression for each group.

### ➗ 9. Calculate Difference Ratios (DAvs_X_Y)
- Calculate normalized difference ratios between each group pair based on mean expression.

### 📈 10. Perform MANOVA
- Perform MANOVA analysis on normalized enrichment data to detect differentially enriched terms.

### 📉 11. Extract Significant Terms
- Select significant biological terms with p-value < 0.05.

### 📝 12. Merge and Annotate Results
- Merge functional annotation, gene lists, and significance results into one table.

### 💾 13. Save Final Results
- Save the final merged table into an Excel file (`enrichment_analysis_results.xlsx`).

---

# 🚀 Example Commands

_(Inside R, or using Rscript if saved as scripts)_

### Run Differential Protein Abundance Analysis:

```r
# Inside R
source("run_protein_abundance_MANOVA.R")

Here is an example to build a dendrogram plot [here](https://github.com/isky20/PLOT_enrichment_TREE2/tree/main).
