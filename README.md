The main steps in enricment analysis 

```
- Load CSV file and keep the max intensity per protein.
- Select proteins with intensity > 0 for each condition.
- Run STRING enrichment analysis for each condition.
- Combine enrichment results by matching term and category.
- Collect all genes involved in each term across conditions.
- Calculate average intensity for each protein per condition.
- Compute DAv (difference in abundance) between conditions.
- Use MANOVA to find terms with significant changes.
- Keep significant terms and add related intensity, DAv, and gene info.
- Save everything to an Excel file.
```
 MSstats Differential Expression Workflow
``` 
- Remove duplicates, keeping the row with the highest intensity per protein.
- Remove rows with all zero or NA values.
- Convert to long format using pivot_longer().
- Add MSstats columns (e.g., Condition, BioReplicate).
- Run dataProcess() to clean and normalize data.
- Set up comparison (e.g., F vs E).
- Run groupComparison() for differential analysis.
- Export results to CSV.
- Create volcano plot to visualize differences.
```
