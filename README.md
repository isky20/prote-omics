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
Step 1: Remove duplicate protein IDs by keeping the row with the highest intensity per protein.

Step 3: Remove rows where all intensity values are zero or NA.

Step 4: Reshape the data from wide to long format using pivot_longer().

Step 5: Add required MSstats columns.

Step 6: Run dataProcess() to normalize and summarize the data.

Step 7: Define the comparison matrix for conditions (e.g., F vs E).

Step 8: Run groupComparison() to perform differential analysis.

Step 9: Save the results to CSV and generate a volcano plot.
```
