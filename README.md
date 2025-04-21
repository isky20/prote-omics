The main steps in enricment analysis 

```
1- Load the input CSV file and aggregate intensities by protein using the maximum value.

2- Run STRING enrichment analysis on protein (psm per conditions >0). 

3- Merge enrichment results across all conditions by term, description, and category.

4- Create a unique combined gene list for each enriched term across conditions.

5- Calculate the mean intensity values for each experimental condition.

6- Compute Differential Abundance values (DAvs) between each pair of conditions.

7- Perform MANOVA to identify enriched terms significantly differing across groups.

8- Filter significant terms and merge them with intensity, DAv, category, and gene data.

9- Export the final combined results as an Excel file. 
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
