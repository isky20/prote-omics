the main steps in enricment analysis 

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
