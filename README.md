## Distributional Data Analysis Uncovers Hundreds of Novel and Heritable Phenomic Features from Temporal Cotton and Maize Drone Imagery

This repository serves as a supplementary data source for our manuscript that introduces the concept of distributional data analysis for agricultural research and will be updated throughout the preprint and publication process.
All analyses mentioned in the paper are included below as RMarkdown documents. Files needed to run the RMarkdown scripts are included in the `R files` folder. Files needed to run the python scripts are included in the `Python files` folder.

Our preprint can be found at [this biorxiv link](https://insert_link_here).

### Notes regarding raw data
In total, this manuscript relied upon several terabytes of files, between raw drone images, orthomosaics, and CSV files that contain the distributional quantile data (the raw quantile data from all G2F locations was a 58 gigabyte file alone, excluding columns for information like genotype, range, row, replicate)! As such, we simply cannot store all of this data on GitHub. The scripts provided below will allow readers to reproduce the analyses in the paper and generate the CSV files that are required to produce the figures from which conclusions in the paper are drawn. For the section that details extracting distributional data from orthomosaics: the orthomosaics, shapefiles, and raw drone images are in the process of being uploaded to the [Data 2 Science website](https://ps2.d2s.org/) hosted by Purdue University and created by [Dr. Jinha Jung](https://engineering.purdue.edu/CCE/People/ptProfile?resource_id=222078).

### Extraction of distributional data from orthomosaics

### Cotton Distributional Data Analysis - R Scripts
"AM" indicates "ANOVA Model", and "HM" indicates "heritability model" (see Table 1 in the paper)
- $\textbf{ANOVA to quantify Genotype × Time interaction with median VIs}$: [AM1 / HM1](https://rpubs.com/ajdesalvio/cotton_maize_anova1)
- $\textbf{ANOVA to quantify heritability and Genotype variance for each time point with median VIs}$: [AM2 / HM2](https://rpubs.com/ajdesalvio/cotton_maize_anova2)
- $\textbf{ANOVA to quantify heritability and Genotype variance for each quantile WITHIN time point}$: [AM3 / HM3](https://rpubs.com/ajdesalvio/cotton_maize_anova3)
- $\textbf{ANOVA to quantify Genotype × Time interaction at each quantile level}$: [AM4 / HM4](https://rpubs.com/ajdesalvio/cotton_maize_anova4)
- $\textbf{ANOVA to quantify Genotype × Quantile interaction at each time point}$: [AM5 / HM5](https://rpubs.com/ajdesalvio/cotton_maize_anova5)

### Distributional Data Analysis - Python Scripts

See the `RMarkdown` folder for the RMarkdown files used to generate each RPubs document.
