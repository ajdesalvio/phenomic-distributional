## Distributional Data Analysis Uncovers Hundreds of Novel and Heritable Phenomic Features from Temporal Cotton and Maize Drone Imagery

This repository serves as a supplementary data source for our manuscript that introduces the concept of distributional data analysis for agricultural research and will be updated throughout the preprint and publication process.
All analyses mentioned in the paper are included below as RMarkdown documents. Files needed to run the RMarkdown scripts are included in the `R_files` folder. Files needed to run the python scripts are included in the `Python_files` folder.

Our preprint can be found at [this biorxiv link](https://www.biorxiv.org/content/10.1101/2025.09.05.674557v1).

### Notes regarding raw data
In total, this manuscript relied upon several terabytes of files, between raw drone images, orthomosaics, and CSV files that contain the distributional quantile data (the raw quantile data from all G2F locations was a 58 gigabyte file alone, excluding columns for information like genotype, range, row, replicate)! As such, we simply cannot store all of this data on GitHub. The scripts provided below will allow readers to reproduce the analyses in the paper and generate the CSV files that are required to produce the figures from which conclusions in the paper are drawn. For the section that details extracting distributional data from orthomosaics: the orthomosaics, shapefiles, and raw drone images are in the process of being uploaded to the [Data 2 Science website](https://ps2.d2s.org/) hosted by Purdue University and created by [Dr. Jinha Jung](https://engineering.purdue.edu/CCE/People/ptProfile?resource_id=222078).

## Menu
[1 - Python setup - creating a conda virtual environment from a .yml file](#installation)

[2 - Shapefile creation from a CSV file](#p1)

[3 - Image extraction from orthomosaics](#p2)

[4 - Soil masking](#p3)

[5 - Extraction of standard VIs and distributional data](#p4)

[6 - Mixed function-on-scalar distributional regression model for phenomic data](#p5)
     
[7 - Yield prediction and feature importance quantificaiton with maize distributional data](#p6)
     
[8 - Correlation analysis of maize quantile phenomic data](#p7)

[9 - Visualization of raw data](#p8)

[10 - Calculation and visualization of heritability](#p9)

<br />

<div id="installation" />

### 1 - Python setup - creating a conda virtual environment from a .yml file
All of the libraries required to use our python scripts are included in the file `Python_Phenomics.yaml` stored in the `Python_files` directory. First, download Anaconda Navigator if you haven't already done so, and then launch Anaconda Prompt. From the console, navigate to the directory where you downloaded the yaml file and set up an Anaconda virtual environment using the file. An example (Windows) would look like this:

`cd C:\Users\aaron.desalvio\Downloads`

`conda env create -f Python_Phenomics.yaml`

Once the installation is complete, activate the environment. If you want to use the Spyder IDE, simply type `spyder` after activating the environment and it will launch.

`conda activate Python_Phenomics`

<br />

<div id="p1" />

### 2 - Shapefile creation from a CSV file
To create a shapefile, use the `create_shapefile_from_csv()` function within the `Create_Shapefile.py` script. A CSV file with unique plot IDs (or single plant IDs) must be provided to the function for it to work. An example usage of the function is as follows:

```python
# Path to CSV
csv_path = 'C:/Users/aaron.desalvio/Documents/TESTING/Python_Phenomics/RGB_G2F_2024'
output_path = 'C:/Users/aaron.desalvio/Documents/TESTING/Python_Phenomics/RGB_G2F_2024/Shapefiles'

gdf = create_shapefile_from_csv(
    csv_file=os.path.join(csv_path, 'CS24_G2F_Plots.csv'),
    shapefile_output=os.path.join(output_path, '20240501_G2F_V2.shp'),
    A=(746593.958,3381859.705),
    B=(746723.366,3381741.565),
    ft_toggle=0.3048,
    rangelength=25,
    rangebuffer=2,
    rowwidth=5.04,
    rowbuffer=0.3,
    utm_epsg=32614,
    plot_id_colname='Range_Row',
    range_colname='Range',
    row_colname='Row',
    remap_rows=True,
    shift_rows_to_one=True
)
```
`csv_path` and `output_path` indicate the paths where the initial CSV file is stored and the location where you'd like to save the resulting shapefile, respectively.
"A" and "B" indcate the coordinates of the "AB line", which provides a line that the function will use to specify the angle of each component of the shapefile. Examples of where to obtain the AB coordinates in an orthomosaic could look like:


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
