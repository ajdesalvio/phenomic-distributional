# Software

We tested the workflows on Windows with R 4.5.2 and Python 3.12.7. Add `Rscript` and `python` to your PATH, then run the commands from the repository root. Install Arial to match the prediction figures. Pixel comparisons were done on Windows.

The first command installs missing R packages from CRAN. The second compares your installed versions with those in `R_versions.csv`.

```sh
Rscript -e "x <- read.csv('environment/R_versions.csv'); install.packages(setdiff(x$Package[x$Package != 'R'], rownames(installed.packages())), repos='https://cloud.r-project.org')"
Rscript environment/check.R
```

The table records the versions we tested for direct dependencies. Use those versions when reproducing the numerical results. The check fails if a required package is missing, and the cotton functional workflow specifically requires `clubSandwich` 0.7.0. If you install packages in a custom library, set `R_LIBS_USER` before opening R. We used a separate local library for `clubSandwich` when validating the manuscript results.

`RhpcBLASctl` is optional for limiting BLAS threads in cotton fits.

Only HAM requires third-party Python packages:

```sh
python -m pip install -r environment/requirements-ham.txt
```

The projection partitioner and prediction-table exporter use the Python standard library. Use PowerPoint or compatible software to edit and export Figures 1 and 2.

Run these optional checks on synthetic data.

```sh
Rscript validation/check_preprocessing.R output/check_preprocessing
Rscript validation/check_loeo.R
Rscript validation/check_nested_cv.R output/check_nested_cv
```

The nested CV check fits small synthetic models, resumes from checkpoints, and confirms that changed inputs are rejected. Choose new output directories when repeating the checks.
