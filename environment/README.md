# Software

Tested on Windows with R 4.5.2 and Python 3.12.7. Put `Rscript` and `python` on PATH and run from the package root. Arial is required for matching prediction-figure typography. Cross-platform pixel identity has not been tested.

Install the R packages listed in `R_versions.csv` (this installs missing packages from CRAN, not historical versions):

```sh
Rscript -e "x <- read.csv('environment/R_versions.csv'); install.packages(setdiff(x$Package[x$Package != 'R'], rownames(installed.packages())), repos='https://cloud.r-project.org')"
Rscript environment/check.R
```

The version table records the tested direct dependencies; it is not a complete dependency lockfile. The check reports installed versus tested versions and fails if a required package is missing. Prefer the recorded versions for numerical reproduction; `clubSandwich` 0.7.0 is required by the functional workflow. If dependencies live in a custom library, set `R_LIBS_USER` to that library before starting Rscript. The local manuscript validation used a separate library for `clubSandwich`; it is not bundled here.

`RhpcBLASctl` is optional for limiting BLAS threads in cotton fits.

Only HAM requires third-party Python packages:

```sh
python -m pip install -r environment/requirements-ham.txt
```

The projection partitioner and prediction-table exporter use the Python standard library. PowerPoint or compatible presentation software is needed for authored Figure 1/2 editing and export.

Optional small validation checks (synthetic data, no manuscript refitting):

```sh
Rscript validation/check_preprocessing.R output/check_preprocessing
Rscript validation/check_loeo.R
Rscript validation/check_nested_cv.R output/check_nested_cv
```

The last check performs small synthetic fits, tests checkpoint resume, and rejects changed inputs. Repeated checks need new output directories.
