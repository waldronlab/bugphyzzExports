Hourly export status: ![hourly export](https://github.com/waldronlab/bugphyzzExports/actions/workflows/export-bugphyzz.yml/badge.svg)

# bugphyzzExports

This repository contains the devel
[bugphyzz](https://github.com/waldronlab/bugphyzz) signature files and a full
data dump. The devel files are generated weekly.

The release version of these files will soon be available on Zenodo.

## File creation

If desired, you can generate the signatures and full dump file.

### Requirements

In addition to installing this repository, you will need to install

* [BiocFileCache](https://www.bioconductor.org/packages/BiocFileCache)
* [bugphyzz](https://github.com/waldronlab/bugphyzz)
* [castor](https://cran.r-project.org/web/packages/castor/)
* [logr](https://cran.r-project.org/web/packages/logr/)
* [phytools](https://cran.r-project.org/web/packages/phytools/)
* [purrr](https://cran.r-project.org/web/packages/purrr)
* [readr](https://readr.tidyverse.org/)
* [tibble](https://cran.r-project.org/web/packages/tibble/)
* [tidyr](https://cran.r-project.org/web/packages/tidyr/)

Install using `BiocManager`:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
dependencies <- c("BiocFileCache",
                  "waldronlab/bugphyzz",
                  "castor",
                  "logr",
                  "phytools",
                  "purrr",
                  "readr",
                  "tibble",
                  "tidyr",
                  "waldronlab/bugphyzzExports")
BiocManager::install(dependencies)
```
### Run export_bugphyzz.R

Run the script, which will produce the files in the directory where the script
is run.

```
Rscript bugphyzzExports\inst\script\export_bugphyzz.R
```
#### For Internal Use

When on supermicro

```
/usr/bin/Rscript --vanilla inst/scripts/export_bugphyzz.R
```
