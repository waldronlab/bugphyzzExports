Weekly export status: ![weekly export](https://github.com/waldronlab/bugphyzzExports/actions/workflows/export-bugphyzz.yml/badge.svg)

# bugphyzzExports

This repository contains the devel
[bugphyzz](https://github.com/waldronlab/bugphyzz) signature files and a full
data dump. The devel files are generated weekly.

## File creation

If desired, you can generate the signatures and full dump file.

### Requirements

* [bugphyzz](https://github.com/waldronlab/bugphyzz)
* [bugsigdbr](https://bioconductor.org/packages/release/bioc/html/bugsigdbr.html)
* [castor](https://cran.r-project.org/web/packages/castor/)
* [dplyr](https://cran.r-project.org/web/packages/dplyr)
* [logr](https://cran.r-project.org/web/packages/logr/)
* [phytools](https://cran.r-project.org/web/packages/phytools/)
* [purrr](https://cran.r-project.org/web/packages/purrr)
* [rlang](https://cran.r-project.org/web/packages/rlang)
* [sessioninfo](https://cran.r-project.org/web/packages/sessioninfo)
* [stringr](https://cran.r-project.org/web/packages/stringr)
* [taxPPro](https://github.com/waldronlab/taxPPro)
* [tibble](https://cran.r-project.org/web/packages/tibble/)
* [tidyr](https://cran.r-project.org/web/packages/tidyr/)

Install using `BiocManager`:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
dependencies <- c(
    "waldronlab/bugphyzz",
    "bugsigdbr",
    "castor",
    "dplyr",
    "logr",
    "phytools"
    "purrr",
    "rlang",
    "sessioninfo",
    "stringr",
    "waldronlab/taxPPro",
    "tibble",
    "tidyr"
)
BiocManager::install(dependencies)
```
### Run export_bugphyzz.R

Run the script, which will produce the files in the directory where the script
is run. Preferably run inside the project main directory.

On a linux-like terminal:

```
Rscript inst/script/export_bugphyzz.R
```
#### For Internal Use

When on supermicro:

```
/usr/bin/Rscript --vanilla inst/scripts/export_bugphyzz.R
```

### LICENSE

The files are available under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode).

