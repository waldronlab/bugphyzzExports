# bugphyzzExports

This repository contains the devel
[bugphyzz](https://github.com/waldronlab/bugphyzz) signature files and a full
data dump. The devel files are generated weekly.

The release version of these files will soon be available on Zenodo.

## File creation

If desired, you can generate the signatures and full dump file.

### Requirements

You will need to install

* [bugphyzz](https://github.com/waldronlab/bugphyzz)
* [bugsigdbr](https://www.bioconductor.org/packages/bugsigdbr/)
* [readr](https://readr.tidyverse.org/)

Install using `BiocManager`:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("waldronlab/bugphyzz", "bugsigdbr", "readr"))
```
### Run export_bugphyzz.R

Run the script, which will produce the files in the directory where the script
is run.

```
Rscript bugphyzzExports\inst\script\export_bugphyzz.R
```
