Weekly export status: ![weekly export](https://github.com/waldronlab/bugphyzzExports/actions/workflows/export-bugphyzz.yml/badge.svg)

# bugphyzzExports

bugphyzz is a database that harmonizes physiological and other microbial
trait annotations from different sources using a controlled vocabulary and
ontology terms. Furthermore, these annotations are propagated to
uncharacterized microbes through Ancestral State Reconstruction (ASR).

You can learn more about this project [here](https://github.com/waldronlab/bugphyzz).

This repository contains the code for resolving conflicting annotations
and run the ASR step. It also contains the devel version of the annotations
(before being released on Zenodo) distributed across different text files.
The *.csv files contain the data in tabular format and are imported through
the `bugphyzz::importBugphyzz` function in [R](https://github.com/waldronlab/bugphyzz).
The *gmt files contain lists of microbial signatures in GMT format
created with the `bugphyzz::makeSignatures` function.
Annotations are at genus, species and strain level.

The data schema is described [here](https://github.com/waldronlab/bugphyzz)

The devel files are generated weekly.

## File creation

If desired, anyone can generate the *.csv and *.gmt files.

### 1. Install required R packages

The first step is downloading the repo:

```bash
git clone https://github.com/waldronlab/bugphyzzExports.git
cd bugphyzzExports
```
The following packages need to be installed in the R environment:

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

This could be accomplished for example with:

```r
## Inside an R session
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
Or running `devtools::install_deps(dependencies=TRUE)` in an R session within
the main directory.

### 2. Run the inst/scripts/export_bugphyzz.R script

Run the script, which will produce the files in the directory where the script
is run. Preferably run inside the main directory of the project.

On a linux-like terminal:

```
Rscript inst/scripts/export_bugphyzz.R
```

On supermicro (for internal use):

```
/usr/bin/Rscript --vanilla inst/scripts/export_bugphyzz.R
```

## License

The files are available under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode).

## Zenodo

Find this dataset on Zenodo (latest realease version): https://zenodo.org/doi/10.5281/zenodo.10980653 

A release deposit holds the three `bugphyzz_*.csv` files and the GMT files. `bugphyzz_asr_dropped.csv.gz` stays on `devel`: it is the imputed annotations the validation thresholds removed, kept for threshold experiments (`bugphyzz::importBugphyzz(includeDroppedAsr = TRUE)`), not release data.

## Versioning (recommended)

Some recommendations about versioning for relase.

Format: x.y.z  
Example: 1.0.2

The third digit (z) should be used to fix typos or any other minor
trouble with the annotations. Essentially these are the same annotations,
but with minor adjustments.

The second digit (y) should be used for major adjustments such as fixing the
way conflicting annotations are handled or
adjusting ASR methods/parameters, say choosing a different phylogenetic tree
or using a different package for running ASR.

The first digit (x) should be reserved for major changes, such as adding
new datasets or using a completely different approach for propagating
annotations, etc.

## Validation of the ASR method

A 10-fold cross-validation (https://github.com/waldronlab/taxPProValidation/)
scores the ASR method per attribute. The export keeps an attribute's imputed
annotations all together or not at all ([ADR 0006](docs/adr/0006-impute-whole-attributes-by-distance-matched-kappa.md)):
they are kept when the lower 95% bound of the attribute's Cohen's kappa
(weighted for ordered states; R2 for numeric attributes), with held-out taxa
weighted to match the phylogenetic distance of the taxa being imputed, is at
least 0.7 (0.5 for R2), and only rows with `Score` >= 0.5 are kept. The
scores come from `validation_attributes.tsv` in the installed bugphyzz
package, which `importBugphyzz()` also attaches as `Validation`. The imputed
annotations the export leaves out are in `bugphyzz_asr_dropped.csv.gz`
(`devel` only), and `importBugphyzz()` reads them when any of its thresholds
is set below the export's.

