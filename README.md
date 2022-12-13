# bugphyzzExports

This repository contains data exports from/for the bugphyzz project.

## Data sources

Data files come from different sources:

1. bugphyzz spreadsheets hosted on google drive (imported through the `bugphyzz::physiologies` function).
2. BacDive file hosted on google drive.
3. Madin et al data hosted on google drive.
4. PATRIC data hosted on google drive.

## Exports workflow 

All of this workflow happens in an Rscript (see `inst/scripts/dump_release.R`.

1. Import physiologies from google spreadsheets with `bugphyzz::physiologies`. This function might become a hidden function.
2. Import data from BacDive, Madin et al, and PATRIC (whcih are in a google drive).
3. Combine datasets.
4. Solve duplicates, conflicts, and agreements (and other possible cases).
5. Propagate annotations (run ASR and inheritance algorithms).
6. Merge data into a single object and export as tsv or csv (could be very large).
7. Use the single object and `bugphyzz::getSignatures` to create signatures.
8. Export .gmt files.

## Zenodo version

The entire repository will be used for Zenodo.

## Intended usage

+ The exported data, character delimited text (csv or tsv) and .gmt files, could be used as is, i.e.,
outside of the R environment.
+ The dump release could be also be imported back into R through the bugphyzz
package with the `bugphyzz::importBugphyzz` function. This function should have
two arguments, version and cache, mirroring `bugsigdbr::imoprtBugSigDB`.

