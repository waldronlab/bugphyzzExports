## Steps
##
## At some point remove row without parent information.
##
## Import with physiologies
## Group by and combine taxon names by taxids (combine synonyms).
## Use a first ASR step to get the NCBI_ID of taxa with missing taxid.
## Filter, keeping only row wit taxid.
## Resolve any conflicts or agreements.
## Run ASR/Inheritance  .
## Create new data.frames.
## Combine data.frames and export the full_dump file.

## Some pacakges installed from github with remotes
## BiocManager::install('waldronlab/bugphyzz', force = TRUE)
## BiocManager::install('sdgamboa/taxPPro', force = TRUE)

library(bugphyzz)
library(taxPPro)
library(purrr)
library(dplyr)

## Import of physiologies
phys <- physiologies(full_source = FALSE, remove_false = TRUE)
x <- phys$aerophilicity
y <- phys$width

x$NCBI_ID[which(is.na(x$NCBI_ID))] <- 'unknown'
x$NCBI_ID[which(is.na(x$Parent_NCBI_ID))] <- 'unknown'

y$NCBI_ID[which(is.na(y$NCBI_ID))] <- 'unknown'
y$NCBI_ID[which(is.na(y$Parent_NCBI_ID))] <- 'unknown'

## Remove taxa without parent NCBI_ID
x <- x[x$NCBI_ID != 'unknown',]
y <- y[y$NCBI_ID != 'unknown',]

## Separate taxa without NCBI ID from taxa with NCBI ID
which(x$NCBI_ID == 'unknown')
which(y$NCBI_ID == 'unknown')





## Code for filtering

## Code for resolving conflicts

## Code for propagation

## Code for exporting
full_dump <- reduce(phys, bind_rows)

##
fname <- paste0("full_dump_bugphyzz_", Sys.Date(), ".csv.bz2")
unlink(fname)
con <- bzfile(fname, "w")
write.csv(full_dump, file=con, quote = TRUE)
close(con)

## Code for creating bugphyzz signatures and exporting to gmt
