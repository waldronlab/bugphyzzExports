## Steps
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





