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
source('functions.R')

phys_names <- c('aerophilicity', 'growth temperature')
phys <- physiologies(phys_names, full_source = FALSE, remove_false = TRUE)

data_ready <- phys |>
    map(~ tryCatch(error = function(e) e, getDataReadyForPropagation(.x)))
# any(map_int(data_ready, rlang::is_error))

aer <- phys$aerophilicity
gt <- phys$`growth temperature`

x <- getDataReadyForPropagation(gt)
y <- getDuplicates(x)


chooseColVal(y)
## Code for propagation
# data('tree_list')
# tree <- data.tree::as.Node(tree_list)
# propagated <- data_ready |>
    # map(~ tryCatch(error = function(e) e, propagate(tree, .x)))



## Code for propagation
# data('tree_list')
# tree <- data.tree::as.Node(tree_list)
#
# propagated <- data_ready |>
#     map(~ tryCatch(error = function(e) e, propagate(tree, .x)))

## Code for exporting
# full_dump <- reduce(output, bind_rows)
# fname <- paste0("full_dump_bugphyzz_", Sys.Date(), ".csv.bz2")
# unlink(fname)
# con <- bzfile(fname, "w")
# write.csv(full_dump, file=con, quote = TRUE)
# close(con)

## Code for creating bugphyzz signatures and exporting to gmt
