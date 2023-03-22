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

phys <- physiologies(full_source = FALSE, remove_false = TRUE)
output <- vector('list', length(phys))
for (i in seq_along(output)) {
    names(output)[i] <- names(phys)[i]
    output[[i]] <- tryCatch(
        error = function(e) e,
        {
            getDataReady(phys[[i]])
        }
    )
}
output <- discard(output, ~ rlang::is_error(.x))

## Code for propagation
# data('tree_list')
# tree <- data.tree::as.Node(tree_list)



# x_propagated <- propagate(data_tree = tree, df = resolvedConflicts)
# y_propagated <- propagate(data_tree = tree, df = resolvedConflicts)

## Code for exporting
# full_dump <- reduce(output, bind_rows)
# fname <- paste0("full_dump_bugphyzz_", Sys.Date(), ".csv.bz2")
# unlink(fname)
# con <- bzfile(fname, "w")
# write.csv(full_dump, file=con, quote = TRUE)
# close(con)

## Code for creating bugphyzz signatures and exporting to gmt
