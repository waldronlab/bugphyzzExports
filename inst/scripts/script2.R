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

## Import of physiologies
# phys <- physiologies(full_source = FALSE, remove_false = TRUE)
phys_names <- c('aerophilicity', 'growth temperature')
phys <- physiologies(phys_names, full_source = FALSE, remove_false = TRUE)
x <- phys[['aerophilicity']]
y <- phys[['growth temperature']]

## replace NAs with unknowns in NCBI_ID and Parent_NCBI_ID
x$NCBI_ID[which(is.na(x$NCBI_ID))] <- 'unknown'
x$Parent_NCBI_ID[which(is.na(x$Parent_NCBI_ID))] <- 'unknown'

y$NCBI_ID[which(is.na(y$NCBI_ID))] <- 'unknown'
y$Parent_NCBI_ID[which(is.na(y$Parent_NCBI_ID))] <- 'unknown'

## Remove taxa without parent NCBI_ID
x <- x[x$Parent_NCBI_ID != 'unknown',]
y <- y[y$Parent_NCBI_ID != 'unknown',]

## Remove taxa for which we do not have enough information
x <- x[!is.na(x$Rank), ]
x <- x[!is.na(x$Evidence), ]
x <- x[!is.na(x$Frequency), ]
x <- x[!is.na(x$Confidence_in_curation), ]
x <- unique(x)

y <- y[!is.na(y$Rank), ]
y <- y[!is.na(y$Evidence), ]
y <- y[!is.na(y$Frequency), ]
y <- y[!is.na(y$Confidence_in_curation), ]
y <- unique(y)

## Remove indirect parents
# x <- filterDirectParents(x)
# y <- filterDirectParents(y)

##  Add frequency scores
x <- freq2Scores(x)
y <- freq2Scores(y)

## Separate taxa without NCBI ID from taxa with NCBI ID
x_yesid <- x[which(x$NCBI_ID != 'unknown'),]
x_noid <- x[which(x$NCBI_ID == 'unknown'),]

y_yesid <- y[which(y$NCBI_ID != 'unknown'),]
y_noid <- y[which(y$NCBI_ID == 'unknown'),]

## Let's get a small ASR code
x_noid_asr <- calcParentScores(x_noid)
y_noid_asr <- calcParentScores(y_noid)

## let's combine and get a new dataset
x_new <- dplyr::bind_rows(x_yesid, x_noid_asr)
cols <- c(
    'NCBI_ID', 'Rank',
    'Attribute', 'Attribute_value', 'Attribute_source',
    'Evidence', 'Frequency',
    'Attribute_type', 'Attribute_group',
    'Confidence_in_curation', 'Score'
)
x_new <- unique(x_new[,cols])

y_new <- dplyr::bind_rows(y_yesid, y_noid_asr)
cols2 <- c(
    'NCBI_ID', 'Rank',
    'Attribute', 'Attribute_value_min', 'Attribute_value_max',
    'Attribute_source',
    'Evidence', 'Frequency',
    'Attribute_type', 'Attribute_group',
    'Confidence_in_curation', 'Score'
)
y_new <- unique(y_new[,cols2])

x_ready <- prepareData2(x_new)
y_ready <- prepareData2(y_new)

## Code for resolving conflicts (this only applies for numeric values)
## This only applies to character/numerical values
resolvedAgreements <- resolveAgreements(x_ready)
resolvedConflicts <- resolveConflicts(resolvedAgreements)

## Code for propagation
data('tree_list')
tree <- data.tree::as.Node(tree_list)
x_propagated <- propagate(data_tree = tree, df = x_ready)
y_propagated <- propagate(data_tree = tree, df = y_ready)

## Code for exporting
# full_dump <- reduce(phys, bind_rows)

x_ready

##
# fname <- paste0("full_dump_bugphyzz_", Sys.Date(), ".csv.bz2")
# unlink(fname)
# con <- bzfile(fname, "w")
# write.csv(full_dump, file=con, quote = TRUE)
# close(con)

## Code for creating bugphyzz signatures and exporting to gmt
