#!/usr/bin/env Rscript
## Fixture test for import.R: run with `Rscript test_import.R` from this directory.
suppressPackageStartupMessages(library(dplyr))
fx <- tempfile(fileext = ".csv")
write.csv(data.frame(
    tax_id = c(1, 1, 1, 2, 2, 3, 3), species_tax_id = c(1, 1, 1, 2, 2, 3, 3), species = c("A", "A", "A", "B", "B", "C", "C"),
    genus = "G", data_source = c("s1", "s2", "s3", "s1", "s2", "s1", "s2"),
    gram_stain = c("positive", "positive", "negative", "positive", NA, NA, NA),
    metabolism = c("aerobic", "aerobic", "aerobic", "facultative", "aerobic", NA, NA),
    sporulation = c("yes", "no", NA, NA, NA, NA, NA), motility = c("yes", "yes", "no", NA, NA, NA, NA),
    cell_shape = c("vibrio", "vibrio", NA, "bacillus", "bacillus", NA, NA),
    genome_size = c(100, 200, NA, NA, NA, 5, NA), gc_content = c(50, 52, NA, NA, NA, 90, NA),
    coding_genes = NA, optimum_tmp = NA, optimum_ph = NA, org_name = NA, family = NA, order = NA, class = NA, phylum = NA, superkingdom = NA), fx, row.names = FALSE)
res <- tempfile(fileext = ".rds")
Sys.setenv(MADIN_CSV = fx, MADIN_TEST = res)
source("import.R")
out <- readRDS(res)
g <- out[["gram stain"]]; stopifnot(nrow(g) == 3)                      # A: positive 2/3 + negative 1/3 (>= 0.2) ; B: positive 1/1
stopifnot(g$Frequency[g$NCBI_ID == "1" & g$Attribute == "gram stain positive"] == "usually")
stopifnot(g$Frequency[g$NCBI_ID == "1" & g$Attribute == "gram stain negative"] == "sometimes")
stopifnot(g$Frequency[g$NCBI_ID == "2"] == "always")
a <- out$aerophilicity; stopifnot(setequal(a$Attribute[a$NCBI_ID == "2"], c("facultatively anaerobic", "aerobic")))  # 50/50 kept as two rows
stopifnot(nrow(out[["spore formation"]][out[["spore formation"]]$NCBI_ID == "1", ]) == 0)  # binary tie dropped
stopifnot(out$motility$Attribute_value[out$motility$NCBI_ID == "1"] == "TRUE")             # binary majority
stopifnot(out$shape$Attribute[out$shape$NCBI_ID == "1"] == "vibrio")
stopifnot(out[["genome size"]]$Attribute_value[out[["genome size"]]$NCBI_ID == "1"] == 150) # median
stopifnot(!"3" %in% out[["gc content"]]$NCBI_ID, out[["gc content"]]$Attribute_value[out[["gc content"]]$NCBI_ID == "1"] == "50-52")
cat("import.R fixture test: OK\n")
