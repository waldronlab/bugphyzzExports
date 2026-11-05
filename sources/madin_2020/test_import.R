#!/usr/bin/env Rscript
## Fixture test for import.R: run with `Rscript test_import.R` from this directory.
suppressPackageStartupMessages(library(dplyr))
fx <- tempfile(fileext = ".csv")
write.csv(data.frame(
    tax_id = c(1, 1, 1, 2, 2, 3, 3, 1, 2), species_tax_id = c(1, 1, 1, 2, 2, 3, 3, 1, 2), species = c("A", "A", "A", "B", "B", "C", "C", "A", "B"),
    genus = "G", data_source = c("s1", "s2", "s3", "s1", "s2", "s1", "s2", "protraits", "microbe-directory"),
    gram_stain = c("positive", "positive", "negative", "positive", NA, NA, NA, "negative", "negative"),
    metabolism = c("aerobic", "aerobic", "aerobic", "facultative", "aerobic", NA, NA, "anaerobic", "anaerobic"),
    sporulation = c("yes", "no", NA, NA, NA, NA, NA, "yes", NA), motility = c("yes", "yes", "no", NA, NA, NA, NA, "no", NA),
    cell_shape = c("vibrio", "vibrio", NA, "bacillus", "bacillus", NA, NA, "coccus", NA),
    genome_size = c(100, 200, NA, NA, NA, 5, NA, 9000, NA), gc_content = c(50, 52, NA, NA, NA, 90, NA, 10, NA),
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
## direct-source records excluded: the protraits and microbe-directory rows carry conflicting values
## (A gram negative, anaerobic, coccus, genome 9000; B gram negative) that must not appear
ds <- unlist(lapply(out, function(x) x$Data_source)); stopifnot(!any(grepl("protraits|microbe-directory", ds)))
stopifnot(!"gram stain negative" %in% out[["gram stain"]]$Attribute[out[["gram stain"]]$NCBI_ID == "2"])
stopifnot(!"anaerobic" %in% out$aerophilicity$Attribute[out$aerophilicity$NCBI_ID == "1"])
stopifnot(!"coccus" %in% out$shape$Attribute[out$shape$NCBI_ID == "1"])
stopifnot(out[["gram stain"]]$Data_source[out[["gram stain"]]$NCBI_ID == "1" & out[["gram stain"]]$Attribute == "gram stain positive"] == "s1;s2")
cat("import.R fixture test: OK\n")
