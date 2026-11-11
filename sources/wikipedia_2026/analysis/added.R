## What Wikipedia_2026 adds to the curated data, before conflict resolution and ASR.
## Usage: Rscript added.R <physiologies.rds> <wikipedia_tsv_dir> <records.jsonl> <out_dir>
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 250)
args <- commandArgs(trailingOnly = TRUE)
ph <- readRDS(args[1]); wdir <- args[2]; outdir <- args[4]
rec <- jsonlite::stream_in(file(args[3]), verbose = FALSE)[, c("taxid", "rank", "name")]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rank3 <- function(r) ifelse(r %in% c("genus", "species"), r, "below species")
files <- list.files(wdir, pattern = "\\.tsv$"); files <- setdiff(files, "overrides.tsv")
## every taxon bugphyzz has a curated row for, in any physiology; and the species its strains belong to
ids_any <- unique(unlist(lapply(ph, \(d) d$NCBI_ID)))
ids_parent <- unique(unlist(lapply(ph, \(d) d$Parent_NCBI_ID[d$Rank %in% "strain" & d$Parent_rank %in% "species"])))
sink(file.path(outdir, "added.txt"), split = TRUE)
res <- list(); new_rows <- list(); wk_ids <- character()
for (f in files) {
    phys <- gsub("_", " ", sub("\\.tsv$", "", f))
    w <- read.delim(file.path(wdir, f), quote = "", colClasses = "character", check.names = FALSE)
    w$rank <- rank3(rec$rank[match(w$NCBI_ID, rec$taxid)]); wk_ids <- union(wk_ids, w$NCBI_ID)
    d <- ph[[phys]]
    if (is.null(d)) d <- data.frame(NCBI_ID = character(), Parent_NCBI_ID = character(), Rank = character(), Parent_rank = character(), Attribute = character(), Attribute_value = logical())
    d <- d[is.na(d$Attribute_source) | d$Attribute_source != "Bergey's Manual", ]
    have <- unique(d$NCBI_ID)
    have_strain <- unique(d$Parent_NCBI_ID[d$Rank %in% "strain" & d$Parent_rank %in% "species"])
    w$status <- factor(ifelse(w$NCBI_ID %in% have, "has a curated value", ifelse(w$NCBI_ID %in% have_strain, "only its strains have one", "no curated value")),
                       levels = c("has a curated value", "only its strains have one", "no curated value"))
    if (phys == "habitat") {   # union attribute: a row is new when that habitat term is not recorded for the taxon
        k <- paste(d$NCBI_ID, d$Attribute)[d$Attribute_value %in% TRUE]
        w$new_term <- !paste(w$NCBI_ID, w$Attribute) %in% k
    }
    tx <- distinct(w, NCBI_ID, rank, status)
    res[[phys]] <- tx |> count(physiology = phys, rank, status, .drop = FALSE) |> pivot_wider(names_from = status, values_from = n, values_fill = 0)
    res[[phys]]$rows <- as.integer(table(factor(w$rank, res[[phys]]$rank)))
    if (phys == "habitat") res[[phys]]$new_habitat_terms <- as.integer(tapply(w$new_term, factor(w$rank, res[[phys]]$rank), sum))
    new_rows[[phys]] <- w[w$status != "has a curated value", ] |> mutate(physiology = phys)
}
res <- bind_rows(res)
res <- res |> mutate(taxa = `has a curated value` + `only its strains have one` + `no curated value`) |>
    select(physiology, rank, rows, taxa, `has a curated value`, `only its strains have one`, `no curated value`, any_of("new_habitat_terms"))
cat("Wikipedia_2026 taxa per physiology, by whether bugphyzz already has a curated value of that physiology for the taxon\n")
print(as.data.frame(res), row.names = FALSE)
tot <- res |> group_by(physiology) |> summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
cat("\nAll ranks:\n"); print(as.data.frame(tot), row.names = FALSE)
write.table(res, file.path(outdir, "added_by_physiology.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
wk <- data.frame(NCBI_ID = wk_ids); wk$rank <- rank3(rec$rank[match(wk$NCBI_ID, rec$taxid)])
wk$in_bugphyzz <- ifelse(wk$NCBI_ID %in% ids_any, "taxon has curated rows", ifelse(wk$NCBI_ID %in% ids_parent, "only its strains do", "taxon is new to bugphyzz"))
cat("\nWikipedia_2026 taxa with at least one row, against all curated physiologies:\n")
print(as.data.frame(count(wk, rank, in_bugphyzz) |> pivot_wider(names_from = in_bugphyzz, values_from = n, values_fill = 0)), row.names = FALSE)
nr <- bind_rows(new_rows) |> select(physiology, NCBI_ID, Taxon_name, rank, Attribute, Attribute_value, Frequency, status, Note)
write.table(nr, file.path(outdir, "new_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n%d rows for taxa without a curated value of that physiology written to new_annotations.tsv\n", nrow(nr)))
sink()
