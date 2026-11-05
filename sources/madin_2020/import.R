#!/usr/bin/env Rscript
## Import Madin et al. 2020 (Sci Data 7:170) condensed species traits into
## bugphyzz curation-sheet rows, one TSV per physiology.
##
## Source: https://github.com/bacteria-archaea-traits/bacteria-archaea-traits,
## output/condensed_traits_NCBI.csv at the commit below. One row per
## (taxon, data source); we take the majority value per NCBI species taxid.
## Everything with a bugphyzz attribute is imported, including aerophilicity
## (`metabolism`), which duplicates the Madin_2020 rows pasted into the
## aerophilicity Google Sheet in 2022; those sheet rows should be removed
## once this file is in use. Keep all data; decide at conflict resolution.
##
## Usage: Rscript import.R [out_dir]
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr)})
commit <- "396dcd4d98ca20b613ba280f5514611651b6fb91" # 2020-10-16; update when refreshing
url <- sprintf("https://raw.githubusercontent.com/bacteria-archaea-traits/bacteria-archaea-traits/%s/output/condensed_traits_NCBI.csv", commit)
out_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
local <- Sys.getenv("MADIN_CSV", unset = "")
m <- read_csv(if (nzchar(local)) local else url, show_col_types = FALSE, guess_max = 2e5)
SOURCE <- "Madin_2020"
cols <- c("NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency", "Note")
freq_of <- function(frac) ifelse(frac >= 0.9, "always", ifelse(frac >= 0.5, "usually", "sometimes"))
species <- m |> filter(!is.na(species_tax_id)) |> mutate(NCBI_ID = as.character(species_tax_id), Taxon_name = species)

## Discrete traits: majority value per species, Frequency from its share of records.
discrete <- function(col, phys, map) {
    d <- species |> filter(!is.na(.data[[col]])) |> mutate(val = unname(map[.data[[col]]])) |> filter(!is.na(val)) |>
        distinct(NCBI_ID, Taxon_name, data_source, val) |> count(NCBI_ID, Taxon_name, val) |>
        group_by(NCBI_ID) |> mutate(frac = n / sum(n)) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup()
    tibble(NCBI_ID = d$NCBI_ID, Genome_ID = "", Accession_ID = "", Taxon_name = d$Taxon_name,
           Attribute = if (phys %in% c("motility", "spore formation")) phys else d$val,
           Attribute_value = if (phys %in% c("motility", "spore formation")) d$val else "TRUE",
           Attribute_source = SOURCE, Frequency = freq_of(d$frac), Note = sprintf("%d of %d records", d$n, round(d$n / d$frac)))
}
out <- list(
    "aerophilicity" = discrete("metabolism", "aerophilicity", c(aerobic = "aerobic", `obligate aerobic` = "aerobic", anaerobic = "anaerobic", `obligate anaerobic` = "anaerobic", `strictly anaerobic` = "anaerobic", facultative = "facultatively anaerobic", microaerophilic = "microaerophilic")),
    "gram stain" = discrete("gram_stain", "gram stain", c(positive = "gram stain positive", negative = "gram stain negative")),
    "spore formation" = discrete("sporulation", "spore formation", c(yes = "TRUE", no = "FALSE")),
    "motility" = discrete("motility", "motility", c(yes = "TRUE", flagella = "TRUE", gliding = "TRUE", `axial filament` = "TRUE", no = "FALSE")),
    "shape" = discrete("cell_shape", "shape", c(bacillus = "bacillus", coccus = "coccus", coccobacillus = "coccobacillus", spiral = "spirillum", spirochete = "spirillum", filament = "filamentous", pleomorphic = "pleiomorphic cell", vibrio = "bacillus"))
)
## Numeric traits: median per species.
numeric <- function(col, phys, attr, unit) {
    d <- species |> filter(!is.na(.data[[col]])) |> group_by(NCBI_ID, Taxon_name) |> summarise(v = median(.data[[col]]), n = n(), .groups = "drop")
    tibble(NCBI_ID = d$NCBI_ID, Genome_ID = "", Accession_ID = "", Taxon_name = d$Taxon_name, Attribute = attr, Attribute_value = d$v,
           Attribute_source = SOURCE, Frequency = "always", Unit = unit, Note = sprintf("median of %d records", d$n))
}
out[["genome size"]] <- numeric("genome_size", "genome size", "Genome size", "base pair")
out[["coding genes"]] <- numeric("coding_genes", "coding genes", "genes", "genes")
out[["growth temperature"]] <- numeric("optimum_tmp", "growth temperature", "growth temperature", "Celcius")
out[["optimal ph"]] <- numeric("optimum_ph", "optimal ph", "optimal ph", "pH")
for (phys in names(out)) {
    f <- file.path(out_dir, paste0(gsub(" ", "_", phys), ".tsv"))
    write_tsv(out[[phys]], f, na = "")
    cat(sprintf("%s: %d rows, %d species\n", phys, nrow(out[[phys]]), n_distinct(out[[phys]]$NCBI_ID)))
}
