## Change in the export (curated values, inheritance and ASR) when Wikipedia_2026 is added.
## Usage: Rscript compare_exports.R <baseline_export_dir> <with_wikipedia_export_dir> <out_dir>
## Each directory holds bugphyzz_binary.csv, bugphyzz_multistate.csv, bugphyzz_numeric.csv
## and bugphyzz_asr_dropped.csv.gz written by inst/scripts/export_bugphyzz.R.
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 250)
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[3]; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rd <- function(dir) {
    x <- lapply(c("binary", "multistate", "numeric"), \(t)
        read.csv(file.path(dir, paste0("bugphyzz_", t, ".csv")), comment.char = "#", colClasses = "character"))
    x <- bind_rows(x)
    x$Score <- as.numeric(x$Score)
    x$kind <- ifelse(x$Evidence %in% c("asr", "inh", "tax"), x$Evidence, "curated")
    x$group <- ifelse(x$Attribute_type == "multistate-union", "habitat and other union attributes",
               ifelse(x$Attribute_type == "multistate-intersection", sub(":.*", "", x$Attribute), x$Attribute))
    x
}
b <- rd(args[1]); w <- rd(args[2])
## the multistate export has the state in Attribute_value and the physiology in Attribute
sink(file.path(outdir, "compare_exports.txt"), split = TRUE)
cat("Rows in the export by how the value was obtained (curated source, inh = inherited, tax = taxonomic pooling, asr = ancestral state reconstruction)\n")
tab <- function(x, lab) count(x, Attribute, kind) |> mutate(run = lab)
t1 <- bind_rows(tab(b, "baseline"), tab(w, "wikipedia")) |> pivot_wider(names_from = c(kind, run), values_from = n, values_fill = 0, names_sort = TRUE)
print(as.data.frame(t1), row.names = FALSE)
cat("\nTotals:\n"); print(as.data.frame(bind_rows(tab(b, "baseline"), tab(w, "wikipedia")) |> count(run, kind, wt = n) |> pivot_wider(names_from = run, values_from = n)), row.names = FALSE)

cat("\nRows carrying Attribute_source Wikipedia_2026 in the export (the ones that survived conflict resolution):\n")
print(as.data.frame(w |> filter(Attribute_source %in% "Wikipedia_2026") |> count(Attribute, Rank) |> pivot_wider(names_from = Rank, values_from = n, values_fill = 0)), row.names = FALSE)

## taxon by attribute coverage
cov <- function(x) distinct(x, Attribute, NCBI_ID, Rank)
cb <- cov(b); cw <- cov(w)
gained <- anti_join(cw, cb, by = c("Attribute", "NCBI_ID")); lost <- anti_join(cb, cw, by = c("Attribute", "NCBI_ID"))
cat("\nTaxa annotated per attribute: baseline, with Wikipedia, gained, lost\n")
t2 <- full_join(count(cb, Attribute, name = "baseline"), count(cw, Attribute, name = "wikipedia"), by = "Attribute") |>
    left_join(count(gained, Attribute, name = "gained"), by = "Attribute") |> left_join(count(lost, Attribute, name = "lost"), by = "Attribute") |>
    mutate(across(where(is.numeric), \(v) coalesce(v, 0L)))
print(as.data.frame(t2), row.names = FALSE)
cat(sprintf("\nDistinct taxa in the export: baseline %d, with Wikipedia %d\n", n_distinct(b$NCBI_ID), n_distinct(w$NCBI_ID)))
g2 <- left_join(gained, distinct(w, Attribute, NCBI_ID, kind), by = c("Attribute", "NCBI_ID"))
cat("\nGained taxon by attribute pairs, by rank and by how the new value was obtained:\n")
print(as.data.frame(count(g2, Rank, kind) |> pivot_wider(names_from = kind, values_from = n, values_fill = 0)), row.names = FALSE)

## values that changed for taxon by attribute pairs present in both runs
top <- function(x) x |> filter(Attribute_type != "numeric", Attribute_type != "multistate-union") |>
    group_by(Attribute, NCBI_ID) |> slice_max(Score, n = 1, with_ties = FALSE) |> ungroup() |>
    select(Attribute, NCBI_ID, Taxon_name, Rank, value = Attribute_value, kind, Score, Attribute_source)
j <- inner_join(top(b), top(w), by = c("Attribute", "NCBI_ID"), suffix = c(".base", ".wiki"))
ch <- j[j$value.base != j$value.wiki, ]
cat("\nDiscrete attributes: taxa whose top-scoring state changed, by how the value was obtained before and after\n")
print(as.data.frame(count(ch, Attribute, kind.base, kind.wiki) |> arrange(Attribute, desc(n))), row.names = FALSE)
cat(sprintf("Changed: %d of %d taxon by attribute pairs present in both runs (%.2f%%)\n", nrow(ch), nrow(j), 100 * nrow(ch) / nrow(j)))
asr <- j[j$kind.base == "asr" & j$kind.wiki == "asr", ]
cat("\nASR-only pairs in both runs: same top state, and change in its score\n")
print(as.data.frame(asr |> group_by(Attribute) |> summarise(pairs = n(), same_state = round(100 * mean(value.base == value.wiki), 2),
    median_abs_score_change = round(median(abs(Score.wiki - Score.base)), 4), .groups = "drop")), row.names = FALSE)
write.table(ch, file.path(outdir, "changed_states.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

num <- function(x) x |> filter(Attribute_type == "numeric") |> mutate(v = as.numeric(Attribute_value)) |>
    group_by(Attribute, NCBI_ID) |> summarise(v = mean(v), kind = kind[1], .groups = "drop")
jn <- inner_join(num(b), num(w), by = c("Attribute", "NCBI_ID"), suffix = c(".base", ".wiki"))
cat("\nNumeric attributes, taxa in both runs: change in the exported value\n")
print(as.data.frame(jn |> group_by(Attribute, kind.wiki) |> summarise(taxa = n(), changed = sum(abs(v.wiki - v.base) > 1e-9),
    median_abs_change_when_changed = round(median(abs(v.wiki - v.base)[abs(v.wiki - v.base) > 1e-9]), 3), .groups = "drop")), row.names = FALSE)

## attributes that pass or fail validation
dr <- function(dir) { x <- read.csv(gzfile(file.path(dir, "bugphyzz_asr_dropped.csv.gz")), comment.char = "#", colClasses = "character"); count(x, Attribute, name = "n") }
d <- full_join(dr(args[1]), dr(args[2]), by = "Attribute", suffix = c(".base", ".wiki"))
cat("\nASR annotations dropped (attribute below the validation threshold, or Score < 0.5):\n"); print(as.data.frame(d), row.names = FALSE)
sink()
