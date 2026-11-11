## Agreement of Wikipedia_2026 with the other curated sources, and the same
## measure between every other pair of sources, so that Wikipedia can be placed
## among sources whose confidence is already set.
## Usage: Rscript compare_sources.R <physiologies.rds> <wikipedia_tsv_dir> <records.jsonl> <out_dir>
## <physiologies.rds> is saveRDS(bugphyzz::physiologies("all")) without the
## Wikipedia source; <records.jsonl> is collect.py's output (for the ranks).
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 250)
args <- commandArgs(trailingOnly = TRUE)
ph <- readRDS(args[1]); wdir <- args[2]; outdir <- args[4]
rec <- jsonlite::stream_in(file(args[3]), verbose = FALSE)[, c("taxid", "rank")]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
WIKI <- "Wikipedia_2026"
multistate <- c("aerophilicity", "gram stain", "shape", "arrangement", "spore shape")
binary <- c("motility", "spore formation", "animal pathogen", "plant pathogenicity", "host-associated")
numeric <- c("growth temperature", "optimal ph", "gc content", "length", "width")
conf <- ph |> lapply(\(d) distinct(d, Attribute_source, Confidence_in_curation, Evidence)) |> bind_rows() |>
    filter(!is.na(Attribute_source)) |> distinct(Attribute_source, .keep_all = TRUE)

read_wiki <- function(phys) {
    f <- file.path(wdir, paste0(gsub(" ", "_", phys), ".tsv"))
    if (!file.exists(f)) return(NULL)
    d <- read.delim(f, quote = "", colClasses = "character", check.names = FALSE)
    d$rank <- rec$rank[match(d$NCBI_ID, rec$taxid)]
    d
}
## One key per row: a strain under a species counts for the species, so that
## BacDive (strain records) can be compared with sources that annotate species.
key_of <- function(d) ifelse(d$Rank %in% "strain" & d$Parent_rank %in% "species", d$Parent_NCBI_ID, d$NCBI_ID)
norm_state <- function(x) recode(x, "obligately aerobic" = "aerobic", "obligately anaerobic" = "anaerobic", "aerotolerant" = "anaerobic")

## long table of discrete states: source, key, state, n records
states <- function(phys) {
    d <- ph[[phys]]
    d <- d[!is.na(d$Attribute_source) & d$Attribute_source != "Bergey's Manual", ]  # superseded by Bergeys_2026
    d$key <- key_of(d)
    if (phys %in% binary) {
        d <- d[!is.na(d$Attribute_value), ]; d$state <- as.character(d$Attribute_value)
    } else {
        d <- d[d$Attribute_value %in% TRUE, ]; d$state <- norm_state(d$Attribute)
    }
    ref <- count(d, source = Attribute_source, key, state)
    w <- read_wiki(phys)
    if (!is.null(w)) {
        w$state <- if (phys %in% binary) w$Attribute_value else norm_state(w$Attribute)
        if (!phys %in% binary) w <- w[w$Attribute_value == "TRUE", ]
        ref <- bind_rows(ref, count(w, source = Attribute_source, key = NCBI_ID, state))
    }
    ref |> group_by(source, key) |> mutate(frac = n / sum(n)) |> ungroup()
}
## per source and key: the states held by at least 20% of its records, and the majority state
summarise_states <- function(s) {
    s |> group_by(source, key) |>
        summarise(set = list(sort(state[frac >= 0.2])), major = state[which.max(n)], nrec = sum(n), .groups = "drop")
}
kappa2 <- function(a, b) {
    if (length(a) < 2) return(NA_real_)
    lv <- union(a, b); t <- table(factor(a, lv), factor(b, lv)); n <- sum(t)
    po <- sum(diag(t)) / n; pe <- sum(rowSums(t) * colSums(t)) / n^2
    if (pe == 1) NA_real_ else (po - pe) / (1 - pe)
}
agree <- function(ss, a, b, keys = NULL) {
    x <- ss[ss$source == a, ]; y <- ss[ss$source == b, ]
    if (!is.null(keys)) x <- x[x$key %in% keys, ]
    ## only states both sources can express (BacDive here has no "microaerophilic" or "coccobacillus")
    vocab <- intersect(unlist(ss$set[ss$source == a]), unlist(ss$set[ss$source == b]))
    x$set <- lapply(x$set, intersect, vocab); y$set <- lapply(y$set, intersect, vocab)
    x <- x[lengths(x$set) > 0, ]; y <- y[lengths(y$set) > 0, ]
    j <- inner_join(x, y, by = "key", suffix = c(".a", ".b"))
    if (!nrow(j)) return(NULL)
    share <- mapply(\(p, q) length(intersect(p, q)) > 0, j$set.a, j$set.b)
    single <- lengths(j$set.a) == 1 & lengths(j$set.b) == 1
    sa <- vapply(j$set.a[single], `[`, "", 1); sb <- vapply(j$set.b[single], `[`, "", 1)
    tibble(source = a, reference = b, taxa = nrow(j), share_a_state = round(100 * mean(share), 1),
           single_state_taxa = sum(single), same_state = round(100 * mean(sa == sb), 1), kappa = round(kappa2(sa, sb), 3))
}

sink(file.path(outdir, "compare_sources.txt"), split = TRUE)
cat("Agreement between sources on the same taxon. Strain records count for their species.\n",
    "share_a_state: % of taxa where the two sources have a state in common (states held by >= 20% of a source's records for the taxon).\n",
    "same_state, kappa: among taxa where each source gives exactly one state.\n",
    "Only states that both sources of a pair use are compared.\n", sep = "")
all_pairs <- list(); disagreements <- list()
for (phys in c(multistate, binary)) {
    s <- states(phys); ss <- summarise_states(s)
    srcs <- setdiff(unique(ss$source), WIKI)
    if (!WIKI %in% ss$source) { cat("\n===", phys, ": no Wikipedia rows\n"); next }
    wk <- ss[ss$source == WIKI, ]
    wk$rank <- rec$rank[match(wk$key, rec$taxid)]
    cat(sprintf("\n=== %s: Wikipedia taxa %d (genus %d, species %d, below species %d)\n", phys, nrow(wk),
                sum(wk$rank == "genus"), sum(wk$rank == "species"), sum(!wk$rank %in% c("genus", "species"))))
    res <- list()
    for (r in c("species", "genus")) {
        keys <- wk$key[wk$rank == r]
        for (b in srcs) { a1 <- agree(ss, WIKI, b, keys); if (!is.null(a1)) res[[length(res) + 1]] <- mutate(a1, rank = r) }
    }
    res <- bind_rows(res)
    if (nrow(res)) res <- left_join(res, conf, by = c("reference" = "Attribute_source")) |>
        rename(reference_confidence = Confidence_in_curation, reference_evidence = Evidence)
    cat("Wikipedia against each source:\n"); print(as.data.frame(filter(res, taxa >= 10)), row.names = FALSE)
    ## the same measure between the other sources, on the keys Wikipedia could be compared on
    ## (all keys, not only Wikipedia's, so the numbers are stable)
    oth <- list()
    for (a in srcs) for (b in srcs) if (a < b) { a1 <- agree(ss, a, b); if (!is.null(a1)) oth[[length(oth) + 1]] <- a1 }
    oth <- bind_rows(oth)
    if (nrow(oth)) {
        cat("Other pairs of sources, all taxa they share:\n"); print(as.data.frame(filter(oth, taxa >= 30)), row.names = FALSE)
    }
    all_pairs[[phys]] <- bind_rows(mutate(res, physiology = phys), mutate(oth, physiology = phys, rank = "any"))
    ## taxa where Wikipedia gives one state and a high-confidence source gives another, unanimously
    for (b in intersect(c("BacDive", "Bergeys_2026"), srcs)) {
        rf <- ss[ss$source == b, ]; vocab <- unlist(rf$set)
        j <- inner_join(wk, rf, by = "key", suffix = c(".w", ".r"))
        j <- j[vapply(j$set.w, \(p) any(p %in% vocab), TRUE), ]
        bad <- j[!mapply(\(p, q) length(intersect(p, q)) > 0, j$set.w, j$set.r), ]
        if (nrow(bad)) disagreements[[length(disagreements) + 1]] <- tibble(
            physiology = phys, NCBI_ID = bad$key, rank = bad$rank, wikipedia = vapply(bad$set.w, paste, "", collapse = "|"),
            reference = b, reference_states = vapply(bad$set.r, paste, "", collapse = "|"), reference_records = bad$nrec.r)
    }
}
pairs <- bind_rows(all_pairs)
write.table(pairs, file.path(outdir, "agreement_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dis <- bind_rows(disagreements)
if (nrow(dis)) {
    nm <- lapply(c(multistate, binary), read_wiki) |> bind_rows() |> distinct(NCBI_ID, Taxon_name, Note)
    dis <- left_join(dis, distinct(nm, NCBI_ID, .keep_all = TRUE), by = "NCBI_ID")
    write.table(dis, file.path(outdir, "disagreements_with_high_confidence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\n=== Wikipedia states that share nothing with a high-confidence source (full list in disagreements_with_high_confidence.tsv)\n")
    print(as.data.frame(count(dis, physiology, reference, wikipedia, reference_states, sort = TRUE) |> head(40)), row.names = FALSE)
}

## Summary: each source against BacDive (strain-level laboratory records, confidence high), species only
cat("\n=== Each source against BacDive, taxa where both give one state (same_state %, kappa, n)\n")
summ <- pairs |> filter((reference == "BacDive" & (rank %in% c("species", "any"))) | (source == "BacDive" & rank == "any")) |>
    mutate(other = ifelse(source == "BacDive", reference, source)) |>
    left_join(conf, by = c("other" = "Attribute_source")) |>
    transmute(physiology, source = other, confidence = ifelse(other == WIKI, "?", Confidence_in_curation),
              taxa, share_a_state, single_state_taxa, same_state, kappa) |> arrange(physiology, desc(same_state))
print(as.data.frame(summ), row.names = FALSE)
write.table(summ, file.path(outdir, "agreement_with_bacdive.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Numeric attributes: Wikipedia's value against the mean of each source's values for the taxon
cat("\n=== Numeric attributes: Wikipedia against each source (value = midpoint of the range)\n")
num <- list()
for (phys in numeric) {
    w <- read_wiki(phys); if (is.null(w)) next
    v <- strsplit(w$Attribute_value, "(?<=[0-9])-", perl = TRUE)
    w$val <- vapply(v, \(z) mean(as.numeric(z)), 0)
    d <- ph[[phys]]; d <- d[!is.na(d$Attribute_source) & d$Attribute_source != "Bergey's Manual", ]; d$key <- key_of(d)
    if ("Attribute_value_min" %in% names(d)) d$val <- (as.numeric(d$Attribute_value_min) + as.numeric(d$Attribute_value_max)) / 2
    else d$val <- as.numeric(d$Attribute_value)
    r <- d |> filter(is.finite(val)) |> group_by(source = Attribute_source, key) |> summarise(ref = mean(val), .groups = "drop")
    j <- inner_join(select(w, key = NCBI_ID, rank, val), r, by = "key")
    tol <- c("growth temperature" = 5, "optimal ph" = 0.5, "gc content" = 2, "length" = 1, "width" = 0.25)[[phys]]
    cat(sprintf("%s: Wikipedia rows %d\n", phys, nrow(w)))
    if (!nrow(j)) next
    num[[phys]] <- j |> group_by(physiology = phys, reference = source) |>
        summarise(taxa = n(), median_abs_diff = round(median(abs(val - ref)), 2), within_tolerance = round(100 * mean(abs(val - ref) <= tol), 1),
                  tolerance = tol, spearman = if (n() >= 5) round(suppressWarnings(cor(val, ref, method = "spearman")), 2) else NA_real_, .groups = "drop")
}
num <- bind_rows(num); print(as.data.frame(num), row.names = FALSE)
write.table(num, file.path(outdir, "agreement_numeric.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink()
