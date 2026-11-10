## Run from the repo root: Rscript inst/scripts/test_filter_asr.R
suppressPackageStartupMessages(library(dplyr))
source("inst/scripts/filter_asr.R")
## one row per attribute, as export_bugphyzz.R builds it from
## bugphyzz's validation_attributes.tsv (rank "all", score_lcb)
validation <- data.frame(
    .attr_key = c("aerophilicity", "motility", "length", "gc content"),
    value = c(0.786, 0.586, 0.22, 0.942)
)
aero <- data.frame(
    Attribute = "aerophilicity", Attribute_type = "multistate-intersection",
    Attribute_value = c("aerobic", "aerobic", "anaerobic", "anaerobic",
                        "facultatively anaerobic", "aerobic"),
    Evidence = c("exp", "asr", "asr", "asr", "asr", NA),
    Frequency = c("always", "usually", "usually", "never", "sometimes", "always"),
    Score = c(1, 0.95, 0.9, 0, 0.6, 1),
    stringsAsFactors = FALSE
)
res <- filterAsr(aero, validation, 0.7, 0.5, 0.5)
## the attribute passes (0.786 >= 0.7): every imputed state is kept, whatever
## the state, except the row with Score 0 (below 0.5); curated rows and the
## NA-evidence row stay
stopifnot(
    nrow(res$kept) == 5, nrow(res$dropped) == 1,
    nrow(res$kept) + nrow(res$dropped) == nrow(aero),
    res$n_before == 4, res$n_after == 3,
    !"value" %in% names(res$kept), "Validation" %in% names(res$dropped),
    res$dropped$Attribute_value == "anaerobic", res$dropped$Score == 0,
    res$dropped$Validation == 0.786,
    "facultatively anaerobic" %in% res$kept$Attribute_value,
    !any(is.na(res$kept$Attribute_value))
)
## a higher Score threshold removes imputed rows within a passing attribute
res <- filterAsr(aero, validation, 0.7, 0.5, 0.8)
stopifnot(nrow(res$kept) == 4, identical(sort(res$dropped$Score), c(0, 0.6)))
## binary: the attribute fails (0.586 < 0.7), so every imputed row drops
mot <- data.frame(Attribute = "motility", Attribute_type = "binary",
    Attribute_value = c("TRUE", "FALSE", "TRUE"), Evidence = c("asr", "exp", "asr"),
    Frequency = "always", Score = c(0.99, 1, 0.7), stringsAsFactors = FALSE)
res <- filterAsr(mot, validation, 0.7, 0.5, 0.5)
stopifnot(nrow(res$kept) == 1, res$kept$Evidence == "exp",
          nrow(res$dropped) == 2, all(res$dropped$Validation == 0.586))
res <- filterAsr(mot, validation, 0.5, 0.5, 0.5)
stopifnot(nrow(res$kept) == 3, nrow(res$dropped) == 0)
## numeric: R2 threshold, "range" counts as numeric
len <- data.frame(Attribute = "length", Attribute_type = "range",
    Attribute_value = "1-2", Evidence = "asr", Frequency = "always", Score = 0.9,
    stringsAsFactors = FALSE)
stopifnot(nrow(filterAsr(len, validation, 0.7, 0.5, 0.5)$kept) == 0,
          nrow(filterAsr(len, validation, 0.7, 0.2, 0.5)$kept) == 1)
gc <- data.frame(Attribute = "gc content", Attribute_type = "numeric",
    Attribute_value = 50, Evidence = "asr", Frequency = "always", Score = 0.5,
    stringsAsFactors = FALSE)
stopifnot(nrow(filterAsr(gc, validation, 0.7, 0.5, 0.5)$kept) == 1)
## an attribute with no validation row loses its imputed rows
hab <- data.frame(Attribute = "habitat", Attribute_type = "multistate-union",
    Attribute_value = c("soil", "soil"), Evidence = c("asr", "exp"),
    Frequency = "always", Score = c(1, 1), stringsAsFactors = FALSE)
res <- filterAsr(hab, validation, 0.7, 0.5, 0.5)
stopifnot(nrow(res$kept) == 1, is.na(res$dropped$Validation))
## an imputed row with NA Score fails the Score rule and goes to dropped
na <- aero[2, ]; na$Score <- NA_real_
res <- filterAsr(rbind(aero, na), validation, 0.7, 0.5, 0.5)
stopifnot(nrow(res$dropped) == 2, any(is.na(res$dropped$Score)))
## empty input
e <- filterAsr(aero[0, ], validation, 0.7, 0.5, 0.5)
stopifnot(nrow(e$kept) == 0, nrow(e$dropped) == 0, e$n_before == 0)
## round trip through the gz side file as export_bugphyzz.R writes it
tmp <- tempfile(fileext = ".csv.gz")
d <- filterAsr(mot, validation, 0.7, 0.5, 0.5)$dropped
con <- gzfile(tmp, "w"); cat("# header\n", file = con)
suppressWarnings(write.table(d, con, quote = TRUE, sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE)); close(con)
back <- read.csv(gzfile(tmp), skip = 1, colClasses = "character")
stopifnot(nrow(back) == nrow(d), identical(names(back), names(d)))
cat("ok\n")
